//! The original's Python pipeline (`scripts/*_pipeline.py`), which its pipeline script
//! (`run_*.py`) runs for `--platform ilmn`, and for ONT and HiFi with `--disable_c_impl`.
//!
//! It differs from the C-implementation pipeline in three places: pileup tensors come from
//! `samtools mpileup` text (`CreateTensorPileup.py`), each contig's reads are haplotagged into
//! `tmp/phase_output/phase_bam/<contig>.bam`, and full-alignment tensors come from those reads
//! (`CreateTensorFullAlignment.py`), for Illumina after short-read realignment
//! (`RealignReads.py`). Illumina also merges by candidate region.

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use rayon::prelude::*;
use snv_call::vcf_header;
use snv_infer::{
    FullAlignmentModel, InferenceDevice, PileupModel, Precision, load_full_alignment_model, load_pileup_model,
};
use snv_merge::{
    MergeMode, MergeOptions, Postprocess, SortOptions, bgzip, compress_and_index, merge_vcf, sort_vcf, tabix_index,
};
use snv_mpileup::{FastMode, PileupJob};
use snv_phase::{PhaseJob, haplotag_contig, haplotag_output, phase_contig, phase_output};
use snv_pileup::{CandidateAlleles, ChunkRegion, Fasta, TENSOR_LEN, known_sites_chunk};
use snv_select::{
    CandidateOptions, CandidateVcf, CandidateWindow, FullQualOptions, HetSnpOptions, open_vcf,
    select_candidates_for_full_alignment, select_full_qual, select_het_snps_for_phasing, select_phase_qual,
};

use crate::chunk::{self, VcfOutput};
use crate::contigs::{ContigRequest, check_envs};
use crate::genotyping::add_back_missing_variants;
use crate::parallel::{self, Batching, MemoryTrimmer, StageTimes};
use crate::pipeline::{
    CPU_BATCH, FA_MODEL, PILEUP_MODEL, TRIM_INTERVAL, concurrent_phasing, ensure_merge_output, for_each_contig,
    link_unphased_bams, log_times, longphase, remove_intermediate,
};
use crate::settings::{Phaser, Platform, RunArgs, Settings};

pub(crate) fn run(args: &RunArgs) -> Result<()> {
    let started = Instant::now();
    let settings = Settings::resolve(args)?;
    if settings.device != InferenceDevice::Cpu {
        // The original's Python pipeline has no GPU path: its calling jobs run on the CPU.
        log::warn!("the Python pipeline calls on the CPU; ignoring the GPU device");
    }
    let platform = match args.platform {
        Platform::Ont => snv_realign::Platform::Ont,
        Platform::Hifi => snv_realign::Platform::Hifi,
        Platform::Ilmn => snv_realign::Platform::Ilmn,
    };
    let threads = args.threads.clamp(1, std::thread::available_parallelism().map_or(usize::MAX, |n| n.get()));
    let bed = crate::pipeline::check_inputs(args, &settings)?;
    crate::pipeline::clear_previous_run(&args.output)?;
    let tmp = args.output.join("tmp");
    let pileup_dir = tmp.join("pileup_output");
    let phase_dir = tmp.join("phase_output").join("phase_vcf");
    let phase_bam_dir = tmp.join("phase_output").join("phase_bam");
    let merge_dir = tmp.join("merge_output");
    let fa_dir = tmp.join("full_alignment_output");
    let bed_dir = fa_dir.join("candidate_bed");
    for dir in [&pileup_dir, &phase_dir, &phase_bam_dir, &merge_dir, &bed_dir, &tmp.join("gvcf_tmp_output")] {
        fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    }
    let gvcf = args.gvcf.steps(&tmp, bed.as_deref());
    let cmd_file = tmp.join("CMD");
    fs::write(&cmd_file, std::env::args().collect::<Vec<_>>().join(" "))?;
    let contigs = check_envs(
        &ContigRequest {
            bam: &args.bam,
            reference: &args.reference,
            bed: bed.as_deref(),
            vcf: args.vcf_fn.as_deref(),
            ctg_name: args.ctg_name.as_deref(),
            include_all_ctgs: args.include_all_ctgs,
            min_contig_size: settings.min_contig_size,
            chunk_size: args.chunk_size,
            chunk_num: args.chunk_num,
            gender: settings.gender,
        },
        &tmp,
    )?;
    let contigs_file = tmp.join("CONTIGS");
    let names = contigs.names();
    if names.is_empty() {
        log::warn!("No contig intersection found, output header only in merge_output.vcf.gz");
        let header = vcf_header(Some(&args.reference), Some(&cmd_file), &args.sample_name)?;
        let merged = args.output.join("merge_output.vcf");
        fs::write(&merged, header.trim_end_matches('\n'))?;
        compress_and_index(&merged)?;
        return Ok(());
    }

    // Stage 1: one `CallVarBam --pileup` job per chunk: `CreateTensorPileup` tensors, scored on
    // one thread in batches of 200.
    let jobs: Vec<(&str, i64, i64)> =
        contigs.chunks.iter().flat_map(|(c, n)| (1..=*n).map(move |id| (c.as_str(), id, *n))).collect();
    log::info!("[1/7] pileup calling: {} chunks on {threads} threads", jobs.len());
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build()?;
    let output = VcfOutput {
        reference: &args.reference,
        cmd_file: Some(&cmd_file),
        sample_name: &args.sample_name,
        // The original's pipeline does not pass --qual to its calling jobs.
        qual: 2,
        iupac: settings.iupac,
        max_indel_length: settings.candidates.max_indel_length,
        likelihoods: args.gvcf.likelihoods(),
    };
    let (min_snp_af, min_indel_af) = settings.min_af;
    let pileup_job = |&(contig, id, count): &(&str, i64, i64)| -> Result<Vec<snv_pileup::Candidate>> {
        let mut bam = snv_mpileup::AlignmentFile::open(&args.bam, &args.reference)?;
        let fasta = Fasta::open(&args.reference)?;
        let extend_bed = contigs.extend_bed(contig);
        let chunk = match (&contigs.known_sites, &extend_bed) {
            (Some(known), _) => {
                let sites = known.get(contig).map_or(&[][..], Vec::as_slice);
                ChunkRegion::KnownSites { sites: known_sites_chunk(sites, id - 1, count).to_vec() }
            }
            (None, Some(_)) => ChunkRegion::BedChunk { id: id - 1, count },
            (None, None) => ChunkRegion::FaiChunk { id: id - 1, count },
        };
        let job = PileupJob {
            contig,
            chunk,
            extend_bed: extend_bed.as_deref(),
            confident_bed: args.bed_fn.as_deref(),
            min_snp_af,
            min_indel_af,
            min_coverage: settings.min_coverage,
            min_mapq: settings.candidates.min_mapq,
            alleles: if args.call_snp_only { CandidateAlleles::SnpsOnly } else { CandidateAlleles::SnpsAndIndels },
            flanks: settings.candidates.flanks,
            fast_mode: if args.fast_mode && args.platform == Platform::Ont { FastMode::On } else { FastMode::Off },
        };
        let context = || format!("pileup chunk {contig} {id}/{count}");
        let Some(gvcf) = &gvcf else {
            return snv_mpileup::pileup_candidates(&mut bam, &fasta, &job).with_context(context);
        };
        // `--gvcf` runs `samtools mpileup -a`, which also changes the candidates.
        let (candidates, chunk) =
            snv_mpileup::pileup_candidates_with_sites(&mut bam, &fasta, &job).with_context(context)?;
        if let Some(chunk) = chunk {
            gvcf.write_mpileup_chunk(&args.reference, &args.sample_name, &chunk).with_context(context)?;
        }
        Ok(candidates)
    };
    let pileup_model = &settings.pileup_model;
    let times = StageTimes::default();
    parallel::score_in_batches(
        &pool,
        &jobs,
        Batching { tensor_len: TENSOR_LEN, batch: CPU_BATCH },
        &times,
        |job| {
            let candidates = pileup_job(job)?;
            let inputs = chunk::model_inputs(&candidates)?;
            Ok((candidates, inputs))
        },
        |batch| {
            // The 8-bit pileup network's QUALs need the main pipeline's rescoring near cut-offs.
            let config = settings.model_config(InferenceDevice::Cpu, Precision::Float32);
            PILEUP_MODEL.with(pileup_model, |p| Ok(load_pileup_model(p, config)?), |m| Ok(m.predict(batch)?))
        },
        |(contig, id, _), candidates, probs| {
            let path = pileup_dir.join(format!("pileup_{contig}_{id}.vcf"));
            // Like `CallVariants`, this removes a chunk VCF without records.
            output.write(&path, &candidates, &probs)
        },
    )?;
    log_times("pileup", &times);
    log::info!("pileup calling done in {:.1?}", started.elapsed());

    let qual = i32::try_from(args.qual).context("--qual out of range")?;
    let sorted = |input_dir: &Path, prefix: &str, output: PathBuf| {
        let mut sort = SortOptions::new(input_dir, output);
        sort.prefix = Some(prefix.to_owned());
        sort.sample_name.clone_from(&args.sample_name);
        sort.reference = Some(args.reference.clone());
        sort.contigs = Some(contigs_file.clone());
        sort.cmd_file = Some(cmd_file.clone());
        sort
    };
    let mut sort = sorted(&pileup_dir, "pileup", args.output.join("pileup.vcf"));
    if args.pileup_only {
        sort.postprocess = Postprocess::PileupOnly;
    }
    sort.ref_calls = settings.ref_calls;
    sort.haploid = settings.haploid;
    sort.gender = settings.gender;
    sort.par_regions.clone_from(&settings.par_regions);
    sort.qual = qual;
    let summary = sort_vcf(&sort)?;
    log::info!("merged {} pileup records into {}", summary.records, summary.vcf.display());
    if summary.records == 0 {
        log::info!("Exit in pileup variant calling");
        return ensure_merge_output(&args.output);
    }
    if args.pileup_only {
        remove_intermediate(args, &tmp)?;
        log::info!("Only call pileup output with --pileup_only, output file: {}", summary.vcf.display());
        return Ok(());
    }

    let _trimmer = MemoryTrimmer::start(TRIM_INTERVAL);
    if args.no_phasing_for_fa {
        log::info!("[2/7] No phasing for full alignment calling");
        link_unphased_bams(args, &phase_bam_dir, &names)?;
    } else {
        // Stage 2: het SNPs above the phasing quality cut-off, per contig.
        log::info!("[2/7] Select heterozygous SNP variants for Whatshap phasing and haplotagging");
        let phase_qual = select_phase_qual(open_vcf(&summary.vcf)?, settings.var_pct_phasing)?;
        for warning in phase_qual.warnings() {
            log::warn!("{warning}");
        }
        log::info!("{}", phase_qual.info_message());
        phase_qual.write_in(&phase_dir)?;
        let selected = pool.install(|| {
            names
                .par_iter()
                .map(|contig| {
                    let options = HetSnpOptions { contig: Some((*contig).to_owned()), ..HetSnpOptions::default() };
                    Ok(select_het_snps_for_phasing(&summary.vcf, &phase_dir, &options)?.info_message())
                })
                .collect::<Result<Vec<String>>>()
        })?;
        for message in selected {
            log::info!("{message}");
        }

        // Stages 3 and 4: phase each contig's het SNPs, then tag its reads by haplotype.
        let phaser = &settings.intermediate_phasing;
        let tool = if *phaser == Phaser::WhatsHap { "lrpipe-phase" } else { "LongPhase" };
        log::info!("[3/7] Phase VCF file using {tool}, and [4/7] haplotag input BAM file");
        for_each_contig(&contigs, threads, concurrent_phasing(phaser, threads), |contig, per_job| {
            let phased = phase_dir.join(format!("phased_{contig}.vcf.gz"));
            let variants = phase_dir.join(format!("{contig}.vcf"));
            match phaser {
                Phaser::WhatsHap => phase_contig(&PhaseJob {
                    variants: &variants,
                    alignments: &args.bam,
                    reference: &args.reference,
                    contig,
                    output: &phased,
                    threads: per_job,
                })?,
                Phaser::LongPhase(bin) => {
                    longphase(bin, &variants, args, per_job, &phase_dir.join(format!("phased_{contig}")))?;
                    bgzip(&phase_dir.join(format!("phased_{contig}.vcf")))?;
                }
            }
            // `tabix -f -p vcf phased_<contig>.vcf.gz`, which `whatshap haplotag` needs.
            if !tabix_index(&phased)? {
                anyhow::bail!("tabix could not index {}", phased.display());
            }
            let bam = phase_bam_dir.join(format!("{contig}.bam"));
            haplotag_contig(&phased, &args.bam, &args.reference, contig, &bam, per_job)?;
            Ok(())
        })?;
        log::info!("phasing and haplotagging done in {:.1?}", started.elapsed());
    }

    // Stage 5: full-alignment candidates, the low-quality pileup calls.
    log::info!("[5/7] Select candidates for full-alignment calling");
    let full_options = FullQualOptions {
        var_pct_full: settings.var_pct_full,
        ref_pct_full: settings.ref_pct_full,
        candidate_vcf: if args.vcf_fn.is_some() { CandidateVcf::Present } else { CandidateVcf::Absent },
    };
    let full_qual = select_full_qual(open_vcf(&summary.vcf)?, &full_options)?;
    for warning in full_qual.warnings() {
        log::warn!("{warning}");
    }
    for message in full_qual.info_messages() {
        log::info!("{message}");
    }
    full_qual.write_in(&bed_dir)?;
    let per_contig = pool.install(|| {
        names
            .par_iter()
            .map(|contig| {
                let options = CandidateOptions {
                    contig: Some((*contig).to_owned()),
                    var_pct_full: settings.var_pct_full,
                    ref_pct_full: settings.ref_pct_full,
                    window: CandidateWindow::from_platform(args.platform.as_str()),
                    ..CandidateOptions::default()
                };
                let beds = select_candidates_for_full_alignment(&summary.vcf, &bed_dir, &options)?;
                Ok((beds.messages().to_vec(), beds.full_aln_file().map(Path::to_owned)))
            })
            .collect::<Result<Vec<_>>>()
    })?;
    let mut lists: Vec<PathBuf> = Vec::new();
    for (messages, list) in per_contig {
        for message in messages {
            log::info!("{message}");
        }
        lists.extend(list);
    }
    // `sorted(candidate_bed_path.glob("FULL_ALN_FILE_*"))`, concatenated into FULL_ALN_FILES.
    lists.sort();
    if lists.is_empty() {
        log::info!("No Candidate found! Exit in selecting full-alignment candidates");
        return ensure_merge_output(&args.output);
    }
    let mut all_lists = String::new();
    for list in &lists {
        let text = fs::read_to_string(list).with_context(|| format!("reading {}", list.display()))?;
        for line in text.lines() {
            all_lists.push_str(line);
            all_lists.push('\n');
        }
    }
    fs::write(bed_dir.join("FULL_ALN_FILES"), &all_lists)?;
    let fa_jobs: Vec<PathBuf> = all_lists.lines().filter(|l| !l.is_empty()).map(PathBuf::from).collect();

    // Stage 6: one `CallVarBam` job per candidate BED: the reads of the haplotagged BAM (for
    // Illumina realigned), `CreateTensorFullAlignment` tensors, scored on one thread in
    // batches of 200.
    log::info!("[6/7] Call low-quality variants using full-alignment model: {} jobs", fa_jobs.len());
    let fa_model = &settings.fa_model;
    let min_mapq = settings.candidates.min_mapq;
    let fa_depth = platform.matrix_depth();
    let times = StageTimes::default();
    parallel::score_in_batches(
        &pool,
        &fa_jobs,
        Batching {
            tensor_len: snv_infer::full_alignment_tensor_len(fa_depth, snv_infer::FULL_ALIGNMENT_CHANNELS),
            batch: CPU_BATCH,
        },
        &times,
        |bed| {
            let name = bed.file_name().and_then(|n| n.to_str()).context("candidate BED without a name")?;
            // GNU parallel's `{1/.}`: the BED name without its last extension.
            let contig = Path::new(name).file_stem().and_then(|n| n.to_str()).unwrap_or(name);
            let bam = phase_bam_dir.join(format!("{contig}.bam"));
            let job = snv_realign::FullAlignmentJob {
                platform,
                bam: &bam,
                reference: &args.reference,
                contig,
                regions: bed,
                min_mapq,
                min_coverage: settings.min_coverage,
                known_sites: contigs.known_sites.as_ref().and_then(|k| k.get(contig)).map(Vec::as_slice),
            };
            let candidates =
                snv_realign::full_alignment_candidates(&job).with_context(|| format!("full-alignment job {name}"))?;
            let inputs = snv_realign::padded_inputs(&candidates, fa_depth);
            let (positions, alt_info) = candidates.into_iter().map(|c| (c.position, c.alt_info)).unzip();
            Ok(((name.to_owned(), positions, alt_info), inputs))
        },
        |batch| {
            FA_MODEL.with(
                fa_model,
                |p| {
                    let config = settings.model_config(InferenceDevice::Cpu, settings.precision);
                    Ok(load_full_alignment_model(p, config, fa_depth)?)
                },
                |m| Ok(m.predict(batch)?),
            )
        },
        |_, (name, positions, alt_info): (String, Vec<String>, Vec<String>), probs| {
            let path = fa_dir.join(format!("full_alignment_{name}.vcf"));
            output.write_full_alignment(&path, &positions, &alt_info, &probs)
        },
    )?;
    log_times("full-alignment", &times);
    log::info!("full-alignment calling done in {:.1?}", started.elapsed());

    // The original's Python pipeline does not pass --cmd_fn when sorting the full-alignment calls.
    let mut sort = sorted(&fa_dir, "full_alignment", args.output.join("full_alignment.vcf"));
    sort.cmd_file = None;
    let fa_summary = sort_vcf(&sort)?;
    log::info!("merged {} full-alignment records into {}", fa_summary.records, fa_summary.vcf.display());
    if fa_summary.records == 0 {
        log::info!("Exit in full-alignment variant calling");
        return ensure_merge_output(&args.output);
    }
    if let Some(gvcf) = &gvcf {
        gvcf.sort_chunks(&args.reference, &contigs_file, &cmd_file)?;
    }

    // Stage 7: per contig, `MergeVcf`; for Illumina `MergeVcf_illumina`, where full-alignment
    // calls replace the pileup calls inside the candidate regions.
    log::info!("[7/7] Merge pileup VCF and full-alignment VCF");
    pool.install(|| {
        names.par_iter().try_for_each(|contig| -> Result<()> {
            let merged = merge_dir.join(format!("merge_{contig}.vcf"));
            let mut merge = MergeOptions::new(&summary.vcf, &fa_summary.vcf, &merged);
            merge.contig = Some((*contig).to_owned());
            if let Some(gvcf) = &gvcf {
                merge.gvcf = gvcf.merge_output(&merged, &args.reference);
            }
            if args.platform == Platform::Ilmn {
                merge.mode = MergeMode::Regions { candidate_beds: bed_dir.clone() };
            }
            merge.ref_calls = settings.ref_calls;
            merge.haploid = settings.haploid;
            merge.gender = settings.gender;
            merge.par_regions.clone_from(&settings.par_regions);
            merge.qual = qual;
            merge_vcf(&merge)?;
            Ok(())
        })
    })?;
    let final_sort = sorted(&merge_dir, "merge", args.output.join("merge_output.vcf"));
    let final_summary = sort_vcf(&final_sort)?;
    log::info!("{} calls in {} ({:.1?})", final_summary.records, final_summary.vcf.display(), started.elapsed());
    if final_summary.records == 0 {
        log::info!("Exit in variant merging");
        return Ok(());
    }
    if let Some(gvcf) = &gvcf {
        gvcf.sort_merged(&merge_dir, &args.output, &final_sort)?;
    }

    if let Some(phaser) = &settings.final_phasing {
        let tool = if *phaser == Phaser::WhatsHap { "lrpipe-phase" } else { "LongPhase" };
        log::info!("[7/7] Phasing VCF output in parallel using {tool}");
        for_each_contig(&contigs, threads, concurrent_phasing(phaser, threads), |contig, per_job| {
            let variants = merge_dir.join(format!("merge_{contig}.vcf"));
            match phaser {
                Phaser::WhatsHap => Ok(phase_output(&PhaseJob {
                    variants: &variants,
                    alignments: &args.bam,
                    reference: &args.reference,
                    contig,
                    output: &merge_dir.join(format!("phased_merge_{contig}.vcf")),
                    threads: per_job,
                })?),
                Phaser::LongPhase(bin) => {
                    longphase(bin, &variants, args, per_job, &merge_dir.join(format!("phased_merge_{contig}")))
                }
            }
        })?;
        let phased = sort_vcf(&sorted(&merge_dir, "phased_merge", args.output.join("phased_merge_output.vcf")))?;
        log::info!("{} phased calls in {}", phased.records, phased.vcf.display());
        if settings.haplotag {
            log::info!("[4/7] Haplotag input BAM file using lrpipe-phase");
            let bam = args.output.join("phased_output.bam");
            haplotag_output(&phased.vcf, &args.bam, &args.reference, &bam, threads)?;
            log::info!("Finish calling, phased output BAM file: {}", bam.display());
        }
    }
    if let Some(vcf) = &args.vcf_fn {
        log::info!("Double check re-genotyping variants");
        let added = add_back_missing_variants(vcf, &final_summary.vcf, &args.output.join("merge_output.vcf"))?;
        log::info!("{added}");
    }
    remove_intermediate(args, &tmp)?;
    log::info!("Finish calling, output file: {}", final_summary.vcf.display());
    Ok(())
}

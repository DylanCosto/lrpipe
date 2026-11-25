//! `lrpipe-snv run`: the original's pipeline (`scripts/*_c_impl_pipeline.py`) from the BAM to
//! `merge_output.vcf.gz`, writing its output layout so each step can be compared with it.

use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, Read};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, mpsc};
use std::time::{Duration, Instant};

use anyhow::{Context, Result, anyhow, bail};
use pasta::bgzf::{self, VirtualPosition};
use rayon::prelude::*;
use snv_call::{Haploid, vcf_header};
use snv_fullalign::{FeatureOptions, FullAlignment, FullAlignmentJob};
use snv_infer::{
    FullAlignmentModel, InferenceDevice, ModelConfig, PileupModel, Precision, full_alignment_tensor_len,
    load_full_alignment_model, load_pileup_model,
};
use snv_merge::{
    ChunkLayout, Gender, MergeMode, MergeOptions, Postprocess, SortOptions, bgzip, compress_and_index, merge_vcf,
    sort_vcf,
};
use snv_phase::{PhaseJob, haplotag_output, phase_contig, phase_output};
use snv_pileup::{
    AlignmentFile, Candidate, ChunkRegion, Fasta, TENSOR_LEN, candidate_tensors, candidate_tensors_with_sites,
    known_sites_chunk,
};
use snv_select::{
    CandidateOptions, CandidateVcf, CandidateWindow, FullQualOptions, HetSnpOptions, open_vcf,
    select_candidates_for_full_alignment_from, select_full_qual, select_het_snps_for_phasing_from, select_phase_qual,
};

use lrpipe_reads::progress::Progress;

use crate::chunk::{self, VcfOutput};
use crate::contigs::{ContigRequest, Contigs, check_envs};
use crate::genotyping::add_back_missing_variants;
use crate::gvcf;
use crate::parallel::{self, Batching, MemoryTrimmer, StageTimes, ThreadModel};
use crate::rescore::{self, PackedInputs, ScoredChunk};
use crate::settings::{GpuWait, Phaser, Platform, ReadOrder, RunArgs, Settings};

/// Candidates per batch in the original's GPU calling (`param.predictBatchSize * 5`).
const GPU_BATCH: usize = 1000;
/// Candidates per batch in the original's CPU jobs (`param.predictBatchSize * threads`, one
/// thread).
pub(crate) const CPU_BATCH: usize = 200;
/// How often freed heap memory goes back to the system from phasing on ([`MemoryTrimmer`]).
pub(crate) const TRIM_INTERVAL: Duration = Duration::from_secs(5);

/// Builds model inputs for each `(index, job)` on the worker pool and scores them on one GPU
/// thread as they arrive; results come back in index order. `jobs` is drawn on the calling
/// thread, so it may block until more jobs are ready.
fn score_on_gpu<'a, J: Sync + 'a, P: Send, R: Send>(
    pool: &rayon::ThreadPool,
    jobs: impl IntoIterator<Item = (usize, &'a J)>,
    produce: impl Fn(&J) -> Result<P> + Sync,
    mut score: impl FnMut(P) -> Result<R> + Send,
) -> Result<Vec<R>> {
    // A short queue keeps only a few jobs' tensors in memory while the GPU catches up.
    let (tx, rx) = mpsc::sync_channel::<(usize, P)>(2);
    let failure: Mutex<Option<anyhow::Error>> = Mutex::new(None);
    std::thread::scope(|scope| {
        let gpu = scope.spawn(move || -> Result<Vec<(usize, R)>> {
            rx.into_iter().map(|(i, input)| Ok((i, score(input)?))).collect()
        });
        let (produce, failure) = (&produce, &failure);
        pool.in_place_scope(|s| {
            for (i, job) in jobs {
                let tx = tx.clone();
                s.spawn(move |_| {
                    if failure.lock().expect("lock poisoned").is_some() {
                        return;
                    }
                    let sent =
                        produce(job).and_then(|input| tx.send((i, input)).map_err(|_| anyhow!("GPU scoring stopped")));
                    if let Err(e) = sent {
                        failure.lock().expect("lock poisoned").get_or_insert(e);
                    }
                });
            }
        });
        drop(tx);
        let scored = gpu.join().map_err(|_| anyhow!("GPU thread panicked"))?;
        // A GPU error closes the queue, so report it rather than the send failure it causes.
        let mut scored = scored?;
        if let Some(e) = failure.lock().expect("lock poisoned").take() {
            return Err(e);
        }
        scored.sort_by_key(|(i, _)| *i);
        Ok(scored.into_iter().map(|(_, r)| r).collect())
    })
}

/// A full-alignment job's candidates' positions, alt info and probabilities.
type ScoredJob = (Vec<String>, Vec<String>, Vec<f32>);

/// What a full-alignment producer hands the GPU: a batch of a job's matrices, or the end of
/// the job with its candidates' positions and alt info.
enum FaMessage {
    Batch(usize, Vec<i8>),
    Done(usize, Vec<String>, Vec<String>),
}

/// Full-alignment jobs on the GPU: `produce` builds a job's matrices on `pool` and hands them
/// to its sink a batch at a time, and the GPU scores each batch as it arrives, so a job's
/// matrices are never all in memory. Each job's positions, alt info and probabilities go to
/// `done` as the job finishes, on the GPU thread.
fn score_full_alignment_on_gpu<'a, J: Sync + 'a>(
    pool: &rayon::ThreadPool,
    jobs: impl IntoIterator<Item = (usize, &'a J)>,
    produce: impl Fn(&J, &mut dyn FnMut(&[i8])) -> Result<(Vec<String>, Vec<String>)> + Sync,
    mut score: impl FnMut(&[i8]) -> Result<Vec<f32>> + Send,
    mut done: impl FnMut(usize, ScoredJob) -> Result<()> + Send,
) -> Result<()> {
    // A short queue keeps a few batches in memory while the GPU catches up.
    let (tx, rx) = mpsc::sync_channel::<FaMessage>(4);
    let failure: Mutex<Option<anyhow::Error>> = Mutex::new(None);
    std::thread::scope(|scope| {
        let gpu = scope.spawn(move || -> Result<()> {
            let mut probs: HashMap<usize, Vec<f32>> = HashMap::new();
            for message in rx {
                match message {
                    FaMessage::Batch(i, tensors) => probs.entry(i).or_default().extend(score(&tensors)?),
                    FaMessage::Done(i, positions, alt_info) => {
                        done(i, (positions, alt_info, probs.remove(&i).unwrap_or_default()))?;
                    }
                }
            }
            Ok(())
        });
        let (produce, failure) = (&produce, &failure);
        pool.in_place_scope(|s| {
            for (i, job) in jobs {
                let tx = tx.clone();
                s.spawn(move |_| {
                    if failure.lock().expect("lock poisoned").is_some() {
                        return;
                    }
                    let mut stopped = false;
                    let mut sink = |batch: &[i8]| {
                        stopped |= tx.send(FaMessage::Batch(i, batch.to_vec())).is_err();
                    };
                    let sent = produce(job, &mut sink).and_then(|(positions, alt_info)| {
                        if stopped {
                            bail!("GPU scoring stopped");
                        }
                        tx.send(FaMessage::Done(i, positions, alt_info)).map_err(|_| anyhow!("GPU scoring stopped"))
                    });
                    if let Err(e) = sent {
                        failure.lock().expect("lock poisoned").get_or_insert(e);
                    }
                });
            }
        });
        drop(tx);
        let scored = gpu.join().map_err(|_| anyhow!("GPU thread panicked"))?;
        // A GPU error closes the queue, so report it rather than the send failure it causes.
        scored?;
        if let Some(e) = failure.lock().expect("lock poisoned").take() {
            return Err(e);
        }
        Ok(())
    })
}

pub(crate) fn run(args: &RunArgs) -> Result<()> {
    if args.platform == Platform::Ilmn {
        log::warn!("Illumina platform will disable C implement to support short read realignment process!");
        return crate::python::run(args);
    }
    if args.disable_c_impl {
        if !(args.enable_dwell_time && args.platform == Platform::Ont) {
            return crate::python::run(args);
        }
        log::warn!("Dwell time requires the C implementation; re-enabling C implementation workflow.");
    }
    let started = Instant::now();
    let settings = Settings::resolve(args)?;
    #[cfg(feature = "gpu")]
    if let (InferenceDevice::Cuda(_), Some(nvrtc)) = (settings.device, &settings.nvrtc) {
        snv_infer::use_nvrtc(nvrtc)?;
    }
    if let (InferenceDevice::Cuda(i), GpuWait::Block) = (settings.device, settings.gpu_wait)
        && !snv_infer::block_while_waiting(i)
    {
        log::warn!("--gpu-wait block: could not set the CUDA runtime's blocking wait; waiting on a core");
    }
    let threads = args.threads.clamp(1, std::thread::available_parallelism().map_or(usize::MAX, |n| n.get()));
    let platform = args.platform.as_str();
    let canonical_read_order = args.read_order == ReadOrder::Canonical;
    let bed = check_inputs(args, &settings)?;
    clear_previous_run(&args.output)?;
    let tmp = args.output.join("tmp");
    let pileup_dir = tmp.join("pileup_output");
    let phase_dir = tmp.join("phase_output").join("phase_vcf");
    let merge_dir = tmp.join("merge_output");
    for dir in [&pileup_dir, &phase_dir, &merge_dir] {
        fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    }
    let cmd_file = tmp.join("CMD");
    fs::write(&cmd_file, std::env::args().collect::<Vec<_>>().join(" "))?;
    gvcf::create_tmp_dir(&tmp)?;
    let gvcf = args.gvcf.steps(&tmp, bed.as_deref());
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
        fs::write(&merged, header)?;
        compress_and_index(&merged)?;
        if gvcf.is_some() {
            let empty = args.output.join("merge_output.gvcf");
            fs::write(&empty, snv_call::gvcf_header(&args.reference, &args.sample_name)?.0)?;
            compress_and_index(&empty)?;
        }
        finish_without_calls(&args.output, settings.final_phasing.is_some())?;
        return remove_intermediate(args, &tmp);
    }

    // Stage 1: one pileup job per chunk, each on one thread with the original's batch
    // size, so every chunk VCF matches the corresponding `CallVariantsFromCffi --pileup` job.
    let jobs: Vec<(&str, i64, i64)> =
        contigs.chunks.iter().flat_map(|(c, n)| (1..=*n).map(move |id| (c.as_str(), id, *n))).collect();
    log::info!("[1/7] pileup calling: {} chunks on {threads} threads", jobs.len());
    // How far the run has got, for `lrpipe run`'s display: each phase's share of a whole genome's
    // time (HG002 ONT: pileup 38%, full-alignment 49%).
    let progress = Progress::from_env();
    let pileup_done = AtomicUsize::new(0);
    let pileup_jobs = jobs.len();
    progress.mark(0.0, &format!("pileup calling: 0 of {pileup_jobs} chunks"));
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
    let region_candidates =
        |contig: &str, region: &ChunkRegion, extend_bed: Option<&Path>, gvcf: Option<&gvcf::Gvcf>| -> Result<_> {
            let mut bam =
                AlignmentFile::open(&args.bam, &args.reference)?.with_canonical_read_order(canonical_read_order);
            let fasta = Fasta::open(&args.reference)?;
            let params = &settings.candidates;
            let Some(gvcf) = gvcf else {
                return Ok(candidate_tensors(&mut bam, &fasta, contig, region, extend_bed, params)?);
            };
            let (candidates, sites) =
                candidate_tensors_with_sites(&mut bam, &fasta, contig, region, extend_bed, params)?;
            if let Some(sites) = sites {
                gvcf.write_chunk(&args.reference, &args.sample_name, &sites)?;
            }
            Ok(candidates)
        };
    let call_region = |contig: &str, region: &ChunkRegion, extend_bed: Option<&Path>| {
        region_candidates(contig, region, extend_bed, gvcf.as_ref())
    };
    let pileup_candidates = |&(contig, id, count): &(&str, i64, i64)| -> Result<Vec<Candidate>> {
        let extend_bed = contigs.extend_bed(contig);
        let region = match (&contigs.known_sites, &extend_bed) {
            (Some(known), _) => {
                let sites = known.get(contig).map_or(&[][..], Vec::as_slice);
                ChunkRegion::KnownSites { sites: known_sites_chunk(sites, id - 1, count).to_vec() }
            }
            (None, Some(_)) => ChunkRegion::BedChunk { id: id - 1, count },
            (None, None) => ChunkRegion::FaiChunk { id: id - 1, count },
        };
        let candidates = call_region(contig, &region, extend_bed.as_deref())
            .with_context(|| format!("pileup chunk {contig} {id}/{count}"))?;
        let done = pileup_done.fetch_add(1, Ordering::Relaxed) + 1;
        let share = done as f64 / pileup_jobs as f64;
        let what = format!("pileup calling: {done} of {pileup_jobs} chunks");
        if done == pileup_jobs {
            progress.mark(PILEUP_SHARE, &what)
        } else {
            progress.report(PILEUP_SHARE * share, &what)
        }
        Ok(candidates)
    };
    let pileup_model = &settings.pileup_model;
    // `--approximate` on the CPU uses the 8-bit pileup network only if it keeps this run's own
    // candidates within what the rescoring catches.
    let mut pileup_precision = settings.precision;
    if settings.device == InferenceDevice::Cpu && settings.precision == Precision::Mixed {
        let checked = if contigs.known_sites.is_some() {
            Err("the run calls known sites only".to_owned())
        } else {
            // The sampled intervals overlap the calling chunks. Only calling writes gVCF blocks.
            let probe = |contig: &str, region: &ChunkRegion, extend_bed: Option<&Path>| {
                region_candidates(contig, region, extend_bed, None)
            };
            check_eight_bit_pileup(&pool, &jobs, &contigs, &probe, &output, &settings)?
        };
        if let Err(reason) = checked {
            log::info!("pileup calling in float32: the 8-bit pileup network is not used for this run, as {reason}");
            pileup_precision = Precision::Float32;
        }
    }
    let cpu_pileup = settings.model_config(InferenceDevice::Cpu, pileup_precision);
    match settings.device {
        InferenceDevice::Cpu => {
            let times = StageTimes::default();
            // With `--approximate` the chunks are kept, their inputs packed, until the candidates
            // near a cut-off have been scored again in `float32`; then their VCFs are written.
            let approximate = pileup_precision == Precision::Mixed;
            let scored: Mutex<Vec<ScoredChunk>> = Mutex::new(Vec::new());
            parallel::score_in_batches(
                &pool,
                &jobs,
                Batching { tensor_len: TENSOR_LEN, batch: CPU_BATCH },
                &times,
                |job| {
                    let mut candidates = pileup_candidates(job)?;
                    let inputs = chunk::model_inputs(&candidates)?;
                    let packed = if approximate {
                        for c in &mut candidates {
                            c.tensor = Vec::new();
                        }
                        inputs.chunks(TENSOR_LEN).map(PackedInputs::new).collect()
                    } else {
                        Vec::new()
                    };
                    Ok(((candidates, packed), inputs))
                },
                |batch| {
                    PILEUP_MODEL.with(
                        pileup_model,
                        |p| Ok(load_pileup_model(p, cpu_pileup)?),
                        |m| Ok(m.predict(batch)?),
                    )
                },
                |&(contig, id, _), (candidates, packed), probs| {
                    if approximate {
                        let chunk = ScoredChunk { contig: contig.to_owned(), id, candidates, probs, packed };
                        scored.lock().expect("lock poisoned").push(chunk);
                        return Ok(());
                    }
                    output.write(&pileup_dir.join(format!("pileup_{contig}_{id}.vcf")), &candidates, &probs)
                },
            )?;
            log_times("pileup", &times);
            if approximate {
                let mut scored = scored.into_inner().expect("lock poisoned");
                let exact = settings.model_config(InferenceDevice::Cpu, Precision::Float32);
                let load = || Ok(load_pileup_model(pileup_model, exact)?);
                let started = Instant::now();
                let n = rescore::rescore_near_cutoffs(
                    &pool,
                    &mut scored,
                    &output,
                    &settings,
                    args.qual as f64,
                    load,
                    CPU_BATCH,
                )?;
                log::info!(
                    "rescored {n} pileup candidates near a quality cut-off in float32 in {:.1?}",
                    started.elapsed()
                );
                pool.install(|| {
                    scored.par_iter().try_for_each(|c| {
                        let path = pileup_dir.join(format!("pileup_{}_{}.vcf", c.contig, c.id));
                        output.write(&path, &c.candidates, &c.probs)
                    })
                })?;
            }
        }
        device => {
            let mut model = load_pileup_model(pileup_model, settings.model_config(device, settings.precision))?;
            let calls = score_on_gpu(&pool, jobs.iter().enumerate(), pileup_candidates, |mut candidates| {
                let inputs = gpu_pileup_inputs(&candidates)?;
                let probs = chunk::score(&mut model, &inputs, GPU_BATCH)?;
                // Only the strings are needed from here on; a genome's tensors would not fit in
                // memory.
                for c in &mut candidates {
                    c.tensor = Vec::new();
                }
                Ok((candidates, probs))
            })?;
            let jobs: Vec<_> = calls.iter().map(|(c, p)| (c.as_slice(), p.as_slice())).collect();
            pool.install(|| output.write_jobs(&pileup_dir.join("pileup_0.vcf"), &jobs))?;
        }
    }
    if !contigs.unchunked.is_empty() {
        unchunked_pileup(&pool, &contigs.unchunked, &call_region, &output, (&settings, cpu_pileup), &pileup_dir)?;
    }
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
    sort.layout = chunk_layout(settings.device);
    let summary = sort_vcf(&sort)?;
    log::info!("merged {} pileup records into {}", summary.records, summary.vcf.display());
    if summary.records == 0 {
        log::info!("Exit in pileup variant calling");
        if let Some(gvcf) = &gvcf {
            gvcf.sort_chunks(&args.reference, &contigs_file, &cmd_file)?;
            for contig in &names {
                let merged = merge_dir.join(format!("merge_{contig}.vcf"));
                let mut merge = MergeOptions::new(&summary.vcf, &summary.vcf, &merged);
                merge.contig = Some((*contig).to_owned());
                merge.gvcf = gvcf.merge_output(&merged, &args.reference);
                merge.haploid = settings.haploid;
                merge.gender = settings.gender;
                merge.par_regions.clone_from(&settings.par_regions);
                merge_vcf(&merge)?;
            }
            // A fresh sort: the pileup's post-processes rows (`--pileup_only`) as calls, which
            // would mark every reference block LowQual and drop haploid ones.
            gvcf.sort_merged(
                &merge_dir,
                &args.output,
                &sorted(&merge_dir, "merge", args.output.join("merge_output.vcf")),
            )?;
        }
        finish_without_calls(&args.output, settings.final_phasing.is_some())?;
        return remove_intermediate(args, &tmp);
    }
    if args.pileup_only {
        remove_intermediate(args, &tmp)?;
        log::info!("Only call pileup output with --pileup_only, output file: {}", summary.vcf.display());
        return Ok(());
    }

    // Stages 2 and 5 select from the pileup calls: find each contig's records once, and each
    // contig's selection decompresses only those.
    let pileup_vcf = PileupVcf::index(&summary.vcf)?;
    progress.mark(PILEUP_SHARE, "selecting heterozygous SNPs");
    if args.no_phasing_for_fa {
        log::info!("[2/7] No phasing for full alignment calling");
        link_unphased_bams(args, &tmp.join("phase_output").join("phase_bam"), &names)?;
    } else {
        // Stage 2: het SNPs above the phasing quality cut-off, per contig.
        log::info!("[2/7] selecting heterozygous SNPs for phasing");
        let phase_qual = select_phase_qual(pileup_vcf.all()?, settings.var_pct_phasing)?;
        for warning in phase_qual.warnings() {
            log::warn!("{warning}");
        }
        log::info!("{}", phase_qual.info_message());
        phase_qual.write_in(&phase_dir)?;
        // Each contig's selection only reads the shared pileup VCF (the original runs them with
        // `parallel`).
        let selected = pool.install(|| {
            names
                .par_iter()
                .map(|contig| {
                    let options = HetSnpOptions { contig: Some((*contig).to_owned()), ..HetSnpOptions::default() };
                    let records = pileup_vcf.contig(contig)?;
                    Ok(select_het_snps_for_phasing_from(records, &phase_dir, &options)?.info_message())
                })
                .collect::<Result<Vec<String>>>()
        })?;
        for message in selected {
            log::info!("{message}");
        }
    }
    // Stage 3, phasing each contig, longest first and several at once. It runs alongside stage 6:
    // a contig's full-alignment jobs start once it is phased, while its reads are still cached.
    let phasing_jobs = concurrent_phasing(&settings.intermediate_phasing, threads);
    let phase = |contig: &str, per_job: usize| -> Result<()> {
        let variants = phase_dir.join(format!("{contig}.vcf"));
        match &settings.intermediate_phasing {
            Phaser::WhatsHap => Ok(phase_contig(&PhaseJob {
                variants: &variants,
                alignments: &args.bam,
                reference: &args.reference,
                contig,
                output: &phase_dir.join(format!("phased_{contig}.vcf.gz")),
                threads: per_job,
            })?),
            Phaser::LongPhase(bin) => {
                let prefix = phase_dir.join(format!("phased_{contig}"));
                longphase(bin, &variants, args, per_job, &prefix)?;
                bgzip(&phase_dir.join(format!("phased_{contig}.vcf")))?;
                Ok(())
            }
        }
    };

    // Stage 5: full-alignment candidates, the low-quality pileup calls. They depend on the pileup
    // calls only, so they are chosen before phasing.
    log::info!("[5/7] selecting candidates for full-alignment calling");
    progress.mark(FA_START - 0.03, "selecting full-alignment candidates");
    let fa_dir = tmp.join("full_alignment_output");
    let bed_dir = fa_dir.join("candidate_bed");
    fs::create_dir_all(&bed_dir).with_context(|| format!("creating {}", bed_dir.display()))?;
    let full_options = FullQualOptions {
        var_pct_full: settings.var_pct_full,
        ref_pct_full: settings.ref_pct_full,
        candidate_vcf: if args.vcf_fn.is_some() { CandidateVcf::Present } else { CandidateVcf::Absent },
    };
    let full_qual = select_full_qual(pileup_vcf.all()?, &full_options)?;
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
                    window: CandidateWindow::from_platform(platform),
                    ..CandidateOptions::default()
                };
                let records = pileup_vcf.contig(contig)?;
                let beds = select_candidates_for_full_alignment_from(records, &summary.vcf, &bed_dir, &options)?;
                Ok((beds.messages().to_vec(), beds.full_aln_file().map(Path::to_owned)))
            })
            .collect::<Result<Vec<_>>>()
    })?;
    drop(pileup_vcf);
    let mut lists: Vec<PathBuf> = Vec::new();
    for (messages, list) in per_contig {
        for message in messages {
            log::info!("{message}");
        }
        lists.extend(list);
    }
    // The final calls of one contig phased (`--enable_phasing`).
    let phase_final = |phaser: &Phaser, contig: &str, per_job: usize| -> Result<()> {
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
    };
    // The final phasing's end: each contig's merged calls phased (unless that was done as its
    // calls were scored), then sorted into phased_merge_output.vcf.gz, and the reads haplotagged.
    let finish_phasing = |phaser: &Phaser, phased_already: bool| -> Result<()> {
        if !phased_already {
            let tool = if *phaser == Phaser::WhatsHap { "lrpipe-phase" } else { "LongPhase" };
            log::info!("[7/7] Phasing VCF output in parallel using {tool}");
            for_each_contig(&contigs, threads, concurrent_phasing(phaser, threads), |contig, per_job| {
                phase_final(phaser, contig, per_job)
            })?;
        }
        let phased = sort_vcf(&sorted(&merge_dir, "phased_merge", args.output.join("phased_merge_output.vcf")))?;
        log::info!("{} phased calls in {}", phased.records, phased.vcf.display());
        if settings.haplotag {
            log::info!("[4/7] Haplotag input BAM file using lrpipe-phase");
            let bam = args.output.join("phased_output.bam");
            haplotag_output(&phased.vcf, &args.bam, &args.reference, &bam, threads)?;
            log::info!("Finish calling, phased output BAM file: {}", bam.display());
        }
        Ok(())
    };
    // Stage 7 per contig: full-alignment calls replace the pileup calls they re-called.
    let merge_contig = |contig: &str, full_alignment: &Path| -> Result<()> {
        let merged = merge_dir.join(format!("merge_{contig}.vcf"));
        let mut merge = MergeOptions::new(&summary.vcf, full_alignment, &merged);
        merge.contig = Some(contig.to_owned());
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
    };
    let finish_merging = |full_alignment: &Path, already_merged: bool| -> Result<_> {
        if let Some(gvcf) = &gvcf {
            gvcf.sort_chunks(&args.reference, &contigs_file, &cmd_file)?;
        }
        if !already_merged {
            log::info!("[7/7] merging pileup and full-alignment calls");
            pool.install(|| names.par_iter().try_for_each(|contig| merge_contig(contig, full_alignment)))?;
        }
        progress.mark(FA_START + FA_SHARE, "merging the calls");
        let final_sort = sorted(&merge_dir, "merge", args.output.join("merge_output.vcf"));
        let final_summary = sort_vcf(&final_sort)?;
        log::info!("{} calls in {} ({:.1?})", final_summary.records, final_summary.vcf.display(), started.elapsed());
        if let Some(gvcf) = &gvcf {
            gvcf.sort_merged(&merge_dir, &args.output, &final_sort)?;
        }
        Ok(final_summary)
    };
    // Without full-alignment calls, the pileup supplies every call. Haploid genotypes (for every
    // contig, or a male sample's chrX and chrY outside the PARs) and reference blocks still need
    // final merging. Plain diploid VCFs keep the pileup records.
    let finish_pileup_calls = || -> Result<()> {
        if gvcf.is_some() || settings.haploid != Haploid::Off || settings.gender == Gender::Male {
            let empty = fa_dir.join("no_full_alignment.vcf");
            let header = vcf_header(Some(&args.reference), Some(&cmd_file), &args.sample_name)?;
            fs::write(&empty, header).with_context(|| format!("writing {}", empty.display()))?;
            finish_merging(&empty, false)?;
        } else {
            finish_without_calls(&args.output, false)?;
            if settings.final_phasing.is_some() {
                split_by_contig(&args.output.join("merge_output.vcf.gz"), &merge_dir, &names, "merge")?;
            }
        }
        if let Some(phaser) = &settings.final_phasing {
            finish_phasing(phaser, false)?;
        }
        Ok(())
    };
    // `sorted(candidate_bed_path.glob("FULL_ALN_FILE_*"))`, concatenated into FULL_ALN_FILES.
    lists.sort();
    if lists.is_empty() {
        if !args.no_phasing_for_fa {
            log::info!("[3/7] phasing");
            for_each_contig(&contigs, threads, phasing_jobs, phase)?;
        }
        log::info!("No Candidate found! Exit in selecting full-alignment candidates");
        finish_pileup_calls()?;
        return remove_intermediate(args, &tmp);
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
    // GNU parallel's `{/.}`: the BED name without its last extension (`chr20.4_13` -> `chr20`).
    let job_contig = |bed: &Path| -> Result<(String, String)> {
        let name = bed.file_name().and_then(|n| n.to_str()).context("candidate BED without a name")?;
        let contig = Path::new(name).file_stem().and_then(|n| n.to_str()).unwrap_or(name);
        Ok((name.to_owned(), contig.to_owned()))
    };
    let mut jobs_of: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, bed) in fa_jobs.iter().enumerate() {
        jobs_of.entry(job_contig(bed)?.1).or_default().push(i);
    }

    // Stage 6: one full-alignment job per candidate BED, as `CallVariantsFromCffi
    // --full_aln_regions` jobs: one thread each and batches of 200.
    log::info!("[3/7] phasing, and [6/7] full-alignment calling: {} jobs on {threads} threads", fa_jobs.len());
    let fa_platform = args.platform.full_alignment();
    let mut features = FeatureOptions::new();
    features.min_mapq = u64::try_from(settings.min_mq)?;
    features.platform = fa_platform;
    features.indel_limit = settings.indel_limit;
    features.dwell_time = settings.dwell_time;
    features.canonical_read_order = canonical_read_order;
    let fa_job = |bed: &PathBuf| -> Result<(String, FullAlignmentJob)> {
        let (name, contig) = job_contig(bed)?;
        let mut job = FullAlignmentJob::new(&args.bam, &args.reference, &contig, bed).options(features.clone());
        if !args.no_phasing_for_fa {
            job = job.phased_vcf(phase_dir.join(format!("phased_{contig}.vcf.gz")));
        }
        Ok((name, job))
    };
    let fa_done = AtomicUsize::new(0);
    let fa_count = fa_jobs.len();
    progress.mark(FA_START, &format!("full-alignment calling: 0 of {fa_count} jobs"));
    // A job's tensors are built, on the CPU or the GPU.
    let fa_built = || {
        let done = fa_done.fetch_add(1, Ordering::Relaxed) + 1;
        let share = done as f64 / fa_count as f64;
        let what = format!("full-alignment calling: {done} of {fa_count} jobs");
        if done == fa_count {
            progress.mark(FA_START + FA_SHARE, &what);
        } else {
            progress.report(FA_START + FA_SHARE * share, &what);
        }
    };
    let fa_features = |bed: &PathBuf| -> Result<(String, FullAlignment)> {
        let (name, job) = fa_job(bed)?;
        let features = job.run().with_context(|| format!("full-alignment tensors for {name}"))?;
        fa_built();
        Ok((name, features))
    };
    let fa_model = &settings.fa_model;
    let depth = fa_platform.matrix_depth();
    let channels = settings.dwell_time.channels();
    let load_fa = |path: &Path, device| -> Result<FaModel> {
        let model = load_full_alignment_model(path, settings.model_config(device, settings.precision), depth)?;
        // The original builds the network for the channels the option asks for, and loading a
        // checkpoint of another width fails.
        if model.channels() != channels {
            bail!(
                "{} takes {} input channels, but {} tensors have {channels}{}",
                path.display(),
                model.channels(),
                if channels > snv_fullalign::CHANNELS { "dwell-time" } else { "full-alignment" },
                if channels > snv_fullalign::CHANNELS { "" } else { " (a _with_mv model needs --enable_dwell_time)" },
            );
        }
        Ok(model)
    };
    // On the GPU, with lrpipe-phase's final phasing, each contig is merged and phased as soon as
    // its full-alignment jobs are scored, alongside the other contigs' calling, from its own
    // full-alignment file; the files and records are those of merging and phasing at the end.
    let early_final = settings.device != InferenceDevice::Cpu
        && matches!(settings.final_phasing, Some(Phaser::WhatsHap))
        && gvcf.is_none()
        && args.platform != Platform::Ilmn;
    let contig_of: Vec<String> = fa_jobs.iter().map(|bed| job_contig(bed).map(|(_, c)| c)).collect::<Result<_>>()?;
    let finish_contig = |contig: &str, mut jobs: Vec<(usize, ScoredJob)>, per_job: usize| -> Result<()> {
        jobs.sort_by_key(|(i, _)| *i);
        let calls: Vec<_> = jobs.iter().map(|(_, (p, a, y))| (p.as_slice(), a.as_slice(), y.as_slice())).collect();
        pool.install(|| {
            output.write_full_alignment_jobs(&fa_dir.join(format!("full_alignment_{contig}.vcf")), &calls)
        })?;
        drop(jobs);
        // The prefix ends with the dot so that chr1's does not take chr10's file.
        let mut sort =
            sorted(&fa_dir, &format!("full_alignment_{contig}."), fa_dir.join(format!("contig_{contig}.vcf")));
        sort.layout = chunk_layout(settings.device);
        let full_alignment = sort_vcf(&sort)?;
        merge_contig(contig, &full_alignment.vcf)?;
        phase_final(&Phaser::WhatsHap, contig, per_job)
    };

    // Contigs merged and phased at once, as many as phasing runs, with as many threads. They run
    // beside it, so the peak memory grows (HG003 ONT: 13.5 GiB against 9.8); half as many jobs
    // with twice the threads peaked as high and took longer (679 s against 665 s in all).
    let early_jobs = if early_final { phasing_jobs } else { 0 };
    let per_early_job = threads / early_jobs.max(1);
    // Contigs whose full-alignment jobs are all scored, with them, for merging and phasing.
    let (early_tx, early_rx) = mpsc::channel::<(String, Vec<(usize, ScoredJob)>)>();
    let early_rx = Mutex::new(early_rx);

    let _trimmer = MemoryTrimmer::start(TRIM_INTERVAL);
    std::thread::scope(|scope| -> Result<()> {
        let early_workers: Vec<_> = (0..early_jobs)
            .map(|_| {
                scope.spawn(|| -> Result<()> {
                    loop {
                        // Its own statement, so the guard drops before the contig is phased.
                        let next = early_rx.lock().expect("lock poisoned").recv();
                        let Ok((contig, jobs)) = next else { return Ok(()) };
                        finish_contig(&contig, jobs, per_early_job.max(1))
                            .with_context(|| format!("merging and phasing {contig}"))?;
                    }
                })
            })
            .collect();
        let mut early_tx = early_final.then_some(early_tx);
        if let Some(tx) = &early_tx {
            log::info!("[7/7] merging and phasing each contig with lrpipe-phase as its calls finish");
            // Contigs without full-alignment jobs are ready now.
            for contig in names.iter().filter(|c| !jobs_of.contains_key(**c)) {
                let _ = tx.send(((*contig).to_owned(), Vec::new()));
            }
        }
        // Contigs whose phased VCF is written, in that order.
        let (phased_tx, phased_rx) = mpsc::channel::<String>();
        let phasing = scope.spawn(|| -> Result<()> {
            let phased_tx = phased_tx;
            if args.no_phasing_for_fa {
                for contig in &names {
                    // A closed channel means full alignment failed; its error is reported.
                    let _ = phased_tx.send((*contig).to_owned());
                }
                return Ok(());
            }
            for_each_contig(&contigs, threads, phasing_jobs, |contig, per_job| {
                phase(contig, per_job)?;
                let _ = phased_tx.send(contig.to_owned());
                Ok(())
            })?;
            log::info!("phasing done in {:.1?}", started.elapsed());
            Ok(())
        });
        let ready = phased_rx.iter().flat_map(|contig| jobs_of.get(&contig).cloned().unwrap_or_default());
        let calling = match settings.device {
            InferenceDevice::Cpu => {
                let times = StageTimes::default();
                let done = parallel::score_in_batches(
                    &pool,
                    ready.map(|i| &fa_jobs[i]),
                    Batching { tensor_len: full_alignment_tensor_len(depth, channels), batch: CPU_BATCH },
                    &times,
                    |bed| {
                        let (name, features) = fa_features(bed)?;
                        let (tensors, positions, alt_info) = features.into_parts();
                        Ok(((name, positions, alt_info), tensors))
                    },
                    |batch| FA_MODEL.with(fa_model, |p| load_fa(p, InferenceDevice::Cpu), |m| Ok(m.predict(batch)?)),
                    |_, (name, positions, alt_info), probs| {
                        let path = fa_dir.join(format!("full_alignment_{name}.vcf"));
                        output.write_full_alignment(&path, &positions, &alt_info, &probs)
                    },
                );
                log_times("full-alignment", &times);
                done
            }
            device => load_fa(fa_model, device).and_then(|mut model| {
                let mut scored: Vec<(usize, ScoredJob)> = Vec::new();
                // With early_tx: the jobs of contigs not yet complete.
                let mut pending: HashMap<&str, Vec<(usize, ScoredJob)>> = HashMap::new();
                // Batches of GPU_BATCH candidates from each job's start, as the whole job's
                // tensors would be split.
                score_full_alignment_on_gpu(
                    &pool,
                    ready.map(|i| (i, &fa_jobs[i])),
                    |bed, sink| {
                        let (name, job) = fa_job(bed)?;
                        let built = job
                            .run_batched(GPU_BATCH, sink)
                            .with_context(|| format!("full-alignment tensors for {name}"))?;
                        fa_built();
                        Ok(built)
                    },
                    |batch| Ok(model.predict(batch)?),
                    |i, job| {
                        let Some(tx) = &early_tx else {
                            scored.push((i, job));
                            return Ok(());
                        };
                        let contig = contig_of[i].as_str();
                        let jobs = pending.entry(contig).or_default();
                        jobs.push((i, job));
                        if jobs.len() == jobs_of.get(contig).map_or(0, Vec::len) {
                            let jobs = pending.remove(contig).unwrap_or_default();
                            tx.send((contig.to_owned(), jobs)).map_err(|_| anyhow!("final phasing stopped"))?;
                        }
                        Ok(())
                    },
                )?;
                if early_tx.is_some() {
                    return Ok(());
                }
                scored.sort_by_key(|(i, _)| *i);
                let jobs: Vec<_> =
                    scored.iter().map(|(_, (p, a, y))| (p.as_slice(), a.as_slice(), y.as_slice())).collect();
                pool.install(|| output.write_full_alignment_jobs(&fa_dir.join("full_alignment_0.vcf"), &jobs))
                    .map(|_| ())
            }),
        };
        // No more contigs: the workers finish those they have, and their errors come first, as
        // calling stops when one fails.
        early_tx.take();
        for worker in early_workers {
            worker.join().map_err(|_| anyhow!("final phasing thread panicked"))??;
        }
        // A phasing error stops the jobs early, so report it before what that causes.
        phasing.join().map_err(|_| anyhow!("phasing thread panicked"))??;
        calling
    })?;
    if early_final {
        log::info!("full-alignment calling, and merging and phasing each contig, done in {:.1?}", started.elapsed());
    } else {
        log::info!("full-alignment calling done in {:.1?}", started.elapsed());
    }

    let mut sort = sorted(&fa_dir, "full_alignment", args.output.join("full_alignment.vcf"));
    sort.layout = chunk_layout(settings.device);
    let fa_summary = sort_vcf(&sort)?;
    log::info!("merged {} full-alignment records into {}", fa_summary.records, fa_summary.vcf.display());
    if fa_summary.records == 0 {
        log::info!("Exit in full-alignment variant calling");
        finish_pileup_calls()?;
        return remove_intermediate(args, &tmp);
    }
    let final_summary = finish_merging(&fa_summary.vcf, early_final)?;
    if final_summary.records == 0 {
        log::info!("Exit in variant merging");
        finish_without_calls(&args.output, settings.final_phasing.is_some())?;
        return remove_intermediate(args, &tmp);
    }

    if let Some(phaser) = &settings.final_phasing {
        finish_phasing(phaser, early_final)?;
    }
    if let Some(vcf) = &args.vcf_fn {
        log::info!("Double check re-genotyping variants");
        let added = add_back_missing_variants(vcf, &final_summary.vcf, &args.output.join("merge_output.vcf"))?;
        log::info!("{added}");
    }
    remove_intermediate(args, &tmp)?;
    log::info!("Finish calling, output file: {}", final_summary.vcf.display());
    progress.report(1.0, "done");
    Ok(())
}

/// Pileup calling's share of a run, for its progress: 38% of HG002 ONT's whole genome.
const PILEUP_SHARE: f64 = 0.4;

/// Where full-alignment calling starts in a run's progress, and its share: 49% of HG002 ONT's
/// whole genome. Selecting and phasing the SNPs takes what is between pileup calling and it.
const FA_START: f64 = 0.48;
const FA_SHARE: f64 = 0.48;

/// Contigs phased at once. As lrpipe-phase itself does for several chromosomes, up to threads/2
/// run with 2 threads each (its output does not depend on the thread count); the original
/// runs LongPhase as `max(1, min(4, threads // 2))` jobs.
pub(crate) fn concurrent_phasing(phaser: &Phaser, threads: usize) -> usize {
    match phaser {
        Phaser::WhatsHap => (threads / 2).max(1),
        Phaser::LongPhase(_) => (threads / 2).clamp(1, 4),
    }
}

/// Makes the program `command` starts end when the thread starting it does, lrpipe-snv's
/// ending included (Linux's `PR_SET_PDEATHSIG`): a LongPhase run outlived an lrpipe-snv that was
/// killed, as `lrpipe run` kills a stage it stops, and held the machine. Elsewhere it does
/// nothing.
fn ends_with_this_thread(command: &mut std::process::Command) {
    #[cfg(target_os = "linux")]
    {
        use std::os::unix::process::CommandExt;
        // SAFETY: getpid is async-signal-safe.
        let parent = unsafe { libc::getpid() };
        let hook = move || {
            // SAFETY: prctl and getppid are async-signal-safe and touch no memory of ours.
            unsafe {
                if libc::prctl(libc::PR_SET_PDEATHSIG, libc::SIGKILL) != 0 {
                    return Err(std::io::Error::last_os_error());
                }
                // lrpipe-snv may have gone before the request took hold.
                if libc::getppid() != parent {
                    return Err(std::io::Error::from_raw_os_error(libc::ESRCH));
                }
            }
            Ok(())
        };
        // SAFETY: the hook runs between fork and exec, calling only async-signal-safe functions
        // and allocating nothing (its errors are OS error codes).
        unsafe {
            command.pre_exec(hook);
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = command;
}

/// `longphase phase -s VARIANTS -b BAM -r REF -t THREADS -o PREFIX --ont|--pb`, which writes
/// `PREFIX.vcf`.
pub(crate) fn longphase(bin: &Path, variants: &Path, args: &RunArgs, threads: usize, prefix: &Path) -> Result<()> {
    let platform = if args.platform == Platform::Ont { "--ont" } else { "--pb" };
    let mut command = std::process::Command::new(bin);
    ends_with_this_thread(&mut command);
    let output = command
        .arg("phase")
        .arg("-s")
        .arg(variants)
        .arg("-b")
        .arg(&args.bam)
        .arg("-r")
        .arg(&args.reference)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-o")
        .arg(prefix)
        .arg(platform)
        .output()
        .with_context(|| format!("running {}", bin.display()))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("LongPhase failed on {} ({}): {}", variants.display(), output.status, stderr.trim_end());
    }
    Ok(())
}

/// Runs `job` on every contig, longest first, `concurrent` at a time, each with an equal share
/// of `threads`.
pub(crate) fn for_each_contig(
    contigs: &Contigs,
    threads: usize,
    concurrent: usize,
    job: impl Fn(&str, usize) -> Result<()> + Sync,
) -> Result<()> {
    let mut by_length: Vec<&(String, i64)> = contigs.chunks.iter().collect();
    by_length.sort_by_key(|(_, chunks)| std::cmp::Reverse(*chunks));
    let concurrent = concurrent.clamp(1, by_length.len().max(1));
    let per_job = (threads / concurrent).max(1);
    let queue = Mutex::new(by_length.into_iter());
    std::thread::scope(|scope| -> Result<()> {
        let workers: Vec<_> = (0..concurrent)
            .map(|_| {
                scope.spawn(|| -> Result<()> {
                    loop {
                        // Its own statement, so the guard drops before the job; a `while let`
                        // would hold it for the whole body and run the contigs one at a time.
                        let next = queue.lock().expect("lock poisoned").next();
                        let Some((contig, _)) = next else { break };
                        job(contig, per_job)?;
                    }
                    Ok(())
                })
            })
            .collect();
        for worker in workers {
            worker.join().map_err(|_| anyhow!("worker thread panicked"))??;
        }
        Ok(())
    })
}

/// `--chunk_num -1` with `--bed_fn`: the original's one pileup job, for a contig named `None`,
/// calls each of the BED's `regions` in turn into `pileup_None_0.vcf` (and their gVCF blocks),
/// which the sort step then ignores, as it does not name a contig. The original's GPU calling
/// builds no tensors for that job, so it is written on the CPU only.
fn unchunked_pileup(
    pool: &rayon::ThreadPool,
    regions: &[(String, i64, i64)],
    region_candidates: &(impl Fn(&str, &ChunkRegion, Option<&Path>) -> Result<Vec<Candidate>> + Sync),
    output: &VcfOutput<'_>,
    (settings, config): (&Settings, ModelConfig),
    pileup_dir: &Path,
) -> Result<()> {
    if settings.device != InferenceDevice::Cpu {
        return Ok(());
    }
    log::info!("calling the {} BED regions of the unchunked job", regions.len());
    let calls: Vec<(Vec<Candidate>, Vec<f32>)> = pool.install(|| {
        regions
            .par_iter()
            .map(|(contig, start, end)| {
                let region = ChunkRegion::Range { start: *start, end: *end };
                let candidates = region_candidates(contig, &region, None)
                    .with_context(|| format!("pileup region {contig}:{start}-{end}"))?;
                let inputs = chunk::model_inputs(&candidates)?;
                let probs = PILEUP_MODEL.with(
                    &settings.pileup_model,
                    |p| Ok(load_pileup_model(p, config)?),
                    |m| chunk::score(m, &inputs, CPU_BATCH),
                )?;
                Ok((candidates, probs))
            })
            .collect::<Result<_>>()
    })?;
    let jobs: Vec<_> = calls.iter().map(|(c, p)| (c.as_slice(), p.as_slice())).collect();
    pool.install(|| output.write_each(&pileup_dir.join("pileup_None_0.vcf"), &jobs))
}

/// `--no_phasing_for_fa`: `phase_bam/<contig>.bam` (and its index) link to the input BAM.
///
/// The original links only a `.bai`, which leaves a CRAM (or a BAM indexed with `.csi`) without an
/// index there; its `.crai` or `.csi` is linked as well, under the name htslib looks for.
pub(crate) fn link_unphased_bams(args: &RunArgs, phase_bam: &Path, contigs: &[&str]) -> Result<()> {
    fs::create_dir_all(phase_bam).with_context(|| format!("creating {}", phase_bam.display()))?;
    let index = ["bai", "csi", "crai"].into_iter().find_map(|extension| {
        let mut appended = args.bam.as_os_str().to_owned();
        appended.push(format!(".{extension}"));
        [PathBuf::from(appended), args.bam.with_extension(extension)]
            .into_iter()
            .find(|p| p.exists())
            .map(|path| (path, extension))
    });
    let link = |target: &Path, link: PathBuf| -> Result<()> {
        if link.symlink_metadata().is_ok() {
            fs::remove_file(&link)?;
        }
        std::os::unix::fs::symlink(target, &link).with_context(|| format!("linking {}", link.display()))
    };
    for contig in contigs {
        link(&args.bam, phase_bam.join(format!("{contig}.bam")))?;
        if let Some((index, extension)) = &index {
            link(index, phase_bam.join(format!("{contig}.bam.{extension}")))?;
        }
    }
    Ok(())
}

/// `--remove_intermediate_dir`.
pub(crate) fn remove_intermediate(args: &RunArgs, tmp: &Path) -> Result<()> {
    if args.remove_intermediate_dir {
        log::info!("Removing intermediate files in {}", tmp.display());
        // `shutil.rmtree(tmp_path, ignore_errors=True)`.
        let _ = fs::remove_dir_all(tmp);
    }
    Ok(())
}

/// Sub-chunks of a chunk the 8-bit check takes one from.
const CHECK_SPLIT: i64 = 64;

/// Chunks the 8-bit check samples, spread over the run's.
const CHECK_CHUNKS: usize = 4;

/// [`rescore::check_eight_bit`] on candidates from the middle sixty-fourth of a few chunks spread
/// over `jobs`: `Err` says why the 8-bit pileup network is not used.
fn check_eight_bit_pileup(
    pool: &rayon::ThreadPool,
    jobs: &[(&str, i64, i64)],
    contigs: &Contigs,
    region_candidates: &(impl Fn(&str, &ChunkRegion, Option<&Path>) -> Result<Vec<Candidate>> + Sync),
    output: &VcfOutput<'_>,
    settings: &Settings,
) -> Result<std::result::Result<(), String>> {
    let step = jobs.len().div_ceil(CHECK_CHUNKS).max(1);
    let sampled: Vec<&(&str, i64, i64)> = jobs.iter().skip(step / 2).step_by(step).collect();
    let parts = pool.install(|| {
        sampled
            .par_iter()
            .map(|&&(contig, id, count)| {
                let extend_bed = contigs.extend_bed(contig);
                let (id, count) = ((id - 1) * CHECK_SPLIT + CHECK_SPLIT / 2, count * CHECK_SPLIT);
                let region = match extend_bed {
                    Some(_) => ChunkRegion::BedChunk { id, count },
                    None => ChunkRegion::FaiChunk { id, count },
                };
                region_candidates(contig, &region, extend_bed.as_deref())
            })
            .collect::<Result<Vec<Vec<Candidate>>>>()
    })?;
    let candidates: Vec<Candidate> = parts.into_iter().flatten().collect();
    if candidates.is_empty() {
        return Ok(Err("its sampled regions held no candidates".to_owned()));
    }
    let inputs = chunk::model_inputs(&candidates)?;
    let config = |precision| settings.model_config(InferenceDevice::Cpu, precision);
    let mut mixed = load_pileup_model(&settings.pileup_model, config(Precision::Mixed))?;
    let mut exact = load_pileup_model(&settings.pileup_model, config(Precision::Float32))?;
    let started = Instant::now();
    let checked = rescore::check_eight_bit(&candidates, &inputs, output, &mut mixed, &mut exact, CPU_BATCH)?;
    log::info!("checked the 8-bit pileup network on {} candidates in {:.1?}", candidates.len(), started.elapsed());
    Ok(checked)
}

/// A pileup network on whichever backend the run uses.
pub(crate) type PileupNet = Box<dyn PileupModel + Send>;
/// A full-alignment network on whichever backend the run uses.
pub(crate) type FaModel = Box<dyn FullAlignmentModel + Send>;

thread_local! {
    static PILEUP_CELL: RefCell<Option<(PathBuf, PileupNet)>> = const { RefCell::new(None) };
    static FA_CELL: RefCell<Option<(PathBuf, FaModel)>> = const { RefCell::new(None) };
}
pub(crate) static PILEUP_MODEL: ThreadModel<PileupNet> = ThreadModel::new(&PILEUP_CELL);
pub(crate) static FA_MODEL: ThreadModel<FaModel> = ThreadModel::new(&FA_CELL);

pub(crate) fn log_times(stage: &str, times: &StageTimes) {
    let (build, score) = times.totals();
    log::info!("{stage}: {:.1} CPU-s building inputs, {:.1} CPU-s inference", build.as_secs_f64(), score.as_secs_f64());
}

/// One contig's records in a bgzipped VCF: where they start, and their uncompressed byte range.
struct Span {
    start: VirtualPosition,
    bytes: Range<u64>,
}

/// A bgzipped pileup VCF, its header and where each contig's records are in it. The records stay
/// on disk: noisy reads give a whole genome tens of GB of them (HG002 ONT at 10x: 265M records,
/// 18.6 GB, nearly all RefCall), more than the memory of the machines that call them.
struct PileupVcf {
    path: PathBuf,
    header: Vec<u8>,
    /// Each contig's records; one span per contig, as the file is sorted.
    spans: HashMap<String, Vec<Span>>,
}

impl PileupVcf {
    fn index(path: &Path) -> Result<Self> {
        let context = || format!("reading {}", path.display());
        let mut reader = bgzf::io::Reader::new(File::open(path).with_context(context)?);
        let mut spans: HashMap<String, Vec<Span>> = HashMap::new();
        let (mut header, mut line, mut last) = (Vec::new(), Vec::new(), None::<Vec<u8>>);
        let mut end = 0;
        loop {
            let start = reader.virtual_position();
            line.clear();
            let read = reader.read_until(b'\n', &mut line).with_context(context)?;
            if read == 0 {
                break;
            }
            let begin = end;
            end += read as u64;
            if line.starts_with(b"#") {
                if last.is_none() {
                    header.extend_from_slice(&line);
                }
                continue;
            }
            let contig = line.split(|&b| b == b'\t').next().unwrap_or_default();
            let contig_spans = spans.entry(String::from_utf8_lossy(contig).into_owned()).or_default();
            match contig_spans.last_mut() {
                Some(span) if last.as_deref() == Some(contig) => span.bytes.end = end,
                _ => contig_spans.push(Span { start, bytes: begin..end }),
            }
            if last.as_deref() != Some(contig) {
                last = Some(contig.to_vec());
            }
        }
        Ok(Self { path: path.to_owned(), header, spans })
    }

    fn all(&self) -> Result<impl BufRead + '_> {
        Ok(open_vcf(&self.path)?)
    }

    /// The header and the records of `contig`, or the whole file when its records are not in
    /// one piece: the selections filter by contig either way.
    fn contig(&self, contig: &str) -> Result<Box<dyn BufRead + '_>> {
        let header = self.header.as_slice();
        Ok(match self.spans.get(contig).map(Vec::as_slice) {
            None => Box::new(header),
            Some([span]) => {
                let file = File::open(&self.path).with_context(|| format!("reading {}", self.path.display()))?;
                let mut reader = bgzf::io::Reader::new(file);
                reader.seek(span.start)?;
                Box::new(header.chain(reader.take(span.bytes.end - span.bytes.start)))
            }
            Some(_) => Box::new(self.all()?),
        })
    }
}

/// The pileup network's inputs on a GPU: the CPU path's, with its high-depth scaling.
///
/// Clair3's own GPU calling scores the tensors it saved as `np.int8`, so counts past 127 wrap,
/// and skips that scaling. Doing as it did made CUDA calls collapse at high depth (K-12 HiFi at
/// 295x: SNP F1 42.5% against the CPU's 96.4%, issue #45), so every device takes these inputs.
fn gpu_pileup_inputs(candidates: &[Candidate]) -> Result<Vec<i32>> {
    chunk::model_inputs(candidates)
}

/// `SortVcf --use_gpu`: the original's GPU calling writes one VCF per GPU, not one per chunk.
fn chunk_layout(device: InferenceDevice) -> ChunkLayout {
    match device {
        InferenceDevice::Cpu => ChunkLayout::ByContig,
        _ => ChunkLayout::Unordered,
    }
}

/// `ensure_merge_output`: when a stage ends the run early, the pileup calls become the final
/// `merge_output.vcf.gz`.
pub(crate) fn ensure_merge_output(output: &Path) -> Result<()> {
    // `clear_previous_run` removed any earlier run's merged VCF, so one here is this run's.
    let merged = output.join("merge_output.vcf.gz");
    let pileup = output.join("pileup.vcf.gz");
    if !merged.exists() && pileup.exists() {
        log::info!("No merged VCF produced, finalizing pileup output as {}", merged.display());
        fs::copy(&pileup, &merged)?;
        // `tabix -f -p vcf` on a copy of an indexed file writes the same index.
        let pileup_index = output.join("pileup.vcf.gz.tbi");
        if pileup_index.exists() {
            fs::copy(pileup_index, output.join("merge_output.vcf.gz.tbi"))?;
        }
    }
    Ok(())
}

/// Ends a run that has no calls with every file it was asked for: the merged VCF and, with final
/// phasing, the phased one (`phased_merge_output.vcf.gz`), both without records and indexed.
/// Clair3 writes no phased VCF there, so a run asked to phase ended without the file it
/// promised (LRP-020).
fn finish_without_calls(output: &Path, phasing: bool) -> Result<()> {
    ensure_merge_output(output)?;
    let merged = output.join("merge_output.vcf.gz");
    let index = output.join("merge_output.vcf.gz.tbi");
    if merged.exists() && !index.exists() && !snv_merge::tabix_index(&merged)? {
        bail!("could not index {}", merged.display());
    }
    if phasing && merged.exists() {
        fs::copy(&merged, output.join("phased_merge_output.vcf.gz"))?;
        if index.exists() {
            fs::copy(index, output.join("phased_merge_output.vcf.gz.tbi"))?;
        }
    }
    Ok(())
}

/// Writes `vcf`'s records of each of `contigs` to `dir/{prefix}_{contig}.vcf` under its header,
/// a file for each contig, with or without records.
fn split_by_contig(vcf: &Path, dir: &Path, contigs: &[&str], prefix: &str) -> Result<()> {
    let mut header = String::new();
    let mut records: HashMap<&str, String> = contigs.iter().map(|c| (*c, String::new())).collect();
    for line in open_vcf(vcf)?.lines() {
        let line = line.with_context(|| format!("reading {}", vcf.display()))?;
        if line.starts_with('#') {
            header.push_str(&line);
            header.push('\n');
        } else if let Some(out) = line.split('\t').next().and_then(|c| records.get_mut(c)) {
            out.push_str(&line);
            out.push('\n');
        }
    }
    fs::create_dir_all(dir)?;
    for (contig, text) in records {
        fs::write(dir.join(format!("{prefix}_{contig}.vcf")), format!("{header}{text}"))?;
    }
    Ok(())
}

/// The files a run writes in its output directory besides `tmp/`, with their indexes.
const RUN_OUTPUTS: [&str; 18] = [
    "pileup.vcf",
    "pileup.vcf.gz",
    "pileup.vcf.gz.tbi",
    "full_alignment.vcf",
    "full_alignment.vcf.gz",
    "full_alignment.vcf.gz.tbi",
    "merge_output.vcf",
    "merge_output.vcf.gz",
    "merge_output.vcf.gz.tbi",
    "phased_merge_output.vcf",
    "phased_merge_output.vcf.gz",
    "phased_merge_output.vcf.gz.tbi",
    "merge_output.gvcf",
    "merge_output.gvcf.gz",
    "merge_output.gvcf.gz.tbi",
    "phased_output.bam",
    "phased_output.bam.bai",
    "phased_output.bam.csi",
];

/// The file that says lrpipe-snv made what is beside it, and so may remove it: in the output
/// directory it covers the [`RUN_OUTPUTS`] there, and in `tmp/` that `tmp/`. A run writes both
/// before it writes anything else. `tmp/` has its own, so a `tmp/` made later in a directory
/// lrpipe-snv once wrote to (a project's, with `--output .`) is still someone else's.
pub(crate) const OWNER_MARKER: &str = ".lrpipe-snv";

/// Checks that each input a run reads can be opened (the BAM or CRAM and its index, the FASTA
/// and its `.fai`, a `--bed_fn` or `--vcf_fn`, and the model's checkpoints) before the run
/// changes anything in its output directory: a mistyped path must not cost an earlier run's
/// results (LSNV-001).
///
/// # Errors
///
/// Naming the input that cannot be read.
///
/// Returns `--bed-fn`'s rows, read once here.
pub(crate) fn check_inputs(args: &RunArgs, settings: &Settings) -> Result<Option<Vec<(String, i64, i64)>>> {
    AlignmentFile::open(&args.bam, &args.reference)
        .with_context(|| format!("checking --bam {} (and its index)", args.bam.display()))?;
    Fasta::open(&args.reference)
        .with_context(|| format!("checking --ref {} (and its .fai)", args.reference.display()))?;
    if let Some(index) = lrpipe_reads::stale_index(&args.bam) {
        log::warn!(
            "the index {} is older than --bam {}: if the reads were written again since, index them again (`samtools \
             index`), as regions read through an old index can come out empty",
            index.display(),
            args.bam.display()
        );
    }
    lrpipe_reads::check_fai(&args.reference)
        .map_err(lrpipe_reads::OpenError::into_io_error)
        .with_context(|| format!("checking --ref {} (and its .fai)", args.reference.display()))?;
    for (option, path) in [("--bed-fn", &args.bed_fn), ("--vcf-fn", &args.vcf_fn)] {
        if let Some(path) = path {
            readable_file(path).with_context(|| format!("checking {option} {}", path.display()))?;
        }
    }
    let mut checkpoints = vec![&settings.pileup_model];
    if !args.pileup_only {
        checkpoints.push(&settings.fa_model);
    }
    for checkpoint in checkpoints {
        readable_file(checkpoint)
            .with_context(|| format!("checking the model in --model-path {}", args.model_path.display()))?;
    }
    // Everything else a run reads before it writes, checked here so that a mistake fails with the
    // previous run's results intact: they are removed next.
    let bed = args
        .bed_fn
        .as_deref()
        .map(|bed| crate::contigs::bed_rows(bed).with_context(|| format!("checking --bed-fn {}", bed.display())))
        .transpose()?;
    if args.sample_name.is_empty() || args.sample_name.contains(['\t', '\n', '\r']) {
        bail!("--sample-name {:?}: a sample name must be non-empty and have no tabs or line breaks", args.sample_name);
    }
    let cpu = settings.model_config(InferenceDevice::Cpu, Precision::Float32);
    load_pileup_model(&settings.pileup_model, cpu)
        .with_context(|| format!("checking the pileup model {}", settings.pileup_model.display()))?;
    if !args.pileup_only {
        load_full_alignment_model(&settings.fa_model, cpu, snv_infer::FULL_ALIGNMENT_DEPTH)
            .with_context(|| format!("checking the full-alignment model {}", settings.fa_model.display()))?;
    }
    // A device this build cannot use (CUDA without `gpu`, Metal without `metal` on macOS).
    settings.model_config(settings.device, settings.precision).check().context("checking --device")?;
    #[cfg(feature = "gpu")]
    if let InferenceDevice::Cuda(device) = settings.device {
        snv_infer::check_cuda_device(device).context("checking --device")?;
    }
    Ok(bed)
}

/// Fails unless `path` is a file this process can open for reading.
fn readable_file(path: &Path) -> Result<()> {
    let file = File::open(path).with_context(|| format!("cannot read {}", path.display()))?;
    if !file.metadata().with_context(|| format!("cannot read {}", path.display()))?.is_file() {
        bail!("{} is not a file", path.display());
    }
    Ok(())
}

/// Removes what an earlier run left in `output`: `tmp/` and [`RUN_OUTPUTS`]. Each run's chunks
/// are found by listing `tmp/`, and its early exits keep a merged VCF they find, so without this
/// an interrupted run's chunks, or a final VCF from other settings, could join this run's
/// results or stay beside them (LRP-002). Clair3 leaves them; lrpipe-snv does not.
///
/// Only what an [`OWNER_MARKER`] covers is removed: `tmp/` with one inside it, the results with
/// one in `output`. Any other `tmp/`, or file named like a result, is someone else's (a project's
/// own `tmp/` with `--output .`, an earlier Clair3 run), and the run is refused with nothing
/// changed (LRUN-09). Then `output` and a new `tmp/` get their markers, before the run writes
/// anything else there.
///
/// # Errors
///
/// If `output` holds such files, or they cannot be removed, or the markers cannot be written.
pub(crate) fn clear_previous_run(output: &Path) -> Result<()> {
    let tmp = output.join("tmp");
    let marker = output.join(OWNER_MARKER);
    let tmp_marker = tmp.join(OWNER_MARKER);
    let mut found = Vec::new();
    // An empty `tmp/` holds nothing to lose: a run killed between making its `tmp/` and marking
    // it leaves one, which must not turn every later run away.
    let empty_dir = |dir: &Path| fs::read_dir(dir).is_ok_and(|mut entries| entries.next().is_none());
    if tmp.symlink_metadata().is_ok() && !tmp_marker.is_file() && !empty_dir(&tmp) {
        found.push("tmp");
    }
    if !marker.is_file() {
        found.extend(RUN_OUTPUTS.into_iter().filter(|name| output.join(name).symlink_metadata().is_ok()));
    }
    if !found.is_empty() {
        bail!(
            "{} already holds {}, which lrpipe-snv did not write (there is no {OWNER_MARKER} beside them), and a \
             run would replace them; give another --output, or move them away",
            output.display(),
            found.join(", ")
        );
    }
    if tmp.symlink_metadata().is_ok() {
        fs::remove_dir_all(&tmp).with_context(|| format!("removing the previous run's {}", tmp.display()))?;
    }
    for name in RUN_OUTPUTS {
        let path = output.join(name);
        match fs::remove_file(&path) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                return Err(e).with_context(|| format!("removing the previous run's {}", path.display()));
            }
            _ => {}
        }
    }
    fs::create_dir_all(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
    for (path, covers) in [(&marker, "its results"), (&tmp_marker, "this tmp/")] {
        fs::write(path, format!("lrpipe-snv wrote {covers} here, and a run into this directory replaces them.\n"))
            .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Whether process `pid` is gone or a zombie (killed, not yet reaped by its new parent).
    #[cfg(target_os = "linux")]
    fn ended(pid: u32) -> bool {
        fs::read_to_string(format!("/proc/{pid}/stat"))
            .map_or(true, |s| s.rsplit_once(") ").is_some_and(|(_, rest)| rest.starts_with('Z')))
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn longphase_ends_with_the_thread_that_started_it() {
        // A LongPhase run outlived an lrpipe-snv that was killed: its starting thread's end
        // ends it now.
        let pid = std::thread::spawn(|| {
            let mut command = std::process::Command::new("sleep");
            command.arg("30");
            ends_with_this_thread(&mut command);
            command.spawn().expect("sleep runs").id()
        })
        .join()
        .expect("the thread ran");
        let started = Instant::now();
        while !ended(pid) {
            assert!(started.elapsed() < Duration::from_secs(5), "process {pid} outlived its thread");
            std::thread::sleep(Duration::from_millis(20));
        }
    }

    #[test]
    fn gpu_pileup_inputs_are_the_cpus_at_high_depth() {
        // At 295x a CUDA run wrapped counts past 127 to `int8` and skipped the high-depth
        // scaling, and its calls collapsed (#45): every device takes the CPU path's inputs.
        let tensor: Vec<i32> = (0..TENSOR_LEN).map(|i| [0, 12, 140, 295, -150][i % 5]).collect();
        let candidates = vec![
            Candidate { tensor: tensor.clone(), position: "c:1001:A".to_owned(), alt_info: "295-RA".to_owned() },
            Candidate { tensor, position: "c:1002:A".to_owned(), alt_info: "59-RA".to_owned() },
        ];
        let gpu = gpu_pileup_inputs(&candidates).expect("inputs");
        assert_eq!(gpu[..5], [0, 5, 68, 144, -73], "the deep candidate is scaled, not wrapped");
        assert_eq!(gpu[TENSOR_LEN..TENSOR_LEN + 5], [0, 12, 140, 295, -150], "the shallow one is left as it is");
        assert_eq!(gpu, chunk::model_inputs(&candidates).expect("inputs"));
    }

    #[test]
    fn a_run_keeps_nothing_an_earlier_run_left() {
        // LRP-002: an earlier run's chunk and its final VCF and index are there. After clearing,
        // a pileup-only exit finalizes this run's pileup VCF with its own index, and no chunk
        // of the earlier run is left for this run's sorting to find.
        let out = std::env::temp_dir().join(format!("lrpipe-snv-clear-{}", std::process::id()));
        let _ = fs::remove_dir_all(&out);
        let chunks = out.join("tmp/pileup_output");
        fs::create_dir_all(&chunks).expect("dirs");
        fs::write(chunks.join("pileup_chr1.1_2.vcf"), "old chunk").expect("chunk");
        for (name, bytes) in [("merge_output.vcf.gz", "old vcf"), ("merge_output.vcf.gz.tbi", "old index")] {
            fs::write(out.join(name), bytes).expect("old output");
        }
        fs::write(out.join("notes.txt"), "the user's").expect("other file");
        fs::write(out.join(OWNER_MARKER), "").expect("the earlier run's marker");
        fs::write(out.join("tmp").join(OWNER_MARKER), "").expect("the earlier run's tmp/ marker");

        clear_previous_run(&out).expect("cleared");
        assert!(out.join(OWNER_MARKER).is_file(), "this run's marker");
        assert!(out.join("tmp").join(OWNER_MARKER).is_file(), "this run's tmp/ marker");
        assert!(!chunks.exists(), "the earlier run's chunks are gone");
        assert!(!out.join("merge_output.vcf.gz").exists() && !out.join("merge_output.vcf.gz.tbi").exists());
        assert!(out.join("notes.txt").exists(), "other files stay");

        fs::write(out.join("pileup.vcf.gz"), "new vcf").expect("pileup");
        fs::write(out.join("pileup.vcf.gz.tbi"), "new index").expect("pileup index");
        ensure_merge_output(&out).expect("finalized");
        assert_eq!(fs::read(out.join("merge_output.vcf.gz")).expect("merged"), b"new vcf");
        assert_eq!(fs::read(out.join("merge_output.vcf.gz.tbi")).expect("index"), b"new index");
        let _ = fs::remove_dir_all(&out);
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("lrpipe-snv-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).expect("dir");
        dir
    }

    #[test]
    fn a_tmp_lrpipe_snv_did_not_make_is_refused_untouched() {
        // LRUN-09: `--output .` in a project with its own tmp/ removed it, and files named like
        // Clair3's results, before checking anything.
        let out = scratch("foreign-tmp");
        fs::create_dir_all(out.join("tmp/important")).expect("dirs");
        fs::write(out.join("tmp/important/thesis.txt"), "precious").expect("user file");
        let err = clear_previous_run(&out).expect_err("refused");
        assert!(format!("{err:#}").contains("lrpipe-snv did not write"), "{err:#}");
        assert_eq!(fs::read_to_string(out.join("tmp/important/thesis.txt")).expect("kept"), "precious");
        assert!(!out.join(OWNER_MARKER).exists(), "the directory is not claimed");

        // A result of another tool's (an earlier Clair3 run's) without a tmp/.
        fs::remove_dir_all(out.join("tmp")).expect("tmp");
        fs::write(out.join("merge_output.vcf.gz"), "Clair3's calls").expect("calls");
        let err = clear_previous_run(&out).expect_err("refused");
        assert!(format!("{err:#}").contains("merge_output.vcf.gz"), "{err:#}");
        assert_eq!(fs::read_to_string(out.join("merge_output.vcf.gz")).expect("kept"), "Clair3's calls");

        // A tmp/ made after a run that left none (`--remove_intermediate_dir`) in a directory
        // lrpipe-snv wrote to is still someone else's.
        fs::remove_file(out.join("merge_output.vcf.gz")).expect("calls");
        clear_previous_run(&out).expect("an empty directory is cleared");
        fs::remove_dir_all(out.join("tmp")).expect("the run removes its tmp/");
        fs::create_dir_all(out.join("tmp/important")).expect("dirs");
        fs::write(out.join("tmp/important/thesis.txt"), "precious").expect("user file");
        let err = clear_previous_run(&out).expect_err("refused");
        assert!(format!("{err:#}").contains("tmp"), "{err:#}");
        assert_eq!(fs::read_to_string(out.join("tmp/important/thesis.txt")).expect("kept"), "precious");
        let _ = fs::remove_dir_all(&out);
    }

    #[test]
    fn an_empty_unmarked_tmp_does_not_block_later_runs() {
        // A run killed between making tmp/ and marking it leaves an empty, unmarked tmp/; it
        // holds nothing to lose, so the next run goes ahead and marks it.
        let out = scratch("empty-tmp");
        fs::create_dir_all(out.join("tmp")).expect("tmp");
        fs::write(out.join(OWNER_MARKER), "").expect("the killed run's marker");
        clear_previous_run(&out).expect("an empty tmp/ is no one's");
        assert!(out.join("tmp").join(OWNER_MARKER).is_file(), "marked");
        let _ = fs::remove_dir_all(&out);
    }

    fn run_args(line: &[&std::ffi::OsStr]) -> RunArgs {
        use clap::Parser;
        #[derive(Debug, Parser)]
        struct Cli {
            #[command(flatten)]
            run: RunArgs,
        }
        Cli::try_parse_from(std::iter::once(std::ffi::OsStr::new("lrpipe-snv")).chain(line.iter().copied()))
            .expect("parses")
            .run
    }

    #[test]
    fn a_mistyped_input_keeps_the_previous_runs_results() {
        // LSNV-001: a rerun with a mistyped --model-path, --bam or --ref removed the earlier run's
        // results before the input was found missing.
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../ci/fixture");
        let out = scratch("mistyped");
        let model = out.join("model");
        fs::create_dir_all(&model).expect("model dir");
        for name in ["pileup.pt", "full_alignment.pt"] {
            fs::write(model.join(name), "").expect("checkpoint");
        }
        let results = out.join("results");
        fs::create_dir_all(results.join("tmp/pileup_output")).expect("dirs");
        fs::write(results.join("tmp/pileup_output/pileup_chr20.1_1.vcf"), "chunk").expect("chunk");
        fs::write(results.join("merge_output.vcf.gz"), "earlier calls").expect("calls");
        fs::write(results.join(OWNER_MARKER), "").expect("marker");
        fs::write(results.join("tmp").join(OWNER_MARKER), "").expect("tmp/ marker");
        let (bam, reference) = (fixture.join("ont.bam"), fixture.join("chr20_window.fa.gz"));
        let missing = out.join("typo");
        let cases: [(&Path, &Path, &Path, &str); 4] = [
            (&bam, &reference, &missing, "--model-path"),
            (&bam, &reference, &model.join("pileup.pt"), "--model-path"),
            (&missing, &reference, &model, "--bam"),
            (&bam, &missing, &model, "--ref"),
        ];
        for (bam, reference, model, option) in cases {
            let args = run_args(&[
                "--bam".as_ref(),
                bam.as_os_str(),
                "--ref".as_ref(),
                reference.as_os_str(),
                "--model-path".as_ref(),
                model.as_os_str(),
                "--output".as_ref(),
                results.as_os_str(),
                "--threads".as_ref(),
                "1".as_ref(),
            ]);
            let err = run(&args).expect_err("an input is missing");
            assert!(format!("{err:#}").contains(option), "{option}: {err:#}");
            assert_eq!(fs::read_to_string(results.join("merge_output.vcf.gz")).expect("kept"), "earlier calls");
            assert!(results.join("tmp/pileup_output/pileup_chr20.1_1.vcf").exists(), "the chunks are kept");
        }
        let _ = fs::remove_dir_all(&out);
    }

    #[test]
    fn a_merged_vcf_this_run_wrote_keeps_its_own_index() {
        let out = std::env::temp_dir().join(format!("lrpipe-snv-keep-{}", std::process::id()));
        let _ = fs::remove_dir_all(&out);
        fs::create_dir_all(&out).expect("dir");
        for (name, bytes) in [
            ("merge_output.vcf.gz", "merged"),
            ("merge_output.vcf.gz.tbi", "merged index"),
            ("pileup.vcf.gz", "pileup"),
            ("pileup.vcf.gz.tbi", "pileup index"),
        ] {
            fs::write(out.join(name), bytes).expect("file");
        }
        ensure_merge_output(&out).expect("kept");
        assert_eq!(fs::read(out.join("merge_output.vcf.gz.tbi")).expect("index"), b"merged index");
        let _ = fs::remove_dir_all(&out);
    }

    fn read_all(mut reader: impl BufRead) -> String {
        let mut text = String::new();
        reader.read_to_string(&mut text).unwrap();
        text
    }

    #[test]
    fn pileup_vcf_reads_each_contig_from_its_blocks() {
        let dir = std::env::temp_dir().join(format!("lrpipe-snv-pileup-vcf-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\n";
        // chr1 spans several BGZF blocks; chr3's records are in two pieces.
        let rows = |contig: &str, n: usize| {
            (0..n).map(|i| format!("{contig}\t{i}\t.\tA\tC\t5\tRefCall\n")).collect::<String>()
        };
        let (chr1, chr2) = (rows("chr1", 20_000), rows("chr2", 3));
        let (chr3a, chr3b) = (rows("chr3", 2), rows("chr3", 1));
        let text = format!("{header}{chr1}{chr3a}{chr2}{chr3b}");
        let vcf = dir.join("pileup.vcf");
        fs::write(&vcf, &text).unwrap();
        let pileup = PileupVcf::index(&bgzip(&vcf).unwrap()).unwrap();

        assert_eq!(read_all(pileup.all().unwrap()), text);
        assert_eq!(read_all(pileup.contig("chr1").unwrap()), format!("{header}{chr1}"));
        assert_eq!(read_all(pileup.contig("chr2").unwrap()), format!("{header}{chr2}"));
        assert_eq!(read_all(pileup.contig("chr3").unwrap()), text);
        assert_eq!(read_all(pileup.contig("chrX").unwrap()), header);
        fs::remove_dir_all(&dir).unwrap();
    }
}

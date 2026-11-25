//! Spectre's `CNVCaller` end to end (spectreCNV.py:59-134): with SNVs, estimate the
//! copy-neutral depth; normalize, call, score, refine, score again; with SNVs, call LoH; write
//! the `.spc`, BED and VCF.

use std::path::{Path, PathBuf};

use crate::Error;
use crate::breakpoints;
use crate::call::{RawCalls, call, new_ids};
use crate::candidate::{Arena, CandidateId};
use crate::coverage::{ContigDepth, read_regions};
use crate::extensions::{Caller, Extensions, Prepared, finish_copy_numbers, prepare};
use crate::genome::Genome;
use crate::hmm;
use crate::loh::{LohCalls, LohOptions, call_loh};
use crate::merged::{Call, merge};
use crate::metadata::{Metadata, check_mdr, write_mdr};
use crate::normalize::{NormalizeOptions, Normalized, normalize};
use crate::output::{bed, spc, vcf};
use crate::refine::{ContigContext, RefineOptions, refine};
use crate::score::ZScore;
use crate::snv::{CnNeutral, SnvVcf, cn_neutral};
use crate::thresholds::{Clamping, Thresholds, WindowSample, seed_rng};

/// `CNVCaller`'s options.
#[derive(Clone, Debug)]
pub struct Options {
    /// `-c`: mosdepth `regions.bed.gz` (1 kb windows from `mosdepth -x -b 1000 -Q 20`); not read
    /// by [`call_cnvs_with_depths`].
    pub coverage: PathBuf,
    /// `-s`.
    pub sample_id: String,
    /// `-d`.
    pub output_dir: PathBuf,
    /// `-r`: the reference FASTA, with a `.fai` index.
    pub reference: PathBuf,
    /// `-m`: the `.mdr`; made from the reference if it does not exist. Default
    /// `<output_dir>/metadata.mdr`.
    pub metadata: Option<PathBuf>,
    /// `-l`.
    pub blacklist: Option<PathBuf>,
    /// `-v`: a single-sample SNV VCF (Clair3 or lrpipe-snv), for the copy-neutral depth and LoH.
    pub snv: Option<PathBuf>,
    /// Seed of numpy's generator for the dip test, which Spectre leaves unseeded (the
    /// reference runs seed it with 0). It only matters when heterozygous depths are not
    /// clearly unimodal.
    pub dip_seed: u32,
    pub loh: LohOptions,
    /// `-o`: contigs to call; default all kept contigs.
    pub only_chr: Option<Vec<String>>,
    /// `-p`.
    pub ploidy: i64,
    /// `-pc`.
    pub ploidy_chr: Vec<(String, i64)>,
    /// `-n`: N run length for `.mdr` creation.
    pub n_size: i64,
    /// `-mchl`.
    pub min_chr_len: u64,
    /// `-a`: no threshold clamping.
    pub clamping: Clamping,
    /// `-02`: written to the `.spc` only.
    pub mosdepth_cov_genome_chr_diff: f64,
    /// `-06`: written to the `.spc` only.
    pub dist_proportion: f64,
    pub normalize: NormalizeOptions,
    pub refine: RefineOptions,
    /// Changes to Spectre's method; none by default.
    pub extensions: Extensions,
    /// How the VCF is written; Spectre's way by default.
    pub vcf_style: vcf::Style,
}

/// `-mchl`'s default: contigs of this length or shorter are not called.
pub const DEFAULT_MIN_CHR_LEN: u64 = 1_000_000;

impl Options {
    /// The defaults for everything but the required paths and sample ID.
    pub fn new(coverage: PathBuf, sample_id: String, output_dir: PathBuf, reference: PathBuf) -> Self {
        Options {
            coverage,
            sample_id,
            output_dir,
            reference,
            metadata: None,
            blacklist: None,
            snv: None,
            dip_seed: 0,
            loh: LohOptions::default(),
            only_chr: None,
            ploidy: 2,
            ploidy_chr: Vec::new(),
            n_size: 5,
            min_chr_len: DEFAULT_MIN_CHR_LEN,
            clamping: Clamping::PloidyBounds,
            mosdepth_cov_genome_chr_diff: 0.10,
            dist_proportion: 0.25,
            normalize: NormalizeOptions::default(),
            refine: RefineOptions::default(),
            extensions: Extensions::default(),
            vcf_style: vcf::Style::Spectre,
        }
    }

    /// Checks the reference index without starting a call.
    ///
    /// # Errors
    /// If the reference or index cannot be read or the index does not fit the reference.
    pub fn check_reference(&self) -> Result<(), Error> {
        lrpipe_reads::check_fai(&self.reference).map_err(Error::from)
    }
}

/// Everything a run computed.
#[derive(Debug)]
pub struct Run {
    pub genome: Genome,
    /// What the extensions changed before normalization.
    pub prepared: Prepared,
    pub norm: Normalized,
    pub sample: WindowSample,
    pub thresholds: Thresholds,
    pub strict: Thresholds,
    pub arena: Arena,
    pub raw: RawCalls,
    pub calls: Vec<(String, Vec<CandidateId>)>,
    /// With SNVs: the copy-neutral depth and the LoH calls.
    pub cn_neutral: Option<CnNeutral>,
    pub loh: Option<LohCalls>,
    /// CNV and passing LoH calls together, as written to the BED and VCF.
    pub merged: Vec<(String, Vec<Call>)>,
    /// With SV breakpoints: how many calls took them.
    pub snapped: Option<usize>,
}

/// Runs the calling stages on the windows of `opts.coverage`.
///
/// # Errors
///
/// If `--min-cnv-len` is outside what Spectre accepts, or an input cannot be read or used.
pub fn call_cnvs(opts: &Options) -> Result<Run, Error> {
    call_cnvs_with_depths(opts, read_regions(&opts.coverage)?)
}

/// Runs the calling stages on window depths already in memory (as `read_regions` gives them,
/// such as from lrpipe-depth), ignoring `opts.coverage`.
///
/// # Errors
///
/// As [`call_cnvs`].
pub fn call_cnvs_with_depths(opts: &Options, depths: Vec<ContigDepth>) -> Result<Run, Error> {
    // analysis.py:86-101 (bug B1): Spectre exits, without a message, unless this holds.
    let mcl = opts.refine.min_cnv_len;
    if !(10_000 < mcl && mcl <= 100_000) {
        return Err(Error::Input(format!(
            "--min-cnv-len {mcl}: Spectre 0.2.1 only runs with values above 10000 and up to 100000"
        )));
    }
    opts.check_reference()?;
    let genome = Genome::from_fasta_index(&opts.reference, opts.min_chr_len, opts.ploidy, &opts.ploidy_chr)?;
    let mut only_chr = opts.only_chr.clone().unwrap_or_else(|| genome.contigs.iter().map(|c| c.name.clone()).collect());
    let mdr = opts.metadata.clone().unwrap_or_else(|| opts.output_dir.join("metadata.mdr"));
    // Metadata that is there is checked against the reference, whatever the output's style: one
    // made from another reference masks that reference's N runs. The check reads the whole
    // FASTA, so what it checked is noted beside the metadata and not checked again while
    // neither file changes.
    let checked = mdr.with_added_extension("checked");
    if !mdr.exists() {
        write_mdr(&opts.reference, opts.n_size, &mdr)?;
        note_checked(&checked, &opts.reference, opts.n_size, &mdr);
    } else if checked_key(&opts.reference, opts.n_size, &mdr)
        .is_none_or(|key| std::fs::read_to_string(&checked).ok().as_deref() != Some(key.as_str()))
    {
        check_mdr(&opts.reference, opts.n_size, &mdr)?;
        note_checked(&checked, &opts.reference, opts.n_size, &mdr);
    }
    let snv = opts.snv.as_deref().map(SnvVcf::read).transpose()?;
    let cn = snv
        .as_ref()
        .map(|vcf| cn_neutral(vcf, &genome, &mut pynum::numpy_random::RandomState::new(opts.dip_seed)))
        .transpose()?;
    let mut norm_opts = opts.normalize.clone();
    norm_opts.snv_cn_neutral = cn.as_ref().map(|c| c.median);
    let mut depths = depths;
    let prepared = prepare(&mut depths, &opts.extensions, &genome, &opts.reference)?;
    only_chr.retain(|c| !prepared.dropped.contains(c));
    let (norm, metadata) =
        normalize(&genome, depths, &only_chr, &norm_opts, |bin| Metadata::load(&mdr, opts.blacklist.as_deref(), bin))?;

    let mut rng = seed_rng(&opts.sample_id);
    let sample = WindowSample::draw(&norm, &opts.sample_id, &mut rng);
    let thresholds = Thresholds::from_quantiles(&sample, opts.ploidy, opts.clamping);
    let mut arena = Arena::default();
    let (raw, strict, mut calls) = match &opts.extensions.caller {
        Caller::Thresholds => {
            let raw = call(&norm, &thresholds, &mut rng, &mut arena);
            // `get_cnv_metrics()`: every raw candidate gets zero scores. (It also redraws the
            // window sample, which comes out the same: quirk Q-EXCL.)
            for (_, ids) in &raw.contigs {
                for &id in ids {
                    arena.get_mut(id).statistics = Some(ZScore::unrefined());
                }
            }
            // analysis.py:805-807: with no analysed contig the strict thresholds are never
            // computed and the quantile ones are recorded instead.
            let strict = if raw.contigs.is_empty() {
                thresholds.clone()
            } else {
                Thresholds::from_sd(&sample, 1.85, opts.ploidy, opts.clamping)
            };
            let mut calls = Vec::new();
            for (contig, ids) in &raw.contigs {
                let coverage =
                    norm.contigs.iter().find(|c| &c.name == contig).expect("raw calls come from analysed contigs");
                let ctx = ContigContext {
                    coverage,
                    metadata: &metadata,
                    bin_size: norm.bin_size,
                    thresholds: &thresholds,
                    strict: &strict,
                };
                calls.push((contig.clone(), refine(ids, &mut arena, &ctx, &opts.refine)?));
            }
            (raw, strict, calls)
        }
        Caller::Hmm(hmm_opts) => {
            let strict = Thresholds::from_sd(&sample, 1.85, opts.ploidy, opts.clamping);
            let mut raw = RawCalls::default();
            let min_len = opts.refine.min_cnv_len;
            for (contig, cands) in hmm::call(&norm, opts.ploidy, snv.as_ref(), min_len, &strict, hmm_opts) {
                let ids = new_ids(&mut rng, &mut raw.existing_ids, cands.len());
                let ids = cands
                    .into_iter()
                    .zip(ids)
                    .map(|(mut cand, id)| {
                        cand.id = format!("Spectre.{}.{id}", cand.cnv_type.as_str());
                        arena.push(cand)
                    })
                    .collect();
                raw.contigs.push((contig, ids));
            }
            let calls = raw.contigs.clone();
            (raw, strict, calls)
        }
    };
    for (_, ids) in &calls {
        for &id in ids {
            let z = ZScore::of(&arena.get(id).cov, &sample);
            arena.get_mut(id).statistics = Some(z);
        }
    }
    finish_copy_numbers(&mut arena, &mut calls, &opts.extensions, &prepared, genome.ploidy);
    let snapped = match &opts.extensions.sv_breakpoints {
        Some(path) => {
            Some(breakpoints::snap(&mut arena, &calls, &breakpoints::read(path, opts.refine.min_cnv_len / 2)?))
        }
        None => None,
    };
    // LoH IDs continue Python's `random` from the last redraw of the window sample.
    let loh = match (&snv, &cn) {
        (Some(vcf), Some(cn)) => {
            let mut rng = sample.redrawn_rng(&opts.sample_id);
            let mut existing = raw.existing_ids.clone();
            // A one-copy chromosome has no heterozygous SNVs to lose.
            let loh_chr: Vec<String> =
                only_chr.iter().filter(|c| !prepared.haploid.iter().any(|(h, _)| h == *c)).cloned().collect();
            Some(call_loh(vcf, &loh_chr, cn.median, &opts.loh, &mut rng, &mut existing))
        }
        _ => None,
    };
    let merged = merge(&arena, &calls, loh.as_ref());
    Ok(Run {
        genome,
        prepared,
        norm,
        sample,
        thresholds,
        strict,
        arena,
        raw,
        calls,
        cn_neutral: cn,
        loh,
        merged,
        snapped,
    })
}

/// The output paths Spectre uses under `dir`.
#[derive(Clone, Debug)]
pub struct OutputPaths {
    pub spc: PathBuf,
    pub bed: PathBuf,
    pub vcf: PathBuf,
}

impl OutputPaths {
    pub fn new(dir: &Path, sample_id: &str) -> Self {
        OutputPaths {
            spc: dir.join(format!("{sample_id}.spc.gz")),
            bed: dir.join(format!("{sample_id}_cnv.bed.gz")),
            vcf: dir.join(format!("{sample_id}.vcf.gz")),
        }
    }
}

/// What a check of `mdr` against `reference` holds for: both files' sizes and modification
/// times, the reference's path and the N-run size. `None` when either cannot be read.
fn checked_key(reference: &Path, n_size: i64, mdr: &Path) -> Option<String> {
    let stat = |path: &Path| -> Option<String> {
        let meta = std::fs::metadata(path).ok()?;
        let modified = meta.modified().ok()?.duration_since(std::time::UNIX_EPOCH).ok()?.as_nanos();
        Some(format!("{}\t{modified}", meta.len()))
    };
    let path = reference.canonicalize().ok()?;
    Some(format!("{}\t{}\t{n_size}\t{}\n", path.display(), stat(reference)?, stat(mdr)?))
}

/// Notes that `mdr` was checked against (or made from) `reference`. A note that cannot be
/// written (a metadata file shared read-only) only means the check is made again.
fn note_checked(note: &Path, reference: &Path, n_size: i64, mdr: &Path) {
    if let Some(key) = checked_key(reference, n_size, mdr) {
        let _ = std::fs::write(note, key);
    }
}

/// Writes the `.spc`, BED and VCF of a run.
///
/// # Errors
///
/// If a file cannot be written.
pub fn write_outputs(run: &Run, opts: &Options) -> Result<OutputPaths, Error> {
    std::fs::create_dir_all(&opts.output_dir).map_err(Error::io(&opts.output_dir))?;
    let paths = OutputPaths::new(&opts.output_dir, &opts.sample_id);
    let metrics = spc::Metrics {
        bin_size: run.norm.bin_size,
        min_cnv_len: opts.refine.min_cnv_len,
        normalization_value: run.norm.normalization_value,
        max_std_outlier_rm: opts.normalize.max_std_outlier_rm as i64,
        mosdepth_cov_genome_chr_diff: opts.mosdepth_cov_genome_chr_diff,
        cov_diff_threshold: opts.refine.cov_diff_threshold,
        dist_proportion: opts.dist_proportion,
        dist_min_overwrite: opts.refine.dist_min_overwrite as i64,
        thresholds: &run.thresholds,
        strict: &run.strict,
        sample: &run.sample,
    };
    let doc = spc::document(
        &opts.sample_id,
        &run.genome,
        &run.arena,
        &run.raw.contigs,
        &run.calls,
        run.loh.as_ref(),
        &metrics,
    );
    spc::write(&paths.spc, &doc)?;
    bed::write(&paths.bed, &run.arena, &run.merged)?;
    vcf::write(&paths.vcf, &opts.reference, &run.genome, &opts.sample_id, &run.arena, &run.merged, opts.vcf_style)?;
    Ok(paths)
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    #[test]
    fn a_foreign_fai_is_refused_before_creating_metadata() {
        let dir = std::env::temp_dir().join(format!("lrpipe-cnv-foreign-fai-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");
        let reference = dir.join("ref.fa");
        fs::write(&reference, ">chr1\nACGT\nACGT\n").expect("write reference");
        fs::write(reference.with_added_extension("fai"), "chr1\t8\t7\t4\t5\n").expect("write a foreign index");
        let opts = Options::new(dir.join("unused.bed"), "sample".to_owned(), dir.join("out"), reference);

        let err = call_cnvs_with_depths(&opts, Vec::new()).expect_err("reject foreign index");
        let Error::Reference(err) = err else {
            panic!("expected an invalid reference error, got {err}");
        };
        assert!(err.into_io_error().to_string().contains("its .fai is not this FASTA's index"));
        assert!(!opts.output_dir.exists());
        fs::remove_dir_all(dir).expect("remove scratch directory");
    }
}

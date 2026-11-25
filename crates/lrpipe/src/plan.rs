//! The stages of a sample's run: which program each runs, with what arguments, reading and
//! writing which files, after which other stages.
//!
//! lrpipe-snv (SNVs and indels) and `lrpipe sv`'s single read pass (SVs, depth, tandem repeats
//! and bedMethyl methylation) read the alignments at the same time; lrpipe-cnv calls CNVs once
//! both have finished, from the pass's depth windows, lrpipe-snv's heterozygous SNVs and the SVs'
//! breakpoints. pb-CpG-tools' methylation tables (lrpipe-methyl's `cpg`) are a stage of their
//! own, reading the alignments alongside the others. Then `lrpipe merge` writes the variant VCFs
//! as one, and last `lrpipe qc` summarises the run from its files ([`crate::qc`]). [`Analyses`]
//! chooses what runs: each analysis and the stages it needs.
//!
//! When it phases the SNVs, the run haplotags the reads from them (`lrpipe haplotag`, into a
//! haplotag list rather than a second BAM), for the SVs' phases, the tandem repeats' alleles and
//! the CpG tables per haplotype. The read pass still runs beside lrpipe-snv: it calls the SVs
//! unphased and writes what phasing them needs, and `lrpipe sv phase` phases them from the list.
//! The tandem repeats (`lrpipe tr`, then a pass of their own) and the CpG tables wait for the
//! list.
//!
//! Reads that are not one sorted, indexed BAM or CRAM are made one first, by a stage every stage
//! reading them follows: several BAMs or CRAMs of the sample ([`Merged`]) are merged by
//! `merge-reads` into one BAM ([`Layout::merged_reads`]), and unaligned reads (dorado's or
//! PacBio's unaligned BAM, FASTQ) mapped by bonefish, in an align stage, into one
//! ([`Layout::aligned_bam`]). The stages then read that BAM as they read aligned input.

use std::ffi::OsString;
use std::fmt;
use std::fs;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use crate::config::{Estimate, MethylationReport, Site};
use crate::provenance::{self, hex};
use lrpipe_reads::regions::Regions;

/// The pipeline's stages, in the order they are listed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StageName {
    /// Unaligned reads mapped into a sorted, indexed BAM (bonefish).
    Align,
    /// Several BAMs or CRAMs merged into one sorted, indexed BAM (`lrpipe merge-reads`), which
    /// every other stage reads.
    MergeReads,
    /// SNVs and indels (lrpipe-snv).
    Snv,
    /// SVs, depth, tandem repeats and bedMethyl methylation in one pass over the reads
    /// (`lrpipe sv`).
    Reads,
    /// The reads' haplotypes from the phased SNVs, as a haplotag list (`lrpipe haplotag`).
    Haplotag,
    /// Tandem repeats genotyped with the haplotag list (`lrpipe tr`), when the run phases; else
    /// the read pass genotypes them.
    Repeats,
    /// The SVs phased from the haplotag list (`lrpipe sv phase`).
    SvPhase,
    /// Copy number variants (lrpipe-cnv).
    Cnv,
    /// pb-CpG-tools' methylation tables (lrpipe-methyl); bedMethyl comes from the read pass.
    Methylation,
    /// The variant VCFs as one (`lrpipe merge`).
    Merge,
    /// The run's QC summary from the other stages' outputs (`lrpipe qc`).
    Qc,
}

impl StageName {
    /// Every stage, in the order they are listed: what a stage's name is read back from.
    pub const ALL: [StageName; 11] = [
        Self::Align,
        Self::MergeReads,
        Self::Snv,
        Self::Reads,
        Self::Haplotag,
        Self::Repeats,
        Self::SvPhase,
        Self::Cnv,
        Self::Methylation,
        Self::Merge,
        Self::Qc,
    ];

    /// The stage named `name` ([`Self::as_str`]).
    #[must_use]
    pub fn named(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|s| s.as_str() == name)
    }

    /// Whether the stage failing fails the run: every stage's does but the QC summary's, which
    /// the run's results do not need (its failure is recorded and warned of).
    #[must_use]
    pub fn fails_the_run(self) -> bool {
        self != Self::Qc
    }

    /// The stage's name in logs, the manifest and file names.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Align => "align",
            Self::MergeReads => "merge-reads",
            Self::Snv => "snv",
            Self::Reads => "reads",
            Self::Haplotag => "haplotag",
            Self::SvPhase => "sv-phase",
            Self::Repeats => "repeats",
            Self::Cnv => "cnv",
            Self::Methylation => "methylation",
            Self::Merge => "merge",
            Self::Qc => "qc",
        }
    }
}

impl fmt::Display for StageName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// What a run can produce.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Analysis {
    /// SNVs and indels (lrpipe-snv).
    Snv,
    /// The SNVs and indels phased (Clair3's `--enable_phasing`: lrpipe-phase, as WhatsHap);
    /// on HG002 HiFi it adds about 100 s to lrpipe-snv.
    Phase,
    /// Structural variants (the read pass).
    Sv,
    /// Copy number variants and loss of heterozygosity (lrpipe-cnv). It needs the SNVs, and the
    /// read pass's depth and SV breakpoints, so those run too.
    Cnv,
    /// Tandem repeat genotypes (the read pass), with a repeat catalogue configured.
    Tr,
    /// Read depth in 1 kb windows (the read pass).
    Depth,
    /// CpG methylation: bedMethyl from the read pass, or pb-CpG-tools' tables
    /// ([`MethylationReport`]).
    Methyl,
}

impl Analysis {
    /// Every analysis, in the order they are listed.
    pub const ALL: [Analysis; 7] =
        [Analysis::Snv, Analysis::Phase, Analysis::Sv, Analysis::Cnv, Analysis::Tr, Analysis::Depth, Analysis::Methyl];

    /// The name `--analyses` takes, as the lrpipe subcommand.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Snv => "snv",
            Self::Phase => "phase",
            Self::Sv => "sv",
            Self::Cnv => "cnv",
            Self::Tr => "tr",
            Self::Depth => "depth",
            Self::Methyl => "methyl",
        }
    }

    fn bit(self) -> u8 {
        1 << self as u8
    }
}

/// A set of [`Analysis`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Analyses(u8);

impl Analyses {
    /// Every analysis.
    #[must_use]
    pub fn all() -> Self {
        Self::of(&Analysis::ALL)
    }

    /// The listed analyses.
    #[must_use]
    pub fn of(analyses: &[Analysis]) -> Self {
        Self(analyses.iter().fold(0, |set, a| set | a.bit()))
    }

    /// Whether `analysis` is chosen.
    #[must_use]
    pub fn contains(self, analysis: Analysis) -> bool {
        self.0 & analysis.bit() != 0
    }

    /// Whether lrpipe-snv runs: for the SNVs, their phasing, or the CNVs, which are called from
    /// them.
    #[must_use]
    pub fn calls_snvs(self) -> bool {
        [Analysis::Snv, Analysis::Phase, Analysis::Cnv].into_iter().any(|a| self.contains(a))
    }
}

/// Whether the run writes the merged VCF.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MergedVcf {
    /// `SAMPLE.vcf.gz`: the chosen variant analyses' VCFs as one.
    #[default]
    On,
    /// Only each tool's own files.
    Off,
}

/// Whether lrpipe-snv also writes a gVCF, for joint genotyping (GLnexus).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Gvcf {
    /// Only the VCFs.
    #[default]
    Off,
    /// `snv/merge_output.gvcf.gz` too (Clair3's `--gvcf`), with its index.
    On,
}

/// Whether the reads carry haplotypes, as the `HP` tags of haplotagging, which the SVs' phases
/// and the CpG tables per haplotype are read from when the run does not phase the SNVs. When it
/// does, it haplotags the reads from them itself, whatever tags they carry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Haplotypes {
    /// The reads carry `HP` tags.
    Tagged,
    /// They do not: the SVs are unphased, and no CpG tables per haplotype are written (they would
    /// be empty).
    Untagged,
}

/// Whether a run that phases its SNVs haplotags the reads from them (`--haplotag`), for the SVs'
/// phases, the tandem repeats' alleles and the CpG tables per haplotype.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Haplotag {
    /// It does: those follow the run's phasing.
    #[default]
    On,
    /// It does not: they follow the reads' own `HP` tags, if any, and run beside lrpipe-snv
    /// rather than after its phasing.
    Off,
}

/// How the read pass calls SVs (`--sv-mode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SvMode {
    /// lrpipe-sv's `--accuracy`: a repeat expansion's two alleles called apart, and a read's
    /// indels joined only within 50 bp. On HG002 against GIAB v5.0q it finds about 1,600 more
    /// SVs on HiFi and on ONT, F1 0.835 and 0.836 where Sniffles2's calls score 0.807 and 0.810.
    #[default]
    Accuracy,
    /// Sniffles2 2.8.1's calls, record for record.
    Sniffles2,
}

/// Whether the run writes the reads haplotagged, as a BAM (`phase/SAMPLE.haplotagged.bam`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HaplotaggedReads {
    /// Only the haplotag list.
    #[default]
    Off,
    /// The BAM too, with its index, for IGV and other tools (`--haplotagged-reads`).
    On,
}

/// Sequencing platform.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Platform {
    /// PacBio HiFi.
    Hifi,
    /// Oxford Nanopore.
    Ont,
}

impl Platform {
    /// The name every tool takes for it.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Hifi => "hifi",
            Self::Ont => "ont",
        }
    }
}

/// The sample's sex chromosomes, or a haploid bacterial genome.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    /// XX.
    Female,
    /// XY: one copy of chrX outside the pseudo-autosomal regions, and chrY.
    Male,
    /// Haploid across the reference, without human sex chromosomes.
    NotApplicable,
}

impl Sex {
    /// The value recorded by the run and used by QC.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Female => "female",
            Self::Male => "male",
            Self::NotApplicable => "not-applicable",
        }
    }

    fn karyotype(self) -> &'static str {
        match self {
            Self::Female => "xx",
            Self::Male => "xy",
            Self::NotApplicable => unreachable!("only human analyses require a karyotype"),
        }
    }
}

/// The sample to run.
#[derive(Clone, Debug)]
pub struct Sample {
    /// Name in every output's sample column and file names.
    pub name: String,
    /// The sorted, indexed BAM or CRAM the stages read: the input; for several aligned files the
    /// merge-reads stage's [`Layout::merged_reads`]; for unaligned input the align stage's
    /// [`Layout::aligned_bam`]. The stages naming it as an input follow the stage writing it.
    pub reads: PathBuf,
    /// Unaligned reads (unaligned BAM, FASTQ or FASTQ.gz), which the align stage maps into
    /// `reads`, one file after another; none for aligned input.
    pub unaligned: Vec<PathBuf>,
    /// Whether `unaligned` is FASTQ, which has no read groups of its own.
    pub unaligned_fastq: bool,
    /// Several aligned BAMs and CRAMs, which the merge-reads stage merges into `reads`; none for
    /// one file or for unaligned reads.
    pub merged: Option<Merged>,
    /// Platform the reads come from.
    pub platform: Platform,
    /// Sex chromosomes for human calling, or no sex chromosomes for haploid bacterial calling.
    pub sex: Sex,
    /// The Clair3 model directory lrpipe-snv calls the reads with ([`crate::snv_model::choose`]);
    /// without one, the site's for the platform.
    pub model: Option<PathBuf>,
    /// The regions the analyses are restricted to (`lrpipe run --regions`); none for the whole
    /// genome.
    pub regions: Option<Targets>,
}

/// The regions a run's analyses are restricted to, which [`Layout::regions_bed`] holds (sorted
/// in the reference's order and merged), with the depth's windows in them in
/// [`Layout::region_windows`]. Each tool is restricted by its own region option, with its own
/// semantics at the regions' edges; the align and merge-reads stages, which make the reads the
/// others read, take all of them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Targets {
    /// Every region is a whole contig: these, in the reference's order. The tools that can be
    /// restricted to contigs are (lrpipe-snv's `--ctg-name`, the read pass's `--contig`). The
    /// depth, SVs and repeats' loci are then a run of the whole genome's records of those contigs;
    /// lrpipe-snv's are not quite, as it takes its quality cut-offs from every contig it calls
    /// (Clair3's), and nor is what follows from them (the phasing, haplotags and haplotypes).
    Contigs(Vec<String>),
    /// Regions that are not all whole contigs, restricted to by the BED.
    Bed,
}

impl Targets {
    /// The targets of `regions` (as [`crate::checks::regions`] gives them) on a reference of
    /// `contigs` (each name and length).
    #[must_use]
    pub fn of(regions: &Regions, contigs: &[(String, u64)]) -> Self {
        let whole = regions.iter().all(|(contig, intervals)| {
            let length = contigs.iter().find(|(c, _)| c == contig).map(|&(_, l)| l);
            intervals.len() == 1 && length.is_some_and(|l| intervals[0] == (0, l as i64))
        });
        if whole { Self::Contigs(regions.iter().map(|(c, _)| c.to_owned()).collect()) } else { Self::Bed }
    }
}

/// Whether a run of the whole genome calls anything on contig `name` of `length` bases: SNVs on
/// lrpipe-snv's contigs (Clair3's: chr1-22, X and Y, or 1-22, X and Y), SVs on those of 1 Mbp or
/// more (Sniffles2's). The alt, decoy, HLA, EBV and mitochondrial contigs of a human reference
/// are neither.
#[must_use]
pub fn called(name: &str, length: u64) -> bool {
    lrpipe_snv::is_default_contig(name) || length >= lrpipe_sv::cli::MIN_CONTIG_LEN.unsigned_abs()
}

/// `--regions` as a run takes them ([`restrict`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Restriction {
    /// The regions on contigs a run of the whole genome calls ([`called`]).
    pub regions: Regions,
    /// How the stages are restricted to them; none when they cover every contig a run of the
    /// whole genome calls, whole: the run is the whole genome's.
    pub targets: Option<Targets>,
    /// The contigs of regions left out, which a run of the whole genome does not call.
    pub left_out: Vec<String>,
}

/// `regions` (as [`crate::checks::regions`] gives them) on a reference of `contigs`, as a run
/// takes them: without those on contigs a run of the whole genome calls nothing on (a BED of
/// every contig of the `.fai` names the alt and decoy contigs too), and as no restriction when
/// they cover every contig it calls (rather than a list of them all).
#[must_use]
pub fn restrict(regions: &Regions, contigs: &[(String, u64)]) -> Restriction {
    restrict_contigs(regions, contigs, called)
}

/// Restrict to regions on contigs eligible for the selected calling mode.
pub(crate) fn restrict_contigs(
    regions: &Regions,
    contigs: &[(String, u64)],
    called: impl Fn(&str, u64) -> bool,
) -> Restriction {
    let length = |name: &str| contigs.iter().find(|(c, _)| c == name).map_or(0, |&(_, l)| l);
    type Contig<'a> = (&'a str, &'a [(i64, i64)]);
    let (kept, left_out): (Vec<Contig<'_>>, Vec<Contig<'_>>) = regions.iter().partition(|&(c, _)| called(c, length(c)));
    let kept = Regions::from_intervals(
        kept.iter().flat_map(|&(c, intervals)| intervals.iter().map(move |&(s, e)| (c.to_owned(), s, e))),
    );
    let targets = match Targets::of(&kept, contigs) {
        Targets::Contigs(named) if contigs.iter().filter(|&(c, l)| called(c, *l)).all(|(c, _)| named.contains(c)) => {
            None
        }
        targets => Some(targets),
    };
    Restriction { regions: kept, targets, left_out: left_out.iter().map(|&(c, _)| c.to_owned()).collect() }
}

/// The window the depth is reported in, and the CNVs called from (Spectre's 1 kb).
pub const DEPTH_WINDOW: u32 = 1000;

/// The depth's windows in `regions`, as a BED ([`Layout::region_windows`]): those of a run of
/// the whole genome (every [`DEPTH_WINDOW`] bases from each contig's start) that a region reaches,
/// cut to it. A window wholly in a region is the whole genome's, and its depth too.
#[must_use]
pub fn windows_in(regions: &Regions) -> String {
    let size = i64::from(DEPTH_WINDOW);
    let mut out = String::new();
    for (contig, intervals) in regions.iter() {
        for &(start, end) in intervals {
            let mut window = start - start % size;
            while window < end {
                out.push_str(&format!("{contig}\t{}\t{}\n", window.max(start), (window + size).min(end)));
                window += size;
            }
        }
    }
    out
}

/// Several aligned BAMs and CRAMs of a sample (with the same reference sequences), which the
/// `merge-reads` stage merges into one BAM, [`Layout::merged_reads`] (classified by
/// [`crate::input::classify`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Merged {
    /// The files, in the order given.
    pub inputs: Vec<PathBuf>,
    /// The merged BAM's index.
    pub index: BamIndex,
}

/// A BAM index: a BAI, or a CSI where a reference sequence is too long for a BAI (over 2^29
/// bases), as `samtools index` chooses ([`crate::merge_reads::bam_index_for`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BamIndex {
    /// `.bai`.
    Bai,
    /// `.csi`.
    Csi,
}

impl BamIndex {
    /// The index's file extension.
    #[must_use]
    pub fn extension(self) -> &'static str {
        match self {
            Self::Bai => "bai",
            Self::Csi => "csi",
        }
    }
}

/// Whether `reads` is CRAM (else BAM), by its name, as every stage takes it.
#[must_use]
pub fn is_cram(reads: &Path) -> bool {
    reads.as_os_str().as_encoded_bytes().ends_with(b".cram")
}

/// The index beside `reads`, as lrpipe's tools (pasta) look for it: `x.bam.bai` or `x.bam.csi`
/// for BAM, `x.cram.crai` for CRAM. (htslib also takes `x.bai`, which the tools do not.)
#[must_use]
pub fn index_of(reads: &Path) -> Option<PathBuf> {
    let extensions: &[&str] = if is_cram(reads) { &["crai"] } else { &["bai", "csi"] };
    extensions.iter().map(|e| Layout::index(reads, e)).find(|p| p.exists())
}

/// Whether lrpipe-snv's intermediate files stay after it finishes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Intermediates {
    /// Delete them.
    #[default]
    Remove,
    /// Keep them.
    Keep,
}

/// Threads given to the two stages that read the alignments at the same time, and to the align
/// stage, which runs alone before them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Threads {
    /// lrpipe-snv's.
    pub snv: usize,
    /// `lrpipe sv`'s.
    pub reads: usize,
    /// bonefish's: every core.
    pub align: usize,
}

impl Threads {
    /// `n` threads for each, on `n` cores they share. lrpipe-snv keeps fewer than that busy
    /// while it waits on the GPU, and runs alone with all of them once the shorter read pass is
    /// done; on HG002 HiFi with 16 cores this finished in 544 s, against 568 s with 8 threads
    /// each.
    #[must_use]
    pub fn each(n: usize) -> Self {
        let n = n.max(1);
        Self { snv: n, reads: n, align: n }
    }

    /// `n` threads between them, for when they cannot share `n` cores: lrpipe-snv takes the
    /// larger half (each has at least one).
    #[must_use]
    pub fn split(n: usize) -> Self {
        Self { snv: n.max(1) - n / 2, reads: (n / 2).max(1), align: n.max(1) }
    }
}

/// Run-wide choices besides the sample.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Thread counts.
    pub threads: Threads,
    /// What happens to lrpipe-snv's intermediate files.
    pub intermediates: Intermediates,
    /// What the run produces.
    pub analyses: Analyses,
    /// Whether it writes the merged VCF.
    pub merged_vcf: MergedVcf,
    /// Whether lrpipe-snv writes a gVCF.
    pub gvcf: Gvcf,
    /// Whether the reads are haplotagged.
    pub haplotypes: Haplotypes,
    /// Whether the run haplotags the reads from its phasing.
    pub haplotag: Haplotag,
    /// Whether the run writes them haplotagged from its phasing.
    pub haplotagged_reads: HaplotaggedReads,
    /// How the read pass calls SVs.
    pub sv_mode: SvMode,
    /// Whether the reference's pseudo-autosomal regions are known.
    pub pars: Pars,
    /// GiB the stages running together may use, which the align stage's settings fit
    /// ([`align_stage`]).
    pub memory_gib: f64,
}

/// The reference's pseudo-autosomal regions (PARs), which calling an XY sample's chrX and chrY
/// haploid outside them needs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Pars {
    /// Known (GRCh38, GRCh37 or CHM13), written to [`Layout::par_bed`] before the run.
    Known,
    /// Another reference: SNVs and indels on chrX and chrY are called diploid in XY samples, and
    /// tandem repeats and copy number count all of chrX as one copy (lrpipe#44).
    #[default]
    Unknown,
}

/// Where a run's files go, under one output directory.
#[derive(Clone, Debug)]
pub struct Layout {
    root: PathBuf,
    sample: String,
}

impl Layout {
    /// The layout of `sample`'s run in `root`.
    pub fn new(root: impl Into<PathBuf>, sample: &str) -> Self {
        Self { root: root.into(), sample: sample.to_owned() }
    }

    /// The output directory.
    #[must_use]
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// The sample.
    #[must_use]
    pub fn sample(&self) -> &str {
        &self.sample
    }

    /// The sample whose runs the output directory holds, for refusing another's.
    #[must_use]
    pub fn owner(&self) -> PathBuf {
        self.root.join(".lrpipe").join("sample")
    }

    fn file(&self, dir: &str, suffix: &str) -> PathBuf {
        self.root.join(dir).join(format!("{}{suffix}", self.sample))
    }

    /// The reads' haplotypes from the run's phasing: a haplotag list (WhatsHap's
    /// `--output-haplotag-list`: read name, haplotype, phase set, chromosome).
    #[must_use]
    pub fn haplotags(&self) -> PathBuf {
        self.file("phase", ".haplotags.tsv.gz")
    }

    /// The reads haplotagged from the run's phasing (`--haplotagged-reads`).
    #[must_use]
    pub fn haplotagged_reads(&self) -> PathBuf {
        self.file("phase", ".haplotagged.bam")
    }

    /// The read pass's SVs before they are phased from the haplotag list.
    #[must_use]
    pub fn unphased_sv_vcf(&self) -> PathBuf {
        self.file("sv", ".unphased.sv.vcf")
    }

    /// The read pass's SV candidates before they are phased.
    #[must_use]
    pub fn unphased_snf(&self) -> PathBuf {
        self.file("sv", ".unphased.snf")
    }

    /// What phasing the read pass's SVs needs (`lrpipe sv --phasing-evidence`).
    #[must_use]
    pub fn sv_phasing_evidence(&self) -> PathBuf {
        self.file("sv", ".phasing-evidence.gz")
    }

    /// The sample's several BAMs and CRAMs merged into one BAM ([`Merged`]).
    #[must_use]
    pub fn merged_reads(&self) -> PathBuf {
        self.file("reads", ".bam")
    }

    /// The align stage's directory.
    #[must_use]
    pub fn align_dir(&self) -> PathBuf {
        self.root.join("align")
    }

    /// The unaligned reads mapped by bonefish: a coordinate-sorted BAM (bonefish writes no
    /// CRAM), with its index beside it.
    #[must_use]
    pub fn aligned_bam(&self) -> PathBuf {
        self.file("align", ".bam")
    }

    /// The align stage's scratch directory, for the sorter's runs; removed once it succeeds.
    #[must_use]
    pub fn align_scratch(&self) -> PathBuf {
        self.align_dir().join("tmp")
    }

    /// lrpipe-snv's output directory.
    #[must_use]
    pub fn snv_dir(&self) -> PathBuf {
        self.root.join("snv")
    }

    /// SNVs and indels.
    #[must_use]
    pub fn snv_vcf(&self) -> PathBuf {
        self.snv_dir().join("merge_output.vcf.gz")
    }

    /// SNVs and indels, phased.
    #[must_use]
    pub fn phased_snv_vcf(&self) -> PathBuf {
        self.snv_dir().join("phased_merge_output.vcf.gz")
    }

    /// SNVs and indels with the reference blocks between them, for joint genotyping.
    #[must_use]
    pub fn snv_gvcf(&self) -> PathBuf {
        self.snv_dir().join("merge_output.gvcf.gz")
    }

    /// Every chosen variant analysis' calls in one VCF.
    #[must_use]
    pub fn merged_vcf(&self) -> PathBuf {
        self.root.join(format!("{}.vcf.gz", self.sample))
    }

    /// Structural variants.
    #[must_use]
    pub fn sv_vcf(&self) -> PathBuf {
        self.file("sv", ".sv.vcf")
    }

    /// The SV candidates, for multi-sample calling (`lrpipe sv -i a.snf b.snf ...`).
    #[must_use]
    pub fn snf(&self) -> PathBuf {
        self.file("sv", ".snf")
    }

    /// Prefix of the depth files.
    #[must_use]
    pub fn depth_prefix(&self) -> PathBuf {
        self.file("depth", "")
    }

    /// Depth in 1 kb windows (mosdepth `-x -Q 20 -b 1000`).
    #[must_use]
    pub fn depth_windows(&self) -> PathBuf {
        self.file("depth", ".regions.bed.gz")
    }

    /// Depth summary per contig.
    #[must_use]
    pub fn depth_summary(&self) -> PathBuf {
        self.file("depth", ".mosdepth.summary.txt")
    }

    /// The reads' statistics, counted by the read pass (`lrpipe sv --read-stats`).
    #[must_use]
    pub fn read_stats(&self) -> PathBuf {
        self.file("qc", ".read-stats.json")
    }

    /// The QC summary, as JSON.
    #[must_use]
    pub fn qc_json(&self) -> PathBuf {
        self.file("qc", ".qc.json")
    }

    /// The QC summary, as a table of its values.
    #[must_use]
    pub fn qc_tsv(&self) -> PathBuf {
        self.file("qc", ".qc.tsv")
    }

    /// Tandem repeat genotypes.
    #[must_use]
    pub fn repeats_vcf(&self) -> PathBuf {
        self.file("repeats", ".repeats.vcf.gz")
    }

    /// Methylation as modkit's bedMethyl.
    #[must_use]
    pub fn bedmethyl(&self) -> PathBuf {
        self.file("methylation", ".bedmethyl.gz")
    }

    /// Prefix of the CpG tables.
    #[must_use]
    pub fn cpg_prefix(&self) -> PathBuf {
        self.file("methylation", ".cpg")
    }

    /// The CpG table of `reads`: `combined`, `hap1` or `hap2`.
    #[must_use]
    pub fn cpg_table(&self, reads: &str) -> PathBuf {
        self.file("methylation", &format!(".cpg.{reads}.bed.gz"))
    }

    /// lrpipe-cnv's output directory.
    #[must_use]
    pub fn cnv_dir(&self) -> PathBuf {
        self.root.join("cnv")
    }

    /// Copy number variants.
    #[must_use]
    pub fn cnv_vcf(&self) -> PathBuf {
        self.file("cnv", ".vcf.gz")
    }

    /// Copy number segments.
    #[must_use]
    pub fn cnv_bed(&self) -> PathBuf {
        self.file("cnv", "_cnv.bed.gz")
    }

    /// The `.mdr` of N runs that lrpipe-cnv makes from the reference whose [`reference_key`] is
    /// `key`, when the configuration names none.
    #[must_use]
    pub fn cnv_metadata(&self, key: &str) -> PathBuf {
        self.cnv_dir().join(format!("reference-{key}.mdr"))
    }

    /// One log per stage.
    #[must_use]
    pub fn log(&self, stage: StageName) -> PathBuf {
        self.root.join("logs").join(format!("{stage}.log"))
    }

    /// `path`'s index: `path` with `suffix` (`tbi`, `csi`) appended.
    #[must_use]
    pub fn index(path: &Path, suffix: &str) -> PathBuf {
        let mut name = path.as_os_str().to_owned();
        name.push(".");
        name.push(suffix);
        PathBuf::from(name)
    }

    /// The regions of `lrpipe run --regions`, as a BED sorted in the reference's order and merged
    /// ([`Targets`]).
    #[must_use]
    pub fn regions_bed(&self) -> PathBuf {
        self.root.join(".lrpipe").join("regions.bed")
    }

    /// The depth's windows in the regions ([`Targets::Bed`]): the [`DEPTH_WINDOW`] windows of a
    /// run of the whole genome, each cut to the regions, as a BED for the read pass's
    /// `--depth-by`.
    #[must_use]
    pub fn region_windows(&self) -> PathBuf {
        self.root.join(".lrpipe").join("regions.windows.bed")
    }

    /// The reference's PARs as a BED, for lrpipe-snv's `--par-regions-bed`.
    #[must_use]
    pub fn par_bed(&self) -> PathBuf {
        self.root.join(".lrpipe").join("par.bed")
    }

    /// Record of a stage's last successful run, for resuming.
    #[must_use]
    pub fn stamp(&self, stage: StageName) -> PathBuf {
        self.root.join(".lrpipe").join(format!("{stage}.done"))
    }

    /// Where a running stage writes how far it has got, for the run's display
    /// ([`lrpipe_reads::progress`]).
    #[must_use]
    pub fn progress(&self, stage: StageName) -> PathBuf {
        self.root.join(".lrpipe").join("progress").join(stage.as_str())
    }

    /// The run's manifest.
    #[must_use]
    pub fn manifest(&self) -> PathBuf {
        self.root.join("manifest.json")
    }
}

/// One program run.
#[derive(Clone, Debug, PartialEq)]
pub struct Stage {
    /// Which stage.
    pub name: StageName,
    /// The program.
    pub program: PathBuf,
    /// Its arguments.
    pub args: Vec<OsString>,
    /// Files it reads; a change to any reruns it on resume.
    pub inputs: Vec<PathBuf>,
    /// Files it must write.
    pub outputs: Vec<PathBuf>,
    /// Stages that must finish first.
    pub after: Vec<StageName>,
    /// Expected peak memory in GiB, for keeping stages that run together within the
    /// machine's memory.
    pub memory_gib: f64,
    /// Whether the run waits on it ([`Priority::Critical`]: lrpipe-snv, then the CNVs) or it
    /// fills the CPU the critical stages leave idle ([`Priority::Background`]).
    pub priority: Priority,
    /// A directory of its own scratch files: made empty before it starts, and removed once it
    /// succeeds.
    pub scratch: Option<PathBuf>,
}

/// How a stage shares the CPU with the others.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Priority {
    /// On the run's longest path: the normal CPU priority.
    Critical,
    /// Off it: a lower CPU priority ([`BACKGROUND_NICE`]), so it takes the CPU the critical
    /// stages leave idle (lrpipe-snv waiting on the GPU, say) without slowing them.
    Background,
}

/// The `nice` value of [`Priority::Background`] stages.
pub const BACKGROUND_NICE: i32 = 10;

/// Peak memory in GiB with 16 threads, rounded up with room for run-to-run variation (lrpipe-snv's
/// peak moves by 2 GiB between runs, as phasing and full alignment overlap differently). Measured
/// 2026-11-02 on whole genomes: lrpipe-snv 8.9-11.7 GiB on HG002 HiFi and 10.9 on the HG002 ONT
/// CRAM, and 13.4-13.5 on the HG003 ONT BAM once it phased the final calls of each contig beside
/// its calling; the read pass 6.3-6.9 GiB on ONT with bedMethyl and 3.7 on
/// HiFi; lrpipe-cnv 1.0; methylation 1.0 on HiFi (ONT's is an estimate). lrpipe-snv's grows with
/// its threads.
pub(crate) fn expected_memory(name: StageName, platform: Platform) -> f64 {
    match (name, platform) {
        // Set by the settings it runs with (ALIGN_SETTINGS).
        (StageName::Align, _) => ALIGN_SETTINGS[0].peak_gib(16, false),
        // It streams its inputs, a few MiB of each in flight and the BGZF blocks being
        // compressed: 0.13 GiB merging two BAMs of HG003 ONT chr20 (2.5 GB). Each CRAM input
        // adds [`CRAM_INPUT_MEMORY`].
        (StageName::MergeReads, _) => 0.5,
        (StageName::Snv, Platform::Hifi) => 13.0,
        (StageName::Snv, Platform::Ont) => 14.5,
        (StageName::Reads, Platform::Hifi) => 5.0,
        (StageName::Reads, Platform::Ont) => 6.5,
        // The list alone, whole-genome HG002 with 16 threads (2026-11-03): 2.8 GiB from the HiFi
        // BAM and 5.6 from the ONT CRAM ([`HAPLOTAGGED_READS_MEMORY`] more with the BAM).
        (StageName::Haplotag, Platform::Hifi) => 3.5,
        (StageName::Haplotag, Platform::Ont) => 6.5,
        // The evidence and the read pass's VCF and SNF besides the list: on HG002 chr20 (2.1% of
        // the genome) 16 MB (ONT) and 13 MB (HiFi), so some 0.7 GiB for the genome.
        (StageName::SvPhase, _) => 0.75 + HAPLOTAG_LIST_MEMORY,
        // `lrpipe tr --haplotags` on whole-genome HG002, 16 threads (2026-11-03), the list
        // included: 3.9 GiB from the HiFi BAM, 6.7 from the ONT CRAM.
        (StageName::Repeats, Platform::Hifi) => 4.5,
        (StageName::Repeats, Platform::Ont) => 7.5,
        (StageName::Cnv, _) => 1.5,
        (StageName::Methylation, Platform::Hifi) => 2.5,
        (StageName::Methylation, Platform::Ont) => 2.5,
        // It streams its inputs.
        (StageName::Merge, _) => 0.5,
        // It streams the VCFs and tables, but for the phased SNVs, which it reads a contig at a
        // time.
        (StageName::Qc, _) => 1.0,
    }
}

/// What each CRAM input adds to the merge-reads stage's peak memory, in GiB: its containers being
/// decoded, and the reference it is decoded against. That is the bases of its slices' spans,
/// or a whole contig for a slice spanning half of one (a short contig), only the contig being
/// read kept.
const CRAM_INPUT_MEMORY: f64 = 0.25;

/// Settings of bonefish that bound its memory, and its peak with them.
#[derive(Clone, Copy, Debug)]
struct AlignSettings {
    /// Its arguments.
    args: &'static [&'static str],
    /// The most threads it runs with, if fewer than the run's.
    max_threads: Option<usize>,
    /// Peak memory in GiB mapping a human genome with `threads` threads.
    peak_gib: f64,
    /// The threads `peak_gib` was measured with.
    threads: usize,
}

impl AlignSettings {
    /// The threads it runs with, of the run's `threads`.
    fn threads(&self, threads: usize) -> usize {
        self.max_threads.map_or(threads, |most| threads.min(most)).max(1)
    }

    /// Its expected peak with the run's `threads`, mapping to an index built from the reference
    /// first (`built`) or to one built before.
    fn peak_gib(&self, threads: usize, built: bool) -> f64 {
        let threads = self.threads(threads) as f64;
        let mapping = self.peak_gib + ALIGN_THREAD_GIB * (threads - self.threads as f64);
        if built { mapping.max(INDEX_BUILD_GIB) } else { mapping }
    }
}

/// bonefish's settings, from the most memory to the least, of which the align stage takes the
/// first whose peak fits the run's memory, else the last. bonefish 0.5.0's documentation
/// (docs/usage.md) gives its peaks with 16 threads mapping a human genome to a sorted BAM with
/// `-y` carrying methylation tags: about 23 GB with its defaults, and about 19 GB with
/// `--sort-mem 1G` (the sorter's batches, of which two are alive at a time, at 1 GB in place of
/// 4). With `-t 4 -K 10M` (four threads and batches of 10 Mbp of reads in place of 100 Mbp) it
/// maps a human genome in under 9 GB, to which the two sort batches add 2. Other thread counts
/// move each by [`ALIGN_THREAD_GIB`] a thread.
///
/// None of these changes what bonefish writes, so they are left out of the stage's fingerprint
/// ([`crate::run`]), and a run resumed with other threads or less memory does not realign: its
/// output is the same at every thread count and batch size (`-K`; docs/usage.md and
/// docs/validation/2026-10-14-scaling-memory), and the sorter's batches (`--sort-mem`) are sorted
/// stably and merged with ties to the earlier batch, so they keep the mapper's order whatever their
/// size (bonefish-core/src/sorted_bam.rs: `Batch::sort`, `merge_and_write`).
const ALIGN_SETTINGS: [AlignSettings; 3] = [
    AlignSettings { args: &[], max_threads: None, peak_gib: 23.0, threads: 16 },
    AlignSettings { args: &["--sort-mem", "1G"], max_threads: None, peak_gib: 19.0, threads: 16 },
    AlignSettings { args: &["--sort-mem", "1G", "-K", "10M"], max_threads: Some(4), peak_gib: 11.0, threads: 4 },
];

/// What each thread adds to bonefish's peak, in GiB: with 10 Mbp batches its SAM output peaked at
/// 8.6 GB with 4 threads and 10.4 GB with 16, and with 500 Mbp batches at 13.5 GB with 1 and 15.0
/// with 16 (bonefish's docs/validation/2026-10-14-scaling-memory).
const ALIGN_THREAD_GIB: f64 = 0.15;

/// bonefish's peak memory in GiB building its index of a human genome (GRCh38: 10.5 GB, in its
/// docs/validation/2026-10-14-aligner-comparison), which it does before mapping when it is given
/// the reference rather than an index: the stage's peak is the larger of the two.
const INDEX_BUILD_GIB: f64 = 10.5;

/// What counting bedMethyl adds to the read pass's peak memory, in GiB (ONT 2026-11-01: 6.2
/// GiB with it against 5.7): the CpG counts of the contig being read and of the one being
/// written.
const BEDMETHYL_MEMORY: f64 = 0.5;

/// What a haplotag list of a 30x genome takes in memory, in GiB, at its peak while it is read:
/// about 24 bytes and the name per read it names ([`lrpipe_reads::haplotags::Haplotags`]).
/// Measured 2026-11-03 on whole-genome HG002: 6.53 million HiFi reads listed in 378 MiB, 400 MiB
/// at the peak of reading it; 8.16 million ONT reads in 467 MiB, 470 at the peak.
const HAPLOTAG_LIST_MEMORY: f64 = 0.75;

/// What writing the haplotagged reads (`--haplotagged-reads`) adds to the haplotag stage, in GiB:
/// with the BAM, HiFi took 5.4 GiB (2026-11-03), against 2.8 with the list alone.
const HAPLOTAGGED_READS_MEMORY: f64 = 3.0;

/// What `--phasing-evidence` adds to the read pass, in GiB: the read names of the two contigs in
/// memory at a time (HG002 chr1, 16 threads: 19 MB more on HiFi, 17 MB on ONT).
const PHASING_EVIDENCE_MEMORY: f64 = 0.1;

/// What tells a reference apart for the `.mdr` of N runs lrpipe-cnv makes from it: its path (as
/// [`provenance::canonical`] spells it), size, modification time and `.fai`. A run in the same
/// output directory with another reference (GRCh38, then CHM13) makes its own `.mdr` instead of
/// masking the earlier reference's gaps (LRUN-08). What cannot be read is left out.
#[must_use]
pub fn reference_key(reference: &Path) -> String {
    let mut h = Sha256::new();
    h.update(provenance::canonical(reference).as_os_str().as_bytes());
    if let Ok(m) = fs::metadata(reference) {
        h.update(format!("\0{}\0{}.{}", m.len(), m.mtime(), m.mtime_nsec()));
    }
    h.update(b"\0");
    if let Ok(fai) = fs::read(Layout::index(reference, "fai")) {
        h.update(fai);
    }
    hex(&h.finalize()[..8])
}

/// Why a run could not be planned.
#[derive(Debug, thiserror::Error)]
pub enum PlanError {
    /// The configuration names no lrpipe-snv model for the platform.
    #[error("the site configuration has no snv.{0}_model")]
    NoModel(&'static str),
    /// Tandem repeats are asked for without a catalogue.
    #[error("tandem repeats (tr) need a [repeats] catalog in the site configuration")]
    NoRepeatCatalog,
    /// Nothing to run.
    #[error("no analyses chosen")]
    NoAnalyses,
    /// Haploid bacterial calling supports the small-variant stage.
    #[error("--bacterial supports only the snv analysis (SNVs and small indels)")]
    BacterialAnalyses,
    /// CNVs are asked for in regions.
    #[error(
        "cnv cannot be called in --regions: lrpipe-cnv normalises each window's depth against the whole \
         genome's, which the regions alone do not give"
    )]
    CnvInRegions,
}

/// A tool binary: the configured path, else its usual name on `PATH`.
fn tool(configured: Option<&Path>, name: &str) -> PathBuf {
    configured.map_or_else(|| PathBuf::from(name), Path::to_path_buf)
}

/// Builds argument lists from strings and paths alike.
#[derive(Debug, Default)]
struct Args(Vec<OsString>);

impl Args {
    fn arg(mut self, a: impl Into<OsString>) -> Self {
        self.0.push(a.into());
        self
    }

    fn opt(self, key: &str, value: impl Into<OsString>) -> Self {
        self.arg(key).arg(value)
    }
}

/// The stages of `sample`'s run that its chosen analyses need, in dependency order.
///
/// # Errors
///
/// If no analyses are chosen, the configuration has no lrpipe-snv model for the sample's
/// platform and the SNVs are needed, tandem repeats are chosen without a catalogue, or a
/// bacterial run requests analyses other than SNVs and small indels.
pub fn plan(site: &Site, sample: &Sample, options: Options, layout: &Layout) -> Result<Vec<Stage>, PlanError> {
    let mut stages = plan_calling(site, sample, options, layout)?;
    // The stage writing `sample.reads`, when the input is not them already.
    let first = if !sample.unaligned.is_empty() {
        Some(align_stage(site, sample, options, layout))
    } else {
        sample.merged.as_ref().map(|merged| merge_reads_stage(site, sample, merged, options, layout))
    };
    if let Some(first) = first {
        // Every stage reading the reads follows the one writing them.
        for stage in stages.iter_mut().filter(|s| s.inputs.contains(&sample.reads)) {
            stage.after.insert(0, first.name);
        }
        stages.insert(0, first);
    }
    // Last, after every stage writing what it reads.
    let qc = crate::qc::stage(site, sample, layout, &stages);
    stages.push(qc);
    Ok(stages)
}

/// The align stage: bonefish maps the sample's unaligned reads, one file after another, into
/// [`Layout::aligned_bam`] and its index, with the preset of its platform (`-x dna` for Oxford
/// Nanopore, `-x hifi` for PacBio HiFi) and `-y`, which keeps each read's tags (the `MM`/`ML`/`MN`
/// of base modifications, `RG`) and the input's `@RG` lines, which name the basecall model.
/// FASTQ has no `@RG` lines, so its reads get one read group named after the sample (`-R`), in
/// place of any `RG` tag in their comments, which would name a read group the BAM lacks.
///
/// It maps to the site's bonefish index when it has one, else to the reference, and takes the
/// first of [`ALIGN_SETTINGS`] that fits the run's memory. Every core is its: it runs alone, as
/// every other stage follows it. The reference is its input either way: an index is of the
/// reference, which the other stages call against, so another reference realigns. The BAM's
/// index is a BAI, or a CSI for a contig past 2^29 bases (htslib's rule, as bonefish follows
/// it), by the reference's `.fai`, which must list the contigs bonefish maps to (`lrpipe run`
/// checks it lists the index's).
#[must_use]
pub fn align_stage(site: &Site, sample: &Sample, options: Options, layout: &Layout) -> Stage {
    let built = site.align.index.is_none();
    let threads = options.threads.align;
    let settings = ALIGN_SETTINGS
        .iter()
        .find(|s| s.peak_gib(threads, built) <= options.memory_gib)
        .unwrap_or(&ALIGN_SETTINGS[ALIGN_SETTINGS.len() - 1]);
    let target = site.align.index.as_ref().unwrap_or(&site.reference);
    let bam = layout.aligned_bam();
    let preset = match sample.platform {
        Platform::Hifi => "hifi",
        Platform::Ont => "dna",
    };
    let mut args = Args::default().opt("-t", settings.threads(threads).to_string()).opt("-x", preset).arg("-y");
    if sample.unaligned_fastq {
        args = args.opt("-R", format!("@RG\\tID:{0}\\tSM:{0}", sample.name));
    }
    for arg in settings.args {
        args = args.arg(*arg);
    }
    args = args.opt("-o", &bam).opt("--sort-tmp", layout.align_scratch().join(layout.sample())).arg(target);
    for r in &sample.unaligned {
        args = args.arg(r);
    }
    let index =
        crate::merge_reads::bam_index_for(usize::try_from(longest_contig(&site.reference)).unwrap_or(usize::MAX))
            .extension();
    let mut inputs = sample.unaligned.clone();
    inputs.push(site.reference.clone());
    if let Some(index) = &site.align.index {
        inputs.push(index.clone());
    }
    Stage {
        name: StageName::Align,
        program: tool(site.tools.bonefish.as_deref(), "bonefish"),
        args: args.0,
        inputs,
        outputs: vec![bam.clone(), Layout::index(&bam, index)],
        after: Vec::new(),
        memory_gib: settings.peak_gib(threads, built),
        priority: Priority::Critical,
        scratch: Some(layout.align_scratch()),
    }
}

/// The repeat VCF and its index, which the reads pass and lrpipe-tr both write: a `.tbi`, or a
/// `.csi` for a contig longer than tabix allows (2^29 bases). Declaring the index makes a resume
/// rerun the stage when it is missing.
fn repeats_outputs(reference: &Path, layout: &Layout) -> [PathBuf; 2] {
    let vcf = layout.repeats_vcf();
    let suffix = if longest_contig(reference) > 1 << 29 { "csi" } else { "tbi" };
    let index = Layout::index(&vcf, suffix);
    [vcf, index]
}

/// The length of the reference's longest contig, by its `.fai`; 0 if it cannot be read.
fn longest_contig(reference: &Path) -> u64 {
    let fai = fs::read_to_string(Layout::index(reference, "fai")).unwrap_or_default();
    fai.lines().filter_map(|l| l.split('\t').nth(1)?.parse::<u64>().ok()).max().unwrap_or(0)
}

/// The stages calling the sample's aligned reads, without the align stage.
fn plan_calling(site: &Site, sample: &Sample, options: Options, layout: &Layout) -> Result<Vec<Stage>, PlanError> {
    let reference = &site.reference;
    let chosen = |a| options.analyses.contains(a);
    if Analysis::ALL.iter().all(|&a| !chosen(a)) {
        return Err(PlanError::NoAnalyses);
    }
    if sample.sex == Sex::NotApplicable && options.analyses != Analyses::of(&[Analysis::Snv]) {
        return Err(PlanError::BacterialAnalyses);
    }
    if chosen(Analysis::Tr) && site.repeats.is_none() {
        return Err(PlanError::NoRepeatCatalog);
    }
    if chosen(Analysis::Cnv) && sample.regions.is_some() {
        return Err(PlanError::CnvInRegions);
    }
    let methylation = if chosen(Analysis::Methyl) { methylation_report(site, sample) } else { MethylationReport::None };
    // CNVs are called from the SNVs and the read pass's depth and SV breakpoints.
    let cnv = chosen(Analysis::Cnv);
    let phase = chosen(Analysis::Phase);
    let snv = options.analyses.calls_snvs();
    let svs = chosen(Analysis::Sv) || cnv;
    let depth = chosen(Analysis::Depth) || cnv;
    let repeats = site.repeats.as_ref().filter(|_| chosen(Analysis::Tr));
    let bedmethyl = methylation == MethylationReport::Bedmethyl;
    let cpg = methylation == MethylationReport::Cpg;
    // With the SNVs phased, the reads are haplotagged from them for what reads haplotypes: the
    // SVs, which are phased after the read pass, the tandem repeats and the CpG tables per
    // haplotype. CNVs alone need the SVs' breakpoints, not their phases.
    let haplotag = phase
        && options.haplotag == Haplotag::On
        && (chosen(Analysis::Sv) || cpg || repeats.is_some() || options.haplotagged_reads == HaplotaggedReads::On);
    let phase_svs = haplotag && chosen(Analysis::Sv);
    // The repeats in the read pass, or in a stage of their own after the haplotag list.
    let (pass_repeats, listed_repeats) = if haplotag { (None, repeats) } else { (repeats, None) };
    // The stage that writes the SVs' VCF.
    let sv_writer = if phase_svs { StageName::SvPhase } else { StageName::Reads };
    let mut stages = Vec::new();
    let reads = sample.reads.as_path();

    if snv {
        stages.push(snv_stage(site, sample, reads, options, layout, phase)?);
    }

    if svs || depth || pass_repeats.is_some() || bedmethyl {
        let mut args = Args::default()
            .arg("sv")
            .opt("--input", reads)
            .opt("--reference", reference)
            .opt("--sample-id", &sample.name)
            .opt("--threads", options.threads.reads.to_string())
            .arg("--allow-overwrite");
        let mut inputs = vec![reads.to_path_buf(), reference.clone()];
        // The reads' statistics for the QC summary, counted from the records the pass reads
        // anyway (in regions, of the regions' reads alone).
        args = args.opt("--read-stats", layout.read_stats());
        let mut outputs = vec![layout.read_stats()];
        // Whole contigs by `--contig`, which restricts the depth, repeats and methylation too, all
        // counted in the pass; other regions by Sniffles2's `--regions`, and the rest by options
        // of their own below, each read apart from the pass.
        let bed = match &sample.regions {
            None => None,
            Some(Targets::Contigs(contigs)) => {
                for contig in contigs {
                    args = args.opt("--contig", contig);
                }
                None
            }
            Some(Targets::Bed) => {
                args = args.opt("--regions", layout.regions_bed());
                inputs.push(layout.regions_bed());
                Some(layout.regions_bed())
            }
        };
        if phase_svs {
            // Unphased, with what phasing them from the haplotag list needs.
            args = args
                .opt("--vcf", layout.unphased_sv_vcf())
                .opt("--snf", layout.unphased_snf())
                .opt("--phasing-evidence", layout.sv_phasing_evidence());
            outputs.extend([layout.unphased_sv_vcf(), layout.unphased_snf(), layout.sv_phasing_evidence()]);
        } else if svs {
            // The candidates too, for multi-sample calling, which cost next to nothing: on HG003
            // ONT chr20 (a 2.5 GB BAM, 16 threads) the pass took 1.1-1.2 s with them and 1.3-1.8 s
            // without, and 20 MB more memory, for a 0.9 MB `.snf` (some 45 MB for a whole genome).
            args = args.opt("--vcf", layout.sv_vcf()).opt("--snf", layout.snf());
            outputs.extend([layout.sv_vcf(), layout.snf()]);
        }
        if svs && options.sv_mode == SvMode::Accuracy {
            args = args.arg("--accuracy");
        }
        if svs && let Some(trf) = &site.sv.tandem_repeats {
            args = args.opt("--tandem-repeats", trf);
            inputs.push(trf.clone());
        }
        if depth {
            // In regions, mosdepth's `--by` BED of the windows in them, counted from the regions'
            // reads alone, in the pass, rather than from every read of their contigs.
            let by: OsString = if bed.is_some() {
                inputs.push(layout.region_windows());
                layout.region_windows().into()
            } else {
                DEPTH_WINDOW.to_string().into()
            };
            args = args
                .opt("--depth", layout.depth_prefix())
                .opt("--depth-by", by)
                .arg("--depth-fast")
                .opt("--depth-mapq", "20");
            if bed.is_some() {
                args = args.arg("--depth-in-regions");
            }
            outputs.extend([
                layout.depth_windows(),
                Layout::index(&layout.depth_windows(), "csi"),
                layout.depth_summary(),
            ]);
        }
        if let Some(repeats) = pass_repeats {
            args = args
                .opt("--repeats", &repeats.catalog)
                .opt("--repeats-vcf", layout.repeats_vcf())
                .opt("--repeats-platform", sample.platform.as_str())
                .opt("--karyotype", sample.sex.karyotype());
            if let Some(bed) = &bed {
                args = args.opt("--repeats-regions", bed);
            }
            inputs.push(repeats.catalog.clone());
            outputs.extend(repeats_outputs(&site.reference, layout));
        }
        if bedmethyl {
            // modkit's CpG recipe, counted from the pass's reads.
            args = args.opt("--bedmethyl", layout.bedmethyl());
            if let Some(bed) = &bed {
                args = args.opt("--bedmethyl-regions", bed);
            }
            outputs.push(layout.bedmethyl());
        }
        stages.push(Stage {
            name: StageName::Reads,
            program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
            args: args.0,
            inputs,
            outputs,
            after: Vec::new(),
            memory_gib: expected_memory(StageName::Reads, sample.platform)
                + if bedmethyl { BEDMETHYL_MEMORY } else { 0.0 }
                + if phase_svs { PHASING_EVIDENCE_MEMORY } else { 0.0 },
            // Even on ONT, where the pass ends after lrpipe-snv, sharing the CPU equally slows
            // lrpipe-snv more than it speeds the pass (1,315 s against 1,237 s in all).
            priority: Priority::Background,
            scratch: None,
        });
    }

    if haplotag {
        stages.push(haplotag_stage(site, sample, reads, options, layout));
    }

    if let Some(repeats) = listed_repeats {
        let mut args = Args::default()
            .arg("tr")
            .opt("--reads", reads)
            .opt("--genome", reference)
            .opt("--repeats", &repeats.catalog)
            .opt("--output", layout.repeats_vcf())
            .opt("--sample", &sample.name)
            .opt("--karyotype", sample.sex.karyotype())
            .opt("--platform", sample.platform.as_str())
            .opt("--threads", options.threads.reads.to_string())
            .opt("--haplotags", layout.haplotags());
        let mut inputs = vec![reads.to_path_buf(), reference.clone(), repeats.catalog.clone(), layout.haplotags()];
        // The catalogue's loci in the regions.
        if sample.regions.is_some() {
            args = args.opt("--regions", layout.regions_bed());
            inputs.push(layout.regions_bed());
        }
        stages.push(Stage {
            name: StageName::Repeats,
            program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
            args: args.0,
            inputs,
            outputs: repeats_outputs(reference, layout).into(),
            after: vec![StageName::Haplotag],
            memory_gib: expected_memory(StageName::Repeats, sample.platform),
            // The merge waits for it.
            priority: Priority::Critical,
            scratch: None,
        });
    }

    if phase_svs {
        let args = Args::default()
            .arg("sv")
            .arg("phase")
            .opt("--evidence", layout.sv_phasing_evidence())
            .opt("--haplotags", layout.haplotags())
            .opt("--input-vcf", layout.unphased_sv_vcf())
            .opt("--vcf", layout.sv_vcf())
            .opt("--input-snf", layout.unphased_snf())
            .opt("--snf", layout.snf());
        stages.push(Stage {
            name: StageName::SvPhase,
            program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
            args: args.0,
            inputs: vec![
                layout.sv_phasing_evidence(),
                layout.haplotags(),
                layout.unphased_sv_vcf(),
                layout.unphased_snf(),
            ],
            outputs: vec![layout.sv_vcf(), layout.snf()],
            after: vec![StageName::Reads, StageName::Haplotag],
            memory_gib: expected_memory(StageName::SvPhase, sample.platform),
            priority: Priority::Critical,
            scratch: None,
        });
    }

    if cnv {
        let mut args = Args::default()
            .arg("cnv")
            .opt("--coverage", layout.depth_windows())
            .opt("--sample-id", &sample.name)
            .opt("--output-dir", layout.cnv_dir())
            .opt("--reference", reference)
            .opt("--snv", layout.snv_vcf())
            .arg("--improved")
            .opt("--sex", sample.sex.as_str())
            .opt("--sv-breakpoints", layout.sv_vcf())
            .opt("--vcf-style", "standard");
        let mut inputs = vec![layout.depth_windows(), layout.snv_vcf(), layout.sv_vcf(), reference.clone()];
        if let Some(blacklist) = &site.cnv.blacklist {
            args = args.opt("--blacklist", blacklist);
            inputs.push(blacklist.clone());
        }
        // The reference's N runs: the configured `.mdr` (made before the run if it is not
        // there yet), else one that lrpipe-cnv makes for this reference when it is not there,
        // never one left by a run with another (LRUN-08).
        if let Some(metadata) = &site.cnv.metadata {
            args = args.opt("--metadata", metadata);
            inputs.push(metadata.clone());
        } else {
            args = args.opt("--metadata", layout.cnv_metadata(&reference_key(reference)));
        }
        stages.push(Stage {
            name: StageName::Cnv,
            program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
            args: args.0,
            inputs,
            outputs: [layout.cnv_vcf(), layout.cnv_bed()].into_iter().flat_map(with_tabix_index).collect(),
            // The SVs' breakpoints are the phased VCF's, which has the calls the phasing rescue
            // passed too.
            after: [StageName::Snv, StageName::Reads]
                .into_iter()
                .chain(phase_svs.then_some(StageName::SvPhase))
                .collect(),
            memory_gib: expected_memory(StageName::Cnv, sample.platform),
            priority: Priority::Critical,
            scratch: None,
        });
    }

    if cpg {
        stages.push(cpg_tables(site, sample, reads, options, layout, haplotag));
    }

    if options.merged_vcf == MergedVcf::On {
        let vcfs: Vec<(&str, PathBuf, StageName)> = [
            (
                chosen(Analysis::Snv) || phase,
                "--snv",
                if phase { layout.phased_snv_vcf() } else { layout.snv_vcf() },
                StageName::Snv,
            ),
            (chosen(Analysis::Sv), "--sv", layout.sv_vcf(), sv_writer),
            (cnv, "--cnv", layout.cnv_vcf(), StageName::Cnv),
            (
                repeats.is_some(),
                "--tr",
                layout.repeats_vcf(),
                if listed_repeats.is_some() { StageName::Repeats } else { StageName::Reads },
            ),
        ]
        .into_iter()
        .filter(|(wanted, ..)| *wanted)
        .map(|(_, flag, path, stage)| (flag, path, stage))
        .collect();
        if !vcfs.is_empty() {
            stages.push(merge_stage(site, sample, layout, &vcfs));
        }
    }
    Ok(stages)
}

/// The `merge-reads` stage: the sample's several BAMs and CRAMs merged into one sorted BAM,
/// [`Layout::merged_reads`], with its index. It streams them, merging their records by position,
/// and runs before every other stage, alone, with all the threads. The merged BAM takes about the
/// inputs' size in disk space (more when they are CRAM, which is smaller).
fn merge_reads_stage(site: &Site, sample: &Sample, merged: &Merged, options: Options, layout: &Layout) -> Stage {
    let (inputs, index) = (&merged.inputs, merged.index);
    let merged = layout.merged_reads();
    let mut args =
        Args::default().arg("merge-reads").opt("--output", &merged).opt("--threads", options.threads.align.to_string());
    // Each input and its index: either replaced reruns the merge, and so every stage after it.
    let mut read = Vec::new();
    for input in inputs {
        read.push(input.clone());
        read.extend(index_of(input));
    }
    if inputs.iter().any(|i| is_cram(i)) {
        args = args.opt("--reference", &site.reference);
        read.push(site.reference.clone());
    }
    for input in inputs {
        args = args.arg(input);
    }
    Stage {
        name: StageName::MergeReads,
        program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
        args: args.0,
        inputs: read,
        outputs: vec![merged.clone(), Layout::index(&merged, index.extension())],
        after: Vec::new(),
        memory_gib: expected_memory(StageName::MergeReads, sample.platform)
            + CRAM_INPUT_MEMORY * inputs.iter().filter(|i| is_cram(i)).count() as f64,
        priority: Priority::Critical,
        scratch: None,
    }
}

/// lrpipe-snv's stage: its SNVs and indels, and with `phase` the same phased.
fn snv_stage(
    site: &Site,
    sample: &Sample,
    reads: &Path,
    options: Options,
    layout: &Layout,
    phase: bool,
) -> Result<Stage, PlanError> {
    let configured = match sample.platform {
        Platform::Hifi => site.snv.hifi_model.as_ref(),
        Platform::Ont => site.snv.ont_model.as_ref(),
    };
    let model = sample.model.as_ref().or(configured).ok_or(PlanError::NoModel(sample.platform.as_str()))?;
    let mut args = Args::default()
        .arg("snv")
        .arg("run")
        .opt("--bam", reads)
        .opt("--ref", &site.reference)
        .opt("--model-path", model)
        .opt("--platform", sample.platform.as_str())
        .opt("--sample-name", &sample.name)
        .opt("--threads", options.threads.snv.to_string())
        .opt("--device", &site.snv.device)
        // Include read order in the stage fingerprint so resume refreshes calls made in file order.
        .opt("--read-order", "canonical")
        .opt("--output", layout.snv_dir());
    // Clair3's: whole contigs as `--ctg_name`, which calls them as a run of the whole genome
    // does, other regions as `--bed_fn`, whose candidates must overlap them.
    let mut bed = None;
    match &sample.regions {
        None => {}
        Some(Targets::Contigs(contigs)) => args = args.opt("--ctg-name", contigs.join(",")),
        Some(Targets::Bed) => {
            args = args.opt("--bed-fn", layout.regions_bed());
            bed = Some(layout.regions_bed());
        }
    }
    if let Some(nvrtc) = &site.snv.nvrtc {
        args = args.opt("--nvrtc", nvrtc);
    }
    // The sample's sex: chrY is not called in XX samples, and chrX and chrY are haploid outside
    // the PARs in XY ones. Without it, as before, both were called diploid (LRP-003).
    let mut par_bed = None;
    match (sample.sex, options.pars) {
        (Sex::Female, _) => args = args.opt("--gender", "female"),
        (Sex::Male, Pars::Known) => {
            args = args.opt("--gender", "male").opt("--par-regions-bed", layout.par_bed());
            par_bed = Some(layout.par_bed());
        }
        (Sex::Male, Pars::Unknown) => {}
        (Sex::NotApplicable, _) => {
            args = args.arg("--include-all-ctgs").arg("--haploid-precise").arg("--no-phasing-for-fa");
        }
    }
    // `--approximate` changes only CPU calling, which it makes more than twice as fast with 8-bit
    // networks (HG002 ONT chr20: 60 s against 139 s, 99.2% of the PASS calls the same).
    args = args.arg("--approximate");
    // Clair3's `--enable_phasing`: lrpipe-phase phases the final calls too, into
    // phased_merge_output.vcf.gz.
    let mut outputs = with_tabix_index(layout.snv_vcf()).to_vec();
    if phase {
        args = args.arg("--enable-phasing=true");
        outputs.extend(with_tabix_index(layout.phased_snv_vcf()));
    }
    if options.gvcf == Gvcf::On {
        args = args.arg("--gvcf=true");
        outputs.extend(with_tabix_index(layout.snv_gvcf()));
    }
    if options.intermediates == Intermediates::Remove {
        args = args.arg("--remove-intermediate-dir=true");
    }
    Ok(Stage {
        name: StageName::Snv,
        program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
        args: args.0,
        // The networks' weights too: a model replaced in the same directory reruns calling.
        inputs: [
            reads.to_path_buf(),
            site.reference.clone(),
            crate::snv_model::network(model, "pileup"),
            crate::snv_model::network(model, "full_alignment"),
        ]
        .into_iter()
        .chain(par_bed)
        .chain(bed)
        .collect(),
        outputs,
        after: Vec::new(),
        memory_gib: expected_memory(StageName::Snv, sample.platform),
        priority: Priority::Critical,
        scratch: None,
    })
}

/// A bgzipped file and its tabix index, both of which a stage must write.
fn with_tabix_index(path: PathBuf) -> [PathBuf; 2] {
    let index = Layout::index(&path, "tbi");
    [path, index]
}

/// `lrpipe merge`'s stage: `vcfs` (its flag, the file, the stage writing it) as one VCF.
fn merge_stage(site: &Site, sample: &Sample, layout: &Layout, vcfs: &[(&str, PathBuf, StageName)]) -> Stage {
    let mut args = Args::default()
        .arg("merge")
        .opt("--reference", &site.reference)
        .opt("--sample", &sample.name)
        .opt("--output", layout.merged_vcf());
    let mut after: Vec<StageName> = Vec::new();
    for (flag, path, stage) in vcfs {
        args = args.opt(flag, path);
        if !after.contains(stage) {
            after.push(*stage);
        }
    }
    let mut inputs: Vec<PathBuf> = vcfs.iter().map(|(_, path, _)| path.clone()).collect();
    inputs.push(PathBuf::from(format!("{}.fai", site.reference.display())));
    Stage {
        name: StageName::Merge,
        program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
        args: args.0,
        inputs,
        outputs: with_tabix_index(layout.merged_vcf()).to_vec(),
        after,
        memory_gib: expected_memory(StageName::Merge, sample.platform),
        priority: Priority::Critical,
        scratch: None,
    }
}

/// The methylation the configuration asks for on the sample's platform.
fn methylation_report(site: &Site, sample: &Sample) -> MethylationReport {
    match (site.methylation.report, sample.platform) {
        (MethylationReport::Auto, Platform::Ont) => MethylationReport::Bedmethyl,
        (MethylationReport::Auto, Platform::Hifi) => MethylationReport::Cpg,
        (report, _) => report,
    }
}

/// The longest contig a BAI indexes, 2^29 bases; a BAM with a longer one gets a CSI.
const BAI_MAX_LEN: u64 = 1 << 29;

/// The index `lrpipe haplotag --write-index` writes for a BAM aligned to `reference`: `bai`, or
/// `csi` when the reference's `.fai` has a contig longer than a BAI allows (or `bai` if it cannot
/// be read: the reads' header holds the same contigs).
fn bam_index_suffix(reference: &Path) -> &'static str {
    let longest = fs::read_to_string(Layout::index(reference, "fai"))
        .ok()
        .and_then(|fai| fai.lines().filter_map(|l| l.split('\t').nth(1)?.parse::<u64>().ok()).max());
    if longest.is_some_and(|l| l > BAI_MAX_LEN) { "csi" } else { "bai" }
}

/// `lrpipe haplotag`'s stage: the reads' haplotypes from the phased SNVs, as a haplotag list
/// (and with [`HaplotaggedReads::On`] the reads haplotagged, indexed), as WhatsHap's `haplotag`
/// with Clair3's options (`--ignore-read-groups`, the reference). Reads on contigs the phased VCF
/// does not name are untagged.
fn haplotag_stage(site: &Site, sample: &Sample, reads: &Path, options: Options, layout: &Layout) -> Stage {
    let mut args = Args::default()
        .arg("haplotag")
        .opt("--reference", &site.reference)
        .arg("--ignore-read-groups")
        .arg("--untag-missing-contigs")
        .opt("--threads", options.threads.snv.to_string())
        .opt("--output-haplotag-list", layout.haplotags());
    let mut inputs = vec![layout.phased_snv_vcf(), reads.to_path_buf(), site.reference.clone()];
    // The reads in the regions: whole contigs as WhatsHap's `--regions`, other regions from the
    // BED, where a read reaching two is listed once.
    match &sample.regions {
        None => {}
        Some(Targets::Contigs(contigs)) => {
            for contig in contigs {
                args = args.opt("--regions", contig);
            }
        }
        Some(Targets::Bed) => {
            args = args.opt("--regions-bed", layout.regions_bed());
            inputs.push(layout.regions_bed());
        }
    }
    let mut outputs = vec![layout.haplotags()];
    match options.haplotagged_reads {
        HaplotaggedReads::Off => args = args.arg("--list-only"),
        HaplotaggedReads::On => {
            args = args.opt("--output", layout.haplotagged_reads()).arg("--write-index");
            let index = Layout::index(&layout.haplotagged_reads(), bam_index_suffix(&site.reference));
            outputs.extend([layout.haplotagged_reads(), index]);
        }
    }
    let args = args.arg(layout.phased_snv_vcf()).arg(reads);
    Stage {
        name: StageName::Haplotag,
        program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
        args: args.0,
        inputs,
        outputs,
        after: vec![StageName::Snv],
        memory_gib: expected_memory(StageName::Haplotag, sample.platform)
            + if options.haplotagged_reads == HaplotaggedReads::On { HAPLOTAGGED_READS_MEMORY } else { 0.0 },
        // The SVs' phases and the CpG tables wait for it.
        priority: Priority::Critical,
        scratch: None,
    }
}

/// lrpipe-methyl's stage for pb-CpG-tools' tables: all reads', and haplotype 1's and 2's when the
/// reads are haplotagged, by the run (`haplotag`: from the haplotag list, after it is made) or
/// before it. It takes the read pass's threads.
fn cpg_tables(site: &Site, sample: &Sample, reads: &Path, options: Options, layout: &Layout, haplotag: bool) -> Stage {
    let reference = &site.reference;
    let estimate = match site.methylation.estimate {
        Estimate::Count => "count",
        Estimate::Model => "model",
    };
    let mut args = Args::default()
        .arg("methyl")
        .arg("cpg")
        .opt("--bam", reads)
        .opt("--ref", reference)
        .opt("--output-prefix", layout.cpg_prefix())
        .opt("--pileup-mode", estimate)
        .opt("--threads", options.threads.reads.to_string());
    let mut inputs = vec![reads.to_path_buf(), reference.clone()];
    // The CpG sites in the regions.
    if sample.regions.is_some() {
        args = args.opt("--include-bed", layout.regions_bed());
        inputs.push(layout.regions_bed());
    }
    let tables: &[&str] = match (haplotag, options.haplotypes) {
        (true, _) => {
            args = args.opt("--haplotags", layout.haplotags());
            inputs.push(layout.haplotags());
            &["combined", "hap1", "hap2"]
        }
        (false, Haplotypes::Tagged) => &["combined", "hap1", "hap2"],
        (false, Haplotypes::Untagged) => {
            args = args.arg("--combined-only");
            &["combined"]
        }
    };
    Stage {
        name: StageName::Methylation,
        program: tool(site.tools.lrpipe.as_deref(), "lrpipe"),
        args: args.0,
        inputs,
        outputs: tables.iter().map(|t| layout.cpg_table(t)).collect(),
        // Without the haplotag list, methylation reads only the alignments, so it runs from the
        // start alongside the rest.
        after: if haplotag { vec![StageName::Haplotag] } else { Vec::new() },
        memory_gib: expected_memory(StageName::Methylation, sample.platform)
            + if haplotag { HAPLOTAG_LIST_MEMORY } else { 0.0 },
        priority: Priority::Background,
        scratch: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(extra: &str) -> Site {
        let text = format!(
            "reference = \"/ref/genome.fa\"\n[snv]\nhifi_model = \"/models/hifi\"\ndevice = \"cuda:0\"\n{extra}"
        );
        Site::parse(&text, Path::new("/")).expect("valid configuration")
    }

    fn sample(platform: Platform) -> Sample {
        Sample {
            name: "HG002".into(),
            reads: "/reads/hg002.bam".into(),
            unaligned: Vec::new(),
            unaligned_fastq: false,
            merged: None,
            platform,
            sex: Sex::Male,
            model: None,
            regions: None,
        }
    }

    const OPTIONS: Options = Options {
        threads: Threads { snv: 8, reads: 8, align: 8 },
        intermediates: Intermediates::Remove,
        // All but tandem repeats.
        analyses: Analyses(0b110_1111),
        merged_vcf: MergedVcf::On,
        gvcf: Gvcf::Off,
        haplotypes: Haplotypes::Tagged,
        haplotag: Haplotag::On,
        haplotagged_reads: HaplotaggedReads::Off,
        sv_mode: SvMode::Accuracy,
        pars: Pars::Known,
        memory_gib: 64.0,
    };

    fn with(analyses: &[Analysis]) -> Options {
        Options { analyses: Analyses::of(analyses), ..OPTIONS }
    }

    fn words(stage: &Stage) -> String {
        stage.args.iter().map(|a| a.to_string_lossy().into_owned()).collect::<Vec<_>>().join(" ")
    }

    #[test]
    fn full_configuration_runs_every_stage() {
        let site = site("[repeats]\ncatalog = \"/cat.bed\"\n[cnv]\nblacklist = \"/bl.bed\"\n");
        let layout = Layout::new("/out", "HG002");
        let stages = plan(&site, &sample(Platform::Hifi), with(&Analysis::ALL), &layout).expect("plannable");
        let names: Vec<StageName> = stages.iter().map(|s| s.name).collect();
        assert_eq!(names, StageName::ALL[2..], "all but align and merge-reads, for one aligned file");
        let named = |name: StageName| stages.iter().find(|s| s.name == name).expect("planned");
        assert!(words(&stages[0]).contains("--enable-phasing=true"), "{}", words(&stages[0]));
        for weights in ["/models/hifi/pileup.pt", "/models/hifi/full_alignment.pt"] {
            assert!(stages[0].inputs.contains(&PathBuf::from(weights)), "a new model reruns calling: {weights}");
        }
        let merge = named(StageName::Merge);
        assert!(merge.outputs.contains(&PathBuf::from("/out/HG002.vcf.gz.tbi")), "resume needs the index");
        let merge = words(merge);
        assert!(
            merge.ends_with(
                "--output /out/HG002.vcf.gz --snv /out/snv/phased_merge_output.vcf.gz --sv /out/sv/HG002.sv.vcf \
                 --cnv /out/cnv/HG002.vcf.gz --tr /out/repeats/HG002.repeats.vcf.gz"
            ),
            "{merge}"
        );

        let reads = words(named(StageName::Reads));
        assert!(reads.contains("--sample-id HG002"), "{reads}");
        // With phasing the repeats wait for the haplotag list, in a stage of their own.
        assert!(!reads.contains("--repeats"), "{reads}");
        let repeats = named(StageName::Repeats);
        assert_eq!(
            words(repeats),
            "tr --reads /reads/hg002.bam --genome /ref/genome.fa --repeats /cat.bed --output \
             /out/repeats/HG002.repeats.vcf.gz --sample HG002 --karyotype xy --platform hifi --threads 8 --haplotags \
             /out/phase/HG002.haplotags.tsv.gz"
        );
        assert_eq!(repeats.after, [StageName::Haplotag]);
        // A missing index reruns the stage, which writes it (lrpipe#42).
        assert_eq!(repeats.outputs, [layout.repeats_vcf(), Layout::index(&layout.repeats_vcf(), "tbi")]);
        assert!(reads.contains("--depth /out/depth/HG002 --depth-by 1000 --depth-fast --depth-mapq 20"), "{reads}");

        let cnv = named(StageName::Cnv);
        assert!(cnv.inputs.contains(&layout.snv_vcf()), "CNVs use lrpipe-snv's SNVs");
        assert!(cnv.inputs.contains(&layout.depth_windows()), "CNVs use the pass's depth");
        assert!(words(cnv).contains("--sex male"));
        assert_eq!(
            cnv.outputs,
            [
                "/out/cnv/HG002.vcf.gz",
                "/out/cnv/HG002.vcf.gz.tbi",
                "/out/cnv/HG002_cnv.bed.gz",
                "/out/cnv/HG002_cnv.bed.gz.tbi"
            ]
            .map(PathBuf::from)
        );
        let mdr = layout.cnv_metadata(&reference_key(Path::new("/ref/genome.fa")));
        assert!(words(cnv).contains(&format!("--metadata {}", mdr.display())), "{}", words(cnv));

        // Every stage's inputs are either outside the run, the PAR BED written before it, or
        // written by a stage it follows.
        assert!(stages[0].inputs.contains(&layout.par_bed()), "an XY sample's SNVs use the PARs");
        assert!(words(&stages[0]).contains("--gender male --par-regions-bed /out/.lrpipe/par.bed"));
        for stage in &stages {
            for input in stage.inputs.iter().filter(|i| i.starts_with("/out") && **i != layout.par_bed()) {
                let writer =
                    stages.iter().find(|s| s.outputs.contains(input)).unwrap_or_else(|| panic!("{input:?} is written"));
                assert!(
                    stage.after.contains(&writer.name),
                    "{} reads {input:?} before {} writes it",
                    stage.name,
                    writer.name
                );
            }
        }
    }

    #[test]
    fn methylation_follows_the_platform() {
        let layout = Layout::new("/out", "S");
        let site = site("ont_model = \"/models/ont\"\n");
        let hifi = plan(&site, &sample(Platform::Hifi), OPTIONS, &layout).expect("plannable");
        let cpg = hifi.iter().find(|s| s.name == StageName::Methylation).expect("a methylation stage");
        assert!(words(cpg).starts_with("methyl cpg --bam /reads/hg002.bam"), "{}", words(cpg));
        assert_eq!(cpg.program, PathBuf::from("lrpipe"));
        assert!(words(cpg).contains("--output-prefix /out/methylation/S.cpg --pileup-mode count"));
        assert_eq!(cpg.outputs[2], PathBuf::from("/out/methylation/S.cpg.hap2.bed.gz"));
        assert_eq!(cpg.after, [StageName::Haplotag], "the haplotag list of the run's phasing");
        assert_eq!(cpg.priority, Priority::Background);
        // Without phasing, from the start, from the reads' tags.
        let unphased = Options { analyses: Analyses(OPTIONS.analyses.0 & !Analysis::Phase.bit()), ..OPTIONS };
        let hifi = plan(&site, &sample(Platform::Hifi), unphased, &layout).expect("plannable");
        let cpg = hifi.iter().find(|s| s.name == StageName::Methylation).expect("a methylation stage");
        assert!(cpg.after.is_empty() && !words(cpg).contains("--haplotags"), "{}", words(cpg));

        // bedMethyl comes from the read pass.
        let ont = plan(&site, &sample(Platform::Ont), OPTIONS, &layout).expect("plannable");
        assert!(ont.iter().all(|s| s.name != StageName::Methylation));
        let reads = ont.iter().find(|s| s.name == StageName::Reads).expect("a read pass");
        assert!(words(reads).contains("--bedmethyl /out/methylation/S.bedmethyl.gz"), "{}", words(reads));
        assert!(reads.outputs.contains(&layout.bedmethyl()));

        // Reads without HP tags, not phased: the combined table alone, not empty ones per
        // haplotype.
        let untagged = Options { haplotypes: Haplotypes::Untagged, ..unphased };
        let hifi = plan(&site, &sample(Platform::Hifi), untagged, &layout).expect("plannable");
        let cpg = hifi.iter().find(|s| s.name == StageName::Methylation).expect("a methylation stage");
        assert!(words(cpg).contains("--combined-only"), "{}", words(cpg));
        assert_eq!(cpg.outputs, [layout.cpg_table("combined")]);

        let none = self::site("ont_model = \"/models/ont\"\n[methylation]\nreport = \"none\"\n");
        let stages = plan(&none, &sample(Platform::Ont), OPTIONS, &layout).expect("plannable");
        assert!(stages.iter().all(|s| s.name != StageName::Methylation));
    }

    #[test]
    fn phasing_haplotags_the_reads_for_the_svs_and_the_cpg_tables() {
        // F1b: the phased SNVs give a haplotag list, which phases the SVs the read pass called
        // beside lrpipe-snv and gives the CpG tables, whatever tags the reads carry.
        let layout = Layout::new("/out", "S");
        let site = site("ont_model = \"/models/ont\"\n");
        for haplotypes in [Haplotypes::Tagged, Haplotypes::Untagged] {
            let options = Options { haplotypes, ..with(&[Analysis::Phase, Analysis::Sv, Analysis::Methyl]) };
            let stages = plan(&site, &sample(Platform::Hifi), options, &layout).expect("plannable");
            let names: Vec<StageName> = stages.iter().map(|s| s.name).collect();
            use StageName::{Haplotag, Merge, Methylation, Qc, Reads, Snv, SvPhase};
            assert_eq!(names, [Snv, Reads, Haplotag, SvPhase, Methylation, Merge, Qc]);
            let named = |name: StageName| stages.iter().find(|s| s.name == name).expect("planned");

            let reads = named(Reads);
            assert!(words(reads).contains(
                "--vcf /out/sv/S.unphased.sv.vcf --snf /out/sv/S.unphased.snf --phasing-evidence /out/sv/S.phasing-evidence.gz"
            ), "{}", words(reads));
            assert!(reads.after.is_empty(), "the read pass runs beside lrpipe-snv");

            let haplotag = named(Haplotag);
            assert_eq!(
                words(haplotag),
                "haplotag --reference /ref/genome.fa --ignore-read-groups --untag-missing-contigs --threads 8 \
                 --output-haplotag-list /out/phase/S.haplotags.tsv.gz --list-only /out/snv/phased_merge_output.vcf.gz \
                 /reads/hg002.bam"
            );
            assert_eq!(haplotag.after, [Snv]);
            assert_eq!(haplotag.outputs, [layout.haplotags()]);
            assert!(haplotag.inputs.contains(&layout.phased_snv_vcf()), "a new phasing reruns it");

            let sv_phase = named(SvPhase);
            assert_eq!(
                words(sv_phase),
                "sv phase --evidence /out/sv/S.phasing-evidence.gz --haplotags /out/phase/S.haplotags.tsv.gz \
                 --input-vcf /out/sv/S.unphased.sv.vcf --vcf /out/sv/S.sv.vcf --input-snf /out/sv/S.unphased.snf \
                 --snf /out/sv/S.snf"
            );
            assert_eq!(sv_phase.outputs, [layout.sv_vcf(), layout.snf()]);

            let cpg = named(Methylation);
            assert!(words(cpg).ends_with("--haplotags /out/phase/S.haplotags.tsv.gz"), "{}", words(cpg));
            assert_eq!(cpg.outputs, ["combined", "hap1", "hap2"].map(|t| layout.cpg_table(t)));
            assert_eq!(cpg.after, [Haplotag]);
            assert!(cpg.memory_gib > expected_memory(Methylation, Platform::Hifi), "it holds the list");
            assert_eq!(named(Merge).after, [Snv, SvPhase]);
            for stage in &stages {
                for input in stage.inputs.iter().filter(|i| i.starts_with("/out") && **i != layout.par_bed()) {
                    let writer = stages.iter().find(|s| s.outputs.contains(input)).expect("written by a stage");
                    assert!(
                        stage.after.contains(&writer.name),
                        "{} reads {input:?} before {}",
                        stage.name,
                        writer.name
                    );
                }
            }
        }

        // The CNVs take the phased SVs' breakpoints.
        let stages_of = |options| plan(&site, &sample(Platform::Ont), options, &layout).expect("plannable");
        let stages = stages_of(with(&[Analysis::Phase, Analysis::Sv, Analysis::Cnv]));
        let cnv = stages.iter().find(|s| s.name == StageName::Cnv).expect("a cnv stage");
        assert_eq!(cnv.after, [StageName::Snv, StageName::Reads, StageName::SvPhase]);
        // CNVs alone need the SVs' breakpoints, not their phases: nothing is haplotagged.
        let stages = stages_of(with(&[Analysis::Phase, Analysis::Cnv]));
        let names: Vec<StageName> = stages.iter().map(|s| s.name).collect();
        assert_eq!(names, [StageName::Snv, StageName::Reads, StageName::Cnv, StageName::Merge, StageName::Qc]);
        assert!(words(&stages[1]).contains("--vcf /out/sv/S.sv.vcf"), "{}", words(&stages[1]));

        // Tandem repeats are genotyped with the list, after it; the read pass has none.
        let tr_site = self::site("ont_model = \"/models/ont\"\n[repeats]\ncatalog = \"/cat.bed\"\n");
        let stages =
            plan(&tr_site, &sample(Platform::Ont), with(&[Analysis::Phase, Analysis::Tr]), &layout).expect("plannable");
        let names: Vec<StageName> = stages.iter().map(|s| s.name).collect();
        assert_eq!(names, [StageName::Snv, StageName::Haplotag, StageName::Repeats, StageName::Merge, StageName::Qc]);
        assert!(words(&stages[2]).contains("--platform ont"), "{}", words(&stages[2]));
        assert_eq!(stages[3].after, [StageName::Snv, StageName::Repeats]);
        // Without phasing, in the read pass, from the reads' tags.
        let stages =
            plan(&tr_site, &sample(Platform::Ont), with(&[Analysis::Snv, Analysis::Tr]), &layout).expect("plannable");
        assert!(stages.iter().all(|s| s.name != StageName::Repeats && s.name != StageName::Haplotag));
        assert!(words(&stages[1]).contains("--repeats /cat.bed"), "{}", words(&stages[1]));

        // Nothing reads haplotypes: no haplotagging (ONT's bedMethyl is not per haplotype).
        let stages =
            plan(&site, &sample(Platform::Ont), with(&[Analysis::Phase, Analysis::Methyl, Analysis::Depth]), &layout)
                .expect("plannable");
        assert!(stages.iter().all(|s| s.name != StageName::Haplotag), "{stages:?}");
        // Unless the haplotagged reads are asked for: the BAM and its index, beside the list.
        let bam = Options { haplotagged_reads: HaplotaggedReads::On, ..with(&[Analysis::Phase]) };
        let stages = plan(&site, &sample(Platform::Ont), bam, &layout).expect("plannable");
        let haplotag = stages.iter().find(|s| s.name == StageName::Haplotag).expect("a haplotag stage");
        assert!(
            words(haplotag).contains("--output /out/phase/S.haplotagged.bam --write-index /out/snv/"),
            "{}",
            words(haplotag)
        );
        assert!(!words(haplotag).contains("--list-only"));
        let bai = Layout::index(&layout.haplotagged_reads(), "bai");
        assert_eq!(haplotag.outputs, [layout.haplotags(), layout.haplotagged_reads(), bai]);
        let list_only = stages_of(with(&[Analysis::Phase, Analysis::Sv]));
        let list_only = list_only.iter().find(|s| s.name == StageName::Haplotag).expect("a haplotag stage");
        assert!(haplotag.memory_gib > list_only.memory_gib, "the BAM takes more memory");
    }

    #[test]
    fn without_phasing_the_reads_own_tags_are_their_haplotypes() {
        // F1a as before: the SVs called and phased in the read pass, from the reads' tags, and the
        // CpG tables per haplotype only for haplotagged reads.
        let layout = Layout::new("/out", "S");
        let site = site("");
        let options = with(&[Analysis::Snv, Analysis::Sv, Analysis::Methyl]);
        let stages = plan(&site, &sample(Platform::Hifi), options, &layout).expect("plannable");
        let names: Vec<StageName> = stages.iter().map(|s| s.name).collect();
        assert_eq!(names, [StageName::Snv, StageName::Reads, StageName::Methylation, StageName::Merge, StageName::Qc]);
        assert!(words(&stages[1]).contains("--vcf /out/sv/S.sv.vcf --snf /out/sv/S.snf"), "{}", words(&stages[1]));
        assert!(!words(&stages[1]).contains("--phasing-evidence"));
        assert!(stages[2].after.is_empty() && !words(&stages[2]).contains("--haplotags"));
        let untagged = Options { haplotypes: Haplotypes::Untagged, ..options };
        let stages = plan(&site, &sample(Platform::Hifi), untagged, &layout).expect("plannable");
        assert!(words(&stages[2]).contains("--combined-only"));
    }

    #[test]
    fn svs_are_called_for_accuracy_unless_sniffles2_is_asked_for() {
        let layout = Layout::new("/out", "S");
        let site = site("");
        let reads_stage = |options: Options| {
            let stages = plan(&site, &sample(Platform::Hifi), options, &layout).expect("plannable");
            words(stages.iter().find(|s| s.name == StageName::Reads).expect("a read pass"))
        };
        assert!(reads_stage(OPTIONS).contains(" --accuracy"));
        assert!(!reads_stage(Options { sv_mode: SvMode::Sniffles2, ..OPTIONS }).contains("--accuracy"));
        // A read pass without SVs (depth alone) has no SV options.
        assert!(!reads_stage(with(&[Analysis::Depth])).contains("--accuracy"));
    }

    #[test]
    fn the_haplotagged_reads_index_follows_the_longest_contig() {
        // `lrpipe haplotag --write-index` writes a CSI for contigs longer than a BAI allows.
        let dir = std::env::temp_dir().join(format!("lrpipe-bai-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");
        let reference = dir.join("ref.fa");
        for (len, suffix) in [(536_870_912_u64, "bai"), (536_870_913, "csi")] {
            fs::write(dir.join("ref.fa.fai"), format!("chr1\t1000\t6\t60\t61\nbig\t{len}\t9\t60\t61\n")).expect("fai");
            assert_eq!(bam_index_suffix(&reference), suffix, "{len}");
        }
        assert_eq!(bam_index_suffix(&dir.join("none.fa")), "bai");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn minimal_configuration_skips_repeats() {
        let stages = plan(&site(""), &sample(Platform::Hifi), OPTIONS, &Layout::new("/out", "S")).expect("plannable");
        assert_eq!(stages.len(), 8);
        assert!(!words(&stages[1]).contains("--repeats"));
        assert!(!words(&stages[6]).contains("--tr"));
        assert!(words(&stages[0]).contains("--device cuda:0"));
    }

    #[test]
    fn only_the_chosen_analyses_and_what_they_need_run() {
        let layout = Layout::new("/out", "S");
        let site = site("ont_model = \"/models/ont\"\n[repeats]\ncatalog = \"/cat.bed\"\n");
        let names = |stages: &[Stage]| stages.iter().map(|s| s.name).collect::<Vec<_>>();

        // SNVs alone: lrpipe-snv, unphased, and the merge of its calls; with phasing, phased.
        let snv = plan(&site, &sample(Platform::Ont), with(&[Analysis::Snv]), &layout).expect("plannable");
        assert_eq!(names(&snv), [StageName::Snv, StageName::Merge, StageName::Qc]);
        assert!(words(&snv[0]).contains("--read-order canonical"), "{}", words(&snv[0]));
        assert!(!words(&snv[0]).contains("--enable-phasing"), "{}", words(&snv[0]));
        assert!(words(&snv[1]).ends_with("--snv /out/snv/merge_output.vcf.gz"), "{}", words(&snv[1]));
        let phased = plan(&site, &sample(Platform::Ont), with(&[Analysis::Phase]), &layout).expect("plannable");
        assert_eq!(names(&phased), [StageName::Snv, StageName::Merge, StageName::Qc]);
        assert!(words(&phased[0]).contains("--enable-phasing=true"), "{}", words(&phased[0]));
        let tbi = |p: PathBuf| Layout::index(&p, "tbi");
        assert_eq!(
            phased[0].outputs,
            [layout.snv_vcf(), tbi(layout.snv_vcf()), layout.phased_snv_vcf(), tbi(layout.phased_snv_vcf())]
        );
        assert!(words(&phased[1]).ends_with("--snv /out/snv/phased_merge_output.vcf.gz"), "{}", words(&phased[1]));

        // SVs alone: the read pass without depth, repeats or methylation.
        let sv = plan(&site, &sample(Platform::Ont), with(&[Analysis::Sv]), &layout).expect("plannable");
        assert_eq!(names(&sv), [StageName::Reads, StageName::Merge, StageName::Qc]);
        let reads = words(&sv[0]);
        assert!(reads.contains("--vcf /out/sv/S.sv.vcf"), "{reads}");
        assert!(
            !reads.contains("--depth") && !reads.contains("--repeats") && !reads.contains("--bedmethyl"),
            "{reads}"
        );
        assert!(reads.contains("--snf /out/sv/S.snf"), "{reads}");
        assert_eq!(sv[0].outputs, [layout.read_stats(), layout.sv_vcf(), layout.snf()]);

        // Repeats and depth without SVs: no SV calls, and the merge holds only the repeats.
        let tr =
            plan(&site, &sample(Platform::Ont), with(&[Analysis::Tr, Analysis::Depth]), &layout).expect("plannable");
        assert_eq!(names(&tr), [StageName::Reads, StageName::Merge, StageName::Qc]);
        assert!(!words(&tr[0]).contains("--vcf") && !words(&tr[0]).contains("--snf"), "{}", words(&tr[0]));
        assert!(words(&tr[1]).ends_with("--tr /out/repeats/S.repeats.vcf.gz"), "{}", words(&tr[1]));
        // A missing index reruns the pass, which writes it (lrpipe#42).
        let repeats_tbi = Layout::index(&layout.repeats_vcf(), "tbi");
        assert!(tr[0].outputs.contains(&repeats_tbi), "{:?}", tr[0].outputs);

        // CNVs need the SNVs, depth and SV breakpoints, but the merge holds only the CNVs.
        let cnv = plan(&site, &sample(Platform::Ont), with(&[Analysis::Cnv]), &layout).expect("plannable");
        assert_eq!(names(&cnv), [StageName::Snv, StageName::Reads, StageName::Cnv, StageName::Merge, StageName::Qc]);
        assert!(words(&cnv[1]).contains("--vcf") && words(&cnv[1]).contains("--depth"));
        assert!(words(&cnv[3]).ends_with("--output /out/S.vcf.gz --cnv /out/cnv/S.vcf.gz"), "{}", words(&cnv[3]));
        assert_eq!(cnv[3].after, [StageName::Cnv]);

        // Methylation and depth alone on ONT: a pass for bedMethyl and depth, nothing to merge.
        let methyl = plan(&site, &sample(Platform::Ont), with(&[Analysis::Methyl, Analysis::Depth]), &layout)
            .expect("plannable");
        assert_eq!(names(&methyl), [StageName::Reads, StageName::Qc]);
        assert!(words(&methyl[0]).contains("--bedmethyl"));

        // No merged VCF.
        let off = Options { merged_vcf: MergedVcf::Off, ..OPTIONS };
        assert!(
            plan(&site, &sample(Platform::Ont), off, &layout)
                .expect("plannable")
                .iter()
                .all(|s| s.name != StageName::Merge)
        );

        // Tandem repeats without a catalogue, and nothing at all, are errors.
        let bare = self::site("ont_model = \"/models/ont\"\n");
        assert!(matches!(
            plan(&bare, &sample(Platform::Ont), with(&[Analysis::Tr]), &layout),
            Err(PlanError::NoRepeatCatalog)
        ));
        assert!(matches!(plan(&site, &sample(Platform::Ont), with(&[]), &layout), Err(PlanError::NoAnalyses)));
    }

    #[test]
    fn a_gvcf_is_written_only_when_asked_for() {
        let layout = Layout::new("/out", "S");
        let site = site("ont_model = \"/models/ont\"\n");
        let snv = |options| plan(&site, &sample(Platform::Ont), options, &layout).expect("plannable").remove(0);
        let off = snv(OPTIONS);
        assert!(!words(&off).contains("--gvcf"), "{}", words(&off));
        assert!(!off.outputs.contains(&layout.snv_gvcf()));
        let on = snv(Options { gvcf: Gvcf::On, ..OPTIONS });
        assert!(words(&on).contains("--gvcf=true"), "{}", words(&on));
        let gvcf = PathBuf::from("/out/snv/merge_output.gvcf.gz");
        assert_eq!(on.outputs[on.outputs.len() - 2..], [gvcf.clone(), Layout::index(&gvcf, "tbi")]);
        // Only lrpipe-snv's command changes: the other stages are as they were.
        let others = |options| plan(&site, &sample(Platform::Ont), options, &layout).expect("plannable").split_off(1);
        assert_eq!(others(OPTIONS), others(Options { gvcf: Gvcf::On, ..OPTIONS }));
    }

    #[test]
    fn the_cnv_metadata_follows_the_reference() {
        // LRUN-08: lrpipe-cnv makes the `.mdr` only when it is not there, so a run with another
        // reference in the same output directory must name another one.
        let dir = std::env::temp_dir().join(format!("lrpipe-mdr-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");
        let layout = Layout::new("/out", "S");
        let metadata = |site: &Site| {
            let stages = plan(site, &sample(Platform::Hifi), with(&[Analysis::Cnv]), &layout).expect("plannable");
            let cnv = stages.into_iter().find(|s| s.name == StageName::Cnv).expect("a cnv stage");
            let at = cnv.args.iter().position(|a| a == "--metadata").expect("--metadata");
            (PathBuf::from(&cnv.args[at + 1]), cnv)
        };
        let reference = |name: &str, fai: &str| {
            let text = format!("reference = \"{name}\"\n[snv]\nhifi_model = \"/models/hifi\"\n");
            fs::write(dir.join(name), ">chr1\nACGTNNNNNNNNACGT\n").expect("reference");
            fs::write(dir.join(format!("{name}.fai")), fai).expect("index");
            Site::parse(&text, &dir).expect("valid configuration")
        };
        let grch38 = reference("grch38.fa", "chr1\t248956422\t6\t60\t61\n");
        let (mdr, cnv) = metadata(&grch38);
        assert!(mdr.starts_with(layout.cnv_dir()), "{}", mdr.display());
        assert!(!cnv.inputs.contains(&mdr), "lrpipe-cnv makes it from the reference, an input");
        assert_eq!(metadata(&grch38).0, mdr, "the same reference keeps its .mdr");
        let chm13 = reference("chm13.fa", "chr1\t248387328\t6\t60\t61\n");
        assert_ne!(metadata(&chm13).0, mdr, "another reference gets another .mdr");
        // The same path with another build written over it.
        fs::write(dir.join("grch38.fa.fai"), "chr1\t248387328\t6\t60\t61\n").expect("index");
        assert_ne!(metadata(&grch38).0, mdr, "a reference replaced in place gets another .mdr");

        // A configured one is the stage's input, whatever the reference.
        let site = site("[cnv]\nmetadata = \"/ref/genome.mdr\"\n");
        let (mdr, cnv) = metadata(&site);
        assert_eq!(mdr, PathBuf::from("/ref/genome.mdr"));
        assert!(cnv.inputs.contains(&mdr) && !cnv.outputs.contains(&mdr));
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn missing_model_for_platform_is_an_error() {
        let err = plan(&site(""), &sample(Platform::Ont), OPTIONS, &Layout::new("/out", "S"));
        assert!(matches!(err, Err(PlanError::NoModel("ont"))));
    }

    #[test]
    fn the_chosen_model_is_the_snv_stages_and_its_weights_are_inputs() {
        let chosen = Sample { model: Some("/models/r1041_e82_400bps_sup_v500".into()), ..sample(Platform::Ont) };
        let stages = plan(&site(""), &chosen, with(&[Analysis::Snv]), &Layout::new("/out", "S")).expect("plannable");
        assert!(words(&stages[0]).contains("--model-path /models/r1041_e82_400bps_sup_v500 "), "{}", words(&stages[0]));
        for weights in ["pileup.pt", "full_alignment.pt"] {
            let path = Path::new("/models/r1041_e82_400bps_sup_v500").join(weights);
            assert!(stages[0].inputs.contains(&path), "another model reruns calling: {}", path.display());
        }
        // Over the site's for the platform.
        let hifi = Sample { model: Some("/models/hifi_revio".into()), ..sample(Platform::Hifi) };
        let stages = plan(&site(""), &hifi, with(&[Analysis::Snv]), &Layout::new("/out", "S")).expect("plannable");
        assert!(words(&stages[0]).contains("--model-path /models/hifi_revio "), "{}", words(&stages[0]));
        assert!(Analyses::of(&[Analysis::Cnv]).calls_snvs() && !Analyses::of(&[Analysis::Sv]).calls_snvs());
    }

    /// F5: several BAMs and CRAMs are merged first, by a stage every other one follows and reads
    /// the output of; each input and its index are the merge's inputs; one file is read as it is.
    #[test]
    fn several_reads_are_merged_first() {
        let dir = std::env::temp_dir().join(format!("lrpipe-merge-reads-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");
        let (bam, cram) = (dir.join("a.bam"), dir.join("b.cram"));
        for file in [&bam, &Layout::index(&bam, "bai"), &cram, &Layout::index(&cram, "crai")] {
            fs::write(file, "x").expect("written");
        }
        let site = site("ont_model = \"/models/ont\"\n[repeats]\ncatalog = \"/cat.bed\"\n");
        let layout = Layout::new("/out", "S");
        let several = Sample {
            reads: layout.merged_reads(),
            merged: Some(Merged { inputs: vec![bam.clone(), cram.clone()], index: BamIndex::Bai }),
            ..sample(Platform::Ont)
        };
        let stages = plan(&site, &several, with(&Analysis::ALL), &layout).expect("plannable");
        let merge = &stages[0];
        assert_eq!(merge.name, StageName::MergeReads);
        assert_eq!(
            words(merge),
            format!(
                "merge-reads --output /out/reads/S.bam --threads 8 --reference /ref/genome.fa {} {}",
                bam.display(),
                cram.display()
            )
        );
        assert_eq!(
            merge.inputs,
            [
                bam.clone(),
                Layout::index(&bam, "bai"),
                cram.clone(),
                Layout::index(&cram, "crai"),
                PathBuf::from("/ref/genome.fa")
            ]
        );
        assert_eq!(merge.outputs, ["/out/reads/S.bam", "/out/reads/S.bam.bai"].map(PathBuf::from));
        assert_eq!(merge.memory_gib, 0.75, "a CRAM input adds to it");
        assert!(merge.after.is_empty());
        let merged = layout.merged_reads();
        // The stages reading the reads follow it (as they follow the align stage), and read what
        // it writes.
        for stage in &stages[1..] {
            assert_eq!(stage.inputs.contains(&merged), stage.after.contains(&StageName::MergeReads), "{}", stage.name);
            // But the QC, which reads their indexes for the sex, as the run estimated it.
            let named = words(stage).contains("a.bam") || words(stage).contains("b.cram");
            assert_eq!(named, stage.name == StageName::Qc, "{}", words(stage));
        }
        let cnv = stages.iter().find(|s| s.name == StageName::Cnv).expect("planned");
        assert_eq!(cnv.after, [StageName::Snv, StageName::Reads, StageName::SvPhase], "cnv reads none of the reads");
        let haplotag = stages.iter().find(|s| s.name == StageName::Haplotag).expect("planned");
        assert!(haplotag.inputs.contains(&merged) && words(haplotag).ends_with(" /out/reads/S.bam"));
        assert_eq!(haplotag.after, [StageName::MergeReads, StageName::Snv]);
        for name in [StageName::Snv, StageName::Reads] {
            let stage = stages.iter().find(|s| s.name == name).expect("planned");
            assert!(stage.inputs.contains(&merged), "{name} reads the merged BAM: {:?}", stage.inputs);
            assert!(words(stage).contains(" /out/reads/S.bam "), "{}", words(stage));
        }

        // BAMs alone need no reference, and a long reference sequence takes a CSI.
        let bams = Sample {
            reads: layout.merged_reads(),
            merged: Some(Merged { inputs: vec![bam.clone(), bam.clone()], index: BamIndex::Csi }),
            ..sample(Platform::Ont)
        };
        let merge = plan(&site, &bams, with(&[Analysis::Sv]), &layout).expect("plannable").remove(0);
        assert!(!words(&merge).contains("--reference"), "{}", words(&merge));
        assert_eq!(merge.outputs[1], PathBuf::from("/out/reads/S.bam.csi"));

        // One file: no merge, as before.
        let one = plan(&site, &sample(Platform::Ont), with(&Analysis::ALL), &layout).expect("plannable");
        assert!(one.iter().all(|s| s.name != StageName::MergeReads && !s.after.contains(&StageName::MergeReads)));
        assert!(words(&one[0]).contains("--bam /reads/hg002.bam"), "{}", words(&one[0]));
        let _ = fs::remove_dir_all(dir);
    }

    /// Every stage is in [`StageName::ALL`], in order, and is named back from its name: a stage
    /// left out could not be read from a manifest (a new variant fails to compile here until it
    /// is placed).
    #[test]
    fn every_stage_is_listed_and_named_back() {
        let place = |s: StageName| match s {
            StageName::Align => 0,
            StageName::MergeReads => 1,
            StageName::Snv => 2,
            StageName::Reads => 3,
            StageName::Haplotag => 4,
            StageName::Repeats => 5,
            StageName::SvPhase => 6,
            StageName::Cnv => 7,
            StageName::Methylation => 8,
            StageName::Merge => 9,
            StageName::Qc => 10,
        };
        for (i, stage) in StageName::ALL.into_iter().enumerate() {
            assert_eq!(place(stage), i, "{stage}");
            assert_eq!(StageName::named(stage.as_str()), Some(stage));
        }
        assert_eq!(StageName::named("sv-phase"), Some(StageName::SvPhase));
        assert_eq!(StageName::named("phase"), None);
    }

    /// The HG002 sample as unaligned reads: two files, which the align stage maps into the run's BAM.
    fn unaligned(platform: Platform, layout: &Layout) -> Sample {
        Sample {
            reads: layout.aligned_bam(),
            unaligned: vec!["/reads/flowcell1.bam".into(), "/reads/flowcell2.bam".into()],
            ..sample(platform)
        }
    }

    #[test]
    fn unaligned_reads_are_aligned_first_and_every_stage_reads_the_alignments() {
        let layout = Layout::new("/out", "S");
        let site = site("ont_model = \"/models/ont\"\n[tools]\nbonefish = \"/opt/bin/bonefish\"\n");
        let stages = plan(&site, &unaligned(Platform::Ont, &layout), OPTIONS, &layout).expect("plannable");
        let names: Vec<StageName> = stages.iter().map(|s| s.name).collect();
        use StageName::{Align, Cnv, Haplotag, Merge, Qc, Reads, Snv, SvPhase};
        assert_eq!(names, [Align, Snv, Reads, Haplotag, SvPhase, Cnv, Merge, Qc]);
        let align = &stages[0];
        assert_eq!(align.program, PathBuf::from("/opt/bin/bonefish"));
        assert_eq!(
            words(align),
            "-t 8 -x dna -y -o /out/align/S.bam --sort-tmp /out/align/tmp/S /ref/genome.fa /reads/flowcell1.bam \
             /reads/flowcell2.bam"
        );
        let inputs = ["/reads/flowcell1.bam", "/reads/flowcell2.bam", "/ref/genome.fa"].map(PathBuf::from);
        assert_eq!(align.inputs, inputs, "new reads or another reference realign");
        // The reference's .fai cannot be read: its contigs are taken to fit a BAI.
        assert_eq!(align.outputs, ["/out/align/S.bam", "/out/align/S.bam.bai"].map(PathBuf::from));
        assert_eq!(align.scratch, Some(PathBuf::from("/out/align/tmp")));
        assert!(align.after.is_empty());
        // The stages reading the alignments follow it, and read what it writes.
        for stage in &stages[1..] {
            let reads = stage.inputs.contains(&layout.aligned_bam());
            assert_eq!(reads, stage.after.contains(&StageName::Align), "{}", stage.name);
            assert!(!stage.inputs.iter().any(|i| i.starts_with("/reads")), "{} reads the unaligned reads", stage.name);
        }
        assert!(words(&stages[1]).contains("--bam /out/align/S.bam "), "{}", words(&stages[1]));
        assert!(words(&stages[2]).contains("--input /out/align/S.bam "), "{}", words(&stages[2]));
        // The haplotagging reads the alignments too, from the phased SNVs; the CNVs none.
        assert!(
            words(&stages[3]).ends_with("/out/snv/phased_merge_output.vcf.gz /out/align/S.bam"),
            "{}",
            words(&stages[3])
        );
        assert_eq!(stages[3].after, [Align, Snv]);
        assert_eq!(stages[5].after, [Snv, Reads, SvPhase], "cnv reads none of the alignments");

        // The same stages, but for the align stage, as for aligned reads at that path.
        let aligned = Sample { reads: layout.aligned_bam(), ..sample(Platform::Ont) };
        let mut others = plan(&site, &aligned, OPTIONS, &layout).expect("plannable");
        for stage in &mut others {
            if stage.inputs.contains(&layout.aligned_bam()) {
                stage.after.insert(0, StageName::Align);
            }
            // The QC reads the index the align stage writes, which is not on disk here.
            if stage.name == StageName::Qc {
                stage.inputs.insert(1, Layout::index(&layout.aligned_bam(), "bai"));
            }
        }
        assert_eq!(stages[1..], others[..]);
        // HiFi reads take bonefish's hifi preset, and its HiFi CpG tables follow it too.
        let hifi = plan(&site, &unaligned(Platform::Hifi, &layout), OPTIONS, &layout).expect("plannable");
        assert!(words(&hifi[0]).starts_with("-t 8 -x hifi -y "), "{}", words(&hifi[0]));
        let cpg = hifi.iter().find(|s| s.name == StageName::Methylation).expect("CpG tables");
        assert_eq!(cpg.after, [StageName::Align, StageName::Haplotag], "the reads, and the haplotag list of them");
        // FASTQ, which has no read groups, gets one named after the sample.
        let fastq = Sample {
            unaligned: vec!["/reads/r.fastq.gz".into()],
            unaligned_fastq: true,
            ..unaligned(Platform::Ont, &layout)
        };
        let stage = align_stage(&site, &fastq, OPTIONS, &layout);
        assert!(words(&stage).starts_with("-t 8 -x dna -y -R @RG\\tID:HG002\\tSM:HG002 -o "), "{}", words(&stage));
    }

    #[test]
    fn the_align_stage_maps_to_the_index_and_fits_the_memory() {
        let layout = Layout::new("/out", "S");
        let indexed = site("[align]\nindex = \"/ref/genome.raw\"\n");
        let sample = unaligned(Platform::Ont, &layout);
        let align = |site: &Site, memory_gib, threads| {
            let options = Options { memory_gib, threads: Threads::each(threads), ..OPTIONS };
            align_stage(site, &sample, options, &layout)
        };
        let full = align(&indexed, 64.0, 16);
        assert!(
            words(&full).ends_with("/ref/genome.raw /reads/flowcell1.bam /reads/flowcell2.bam"),
            "{}",
            words(&full)
        );
        // The index and the reference it is of: another of either realigns.
        assert!(full.inputs.ends_with(&["/ref/genome.fa", "/ref/genome.raw"].map(PathBuf::from)), "{:?}", full.inputs);
        // bonefish's documented peaks: with 16 threads 23 GB with its defaults and 19 with
        // --sort-mem 1G, and with -t 4 -K 10M under 9, to which the two 1 GB sort batches add 2;
        // 0.15 GiB more or less for each thread more or fewer.
        for (memory, threads, peak, settings) in [
            (64.0, 16, 23.0, "-t 16 -x dna -y -o"),
            (23.0, 16, 23.0, "-t 16 -x dna -y -o"),
            (22.9, 16, 19.0, "-t 16 -x dna -y --sort-mem 1G -o"),
            (16.0, 16, 11.0, "-t 4 -x dna -y --sort-mem 1G -K 10M -o"),
            (16.0, 2, 10.7, "-t 2 -x dna -y --sort-mem 1G -K 10M -o"),
            // Eight threads need 1.2 GiB less than sixteen; 64 need 7.2 more.
            (21.8, 8, 21.8, "-t 8 -x dna -y -o"),
            (29.0, 64, 26.2, "-t 64 -x dna -y --sort-mem 1G -o"),
            (26.1, 64, 11.0, "-t 4 -x dna -y --sort-mem 1G -K 10M -o"),
            // Less than its smallest peak: the smallest, which runs alone.
            (4.0, 16, 11.0, "-t 4 -x dna -y --sort-mem 1G -K 10M -o"),
        ] {
            let stage = align(&indexed, memory, threads);
            assert!(words(&stage).starts_with(settings), "{memory} GiB, {threads} threads: {}", words(&stage));
            assert!((stage.memory_gib - peak).abs() < 1e-9, "{memory} GiB, {threads} threads: {}", stage.memory_gib);
        }
        // Without an index bonefish builds one first, which peaks at 10.5 GB: more than mapping
        // with one thread and 10 Mbp batches.
        let from_reference = site("");
        assert!((align(&indexed, 4.0, 1).memory_gib - 10.55).abs() < 1e-9);
        assert!((align(&from_reference, 4.0, 1).memory_gib - 10.55).abs() < 1e-9);
        assert!(
            words(&align(&from_reference, 64.0, 16))
                .ends_with("/ref/genome.fa /reads/flowcell1.bam /reads/flowcell2.bam")
        );
        assert!(ALIGN_SETTINGS.iter().all(|s| s.peak_gib(1, true) >= INDEX_BUILD_GIB));
        let small = AlignSettings { args: &[], max_threads: None, peak_gib: 8.0, threads: 1 };
        assert_eq!((small.peak_gib(1, false), small.peak_gib(1, true)), (8.0, INDEX_BUILD_GIB));

        // A contig past 2^29 bases needs a CSI, which bonefish writes in place of a BAI.
        let dir = std::env::temp_dir().join(format!("lrpipe-csi-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");
        for (fai, index) in
            [("chr1\t248956422\t6\t60\t61\n", "bai"), ("chr1\t248956422\t6\t60\t61\n1H\t536870913\t9\t60\t61\n", "csi")]
        {
            fs::write(dir.join("ref.fa.fai"), fai).expect("index");
            let site = Site::parse("reference = \"ref.fa\"\n[snv]\n", &dir).expect("valid configuration");
            assert_eq!(align(&site, 64.0, 1).outputs[1], PathBuf::from(format!("/out/align/S.bam.{index}")));
        }
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn every_stage_gets_at_least_one_thread() {
        assert_eq!(Threads::each(16), Threads { snv: 16, reads: 16, align: 16 });
        assert_eq!(Threads::each(0), Threads { snv: 1, reads: 1, align: 1 });
        // The align stage runs alone, with every core.
        assert_eq!(Threads::split(7), Threads { snv: 4, reads: 3, align: 7 });
        assert_eq!(Threads::split(1), Threads { snv: 1, reads: 1, align: 1 });
    }

    #[test]
    fn regions_restrict_each_stage_by_its_tools_own_option() {
        let site = site("ont_model = \"/models/ont\"\n[repeats]\ncatalog = \"/cat.bed\"\n");
        let layout = Layout::new("/out", "S");
        let (bed, windows) = ("/out/.lrpipe/regions.bed", "/out/.lrpipe/regions.windows.bed");
        let analyses = [Analysis::Snv, Analysis::Phase, Analysis::Sv, Analysis::Tr, Analysis::Depth, Analysis::Methyl];
        let in_regions = |platform, targets| Sample { regions: Some(targets), ..sample(platform) };
        let options = Options { haplotagged_reads: HaplotaggedReads::On, ..with(&analyses) };
        let stages = plan(&site, &in_regions(Platform::Hifi, Targets::Bed), options, &layout).expect("plannable");
        let named = |name: StageName| stages.iter().find(|s| s.name == name).expect("planned");
        let regions_read = |name: StageName| named(name).inputs.contains(&PathBuf::from(bed));
        let snv = words(named(StageName::Snv));
        assert!(snv.contains(&format!("--bed-fn {bed}")) && !snv.contains("--ctg-name"), "{snv}");
        let reads = words(named(StageName::Reads));
        for wanted in [
            format!("--regions {bed}"),
            format!("--depth-by {windows} --depth-fast --depth-mapq 20 --depth-in-regions"),
        ] {
            assert!(reads.contains(&wanted), "{wanted}: {reads}");
        }
        assert!(named(StageName::Reads).inputs.contains(&PathBuf::from(windows)), "new windows rerun the depth");
        assert!(words(named(StageName::Repeats)).ends_with(&format!("--regions {bed}")));
        assert!(words(named(StageName::Haplotag)).contains(&format!("--regions-bed {bed}")));
        assert!(words(named(StageName::Methylation)).contains(&format!("--include-bed {bed}")));
        for name in [StageName::Snv, StageName::Reads, StageName::Haplotag, StageName::Repeats, StageName::Methylation]
        {
            assert!(regions_read(name), "{name}: other regions rerun it");
        }
        // The SVs' phases, read from the pass and the list, and the merge follow from those.
        assert!(!words(named(StageName::SvPhase)).contains("regions"));

        // ONT: the repeats and bedMethyl of the read pass, without the list.
        let options = Options { analyses: Analyses::of(&[Analysis::Tr, Analysis::Methyl]), ..OPTIONS };
        let ont = plan(&site, &in_regions(Platform::Ont, Targets::Bed), options, &layout).expect("plannable");
        let reads = words(&ont[0]);
        for wanted in ["--regions", "--repeats-regions", "--bedmethyl-regions"] {
            assert!(reads.contains(&format!("{wanted} {bed}")), "{wanted}: {reads}");
        }
        assert!(!reads.contains("--depth"), "{reads}");

        // Whole contigs by contig: as the whole genome's runs of them.
        let contigs = Targets::Contigs(vec!["chr20".into(), "chr21".into()]);
        let stages = plan(&site, &in_regions(Platform::Hifi, contigs), with(&analyses), &layout).expect("plannable");
        let named = |name: StageName| stages.iter().find(|s| s.name == name).expect("planned");
        let snv = words(named(StageName::Snv));
        assert!(snv.contains("--ctg-name chr20,chr21") && !snv.contains("--bed-fn"), "{snv}");
        let reads = words(named(StageName::Reads));
        assert!(reads.contains("--contig chr20 --contig chr21") && !reads.contains("--regions"), "{reads}");
        assert!(reads.contains("--depth-by 1000") && !reads.contains("--depth-in-regions"), "{reads}");
        assert!(words(named(StageName::Haplotag)).contains("--regions chr20 --regions chr21"));
        assert!(words(named(StageName::Repeats)).ends_with(&format!("--regions {bed}")));

        // No CNVs in regions, and no regions for the align stage.
        let all = plan(&site, &in_regions(Platform::Hifi, Targets::Bed), with(&Analysis::ALL), &layout);
        assert!(matches!(all, Err(PlanError::CnvInRegions)), "{all:?}");
        let unaligned = Sample { regions: Some(Targets::Bed), ..unaligned(Platform::Ont, &layout) };
        let site = site_with_bonefish();
        let stages = plan(&site, &unaligned, with(&[Analysis::Snv]), &layout).expect("plannable");
        assert_eq!(stages[0].name, StageName::Align);
        assert!(!words(&stages[0]).contains("regions") && !stages[0].inputs.contains(&PathBuf::from(bed)));
    }

    fn site_with_bonefish() -> Site {
        site("ont_model = \"/models/ont\"\n[tools]\nbonefish = \"/opt/bin/bonefish\"\n")
    }

    #[test]
    fn whole_contigs_are_told_apart_and_windows_follow_the_genomes() {
        let contigs: Vec<(String, u64)> = vec![("chr1".into(), 2500), ("chr2".into(), 900)];
        let regions = |intervals: &[(&str, i64, i64)]| {
            Regions::from_intervals(intervals.iter().map(|&(c, s, e)| (c.to_owned(), s, e)))
        };
        let whole = regions(&[("chr1", 0, 2500), ("chr2", 0, 900)]);
        assert_eq!(Targets::of(&whole, &contigs), Targets::Contigs(vec!["chr1".into(), "chr2".into()]));
        assert_eq!(Targets::of(&regions(&[("chr1", 0, 2499)]), &contigs), Targets::Bed);
        assert_eq!(Targets::of(&regions(&[("chr1", 0, 1000), ("chr1", 1001, 2500)]), &contigs), Targets::Bed);
        // The whole genome's 1 kb windows, cut to each region.
        assert_eq!(windows_in(&whole), "chr1\t0\t1000\nchr1\t1000\t2000\nchr1\t2000\t2500\nchr2\t0\t900\n");
        assert_eq!(
            windows_in(&regions(&[("chr1", 500, 2100)])),
            "chr1\t500\t1000\nchr1\t1000\t2000\nchr1\t2000\t2100\n"
        );
        assert_eq!(windows_in(&regions(&[("chr1", 1000, 1001)])), "chr1\t1000\t1001\n");

        // Contigs a run of the whole genome calls nothing on are left out, and every contig it
        // calls, whole, is the whole genome.
        let contigs: Vec<(String, u64)> =
            [("chr1", 2_500_000), ("chr2", 900), ("chrM", 16_569), ("chrUn_x", 1_200_000)]
                .into_iter()
                .map(|(c, l)| (c.to_owned(), l))
                .collect();
        let every =
            regions(&[("chr1", 0, 2_500_000), ("chr2", 0, 900), ("chrM", 0, 16_569), ("chrUn_x", 0, 1_200_000)]);
        let restriction = restrict(&every, &contigs);
        assert_eq!((restriction.targets, restriction.left_out), (None, vec!["chrM".to_owned()]));
        let some = restrict(&regions(&[("chr2", 0, 900), ("chrM", 5, 10)]), &contigs);
        assert_eq!(some.targets, Some(Targets::Contigs(vec!["chr2".into()])));
        assert_eq!((some.regions.bed(), some.left_out), ("chr2\t0\t900\n".to_owned(), vec!["chrM".to_owned()]));
        let none = restrict(&regions(&[("chrM", 5, 10)]), &contigs);
        assert!(none.regions.is_empty(), "{none:?}");
    }
}

//! VCF merging: `SortVcf` merges per-chunk VCFs into one sorted, bgzipped and
//! tabix-indexed VCF, and `MergeVcf` merges the pileup and full-alignment calls of a contig.
//!
//! Ports of `preprocess/SortVcf.py` (`sort_vcf_from`) and `preprocess/MergeVcf.py` from the
//! original (v2.0.3), producing the same bytes, including their gVCF paths: lz4-compressed
//! intermediate gVCFs (`--vcf_fn_suffix .tmp.gvcf`, `MergeVcf --gvcf`) and the final bgzipped gVCF.
//! The original shells out to `pigz`, `lz4`, `bgzip -f` and `tabix -f -p vcf`; this crate does all
//! of them in-process. [`sort_vcf`] and [`merge_vcf`] are the entry points and
//! [`compress_and_index`] exposes the compression step on its own.
//!
//! The stdin variant of `SortVcf` (no `--input_dir`) is not ported.

mod bed;
mod gvcf;
mod hts;
mod malformed;
mod merge_vcf;
mod postprocess;
mod pytext;
mod sort;

use std::path::{Path, PathBuf};

pub use merge_vcf::{GvcfOutput, MergeMode, MergeOptions, MergeSummary, merge_vcf};
pub use snv_call::{Haploid, RefCalls};
pub use sort::sort_vcf;

/// Errors from merging VCFs.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// [`MergeMode::Regions`] needs a contig, where the original raises a `TypeError`
    /// (`None + '.'`).
    #[error("merging by candidate regions (--platform ilmn) needs a contig name")]
    ContigUnset,
    /// A candidate BED row is malformed, where the original raises or exits.
    #[error("invalid candidate BED row {line} in {path:?}: {reason}")]
    CandidateBed {
        /// The concatenated BED file (`full_aln_regions_<contig>`).
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// What is wrong with it.
        #[source]
        reason: RowError,
    },
    /// `--input_dir` (or `--bed_fn_prefix`) does not exist.
    #[error("input directory {path:?} does not exist")]
    InputDir {
        /// The directory.
        path: PathBuf,
    },
    /// `--contigs_fn` is missing or does not exist.
    #[error("cannot find contig file {path:?}")]
    ContigsFile {
        /// The file, if one was given.
        path: Option<PathBuf>,
    },
    /// The merged header has no `##cmdline` and no `--cmd_fn` was given, where the original
    /// raises a `TypeError` (`os.path.exists(None)`).
    #[error("the merged header has no ##cmdline and no command file was given")]
    CmdFileUnset,
    /// A line of an input VCF cannot be processed, where the original raises.
    #[error("malformed line {line} of {path:?}: {reason}")]
    Line {
        /// The input VCF.
        path: PathBuf,
        /// 1-based line number.
        line: u64,
        /// What is wrong with it.
        #[source]
        reason: RowError,
    },
    /// A merged record cannot be post-processed (`--pileup_only`), where the original raises.
    #[error("cannot post-process record {row:?}: {reason}")]
    Record {
        /// The record.
        row: String,
        /// What is wrong with it.
        #[source]
        reason: RowError,
    },
    /// A PAR BED row is malformed, where the original raises or exits.
    #[error("invalid PAR BED row {line} in {path:?}: {reason}")]
    ParBed {
        /// The BED file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// What is wrong with it.
        #[source]
        reason: RowError,
    },
    /// The header for an empty output cannot be built.
    #[error("cannot build the VCF header")]
    Header(#[from] snv_call::Error),
    /// A path is not valid UTF-8 or contains a NUL byte, so htslib cannot open it.
    #[error("unsupported path {path:?}")]
    Path {
        /// The path.
        path: PathBuf,
    },
    /// htslib could not write the BGZF file.
    #[error("cannot write BGZF file {path:?}")]
    Compress {
        /// The output file.
        path: PathBuf,
    },
    /// A file could not be read or written.
    #[error("I/O error on {path:?}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
}

/// Why a VCF or BED line is rejected: the Python exception the original raises on it.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum RowError {
    /// A column is missing (`IndexError`).
    #[error("column {index} is missing")]
    MissingColumn {
        /// 0-based column index.
        index: usize,
    },
    /// A column is not an integer (`ValueError` from `int()`), or is outside `i64`.
    #[error("{value:?} is not an integer")]
    InvalidInteger {
        /// The column.
        value: String,
    },
    /// A column is not a float (`ValueError` from `float()`).
    #[error("{value:?} is not a number")]
    InvalidFloat {
        /// The column.
        value: String,
    },
    /// A BED interval ends before it starts or is negative (`sys.exit` in `bed_tree_from`).
    #[error("invalid interval {start}-{end}")]
    InvalidInterval {
        /// Start column.
        start: i64,
        /// End column.
        end: i64,
    },
    /// A BED interval is empty once padded (`ValueError` from `IntervalTree.addi`).
    #[error("null interval {start}-{end} after padding")]
    NullInterval {
        /// Padded start.
        start: i64,
        /// Padded end.
        end: i64,
    },
    /// A line is not valid UTF-8 (`UnicodeDecodeError`).
    #[error("invalid UTF-8")]
    InvalidUtf8,
}

/// Whether merged records are post-processed (`--pileup_only`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Postprocess {
    /// Records are copied unchanged (the default).
    #[default]
    Off,
    /// Apply `postprocess_row_with_params`: drop RefCalls unless [`SortOptions::ref_calls`]
    /// shows them, rewrite haploid genotypes and mark `LowQual`. The original does the last two
    /// on [`ChunkLayout::Unordered`] only, writing the rows per contig unchanged.
    PileupOnly,
}

/// Sample sex for haploid calling on the sex chromosomes (`--gender`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Gender {
    /// No sex-aware calling (the default).
    #[default]
    Unknown,
    /// Haploid-precise calling on chrX/chrY (or X/Y), outside the PAR regions.
    Male,
    /// No sex-aware calling.
    Female,
}

/// Whether fused lines are repaired (`--check_malformed_records`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MalformedRecords {
    /// Split fused header and record lines (the default).
    #[default]
    Repair,
    /// Copy lines as they are.
    Keep,
}

/// How chunk files map to contigs (`--use_gpu`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChunkLayout {
    /// Each file holds one contig whose name is in the file name (the default): files are read
    /// per contig and a file stops at its first record of another contig.
    #[default]
    ByContig,
    /// Files mix contigs (the original's GPU calling): every file is read once and records are
    /// grouped by contig.
    Unordered,
}

/// Which `##contig` lines the final gVCF (`--vcf_fn_suffix .gvcf`) keeps
/// (`--output_all_contigs_in_gvcf_header`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GvcfContigs {
    /// Only the contigs of [`SortOptions::contigs`] (the default, as GATK does).
    #[default]
    Processed,
    /// Every contig of the inputs' headers.
    All,
}

/// What [`sort_vcf`] reads and writes, mirroring `SortVcf.py`'s arguments.
///
/// A [`suffix`](Self::suffix) containing `gvcf` reads lz4-compressed inputs (as the original
/// does when `lz4` is installed). With `.gvcf` exactly the output is a final gVCF, bgzipped and
/// indexed like a VCF; any other gVCF suffix (`.tmp.gvcf`) writes `output` itself as an lz4
/// frame, neither bgzipped nor indexed, and never replaces it by an empty VCF.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct SortOptions {
    /// Directory holding the chunk VCFs (`--input_dir`).
    pub input_dir: PathBuf,
    /// Uncompressed output path (`--output_fn`). It is removed once compressed to `.gz`, and
    /// the index is written next to that.
    pub output: PathBuf,
    /// Only file names starting with this are read (`--vcf_fn_prefix`).
    pub prefix: Option<String>,
    /// Only file names ending with this are read (`--vcf_fn_suffix`, default `.vcf`).
    pub suffix: String,
    /// Sample name in the header of an empty output (`--sampleName`, default `SAMPLE`).
    pub sample_name: String,
    /// Reference FASTA (`--ref_fn`), whose `.fai` gives the header of an empty output.
    pub reference: Option<PathBuf>,
    /// File listing the contigs to output, one per line (`--contigs_fn`); required unless no
    /// input file matches.
    pub contigs: Option<PathBuf>,
    /// File holding the command line for `##cmdline` (`--cmd_fn`).
    pub cmd_file: Option<PathBuf>,
    /// Record post-processing (`--pileup_only`).
    pub postprocess: Postprocess,
    /// RefCall records under [`Postprocess::PileupOnly`] (`--print_ref_calls`, default hidden).
    pub ref_calls: RefCalls,
    /// Global haploid mode (`--haploid_precise` takes precedence over `--haploid_sensitive`).
    pub haploid: Haploid,
    /// Sample sex (`--gender`).
    pub gender: Gender,
    /// BED of pseudo-autosomal regions kept diploid for [`Gender::Male`] (`--par_regions_bed`;
    /// `EMPTY` is `None`). It is read even when unused, as the original does.
    pub par_regions: Option<PathBuf>,
    /// Variants with QUAL at or below this are marked `LowQual` (`--qual`, default 2); 0
    /// disables marking, as Python's truthiness does.
    pub qual: i32,
    /// Repair of fused lines (`--check_malformed_records`).
    pub malformed: MalformedRecords,
    /// How chunk files map to contigs (`--use_gpu`).
    pub layout: ChunkLayout,
    /// The `##contig` lines of a final gVCF (`--output_all_contigs_in_gvcf_header`).
    pub gvcf_contigs: GvcfContigs,
}

impl SortOptions {
    /// Options with `SortVcf.py`'s defaults, reading `input_dir` and writing `output`.
    pub fn new(input_dir: impl Into<PathBuf>, output: impl Into<PathBuf>) -> Self {
        Self {
            input_dir: input_dir.into(),
            output: output.into(),
            prefix: None,
            suffix: ".vcf".to_owned(),
            sample_name: "SAMPLE".to_owned(),
            reference: None,
            contigs: None,
            cmd_file: None,
            postprocess: Postprocess::Off,
            ref_calls: RefCalls::Hide,
            haploid: Haploid::Off,
            gender: Gender::Unknown,
            par_regions: None,
            qual: 2,
            malformed: MalformedRecords::Repair,
            layout: ChunkLayout::ByContig,
            gvcf_contigs: GvcfContigs::Processed,
        }
    }
}

/// Why [`sort_vcf`] wrote a header-only VCF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum EmptyReason {
    /// No file name starts with the prefix.
    NoPrefixMatch,
    /// No file name ends with the suffix.
    NoSuffixMatch,
    /// The matching files hold no lines.
    NoLines,
    /// The matching files hold no records for the listed contigs.
    NoRecords,
}

/// Whether `tabix` indexed the output. The original ignores a failure; so does [`sort_vcf`], which
/// leaves any older index in place as `tabix` does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexStatus {
    /// `.tbi` written.
    Built,
    /// htslib could not index the file (for example, records out of order).
    Failed,
    /// The output is an lz4-compressed intermediate gVCF, which the original does not index.
    NotIndexed,
}

/// What [`sort_vcf`] wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Summary {
    /// The bgzipped VCF (`<output>.gz`), or the output itself for an lz4-compressed gVCF.
    pub vcf: PathBuf,
    /// Lines written after the header.
    pub records: u64,
    /// Set when the output is a header-only VCF. Never set for an lz4-compressed gVCF.
    pub empty: Option<EmptyReason>,
    /// Whether the `.tbi` index was built.
    pub index: IndexStatus,
    /// Where an empty `pileup.vcf.gz` was also copied (`merge_output.vcf.gz` next to it, as
    /// `print_calling_step` does), if it was.
    pub copied_to: Option<PathBuf>,
}

/// `bgzip -f path`: compresses `path` to `<path>.gz` as htslib's `bgzip` writes it, and removes
/// `path`.
///
/// # Errors
///
/// Fails if `path` cannot be read or the `.gz` cannot be written.
pub fn bgzip(path: &Path) -> Result<PathBuf, Error> {
    hts::bgzip(path)
}

/// `tabix -f -p vcf GZ`: indexes a bgzipped VCF, returning whether the index was built.
///
/// # Errors
///
/// Fails if a path is not UTF-8.
pub fn tabix_index(gz: &Path) -> Result<bool, Error> {
    hts::tabix_vcf(gz)
}

/// `compress_index_vcf`: bgzips `vcf` to `<vcf>.gz` as `bgzip -f` does (removing `vcf`) and
/// indexes it as `tabix -f -p vcf` does.
///
/// Returns the compressed path and whether the index was built; like the original, an indexing
/// failure is not an error.
///
/// # Errors
///
/// Fails if `vcf` cannot be read, the `.gz` cannot be written, or a path is not UTF-8.
pub fn compress_and_index(vcf: &Path) -> Result<(PathBuf, IndexStatus), Error> {
    let parsable = hts::tabix_parses_all(vcf)?;
    let gz = hts::bgzip(vcf)?;
    let index = if parsable && hts::tabix_vcf(&gz)? {
        IndexStatus::Built
    } else {
        log::warn!("tabix could not index {}", gz.display());
        IndexStatus::Failed
    };
    Ok((gz, index))
}

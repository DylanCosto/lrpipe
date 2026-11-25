//! The original's selection steps between pileup and full-alignment calling.
//!
//! A port of `preprocess/SelectQual.py` (both modes) and the phasing path of
//! `preprocess/SelectHetSnp.py` (`FiterHeteSnpPhasing`), version 2.0.3, producing the same
//! bytes:
//!
//! - [`select_phase_qual`] (`SelectQual --phase`): the QUAL cut-off below which heterozygous
//!   SNPs are left out of WhatsHap phasing, written to `phase_qual`.
//! - [`select_het_snps`] / [`select_het_snps_for_phasing`] (`SelectHetSnp`): one contig's
//!   heterozygous SNPs at or above that cut-off, written to `{contig}.vcf`.
//! - [`select_full_qual`] (`SelectQual` without `--phase`): the variant and reference-call QUAL
//!   cut-offs for full-alignment calling, written to `qual`.
//! - [`select_candidates`] / [`select_candidates_for_full_alignment`] (`SelectCandidates`): one
//!   contig's calls below those cut-offs, written as BED chunks (`{contig}.{i}_{n}`) listed in
//!   `FULL_ALN_FILE_{contig}`.
//!
//! The scripts read the pileup VCF through `pigz -fdc`; [`open_vcf`] does the same in process.

mod candidates;
mod het_snp;
mod lines;
mod py;
mod select_qual;

use std::path::PathBuf;

pub use candidates::{
    BedChunk, CandidateBeds, CandidateFilter, CandidateOptions, CandidateWindow, Candidates, DEFAULT_SPLIT_BED_SIZE,
    FULL_ALN_FILE_PREFIX, PhasingInfoFlag, QualCutoffs, SeqEntropyCalling, read_qual_cutoffs, select_candidates,
    select_candidates_for_full_alignment, select_candidates_for_full_alignment_from,
};
pub use het_snp::{
    HetSnpFilter, HetSnpMode, HetSnpOptions, HetSnps, PhaseFlag, read_phase_qual_cutoff, select_het_snps,
    select_het_snps_for_phasing, select_het_snps_for_phasing_from,
};
pub use lines::open_vcf;
pub use select_qual::{
    CandidateVcf, DEFAULT_VAR_PCT_PHASING, FULL_QUAL_FILE, FullQual, FullQualOptions, PHASE_QUAL_FILE, PhaseQual,
    select_full_qual, select_phase_qual,
};

/// ANSI colour codes `shared.utils.log_warning` wraps warnings in.
const WARNING_START: &str = "\x1b[93m";
const WARNING_END: &str = "\x1b[0m";

/// `log_warning(message)`: the message in yellow, as the scripts print it to stderr.
fn log_warning(message: &str) -> String {
    format!("{WARNING_START}{message}{WARNING_END}")
}

/// Errors from reading the pileup VCF, computing cut-offs or writing the outputs.
///
/// Most input errors are places where the Python script would raise (an `IndexError`,
/// `ValueError` or `OverflowError`) and stop.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A file could not be opened, read or written.
    #[error("cannot access {path:?}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// The VCF stream could not be read (for example, corrupt gzip data).
    #[error("cannot read the VCF stream")]
    Read {
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// A line is not valid UTF-8.
    #[error("line {line}: invalid UTF-8")]
    InvalidUtf8 {
        /// The 1-based line number.
        line: u64,
    },
    /// A line is empty or only whitespace, where `SelectHetSnp` indexes its first character.
    #[error("line {line}: empty line")]
    EmptyLine {
        /// The 1-based line number.
        line: u64,
    },
    /// A record has fewer whitespace-separated columns than the script indexes.
    #[error("line {line}: no column {column} (0-based)")]
    MissingColumn {
        /// The 1-based line number.
        line: u64,
        /// The 0-based column Python indexes.
        column: usize,
    },
    /// A QUAL or POS column is not a number Python's `float()` or `int()` accepts.
    #[error("line {line}: invalid number {value:?}")]
    InvalidNumber {
        /// The 1-based line number.
        line: u64,
        /// The column's text.
        value: String,
    },
    /// A QUAL is NaN, which Python's sort orders unpredictably, so the selection cannot be
    /// reproduced.
    #[error("line {line}: NaN quality cannot be sorted reproducibly")]
    NanQuality {
        /// The 1-based line number.
        line: u64,
    },
    /// The FORMAT column contains `PL` as a substring but not as a field (`list.index` raises).
    #[error("line {line}: FORMAT {format:?} contains \"PL\" but has no PL field")]
    PlNotInFormat {
        /// The 1-based line number.
        line: u64,
        /// The FORMAT column.
        format: String,
    },
    /// The sample column has no value at the PL field's index.
    #[error("line {line}: sample column has no value at PL index {index}")]
    PlMissingInSample {
        /// The 1-based line number.
        line: u64,
        /// The PL field's index in FORMAT.
        index: usize,
    },
    /// `int(pct * count)` is not finite, so `int()` raises.
    #[error("proportion {pct} of {count} values is not a finite number")]
    InvalidProportion {
        /// The proportion.
        pct: f64,
        /// The number of values it applies to.
        count: usize,
    },
    /// The phasing cut-off is infinite and `round()` cannot make it an integer.
    #[error("quality cut-off {cutoff} cannot be rounded to an integer")]
    NonFiniteCutoff {
        /// The cut-off.
        cutoff: f64,
    },
    /// The `phase_qual` file does not hold a number.
    #[error("{path:?} does not hold a quality cut-off: {value:?}")]
    InvalidCutoffFile {
        /// The file.
        path: PathBuf,
        /// Its contents.
        value: String,
    },
    /// A `SelectCandidates` option that is not ported was given.
    #[error("SelectCandidates option {option} is not ported")]
    UnsupportedOption {
        /// The option, as on the command line.
        option: &'static str,
    },
    /// Candidate BED lines need a contig name, and none was given (`'\t'.join` raises on
    /// `None`).
    #[error("candidates cannot be written without a contig name (--ctgName)")]
    MissingContigName,
    /// A `SelectHetSnp` mode that is not ported was requested.
    #[error("SelectHetSnp mode {function} is not ported")]
    UnsupportedMode {
        /// The Python function `main()` would call.
        function: &'static str,
    },
}

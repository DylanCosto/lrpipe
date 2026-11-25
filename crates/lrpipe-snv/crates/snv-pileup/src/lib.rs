//! Pileup features and candidate tensors.
//!
//! A port of the original's C pileup (`calculate_*_pileup`) and `CreateTensorPileup` (Python),
//! version 2.0.3, that produces the same tensors and candidate strings. [`candidate_tensors`] is
//! the entry point; [`pileup_counts`] exposes the per-position counts underneath it, and
//! [`candidate_tensors_with_sites`] also returns the positions gVCF output needs.

mod bed;
mod candidates;
mod counts;
mod gvcf;
mod hts;
mod khash;

use std::path::PathBuf;

pub use bed::{Bed, Interval};
pub use candidates::{
    Candidate, ChunkRegion, FLANKING_BASES, POSITIONS, TENSOR_LEN, alt_info_depth, candidate_tensors,
    candidate_tensors_with_sites, known_sites_chunk, scale_for_depth,
};
pub use counts::{
    CandidateAlleles, CandidateParams, FEATURES, FlankRequirement, PileupCounts, SiteDepth, pileup_counts,
};
pub use gvcf::{ChunkSites, Site};
pub use hts::{AlignmentFile, Column, Entry, Fasta, Pileup, Read, Reads};

/// Errors from reading inputs or building features.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A path is not valid UTF-8 or contains a NUL byte.
    #[error("unsupported path {path:?}")]
    Path {
        /// The path.
        path: PathBuf,
    },
    /// htslib could not open a file.
    #[error("cannot open {path:?}")]
    Open {
        /// The file.
        path: PathBuf,
    },
    /// htslib could not load an alignment index.
    #[error("cannot load the index of {path:?}")]
    Index {
        /// The alignment file.
        path: PathBuf,
    },
    /// A region string could not be parsed or queried.
    #[error("invalid region {region:?}")]
    Region {
        /// The region string.
        region: String,
    },
    /// A contig is missing from the reference.
    #[error("contig {contig:?} is not in the reference")]
    Contig {
        /// The contig name.
        contig: String,
    },
    /// htslib could not start or advance a pileup.
    #[error("pileup failed")]
    Pileup,
    /// An alignment record could not be decoded.
    #[error("failed to read alignment record (BAM/CRAM decode error)")]
    Decode,
    /// A file could not be read.
    #[error("cannot read {path:?}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// A BED row is malformed.
    #[error("invalid BED row {line} in {path:?}")]
    Bed {
        /// The BED file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
    },
    /// Chunking by BED was requested without a BED file.
    #[error("chunking by BED needs an --extend_bed file")]
    MissingBed,
    /// The BED file has no rows on the contig.
    #[error("the BED file has no intervals on {contig:?}")]
    EmptyBed {
        /// The contig name.
        contig: String,
    },
}

//! The original's short-read realignment and Illumina full-alignment tensors.
//!
//! For `--platform ilmn`, each full-alignment job of the original (2.0.3) pipes `RealignReads.py`
//! (de Bruijn graph haplotypes, `debruijn_graph.cpp`, and Smith-Waterman realignment,
//! `realigner.cpp` with SSW) into `CreateTensorFullAlignment.py` through `samtools mpileup -`.
//! [`full_alignment_candidates`] runs that chain in-process and returns the tensors.

mod cppsort;
mod dbg;
mod full_alignment;
mod pyrandom;
mod reads;
mod realigner;
mod ssw;

use std::path::PathBuf;

pub use full_alignment::{
    FullAlignmentCandidate, FullAlignmentJob, Platform, full_alignment_candidates, padded_inputs,
};
pub use reads::{RealignInput, RealignedReads};

/// Errors from realignment and tensor building.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Reading alignments or building the pileup failed.
    #[error(transparent)]
    Mpileup(#[from] snv_mpileup::Error),
    /// The reference could not be read.
    #[error(transparent)]
    Reference(#[from] snv_pileup::Error),
    /// A file could not be read.
    #[error("cannot read {path:?}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// The original exits (or Python raises) here: the message.
    #[error("{0}")]
    Exit(String),
}

//! The original's short-read (Illumina) tensors, built from `samtools mpileup` columns.
//!
//! For `--platform ilmn` the original (2.0.3) does not use its C pileup code:
//! `CreateTensorPileup.py` and `CreateTensorFullAlignment.py` parse the text of `samtools mpileup`
//! instead. This crate reproduces that text column by column with htslib's pileup engine
//! ([`text::read_bases`]) and ports the Python that turns it into tensors: [`pileup_candidates`]
//! for the pileup model.

mod bed;
mod hts;
mod pileup;
pub mod text;

use std::path::{Path, PathBuf};

pub use bed::Regions;
pub use hts::{AlignmentFile, Column, Entry, Mpileup, Positions, ReadFilter, Record, RecordView, View};
pub use pileup::{FastMode, GvcfChunk, GvcfSites, PileupJob, pileup_candidates, pileup_candidates_with_sites};

/// `param.SAMTOOLS_VIEW_FILTER_FLAG`: unmapped, mate unmapped, secondary and supplementary.
pub const EXCLUDED_FLAGS: u16 = 2316;

/// Errors from reading inputs or building tensors.
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
    /// htslib could not start or advance a pileup (for example on unsorted input).
    #[error("samtools mpileup failed")]
    Pileup,
    /// An alignment record could not be decoded.
    #[error("failed to read alignment record (BAM/CRAM decode error)")]
    Decode,
    /// A SAM line could not be parsed.
    #[error("cannot parse SAM line {line:?}")]
    Parse {
        /// The line.
        line: String,
    },
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
    /// The reference or a BED file could not be read.
    #[error(transparent)]
    Reference(#[from] snv_pileup::Error),
    /// The original exits here: the message it prints.
    #[error("{0}")]
    Exit(String),
}

/// Reads a plain or gzipped text file, as `pigz -fdc` does.
fn read_text(path: &Path) -> Result<String, Error> {
    use std::io::Read as _;
    let io_err = |source| Error::Io { path: path.to_owned(), source };
    let bytes = std::fs::read(path).map_err(io_err)?;
    let mut text = String::new();
    if bytes.starts_with(&[0x1f, 0x8b]) {
        flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut text).map_err(io_err)?;
    } else {
        text = String::from_utf8_lossy(&bytes).into_owned();
    }
    Ok(text)
}

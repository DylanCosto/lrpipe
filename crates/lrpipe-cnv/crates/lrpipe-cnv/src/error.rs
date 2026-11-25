//! Errors from reading inputs and calling.

use std::io;
use std::path::PathBuf;

/// Why a run could not proceed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// The reference or its indexes could not be opened or validated.
    #[error(transparent)]
    Reference(#[from] lrpipe_reads::OpenError),
    /// A file could not be read.
    #[error("reading {path}")]
    Io {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    /// A line of an input file is not in the expected format.
    #[error("{path}:{line}: {message}")]
    Format { path: PathBuf, line: usize, message: String },
    /// The inputs cannot be analysed together.
    #[error("{0}")]
    Input(String),
}

impl Error {
    pub(crate) fn io(path: impl Into<PathBuf>) -> impl FnOnce(io::Error) -> Error {
        let path = path.into();
        move |source| Error::Io { path, source }
    }

    pub(crate) fn format(path: impl Into<PathBuf>, line: usize, message: impl Into<String>) -> Error {
        Error::Format { path: path.into(), line, message: message.into() }
    }
}

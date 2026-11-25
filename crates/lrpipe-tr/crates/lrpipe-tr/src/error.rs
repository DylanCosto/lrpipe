//! Errors from reading inputs and genotyping.

use std::io;
use std::path::PathBuf;

/// Why a run could not proceed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A file could not be read or written.
    #[error("{path}")]
    Io {
        /// The file.
        path: PathBuf,
        /// What went wrong.
        #[source]
        source: io::Error,
    },
    /// A line of an input file is not in the expected format.
    #[error("{path}:{line}: {message}")]
    Format {
        /// The file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// What is wrong with it.
        message: String,
    },
    /// The inputs cannot be used together.
    #[error("{0}")]
    Input(String),
}

impl Error {
    pub(crate) fn io(path: impl Into<PathBuf>) -> impl FnOnce(io::Error) -> Error {
        let path = path.into();
        move |source| Error::Io { path, source }
    }
}

//! Output files: plain text, or BGZF when the name ends in `.gz`; `-` and `stdout` are the
//! standard output, as modkit takes them.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use pasta::bgzf::io::MultithreadedWriter;

/// Threads compressing a BGZF output.
const BGZF_THREADS: usize = 4;

/// An output file.
pub enum Output {
    /// Plain text.
    Plain(BufWriter<File>),
    /// BGZF, compressed on worker threads.
    Bgzf(MultithreadedWriter<File>),
    /// Plain text on the standard output.
    Stdout(BufWriter<io::Stdout>),
}

impl std::fmt::Debug for Output {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Output::Plain(_) => "Output::Plain",
            Output::Bgzf(_) => "Output::Bgzf",
            Output::Stdout(_) => "Output::Stdout",
        })
    }
}

impl Output {
    /// Creates `path`, BGZF-compressed when it ends in `.gz`; [`is_stdout`] names the standard
    /// output.
    ///
    /// # Errors
    ///
    /// If the file cannot be created.
    pub fn create(path: &Path) -> io::Result<Self> {
        if is_stdout(path) {
            return Ok(Output::Stdout(BufWriter::with_capacity(1 << 20, io::stdout())));
        }
        let file = File::create(path)?;
        Ok(if path.extension().is_some_and(|e| e == "gz") {
            let threads = NonZero::new(BGZF_THREADS).expect("a positive count");
            Output::Bgzf(MultithreadedWriter::with_worker_count(threads, file))
        } else {
            Output::Plain(BufWriter::with_capacity(1 << 20, file))
        })
    }

    /// Flushes everything and, for BGZF, writes the end-of-file block.
    ///
    /// # Errors
    ///
    /// If writing fails.
    pub fn finish(self) -> io::Result<()> {
        match self {
            Output::Plain(mut w) => w.flush(),
            Output::Bgzf(mut w) => w.finish().map(drop),
            Output::Stdout(mut w) => w.flush(),
        }
    }
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Output::Plain(w) => w.write(buf),
            Output::Bgzf(w) => w.write(buf),
            Output::Stdout(w) => w.write(buf),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Output::Plain(w) => w.write_all(buf),
            Output::Bgzf(w) => w.write_all(buf),
            Output::Stdout(w) => w.write_all(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Output::Plain(w) => w.flush(),
            Output::Bgzf(w) => w.flush(),
            Output::Stdout(w) => w.flush(),
        }
    }
}

/// Whether `path` names the standard output: `-` or `stdout`, as modkit's `pileup` takes them.
#[must_use]
pub fn is_stdout(path: &Path) -> bool {
    path == Path::new("-") || path == Path::new("stdout")
}

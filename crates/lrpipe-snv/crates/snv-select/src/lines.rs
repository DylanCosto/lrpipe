//! Reading text lines the way the two Python scripts see them, and opening a VCF the way
//! `pigz -fdc` does.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::Error;

/// How a text stream is split into lines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Newlines {
    /// `sys.stdin` on POSIX: lines end at `\n` only and keep their terminator.
    Lf,
    /// A `Popen(..., universal_newlines=True)` pipe: `\n`, `\r\n` and a lone `\r` all end a
    /// line, and the terminator is dropped (Python turns it into `\n`; every caller strips it).
    Universal,
}

/// The decoded lines of a stream, numbered from 1, failing on invalid UTF-8 as Python's strict
/// decoder does.
pub(crate) struct Lines<R> {
    reader: R,
    newlines: Newlines,
    buf: Vec<u8>,
    pending: VecDeque<String>,
    line_no: u64,
}

impl<R: BufRead> Lines<R> {
    pub(crate) fn new(reader: R, newlines: Newlines) -> Self {
        Self { reader, newlines, buf: Vec::new(), pending: VecDeque::new(), line_no: 0 }
    }

    /// The next line and its number, or `None` at the end of the stream.
    pub(crate) fn next_line(&mut self) -> Result<Option<(u64, String)>, Error> {
        loop {
            if let Some(line) = self.pending.pop_front() {
                self.line_no += 1;
                return Ok(Some((self.line_no, line)));
            }
            self.buf.clear();
            if self.reader.read_until(b'\n', &mut self.buf).map_err(|source| Error::Read { source })? == 0 {
                return Ok(None);
            }
            let chunk = std::str::from_utf8(&self.buf).map_err(|_| Error::InvalidUtf8 { line: self.line_no + 1 })?;
            match self.newlines {
                Newlines::Lf => self.pending.push_back(chunk.to_owned()),
                Newlines::Universal => {
                    // The chunk ends at `\n` or at the end of the stream, so every `\r` in it
                    // other than one right before that `\n` (or at the very end) is a line break.
                    let body = chunk.strip_suffix('\n').unwrap_or(chunk);
                    let body = body.strip_suffix('\r').unwrap_or(body);
                    self.pending.extend(body.split('\r').map(str::to_owned));
                }
            }
        }
    }
}

/// Opens a VCF as `pigz -fdc` would read it: gzip (including multi-member BGZF) is
/// decompressed and any other file is passed through unchanged.
///
/// # Errors
///
/// Fails if the file cannot be opened or read.
pub fn open_vcf(path: &Path) -> Result<impl BufRead, Error> {
    let io_err = |source| Error::Io { path: path.to_owned(), source };
    let mut file = BufReader::new(File::open(path).map_err(io_err)?);
    let is_gzip = file.fill_buf().map_err(io_err)?.starts_with(&[0x1f, 0x8b]);
    let inner: Box<dyn Read + Send> = if is_gzip { Box::new(MultiGzDecoder::new(file)) } else { Box::new(file) };
    Ok(BufReader::with_capacity(1 << 20, inner))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(input: &str, newlines: Newlines) -> Vec<String> {
        let mut lines = Lines::new(input.as_bytes(), newlines);
        let mut out = Vec::new();
        while let Some((_, line)) = lines.next_line().expect("valid input") {
            out.push(line);
        }
        out
    }

    #[test]
    fn lf_keeps_carriage_returns_and_terminators() {
        // Python: ['a\rb\n', 'c\r\n', 'd']
        assert_eq!(all("a\rb\nc\r\nd", Newlines::Lf), ["a\rb\n", "c\r\n", "d"]);
    }

    #[test]
    fn universal_splits_at_every_line_ending() {
        // Python: ['a\n', 'b\n', 'c\n', 'd\n', '\n', 'e']
        assert_eq!(all("a\rb\nc\r\nd\r\re", Newlines::Universal), ["a", "b", "c", "d", "", "e"]);
        assert_eq!(all("a\r", Newlines::Universal), ["a"]);
        assert_eq!(all("\n", Newlines::Universal), [""]);
        assert!(all("", Newlines::Universal).is_empty());
    }

    #[test]
    fn invalid_utf8_is_an_error() {
        let mut lines = Lines::new(&b"ok\n\xff\n"[..], Newlines::Lf);
        assert!(lines.next_line().is_ok());
        assert!(matches!(lines.next_line(), Err(Error::InvalidUtf8 { line: 2 })));
    }
}

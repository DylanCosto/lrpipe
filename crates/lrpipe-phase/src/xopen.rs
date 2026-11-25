//! Files opened as Python's `xopen` (2.x) opens them: compressed with gzip,
//! bzip2 or xz according to the file name extension, or (when reading a file
//! without such an extension) according to its first bytes.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result, bail};

#[derive(Clone, Copy, PartialEq, Debug)]
enum Format {
    Gz,
    Bz2,
    Xz,
    Zst,
}

/// `_detect_format_from_extension`
fn format_from_extension(path: &str) -> Option<Format> {
    [("bz2", Format::Bz2), ("xz", Format::Xz), ("gz", Format::Gz), ("zst", Format::Zst)]
        .into_iter()
        .find(|(ext, _)| path.ends_with(&format!(".{ext}")))
        .map(|(_, f)| f)
}

/// `_detect_format_from_content`
fn format_from_content(bs: &[u8]) -> Option<Format> {
    if bs.starts_with(b"\x1f\x8b") {
        Some(Format::Gz)
    } else if bs.starts_with(b"BZh") {
        Some(Format::Bz2)
    } else if bs.starts_with(b"\xfd7zXZ\x00") {
        Some(Format::Xz)
    } else if bs.starts_with(b"\x28\xb5\x2f\xfd") {
        Some(Format::Zst)
    } else {
        None
    }
}

/// Opens `path` for reading (binary), decompressing as xopen does.
pub fn open(path: &str) -> Result<Box<dyn BufRead>> {
    let mut file = BufReader::new(File::open(path).with_context(|| format!("cannot open {path}"))?);
    let format = match format_from_extension(path) {
        Some(f) => Some(f),
        None => format_from_content(file.fill_buf()?),
    };
    Ok(match format {
        None => Box::new(file),
        Some(Format::Gz) => Box::new(BufReader::new(flate2::bufread::MultiGzDecoder::new(file))),
        Some(Format::Bz2) => Box::new(BufReader::new(bzip2::bufread::MultiBzDecoder::new(file))),
        Some(Format::Xz) => Box::new(BufReader::new(lzma_rust2::XzReader::new(file, true))),
        Some(Format::Zst) => bail!("cannot read {path}: zstandard is not supported"),
    })
}

/// A file opened for writing by `create`
pub enum Writer {
    Plain(BufWriter<File>),
    Stdout(io::Stdout),
    Gz(flate2::write::GzEncoder<BufWriter<File>>),
    Bz2(bzip2::write::BzEncoder<BufWriter<File>>),
    Xz(lzma_rust2::XzWriter<BufWriter<File>>),
}

/// Opens `path` for writing ("-": standard output), compressed according to
/// its extension. `compresslevel` is used for gzip and xz (default 6).
pub fn create(path: &str, compresslevel: Option<u32>) -> Result<Writer> {
    if path == "-" {
        return Ok(Writer::Stdout(io::stdout()));
    }
    let format = format_from_extension(path);
    if format == Some(Format::Zst) {
        bail!("cannot write {path}: zstandard is not supported");
    }
    let file = BufWriter::new(File::create(path).with_context(|| format!("cannot create {path}"))?);
    let level = compresslevel.unwrap_or(6);
    Ok(match format {
        None | Some(Format::Zst) => Writer::Plain(file),
        Some(Format::Gz) => Writer::Gz(flate2::write::GzEncoder::new(file, flate2::Compression::new(level))),
        // Python's bz2 module (xopen ignores the level for bzip2)
        Some(Format::Bz2) => Writer::Bz2(bzip2::write::BzEncoder::new(file, bzip2::Compression::best())),
        Some(Format::Xz) => Writer::Xz(
            lzma_rust2::XzWriter::new(file, lzma_rust2::XzOptions::with_preset(level)).context("xz compression")?,
        ),
    })
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Writer::Plain(w) => w.write(buf),
            Writer::Stdout(w) => w.write(buf),
            Writer::Gz(w) => w.write(buf),
            Writer::Bz2(w) => w.write(buf),
            Writer::Xz(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Writer::Plain(w) => w.flush(),
            Writer::Stdout(w) => w.flush(),
            Writer::Gz(w) => w.flush(),
            Writer::Bz2(w) => w.flush(),
            Writer::Xz(w) => w.flush(),
        }
    }
}

impl Writer {
    /// Writes the end of the compressed stream and flushes the file.
    pub fn finish(self) -> Result<()> {
        match self {
            Writer::Plain(mut w) => w.flush()?,
            Writer::Stdout(mut w) => w.flush()?,
            Writer::Gz(w) => w.finish()?.flush()?,
            Writer::Bz2(w) => w.finish()?.flush()?,
            Writer::Xz(w) => w.finish()?.flush()?,
        }
        Ok(())
    }
}

/// The lines of a UTF-8 text with Python's universal newlines ("\n",
/// "\r\n" and "\r" end a line, returned as "\n"), like iterating over a
/// file opened in text mode.
pub struct TextLines<R: BufRead> {
    r: R,
    path: String,
    pending: std::collections::VecDeque<String>,
    raw: Vec<u8>,
}

impl<R: BufRead> TextLines<R> {
    pub fn new(r: R, path: &str) -> Self {
        TextLines { r, path: path.to_string(), pending: Default::default(), raw: Vec::new() }
    }
}

impl<R: BufRead> Iterator for TextLines<R> {
    type Item = Result<String>;

    fn next(&mut self) -> Option<Result<String>> {
        while self.pending.is_empty() {
            self.raw.clear();
            match self.r.read_until(b'\n', &mut self.raw) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(anyhow::Error::new(e).context(format!("reading {}", self.path)))),
            }
            let ends_line = self.raw.last() == Some(&b'\n');
            let mut body = &self.raw[..self.raw.len() - usize::from(ends_line)];
            if ends_line && body.last() == Some(&b'\r') {
                body = &body[..body.len() - 1];
            }
            let text = match std::str::from_utf8(body) {
                Ok(t) => t,
                Err(_) => return Some(Err(anyhow::anyhow!("UnicodeDecodeError: {} is not UTF-8", self.path))),
            };
            let parts: Vec<&str> = text.split('\r').collect();
            let last = parts.len() - 1;
            for (k, part) in parts.into_iter().enumerate() {
                if k < last || ends_line {
                    self.pending.push_back(format!("{part}\n"));
                } else if !part.is_empty() {
                    self.pending.push_back(part.to_string());
                }
            }
        }
        self.pending.pop_front().map(Ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn universal_newlines() {
        let lines: Vec<String> = TextLines::new(&b"a\r\nb\rc\n\nd"[..], "x").map(|l| l.unwrap()).collect();
        assert_eq!(lines, vec!["a\n", "b\n", "c\n", "\n", "d"]);
    }
}

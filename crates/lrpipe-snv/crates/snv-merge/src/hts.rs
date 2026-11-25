//! `bgzip -f` and `tabix -f -p vcf` in-process: htslib's block layout and indexing, on pasta.

use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write as _};
use std::path::{Path, PathBuf};

use pasta::bgzf::io::raw::{self, EOF_MARKER};
use pasta::bgzf::io::writer::CompressionLevel;
use rayon::prelude::*;

use crate::Error;

/// `WINDOW_SIZE` in `bgzip.c`: `BGZF_BLOCK_SIZE`.
const WINDOW_SIZE: usize = 0xff00;

/// `BGZF_BLOCK_SIZE`: the uncompressed bytes of a full block.
const BLOCK_SIZE: usize = 0xff00;

/// The `bgzf_flush_try` size bgzip uses to force a block boundary after a line.
const FORCE_FLUSH: usize = 65536;

/// Blocks compressed together, on as many threads: 16 MB of text.
const BLOCK_BATCH: usize = 256;

/// `path` with `suffix` appended to its file name.
pub(crate) fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// Whether bgzip treats `path` as text, as `hts_detect_format` decides for the uncompressed
/// files lrpipe-snv writes: not empty, not gzip, and its first kilobyte printable or whitespace.
fn is_textual(path: &Path) -> Result<bool, Error> {
    let io_err = |source| Error::Io { path: path.to_owned(), source };
    let mut head = Vec::with_capacity(1024);
    File::open(path).map_err(io_err)?.take(1024).read_to_end(&mut head).map_err(io_err)?;
    if head.is_empty() || head.starts_with(&[0x1f, 0x8b]) {
        return Ok(false);
    }
    Ok(head.iter().all(|&b| b.is_ascii_graphic() || b.is_ascii_whitespace()))
}

/// A BGZF file being written as htslib's `bgzf_write` writes it: full blocks of `BLOCK_SIZE`
/// bytes, and shorter ones where `flush_try` ends them. The blocks are cut in order and
/// compressed [`BLOCK_BATCH`] at a time in parallel, each on its own as bgzip compresses it, so
/// the file is the same.
struct BgzfWriter {
    out: BufWriter<File>,
    block: Vec<u8>,
    /// Blocks cut but not yet compressed.
    pending: Vec<Vec<u8>>,
}

impl BgzfWriter {
    fn create(path: &Path) -> io::Result<Self> {
        Ok(Self {
            out: BufWriter::new(File::create(path)?),
            block: Vec::with_capacity(BLOCK_SIZE),
            pending: Vec::with_capacity(BLOCK_BATCH),
        })
    }

    /// Compresses the pending blocks and writes them in order.
    fn write_pending(&mut self) -> io::Result<()> {
        let frames = self
            .pending
            .par_iter()
            .map(|block| {
                let mut frame = Vec::new();
                raw::deflate_frame(block, CompressionLevel::default(), &mut frame).map(|()| frame)
            })
            .collect::<io::Result<Vec<Vec<u8>>>>()?;
        for frame in frames {
            self.out.write_all(&frame)?;
        }
        self.pending.clear();
        Ok(())
    }

    fn write(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            let n = (BLOCK_SIZE - self.block.len()).min(data.len());
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.block.len() == BLOCK_SIZE {
                self.flush()?;
            }
        }
        Ok(())
    }

    /// `bgzf_flush`: ends the current block, if it holds anything.
    fn flush(&mut self) -> io::Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        let block = std::mem::replace(&mut self.block, Vec::with_capacity(BLOCK_SIZE));
        self.pending.push(block);
        if self.pending.len() == BLOCK_BATCH {
            self.write_pending()?;
        }
        Ok(())
    }

    /// `bgzf_flush_try(fp, size)`: ends the current block if `size` more bytes would not fit.
    fn flush_try(&mut self, size: usize) -> io::Result<()> {
        if self.block.len() + size > BLOCK_SIZE { self.flush() } else { Ok(()) }
    }

    /// `bgzf_close`: writes the last block and the EOF marker.
    fn close(mut self) -> io::Result<()> {
        self.flush()?;
        self.write_pending()?;
        self.out.write_all(&EOF_MARKER)?;
        self.out.flush()
    }
}

/// Reads until `buf` is full or the file ends, as `hread` does.
fn fill(src: &mut File, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match src.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

/// `bgzip -f path`: compresses `path` to `path.gz` and removes `path`.
///
/// Text input is written as htslib 1.19's `bgzip` writes it, with a block boundary after the
/// header and after the last complete line of each 0xff00-byte window.
pub(crate) fn bgzip(path: &Path) -> Result<PathBuf, Error> {
    let gz = with_suffix(path, ".gz");
    let io_err = |source| Error::Io { path: path.to_owned(), source };
    let textual = is_textual(path)?;
    let mut src = File::open(path).map_err(io_err)?;
    let compress_err = |_: io::Error| Error::Compress { path: gz.clone() };
    let mut out = BgzfWriter::create(&gz).map_err(compress_err)?;
    let mut buffer = vec![0u8; WINDOW_SIZE];
    if textual {
        let (mut in_header, mut n, mut long_line) = (true, 0usize, false);
        loop {
            let c = fill(&mut src, &mut buffer[n..]).map_err(io_err)?;
            if c == 0 {
                break;
            }
            let c2 = c + n;
            let mut flush = false;
            if in_header && (long_line || buffer[0] == b'@' || buffer[0] == b'#') {
                // Scan forward to the start of the last header line.
                let mut last_start = 0;
                n = 0;
                while n < c2 {
                    let byte = buffer[n];
                    n += 1;
                    if byte != b'\n' {
                        continue;
                    }
                    last_start = n;
                    if n < c2 && !(buffer[n] == b'@' || buffer[n] == b'#') {
                        in_header = false;
                        break;
                    }
                }
                if last_start == 0 {
                    n = c2;
                    long_line = true;
                } else {
                    n = last_start;
                    flush = true;
                    long_line = false;
                }
            } else {
                // Scan backwards to the last newline.
                match buffer[..c2].iter().rposition(|&b| b == b'\n') {
                    Some(i) => {
                        flush = true;
                        n = i + 1;
                    }
                    None => n = c2,
                }
            }
            out.write(&buffer[..n]).map_err(compress_err)?;
            if flush {
                out.flush_try(FORCE_FLUSH).map_err(compress_err)?;
            }
            buffer.copy_within(n..c2, 0);
            n = c2 - n;
        }
        out.write(&buffer[..n]).map_err(compress_err)?;
    } else {
        loop {
            let c = fill(&mut src, &mut buffer).map_err(io_err)?;
            if c == 0 {
                break;
            }
            out.write(&buffer[..c]).map_err(compress_err)?;
        }
    }
    out.close().map_err(compress_err)?;
    drop(src);
    fs::remove_file(path).map_err(io_err)?;
    Ok(gz)
}

/// Whether `tbx_parse1` (VCF preset) can parse `line`, a `bgzf_getline` line without its `\n`:
/// meta lines are skipped, and a data line needs a second column (after a tab or NUL) where
/// `strtoll(s, &end, 0)` reads at least one digit. `strtoll` is given the rest of the line,
/// so it skips whitespace (tabs included) past the column, and stops at a NUL.
fn tabix_parses(line: &[u8]) -> bool {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    if line.first() == Some(&b'#') {
        return true;
    }
    let Some(first_delimiter) = line.iter().position(|&b| b == b'\t' || b == 0) else { return false };
    let rest = &line[first_delimiter + 1..];
    let rest = rest.split(|&b| b == 0).next().unwrap_or_default();
    // C `isspace`, which unlike `u8::is_ascii_whitespace` includes vertical tab.
    let skip = rest.iter().take_while(|&&b| matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')).count();
    let rest = &rest[skip..];
    let rest = rest.strip_prefix(b"+").or_else(|| rest.strip_prefix(b"-")).unwrap_or(rest);
    rest.first().is_some_and(u8::is_ascii_digit)
}

/// Whether every line of the uncompressed VCF at `path` is one htslib's indexer accepts.
///
/// The original's `tabix` (htslib 1.24) fails on a line it cannot parse, as pasta's port of it
/// does; checking before compressing says so without building a partial index.
pub(crate) fn tabix_parses_all(path: &Path) -> Result<bool, Error> {
    let io_err = |source| Error::Io { path: path.to_owned(), source };
    let mut reader = io::BufReader::new(File::open(path).map_err(io_err)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if io::BufRead::read_until(&mut reader, b'\n', &mut line).map_err(io_err)? == 0 {
            return Ok(true);
        }
        if !tabix_parses(line.strip_suffix(b"\n").unwrap_or(&line)) {
            return Ok(false);
        }
    }
}

/// `tabix -f -p vcf path`: writes `path.tbi`. Returns whether the index was built.
pub(crate) fn tabix_vcf(path: &Path) -> Result<bool, Error> {
    Ok(snv_hts::index_vcf(path).is_ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tabix_parse_failures() {
        for ok in [
            &b"#CHROM\tPOS"[..],
            b"chr20\t200",
            b"chr20\t200x\t.",
            b"chr20\t-5",
            b"c\t 0x\r",
            b"chr\0 7",
            b"c\t\x0b9",
            b" #x\t1x",
            b"chr20\t\t150",
        ] {
            assert!(tabix_parses(ok), "{:?}", String::from_utf8_lossy(ok));
        }
        for bad in [&b"chr2"[..], b"", b"\r", b"chr20\tabc", b"chr20\t\tx1", b"chr20\t+", b"c\t\x001", b"c\t"] {
            assert!(!tabix_parses(bad), "{:?}", String::from_utf8_lossy(bad));
        }
    }
}

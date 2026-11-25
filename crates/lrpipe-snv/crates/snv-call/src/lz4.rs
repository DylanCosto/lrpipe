//! The lz4 frames the original's intermediate gVCF files are written and read through.
//!
//! The original pipes them through `lz4 -c` and `lz4 -fdc` (`compressReaderWriter`,
//! `preprocess/utils.py`) when `lz4` is on the `PATH`, as it is in the original's conda
//! environment. [`Lz4Writer`] writes frames laid out as lz4 1.10's command-line tool lays them
//! out, with blocks compressed by `lz4_flex` (in Rust, so the bytes of compressed blocks differ
//! from the tool's), and [`open_lz4`] reads any lz4 frame, the tool's included.

use std::fs::File;
use std::hash::Hasher as _;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use lz4_flex::block::{compress_into, get_maximum_output_size};
use lz4_flex::frame::FrameDecoder;
use twox_hash::XxHash32;

/// `LZ4IO_compressFilename_extRess_MT` reads its input in 4 MiB chunks; each becomes one block.
const CHUNK: usize = 4 << 20;
/// Frame magic number, little-endian.
const MAGIC: [u8; 4] = [0x04, 0x22, 0x4d, 0x18];
/// FLG byte: version 01, independent blocks, no block checksums, no content size, content
/// checksum (the `lz4` defaults).
const FLG: u8 = 0x64;
/// Set in a block's size field when the block is stored uncompressed.
const UNCOMPRESSED: u32 = 0x8000_0000;

/// Writes an lz4 frame as `lz4 -c` (lz4 1.10.0) lays one out for the same input.
///
/// The tool reads up to 4 MiB first. A shorter input becomes a single-block frame
/// (`LZ4F_compressFrame`), whose header names the smallest block size that holds it; a longer
/// one becomes a 4 MiB-block frame with each 4 MiB chunk compressed as one independent block.
/// Blocks are compressed by `lz4_flex`, or stored when that does not save a byte, and the frame
/// ends with the content's xxHash32.
///
/// Nothing is written to the inner writer until 4 MiB have been buffered or [`finish`] is
/// called, so dropping the writer without finishing it loses the frame.
///
/// [`finish`]: Lz4Writer::finish
pub struct Lz4Writer<W: Write> {
    out: W,
    chunk: Vec<u8>,
    hasher: XxHash32,
    /// Whether the 4 MiB-block frame header has been written.
    multi_block: bool,
}

impl<W: Write> std::fmt::Debug for Lz4Writer<W> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Lz4Writer")
            .field("buffered", &self.chunk.len())
            .field("multi_block", &self.multi_block)
            .finish_non_exhaustive()
    }
}

impl<W: Write> Lz4Writer<W> {
    /// A writer compressing into `out`.
    pub fn new(out: W) -> Self {
        Self { out, chunk: Vec::with_capacity(CHUNK), hasher: XxHash32::with_seed(0), multi_block: false }
    }

    /// Writes the rest of the frame and returns the inner writer, flushed.
    ///
    /// # Errors
    ///
    /// Fails if the inner writer does.
    pub fn finish(mut self) -> io::Result<W> {
        if !self.multi_block {
            write_header(&mut self.out, block_size_id(self.chunk.len()))?;
        }
        if !self.chunk.is_empty() {
            write_block(&mut self.out, &self.chunk)?;
        }
        self.out.write_all(&[0; 4])?;
        self.out.write_all(&self.hasher.finish_32().to_le_bytes())?;
        self.out.flush()?;
        Ok(self.out)
    }

    fn flush_chunk(&mut self) -> io::Result<()> {
        if !self.multi_block {
            write_header(&mut self.out, 7)?;
            self.multi_block = true;
        }
        write_block(&mut self.out, &self.chunk)?;
        self.chunk.clear();
        Ok(())
    }
}

impl<W: Write> Write for Lz4Writer<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = buf.len().min(CHUNK - self.chunk.len());
        self.chunk.extend_from_slice(&buf[..n]);
        self.hasher.write(&buf[..n]);
        if self.chunk.len() == CHUNK {
            self.flush_chunk()?;
        }
        Ok(n)
    }

    /// Buffered data stays buffered: a partial block is only written by [`Lz4Writer::finish`].
    fn flush(&mut self) -> io::Result<()> {
        self.out.flush()
    }
}

/// `LZ4F_optimalBSID(max4MB, size)`: 64 KB (4), 256 KB (5), 1 MB (6) or 4 MB (7).
fn block_size_id(size: usize) -> u8 {
    let mut id = 4;
    let mut max = 64 << 10;
    while id < 7 && size > max {
        id += 1;
        max <<= 2;
    }
    id
}

fn write_header(out: &mut impl Write, block_size_id: u8) -> io::Result<()> {
    let descriptor = [FLG, block_size_id << 4];
    let checksum = (XxHash32::oneshot(0, &descriptor) >> 8) as u8;
    out.write_all(&MAGIC)?;
    out.write_all(&descriptor)?;
    out.write_all(&[checksum])
}

/// As `LZ4F_makeBlock`: compressed when that saves a byte, else stored.
fn write_block(out: &mut impl Write, data: &[u8]) -> io::Result<()> {
    let mut compressed = vec![0u8; get_maximum_output_size(data.len())];
    let n = compress_into(data, &mut compressed).map_err(io::Error::other)?;
    if n < data.len() {
        out.write_all(&(n as u32).to_le_bytes())?;
        out.write_all(&compressed[..n])
    } else {
        out.write_all(&(data.len() as u32 | UNCOMPRESSED).to_le_bytes())?;
        out.write_all(data)
    }
}

/// Opens `path` as `lz4 -fdc` streams it: an lz4 frame is decompressed and anything else is
/// passed through unchanged.
///
/// # Errors
///
/// Fails if the file cannot be opened or read; the original ignores `lz4` failing, so callers that
/// follow it read such a file as empty.
pub fn open_lz4(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut reader = BufReader::with_capacity(1 << 20, File::open(path)?);
    let is_lz4 = reader.fill_buf()?.starts_with(&MAGIC);
    Ok(if is_lz4 { Box::new(BufReader::with_capacity(1 << 20, FrameDecoder::new(reader))) } else { Box::new(reader) })
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;

    fn compress(data: &[u8], piece: usize) -> Vec<u8> {
        let mut w = Lz4Writer::new(Vec::new());
        for part in data.chunks(piece.max(1)) {
            w.write_all(part).unwrap();
        }
        w.finish().unwrap()
    }

    #[test]
    fn frames_match_the_lz4_tool() {
        // `printf '' | lz4 -c` and `printf 'hello\n' | lz4 -c` (lz4 1.10.0): too short to
        // compress, so the blocks are stored and the frames are the tool's byte for byte.
        assert_eq!(compress(b"", 1), hex("04224d18 6440a7 00000000 055dcc02"));
        assert_eq!(compress(b"hello\n", 1), hex("04224d18 6440a7 06000080 68656c6c6f0a 00000000 f95b6b94"));
    }

    #[test]
    fn reads_the_lz4_tools_frames() {
        // `printf 'abcabcabcabcabcabcabcabcabcabc\n' | lz4 -c` (lz4 1.10.0): a compressed block.
        let frame = hex("04224d18 6440a7 0d000000 3f61626303000450636162630a 00000000 e7e1cdca");
        let mut decoded = Vec::new();
        FrameDecoder::new(frame.as_slice()).read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, b"abcabcabcabcabcabcabcabcabcabc\n");
    }

    fn hex(s: &str) -> Vec<u8> {
        let s: String = s.chars().filter(|c| !c.is_whitespace()).collect();
        (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap()).collect()
    }

    #[test]
    fn block_size_follows_the_input() {
        assert_eq!([0, 65536, 65537, 1 << 20, (1 << 20) + 1].map(block_size_id), [4, 4, 5, 6, 7]);
    }

    #[test]
    fn round_trips_across_chunks() {
        let data: Vec<u8> = (0..(CHUNK * 2 + 12345)).map(|i| (i % 251) as u8 ^ (i / 7919) as u8).collect();
        let frame = compress(&data, 100_003);
        // The header names 4 MB blocks once the input reaches a full chunk.
        assert_eq!(frame[5], 0x70);
        let mut decoded = Vec::new();
        FrameDecoder::new(frame.as_slice()).read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, data);
    }
}

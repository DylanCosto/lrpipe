//! BAM files: records read on top of pasta::bgzf (records keep their
//! bytes), region queries with a BAI/CSI index, and writing.

use std::fs::File;
use std::io::{self, Read, Write};
use std::num::NonZero;
use std::path::Path;

use anyhow::{Context, Result, bail};
use pasta::bgzf::{self, VirtualPosition};
use pasta::core::Position;
use pasta::csi::BinningIndex;
use pasta::csi::binning_index::index::reference_sequence::bin::Chunk;

use super::header::Header;
use super::record::Record;

/// A BGZF reader, decompressing with worker threads if requested.
enum Bgzf {
    Single(bgzf::io::Reader<File>),
    Multi(bgzf::io::MultithreadedReader<File>),
}

impl Bgzf {
    fn open(path: &Path, threads: usize) -> io::Result<Bgzf> {
        let file = File::open(path)?;
        Ok(match NonZero::new(threads).filter(|t| t.get() > 1) {
            Some(t) => Bgzf::Multi(bgzf::io::MultithreadedReader::with_worker_count(t, file)),
            None => Bgzf::Single(bgzf::io::Reader::new(file)),
        })
    }

    fn virtual_position(&self) -> VirtualPosition {
        match self {
            Bgzf::Single(r) => r.virtual_position(),
            Bgzf::Multi(r) => r.virtual_position(),
        }
    }

    fn seek(&mut self, pos: VirtualPosition) -> io::Result<()> {
        use pasta::bgzf::io::Seek as _;
        match self {
            Bgzf::Single(r) => r.seek(pos).map(|_| ()),
            Bgzf::Multi(r) => r.seek_to_virtual_position(pos).map(|_| ()),
        }
    }
}

impl Read for Bgzf {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Bgzf::Single(r) => r.read(buf),
            Bgzf::Multi(r) => r.read(buf),
        }
    }
}

/// Reads up to `buf.len()` bytes; Ok(false) at a clean end of file.
fn read_exact_or_eof(r: &mut impl Read, buf: &mut [u8]) -> io::Result<bool> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    if n == 0 {
        Ok(false)
    } else if n < buf.len() {
        Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated file"))
    } else {
        Ok(true)
    }
}

/// The region being iterated (htslib's `hts_itr_t` for BAM files).
enum Iter {
    /// All records from the current position to the end of the file
    Rest,
    /// Records of `tid` overlapping [beg, end), from the index chunks
    Region {
        tid: i32,
        beg: i64,
        end: i64,
        chunks: Vec<Chunk>,
        i: usize,
        curr_off: Option<VirtualPosition>,
    },
    /// Unplaced unmapped records (`HTS_IDX_NOCOOR`)
    Unmapped,
    Done,
}

pub struct BamReader {
    bgzf: Bgzf,
    pub header: Header,
    /// The header text as stored in the file (`l_text` bytes)
    pub raw_text: Vec<u8>,
    index: Option<Box<dyn BinningIndex + Send>>,
    first_record: VirtualPosition,
    iter: Iter,
    buf: Vec<u8>,
}

/// The index of `path` as htslib looks for it: `X.bam.csi`, `X.csi`,
/// `X.bam.bai`, `X.bai` (only the first two and `.tbi`/`.crai` variants for
/// other formats).
pub fn find_index(path: &Path, exts: &[&str]) -> Option<std::path::PathBuf> {
    let s = path.as_os_str().to_string_lossy();
    for ext in exts {
        let full = format!("{s}{ext}");
        if Path::new(&full).exists() {
            return Some(full.into());
        }
        if let Some(dot) = s.rfind('.')
            && dot > 0
        {
            let replaced = format!("{}{ext}", &s[..dot]);
            if Path::new(&replaced).exists() {
                return Some(replaced.into());
            }
        }
    }
    None
}

impl BamReader {
    pub fn open(path: &Path, threads: usize) -> Result<BamReader> {
        let mut bgzf = Bgzf::open(path, threads).with_context(|| format!("opening {}", path.display()))?;
        let (header, raw_text) =
            read_header(&mut bgzf).with_context(|| format!("reading the header of {}", path.display()))?;
        let first_record = bgzf.virtual_position();
        Ok(BamReader { bgzf, header, raw_text, index: None, first_record, iter: Iter::Rest, buf: Vec::new() })
    }

    /// Loads the BAI/CSI index (Ok(false) if there is none).
    pub fn load_index(&mut self, path: &Path) -> Result<bool> {
        let Some(index_path) = find_index(path, &[".csi", ".bai"]) else { return Ok(false) };
        let csi = index_path.extension().is_some_and(|e| e == "csi");
        self.index = Some(if csi {
            Box::new(pasta::csi::fs::read(&index_path)?)
        } else {
            Box::new(pasta::bam::bai::fs::read(&index_path)?)
        });
        Ok(true)
    }

    pub fn has_index(&self) -> bool {
        self.index.is_some()
    }

    /// Iterate over all records (`HTS_IDX_START`).
    pub fn fetch_all(&mut self) -> Result<()> {
        self.bgzf.seek(self.first_record)?;
        self.iter = Iter::Rest;
        Ok(())
    }

    /// Iterate over the records of `tid` overlapping [beg, end).
    pub fn fetch(&mut self, tid: usize, beg: i64, end: i64) -> Result<()> {
        let Some(index) = &self.index else { bail!("no index available") };
        let max = 1i64 << (index.min_shift() as u32 + 3 * index.depth() as u32);
        let chunks = if beg >= max || end <= beg {
            Vec::new()
        } else {
            let start = Position::try_from((beg + 1) as usize)?;
            let stop = Position::try_from(end.min(max) as usize)?;
            index.query(tid, (start..=stop).into())?
        };
        self.iter = Iter::Region { tid: tid as i32, beg, end, chunks, i: 0, curr_off: None };
        Ok(())
    }

    /// Iterate over the unplaced unmapped records at the end of the file.
    pub fn fetch_unmapped(&mut self) -> Result<()> {
        let Some(index) = &self.index else { bail!("no index available") };
        let start = index.last_first_record_start_position().unwrap_or(self.first_record);
        self.bgzf.seek(start)?;
        self.iter = Iter::Unmapped;
        Ok(())
    }

    fn read_record(&mut self) -> Result<Option<Record>> {
        let mut len = [0u8; 4];
        if !read_exact_or_eof(&mut self.bgzf, &mut len)? {
            return Ok(None);
        }
        let n = i32::from_le_bytes(len);
        if n < 32 {
            bail!("invalid BAM record (block size {n})");
        }
        self.buf.resize(n as usize, 0);
        if !read_exact_or_eof(&mut self.bgzf, &mut self.buf)? {
            bail!("truncated file");
        }
        Ok(Some(Record::from_bam(&self.buf)?))
    }

    /// The next record of the current iteration (`hts_itr_next`).
    pub fn read(&mut self) -> Result<Option<Record>> {
        loop {
            match &mut self.iter {
                Iter::Done => return Ok(None),
                Iter::Rest => return self.read_record(),
                Iter::Unmapped => match self.read_record()? {
                    None => return Ok(None),
                    Some(r) if r.tid < 0 => return Ok(Some(r)),
                    Some(_) => {}
                },
                Iter::Region { chunks, i, curr_off, .. } => {
                    let at_chunk_end = match *curr_off {
                        None => true,
                        Some(off) => off >= chunks[*i].end(),
                    };
                    if at_chunk_end {
                        let next = if curr_off.is_none() { 0 } else { *i + 1 };
                        if next >= chunks.len() {
                            self.iter = Iter::Done;
                            continue;
                        }
                        let adjacent = curr_off.is_some() && chunks[*i].end() == chunks[next].start();
                        let start = chunks[next].start();
                        *i = next;
                        if !adjacent {
                            self.bgzf.seek(start)?;
                        }
                        *curr_off = Some(self.bgzf.virtual_position());
                    }
                    let rec = self.read_record()?;
                    let Iter::Region { tid, beg, end, curr_off, .. } = &mut self.iter else { unreachable!() };
                    let Some(rec) = rec else {
                        self.iter = Iter::Done;
                        return Ok(None);
                    };
                    *curr_off = Some(self.bgzf.virtual_position());
                    let (rbeg, rend) = (rec.pos(), rec.end_pos());
                    if rec.tid != *tid || rbeg >= *end {
                        self.iter = Iter::Done;
                        return Ok(None);
                    }
                    if rend > *beg && *end > rbeg {
                        return Ok(Some(rec));
                    }
                }
            }
        }
    }
}

fn read_i32(r: &mut impl Read) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

/// The parsed header and the raw header text
fn read_header(r: &mut impl Read) -> Result<(Header, Vec<u8>)> {
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != b"BAM\x01" {
        bail!("invalid BAM header");
    }
    let l_text = read_i32(r)?;
    let mut text = vec![0u8; l_text.max(0) as usize];
    r.read_exact(&mut text)?;
    let raw_text = text.clone();
    // htslib stops at the first NUL
    if let Some(end) = text.iter().position(|&c| c == 0) {
        text.truncate(end);
    }
    let n_ref = read_i32(r)?;
    let mut refs = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref {
        let l_name = read_i32(r)?;
        let mut name = vec![0u8; l_name.max(0) as usize];
        r.read_exact(&mut name)?;
        if let Some(end) = name.iter().position(|&c| c == 0) {
            name.truncate(end);
        }
        let len = read_i32(r)? as u32;
        refs.push((String::from_utf8_lossy(&name).into_owned(), len));
    }
    Ok((Header { text: String::from_utf8_lossy(&text).into_owned(), refs }, raw_text))
}

/// htslib's `sam_hdr_sanitise` on the text of a BAM header: rejects lines
/// that do not start with '@' and adds a missing final newline. NUL padding
/// is kept.
pub fn sanitise_header_text(text: &[u8]) -> Result<Vec<u8>> {
    let mut text = text.to_vec();
    if text.is_empty() {
        return Ok(text);
    }
    let mut lnum = 0;
    let mut last = b'\n';
    let mut i = 0;
    while i < text.len() && text[i] != 0 {
        if last == b'\n' {
            lnum += 1;
            if text[i] != b'@' {
                bail!("Malformed SAM header at line {lnum}");
            }
        }
        last = text[i];
        i += 1;
    }
    if last != b'\n' {
        if i < text.len() {
            text[i] = b'\n';
        } else {
            text.push(b'\n');
        }
    }
    Ok(text)
}

/// A BGZF writer, compressing with worker threads if requested.
enum BgzfWriter {
    Single(bgzf::io::Writer<Box<dyn Write + Send>>),
    Multi(bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Write for BgzfWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            BgzfWriter::Single(w) => w.write(buf),
            BgzfWriter::Multi(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            BgzfWriter::Single(w) => w.flush(),
            BgzfWriter::Multi(w) => w.flush(),
        }
    }
}

/// A BAM writer (`bam_hdr_write` + `bam_write1`).
pub struct BamWriter {
    inner: BgzfWriter,
    buf: Vec<u8>,
}

impl BamWriter {
    /// `level`: compression level 0-9 (None: the default, 6 as in htslib).
    pub fn create(
        out: Box<dyn Write + Send>,
        header: &Header,
        level: Option<u32>,
        threads: usize,
    ) -> Result<BamWriter> {
        Self::create_with_text(out, header.text.as_bytes(), &header.refs, level, threads)
    }

    /// Like `create`, with the header text given as bytes (written as is).
    pub fn create_with_text(
        out: Box<dyn Write + Send>,
        text: &[u8],
        refs: &[(String, u32)],
        level: Option<u32>,
        threads: usize,
    ) -> Result<BamWriter> {
        let level = bgzf::io::writer::CompressionLevel::new(level.unwrap_or(6) as u8)
            .ok_or_else(|| anyhow::anyhow!("invalid compression level"))?;
        let inner = match NonZero::new(threads).filter(|t| t.get() > 1) {
            Some(t) => BgzfWriter::Multi(
                bgzf::io::multithreaded_writer::Builder::default()
                    .set_compression_level(level)
                    .set_worker_count(t)
                    .build_from_writer(out),
            ),
            None => BgzfWriter::Single(
                bgzf::io::writer::Builder::default().set_compression_level(level).build_from_writer(out),
            ),
        };
        let mut w = BamWriter { inner, buf: Vec::new() };
        let mut h = Vec::new();
        h.extend_from_slice(b"BAM\x01");
        h.extend_from_slice(&(text.len() as i32).to_le_bytes());
        h.extend_from_slice(text);
        h.extend_from_slice(&(refs.len() as i32).to_le_bytes());
        for (name, len) in refs {
            h.extend_from_slice(&(name.len() as i32 + 1).to_le_bytes());
            h.extend_from_slice(name.as_bytes());
            h.push(0);
            h.extend_from_slice(&len.to_le_bytes());
        }
        w.inner.write_all(&h)?;
        Ok(w)
    }

    pub fn write(&mut self, rec: &Record) -> Result<()> {
        self.buf.clear();
        rec.write_bam(&mut self.buf)?;
        self.inner.write_all(&self.buf)?;
        Ok(())
    }

    /// Writes the remaining data and the end-of-file block.
    pub fn finish(self) -> Result<()> {
        let mut out = match self.inner {
            BgzfWriter::Single(w) => w.finish()?,
            BgzfWriter::Multi(mut w) => w.finish()?,
        };
        out.flush()?;
        Ok(())
    }
}

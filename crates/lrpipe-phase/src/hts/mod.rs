//! VCF/BCF reading and writing as WhatsHap does it through pysam, which uses
//! htslib 1.23: header and record handling are ports of htslib's (`header`,
//! `record`), compression and indexes come from pasta. The API mirrors the
//! pysam calls WhatsHap makes.

mod header;
pub mod num;
mod record;

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use pasta::bgzf::{self, VirtualPosition};
use pasta::core::Position;
use pasta::csi::BinningIndex;
use pasta::csi::binning_index::index::reference_sequence::bin::Chunk;

pub use header::{
    BCF_DT_CTG, BCF_DT_ID, BCF_DT_SAMPLE, BCF_HL_FLT, BCF_HL_FMT, BCF_HL_GEN, BCF_HL_INFO, BCF_HT_FLAG, BCF_HT_INT,
    BCF_HT_REAL, BCF_HT_STR, BCF_VL_FIXED, BCF_VL_VAR,
};
use header::{BCF_VL_A, BCF_VL_G, BCF_VL_R, HeaderCore, Hrec};
pub use record::{
    BCF_BT_CHAR, BCF_BT_FLOAT, BCF_BT_INT8, BCF_BT_INT16, BCF_BT_INT32, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING,
    INT16_VECTOR_END, INT32_MISSING, INT32_VECTOR_END, Record,
};
use record::{BCF_ERR_LIMITS, BcfCore, FLOAT_MISSING, FLOAT_VECTOR_END, FmtValues, InfoValues};

use crate::alignment::bam::find_index;

/// Metadata of a header ID for one header line type (INFO/FORMAT).
#[derive(Clone, Debug, PartialEq)]
pub struct IdInfo {
    pub number: String,
    pub typ: String,
}

#[derive(Clone, Debug)]
pub struct Header {
    pub(crate) core: HeaderCore,
}

impl Header {
    pub fn dup(&self) -> Header {
        Header { core: self.core.dup().expect("a parsed header can be parsed again") }
    }

    pub fn samples(&self) -> Vec<String> {
        self.core.samples.clone()
    }

    pub fn contigs(&self) -> Vec<String> {
        self.core.ids(BCF_DT_CTG).iter().flatten().cloned().collect()
    }

    /// Contigs with their lengths as pysam's `VariantContig.length` reports
    /// them (the low 32 bits of the header's length; 0 means unknown).
    pub fn contig_lengths(&self) -> Vec<(String, u32)> {
        let ids = self.core.ids(BCF_DT_CTG);
        (0..ids.len())
            .filter_map(|rid| ids[rid].clone().map(|name| (name, self.core.contig_length(rid as i32) as u32)))
            .collect()
    }

    /// IDs defined for the given header line type (BCF_HL_INFO or BCF_HL_FMT),
    /// in dictionary order, with Number and Type as pysam reports them.
    pub fn defined(&self, hl: i32) -> Vec<(String, IdInfo)> {
        let mut out = Vec::new();
        for (id, key) in self.core.ids(BCF_DT_ID).iter().enumerate() {
            let Some(key) = key else { continue };
            let info = self.core.info(hl, id as i32);
            if info & 0xf == 0xf {
                continue;
            }
            let number = match (info >> 8) & 0xf {
                BCF_VL_FIXED => (info >> 12).to_string(),
                BCF_VL_VAR => ".".to_string(),
                2 => "A".to_string(),
                3 => "G".to_string(),
                4 => "R".to_string(),
                _ => "?".to_string(),
            };
            let typ = match (info >> 4) & 0xf {
                BCF_HT_FLAG => "Flag",
                BCF_HT_INT => "Integer",
                BCF_HT_REAL => "Float",
                BCF_HT_STR => "String",
                _ => "?",
            }
            .to_string();
            out.push((key.clone(), IdInfo { number, typ }));
        }
        out
    }

    pub fn id(&self, key: &str) -> i32 {
        self.core.id2int(BCF_DT_ID, key)
    }

    pub fn is_defined(&self, hl: i32, key: &str) -> bool {
        let id = self.id(key);
        id >= 0 && self.core.info(hl, id) & 0xf != 0xf
    }

    pub fn id2type(&self, hl: i32, id: i32) -> u32 {
        self.core.id2type(hl, id)
    }

    fn sync(&mut self) {
        if self.core.dirty {
            self.core.sync();
        }
    }

    /// pysam `VariantHeader.add_line`
    pub fn add_line(&mut self, line: &str) -> Result<()> {
        self.core.append(line)?;
        self.sync();
        Ok(())
    }

    /// pysam `VariantHeader.add_meta(key, value=...)`
    pub fn add_meta_value(&mut self, key: &str, value: &str) -> Result<()> {
        let mut hrec = Hrec::new(key);
        hrec.value = Some(value.to_string());
        self.core.add_hrec(hrec)?;
        self.sync();
        Ok(())
    }

    /// pysam `VariantHeader.add_meta(key, items=[(k, v, quoted), ...])`
    pub fn add_meta_items(&mut self, key: &str, items: &[(&str, &str, bool)]) -> Result<()> {
        let mut hrec = Hrec::new(key);
        for (i, (k, v, quoted)) in items.iter().enumerate() {
            hrec.add_key(k);
            hrec.set_val(i, v, *quoted);
        }
        self.core.add_hrec(hrec)?;
        self.sync();
        Ok(())
    }

    /// `bcf_hdr_remove` (pysam `remove_header()` / `VariantHeaderRecord.remove()`)
    pub fn remove(&mut self, hl: i32, key: &str) {
        self.core.remove(hl, key);
    }

    /// pysam's `for hr in header.records: if hr.key == key: hr.remove(); break`:
    /// a generic line is removed with all generic lines of that key, a
    /// structured one (whose value pysam passes as a NULL key) with all
    /// structured lines.
    pub fn remove_first_record_with_key(&mut self, key: &str) {
        let Some(h) = self.core.hrecs.iter().find(|h| h.key == key) else { return };
        if h.typ == BCF_HL_GEN {
            self.core.remove(BCF_HL_GEN, key);
        } else {
            self.core.remove_all(h.typ);
        }
    }

    /// Keys of all header records (pysam `header.records`, `.key`).
    pub fn record_keys(&self) -> Vec<String> {
        self.core.hrecs.iter().map(|h| h.key.clone()).collect()
    }

    pub fn nsamples(&self) -> usize {
        self.core.samples.len()
    }
}

// --- input ------------------------------------------------------------------

/// A decompressed input stream (htslib detects the compression).
enum Input {
    Plain(BufReader<File>),
    Bgzf(bgzf::io::Reader<File>),
    Gzip(Box<BufReader<ReadFull<flate2::read::MultiGzDecoder<File>>>>),
    /// standard input ("-"), decompressed if it is gzip/BGZF
    Stream(Box<dyn BufRead + Send>),
}

impl Read for Input {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Input::Plain(r) => r.read(buf),
            Input::Bgzf(r) => r.read(buf),
            Input::Gzip(r) => r.read(buf),
            Input::Stream(r) => r.read(buf),
        }
    }
}

impl BufRead for Input {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Input::Plain(r) => r.fill_buf(),
            Input::Bgzf(r) => r.fill_buf(),
            Input::Gzip(r) => r.fill_buf(),
            Input::Stream(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, n: usize) {
        match self {
            Input::Plain(r) => r.consume(n),
            Input::Bgzf(r) => r.consume(n),
            Input::Gzip(r) => r.consume(n),
            Input::Stream(r) => r.consume(n),
        }
    }
}

impl Input {
    fn open(path: &str) -> Result<Input> {
        if path == "-" {
            let mut stdin = BufReader::new(ReadFull(io::stdin()));
            let gzip = stdin.fill_buf()?.starts_with(b"\x1f\x8b");
            return Ok(if gzip {
                Input::Stream(Box::new(BufReader::new(ReadFull(flate2::bufread::MultiGzDecoder::new(stdin)))))
            } else {
                Input::Stream(Box::new(stdin))
            });
        }
        let mut file = File::open(path).with_context(|| format!("could not open variant file {path}"))?;
        let mut magic = [0u8; 16];
        let n = read_up_to(&mut file, &mut magic)?;
        drop(file);
        let file = File::open(path)?;
        let gzip = n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b;
        // BGZF: FEXTRA with a "BC" subfield first
        let bgzf = gzip && n >= 14 && magic[3] & 4 != 0 && magic[12] == b'B' && magic[13] == b'C';
        Ok(if bgzf {
            Input::Bgzf(bgzf::io::Reader::new(file))
        } else if gzip {
            Input::Gzip(Box::new(BufReader::new(ReadFull(flate2::read::MultiGzDecoder::new(file)))))
        } else {
            Input::Plain(BufReader::new(file))
        })
    }

    /// `hts_getline`: the next line without its newline (and a `\r` before
    /// it); Ok(false) at the end of the file.
    fn getline(&mut self, line: &mut Vec<u8>) -> io::Result<bool> {
        getline(self, line)
    }
}

fn getline(r: &mut impl BufRead, line: &mut Vec<u8>) -> io::Result<bool> {
    line.clear();
    if r.read_until(b'\n', line)? == 0 {
        return Ok(false);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.len() > 1 && line.last() == Some(&b'\r') {
        line.pop();
    }
    Ok(true)
}

/// A reader whose reads fill the buffer unless the stream ends, so that a
/// `BufReader` over it sees enough bytes in one `fill_buf` to detect the format
/// (pipes and decompressors can return short reads).
struct ReadFull<R>(R);

impl<R: Read> Read for ReadFull<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        read_up_to(&mut self.0, buf)
    }
}

fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// pysam's error check of a record that was read (`makeVariantRecord`).
fn check_record_errcode(rec: &Record) -> Result<()> {
    let mut msg = Vec::new();
    for (bit, text) in [
        (4, "invalid number of columns"),
        (8, "limits violated"),
        (16, "invalid character found"),
        (32, "invalid contig"),
        (64, "invalid tag"),
    ] {
        if rec.errcode & bit != 0 {
            msg.push(text);
        }
    }
    if !msg.is_empty() {
        bail!("Error(s) reading record: {}", msg.join(", "));
    }
    Ok(())
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Format {
    Vcf,
    Bcf,
}

/// The index of a VCF/BCF file, as htslib looks for it.
fn load_index(path: &str, bcf: bool) -> Result<Option<Box<dyn BinningIndex + Send>>> {
    let exts: &[&str] = if bcf { &[".csi"] } else { &[".csi", ".tbi"] };
    let Some(index_path) = find_index(Path::new(path), exts) else { return Ok(None) };
    let index: Box<dyn BinningIndex + Send> = if index_path.extension().is_some_and(|e| e == "csi") {
        Box::new(pasta::csi::fs::read(&index_path)?)
    } else {
        Box::new(pasta::tabix::fs::read(&index_path)?)
    };
    Ok(Some(index))
}

/// Sequence names of a tabix (or CSI with tabix metadata) index.
fn index_names(index: &dyn BinningIndex) -> Vec<String> {
    index
        .header()
        .map(|h| h.reference_sequence_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect())
        .unwrap_or_default()
}

pub struct VcfFile {
    input: Input,
    format: Format,
    pub header: Header,
    line: Vec<u8>,
    buf: Vec<u8>,
}

impl VcfFile {
    pub fn open(path: &str) -> Result<VcfFile> {
        let mut input = Input::open(path)?;
        let is_bcf = input.fill_buf()?.starts_with(b"BCF");
        let mut header = HeaderCore::new();
        let mut line = Vec::new();
        if is_bcf {
            let mut magic = [0u8; 5];
            input.read_exact(&mut magic).context("Failed to read the header (reading BCF in text mode?)")?;
            if &magic != b"BCF\x02\x02" {
                bail!("Invalid BCF2 magic string: only BCFv2.2 is supported");
            }
            let mut len = [0u8; 4];
            input.read_exact(&mut len)?;
            let mut text = vec![0u8; u32::from_le_bytes(len) as usize];
            input.read_exact(&mut text).context("Failed to read BCF header")?;
            let end = text.iter().position(|&c| c == 0).unwrap_or(text.len());
            header.parse(&text[..end]).context("Failed to read BCF header")?;
        } else {
            // vcf_hdr_read
            let mut text = Vec::new();
            while input.getline(&mut line)? {
                if line.is_empty() {
                    continue;
                }
                if line[0] != b'#' {
                    bail!("No sample line in {path}");
                }
                text.extend_from_slice(&line);
                text.push(b'\n');
                if line.get(1) != Some(&b'#') {
                    break;
                }
            }
            if text.is_empty() {
                bail!("could not read VCF header from {path}");
            }
            header.parse(&text).with_context(|| format!("could not read VCF header from {path}"))?;
            // contigs in the tabix index but not in the header
            if let Ok(Some(index)) = load_index(path, false) {
                for name in index_names(index.as_ref()) {
                    if !header.contig_has_hrec(&name) {
                        header.add_missing_contig(&name)?;
                    }
                }
                if header.dirty {
                    header.sync();
                }
            }
        }
        // pysam's VariantFile.open calls tell() after the header, which fails
        // for gzip that is not BGZF
        if matches!(input, Input::Gzip(_)) {
            bail!("seek not implemented in files compressed by method 1");
        }
        let format = if is_bcf { Format::Bcf } else { Format::Vcf };
        Ok(VcfFile { input, format, header: Header { core: header }, line, buf: Vec::new() })
    }

    /// The next BCF record: Ok(None) at the end of the file.
    fn read_bcf(&mut self) -> Result<Option<Record>> {
        let mut x = [0u8; 32];
        let n = read_up_to(&mut self.input, &mut x)?;
        if n == 0 {
            return Ok(None);
        }
        if n < 32 {
            bail!("truncated file");
        }
        let u32_at = |i: usize| u32::from_le_bytes(x[i..i + 4].try_into().unwrap());
        let shared_len = u32_at(0);
        if shared_len < 24 {
            bail!("truncated file");
        }
        let shared_len = (shared_len - 24) as usize;
        let indiv_len = u32_at(4) as usize;
        let pos = u32_at(12);
        let mut core = BcfCore {
            rid: u32_at(8) as i32,
            pos: if pos == u32::MAX { -1 } else { pos as i64 },
            rlen: u32_at(16) as i32 as i64,
            qual: u32_at(20),
            n_info: u16::from_le_bytes([x[24], x[25]]) as usize,
            n_allele: u16::from_le_bytes([x[26], x[27]]) as usize,
            n_sample: u32_at(28) & 0xff_ffff,
            n_fmt: x[31] as usize,
        };
        if indiv_len == 0 || core.n_sample == 0 {
            core.n_fmt = 0;
        }
        self.buf.resize(shared_len + indiv_len, 0);
        if read_up_to(&mut self.input, &mut self.buf)? != self.buf.len() {
            bail!("truncated file");
        }
        let (shared, indiv) = self.buf.split_at(shared_len);
        let mut rec = Record::from_bcf(&core, shared, indiv).map_err(|_| anyhow!("unable to parse next record"))?;
        if rec.check_bcf(&self.header.core).is_err() {
            bail!("unable to parse next record");
        }
        Ok(Some(rec))
    }

    /// Read next record; Ok(false) at EOF.
    pub fn read(&mut self, rec: &mut Record) -> Result<bool> {
        // Same error handling as pysam's VariantFile.__next__ and
        // makeVariantRecord: undefined contigs/tags are tolerated.
        match self.format {
            Format::Bcf => match self.read_bcf()? {
                None => return Ok(false),
                Some(r) => *rec = r,
            },
            Format::Vcf => {
                if !self.input.getline(&mut self.line)? {
                    return Ok(false);
                }
                let (r, ret) = Record::parse_vcf(&self.line, &mut self.header.core);
                if ret < 0 {
                    if r.errcode != 0 {
                        bail!("unable to parse next record");
                    }
                    bail!("truncated file");
                }
                *rec = r;
            }
        }
        check_record_errcode(rec)?;
        Ok(true)
    }
}

// --- output -----------------------------------------------------------------

/// BGZF block size of htslib (`BGZF_BLOCK_SIZE`)
const BGZF_BLOCK_SIZE: usize = 0xff00;

enum Output {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Bgzf(bgzf::io::Writer<Box<dyn Write + Send>>),
}

impl Output {
    fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        match self {
            Output::Plain(w) => w.write_all(data),
            Output::Bgzf(w) => w.write_all(data),
        }
    }

    /// `bgzf_flush_try`: starts a new block unless `len` more bytes fit.
    fn flush_try(&mut self, len: usize) -> io::Result<()> {
        if let Output::Bgzf(w) = self
            && usize::from(w.virtual_position().uncompressed()) + len > BGZF_BLOCK_SIZE
        {
            w.flush()?;
        }
        Ok(())
    }

    /// `bgzf_flush`
    fn flush_block(&mut self) -> io::Result<()> {
        match self {
            Output::Plain(_) => Ok(()),
            Output::Bgzf(w) => w.flush(),
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Output::Plain(mut w) => w.flush(),
            Output::Bgzf(w) => w.finish()?.flush(),
        }
    }
}

pub struct VcfWriter {
    out: Option<Output>,
    format: Format,
    pub header: Header,
    header_written: bool,
    line: Vec<u8>,
}

impl VcfWriter {
    /// Like pysam `VariantFile(path, "w", header=header)`; "-" is stdout.
    pub fn create(path: &str, header: &Header) -> Result<VcfWriter> {
        Self::create_as(path, header, false)
    }

    /// Like pysam `VariantFile(open(path, "w"), "w", header=header)`:
    /// uncompressed VCF whatever the file name.
    pub fn create_plain(path: &str, header: &Header) -> Result<VcfWriter> {
        Self::create_as(path, header, true)
    }

    fn create_as(path: &str, header: &Header, plain: bool) -> Result<VcfWriter> {
        let (format, compressed) = if plain {
            (Format::Vcf, false)
        } else if path.ends_with(".gz") {
            (Format::Vcf, true)
        } else if path.ends_with(".bcf") {
            (Format::Bcf, true)
        } else {
            (Format::Vcf, false)
        };
        let sink: Box<dyn Write + Send> = if path == "-" {
            Box::new(io::stdout())
        } else {
            Box::new(File::create(path).with_context(|| format!("could not open {path} for writing"))?)
        };
        let out = if compressed {
            // htslib's default level (zlib's 6)
            let level = bgzf::io::writer::CompressionLevel::new(6).expect("valid level");
            Output::Bgzf(bgzf::io::writer::Builder::default().set_compression_level(level).build_from_writer(sink))
        } else {
            Output::Plain(BufWriter::with_capacity(1 << 16, sink))
        };
        Ok(VcfWriter { out: Some(out), format, header: header.dup(), header_written: false, line: Vec::new() })
    }

    fn out(&mut self) -> &mut Output {
        self.out.as_mut().expect("writer is open")
    }

    fn write_header(&mut self) -> Result<()> {
        if self.header_written {
            return Ok(());
        }
        self.header_written = true;
        self.header.sync();
        match self.format {
            Format::Vcf => {
                let text = self.header.core.format(false);
                let out = self.out();
                out.write_all(text.as_bytes())?;
                out.flush_block()?;
            }
            Format::Bcf => {
                let mut text = self.header.core.format(true).into_bytes();
                text.push(0);
                let out = self.out();
                out.write_all(b"BCF\x02\x02")?;
                out.write_all(&(text.len() as u32).to_le_bytes())?;
                out.write_all(&text)?;
                out.flush_block()?;
            }
        }
        Ok(())
    }

    /// pysam `VariantFile.write(record)`; `record_header` is the header the
    /// record was read with.
    pub fn write(&mut self, record: &mut Record, record_header: &Header) -> Result<()> {
        self.write_header()?;
        if record.n_sample as usize != self.header.nsamples() {
            bail!("Invalid VariantRecord. Number of samples does not match header");
        }
        record.sync_end(record_header)?;
        // bcf_write
        self.header.sync();
        let mut line = std::mem::take(&mut self.line);
        line.clear();
        match self.format {
            Format::Vcf => record.format_vcf(&self.header.core, &mut line).context("Can't write record")?,
            Format::Bcf => {
                if record.errcode & !BCF_ERR_LIMITS != 0 {
                    bail!("Can't write record: unchecked error {}", record.errcode);
                }
                record.encode_bcf(&mut line).context("Can't write record")?;
            }
        }
        let vcf = self.format == Format::Vcf;
        let out = self.out();
        if vcf {
            out.flush_try(line.len())?;
        }
        out.write_all(&line)?;
        self.line = line;
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        self.write_header()?;
        self.out.take().expect("writer is open").finish().context("error closing output VCF")?;
        Ok(())
    }
}

impl Drop for VcfWriter {
    fn drop(&mut self) {
        if let Some(out) = self.out.take() {
            let _ = out.finish();
        }
    }
}

// --- region queries -----------------------------------------------------------

/// A Python value assigned to a FORMAT field of one sample
#[derive(Clone, Debug)]
pub enum PyFormatValue {
    Ints(Vec<Option<i64>>),
    Floats(Vec<Option<f64>>),
    /// GT: allele indices (None: missing)
    Alleles(Vec<Option<u32>>),
    /// n times None (a scalar None is 1)
    Nulls(usize),
}

/// Errors from `IndexedVcf::fetch`, mirroring the exceptions pysam raises.
#[derive(Debug)]
pub enum VcfFetchError {
    /// "fetch requires an index"
    IndexMissing,
    /// "invalid contig" (contig not in the VCF header)
    InvalidContig,
    Other(anyhow::Error),
}

/// A VCF/BCF opened for region queries, like pysam's `VariantFile.fetch`.
pub struct IndexedVcf {
    file: VcfFile,
    index: Option<Box<dyn BinningIndex + Send>>,
    /// a second reader of the file for the queries
    bgzf: Option<bgzf::io::Reader<File>>,
    /// pysam's `index.refmap`: sequence name -> position in the index's
    /// sequence-name list (used as the tid of the query, as pysam does)
    refmap: HashMap<String, usize>,
    /// the index's sequence names (tabix: `tbx_name2id`)
    tbx_names: HashMap<String, usize>,
    line: Vec<u8>,
}

impl IndexedVcf {
    pub fn open(path: &str) -> Result<IndexedVcf> {
        let file = VcfFile::open(path)?;
        let mut refmap = HashMap::new();
        let mut tbx_names = HashMap::new();
        let mut index = None;
        let mut reader = None;
        let bgzf = matches!(file.input, Input::Bgzf(_));
        if bgzf && let Some(idx) = load_index(path, file.format == Format::Bcf)? {
            if file.format == Format::Bcf {
                // hts_idx_seqnames with bcf_hdr_id2name: every sequence of the index
                let n = idx.reference_sequences().count();
                for i in 0..n {
                    if let Some(name) = file.header.core.id_key(BCF_DT_CTG, i as i32) {
                        refmap.insert(name.to_string(), i);
                    }
                }
            } else {
                for (i, name) in index_names(idx.as_ref()).into_iter().enumerate() {
                    refmap.insert(name.clone(), i);
                    tbx_names.insert(name, i);
                }
            }
            index = Some(idx);
            reader = Some(bgzf::io::Reader::new(File::open(path)?));
        }
        Ok(IndexedVcf { file, index, bgzf: reader, refmap, tbx_names, line: Vec::new() })
    }

    pub fn header(&self) -> &Header {
        &self.file.header
    }

    /// Whether an index was found (pysam's `VariantFile.index is not None`).
    pub fn has_index(&self) -> bool {
        self.index.is_some()
    }

    /// Records overlapping [start, stop) on `contig` (stop None: to the end).
    pub fn fetch(&mut self, contig: &str, start: i64, stop: Option<i64>) -> Result<Vec<Record>, VcfFetchError> {
        let Some(index) = &self.index else { return Err(VcfFetchError::IndexMissing) };
        if self.file.header.core.id2int(BCF_DT_CTG, contig) < 0 {
            return Err(VcfFetchError::InvalidContig);
        }
        let Some(&tid) = self.refmap.get(contig) else { return Ok(Vec::new()) };
        // MAX_POS in pysam
        let stop = stop.unwrap_or((1 << 31) - 1);
        let beg = start.max(0);
        let max = 1i64 << (index.min_shift() as u32 + 3 * index.depth() as u32);
        let chunks = if beg >= max || stop <= beg {
            Vec::new()
        } else {
            let query = || -> Result<Vec<Chunk>> {
                let s = Position::try_from((beg + 1) as usize)?;
                let e = Position::try_from(stop.min(max) as usize)?;
                Ok(index.query(tid, (s..=e).into())?)
            };
            query().map_err(|e| {
                VcfFetchError::Other(e.context(format!("unable to fetch {contig}:{}-{stop}", start + 1)))
            })?
        };
        self.read_chunks(tid as i32, beg, stop, &chunks).map_err(VcfFetchError::Other)
    }

    /// `hts_itr_next` over the chunks, keeping records of `tid` overlapping
    /// [beg, end).
    fn read_chunks(&mut self, tid: i32, beg: i64, end: i64, chunks: &[Chunk]) -> Result<Vec<Record>> {
        let mut out = Vec::new();
        let bgzf = self.bgzf.as_mut().expect("indexed files have a reader");
        let mut i = 0;
        let mut curr_off: Option<VirtualPosition> = None;
        loop {
            let at_chunk_end = curr_off.is_none_or(|off| off >= chunks[i].end());
            if at_chunk_end {
                let next = if curr_off.is_none() { 0 } else { i + 1 };
                if next >= chunks.len() {
                    break;
                }
                let adjacent = curr_off.is_some() && chunks[i].end() == chunks[next].start();
                i = next;
                if !adjacent {
                    bgzf.seek(chunks[i].start())?;
                }
            }
            let (rec, rtid, rbeg, rend) = match self.file.format {
                Format::Vcf => {
                    // tbx_readrec: skip meta lines
                    loop {
                        if !getline(bgzf, &mut self.line)? {
                            return Ok(out);
                        }
                        if self.line.first() != Some(&b'#') {
                            break;
                        }
                    }
                    let Some((name, rbeg, rend)) = tbx_parse_vcf(&self.line) else {
                        bail!("unable to fetch next record");
                    };
                    let Some(&rtid) = self.tbx_names.get(&name) else { bail!("unable to fetch next record") };
                    (None, rtid as i32, rbeg, rend)
                }
                Format::Bcf => {
                    let rec = read_bcf_record(bgzf, &self.file.header.core)?;
                    let Some(rec) = rec else { return Ok(out) };
                    let (t, b, e) = (rec.rid, rec.pos, rec.pos + rec.rlen);
                    (Some(rec), t, b, e)
                }
            };
            curr_off = Some(bgzf.virtual_position());
            if rtid != tid || rbeg >= end {
                break;
            }
            if rend > beg && end > rbeg {
                let rec = match rec {
                    Some(r) => r,
                    None => {
                        let (r, ret) = Record::parse_vcf(&self.line, &mut self.file.header.core);
                        if ret < 0 {
                            bail!("error in vcf_parse");
                        }
                        r
                    }
                };
                check_record_errcode(&rec)?;
                out.push(rec);
            }
        }
        Ok(out)
    }
}

/// A BCF record from an indexed reader (`bcf_readrec`).
fn read_bcf_record(r: &mut bgzf::io::Reader<File>, h: &HeaderCore) -> Result<Option<Record>> {
    let mut x = [0u8; 32];
    let n = read_up_to(r, &mut x)?;
    if n == 0 {
        return Ok(None);
    }
    if n < 32 {
        bail!("truncated file");
    }
    let u32_at = |i: usize| u32::from_le_bytes(x[i..i + 4].try_into().unwrap());
    let shared_len = u32_at(0).checked_sub(24).ok_or_else(|| anyhow!("truncated file"))? as usize;
    let indiv_len = u32_at(4) as usize;
    let pos = u32_at(12);
    let mut core = BcfCore {
        rid: u32_at(8) as i32,
        pos: if pos == u32::MAX { -1 } else { pos as i64 },
        rlen: u32_at(16) as i32 as i64,
        qual: u32_at(20),
        n_info: u16::from_le_bytes([x[24], x[25]]) as usize,
        n_allele: u16::from_le_bytes([x[26], x[27]]) as usize,
        n_sample: u32_at(28) & 0xff_ffff,
        n_fmt: x[31] as usize,
    };
    if indiv_len == 0 || core.n_sample == 0 {
        core.n_fmt = 0;
    }
    let mut buf = vec![0u8; shared_len + indiv_len];
    if read_up_to(r, &mut buf)? != buf.len() {
        bail!("truncated file");
    }
    let (shared, indiv) = buf.split_at(shared_len);
    let mut rec = Record::from_bcf(&core, shared, indiv).map_err(|_| anyhow!("unable to fetch next record"))?;
    rec.check_bcf(h).map_err(|_| anyhow!("unable to fetch next record"))?;
    Ok(Some(rec))
}

/// C `strtoll(s, &end, base)` for base 10 or 0: (value, bytes used).
fn strtoll(s: &[u8], base0: bool) -> (i64, usize) {
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let mut i = 0;
    while matches!(at(i), b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r') {
        i += 1;
    }
    let neg = at(i) == b'-';
    if matches!(at(i), b'-' | b'+') {
        i += 1;
    }
    let mut base = 10u32;
    if base0 && at(i) == b'0' {
        if matches!(at(i + 1), b'x' | b'X') && at(i + 2).is_ascii_hexdigit() {
            base = 16;
            i += 2;
        } else {
            base = 8;
        }
    }
    let start = i;
    let mut v: i128 = 0;
    while let Some(d) = (at(i) as char).to_digit(base) {
        v = (v * base as i128 + d as i128).min(1 << 64);
        i += 1;
    }
    if i == start {
        return (0, 0);
    }
    let v = if neg { -v } else { v };
    (v.clamp(i64::MIN as i128, i64::MAX as i128) as i64, i)
}

/// `tbx_parse1` for the VCF preset (htslib 1.23): contig name, 0-based
/// begin and end. The end is the longest of REF, INFO/END (if after the
/// begin), the largest |SVLEN| of <DEL>/<DUP>/<CNV>/<INV> alleles and, with
/// gVCF alleles, the largest FORMAT/LEN. Quirks kept: SVLEN values are read
/// up to the end of INFO (commas of later fields count), END=/SVLEN= are
/// found by substring search, and LEN is only looked for if FORMAT has it.
fn tbx_parse_vcf(line: &[u8]) -> Option<(String, i64, i64)> {
    let mut name = None;
    let (mut beg, mut end) = (-1i64, -1i64);
    let mut svlenals: Vec<bool> = Vec::new();
    let (mut use_svlen, mut getlen) = (false, false);
    let (mut reflen, mut svlen, mut fmtlen) = (0i64, 0i64, 0i64);
    let mut lenpos: Option<usize> = None;
    let find = |hay: &[u8], needle: &[u8]| hay.windows(needle.len()).position(|w| w == needle);
    for (k, field) in line.split(|&c| c == b'\t').enumerate() {
        let id = k + 1;
        match id {
            1 => name = Some(String::from_utf8_lossy(field).into_owned()),
            2 => {
                let (v, used) = strtoll(field, true);
                if used == 0 {
                    return None;
                }
                beg = (v - 1).max(0);
                end = 1;
            }
            4 => {
                if !field.is_empty() {
                    end = beg + field.len() as i64;
                }
                reflen = field.len() as i64;
                svlenals.push(false);
            }
            5 => {
                for a in field.split(|&c| c == b',').take(65535) {
                    let sv = svlen_on_ref(a);
                    use_svlen |= sv;
                    getlen |= !sv && (a == b"<*>" || a == b"<NON_REF>");
                    svlenals.push(sv);
                }
            }
            8 => {
                let s = if field.starts_with(b"END=") {
                    Some(4)
                } else if find(field, b"END=").is_some() {
                    find(field, b";END=").map(|p| p + 5)
                } else {
                    None
                };
                if let Some(s) = s
                    && field.get(s) != Some(&b'.')
                {
                    let (e, _) = strtoll(&field[s..], true);
                    if e > beg {
                        end = e;
                    } else {
                        log::warn!("VCF INFO/END={e} is smaller than POS; this tag will be ignored");
                    }
                }
                let s = if field.starts_with(b"SVLEN=") {
                    Some(6)
                } else if find(field, b"SVLEN=").is_some() {
                    find(field, b";SVLEN=").map(|p| p + 7)
                } else {
                    None
                };
                if let Some(mut s) = s {
                    let mut d = 1;
                    while d < svlenals.len() {
                        let tmp = if use_svlen && svlenals[d] { strtoll(&field[s..], false).0.abs() } else { 1 };
                        svlen = svlen.max(tmp);
                        d += 1;
                        match field[s..].iter().position(|&c| c == b',') {
                            Some(p) => s += p + 1,
                            None => break,
                        }
                    }
                }
            }
            9 if getlen => {
                lenpos = field.split(|&c| c == b':').position(|f| f == b"LEN");
                if lenpos.is_none() {
                    break;
                }
            }
            _ if id > 9 && getlen => {
                if let Some(lp) = lenpos {
                    let tmp = field.split(|&c| c == b':').nth(lp).map_or(0, |f| strtoll(f, false).0);
                    fmtlen = fmtlen.max(tmp);
                }
            }
            _ => {}
        }
    }
    let name = name?;
    if beg < 0 {
        return None;
    }
    end = end.max(beg + reflen.max(svlen).max(fmtlen));
    Some((name, beg, end))
}

fn svlen_on_ref(alt: &[u8]) -> bool {
    record::svlen_on_ref_for_vcf_alt(alt)
}

// --- pysam record API -----------------------------------------------------------

/// A scalar FORMAT value as pysam returns it.
#[derive(Clone, Debug, PartialEq)]
pub enum FormatScalar {
    Int(i64),
    Float(f32),
}

/// Value `i` of a FORMAT/INFO vector as (value, missing, vector end).
fn value_at(typ: u8, p: &[u8], i: usize) -> Option<(i64, bool, bool)> {
    Some(match typ {
        BCF_BT_INT8 => {
            let v = p[i] as i8;
            (v as i64, v == INT8_MISSING, v == INT8_VECTOR_END)
        }
        BCF_BT_INT16 => {
            let v = i16::from_le_bytes([p[2 * i], p[2 * i + 1]]);
            (v as i64, v == INT16_MISSING, v == INT16_VECTOR_END)
        }
        BCF_BT_INT32 => {
            let v = i32::from_le_bytes(p[4 * i..4 * i + 4].try_into().unwrap());
            (v as i64, v == INT32_MISSING, v == INT32_VECTOR_END)
        }
        _ => return None,
    })
}

fn float_at(p: &[u8], i: usize) -> u32 {
    u32::from_le_bytes(p[4 * i..4 * i + 4].try_into().unwrap())
}

impl Record {
    pub fn rid(&self) -> i32 {
        self.rid
    }

    pub fn chrom(&self, header: &Header) -> String {
        header.core.id_key(BCF_DT_CTG, self.rid).unwrap_or("").to_string()
    }

    /// 0-based start
    pub fn pos(&self) -> i64 {
        self.pos
    }

    pub fn n_sample(&self) -> usize {
        self.n_sample as usize
    }

    pub fn alleles(&self) -> Result<Vec<String>> {
        Ok(self.alleles.iter().map(|a| String::from_utf8_lossy(a).into_owned()).collect())
    }

    /// pysam `record.alts` (None if fewer than two alleles)
    pub fn alts(&self) -> Result<Option<Vec<String>>> {
        let a = self.alleles()?;
        if a.len() < 2 {
            return Ok(None);
        }
        Ok(Some(a[1..].to_vec()))
    }

    /// pysam `record.ref`
    pub fn reference(&self) -> Result<Option<String>> {
        Ok(self.alleles()?.into_iter().next())
    }

    fn key(header: &Header, id: i32) -> String {
        header.core.id_key(BCF_DT_ID, id).unwrap_or("").to_string()
    }

    /// INFO keys present (pysam `record.info` iteration; skips END)
    pub fn info_keys(&self, header: &Header) -> Result<Vec<String>> {
        Ok(self.info.iter().filter(|i| !i.deleted).map(|i| Self::key(header, i.key)).filter(|k| k != "END").collect())
    }

    /// FORMAT keys present (pysam `record.format` iteration)
    pub fn format_keys(&self, header: &Header) -> Result<Vec<String>> {
        Ok(self.fmt.iter().filter(|f| !f.deleted).map(|f| Self::key(header, f.id)).collect())
    }

    fn fmt_of(&self, header: &Header, key: &str) -> Option<usize> {
        self.get_fmt(&header.core, key)
    }

    /// pysam `del record.format[key]` (`bcf_update_format` with no values)
    pub fn format_delete(&mut self, header: &Header, key: &str) -> Result<()> {
        self.update_format(&header.core, key, None)
    }

    /// pysam `record.qual = None`
    pub fn set_qual_missing(&mut self) {
        self.qual = FLOAT_MISSING;
    }

    /// pysam `call[key] = value` (`bcf_format_set_value`); deleting a value
    /// of one sample (`del call[key]`) assigns `PyFormatValue::Nulls` with
    /// the count from `format_null_count`.
    pub fn format_set_value(&mut self, header: &Header, sample: usize, key: &str, value: &PyFormatValue) -> Result<()> {
        // bcf_get_fmt also returns a deleted field (its data pointer is NULL)
        let fi = self.fmt_of(header, key);
        let deleted = fi.is_some_and(|i| self.fmt[i].deleted);
        let fmt_id = match fi {
            Some(i) => self.fmt[i].id,
            None => {
                if !header.is_defined(BCF_HL_FMT, key) {
                    bail!("KeyError: 'unknown format: {key}'");
                }
                header.id(key)
            }
        };
        let mut fmt_type = header.id2type(BCF_HL_FMT, fmt_id);
        if fmt_type == BCF_HT_FLAG {
            bail!("Flag types are not allowed on FORMATs");
        }
        let is_gt = fmt_id == header.id("GT");
        // the values as pysam sees them after bcf_check_values
        enum Values {
            Ints(Vec<Option<i64>>),
            Floats(Vec<Option<f64>>),
            Bytes(Vec<u8>),
        }
        let n_allele = self.n_allele() as u32;
        let encoded;
        let value = if is_gt {
            fmt_type = BCF_HT_INT;
            let alleles: Vec<Option<u32>> = match value {
                PyFormatValue::Alleles(a) => a.clone(),
                PyFormatValue::Nulls(n) => vec![None; *n],
                _ => bail!("TypeError: invalid GT value"),
            };
            let mut ints = Vec::with_capacity(alleles.len());
            for a in alleles {
                ints.push(Some(match a {
                    None => 0,
                    Some(a) if a < n_allele => i64::from((a + 1) << 1),
                    Some(_) => bail!("Invalid allele index"),
                }));
            }
            encoded = PyFormatValue::Ints(ints);
            &encoded
        } else {
            value
        };
        let (mut count, scalar) = self.value_count(header, fmt_id, sample)?;
        if is_gt {
            count = -1;
        }
        let given = match value {
            PyFormatValue::Ints(v) => v.len(),
            PyFormatValue::Floats(v) => v.len(),
            PyFormatValue::Nulls(n) => *n,
            PyFormatValue::Alleles(v) => v.len(),
        } as i64;
        if count != -1 && count != given {
            if scalar {
                bail!("TypeError: value expected to be scalar, given len={given}");
            }
            bail!("TypeError: values expected to be {count}-tuple, given len={given}");
        }
        let values = match (fmt_type, value) {
            (BCF_HT_INT, PyFormatValue::Ints(v)) => {
                if v.iter().flatten().any(|&x| x <= INT32_MISSING as i64 || x > i32::MAX as i64) {
                    bail!("Integer value too small/large to store in VCF/BCF");
                }
                Values::Ints(v.clone())
            }
            (BCF_HT_INT, PyFormatValue::Nulls(n)) => Values::Ints(vec![None; *n]),
            (BCF_HT_REAL, PyFormatValue::Floats(v)) => Values::Floats(v.clone()),
            (BCF_HT_REAL, PyFormatValue::Ints(v)) => Values::Floats(v.iter().map(|x| x.map(|x| x as f64)).collect()),
            (BCF_HT_REAL, PyFormatValue::Nulls(n)) => Values::Floats(vec![None; *n]),
            (BCF_HT_STR, PyFormatValue::Nulls(n)) => Values::Bytes(vec![b','; n.saturating_sub(1)]),
            _ => bail!("TypeError: unsupported value for FORMAT/{key}"),
        };
        let len = match &values {
            Values::Ints(v) => v.len(),
            Values::Floats(v) => v.len(),
            Values::Bytes(v) => v.len(),
        };
        let (bt_type, bt_len) = fi.map_or((u8::MAX, -1i64), |i| (self.fmt[i].typ, self.fmt[i].n as i64));
        let mut realloc = len as i64 > bt_len;
        if !realloc && let Values::Ints(v) = &values {
            realloc = v.iter().flatten().any(|&x| match bt_type {
                BCF_BT_INT8 => !((INT8_MISSING as i64) < x && x <= i8::MAX as i64),
                BCF_BT_INT16 => !((INT16_MISSING as i64) < x && x <= i16::MAX as i64),
                _ => false,
            });
        }
        let vlen = count < 0;
        // bcf_object_to_array
        let write = |values: &Values, data: &mut [u8], typ: u8, n: usize| -> Result<()> {
            match (values, typ) {
                (Values::Bytes(b), BCF_BT_CHAR) => {
                    data[..b.len()].copy_from_slice(b);
                    data[b.len()..n].fill(0);
                }
                (Values::Ints(v), BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32) => write_int_array(v, data, typ, n),
                (Values::Floats(v), BCF_BT_FLOAT) => {
                    for i in 0..n {
                        let bits = match v.get(i) {
                            Some(Some(x)) => (*x as f32).to_bits(),
                            Some(None) => FLOAT_MISSING,
                            None => FLOAT_VECTOR_END,
                        };
                        data[4 * i..4 * i + 4].copy_from_slice(&bits.to_le_bytes());
                    }
                }
                _ => bail!("TypeError: unsupported type"),
            }
            Ok(())
        };
        let nsamples = self.n_sample();
        if deleted && (!realloc || nsamples > 1) {
            // pysam would write to or copy from the NULL data of the field
            bail!("FORMAT/{key} was deleted from the record and cannot be set again (pysam crashes)");
        }
        if let Some(i) = fi
            && !realloc
        {
            let f = &mut self.fmt[i];
            let (typ, n) = (f.typ, f.n as usize);
            return write(&values, f.sample_mut(sample), typ, n);
        }
        let mut alloc_len = len.max(1);
        if let Some(i) = fi {
            alloc_len = alloc_len.max(self.fmt[i].n as usize);
        }
        let (dst_type, width) = match fmt_type {
            BCF_HT_INT => (BCF_BT_INT32, 4),
            BCF_HT_REAL => (BCF_BT_FLOAT, 4),
            BCF_HT_STR => (BCF_BT_CHAR, 1),
            _ => bail!("Unsupported FORMAT type"),
        };
        // bcf_empty_array
        let mut new_values = vec![0u8; nsamples * alloc_len * width];
        let empty: Vec<u8> = match fmt_type {
            BCF_HT_INT => (if vlen { INT32_VECTOR_END } else { INT32_MISSING }).to_le_bytes().to_vec(),
            BCF_HT_REAL => (if vlen { FLOAT_VECTOR_END } else { FLOAT_MISSING }).to_le_bytes().to_vec(),
            _ => vec![if vlen { 0 } else { 0x07 }],
        };
        for chunk in new_values.chunks_mut(width) {
            chunk.copy_from_slice(&empty);
        }
        let dst_size = alloc_len * width;
        if let Some(i) = fi
            && nsamples > 1
        {
            let f = &self.fmt[i];
            for s in 0..nsamples {
                let dst = &mut new_values[s * dst_size..(s + 1) * dst_size];
                let (src, n) = (f.sample(s), f.n as usize);
                match (f.typ, dst_type) {
                    (BCF_BT_CHAR, BCF_BT_CHAR) => {
                        dst[..n].copy_from_slice(&src[..n]);
                        dst[n..].fill(0);
                    }
                    (BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32, BCF_BT_INT32) => {
                        let mut ints = vec![0i32; alloc_len];
                        copy_expand_int(src, f.typ, n, &mut ints, vlen);
                        for (k, v) in ints.iter().enumerate() {
                            dst[4 * k..4 * k + 4].copy_from_slice(&v.to_le_bytes());
                        }
                    }
                    (BCF_BT_FLOAT, BCF_BT_FLOAT) => {
                        dst[..4 * n].copy_from_slice(&src[..4 * n]);
                        let pad = if vlen { FLOAT_VECTOR_END } else { FLOAT_MISSING };
                        for k in n..alloc_len {
                            dst[4 * k..4 * k + 4].copy_from_slice(&pad.to_le_bytes());
                        }
                    }
                    _ => bail!("TypeError: unsupported types"),
                }
            }
        }
        write(&values, &mut new_values[sample * dst_size..(sample + 1) * dst_size], dst_type, alloc_len)?;
        let result = match dst_type {
            BCF_BT_INT32 => {
                let ints: Vec<i32> = new_values.as_chunks::<4>().0.iter().map(|c| i32::from_le_bytes(*c)).collect();
                self.update_format(&header.core, key, Some(FmtValues::Int(&ints)))
            }
            BCF_BT_FLOAT => {
                let bits: Vec<u32> = new_values.as_chunks::<4>().0.iter().map(|c| u32::from_le_bytes(*c)).collect();
                self.update_format(&header.core, key, Some(FmtValues::Float(&bits)))
            }
            _ => self.update_format(&header.core, key, Some(FmtValues::Str(&new_values))),
        };
        result.map_err(|e| e.context("Unable to update format values"))
    }

    /// pysam `bcf_get_value_count` for a FORMAT field: (count, scalar);
    /// count -1 for a variable number
    fn value_count(&self, header: &Header, id: i32, sample: usize) -> Result<(i64, bool)> {
        let length = header.core.id2length(BCF_HL_FMT, id);
        let number = header.core.id2number(BCF_HL_FMT, id) as i64;
        Ok(if id == header.id("GT") {
            (number, false)
        } else if length == BCF_VL_FIXED {
            (number, number == 1)
        } else if length == BCF_VL_R {
            (self.n_allele() as i64, false)
        } else if length == BCF_VL_A {
            (self.n_allele() as i64 - 1, false)
        } else if length == BCF_VL_G {
            (self.genotype_count(header, sample), false)
        } else if length == BCF_VL_VAR {
            (-1, false)
        } else {
            bail!("Unknown format length");
        })
    }

    /// pysam `bcf_genotype_count`: genotypes for the ploidy of the sample's GT
    fn genotype_count(&self, header: &Header, sample: usize) -> i64 {
        let gt = header.id("GT");
        let Some(f) = self.fmt.iter().find(|f| f.id == gt && !f.deleted) else { return 0 };
        if f.n == 0 || self.n_sample() == 0 {
            return 0;
        }
        let base = f.sample(sample);
        let ploidy = (0..f.n as usize)
            .take_while(|&i| value_at(f.typ, base, i).is_some_and(|(_, _, vector_end)| !vector_end))
            .count() as i64;
        // comb(alleles + ploidy - 1, ploidy)
        let n = self.n_allele() as i64 + ploidy - 1;
        (0..ploidy).fold(1i64, |acc, i| acc * (n - i) / (i + 1))
    }

    /// The value pysam assigns when deleting FORMAT `key` of one sample
    /// (`bcf_format_del_value`)
    pub fn format_null_value(&self, header: &Header, sample: usize, key: &str) -> Result<PyFormatValue> {
        let id = match self.fmt_of(header, key).filter(|&i| !self.fmt[i].deleted) {
            Some(i) => self.fmt[i].id,
            None => bail!("KeyError: '{key}'"),
        };
        let (count, scalar) = self.value_count(header, id, sample)?;
        Ok(PyFormatValue::Nulls(if count <= 0 {
            0
        } else if scalar {
            1
        } else {
            count as usize
        }))
    }

    /// pysam `call.keys()`: FORMAT keys with data
    pub fn format_keys_present(&self, header: &Header) -> Vec<String> {
        self.format_keys(header).unwrap_or_default()
    }

    /// pysam `key in record.samples[i]`
    pub fn has_format(&self, header: &Header, key: &str) -> Result<bool> {
        Ok(self.fmt_of(header, key).is_some_and(|i| !self.fmt[i].deleted))
    }

    /// Whether the sample has a non-missing first value for FORMAT `key`.
    pub fn format_value_present(&self, header: &Header, sample: usize, key: &str) -> Result<bool> {
        let Some(i) = self.fmt_of(header, key) else { return Ok(false) };
        let f = &self.fmt[i];
        if f.deleted || f.n == 0 {
            return Ok(false);
        }
        let base = f.sample(sample);
        Ok(match f.typ {
            BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 => {
                let (_, missing, vector_end) = value_at(f.typ, base, 0).expect("integer type");
                !missing && !vector_end
            }
            BCF_BT_FLOAT => {
                let b = float_at(base, 0);
                b != FLOAT_MISSING && b != FLOAT_VECTOR_END
            }
            BCF_BT_CHAR => base[0] != 0 && !(base[0] == b'.' && (f.size == 1 || base[1] == 0)),
            _ => true,
        })
    }

    /// The first FORMAT field if it is a non-empty GT.
    fn first_fmt_gt(&self, header: &Header) -> Option<usize> {
        let f = self.fmt.first()?;
        (Self::key(header, f.id) == "GT" && f.n != 0).then_some(0)
    }

    /// pysam `record.samples[i].get("GT")`: None if there is no GT field,
    /// otherwise the allele indices (None for missing). Like pysam, this
    /// reads the *first* FORMAT field and returns () if that is not GT.
    pub fn gt(&self, header: &Header, sample: usize) -> Result<Option<Vec<Option<u32>>>> {
        if !self.has_format(header, "GT")? {
            return Ok(None);
        }
        let n_allele = self.n_allele() as i64;
        let Some(fi) = self.first_fmt_gt(header) else { return Ok(Some(Vec::new())) };
        if sample >= self.n_sample() {
            return Ok(Some(Vec::new()));
        }
        let f = &self.fmt[fi];
        let base = f.sample(sample);
        let mut alleles = Vec::new();
        for i in 0..f.n as usize {
            let Some((val, _, vector_end)) = value_at(f.typ, base, i) else { bail!("unsupported GT type") };
            if vector_end {
                break;
            }
            let a = if val == 0 { -1 } else { (val >> 1) - 1 };
            alleles.push(if 0 <= a && a < n_allele { Some(a as u32) } else { None });
        }
        Ok(Some(alleles))
    }

    /// pysam `record.samples[i].phased`
    pub fn sample_phased(&self, header: &Header, sample: usize) -> Result<bool> {
        if sample >= self.n_sample() || self.fmt.is_empty() {
            return Ok(false);
        }
        let Some(fi) = self.first_fmt_gt(header) else { return Ok(false) };
        let f = &self.fmt[fi];
        let base = f.sample(sample);
        let mut phased = false;
        for i in 0..f.n as usize {
            let Some((v, missing, vector_end)) = value_at(f.typ, base, i) else { return Ok(false) };
            if vector_end {
                break;
            } else if missing {
                continue;
            } else if i > 0 && v & 1 == 0 {
                return Ok(false);
            } else {
                phased = true;
            }
        }
        Ok(phased)
    }

    /// pysam `record.samples[i].phased = value` (`bcf_sample_set_phased`)
    pub fn set_phased(&mut self, header: &Header, sample: usize, phased: bool) -> Result<()> {
        if sample >= self.n_sample() || self.fmt.is_empty() {
            return Ok(());
        }
        let Some(fi) = self.first_fmt_gt(header) else { bail!("Cannot set phased before genotype is set") };
        let f = &mut self.fmt[fi];
        let (typ, n) = (f.typ, f.n as usize);
        let base = f.sample_mut(sample);
        let ph = phased as u8;
        for i in 0..n {
            let Some((_, missing, vector_end)) = value_at(typ, base, i) else { continue };
            if vector_end {
                break;
            } else if missing {
                continue;
            } else if i > 0 {
                // the phase bit is in the lowest byte (little endian)
                let w = 1 << (typ - 1);
                base[i * w] = (base[i * w] & 0xFE) | ph;
            }
        }
        Ok(())
    }

    /// pysam `record.samples[i]["GT"] = alleles` (alleles are allele indices)
    pub fn set_gt(&mut self, header: &Header, sample: usize, alleles: &[u32]) -> Result<()> {
        let n_allele = self.n_allele() as u32;
        let mut values = Vec::new();
        for &a in alleles {
            if a >= n_allele {
                bail!("Invalid allele index");
            }
            values.push(Some(((a as i32 + 1) << 1) as i64));
        }
        self.format_set_int(header, sample, "GT", &values)
    }

    /// pysam `record.samples[i][key] = value` for integer-valued FORMAT
    /// fields (including GT, already encoded). `None` entries are missing
    /// values; for a scalar field pass a single-element slice.
    pub fn format_set_int(&mut self, header: &Header, sample: usize, key: &str, values: &[Option<i64>]) -> Result<()> {
        let fi = self.fmt_of(header, key);
        let fmt_id = fi.map_or_else(|| header.id(key), |i| self.fmt[i].id);
        if fmt_id < 0 || !header.is_defined(BCF_HL_FMT, key) {
            bail!("unknown format: {key}");
        }
        let is_gt = key == "GT";
        let fmt_type = if is_gt { BCF_HT_INT } else { header.id2type(BCF_HL_FMT, fmt_id) };
        if fmt_type != BCF_HT_INT && fmt_type != BCF_HT_REAL {
            bail!("unsupported FORMAT type for {key}");
        }
        if values.iter().flatten().any(|&v| v <= INT32_MISSING as i64 || v > i32::MAX as i64) {
            bail!("Integer value too small/large to store in VCF/BCF");
        }
        let vlen = is_gt || (header.core.id2length(BCF_HL_FMT, fmt_id) != BCF_VL_FIXED);
        if fmt_type == BCF_HT_REAL {
            // Float FORMAT with integer values (PS declared as Float): pysam
            // writes floats.
            return self.format_set_float_from_int(header, sample, key, fi, values);
        }
        let (bt_type, bt_len) = fi.map_or((u8::MAX, -1), |i| (self.fmt[i].typ, self.fmt[i].n as i64));
        let mut realloc = values.len() as i64 > bt_len;
        if !realloc {
            for v in values.iter().flatten() {
                let fits = match bt_type {
                    BCF_BT_INT8 => (INT8_MISSING as i64) < *v && *v <= i8::MAX as i64,
                    BCF_BT_INT16 => (INT16_MISSING as i64) < *v && *v <= i16::MAX as i64,
                    _ => true,
                };
                if !fits {
                    realloc = true;
                    break;
                }
            }
        }
        if let Some(i) = fi
            && !realloc
        {
            let f = &mut self.fmt[i];
            let (typ, n) = (f.typ, f.n as usize);
            write_int_array(values, f.sample_mut(sample), typ, n);
            return Ok(());
        }
        let mut alloc_len = values.len().max(1);
        if let Some(i) = fi {
            alloc_len = alloc_len.max(self.fmt[i].n as usize);
        }
        let nsamples = self.n_sample();
        let fill = if vlen { INT32_VECTOR_END } else { INT32_MISSING };
        let mut new_values = vec![fill; nsamples * alloc_len];
        if let Some(i) = fi
            && nsamples > 1
        {
            let f = &self.fmt[i];
            for s in 0..nsamples {
                copy_expand_int(
                    f.sample(s),
                    f.typ,
                    f.n as usize,
                    &mut new_values[s * alloc_len..(s + 1) * alloc_len],
                    vlen,
                );
            }
        }
        let mut sample_bytes = vec![0u8; alloc_len * 4];
        write_int_array(values, &mut sample_bytes, BCF_BT_INT32, alloc_len);
        for (k, c) in sample_bytes.as_chunks::<4>().0.iter().enumerate() {
            new_values[sample * alloc_len + k] = i32::from_le_bytes(*c);
        }
        self.update_format(&header.core, key, Some(FmtValues::Int(&new_values)))
            .map_err(|e| e.context("Unable to update format values"))
    }

    fn format_set_float_from_int(
        &mut self,
        header: &Header,
        sample: usize,
        key: &str,
        fi: Option<usize>,
        values: &[Option<i64>],
    ) -> Result<()> {
        let nsamples = self.n_sample();
        let mut alloc_len = values.len().max(1);
        if let Some(i) = fi {
            alloc_len = alloc_len.max(self.fmt[i].n as usize);
        }
        let mut new_values = vec![FLOAT_MISSING; nsamples * alloc_len];
        if let Some(i) = fi {
            let f = &self.fmt[i];
            if f.typ == BCF_BT_FLOAT {
                for s in 0..nsamples {
                    for j in 0..f.n as usize {
                        new_values[s * alloc_len + j] = float_at(f.sample(s), j);
                    }
                }
            }
        }
        for j in 0..alloc_len {
            new_values[sample * alloc_len + j] = match values.get(j) {
                Some(Some(v)) => (*v as f32).to_bits(),
                Some(None) => FLOAT_MISSING,
                None => FLOAT_VECTOR_END,
            };
        }
        self.update_format(&header.core, key, Some(FmtValues::Float(&new_values)))
            .map_err(|e| e.context("Unable to update format values"))
    }

    /// pysam's `bcf_sync_end`, called before every write.
    pub fn sync_end(&mut self, header: &Header) -> Result<()> {
        let ref_len = self.alleles.first().map_or(0, |a| a.len()) as i64;
        let symbolic = self.alleles.iter().skip(1).any(|a| a.starts_with(b"<") && a.ends_with(b">"));
        if !symbolic && (self.alleles.is_empty() || self.rlen == ref_len) {
            if header.is_defined(BCF_HL_INFO, "END")
                && let Some(info) = self.get_info(&header.core, "END")
                && !info.deleted
            {
                self.update_info(&header.core, "END", None).map_err(|e| e.context("Unable to delete END"))?;
            }
        } else {
            if !header.is_defined(BCF_HL_INFO, "END") {
                bail!("END not defined in header");
            }
            // bcf_info_set_value(record, "END", pos + rlen)
            let id = header.id("END");
            let end = self.pos + self.rlen;
            let n_allele = self.alleles.len() as i64;
            let count = match header.core.id2length(BCF_HL_INFO, id) {
                BCF_VL_FIXED => header.core.id2number(BCF_HL_INFO, id) as i64,
                BCF_VL_VAR => -1,
                2 => n_allele - 1,
                3 => n_allele * (n_allele + 1) / 2,
                4 => n_allele,
                _ => bail!("Unknown format length"),
            };
            let typ = header.id2type(BCF_HL_INFO, id);
            if typ != BCF_HT_FLAG && count != -1 && count != 1 {
                bail!("values expected to be {count}-tuple, given len=1");
            }
            match typ {
                BCF_HT_INT => {
                    let v = i32::try_from(end)
                        .ok()
                        .filter(|&v| v != INT32_MISSING)
                        .ok_or_else(|| anyhow!("Integer value too small/large to store in VCF/BCF"))?;
                    // the old values are padded with missing values
                    if self.get_info(&header.core, "END").is_some_and(|i| i.len > 1) {
                        bail!("Unable to update INFO values");
                    }
                    self.update_info(&header.core, "END", Some(InfoValues::Int(&[v])))
                        .map_err(|e| e.context("Unable to update INFO values"))?;
                }
                BCF_HT_FLAG if end == 0 => self.update_info(&header.core, "END", None)?,
                BCF_HT_FLAG if end == 1 => bail!("Unable to update INFO values"),
                BCF_HT_FLAG => bail!("Flag values must be: True, False, None, 1, 0"),
                BCF_HT_STR => bail!("Argument must be string, bytes or unicode."),
                // htslib refuses non-integer END values
                _ => bail!("Unable to update INFO values"),
            }
        }
        Ok(())
    }

    /// Sets rlen as htslib 1.23.1 does after parsing a VCF line (`get_rlen`).
    pub fn set_rlen_from_vcf(&mut self, header: &Header) -> Result<()> {
        self.rlen = self.get_rlen(&header.core, self.alleles.first().map_or(0, |a| a.len() as i64));
        Ok(())
    }

    /// pysam `record.samples[i].get(key)` for a scalar numeric field:
    /// Ok(None) if the field is absent (KeyError), Ok(Some(None)) if missing.
    pub fn format_scalar(&self, header: &Header, sample: usize, key: &str) -> Result<Option<Option<FormatScalar>>> {
        let Some(i) = self.fmt_of(header, key).filter(|&i| !self.fmt[i].deleted) else { return Ok(None) };
        let f = &self.fmt[i];
        if f.n == 0 || f.size == 0 {
            return Ok(Some(None));
        }
        let base = f.sample(sample);
        Ok(Some(match f.typ {
            BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 => {
                let (v, missing, vector_end) = value_at(f.typ, base, 0).expect("integer type");
                (!missing && !vector_end).then_some(FormatScalar::Int(v))
            }
            BCF_BT_FLOAT => {
                let b = float_at(base, 0);
                (b != FLOAT_MISSING && b != FLOAT_VECTOR_END).then_some(FormatScalar::Float(f32::from_bits(b)))
            }
            _ => bail!("unsupported type for FORMAT {key}"),
        }))
    }

    /// pysam `record.samples[i].get(key)` for a string field: None if the
    /// field is absent, otherwise the comma-separated values (None for
    /// empty ones).
    pub fn format_strings(&self, header: &Header, sample: usize, key: &str) -> Result<Option<Vec<Option<String>>>> {
        let Some(i) = self.fmt_of(header, key).filter(|&i| !self.fmt[i].deleted) else { return Ok(None) };
        let f = &self.fmt[i];
        if f.typ != BCF_BT_CHAR {
            bail!("FORMAT {key} is not a string");
        }
        let n = f.size;
        if n == 0 {
            return Ok(Some(Vec::new()));
        }
        let data = f.sample(sample);
        let bytes: &[u8] = if data[n - 1] == 0 {
            let end = data.iter().position(|&b| b == 0).unwrap();
            &data[..end]
        } else {
            data
        };
        Ok(Some(
            bytes
                .split(|&b| b == b',')
                .map(|v| if v.is_empty() || v == [7u8] { None } else { Some(String::from_utf8_lossy(v).into_owned()) })
                .collect(),
        ))
    }

    /// pysam `record.samples[i].get(key)` for a numeric array field: None if
    /// the field is absent, otherwise the values up to the vector end (None
    /// for missing values). Floats are returned as f32 widened to f64, like
    /// pysam does.
    pub fn format_numbers(&self, header: &Header, sample: usize, key: &str) -> Result<Option<Vec<Option<f64>>>> {
        let Some(i) = self.fmt_of(header, key).filter(|&i| !self.fmt[i].deleted) else { return Ok(None) };
        let f = &self.fmt[i];
        let base = f.sample(sample);
        let mut out = Vec::new();
        for i in 0..f.n as usize {
            let v = match f.typ {
                BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 => {
                    let (v, missing, vector_end) = value_at(f.typ, base, i).expect("integer type");
                    if vector_end {
                        break;
                    }
                    (!missing).then_some(v as f64)
                }
                BCF_BT_FLOAT => {
                    let b = float_at(base, i);
                    if b == FLOAT_VECTOR_END {
                        break;
                    }
                    (b != FLOAT_MISSING).then_some(f32::from_bits(b) as f64)
                }
                _ => bail!("unsupported type for FORMAT {key}"),
            };
            out.push(v);
        }
        Ok(Some(out))
    }

    /// pysam `record.samples[i][key] = value` for a String FORMAT field
    /// (`None` sets it to missing).
    pub fn format_set_str(&mut self, header: &Header, sample: usize, key: &str, value: Option<&str>) -> Result<()> {
        let fi = self.fmt_of(header, key);
        if fi.is_none() && !header.is_defined(BCF_HL_FMT, key) {
            bail!("unknown format: {key}");
        }
        let values: &[u8] = value.map(|v| v.as_bytes()).unwrap_or(b"");
        let (bt_len, bt_n) = fi.map_or((-1i64, 0usize), |i| (self.fmt[i].n as i64, self.fmt[i].n as usize));
        let realloc = values.len() as i64 > bt_len;
        if let Some(i) = fi
            && !realloc
        {
            let dst = self.fmt[i].sample_mut(sample);
            dst[..values.len()].copy_from_slice(values);
            dst[values.len()..bt_n].fill(0);
            return Ok(());
        }
        let alloc_len = values.len().max(1).max(bt_n);
        let nsamples = self.n_sample();
        // bcf_empty_array for a vlen string: all vector_end (0)
        let mut new_values = vec![0u8; nsamples * alloc_len];
        if let Some(i) = fi
            && nsamples > 1
        {
            let f = &self.fmt[i];
            for s in 0..nsamples {
                let src = &f.sample(s)[..f.n as usize];
                new_values[s * alloc_len..s * alloc_len + src.len()].copy_from_slice(src);
            }
        }
        let dst = &mut new_values[sample * alloc_len..(sample + 1) * alloc_len];
        dst[..values.len()].copy_from_slice(values);
        dst[values.len()..].fill(0);
        self.update_format(&header.core, key, Some(FmtValues::Str(&new_values)))
            .map_err(|e| e.context("Unable to update format values"))
    }
}

/// pysam `bcf_object_to_array` for integer types.
fn write_int_array(values: &[Option<i64>], data: &mut [u8], bt_type: u8, n: usize) {
    for i in 0..n {
        let v = values.get(i);
        match bt_type {
            BCF_BT_INT8 => {
                data[i] = match v {
                    Some(Some(v)) => *v as i8,
                    Some(None) => INT8_MISSING,
                    None => INT8_VECTOR_END,
                } as u8
            }
            BCF_BT_INT16 => {
                let x = match v {
                    Some(Some(v)) => *v as i16,
                    Some(None) => INT16_MISSING,
                    None => INT16_VECTOR_END,
                };
                data[2 * i..2 * i + 2].copy_from_slice(&x.to_le_bytes());
            }
            BCF_BT_INT32 => {
                let x = match v {
                    Some(Some(v)) => *v as i32,
                    Some(None) => INT32_MISSING,
                    None => INT32_VECTOR_END,
                };
                data[4 * i..4 * i + 4].copy_from_slice(&x.to_le_bytes());
            }
            _ => panic!("unsupported type"),
        }
    }
}

/// pysam `bcf_copy_expand_array` (integer source types to int32).
fn copy_expand_int(src: &[u8], src_type: u8, src_values: usize, dst: &mut [i32], vlen: bool) {
    for (i, d) in dst.iter_mut().enumerate().take(src_values) {
        let (v, missing, vector_end) = value_at(src_type, src, i).expect("integer type");
        *d = if missing {
            INT32_MISSING
        } else if vector_end {
            INT32_VECTOR_END
        } else {
            v as i32
        };
    }
    for d in dst.iter_mut().skip(src_values) {
        *d = if vlen { INT32_VECTOR_END } else { INT32_MISSING };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tabix_end() {
        assert_eq!(tbx_parse_vcf(b"1\t100\t.\tACG\tA\t.\t.\tEND=500"), Some(("1".into(), 99, 500)));
        assert_eq!(tbx_parse_vcf(b"1\t100\t.\tA\t<DEL>\t.\t.\tSVLEN=-300"), Some(("1".into(), 99, 399)));
        assert_eq!(tbx_parse_vcf(b"1\t100\t.\tA\t<*>\t.\t.\t.\tGT:LEN\t0/0:50"), Some(("1".into(), 99, 149)));
        assert_eq!(tbx_parse_vcf(b"1\t100\t.\tA\tC\t.\t.\tXEND=900"), Some(("1".into(), 99, 100)));
    }
}

//! BAM/CRAM input and output on pasta.

pub mod bam;
pub mod cram;
pub mod header;
pub mod record;

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use pasta::sam::alignment::io::Write as _;
use pasta::{fasta, sam};

use self::bam::{BamReader, BamWriter};
use self::cram::CramReader;
pub use self::cram::Decoding;
pub use self::header::Header;
pub use self::record::{Aux, Record};

/// What to iterate over after `AlignmentFile::fetch`.
#[derive(Clone, Copy, Debug)]
pub enum Fetch {
    /// All records from the start (`HTS_IDX_START`)
    All,
    /// Unplaced unmapped records (`HTS_IDX_NOCOOR`)
    Unmapped,
    /// Records of a reference overlapping [start, end)
    Region(usize, i64, i64),
}

enum Inner {
    Bam(Box<BamReader>),
    Cram(Box<CramReader>),
}

/// An indexed BAM or CRAM file (format detected from its contents).
pub struct AlignmentFile {
    inner: Inner,
    path: PathBuf,
}

impl AlignmentFile {
    /// Opens `path`; `reference` is needed to decode CRAM files. `threads`
    /// decompress BAM files, and decode CRAM regions read for alleles
    /// ([`Decoding::Alleles`]).
    pub fn open(path: &str, reference: Option<&Path>, threads: usize) -> Result<AlignmentFile> {
        let p = Path::new(path);
        let mut magic = [0u8; 4];
        let n = File::open(p).with_context(|| format!("opening {path}"))?.read(&mut magic)?;
        let inner = if n >= 4 && &magic == b"CRAM" {
            Inner::Cram(Box::new(CramReader::open(p, reference, threads)?))
        } else if n >= 2 && magic[..2] == [0x1f, 0x8b] {
            Inner::Bam(Box::new(BamReader::open(p, threads)?))
        } else {
            bail!("{path} is not a BAM or CRAM file");
        };
        Ok(AlignmentFile { inner, path: p.to_path_buf() })
    }

    pub fn header(&self) -> &Header {
        match &self.inner {
            Inner::Bam(r) => &r.header,
            Inner::Cram(r) => &r.header,
        }
    }

    pub fn is_cram(&self) -> bool {
        matches!(self.inner, Inner::Cram(_))
    }

    /// What decoding gives CRAM records; BAM records are read as stored.
    pub fn set_decoding(&mut self, decoding: Decoding) {
        if let Inner::Cram(r) = &mut self.inner {
            r.set_decoding(decoding);
        }
    }

    /// See [`cram::CramReader::share_reference`]; BAM files need no reference.
    pub fn share_reference(&mut self, repository: &pasta::fasta::Repository) {
        if let Inner::Cram(r) = &mut self.inner {
            r.share_reference(repository.clone());
        }
    }

    /// Loads the index next to the file (Ok(false) if there is none).
    pub fn load_index(&mut self) -> Result<bool> {
        match &mut self.inner {
            Inner::Bam(r) => r.load_index(&self.path),
            Inner::Cram(r) => r.load_index(&self.path),
        }
    }

    pub fn fetch(&mut self, what: Fetch) -> Result<()> {
        match (&mut self.inner, what) {
            (Inner::Bam(r), Fetch::All) => r.fetch_all(),
            (Inner::Bam(r), Fetch::Unmapped) => r.fetch_unmapped(),
            (Inner::Bam(r), Fetch::Region(tid, beg, end)) => r.fetch(tid, beg, end),
            (Inner::Cram(r), Fetch::All) => r.fetch_all(),
            (Inner::Cram(r), Fetch::Unmapped) => r.fetch_unmapped(),
            (Inner::Cram(r), Fetch::Region(tid, beg, end)) => r.fetch(tid, beg, end),
        }
    }

    /// Whether a fetch of all of `tid` gives a record. From an indexed CRAM the index says,
    /// without decoding; otherwise this fetches and reads one record (ending the current fetch).
    pub fn has_records(&mut self, tid: usize) -> Result<bool> {
        if let Inner::Cram(r) = &self.inner
            && let Some(has) = r.has_records(tid)
        {
            return Ok(has);
        }
        self.fetch(Fetch::Region(tid, 0, i64::MAX >> 1))?;
        Ok(self.read()?.is_some())
    }

    /// The next record of the current fetch.
    pub fn read(&mut self) -> Result<Option<Record>> {
        match &mut self.inner {
            Inner::Bam(r) => r.read(),
            Inner::Cram(r) => r.read(),
        }
    }
}

/// A BAM or CRAM output file.
pub enum AlignmentWriter {
    Bam(BamWriter),
    Cram(Box<CramWriter>),
}

impl AlignmentWriter {
    pub fn write(&mut self, rec: &Record) -> Result<()> {
        match self {
            AlignmentWriter::Bam(w) => w.write(rec),
            AlignmentWriter::Cram(w) => w.write(rec),
        }
    }

    pub fn finish(self) -> Result<()> {
        match self {
            AlignmentWriter::Bam(w) => w.finish(),
            AlignmentWriter::Cram(w) => w.finish(),
        }
    }
}

/// A CRAM writer: records are converted from BAM layout to pasta records.
pub struct CramWriter {
    writer: pasta::cram::io::Writer<File>,
    header: sam::Header,
    buf: Vec<u8>,
}

impl CramWriter {
    pub fn create(path: &Path, header: &Header, reference: &Path) -> Result<CramWriter> {
        let fasta_reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(reference)
            .with_context(|| format!("opening reference {}", reference.display()))?;
        let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta_reader));
        let sam_header: sam::Header = header.text.parse().context("parsing the output header")?;
        let mut writer = pasta::cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(path)
            .with_context(|| format!("creating {}", path.display()))?;
        writer.write_header(&sam_header)?;
        Ok(CramWriter { writer, header: sam_header, buf: Vec::new() })
    }

    pub fn write(&mut self, rec: &Record) -> Result<()> {
        self.buf.clear();
        rec.write_bam(&mut self.buf)?;
        let mut reader = pasta::bam::io::Reader::from(&self.buf[..]);
        let mut record = pasta::bam::Record::default();
        reader.read_record(&mut record)?;
        let buf = sam::alignment::RecordBuf::try_from_alignment_record(&self.header, &record)?;
        self.writer.write_alignment_record(&self.header, &buf)?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.writer.try_finish(&self.header)?;
        self.writer.get_mut().flush()?;
        Ok(())
    }
}

//! The alignment and reference reading the original's C library does with htslib, on pasta
//! (through `snv-hts`).
//!
//! Reads are filtered exactly like `read_bam` in `medaka_bamiter.c`, both for [`Reads`] (which
//! the pileup counts are built from) and for the `bam_plp` engine ([`Pileup`], which the counts
//! are tested against), so they match the original's.

use std::path::Path;

use snv_hts::{Alignments, Overlaps, ReadSource as _, Records};

use crate::Error;

/// `BAM_FREVERSE`.
const FLAG_REVERSE: u16 = 0x10;

/// Reads `medaka_bamiter.c::read_bam` drops before they reach the pileup: unmapped, secondary,
/// supplementary, QC-failed and duplicate.
const EXCLUDED_FLAGS: u16 = 0x4 | 0x100 | 0x800 | 0x200 | 0x400;

/// `bam_mplp_set_maxcnt(mplp, 1 << 20)` in the original's C pileup (`calculate_*_pileup`).
const MAX_PILEUP_DEPTH: usize = 1 << 20;

impl From<snv_hts::Error> for Error {
    fn from(e: snv_hts::Error) -> Self {
        match e {
            snv_hts::Error::Open { path, .. } | snv_hts::Error::Fasta { path, .. } => Error::Open { path },
            snv_hts::Error::Index { path, .. } => Error::Index { path },
            snv_hts::Error::Region { region } => Error::Region { region },
            snv_hts::Error::Contig { contig } => Error::Contig { contig },
            snv_hts::Error::Pileup(_) => Error::Pileup,
            snv_hts::Error::Decode(_) | snv_hts::Error::Parse { .. } => Error::Decode,
        }
    }
}

/// `read_bam`'s filter: whether a read reaches the pileup.
fn kept(read: &snv_hts::Read, min_mapq: u8) -> bool {
    read.flags() & EXCLUDED_FLAGS == 0 && read.mapping_quality() >= min_mapq
}

/// An indexed BAM/CRAM file with its header (`bam_fset` in the original).
#[derive(Debug)]
pub struct AlignmentFile(Alignments);

impl AlignmentFile {
    /// Opens `path` and its index, using `reference` to decode CRAM.
    ///
    /// # Errors
    ///
    /// Fails if the file, its index or its header cannot be read.
    pub fn open(path: impl AsRef<Path>, reference: impl AsRef<Path>) -> Result<Self, Error> {
        Ok(Self(Alignments::open(path, reference)?))
    }

    /// Order alignments that start at the same coordinate before building the pileup.
    #[must_use]
    pub fn with_canonical_read_order(self, enabled: bool) -> Self {
        Self(self.0.with_canonical_read_order(enabled))
    }

    /// Target id of `contig` in the header, if present.
    pub fn tid(&self, contig: &str) -> Option<i32> {
        self.0.tid(contig)
    }

    /// Mapped reads on `contig` as the index counts them (`samtools idxstats`); `None` if the
    /// contig is not in the header or the index keeps no counts (as CRAM indexes do not).
    pub fn mapped_reads(&self, contig: &str) -> Option<u64> {
        self.0.mapped_reads(contig)
    }

    /// Starts a pileup over the reads `sam_itr_querys` returns for `region` (a samtools-style
    /// region string, passed through unchanged as the original does).
    ///
    /// # Errors
    ///
    /// Fails if the region cannot be parsed or queried.
    pub fn pileup(&mut self, region: &str, min_mapq: u8) -> Result<Pileup<'_>, Error> {
        let records = self.0.query(region)?;
        let filter: Box<dyn FnMut(&mut snv_hts::Read) -> bool> = Box::new(move |r| kept(r, min_mapq));
        Ok(Pileup(snv_hts::Pileup::new(records, MAX_PILEUP_DEPTH, Overlaps::Ignore, filter)))
    }

    /// The reads `sam_itr_querys` returns for `region` that the pileup engine would see (the
    /// same filter as [`pileup`](Self::pileup)), in coordinate order.
    ///
    /// # Errors
    ///
    /// Fails if the region cannot be parsed or queried.
    pub fn reads(&mut self, region: &str, min_mapq: u8) -> Result<Reads<'_>, Error> {
        Ok(Reads { records: self.0.query(region)?, min_mapq, current: None })
    }
}

/// The filtered reads of a region, one at a time.
#[derive(Debug)]
pub struct Reads<'file> {
    records: Records<'file>,
    min_mapq: u8,
    current: Option<snv_hts::Read>,
}

impl Reads<'_> {
    /// The next read, valid until the next call.
    ///
    /// # Errors
    ///
    /// Fails if a record cannot be decoded.
    pub fn next_read(&mut self) -> Result<Option<Read<'_>>, Error> {
        loop {
            match self.records.next_read()? {
                None => return Ok(None),
                Some(read) if kept(&read, self.min_mapq) => {
                    return Ok(Some(Read(self.current.insert(read))));
                }
                Some(_) => {}
            }
        }
    }
}

/// One alignment record.
#[derive(Clone, Copy)]
pub struct Read<'r>(&'r snv_hts::Read);

impl std::fmt::Debug for Read<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Read").field("tid", &self.tid()).field("pos", &self.pos()).finish_non_exhaustive()
    }
}

impl Read<'_> {
    /// Target id of the contig.
    pub fn tid(&self) -> i32 {
        self.0.reference_sequence_id()
    }

    /// 0-based leftmost reference position.
    pub fn pos(&self) -> i64 {
        self.0.position()
    }

    /// The read maps to the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.0.flags() & FLAG_REVERSE != 0
    }

    /// The CIGAR operations, each `length << 4 | op`.
    pub fn cigar(&self) -> &[u32] {
        self.0.cigar()
    }

    /// Number of query bases.
    pub fn query_len(&self) -> i32 {
        i32::try_from(self.0.sequence_len()).expect("read length fits i32")
    }

    /// The query bases as `seq_nt16` codes, two per byte, high nibble first (`bam_get_seq`).
    pub fn packed_seq(&self) -> &[u8] {
        let raw = self.0.as_bytes();
        let name_len = usize::from(raw[8]);
        let n_cigar = usize::from(u16::from_le_bytes([raw[12], raw[13]]));
        let start = 32 + name_len + 4 * n_cigar;
        &raw[start..start + self.0.sequence_len().div_ceil(2)]
    }

    /// 4-bit `seq_nt16` code of query base `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is outside the read.
    pub fn base_code(&self, i: i32) -> u8 {
        let n = self.0.sequence_len();
        let at = usize::try_from(i).ok().filter(|&i| i < n);
        let Some(i) = at else { panic!("query position {i} outside read of {n} bases") };
        self.0.base(i)
    }
}

type Filter = Box<dyn FnMut(&mut snv_hts::Read) -> bool>;

/// A running pileup over one region.
#[derive(Debug)]
pub struct Pileup<'file>(snv_hts::Pileup<Records<'file>, Filter>);

impl Pileup<'_> {
    /// Advances to the next pileup column.
    ///
    /// # Errors
    ///
    /// Fails if a record cannot be decoded (the original exits in that case).
    pub fn next_column(&mut self) -> Result<Option<Column<'_>>, Error> {
        Ok(self.0.next_column()?.map(|c| Column {
            tid: c.tid,
            pos: i32::try_from(c.pos).expect("position fits i32"),
            column: c,
        }))
    }
}

/// One pileup column: every filtered read overlapping reference position `pos`.
pub struct Column<'plp> {
    /// Target id of the contig.
    pub tid: i32,
    /// 0-based reference position.
    pub pos: i32,
    column: snv_hts::Column<'plp>,
}

impl std::fmt::Debug for Column<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Column")
            .field("tid", &self.tid)
            .field("pos", &self.pos)
            .field("depth", &self.column.entries().len())
            .finish()
    }
}

impl<'plp> Column<'plp> {
    /// The reads in this column, in htslib's order.
    pub fn entries(&self) -> impl Iterator<Item = Entry<'plp>> + 'plp {
        let column = self.column;
        column.entries().iter().map(move |e| Entry { entry: e, read: column.read(e) })
    }
}

/// One read in a pileup column.
#[derive(Clone, Copy)]
pub struct Entry<'plp> {
    entry: &'plp snv_hts::Entry,
    read: &'plp snv_hts::Read,
}

impl std::fmt::Debug for Entry<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Entry")
            .field("qpos", &self.qpos())
            .field("indel", &self.indel())
            .field("is_del", &self.is_del())
            .field("is_refskip", &self.is_refskip())
            .finish()
    }
}

impl Entry<'_> {
    /// The read skips this position with an `N` CIGAR operation.
    pub fn is_refskip(&self) -> bool {
        self.entry.is_refskip
    }

    /// The read has a deletion at this position.
    pub fn is_del(&self) -> bool {
        self.entry.is_del
    }

    /// Length of the indel after this position: positive for insertions, negative for deletions.
    pub fn indel(&self) -> i32 {
        self.entry.indel
    }

    /// Query position of this column in the read.
    pub fn qpos(&self) -> i32 {
        i32::try_from(self.entry.qpos).expect("query position fits i32")
    }

    /// The read maps to the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.read.flags() & FLAG_REVERSE != 0
    }

    /// 4-bit `seq_nt16` code of query base `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is outside the read.
    pub fn base_code(&self, i: i32) -> u8 {
        Read(self.read).base_code(i)
    }
}

/// An indexed FASTA file (`faidx_t`).
#[derive(Debug)]
pub struct Fasta(snv_hts::Fasta);

impl Fasta {
    /// Opens `path` and its `.fai` index.
    ///
    /// # Errors
    ///
    /// Fails if the FASTA or its index cannot be loaded.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        Ok(Self(snv_hts::Fasta::open(path)?))
    }

    /// Length of `contig` from the `.fai` index.
    ///
    /// # Errors
    ///
    /// Fails if the contig is not in the index.
    pub fn contig_len(&self, contig: &str) -> Result<i64, Error> {
        Ok(self.0.contig_len(contig)?)
    }

    /// `faidx_fetch_seq(fai, contig, beg, end)`: bases `beg..=end` (0-based), with htslib's
    /// clamping of `beg < 0` to 0 and of `end` to the contig end.
    ///
    /// # Errors
    ///
    /// Fails if the contig is not in the index.
    pub fn fetch(&self, contig: &str, beg: i32, end: i32) -> Result<Vec<u8>, Error> {
        Ok(self.0.fetch(contig, i64::from(beg), i64::from(end))?)
    }
}

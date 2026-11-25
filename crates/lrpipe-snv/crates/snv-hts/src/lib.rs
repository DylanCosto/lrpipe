//! htslib's alignment reading as the original uses it, on [pasta]: indexed BAM and CRAM
//! region queries (`sam_itr_querys`), records as htslib holds them (`bam1_t`, here pasta's encoded
//! [`Read`]), the `bam_plp` pileup engine, SAM text in and out (`sam_parse1`, `sam_format1`),
//! and FASTA access (`faidx_fetch_seq`).
//!
//! pasta's pileup is a port of htslib's, and its queries return the reads htslib's do (a placed
//! unmapped read spans one base, as with `bam_endpos`), so callers see the same reads in the
//! same order by default.

use std::cell::RefCell;
use std::fs::File;
use std::io::{self, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use pasta::bam;
use pasta::bgzf;
use pasta::core::region::Interval;
use pasta::core::{Position, Region};
use pasta::cram;
use pasta::cram::io::reader::{Fields, MdNm, Tags};
use pasta::fasta;
use pasta::sam;
use pasta::sam::alignment::RecordBuf;
use pasta::sam::alignment::io::Write as _;
use pasta::sam::alignment::record::data::field::Tag;

pub use pasta::pileup::{Entry, Overlaps, Read};

/// Why reading failed.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A file could not be opened or its header read.
    #[error("cannot open {path:?}")]
    Open {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// An alignment index could not be loaded.
    #[error("cannot load the index of {path:?}")]
    Index {
        /// The alignment file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A region string names no contig of the header, or has an invalid range.
    #[error("invalid region {region:?}")]
    Region {
        /// The region string.
        region: String,
    },
    /// A record could not be read or decoded.
    #[error("failed to read alignment record (BAM/CRAM decode error)")]
    Decode(#[source] io::Error),
    /// A SAM line could not be parsed.
    #[error("cannot parse SAM line {line:?}")]
    Parse {
        /// The line.
        line: String,
    },
    /// The pileup engine was given reads out of order.
    #[error("pileup failed")]
    Pileup(#[source] pasta::pileup::Error),
    /// A contig is missing from the reference.
    #[error("contig {contig:?} is not in the reference")]
    Contig {
        /// The contig name.
        contig: String,
    },
    /// A FASTA could not be read.
    #[error("cannot read {path:?}")]
    Fasta {
        /// The FASTA.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
}

/// `BAM_FUNMAP`.
pub const FLAG_UNMAPPED: u16 = 0x4;

/// One past the last reference position the read covers (`bam_endpos`): a read without
/// reference-consuming operations, or unmapped, covers its first position.
#[must_use]
pub fn end_pos(read: &Read) -> i64 {
    if read.flags() & FLAG_UNMAPPED != 0 {
        return read.position() + 1;
    }
    let span: i64 =
        read.cigar().iter().filter(|&&op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8)).map(|&op| i64::from(op >> 4)).sum();
    read.position() + span.max(1)
}

enum Inner {
    Bam(Box<bam::io::IndexedReader<bgzf::io::Reader<File>>>),
    Cram(Box<cram::io::IndexedReader<File>>),
}

/// An indexed BAM or CRAM file with its header.
pub struct Alignments {
    path: PathBuf,
    /// The FASTA CRAM is decoded with.
    reference: PathBuf,
    header: sam::Header,
    names: Vec<String>,
    lens: Vec<i64>,
    inner: Inner,
    canonical_read_order: bool,
}

impl std::fmt::Debug for Alignments {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Alignments").field("path", &self.path).finish_non_exhaustive()
    }
}

fn is_cram(path: &Path) -> bool {
    path.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"))
}

/// A FASTA as a CRAM reference: loaded a contig at a time.
fn repository(reference: &Path) -> io::Result<fasta::Repository> {
    let reader = fasta::io::indexed_reader::Builder::default().build_from_path(reference)?;
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader)))
}

/// Reference repositories shared by the CRAM queries of a process, one for each FASTA and
/// contig queried. lrpipe-snv opens a reader for each chunk, and a repository of each reader's
/// own would load the whole contig again for each, and hash it again for pasta's check of
/// slices without a reference checksum against `@SQ` M5. A contig's repository is kept while a
/// query holds it and for the last [`IDLE`](shared::IDLE) contigs no query holds, then dropped
/// with its sequences.
mod shared {
    use std::collections::HashMap;
    use std::io;
    use std::path::{Path, PathBuf};
    use std::sync::{LazyLock, Mutex, MutexGuard, PoisonError};

    use pasta::fasta;

    /// Contigs kept loaded with no query holding them, for the next chunks of the same contigs.
    pub(super) const IDLE: usize = 2;

    #[derive(Default)]
    struct Repositories {
        entries: HashMap<(PathBuf, String), Entry>,
        /// Counts releases, to find the contig held least recently.
        clock: u64,
    }

    struct Entry {
        repository: fasta::Repository,
        holders: usize,
        released: u64,
    }

    static REPOSITORIES: LazyLock<Mutex<Repositories>> = LazyLock::new(Mutex::default);

    // The map stays whole whatever panicked while it was locked.
    fn repositories() -> MutexGuard<'static, Repositories> {
        REPOSITORIES.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// A query's hold on the shared repository of `reference` for one contig.
    pub(super) struct Lease {
        key: (PathBuf, String),
        repository: fasta::Repository,
    }

    impl Lease {
        /// Holds the repository of `reference` for `contig`, opening it if none is kept.
        pub(super) fn new(reference: &Path, contig: &str) -> io::Result<Self> {
            let key = (reference.to_owned(), contig.to_owned());
            let mut repositories = repositories();
            let entry = match repositories.entries.entry(key.clone()) {
                std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(Entry { repository: super::repository(reference)?, holders: 0, released: 0 })
                }
            };
            entry.holders += 1;
            Ok(Self { key, repository: entry.repository.clone() })
        }

        pub(super) fn repository(&self) -> &fasta::Repository {
            &self.repository
        }
    }

    impl Drop for Lease {
        fn drop(&mut self) {
            let mut repositories = repositories();
            repositories.clock += 1;
            let clock = repositories.clock;
            if let Some(entry) = repositories.entries.get_mut(&self.key) {
                entry.holders -= 1;
                entry.released = clock;
            }
            let mut idle: Vec<(u64, (PathBuf, String))> = repositories
                .entries
                .iter()
                .filter(|(_, e)| e.holders == 0)
                .map(|(k, e)| (e.released, k.clone()))
                .collect();
            if idle.len() > IDLE {
                idle.sort_unstable();
                for (_, key) in &idle[..idle.len() - IDLE] {
                    repositories.entries.remove(key);
                }
            }
        }
    }

    /// How many contigs' repositories are kept (for tests).
    #[cfg(test)]
    pub(super) fn kept(reference: &Path) -> usize {
        repositories().entries.keys().filter(|(r, _)| r == reference).count()
    }
}

impl Alignments {
    /// Opens `path` (BAM with `.bai` or `.csi`, CRAM with `.crai`) and its index, decoding CRAM
    /// with `reference`.
    ///
    /// # Errors
    ///
    /// If the file, its index or its header cannot be read.
    pub fn open(path: impl AsRef<Path>, reference: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let open = |source| Error::Open { path: path.to_owned(), source };
        let index = |source| Error::Index { path: path.to_owned(), source };
        let (header, inner) = if is_cram(path) {
            let index = cram::fs::read_associated_index(path).map_err(index)?;
            let mut r = cram::io::indexed_reader::Builder::default()
                .set_index(index)
                .set_reference_sequence_repository(repository(reference.as_ref()).map_err(open)?)
                .build_from_path(path)
                .map_err(open)?;
            (r.read_header().map_err(open)?, Inner::Cram(Box::new(r)))
        } else {
            let index = bam::fs::read_associated_index(path).map_err(index)?;
            let mut r =
                bam::io::indexed_reader::Builder::default().set_index(index).build_from_path(path).map_err(open)?;
            (r.read_header().map_err(open)?, Inner::Bam(Box::new(r)))
        };
        let names = header.reference_sequences().keys().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
        let lens = header
            .reference_sequences()
            .values()
            .map(|m| i64::try_from(usize::from(m.length())).unwrap_or(i64::MAX))
            .collect();
        Ok(Self {
            path: path.to_owned(),
            reference: reference.as_ref().to_owned(),
            header,
            names,
            lens,
            inner,
            canonical_read_order: false,
        })
    }

    /// Order alignments that start at the same coordinate before returning them from queries.
    #[must_use]
    pub fn with_canonical_read_order(mut self, enabled: bool) -> Self {
        self.canonical_read_order = enabled;
        self
    }

    /// The file's header.
    #[must_use]
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Target id of `contig` in the header, if present.
    #[must_use]
    pub fn tid(&self, contig: &str) -> Option<i32> {
        self.header.reference_sequences().get_index_of(contig.as_bytes()).and_then(|i| i32::try_from(i).ok())
    }

    /// The name of target `tid`.
    #[must_use]
    pub fn target_name(&self, tid: i32) -> Option<&str> {
        usize::try_from(tid).ok().and_then(|i| self.names.get(i)).map(String::as_str)
    }

    /// The contig names in target-id order.
    #[must_use]
    pub fn target_names(&self) -> &[String] {
        &self.names
    }

    /// The contig lengths in target-id order.
    #[must_use]
    pub fn target_lens(&self) -> &[i64] {
        &self.lens
    }

    /// Mapped reads on `contig` as the index counts them (`samtools idxstats`): the contig's own
    /// count; 0 for a contig a BAI or CSI has no bins for, as every read is in one; `None` for a
    /// contig not in the header, or with bins but without its optional count (LRP-005), and for
    /// a CRAM index, which has no counts. A count elsewhere in the index says nothing of this
    /// contig's: the counts are optional per contig.
    #[must_use]
    pub fn mapped_reads(&self, contig: &str) -> Option<u64> {
        let tid = usize::try_from(self.tid(contig)?).ok()?;
        match &self.inner {
            Inner::Cram(_) => None,
            Inner::Bam(r) => {
                let index = r.index();
                match index.reference_sequences().nth(tid)?.metadata() {
                    Some(m) => Some(m.mapped_record_count()),
                    // A contig whose reads are in no chunk has none.
                    None if index.query(tid, Interval::from(..)).ok()?.is_empty() => Some(0),
                    None => None,
                }
            }
        }
    }

    /// The reads `sam_itr_querys` returns for `region`, a samtools-style region string
    /// (`contig`, `contig:beg`, `contig:beg-end`, 1-based and inclusive, `,` allowed in numbers).
    ///
    /// # Errors
    ///
    /// If the region names no contig, has an invalid range, or cannot be queried.
    pub fn query(&mut self, region: &str) -> Result<Records<'_>, Error> {
        let bad = || Error::Region { region: region.to_owned() };
        let (name, beg, end) = parse_region(region, &self.header).ok_or_else(bad)?;
        self.query_interval_named(&name, beg, end, region)
    }

    /// The reads overlapping `[beg, end)` (0-based) of `contig`, as `sam_itr_queryi` returns
    /// them (rust-htslib's `fetch((tid, beg, end))`).
    ///
    /// # Errors
    ///
    /// If the contig is not in the header, or the range cannot be queried.
    pub fn query_interval(&mut self, contig: &str, beg: i64, end: i64) -> Result<Records<'_>, Error> {
        let region = format!("{contig}:{}-{end}", beg + 1);
        if self.tid(contig).is_none() || end <= beg.max(0) {
            return Err(Error::Region { region });
        }
        self.query_interval_named(contig, beg.max(0), Some(end), &region)
    }

    fn query_interval_named(
        &mut self,
        name: &str,
        beg: i64,
        end: Option<i64>,
        region: &str,
    ) -> Result<Records<'_>, Error> {
        let bad = || Error::Region { region: region.to_owned() };
        let start = Position::new(usize::try_from(beg).map_err(|_| bad())? + 1).ok_or_else(bad)?;
        let query_region = match end {
            Some(end) => {
                Region::new(name, start..=Position::new(usize::try_from(end).map_err(|_| bad())?).ok_or_else(bad)?)
            }
            None => Region::new(name, start..),
        };
        let Self { header, inner, names, lens, reference, canonical_read_order, .. } = self;
        let inner = match inner {
            Inner::Bam(r) => {
                RecordsInner::Bam(r.query(header, &query_region).map_err(Error::Decode)?, bam::Record::default())
            }
            Inner::Cram(r) => {
                let id = header.reference_sequences().get_index_of(name.as_bytes()).ok_or_else(bad)?;
                RecordsInner::Cram(CramQuery {
                    reader: r,
                    reference: shared::Lease::new(reference, name).map_err(Error::Decode)?,
                    header,
                    id,
                    interval: query_region.interval(),
                    entry: 0,
                    pending: Vec::new().into_iter(),
                })
            }
        };
        Ok(Records { inner, beg, end, names, lens, tied_reads: canonical_read_order.then(TiedReads::default) })
    }
}

impl Alignments {
    /// The file's first `n` records from its start, in file order, whatever their contig
    /// (`sam_read1` from just after the header).
    ///
    /// # Errors
    ///
    /// If a record cannot be read or decoded.
    pub fn first_reads(&mut self, n: usize) -> Result<Vec<Read>, Error> {
        let mut reads = Vec::new();
        match &mut self.inner {
            Inner::Bam(r) => {
                for record in r.records().take(n) {
                    reads.push(Read::new(record.map_err(Error::Decode)?.as_ref().to_vec()).map_err(Error::Decode)?);
                }
            }
            Inner::Cram(r) => {
                let mut encoder = bam::io::Writer::from(Vec::new());
                for record in r.records(&self.header).take(n) {
                    let raw =
                        encode(&mut encoder, &self.header, &record.map_err(Error::Decode)?).map_err(Error::Decode)?;
                    reads.push(Read::new(raw).map_err(Error::Decode)?);
                }
            }
        }
        Ok(reads)
    }
}

/// `hts_parse_region` for the names in `header`: the contig, the 0-based start and the
/// exclusive end (`None` for the contig's end).
fn parse_region(region: &str, header: &sam::Header) -> Option<(String, i64, Option<i64>)> {
    let contigs = header.reference_sequences();
    if contigs.contains_key(region.as_bytes()) {
        return Some((region.to_owned(), 0, None));
    }
    let (name, range) = region.rsplit_once(':')?;
    if !contigs.contains_key(name.as_bytes()) {
        return None;
    }
    let number = |s: &str| -> Option<i64> { s.replace(',', "").parse().ok() };
    let (beg, end) = match range.split_once('-') {
        None => (number(range)?, None),
        Some((b, "")) => (number(b)?, None),
        Some((b, e)) => (if b.is_empty() { 1 } else { number(b)? }, Some(number(e)?)),
    };
    let beg = (beg - 1).max(0);
    if end.is_some_and(|e| e <= beg) {
        return None;
    }
    Some((name.to_owned(), beg, end))
}

enum RecordsInner<'a> {
    Bam(bam::io::reader::Query<'a, bgzf::io::Reader<File>>, bam::Record),
    Cram(CramQuery<'a>),
}

/// A CRAM region query that encodes the records straight to BAM: pasta's query would build a
/// `RecordBuf` for each, to be written out as BAM again. It gives the records pasta's query
/// does, in the same order, with only what lrpipe-snv reads of them (see `read_slice`).
struct CramQuery<'a> {
    reader: &'a mut cram::io::IndexedReader<File>,
    /// The contig's reference, shared with the process's other queries of it.
    reference: shared::Lease,
    header: &'a sam::Header,
    id: usize,
    interval: Interval,
    /// The next index entry to look at.
    entry: usize,
    /// The current slice's records in the query, encoded.
    pending: std::vec::IntoIter<Vec<u8>>,
}

impl CramQuery<'_> {
    /// The next record, encoded.
    fn next(&mut self) -> io::Result<Option<Vec<u8>>> {
        loop {
            if let Some(raw) = self.pending.next() {
                return Ok(Some(raw));
            }
            if !self.read_slice()? {
                return Ok(None);
            }
        }
    }

    /// Decodes the next slice the index has in the query into `pending`; false at the end.
    fn read_slice(&mut self) -> io::Result<bool> {
        let index = self.reader.index();
        let Some(offset) = index[self.entry..].iter().position(|x| slice_overlaps(x, self.id, self.interval)) else {
            self.entry = index.len();
            return Ok(false);
        };
        self.entry += offset;
        let (position, landmark) = (index[self.entry].offset(), index[self.entry].landmark());
        self.entry += 1;
        self.reader.get_mut().seek(SeekFrom::Start(position))?;
        let mut container = cram::io::reader::Container::default();
        if self.reader.read_container(&mut container)? == 0 {
            self.entry = self.reader.index().len();
            return Ok(false);
        }
        let compression_header = container.compression_header()?;
        let slice = container
            .header()
            .landmarks()
            .iter()
            .position(|&n| n as u64 == landmark)
            .and_then(|i| container.slices().nth(i))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid slice landmark"))??;
        // Only what the callers read: every field but the read group, and of the tags only the
        // move table (`mv`, for dwell time). Nothing reads MD or NM, and leaving out the other
        // tags skips most of an ONT CRAM, its `MM`/`ML` methylation calls.
        let fields = Fields::all() - Fields::READ_GROUP;
        let tags = Tags::Only(vec![Tag::new(b'm', b'v')]);
        let (core, external) = slice.decode_blocks_for(&compression_header, fields, &tags)?;
        let records = slice.records_with_fields(
            self.reference.repository().clone(),
            self.header,
            &compression_header,
            &core,
            &external,
            MdNm::AsStored,
            fields,
            &tags,
        )?;
        let mut pending = Vec::with_capacity(records.len());
        for record in &records {
            let mut raw = Vec::new();
            record.encode_bam_record(self.header, &mut raw)?;
            let read = Read::new(raw)?;
            if read_overlaps(&read, self.id, self.interval) {
                pending.push(read.into_bytes());
            }
        }
        self.pending = pending.into_iter();
        Ok(true)
    }
}

/// Whether an index entry's slice can hold records of `id` overlapping `interval` (pasta's
/// query's test).
fn slice_overlaps(entry: &cram::crai::Record, id: usize, interval: Interval) -> bool {
    if entry.reference_sequence_id() != Some(id) {
        return false;
    }
    if interval.start().is_none() && interval.end().is_none() {
        return true;
    }
    let Some(start) = entry.alignment_start() else { return false };
    let span = entry.alignment_span();
    if span == 0 {
        return false;
    }
    let end = start.checked_add(span - 1).expect("slice end overflows");
    interval.intersects((start..=end).into())
}

/// Whether a record lies on `id` and overlaps `interval` (pasta's query's test: an alignment
/// that consumes no reference covers its start).
fn read_overlaps(read: &Read, id: usize, interval: Interval) -> bool {
    if usize::try_from(read.reference_sequence_id()).ok() != Some(id) {
        return false;
    }
    if interval.start().is_none() && interval.end().is_none() {
        return true;
    }
    let Some(start) = usize::try_from(read.position() + 1).ok().and_then(Position::new) else { return false };
    let span: usize =
        read.cigar().iter().filter(|&&op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8)).map(|&op| (op >> 4) as usize).sum();
    let end = start.checked_add(span.max(1) - 1).expect("alignment end overflows");
    interval.intersects((start..=end).into())
}

/// The reads of a region query, in coordinate order.
pub struct Records<'a> {
    inner: RecordsInner<'a>,
    beg: i64,
    end: Option<i64>,
    names: &'a [String],
    lens: &'a [i64],
    tied_reads: Option<TiedReads>,
}

impl std::fmt::Debug for Records<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Records").field("beg", &self.beg).field("end", &self.end).finish_non_exhaustive()
    }
}

impl<'a> Records<'a> {
    /// The query's 0-based start and exclusive end (`None` for the contig's end), as
    /// `hts_itr_t`'s `beg` and `end`.
    #[must_use]
    pub fn range(&self) -> (i64, Option<i64>) {
        (self.beg, self.end)
    }

    /// The file's contig names in target-id order.
    #[must_use]
    pub fn target_names(&self) -> &'a [String] {
        self.names
    }

    /// The file's contig lengths in target-id order.
    #[must_use]
    pub fn target_lens(&self) -> &'a [i64] {
        self.lens
    }
}

/// A source of reads in coordinate order, for [`Pileup`].
pub trait ReadSource {
    /// The next read, or `None` at the end.
    ///
    /// # Errors
    ///
    /// If a read cannot be read or decoded.
    fn next_read(&mut self) -> Result<Option<Read>, Error>;
}

impl ReadSource for Records<'_> {
    fn next_read(&mut self) -> Result<Option<Read>, Error> {
        match &mut self.tied_reads {
            Some(tied_reads) => tied_reads.next_read(&mut self.inner),
            None => self.inner.next_read(),
        }
    }
}

impl ReadSource for RecordsInner<'_> {
    fn next_read(&mut self) -> Result<Option<Read>, Error> {
        let raw = match self {
            Self::Bam(query, record) => {
                if query.read_record(record).map_err(Error::Decode)? == 0 {
                    return Ok(None);
                }
                record.as_ref().to_vec()
            }
            Self::Cram(query) => match query.next().map_err(Error::Decode)? {
                Some(raw) => raw,
                None => return Ok(None),
            },
        };
        Read::new(raw).map(Some).map_err(Error::Decode)
    }
}

/// Buffer one coordinate's reads so a different order of ties cannot change the caller's
/// read deduplication or full-alignment sampling.
#[derive(Default)]
struct TiedReads {
    reads: Vec<Read>,
    lookahead: Option<Read>,
}

impl TiedReads {
    fn next_read(&mut self, source: &mut impl ReadSource) -> Result<Option<Read>, Error> {
        if let Some(read) = self.reads.pop() {
            return Ok(Some(read));
        }
        let first = match self.lookahead.take() {
            Some(first) => first,
            None => match source.next_read()? {
                Some(first) => first,
                None => return Ok(None),
            },
        };
        let coordinate = (first.reference_sequence_id(), first.position());
        self.reads.push(first);
        loop {
            match source.next_read() {
                Ok(Some(read)) if (read.reference_sequence_id(), read.position()) == coordinate => {
                    self.reads.push(read)
                }
                Ok(next) => {
                    self.lookahead = next;
                    break;
                }
                Err(error) => {
                    self.reads.clear();
                    return Err(error);
                }
            }
        }
        // Popped from the end, so sorted in reverse.
        if self.reads.len() > LARGE_TIE {
            self.reads.sort_unstable_by(|a, b| compare_many_tied_reads(b, a));
        } else {
            self.reads.sort_unstable_by(|a, b| compare_tied_reads(b, a));
        }
        Ok(self.reads.pop())
    }
}

/// More reads than this starting at one position are an amplicon's or a targeted capture's, not
/// a whole genome's chance ties (two or three reads), and are ordered by
/// [`compare_many_tied_reads`].
const LARGE_TIE: usize = 64;

/// A read name's hash (FNV-1a): an order of names unrelated to what they say.
fn name_hash(name: &[u8]) -> u64 {
    name.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3))
}

/// The order of a large tie: by a hash of the read name first, then as [`compare_tied_reads`].
/// The stages that take the first reads up to a depth then take them as a coordinate-sorted
/// file has them, in no order of strand, flow cell or time; by strand and name they would take
/// the forward strand's reads alone, the earliest named first.
fn compare_many_tied_reads(a: &Read, b: &Read) -> std::cmp::Ordering {
    name_hash(a.name()).cmp(&name_hash(b.name())).then_with(|| compare_tied_reads(a, b))
}

fn compare_tied_reads(a: &Read, b: &Read) -> std::cmp::Ordering {
    // Sort by strand and read name first. Names can be shared by multiple alignments. Use
    // the fields that can distinguish their features, avoiding format-specific auxiliary tags.
    (a.is_reverse_complemented(), a.name())
        .cmp(&(b.is_reverse_complemented(), b.name()))
        .then_with(|| a.flags().cmp(&b.flags()))
        .then_with(|| b.mapping_quality().cmp(&a.mapping_quality()))
        .then_with(|| a.cigar().cmp(b.cigar()))
        .then_with(|| a.quality_scores().cmp(b.quality_scores()))
        .then_with(|| a.tag(*b"mv").cmp(&b.tag(*b"mv")))
        .then_with(|| packed_seq(a).cmp(packed_seq(b)))
}

/// `record` as an encoded BAM record, without its `block_size`.
fn encode(encoder: &mut bam::io::Writer<Vec<u8>>, header: &sam::Header, record: &RecordBuf) -> io::Result<Vec<u8>> {
    encoder.get_mut().clear();
    encoder.write_alignment_record(header, record)?;
    Ok(encoder.get_ref()[size_of::<u32>()..].to_vec())
}

/// Where the stored CIGAR, sequence and auxiliary fields of an encoded record start.
fn layout(read: &Read) -> (usize, usize) {
    let raw = read.as_bytes();
    let name_len = usize::from(raw[8]);
    let n_cigar = usize::from(u16::from_le_bytes([raw[12], raw[13]]));
    let seq = 32 + name_len + 4 * n_cigar;
    let aux = seq + read.sequence_len().div_ceil(2) + read.sequence_len();
    (seq, aux)
}

/// The query bases as `seq_nt16` codes, two per byte, high nibble first (`bam_get_seq`).
#[must_use]
pub fn packed_seq(read: &Read) -> &[u8] {
    let (seq, _) = layout(read);
    &read.as_bytes()[seq..seq + read.sequence_len().div_ceil(2)]
}

/// An auxiliary field's value.
#[derive(Clone, Debug, PartialEq)]
pub enum Aux<'a> {
    /// `A`.
    Char(u8),
    /// `c`, `C`, `s`, `S`, `i` or `I`.
    Int(i64),
    /// `f`.
    Float(f32),
    /// `d`, HTSlib's non-standard double, which samtools writes to BAM.
    Double(f64),
    /// `Z`, without its NUL.
    Text(&'a [u8]),
    /// `H`, without its NUL.
    Hex(&'a [u8]),
    /// `B` with an integer subtype.
    Ints(Vec<i64>),
    /// `B:f`.
    Floats(Vec<f32>),
}

/// The read's auxiliary fields in order (`bam_aux_first`/`bam_aux_next`); parsing stops at a
/// malformed field.
#[must_use]
pub fn aux_fields(read: &Read) -> Vec<([u8; 2], Aux<'_>)> {
    let raw = read.as_bytes();
    let (_, mut at) = layout(read);
    let mut fields = Vec::new();
    while at + 3 <= raw.len() {
        let tag = [raw[at], raw[at + 1]];
        let ty = raw[at + 2];
        let v = &raw[at + 3..];
        let int = |size: usize, signed: bool, b: &[u8]| -> Option<i64> {
            let b = b.get(..size)?;
            Some(match (size, signed) {
                (1, true) => i64::from(b[0].cast_signed()),
                (1, false) => i64::from(b[0]),
                (2, true) => i64::from(i16::from_le_bytes([b[0], b[1]])),
                (2, false) => i64::from(u16::from_le_bytes([b[0], b[1]])),
                (4, true) => i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])),
                _ => i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
            })
        };
        let width = |t: u8| match t {
            b'c' | b'C' | b'A' => Some((1, t == b'c')),
            b's' | b'S' => Some((2, t == b's')),
            b'i' | b'I' | b'f' => Some((4, t == b'i')),
            _ => None,
        };
        let (value, len) = match ty {
            b'A' => match v.first() {
                Some(&c) => (Aux::Char(c), 1),
                None => break,
            },
            b'f' => match v.get(..4) {
                Some(b) => (Aux::Float(f32::from_le_bytes([b[0], b[1], b[2], b[3]])), 4),
                None => break,
            },
            b'd' => match v.first_chunk::<8>() {
                Some(b) => (Aux::Double(f64::from_le_bytes(*b)), 8),
                None => break,
            },
            b'Z' | b'H' => {
                let Some(end) = v.iter().position(|&b| b == 0) else { break };
                (if ty == b'Z' { Aux::Text(&v[..end]) } else { Aux::Hex(&v[..end]) }, end + 1)
            }
            b'B' => {
                let (Some(&sub), Some(count)) = (v.first(), v.get(1..5)) else { break };
                let count = u32::from_le_bytes([count[0], count[1], count[2], count[3]]) as usize;
                let Some((size, signed)) = width(sub) else { break };
                let Some(data) = v.get(5..5 + size * count) else { break };
                let value = if sub == b'f' {
                    Aux::Floats(data.as_chunks::<4>().0.iter().map(|b| f32::from_le_bytes(*b)).collect())
                } else {
                    Aux::Ints(data.chunks_exact(size).filter_map(|b| int(size, signed, b)).collect())
                };
                (value, 5 + size * count)
            }
            t => match width(t).and_then(|(size, signed)| int(size, signed, v).map(|x| (x, size))) {
                Some((x, size)) => (Aux::Int(x), size),
                None => break,
            },
        };
        fields.push((tag, value));
        at += 3 + len;
    }
    fields
}

/// The first auxiliary field named `tag` (`bam_aux_get`).
#[must_use]
pub fn aux<'r>(read: &'r Read, tag: [u8; 2]) -> Option<Aux<'r>> {
    aux_fields(read).into_iter().find(|(t, _)| *t == tag).map(|(_, v)| v)
}

/// `sam_parse1`: a SAM line as a read, with `header`'s contigs.
///
/// # Errors
///
/// If the line is not a valid SAM record for `header`.
pub fn parse_sam(header: &sam::Header, line: &str) -> Result<Read, Error> {
    let bad = || Error::Parse { line: line.to_owned() };
    let record = sam::Record::try_from(line.as_bytes()).map_err(|_| bad())?;
    let record = RecordBuf::try_from_alignment_record(header, &record).map_err(|_| bad())?;
    let raw = encode(&mut bam::io::Writer::from(Vec::new()), header, &record).map_err(|_| bad())?;
    Read::new(raw).map_err(|_| bad())
}

/// `sam_format1`: `read` as `samtools view` prints it, without the newline.
///
/// # Errors
///
/// If the read cannot be decoded.
pub fn format_sam(header: &sam::Header, read: &Read) -> Result<String, Error> {
    let record = bam::Record::try_from(read.as_bytes().to_vec()).map_err(Error::Decode)?;
    let mut writer = sam::io::Writer::new(Vec::new());
    writer.write_alignment_record(header, &record).map_err(Error::Decode)?;
    let mut text = writer.into_inner();
    if text.last() == Some(&b'\n') {
        text.pop();
    }
    String::from_utf8(text).map_err(|e| Error::Decode(io::Error::new(io::ErrorKind::InvalidData, e)))
}

/// A running pileup over the reads of a [`ReadSource`] that `filter` keeps (htslib's `bam_plp`
/// with a read callback).
pub struct Pileup<S, F> {
    source: S,
    filter: F,
    engine: pasta::pileup::Pileup,
    /// Whether the engine may have positions ready since the last push.
    draining: bool,
}

impl<S, F> std::fmt::Debug for Pileup<S, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Pileup").finish_non_exhaustive()
    }
}

impl<S: ReadSource, F: FnMut(&mut Read) -> bool> Pileup<S, F> {
    /// A pileup of `source`'s reads that `filter` keeps (it may also change them), holding at
    /// most `max_depth` reads that start at one position (`bam_plp_set_maxcnt`), with mates'
    /// overlaps adjusted or not (`bam_plp_init_overlaps`).
    pub fn new(source: S, max_depth: usize, overlaps: Overlaps, filter: F) -> Self {
        Self { source, filter, engine: pasta::pileup::Pileup::new(max_depth, overlaps), draining: false }
    }

    /// The source, to read its state between columns.
    pub fn source(&self) -> &S {
        &self.source
    }

    /// The column [`next_column`](Self::next_column) last returned, again.
    pub fn current(&self) -> Column<'_> {
        let (tid, pos) = self.engine.position();
        Column { tid, pos, engine: &self.engine }
    }

    /// Advances to the next column (`bam_plp_auto`).
    ///
    /// # Errors
    ///
    /// If a read cannot be read, or the reads are out of order.
    pub fn next_column(&mut self) -> Result<Option<Column<'_>>, Error> {
        loop {
            if self.draining {
                if self.engine.next_position().map_err(Error::Pileup)? {
                    let (tid, pos) = self.engine.position();
                    return Ok(Some(Column { tid, pos, engine: &self.engine }));
                }
                self.draining = false;
            }
            if self.engine.is_eof() {
                return Ok(None);
            }
            let read = loop {
                match self.source.next_read()? {
                    None => break None,
                    Some(mut read) => {
                        if (self.filter)(&mut read) {
                            break Some(read);
                        }
                    }
                }
            };
            self.engine.push(read).map_err(Error::Pileup)?;
            self.draining = true;
        }
    }
}

/// One pileup column: every kept read overlapping reference position `pos`.
#[derive(Clone, Copy)]
pub struct Column<'p> {
    /// Target id of the contig.
    pub tid: i32,
    /// 0-based reference position.
    pub pos: i64,
    engine: &'p pasta::pileup::Pileup,
}

impl std::fmt::Debug for Column<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Column").field("tid", &self.tid).field("pos", &self.pos).finish_non_exhaustive()
    }
}

impl<'p> Column<'p> {
    /// The reads in this column, in htslib's order.
    #[must_use]
    pub fn entries(&self) -> &'p [Entry] {
        self.engine.entries()
    }

    /// The read of an entry of this column.
    #[must_use]
    pub fn read(&self, entry: &Entry) -> &'p Read {
        self.engine.read(entry)
    }
}

/// An indexed FASTA file (`faidx_t`), plain or bgzipped (with `.gzi`).
pub struct Fasta {
    path: PathBuf,
    index: fasta::fai::Index,
    reader: RefCell<fasta::io::IndexedReader<fasta::io::BufReader<File>>>,
}

impl std::fmt::Debug for Fasta {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Fasta").field("path", &self.path).finish_non_exhaustive()
    }
}

impl Fasta {
    /// Opens `path` and its `.fai` index.
    ///
    /// # Errors
    ///
    /// If the FASTA or its index cannot be loaded.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|source| Error::Open { path: path.to_owned(), source })?;
        let index = reader.index().clone();
        Ok(Self { path: path.to_owned(), index, reader: RefCell::new(reader) })
    }

    /// The contigs' names and lengths, in index order.
    pub fn contigs(&self) -> impl Iterator<Item = (&[u8], u64)> {
        self.index.as_ref().iter().map(|r| (r.name().as_ref(), r.length()))
    }

    /// Length of `contig` from the `.fai` index.
    ///
    /// # Errors
    ///
    /// If the contig is not in the index.
    pub fn contig_len(&self, contig: &str) -> Result<i64, Error> {
        self.index
            .as_ref()
            .iter()
            .find(|r| AsRef::<[u8]>::as_ref(r.name()) == contig.as_bytes())
            .and_then(|r| i64::try_from(r.length()).ok())
            .ok_or_else(|| Error::Contig { contig: contig.to_owned() })
    }

    /// `faidx_fetch_seq(fai, contig, beg, end)`: bases `beg..=end` (0-based) as stored, with
    /// htslib's clamping (`end < beg` fetches from `end`; `beg` to `0..=len`, `end` to
    /// `0..len`).
    ///
    /// # Errors
    ///
    /// If the contig is not in the index or cannot be read.
    pub fn fetch(&self, contig: &str, beg: i64, end: i64) -> Result<Vec<u8>, Error> {
        let len = self.contig_len(contig)?;
        let mut beg = if end < beg { end } else { beg };
        let mut end = end;
        beg = if beg < 0 { 0 } else { beg.min(len) };
        end = if end < 0 {
            0
        } else if len <= end {
            len - 1
        } else {
            end
        };
        // `fai_retrieve(beg, end + 1)`: nothing when the range is empty.
        if end < beg {
            return Ok(Vec::new());
        }
        let start = Position::new(usize::try_from(beg + 1).expect("positive")).expect("positive");
        let stop = Position::new(usize::try_from(end + 1).expect("positive")).expect("positive");
        let read = |source| Error::Fasta { path: self.path.clone(), source };
        let record = self.reader.borrow_mut().query(&Region::new(contig, start..=stop)).map_err(read)?;
        Ok(record.sequence().as_ref().to_vec())
    }
}

/// `samtools index -b`: writes `path.bai` for the coordinate-sorted BAM at `path`. As with
/// samtools, the header need not say it is sorted; out-of-order records fail.
///
/// # Errors
///
/// If the BAM cannot be read, is not sorted, or the index cannot be written.
pub fn index_bam(path: &Path) -> io::Result<()> {
    use pasta::csi::binning_index::index::reference_sequence::bin::Chunk;
    use pasta::sam::alignment::Record as _;
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    let mut indexer =
        bam::index::Indexer::builder().build().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut record = bam::Record::default();
    let mut start = reader.get_ref().virtual_position();
    while reader.read_record(&mut record)? != 0 {
        let end = reader.get_ref().virtual_position();
        let context = match (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(id), Some(first), Some(last)) => Some((id, first, last, !record.flags().is_unmapped())),
            _ => None,
        };
        indexer
            .add_record(context, Chunk::new(start, end))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        start = end;
    }
    match indexer.build(header.reference_sequences().len()) {
        bam::Index::Bai(index) => bam::bai::fs::write(path.with_added_extension("bai"), &index),
        bam::Index::Csi(_) => Err(io::Error::other("contigs too long for a BAI")),
    }
}

/// `samtools faidx`: writes `path.fai` for the FASTA at `path`.
///
/// # Errors
///
/// If the FASTA cannot be read or the index written.
pub fn index_fasta(path: &Path) -> io::Result<()> {
    let index = fasta::fs::index(path)?;
    fasta::fai::fs::write(path.with_added_extension("fai"), &index)
}

/// `tabix -p vcf`: writes `path.tbi` for the bgzipped, coordinate-sorted VCF at `path`, as
/// htslib 1.24's `tbx_index_build` does (pasta's port of it): only the contig and position are
/// required, and the end is the longest of the reference allele, `INFO/END` and `SVLEN`.
///
/// # Errors
///
/// If the VCF cannot be read, a data line cannot be parsed or is out of order (htslib 1.24's
/// `tabix` fails there too), or the index cannot be written.
pub fn index_vcf(path: &Path) -> io::Result<()> {
    use io::BufRead as _;
    use pasta::variant::index::IndexFormat;
    use pasta::variant::index::tabix::{TabixBuilder, TabixConfig};
    let mut reader = bgzf::io::Reader::new(File::open(path)?);
    let mut builder = TabixBuilder::new(TabixConfig::VCF, 0);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
            if line.last() == Some(&b'\r') {
                line.pop();
            }
        }
        builder.push_line(&line, u64::from(reader.virtual_position()))?;
    }
    let tbx = builder.finish(u64::from(reader.virtual_position()));
    tbx.index.save(&path.with_added_extension("tbi"), IndexFormat::Tbi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pasta::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZero;

    fn header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence("chr20", Map::<ReferenceSequence>::new(NonZero::new(1000).expect("non-zero")))
            .add_reference_sequence("chr1:2", Map::<ReferenceSequence>::new(NonZero::new(500).expect("non-zero")))
            .build()
    }

    #[test]
    fn regions_parse_as_htslib_does() {
        let h = header();
        assert_eq!(parse_region("chr20", &h), Some(("chr20".into(), 0, None)));
        assert_eq!(parse_region("chr20:101-200", &h), Some(("chr20".into(), 100, Some(200))));
        assert_eq!(parse_region("chr20:1,001-2,000", &h), Some(("chr20".into(), 1000, Some(2000))));
        assert_eq!(parse_region("chr20:5", &h), Some(("chr20".into(), 4, None)));
        assert_eq!(parse_region("chr1:2", &h), Some(("chr1:2".into(), 0, None)));
        assert_eq!(parse_region("chr1:2:10-20", &h), Some(("chr1:2".into(), 9, Some(20))));
        assert_eq!(parse_region("chrX:1-10", &h), None);
        assert_eq!(parse_region("chr20:20-10", &h), None);
    }

    #[test]
    fn sam_lines_round_trip() {
        let h = header();
        let line = "r1\t0\tchr20\t11\t60\t3M1I2M\t*\t0\t0\tACGTAC\tIIIIII\tHP:i:1";
        let read = parse_sam(&h, line).expect("parses");
        assert_eq!((read.reference_sequence_id(), read.position(), read.mapping_quality()), (0, 10, 60));
        assert_eq!(read.integer_tag(*b"HP"), Some(1));
        assert_eq!(end_pos(&read), 15);
        assert_eq!(format_sam(&h, &read).expect("formats"), line);
        assert_eq!(packed_seq(&read), &[0x12, 0x48, 0x12]);
        let tagged = parse_sam(&h, "r2\t16\tchr20\t5\t9\t2M\t*\t0\t0\tAC\tII\tmv:B:c,5,1,0,1\tXf:f:1.5\tZZ:Z:x y")
            .expect("parses");
        assert_eq!(aux(&tagged, *b"mv"), Some(Aux::Ints(vec![5, 1, 0, 1])));
        // A `d` field is read, and does not hide the fields after it.
        let double = parse_sam(&h, "r3\t0\tchr20\t5\t9\t2M\t*\t0\t0\tAC\tII\tXD:d:2.5\tHP:i:2").expect("parses");
        assert_eq!(aux(&double, *b"XD"), Some(Aux::Double(2.5)));
        assert_eq!(aux(&double, *b"HP"), Some(Aux::Int(2)));
        assert_eq!(aux(&tagged, *b"Xf"), Some(Aux::Float(1.5)));
        assert_eq!(aux(&tagged, *b"ZZ"), Some(Aux::Text(b"x y")));
        assert_eq!(aux(&tagged, *b"HP"), None);
    }

    #[test]
    fn coordinate_ties_are_independent_of_file_order() {
        let h = header();
        let canonical = [
            "early\t16\tchr20\t9\t60\t2M\t*\t0\t0\tAC\tII",
            "dup\t0\tchr20\t10\t60\t2M\t*\t0\t0\tAC\tII\tmv:B:c,5,1,0",
            "dup\t0\tchr20\t10\t60\t2M\t*\t0\t0\tAG\tII\tmv:B:c,5,1,0",
            "dup\t0\tchr20\t10\t60\t2M\t*\t0\t0\tAC\tII\tmv:B:c,5,1,1",
            "dup\t0\tchr20\t10\t40\t2M\t*\t0\t0\tAC\tII",
            "a\t16\tchr20\t10\t60\t2M\t*\t0\t0\tAC\tII",
            "late\t0\tchr20\t11\t60\t2M\t*\t0\t0\tAC\tII",
        ];
        let file_order = [0, 5, 4, 3, 2, 1, 6].map(|i| canonical[i]);
        let dir = std::env::temp_dir().join(format!("snv-hts-coordinate-ties-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let path = dir.join("ties.bam");
        let mut writer = bam::io::Writer::new(File::create(&path).expect("BAM"));
        writer.write_header(&h).expect("header");
        for line in file_order {
            let record = sam::Record::try_from(line.as_bytes()).expect("record");
            writer.write_alignment_record(&h, &record).expect("write");
        }
        writer.try_finish().expect("finish BAM");
        drop(writer);
        index_bam(&path).expect("index");

        for (enabled, expected) in [(false, file_order), (true, canonical)] {
            let mut source = Alignments::open(&path, "").expect("open BAM").with_canonical_read_order(enabled);
            let mut records = source.query("chr20:9-11").expect("region");
            for expected in expected {
                let read = records.next_read().expect("valid source").expect("read");
                assert_eq!(format_sam(&h, &read).expect("SAM"), expected);
            }
            assert!(records.next_read().expect("valid source").is_none());
        }
        std::fs::remove_dir_all(&dir).expect("removed");
    }

    #[test]
    fn a_large_tie_is_not_taken_one_strand_first() {
        // An amplicon's reads all start at one position. By strand and name the first half would
        // be the forward strand's alone; by the hash of their names both strands come throughout,
        // in the same order whatever the file's.
        let h = header();
        let lines: Vec<String> = (0..2 * LARGE_TIE)
            .map(|i| format!("r{i:03}\t{}\tchr20\t10\t60\t2M\t*\t0\t0\tAC\tII", if i % 2 == 0 { 0 } else { 16 }))
            .collect();
        let dir = std::env::temp_dir().join(format!("snv-hts-large-tie-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let read_all = |name: &str, order: &mut dyn Iterator<Item = &String>| -> Vec<String> {
            let path = dir.join(name);
            let mut writer = bam::io::Writer::new(File::create(&path).expect("BAM"));
            writer.write_header(&h).expect("header");
            for line in order {
                let record = sam::Record::try_from(line.as_bytes()).expect("record");
                writer.write_alignment_record(&h, &record).expect("write");
            }
            writer.try_finish().expect("finish BAM");
            drop(writer);
            index_bam(&path).expect("index");
            let mut source = Alignments::open(&path, "").expect("open BAM").with_canonical_read_order(true);
            let mut records = source.query("chr20:10-10").expect("region");
            let mut out = Vec::new();
            while let Some(read) = records.next_read().expect("valid source") {
                out.push(format_sam(&h, &read).expect("SAM"));
            }
            out
        };
        let forward = read_all("forward.bam", &mut lines.iter());
        let reversed = read_all("reversed.bam", &mut lines.iter().rev());
        assert_eq!(forward.len(), 2 * LARGE_TIE);
        assert_eq!(forward, reversed, "the order does not follow the file's");
        let reverse_strand = forward[..LARGE_TIE].iter().filter(|l| l.split('\t').nth(1) == Some("16")).count();
        assert!((LARGE_TIE / 4..=3 * LARGE_TIE / 4).contains(&reverse_strand), "{reverse_strand} of the first half");
        std::fs::remove_dir_all(&dir).expect("removed");
    }

    /// Queries of one contig share its reference, and a contig no query holds is kept for the
    /// next ones only among the last few released.
    #[test]
    fn queries_share_a_contigs_reference_and_keep_few_idle() {
        let dir = std::env::temp_dir().join(format!("snv-hts-shared-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let fasta = dir.join("ref.fa");
        let (mut fa, mut fai) = (String::new(), String::new());
        for name in ["c1", "c2", "c3", "c4"] {
            fa.push_str(&format!(">{name}\n"));
            fai.push_str(&format!("{name}\t8\t{}\t8\t9\n", fa.len()));
            fa.push_str("ACGTACGT\n");
        }
        std::fs::write(&fasta, fa).expect("FASTA");
        std::fs::write(dir.join("ref.fa.fai"), fai).expect("FASTA index");

        let first = shared::Lease::new(&fasta, "c1").expect("a lease");
        let second = shared::Lease::new(&fasta, "c1").expect("a lease");
        first.repository().get(b"c1").expect("c1").expect("read");
        assert_eq!(second.repository().len(), 1, "the second query sees the contig the first loaded");
        assert_eq!(shared::kept(&fasta), 1);
        drop((first, second));
        for name in ["c2", "c3", "c4"] {
            drop(shared::Lease::new(&fasta, name).expect("a lease"));
            assert!(shared::kept(&fasta) <= shared::IDLE);
        }
        let held = shared::Lease::new(&fasta, "c1").expect("a lease");
        drop(shared::Lease::new(&fasta, "c2").expect("a lease"));
        assert_eq!(shared::kept(&fasta), shared::IDLE + 1, "a held contig is kept besides the idle ones");
        drop(held);
        assert_eq!(shared::kept(&fasta), shared::IDLE);
        std::fs::remove_dir_all(&dir).expect("removed");
    }
}

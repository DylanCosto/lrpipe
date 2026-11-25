//! The reading `samtools mpileup` and `samtools view` do, on pasta (through `snv-hts`).
//!
//! Columns come from the `bam_mplp` engine (pasta's port of htslib's), set up as `samtools
//! mpileup` sets it up (overlap detection, a depth cap of 8000), and fed through a filter like
//! `mplp_func` in samtools 1.24's `bam_plcmd.c`. Reads come from an indexed file (`-r REGION`)
//! or from SAM lines (`samtools mpileup -`).

use std::collections::VecDeque;
use std::ffi::CStr;
use std::path::Path;
use std::rc::Rc;

use snv_hts::{Alignments, Overlaps, Read, ReadSource, Records};

use crate::Error;
use crate::bed::Regions;

/// `--max-depth` default of `samtools mpileup`.
const MAX_DEPTH: usize = 8000;

/// `BAM_FREVERSE`, `BAM_FUNMAP`, `BAM_FPAIRED` and `BAM_FPROPER_PAIR`.
const FLAG_REVERSE: u16 = 0x10;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;

impl From<snv_hts::Error> for Error {
    fn from(e: snv_hts::Error) -> Self {
        match e {
            snv_hts::Error::Open { path, .. } | snv_hts::Error::Fasta { path, .. } => Error::Open { path },
            snv_hts::Error::Index { path, .. } => Error::Index { path },
            snv_hts::Error::Region { region } => Error::Region { region },
            snv_hts::Error::Parse { line } => Error::Parse { line },
            snv_hts::Error::Pileup(_) => Error::Pileup,
            snv_hts::Error::Decode(_) | snv_hts::Error::Contig { .. } => Error::Decode,
        }
    }
}

/// An indexed BAM/CRAM file with its header.
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

    /// `samtools mpileup FILE -r REGION` with `filter`: a pileup over the reads the index
    /// returns for `region`, keeping only columns inside it.
    ///
    /// # Errors
    ///
    /// Fails if the region cannot be parsed or queried.
    pub fn pileup(&mut self, region: &str, filter: ReadFilter) -> Result<Mpileup<'_>, Error> {
        let records = self.0.query(region)?;
        let (beg, end) = records.range();
        let (names, lens) = (records.target_names(), records.target_lens());
        Ok(Mpileup::new(Source::Region(Box::new(records)), names, lens, filter, Some((beg, end.unwrap_or(i64::MAX)))))
    }

    /// `samtools view FILE REGION -q MIN_MAPQ [-L BED]` without the text: the records of
    /// `region`, in file order, with a mapping quality of at least `min_mapq` and overlapping
    /// `bed` if given.
    ///
    /// # Errors
    ///
    /// Fails if the region cannot be parsed or queried.
    pub fn view(&mut self, region: &str, min_mapq: u8, bed: Option<Regions>) -> Result<View<'_>, Error> {
        let bed = bed.map(|b| b.with_target_names(self.0.target_names().to_vec()));
        Ok(View { records: self.0.query(region)?, min_mapq, bed, done: false })
    }

    /// The record as `samtools view` prints it (`sam_format1`), without the newline.
    ///
    /// # Errors
    ///
    /// Fails if the record cannot be formatted.
    pub fn format(&self, record: &Record) -> Result<String, Error> {
        Ok(snv_hts::format_sam(self.0.header(), &record.0)?)
    }

    /// The name of target `tid`.
    pub fn target_name(&self, tid: i32) -> Option<&str> {
        self.0.target_name(tid)
    }

    /// Target id of `contig` in the header, if present.
    pub fn tid(&self, contig: &str) -> Option<i32> {
        self.0.tid(contig)
    }

    /// `samtools mpileup -` over SAM text: a pileup over the records `lines` parse to with this
    /// file's header, as the original pipes `RealignReads` into `samtools mpileup -`. Every
    /// column is kept.
    ///
    /// # Errors
    ///
    /// Parse errors surface from [`Mpileup::next_column`]; this does not fail.
    pub fn pileup_sam_lines<'a, I>(&'a self, lines: I, filter: ReadFilter) -> Result<Mpileup<'a>, Error>
    where
        I: Iterator<Item = Result<String, Error>> + 'a,
    {
        let source = Source::Lines { lines: Box::new(lines), alignments: &self.0 };
        Ok(Mpileup::new(source, self.0.target_names(), self.0.target_lens(), filter, None))
    }
}

/// The records of a `samtools view` region query (see [`AlignmentFile::view`]).
#[derive(Debug)]
pub struct View<'file> {
    records: Records<'file>,
    min_mapq: u8,
    bed: Option<Regions>,
    done: bool,
}

impl Iterator for View<'_> {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.done {
            let read = match self.records.next_read() {
                Ok(Some(read)) => read,
                Ok(None) => {
                    self.done = true;
                    return None;
                }
                Err(e) => {
                    self.done = true;
                    return Some(Err(e.into()));
                }
            };
            if read.mapping_quality() < self.min_mapq {
                continue;
            }
            // `samtools view -L`: reads overlapping a BED interval.
            if let Some(bed) = &self.bed
                && !bed.overlaps(read.reference_sequence_id(), read.position(), snv_hts::end_pos(&read))
            {
                continue;
            }
            return Some(Ok(Record(read)));
        }
        None
    }
}

/// The read filters of `samtools mpileup` (`mplp_func`) with the options the original passes.
#[derive(Debug, Clone, Default)]
pub struct ReadFilter {
    /// `--excl-flags`.
    pub excluded_flags: u16,
    /// `--min-MQ`.
    pub min_mapq: u8,
    /// `-l BED`: the columns written must overlap it, and without `-a` so must the reads.
    pub bed: Option<Regions>,
    /// `-a`: whether positions without reads get columns too.
    pub positions: Positions,
}

impl ReadFilter {
    /// `mplp_func`: whether a read reaches the pileup.
    fn keeps(&self, read: &Read) -> bool {
        let flag = read.flags();
        if read.reference_sequence_id() < 0 || flag & FLAG_UNMAPPED != 0 || flag & self.excluded_flags != 0 {
            return false;
        }
        // With `-a`, `mplp_func` keeps reads outside the BED.
        if let Some(bed) = &self.bed
            && self.positions == Positions::Covered
            && !bed.overlaps(read.reference_sequence_id(), read.position(), snv_hts::end_pos(read))
        {
            return false;
        }
        if read.mapping_quality() < self.min_mapq {
            return false;
        }
        // Anomalous read pairs are skipped unless `-A`.
        !(flag & FLAG_PAIRED != 0 && flag & FLAG_PROPER_PAIR == 0)
    }
}

/// Which positions `samtools mpileup` writes a column for.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Positions {
    /// Positions with reads (the default).
    #[default]
    Covered,
    /// `-a`: also every position without reads between the region's start and the end of the
    /// region (or contig), once any read has a column in the region. These columns have no
    /// entries.
    All,
}

/// An owned alignment record (`bam1_t`).
#[derive(Debug, Clone)]
pub struct Record(Read);

impl Record {
    /// A read-only view of the record.
    pub fn view(&self) -> RecordView<'_> {
        RecordView(&self.0)
    }

    /// 0-based leftmost position.
    pub fn pos(&self) -> i64 {
        self.0.position()
    }

    /// One past the last reference position (`bam_endpos`).
    pub fn end_pos(&self) -> i64 {
        snv_hts::end_pos(&self.0)
    }

    /// Target id.
    pub fn tid(&self) -> i32 {
        self.0.reference_sequence_id()
    }

    /// Mapping quality.
    pub fn mapq(&self) -> u8 {
        self.0.mapping_quality()
    }
}

/// Where the pileup's reads come from.
enum Source<'a> {
    /// An index query (`samtools mpileup FILE -r REGION`).
    Region(Box<Records<'a>>),
    /// SAM lines (`samtools mpileup -`).
    Lines { lines: Box<dyn Iterator<Item = Result<String, Error>> + 'a>, alignments: &'a Alignments },
}

impl ReadSource for Source<'_> {
    fn next_read(&mut self) -> Result<Option<Read>, snv_hts::Error> {
        match self {
            Source::Region(records) => records.next_read(),
            Source::Lines { lines, alignments } => match lines.next() {
                None => Ok(None),
                // The line iterator's own errors (reading the realigner's output) are kept as
                // parse failures of a line, as htslib would report them.
                Some(Err(e)) => Err(snv_hts::Error::Parse { line: e.to_string() }),
                Some(Ok(line)) => snv_hts::parse_sam(alignments.header(), &line).map(Some),
            },
        }
    }
}

type Filter<'a> = Box<dyn FnMut(&mut Read) -> bool + 'a>;

/// A running pileup: `samtools mpileup`'s column loop.
pub struct Mpileup<'file> {
    pileup: snv_hts::Pileup<Source<'file>, Filter<'file>>,
    filter: Rc<ReadFilter>,
    names: &'file [String],
    lens: &'file [i64],
    /// `-r`: columns outside `[beg, end)` are skipped.
    region: Option<(i64, i64)>,
    /// With `-a`: the last column of the engine inside the region (`last_tid`, `last_pos`).
    last: Option<(i32, i64)>,
    /// With `-a`: positions without reads still to write, before `held`.
    gaps: VecDeque<Gap>,
    /// The engine's column that waits for the gaps before it: tid and position.
    held: Option<(i32, i64)>,
    /// The engine has no more columns.
    done: bool,
}

/// Positions `next..end` of target `tid`, written without reads under `-a`.
#[derive(Debug, Clone, Copy)]
struct Gap {
    tid: i32,
    next: i64,
    end: i64,
}

impl std::fmt::Debug for Mpileup<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mpileup").field("region", &self.region).finish_non_exhaustive()
    }
}

impl<'file> Mpileup<'file> {
    fn new(
        source: Source<'file>,
        names: &'file [String],
        lens: &'file [i64],
        mut filter: ReadFilter,
        region: Option<(i64, i64)>,
    ) -> Self {
        if let Some(bed) = filter.bed.take() {
            filter.bed = Some(bed.with_target_names(names.to_vec()));
        }
        let filter = Rc::new(filter);
        let keeps = Rc::clone(&filter);
        let engine_filter: Filter<'file> = Box::new(move |read| keeps.keeps(read));
        Self {
            pileup: snv_hts::Pileup::new(source, MAX_DEPTH, Overlaps::Adjust, engine_filter),
            filter,
            names,
            lens,
            region,
            last: None,
            gaps: VecDeque::new(),
            held: None,
            done: false,
        }
    }

    /// Advances to the next column `samtools mpileup` writes.
    ///
    /// # Errors
    ///
    /// Fails if a record cannot be read or parsed, or the input is not sorted.
    pub fn next_column(&mut self) -> Result<Option<Column<'_>>, Error> {
        loop {
            if let Some(gap) = self.gaps.front_mut() {
                let (tid, pos) = (gap.tid, gap.next);
                if pos >= gap.end {
                    self.gaps.pop_front();
                    continue;
                }
                gap.next += 1;
                if self.overlaps_bed(tid, pos) {
                    return Ok(Some(Column { contig: self.target_name(tid), pos, column: None }));
                }
                continue;
            }
            if let Some((tid, pos)) = self.held.take() {
                if !self.overlaps_bed(tid, pos) {
                    continue;
                }
                let contig = self.target_name(tid);
                return Ok(Some(Column { contig, pos, column: Some(self.pileup.current()) }));
            }
            if self.done {
                return Ok(None);
            }
            let Some(column) = self.pileup.next_column()? else {
                self.done = true;
                // The rest of the last contig with a column, up to the region's end.
                if self.filter.positions == Positions::All
                    && let Some((last_tid, last_pos)) = self.last
                {
                    let end = self.region.map_or(i64::MAX, |(_, end)| end).min(self.target_len(last_tid));
                    self.gaps.push_back(Gap { tid: last_tid, next: last_pos + 1, end });
                }
                continue;
            };
            let (tid, pos) = (column.tid, column.pos);
            if self.region.is_some_and(|(beg, end)| pos < beg || pos >= end) {
                continue;
            }
            if self.filter.positions == Positions::All {
                let beg = self.region.map_or(0, |(beg, _)| beg);
                let next = match self.last {
                    Some((last_tid, last_pos)) if last_tid == tid => last_pos + 1,
                    Some((last_tid, last_pos)) => {
                        // Without a region, the rest of the previous contig comes first.
                        if self.region.is_none() {
                            let end = self.target_len(last_tid);
                            self.gaps.push_back(Gap { tid: last_tid, next: last_pos + 1, end });
                        }
                        beg
                    }
                    None => beg,
                };
                self.gaps.push_back(Gap { tid, next: next.max(beg), end: pos });
                self.last = Some((tid, pos));
            }
            self.held = Some((tid, pos));
        }
    }

    /// Whether the `-l` BED (if any) overlaps 0-based position `pos` of target `tid`.
    fn overlaps_bed(&self, tid: i32, pos: i64) -> bool {
        self.filter.bed.as_ref().is_none_or(|bed| bed.overlaps(tid, pos, pos + 1))
    }

    /// The name of target `tid`.
    fn target_name(&self, tid: i32) -> &'file str {
        usize::try_from(tid).ok().and_then(|i| self.names.get(i)).map_or("", String::as_str)
    }

    /// The length of target `tid`.
    fn target_len(&self, tid: i32) -> i64 {
        usize::try_from(tid).ok().and_then(|i| self.lens.get(i)).copied().unwrap_or(0)
    }
}

/// One pileup column: every filtered read overlapping reference position `pos`.
pub struct Column<'plp> {
    /// Contig name.
    pub contig: &'plp str,
    /// 0-based reference position.
    pub pos: i64,
    /// The engine's column; `None` for a position without reads (`-a`).
    column: Option<snv_hts::Column<'plp>>,
}

impl std::fmt::Debug for Column<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Column").field("pos", &self.pos).field("depth", &self.entries().len()).finish()
    }
}

impl<'plp> Column<'plp> {
    /// The reads in this column, in htslib's order.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = Entry<'plp>> + 'plp {
        let entries = self.column.map_or(&[][..], |c| c.entries());
        let column = self.column;
        entries
            .iter()
            .map(move |e| Entry { entry: e, read: column.expect("entries come from an engine column").read(e) })
    }
}

/// One read in a pileup column.
#[derive(Clone, Copy)]
pub struct Entry<'plp> {
    entry: &'plp snv_hts::Entry,
    read: &'plp Read,
}

impl std::fmt::Debug for Entry<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Entry").field("qpos", &self.qpos()).field("indel", &self.indel()).finish()
    }
}

impl<'plp> Entry<'plp> {
    /// The record.
    pub fn record(&self) -> RecordView<'plp> {
        RecordView(self.read)
    }

    /// This column is the read's first base.
    pub fn is_head(&self) -> bool {
        self.entry.is_head
    }

    /// This column is the read's last base.
    pub fn is_tail(&self) -> bool {
        self.entry.is_tail
    }

    /// The read has a deletion (or a reference skip) here.
    pub fn is_del(&self) -> bool {
        self.entry.is_del
    }

    /// The read skips this position with an `N` CIGAR operation.
    pub fn is_refskip(&self) -> bool {
        self.entry.is_refskip
    }

    /// Indel after this position: positive for insertions, negative for deletions.
    pub fn indel(&self) -> i32 {
        self.entry.indel
    }

    /// Query position of this column.
    pub fn qpos(&self) -> i32 {
        i32::try_from(self.entry.qpos).expect("query position fits i32")
    }

    /// Index of the CIGAR operation holding this column.
    pub fn cigar_index(&self) -> usize {
        usize::try_from(self.entry.cigar_ind).expect("CIGAR index is not negative")
    }
}

/// A borrowed alignment record.
#[derive(Clone, Copy)]
pub struct RecordView<'a>(&'a Read);

impl std::fmt::Debug for RecordView<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RecordView").field("qname", &self.qname()).finish()
    }
}

impl<'a> RecordView<'a> {
    /// Target id.
    pub fn tid(&self) -> i32 {
        self.0.reference_sequence_id()
    }

    /// 0-based leftmost position.
    pub fn pos(&self) -> i64 {
        self.0.position()
    }

    /// SAM flag.
    pub fn flag(&self) -> u16 {
        self.0.flags()
    }

    /// Mapping quality.
    pub fn mapq(&self) -> u8 {
        self.0.mapping_quality()
    }

    /// The read maps to the reverse strand.
    pub fn is_reverse(&self) -> bool {
        self.flag() & FLAG_REVERSE != 0
    }

    /// Number of query bases.
    pub fn len(&self) -> i32 {
        i32::try_from(self.0.sequence_len()).expect("read length fits i32")
    }

    /// Whether the read has no bases.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The read name.
    pub fn qname(&self) -> &'a str {
        std::str::from_utf8(self.0.name()).unwrap_or("")
    }

    /// The CIGAR operations as `(op, length)` with htslib's op codes.
    pub fn cigar(&self) -> impl ExactSizeIterator<Item = (u32, u32)> + 'a {
        self.0.cigar().iter().map(|&c| (c & 0xf, c >> 4))
    }

    /// 4-bit `seq_nt16` code of query base `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is outside the read.
    pub fn base_code(&self, i: i32) -> u8 {
        self.0.base(self.index(i))
    }

    /// Base quality of query base `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is outside the read.
    pub fn qual(&self, i: i32) -> u8 {
        self.0.quality_scores()[self.index(i)]
    }

    fn index(&self, i: i32) -> usize {
        let n = self.0.sequence_len();
        match usize::try_from(i) {
            Ok(i) if i < n => i,
            _ => panic!("query position {i} outside read of {n} bases"),
        }
    }

    /// The value of aux tag `tag` if it is an integer (`[cCsSiI]`), without formatting it.
    pub fn aux_integer(&self, tag: &CStr) -> Option<i64> {
        let (ty, _) = self.0.tag(two(tag)?)?;
        if !b"IiCcSs".contains(&ty) {
            return None;
        }
        self.0.integer_tag(two(tag)?)
    }

    /// The value of aux tag `tag` as `samtools mpileup --output-extra` prints it, or `None`
    /// when the read does not have it.
    pub fn aux_text(&self, tag: &CStr) -> Option<String> {
        let tag = two(tag)?;
        let (ty, value) = self.0.tag(tag)?;
        Some(match ty {
            b'Z' | b'H' => String::from_utf8_lossy(self.0.string_tag(tag).unwrap_or_default()).into_owned(),
            b'I' | b'i' | b'C' | b'c' | b'S' | b's' => self.0.integer_tag(tag).unwrap_or_default().to_string(),
            b'A' => char::from(value[0]).to_string(),
            // Floats (`kputd`) and arrays are not ported; the original only asks for `HP`.
            _ => "*".to_owned(),
        })
    }
}

/// A two-letter tag name.
fn two(tag: &CStr) -> Option<[u8; 2]> {
    tag.to_bytes().try_into().ok()
}

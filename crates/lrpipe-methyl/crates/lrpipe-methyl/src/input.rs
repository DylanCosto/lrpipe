//! Reads over a region from an indexed BAM or CRAM, as [`ModRead`]s.

use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use pasta::bam;
use pasta::sam::{
    self,
    alignment::record::{
        Flags,
        cigar::op::Kind,
        data::field::{Tag, Value, value::Array},
    },
};

use crate::mods::{Bases, decode};
use crate::pileup::RefWindow;
use crate::read::{ModRead, Op};

pub use lrpipe_reads::haplotags::Haplotags;
pub use lrpipe_reads::is_cram;
use lrpipe_reads::{Decoding, Fields, IndexedReads, MdNm, OpenError, Tags};

/// Why reads could not be read.
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    /// A file could not be opened or read.
    #[error("reading {}", path.display())]
    Io {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// CRAM was given without a reference to decode it.
    #[error("{} is CRAM, which needs a reference", .0.display())]
    CramWithoutReference(PathBuf),
}

impl From<OpenError> for InputError {
    fn from(e: OpenError) -> Self {
        match e {
            OpenError::Io { path, source } => InputError::Io { path, source },
            OpenError::CramWithoutReference(path) => InputError::CramWithoutReference(path),
        }
    }
}

/// Which alignments count, and which of their calls are decoded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadFilter {
    /// Alignments with any of these flags are skipped.
    pub skip: Flags,
    /// Alignments with a lower mapping quality are skipped.
    pub min_mapq: u8,
    /// Canonical bases whose calls are decoded.
    pub bases: Bases,
}

impl ReadFilter {
    /// Primary alignments only, any mapping quality, every call.
    pub const PRIMARY: ReadFilter = ReadFilter {
        skip: Flags::UNMAPPED.union(Flags::SECONDARY).union(Flags::SUPPLEMENTARY),
        min_mapq: 0,
        bases: Bases::ALL,
    };

    fn passes(self, flags: Flags, mapq: u8) -> bool {
        !flags.intersects(self.skip) && mapq >= self.min_mapq
    }
}

impl Default for ReadFilter {
    fn default() -> Self {
        ReadFilter::PRIMARY
    }
}

/// Reference sequences for decoding CRAM, loaded a contig at a time and shared by every
/// reader opened with it.
#[derive(Clone, Debug)]
pub struct SharedReference(lrpipe_reads::Reference);

impl SharedReference {
    /// Opens a FASTA with a `.fai` index.
    ///
    /// # Errors
    ///
    /// If the FASTA or its index cannot be read.
    pub fn open(path: &Path) -> Result<Self, InputError> {
        Ok(SharedReference(lrpipe_reads::Reference::open(path)?))
    }

    /// Drops every contig read so far (they are kept until then, and a genome's worth would
    /// stay in memory); a contig needed again is read again.
    pub fn clear(&self) {
        self.0.clear();
    }

    /// `contig`'s length by the `.fai` index; `None` when the reference has no such contig.
    #[must_use]
    pub fn length(&self, contig: &str) -> Option<u64> {
        self.0.length(contig)
    }

    /// Bases `[start, end)` of `contig`, uppercase, clipped to the contig; `None` when the
    /// reference has no such contig.
    ///
    /// # Errors
    ///
    /// If the sequence cannot be read.
    pub fn window(&self, contig: &str, start: i64, end: i64) -> io::Result<Option<RefWindow>> {
        let Some(seq) = self.0.sequence(contig)? else { return Ok(None) };
        let bases: &[u8] = seq.as_ref().as_ref();
        let s = start.clamp(0, bases.len() as i64) as usize;
        let e = end.clamp(0, bases.len() as i64) as usize;
        Ok(Some(RefWindow { start: s as i64, seq: bases[s..e.max(s)].to_ascii_uppercase() }))
    }
}

/// The tags methylation reads: the modification calls (`MM`, `ML`, `MN`, and the older `Mm` and
/// `Ml`) and the haplotype (`HP`).
pub const TAGS: [Tag; 6] = [
    Tag::BASE_MODIFICATIONS,
    Tag::BASE_MODIFICATION_PROBABILITIES,
    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
    Tag::new(b'M', b'm'),
    Tag::new(b'M', b'l'),
    Tag::new(b'H', b'P'),
];

/// An open, indexed alignment file.
#[derive(Debug)]
pub struct Reads {
    reads: IndexedReads,
    /// The haplotypes the alignments are given in place of their `HP` tags.
    haplotags: Option<Arc<Haplotags>>,
}

impl Reads {
    /// Opens `path` (BAM with `.bai`/`.csi`, CRAM with `.crai`, decoded with `reference`).
    ///
    /// # Errors
    ///
    /// If the file or its index cannot be opened, or CRAM is given without a reference.
    pub fn open(path: &Path, reference: Option<&SharedReference>) -> Result<Self, InputError> {
        let mut reads = IndexedReads::open(path, reference.map(|r| &r.0))?;
        reads.set_decoding(Decoding {
            // Nothing here reads MD or NM.
            md_nm: MdNm::AsStored,
            // Nor base qualities or read groups.
            fields: Fields::all() - Fields::QUALITY_SCORES - Fields::READ_GROUP,
            tags: Tags::Only(TAGS.to_vec()),
        });
        Ok(Reads { reads, haplotags: None })
    }

    /// Gives the alignments the haplotypes `haplotags` assigns their reads, in place of their
    /// `HP` tags, as `haplotag` would have tagged them: a read's primary alignment its
    /// haplotype, any other none.
    pub fn set_haplotags(&mut self, haplotags: Arc<Haplotags>) {
        self.haplotags = Some(haplotags);
    }

    /// The file's header.
    #[must_use]
    pub fn header(&self) -> &sam::Header {
        self.reads.header()
    }

    /// The compressed bytes the index says a fetch of `[start, end)` of `contig` reads: BAM
    /// chunks' spans, or CRAM slices' lengths. A measure of how much data the range holds,
    /// which grows with depth; 0 for a contig the reads do not name.
    ///
    /// # Errors
    ///
    /// If the BAM index cannot be queried.
    pub fn indexed_bytes(&self, contig: &str, start: i64, end: i64) -> Result<u64, InputError> {
        self.reads.indexed_bytes(contig, start, end).map_err(|source| self.error(source))
    }

    /// Appends to `out` the alignments overlapping `[start, end)` of `contig` that pass
    /// `filter` and carry modification calls, in file order. Alignments without `MM` tags, and
    /// those whose tags do not decode or whose hard clips make `MM` positions unplaceable, are
    /// left out and counted.
    ///
    /// # Errors
    ///
    /// If records cannot be read.
    pub fn fetch(
        &mut self,
        contig: &str,
        start: i64,
        end: i64,
        filter: ReadFilter,
        out: &mut Vec<ModRead>,
    ) -> Result<Fetched, InputError> {
        self.fetch_each(contig, start, end, filter, |read| out.push(read))
    }

    /// Calls `f` with each alignment [`fetch`](Self::fetch) would append, as it is read,
    /// so the alignments of a range need not be in memory together.
    ///
    /// # Errors
    ///
    /// If records cannot be read.
    pub fn fetch_each(
        &mut self,
        contig: &str,
        start: i64,
        end: i64,
        filter: ReadFilter,
        mut f: impl FnMut(ModRead),
    ) -> Result<Fetched, InputError> {
        let mut fetched = Fetched::default();
        let list = self.haplotags.as_deref();
        let listed = list.and_then(|h| h.contig(contig));
        let result = self.reads.query(contig, start, end, |record| {
            match convert(record, filter)? {
                Converted::Read(mut r) => {
                    fetched.kept += 1;
                    if let Some(list) = list {
                        let name = record.name().map_or(&b""[..], |n| n.as_ref());
                        r.haplotype = list
                            .of_contig_alignment(listed, name, r.flags)
                            .and_then(|a| matches!(a.haplotype, 1 | 2).then_some(a.haplotype));
                    }
                    f(r);
                }
                Converted::Untagged => fetched.untagged += 1,
                Converted::Unplaceable => fetched.unplaceable += 1,
                Converted::Filtered => {}
            }
            Ok(())
        });
        result.map_err(|source| self.error(source))?;
        Ok(fetched)
    }

    fn error(&self, source: io::Error) -> InputError {
        InputError::Io { path: self.reads.path().into(), source }
    }
}

/// What a fetch took and left out. An alignment is counted again by every window it reaches.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Fetched {
    /// Alignments with modification calls.
    pub kept: usize,
    /// Alignments the filter let through without an `MM` tag.
    pub untagged: usize,
    /// Alignments whose `MM`/`ML` tags do not fit their sequence.
    pub unplaceable: usize,
}

impl std::ops::AddAssign for Fetched {
    fn add_assign(&mut self, other: Self) {
        self.kept += other.kept;
        self.untagged += other.untagged;
        self.unplaceable += other.unplaceable;
    }
}

/// What became of one record.
pub(crate) enum Converted {
    Read(ModRead),
    /// The filter dropped it, or it is unmapped.
    Filtered,
    Untagged,
    Unplaceable,
}

/// A record as a [`ModRead`], when it passes `filter` and its tags fit its sequence.
pub(crate) fn convert(r: &bam::Record, filter: ReadFilter) -> io::Result<Converted> {
    let flags = r.flags();
    let mapq = r.mapping_quality().map_or(255, u8::from);
    if !filter.passes(flags, mapq) {
        return Ok(Converted::Filtered);
    }
    let Some(start) = r.alignment_start().transpose()? else { return Ok(Converted::Filtered) };
    let data = r.data();
    let tag = |a: Tag, legacy: [u8; 2]| -> io::Result<Option<Value<'_>>> {
        Ok(match data.get(&a).transpose()? {
            Some(v) => Some(v),
            None => data.get(&Tag::new(legacy[0], legacy[1])).transpose()?,
        })
    };
    let Some(Value::String(mm)) = tag(Tag::BASE_MODIFICATIONS, *b"Mm")? else { return Ok(Converted::Untagged) };
    let mm = mm.to_vec();
    let ml: Vec<u8> = match tag(Tag::BASE_MODIFICATION_PROBABILITIES, *b"Ml")? {
        Some(Value::Array(Array::UInt8(v))) => v.iter().collect::<io::Result<_>>()?,
        _ => Vec::new(),
    };
    let mn = match data.get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH).transpose()? {
        Some(v) => v.as_int(),
        None => None,
    };
    let haplotype = match data.get(&Tag::new(b'H', b'P')).transpose()? {
        Some(v) => v.as_int().and_then(|h| u8::try_from(h).ok()).filter(|h| matches!(h, 1 | 2)),
        None => None,
    };
    let mut cigar = Vec::new();
    let mut hard_clipped = false;
    for op in r.cigar().iter() {
        let op = op?;
        let len = op.len() as u32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => cigar.push((Op::Match, len)),
            Kind::Insertion | Kind::SoftClip => cigar.push((Op::Read, len)),
            Kind::Deletion => cigar.push((Op::Deletion, len)),
            Kind::Skip => cigar.push((Op::Skip, len)),
            Kind::HardClip => hard_clipped = true,
            Kind::Pad => {}
        }
    }
    // BAM's base alphabet (`=ACMGRSVTWYHKDBN`) is already uppercase.
    let seq: Vec<u8> = r.sequence().iter().collect();
    // MM counts bases of the whole read; with hard clips SEQ is only part of it, which MN
    // (the length MM was written for) must then confirm.
    if mn.is_some_and(|n| n != seq.len() as i64) || (hard_clipped && mn.is_none()) {
        return Ok(Converted::Unplaceable);
    }
    let reverse = flags.is_reverse_complemented();
    let Ok(groups) = decode(&mm, &ml, &seq, reverse, filter.bases) else { return Ok(Converted::Unplaceable) };
    let start = usize::from(start) as i64 - 1;
    let mut read = ModRead::new(start, reverse, cigar, seq, haplotype, groups);
    read.flags = flags.bits();
    read.mapq = mapq;
    Ok(Converted::Read(read))
}

//! Reads over a region from an indexed BAM or CRAM: alignment position, CIGAR, bases and
//! haplotag.

use std::io;
use std::path::Path;

use pasta::bam;
use pasta::sam::alignment::record::{
    Flags,
    cigar::op::Kind,
    data::field::{Tag, Value},
};

pub use lrpipe_reads::is_cram;
use std::sync::Arc;

use lrpipe_reads::haplotags::Haplotags;
use lrpipe_reads::{Decoding, Fields, IndexedReads, MdNm, OpenError, Tags};

use crate::Error;

/// CIGAR operations, reduced to what extraction needs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// Aligned bases (M, =, X).
    Match,
    /// Inserted read bases (I).
    Ins,
    /// Deleted or skipped reference bases (D, N).
    Del,
    /// Soft-clipped read bases (S).
    SoftClip,
}

/// One alignment.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Read {
    /// 0-based start on the reference.
    pub start: i64,
    /// Exclusive end on the reference.
    pub end: i64,
    /// Operations and their lengths; hard clips and padding left out.
    pub cigar: Vec<(Op, u32)>,
    /// Uppercase bases of the whole read, soft clips included.
    pub seq: Vec<u8>,
    /// The `HP` tag (1 or 2) of a haplotagged read.
    pub haplotype: Option<u8>,
    /// The `PS` tag: the phase block the haplotype belongs to. Haplotypes of different phase
    /// sets are unrelated.
    pub phase_set: Option<i64>,
    /// A hash of the read's name ([`name_hash`]), by which two alignments of one read are told
    /// to be one molecule; `None` when the name is unknown.
    pub name: Option<u64>,
    /// Whether this is a supplementary alignment.
    pub supplementary: bool,
    /// Mapping quality (255 when unavailable).
    pub mapq: u8,
}

/// A read name's hash (FNV-1a), for [`Read::name`].
#[must_use]
pub fn name_hash(name: &[u8]) -> u64 {
    name.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3))
}

/// An open, indexed alignment file.
#[derive(Debug)]
pub struct Reads {
    reads: IndexedReads,
    /// The haplotypes the alignments are given in place of their `HP` tags.
    haplotags: Option<Arc<Haplotags>>,
}

/// Records with any of these flags are skipped: unmapped, secondary, QC fail, duplicate.
const SKIP: Flags = Flags::UNMAPPED.union(Flags::SECONDARY).union(Flags::QC_FAIL).union(Flags::DUPLICATE);

/// Reference sequences for decoding CRAM, loaded a contig at a time on first use and shared
/// by every reader opened with it (cloning shares the cache), so threads hold one copy.
#[derive(Clone, Debug)]
pub struct SharedReference(lrpipe_reads::Reference);

impl SharedReference {
    /// Opens a FASTA with a `.fai` index.
    ///
    /// # Errors
    ///
    /// If the FASTA or its index cannot be read.
    pub fn open(path: &Path) -> Result<Self, Error> {
        Ok(SharedReference(lrpipe_reads::Reference::open(path).map_err(open_error)?))
    }
}

fn open_error(e: OpenError) -> Error {
    match e {
        OpenError::Io { path, source } => Error::Io { path, source },
        OpenError::CramWithoutReference(path) => Error::Input(format!("{}: CRAM needs a reference", path.display())),
    }
}

impl Reads {
    /// Opens `path` (BAM with `.bai`/`.csi`, CRAM with `.crai`); CRAM is decoded with
    /// `reference`.
    ///
    /// # Errors
    ///
    /// If the file or its index cannot be opened, or CRAM is given without a reference.
    pub fn open(path: &Path, reference: Option<&SharedReference>) -> Result<Self, Error> {
        let mut reads = IndexedReads::open(path, reference.map(|r| &r.0)).map_err(open_error)?;
        reads.set_decoding(Decoding {
            // Nothing here reads MD or NM, base qualities, read groups or tags but the haplotype and its phase set.
            md_nm: MdNm::AsStored,
            fields: Fields::all() - Fields::QUALITY_SCORES - Fields::READ_GROUP,
            tags: Tags::Only(vec![Tag::new(b'H', b'P'), Tag::new(b'P', b'S')]),
        });
        Ok(Reads { reads, haplotags: None })
    }

    /// Gives the alignments the haplotypes `haplotags` assigns their reads, in place of their
    /// `HP` tags, as `haplotag` would have tagged them: a read's primary alignment its
    /// haplotype, any other none.
    pub fn set_haplotags(&mut self, haplotags: Arc<Haplotags>) {
        self.haplotags = Some(haplotags);
    }

    /// Appends to `out` the alignments overlapping `[start, end)` of `contig` with MAPQ of at
    /// least `min_mapq`, in file order. A contig the file does not have gives none.
    ///
    /// # Errors
    ///
    /// If records cannot be read.
    pub fn fetch(
        &mut self,
        contig: &str,
        start: i64,
        end: i64,
        min_mapq: u8,
        out: &mut Vec<Read>,
    ) -> Result<(), Error> {
        let err = |e: io::Error| Error::Input(format!("reading {contig}:{start}-{end}: {e}"));
        let list = self.haplotags.as_deref();
        let listed = list.and_then(|h| h.contig(contig));
        self.reads
            .query(contig, start, end, |record| {
                let Some(mut read) = convert(record, min_mapq)? else { return Ok(()) };
                if let Some(list) = list {
                    let name = record.name().map_or(&b""[..], |n| n.as_ref());
                    let assigned = list.of_contig_alignment(listed, name, record.flags().bits());
                    read.haplotype = assigned.map(|a| a.haplotype);
                    read.phase_set = assigned.map(|a| a.phase_set);
                }
                out.push(read);
                Ok(())
            })
            .map_err(err)
    }
}

/// An alignment record's fields, as a program that decodes records itself has them.
#[derive(Clone, Copy, Debug)]
pub struct RecordFields<'r> {
    /// SAM flags.
    pub flags: u16,
    /// 0-based alignment start (`None` when unplaced).
    pub start: Option<i64>,
    /// Mapping quality (255 when unavailable).
    pub mapq: u8,
    /// CIGAR operations as BAM codes (0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X) and lengths.
    pub cigar: &'r [(u8, u32)],
    /// The `HP` tag, when an integer.
    pub haplotype: Option<u8>,
    /// The `PS` tag, when an integer.
    pub phase_set: Option<i64>,
    /// The read's name (empty when unknown).
    pub name: &'r [u8],
}

impl Read {
    /// A read from a record's fields and its bases (any case), or `None` when lrpipe-tr does not
    /// use it: unmapped, secondary, QC fail or duplicate, MAPQ below `min_mapq`, or no bases.
    pub fn from_record(rec: RecordFields<'_>, seq: impl IntoIterator<Item = u8>, min_mapq: u8) -> Option<Read> {
        if !passes(rec.flags, rec.mapq, min_mapq) {
            return None;
        }
        let start = rec.start?;
        let mut cigar = Vec::with_capacity(rec.cigar.len());
        let mut span = 0i64;
        for &(code, len) in rec.cigar {
            let op = match code {
                0 | 7 | 8 => Op::Match,
                1 => Op::Ins,
                2 | 3 => Op::Del,
                4 => Op::SoftClip,
                _ => continue,
            };
            if matches!(op, Op::Match | Op::Del) {
                span += i64::from(len);
            }
            cigar.push((op, len));
        }
        let seq: Vec<u8> = seq.into_iter().map(|b| b.to_ascii_uppercase()).collect();
        if seq.is_empty() {
            return None;
        }
        Some(Read {
            start,
            end: start + span,
            cigar,
            seq,
            haplotype: rec.haplotype,
            phase_set: rec.phase_set,
            name: (!rec.name.is_empty()).then(|| name_hash(rec.name)),
            supplementary: Flags::from_bits_retain(rec.flags).contains(Flags::SUPPLEMENTARY),
            mapq: rec.mapq,
        })
    }
}

/// Whether a record with these flags and mapping quality is used.
pub fn passes(flags: u16, mapq: u8, min_mapq: u8) -> bool {
    !Flags::from_bits_retain(flags).intersects(SKIP) && mapq >= min_mapq
}

fn convert(r: &bam::Record, min_mapq: u8) -> io::Result<Option<Read>> {
    let flags = r.flags().bits();
    let mapq = r.mapping_quality().map_or(255, u8::from);
    if !passes(flags, mapq, min_mapq) {
        return Ok(None);
    }
    let start = r.alignment_start().transpose()?.map(|p| usize::from(p) as i64 - 1);
    let mut cigar = Vec::new();
    for op in r.cigar().iter() {
        let op = op?;
        let code = match op.kind() {
            Kind::Match => 0,
            Kind::Insertion => 1,
            Kind::Deletion => 2,
            Kind::Skip => 3,
            Kind::SoftClip => 4,
            Kind::HardClip => 5,
            Kind::Pad => 6,
            Kind::SequenceMatch => 7,
            Kind::SequenceMismatch => 8,
        };
        cigar.push((code, op.len() as u32));
    }
    let haplotype = match r.data().get(&Tag::new(b'H', b'P')).transpose()? {
        Some(v @ (Value::Int8(_) | Value::UInt8(_) | Value::Int16(_) | Value::UInt16(_) | Value::Int32(_))) => {
            v.as_int().and_then(|h| u8::try_from(h).ok())
        }
        _ => None,
    };
    let phase_set = match r.data().get(&Tag::new(b'P', b'S')).transpose()? {
        Some(
            v @ (Value::Int8(_)
            | Value::UInt8(_)
            | Value::Int16(_)
            | Value::UInt16(_)
            | Value::Int32(_)
            | Value::UInt32(_)),
        ) => v.as_int(),
        _ => None,
    };
    let name = r.name().map_or(&b""[..], |n| n.as_ref());
    let fields = RecordFields { flags, start, mapq, cigar: &cigar, haplotype, phase_set, name };
    Ok(Read::from_record(fields, r.sequence().iter(), min_mapq))
}

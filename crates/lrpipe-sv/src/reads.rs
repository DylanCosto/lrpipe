//! Reading alignments and extracting SV evidence ("leads"), as in Sniffles2's `leadprov.py`.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};
use pasta::bam;
use pasta::core::{Position, Region};
use pasta::cram;
use pasta::csi::{self, BinningIndex};
use pasta::fasta;
use pasta::sam;

use lrpipe_reads::haplotags::ContigHaplotags;
use lrpipe_reads::{CramSource, Decoding, Fields, MdNm};
use pasta::sam::alignment::record::data::field::Tag;

use crate::config::{Config, Haplotypes};
use crate::read_stats::ReadStats;

pub const CINS: u8 = 1;
pub const CDEL: u8 = 2;
pub const CSOFT_CLIP: u8 = 4;
pub const CHARD_CLIP: u8 = 5;

const FLAG_REVERSE: u16 = 0x10;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// SV types that produce VCF records, in Sniffles2's processing order (`sv.TYPES`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SvType {
    Ins,
    Del,
    Dup,
    Inv,
    Bnd,
}

impl SvType {
    pub const ALL: [SvType; 5] = [SvType::Ins, SvType::Del, SvType::Dup, SvType::Inv, SvType::Bnd];

    pub fn as_str(self) -> &'static str {
        match self {
            SvType::Ins => "INS",
            SvType::Del => "DEL",
            SvType::Dup => "DUP",
            SvType::Inv => "INV",
            SvType::Bnd => "BND",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Inline,
    BndSa,
    SplitPrim,
    SplitSup,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BndInfo {
    pub mate_contig: Arc<str>,
    pub mate_ref_start: i64,
    /// True for N..., false for ...N
    pub is_first: bool,
    /// True for ]...], false for [...[
    pub is_reverse: bool,
}

/// One piece of evidence for an SV from one alignment (Sniffles2's `Lead`).
#[derive(Debug, Clone)]
pub struct Lead {
    pub read_id: u64,
    pub qname: Arc<[u8]>,
    pub ref_start: i64,
    pub ref_end: i64,
    pub qry_start: i64,
    pub qry_end: i64,
    pub reverse: bool,
    pub mapq: u8,
    /// Length-adjusted mismatch rate of the read (-1 if unknown).
    pub nm: f64,
    pub source: Source,
    /// `None` for soft-clip based (long) insertions.
    pub svlen: Option<i64>,
    pub seq: Option<Vec<u8>>,
    pub bnd: Option<Box<BndInfo>>,
    /// HP tag (0 if absent).
    pub hap: u8,
    /// PS tag as Python's `str()` of it, `"NULL"` if absent; `None` for BND leads.
    pub phase_set: Option<Arc<str>>,
    pub is_sa: bool,
    pub read_len: i64,
}

/// Evidence and coverage collected from one stretch of a contig.
#[derive(Debug, Default)]
pub struct ChunkLeads {
    /// Leads per [`SvType`], with their leadtab bin, in read order.
    pub leads: [Vec<(i64, Lead)>; 5],
    /// SINGLE_LEFT and SINGLE_RIGHT break leads (soft clips), with their bin; only for `--snf`.
    pub single: [Vec<(i64, Lead)>; 2],
    /// Reference span of each read used, for coverage.
    pub spans: Vec<(i64, i64)>,
    /// (first bin, end bin, HP) per read, for the reference haplotype counts.
    pub hap_bins: Vec<(i64, i64, u8)>,
    /// With `--phasing-evidence`, beside each of `hap_bins`: the read's name when the alignment
    /// is its primary one (which a haplotag list tags), else `None`.
    pub hap_reads: Vec<Option<Arc<[u8]>>>,
    /// Length-adjusted mismatch rate of each read with an NM tag, in read order.
    pub nms: Vec<f64>,
    pub read_count: u64,
    /// `--depth`: depth of the reads that start in the chunk's depth range.
    pub depth: Option<lrpipe_depth::ContigCounter>,
    /// `--repeats`: observations of the contig's repeat loci by the reads starting in the
    /// chunk, in read order.
    pub repeats: Vec<(u32, lrpipe_tr::genotype::Observation)>,
    /// `--read-stats`: the records starting in the chunk's statistics range.
    pub read_stats: Option<ReadStats>,
}

/// `--repeats`: the contig's loci and reference, for genotyping repeats from the same reads.
#[derive(Clone, Copy, Debug)]
pub struct RepeatPiece<'a> {
    pub loci: &'a [lrpipe_tr::catalog::Locus],
    pub reference: &'a lrpipe_tr::extract::RefWindow,
    pub opts: &'a lrpipe_tr::run::RunOptions,
}

/// Where one contig's reads take their haplotypes from ([`Haplotypes`]).
#[derive(Clone, Copy, Debug)]
enum ContigHaplotypes<'a> {
    /// The HP and PS tags.
    Tags,
    /// A haplotag list, and its reads of the contig if it has any.
    List(&'a lrpipe_reads::haplotags::Haplotags, Option<ContigHaplotags<'a>>),
    /// Nowhere (`--phasing-evidence`): the SVs are phased later.
    Later,
}

impl ContigHaplotypes<'_> {
    /// What a haplotag list assigns an alignment: its read's haplotype and phase set for a
    /// primary alignment (`flags` without the secondary and supplementary bits), none for
    /// another, or none for a read the list leaves out. Counts the primary ones looked up.
    fn listed(
        list: &lrpipe_reads::haplotags::Haplotags,
        listed: Option<ContigHaplotags<'_>>,
        name: &[u8],
        flags: u16,
    ) -> Option<(i64, String)> {
        let a = list.of_contig_alignment(listed, name, flags)?;
        Some((i64::from(a.haplotype), a.phase_set.to_string()))
    }
}

/// Offers a record starting at `pos` to the repeat collector: decoded into a lrpipe-tr read only
/// when lrpipe-tr uses it and it spans a locus. A haplotag list (`--haplotags`) gives the read's
/// haplotype in place of its HP tag.
fn offer_repeat(
    collector: &mut lrpipe_tr::stream::Collector<'_>,
    view: &RecordView<'_>,
    cigar: &[(u8, i64)],
    pos: i64,
    min_mapq: u8,
    haplotypes: ContigHaplotypes<'_>,
) -> Result<()> {
    if !lrpipe_tr::reads::passes(view.flags, view.mapq, min_mapq) {
        return Ok(());
    }
    let span: i64 = cigar.iter().filter(|(op, _)| matches!(op, 0 | 2 | 3 | 7 | 8)).map(|&(_, l)| l).sum();
    // Soft clips (op 4) at either end: an expansion the aligner did not align through.
    let clip =
        |op: Option<&(u8, i64)>| op.filter(|(o, _)| *o == 4).map_or(0, |&(_, l)| usize::try_from(l).unwrap_or(0));
    let clips = [clip(cigar.first()), clip(cigar.last())];
    if !collector.spans_any(pos, pos + span, clips) {
        return Ok(());
    }
    let ops: Vec<(u8, u32)> = cigar.iter().map(|&(op, len)| (op, len as u32)).collect();
    let (haplotype, phase_set) = match haplotypes {
        // Not counted: the SV leads count each alignment.
        ContigHaplotypes::List(_, listed) => listed
            .filter(|_| view.flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) == 0)
            .and_then(|l| l.get(view.name))
            .map_or((None, None), |a| (Some(i64::from(a.haplotype)), Some(a.phase_set))),
        ContigHaplotypes::Tags | ContigHaplotypes::Later => {
            let tags = parse_tags(view.tags)?;
            (tags.hp, tags.ps.as_deref().and_then(|ps| ps.parse().ok()))
        }
    };
    let haplotype = haplotype.and_then(|h| u8::try_from(h).ok());
    let fields = lrpipe_tr::reads::RecordFields {
        flags: view.flags,
        start: Some(pos),
        mapq: view.mapq,
        cigar: &ops,
        haplotype,
        phase_set,
        name: view.name,
    };
    if let Some(read) = lrpipe_tr::reads::Read::from_record(fields, view.seq.slice(0, view.seq.len() as i64), min_mapq)
    {
        collector.push(&read);
    }
    Ok(())
}

/// An open, indexed alignment file: BAM, or CRAM decoded against the reference.
pub struct Alignments {
    path: String,
    format: Format,
    pub contigs: Vec<(String, i64)>,
}

enum Format {
    Bam { header: sam::Header, index: Box<dyn BinningIndex + Send + Sync> },
    Cram { header: sam::Header, index: cram::crai::Index, reference: fasta::Repository, decoding: Decoding },
}

/// The tags the read pass reads besides methylation's: the SV caller's `NM`, `SA`, `HP` and
/// `PS` (repeats read `HP` too).
const PASS_TAGS: [Tag; 4] = [Tag::EDIT_DISTANCE, Tag::OTHER_ALIGNMENTS, Tag::new(b'H', b'P'), Tag::new(b'P', b'S')];

/// What CRAM decoding gives the read pass: every field but the base qualities and the read
/// group, which none of its parts reads (`--read-stats` takes the qualities from BAM alone), the
/// `NM` tags htslib adds (Sniffles2 reads `NM`), and only the tags it reads.
fn pass_decoding() -> Decoding {
    Decoding {
        md_nm: MdNm::Regenerate,
        fields: Fields::all() - Fields::QUALITY_SCORES - Fields::READ_GROUP,
        tags: lrpipe_reads::Tags::Only(PASS_TAGS.to_vec()),
    }
}

impl std::fmt::Debug for Alignments {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Alignments").field("path", &self.path).finish_non_exhaustive()
    }
}

const INDEX_HINT: &str = "Please verify that your input file is sorted + indexed and that the index file is valid \
     and in the right location.";

impl Alignments {
    /// Opens a BAM, or a CRAM (by its `.cram` extension, as Sniffles2 decides), with its index.
    ///
    /// # Errors
    ///
    /// Fails if the file or its index cannot be read, or for CRAM without a reference.
    pub fn open(path: &str, reference: Option<&str>) -> Result<Self> {
        let is_cram = Path::new(path).extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        let format = if is_cram {
            let Some(reference) = reference else { bail!("CRAM input requires --reference") };
            let reference = lrpipe_reads::Reference::open(Path::new(reference))
                .with_context(|| format!("opening {reference} (and its .fai)"))?
                .repository()
                .clone();
            let index = cram::fs::read_associated_index(path)
                .with_context(|| format!("Unable to load index for input file '{path}'. {INDEX_HINT}"))?;
            let mut reader = File::open(path).map(cram::io::Reader::new).with_context(|| format!("opening {path}"))?;
            let header = reader.read_header().with_context(|| format!("reading header of {path}"))?;
            Format::Cram { header, index, reference, decoding: pass_decoding() }
        } else {
            let mut reader = File::open(path).map(bam::io::Reader::new).with_context(|| format!("opening {path}"))?;
            let header = reader.read_header().with_context(|| format!("reading header of {path}"))?;
            let index = read_index(path)
                .with_context(|| format!("Unable to load index for input file '{path}'. {INDEX_HINT}"))?;
            Format::Bam { header, index }
        };
        let header = match &format {
            Format::Bam { header, .. } | Format::Cram { header, .. } => header,
        };
        let contigs = header
            .reference_sequences()
            .iter()
            .map(|(name, rs)| (name.to_string(), usize::from(rs.length()) as i64))
            .collect();
        Ok(Alignments { path: path.to_string(), format, contigs })
    }

    /// Decodes `tags` from CRAM too, for a part of the pass that reads them (methylation).
    pub fn decode_tags(&mut self, tags: &[Tag]) {
        if let Format::Cram { decoding: Decoding { tags: lrpipe_reads::Tags::Only(only), .. }, .. } = &mut self.format {
            for &tag in tags {
                if !only.contains(&tag) {
                    only.push(tag);
                }
            }
        }
    }

    /// Drops cached reference sequences (CRAM), so that only the contig being read is kept.
    pub fn start_contig(&self) {
        if let Format::Cram { reference, .. } = &self.format {
            reference.clear();
        }
    }

    /// Extracts leads from reads on contig `tid` whose start lies in [start, end).
    /// Leads outside `region` (the Sniffles2 region being read, [start, end)) are dropped.
    /// With `depth` = (depth start, depth end, counter), every record starting in [depth start,
    /// depth end) that the query from `start` returns (that reaches `start` if it starts before),
    /// filtered or not, is also pushed to the counter (returned in [`ChunkLeads::depth`]). With
    /// `repeats`, records starting in [start, end) are offered to a repeat collector (returned
    /// in [`ChunkLeads::repeats`]). With `methylation` = (its end, piece), every record starting
    /// in [start, its end) is counted in the piece. With `read_stats` = (its start, its end),
    /// every record starting in that range that the query from `start` returns is counted in
    /// [`ChunkLeads::read_stats`].
    #[allow(clippy::too_many_arguments)]
    pub fn chunk_leads(
        &self,
        cfg: &Config,
        tid: usize,
        start: i64,
        end: i64,
        region: (i64, i64),
        read_id_base: u64,
        depth: Option<(i64, i64, lrpipe_depth::ContigCounter)>,
        repeats: Option<RepeatPiece<'_>>,
        methylation: Option<(i64, &mut lrpipe_methyl::stream::Piece<'_>)>,
        read_stats: Option<(i64, i64)>,
    ) -> Result<ChunkLeads> {
        let read = self.read_chunk(cfg, tid, start, end, region, read_id_base, depth, repeats, methylation, read_stats);
        read.with_context(|| {
            format!(
                "reading {}:{}-{end} from '{}' (truncated or corrupt file?)",
                self.contigs[tid].0,
                start + 1,
                self.path
            )
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn read_chunk(
        &self,
        cfg: &Config,
        tid: usize,
        start: i64,
        end: i64,
        region: (i64, i64),
        read_id_base: u64,
        depth: Option<(i64, i64, lrpipe_depth::ContigCounter)>,
        repeats: Option<RepeatPiece<'_>>,
        mut methylation: Option<(i64, &mut lrpipe_methyl::stream::Piece<'_>)>,
        read_stats: Option<(i64, i64)>,
    ) -> Result<ChunkLeads> {
        let mut out = ChunkLeads::default();
        let mut stats = read_stats.map(|range| (range, ReadStats::default()));
        let mut collector =
            repeats.map(|r| (lrpipe_tr::stream::Collector::new(r.loci, r.reference, r.opts), r.opts.min_mapq));
        let in_chunk = |pos: Option<i64>| pos.filter(|&p| p >= start && p < end);
        // Depth may need reads past `end`; lead extraction drops them by their start.
        let (depth_start, depth_end, mut depth) = match depth {
            Some((s, e, counter)) => (s.min(start), e.max(end), Some(counter)),
            None => (start, end, None),
        };
        // Methylation and the statistics too may need reads past `end`.
        let query_end = methylation.as_ref().map_or(depth_end, |&(e, _)| e.max(depth_end));
        let query_end = stats.as_ref().map_or(query_end, |&((_, e), _)| e.max(query_end));
        if query_end <= start {
            out.depth = depth;
            out.read_stats = stats.map(|(_, s)| s);
            return Ok(out);
        }
        let contig = &self.contigs[tid].0;
        let query_region = Region::new(
            contig.as_str(),
            Position::try_from(start as usize + 1)?..=Position::try_from(query_end as usize)?,
        );
        let in_depth_range = |pos: Option<i64>| pos.filter(|&p| p >= depth_start && p < depth_end);
        let haplotypes = match &cfg.haplotypes {
            Haplotypes::Tags => ContigHaplotypes::Tags,
            Haplotypes::List(list) => ContigHaplotypes::List(list, list.contig(contig)),
            Haplotypes::Later => ContigHaplotypes::Later,
        };
        let mut ex = Extractor {
            cfg,
            contig,
            haplotypes,
            region,
            read_id: read_id_base,
            out: &mut out,
            ins_bin_counts: HashMap::new(),
            cigar_buf: Vec::new(),
        };
        // Fetched on first use: a contig without reads need not be in the reference.
        let mut ref_seq = None;
        let mut handle = |record: &bam::Record, format: &Format| -> Result<()> {
            let mut cigar = std::mem::take(&mut ex.cigar_buf);
            decode_cigar_into(record.cigar().as_bytes(), &mut cigar);
            let aln_start = record.alignment_start().transpose()?.map(|p| usize::from(p) as i64 - 1);
            let seq = record.sequence();
            // NM is only for SV leads.
            let computed_nm = match (format, aln_start) {
                (Format::Cram { reference, .. }, Some(s))
                    if cfg.calls_svs()
                        && record.data().get(&sam::alignment::record::data::field::Tag::EDIT_DISTANCE).is_none() =>
                {
                    let r = match &mut ref_seq {
                        Some(r) => r,
                        slot @ None => slot.insert(reference.get(contig.as_bytes()).transpose()?),
                    };
                    let bases: Vec<u8> = seq.iter().collect();
                    r.as_ref().map(|r| edit_distance_to_ref(&cigar, &bases, (**r).as_ref(), s))
                }
                _ => None,
            };
            let view = RecordView {
                flags: record.flags().bits(),
                start: aln_start,
                mapq: record.mapping_quality().map_or(255, u8::from),
                name: record.name().map_or(&b""[..], |n| n.as_ref()),
                seq: Seq { bytes: seq.as_bytes(), len: seq.len() },
                tags: record.data().as_bytes(),
                computed_nm,
            };
            if let Some(((stats_start, stats_end), stats)) = stats.as_mut()
                && view.start.is_some_and(|p| p >= *stats_start && p < *stats_end)
            {
                stats.add(view.flags, view.seq.len, qualities(record, format));
            }
            if let Some(counter) = depth.as_mut()
                && let Some(pos) = in_depth_range(view.start)
            {
                counter.push(&lrpipe_depth::Alignment {
                    tid,
                    pos,
                    flags: view.flags,
                    // CRAM stores no mapping quality for unmapped reads; htslib reports 0.
                    mapq: match format {
                        Format::Cram { .. } if view.flags & FLAG_UNMAPPED != 0 && view.mapq == 255 => 0,
                        _ => view.mapq,
                    },
                    mate_tid: record.mate_reference_sequence_id().transpose()?,
                    mate_pos: record.mate_alignment_start().transpose()?.map_or(-1, |p| usize::from(p) as i64 - 1),
                    template_len: record.template_length(),
                    name: view.name,
                    cigar: lrpipe_depth::Cigar::Packed(record.cigar().as_bytes()),
                    read_group: None,
                });
            }
            if let (Some((c, q)), Some(pos)) = (collector.as_mut(), in_chunk(view.start)) {
                offer_repeat(c, &view, &cigar, pos, *q, haplotypes)?;
            }
            if let Some((methyl_end, piece)) = methylation.as_mut()
                && view.start.is_some_and(|p| p >= start && p < *methyl_end)
            {
                piece.add(record)?;
            }
            let r = if cfg.calls_svs() { ex.process(&view, &cigar, start, end) } else { Ok(()) };
            ex.cigar_buf = cigar;
            r
        };
        self.each_record(tid, &query_region, &mut handle)?;
        out.depth = depth;
        out.read_stats = stats.map(|(_, s)| s);
        if let Some((c, _)) = collector {
            out.repeats = c.into_observations();
        }
        Ok(out)
    }

    /// Counts every record of contigs `tids` in `stats` (`--read-stats` for the contigs the pass
    /// does not otherwise read: the short ones, a whole genome's hundreds or thousands).
    ///
    /// From BAM, a query of each. From CRAM, one pass over the slices holding their records, in
    /// file order, each slice decoded once though it holds many contigs, and only for its
    /// records' flags, reference and CIGAR (their lengths): no reference sequence is read, and
    /// none is kept.
    ///
    /// # Errors
    ///
    /// If the reads cannot be read.
    pub fn count_contigs(&self, tids: &[usize], stats: &mut ReadStats) -> Result<()> {
        match &self.format {
            Format::Bam { .. } => tids.iter().try_for_each(|&tid| self.count_contig(tid, stats)),
            Format::Cram { header, index, .. } => count_cram_contigs(&self.path, header, index, tids, stats)
                .with_context(|| format!("reading the short contigs of '{}' (truncated or corrupt file?)", self.path)),
        }
    }

    /// Counts every record of contig `tid` of a BAM in `stats`.
    fn count_contig(&self, tid: usize, stats: &mut ReadStats) -> Result<()> {
        let (name, len) = &self.contigs[tid];
        if *len < 1 {
            return Ok(());
        }
        let region = Region::new(name.as_str(), Position::MIN..=Position::try_from(*len as usize)?);
        self.each_record(tid, &region, |record, format| {
            stats.add(record.flags().bits(), record.sequence().len(), qualities(record, format));
            Ok(())
        })
        .with_context(|| format!("reading {name} from '{}' (truncated or corrupt file?)", self.path))
    }

    /// Calls `handle` on each record of contig `tid` overlapping `region`, in order.
    fn each_record(
        &self,
        tid: usize,
        region: &Region,
        mut handle: impl FnMut(&bam::Record, &Format) -> Result<()>,
    ) -> Result<()> {
        match &self.format {
            format @ Format::Bam { header, index } => {
                let mut reader = File::open(&self.path).map(bam::io::Reader::new)?;
                let mut query = reader.query(header, index, region)?;
                let mut record = bam::Record::default();
                while query.read_record(&mut record)? != 0 {
                    handle(&record, format)?;
                }
            }
            format @ Format::Cram { header, index, reference, decoding } => {
                // Records come encoded as BAM, with what the pass reads.
                let source = CramSource { header, index, reference, decoding };
                let mut reader = File::open(&self.path).map(cram::io::Reader::new)?;
                let mut failure = None;
                source
                    .query(&mut reader, tid, region.interval(), |record| {
                        handle(record, format).map_err(|e| {
                            failure = Some(e);
                            io::Error::other("stopped")
                        })
                    })
                    .or_else(|e| failure.take().map_or(Err(e.into()), Err))?;
            }
        }
        Ok(())
    }
}

/// [`Alignments::count_contigs`] of a CRAM: the slices its index lists for `tids`, each once,
/// their records decoded for their flags, reference and CIGAR alone.
fn count_cram_contigs(
    path: &str,
    header: &sam::Header,
    index: &cram::crai::Index,
    tids: &[usize],
    stats: &mut ReadStats,
) -> Result<()> {
    use pasta::sam::alignment::Record as _;
    let wanted: std::collections::HashSet<usize> = tids.iter().copied().collect();
    // (container offset, slice landmark): a multi-reference slice has a line for each contig.
    let mut slices: Vec<(u64, u64)> = index
        .iter()
        .filter(|r| r.reference_sequence_id().is_some_and(|t| wanted.contains(&t)))
        .map(|r| (r.offset(), r.landmark()))
        .collect();
    slices.sort_unstable();
    slices.dedup();
    let fields = Fields::FLAGS | Fields::REFERENCE_SEQUENCE_ID | Fields::CIGAR;
    let tags = lrpipe_reads::Tags::Only(Vec::new());
    let mut reader = File::open(path).map(cram::io::Reader::new)?;
    let mut container = cram::io::reader::Container::default();
    for group in slices.chunk_by(|a, b| a.0 == b.0) {
        reader.seek(io::SeekFrom::Start(group[0].0))?;
        if reader.read_container(&mut container)? == 0 {
            bail!("no container at byte {} of the index", group[0].0);
        }
        let compression = container.compression_header()?;
        let landmarks = container.header().landmarks().to_vec();
        for (slice, landmark) in container.slices().zip(landmarks) {
            let slice = slice?;
            if !group.iter().any(|&(_, l)| l == landmark as u64) {
                continue;
            }
            let (core, external) = slice.decode_blocks_for(&compression, fields, &tags)?;
            let records = slice.records_with_fields(
                fasta::Repository::default(),
                header,
                &compression,
                &core,
                &external,
                MdNm::AsStored,
                fields,
                &tags,
            )?;
            for record in &records {
                if !record.reference_sequence_id(header).transpose()?.is_some_and(|t| wanted.contains(&t)) {
                    continue;
                }
                // The read's length: the bases its CIGAR consumes of it (SEQ's, by SAM's rules).
                let mut length = 0;
                for op in record.cigar().iter() {
                    let op = op?;
                    if op.kind().consumes_read() {
                        length += op.len();
                    }
                }
                stats.add(record.flags()?.bits(), length, &[]);
            }
        }
    }
    Ok(())
}

/// A record's base qualities for `--read-stats`: none from CRAM, whose qualities the pass does
/// not decode (reading the placeholder its records hold took 8% of the pass's CPU on HG002 ONT
/// chr20).
fn qualities<'r>(record: &'r bam::Record, format: &Format) -> &'r [u8] {
    match format {
        Format::Bam { .. } => record.quality_scores().as_bytes(),
        Format::Cram { .. } => &[],
    }
}

/// The fields of an alignment record that lead extraction reads.
struct RecordView<'a> {
    flags: u16,
    /// 0-based alignment start.
    start: Option<i64>,
    mapq: u8,
    name: &'a [u8],
    seq: Seq<'a>,
    /// BAM auxiliary fields.
    tags: &'a [u8],
    /// NM computed from the reference, for CRAM records whose NM was not stored (htslib
    /// regenerates it when decoding, so Sniffles2 sees it).
    computed_nm: Option<i64>,
}

/// Read bases, 4-bit packed as BAM stores them.
#[derive(Clone, Copy)]
struct Seq<'a> {
    bytes: &'a [u8],
    len: usize,
}

impl Seq<'_> {
    fn len(&self) -> usize {
        self.len
    }

    /// Python-style slice `seq[a:b]` (bounds clipped).
    fn slice(&self, a: i64, b: i64) -> Vec<u8> {
        seq_slice(self.bytes, self.len, a, b)
    }
}

/// The alignment's edit distance (NM) as htslib derives it when decoding CRAM: aligned bases that
/// differ from the reference (case-insensitive), plus inserted and deleted bases.
fn edit_distance_to_ref(cigar: &[(u8, i64)], seq: &[u8], reference: &[u8], start: i64) -> i64 {
    let (mut q, mut r, mut nm) = (0usize, start.max(0) as usize, 0i64);
    for &(op, len) in cigar {
        let n = len as usize;
        match op {
            0 | 7 | 8 => {
                let read = seq.get(q..q + n).unwrap_or_default();
                let refs = reference.get(r..(r + n).min(reference.len())).unwrap_or_default();
                nm += read.iter().zip(refs).filter(|(a, b)| !a.eq_ignore_ascii_case(b)).count() as i64;
                q += n;
                r += n;
            }
            CINS => {
                nm += len;
                q += n;
            }
            CDEL => {
                nm += len;
                r += n;
            }
            3 => r += n,
            CSOFT_CLIP => q += n,
            _ => {}
        }
    }
    nm
}

/// The BAM's index: `PATH.bai`, `PATH` with `.bai` for `.bam`, then the same with `.csi`, as
/// `lrpipe run` and htslib accept them. Sniffles2 (pysam) reads a CSI too; this reader read
/// only a BAI (LRP-009).
fn read_index(path: &str) -> io::Result<Box<dyn BinningIndex + Send + Sync>> {
    let candidates = |ext: &str| [PathBuf::from(format!("{path}.{ext}")), Path::new(path).with_extension(ext)];
    if let Some(bai) = candidates("bai").into_iter().find(|p| p.exists()) {
        return Ok(Box::new(bam::bai::fs::read(bai)?));
    }
    if let Some(csi) = candidates("csi").into_iter().find(|p| p.exists()) {
        return Ok(Box::new(csi::fs::read(csi)?));
    }
    Err(io::Error::new(io::ErrorKind::NotFound, format!("no {path}.bai or {path}.csi")))
}

/// Decodes a packed BAM CIGAR into (operation, length) pairs, reusing `out`.
fn decode_cigar_into(raw: &[u8], out: &mut Vec<(u8, i64)>) {
    out.clear();
    out.extend(raw.as_chunks::<4>().0.iter().map(|&c| {
        let v = u32::from_le_bytes(c);
        ((v & 0xF) as u8, i64::from(v >> 4))
    }));
}

/// pysam's `query_alignment_start` and `query_alignment_end`: the read range without soft clips
/// (hard clips are skipped over).
fn soft_clipped_range(cigar: &[(u8, i64)], seq_len: usize) -> (i64, i64) {
    let clip = |ops: &mut dyn Iterator<Item = &(u8, i64)>| {
        let mut n = 0;
        for &(op, len) in ops {
            match op {
                CSOFT_CLIP => n += len,
                CHARD_CLIP => continue,
                _ => break,
            }
        }
        n
    };
    (clip(&mut cigar.iter()), seq_len as i64 - clip(&mut cigar.iter().rev()))
}

/// Reference bases covered by an alignment (M, D, N, =, X).
fn ref_span(cigar: &[(u8, i64)]) -> i64 {
    cigar.iter().filter(|(op, _)| matches!(op, 0 | 2 | 3 | 7 | 8)).map(|&(_, l)| l).sum()
}

/// Sniffles2's per-read mismatch rate: NM without indels longer than 10 bp, per aligned base.
fn adjusted_nm(nm_raw: i64, cigar: &[(u8, i64)], alen: i64) -> f64 {
    let large: i64 = cigar.iter().filter(|&&(op, len)| (op == CINS || op == CDEL) && len > 10).map(|&(_, l)| l).sum();
    (nm_raw - large) as f64 / (alen + 1) as f64
}

/// Tags of a record that Sniffles2 looks at.
#[derive(Default)]
struct Tags<'a> {
    sa: Option<&'a [u8]>,
    nm: Option<i64>,
    hp: Option<i64>,
    ps: Option<String>,
}

fn parse_tags(mut d: &[u8]) -> Result<Tags<'_>> {
    let mut t = Tags::default();
    while d.len() >= 3 {
        let tag = [d[0], d[1]];
        let ty = d[2];
        d = &d[3..];
        // Integer value of a tag whose bytes `d` were cut to the type's size.
        let int_val = |d: &[u8]| -> i64 {
            match ty {
                b'c' => i64::from(d[0] as i8),
                b'C' => i64::from(d[0]),
                b's' => i64::from(i16::from_le_bytes([d[0], d[1]])),
                b'S' => i64::from(u16::from_le_bytes([d[0], d[1]])),
                b'i' => i64::from(i32::from_le_bytes([d[0], d[1], d[2], d[3]])),
                _ => i64::from(u32::from_le_bytes([d[0], d[1], d[2], d[3]])),
            }
        };
        let size = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            // HTSlib's non-standard double, which samtools writes to BAM.
            b'd' => 8,
            b'Z' | b'H' => d.iter().position(|&b| b == 0).context("unterminated string tag")? + 1,
            b'B' => {
                let (&sub, n) = d.split_first().context("truncated array tag")?;
                let n = n.first_chunk::<4>().context("truncated array tag")?;
                let w = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => bail!("unknown array tag subtype {}", sub as char),
                };
                (u32::from_le_bytes(*n) as usize)
                    .checked_mul(w)
                    .and_then(|n| n.checked_add(5))
                    .context("array tag too long")?
            }
            _ => bail!("unknown tag type {}", ty as char),
        };
        let v = d.get(..size).with_context(|| format!("truncated {} tag", String::from_utf8_lossy(&tag)))?;
        let is_int = matches!(ty, b'c' | b'C' | b's' | b'S' | b'i' | b'I');
        match &tag {
            b"SA" if ty == b'Z' => t.sa = Some(&v[..size - 1]),
            b"NM" if is_int => t.nm = Some(int_val(v)),
            b"HP" if is_int => t.hp = Some(int_val(v)),
            b"PS" if is_int => t.ps = Some(int_val(v).to_string()),
            b"PS" if ty == b'Z' => t.ps = Some(String::from_utf8_lossy(&v[..size - 1]).into_owned()),
            _ => {}
        }
        d = &d[size..];
    }
    Ok(t)
}

const SEQ_CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Python-style slice `seq[a:b]` of a 4-bit packed sequence.
fn seq_slice(packed: &[u8], len: usize, a: i64, b: i64) -> Vec<u8> {
    let a = a.clamp(0, len as i64) as usize;
    let b = b.clamp(0, len as i64) as usize;
    (a..b.max(a)).map(|i| SEQ_CODES[((packed[i / 2] >> (4 * (1 - i % 2))) & 0xF) as usize]).collect()
}

/// `CIGAR_analyze`: (clip at start, clip at end, refspan, readspan) of an SA CIGAR string.
fn cigar_analyze(cigar: &[u8]) -> Option<(i64, i64, i64, i64)> {
    let mut num: Option<i64> = None;
    let (mut readspan, mut refspan, mut clip) = (0, 0, 0);
    let mut clip_start: Option<i64> = None;
    for &c in cigar {
        if c.is_ascii_digit() {
            num = Some(num.unwrap_or(0) * 10 + (c - b'0') as i64);
            continue;
        }
        let oplen = num.take()?;
        let mut h = false;
        if b"MIX=".contains(&c) {
            readspan += oplen;
            h = true;
        }
        if b"MDX=N".contains(&c) {
            refspan += oplen;
            h = true;
        }
        if !h {
            if b"SH".contains(&c) {
                if clip_start.is_none() && readspan + refspan > 0 {
                    clip_start = Some(clip);
                }
                clip += oplen;
            } else {
                return None;
            }
        }
    }
    let clip_start = clip_start.unwrap_or(clip);
    Some((clip_start, clip - clip_start, refspan, readspan))
}

struct SaEntry<'a> {
    refname: &'a str,
    pos: &'a str,
    strand: &'a str,
    cigar: &'a str,
    mapq: &'a str,
    nm: &'a str,
}

fn parse_sa(sa: &[u8]) -> Result<Vec<SaEntry<'_>>> {
    let s = std::str::from_utf8(sa)?;
    s.split(';')
        .filter(|p| !p.is_empty())
        .map(|p| {
            let f: Vec<&str> = p.split(',').collect();
            if f.len() != 6 {
                bail!("malformed SA tag entry '{p}'");
            }
            Ok(SaEntry { refname: f[0], pos: f[1], strand: f[2], cigar: f[3], mapq: f[4], nm: f[5] })
        })
        .collect()
}

/// A split segment of a read (`SPLIT_PRIM` / `SPLIT_SUP` lead) during split classification.
#[derive(Debug, Clone)]
struct Split {
    contig: Arc<str>,
    ref_start: i64,
    ref_end: i64,
    qry_start: i64,
    qry_end: i64,
    reverse: bool,
    mapq: i64,
    source: Source,
    seq: Option<Vec<u8>>,
    svtypes: Vec<(SvType, i64, Option<i64>)>,
}

/// The parts of one alignment record that lead extraction needs.
struct Aln<'r> {
    cigar: &'r [(u8, i64)],
    ref_start: i64,
    ref_end: i64,
    /// Query range of the alignment, excluding soft clips (`query_alignment_start/end`).
    qstart: i64,
    qend: i64,
    reverse: bool,
    is_supp: bool,
    mapq: u8,
    seq: Seq<'r>,
    tags: &'r Tags<'r>,
}

struct Extractor<'a> {
    cfg: &'a Config,
    contig: &'a str,
    /// Where the reads' haplotypes come from.
    haplotypes: ContigHaplotypes<'a>,
    /// Sniffles2 region being read, [start, end): leads outside it are "externals".
    region: (i64, i64),
    read_id: u64,
    out: &'a mut ChunkLeads,
    /// INS leads per bin so far in this chunk.
    ins_bin_counts: HashMap<i64, usize>,
    /// Reused CIGAR storage.
    cigar_buf: Vec<(u8, i64)>,
}

impl Extractor<'_> {
    fn emit(&mut self, svtype: SvType, contig: &str, lead: Lead) {
        // Leads outside the task region are "externals" in Sniffles2 and never used.
        if contig != self.contig || lead.ref_start < self.region.0 || lead.ref_start >= self.region.1 {
            return;
        }
        let bin = lead.ref_start / self.cfg.args.cluster_binsize * self.cfg.args.cluster_binsize;
        let mut lead = lead;
        if svtype == SvType::Ins {
            // Sniffles2 keeps the sequence of only the first leads of a bin. A lead's index within
            // this chunk's part of the bin is a lower bound of its index in the whole bin, so the
            // sequence can be dropped now (Leadtab::new applies the exact rule).
            let n = self.ins_bin_counts.entry(bin).or_insert(0);
            if *n >= self.cfg.args.consensus_max_reads_bin {
                lead.seq = None;
            }
            *n += 1;
        }
        self.out.leads[svtype as usize].push((bin, lead));
    }

    /// Records a single break lead (0: left end of the alignment, 1: right end).
    fn emit_single(&mut self, side: usize, lead: Lead) {
        if lead.ref_start < self.region.0 || lead.ref_start >= self.region.1 {
            return;
        }
        let bin = lead.ref_start / self.cfg.args.cluster_binsize * self.cfg.args.cluster_binsize;
        self.out.single[side].push((bin, lead));
    }

    fn process(&mut self, rec: &RecordView<'_>, cigar: &[(u8, i64)], start: i64, end: i64) -> Result<()> {
        let cfg = self.cfg;
        let flags = rec.flags;
        if flags & FLAG_UNMAPPED != 0 || cigar.is_empty() {
            return Ok(());
        }
        let mapq = rec.mapq;
        let seq_len = rec.seq.len();
        let (qstart, qend) = soft_clipped_range(cigar, seq_len);
        let alen = qend - qstart;
        if mapq < cfg.args.mapq || flags & FLAG_SECONDARY != 0 || alen < cfg.args.min_alignment_length {
            return Ok(());
        }
        if let Some(ex) = cfg.args.exclude_flags
            && flags & ex != 0
        {
            return Ok(());
        }
        let Some(ref_start) = rec.start else { return Ok(()) };
        if ref_start < start || ref_start >= end {
            return Ok(());
        }
        let ref_end = ref_start + ref_span(cigar);
        let mut tags = parse_tags(rec.tags)?;
        if tags.nm.is_none() {
            tags.nm = rec.computed_nm;
        }
        match self.haplotypes {
            ContigHaplotypes::Tags => {}
            ContigHaplotypes::List(list, listed) => {
                (tags.hp, tags.ps) = ContigHaplotypes::listed(list, listed, rec.name, flags).unzip();
            }
            ContigHaplotypes::Later => (tags.hp, tags.ps) = (None, None),
        }
        let hp = tags.hp.unwrap_or(0);
        if !(0..=2).contains(&hp) {
            bail!("HP tag value {hp} is not supported (Sniffles2 expects 0, 1 or 2)");
        }
        let hp = hp as u8;
        let ps: Arc<str> = Arc::from(tags.ps.as_deref().unwrap_or("NULL"));

        self.read_id += 1;
        self.out.read_count += 1;
        self.out.spans.push((ref_start, ref_end));

        // Always computed: Sniffles2 phases by default (advanced tags).
        let nm = match tags.nm {
            Some(nm_raw) => {
                let nm = adjusted_nm(nm_raw, cigar, alen);
                self.out.nms.push(nm);
                nm
            }
            None => -1.0,
        };
        let aln = Aln {
            cigar,
            ref_start,
            ref_end,
            qstart,
            qend,
            reverse: flags & FLAG_REVERSE != 0,
            is_supp: flags & FLAG_SUPPLEMENTARY != 0,
            mapq,
            seq: rec.seq,
            tags: &tags,
        };
        let base = Lead {
            read_id: self.read_id,
            qname: Arc::from(rec.name),
            ref_start: 0,
            ref_end: 0,
            qry_start: 0,
            qry_end: 0,
            reverse: aln.reverse,
            mapq,
            nm,
            source: Source::Inline,
            svlen: None,
            seq: None,
            bnd: None,
            hap: hp,
            phase_set: Some(ps),
            is_sa: aln.is_supp,
            read_len: alen,
        };

        self.read_iterindels(&base, &aln);
        if let Some(sa) = aln.tags.sa {
            let sa = parse_sa(sa)?;
            if let Some(lead) = lead_for_bnd(&base, &sa, &aln) {
                self.emit(SvType::Bnd, self.contig, lead);
            }
            if !aln.is_supp {
                self.read_itersplits(&base, &sa, &aln)?;
            }
        }

        let bs = cfg.args.cluster_binsize;
        self.out.hap_bins.push((ref_start / bs * bs, ref_end / bs * bs, hp));
        if matches!(self.haplotypes, ContigHaplotypes::Later) {
            self.out.hap_reads.push((!aln.is_supp).then(|| Arc::clone(&base.qname)));
        }
        Ok(())
    }

    /// `read_iterindels`: INS and DEL leads from the CIGAR, and long-insertion leads from soft
    /// clips of reads without split alignments.
    fn read_iterindels(&mut self, base: &Lead, aln: &Aln<'_>) {
        let cfg = self.cfg;
        let use_clips = cfg.args.detect_large_ins && !aln.is_supp && aln.tags.sa.is_none();
        let minsvlen = cfg.minsvlen_screen;
        let longinslen = cfg.args.long_ins_length as f64 / 2.0;
        let (mut pos_read, mut pos_ref) = (0i64, aln.ref_start);
        for &(op, len) in aln.cigar {
            let (add_read, add_ref, event) = match op {
                0 | 7 | 8 => (1, 1, false),
                CINS => (1, 0, true),
                CDEL => (0, 1, true),
                3 => (0, 1, false),
                CSOFT_CLIP => (1, 0, true),
                _ => (0, 0, false),
            };
            if event && len >= minsvlen {
                if op == CINS {
                    let s = (len <= cfg.args.dev_seq_cache_maxlen).then(|| aln.seq.slice(pos_read, pos_read + len));
                    let lead = Lead {
                        ref_start: pos_ref,
                        ref_end: pos_ref,
                        qry_start: pos_read,
                        qry_end: pos_read + len,
                        svlen: Some(len),
                        seq: s,
                        ..base.clone()
                    };
                    self.emit(SvType::Ins, self.contig, lead);
                } else if op == CDEL {
                    let lead = Lead {
                        ref_start: pos_ref + len,
                        ref_end: pos_ref,
                        qry_start: pos_read,
                        qry_end: pos_read,
                        svlen: Some(-len),
                        ..base.clone()
                    };
                    self.emit(SvType::Del, self.contig, lead);
                } else if use_clips && op == CSOFT_CLIP && len as f64 >= longinslen {
                    let lead = Lead {
                        ref_start: pos_ref,
                        ref_end: pos_ref,
                        qry_start: pos_read,
                        qry_end: pos_read + len,
                        svlen: None,
                        ..base.clone()
                    };
                    self.emit(SvType::Ins, self.contig, lead);
                } else if op == CSOFT_CLIP && cfg.args.snf.is_some() {
                    // Other soft clips are SINGLE_LEFT/RIGHT breaks: never in the VCF, but their
                    // candidates decide which SNF blocks exist.
                    let lead = Lead {
                        ref_start: pos_ref,
                        ref_end: pos_ref,
                        qry_start: pos_read,
                        qry_end: pos_read + len,
                        svlen: Some(0),
                        ..base.clone()
                    };
                    self.emit_single(usize::from(pos_ref != aln.ref_start), lead);
                }
            }
            pos_read += add_read * len;
            pos_ref += add_ref * len;
        }
    }

    /// `read_itersplits`: INS/DEL/DUP/INV leads from a primary alignment and its supplementaries.
    fn read_itersplits(&mut self, base: &Lead, sa: &[SaEntry<'_>], aln: &Aln<'_>) -> Result<()> {
        let cfg = self.cfg;
        let seq_len = aln.seq.len();
        if sa.len() as f64 > cfg.args.max_splits_base as f64 + cfg.args.max_splits_kb * (seq_len as f64 / 1000.0) {
            return Ok(());
        }
        let qry_start = if aln.reverse { seq_len as i64 - aln.qend } else { aln.qstart };
        let contig: Arc<str> = Arc::from(self.contig);
        let mut splits = vec![Split {
            contig: contig.clone(),
            ref_start: aln.ref_start,
            ref_end: aln.ref_end,
            qry_start,
            qry_end: qry_start + (aln.qend - aln.qstart),
            reverse: aln.reverse,
            mapq: i64::from(aln.mapq),
            source: Source::SplitPrim,
            seq: None,
            svtypes: Vec::new(),
        }];
        for e in sa {
            let Ok(m) = e.mapq.parse::<i64>() else { bail!("malformed SA mapq '{}'", e.mapq) };
            let is_rev = e.strand == "-";
            let Some((fwd, rev, refspan, readspan)) = cigar_analyze(e.cigar.as_bytes()) else {
                report_bad_sa_cigar(e, base);
                return Ok(());
            };
            let pos_zero: i64 = e.pos.parse::<i64>()? - 1;
            let q = if is_rev { rev } else { fwd };
            splits.push(Split {
                contig: if e.refname == self.contig { contig.clone() } else { Arc::from(e.refname) },
                ref_start: pos_zero,
                ref_end: pos_zero + refspan,
                qry_start: q,
                qry_end: q + readspan,
                reverse: is_rev,
                mapq: m,
                source: Source::SplitSup,
                seq: None,
                svtypes: Vec::new(),
            });
        }
        let splits = classify_splits(cfg, &aln.seq, splits);
        for (i, s) in splits.iter().enumerate() {
            for &(svtype, svstart, svlen) in &s.svtypes {
                let min_mapq = s.mapq.min(splits[i.saturating_sub(1)].mapq);
                if !cfg.args.dev_keep_lowqual_splits && min_mapq < i64::from(cfg.args.mapq) {
                    continue;
                }
                let ref_end = match svlen {
                    Some(l) if svtype != SvType::Ins => svstart + l,
                    _ => svstart,
                };
                let lead = Lead {
                    ref_start: svstart,
                    ref_end,
                    qry_start: s.qry_start,
                    qry_end: s.qry_end,
                    reverse: s.reverse,
                    mapq: s.mapq.clamp(0, 255) as u8,
                    source: s.source,
                    svlen,
                    seq: if svtype == SvType::Ins { s.seq.clone() } else { None },
                    is_sa: false,
                    read_len: 0,
                    ..base.clone()
                };
                self.emit(svtype, &s.contig, lead);
            }
        }
        Ok(())
    }
}

fn report_bad_sa_cigar(e: &SaEntry<'_>, base: &Lead) {
    eprintln!(
        "Sniffles2 Error: Malformed CIGAR '{}' with pos {} of read '{}'. Skipping.",
        e.cigar,
        e.pos,
        String::from_utf8_lossy(&base.qname)
    );
}

/// `Lead.for_bnd`: a BND lead from the first SA entry, if it is on the opposite strand.
fn lead_for_bnd(base: &Lead, sa: &[SaEntry<'_>], aln: &Aln<'_>) -> Option<Lead> {
    let clip = |&(op, len): &(u8, i64)| if op == CSOFT_CLIP || op == CHARD_CLIP { len } else { 0 };
    let left = aln.cigar.first().map_or(0, clip);
    let right = aln.cigar.last().map_or(0, clip);
    let (bnd_start, is_first) = if left > right { (aln.ref_start + 1, false) } else { (aln.ref_end, true) };
    let read_strand = if aln.reverse { "-" } else { "+" };
    let e = sa.first()?;
    if e.strand == read_strand {
        return None;
    }
    let pos: i64 = e.pos.parse().ok()?;
    let pos = pos - 1;
    let Some((l, r, refspan, _)) = cigar_analyze(e.cigar.as_bytes()) else {
        report_bad_sa_cigar(e, base);
        return None;
    };
    let is_rev = r > l;
    let mate_ref_start = if is_rev || !is_first && e.strand == read_strand {
        pos + refspan
    } else if is_first {
        pos + 1
    } else {
        pos + 2
    };
    let nm = if aln.tags.nm.is_some() { e.nm.parse::<i64>().map_or(f64::NAN, |v| v as f64) } else { f64::NAN };
    Some(Lead {
        ref_start: bnd_start,
        ref_end: bnd_start,
        qry_start: aln.qstart,
        qry_end: aln.qend,
        nm,
        source: Source::BndSa,
        svlen: Some(0),
        seq: None,
        bnd: Some(Box::new(BndInfo {
            mate_contig: Arc::from(e.refname),
            mate_ref_start,
            is_first,
            is_reverse: is_rev,
        })),
        hap: 0,
        phase_set: None,
        is_sa: false,
        read_len: 0,
        ..base.clone()
    })
}

/// `sv.classify_splits`: SV type hints between consecutive segments of a split read.
fn classify_splits(cfg: &Config, seq: &Seq<'_>, mut leads: Vec<Split>) -> Vec<Split> {
    let minsv = cfg.minsvlen_screen;
    leads.sort_by_key(|l| l.qry_start);
    let mut hints = 0;
    leads[0].svtypes.clear();
    if leads[0].qry_start as f64 >= cfg.args.long_ins_length as f64 * 0.5 {
        let rs = leads[0].ref_start;
        leads[0].svtypes.push((SvType::Ins, rs, None));
    }
    for i in 1..leads.len() {
        let (last, rest) = leads.split_at_mut(i);
        let last = &last[i - 1];
        let curr = &mut rest[0];
        curr.svtypes.clear();
        if curr.contig != last.contig {
            continue;
        }
        let rev = curr.reverse;
        let fwd = !rev;
        let mut push = |curr: &mut Split, t: SvType, start: i64, len: i64| {
            curr.svtypes.push((t, start, Some(len)));
            hints += 1;
        };
        if curr.reverse == last.reverse {
            let qgap = curr.qry_start - last.qry_end;
            if fwd && qgap >= minsv && qgap - (curr.ref_start - last.ref_end) >= minsv {
                let svlen = qgap;
                curr.seq = (svlen <= cfg.args.dev_seq_cache_maxlen).then(|| seq.slice(last.qry_end, curr.qry_start));
                let s = curr.ref_start;
                push(curr, SvType::Ins, s, svlen);
            } else if rev && qgap >= minsv && qgap - (last.ref_start - curr.ref_end) >= minsv {
                let svlen = qgap;
                curr.seq = (svlen <= cfg.args.dev_seq_cache_maxlen).then(|| seq.slice(last.qry_end, curr.qry_start));
                push(curr, SvType::Ins, last.ref_start, svlen);
            } else if fwd && curr.ref_start - last.ref_end >= minsv && (curr.ref_start - last.ref_end) - qgap >= minsv {
                let s = curr.ref_start;
                push(curr, SvType::Del, s, -(curr.ref_start - last.ref_end));
            } else if rev && last.ref_start - curr.ref_end >= minsv && (last.ref_start - curr.ref_end) - qgap >= minsv {
                push(curr, SvType::Del, last.ref_start, -(last.ref_start - curr.ref_end));
            } else if fwd && curr.ref_start <= last.ref_end {
                let svlen = last.ref_end - curr.ref_start;
                if svlen >= minsv {
                    let s = curr.ref_start;
                    push(curr, SvType::Dup, s, svlen);
                }
            } else if rev && last.ref_start <= curr.ref_end {
                let svlen = curr.ref_end - last.ref_start;
                if svlen >= minsv {
                    push(curr, SvType::Dup, last.ref_start, svlen);
                }
            }
        } else if fwd && curr.ref_start <= last.ref_start {
            let svlen = last.ref_start - curr.ref_start;
            if svlen >= minsv {
                let s = curr.ref_start;
                push(curr, SvType::Inv, s, svlen);
            }
        } else if fwd && curr.ref_start > last.ref_start {
            let svlen = curr.ref_start - last.ref_start;
            if svlen >= minsv {
                push(curr, SvType::Inv, last.ref_start, svlen);
            }
        } else if rev && curr.ref_end >= last.ref_end {
            let svlen = curr.ref_end - last.ref_end;
            if svlen >= minsv {
                push(curr, SvType::Inv, last.ref_end, svlen);
            }
        } else if rev && curr.ref_end < last.ref_end {
            let svlen = last.ref_end - curr.ref_end;
            if svlen >= minsv {
                let s = curr.ref_end;
                push(curr, SvType::Inv, s, svlen);
            }
        }
    }
    if hints == 0 && leads.len() > 2 {
        let (c, r) = (leads[0].contig.clone(), leads[0].reverse);
        let filtered: Vec<Split> = leads.into_iter().filter(|l| l.contig == c && l.reverse == r).collect();
        if filtered.len() == 2 {
            return classify_splits(cfg, seq, filtered);
        }
        return filtered;
    }
    leads
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack(ops: &[(u8, u32)]) -> Vec<u8> {
        ops.iter().flat_map(|&(op, len)| ((len << 4) | u32::from(op)).to_le_bytes()).collect()
    }

    #[test]
    fn computes_nm_like_htslib() {
        // 2S, 4M with one mismatch, 2I, 3D, then 2M matching a lowercase reference base
        let reference = b"xxACGTacgtACGT";
        let cigar = [(CSOFT_CLIP, 2), (0, 4), (CINS, 2), (CDEL, 3), (0, 2)];
        let seq = b"NNACTTGGTA";
        assert_eq!(edit_distance_to_ref(&cigar, seq, reference, 2), 1 + 2 + 3);
    }

    #[test]
    fn decodes_cigar_and_derived_ranges() {
        // 5H 10S 100M 60I 50M 70D 30M 20S
        let mut cigar = Vec::new();
        decode_cigar_into(&pack(&[(5, 5), (4, 10), (0, 100), (1, 60), (0, 50), (2, 70), (0, 30), (4, 20)]), &mut cigar);
        assert_eq!(cigar[1], (CSOFT_CLIP, 10));
        let seq_len = 10 + 100 + 60 + 50 + 30 + 20;
        assert_eq!(soft_clipped_range(&cigar, seq_len), (10, seq_len as i64 - 20));
        assert_eq!(ref_span(&cigar), 100 + 50 + 70 + 30);
        // NM 140 minus the two indels over 10 bp, per aligned base + 1.
        assert_eq!(adjusted_nm(140, &cigar, 240), 10.0 / 241.0);
    }

    #[test]
    fn analyzes_sa_cigars() {
        assert_eq!(cigar_analyze(b"100S50M20I30M"), Some((100, 0, 80, 100)));
        assert_eq!(cigar_analyze(b"40M60S"), Some((0, 60, 40, 40)));
        assert_eq!(cigar_analyze(b"10H40M5D10M3H"), Some((10, 3, 55, 50)));
        assert_eq!(cigar_analyze(b"40Q"), None);
        assert_eq!(cigar_analyze(b"M"), None);
    }

    #[test]
    fn slices_packed_sequence_like_python() {
        // ACGTN packed in 4-bit codes: A=1 C=2 G=4 T=8 N=15.
        let packed = [0x12, 0x48, 0xF0];
        assert_eq!(seq_slice(&packed, 5, 0, 5), b"ACGTN");
        assert_eq!(seq_slice(&packed, 5, 1, 3), b"CG");
        assert_eq!(seq_slice(&packed, 5, 3, 99), b"TN");
        assert_eq!(seq_slice(&packed, 5, 4, 2), b"");
        assert_eq!(seq_slice(&packed, 5, -3, 2), b"AC");
    }

    #[test]
    fn parses_relevant_tags() {
        let mut d = Vec::new();
        d.extend_from_slice(b"NMC\x07");
        d.extend_from_slice(b"XYZskip\0");
        d.extend_from_slice(b"HPc\x02");
        d.extend_from_slice(b"PSi");
        d.extend_from_slice(&123456i32.to_le_bytes());
        d.extend_from_slice(b"SAZchr1,100,+,50M,60,0;\0");
        d.extend_from_slice(b"MLBC");
        d.extend_from_slice(&3u32.to_le_bytes());
        d.extend_from_slice(&[1, 2, 3]);
        let t = parse_tags(&d).expect("valid tags");
        assert_eq!(t.nm, Some(7));
        assert_eq!(t.hp, Some(2));
        assert_eq!(t.ps.as_deref(), Some("123456"));
        assert_eq!(t.sa, Some(&b"chr1,100,+,50M,60,0;"[..]));
    }

    #[test]
    fn a_double_tag_is_skipped_and_the_tags_after_it_read() {
        let mut d = Vec::new();
        d.extend_from_slice(b"XDd");
        d.extend_from_slice(&2.5f64.to_le_bytes());
        d.extend_from_slice(b"HPc\x01");
        assert_eq!(parse_tags(&d).expect("valid tags").hp, Some(1));
    }

    #[test]
    fn corrupt_tags_are_errors_not_panics() {
        let array = |sub: u8, n: u32, data: &[u8]| {
            let mut d = vec![b'M', b'L', b'B', sub];
            d.extend_from_slice(&n.to_le_bytes());
            d.extend_from_slice(data);
            d
        };
        for d in [
            b"MLB".to_vec(),             // no subtype
            b"MLBC\x03\x00".to_vec(),    // count cut short
            array(b'C', 3, &[1, 2]),     // fewer elements than the count
            array(b'C', u32::MAX, &[]),  // a count past the record
            array(b'q', 1, &[0; 4]),     // unknown subtype
            b"PSi\x01\x02".to_vec(),     // integer cut short
            b"XDd\x00\x00\x00".to_vec(), // double cut short
        ] {
            assert!(parse_tags(&d).is_err(), "{:?}", String::from_utf8_lossy(&d));
        }
    }
}

//! The original's C full alignment: `calculate_*_full_alignment` from `*_full_alignment_dwell.c`.
//!
//! Two passes, as in C. The first reads every alignment of the region once, haplotags it, and
//! keeps for each read the bases at the "flanking" positions (the 33 positions around every
//! candidate) it covers, while counting alleles at the candidates themselves. The second builds
//! one `matrix_depth x 33 x 8` (9 with dwell time) matrix per candidate from the reads overlapping its window,
//! choosing reads with glibc `rand()` when there are more than `matrix_depth`.
//!
//! The C code visits every reference base of every read and scans all reads for every
//! candidate; this port visits only flanking positions and keeps a sliding window of reads, which
//! touches the same data in the same order.

use std::collections::HashSet;
use std::path::Path;

use snv_hts::{Alignments, Aux, Read, ReadSource as _};

use crate::Error;
use crate::haplotag::{Haplotype, haplotag_read};
use crate::inputs::PhasedSnp;
use crate::khash::KhCounter;
use crate::rand::{GlibcRand, RAND_MAX};
use crate::record::{NT16, ReadView, RefWindow, cigar_op, op, upper_base};
use crate::{DwellTime, FLANKING_BASES, POSITIONS};

/// `expand_reference_region`: reference fetched on each side of the region.
const EXPAND_REFERENCE_REGION: i64 = 2_000_000;
/// `min_haplotag_mq`: reads below this mapping quality are not haplotagged.
const MIN_HAPLOTAG_MQ: u8 = 20;
/// `SAMTOOLS_VIEW_FILTER_FLAG`: unmapped, mate unmapped, secondary, supplementary.
const FILTER_FLAG: u16 = 2316;
/// `BAM_FREVERSE`.
const FLAG_REVERSE: u16 = 0x10;

/// Bytes of the `.rodata` of the original's C library from 16 bytes before `acgt2num` to the end
/// of the 32 bytes after `num2countbase_fa` (2.0.3 conda build, `0xcf0d0..0xcf140`). The
/// original indexes both tables with `char - 'A'` without a bounds check, and `'='` (the base it
/// gives positions skipped by an `N` CIGAR operation) lands before the table, so reproducing its
/// output means reading the neighbouring bytes: zero padding, `acgt2num` itself, `HAP_TYPE` and
/// a string. Indices outside this image (bytes below `'1'` or above `'~'`) read as 0.
const RODATA: [i8; 112] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // padding
    0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, // acgt2num
    0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    100, 0, 25, -100, 0, 0, 75, 0, -50, 0, 0, 0, 0, 100, 0, 0, // num2countbase_fa
    0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    60, 30, 90, 70, 105, 108, 101, 0, 99, 104, 114, 0, 109, 118, 0, 37, // HAP_TYPE, "File", ...
    115, 32, 110, 111, 116, 32, 66, 71, 90, 70, 32, 99, 111, 109, 112, 114,
];
const ACGT2NUM_AT: i32 = 16;
const NUM2COUNTBASE_FA_AT: i32 = 48;

fn rodata(table_at: i32, c: u8) -> i8 {
    // `char` is signed, and `c - 'A'` is computed in `int`.
    let i = table_at + i32::from(c as i8) - i32::from(b'A');
    usize::try_from(i).ok().and_then(|i| RODATA.get(i)).copied().unwrap_or(0)
}

/// `num2countbase_fa[c - 'A']`: the feature value of a base (A 100, C 25, G 75, T 50, N 100,
/// `I` -50, `D` -100).
fn num2countbase(c: u8) -> i8 {
    rodata(NUM2COUNTBASE_FA_AT, c)
}

/// `acgt2num[c - 'A']`: the A/C/G/T count slot of a base; anything else counts as `A`.
fn acgt2num(c: u8) -> usize {
    // Every byte the image covers holds 0..=3 here except `num2countbase_fa`'s, which only
    // lowercase letters reach, and the original upper-cases or never sees those.
    usize::try_from(rodata(ACGT2NUM_AT, c)).unwrap_or(0).min(3)
}

/// `normalize_mq`: `(int)(x < 60 ? 100 * x / 60.0 : 100)`.
fn normalize_mq(x: u8) -> i8 {
    if x < 60 { (100.0 * f64::from(x) / 60.0) as i8 } else { 100 }
}

/// `normalize_bq`: `(int)(x < 40 ? 100 * x / 40.0 : 100)`.
fn normalize_bq(x: u8) -> i8 {
    if x < 40 { (100.0 * f64::from(x) / 40.0) as i8 } else { 100 }
}

/// `normalize_af(count / (float)depth)`. The macro does not parenthesise its argument, so the
/// numerator is `100 * count` in `size_t` divided by the depth as a `float`, not
/// `100 * (count / depth)`.
fn normalize_af(count: u64, depth: u64) -> i8 {
    let depth = depth as f32;
    if f64::from(count as f32 / depth) < 1.0 { (count.wrapping_mul(100) as f32 / depth) as i32 as i8 } else { 100 }
}

/// What a read has at one flanking position (`Pos_info`).
#[derive(Debug, Clone, Copy)]
struct PosInfo {
    /// 4-bit base code, -1 inside a deletion, 0 if never set (inside an `N` operation).
    alt_base: i8,
    bq: i8,
    /// Length of a deletion starting after this position.
    del_length: u32,
    /// Index into the job's insertions of an insertion after this position.
    ins: u32,
    /// With dwell time: the base's signal samples, plus those of an insertion after it.
    signal_length: i32,
}

const NO_INS: u32 = u32::MAX;

impl Default for PosInfo {
    fn default() -> Self {
        Self { alt_base: 0, bq: 0, del_length: 0, ins: NO_INS, signal_length: 0 }
    }
}

/// `compute_signal_lengths_from_mv_tag`: signal samples per base, in alignment order, from the
/// Dorado move table (`mv:B:*`, a stride and then one entry per sample block, non-zero where a
/// new base starts). Blocks before the first base are dropped and counting stops after the
/// last base. `None` without a usable tag.
fn signal_lengths(record: &Read) -> Option<Vec<i32>> {
    let read_len = record.sequence_len();
    // `bam_auxB2i` turns any array element into an integer (floats truncate).
    let moves: Vec<i64> = match snv_hts::aux(record, *b"mv")? {
        Aux::Ints(a) => a,
        Aux::Floats(a) => a.iter().map(|&x| x as i64).collect(),
        _ => return None,
    };
    if moves.len() <= 1 || read_len == 0 {
        return None;
    }
    let mut signals = vec![0i32; read_len];
    let mut base: Option<usize> = None;
    for &movement in &moves[1..] {
        if movement != 0 {
            let next = base.map_or(0, |b| b + 1);
            if next >= read_len {
                break;
            }
            base = Some(next);
        }
        let Some(b) = base else { continue };
        signals[b] += 1;
    }
    if record.flags() & FLAG_REVERSE != 0 {
        signals.reverse();
    }
    Some(signals)
}

/// A read kept for the matrices (`struct Read` after the first pass).
#[derive(Debug)]
struct ReadFeatures {
    start: u64,
    end: u64,
    haplotype: Haplotype,
    strand: i8,
    mq: i8,
    /// Index of the first flanking position the read covers.
    flanking_start: usize,
    /// One entry per covered flanking position, from `flanking_start`.
    pos_info: Vec<PosInfo>,
}

/// Allele counts at one candidate (`Pos_alt_info`).
#[derive(Debug, Default)]
struct PosAltInfo {
    ins_counter: KhCounter<Box<[u8]>>,
    del_counter: KhCounter<i32>,
    acgt_count: [u64; 4],
    depth: u64,
}

/// The centre-position allele of a matrix row (`Alt_info`).
#[derive(Debug, Clone, Copy, Default)]
struct RowAlt {
    ins: Option<u32>,
    alt_base: u8,
    del_length: u32,
    has_alt_info: bool,
}

/// What the C full-alignment function (`calculate_*_full_alignment`) is called with.
#[derive(Debug)]
pub(crate) struct Job<'a> {
    pub(crate) bam: &'a Path,
    pub(crate) reference: &'a Path,
    pub(crate) contig: &'a str,
    /// 0-based start and end of the region, as `hts_parse_reg` returns them.
    pub(crate) beg: i64,
    pub(crate) end: i64,
    /// 0-based candidates, ascending and distinct.
    pub(crate) candidates: &'a [u64],
    /// `None` without haplotagging.
    pub(crate) snps: Option<&'a [PhasedSnp]>,
    pub(crate) min_mapq: u64,
    pub(crate) matrix_depth: usize,
    pub(crate) max_indel_length: u64,
    pub(crate) dwell_time: DwellTime,
    pub(crate) canonical_read_order: bool,
}

/// `fa_data`'s alt-info C strings (the matrices go to the caller's sink as they are built).
#[derive(Debug)]
pub(crate) struct FaData {
    pub(crate) alt_info: Vec<Vec<u8>>,
}

/// The flanking positions of all candidates, sorted, with the candidate (if any) at each.
#[derive(Debug)]
struct Flanking {
    positions: Vec<u64>,
    candidate: Vec<Option<u32>>,
}

impl Flanking {
    /// The C loop inserts `max(0, candidate - 16 + j)` for `j < 33` into a hash in candidate
    /// order; for sorted candidates every new key is larger than all earlier ones, so the
    /// insertion index is the rank.
    fn new(candidates: &[u64]) -> Self {
        let mut positions: Vec<u64> = Vec::new();
        for &c in candidates {
            for j in 0..POSITIONS as u64 {
                // `int result = candidate - flanking_base_num + j`.
                let result = c.wrapping_sub(FLANKING_BASES).wrapping_add(j) as i32;
                let key = u64::try_from(result).unwrap_or(0);
                if positions.last().is_none_or(|&last| key > last) {
                    positions.push(key);
                }
            }
        }
        let mut candidate = vec![None; positions.len()];
        for (i, c) in candidates.iter().enumerate() {
            let k = positions.binary_search(c).expect("every candidate is one of its own flanking positions");
            candidate[k] = Some(u32::try_from(i).expect("fewer than 2^32 candidates"));
        }
        Self { positions, candidate }
    }

    fn len(&self) -> usize {
        self.positions.len()
    }

    /// Index of `pos` among `positions[range]`'s values, as an absolute index.
    fn find(&self, lo: usize, hi: usize, pos: u64) -> Option<usize> {
        self.positions[lo..hi].binary_search(&pos).ok().map(|k| lo + k)
    }

    /// Absolute indices of the positions in `from..to` within `positions[lo..hi]`.
    fn span(&self, lo: usize, hi: usize, from: u64, to: u64) -> std::ops::Range<usize> {
        let slice = &self.positions[lo..hi];
        let a = slice.partition_point(|&p| p < from);
        let b = slice.partition_point(|&p| p < to);
        lo + a..lo + b.max(a)
    }
}

/// Runs the C full-alignment function for one job, handing the matrices to `sink` in order, `batch`
/// candidates at a time (the last batch may be shorter), as soon as each batch is built.
pub(crate) fn calculate_full_alignment(
    job: &Job<'_>,
    batch: usize,
    sink: &mut dyn FnMut(&[i8]),
) -> Result<FaData, Error> {
    let candidates = job.candidates;
    // `max(int, int)`: `start - expand_reference_region` is computed in `size_t` and converted
    // back to `int`, which gives the plain difference.
    let ref_start = (job.beg - EXPAND_REFERENCE_REGION).max(0);
    let ref_end = job.end + EXPAND_REFERENCE_REGION;
    let fasta = snv_pileup::Fasta::open(job.reference).map_err(|source| Error::Reference { source })?;
    let to_c_int = |v: i64| i32::try_from(v).map_err(|_| Error::Position { position: v });
    let seq = fasta
        .fetch(job.contig, to_c_int(ref_start)?, to_c_int(ref_end).unwrap_or(i32::MAX))
        .map_err(|source| Error::Reference { source })?;
    let reference = RefWindow::new(seq, ref_start as u64);

    let alignments_error = |e: snv_hts::Error| Error::Alignments { path: job.bam.to_owned(), source: Box::new(e) };
    let mut alignments = Alignments::open(job.bam, job.reference)
        .map_err(alignments_error)?
        .with_canonical_read_order(job.canonical_read_order);
    if alignments.tid(job.contig).is_none() {
        return Err(Error::Contig { contig: job.contig.to_owned() });
    }
    let mut records = alignments.query_interval(job.contig, job.beg, job.end).map_err(alignments_error)?;

    let flanking = Flanking::new(candidates);
    let snps = job.snps.unwrap_or(&[]);
    let mut pos_alt_info: Vec<PosAltInfo> = candidates.iter().map(|_| PosAltInfo::default()).collect();
    let mut insertions: Vec<Box<[u8]>> = Vec::new();
    let mut reads: Vec<ReadFeatures> = Vec::new();
    let mut read_names: HashSet<Box<[u8]>> = HashSet::new();
    let mut variant_current_pos = 0;
    let mut candidate_current_index = 0;
    let mut builder = Builder {
        candidates,
        flanking: &flanking,
        reference: &reference,
        matrices: Matrices::new(job.matrix_depth, job.dwell_time, batch.clamp(1, candidates.len().max(1))),
        sink,
        window: ReadWindow::default(),
        rng: GlibcRand::new(),
        alt_info: Vec::with_capacity(candidates.len()),
        max_indel_length: job.max_indel_length,
        // Candidates out of order are built once all reads are in, as the C code does.
        streaming: candidates.is_sorted(),
    };

    while let Some(record) = records.next_read().map_err(alignments_error)? {
        if record.flags() & FILTER_FLAG != 0 || u64::from(record.mapping_quality()) < job.min_mapq {
            continue;
        }
        // Later records with an already seen name are skipped.
        if !read_names.insert(record.name().into()) {
            continue;
        }
        let pos = record.position();
        // The records come by position, so the windows that end at or before this one's start
        // have all their reads.
        if let Ok(start) = u64::try_from(pos) {
            builder.build_until(Some(start), &mut reads, &pos_alt_info, &insertions);
        }
        let read = ReadView {
            // `size_t read_start = (int32_t)pos`.
            start: pos as i32 as i64 as u64,
            cigar: record.cigar(),
            seq: snv_hts::packed_seq(&record),
            qual: record.quality_scores(),
            len: record.sequence_len(),
        };
        while variant_current_pos < snps.len() && i64::from(snps[variant_current_pos].position) < pos {
            variant_current_pos += 1;
        }
        while candidate_current_index < flanking.len() && flanking.positions[candidate_current_index] < read.start {
            candidate_current_index += 1;
        }
        let read_end = read.end();
        // `get_overlap_candidate_num`: flanking positions from here that lie before the end.
        let overlap = flanking.positions[candidate_current_index..].partition_point(|&p| p < read_end);
        if overlap == 0 {
            continue;
        }
        let haplotype = if job.snps.is_some() && record.mapping_quality() >= MIN_HAPLOTAG_MQ {
            haplotag_read(snps, variant_current_pos, &read, &reference)
        } else {
            Haplotype::Unphased
        };
        let signals = if job.dwell_time == DwellTime::On { signal_lengths(&record) } else { None };
        let signal = |q: u64| signals.as_ref().and_then(|s| s.get(usize::try_from(q).ok()?)).copied();
        let lo = candidate_current_index;
        let hi = lo + overlap;
        let mut pos_info = vec![PosInfo::default(); overlap];
        let mut ref_pos = read.start;
        let mut query_pos = 0u64;
        for &c in read.cigar {
            let (o, length) = cigar_op(c);
            if op::is_match(o) {
                for k in flanking.span(lo, hi, ref_pos, ref_pos.wrapping_add(length)) {
                    let q = query_pos + (flanking.positions[k] - ref_pos);
                    let info = &mut pos_info[k - lo];
                    let code = read.base_code(q);
                    info.alt_base = code as i8;
                    info.bq = normalize_bq(read.qual(q));
                    if let Some(length) = signal(q) {
                        info.signal_length = length;
                    }
                    if let Some(ci) = flanking.candidate[k] {
                        let alt = &mut pos_alt_info[ci as usize];
                        alt.acgt_count[acgt2num(NT16[usize::from(code)])] += 1;
                        alt.depth += 1;
                    }
                }
                query_pos = query_pos.wrapping_add(length);
                ref_pos = ref_pos.wrapping_add(length);
            } else if o == op::DEL {
                if let Some(k) = flanking.find(lo, hi, ref_pos.wrapping_sub(1)) {
                    pos_info[k - lo].del_length = length as u32;
                    if let Some(ci) = flanking.candidate[k] {
                        pos_alt_info[ci as usize].del_counter.add(length as i32, 1);
                    }
                }
                for k in flanking.span(lo, hi, ref_pos, ref_pos.wrapping_add(length)) {
                    pos_info[k - lo].alt_base = -1;
                    if let Some(ci) = flanking.candidate[k] {
                        pos_alt_info[ci as usize].depth += 1;
                    }
                }
                ref_pos = ref_pos.wrapping_add(length);
            } else if o == op::INS {
                if let Some(k) = flanking.find(lo, hi, ref_pos.wrapping_sub(1)) {
                    let bases: Box<[u8]> =
                        (0..length).map(|i| NT16[usize::from(read.base_code(query_pos + i))]).collect();
                    if let Some(ci) = flanking.candidate[k] {
                        pos_alt_info[ci as usize].ins_counter.add(bases.clone(), 1);
                    }
                    // A second insertion at the same position replaces the first, but its
                    // signal adds to the first's.
                    pos_info[k - lo].ins = u32::try_from(insertions.len()).expect("fewer than 2^32 insertions");
                    insertions.push(bases);
                    if signals.is_some() {
                        let inserted = (0..length).filter_map(|i| signal(query_pos + i)).fold(0i32, i32::wrapping_add);
                        pos_info[k - lo].signal_length = pos_info[k - lo].signal_length.wrapping_add(inserted);
                    }
                }
                query_pos = query_pos.wrapping_add(length);
            } else if o == op::REF_SKIP {
                ref_pos = ref_pos.wrapping_add(length);
            } else if o == op::SOFT_CLIP {
                query_pos = query_pos.wrapping_add(length);
            }
        }
        reads.push(ReadFeatures {
            start: read.start,
            end: read_end,
            haplotype,
            // `normalize_strand(is_reverse)`: 50 for reverse reads, 100 for forward ones.
            strand: if record.flags() & FLAG_REVERSE != 0 { 50 } else { 100 },
            mq: normalize_mq(record.mapping_quality()),
            flanking_start: lo,
            pos_info,
        });
    }

    builder.build_until(None, &mut reads, &pos_alt_info, &insertions);
    builder.flush();
    Ok(FaData { alt_info: builder.alt_info })
}

/// Builds the candidates' matrices in order while the reads come in. A candidate's window takes
/// the reads that start before its end, so it can be built once a read starts at or past that
/// end, and the reads that ended before every later window are dropped: a job then holds the
/// reads around its current position rather than all of its region's.
struct Builder<'a> {
    candidates: &'a [u64],
    flanking: &'a Flanking,
    reference: &'a RefWindow,
    /// The batch being built.
    matrices: Matrices,
    /// Where full batches go.
    sink: &'a mut dyn FnMut(&[i8]),
    window: ReadWindow,
    rng: GlibcRand,
    /// One entry per candidate built so far.
    alt_info: Vec<Vec<u8>>,
    max_indel_length: u64,
    /// Whether candidates may be built before all reads are in (they are sorted).
    streaming: bool,
}

impl Builder<'_> {
    /// Reads dropped at a time, once no window can take them.
    const DROP_BATCH: usize = 4096;

    /// Hands the matrices built since the last batch to the sink and starts the next batch.
    fn flush(&mut self) {
        let built = self.alt_info.len() - self.matrices.first;
        if built == 0 {
            return;
        }
        let m = &mut self.matrices;
        let len = built * m.depth * POSITIONS * m.dwell_time.channels();
        (self.sink)(&m.data[..len]);
        m.data.fill(0);
        m.first += built;
    }

    /// Builds the next candidates, in order, whose windows end at or before `limit`, or all the
    /// rest with `None`; `reads` holds the reads from the window's first on.
    fn build_until(
        &mut self,
        limit: Option<u64>,
        reads: &mut Vec<ReadFeatures>,
        pos_alt_info: &[PosAltInfo],
        insertions: &[Box<[u8]>],
    ) {
        if limit.is_some() && !self.streaming {
            return;
        }
        let depth = self.matrices.depth;
        while let Some(&candidate) = self.candidates.get(self.alt_info.len()) {
            if limit.is_some_and(|l| candidate + FLANKING_BASES + 1 > l) {
                break;
            }
            let i = self.alt_info.len();
            let mut rows = vec![RowAlt::default(); depth];
            // `start_pos = candidate - 16` wraps below 16, and then no read overlaps.
            if candidate >= FLANKING_BASES {
                let overlapping =
                    self.window.advance(reads, candidate - FLANKING_BASES, candidate + FLANKING_BASES + 1);
                let selected = select_reads(overlapping, reads, depth, &mut self.rng);
                let base = self
                    .flanking
                    .find(0, self.flanking.len(), candidate - FLANKING_BASES)
                    .expect("window start is flanking");
                for (d, read) in selected.iter().enumerate() {
                    let Some(r) = read else { continue };
                    let read = &reads[*r as usize];
                    rows[d] = self.matrices.fill_row(
                        i,
                        d,
                        read,
                        base,
                        candidate - FLANKING_BASES,
                        self.reference,
                        insertions,
                    );
                }
            }
            self.matrices.fill_af(i, &rows, &pos_alt_info[i], insertions);
            self.alt_info.push(alt_info_string(candidate, &pos_alt_info[i], self.reference, self.max_indel_length));
            if self.alt_info.len() - self.matrices.first == self.matrices.capacity {
                self.flush();
            }
        }
        // Reads before the window's first are out of every later window (windows only move
        // right); read indices only order reads, so they can shift down.
        let first = self.window.active.first().map_or(self.window.next, |&r| r as usize);
        if first >= Self::DROP_BATCH {
            reads.drain(..first);
            self.window.next -= first;
            let shift = u32::try_from(first).expect("fewer than 2^32 reads");
            for r in &mut self.window.active {
                *r -= shift;
            }
        }
    }
}

/// Reads overlapping the current candidate window, in read order. Windows only move right,
/// so a read that ends before one window never overlaps a later one.
#[derive(Debug, Default)]
struct ReadWindow {
    next: usize,
    active: Vec<u32>,
}

impl ReadWindow {
    /// The reads with `start < end_pos` and `end > start_pos`, as the C scan over all reads
    /// (which stops at the first read starting at or after `end_pos`) finds them.
    fn advance(&mut self, reads: &[ReadFeatures], start_pos: u64, end_pos: u64) -> &[u32] {
        while self.next < reads.len() && reads[self.next].start < end_pos {
            self.active.push(u32::try_from(self.next).expect("fewer than 2^32 reads"));
            self.next += 1;
        }
        self.active.retain(|&r| reads[r as usize].end > start_pos);
        &self.active
    }
}

/// `sort_read_name_by_haplotype`: the read of each matrix row (`None` for padding).
fn select_reads(overlapping: &[u32], reads: &[ReadFeatures], depth: usize, rng: &mut GlibcRand) -> Vec<Option<u32>> {
    let mut haps: Vec<(Haplotype, u32)> = overlapping.iter().map(|&r| (reads[r as usize].haplotype, r)).collect();
    let n = haps.len();
    if n > depth {
        // Shuffle with the process-wide `rand()` stream.
        for i in 0..n - 1 {
            let span = (n - i) as u64;
            let j = i + (u64::from(rng.next_u31()) / (u64::from(RAND_MAX) / span + 1)) as usize;
            haps.swap(i, j);
        }
    }
    let read_num = n.min(depth);
    // `hap_cmp`: haplotype, then read index; the keys are distinct, so any sort agrees with qsort.
    haps[..read_num].sort_unstable();
    let mut rows = vec![None; depth];
    let prefix = if n < depth { (depth - read_num) >> 1 } else { 0 };
    for (k, &(_, r)) in haps[..read_num].iter().enumerate() {
        rows[prefix + k] = Some(r);
    }
    rows
}

/// The output matrices (`int8_t *matrix`) of a batch of candidates: `capacity` of them from
/// candidate `first` on.
#[derive(Debug)]
struct Matrices {
    depth: usize,
    dwell_time: DwellTime,
    first: usize,
    capacity: usize,
    data: Vec<i8>,
}

impl Matrices {
    fn new(depth: usize, dwell_time: DwellTime, capacity: usize) -> Self {
        let data = vec![0; capacity * depth * POSITIONS * dwell_time.channels()];
        Self { depth, dwell_time, first: 0, capacity, data }
    }

    fn at(&mut self, i: usize, d: usize, p: usize, channel: usize) -> &mut i8 {
        let channels = self.dwell_time.channels();
        &mut self.data[(((i - self.first) * self.depth + d) * POSITIONS + p) * channels + channel]
    }

    /// Writes row `d` of candidate `i` from `read`; `base` is the flanking index of
    /// `start_pos`. Returns the row's centre allele.
    #[expect(clippy::too_many_arguments, reason = "mirrors the C loop body's inputs")]
    fn fill_row(
        &mut self,
        i: usize,
        d: usize,
        read: &ReadFeatures,
        base: usize,
        start_pos: u64,
        reference: &RefWindow,
        insertions: &[Box<[u8]>],
    ) -> RowAlt {
        let mut row = RowAlt::default();
        let hap_v = read.haplotype.feature();
        for p in 0..POSITIONS {
            // Positions the read does not cover are skipped (C reads `pos_info` out of bounds
            // first, but skips them either way).
            let Some(info) = (base + p).checked_sub(read.flanking_start).and_then(|o| read.pos_info.get(o)) else {
                continue;
            };
            if info.alt_base < 0 {
                continue;
            }
            let is_center = p == FLANKING_BASES as usize;
            let ref_base = upper_base(reference.at(start_pos + p as u64));
            let ref_v = num2countbase(ref_base);
            let alt_base = NT16[info.alt_base as usize];
            let ins = (info.ins != NO_INS).then(|| &insertions[info.ins as usize]).filter(|b| !b.is_empty());
            let alt_v = if let Some(bases) = ins {
                if p < POSITIONS - 1 {
                    for (k, &b) in bases.iter().take(POSITIONS - p).enumerate() {
                        *self.at(i, d, p + k, 6) = num2countbase(b);
                    }
                }
                if is_center {
                    row = RowAlt { ins: Some(info.ins), alt_base, has_alt_info: true, ..row };
                }
                num2countbase(b'I')
            } else if info.del_length > 0 {
                if is_center {
                    row = RowAlt { del_length: info.del_length, has_alt_info: true, ..row };
                }
                num2countbase(b'D')
            } else if ref_base != alt_base {
                if is_center {
                    row = RowAlt { alt_base, has_alt_info: true, ..row };
                }
                num2countbase(alt_base)
            } else {
                0
            };
            *self.at(i, d, p, 0) = ref_v;
            *self.at(i, d, p, 1) = alt_v;
            *self.at(i, d, p, 2) = read.strand;
            *self.at(i, d, p, 3) = read.mq;
            *self.at(i, d, p, 4) = info.bq;
            *self.at(i, d, p, 7) = hap_v;
            if self.dwell_time == DwellTime::On {
                // `(int8_t)raw_signal`.
                *self.at(i, d, p, 8) = info.signal_length as i8;
            }
        }
        row
    }

    /// The candidate-proportion channel: each row's centre-allele frequency on every position
    /// with a reference feature.
    fn fill_af(&mut self, i: usize, rows: &[RowAlt], counts: &PosAltInfo, insertions: &[Box<[u8]>]) {
        let depth = counts.depth;
        for (j, row) in rows.iter().enumerate() {
            if !row.has_alt_info {
                continue;
            }
            let af = if let Some(ins) = row.ins {
                // `kh_counter_val` is 0 for a missing key.
                let count = counts.ins_counter.get(&insertions[ins as usize]).unwrap_or(0) as i64 as u64;
                if count > 0 { normalize_af(count, depth) } else { 0 }
            } else if row.del_length > 0 {
                // `kh_int_counter_val` is -1 for a missing key, which `size_t` makes huge.
                let count = counts.del_counter.get(&(row.del_length as i32)).unwrap_or(-1) as i64 as u64;
                if count > 0 { normalize_af(count, depth) } else { 0 }
            } else if row.alt_base != 0 {
                normalize_af(counts.acgt_count[acgt2num(row.alt_base)], depth)
            } else {
                0
            };
            if af > 0 {
                for p in 0..POSITIONS {
                    if *self.at(i, j, p, 0) != 0 {
                        *self.at(i, j, p, 5) = af;
                    }
                }
            }
        }
    }
}

/// The alt-info C string of a candidate:
/// `"{pos}-{depth}-{ref}-" + "X{base} {count} "... + "I{ref}{bases} {count} "... +
/// "D{deleted ref} {count} "... + "R{ref} {count} "`, insertions and deletions in khash order.
fn alt_info_string(candidate: u64, counts: &PosAltInfo, reference: &RefWindow, max_indel_length: u64) -> Vec<u8> {
    use std::io::Write as _;

    let center_ref_base = upper_base(reference.at(candidate));
    let ref_slot = acgt2num(center_ref_base);
    // `int ref_count` from a `size_t`.
    let mut ref_count = counts.acgt_count[ref_slot] as i32;
    let mut s: Vec<u8> = Vec::with_capacity(64);
    let w = "writing to a Vec cannot fail";
    write!(s, "{}-{}-", candidate + 1, counts.depth).expect(w);
    s.extend_from_slice(&[center_ref_base, b'-']);
    for (j, &count) in counts.acgt_count.iter().enumerate() {
        if j != ref_slot && count > 0 {
            write!(s, "X{} {count} ", char::from(b"ACGT"[j])).expect(w);
        }
    }
    for (key, val) in counts.ins_counter.iter() {
        ref_count = ref_count.wrapping_sub(val);
        if key.len() as u64 <= max_indel_length {
            s.extend_from_slice(&[b'I', center_ref_base]);
            s.extend_from_slice(key);
            write!(s, " {val} ").expect(w);
        }
    }
    for (&key, val) in counts.del_counter.iter() {
        ref_count = ref_count.wrapping_sub(val);
        // `key <= (int)max_indel_length`.
        if i64::from(key) <= i64::from(max_indel_length as i32) {
            s.push(b'D');
            s.extend_from_slice(reference.c_str(candidate + 1, u64::try_from(key).unwrap_or(0)));
            write!(s, " {val} ").expect(w);
        }
    }
    if ref_count > 0 {
        s.extend_from_slice(&[b'R', center_ref_base]);
        write!(s, " {ref_count} ").expect(w);
    }
    // A C string ends at its first NUL (a NUL reference byte would cut it there).
    if let Some(nul) = s.iter().position(|&b| b == 0) {
        s.truncate(nul);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_lookups() {
        let fa: Vec<i8> = b"ACGTNIDB".iter().map(|&c| num2countbase(c)).collect();
        assert_eq!(fa, [100, 25, 75, 50, 100, -50, -100, 0]);
        let slots: Vec<usize> = b"ACGTN".iter().map(|&c| acgt2num(c)).collect();
        assert_eq!(slots, [0, 1, 2, 3, 0]);
        // '=' indexes before both tables, where the original's C library has zeros.
        assert_eq!(num2countbase(b'='), 0);
        assert_eq!(acgt2num(b'='), 0);
        // 'a' - 'A' = 32 runs past `num2countbase_fa` into `HAP_TYPE`.
        assert_eq!(num2countbase(b'a'), 60);
        assert_eq!(num2countbase(0), 0);
    }

    #[test]
    fn normalizations() {
        assert_eq!(normalize_mq(59), 98);
        assert_eq!(normalize_mq(60), 100);
        assert_eq!(normalize_mq(255), 100);
        assert_eq!(normalize_bq(13), 32);
        assert_eq!(normalize_bq(41), 100);
        // 100 * 1 / 3.0f = 33.333 -> 33; 100 * 2 / 3.0f -> 66.
        assert_eq!(normalize_af(1, 3), 33);
        assert_eq!(normalize_af(2, 3), 66);
        assert_eq!(normalize_af(3, 3), 100);
        assert_eq!(normalize_af(u64::MAX, 5), 100);
        assert_eq!(normalize_af(0, 0), 100);
        // The unparenthesised macro: (100 * 7) / 9.0f = 77.77 -> 77, whereas
        // 100 * (7 / 9.0f) = 77.77777 too; 29/100: 2900/100 = 29 vs 100 * 0.29f = 28.999998.
        assert_eq!(normalize_af(29, 100), 29);
    }

    fn read(start: u64, end: u64, haplotype: Haplotype) -> ReadFeatures {
        ReadFeatures { start, end, haplotype, strand: 100, mq: 100, flanking_start: 0, pos_info: Vec::new() }
    }

    #[test]
    fn window_matches_full_scan() {
        let reads: Vec<ReadFeatures> =
            (0..50).map(|i| read(i * 7, i * 7 + 3 + (i * 37) % 60, Haplotype::Unphased)).collect();
        let mut window = ReadWindow::default();
        for c in (16..400).step_by(5) {
            let (s, e) = (c - 16, c + 17);
            let mut want = Vec::new();
            for (j, r) in reads.iter().enumerate() {
                if r.start >= e {
                    break;
                }
                if r.end > s {
                    want.push(j as u32);
                }
            }
            assert_eq!(window.advance(&reads, s, e), want, "candidate {c}");
        }
    }

    #[test]
    fn selection_pads_and_sorts() {
        let reads: Vec<ReadFeatures> =
            [Haplotype::Second, Haplotype::Unphased, Haplotype::First].into_iter().map(|h| read(0, 10, h)).collect();
        let mut rng = GlibcRand::new();
        let rows = select_reads(&[0, 1, 2], &reads, 6, &mut rng);
        assert_eq!(rows, [None, Some(1), Some(2), Some(0), None, None]);
        // No `rand()` was drawn.
        assert_eq!(rng.next_u31(), 1_804_289_383);
        let rows = select_reads(&[0, 1, 2], &reads, 3, &mut rng);
        assert_eq!(rows, [Some(1), Some(2), Some(0)]);
    }

    #[test]
    fn downsampling_shuffles_with_rand() {
        let reads: Vec<ReadFeatures> = (0..10).map(|_| read(0, 10, Haplotype::Unphased)).collect();
        let ids: Vec<u32> = (0..10).collect();
        let mut rng = GlibcRand::new();
        let rows = select_reads(&ids, &reads, 4, &mut rng);
        // Fisher-Yates as in C with glibc's first nine draws, then the first 4 sorted.
        let mut expect = ids.clone();
        let mut check = GlibcRand::new();
        for i in 0..9 {
            let j = i + (check.next_u31() as usize) / (RAND_MAX as usize / (10 - i) + 1);
            expect.swap(i, j);
        }
        let mut first: Vec<u32> = expect[..4].to_vec();
        first.sort_unstable();
        assert_eq!(rows, first.into_iter().map(Some).collect::<Vec<_>>());
    }

    #[test]
    fn flanking_ranks() {
        let f = Flanking::new(&[5, 40, 50]);
        // 5 -> clamped 0..=21; 40 -> 24..=56; 50 -> up to 66.
        assert_eq!(f.positions[0], 0);
        assert_eq!(f.positions[21], 21);
        assert_eq!(f.positions[22], 24);
        assert_eq!(*f.positions.last().expect("non-empty"), 66);
        assert_eq!(f.candidate[5], Some(0));
        assert_eq!(f.find(0, f.len(), 22), None);
        assert_eq!(f.span(0, f.len(), 20, 26), 20..24);
    }
}

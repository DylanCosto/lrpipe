//! Port of the original's C pileup (`calculate_*_pileup` in `src/*_pileup.c`).
//!
//! For each reference position in a region this builds the 18 pileup features and decides
//! whether the position is a variant candidate, writing the same "alt info" string the C code
//! writes. Where the C code relies on unsigned wrap-around or out-of-bounds access, the
//! behaviour that actually results is reproduced and marked `C quirk`.

use std::fmt::Write as _;

use crate::Error;
#[cfg(test)]
use crate::hts::Column;
use crate::hts::{AlignmentFile, Fasta, Read};
use crate::khash::StrCounter;

/// Features per position (`featlen*` in the C pileup).
pub const FEATURES: usize = 18;

const FWD_INS_ALL: usize = 4;
const FWD_INS_BEST: usize = 5;
const FWD_DEL_ALL: usize = 6;
const FWD_DEL_BEST: usize = 7;
const FWD_DEL: usize = 8;
const REV_INS_ALL: usize = 13;
const REV_INS_BEST: usize = 14;
const REV_DEL_ALL: usize = 15;
const REV_DEL_BEST: usize = 16;
const REV_DEL: usize = 17;
/// Offset of the reverse-strand features (`reverse_pos_start`).
const REVERSE_OFFSET: usize = 9;
/// Reference bases fetched on each side of the region (`mpileup_expand_reference_region`).
const REFERENCE_PADDING: i64 = 1000;
/// Contiguous covered positions required before a candidate (`pileup_flanking_base_num`).
const FLANKING_BASES: usize = 16;

/// Feature column of each `seq_nt16` code, plus 16 for the reverse strand
/// (`num2countbase*` in the C pileup). `None` is the C value -1.
const BASE_FEATURE: [Option<usize>; 32] = {
    let mut t = [None; 32];
    (t[1], t[2], t[4], t[8]) = (Some(0), Some(1), Some(2), Some(3));
    (t[17], t[18], t[20], t[24]) = (Some(9), Some(10), Some(11), Some(12));
    t
};

/// `plp_bases_*` in the C pileup: the base letter of feature columns 0..4.
const FEATURE_BASES: [u8; 4] = *b"ACGT";

/// `seq_nt16_str`.
const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Thresholds for candidate selection (the scalar arguments of the C pileup function).
#[derive(Debug, Clone, PartialEq)]
pub struct CandidateParams {
    /// Minimum depth (`min_depth`).
    pub min_depth: usize,
    /// Minimum SNP allele frequency (`min_snp_af`).
    pub min_snp_af: f32,
    /// Minimum indel allele frequency (`min_indel_af`).
    pub min_indel_af: f32,
    /// Minimum mapping quality of reads (`min_mq`).
    pub min_mapq: u8,
    /// Longest indel written into the alt info (`max_indel_length`).
    pub max_indel_length: usize,
    /// Which alleles make a position a candidate.
    pub alleles: CandidateAlleles,
    /// Whether candidates need 16 covered bases before them.
    pub flanks: FlankRequirement,
}

/// Which alleles can make a position a candidate (`call_snp_only`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateAlleles {
    /// SNPs and indels.
    SnpsAndIndels,
    /// Only SNPs, by allele frequency.
    SnpsOnly,
}

/// Whether a candidate needs fully covered flanks (`call_ht`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlankRequirement {
    /// At least 16 contiguous covered positions before the candidate (the default).
    Contiguous,
    /// No requirement (`--enable_variant_calling_at_sequence_head_and_tail`).
    None,
}

/// Read counts of one position for gVCF output (`pos_ref_count`, `pos_total_count`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SiteDepth {
    /// Reads with the reference base.
    pub reference: u64,
    /// The reference reads plus the deletions, the insertions and `all_alt_count`, which sums
    /// the running maxima of the other bases' counts rather than all of them (`C quirk`).
    pub total: u64,
}

/// Pileup features and candidates of one region (`plp_data`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PileupCounts {
    /// One row of [`FEATURES`] counts per covered position. The reference base's columns hold
    /// minus the strand's base total, as in the original.
    pub rows: Vec<[i32; FEATURES]>,
    /// 0-based reference position of each row (`major`).
    pub positions: Vec<i64>,
    /// gVCF read counts of each row.
    pub depths: Vec<SiteDepth>,
    /// Alt info of each candidate, `"{pos1}-{depth}-{ref}-{alleles}"`, with the trailing space
    /// that Python strips still present.
    pub candidates: Vec<String>,
    /// The region as htslib parsed it: 0-based start and exclusive end. The C code's gVCF
    /// buffers hold `end - start` positions from `start`.
    pub bounds: (i64, i64),
}

/// Region bounds as `hts_parse_reg` reads them: 0-based start, exclusive end.
fn parse_region_bounds(region: &str) -> Result<(String, i64, i64), Error> {
    let bad = || Error::Region { region: region.to_owned() };
    let (contig, bounds) = region.rsplit_once(':').ok_or_else(bad)?;
    let (start, end) = bounds.split_once('-').ok_or_else(bad)?;
    let parse = |s: &str| s.replace(',', "").parse::<i64>().map_err(|_| bad());
    let start = (parse(start)? - 1).max(0);
    let end = parse(end)?;
    if start >= end {
        return Err(bad());
    }
    Ok((contig.to_owned(), start, end))
}

/// The reference around a region, padded by 1000 bases on each side.
struct ReferenceWindow {
    seq: Vec<u8>,
    seq_start: i64,
}

impl ReferenceWindow {
    fn fetch(fasta: &Fasta, contig: &str, start: i64, end: i64) -> Result<Self, Error> {
        // `max` in medaka_common.h is `int max(int, int)`, so these clamp at 0 as written.
        let beg = (start - REFERENCE_PADDING).max(0);
        let end = (end + REFERENCE_PADDING).max(0);
        Ok(Self { seq: fasta.fetch(contig, beg as i32, end as i32)?, seq_start: beg })
    }

    fn offset(&self, pos: i64) -> i64 {
        pos - self.seq_start
    }

    /// Byte at `offset`; past the end the C code reads the NUL terminator.
    fn byte(&self, offset: i64) -> u8 {
        usize::try_from(offset).ok().and_then(|i| self.seq.get(i)).copied().unwrap_or(0)
    }

    /// `%.*s` of `len` bases starting at `offset`, stopping at the end of the buffer.
    fn slice(&self, offset: i64, len: usize) -> &[u8] {
        let start = usize::try_from(offset).unwrap_or(usize::MAX).min(self.seq.len());
        &self.seq[start..(start + len).min(self.seq.len())]
    }
}

/// `upper_base`.
fn upper_base(c: u8) -> u8 {
    c.to_ascii_uppercase()
}

/// `base2index[ref_base - 'A']`: A, C, G, T map to 0..4 and every other byte to 0.
fn reference_feature(ref_base: u8) -> usize {
    match ref_base {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// Computes pileup features and candidates for `region` (a samtools-style region string,
/// 1-based inclusive, passed through to htslib unchanged).
///
/// The columns are what htslib's `bam_plp` engine gives the C code, built read by read: each
/// read's CIGAR is walked once (as `resolve_cigar2` walks it a position at a time) and its bases,
/// deletions and insertions are added to the positions it covers. Reads come in file order, so
/// each position sees them in the order `bam_plp` lists them.
///
/// # Errors
///
/// Fails if the region cannot be parsed or a file cannot be read.
pub fn pileup_counts(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    region: &str,
    params: &CandidateParams,
) -> Result<PileupCounts, Error> {
    let (contig, start, end) = parse_region_bounds(region)?;
    let tid = bam.tid(&contig);
    let reference = ReferenceWindow::fetch(fasta, &contig, start, end)?;
    let mut reads = bam.reads(region, params.min_mapq)?;
    let mut next = || -> Result<Option<ActiveRead>, Error> {
        while let Some(read) = reads.next_read()? {
            if Some(read.tid()) == tid {
                return Ok(Some(ActiveRead::new(&read)));
            }
        }
        Ok(None)
    };
    // Positions are counted a tile at a time, with every read over the tile, so the counts stay
    // in cache; each read keeps its place in its CIGAR from one tile to the next.
    let mut columns = Columns::new((start, end));
    let mut tile = TileCounts::default();
    let mut active: Vec<ActiveRead> = Vec::new();
    let mut pending = next()?;
    let mut tile_start = start;
    while tile_start < end {
        if active.is_empty() {
            let Some(read) = &pending else { break };
            tile_start = tile_start.max(read.x);
            if tile_start >= end {
                break;
            }
        }
        let tile_end = (tile_start + TILE_LENGTH).min(end);
        while let Some(read) = pending.take_if(|read| read.x < tile_end) {
            active.push(read);
            pending = next()?;
        }
        tile.reset(tile_start, tile_end);
        for read in &mut active {
            tile.add(read);
        }
        active.retain(|read| read.end > tile_end);
        columns.add_tile(&mut tile, &reference, params);
        tile_start = tile_end;
    }
    Ok(columns.out)
}

/// Positions counted at once by [`pileup_counts`]: their counts fit in a core's cache.
const TILE_LENGTH: i64 = 2048;

/// CIGAR operation codes (`BAM_C*`).
const CIGAR_MATCH: u32 = 0;
const CIGAR_INS: u32 = 1;
const CIGAR_DEL: u32 = 2;
const CIGAR_REF_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_PAD: u32 = 6;
const CIGAR_EQUAL: u32 = 7;
const CIGAR_DIFF: u32 = 8;

/// The `indel` `resolve_cigar2` gives the last position of CIGAR operation `k`: minus the length
/// of the deletions that follow it (unless it is one), or the length of the insertions that
/// follow it, padding skipped.
fn indel_after(cigar: &[u32], k: usize) -> i32 {
    let (op, rest) = (cigar[k] & 0xf, &cigar[k + 1..]);
    let Some(&next) = rest.first() else { return 0 };
    let (op2, len2) = (next & 0xf, (next >> 4) as i32);
    let len = |c: &u32| (c >> 4) as i32;
    if op2 == CIGAR_DEL && op != CIGAR_DEL {
        -(len2 + rest[1..].iter().take_while(|&&c| c & 0xf == CIGAR_DEL).map(len).sum::<i32>())
    } else if op2 == CIGAR_INS {
        let mut indel = len2;
        for c in &rest[1..] {
            match c & 0xf {
                CIGAR_INS => indel += len(c),
                CIGAR_PAD => {}
                _ => break,
            }
        }
        indel
    } else if op2 == CIGAR_PAD && rest.len() > 1 {
        let mut inserted = 0;
        for c in &rest[1..] {
            match c & 0xf {
                CIGAR_INS => inserted += len(c),
                CIGAR_DEL | CIGAR_MATCH | CIGAR_REF_SKIP | CIGAR_EQUAL | CIGAR_DIFF => break,
                _ => {}
            }
        }
        inserted
    } else {
        0
    }
}

/// A read whose positions are being counted, with its place in its CIGAR.
struct ActiveRead {
    cigar: Vec<u32>,
    /// The query bases, two `seq_nt16` codes per byte.
    seq: Vec<u8>,
    query_len: i32,
    reverse: bool,
    /// The CIGAR operation reached, the reference position it starts at, and the query position.
    k: usize,
    x: i64,
    y: i32,
    /// One past the read's last reference position.
    end: i64,
}

impl ActiveRead {
    fn new(read: &Read<'_>) -> Self {
        let cigar = read.cigar().to_vec();
        let reference_len: i64 = cigar
            .iter()
            .filter(|&&c| matches!(c & 0xf, CIGAR_MATCH | CIGAR_DEL | CIGAR_REF_SKIP | CIGAR_EQUAL | CIGAR_DIFF))
            .map(|&c| i64::from(c >> 4))
            .sum();
        Self {
            cigar,
            seq: read.packed_seq().to_vec(),
            query_len: read.query_len(),
            reverse: read.is_reverse(),
            k: 0,
            x: read.pos(),
            y: 0,
            end: read.pos() + reference_len,
        }
    }

    /// 4-bit `seq_nt16` code of query base `i`.
    fn base_code(&self, i: i32) -> u8 {
        assert!(0 <= i && i < self.query_len, "query position {i} outside read of {} bases", self.query_len);
        (self.seq[i as usize >> 1] >> ((!i & 1) << 2)) & 0xf
    }
}

/// The counts of one tile's positions, added read by read.
#[derive(Default)]
struct TileCounts {
    start: i64,
    end: i64,
    /// Base and deletion features of each position.
    rows: Vec<[i32; FEATURES]>,
    /// Reads at each position, without those skipping it (`N`).
    depth: Vec<u32>,
    /// Whether any read, skipping or not, covers each position (`bam_plp` gives it a column).
    covered: Vec<bool>,
    /// Bases without a feature (`N` and other IUPAC codes) at each position.
    unknown: Vec<u32>,
    /// `(position, length, reverse)` of each deletion after a position, in read order.
    deletions: Vec<(u32, u32, bool)>,
    /// `(position, reverse, bases)` of each insertion after a position, in read order; the bases
    /// are a range of `inserted`.
    insertions: Vec<(u32, bool, std::ops::Range<usize>)>,
    inserted: Vec<u8>,
}

impl TileCounts {
    /// Empties the counts for positions `start..end`.
    fn reset(&mut self, start: i64, end: i64) {
        let len = usize::try_from(end - start).unwrap_or(0);
        (self.start, self.end) = (start, end);
        for (v, zero) in [(&mut self.depth, 0), (&mut self.unknown, 0)] {
            v.clear();
            v.resize(len, zero);
        }
        self.rows.clear();
        self.rows.resize(len, [0; FEATURES]);
        self.covered.clear();
        self.covered.resize(len, false);
        self.deletions.clear();
        self.insertions.clear();
        self.inserted.clear();
    }

    /// Adds the positions of `read` in the tile (as `resolve_cigar2` would see them a position
    /// at a time) and moves the read's place to the tile's end.
    fn add(&mut self, read: &mut ActiveRead) {
        let strand = if read.reverse { 16 } else { 0 };
        while read.k < read.cigar.len() {
            let c = read.cigar[read.k];
            let (op, len) = (c & 0xf, i64::from(c >> 4));
            match op {
                CIGAR_MATCH | CIGAR_EQUAL | CIGAR_DIFF | CIGAR_DEL | CIGAR_REF_SKIP => {
                    let x = read.x;
                    if x >= self.end {
                        return;
                    }
                    let first = (self.start - x).clamp(0, len);
                    let last = (self.end - x).clamp(0, len);
                    for i in first..last {
                        let at = (x + i - self.start) as usize;
                        self.covered[at] = true;
                        if op == CIGAR_REF_SKIP {
                            continue;
                        }
                        let is_del = op == CIGAR_DEL;
                        let qpos = if is_del { read.y } else { read.y + i as i32 };
                        let indel = if i == len - 1 { indel_after(&read.cigar, read.k) } else { 0 };
                        if indel < 0 {
                            self.deletions.push((at as u32, indel.unsigned_abs(), read.reverse));
                        }
                        self.depth[at] += 1;
                        let feature = if is_del {
                            Some(if read.reverse { REV_DEL } else { FWD_DEL })
                        } else {
                            BASE_FEATURE[usize::from(read.base_code(qpos)) + strand]
                        };
                        match feature {
                            Some(f) => self.rows[at][f] += 1,
                            None => self.unknown[at] += 1,
                        }
                        if indel > 0 {
                            let skip = if is_del { 0 } else { 1 };
                            let begin = self.inserted.len();
                            self.inserted
                                .extend((skip..indel + skip).map(|j| NT16[usize::from(read.base_code(qpos + j))]));
                            self.insertions.push((at as u32, read.reverse, begin..self.inserted.len()));
                        }
                    }
                    if x + len > self.end {
                        // The operation goes on in the next tile.
                        return;
                    }
                    read.x += len;
                    if !matches!(op, CIGAR_DEL | CIGAR_REF_SKIP) {
                        read.y += len as i32;
                    }
                }
                CIGAR_INS | CIGAR_SOFT_CLIP => read.y += len as i32,
                _ => {}
            }
            read.k += 1;
        }
    }
}

/// The columns emitted so far, and what carries over from one tile to the next.
struct Columns {
    out: PileupCounts,
    scratch: ColumnScratch,
    pre_pos: i64,
    contiguous_flanking: usize,
}

impl Columns {
    fn new(bounds: (i64, i64)) -> Self {
        Self {
            out: PileupCounts { bounds, ..PileupCounts::default() },
            scratch: ColumnScratch::default(),
            pre_pos: 0,
            contiguous_flanking: 0,
        }
    }

    /// Emits the covered positions of a counted tile, in order, as `pileup_counts` does.
    fn add_tile(&mut self, tile: &mut TileCounts, reference: &ReferenceWindow, params: &CandidateParams) {
        // Stable sorts keep each position's indels in read order, which the insertion hash
        // tables' iteration order depends on.
        tile.deletions.sort_by_key(|&(at, ..)| at);
        tile.insertions.sort_by_key(|&(at, ..)| at);
        let (mut deletions, mut insertions) = (tile.deletions.iter().peekable(), tile.insertions.iter().peekable());
        let out = &mut self.out;
        for at in 0..tile.rows.len() {
            if !tile.covered[at] {
                continue;
            }
            let pos = tile.start + at as i64;
            if self.pre_pos + 1 != pos || self.pre_pos == 0 {
                self.contiguous_flanking = 0;
            } else {
                self.contiguous_flanking += 1;
            }
            self.pre_pos = pos;
            // C quirk: an N (or other IUPAC) base indexes feature -1, which is the last feature
            // of the previous row; for the first row it lands outside the matrix and is lost.
            if tile.unknown[at] > 0
                && let Some(prev) = out.rows.last_mut()
            {
                prev[FEATURES - 1] += tile.unknown[at] as i32;
            }
            let ColumnScratch { dels_f, dels_r } = &mut self.scratch;
            let mut longest_del = 0usize;
            while let Some(&(_, d, reverse)) = deletions.next_if(|&&(a, ..)| a as usize == at) {
                let d = d as usize;
                if d >= dels_f.len() {
                    let new_len = d.max(2 * dels_f.len());
                    dels_f.resize(new_len, 0);
                    dels_r.resize(new_len, 0);
                }
                if reverse {
                    dels_r[d - 1] += 1
                } else {
                    dels_f[d - 1] += 1
                }
                longest_del = longest_del.max(d);
            }
            let mut inserted = Insertions::default();
            while let Some((_, reverse, bases)) = insertions.next_if(|&&(a, ..)| a as usize == at) {
                inserted.add(&tile.inserted[bases.clone()], *reverse);
            }
            let flanks_ok = params.flanks == FlankRequirement::None || self.contiguous_flanking >= FLANKING_BASES;
            let column = ColumnCounts {
                row: tile.rows[at],
                depth: tile.depth[at] as usize,
                dels: (&dels_f[..longest_del], &dels_r[..longest_del]),
                inserted: &inserted,
            };
            let (depth, candidate) = finish_column(&mut out.rows, column, pos, reference, params, flanks_ok);
            dels_f[..longest_del].fill(0);
            dels_r[..longest_del].fill(0);
            out.positions.push(pos);
            out.depths.push(depth);
            if let Some(alt_info) = candidate {
                out.candidates.push(alt_info);
            }
        }
    }
}

/// The inserted sequences after a position, counted per strand and in all
/// (`kh_counter` tables, whose iteration order depends on the order of the increments).
#[derive(Debug, Default)]
struct Insertions {
    forward: StrCounter,
    reverse: StrCounter,
    all: StrCounter,
}

impl Insertions {
    fn add(&mut self, bases: &[u8], reverse: bool) {
        if reverse {
            self.reverse.increment(bases)
        } else {
            self.forward.increment(bases)
        }
        self.all.increment(bases);
    }
}

/// One position's counts before the reference columns are filled in.
struct ColumnCounts<'a> {
    row: [i32; FEATURES],
    depth: usize,
    /// Deletion counts by length - 1, forward and reverse.
    dels: (&'a [usize], &'a [usize]),
    inserted: &'a Insertions,
}

/// The column engine: `pileup_counts` from htslib's `bam_plp` columns, which the read-by-read
/// engine is checked against.
#[cfg(test)]
fn pileup_counts_by_column(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    region: &str,
    params: &CandidateParams,
) -> Result<PileupCounts, Error> {
    let (contig, start, end) = parse_region_bounds(region)?;
    let tid = bam.tid(&contig);
    let reference = ReferenceWindow::fetch(fasta, &contig, start, end)?;
    let mut pileup = bam.pileup(region, params.min_mapq)?;
    let mut out = PileupCounts { bounds: (start, end), ..PileupCounts::default() };
    let mut scratch = ColumnScratch::default();
    let mut pre_pos: i64 = 0;
    let mut contiguous_flanking = 0usize;
    while let Some(column) = pileup.next_column()? {
        if Some(column.tid) != tid {
            continue;
        }
        let pos = i64::from(column.pos);
        if pos < start {
            continue;
        }
        if pos >= end {
            break;
        }
        if pre_pos + 1 != pos || pre_pos == 0 {
            contiguous_flanking = 0;
        } else {
            contiguous_flanking += 1;
        }
        pre_pos = pos;

        let flanks_ok = params.flanks == FlankRequirement::None || contiguous_flanking >= FLANKING_BASES;
        let (depth, candidate) = count_column(&mut out.rows, &mut scratch, &column, pos, &reference, params, flanks_ok);
        out.positions.push(pos);
        out.depths.push(depth);
        if let Some(alt_info) = candidate {
            out.candidates.push(alt_info);
        }
    }
    Ok(out)
}

/// Deletion lengths counted per column before the counts grow (`int del_counts[32]` in C).
const DELETION_SLOTS: usize = 32;

/// Buffers `count_column` reuses from one column to the next, so that a column allocates
/// nothing unless it holds an insertion or a deletion longer than any before.
#[derive(Debug)]
struct ColumnScratch {
    /// Deletion counts by length - 1, forward and reverse strand, all zero between columns.
    /// Zeros past a column's longest deletion count for nothing, so the vectors only grow.
    dels_f: Vec<usize>,
    dels_r: Vec<usize>,
}

impl Default for ColumnScratch {
    fn default() -> Self {
        Self { dels_f: vec![0; DELETION_SLOTS], dels_r: vec![0; DELETION_SLOTS] }
    }
}

/// Counts one `bam_plp` column into a new row of `rows`; returns its gVCF counts and its alt
/// info if it is a candidate.
#[cfg(test)]
fn count_column(
    rows: &mut Vec<[i32; FEATURES]>,
    scratch: &mut ColumnScratch,
    column: &Column<'_>,
    pos: i64,
    reference: &ReferenceWindow,
    params: &CandidateParams,
    flanks_ok: bool,
) -> (SiteDepth, Option<String>) {
    let mut row = [0i32; FEATURES];
    let mut depth = 0usize;
    let ColumnScratch { dels_f, dels_r } = scratch;
    let mut indel_buf = Vec::new();
    // The longest deletion in the column: the counts past it are zero.
    let mut longest_del = 0usize;
    let mut inserted = Insertions::default();

    for entry in column.entries() {
        if entry.is_refskip() {
            continue;
        }
        let reverse = entry.is_reverse();
        if entry.indel() < 0 {
            let d = entry.indel().unsigned_abs() as usize;
            if d >= dels_f.len() {
                let new_len = d.max(2 * dels_f.len());
                dels_f.resize(new_len, 0);
                dels_r.resize(new_len, 0);
            }
            if reverse {
                dels_r[d - 1] += 1
            } else {
                dels_f[d - 1] += 1
            }
            longest_del = longest_del.max(d);
        }

        depth += 1;
        let feature = if entry.is_del() {
            Some(if reverse { REV_DEL } else { FWD_DEL })
        } else {
            BASE_FEATURE[usize::from(entry.base_code(entry.qpos())) + if reverse { 16 } else { 0 }]
        };
        match feature {
            Some(f) => row[f] += 1,
            None => {
                if let Some(prev) = rows.last_mut() {
                    prev[FEATURES - 1] += 1;
                }
            }
        }

        if entry.indel() > 0 {
            let first = if entry.is_del() { 0 } else { 1 };
            indel_buf.clear();
            indel_buf
                .extend((first..entry.indel() + first).map(|j| NT16[usize::from(entry.base_code(entry.qpos() + j))]));
            inserted.add(&indel_buf, reverse);
        }
    }
    let column =
        ColumnCounts { row, depth, dels: (&dels_f[..longest_del], &dels_r[..longest_del]), inserted: &inserted };
    let result = finish_column(rows, column, pos, reference, params, flanks_ok);
    dels_f[..longest_del].fill(0);
    dels_r[..longest_del].fill(0);
    result
}

/// Fills in a column's indel and reference features, pushes its row onto `rows`, and returns
/// its gVCF counts and its alt info if it is a candidate (the C pileup function).
fn finish_column(
    rows: &mut Vec<[i32; FEATURES]>,
    column: ColumnCounts<'_>,
    pos: i64,
    reference: &ReferenceWindow,
    params: &CandidateParams,
    flanks_ok: bool,
) -> (SiteDepth, Option<String>) {
    let ColumnCounts { mut row, depth, dels: (dels_f, dels_r), inserted } = column;
    let mut del_count = 0usize;
    for (dels, all_col, best_col) in [(dels_f, FWD_DEL_ALL, FWD_DEL_BEST), (dels_r, REV_DEL_ALL, REV_DEL_BEST)] {
        let all: usize = dels.iter().sum();
        let best = dels.iter().copied().max().unwrap_or(0);
        row[all_col] = all as i32;
        row[best_col] = best as i32;
        del_count += all;
    }
    let mut ins_count = 0usize;
    for (counter, all_col, best_col) in
        [(&inserted.forward, FWD_INS_ALL, FWD_INS_BEST), (&inserted.reverse, REV_INS_ALL, REV_INS_BEST)]
    {
        let (sum, max) = counter.stats();
        row[all_col] = sum as i32;
        row[best_col] = max as i32;
        ins_count += sum;
    }

    let offset = reference.offset(pos);
    let ref_base = upper_base(reference.byte(offset));
    let ref_fwd = reference_feature(ref_base);
    let mut ref_count = 0usize;
    let mut alt_count = 0usize;
    let mut all_alt_count = 0usize;
    let mut major_alt_base = 0u8;
    let (mut forward_sum, mut reverse_sum) = (0i32, 0i32);
    for (i, &base) in FEATURE_BASES.iter().enumerate() {
        forward_sum += row[i];
        reverse_sum += row[i + REVERSE_OFFSET];
        let count = (row[i] + row[i + REVERSE_OFFSET]) as usize;
        if i == ref_fwd {
            ref_count = count;
        } else if count > alt_count {
            alt_count = count;
            major_alt_base = base;
            // C quirk: only counts that raise the maximum are added.
            all_alt_count += alt_count;
        }
    }
    let site =
        SiteDepth { reference: ref_count as u64, total: (ref_count + all_alt_count + del_count + ins_count) as u64 };
    row[ref_fwd] = -forward_sum;
    row[ref_fwd + REVERSE_OFFSET] = -reverse_sum;
    rows.push(row);

    let depth = depth.max(1);
    let af = |count: usize| count as f32 / depth as f32;
    let pass_af = match params.alleles {
        CandidateAlleles::SnpsOnly => af(alt_count) >= params.min_snp_af,
        CandidateAlleles::SnpsAndIndels => {
            let non_ref_majority = ref_count < alt_count || ref_count < ins_count || ref_count < del_count;
            let tie_with_alt = ref_count > 0 && ref_count == alt_count && ref_base < major_alt_base;
            non_ref_majority
                || tie_with_alt
                || af(alt_count) >= params.min_snp_af
                || af(del_count) >= params.min_indel_af
                || af(ins_count) >= params.min_indel_af
        }
    };
    let is_candidate =
        pass_af && depth >= params.min_depth && matches!(ref_base, b'A' | b'C' | b'G' | b'T') && flanks_ok;
    if !is_candidate {
        return (site, None);
    }

    let row = rows.last().expect("row was just pushed");
    let ref_char = char::from(ref_base);
    let mut alt_info = format!("{}-{depth}-{ref_char}-", pos + 1);
    for (i, &base) in FEATURE_BASES.iter().enumerate() {
        let alt_sum = row[i] + row[i + REVERSE_OFFSET];
        if alt_sum > 0 && i != ref_fwd {
            write!(alt_info, "X{} {alt_sum} ", char::from(base)).expect("writing to a String cannot fail");
        }
    }
    let mut ref_depth = ref_count as i64;
    for (i, (f, r)) in dels_f.iter().zip(dels_r).enumerate() {
        let d = f + r;
        ref_depth -= d as i64;
        if d > 0 && i < params.max_indel_length {
            let deleted = String::from_utf8_lossy(reference.slice(offset + 1, i + 1));
            write!(alt_info, "D{deleted} {d} ").expect("writing to a String cannot fail");
        }
    }
    for (key, val) in inserted.all.iter() {
        ref_depth -= i64::from(val);
        if key.len() <= params.max_indel_length {
            let key = String::from_utf8_lossy(key);
            write!(alt_info, "I{ref_char}{key} {val} ").expect("writing to a String cannot fail");
        }
    }
    if ref_depth > 0 {
        write!(alt_info, "R{ref_char} {ref_depth} ").expect("writing to a String cannot fail");
    }
    (site, Some(alt_info))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_bounds_follow_hts_parse_reg() {
        assert_eq!(parse_region_bounds("chr20:0-100").unwrap(), ("chr20".into(), 0, 100));
        assert_eq!(parse_region_bounds("chr20:1-100").unwrap(), ("chr20".into(), 0, 100));
        assert_eq!(parse_region_bounds("chr20:5,001-6,000").unwrap(), ("chr20".into(), 5000, 6000));
        assert_eq!(parse_region_bounds("HLA-A:10-20").unwrap(), ("HLA-A".into(), 9, 20));
        assert!(parse_region_bounds("chr20").is_err());
    }

    #[test]
    fn reference_window_reads_nul_past_its_end() {
        let window = ReferenceWindow { seq: b"ACGTACGT".to_vec(), seq_start: 100 };
        assert_eq!(window.offset(103), 3);
        assert_eq!(window.byte(3), b'T');
        assert_eq!(window.byte(8), 0);
        assert_eq!(window.slice(6, 5), b"GT");
    }

    #[test]
    fn reference_feature_defaults_to_a() {
        assert_eq!(b"ACGTN".map(reference_feature), [0, 1, 2, 3, 0]);
    }

    /// `x = x * 1103515245 + 12345`, for test data.
    fn next(state: &mut u32) -> u32 {
        *state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        *state >> 16
    }

    /// A read's sequence for `cigar` at 0-based `pos`: the reference with every seventh base
    /// changed and an N at a few, and random inserted and clipped bases.
    fn read_seq(reference: &[u8], pos: usize, cigar: &str, seed: u32) -> String {
        let (mut seq, mut x, mut state, mut number) = (String::new(), pos, seed, 0usize);
        for c in cigar.chars() {
            if let Some(d) = c.to_digit(10) {
                number = number * 10 + d as usize;
                continue;
            }
            for _ in 0..number {
                let random = char::from(b"ACGT"[(next(&mut state) % 4) as usize]);
                match c {
                    'M' | '=' | 'X' => {
                        let r = next(&mut state);
                        let base = match r % 23 {
                            0..=2 => random,
                            3 => 'N',
                            _ => char::from(reference[x]),
                        };
                        seq.push(base);
                        x += 1;
                    }
                    'I' | 'S' => seq.push(random),
                    'D' | 'N' => x += 1,
                    _ => {}
                }
            }
            number = 0;
        }
        seq
    }

    #[test]
    fn reads_give_the_columns_of_htslib_pileup() {
        use pasta::bam;
        use pasta::sam::{self, alignment::io::Write as _};

        let dir = std::env::temp_dir().join(format!("snv-pileup-reads-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temporary directory");
        let mut state = 7u32;
        let mut reference: Vec<u8> = (0..6000).map(|_| b"ACGT"[(next(&mut state) % 4) as usize]).collect();
        reference[150] = b'N';
        let fasta = dir.join("ref.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", String::from_utf8_lossy(&reference))).expect("FASTA");
        snv_hts::index_fasta(&fasta).expect("FASTA index");

        // (1-based POS, flag, MAPQ, CIGAR), sorted by position.
        let reads: &[(usize, u16, u8, &str)] = &[
            (1, 0, 60, "20M"),
            (3, 16, 60, "5M2I8M"),
            (5, 0, 60, "4M3D6M"),
            (5, 0, 60, "3M1D2D5M"),
            (8, 16, 60, "2S6M4N6M"),
            (10, 0, 60, "4M1P2I4M"),
            (10, 16, 60, "3M2I1P1I4M"),
            (12, 0, 60, "5M2D1I4M"),
            (15, 0, 60, "6=2X4M3I"),
            (15, 16, 60, "3H8M2S"),
            (16, 0, 60, "1D9M"),
            (18, 0, 60, "10M"),
            (20, 0, 0, "10M"),
            (20, 1024, 60, "10M"),
            (25, 16, 60, "5M40D5M"),
            (30, 0, 60, "4M2D3M1I5M"),
            (100, 0, 60, "3M20N3M"),
            (140, 0, 60, "5M3I5M"),
            (140, 16, 60, "5M3I5M"),
            (140, 0, 60, "5M2I5M"),
            (140, 16, 60, "5M3I5M"),
            (140, 0, 60, "5M1I5M"),
            (141, 0, 60, "4M4I6M"),
            (145, 16, 60, "10M"),
            (146, 0, 60, "8M"),
            // Across tile boundaries (2,048 positions from the region start or the first read).
            (1500, 16, 60, "3000M"),
            (2000, 0, 60, "100M"),
            (2030, 0, 60, "15M5D20M"),
            (2040, 16, 60, "8M3I10M"),
            (2045, 0, 60, "3M10N12M"),
            (2048, 0, 60, "1M2I5M"),
            (2048, 16, 60, "1M3D4M"),
            (2049, 0, 60, "2D6M"),
            (4040, 0, 60, "6M2I4M30D6M"),
            (4045, 16, 60, "3M1D2D9M"),
        ];
        let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:6000\n".parse().expect("SAM header");
        let path = dir.join("reads.bam");
        {
            let mut writer = bam::io::Writer::new(std::fs::File::create(&path).expect("BAM"));
            writer.write_header(&header).expect("BAM header");
            for (i, &(pos, flag, mapq, cigar)) in reads.iter().enumerate() {
                let seq = read_seq(&reference, pos - 1, cigar, i as u32 + 1);
                let qual = "I".repeat(seq.len());
                let sam = format!("r{i}\t{flag}\tchr1\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}");
                let record = sam::Record::try_from(sam.as_bytes()).expect("SAM line");
                writer.write_alignment_record(&header, &record).expect("write");
            }
            writer.try_finish().expect("BAM end");
        }
        snv_hts::index_bam(&path).expect("BAM index");
        let mut check = bam::io::reader::Builder.build_from_path(&path).expect("BAM reads back");
        check.read_header().expect("BAM header");
        assert_eq!(check.records().count(), reads.len());

        let fasta = Fasta::open(&fasta).expect("FASTA");
        for flanks in [FlankRequirement::Contiguous, FlankRequirement::None] {
            let params = CandidateParams {
                min_depth: 1,
                min_snp_af: 0.0,
                min_indel_af: 0.0,
                min_mapq: 5,
                max_indel_length: 50,
                alleles: CandidateAlleles::SnpsAndIndels,
                flanks,
            };
            for region in ["chr1:1-6000", "chr1:12-145", "chr1:30-31", "chr1:2000-4100", "chr1:2047-2049"] {
                let mut bam = AlignmentFile::open(&path, dir.join("ref.fa")).expect("BAM");
                let by_read = pileup_counts(&mut bam, &fasta, region, &params).expect("by read");
                let by_column = pileup_counts_by_column(&mut bam, &fasta, region, &params).expect("by column");
                assert!(
                    !by_read.candidates.is_empty() || !region.ends_with("-6000") && region != "chr1:2000-4100",
                    "{region}: no candidates"
                );
                assert_eq!(by_read, by_column, "{region} {flanks:?}");
            }
        }
        std::fs::remove_dir_all(&dir).ok();
    }
}

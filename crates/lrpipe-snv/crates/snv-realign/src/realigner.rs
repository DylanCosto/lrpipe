//! Port of `preprocess/realign/realigner.cpp` (`realign_reads`): aligns a window's reads to its
//! candidate haplotypes, then maps each read's best haplotype alignment back to the reference.

use std::collections::{HashMap, VecDeque};

use crate::cppsort::sort_by;
use crate::ssw::Aligner;

/// `kmer_size_`.
const KMER: usize = 32;
/// `max_num_of_mismatches_`.
const MAX_MISMATCHES: i32 = 2;
/// `match_score_` and `mismatch_penalty_` (the fast alignment's scores).
const MATCH: i32 = 4;
const MISMATCH: i32 = 6;
/// `CalculateSswAlignmentScoreThreshold` with `read_size_ = 250` and `similarity_threshold_ =
/// 0.16934`: negative, so 1.
const SSW_SCORE_THRESHOLD: i32 = 1;

/// `Operation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Op {
    Unspecified,
    Match,
    Insert,
    Delete,
    SoftClip,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CigarOp {
    op: Op,
    len: i32,
}

/// `CigarStringToVector`: every `(\d+)([XIDS=])` match, case-insensitively (`M` is skipped).
fn parse_cigar(cigar: &str) -> VecDeque<CigarOp> {
    let bytes = cigar.as_bytes();
    let mut out = VecDeque::new();
    let mut i = 0;
    while i < bytes.len() {
        if !bytes[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let Some(&c) = bytes.get(i) else { break };
        let op = match c.to_ascii_uppercase() {
            b'=' | b'X' => Op::Match,
            b'S' => Op::SoftClip,
            b'D' => Op::Delete,
            b'I' => Op::Insert,
            _ => continue,
        };
        // `atoi` of the digits.
        let len = cigar[start..i].parse::<i64>().map_or(i32::MAX, |v| v.clamp(0, i64::from(i32::MAX)) as i32);
        out.push_back(CigarOp { op, len });
        i += 1;
    }
    out
}

/// `CigarVectorToString`: matches are written `X`.
fn cigar_string(ops: &[CigarOp]) -> String {
    use std::fmt::Write as _;
    let mut s = String::new();
    for op in ops {
        let letter = match op.op {
            Op::Match => "X",
            Op::Insert => "I",
            Op::Delete => "D",
            Op::SoftClip => "S",
            Op::Unspecified => "",
        };
        write!(s, "{}{letter}", op.len).expect("writing to a String cannot fail");
    }
    s
}

/// `ReadAlignment`.
#[derive(Debug, Clone)]
struct ReadAlignment {
    position: i32,
    cigar: String,
    score: i32,
}

impl ReadAlignment {
    const NOT_ALIGNED: i32 = -1;

    fn new() -> Self {
        Self { position: Self::NOT_ALIGNED, cigar: String::new(), score: 0 }
    }
}

/// `HaplotypeReadsAlignment`.
#[derive(Debug, Clone)]
struct HaplotypeAlignment {
    haplotype_index: usize,
    haplotype_score: i32,
    reads: Vec<ReadAlignment>,
    cigar: String,
    cigar_ops: VecDeque<CigarOp>,
    ref_pos: i32,
    positions_map: Vec<i32>,
    is_reference: bool,
}

/// One read handed to `realign_reads`.
#[derive(Debug, Clone)]
pub struct WindowRead<'a> {
    /// Upper-case bases.
    pub seq: &'a [u8],
    /// 0-based alignment start.
    pub position: i32,
    /// CIGAR string.
    pub cigar: &'a str,
}

/// `ReAligner::realign_reads`: each read's `(position, cigar)` after realignment; reads without
/// a better haplotype alignment keep theirs. Realigned CIGARs write matches as `X`.
///
/// `reference` is the window's reference with `ref_prefix` and `ref_suffix` flanking bases,
/// starting at 0-based `ref_start`; every haplotype carries the same flanks.
pub fn realign_reads(
    reads: &[WindowRead<'_>],
    reference: &[u8],
    haplotypes: &[Vec<u8>],
    ref_start: i32,
    ref_prefix: i32,
    ref_suffix: i32,
) -> Vec<(i32, String)> {
    let seqs: Vec<&[u8]> = reads.iter().map(|r| r.seq).collect();
    let mut index: HashMap<&[u8], Vec<(usize, usize)>> = HashMap::new();
    for (id, read) in seqs.iter().enumerate() {
        if read.len() <= KMER {
            continue;
        }
        for i in 0..=read.len() - KMER {
            index.entry(&read[i..i + KMER]).or_default().push((id, i));
        }
    }

    // FastAlignReadsToHaplotypes.
    let mut alignments: Vec<HaplotypeAlignment> = Vec::new();
    for (h, haplotype) in haplotypes.iter().enumerate() {
        let mut scores = vec![ReadAlignment::new(); seqs.len()];
        let score = fast_align(haplotype, reference, &seqs, &index, &mut scores, ref_prefix, ref_suffix);
        if score == 0 {
            scores = vec![ReadAlignment::new(); seqs.len()];
        }
        alignments.push(HaplotypeAlignment {
            haplotype_index: h,
            haplotype_score: score,
            reads: scores,
            cigar: String::new(),
            cigar_ops: VecDeque::new(),
            ref_pos: 0,
            positions_map: Vec::new(),
            is_reference: false,
        });
    }

    // AlignHaplotypesToReference.
    let aligner = Aligner::new(reference);
    for alignment in &mut alignments {
        let haplotype = &haplotypes[alignment.haplotype_index];
        if let Some(a) = aligner.align(haplotype).filter(|a| a.score > 0) {
            alignment.is_reference = a.cigar == format!("{}=", haplotype.len());
            alignment.cigar_ops = parse_cigar(&a.cigar);
            alignment.cigar = a.cigar;
            alignment.ref_pos = a.ref_begin;
        }
    }

    // CalculatePositionMaps.
    for alignment in &mut alignments {
        alignment.positions_map = positions_map(haplotypes[alignment.haplotype_index].len(), &alignment.cigar);
    }

    // SswAlignReadsToHaplotypes.
    for (i, read) in seqs.iter().enumerate() {
        if alignments.iter().any(|a| a.reads[i].score > 0) {
            continue;
        }
        for alignment in &mut alignments {
            if alignment.haplotype_score == 0 {
                continue;
            }
            let hap_aligner = Aligner::new(&haplotypes[alignment.haplotype_index]);
            let Some(a) = hap_aligner.align(read) else { continue };
            let score = i32::from(a.score);
            // `sw_score > 0 && sw_score >= threshold`: the threshold is 1.
            if score >= SSW_SCORE_THRESHOLD && alignment.reads[i].score < score {
                alignment.reads[i] = ReadAlignment { position: a.ref_begin, cigar: a.cigar, score };
            }
        }
    }

    // `sort(read_to_haplotype_alignments_)` by haplotype score, as libstdc++ sorts.
    sort_by(&mut alignments, |a, b| a.haplotype_score < b.haplotype_score);

    // RealignReadsToReference.
    reads
        .iter()
        .enumerate()
        .map(|(i, read)| {
            let Some(best) = best_alignment(&alignments, haplotypes.len(), i) else {
                return (read.position, read.cigar.to_owned());
            };
            let hap = &alignments[best];
            let read_to_hap = &hap.reads[i];
            let map_shift =
                usize::try_from(read_to_hap.position).ok().and_then(|p| hap.positions_map.get(p)).copied().unwrap_or(0);
            let new_position = ref_start + hap.ref_pos + read_to_hap.position + map_shift;
            let ops = read_to_ref(read.seq.len() as i32, read_to_hap, &hap.cigar_ops);
            if ops.is_empty() { (read.position, read.cigar.to_owned()) } else { (new_position, cigar_string(&ops)) }
        })
        .collect()
}

/// `FastAlignReadsToHaplotype`: ungapped placements of reads sharing a 32-mer with the
/// haplotype, with at most two mismatches. Returns the haplotype score, or 0 when a non-reference
/// haplotype has an uncovered base outside the flanks.
fn fast_align(
    haplotype: &[u8],
    reference: &[u8],
    reads: &[&[u8]],
    index: &HashMap<&[u8], Vec<(usize, usize)>>,
    scores: &mut [ReadAlignment],
    ref_prefix: i32,
    ref_suffix: i32,
) -> i32 {
    let is_ref = haplotype == reference;
    let mut coverage = vec![0i32; haplotype.len()];
    let mut haplotype_score = 0;
    if haplotype.len() < KMER {
        return 0;
    }
    for i in 0..=haplotype.len() - KMER {
        let Some(hits) = index.get(&haplotype[i..i + KMER]) else { continue };
        for &(read_id, read_pos) in hits {
            let target_start = (i as i64 - read_pos as i64).max(0) as usize;
            let read = reads[read_id];
            if target_start + read.len() > haplotype.len() {
                continue;
            }
            let alignment = &mut scores[read_id];
            if alignment.position != ReadAlignment::NOT_ALIGNED && alignment.position == target_start as i32 {
                continue;
            }
            let (score, mismatches) = fast_align_strings(&haplotype[target_start..target_start + read.len()], read);
            if mismatches <= MAX_MISMATCHES {
                let old = alignment.score;
                for c in &mut coverage[target_start..target_start + read.len()] {
                    *c += 1;
                }
                if old < score {
                    alignment.score = score;
                    haplotype_score += score - old;
                    alignment.position = target_start as i32;
                    alignment.cigar = format!("{}=", read.len());
                }
            }
        }
        let suffix_start = haplotype.len() as i64 - i64::from(ref_suffix);
        if coverage[i] == 0 && i as i64 >= i64::from(ref_prefix) && (i as i64) < suffix_start && !is_ref {
            return 0;
        }
    }
    haplotype_score
}

/// `FastAlignStrings` with `max_mismatches = 3`: `(score, mismatches)`; `N` matches anything.
fn fast_align_strings(s1: &[u8], s2: &[u8]) -> (i32, i32) {
    let (mut matches, mut mismatches) = (0, 0);
    for (&c1, &c2) in s1.iter().zip(s2) {
        if c1 != c2 && c1 != b'N' && c2 != b'N' {
            mismatches += 1;
            if mismatches == MAX_MISMATCHES + 1 {
                return (0, mismatches);
            }
        } else {
            matches += 1;
        }
    }
    (matches * MATCH - mismatches * MISMATCH, mismatches)
}

/// `SetPositionsMap`: for each haplotype base, the shift to its reference position.
fn positions_map(haplotype_len: usize, cigar: &str) -> Vec<i32> {
    let mut map = vec![0i32; haplotype_len];
    let (mut shift, mut pos) = (0i32, 0usize);
    let bytes = cigar.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if !bytes[i].is_ascii_digit() {
            i += 1;
            continue;
        }
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let Some(&c) = bytes.get(i) else { break };
        let len: usize = cigar[start..i].parse().unwrap_or(0);
        match c.to_ascii_uppercase() {
            b'=' | b'X' => {
                for _ in 0..len {
                    if let Some(slot) = map.get_mut(pos) {
                        *slot = shift;
                    }
                    pos += 1;
                }
            }
            b'S' => {
                shift -= len as i32;
                for _ in 0..len {
                    if let Some(slot) = map.get_mut(pos) {
                        *slot = shift;
                    }
                    pos += 1;
                }
            }
            b'D' => shift += len as i32,
            b'I' => {
                for _ in 0..len {
                    if let Some(slot) = map.get_mut(pos) {
                        *slot = shift;
                    }
                    shift -= 1;
                    pos += 1;
                }
            }
            _ => {}
        }
        i += 1;
    }
    map
}

/// `GetBestReadAlignment`: the best-scoring haplotype (in sorted order); on ties a later
/// non-reference haplotype wins.
fn best_alignment(alignments: &[HaplotypeAlignment], count: usize, read: usize) -> Option<usize> {
    let mut best_score = 0;
    let mut best = None;
    for (h, alignment) in alignments.iter().enumerate().take(count) {
        let score = alignment.reads[read].score;
        if score > best_score || (best_score > 0 && score == best_score && !alignment.is_reference) {
            best_score = score;
            best = Some(h);
        }
    }
    best
}

/// `AlignedLength`: read bases the CIGAR covers.
fn aligned_length(cigar: &[CigarOp]) -> i32 {
    cigar.iter().filter(|op| op.op != Op::Delete).map(|op| op.len).sum()
}

/// `MergeCigarOp`.
fn merge_op(op: CigarOp, read_len: i32, cigar: &mut Vec<CigarOp>) {
    let last = cigar.last().map_or(Op::Unspecified, |c| c.op);
    let aligned = aligned_length(cigar);
    let len = if op.op == Op::Delete { op.len } else { op.len.min(read_len - aligned) };
    if len <= 0 || aligned == read_len {
        return;
    }
    if op.op == last {
        cigar.last_mut().expect("last exists").len += len;
    } else {
        cigar.push(CigarOp { op: op.op, len });
    }
}

/// `LeftTrimHaplotypeToRefAlignment`.
fn left_trim(ops: &VecDeque<CigarOp>, read_to_hap_pos: i32) -> VecDeque<CigarOp> {
    let mut ops = ops.clone();
    let mut cur = 0;
    while cur != read_to_hap_pos {
        let Some(op) = ops.pop_front() else { break };
        if matches!(op.op, Op::Match | Op::SoftClip | Op::Insert) {
            if op.len + cur > read_to_hap_pos {
                ops.push_front(CigarOp { op: op.op, len: op.len - (read_to_hap_pos - cur) });
            }
            cur = (op.len + cur).min(read_to_hap_pos);
        }
    }
    if ops.front().is_some_and(|op| op.op == Op::Delete) {
        ops.pop_front();
    }
    ops
}

fn is_match(op: CigarOp) -> bool {
    matches!(op.op, Op::Match | Op::SoftClip)
}

fn push_front_if_not_empty(op: CigarOp, ops: &mut VecDeque<CigarOp>) {
    if op.len > 0 {
        ops.push_front(op);
    }
}

/// `CalculateReadToRefAlignment`: composes read-to-haplotype and haplotype-to-reference
/// CIGARs; empty when they cannot be composed.
fn read_to_ref(read_len: i32, read_to_hap: &ReadAlignment, hap_to_ref: &VecDeque<CigarOp>) -> Vec<CigarOp> {
    let mut out: Vec<CigarOp> = Vec::new();
    let mut r = parse_cigar(&read_to_hap.cigar);
    let mut h = left_trim(hap_to_ref, read_to_hap.position);
    if let Some(&front) = r.front()
        && front.op == Op::SoftClip
    {
        merge_op(CigarOp { op: Op::SoftClip, len: front.len }, read_len, &mut out);
        r.pop_front();
    }
    while (!r.is_empty() || !h.is_empty()) && aligned_length(&out) < read_len {
        if !r.is_empty() && h.is_empty() {
            merge_op(r.pop_front().expect("non-empty"), read_len, &mut out);
            continue;
        }
        if r.is_empty() && !h.is_empty() {
            break;
        }
        let mut rop = r.pop_front().expect("non-empty");
        let mut hop = h.pop_front().expect("non-empty");
        if is_match(rop) && is_match(hop) {
            let len = rop.len.min(hop.len);
            let op = if rop.op == Op::SoftClip || hop.op == Op::SoftClip { Op::SoftClip } else { Op::Match };
            merge_op(CigarOp { op, len }, read_len, &mut out);
            rop.len -= len;
            push_front_if_not_empty(rop, &mut r);
            hop.len -= len;
            push_front_if_not_empty(hop, &mut h);
        } else if rop.op == Op::Delete && is_match(hop) {
            merge_op(CigarOp { op: Op::Delete, len: rop.len }, read_len, &mut out);
            hop.len -= rop.len;
            push_front_if_not_empty(hop, &mut h);
        } else if hop.op == Op::Delete && is_match(rop) {
            merge_op(CigarOp { op: Op::Delete, len: hop.len }, read_len, &mut out);
            push_front_if_not_empty(rop, &mut r);
        } else if rop.op == Op::Delete && hop.op == Op::Delete {
            merge_op(CigarOp { op: Op::Delete, len: hop.len + rop.len }, read_len, &mut out);
        } else if rop.op == Op::Insert && is_match(hop) {
            rop.len = (read_len - aligned_length(&out)).min(rop.len);
            merge_op(CigarOp { op: Op::Insert, len: rop.len }, read_len, &mut out);
            push_front_if_not_empty(hop, &mut h);
        } else if hop.op == Op::Insert && is_match(rop) {
            hop.len = (read_len - aligned_length(&out)).min(hop.len);
            merge_op(CigarOp { op: Op::Insert, len: hop.len }, read_len, &mut out);
            rop.len = (rop.len - hop.len).max(0);
            push_front_if_not_empty(rop, &mut r);
        } else if hop.op == Op::Insert && rop.op == Op::Insert {
            merge_op(CigarOp { op: Op::Insert, len: hop.len + rop.len }, read_len, &mut out);
        } else {
            return Vec::new();
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigars_round_trip_with_x_for_matches() {
        let ops: Vec<CigarOp> = parse_cigar("5S10=2X3I4D1=").into_iter().collect();
        assert_eq!(cigar_string(&ops), "5S10X2X3I4D1X");
        assert_eq!(aligned_length(&ops), 21);
    }

    #[test]
    fn positions_map_follows_indels() {
        assert_eq!(positions_map(6, "2=1I3="), [0, 0, 0, -1, -1, -1]);
        assert_eq!(positions_map(6, "2=1D4="), [0, 0, 1, 1, 1, 1]);
        assert_eq!(positions_map(4, "1S3="), [-1, -1, -1, -1]);
    }
}

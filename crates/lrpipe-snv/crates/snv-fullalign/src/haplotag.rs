//! The original's C re-implementation of WhatsHap haplotagging (`haplotag_read` and helpers).
//!
//! Each phased SNP a read covers is "realigned": the read's bases around it are compared by
//! edit distance with the reference and with the reference carrying the alternative base. The
//! votes per phase set decide the read's haplotype. Arithmetic is on `size_t`, so subtractions
//! wrap as in C.

use std::collections::HashMap;

use crate::inputs::PhasedSnp;
use crate::record::{NT16, ReadView, RefWindow, cigar_op, op};

/// `overhang`: reference bases taken left of a SNP (one more is taken on the right).
const OVERHANG: u64 = 10;

/// A read's haplotype (`HAP_UNPHASED`, `HAP_1`, `HAP_2`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub(crate) enum Haplotype {
    #[default]
    Unphased = 0,
    First = 1,
    Second = 2,
}

impl Haplotype {
    /// `normalize_hap`: `HAP_TYPE[haplotype]`.
    pub(crate) fn feature(self) -> i8 {
        match self {
            Self::Unphased => 60,
            Self::First => 30,
            Self::Second => 90,
        }
    }
}

/// Direction of a `cigar_prefix_length` walk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Walk {
    Forward,
    Backward,
}

/// `cigar_prefix_length`: walks CIGAR operations `left..right` (backwards from `right - 1`, or
/// forwards from `left`), the first one shortened to `consumed`, until `reference_bases`
/// reference bases are covered, and returns `(ref_bases, query_bases)`. If they never are, the
/// C function leaves its outputs untouched, which are zero in every caller. Clips are skipped.
fn cigar_prefix_length(
    cigar: &[u32],
    reference_bases: u64,
    left: usize,
    right: usize,
    consumed: u64,
    walk: Walk,
) -> (u64, u64) {
    let (mut ref_pos, mut query_pos) = (0u64, 0u64);
    for i in left..right {
        let index = match walk {
            Walk::Backward => left + right - i - 1,
            Walk::Forward => i,
        };
        let (o, oplen) = cigar_op(cigar[index]);
        let length = if i == left { consumed } else { oplen };
        if length == 0 {
            continue;
        }
        if op::is_match(o) {
            query_pos = query_pos.wrapping_add(length);
            ref_pos = ref_pos.wrapping_add(length);
            if ref_pos >= reference_bases {
                return (reference_bases, query_pos.wrapping_add(reference_bases).wrapping_sub(ref_pos));
            }
        } else if o == op::DEL {
            ref_pos = ref_pos.wrapping_add(length);
            if ref_pos >= reference_bases {
                return (reference_bases, query_pos);
            }
        } else if o == op::INS {
            query_pos = query_pos.wrapping_add(length);
        } else if o == op::REF_SKIP {
            return (reference_bases, query_pos);
        }
    }
    (0, 0)
}

/// `levenshtein_n` from `levenshtein.c` (Titus Wormer's single-row edit distance).
pub(crate) fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    if a.is_empty() {
        return b.len();
    }
    if b.is_empty() {
        return a.len();
    }
    let mut cache: Vec<usize> = (1..=a.len()).collect();
    let mut result = 0;
    for (b_index, &code) in b.iter().enumerate() {
        let mut distance = b_index;
        result = b_index;
        for (index, &a_code) in a.iter().enumerate() {
            let b_distance = if code == a_code { distance } else { distance + 1 };
            distance = cache[index];
            result = if distance > result {
                if b_distance > result { result + 1 } else { b_distance }
            } else if b_distance > distance {
                distance + 1
            } else {
                b_distance
            };
            cache[index] = result;
        }
    }
    result
}

/// The prefix of `s` before its first NUL (`strlen`).
fn c_str(s: &[u8]) -> &[u8] {
    s.iter().position(|&b| b == 0).map_or(s, |nul| &s[..nul])
}

/// `realign_read`: which allele the read supports at `snp`, whose position lies in CIGAR
/// operation `i`, `consumed` bases into it, at query position `query_pos`. 0: neither (a tie),
/// 1: reference, 2: alternative.
fn realign_read(
    snp: &PhasedSnp,
    read: &ReadView<'_>,
    i: usize,
    consumed: u64,
    query_pos: u64,
    reference: &RefWindow,
) -> i32 {
    let cigar = read.cigar;
    let (_, middle_length) = cigar_op(cigar[i]);
    let left_consumed = consumed;
    let right_consumed = middle_length.saturating_sub(consumed);
    let (left_ref, left_query) = cigar_prefix_length(cigar, OVERHANG, 0, i + 1, left_consumed, Walk::Backward);
    let (right_ref, right_query) =
        cigar_prefix_length(cigar, OVERHANG + 1, i, cigar.len(), right_consumed, Walk::Forward);

    let qst = query_pos.wrapping_sub(left_query);
    let qen = query_pos.wrapping_add(right_query);
    // `variant->position` is an `int`, converted to `size_t` here; `ref_start` is applied by
    // `RefWindow`.
    let position = snp.position as i64 as u64;
    let rst = position.wrapping_sub(left_ref);
    let ren = position.wrapping_add(right_ref);
    if qen == qst {
        return 0;
    }
    let query: Vec<u8> = (qst..qen).map(|q| NT16[usize::from(read.base_code(q))]).collect();
    let ref_seq = reference.c_str(rst, ren.wrapping_sub(rst));
    // `alt` is a zeroed buffer of `alt_length + 2` bytes holding `ref`, with the alternative
    // base written at `left_ref_bases` (past the end of `ref` when no right flank was found).
    let alt_length = left_ref.wrapping_add(right_ref);
    let mut alt = vec![0u8; usize::try_from(alt_length).unwrap_or(0).saturating_add(2)];
    let copied = ref_seq.len().min(alt.len());
    alt[..copied].copy_from_slice(&ref_seq[..copied]);
    if let Some(slot) = usize::try_from(left_ref).ok().and_then(|l| alt.get_mut(l)) {
        *slot = snp.alt_base;
    }
    let distance_ref = levenshtein(&query, ref_seq);
    let distance_alt = levenshtein(&query, c_str(&alt));
    match distance_ref.cmp(&distance_alt) {
        std::cmp::Ordering::Less => 1,
        std::cmp::Ordering::Greater => 2,
        std::cmp::Ordering::Equal => 0,
    }
}

/// `haplotag_read`: the haplotype of `read` from the phased SNPs `snps[first..]` it covers.
/// `first` is the first SNP at or after the read start (`variant_current_pos`).
pub(crate) fn haplotag_read(snps: &[PhasedSnp], first: usize, read: &ReadView<'_>, reference: &RefWindow) -> Haplotype {
    let n = snps.len();
    let position = |j: usize| snps[j].position as i64 as u64;
    let mut j = first;
    let mut ref_pos = read.start;
    let mut query_pos = 0u64;
    let mut v_position = 0u64;
    // Votes per phase set; only the extremes and whether any set was seen matter.
    let mut cost: HashMap<i32, i32> = HashMap::new();
    let mut vote = |allele: i32, snp: &PhasedSnp| {
        if allele == 0 {
            return;
        }
        let delta = if allele == snp.genotype.code() { 1 } else { -1 };
        *cost.entry(snp.phase_set).or_insert(0) += delta;
    };

    while j < n && position(j) < ref_pos {
        j += 1;
    }
    for (i, &c) in read.cigar.iter().enumerate() {
        let (o, length) = cigar_op(c);
        if j < n {
            v_position = position(j);
        }
        if op::is_match(o) {
            while j < n && v_position < ref_pos.wrapping_add(length) {
                let consumed = v_position.wrapping_sub(ref_pos);
                let allele = realign_read(&snps[j], read, i, consumed, query_pos.wrapping_add(consumed), reference);
                vote(allele, &snps[j]);
                j += 1;
                if j < n {
                    v_position = position(j);
                }
            }
            query_pos = query_pos.wrapping_add(length);
            ref_pos = ref_pos.wrapping_add(length);
        } else if o == op::INS {
            if j < n && v_position == ref_pos {
                let allele = realign_read(&snps[j], read, i, 0, query_pos, reference);
                vote(allele, &snps[j]);
                j += 1;
                if j < n {
                    v_position = position(j);
                }
            }
            query_pos = query_pos.wrapping_add(length);
        } else if o == op::DEL {
            while j < n && v_position < ref_pos.wrapping_add(length) {
                let allele = realign_read(&snps[j], read, i, v_position.wrapping_sub(ref_pos), query_pos, reference);
                vote(allele, &snps[j]);
                j += 1;
                if j < n {
                    v_position = position(j);
                }
            }
            ref_pos = ref_pos.wrapping_add(length);
        } else if o == op::REF_SKIP {
            while j < n && v_position < ref_pos.wrapping_add(length) {
                j += 1;
                if j < n {
                    v_position = position(j);
                }
            }
            ref_pos = ref_pos.wrapping_add(length);
        } else if o == op::SOFT_CLIP {
            query_pos = query_pos.wrapping_add(length);
        }
    }

    let max_value = cost.values().fold(0, |m, &v| m.max(v));
    let min_value = cost.values().fold(0, |m, &v| m.min(v));
    if cost.is_empty() || (max_value == 0 && min_value == 0) {
        Haplotype::Unphased
    } else if max_value > min_value.wrapping_abs() {
        Haplotype::First
    } else {
        Haplotype::Second
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inputs::PhasedGenotype;

    fn cig(ops: &[(u32, u32)]) -> Vec<u32> {
        ops.iter().map(|&(o, l)| (l << 4) | o).collect()
    }

    #[test]
    fn edit_distance() {
        assert_eq!(levenshtein(b"kitten", b"sitting"), 3);
        assert_eq!(levenshtein(b"", b"abc"), 3);
        assert_eq!(levenshtein(b"abc", b""), 3);
        assert_eq!(levenshtein(b"ACGT", b"ACGT"), 0);
        assert_eq!(levenshtein(b"ACGT", b"AGT"), 1);
        assert_eq!(levenshtein(b"flaw", b"lawn"), 2);
    }

    #[test]
    fn prefix_lengths() {
        let c = cig(&[(op::SOFT_CLIP, 5), (op::MATCH, 4), (op::INS, 2), (op::MATCH, 20)]);
        // Forward from the last M, 3 bases into it: 17 left, enough for 11.
        assert_eq!(cigar_prefix_length(&c, 11, 3, 4, 17, Walk::Forward), (11, 11));
        // Backward from the last M with 3 consumed: 3 + insertion 2 + 4 more for 7 reference.
        assert_eq!(cigar_prefix_length(&c, 7, 0, 4, 3, Walk::Backward), (7, 9));
        // Not enough reference before the clip: outputs stay zero.
        assert_eq!(cigar_prefix_length(&c, 10, 0, 4, 3, Walk::Backward), (0, 0));
        let skip = cig(&[(op::MATCH, 4), (op::REF_SKIP, 100), (op::MATCH, 4)]);
        // An `N` ends the walk with the full reference length and the query bases so far.
        assert_eq!(cigar_prefix_length(&skip, 10, 0, 3, 2, Walk::Forward), (10, 2));
    }

    #[test]
    fn votes_follow_genotype() {
        // Reference ACGT... with a read carrying G at position 12 (reference C).
        let reference: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let window = RefWindow::new(reference.clone(), 0);
        let mut read_bases = reference.clone();
        read_bases[13] = b'G';
        let codes: Vec<u8> =
            read_bases.iter().map(|b| NT16.iter().position(|n| n == b).expect("nt16 base") as u8).collect();
        let seq: Vec<u8> = codes.chunks(2).map(|p| (p[0] << 4) | p.get(1).copied().unwrap_or(0)).collect();
        let cigar = cig(&[(op::MATCH, 32)]);
        let qual = vec![30; 32];
        let read = ReadView { start: 0, cigar: &cigar, seq: &seq, qual: &qual, len: 32 };
        let snp = |genotype, alt| PhasedSnp { position: 13, ref_base: b'C', alt_base: alt, genotype, phase_set: 1 };
        let alt_second = [snp(PhasedGenotype::AltOnSecond, b'G')];
        assert_eq!(haplotag_read(&alt_second, 0, &read, &window), Haplotype::Second);
        let alt_first = [snp(PhasedGenotype::AltOnFirst, b'G')];
        assert_eq!(haplotag_read(&alt_first, 0, &read, &window), Haplotype::First);
        // A SNP the read does not support either way (a tie) leaves it unphased.
        let other = [snp(PhasedGenotype::AltOnFirst, b'T')];
        assert_eq!(haplotag_read(&other, 0, &read, &window), Haplotype::Unphased);
    }
}

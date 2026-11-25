//! Port of `align.pyx`: `edit_distance` (without banding),
//! `edit_distance_affine_gap` and the kmer alignment of kmerald.

use std::collections::{HashMap, VecDeque};

/// Unit-cost edit distance between `s` and `t`.
pub fn edit_distance(s: &[u8], t: &[u8]) -> i64 {
    let (mut s, mut t) = (s, t);
    // Skip identical prefixes and suffixes
    while !s.is_empty() && !t.is_empty() && s[0] == t[0] {
        s = &s[1..];
        t = &t[1..];
    }
    while !s.is_empty() && !t.is_empty() && s[s.len() - 1] == t[t.len() - 1] {
        s = &s[..s.len() - 1];
        t = &t[..t.len() - 1];
    }
    let m = s.len();
    let mut costs: Vec<i64> = (0..=m as i64).collect();
    for &tc in t {
        let mut prev = costs[0];
        costs[0] += 1;
        for i in 1..=m {
            let matched = (s[i - 1] == tc) as i64;
            let c = (prev + 1 - matched).min(costs[i] + 1).min(costs[i - 1] + 1);
            prev = costs[i];
            costs[i] = c;
        }
    }
    costs[m]
}

/// `edit_distance_affine_gap` (Gotoh): `mismatch_cost` per query position;
/// the DP runs in single precision as in WhatsHap
pub fn edit_distance_affine_gap(
    query: &[u8],
    reference: &[u8],
    mismatch_cost: f32,
    gap_start: i64,
    gap_extend: i64,
) -> i64 {
    let f = |l: i64| (gap_start + (l - 1) * gap_extend) as f32;
    let (mut sv, mut tv) = (query, reference);
    // Skip identical prefixes and suffixes
    while !sv.is_empty() && !tv.is_empty() && sv[0] == tv[0] {
        sv = &sv[1..];
        tv = &tv[1..];
    }
    while !sv.is_empty() && !tv.is_empty() && sv[sv.len() - 1] == tv[tv.len() - 1] {
        sv = &sv[..sv.len() - 1];
        tv = &tv[..tv.len() - 1];
    }
    let (m, n) = (sv.len(), tv.len());
    let int_max = i32::MAX as f32;
    let (gs, ge) = (gap_start as f32, gap_extend as f32);
    let mut a = vec![int_max; m + 1];
    let mut b: Vec<f32> = (0..=m).map(|i| if i == 0 { 0.0 } else { f(i as i64) }).collect();
    let mut c = vec![int_max; m + 1];
    a[0] = 0.0;
    c[0] = 0.0;
    for j in 1..=n {
        let (mut prev_a, mut prev_b, mut prev_c) = (a[0], b[0], c[0]);
        a[0] = int_max;
        b[0] = int_max;
        c[0] = f(j as i64);
        for i in 1..=m {
            let m_c = if sv[i - 1] == tv[j - 1] { 0.0 } else { mismatch_cost };
            let c_a = prev_a.min(prev_b).min(prev_c) + m_c;
            let c_b = (a[i - 1] + gs).min(b[i - 1] + ge).min(c[i - 1] + gs);
            let c_c = (a[i] + gs).min(b[i] + gs).min(c[i] + ge);
            prev_a = a[i];
            prev_b = b[i];
            prev_c = c[i];
            a[i] = c_a;
            b[i] = c_b;
            c[i] = c_c;
        }
    }
    a[m].min(b[m]).min(c[m]) as i64
}

/// `enumerate_all_kmers`: the 2-bit encoded k-mers (A=0, C=1, G=2, T=3;
/// other characters leave the k-mer unchanged), skipping all-A k-mers
pub fn enumerate_all_kmers(reference: &[u8], k: u32) -> VecDeque<i32> {
    let mask: i32 = ((1i64 << (2 * k)) - 1) as i32;
    let mut h: i32 = 0;
    let mut kmers = VecDeque::new();
    for (i, &c) in reference.iter().enumerate() {
        let code = match c {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        if let Some(code) = code {
            h = (h.wrapping_shl(2) | code) & mask;
        }
        if i as i64 >= k as i64 - 1 && h != 0 {
            kmers.push_back(h);
        }
    }
    kmers
}

/// `kmer_align`: alignment cost of two k-mer sequences with the error model
/// `costs` ((kmer, kmer) -> cost text; (kmer, -5): any other k-mer)
pub fn kmer_align(
    seq1: &VecDeque<i32>,
    seq2: &VecDeque<i32>,
    costs: &HashMap<(i64, i64), f32>,
    gap_penalty: f32,
) -> f32 {
    if seq1 == seq2 {
        return 0.0;
    }
    let (mut m, mut n) = (seq1.len(), seq2.len());
    let mut x = 0;
    while x < m && x < n && seq1[x] == seq2[x] {
        x += 1;
    }
    while m > x && n > x && seq1[m - 1] == seq2[n - 1] {
        m -= 1;
        n -= 1;
    }
    m -= x;
    n -= x;
    let mut score = vec![vec![0f32; n + 1]; m + 1];
    for (i, row) in score.iter_mut().enumerate() {
        row[0] = gap_penalty * i as f32;
    }
    for (j, s) in score[0].iter_mut().enumerate() {
        *s = gap_penalty * j as f32;
    }
    for i in 1..=m {
        for j in 1..=n {
            let (k1, k2) = (seq1[i - 1 + x] as i64, seq2[j - 1 + x] as i64);
            let matched = if k1 == k2 {
                score[i - 1][j - 1]
            } else {
                let mismatching =
                    costs.get(&(k1, k2)).or_else(|| costs.get(&(k1, -5))).copied().unwrap_or(f32::INFINITY);
                score[i - 1][j - 1] + mismatching
            };
            let delete = score[i - 1][j] + gap_penalty;
            let insert = score[i][j - 1] + gap_penalty;
            // Python's min: the first of equal values
            let mut best = matched;
            if delete < best {
                best = delete;
            }
            if insert < best {
                best = insert;
            }
            score[i][j] = best;
        }
    }
    score[m][n]
}

#[cfg(test)]
mod tests {
    use super::{edit_distance, edit_distance_affine_gap};

    #[test]
    fn affine() {
        assert_eq!(edit_distance_affine_gap(b"ACGT", b"ACGT", 15.0, 10, 7), 0);
        assert_eq!(edit_distance_affine_gap(b"ACGT", b"AGGT", 15.0, 10, 7), 15);
        assert_eq!(edit_distance_affine_gap(b"ACGGGT", b"ACT", 15.0, 10, 7), 24);
    }

    #[test]
    fn basic() {
        assert_eq!(edit_distance(b"", b""), 0);
        assert_eq!(edit_distance(b"ACGT", b"ACGT"), 0);
        assert_eq!(edit_distance(b"ACGT", b"AGT"), 1);
        assert_eq!(edit_distance(b"kitten", b"sitting"), 3);
        assert_eq!(edit_distance(b"", b"ACG"), 3);
    }
}

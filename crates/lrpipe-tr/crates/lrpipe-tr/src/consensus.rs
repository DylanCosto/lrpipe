//! An allele's sequence from its reads' sequences.
//!
//! When one sequence is shared by most reads it is the allele. Otherwise every read is aligned
//! to a backbone (banded global alignment by edit distance) and each backbone base, deletion
//! and insertion is put to a vote; the result becomes the backbone of a second round. Ties in
//! the alignment always place gaps the same way (as far right as possible), so reads agree on
//! where an indel in a homopolymer goes and their votes add up.

use std::collections::HashMap;

/// Extra diagonals the band keeps beyond the length difference.
const BAND_SLACK: usize = 16;
/// Voting rounds, each on the previous round's consensus.
const ROUNDS: usize = 2;

/// The consensus of `seqs`, starting from `backbone` when no sequence has a majority.
///
/// # Panics
///
/// If `seqs` is empty.
pub fn consensus(seqs: &[&[u8]], backbone: &[u8]) -> Vec<u8> {
    assert!(!seqs.is_empty(), "consensus of no sequences");
    let mut counts: HashMap<&[u8], usize> = HashMap::new();
    for s in seqs {
        *counts.entry(s).or_default() += 1;
    }
    if let Some((s, _)) = counts.iter().find(|&(_, &n)| n * 2 > seqs.len()) {
        return s.to_vec();
    }
    let mut current = backbone.to_vec();
    for _ in 0..ROUNDS {
        current = vote(seqs, &current);
    }
    current
}

/// Votes at one backbone position: bases A, C, G, T, other, and deletion.
#[derive(Clone, Copy, Default)]
struct Column([u32; 6]);

fn base_index(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// One round: align every read to `backbone`, keep each position's most voted base (or drop
/// it when deletion wins) and insert, before each position, the most common insertion when
/// most reads have one there.
fn vote(seqs: &[&[u8]], backbone: &[u8]) -> Vec<u8> {
    let m = backbone.len();
    let mut columns = vec![Column::default(); m];
    let mut inserts: Vec<HashMap<Vec<u8>, u32>> = vec![HashMap::new(); m + 1];
    let mut with_insert = vec![0u32; m + 1];
    for read in seqs {
        for (gap, ins) in align(backbone, read, &mut columns) {
            with_insert[gap] += 1;
            *inserts[gap].entry(ins).or_default() += 1;
        }
    }
    let n = seqs.len() as u32;
    let mut out = Vec::with_capacity(m + m / 10);
    for gap in 0..=m {
        if with_insert[gap] * 2 > n
            && let Some((ins, _)) = inserts[gap].iter().max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
        {
            out.extend_from_slice(ins);
        }
        if gap < m {
            let c = &columns[gap].0;
            let best =
                (0..6).max_by_key(|&k| (c[k], usize::from(k == base_index(backbone[gap])))).expect("six choices");
            match best {
                0..=3 => out.push(b"ACGT"[best]),
                4 => out.push(backbone[gap]),
                _ => {}
            }
        }
    }
    out
}

/// The cost table of a banded global alignment of a read to a backbone: row i (backbone
/// position) keeps read positions lo[i]..=hi[i] around the diagonal i * n / m, stored from
/// start[i].
struct Band {
    lo: Vec<usize>,
    hi: Vec<usize>,
    start: Vec<usize>,
    cost: Vec<u32>,
}

impl Band {
    fn new(backbone: &[u8], read: &[u8]) -> Self {
        let (m, n) = (backbone.len(), read.len());
        let w = m.abs_diff(n) + BAND_SLACK;
        let center = |i: usize| (i * n).checked_div(m).unwrap_or(0);
        let lo: Vec<usize> = (0..=m).map(|i| center(i).saturating_sub(w)).collect();
        let hi: Vec<usize> = (0..=m).map(|i| (center(i) + w).min(n)).collect();
        let mut start = vec![0usize; m + 2];
        for i in 0..=m {
            start[i + 1] = start[i] + hi[i] - lo[i] + 1;
        }
        const INF: u32 = u32::MAX / 2;
        let mut cost = vec![INF; start[m + 1]];
        // Row 0 (lo[0] is 0): only gaps in the backbone.
        for (j, c) in cost.iter_mut().enumerate().take(hi[0] + 1) {
            *c = j as u32;
        }
        // Each cell from its diagonal, upper and left neighbours where they lie in the band
        // (bands only move right, so the previous row's covers j - 1 and j or ends before).
        for i in 1..=m {
            let (pl, ph, ps) = (lo[i - 1], hi[i - 1], start[i - 1]);
            let (l, h, s) = (lo[i], hi[i], start[i]);
            let base = backbone[i - 1];
            let cell = |cost: &[u32], j: usize| -> u32 {
                let mut best = INF;
                if j > pl && j - 1 <= ph {
                    best = cost[ps + j - 1 - pl] + u32::from(base != read[j - 1]);
                }
                if j >= pl && j <= ph {
                    best = best.min(cost[ps + j - pl] + 1);
                }
                if j > l {
                    best = best.min(cost[s + j - 1 - l] + 1);
                }
                best
            };
            // Where the diagonal, upper and left neighbours all lie in the band, the same
            // minimum without the checks.
            let (mid_lo, mid_hi) = ((l + 1).max(pl + 1), h.min(ph));
            if mid_lo > mid_hi {
                for j in l..=h {
                    cost[s + j - l] = cell(&cost, j);
                }
                continue;
            }
            for j in l..mid_lo {
                cost[s + j - l] = cell(&cost, j);
            }
            let (prev, row) = cost.split_at_mut(s);
            let prev = &prev[ps..];
            let mut left = row[mid_lo - 1 - l];
            for j in mid_lo..=mid_hi {
                let diagonal = prev[j - 1 - pl] + u32::from(base != read[j - 1]);
                let up = prev[j - pl] + 1;
                left = diagonal.min(up).min(left + 1);
                row[j - l] = left;
            }
            for j in mid_hi + 1..=h {
                cost[s + j - l] = cell(&cost, j);
            }
        }
        Self { lo, hi, start, cost }
    }

    /// Where cell (i, j) is stored, if it lies in the band.
    fn at(&self, i: usize, j: usize) -> Option<usize> {
        (j >= self.lo[i] && j <= self.hi[i]).then(|| self.start[i] + j - self.lo[i])
    }
}

/// Aligns `read` to `backbone` (banded global edit distance), adds the read's base or
/// deletion at each backbone position to `columns`, and returns its insertions as (gap before
/// backbone position, inserted bases).
fn align(backbone: &[u8], read: &[u8], columns: &mut [Column]) -> Vec<(usize, Vec<u8>)> {
    let (m, n) = (backbone.len(), read.len());
    let band = Band::new(backbone, read);
    let (cost, at) = (&band.cost, |i, j| band.at(i, j));
    // Traceback from the end, gaps before matches on ties (gaps go right).
    let (mut i, mut j) = (m, n);
    let mut insertions: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut pending: Vec<u8> = Vec::new();
    while i > 0 || j > 0 {
        let here = cost[at(i, j).expect("the path stays in the band")];
        let ins = (j > 0).then(|| at(i, j - 1)).flatten().filter(|&l| cost[l] + 1 == here);
        let del = (i > 0).then(|| at(i - 1, j)).flatten().filter(|&u| cost[u] + 1 == here);
        if ins.is_some() {
            pending.push(read[j - 1]);
            j -= 1;
            continue;
        }
        if !pending.is_empty() {
            pending.reverse();
            insertions.push((i, std::mem::take(&mut pending)));
        }
        if del.is_some() {
            columns[i - 1].0[5] += 1;
            i -= 1;
        } else {
            columns[i - 1].0[base_index(read[j - 1])] += 1;
            i -= 1;
            j -= 1;
        }
    }
    if !pending.is_empty() {
        pending.reverse();
        insertions.push((0, pending));
    }
    insertions
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cost table filled through band lookups for every neighbour, as first written.
    fn reference_costs(backbone: &[u8], read: &[u8]) -> Vec<u32> {
        let (m, n) = (backbone.len(), read.len());
        let w = m.abs_diff(n) + BAND_SLACK;
        let center = |i: usize| (i * n).checked_div(m).unwrap_or(0);
        let lo: Vec<usize> = (0..=m).map(|i| center(i).saturating_sub(w)).collect();
        let hi: Vec<usize> = (0..=m).map(|i| (center(i) + w).min(n)).collect();
        let mut start = vec![0usize; m + 2];
        for i in 0..=m {
            start[i + 1] = start[i] + hi[i] - lo[i] + 1;
        }
        const INF: u32 = u32::MAX / 2;
        let mut cost = vec![INF; start[m + 1]];
        let at = |i: usize, j: usize| (j >= lo[i] && j <= hi[i]).then(|| start[i] + j - lo[i]);
        for i in 0..=m {
            for j in lo[i]..=hi[i] {
                let here = start[i] + j - lo[i];
                if i == 0 && j == 0 {
                    cost[here] = 0;
                    continue;
                }
                let mut best = INF;
                if i > 0
                    && j > 0
                    && let Some(d) = at(i - 1, j - 1)
                {
                    best = best.min(cost[d] + u32::from(backbone[i - 1] != read[j - 1]));
                }
                if i > 0
                    && let Some(u) = at(i - 1, j)
                {
                    best = best.min(cost[u] + 1);
                }
                if j > 0
                    && let Some(l) = at(i, j - 1)
                {
                    best = best.min(cost[l] + 1);
                }
                cost[here] = best;
            }
        }
        cost
    }

    #[test]
    fn direct_band_indexing_changes_no_cost() {
        let seqs = crate::align::tests::sequences(11, 300, 150);
        for pair in seqs.chunks(2) {
            let (a, b) = (&pair[0], &pair[1]);
            assert_eq!(Band::new(a, b).cost, reference_costs(a, b), "{a:?} {b:?}");
        }
    }

    #[test]
    fn majority_sequence_wins() {
        let s: Vec<&[u8]> = vec![b"CAGCAG", b"CAGCAG", b"CAGAG"];
        assert_eq!(consensus(&s, b"CAGAG"), b"CAGCAG");
    }

    #[test]
    fn votes_out_independent_errors() {
        // Every read has one different error; the truth is TTTTTGACGTAAAAAC.
        let truth: &[u8] = b"TTTTTGACGTAAAAAC";
        let reads: Vec<&[u8]> = vec![
            b"TTTTGACGTAAAAAC",   // homopolymer deletion
            b"TTTTTGACGTAAAAAAC", // homopolymer insertion
            b"TTTTTGACCTAAAAAC",  // substitution
            b"TTTTTGACGTAAAAC",   // deletion
            b"TTTTTTGACGTAAAAAC", // insertion
        ];
        assert_eq!(consensus(&reads, reads[0]), truth);
    }

    #[test]
    fn long_alleles_keep_their_length() {
        let truth: Vec<u8> = b"ACGTTGCAGGT".repeat(40);
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for k in 0..7 {
            let mut r = truth.clone();
            // A one-base deletion at a different place in each read.
            r.remove(30 + 50 * k);
            if k % 2 == 0 {
                r.insert(10 + 40 * k, b'A');
            }
            reads.push(r);
        }
        let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        assert_eq!(consensus(&refs, refs[1]), truth);
    }
}

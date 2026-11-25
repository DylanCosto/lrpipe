//! Fitting alignment of a flank into a stretch of read: the whole flank against any substring
//! of the read, by edit distance. Flanks are short (tens of bases) and windows small, so a
//! plain dynamic program is enough. Also the edit distance of two whole sequences, which can
//! be thousands of bases long, bit-parallel.

/// Where `pattern` aligns best inside `text`: the end (exclusive) of the best match in
/// `text`, and its edit distance. Among equally good ends, the one nearest `expected`.
pub fn best_end(pattern: &[u8], text: &[u8], expected: usize) -> (usize, u32) {
    // prev[j]: cost of aligning pattern[..i] ending at text[..j], free start in text.
    let mut prev: Vec<u32> = vec![0; text.len() + 1];
    let mut cur: Vec<u32> = vec![0; text.len() + 1];
    for (i, &p) in pattern.iter().enumerate() {
        cur[0] = i as u32 + 1;
        for j in 1..=text.len() {
            let sub = prev[j - 1] + u32::from(p != text[j - 1]);
            cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    let (j, c) = prev.iter().enumerate().min_by_key(|&(j, &c)| (c, j.abs_diff(expected))).expect("at least one end");
    (j, *c)
}

/// Where `pattern` aligns best inside `text`: the start of the best match, and its edit
/// distance. Among equally good starts, the one nearest `expected`.
pub fn best_start(pattern: &[u8], text: &[u8], expected: usize) -> (usize, u32) {
    let p: Vec<u8> = pattern.iter().rev().copied().collect();
    let t: Vec<u8> = text.iter().rev().copied().collect();
    let (end, cost) = best_end(&p, &t, text.len().saturating_sub(expected));
    (text.len() - end, cost)
}

/// Edit distance between `a` and `b`, by Myers' bit-parallel algorithm in 64-row blocks
/// (Hyyrö's block formulation, as in edlib): about m * n / 64 word steps.
pub fn edit_distance(a: &[u8], b: &[u8]) -> u32 {
    let m = a.len();
    if m == 0 {
        return b.len() as u32;
    }
    let blocks = m.div_ceil(64);
    // Match masks: for each byte of `a`, the rows (bits) it sits in. Rows past `m` match
    // nothing.
    let mut class = [u8::MAX; 256];
    let mut peq: Vec<u64> = Vec::new();
    for (i, &c) in a.iter().enumerate() {
        if class[usize::from(c)] == u8::MAX {
            class[usize::from(c)] = u8::try_from(peq.len() / blocks).expect("at most 255 distinct bytes");
            peq.extend(std::iter::repeat_n(0, blocks));
        }
        peq[usize::from(class[usize::from(c)]) * blocks + i / 64] |= 1 << (i % 64);
    }
    let mut pv = vec![u64::MAX; blocks];
    let mut mv = vec![0u64; blocks];
    // D at the last row of the last (padded) block, which starts at its row count.
    let mut score = (blocks * 64) as i64;
    for &c in b {
        let k = class[usize::from(c)];
        let mut hin: i64 = 1;
        for blk in 0..blocks {
            let eq = if k == u8::MAX { 0 } else { peq[usize::from(k) * blocks + blk] };
            let (p, mm) = (pv[blk], mv[blk]);
            let hin_neg = u64::from(hin < 0);
            let xv = eq | mm;
            let eq = eq | hin_neg;
            let xh = (((eq & p).wrapping_add(p)) ^ p) | eq;
            let mut ph = mm | !(xh | p);
            let mut mh = p & xh;
            let hout = i64::from(ph >> 63 != 0) - i64::from(mh >> 63 != 0);
            ph = (ph << 1) | u64::from(hin > 0);
            mh = (mh << 1) | hin_neg;
            pv[blk] = mh | !(xv | ph);
            mv[blk] = ph & xv;
            hin = hout;
        }
        score += hin;
    }
    // Back up from the padded last row to row m by the vertical differences in between.
    let pad = blocks * 64 - m;
    if pad > 0 {
        let tail = u64::MAX << (64 - pad);
        let (p, mm) = (pv[blocks - 1], mv[blocks - 1]);
        score -= i64::from((p & tail).count_ones()) - i64::from((mm & tail).count_ones());
    }
    score as u32
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    #[test]
    fn fits_flanks_with_errors() {
        assert_eq!(best_end(b"ACGT", b"TTACGTCAG", 0), (6, 0));
        assert_eq!(best_end(b"ACGTA", b"TTACTACAG", 0), (6, 1));
        assert_eq!(best_start(b"GGCA", b"CAGCAGGGCATT", 0), (6, 0));
        // Repetitive text: the match nearest the expected position.
        assert_eq!(best_end(b"CACA", b"CACACACACA", 8), (8, 0));
        assert_eq!(edit_distance(b"CAGCAGCAG", b"CAGCAGCAG"), 0);
        assert_eq!(edit_distance(b"CAGCAGCAG", b"CAGCAG"), 3);
        assert_eq!(edit_distance(b"AAAAT", b"AAAT"), 1);
    }

    /// A small deterministic generator (xorshift), for sequences.
    pub(crate) fn sequences(seed: u64, count: usize, max_len: usize) -> Vec<Vec<u8>> {
        let mut x = seed.max(1);
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        (0..count)
            .map(|_| {
                let len = (next() % (max_len as u64 + 1)) as usize;
                // Mostly a repeat with errors, as reads of a tandem repeat are.
                (0..len).map(|k| if next() % 8 == 0 { b"ACGT"[(next() % 4) as usize] } else { b"CAG"[k % 3] }).collect()
            })
            .collect()
    }

    /// Plain edit distance.
    fn plain_distance(a: &[u8], b: &[u8]) -> u32 {
        let mut prev: Vec<u32> = (0..=b.len() as u32).collect();
        for (i, &x) in a.iter().enumerate() {
            let mut cur = vec![i as u32 + 1; b.len() + 1];
            for (j, &y) in b.iter().enumerate() {
                cur[j + 1] = (prev[j] + u32::from(x != y)).min(prev[j + 1] + 1).min(cur[j] + 1);
            }
            prev = cur;
        }
        prev[b.len()]
    }

    #[test]
    fn bit_parallel_distance_is_the_edit_distance() {
        for (a, b) in [(&b""[..], &b"ACG"[..]), (b"ACG", b""), (b"N", b"A"), (b"ACGTN", b"ACGTN")] {
            assert_eq!(edit_distance(a, b), plain_distance(a, b), "{a:?} {b:?}");
        }
        let seqs = sequences(3, 600, 300);
        for pair in seqs.chunks(2) {
            let (a, b) = (&pair[0], &pair[1]);
            assert_eq!(edit_distance(a, b), plain_distance(a, b), "{} {}", a.len(), b.len());
        }
        // Block edges.
        for len in [63, 64, 65, 127, 128, 129] {
            let a: Vec<u8> = b"CAG".iter().cycle().take(len).copied().collect();
            let mut b = a.clone();
            b.insert(len / 2, b'T');
            b.pop();
            assert_eq!(edit_distance(&a, &b), plain_distance(&a, &b), "{len}");
        }
    }
}

//! A read's sequence of a repeat locus.
//!
//! The locus boundaries are mapped into the read through its alignment: an insertion right
//! at a boundary counts as part of the repeat. Where the alignment has an indel near a
//! boundary on the flank side, where the aligner may have placed repeat bases in the flank or
//! flank bases in the repeat, the boundary is found again by aligning that flank's
//! reference bases into the read.
//!
//! An expansion longer than the aligner aligns through is usually written as an alignment that
//! stops in the repeat, the rest of the read soft-clipped. Such a read counts when its clipped
//! bases hold the far flank: its allele runs from the aligned boundary to where that flank fits.
//!
//! A boundary the locus shares with an adjacent catalog locus is not refined (its flank is the
//! other repeat): an insertion at it, or near it and sliding to it (an aligner places an
//! insertion anywhere along a stretch of its period, as FXN's GAA expansion two bases into the
//! A tail before it), is split between the two loci where it best continues each one's repeat,
//! so each base of it counts in one of them only.

use crate::align::{best_end, best_start};
use crate::catalog::{Locus, SharedBoundary};
use crate::reads::{Op, Read};

/// Settings of extraction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExtractOptions {
    /// Reads must align at least this far past both boundaries.
    pub flank: i64,
    /// An indel within this many reference bases of a boundary, on the flank side, makes the
    /// boundary be found by flank alignment.
    pub check: i64,
    /// Flank bases aligned when refining a boundary.
    pub anchor: usize,
    /// Extra read bases searched on either side of the alignment's boundary.
    pub search: usize,
}

impl Default for ExtractOptions {
    fn default() -> Self {
        ExtractOptions { flank: 50, check: 10, anchor: 30, search: 40 }
    }
}

/// Reference bases `[start, start + seq.len())` of a contig.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RefWindow {
    /// 0-based position of the first base.
    pub start: i64,
    /// Uppercase bases.
    pub seq: Vec<u8>,
}

impl RefWindow {
    /// The whole of `contig` from an indexed FASTA, uppercase.
    ///
    /// # Errors
    ///
    /// If the FASTA or its index cannot be read or lacks the contig.
    pub fn contig(genome: &std::path::Path, contig: &str) -> Result<Self, crate::Error> {
        // Bgzipped or not by its contents, not its name.
        let mut fa = lrpipe_reads::indexed_fasta(genome).map_err(|e| crate::Error::io(genome)(e.into_io_error()))?;
        let rec = fa
            .query(&pasta::core::Region::new(contig, ..))
            .map_err(|e| crate::Error::Input(format!("reference {contig}: {e}")))?;
        Ok(RefWindow { start: 0, seq: rec.sequence().as_ref().to_ascii_uppercase() })
    }

    /// Bases `[from, to)`, clipped to the window.
    pub fn slice(&self, from: i64, to: i64) -> &[u8] {
        let n = self.seq.len() as i64;
        let (a, b) = ((from - self.start).clamp(0, n), (to - self.start).clamp(0, n));
        &self.seq[a as usize..b.max(a) as usize]
    }
}

/// Where reference position `p` falls in the read: the query offset before and after any
/// insertion at `p`, and whether an indel touches the reference positions in `near`.
struct Mapped {
    early: usize,
    late: usize,
}

fn map(read: &Read, p: i64) -> Option<Mapped> {
    let (mut r, mut q) = (read.start, 0usize);
    let (mut early, mut late) = (None, None);
    for &(op, len) in &read.cigar {
        if r == p {
            early.get_or_insert(q);
            late = Some(q);
        }
        let l = i64::from(len);
        match op {
            Op::SoftClip | Op::Ins => q += len as usize,
            Op::Match => {
                if r < p && p < r + l {
                    let at = q + (p - r) as usize;
                    return Some(Mapped { early: at, late: at });
                }
                r += l;
                q += len as usize;
            }
            Op::Del => {
                if r < p && p < r + l {
                    return Some(Mapped { early: q, late: q });
                }
                r += l;
            }
        }
    }
    if r == p {
        early.get_or_insert(q);
        late = Some(q);
    }
    Some(Mapped { early: early?, late: late? })
}

/// Whether the alignment has an indel at reference positions `[from, to]`.
fn indel_near(read: &Read, from: i64, to: i64) -> bool {
    let mut r = read.start;
    for &(op, len) in &read.cigar {
        let l = i64::from(len);
        match op {
            Op::Ins if r >= from && r <= to => return true,
            Op::Del if r <= to && r + l >= from => return true,
            Op::Match | Op::Del => r += l,
            _ => {}
        }
        if r > to {
            break;
        }
    }
    false
}

/// The read's bases between the boundaries of `locus`, or `None` when the read does not span
/// the locus with both flanks (aligned, or one of them in a soft clip) or a flank does not align.
pub fn allele(read: &Read, locus: &Locus, reference: &RefWindow, opts: &ExtractOptions) -> Option<Vec<u8>> {
    let (start, end) = (locus.start, locus.end);
    let (left, right) = (read.start <= start - opts.flank, read.end >= end + opts.flank);
    if !(left && right) {
        return clipped_allele(read, locus, reference, opts);
    }
    let shared = |p: i64, b: &SharedBoundary| -> Option<usize> {
        let at = at_shared(read, p, opts.check)?;
        Some(at.early + left_share(read.seq.get(at.early..at.late)?, p, b, reference))
    };
    let mut qs = match &locus.shared_start {
        Some(b) => shared(start, b)?,
        None => map(read, start)?.early,
    };
    let mut qe = match &locus.shared_end {
        Some(b) => shared(end, b)?,
        None => map(read, end)?.late,
    };
    let (a, w) = (opts.anchor, opts.search);
    let max_cost = (a / 5) as u32;
    // The windows reach from where the checked flank stretch maps in the read, so they span
    // any indel found there, however long.
    if locus.shared_start.is_none() && indel_near(read, start - opts.check, start) {
        let anchor = reference.slice(start - a as i64, start);
        let from = map(read, start - opts.check - a as i64)?.early.saturating_sub(w).min(qs);
        let text = &read.seq[from..(qs + w).min(read.seq.len())];
        let (e, cost) = best_end(anchor, text, qs - from);
        if cost > max_cost {
            return None;
        }
        qs = from + e;
    }
    if locus.shared_end.is_none() && indel_near(read, end, end + opts.check) {
        let anchor = reference.slice(end, end + a as i64);
        let from = qe.saturating_sub(w).max(qs);
        let to = map(read, end + opts.check + a as i64)?.late + w;
        let text = &read.seq[from.min(read.seq.len())..to.max(qe).min(read.seq.len())];
        let (s, cost) = best_start(anchor, text, qe.saturating_sub(from));
        if cost > max_cost {
            return None;
        }
        qe = from + s;
    }
    Some(read.seq.get(qs..qe.max(qs))?.to_vec())
}

/// The soft-clipped bases at the start and at the end of the read's alignment.
pub fn soft_clips(cigar: &[(Op, u32)]) -> [usize; 2] {
    let clip = |op: Option<&(Op, u32)>| match op {
        Some(&(Op::SoftClip, len)) => len as usize,
        _ => 0,
    };
    [clip(cigar.first()), clip(cigar.last())]
}

/// Whether a read aligned over `[start, end)` with `clips` soft-clipped bases could give an allele
/// of locus `[locus_start, locus_end)`: it spans it with `flank` bases on both sides, or spans one
/// flank and reaches the repeat with at least `flank` bases clipped beyond it.
pub fn may_reach(start: i64, end: i64, clips: [usize; 2], locus_start: i64, locus_end: i64, flank: i64) -> bool {
    let (left, right) = (start <= locus_start - flank, end >= locus_end + flank);
    let long = |clip: usize| i64::try_from(clip).is_ok_and(|c| c >= flank);
    (left && right) || (left && long(clips[1]) && end >= locus_start) || (right && long(clips[0]) && start <= locus_end)
}

/// The allele of a read that spans one flank and stops in the repeat, the rest soft-clipped,
/// when the clipped bases hold the other flank. The flank is fitted as `2 * flank` reference
/// bases at an edit distance of at most a fifth of their length: first the bases next to the
/// repeat, then (where the expanded copy differs there, as it often does, its flank part of the
/// repeat's own instability) the same length one and two lengths further out, the allele then
/// taking the skipped bases to be the reference's. A fit this long cannot sit inside the repeat,
/// as a short one of a GC-rich flank can in a CGG repeat.
fn clipped_allele(read: &Read, locus: &Locus, reference: &RefWindow, opts: &ExtractOptions) -> Option<Vec<u8>> {
    let (start, end) = (locus.start, locus.end);
    let clips = soft_clips(&read.cigar);
    if !may_reach(read.start, read.end, clips, start, end, opts.flank) {
        return None;
    }
    let len = 2 * opts.flank;
    let n = usize::try_from(len).ok()?;
    let max_cost = u32::try_from(n / 5).ok()?;
    let w = usize::try_from(opts.flank).ok()?;
    let periods: Vec<usize> = locus.motifs.iter().map(Vec::len).collect();
    if read.start <= start - opts.flank {
        // Anchored on the left flank, the rest clipped: the right flank among the clipped bases,
        // when they start as the repeat (a clip of another kind, a chimera's other part, is not
        // searched: a long flank can fit by chance far into it, in an Alu or L1).
        let qs = map(read, start)?.early;
        let text = read.seq.get(qs..)?;
        let clip_start = read.seq.len() - clips[1];
        if !repeats(read.seq.get(clip_start..clip_start + w)?, &periods) {
            return None;
        }
        for skip in [0, len, 2 * len] {
            let anchor = reference.slice(end + skip, end + skip + len);
            if anchor.len() < n || text.len() < n {
                return None;
            }
            let (s, cost) = best_start(anchor, text, 0);
            if cost <= max_cost {
                let s = s.checked_sub(usize::try_from(skip).ok()?)?;
                return Some(text[..s].to_vec());
            }
        }
        return None;
    }
    // Anchored on the right flank, the start clipped: the left flank among the clipped bases,
    // when they end as the repeat.
    let qe = map(read, end)?.late;
    let text = read.seq.get(..qe)?;
    if !repeats(read.seq.get(clips[0].checked_sub(w)?..clips[0])?, &periods) {
        return None;
    }
    for skip in [0, len, 2 * len] {
        let anchor = reference.slice(start - skip - len, start - skip);
        if anchor.len() < n || text.len() < n {
            return None;
        }
        let (e, cost) = best_end(anchor, text, text.len());
        if cost <= max_cost {
            let e = e + usize::try_from(skip).ok()?;
            return text.get(e..).map(<[u8]>::to_vec);
        }
    }
    None
}

/// Whether `bases` repeat with one of the `periods` (a motif's length): at least three in four
/// of them equal the base one period before.
fn repeats(bases: &[u8], periods: &[usize]) -> bool {
    periods.iter().any(|&k| {
        let pairs = bases.len().saturating_sub(k);
        k > 0 && pairs > 0 && 4 * bases.iter().zip(&bases[k..]).filter(|(a, b)| a == b).count() >= 3 * pairs
    })
}

/// Where shared boundary `p` falls in the read, with an insertion at it, or within `check`
/// reference bases of it and sliding to it, between `early` and `late`. An insertion of `n`
/// bases at `r` slides to `p` when only aligned bases lie between and the read repeats itself
/// with period `n` across them; it is then the same read with the insertion's rotation at `p`.
fn at_shared(read: &Read, p: i64, check: i64) -> Option<Mapped> {
    let at = map(read, p)?;
    if at.late > at.early {
        return Some(at);
    }
    let periodic = |from: usize, len: usize, n: usize| {
        (0..len).all(|k| matches!((read.seq.get(from + k), read.seq.get(from + n + k)), (Some(a), Some(b)) if a == b))
    };
    let (mut r, mut q) = (read.start, 0usize);
    for &(op, len) in &read.cigar {
        let n = len as usize;
        match op {
            Op::Ins if r != p && (r - p).abs() <= check => {
                let between = (r - p).unsigned_abs() as usize;
                if r < p && at.early == q + n + between && periodic(q, between, n) {
                    return Some(Mapped { early: at.early - n, late: at.early });
                }
                if r > p && q == at.early + between && periodic(at.early, between, n) {
                    return Some(Mapped { early: at.early, late: at.early + n });
                }
                q += n;
            }
            Op::SoftClip | Op::Ins => q += n,
            Op::Match => {
                r += i64::from(len);
                q += n;
            }
            Op::Del => r += i64::from(len),
        }
        if r > p + check {
            break;
        }
    }
    Some(at)
}

/// How many leading bases of `ins`, a read's insertion at the shared boundary `p`, belong to
/// the loci ending there: the split at which the most bases continue their side's repeat. The
/// left side's bases are compared with its last reference unit repeated forward, the right
/// side's with its first unit repeated backward, for each motif length of that side's loci.
/// Ties go to the right side, so the split is the same from both loci.
fn left_share(ins: &[u8], p: i64, boundary: &SharedBoundary, reference: &RefWindow) -> usize {
    let n = ins.len();
    // left[k]: of ins[..k], the most bases that continue a left repeat; right[k]: of ins[k..],
    // the most that lead into a right repeat.
    let (mut left, mut right) = (vec![0usize; n + 1], vec![0usize; n + 1]);
    for &period in &boundary.left {
        let unit = reference.slice(p - period as i64, p);
        if unit.len() != period || period == 0 {
            continue;
        }
        let mut matched = 0;
        for (k, &base) in ins.iter().enumerate() {
            matched += usize::from(base == unit[k % period]);
            left[k + 1] = left[k + 1].max(matched);
        }
    }
    for &period in &boundary.right {
        let unit = reference.slice(p, p + period as i64);
        if unit.len() != period || period == 0 {
            continue;
        }
        let mut matched = 0;
        for (k, &base) in ins.iter().enumerate().rev() {
            matched += usize::from(base == unit[period - 1 - (n - 1 - k) % period]);
            right[k] = right[k].max(matched);
        }
    }
    (0..=n).max_by_key(|&k| (left[k] + right[k], std::cmp::Reverse(k))).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locus(start: i64, end: i64) -> Locus {
        Locus {
            id: "t".into(),
            contig: "chr1".into(),
            start,
            end,
            motifs: vec![b"CAG".to_vec()],
            info: String::new(),
            shared_start: None,
            shared_end: None,
        }
    }

    fn read(start: i64, cigar: Vec<(Op, u32)>, seq: &[u8]) -> Read {
        let end = start
            + cigar.iter().filter(|(o, _)| matches!(o, Op::Match | Op::Del)).map(|(_, l)| i64::from(*l)).sum::<i64>();
        Read { start, end, cigar, seq: seq.to_vec(), mapq: 60, ..Read::default() }
    }

    #[test]
    fn insertion_at_the_boundary_belongs_to_the_repeat() {
        // Reference: 60 flank bases, CAGCAG (60..66), 60 flank bases. The read has one more
        // CAG, aligned as an insertion right before the repeat.
        let left: Vec<u8> = (0..60).map(|i| b"ACGTTGCA"[i % 8]).collect();
        let right: Vec<u8> = (0..60).map(|i| b"TTGACCAT"[i % 8]).collect();
        let reference = RefWindow { start: 0, seq: [&left[..], b"CAGCAG", &right[..]].concat() };
        let seq = [&left[..], b"CAGCAGCAG", &right[..]].concat();
        let r = read(0, vec![(Op::Match, 60), (Op::Ins, 3), (Op::Match, 66)], &seq);
        let opts = ExtractOptions { flank: 20, ..ExtractOptions::default() };
        assert_eq!(allele(&r, &locus(60, 66), &reference, &opts).as_deref(), Some(&b"CAGCAGCAG"[..]));
        // Not spanning with the flank.
        let short = read(50, vec![(Op::Match, 76)], &seq[50..126]);
        assert_eq!(allele(&short, &locus(60, 66), &reference, &opts), None);
    }

    #[test]
    fn long_insertion_just_before_the_repeat_is_kept() {
        // A 60-base expansion aligned as an insertion 3 bases into the left flank, as the
        // aligner places it when the flank end resembles the repeat. Flanks without
        // periodicity, so the flank fits in one place only.
        let mut x: u32 = 7;
        let mut random = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                    b"ACGT"[(x >> 16) as usize % 4]
                })
                .collect()
        };
        let (left, right) = (random(60), random(60));
        let reference = RefWindow { start: 0, seq: [&left[..], b"CAGCAG", &right[..]].concat() };
        let expansion = b"CAG".repeat(20);
        let seq = [&left[..57], &expansion[..], &left[57..], b"CAGCAG", &right[..]].concat();
        let r = read(0, vec![(Op::Match, 57), (Op::Ins, 60), (Op::Match, 69)], &seq);
        let opts = ExtractOptions { flank: 20, ..ExtractOptions::default() };
        let got = allele(&r, &locus(60, 66), &reference, &opts).expect("the read is kept");
        // The read has 60 more bases than the reference here, so the allele is 60 longer than
        // the 6-base repeat.
        assert_eq!(got.len(), 66, "allele {}", String::from_utf8_lossy(&got));
    }

    #[test]
    fn misplaced_flank_indel_is_realigned() {
        let left: Vec<u8> = (0..60).map(|i| b"ACGTTGCA"[i % 8]).collect();
        let right: Vec<u8> = (0..60).map(|i| b"TTGACCAT"[i % 8]).collect();
        let reference = RefWindow { start: 0, seq: [&left[..], b"CAGCAG", &right[..]].concat() };
        // Read lacks one CAG, but the aligner put the deletion 3 bases into the left flank's
        // end and shifted: deletion at 57..60, then the read's flank tail aligned to the repeat.
        let seq = [&left[..], b"CAG", &right[..]].concat();
        let r = read(0, vec![(Op::Match, 57), (Op::Del, 3), (Op::Match, 66)], &seq);
        let opts = ExtractOptions { flank: 20, ..ExtractOptions::default() };
        assert_eq!(allele(&r, &locus(60, 66), &reference, &opts).as_deref(), Some(&b"CAG"[..]));
    }

    #[test]
    fn an_insertion_at_a_shared_boundary_counts_once() {
        // Two adjacent loci, as ExpansionHunter's catalog splits ATXN7 and FXN: a left repeat
        // [60, 72) and a right one [72, 81), and a read with an insertion right between them.
        let left_flank: Vec<u8> = (0..60).map(|i| b"ACGTTGCA"[i % 8]).collect();
        let right_flank: Vec<u8> = (0..60).map(|i| b"TTGACCAT"[i % 8]).collect();
        let split = |left: &[u8], right: &[u8], ins: &[u8]| {
            let reference = RefWindow { start: 0, seq: [&left_flank[..], left, right, &right_flank[..]].concat() };
            let seq = [&left_flank[..], left, ins, right, &right_flank[..]].concat();
            let (l, r) = (left.len() as i64, right.len() as i64);
            let at = 60 + l;
            let cigar = vec![(Op::Match, at as u32), (Op::Ins, ins.len() as u32), (Op::Match, (r + 60) as u32)];
            let read = read(0, cigar, &seq);
            let boundary = SharedBoundary { left: vec![l as usize / 4], right: vec![3] };
            let mut a = locus(60, at);
            a.shared_end = Some(boundary.clone());
            let mut b = locus(at, at + r);
            b.shared_start = Some(boundary);
            let opts = ExtractOptions { flank: 20, ..ExtractOptions::default() };
            let got = |l: &Locus| allele(&read, l, &reference, &opts).expect("spans");
            (got(&a), got(&b))
        };
        // ATXN7: GCA and GCC repeats; the read has GCAGCA more of the first and GCC more of the
        // second.
        let (a, b) = split(b"GCAGCAGCAGCA", b"GCCGCCGCC", b"GCAGCAGCC");
        assert_eq!((a.len(), b.len()), (18, 12), "{} | {}", String::from_utf8_lossy(&a), String::from_utf8_lossy(&b));
        // FXN: an A tail, then GAA; three more GAA at the boundary belong to the GAA repeat.
        let (a, b) = split(b"AAAAAAAAAAAA", b"GAAGAAGAA", b"GAAGAAGAA");
        assert_eq!((a.len(), b.len()), (12, 18), "{} | {}", String::from_utf8_lossy(&a), String::from_utf8_lossy(&b));
    }

    #[test]
    fn an_insertion_sliding_to_a_shared_boundary_counts_once() {
        // FXN as aligned: the GAA expansion placed as AAGAAGAAG two bases into the A tail, which
        // is the same read as GAAGAAGAA at the boundary.
        let left_flank: Vec<u8> = (0..60).map(|i| b"ACGTTGCA"[i % 8]).collect();
        let right_flank: Vec<u8> = (0..60).map(|i| b"TTGACCAT"[i % 8]).collect();
        let (tail, gaa) = (b"AAAAAAAAAAAA", b"GAAGAAGAA");
        let reference = RefWindow { start: 0, seq: [&left_flank[..], tail, gaa, &right_flank[..]].concat() };
        let seq = [&left_flank[..], &tail[..10], b"AAGAAGAAG", &tail[10..], gaa, &right_flank[..]].concat();
        let read = read(0, vec![(Op::Match, 70), (Op::Ins, 9), (Op::Match, 71)], &seq);
        let boundary = SharedBoundary { left: vec![1], right: vec![3] };
        let mut a = locus(60, 72);
        a.shared_end = Some(boundary.clone());
        let mut b = locus(72, 81);
        b.shared_start = Some(boundary);
        let opts = ExtractOptions { flank: 20, ..ExtractOptions::default() };
        assert_eq!(allele(&read, &a, &reference, &opts).as_deref(), Some(&tail[..]));
        assert_eq!(allele(&read, &b, &reference, &opts).as_deref(), Some(&b"GAAGAAGAAGAAGAAGAA"[..]));
    }

    fn random(seed: u32, n: usize) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn an_expansion_in_a_soft_clip_is_read_up_to_the_far_flank() {
        // Reference: 400 flank bases, CAG x 10 (400..430), 400 flank bases. The read carries
        // CAG x 300; the aligner aligned the left flank and the first 30 repeat bases and
        // soft-clipped the rest, the expansion and the right flank.
        let (left, right) = (random(7, 400), random(11, 400));
        let reference = RefWindow { start: 0, seq: [&left[..], &b"CAG".repeat(10), &right[..]].concat() };
        let expansion = b"CAG".repeat(300);
        let seq = [&left[..], &expansion[..], &right[..]].concat();
        let opts = ExtractOptions::default();
        let clipped = read(0, vec![(Op::Match, 430), (Op::SoftClip, 1270)], &seq);
        assert_eq!(allele(&clipped, &locus(400, 430), &reference, &opts), Some(expansion.clone()));
        // The same read aligned on its right flank, the start clipped.
        let mut from_right = read(400, vec![(Op::SoftClip, 1270), (Op::Match, 430)], &seq);
        from_right.start = 400;
        assert_eq!(allele(&from_right, &locus(400, 430), &reference, &opts), Some(expansion));
        // A clip that does not start as the repeat (a chimera's other part) is not searched.
        let chimera = [&left[..], &b"CAG".repeat(10), &random(13, 1000)[..]].concat();
        let other = read(0, vec![(Op::Match, 430), (Op::SoftClip, 1000)], &chimera);
        assert_eq!(allele(&other, &locus(400, 430), &reference, &opts), None);
    }
}

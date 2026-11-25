//! A copy-number HMM (an alternative to Spectre's threshold runs): Viterbi segmentation of
//! each contig's normalized depth into copy-number states, with heterozygous SNV allele
//! fractions as a second emission where an SNV VCF is given.
//!
//! States are copy numbers 0 to [`HmmOptions::max_cn`] on the ploidy-scaled depth (so a
//! window at the genome median reads `ploidy`). A state's depth emission is a normal around
//! its copy number, its spread the genome's (estimated robustly from the windows) scaled with
//! the square root of the copy number, mixed with a small uniform share for outliers. A
//! heterozygous SNV in a window adds the likelihood of its allele fraction under the state:
//! a mixture of normals at `m / c` for `m` in `1..c`, or a low constant for fewer than two
//! copies, which should have no heterozygous sites. Runs of masked windows longer than
//! [`HmmOptions::max_gap`] end a segment.

use pynum::sum;

use crate::candidate::{Candidate, CnvType};
use crate::normalize::{ContigCoverage, Normalized};
use crate::snv::{Af, SnvVcf};
use crate::thresholds::Thresholds;

/// Settings of the HMM.
#[derive(Clone, Debug, PartialEq)]
pub struct HmmOptions {
    /// Highest copy number modelled.
    pub max_cn: usize,
    /// Log probability of a change of state between neighbouring windows.
    pub switch_log_prob: f64,
    /// Share of windows whose depth is an outlier under any state.
    pub outlier_share: f64,
    /// Runs of more masked windows than this split a contig into independent pieces.
    pub max_gap: usize,
    /// Heterozygous SNVs need at least this depth to count.
    pub min_snv_depth: i64,
}

impl Default for HmmOptions {
    fn default() -> Self {
        HmmOptions { max_cn: 6, switch_log_prob: -20.0, outlier_share: 0.01, max_gap: 100, min_snv_depth: 8 }
    }
}

/// A segment of one copy-number state: indexes into the contig's windows (`first..=last`, both
/// unmasked) and the state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Segment {
    pub first: usize,
    pub last: usize,
    pub copies: usize,
}

/// The depth spread at the baseline copy number: 1.4826 × the median absolute deviation of
/// windows within one copy of it, over every analysed contig.
pub fn depth_sd(norm: &Normalized, ploidy: i64) -> f64 {
    let base = ploidy as f64;
    let near: Vec<f64> = norm
        .contigs
        .iter()
        .flat_map(|c| c.normalized_ploidy.iter().copied())
        .filter(|v| (v - base).abs() < 1.0)
        .collect();
    if near.is_empty() {
        return 0.3;
    }
    let median = sum::median(&near);
    let dev: Vec<f64> = near.iter().map(|v| (v - median).abs()).collect();
    (1.4826 * sum::median(&dev)).max(0.05)
}

fn log_normal(x: f64, mean: f64, sd: f64) -> f64 {
    let z = (x - mean) / sd;
    -0.5 * z * z - sd.ln() - 0.5 * (2.0 * std::f64::consts::PI).ln()
}

/// Heterozygous SNV allele fractions per window of `c`, with each site's depth.
fn allele_fractions(c: &ContigCoverage, bin_size: i64, snv: &SnvVcf, min_depth: i64) -> Vec<Vec<(f64, i64)>> {
    let mut out = vec![Vec::new(); c.positions.len()];
    for r in snv.records(&c.name) {
        let Af::Value(af) = r.af else { continue };
        if !r.is_het() || r.dp < min_depth || r.qual <= 5.0 {
            continue;
        }
        // Windows are `[p, p + bin_size)` in 0-based coordinates; POS is 1-based.
        let at = r.pos - 1;
        let i = c.positions.partition_point(|&p| p + bin_size <= at);
        if i < c.positions.len() && c.positions[i] <= at {
            out[i].push((af, r.dp));
        }
    }
    out
}

/// Viterbi segmentation of one contig.
pub fn segment(
    c: &ContigCoverage,
    bin_size: i64,
    sd: f64,
    ploidy: i64,
    snv: Option<&SnvVcf>,
    opts: &HmmOptions,
) -> Vec<Segment> {
    let k = opts.max_cn + 1;
    let afs = snv.map(|v| allele_fractions(c, bin_size, v, opts.min_snv_depth));
    let uniform = (opts.outlier_share / (opts.max_cn as f64 + 1.0)).ln();
    let emit = |i: usize, state: usize| -> f64 {
        let x = c.normalized_ploidy[i];
        let spread = sd * (state.max(1) as f64 / ploidy as f64).sqrt();
        let depth = log_normal(x, state as f64, spread) + (1.0 - opts.outlier_share).ln();
        let mut e = depth.max(uniform);
        if let Some(afs) = &afs {
            let mut baf = 0.0;
            for &(af, dp) in &afs[i] {
                let s = (0.25 / dp as f64).sqrt() + 0.02;
                let l = if state < 2 {
                    0.2f64.ln()
                } else {
                    let n = (state - 1) as f64;
                    let p: f64 = (1..state).map(|m| log_normal(af, m as f64 / state as f64, s).exp()).sum::<f64>() / n;
                    (0.95 * p + 0.05).ln()
                };
                baf += l;
            }
            e += baf.clamp(-10.0, 10.0);
        }
        e
    };

    let observed: Vec<usize> = (0..c.normalized_ploidy.len()).filter(|&i| c.normalized_ploidy[i].is_finite()).collect();
    let mut out = Vec::new();
    // Pieces of observed windows without a long masked gap.
    let mut start = 0;
    for j in 1..=observed.len() {
        if j == observed.len() || observed[j] - observed[j - 1] > opts.max_gap {
            out.extend(viterbi(&observed[start..j], k, opts.switch_log_prob, &emit));
            start = j;
        }
    }
    out
}

/// The most likely states of `windows` (indexes), as segments.
fn viterbi(windows: &[usize], k: usize, switch: f64, emit: &impl Fn(usize, usize) -> f64) -> Vec<Segment> {
    if windows.is_empty() {
        return Vec::new();
    }
    let mut score: Vec<f64> = (0..k).map(|s| emit(windows[0], s)).collect();
    let mut back: Vec<Vec<u8>> = Vec::with_capacity(windows.len());
    back.push(vec![0; k]);
    for &w in &windows[1..] {
        let (best, best_score) =
            score.iter().enumerate().fold((0, f64::NEG_INFINITY), |b, (s, &v)| if v > b.1 { (s, v) } else { b });
        let mut from = vec![0u8; k];
        let next: Vec<f64> = (0..k)
            .map(|s| {
                let (prev, v) =
                    if score[s] >= best_score + switch { (s, score[s]) } else { (best, best_score + switch) };
                from[s] = prev as u8;
                v + emit(w, s)
            })
            .collect();
        back.push(from);
        score = next;
    }
    let mut state =
        score.iter().enumerate().fold((0, f64::NEG_INFINITY), |b, (s, &v)| if v > b.1 { (s, v) } else { b }).0;
    let mut states = vec![0usize; windows.len()];
    for t in (0..windows.len()).rev() {
        states[t] = state;
        state = back[t][state] as usize;
    }
    let mut segs: Vec<Segment> = Vec::new();
    for (t, &s) in states.iter().enumerate() {
        match segs.last_mut() {
            Some(seg) if seg.copies == s => seg.last = windows[t],
            _ => segs.push(Segment { first: windows[t], last: windows[t], copies: s }),
        }
    }
    segs
}

/// Candidates from the segments of every analysed contig that are not at `ploidy`, span at
/// least `min_len` bases and whose median depth is beyond `bounds` (below the deletion
/// threshold or above the duplication one), as `(contig, candidates)` in contig order. IDs
/// are left empty.
pub fn call(
    norm: &Normalized,
    ploidy: i64,
    snv: Option<&SnvVcf>,
    min_len: i64,
    bounds: &Thresholds,
    opts: &HmmOptions,
) -> Vec<(String, Vec<Candidate>)> {
    let sd = depth_sd(norm, ploidy);
    norm.contigs
        .iter()
        .map(|c| {
            let cands = segment(c, norm.bin_size, sd, ploidy, snv, opts)
                .into_iter()
                .filter(|s| s.copies as i64 != ploidy && c.positions[s.last] - c.positions[s.first] >= min_len)
                .map(|s| {
                    let t = if (s.copies as i64) < ploidy { CnvType::Del } else { CnvType::Dup };
                    let pos = c.positions[s.first..=s.last].to_vec();
                    let cov = c.normalized_ploidy[s.first..=s.last].to_vec();
                    Candidate::new(&c.name, t, norm.bin_size, pos, cov)
                })
                .filter(|cand| match cand.cnv_type {
                    CnvType::Del => cand.median_cov_norm < bounds.del.value(),
                    CnvType::Dup => cand.median_cov_norm > bounds.dup.value(),
                })
                .collect();
            (c.name.clone(), cands)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::Stats;

    fn contig(depth: Vec<f64>) -> ContigCoverage {
        let stats = Stats { average: 2.0, std_dev: 0.2, median: 2.0, min: 0.0, max: 4.0 };
        ContigCoverage {
            name: "chr1".into(),
            positions: (0..depth.len() as i64).map(|i| i * 1000).collect(),
            raw: depth.clone(),
            normalized: depth.iter().map(|v| v / 2.0).collect(),
            normalized_ploidy: depth,
            stats,
            norm_stats: stats,
        }
    }

    #[test]
    fn finds_a_deletion_and_a_duplication() {
        // Deterministic jitter around the copy number.
        let jitter = |i: usize| ((i * 7919) % 13) as f64 / 13.0 * 0.4 - 0.2;
        let depth: Vec<f64> = (0..1000)
            .map(|i| {
                let cn = if (200..400).contains(&i) {
                    1.0
                } else if (600..750).contains(&i) {
                    3.0
                } else {
                    2.0
                };
                cn + jitter(i)
            })
            .collect();
        let segs = segment(&contig(depth), 1000, 0.15, 2, None, &HmmOptions::default());
        let calls: Vec<(usize, usize, usize)> =
            segs.iter().filter(|s| s.copies != 2).map(|s| (s.first, s.last, s.copies)).collect();
        assert_eq!(calls, vec![(200, 399, 1), (600, 749, 3)]);
    }
}

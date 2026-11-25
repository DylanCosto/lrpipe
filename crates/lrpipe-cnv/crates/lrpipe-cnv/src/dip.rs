//! Hartigan's dip test of unimodality and UniDip, written from their papers: J. A. Hartigan and
//! P. M. Hartigan, "The Dip Test of Unimodality", Annals of Statistics 13 (1985), with the
//! computation of P. M. Hartigan, Algorithm AS 217, Applied Statistics 34 (1985); and S. Maurus
//! and C. Plant, "Skinny-dip: Clustering in a Sea of Noise", KDD 2016.
//!
//! Spectre uses a GPL-3.0 implementation of these (`spectre/libs/unidip`); this one is
//! independent and is checked against that one's results as a black box
//! (`tools/dip_vectors.py`, `tests/dip.rs`).
//!
//! Data are reduced to their sorted distinct values with counts, as Spectre's version does,
//! so the modal interval is reported in distinct-value indexes.

use pynum::numpy_random::RandomState;

/// The dip of a sample and its modal interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dip {
    /// The dip: the largest distance, in CDF units, from the empirical CDF to the closest
    /// unimodal CDF.
    pub dip: f64,
    /// The modal interval, as indexes into the distinct values (`(-1, m - 1)` when there are
    /// at most 4 distinct values).
    pub modal: (isize, isize),
    /// The number of distinct values.
    pub distinct: usize,
}

/// Distinct values in ascending order and their counts.
fn distinct(data: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut v = data.to_vec();
    v.sort_by(f64::total_cmp);
    let (mut xs, mut counts): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    for x in v {
        if xs.last() == Some(&x) {
            *counts.last_mut().expect("pushed with the value") += 1.0;
        } else {
            xs.push(x);
            counts.push(1.0);
        }
    }
    (xs, counts)
}

/// The dip of `data`.
pub fn dip(data: &[f64]) -> Dip {
    let (x, c) = distinct(data);
    let m = x.len();
    if m <= 4 {
        return Dip { dip: 0.0, modal: (-1, m as isize - 1), distinct: m };
    }
    let n: f64 = c.iter().sum();
    // w[k + 1] is the count up to and including x[k]: the CDF jumps from w[k] to w[k + 1] at
    // x[k]. The convex minorant is fitted to the lower points (x[k], w[k]), the concave
    // majorant to the upper points (x[k], w[k + 1]).
    // In CDF units: cumulative counts divided by the total.
    let mut w = vec![0.0; m + 1];
    let mut cum = 0.0;
    for k in 0..m {
        cum += c[k];
        w[k + 1] = cum;
    }
    for v in &mut w {
        *v /= n;
    }
    let lower = |k: usize| w[k];
    let upper = |k: usize| w[k + 1];
    let (mut low, mut high) = (0usize, m - 1);
    let mut d_total = 0.0f64;
    loop {
        let gcm = hull(&x, low, high, &lower, Hull::Convex);
        let lcm = hull(&x, low, high, &upper, Hull::Concave);
        // The largest distance between the two fits, at a vertex of either.
        let mut d = 0.0;
        let (mut new_low, mut new_high) = (low, high);
        let at = |hull: &[usize], f: &dyn Fn(usize) -> f64, k: usize| -> (f64, usize, usize) {
            let s = hull.partition_point(|&v| v <= k).max(1) - 1;
            let (a, b) = (hull[s], hull[(s + 1).min(hull.len() - 1)]);
            let y = if a == b || x[b] == x[a] { f(a) } else { f(a) + (f(b) - f(a)) * (x[k] - x[a]) / (x[b] - x[a]) };
            (y, a, b)
        };
        let mut vertices: Vec<(usize, bool)> =
            gcm.iter().map(|&k| (k, true)).chain(lcm.iter().map(|&k| (k, false))).collect();
        vertices.sort_by_key(|&(k, is_gcm)| (k, !is_gcm));
        for (k, is_gcm) in vertices {
            if is_gcm {
                let (y_lcm, _, right) = at(&lcm, &upper, k);
                let dx = y_lcm - lower(k);
                if dx >= d {
                    d = dx;
                    new_low = k;
                    new_high = right.max(k);
                }
            } else {
                let (y_gcm, left, _) = at(&gcm, &lower, k);
                let dx = upper(k) - y_gcm;
                if dx >= d {
                    d = dx;
                    new_low = left.min(k);
                    new_high = k;
                }
            }
        }
        if d < d_total || (new_low == low && new_high == high && d_total > 0.0) {
            break;
        }
        // The deviations outside the new modal interval: the empirical CDF above the convex
        // minorant left of it, below the concave majorant right of it.
        let mut outer = 0.0f64;
        for pair in gcm.windows(2).filter(|p| p[1] <= new_low) {
            let (a, b) = (pair[0], pair[1]);
            for k in a..=b {
                let line = if x[b] == x[a] {
                    lower(a)
                } else {
                    lower(a) + (lower(b) - lower(a)) * (x[k] - x[a]) / (x[b] - x[a])
                };
                outer = outer.max(upper(k) - line);
            }
        }
        for pair in lcm.windows(2).filter(|p| p[0] >= new_high) {
            let (a, b) = (pair[0], pair[1]);
            for k in a..=b {
                let line = if x[b] == x[a] {
                    upper(a)
                } else {
                    upper(a) + (upper(b) - upper(a)) * (x[k] - x[a]) / (x[b] - x[a])
                };
                outer = outer.max(line - lower(k));
            }
        }
        d_total = d_total.max(outer).max(if new_low == low && new_high == high { d } else { 0.0 });
        if new_low == low && new_high == high {
            break;
        }
        low = new_low;
        high = new_high;
    }
    Dip { dip: d_total / 2.0, modal: (low as isize, high as isize), distinct: m }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Hull {
    Convex,
    Concave,
}

/// Vertex indexes, ascending, of the convex (from below) or concave (from above) hull of the
/// points `(x[k], f(k))` for `k` in `low..=high`.
fn hull(x: &[f64], low: usize, high: usize, f: &dyn Fn(usize) -> f64, kind: Hull) -> Vec<usize> {
    let mut h: Vec<usize> = Vec::new();
    for k in low..=high {
        while h.len() >= 2 {
            let (a, b) = (h[h.len() - 2], h[h.len() - 1]);
            // Cross product of (b - a) and (k - a): drop b when it is not strictly outside.
            let cross = (x[b] - x[a]) * (f(k) - f(a)) - (f(b) - f(a)) * (x[k] - x[a]);
            let drop = match kind {
                Hull::Convex => cross <= 0.0,
                Hull::Concave => cross >= 0.0,
            };
            if drop {
                h.pop()
            } else {
                break;
            };
        }
        h.push(k);
    }
    h
}

/// The dip test: the dip of `data` and its p-value against `trials` uniform samples of as
/// many values as `data` has distinct values, drawn from `rng` (numpy's global generator in
/// Spectre). `None` when every uniform sample has a dip of 0.
pub fn dip_test(data: &[f64], trials: usize, rng: &mut RandomState) -> (Dip, Option<f64>) {
    let d = dip(data);
    let m = d.distinct;
    let mut uniform = vec![0.0; m];
    let mut larger = 0usize;
    let mut total = 0.0;
    for _ in 0..trials {
        for u in &mut uniform {
            *u = rng.uniform();
        }
        let u = dip(&uniform).dip;
        total += u;
        if d.dip < u {
            larger += 1;
        }
    }
    let p = (total != 0.0).then(|| (larger + 1) as f64 / (trials + 1) as f64);
    (d, p)
}

/// UniDip (Maurus and Plant): the modal intervals of `data`, as inclusive index ranges into
/// the sorted data, each a unimodal mode at significance 0.05 with 100 trials per test.
pub fn unidip(data: &[f64], rng: &mut RandomState) -> Vec<(usize, usize)> {
    let mut sorted = data.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut intervals = modes(&sorted, 0, sorted.len(), Model::Yes, rng);
    intervals.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for iv in intervals {
        match merged.last_mut() {
            Some(last) if iv.0 as isize - last.1 as isize <= 1 => last.1 = iv.1,
            _ => merged.push(iv),
        }
    }
    merged
}

const ALPHA: f64 = 0.05;
const TRIALS: usize = 100;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Model {
    Yes,
    No,
}

/// The modes in `sorted[start..end]`.
fn modes(sorted: &[f64], start: usize, end: usize, model: Model, rng: &mut RandomState) -> Vec<(usize, usize)> {
    let dat = &sorted[start..end];
    let (d, p) = dip_test(dat, TRIALS, rng);
    let Some(p) = p else { return Vec::new() };
    if p > ALPHA {
        return match model {
            Model::Yes => vec![(start, end - 1)],
            Model::No => {
                let (a, b) = full_interval(dat, rng);
                vec![(start + a, start + b)]
            }
        };
    }
    // Several modes: recurse into the modal interval, then the parts left and right of it.
    let (m0, m1) = (d.modal.0.max(0) as usize, d.modal.1.max(0) as usize);
    // The modal interval is in distinct-value indexes but taken as data offsets, as in
    // Spectre's version; stop when that does not shrink the range.
    if m1 <= m0 || (m0 == 0 && start + m1 >= end) {
        return vec![(start, end - 1)];
    }
    let inner = modes(sorted, start + m0, (start + m1).min(end), Model::Yes, rng);
    let mut found = inner.clone();
    let mut parts: Vec<(usize, usize)> = inner.iter().map(|&(a, b)| (a - start, b - start + 1)).collect();
    parts.push((m0, m1));
    let l = *parts.iter().min_by_key(|t| t.1).expect("modal interval");
    let h = *parts.iter().max().expect("modal interval");
    if let (_, Some(p)) = dip_test(&dat[..l.1.min(dat.len())], TRIALS, rng)
        && p < ALPHA
    {
        found.extend(modes(sorted, start, start + l.0, Model::No, rng));
    }
    if let (_, Some(p)) = dip_test(&dat[h.0.min(dat.len())..], TRIALS, rng)
        && p < ALPHA
    {
        found.extend(modes(sorted, start + h.1, end, Model::No, rng));
    }
    found
}

/// The full extent of the mode in a non-model part, found by testing the part mirrored about
/// its end (Skinny-dip's mirroring): the modal interval of the mirrored data, mapped back.
fn full_interval(dat: &[f64], rng: &mut RandomState) -> (usize, usize) {
    if dat.is_empty() {
        return (0, 0);
    }
    let last = dat[dat.len() - 1];
    let mut mirrored: Vec<f64> = dat.iter().map(|&v| 2.0 * last - v).chain(dat.iter().copied()).collect();
    mirrored.sort_by(f64::total_cmp);
    let (d, _) = dip_test(&mirrored, TRIALS, rng);
    let n = dat.len();
    let lo = (d.modal.0.max(0) as usize).saturating_sub(n).min(n - 1);
    (lo, n - 1)
}

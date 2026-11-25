//! Refining raw candidates into calls (spec §3.6-3.9; `CNVAnalysis.refine_cnv_calls`,
//! analysis.py:333-715, and `CNV_Post_Processing.cnv_fit`, cnv_post_processing.py).

use pynum::sum;

use crate::Error;
use crate::candidate::{Arena, CandidateId, Side};
use crate::metadata::Metadata;
use crate::normalize::ContigCoverage;
use crate::thresholds::Thresholds;

/// Settings for refinement.
#[derive(Clone, Debug)]
pub struct RefineOptions {
    /// `--cov-diff-threshold` (0.8): the largest median difference of candidates merged.
    pub cov_diff_threshold: f64,
    /// Merges may always bridge this many bases (`dist_min_overwrite`).
    pub dist_min_overwrite: f64,
    /// `--min-cnv-len` (100 kb): calls must span at least this much.
    pub min_cnv_len: i64,
    /// `--qc-noise-allow` (0.5): calls need more than this share of windows outside the
    /// strict thresholds.
    pub qc_noise_allow: f64,
}

impl Default for RefineOptions {
    fn default() -> Self {
        RefineOptions {
            cov_diff_threshold: 0.8,
            dist_min_overwrite: 10_000.0,
            min_cnv_len: 100_000,
            qc_noise_allow: 0.5,
        }
    }
}

/// What refinement needs about one contig.
#[derive(Debug)]
pub struct ContigContext<'a> {
    pub coverage: &'a ContigCoverage,
    pub metadata: &'a Metadata,
    pub bin_size: i64,
    /// The quantile thresholds (used to decide merges).
    pub thresholds: &'a Thresholds,
    /// The strict thresholds (mean ± 1.85 SD of the window sample), for the final filter.
    pub strict: &'a Thresholds,
}

/// `refine_cnv_calls` for one contig: merge, scaffold, grow and trim, resolve overlaps,
/// filter. Candidates are changed in place in `arena`.
///
/// # Errors
///
/// If resolving an overlap needs a window a candidate does not have (a ValueError in Spectre).
pub fn refine(
    candidates: &[CandidateId],
    arena: &mut Arena,
    ctx: &ContigContext<'_>,
    opts: &RefineOptions,
) -> Result<Vec<CandidateId>, Error> {
    let merged = merge(candidates, arena, ctx, opts);
    if merged.len() >= 2 {
        for &id in &merged {
            scaffold(arena, id, ctx);
        }
    }
    for &id in &merged {
        fit(arena, id, ctx);
    }
    let cleaned = clean_by_overlap(&merged, arena)?;
    let kept = cleaned
        .into_iter()
        .filter(|&id| passes_strict_thresholds(arena, id, ctx.strict, opts.qc_noise_allow))
        .filter(|&id| arena.get(id).size as f64 >= opts.min_cnv_len as f64)
        .collect();
    Ok(kept)
}

/// `merge_candidates`: merge passes until none merges, then passes that may also bridge
/// excluded regions, until none merges.
fn merge(
    candidates: &[CandidateId],
    arena: &mut Arena,
    ctx: &ContigContext<'_>,
    opts: &RefineOptions,
) -> Vec<CandidateId> {
    if candidates.len() <= 1 {
        return candidates.to_vec();
    }
    let (mut n, mut merged) = merge_pass(candidates, arena, ctx, opts, Blacklist::Blocks);
    while n > 0 {
        (n, merged) = merge_pass(&merged, arena, ctx, opts, Blacklist::Blocks);
    }
    loop {
        let (n, m) = merge_pass(&merged, arena, ctx, opts, Blacklist::Ignored);
        merged = m;
        if n == 0 {
            break;
        }
    }
    if merged.is_empty() { candidates.to_vec() } else { merged }
}

/// Whether an excluded region between two candidates prevents merging them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Blacklist {
    Blocks,
    Ignored,
}

/// `cnv_candidate_merge`: one left-to-right pass merging neighbouring pairs.
fn merge_pass(
    c: &[CandidateId],
    arena: &mut Arena,
    ctx: &ContigContext<'_>,
    opts: &RefineOptions,
    blacklist: Blacklist,
) -> (usize, Vec<CandidateId>) {
    let (del, dup) = (ctx.thresholds.del.value(), ctx.thresholds.dup.value());
    let mut idx = 0;
    let mut n_merges = 0;
    let mut out = Vec::new();
    let mut cand = c[0];
    let mut last_merged = false;
    while idx < c.len() - 1 {
        let (this, next) = (arena.get(c[idx]), arena.get(c[idx + 1]));
        let (end0, start1) = (this.end, next.start);
        let same_type = this.cnv_type == next.cnv_type;
        let same_status = (this.median_cov_norm - next.median_cov_norm).abs() <= opts.cov_diff_threshold;
        let acc = arena.get(cand);
        let mean_size = sum::mean(&[acc.size as f64, next.size as f64]);
        // `40 * np.e ** (-m / 1e6) + 0.1`: numpy's power, i.e. libm pow.
        let percent = 40.0 * std::f64::consts::E.powf(-mean_size / 1_000_000.0) + 0.1;
        let threshold = mean_size * percent / 100.0;
        let distance = if opts.dist_min_overwrite > threshold { opts.dist_min_overwrite } else { threshold };
        let distance_allowed = acc.end as f64 + (distance + ctx.bin_size as f64) >= next.start as f64;
        let spans_excluded = blacklist == Blacklist::Blocks
            && ctx.metadata.contigs.iter().find(|(name, _)| *name == acc.chromosome).is_some_and(|(_, regions)| {
                regions.iter().any(|r| [r.start, r.end].iter().any(|&v| end0 <= v && v <= start1))
            });
        let outside = !(del <= acc.median_cov_norm && acc.median_cov_norm <= dup);
        if distance_allowed && same_type && same_status && !spans_excluded && outside {
            n_merges += 1;
            let (a, b) = arena.pair_mut(cand, c[idx + 1]);
            let (pos, cov, id, refs) = (b.pos.clone(), b.cov.clone(), b.id.clone(), b.merged_sample_references.clone());
            a.absorb(&pos, &cov, &id, &refs, Side::Back);
            out.push(cand);
            idx += 1;
            if idx >= c.len() - 1 {
                last_merged = true;
            }
        } else {
            arena.get_mut(cand).update_median();
            out.push(cand);
        }
        idx += 1;
        if idx < c.len() {
            cand = c[idx];
        }
    }
    if !last_merged {
        out.push(cand);
    }
    (n_merges, out)
}

/// `scaffold_candidates` for one candidate: gaps between windows are filled with the
/// median of the contig's windows there (the candidate's own window before the gap
/// included and replaced, quirk Q-SCAF), unless one of them is NaN.
fn scaffold(arena: &mut Arena, id: CandidateId, ctx: &ContigContext<'_>) {
    let positions = &ctx.coverage.positions;
    let ploidy = &ctx.coverage.normalized_ploidy;
    let c = arena.get_mut(id);
    let gaps: Vec<usize> =
        c.pos.windows(2).enumerate().filter(|(_, w)| w[1] - w[0] != ctx.bin_size).map(|(k, _)| k).collect();
    for &k in gaps.iter().rev() {
        let index_of = |p: i64| positions.iter().position(|&x| x == p).expect("candidate windows are contig windows");
        let (gs, ge) = (index_of(c.pos[k]), index_of(c.pos[k + 1]));
        let gap = &ploidy[gs..ge];
        if gap.iter().any(|v| v.is_nan()) {
            continue;
        }
        let m = sum::nanmedian(gap);
        c.cov.splice(k..=k, std::iter::repeat_n(m, gap.len()));
        c.pos.splice(k..=k, positions[gs..ge].iter().copied());
    }
    c.update_median();
}

/// `get_coverage_slope`: `v[i] - mean(v[i..=i+w])` at every `w`-th window, held for the `w`
/// windows from there.
fn slope(v: &[f64], w: usize) -> Vec<f64> {
    let mut d = vec![0.0; v.len()];
    for i in (0..v.len()).step_by(w) {
        let s = v[i] - sum::mean(&v[i..(i + w + 1).min(v.len())]);
        for x in &mut d[i..(i + w).min(v.len())] {
            *x = s;
        }
    }
    d
}

/// `cnv_fit` for one candidate: grow each end to the nearest window outside the
/// candidate's depth range with a moderate slope, trim back to the steepest change, and
/// take the windows in between (spec §3.8). The median, CN and GT are not updated.
fn fit(arena: &mut Arena, id: CandidateId, ctx: &ContigContext<'_>) {
    let positions = &ctx.coverage.positions;
    let ploidy = &ctx.coverage.normalized_ploidy;
    let c = arena.get(id);
    let (start, end) = (c.start, c.end);
    let n = (100 * ctx.bin_size) as f64;
    let n = if c.size as f64 * 0.25 > n { c.size as f64 * 0.25 } else { n };

    // The extended region: windows from start - n to end + n, in order.
    let ext: Vec<usize> = (0..positions.len())
        .filter(|&i| {
            let p = positions[i] as f64;
            p >= start as f64 - n && p <= end as f64 + n
        })
        .collect();
    let pos: Vec<i64> = ext.iter().map(|&i| positions[i]).collect();
    let v: Vec<f64> = ext.iter().map(|&i| ploidy[i]).collect();
    let slope5 = slope(&v, 5);
    let slope10 = slope(&v, 10);
    let cnv: Vec<f64> = pos.iter().zip(&v).filter(|(p, _)| **p >= start && **p <= end).map(|(_, &x)| x).collect();
    let cnv_sd = sum::std(&cnv) * 3.0;
    let cnv_mean = sum::mean(&cnv);
    let in_sd = |x: f64| cnv_mean - cnv_sd < x && x < cnv_mean + cnv_sd;
    let (m10, s10) = (sum::series_mean(&slope10), sum::series_std(&slope10));
    let (upper, lower) = (m10 + s10, m10 - s10);
    let moderate = |i: usize| !in_sd(v[i]) && lower < slope10[i] && slope10[i] < upper;
    let new_end = (0..pos.len()).find(|&i| moderate(i) && pos[i] > end).map_or(end, |i| pos[i]);
    let new_start = (0..pos.len()).rev().find(|&i| moderate(i) && pos[i] < start).map_or(start, |i| pos[i]);

    // Trim, over the rows between the new borders.
    let rows: Vec<usize> = (0..pos.len()).filter(|&i| pos[i] >= new_start && pos[i] <= new_end).collect();
    let slope_diff: Vec<f64> = rows
        .iter()
        .enumerate()
        .map(|(k, &i)| if k == 0 { f64::NAN } else { (slope5[i] - slope5[rows[k - 1]]).abs() })
        .collect();
    let steepest = |keep: &dyn Fn(i64) -> bool| -> (Vec<usize>, f64) {
        let ks: Vec<usize> = (0..rows.len()).filter(|&k| keep(pos[rows[k]])).collect();
        let max = ks.iter().map(|&k| slope_diff[k]).filter(|x| !x.is_nan()).fold(f64::NAN, f64::max);
        (ks, max)
    };
    let (right, rmax) = steepest(&|p| end < p && p < new_end);
    let new_end_trim = right.iter().find(|&&k| slope_diff[k] == rmax).map_or(new_end, |&k| pos[rows[k]]);
    let (left, lmax) = steepest(&|p| new_start < p && p < start);
    let new_start_trim = left.iter().rev().find(|&&k| slope_diff[k] == lmax).map_or(new_start, |&k| pos[rows[k]]);

    let keep: Vec<usize> =
        rows.iter().copied().filter(|&i| pos[i] >= new_start_trim && pos[i] <= new_end_trim).collect();
    let c = arena.get_mut(id);
    c.pos = keep.iter().map(|&i| pos[i]).collect();
    c.cov = keep.iter().map(|&i| v[i]).collect();
    c.start = c.pos[0];
    c.end = *c.pos.last().expect("a refined candidate keeps its own windows");
    c.size = (c.end - c.start).abs();
}

/// `clean_by_overlap`, transliterated with its index bugs (spec §3.9): overlapping
/// candidates are joined into the larger one and the other is dropped.
fn clean_by_overlap(l: &[CandidateId], arena: &mut Arena) -> Result<Vec<CandidateId>, Error> {
    let n = l.len();
    let mut remove: Vec<usize> = Vec::new();
    let (mut idx, mut offset) = (0usize, 0usize);
    while idx + offset < n {
        let cur = l[idx + offset];
        let mut early_stop = false;
        while idx + offset < n - 1 && !early_stop {
            let nxt = l[idx + offset + 1];
            let (c, x) = (arena.get(cur), arena.get(nxt));
            if !((x.start < c.start && c.start < x.end) || (x.start < c.end && c.end < x.end)) {
                break;
            }
            let cur_full = c.start <= x.start && x.end <= c.end;
            let nxt_full = x.start <= c.start && c.end <= x.end;
            if cur_full {
                remove.push(idx + offset + 1);
                offset += 1;
            }
            if nxt_full {
                remove.push(idx);
                early_stop = true;
            } else {
                let (master, slave, remove_candidate) =
                    if c.size < x.size { (nxt, cur, idx) } else { (cur, nxt, idx + offset + 1) };
                let (m, s) = arena.pair_mut(master, slave);
                let find = |p: i64| {
                    s.pos
                        .iter()
                        .position(|&q| q == p)
                        .ok_or_else(|| Error::Input(format!("overlap resolution: {} has no window at {p}", s.id)))
                };
                let (pos, cov, side) = if m.end >= s.end {
                    let k = find(m.start)?;
                    (s.pos[..k].to_vec(), s.cov[..k].to_vec(), Side::Front)
                } else {
                    let k = find(m.end)? + 1;
                    (s.pos[k..].to_vec(), s.cov[k..].to_vec(), Side::Back)
                };
                let (id, refs) = (s.id.clone(), s.merged_sample_references.clone());
                m.absorb(&pos, &cov, &id, &refs, side);
                remove.push(remove_candidate);
                offset += 1;
            }
        }
        idx += 1;
    }
    Ok((0..n).filter(|i| !remove.contains(i)).map(|i| l[i]).collect())
}

/// `clean_by_threshold`: the mean depth (NaN-propagating) must lie outside the strict
/// thresholds and enough windows must too.
fn passes_strict_thresholds(arena: &Arena, id: CandidateId, strict: &Thresholds, qc_noise_allow: f64) -> bool {
    let c = arena.get(id);
    let (del, dup) = (strict.del.value(), strict.dup.value());
    let mu = sum::mean(&c.cov);
    let outside = c.cov.iter().filter(|&&v| v < del).count() + c.cov.iter().filter(|&&v| v > dup).count();
    !(del <= mu && mu <= dup) && outside as f64 > c.cov.len() as f64 * qc_noise_allow
}

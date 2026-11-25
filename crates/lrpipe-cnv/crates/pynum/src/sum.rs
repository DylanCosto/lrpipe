//! numpy's float reductions (numpy 1.24) and pandas' (2.0, without bottleneck).
//!
//! numpy adds float64 arrays pairwise in blocks of 8 (`pairwise_sum` in
//! `numpy/core/src/umath/loops_utils.h.src`), so a sum, and every mean and variance built on
//! one, depends on the order of additions. These functions add in the same order.

/// numpy's `pairwise_sum` over `a`.
fn pairwise(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &v in a {
            res += v;
        }
        res
    } else if n <= 128 {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[..8]);
        let blocks = n - n % 8;
        for chunk in a[8..blocks].as_chunks::<8>().0 {
            for (acc, &v) in r.iter_mut().zip(chunk) {
                *acc += v;
            }
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for &v in &a[blocks..] {
            res += v;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise(&a[..n2]) + pairwise(&a[n2..])
    }
}

/// numpy's reduction buffer: a reduction hands its inner loop at most this many elements at
/// a time and carries the total between buffers.
const BUFSIZE: usize = 8192;

/// `np.sum(a)` of a contiguous float64 array: pairwise within each buffer, sequentially
/// across buffers, starting from 0.
pub fn sum(a: &[f64]) -> f64 {
    let mut total = 0.0;
    for chunk in a.chunks(BUFSIZE) {
        total += pairwise(chunk);
    }
    total
}

/// `np.mean(a)`: NaN for an empty array.
pub fn mean(a: &[f64]) -> f64 {
    sum(a) / a.len() as f64
}

/// `np.var(a)` (ddof 0): the mean of squared deviations, each sum pairwise.
pub fn var(a: &[f64]) -> f64 {
    let m = mean(a);
    let sq: Vec<f64> = a.iter().map(|&v| (v - m) * (v - m)).collect();
    sum(&sq) / a.len() as f64
}

/// `np.std(a)` (ddof 0).
pub fn std(a: &[f64]) -> f64 {
    var(a).sqrt()
}

/// A copy of `a` with NaNs replaced by 0, and the count of non-NaN values: numpy's
/// `_replace_nan`. The NaN-aware reductions add over the whole array, zeros included, so
/// their result depends on where the NaNs were.
fn replace_nan(a: &[f64]) -> (Vec<f64>, usize) {
    let mut count = 0;
    let filled = a
        .iter()
        .map(|&v| {
            if v.is_nan() {
                0.0
            } else {
                count += 1;
                v
            }
        })
        .collect();
    (filled, count)
}

/// `np.nanmean(a)`.
pub fn nanmean(a: &[f64]) -> f64 {
    let (filled, count) = replace_nan(a);
    sum(&filled) / count as f64
}

/// `np.nanvar(a, ddof)`: deviations at NaN positions stay 0.
pub fn nanvar(a: &[f64], ddof: usize) -> f64 {
    let (filled, count) = replace_nan(a);
    let avg = sum(&filled) / count as f64;
    let sq: Vec<f64> =
        a.iter().zip(&filled).map(|(&raw, &v)| if raw.is_nan() { 0.0 } else { (v - avg) * (v - avg) }).collect();
    sum(&sq) / (count as f64 - ddof as f64)
}

/// `np.nanstd(a)` (ddof 0).
pub fn nanstd(a: &[f64]) -> f64 {
    nanvar(a, 0).sqrt()
}

/// `np.median(a)` of values without NaN: `np.mean` of the middle value or the two middle
/// values (so a median of -0.0 is 0.0).
///
/// # Panics
///
/// If `a` is empty.
pub fn median(a: &[f64]) -> f64 {
    assert!(!a.is_empty(), "median of an empty array");
    let mut v = a.to_vec();
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 { mean(&v[n / 2..=n / 2]) } else { mean(&v[n / 2 - 1..=n / 2]) }
}

/// `np.nanmedian(a)`: NaN when every value is NaN.
pub fn nanmedian(a: &[f64]) -> f64 {
    let kept: Vec<f64> = a.iter().copied().filter(|v| !v.is_nan()).collect();
    if kept.is_empty() { f64::NAN } else { median(&kept) }
}

/// `np.percentile(a, q)` with the default linear method, for `a` without NaN.
///
/// # Panics
///
/// If `a` is empty.
pub fn percentile(a: &[f64], q: f64) -> f64 {
    assert!(!a.is_empty(), "percentile of an empty array");
    let mut v = a.to_vec();
    v.sort_by(f64::total_cmp);
    let n = v.len();
    let q = q / 100.0;
    let virtual_index = (n - 1) as f64 * q;
    // numpy's `_get_indexes`: past the last value both indexes become -1 (the last value),
    // and the weight is computed against that -1, not against the floor.
    let (lo, hi, previous) = if virtual_index >= (n - 1) as f64 {
        (n - 1, n - 1, -1.0)
    } else {
        let below = virtual_index.floor();
        (below as usize, below as usize + 1, below)
    };
    lerp(v[lo], v[hi], virtual_index - previous)
}

/// numpy's `_lerp`: from `a` for `t < 0.5`, back from `b` otherwise.
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    let diff = b - a;
    if t >= 0.5 { b - diff * (1.0 - t) } else { a + diff * t }
}

/// `pd.Series.mean()` (skipna) of float64 values.
pub fn series_mean(a: &[f64]) -> f64 {
    nanmean(a)
}

/// pandas' `nanvar` (skipna) of float64 values: `((avg - x) ** 2)` with NaNs as 0, summed,
/// over `count - ddof`; NaN when `count <= ddof`. `Series.var(ddof)`, and `np.var(series)`
/// with ddof 0.
pub fn series_var(a: &[f64], ddof: usize) -> f64 {
    let (filled, count) = replace_nan(a);
    if count <= ddof {
        return f64::NAN;
    }
    let avg = sum(&filled) / count as f64;
    let sq: Vec<f64> = a.iter().map(|&v| if v.is_nan() { 0.0 } else { (avg - v) * (avg - v) }).collect();
    sum(&sq) / (count - ddof) as f64
}

/// `pd.Series.std()` (skipna, ddof 1) of float64 values: NaN with fewer than two values.
pub fn series_std(a: &[f64]) -> f64 {
    series_var(a, 1).sqrt()
}

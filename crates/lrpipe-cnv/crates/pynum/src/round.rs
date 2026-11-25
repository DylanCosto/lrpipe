//! numpy's rounding to decimal places (`np.round`, and `round()` of an `np.float64`).

/// `np.round(x, decimals)` for `decimals >= 0`: scale by 10^decimals, round half to even,
/// scale back. Unlike Python's `round(float, d)` this is not correctly rounded: 0.075 rounds
/// to 0.08 because 0.075 * 100 is 7.5000000000000009.
pub fn np_round(x: f64, decimals: u32) -> f64 {
    if decimals == 0 {
        return x.round_ties_even();
    }
    let f = 10f64.powi(decimals as i32);
    (x * f).round_ties_even() / f
}

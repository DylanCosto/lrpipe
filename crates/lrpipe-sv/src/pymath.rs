//! Numeric helpers that reproduce Python's (and numpy's) results bit for bit.
//!
//! Sniffles2 decides filters and formats INFO fields from these values, so an ulp of difference
//! can change a record.

use std::collections::HashMap;
use std::hash::Hash;

use num_bigint::BigUint;

/// `statistics.stdev` of integers (Python 3.11): the correctly rounded square root of the exact
/// sample variance. `None` for fewer than two values.
pub fn stdev(nums: &[i64]) -> Option<f64> {
    let n = nums.len();
    if n < 2 {
        return None;
    }
    let (sx, sxx) = nums.iter().fold((0i128, 0i128), |(s, ss), &x| (s + x as i128, ss + (x as i128) * (x as i128)));
    let n = n as i128;
    // variance = (n*sxx - sx^2) / (n*(n-1))
    let num = n * sxx - sx * sx;
    debug_assert!(num >= 0);
    let den = n * (n - 1);
    Some(float_sqrt_of_frac(BigUint::from(num as u128), BigUint::from(den as u128)))
}

/// `util.stdev` in Sniffles2: `statistics.stdev`, or the integer 0 for fewer than two values.
/// Which of the two it is matters for VCF formatting (`0` vs `0.000`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PyNum {
    Int(i64),
    Float(f64),
}

impl PyNum {
    pub fn as_f64(self) -> f64 {
        match self {
            PyNum::Int(i) => i as f64,
            PyNum::Float(f) => f,
        }
    }
}

pub fn util_stdev(nums: &[i64]) -> PyNum {
    match stdev(nums) {
        Some(s) => PyNum::Float(s),
        None => PyNum::Int(0),
    }
}

const SQRT_BIT_WIDTH: i64 = 2 * 53 + 3;

/// `statistics._float_sqrt_of_frac`: sqrt(n/m) correctly rounded.
fn float_sqrt_of_frac(n: BigUint, m: BigUint) -> f64 {
    let q = (n.bits() as i64 - m.bits() as i64 - SQRT_BIT_WIDTH).div_euclid(2);
    if q >= 0 {
        let numerator = integer_sqrt_of_frac_rto(&n, &(m << (2 * q) as usize));
        scale(numerator, q)
    } else {
        let numerator = integer_sqrt_of_frac_rto(&(n << (-2 * q) as usize), &m);
        scale(numerator, q)
    }
}

/// Square root of n/m rounded to an integer with round-to-odd.
fn integer_sqrt_of_frac_rto(n: &BigUint, m: &BigUint) -> BigUint {
    let a = (n / m).sqrt();
    if &a * &a * m != *n { a | BigUint::from(1u8) } else { a }
}

/// `numerator * 2^q` as a correctly rounded float (`numerator` has at most ~56 significant bits).
fn scale(numerator: BigUint, q: i64) -> f64 {
    let bits = numerator.bits();
    // Keep the integer conversion exact or correctly rounded: Rust's u128 -> f64 cast rounds to
    // nearest, ties to even, like Python's int / int.
    assert!(bits <= 128, "sqrt numerator out of range");
    let digits = numerator.to_u64_digits();
    let v: u128 = digits.iter().rev().fold(0u128, |acc, &d| (acc << 64) | d as u128);
    (v as f64) * 2f64.powi(q as i32)
}

/// Python's `round()` on a float: round half to even.
pub fn py_round(x: f64) -> i64 {
    x.round_ties_even() as i64
}

/// Python's `int()` on a float: truncate toward zero.
pub fn py_int(x: f64) -> i64 {
    x.trunc() as i64
}

/// `util.trim(nums, 25)`: sorted, with a quarter cut from each end.
pub fn trim(nums: &[i64]) -> Vec<i64> {
    let mut v = nums.to_vec();
    v.sort_unstable();
    let trim_n = (v.len() as f64 / 100.0 * 25.0) as usize;
    if trim_n > 0 { v[trim_n..v.len() - trim_n].to_vec() } else { v }
}

/// `util.median_noavg`: the upper median.
pub fn median_noavg(mut nums: Vec<i64>) -> i64 {
    nums.sort_unstable();
    nums[nums.len() / 2]
}

/// `util.median_modes` (`util.center`): upper median of the values whose count is within 2 of the
/// most common value's count.
pub fn center(nums: &[i64]) -> i64 {
    let mut counts: HashMap<i64, u32> = HashMap::new();
    let mut max_count = 0;
    for &n in nums {
        let c = counts.entry(n).or_insert(0);
        *c += 1;
        max_count = max_count.max(*c);
    }
    median_noavg(counts.into_iter().filter(|&(_, c)| max_count - c < 3).map(|(k, _)| k).collect())
}

/// `util.mean` of integers.
pub fn mean_int(nums: &[i64]) -> f64 {
    let s: i128 = nums.iter().map(|&x| x as i128).sum();
    s as f64 / nums.len() as f64
}

/// `util.most_common`: (count, value) pairs sorted descending (by count, then value).
pub fn most_common<T: Ord + Hash + Clone>(items: impl IntoIterator<Item = T>) -> Vec<(usize, T)> {
    let mut counts: HashMap<T, usize> = HashMap::new();
    for item in items {
        *counts.entry(item).or_insert(0) += 1;
    }
    let mut v: Vec<(usize, T)> = counts.into_iter().map(|(k, c)| (c, k)).collect();
    v.sort_unstable_by(|a, b| b.cmp(a));
    v
}

/// `util.most_common_top`: the smallest of the most common values.
pub fn most_common_top<T: Ord + Hash + Clone>(items: impl IntoIterator<Item = T>) -> T {
    let mc = most_common(items);
    let top = mc[0].0;
    mc.into_iter().filter(|(c, _)| *c == top).map(|(_, v)| v).min().expect("non-empty")
}

/// numpy's pairwise summation (`np.add.reduce` on a contiguous float64 array).
pub fn numpy_sum(a: &[f64]) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    // np.add.reduce starts from the first element and adds the pairwise sum of the rest.
    a[0] + pairwise_sum(&a[1..])
}

fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = -0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [0.0; 8];
        r.copy_from_slice(&a[..8]);
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

/// `np.nanmean` of a float list.
pub fn numpy_nanmean(a: &[f64]) -> f64 {
    let vals: Vec<f64> = a.iter().map(|&x| if x.is_nan() { 0.0 } else { x }).collect();
    let cnt = a.iter().filter(|x| !x.is_nan()).count();
    numpy_sum(&vals) / cnt as f64
}

/// Python's `math.log(x, 10)`.
pub fn log10_py(x: f64) -> f64 {
    x.ln() / 10f64.ln()
}

/// Python `f"{v:.3f}"`.
pub fn fmt3(v: f64) -> String {
    format!("{v:.3}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stdev_matches_python() {
        // Values from Python 3.11 statistics.stdev.
        assert_eq!(stdev(&[1, 2, 3, 4]), Some(1.2909944487358056));
        assert_eq!(stdev(&[5, 5, 5]), Some(0.0));
        assert_eq!(stdev(&[149021, 149095, 149100, 149044]), Some(38.738439135652676));
        assert_eq!(stdev(&[7]), None);
    }

    #[test]
    fn trim_uses_float_math() {
        // int(7/100.0*25) == 1 and int(12/100.0*25) == 3 in Python.
        let v: Vec<i64> = (0..7).collect();
        assert_eq!(trim(&v), (1..6).collect::<Vec<_>>());
        let v: Vec<i64> = (0..12).collect();
        assert_eq!(trim(&v), (3..9).collect::<Vec<_>>());
    }

    #[test]
    fn center_is_upper_median_of_modes() {
        assert_eq!(center(&[1, 1, 1, 1, 2, 3]), 1);
        assert_eq!(center(&[1, 1, 1, 2, 2, 3]), 2);
        assert_eq!(center(&[10, 10, 10, 10, 10, 20]), 10);
    }

    #[test]
    fn round_half_even() {
        assert_eq!(py_round(2.5), 2);
        assert_eq!(py_round(3.5), 4);
        assert_eq!(py_round(-0.5), 0);
    }
}

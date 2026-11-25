//! Python, numpy and scipy numerics reproduced bit for bit.
//!
//! A port of a Python tool gets the same output only if it adds, rounds, formats and draws
//! random numbers exactly as the original did. This crate collects those operations, each
//! checked against the reference environment (numpy 1.24, pandas 2.0, scipy 1.10, CPython
//! 3.10) by `tests/vectors.rs`:
//!
//! - [`sum`]: numpy's pairwise summation, and the means, variances, medians and percentiles
//!   built on it, with numpy's and pandas' NaN handling;
//! - [`round`]: numpy's rounding to decimals;
//! - [`repr`]: Python's shortest round-trip float text;
//! - [`random`]: CPython's `random.Random` (Mersenne Twister) seeding, `getrandbits`,
//!   `choice` and `sample`;
//! - [`numpy_random`]: numpy's legacy global `RandomState` uniform draws;
//! - [`special`]: scipy's normal distribution tail (`scipy.stats.norm.sf`, Cephes `ndtr`);
//! - [`json`]: `json.dump` text with an indent.

pub mod json;
pub mod numpy_random;
pub mod random;
pub mod repr;
pub mod round;
pub mod special;
pub mod sum;

mod mt19937;

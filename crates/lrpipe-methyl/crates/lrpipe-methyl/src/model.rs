//! A site's methylation level from its reads' 5mC values, by a fixed model of how methylated and
//! unmethylated reads' values are distributed.
//!
//! Each read at a CpG is methylated or not. Its 5mC value (an `ML` value, binned in sixteenths,
//! or "no call" where the read places no call at the C) is drawn from one distribution for
//! methylated reads and another for unmethylated ones; levels across sites follow a third. A
//! site's level is its posterior mean. Unlike counting values over one half, this uses how
//! informative each value is: a value of 0.6 says little, 0.99 a lot, and a missing call leans
//! unmethylated.
//!
//! The three distributions were fitted once, by expectation-maximisation with each site's level
//! known, on HG002 PacBio Revio reads (aligned to GRCh38) at chr20's CpGs against EM-seq
//! (EpiQC, GSE186383; 721,366 CpGs with at least 10 reads in both), and checked on the other
//! chromosomes only. They are fixed rather than learned from each input: learned per input, the
//! two read distributions are not identifiable (a no-call-heavy "unmethylated" component fits
//! as well as the real one), and the levels moved with the sites sampled. This is not PacBio's
//! trained network, and like it, it describes the chemistry it was fitted on.

use std::sync::LazyLock;

/// Value categories: sixteen bins of `ML` values, and no call.
pub const CATEGORIES: usize = 17;

/// The value standing for no call (one past the largest `ML` value).
pub const NO_CALL: u16 = 256;

/// Levels a site's estimate considers: 0 to 1 in hundredths.
const GRID: usize = 101;

/// How often methylated reads' values fall in each category.
const METHYLATED: [f64; CATEGORIES] = [
    0.00922624, 0.0158128, 0.0191379, 0.0199117, 0.020772, 0.0227244, 0.0241572, 0.0268329, 0.0307083, 0.0361783,
    0.0434482, 0.0515482, 0.06092, 0.0763616, 0.111913, 0.326791, 0.103556,
];

/// How often unmethylated reads' values fall in each category.
const UNMETHYLATED: [f64; CATEGORIES] = [
    0.278018, 0.11487, 0.0810835, 0.0564294, 0.0430801, 0.0362381, 0.0304212, 0.0270665, 0.0245334, 0.0225813,
    0.0207899, 0.0182051, 0.0156848, 0.0134966, 0.0117373, 0.00896345, 0.196802,
];

/// How sites' levels are distributed, at each level of [`GRID`].
const PRIOR: [f64; GRID] = [
    0.100032,
    6.3906e-05,
    0.00122872,
    0.00907988,
    0.0110417,
    0.00733,
    0.00604939,
    0.00578632,
    0.00508086,
    0.0044858,
    0.00435306,
    0.00374368,
    0.0037706,
    0.003585,
    0.00325717,
    0.00336132,
    0.00329071,
    0.00318417,
    0.0031866,
    0.00303212,
    0.00336779,
    0.00304607,
    0.00301771,
    0.00314037,
    0.00309676,
    0.00362402,
    0.00305678,
    0.00328712,
    0.00339732,
    0.00338213,
    0.00380549,
    0.00345289,
    0.00396737,
    0.00400272,
    0.00353906,
    0.00444725,
    0.00419028,
    0.00428535,
    0.00445159,
    0.00423711,
    0.00518433,
    0.00465919,
    0.00481917,
    0.00490256,
    0.00494753,
    0.00504386,
    0.00532964,
    0.00497633,
    0.00655249,
    0.00184937,
    0.0104897,
    0.00185793,
    0.00650698,
    0.00505818,
    0.00538874,
    0.00534038,
    0.00534285,
    0.00529562,
    0.00533791,
    0.00498093,
    0.00591615,
    0.00514652,
    0.00528635,
    0.00537295,
    0.00551755,
    0.00580203,
    0.00475142,
    0.00566629,
    0.00610944,
    0.00560768,
    0.00613618,
    0.00585132,
    0.00631001,
    0.00629572,
    0.00624696,
    0.0080051,
    0.00681829,
    0.00746078,
    0.00771022,
    0.00792978,
    0.0100478,
    0.00900641,
    0.0102836,
    0.0105275,
    0.0120019,
    0.0129092,
    0.0133244,
    0.0152651,
    0.016847,
    0.0180995,
    0.021596,
    0.0228751,
    0.0270381,
    0.02949,
    0.0280042,
    0.0319858,
    0.054207,
    0.0443514,
    0.00490206,
    6.97804e-06,
    0.121688,
];

/// `ln(θ m(v) + (1 - θ) u(v))` for each level θ of the grid and category v, row by row.
static LOG_LIKELIHOOD: LazyLock<Vec<f64>> = LazyLock::new(|| {
    (0..GRID)
        .flat_map(|g| {
            let theta = g as f64 / (GRID - 1) as f64;
            (0..CATEGORIES).map(move |v| (theta * METHYLATED[v] + (1.0 - theta) * UNMETHYLATED[v]).ln())
        })
        .collect()
});

/// The counted reads' values at one site, by category.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Histogram([u16; CATEGORIES]);

impl Histogram {
    /// Adds a read's value: an `ML` value (0 to 255), or [`NO_CALL`].
    ///
    /// # Panics
    ///
    /// If `value` is above [`NO_CALL`].
    pub fn add(&mut self, value: u16) {
        assert!(value <= NO_CALL, "value {value} out of range");
        let bin = if value == NO_CALL { CATEGORIES - 1 } else { usize::from(value) / 16 };
        self.0[bin] = self.0[bin].saturating_add(1);
    }

    /// Whether no value was added.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.iter().all(|&n| n == 0)
    }

    /// The site's methylated fraction (0 to 1): its posterior mean level; 0 with no values.
    #[must_use]
    pub fn level(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let table = &*LOG_LIKELIHOOD;
        let mut log = [0.0; GRID];
        for (g, (l, p)) in log.iter_mut().zip(PRIOR).enumerate() {
            let row = &table[g * CATEGORIES..(g + 1) * CATEGORIES];
            *l = p.ln() + self.0.iter().zip(row).filter(|(n, _)| **n > 0).map(|(&n, l)| f64::from(n) * l).sum::<f64>();
        }
        let top = log.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let (mut weight, mut total) = (0.0, 0.0);
        for (g, l) in log.iter().enumerate() {
            let w = (l - top).exp();
            weight += w * g as f64 / (GRID - 1) as f64;
            total += w;
        }
        weight / total
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(values: &[(u16, usize)]) -> Histogram {
        let mut h = Histogram::default();
        for &(v, n) in values {
            for _ in 0..n {
                h.add(v);
            }
        }
        h
    }

    #[test]
    fn tables_are_distributions() {
        for t in [&METHYLATED[..], &UNMETHYLATED[..], &PRIOR[..]] {
            assert!((t.iter().sum::<f64>() - 1.0).abs() < 1e-4);
            assert!(t.iter().all(|&p| p > 0.0));
        }
    }

    #[test]
    fn confident_reads_give_their_share() {
        assert!(site(&[(250, 20)]).level() > 0.95);
        assert!(site(&[(3, 20)]).level() < 0.05);
        let half = site(&[(250, 10), (3, 10)]).level();
        assert!((0.4..0.6).contains(&half), "{half}");
    }

    #[test]
    fn uncertain_values_say_less_than_confident_ones() {
        // Ten reads just over one half are weaker evidence than ten near one.
        assert!(site(&[(140, 10)]).level() < site(&[(250, 10)]).level());
        // No calls lean unmethylated.
        assert!(site(&[(250, 10), (NO_CALL, 5)]).level() < site(&[(250, 10)]).level());
    }

    #[test]
    fn levels_rise_with_values() {
        let levels: Vec<f64> = (0..=255).step_by(16).map(|v| site(&[(v, 12)]).level()).collect();
        assert!(levels.windows(2).all(|w| w[0] <= w[1]), "{levels:?}");
    }

    #[test]
    fn empty_sites_have_no_level() {
        assert_eq!(Histogram::default().level(), 0.0);
    }
}

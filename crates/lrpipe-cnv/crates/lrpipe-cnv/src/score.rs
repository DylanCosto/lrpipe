//! Quality scores: a z-test of a call's depth against the window sample (spec §3.10;
//! `CNVMetrics.evaluate_cnvs` and `get_z_score`, cnv_metrics.py:160-199, 452-486).

use pynum::{special::norm_sf, sum};

use crate::number::Number;
use crate::thresholds::WindowSample;

/// `statistics["z-score"]` of a candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZScore {
    /// The z-score capped to ±10 (a Python int when capped).
    pub score: Number,
    pub pvalue: Number,
    /// The genotype quality: `min(60, int(-10 log10 p))`.
    pub sample_score: i64,
}

impl ZScore {
    /// What every candidate gets before refinement: zeros.
    pub fn unrefined() -> Self {
        ZScore { score: Number::Int(0), pvalue: Number::Int(0), sample_score: 0 }
    }

    /// The z-test of a refined call's finite depths against the window sample.
    pub fn of(cov: &[f64], sample: &WindowSample) -> Self {
        let x: Vec<f64> = cov.iter().copied().filter(|v| v.is_finite()).collect();
        let (mean, sd) = (sum::mean(&x), sum::std(&x));
        let n = sample.windows.len() as f64;
        let z = (sample.mean - mean) / (sd / n.sqrt());
        // `min(z, 10) if z > 1 else max(-10, z)`: Python keeps the int when it wins, and NaN
        // never wins a comparison.
        let capped = if z > 1.0 {
            if z > 10.0 { Number::Int(10) } else { Number::Float(z) }
        } else if z > -10.0 {
            Number::Float(z)
        } else {
            Number::Int(-10)
        };
        let pvalue = norm_sf(capped.value().abs()) * 2.0;
        let gq = -pvalue.log10() * 10.0;
        let gq = if gq.is_finite() { gq as i64 } else { i64::MAX };
        ZScore { score: capped, pvalue: Number::Float(pvalue), sample_score: gq.min(60) }
    }
}

//! DEL/DUP thresholds from a seeded random sample of windows (spec §3.4; `CNVMetrics`,
//! cnv_metrics.py:30-135, 386-450).

use md5::{Digest, Md5};
use pynum::{random::Random, sum};

use crate::normalize::Normalized;
use crate::number::Number;

/// Whether thresholds are clamped to ploidy ± 1 (`--cancer` turns clamping off).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Clamping {
    #[default]
    PloidyBounds,
    Unclamped,
}

/// One sampled window.
#[derive(Clone, Debug, PartialEq)]
pub struct SampledWindow {
    pub contig: String,
    pub position: i64,
    /// Its normalized, ploidy-scaled depth.
    pub coverage: f64,
}

/// The random sample of windows thresholds and quality scores are computed from.
///
/// Spectre recomputes it whenever CNV calls change, meaning to exclude blacklisted and called
/// regions; because the exclusion mask is taken from an unfiltered copy, nothing is ever
/// excluded and every recomputation gives this same sample (quirk Q-EXCL).
#[derive(Clone, Debug, PartialEq)]
pub struct WindowSample {
    pub windows: Vec<SampledWindow>,
    /// `float(np.mean(sample))`.
    pub mean: f64,
    /// `float(np.std(sample))`: pandas, ddof 0.
    pub sd: f64,
    /// The draw made: `sample(range(population), drawn)`. Spectre redraws it (identically)
    /// before each scoring pass, which leaves Python's `random` where this draw left it.
    pub population: u64,
    pub drawn: usize,
}

impl WindowSample {
    /// Python's `random` as Spectre leaves it after a (re)draw of the window sample: where
    /// LoH IDs are drawn from.
    pub fn redrawn_rng(&self, sample_id: &str) -> Random {
        let mut rng = seed_rng(sample_id);
        rng.sample_range(self.population, self.drawn);
        rng
    }
}

impl WindowSample {
    /// `random_sample_selection`: reseeds `rng` from the sample ID ([`seed_rng`]) and draws
    /// `int((M - 1) * 0.1)` of the first `M - 1` non-NaN, finite windows (the last one is
    /// never chosen), keeping them in genome order. Windows above
    /// [`Normalized::kept_above`] are left out. Later draws (candidate IDs) continue
    /// from `rng`.
    pub fn draw(norm: &Normalized, sample_id: &str, rng: &mut Random) -> Self {
        let rows: Vec<SampledWindow> = norm
            .contigs
            .iter()
            .flat_map(|c| {
                c.positions
                    .iter()
                    .zip(&c.normalized_ploidy)
                    .zip(&c.raw)
                    .filter(|((_, v), raw)| v.is_finite() && norm.kept_above.is_none_or(|cut| **raw <= cut))
                    .map(|((&position, &coverage), _)| SampledWindow { contig: c.name.clone(), position, coverage })
            })
            .collect();
        let last = rows.len().saturating_sub(1);
        let k = (last as f64 * 0.1) as usize;
        *rng = seed_rng(sample_id);
        let mut picks = rng.sample_range(last as u64, k);
        picks.sort_unstable();
        let windows: Vec<SampledWindow> = picks.into_iter().map(|i| rows[i as usize].clone()).collect();
        let cov: Vec<f64> = windows.iter().map(|w| w.coverage).collect();
        WindowSample {
            mean: sum::series_mean(&cov),
            sd: sum::series_var(&cov, 0).sqrt(),
            windows,
            population: last as u64,
            drawn: k,
        }
    }
}

/// `random.seed(hashlib.md5(sample_id.encode()).hexdigest())`.
pub fn seed_rng(sample_id: &str) -> Random {
    let seed: String = Md5::digest(sample_id.as_bytes()).iter().map(|b| format!("{b:02x}")).collect();
    Random::from_str_seed(&seed)
}

/// The calling thresholds on normalized, ploidy-scaled depth.
#[derive(Clone, Debug, PartialEq)]
pub struct Thresholds {
    pub del: Number,
    pub dup: Number,
}

impl Thresholds {
    /// `calculate_del_dup_borders_quantile(0.1, 0.9)`: the 10th and 90th percentiles of the
    /// sampled depths below mean + 2 SD, clamped to `ploidy - 1` and `ploidy + 1`.
    pub fn from_quantiles(sample: &WindowSample, ploidy: i64, clamping: Clamping) -> Self {
        let clip = sample.mean + 2.0 * sample.sd;
        let clipped: Vec<f64> = sample.windows.iter().map(|w| w.coverage).filter(|&v| v < clip).collect();
        // pandas' Series.quantile(q) is np.percentile(values, q * 100.0); NaN when empty.
        let quantile = |q: f64| if clipped.is_empty() { f64::NAN } else { sum::percentile(&clipped, q * 100.0) };
        clamp(Number::Float(quantile(0.1)), Number::Float(quantile(0.9)), ploidy, clamping)
    }

    /// `calculate_del_dup_borders_sd(sd_clip_n)`: mean ∓ `sd_clip_n` SD, clamped the same way.
    /// Spectre's clean-up step uses it with 1.85 (the strict thresholds).
    pub fn from_sd(sample: &WindowSample, sd_clip_n: f64, ploidy: i64, clamping: Clamping) -> Self {
        let del = sample.mean - sample.sd * sd_clip_n;
        let dup = sample.mean + sample.sd * sd_clip_n;
        clamp(Number::Float(del), Number::Float(dup), ploidy, clamping)
    }
}

/// `check_del_dup_borders`: out-of-range thresholds become the Python ints `ploidy ∓ 1`.
fn clamp(del: Number, dup: Number, ploidy: i64, clamping: Clamping) -> Thresholds {
    if clamping == Clamping::Unclamped {
        return Thresholds { del, dup };
    }
    let del = if del.value() < (ploidy - 1) as f64 { Number::Int(ploidy - 1) } else { del };
    let dup = if dup.value() > (ploidy + 1) as f64 { Number::Int(ploidy + 1) } else { dup };
    Thresholds { del, dup }
}

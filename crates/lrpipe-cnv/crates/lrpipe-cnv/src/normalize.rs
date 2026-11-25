//! Masking, genome-wide statistics and per-contig normalization (spec §3.2, §3.3;
//! `CNVAnalysis.data_normalization`, analysis.py:104-291).

use pynum::sum;

use crate::Error;
use crate::coverage::ContigDepth;
use crate::genome::Genome;
use crate::metadata::Metadata;

/// Settings for normalization.
#[derive(Clone, Debug)]
pub struct NormalizeOptions {
    /// `--sample-coverage-overwrite`: replaces the genome-wide median (not used with SNVs).
    pub sample_coverage_overwrite: Option<f64>,
    /// `--dev-max-std-outlier-rm` (default 5): windows deeper than mean + this many standard
    /// deviations are masked.
    pub max_std_outlier_rm: f64,
    /// `--disable-max-coverage`: skip that outlier masking.
    pub outlier_removal: OutlierRemoval,
    /// The copy-number-neutral depth estimated from heterozygous SNVs, which replaces the
    /// genome-wide median as the normalization value (spec §3.1).
    pub snv_cn_neutral: Option<f64>,
}

/// Whether windows far above the genome-wide mean are masked.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OutlierRemoval {
    #[default]
    Enabled,
    Disabled,
    /// Only windows in short runs above the cutoff (fewer than [`MIN_AMPLIFICATION_WINDOWS`]
    /// such windows, NaN windows between them allowed), so that amplifications stay (an
    /// improvement over Spectre, which masks them).
    IsolatedWindows,
}

/// Runs of this many windows above the outlier cutoff are kept by
/// [`OutlierRemoval::IsolatedWindows`].
pub const MIN_AMPLIFICATION_WINDOWS: usize = 25;

impl Default for NormalizeOptions {
    fn default() -> Self {
        NormalizeOptions {
            sample_coverage_overwrite: None,
            max_std_outlier_rm: 5.0,
            outlier_removal: OutlierRemoval::Enabled,
            snv_cn_neutral: None,
        }
    }
}

/// `CoverageStatistics`: summary values of one contig's depths (NaN-aware).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stats {
    pub average: f64,
    pub std_dev: f64,
    pub median: f64,
    pub min: f64,
    pub max: f64,
}

impl Stats {
    fn of(v: &[f64]) -> Self {
        let kept = v.iter().copied().filter(|x| !x.is_nan());
        Stats {
            average: sum::nanmean(v),
            std_dev: sum::nanstd(v),
            median: sum::nanmedian(v),
            min: kept.clone().fold(f64::INFINITY, f64::min),
            max: kept.fold(f64::NEG_INFINITY, f64::max),
        }
    }
}

/// One analysed contig (`coverage_analysis[contig]`).
#[derive(Clone, Debug, PartialEq)]
pub struct ContigCoverage {
    pub name: String,
    pub positions: Vec<i64>,
    /// Depths with N runs, blacklist and outliers set to NaN (`coverage_raw`).
    pub raw: Vec<f64>,
    /// `raw / normalization value` (`normalized_cov`).
    pub normalized: Vec<f64>,
    /// `normalized * ploidy multiplier` (`normalized_cov_ploidy`): what calling uses.
    pub normalized_ploidy: Vec<f64>,
    pub stats: Stats,
    pub norm_stats: Stats,
}

/// The result of normalization.
#[derive(Clone, Debug, PartialEq)]
pub struct Normalized {
    /// Window size, from the first two windows of the first contig with coverage.
    pub bin_size: i64,
    pub genome_wide_median: f64,
    pub genome_wide_mean: f64,
    pub genome_std: f64,
    /// What depths were divided by: the SNV estimate or the genome-wide median.
    pub normalization_value: f64,
    /// The analysed contigs in `--only-chr` order, without contigs that had no usable window.
    pub contigs: Vec<ContigCoverage>,
    /// With [`OutlierRemoval::IsolatedWindows`], the outlier cutoff on raw depth: windows above
    /// it are kept for calling but left out of the threshold sample.
    pub kept_above: Option<f64>,
}

/// Loads, masks and normalizes depth. `load_metadata` is called with the bin size once it
/// is known, as Spectre loads the metadata only then.
///
/// # Errors
///
/// If no genome contig has coverage, the first one has a single window, or an `only_chr`
/// contig has no coverage (a KeyError in Spectre).
pub fn normalize(
    genome: &Genome,
    depths: Vec<ContigDepth>,
    only_chr: &[String],
    opts: &NormalizeOptions,
    load_metadata: impl FnOnce(i64) -> Result<Metadata, Error>,
) -> Result<(Normalized, Metadata), Error> {
    // analysis.py:117-150: contigs in reference order, those without coverage skipped.
    let mut by_name: Vec<ContigDepth> = Vec::new();
    let mut depths = depths;
    for c in &genome.contigs {
        if let Some(i) = depths.iter().position(|d| d.name == c.name) {
            by_name.push(depths.swap_remove(i));
        }
    }
    let first = by_name
        .first()
        .ok_or_else(|| Error::Input("no coverage data found for any chromosome in the reference".into()))?;
    if first.positions.len() < 2 {
        return Err(Error::Input(format!("{}: cannot determine the bin size from one window", first.name)));
    }
    let bin_size = (first.positions[1] - first.positions[0]).abs();
    let metadata = load_metadata(bin_size)?;
    for d in &mut by_name {
        mask(d, &metadata);
    }

    // analysis.py:161-172: statistics over every window of every contig, masked.
    let all: Vec<f64> = by_name.iter().flat_map(|d| d.depth.iter().copied()).collect();
    let genome_wide_mean = sum::nanmean(&all);
    let genome_std = sum::nanstd(&all);
    let mut genome_wide_median = sum::nanmedian(&all);
    if let Some(x) = opts.sample_coverage_overwrite {
        genome_wide_median = x;
    }
    let normalization_value = opts.snv_cn_neutral.unwrap_or(genome_wide_median);
    let cutoff = genome_wide_mean + opts.max_std_outlier_rm * genome_std;

    let mut contigs = Vec::new();
    for name in only_chr {
        let d = by_name
            .iter_mut()
            .find(|d| &d.name == name)
            .ok_or_else(|| Error::Input(format!("no coverage for --only-chr contig {name}")))?;
        let contig = genome.contig(name).ok_or_else(|| Error::Input(format!("{name} is not a reference contig")))?;
        // analysis.py:239-243, in place: later contigs see the masked depths too.
        match opts.outlier_removal {
            OutlierRemoval::Enabled => {
                for v in &mut d.depth {
                    if *v > cutoff {
                        *v = f64::NAN;
                    }
                }
            }
            OutlierRemoval::IsolatedWindows => mask_isolated(&mut d.depth, cutoff),
            OutlierRemoval::Disabled => {}
        }
        if d.depth.is_empty() || d.depth.iter().all(|v| v.is_nan()) {
            continue;
        }
        let normalized: Vec<f64> = d.depth.iter().map(|&v| v / normalization_value).collect();
        let normalized_ploidy: Vec<f64> = normalized.iter().map(|&v| v * contig.ploidy_multiplier as f64).collect();
        contigs.push(ContigCoverage {
            name: name.clone(),
            positions: d.positions.clone(),
            raw: d.depth.clone(),
            stats: Stats::of(&d.depth),
            norm_stats: Stats::of(&normalized_ploidy),
            normalized,
            normalized_ploidy,
        });
    }
    let kept_above = (opts.outlier_removal == OutlierRemoval::IsolatedWindows).then_some(cutoff);
    let normalized = Normalized {
        bin_size,
        genome_wide_median,
        genome_wide_mean,
        genome_std,
        normalization_value,
        contigs,
        kept_above,
    };
    Ok((normalized, metadata))
}

/// Masks windows above `cutoff` in runs of fewer than [`MIN_AMPLIFICATION_WINDOWS`] such
/// windows; NaN windows inside a run do not end it.
fn mask_isolated(depth: &mut [f64], cutoff: f64) {
    let mut run: Vec<usize> = Vec::new();
    let flush = |run: &mut Vec<usize>, depth: &mut [f64]| {
        if run.len() < MIN_AMPLIFICATION_WINDOWS {
            for &i in run.iter() {
                depth[i] = f64::NAN;
            }
        }
        run.clear();
    };
    for i in 0..depth.len() {
        if depth[i] > cutoff {
            run.push(i);
        } else if !depth[i].is_nan() {
            flush(&mut run, depth);
        }
    }
    flush(&mut run, depth);
}

/// analysis.py:293-308: windows whose start lies in `start..=end` of an excluded region.
fn mask(d: &mut ContigDepth, metadata: &Metadata) {
    for r in metadata.regions(&d.name) {
        for (p, v) in d.positions.iter().zip(d.depth.iter_mut()) {
            if *p >= r.start && *p <= r.end {
                *v = f64::NAN;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isolated_outliers_masked_runs_kept() {
        let mut depth = vec![30.0; 100];
        depth[5] = 200.0;
        // A NaN window inside a run does not end it, nor count towards it.
        let run = 40..41 + MIN_AMPLIFICATION_WINDOWS;
        for v in &mut depth[run.clone()] {
            *v = 90.0;
        }
        depth[50] = f64::NAN;
        mask_isolated(&mut depth, 80.0);
        assert!(depth[5].is_nan());
        assert!(depth[run].iter().enumerate().all(|(i, v)| i == 10 || *v == 90.0));
        assert_eq!(depth.iter().filter(|v| v.is_nan()).count(), 2);
    }
}

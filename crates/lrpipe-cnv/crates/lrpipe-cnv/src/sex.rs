//! The sample's sex from sex-chromosome depth (an improvement over Spectre, which treats
//! chrX and chrY like autosomes unless told their ploidy).

use pynum::sum;

use crate::coverage::ContigDepth;

/// The sample's sex, as given or as read from chrX and chrY depth.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    /// One X, one Y.
    Male,
    /// Two X, no Y.
    Female,
}

/// The names of the X and Y chromosomes among `names`, if present.
pub fn sex_contigs<'a>(names: impl IntoIterator<Item = &'a str>) -> (Option<String>, Option<String>) {
    let names: Vec<&str> = names.into_iter().collect();
    let find = |cands: &[&str]| cands.iter().find(|c| names.contains(c)).map(|c| (*c).to_owned());
    (find(&["chrX", "X"]), find(&["chrY", "Y"]))
}

/// Median depth of windows with reads (NaN windows skipped), or 0 without any.
fn covered_median<'a>(depths: impl IntoIterator<Item = &'a f64>) -> f64 {
    let v: Vec<f64> = depths.into_iter().copied().filter(|&x| x > 0.0).collect();
    if v.is_empty() { 0.0 } else { sum::median(&v) }
}

/// Depth of chrX and chrY relative to the autosomes, each a median over windows with reads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SexDepth {
    pub x: f64,
    pub y: f64,
}

impl SexDepth {
    /// Measures the ratios. `x` and `y` name the sex chromosomes; `autosomes` the contigs to
    /// compare with. A missing chromosome has ratio 0; without autosome reads both are NaN.
    pub fn measure(depths: &[ContigDepth], x: Option<&str>, y: Option<&str>, autosomes: &[String]) -> Self {
        let median = covered_median(depths.iter().filter(|d| autosomes.contains(&d.name)).flat_map(|d| &d.depth));
        let ratio = |name: Option<&str>| {
            let d = name.and_then(|n| depths.iter().find(|d| d.name == n));
            d.map_or(0.0, |d| covered_median(&d.depth)) / median
        };
        SexDepth { x: ratio(x), y: ratio(y) }
    }

    /// The sex these ratios show: X near half the autosomes' depth with Y present is male,
    /// X near the autosomes' with (nearly) no Y female; anything else is not clear.
    pub fn sex(self) -> Option<Sex> {
        if (0.3..0.7).contains(&self.x) && self.y > 0.2 {
            Some(Sex::Male)
        } else if self.x > 0.8 && self.y < 0.1 {
            Some(Sex::Female)
        } else {
            None
        }
    }
}

/// The pseudoautosomal regions of chrX (0-based, half-open), recognized by chrX's length
/// (GRCh38, GRCh37 or T2T-CHM13, as every lrpipe tool recognizes them): diploid in males too.
/// lrpipe-cnv had its own table without CHM13, so a male's CHM13 PARs counted as haploid and
/// their depth doubled (LRP-018).
pub fn x_pars(x_len: u64) -> Vec<(i64, i64)> {
    lrpipe_reads::par::Build::from_x_length(x_len).map_or_else(Vec::new, |b| b.x().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(name: &str, depth: f64) -> ContigDepth {
        ContigDepth { name: name.into(), positions: (0..10).map(|i| i * 1000).collect(), depth: vec![depth; 10] }
    }

    #[test]
    fn reads_sex_from_depth_ratios() {
        let autosomes = vec!["chr1".to_owned()];
        let measure = |x: f64, y: f64| {
            let depths = [contig("chr1", 30.0), contig("chrX", x), contig("chrY", y)];
            SexDepth::measure(&depths, Some("chrX"), Some("chrY"), &autosomes).sex()
        };
        assert_eq!(measure(15.0, 14.0), Some(Sex::Male));
        assert_eq!(measure(30.0, 0.5), Some(Sex::Female));
        assert_eq!(measure(22.0, 1.0), None);
    }
}

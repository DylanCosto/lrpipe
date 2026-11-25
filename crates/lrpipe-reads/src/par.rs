//! The pseudoautosomal regions (PARs) of the human sex chromosomes, recognized by the length of
//! chrX (or X) in the reference: diploid in XY samples too. Coordinates are 0-based, half-open.

/// A human reference whose PARs are known.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Build {
    /// GRCh38 (hg38).
    Grch38,
    /// GRCh37 (hg19, b37).
    Grch37,
    /// T2T-CHM13 v2.0 (hs1).
    Chm13,
}

impl Build {
    /// The build whose chrX is `x_len` long, if any.
    #[must_use]
    pub fn from_x_length(x_len: u64) -> Option<Self> {
        match x_len {
            156_040_895 => Some(Self::Grch38),
            155_270_560 => Some(Self::Grch37),
            154_259_566 => Some(Self::Chm13),
            _ => None,
        }
    }

    /// chrX's PARs.
    #[must_use]
    pub fn x(self) -> [(i64, i64); 2] {
        match self {
            Self::Grch38 => [(10_000, 2_781_479), (155_701_382, 156_030_895)],
            Self::Grch37 => [(60_000, 2_699_520), (154_931_043, 155_260_560)],
            Self::Chm13 => [(0, 2_394_410), (153_925_834, 154_259_566)],
        }
    }

    /// chrY's PARs.
    #[must_use]
    pub fn y(self) -> [(i64, i64); 2] {
        match self {
            Self::Grch38 => [(10_000, 2_781_479), (56_887_902, 57_217_415)],
            Self::Grch37 => [(10_000, 2_649_520), (59_034_049, 59_363_566)],
            Self::Chm13 => [(0, 2_458_320), (62_122_809, 62_460_029)],
        }
    }

    /// The PARs as a BED (`chrX` and `chrY`, or `X` and `Y` when `prefix` is empty).
    #[must_use]
    pub fn bed(self, prefix: &str) -> String {
        let mut out = String::new();
        for (contig, pars) in [("X", self.x()), ("Y", self.y())] {
            for (start, end) in pars {
                out.push_str(&format!("{prefix}{contig}\t{start}\t{end}\n"));
            }
        }
        out
    }
}

/// Whether `[start, end)` lies inside one of `pars`.
#[must_use]
pub fn within(pars: &[(i64, i64)], start: i64, end: i64) -> bool {
    pars.iter().any(|&(s, e)| start >= s && end <= e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_by_chrx_length() {
        assert_eq!(Build::from_x_length(156_040_895), Some(Build::Grch38));
        assert_eq!(Build::from_x_length(155_270_560), Some(Build::Grch37));
        assert_eq!(Build::from_x_length(1), None);
        assert!(within(&Build::Grch37.x(), 155_000_000, 155_000_030), "GRCh37 PAR2");
        assert!(!within(&Build::Grch38.x(), 155_000_000, 155_000_030), "not GRCh38's");
        assert_eq!(Build::Grch38.bed("chr").lines().next(), Some("chrX\t10000\t2781479"));
    }
}

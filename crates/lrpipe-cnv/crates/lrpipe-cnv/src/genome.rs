//! The contigs analysed and their ploidy multipliers (spec §2.2, §2.5).

use std::path::Path;

use crate::Error;

/// A reference contig kept for analysis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub len: u64,
    /// The coverage multiplier for this contig: the global ploidy, or `--ploidy-chr`'s
    /// `int(ploidy² / p)` (so `chrX:1` with ploidy 2 multiplies by 4, see spec §2.5).
    pub ploidy_multiplier: i64,
}

/// The reference contigs longer than the minimum length, in FASTA index order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Genome {
    pub contigs: Vec<Contig>,
    pub ploidy: i64,
}

impl Genome {
    /// Reads contig names and lengths from the reference's `.fai` and keeps contigs longer
    /// than `min_chr_len` (strictly, spectre.py:185). `ploidy_chr` holds `--ploidy-chr`
    /// pairs; ones naming a contig not kept are skipped, as Spectre warns and skips them.
    ///
    /// # Errors
    ///
    /// If the index cannot be read.
    pub fn from_fasta_index(
        reference: &Path,
        min_chr_len: u64,
        ploidy: i64,
        ploidy_chr: &[(String, i64)],
    ) -> Result<Self, Error> {
        let fai_path = reference.with_added_extension("fai");
        let index = pasta::fasta::fai::fs::read(&fai_path).map_err(Error::io(&fai_path))?;
        let mut contigs: Vec<Contig> = index
            .as_ref()
            .iter()
            .filter(|r| r.length() > min_chr_len)
            .map(|r| Contig { name: r.name().to_string(), len: r.length(), ploidy_multiplier: ploidy })
            .collect();
        for (name, p) in ploidy_chr {
            if let Some(c) = contigs.iter_mut().find(|c| &c.name == name) {
                // spectre.py:299-300: `int((1 / (int(p) / ploidy)) * ploidy)`.
                c.ploidy_multiplier = ((1.0 / (*p as f64 / ploidy as f64)) * ploidy as f64) as i64;
            }
        }
        Ok(Genome { contigs, ploidy })
    }

    pub fn contig(&self, name: &str) -> Option<&Contig> {
        self.contigs.iter().find(|c| c.name == name)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn ploidy_chr_is_a_multiplier() {
        let m = |p: i64, ploidy: i64| ((1.0 / (p as f64 / ploidy as f64)) * ploidy as f64) as i64;
        assert_eq!([m(1, 2), m(2, 2), m(3, 2), m(4, 2)], [4, 2, 1, 1]);
    }
}

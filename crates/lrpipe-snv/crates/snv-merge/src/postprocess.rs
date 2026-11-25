//! `postprocess_row_with_params` (`SortVcf.py`) and the helpers it takes from `MergeVcf.py`:
//! haploid genotypes, `LowQual` marking, RefCall filtering and the PAR-aware haploid decider.

use std::path::Path;

use crate::bed::{BedError, BedTree, open_pigz};
use crate::pytext::{parse_float, parse_int, split_whitespace, strip};
use crate::{Error, Gender, Haploid, RefCalls, RowError};

/// `SEX_CHROMOSOMES` in `MergeVcf.py`.
const SEX_CHROMOSOMES: [&str; 4] = ["chrX", "X", "chrY", "Y"];

/// The pseudo-autosomal regions `bed_tree_from(par_regions_bed)` loads.
#[derive(Debug, Clone, Default)]
pub(crate) struct ParRegions(BedTree);

impl ParRegions {
    /// Reads a plain or gzipped BED file as `pigz -fdc` feeds it to `bed_tree_from`.
    ///
    /// A file that cannot be opened gives no regions: `pigz` fails, the original reads an empty
    /// stream and never notices.
    pub(crate) fn read(path: &Path) -> Result<Self, Error> {
        let reader = match open_pigz(path) {
            Ok(reader) => reader,
            Err(e) => {
                log::warn!("cannot read PAR regions {}: {e}; no PAR regions applied", path.display());
                return Ok(Self::default());
            }
        };
        BedTree::read(reader, None, 0).map(Self).map_err(|e| match e {
            BedError::Io(source) => Error::Io { path: path.to_owned(), source },
            BedError::Row { line, reason } => Error::ParBed { path: path.to_owned(), line, reason },
        })
    }
}

/// `haploid_mode_decider(args)`: the haploid mode of each record.
#[derive(Debug, Clone)]
pub(crate) struct HaploidDecider {
    haploid: Haploid,
    gender: Gender,
    par: ParRegions,
}

impl HaploidDecider {
    pub(crate) fn new(haploid: Haploid, gender: Gender, par: ParRegions) -> Self {
        Self { haploid, gender, par }
    }

    /// `decide(contig, pos)`: the global modes first, then `--gender male` on the sex
    /// chromosomes outside PAR regions. `pos` is the 1-based VCF position, so a PAR BED row
    /// `[start, end)` covers positions `start + 1..=end`. The original compares the 1-based
    /// position with the 0-based BED, one base off at both ends (LRP-003).
    pub(crate) fn decide(&self, contig: &str, pos: i64) -> Haploid {
        match self.haploid {
            Haploid::Off => {}
            mode => return mode,
        }
        if self.gender == Gender::Male && SEX_CHROMOSOMES.contains(&contig) {
            if !self.par.0.is_empty() && self.par.0.contains(contig, pos - 1) {
                return Haploid::Off;
            }
            return Haploid::Precise;
        }
        Haploid::Off
    }

    /// The first 1-based position after `pos` at which ploidy changes.
    pub(crate) fn next_boundary(&self, contig: &str, pos: i64) -> Option<i64> {
        if self.haploid == Haploid::Off && self.gender == Gender::Male && SEX_CHROMOSOMES.contains(&contig) {
            return self.par.0.next_boundary(contig, pos - 1)?.checked_add(1);
        }
        None
    }
}

/// The `postprocess_row_with_params` settings.
#[derive(Debug, Clone)]
pub(crate) struct RowFilter {
    pub(crate) ref_calls: RefCalls,
    pub(crate) qual: i32,
    pub(crate) decider: HaploidDecider,
}

fn column<'a>(columns: &[&'a str], index: usize) -> Result<&'a str, RowError> {
    columns.get(index).copied().ok_or(RowError::MissingColumn { index })
}

/// `update_haploid_precise_genotype`: `1/1` becomes `1`, `0/0` becomes `0`, anything else
/// gives an empty row.
pub(crate) fn haploid_precise(columns: &[&str]) -> Result<String, RowError> {
    let sample = column(columns, 9)?;
    let mut info = sample.split(':');
    let genotype = info.next().unwrap_or_default().replace('|', "/");
    let genotype = match genotype.as_str() {
        "1/1" => "1",
        "0/0" => "0",
        _ => return Ok(String::new()),
    };
    Ok(with_sample(columns, genotype, info))
}

/// `update_haploid_sensitive_genotype`: any ALT genotype becomes `1`, the rest `0`; a
/// multi-allelic site gives an empty row.
pub(crate) fn haploid_sensitive(columns: &[&str]) -> Result<String, RowError> {
    let sample = column(columns, 9)?;
    let alt_base = column(columns, 4)?;
    if alt_base.contains(',') {
        return Ok(String::new());
    }
    let mut info = sample.split(':');
    let genotype = info.next().unwrap_or_default().replace('|', "/");
    let genotype = if matches!(genotype.as_str(), "0/1" | "1/0" | "1/1") { "1" } else { "0" };
    Ok(with_sample(columns, genotype, info))
}

/// `'\t'.join(columns) + '\n'` with column 9 replaced by `genotype` and the remaining fields.
fn with_sample<'a>(columns: &[&str], genotype: &'a str, info: impl Iterator<Item = &'a str>) -> String {
    let sample: Vec<&str> = std::iter::once(genotype).chain(info).collect();
    let alleles = if columns[4] == "." { 1 } else { 1 + columns[4].split(',').count() };
    let sample = haploid_sample(columns[8], &sample.join(":"), genotype, alleles, Quality::Keep);
    let mut fields: Vec<&str> = columns.to_vec();
    fields[9] = &sample;
    fields.join("\t") + "\n"
}

/// Where a haploid record's GQ comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Quality {
    /// The record's own: a variant call's, which is the same with and without the likelihoods
    /// (`PL` is written only for a gVCF, so a GQ taken from it would differ between the two).
    Keep,
    /// From the projected likelihoods: a gVCF reference block's.
    FromLikelihoods,
}

/// Conditional likelihoods for haploid alleles from the homozygous diploid states. Clair3's
/// model scores diploid genotypes. When its conditional likelihoods disagree with the selected
/// haploid genotype, publish the call without a likelihood (and, with
/// [`Quality::FromLikelihoods`], without a genotype-quality estimate).
pub(crate) fn haploid_sample(format: &str, sample: &str, genotype: &str, alleles: usize, quality: Quality) -> String {
    let mut fields: Vec<String> = sample.split(':').map(str::to_owned).collect();
    fields[0] = genotype.to_owned();
    let keys: Vec<_> = format.split(':').collect();
    if let Some(pl) = keys.iter().position(|&key| key == "PL").filter(|&index| index < fields.len()) {
        let projected = projected_pl(&fields[pl], alleles)
            .filter(|scores| genotype.parse::<usize>().ok().and_then(|index| scores.get(index)) == Some(&0));
        if quality == Quality::FromLikelihoods
            && let Some(gq) = keys.iter().position(|&key| key == "GQ").filter(|&index| index < fields.len())
        {
            fields[gq] = projected.as_ref().filter(|scores| scores.len() >= 2).map_or_else(
                || if genotype == "." { "0".to_owned() } else { ".".to_owned() },
                |scores| {
                    let allele = genotype.parse::<usize>().expect("projected genotype has an allele");
                    scores
                        .iter()
                        .enumerate()
                        .filter(|&(i, _)| i != allele)
                        .map(|(_, &score)| score)
                        .min()
                        .expect("more than one allele")
                        .min(99)
                        .to_string()
                },
            );
        }
        fields[pl] = projected
            .map_or_else(|| ".".to_owned(), |scores| scores.iter().map(u32::to_string).collect::<Vec<_>>().join(","));
    }
    fields.join(":")
}

fn projected_pl(value: &str, alleles: usize) -> Option<Vec<u32>> {
    let scores = value.split(',').map(str::parse::<u32>).collect::<Result<Vec<_>, _>>().ok()?;
    let mut projected = if scores.len() == alleles {
        scores
    } else if scores.len() == alleles * (alleles + 1) / 2 {
        (0..alleles).map(|i| scores[i * (i + 1) / 2 + i]).collect()
    } else {
        return None;
    };
    let min = *projected.iter().min()?;
    projected.iter_mut().for_each(|score| *score -= min);
    Some(projected)
}

/// `MarkLowQual(row, quality_score_for_pass, qual)`: a non-zero threshold marks rows with
/// `qual <= threshold` as `LowQual`.
pub(crate) fn mark_low_qual(row: String, threshold: i32, qual: f64) -> Result<String, RowError> {
    let low = threshold != 0 && qual <= f64::from(threshold);
    if row.is_empty() || !low {
        return Ok(row);
    }
    let mut fields: Vec<&str> = row.split('\t').collect();
    *fields.get_mut(6).ok_or(RowError::MissingColumn { index: 6 })? = "LowQual";
    Ok(fields.join("\t"))
}

impl RowFilter {
    /// `postprocess_row_with_params(args, row, decide)`: `None` drops the row. An empty row is
    /// kept (the original writes it as nothing, or, on the per-contig path, writes the row
    /// unchanged).
    pub(crate) fn apply(&self, row: &str) -> Result<Option<String>, RowError> {
        let columns = split_whitespace(strip(row), None);
        let (ref_base, alt_base) = (column(&columns, 3)?, column(&columns, 4)?);
        let qual = column(&columns, 5)?;
        let qual = parse_float(qual).ok_or_else(|| RowError::InvalidFloat { value: qual.to_owned() })?;
        let is_reference = alt_base == "." || ref_base == alt_base;
        let pos = column(&columns, 1)?;
        let pos = parse_int(pos).ok_or_else(|| RowError::InvalidInteger { value: pos.to_owned() })?;
        let mut row = match self.decider.decide(columns[0], pos) {
            Haploid::Precise => haploid_precise(&columns)?,
            Haploid::Sensitive => haploid_sensitive(&columns)?,
            Haploid::Off => row.to_owned(),
        };
        if !is_reference {
            row = mark_low_qual(row, self.qual, qual)?;
        }
        if is_reference && self.ref_calls == RefCalls::Hide {
            return Ok(None);
        }
        Ok(Some(row))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projects_haploid_likelihoods_without_using_heterozygous_states() {
        let blocks = Quality::FromLikelihoods;
        assert_eq!(haploid_sample("GT:GQ:PL", "1/1:99:66,35,0", "1", 2, blocks), "1:66:66,0");
        assert_eq!(haploid_sample("GT:GQ:PL", "0/1:15:0,4,19", "1", 2, blocks), "1:.:.");
        assert_eq!(haploid_sample("GT:GQ:PL", "./.:0:0,9,90", ".", 2, blocks), ".:0:.");
        // A variant call keeps its GQ, as it has it without the likelihoods.
        assert_eq!(haploid_sample("GT:GQ:PL", "1/1:99:66,35,0", "1", 2, Quality::Keep), "1:99:66,0");
        assert_eq!(haploid_sample("GT:GQ:PL", "0/1:15:0,4,19", "1", 2, Quality::Keep), "1:15:.");
        assert_eq!(haploid_sample("GT:GQ", "1/1:99", "1", 2, Quality::Keep), "1:99");
    }

    #[test]
    fn identifies_both_sides_of_a_par_for_reference_blocks() {
        let bed = b"chrX\t10\t20\nchrX\t20\t30\n";
        let decider =
            HaploidDecider::new(Haploid::Off, Gender::Male, ParRegions(BedTree::read(&bed[..], None, 0).unwrap()));
        assert_eq!(decider.decide("chrX", 10), Haploid::Precise);
        assert_eq!(decider.next_boundary("chrX", 10), Some(11));
        assert_eq!(decider.decide("chrX", 11), Haploid::Off);
        assert_eq!(decider.next_boundary("chrX", 11), Some(31));
        assert_eq!(decider.decide("chrX", 31), Haploid::Precise);
        assert_eq!(decider.next_boundary("chrX", 31), None);
    }

    const HET: &str = "chrX\t100\t.\tA\tC\t5.55\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:5:3:1,2:0.6667\n";
    const HOM: &str = "chrX\t100\t.\tA\tC\t25.5\tPASS\tP\tGT:GQ:DP:AD:AF\t1|1:5:3:1,2:0.6667\n";
    const REF: &str = "chrX\t100\t.\tA\t.\t1.5\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:5:3:3:0.6667\n";

    fn filter(haploid: Haploid, gender: Gender, ref_calls: RefCalls, qual: i32) -> RowFilter {
        RowFilter { ref_calls, qual, decider: HaploidDecider::new(haploid, gender, ParRegions::default()) }
    }

    #[test]
    fn ref_calls_and_low_qual() {
        let f = filter(Haploid::Off, Gender::Unknown, RefCalls::Hide, 10);
        assert_eq!(f.apply(REF), Ok(None));
        assert_eq!(f.apply(HET), Ok(Some(HET.replace("PASS", "LowQual"))));
        assert_eq!(f.apply(HOM), Ok(Some(HOM.to_owned())));
        // --qual 0 disables LowQual marking; RefCalls are never marked.
        let f = filter(Haploid::Off, Gender::Unknown, RefCalls::Show, 0);
        assert_eq!(f.apply(HET), Ok(Some(HET.to_owned())));
        let f = filter(Haploid::Off, Gender::Unknown, RefCalls::Show, 10);
        assert_eq!(f.apply(REF), Ok(Some(REF.to_owned())));
    }

    #[test]
    fn haploid_genotypes() {
        let precise = filter(Haploid::Precise, Gender::Unknown, RefCalls::Show, 2);
        assert_eq!(precise.apply(HET), Ok(Some(String::new())));
        assert_eq!(precise.apply(HOM), Ok(Some(HOM.replace("1|1:", "1:"))));
        assert_eq!(precise.apply(REF), Ok(Some(REF.replace("0/0:", "0:"))));
        let sensitive = filter(Haploid::Sensitive, Gender::Unknown, RefCalls::Show, 2);
        assert_eq!(sensitive.apply(HET), Ok(Some(HET.replace("0/1:", "1:"))));
        assert_eq!(sensitive.apply(&HET.replace("\tC\t", "\tC,G\t")), Ok(Some(String::new())));
        let short = "chrX\t100\t.\tA\tC\t5.55\tPASS\n";
        assert_eq!(precise.apply(short), Err(RowError::MissingColumn { index: 9 }));
        assert_eq!(precise.apply("chrX\t1x\t.\tA\tC\t5\n"), Err(RowError::InvalidInteger { value: "1x".into() }));
    }

    #[test]
    fn male_haploid_outside_par() {
        let par = BedTree::read("chrX\t200\t300\nchrX\t100\t200\n".as_bytes(), None, 0).expect("valid BED");
        let par = ParRegions(par);
        let decider = HaploidDecider::new(Haploid::Off, Gender::Male, par);
        // [100, 200) and [200, 300) are VCF positions 101 to 300.
        assert_eq!(decider.decide("chrX", 100), Haploid::Precise);
        assert_eq!(decider.decide("chrX", 101), Haploid::Off);
        assert_eq!(decider.decide("chrX", 200), Haploid::Off);
        assert_eq!(decider.decide("chrX", 300), Haploid::Off);
        assert_eq!(decider.decide("chrX", 301), Haploid::Precise);
        assert_eq!(decider.decide("chrY", 100), Haploid::Precise);
        assert_eq!(decider.decide("chr1", 100), Haploid::Off);
        let female = HaploidDecider::new(Haploid::Off, Gender::Female, ParRegions::default());
        assert_eq!(female.decide("chrX", 1), Haploid::Off);
        let global = HaploidDecider::new(Haploid::Sensitive, Gender::Male, ParRegions::default());
        assert_eq!(global.decide("chrX", 1), Haploid::Sensitive);
    }
}

//! `preprocess/SelectQual.py`: global QUAL cut-offs from the pileup VCF.

use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use crate::lines::{Lines, Newlines};
use crate::{Error, log_warning, py};

/// The file `SelectQual --phase` writes in its output folder (its `--qual_fn` default).
pub const PHASE_QUAL_FILE: &str = "phase_qual";

/// The file `SelectQual` writes in full-alignment mode (its `--qual_fn` default).
pub const FULL_QUAL_FILE: &str = "qual";

/// `--var_pct_phasing` default: the proportion of heterozygous SNPs, by QUAL, kept for phasing.
pub const DEFAULT_VAR_PCT_PHASING: f64 = 0.7;

/// `--var_pct_full` and `--ref_pct_full` default in `SelectQual.py`. (The original's pipeline
/// script always passes both: 0.7 and 0.1 for ONT, 0.3 and 0.3 otherwise.)
const DEFAULT_PCT_FULL: f64 = 0.3;

/// One record's fields as `row.rstrip().split()` gives them.
struct Record<'a> {
    line: u64,
    columns: Vec<&'a str>,
}

impl<'a> Record<'a> {
    fn column(&self, column: usize) -> Result<&'a str, Error> {
        self.columns.get(column).copied().ok_or(Error::MissingColumn { line: self.line, column })
    }

    /// `float(row[5])`, refusing NaN, which [`lowest_cutoff`] could not order like Python.
    fn qual(&self) -> Result<f64, Error> {
        let text = self.column(5)?;
        let qual =
            py::parse_float(text).ok_or_else(|| Error::InvalidNumber { line: self.line, value: text.to_owned() })?;
        if qual.is_nan() {
            return Err(Error::NanQuality { line: self.line });
        }
        Ok(qual)
    }
}

/// Calls `f` on every non-header record of a stream read like `sys.stdin`.
fn for_each_record(vcf: impl BufRead, mut f: impl FnMut(&Record<'_>) -> Result<(), Error>) -> Result<(), Error> {
    let mut lines = Lines::new(vcf, Newlines::Lf);
    while let Some((line, row)) = lines.next_line()? {
        // `row[0] == '#'` on the raw line, terminator included, so it is never empty.
        if row.starts_with('#') {
            continue;
        }
        f(&Record { line, columns: py::split(&row) })?;
    }
    Ok(())
}

/// `sorted(values)[:int(pct * len(values))][-1]`, or `None` when that slice is empty.
fn lowest_cutoff(mut values: Vec<f64>, pct: f64) -> Result<Option<f64>, Error> {
    // NaN is rejected while reading, so `partial_cmp` is a total preorder; the stable sort then
    // keeps -0.0 and 0.0 in input order, as Python's does.
    values.sort_by(|a, b| a.partial_cmp(b).expect("NaN qualities are rejected while reading"));
    let len = py::prefix_len(pct, values.len()).ok_or(Error::InvalidProportion { pct, count: values.len() })?;
    Ok(len.checked_sub(1).map(|last| values[last]))
}

fn write_file(path: &Path, contents: &str) -> Result<(), Error> {
    fs::write(path, contents).map_err(|source| Error::Io { path: path.to_owned(), source })
}

/// The result of `select_phase_qual_from_stdin`.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseQual {
    cutoff: f64,
    info: String,
    warnings: Vec<String>,
}

impl PhaseQual {
    /// The QUAL cut-off: heterozygous SNPs with a lower QUAL are not used for phasing.
    #[must_use]
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// The `phase_qual` file's contents: `str(cutoff)`, with no newline.
    #[must_use]
    pub fn contents(&self) -> String {
        py::float_repr(self.cutoff)
    }

    /// The `[INFO]` line the script prints to stdout, without its newline.
    #[must_use]
    pub fn info_message(&self) -> &str {
        &self.info
    }

    /// The warning lines (ANSI-coloured, without newlines) the script prints to stderr when no
    /// variant falls below the proportion.
    #[must_use]
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Writes [`contents`](Self::contents) to `dir/phase_qual` and returns that path.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be written.
    pub fn write_in(&self, dir: &Path) -> Result<PathBuf, Error> {
        let path = dir.join(PHASE_QUAL_FILE);
        write_file(&path, &self.contents())?;
        Ok(path)
    }
}

/// `SelectQual --phase --var_pct_phasing <var_pct_phasing>` on a decompressed pileup VCF.
///
/// Among single-base REF/ALT records with genotype exactly `0/1`, the cut-off is the highest
/// QUAL of the lowest `int((1 - var_pct_phasing) * n)`, or `0.0` if there are none.
///
/// # Errors
///
/// Fails on unreadable input, a record with fewer than 10 columns (6 for non-SNPs), an invalid
/// or NaN QUAL, or where the script raises: a non-finite proportion or an infinite cut-off.
pub fn select_phase_qual(vcf: impl BufRead, var_pct_phasing: f64) -> Result<PhaseQual, Error> {
    let low_qual_hete_var_pct = 1.0 - var_pct_phasing;
    let mut quals = Vec::new();
    for_each_record(vcf, |record| {
        let (ref_base, alt_base) = (record.column(3)?, record.column(4)?);
        if ref_base.chars().count() != 1 || alt_base.chars().count() != 1 {
            return Ok(());
        }
        // `qual, gt_info = row[5], row[9]`; QUAL is only converted for 0/1.
        record.column(5)?;
        let gt_info = record.column(9)?;
        if gt_info.split(':').next() == Some("0/1") {
            quals.push(record.qual()?);
        }
        Ok(())
    })?;

    let mut warnings = Vec::new();
    let cutoff = match lowest_cutoff(quals, low_qual_hete_var_pct)? {
        Some(cutoff) => cutoff,
        None => {
            // Reproduced: the count printed is that of the (empty) low-quality slice.
            warnings.push(log_warning(&format!(
                "[WARNING] Cannot find any 0/1 variant in pileup output using variant quality cut-off proportion: {}, total heterozygous variants: 0",
                py::float_repr(low_qual_hete_var_pct)
            )));
            warnings.push(log_warning("[WARNING] Set low variant quality score cut-off to 0.0"));
            0.0
        }
    };
    let rounded = py::round_to_int_repr(cutoff).ok_or(Error::NonFiniteCutoff { cutoff })?;
    let info = format!("[INFO] Select heterozygous pileup variants exceeding phasing quality cutoff {rounded}");
    Ok(PhaseQual { cutoff, info, warnings })
}

/// Whether a known-variants VCF (`--vcf_fn`) was given, which changes how `--ref_pct_full` is
/// applied.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum CandidateVcf {
    /// No file: the pipeline passes `EMPTY` or the path does not exist.
    #[default]
    Absent,
    /// The file exists.
    Present,
}

impl CandidateVcf {
    /// `file_path_from(args.vcf_fn) is not None`: present when `vcf_fn` names an existing
    /// regular file, relative paths resolved against the current directory. (So the pipeline's
    /// `EMPTY` counts as present if a file of that name is in the working directory.)
    #[must_use]
    pub fn from_arg(vcf_fn: &str) -> Self {
        if Path::new(vcf_fn).is_file() { Self::Present } else { Self::Absent }
    }
}

/// The arguments of `SelectQual` in full-alignment mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FullQualOptions {
    /// `--var_pct_full`: the proportion of variant calls, lowest QUAL first, sent to
    /// full-alignment calling.
    pub var_pct_full: f64,
    /// `--ref_pct_full`: the same for reference calls (`0/0`). Without a candidate VCF, zero
    /// means "use `var_pct_full`".
    pub ref_pct_full: f64,
    /// Whether `--vcf_fn` names an existing file.
    pub candidate_vcf: CandidateVcf,
}

impl Default for FullQualOptions {
    /// The script's argparse defaults: 0.3, 0.3 and no candidate VCF.
    fn default() -> Self {
        Self { var_pct_full: DEFAULT_PCT_FULL, ref_pct_full: DEFAULT_PCT_FULL, candidate_vcf: CandidateVcf::Absent }
    }
}

/// The result of `select_qual_from_stdin`.
#[derive(Debug, Clone, PartialEq)]
pub struct FullQual {
    var_cutoff: f64,
    ref_cutoff: f64,
    warnings: Vec<String>,
}

impl FullQual {
    /// The variant QUAL cut-off: variant calls at or below it go to full-alignment calling.
    #[must_use]
    pub fn var_cutoff(&self) -> f64 {
        self.var_cutoff
    }

    /// The reference-call QUAL cut-off.
    #[must_use]
    pub fn ref_cutoff(&self) -> f64 {
        self.ref_cutoff
    }

    /// The `qual` file's contents: `"{var} {ref}"` with Python float formatting, no newline.
    #[must_use]
    pub fn contents(&self) -> String {
        format!("{} {}", py::float_repr(self.var_cutoff), py::float_repr(self.ref_cutoff))
    }

    /// The two `[INFO]` lines the script prints to stdout, without newlines. The cut-offs are
    /// shown as `round(x, 0)`, a float (`23.0`).
    #[must_use]
    pub fn info_messages(&self) -> [String; 2] {
        [
            format!("[INFO] Set variants quality cutoff {}", py::round_to_float_repr(self.var_cutoff)),
            format!("[INFO] Set reference calls quality cutoff {}", py::round_to_float_repr(self.ref_cutoff)),
        ]
    }

    /// The warning lines (ANSI-coloured, without newlines) the script prints to stderr when a
    /// proportion selects nothing.
    #[must_use]
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Writes [`contents`](Self::contents) to `dir/qual` and returns that path.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be written.
    pub fn write_in(&self, dir: &Path) -> Result<PathBuf, Error> {
        let path = dir.join(FULL_QUAL_FILE);
        write_file(&path, &self.contents())?;
        Ok(path)
    }
}

/// `SelectQual` without `--phase` on a decompressed pileup VCF.
///
/// Records with genotype `0/0` are reference calls and every other record (any length, any
/// other genotype) a variant call; each cut-off is the highest QUAL among the lowest
/// `int(pct * n)` of its kind, or `0.0` if there are none.
///
/// # Errors
///
/// Fails on unreadable input, a record with fewer than 10 columns, an invalid or NaN QUAL, or
/// a non-finite proportion.
pub fn select_full_qual(vcf: impl BufRead, options: &FullQualOptions) -> Result<FullQual, Error> {
    let &FullQualOptions { var_pct_full, ref_pct_full: ref_pct_arg, candidate_vcf } = options;
    // `args.ref_pct_full if args.ref_pct_full else var_pct_full`: only zero is falsy (NaN is not).
    let ref_pct_full = if ref_pct_arg == 0.0 { var_pct_full } else { ref_pct_arg };
    let mut ref_quals = Vec::new();
    let mut variant_quals = Vec::new();
    for_each_record(vcf, |record| {
        record.column(5)?;
        let gt_info = record.column(9)?;
        let qual = record.qual()?;
        if gt_info.split(':').next() == Some("0/0") {
            ref_quals.push(qual);
        } else {
            variant_quals.push(qual);
        }
        Ok(())
    })?;

    let mut warnings = Vec::new();
    let variant_count = variant_quals.len();
    let var_cutoff = lowest_cutoff(variant_quals, var_pct_full)?.unwrap_or_else(|| {
        warnings.push(log_warning(&format!(
            "[WARNING] Cannot find any low-quality 0/1 or 1/1 variant in pileup output using variant quality cut-off proportion: {}, total variants: {variant_count}",
            py::float_repr(var_pct_full)
        )));
        warnings.push(log_warning("[WARNING] Set low variant quality score cut-off to 0.0"));
        0.0
    });

    // Reproduced: with a candidate VCF the raw `--ref_pct_full` is used, even zero, while the
    // warning still shows the fallback value.
    let used_ref_pct = match candidate_vcf {
        CandidateVcf::Absent => ref_pct_full,
        CandidateVcf::Present => ref_pct_arg,
    };
    let ref_count = ref_quals.len();
    let ref_cutoff = lowest_cutoff(ref_quals, used_ref_pct)?.unwrap_or_else(|| {
        warnings.push(log_warning(&format!(
            "[WARNING] Cannot find any low-quality 0/0 reference calls in pileup output using reference quality cut-off proportion: {}, total reference calls: {ref_count}",
            py::float_repr(ref_pct_full)
        )));
        warnings.push(log_warning("[WARNING] Set low reference quality score cut-off to 0.0"));
        0.0
    });
    Ok(FullQual { var_cutoff, ref_cutoff, warnings })
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

    fn vcf(records: &[(&str, &str, &str, &str)]) -> String {
        let mut text = HEADER.to_owned();
        for (i, (ref_base, alt, qual, gt)) in records.iter().enumerate() {
            text += &format!("chr1\t{}\t.\t{ref_base}\t{alt}\t{qual}\tPASS\tP\tGT:GQ\t{gt}:5\n", 100 + i);
        }
        text
    }

    #[test]
    fn phase_cutoff_is_the_last_of_the_low_slice() {
        let input = vcf(&[
            ("A", "C", "5.0", "0/1"),
            ("A", "C", "1.5", "0/1"),
            ("A", "C", "9.0", "0/1"),
            ("A", "C", "3.25", "0/1"),
            ("A", "C", "0.5", "1/1"),
            ("AT", "C", "0.1", "0/1"),
            ("A", "C", "0.2", "0|1"),
            ("A", "C,G", "0.3", "0/1"),
            ("A", "C", "4.0", "0/1"),
            ("A", "C", "7.0", "0/1"),
            ("A", "C", "8.0", "0/1"),
            ("A", "C", "2.0", "0/1"),
            ("A", "C", "6.0", "0/1"),
            ("A", "C", "10.0", "0/1"),
        ]);
        // Ten 0/1 SNPs; int(0.30000000000000004 * 10) == 3: 1.5, 2.0, 3.25.
        let phase = select_phase_qual(input.as_bytes(), 0.7).expect("valid input");
        assert_eq!(phase.contents(), "3.25");
        assert_eq!(
            phase.info_message(),
            "[INFO] Select heterozygous pileup variants exceeding phasing quality cutoff 3"
        );
        assert!(phase.warnings().is_empty());
    }

    #[test]
    fn phase_without_het_snps_warns_and_uses_zero() {
        let input = vcf(&[("A", "C", "5.0", "1/1"), ("A", "C", "4.0", "0/1")]);
        // int(0.30000000000000004 * 1) == 0.
        let phase = select_phase_qual(input.as_bytes(), 0.7).expect("valid input");
        assert_eq!(phase.contents(), "0.0");
        assert_eq!(
            phase.warnings(),
            [
                "\x1b[93m[WARNING] Cannot find any 0/1 variant in pileup output using variant quality cut-off proportion: 0.30000000000000004, total heterozygous variants: 0\x1b[0m",
                "\x1b[93m[WARNING] Set low variant quality score cut-off to 0.0\x1b[0m",
            ]
        );
        assert_eq!(
            phase.info_message(),
            "[INFO] Select heterozygous pileup variants exceeding phasing quality cutoff 0"
        );
    }

    #[test]
    fn phase_keeps_negative_zero_order() {
        // Equal under `<`, so the stable sort keeps input order and the last kept one is printed.
        let input = vcf(&[("A", "C", "0.0", "0/1"), ("A", "C", "-0.0", "0/1"), ("A", "C", "1", "0/1")]);
        assert_eq!(select_phase_qual(input.as_bytes(), 1.0 / 3.0).expect("valid").contents(), "-0.0");
    }

    #[test]
    fn phase_infinite_cutoff_fails_like_round() {
        let input = vcf(&[("A", "C", "inf", "0/1")]);
        assert!(matches!(select_phase_qual(input.as_bytes(), -1.0), Err(Error::NonFiniteCutoff { .. })));
    }

    #[test]
    fn full_mode_splits_reference_and_variant_calls() {
        let input = vcf(&[
            ("A", ".", "9.5", "0/0"),
            ("A", ".", "3.5", "0/0"),
            ("A", ".", "1.25", "0/0"),
            ("AT", "A", "2.5", "0/1"),
            ("A", "C", "22.5", "1/1"),
            ("A", "C,G", "4.0", "1/2"),
            ("A", "C", "1e-05", "./."),
        ]);
        let options = FullQualOptions { var_pct_full: 0.7, ref_pct_full: 0.1, candidate_vcf: CandidateVcf::Absent };
        let full = select_full_qual(input.as_bytes(), &options).expect("valid input");
        // Variants: int(0.7 * 4) == 2 -> 2.5; references: int(0.1 * 3) == 0 -> warning, 0.0.
        assert_eq!(full.contents(), "2.5 0.0");
        assert_eq!(
            full.info_messages(),
            ["[INFO] Set variants quality cutoff 2.0", "[INFO] Set reference calls quality cutoff 0.0"]
        );
        assert_eq!(full.warnings().len(), 2);
        assert!(full.warnings()[0].contains("reference quality cut-off proportion: 0.1, total reference calls: 3"));
    }

    #[test]
    fn full_mode_zero_ref_pct_depends_on_candidate_vcf() {
        let input = vcf(&[("A", ".", "9.5", "0/0"), ("A", ".", "3.5", "0/0"), ("A", "C", "8.0", "0/1")]);
        let absent = FullQualOptions { var_pct_full: 1.0, ref_pct_full: 0.0, candidate_vcf: CandidateVcf::Absent };
        assert_eq!(select_full_qual(input.as_bytes(), &absent).expect("valid").contents(), "8.0 9.5");
        let present = FullQualOptions { candidate_vcf: CandidateVcf::Present, ..absent };
        let full = select_full_qual(input.as_bytes(), &present).expect("valid");
        assert_eq!(full.contents(), "8.0 0.0");
        assert!(full.warnings()[0].contains("proportion: 1.0, total reference calls: 2"));
    }

    #[test]
    fn short_records_fail_like_index_errors() {
        let input = format!("{HEADER}chr1\t100\t.\tA\tC\t5.0\n");
        assert!(matches!(
            select_full_qual(input.as_bytes(), &FullQualOptions::default()),
            Err(Error::MissingColumn { line: 3, column: 9 })
        ));
        // Indels are skipped before column 9 is read in phase mode.
        let input = format!("{HEADER}chr1\t100\t.\tAT\tC\t5.0\n");
        assert!(select_phase_qual(input.as_bytes(), 0.7).is_ok());
        let input = format!("{HEADER}\n");
        assert!(matches!(select_phase_qual(input.as_bytes(), 0.7), Err(Error::MissingColumn { line: 3, column: 3 })));
    }

    #[test]
    fn nan_quality_is_refused() {
        let input = vcf(&[("A", "C", "nan", "0/1")]);
        assert!(matches!(select_phase_qual(input.as_bytes(), 0.7), Err(Error::NanQuality { line: 3 })));
        assert!(matches!(select_phase_qual(input.as_bytes(), f64::NAN), Err(Error::NanQuality { .. })));
        let input = vcf(&[("A", "C", "1.0", "0/1")]);
        assert!(matches!(select_phase_qual(input.as_bytes(), f64::NAN), Err(Error::InvalidProportion { .. })));
    }
}

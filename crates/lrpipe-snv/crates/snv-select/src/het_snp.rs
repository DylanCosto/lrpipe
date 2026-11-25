//! `preprocess/SelectHetSnp.py`, `FiterHeteSnpPhasing`: one contig's heterozygous SNPs for
//! WhatsHap.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::lines::{Lines, Newlines, open_vcf};
use crate::{Error, PHASE_QUAL_FILE, py};

/// `--var_pct_full` default in `SelectHetSnp.py`, used only when there is no `phase_qual` file.
const DEFAULT_VAR_PCT_FULL: f64 = 0.3;

/// Whether `--phase` was passed to `SelectHetSnp`. It is a `store_false` flag: passing it
/// turns phasing mode *off*.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum PhaseFlag {
    /// Not passed (as in the pipeline): phasing mode.
    #[default]
    Absent,
    /// Passed.
    Given,
}

/// The function `SelectHetSnp.main()` dispatches to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum HetSnpMode {
    /// `FiterHeteSnpPhasing`, the only one the calling pipeline uses; ported.
    Phasing,
    /// `FiterHeteSnp_FP` (`--phase --calling 1`), a testing function; not ported.
    FalsePositive,
    /// `FiterHeteSnp` (`--phase` otherwise), for training and alt-file input; not ported.
    Training,
}

impl HetSnpMode {
    /// The mode `main()` picks from `--phase` and `--calling`.
    #[must_use]
    pub fn from_args(phase: PhaseFlag, calling: i64) -> Self {
        match (phase, calling) {
            (PhaseFlag::Absent, _) => Self::Phasing,
            (PhaseFlag::Given, 1) => Self::FalsePositive,
            (PhaseFlag::Given, _) => Self::Training,
        }
    }

    /// Succeeds for the ported mode.
    ///
    /// # Errors
    ///
    /// [`Error::UnsupportedMode`] for the others.
    pub fn ensure_ported(self) -> Result<(), Error> {
        match self {
            Self::Phasing => Ok(()),
            Self::FalsePositive => Err(Error::UnsupportedMode { function: "FiterHeteSnp_FP" }),
            Self::Training => Err(Error::UnsupportedMode { function: "FiterHeteSnp" }),
        }
    }
}

/// Which heterozygous SNPs are dropped.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum HetSnpFilter {
    /// Those with QUAL below the cut-off from `phase_qual` (NaN QUALs are kept, as `<` is
    /// false for them).
    BelowCutoff(f64),
    /// Without a `phase_qual` file: the lowest `int(pct * n)` by QUAL, ties broken by first
    /// appearance of the position.
    LowestProportion(f64),
}

/// The arguments of `SelectHetSnp` in phasing mode.
#[derive(Debug, Clone, PartialEq)]
pub struct HetSnpOptions {
    /// `--ctgName`. `None` or empty selects every contig; the output is then named after it
    /// anyway (`None.vcf` or `.vcf`).
    pub contig: Option<String>,
    /// `--var_pct_full`: the proportion dropped when there is no cut-off file.
    pub var_pct_full: f64,
    /// `--qual_fn`: the cut-off file's name in the split folder.
    pub qual_fn: String,
}

impl Default for HetSnpOptions {
    /// The script's argparse defaults: every contig, 0.3 and `phase_qual`.
    fn default() -> Self {
        Self { contig: None, var_pct_full: DEFAULT_VAR_PCT_FULL, qual_fn: PHASE_QUAL_FILE.to_owned() }
    }
}

/// The selected SNPs and the VCF header, ready to write.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HetSnps {
    contig: Option<String>,
    header: String,
    /// Records sorted by position.
    records: Vec<String>,
}

impl HetSnps {
    /// The number of selected SNPs.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether no SNP was selected.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// `str(contig_name)`: the contig, or `None` when none was given.
    fn contig_label(&self) -> &str {
        self.contig.as_deref().unwrap_or("None")
    }

    /// The `[INFO]` line the script prints to stdout, without its newline.
    #[must_use]
    pub fn info_message(&self) -> String {
        format!("[INFO] Total heterozygous SNP positions selected: {}: {}", self.contig_label(), self.len())
    }

    /// The output file's name: `{contig}.vcf`.
    #[must_use]
    pub fn file_name(&self) -> String {
        format!("{}.vcf", self.contig_label())
    }

    /// Writes the VCF: every header line of the input (right-stripped), then the records by
    /// position, each ending in `\n`.
    ///
    /// # Errors
    ///
    /// Fails if writing fails.
    pub fn write_vcf(&self, mut out: impl Write) -> io::Result<()> {
        out.write_all(self.header.as_bytes())?;
        for record in &self.records {
            out.write_all(record.as_bytes())?;
            out.write_all(b"\n")?;
        }
        out.flush()
    }

    /// Writes the VCF to `split_folder/{contig}.vcf` and returns that path.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be written.
    pub fn write_in(&self, split_folder: &Path) -> Result<PathBuf, Error> {
        let path = split_folder.join(self.file_name());
        let io_err = |source| Error::Io { path: path.clone(), source };
        let file = File::create(&path).map_err(io_err)?;
        self.write_vcf(BufWriter::new(file)).map_err(io_err)?;
        Ok(path)
    }
}

/// Reads the cut-off `SelectQual --phase` left in `split_folder/qual_fn`, or `None` if there is
/// no such file.
///
/// # Errors
///
/// Fails if the file exists but cannot be read or does not hold a number.
pub fn read_phase_qual_cutoff(split_folder: &Path, qual_fn: &str) -> Result<Option<f64>, Error> {
    let path = split_folder.join(qual_fn);
    if !path.exists() {
        return Ok(None);
    }
    let text = fs::read_to_string(&path).map_err(|source| Error::Io { path: path.clone(), source })?;
    // float(open(f).read().rstrip())
    match py::parse_float(py::rstrip(&text)) {
        Some(cutoff) => Ok(Some(cutoff)),
        None => Err(Error::InvalidCutoffFile { path, value: text }),
    }
}

/// A heterozygous SNP kept so far, in `dict` insertion order.
struct Snp {
    pos: i64,
    row: String,
    qual: f64,
    /// The line it was last set from, for error messages.
    line: u64,
}

fn column<'a>(columns: &[&'a str], line: u64, column: usize) -> Result<&'a str, Error> {
    columns.get(column).copied().ok_or(Error::MissingColumn { line, column })
}

/// `row` with the PL field removed from FORMAT (the second-to-last column) and from the last
/// column, re-joined with tabs, if FORMAT contains `PL` anywhere; otherwise `row` unchanged.
fn without_pl(row: &str, columns: &[&str], line: u64) -> Result<String, Error> {
    let [.., format, sample] = columns else { unreachable!("a het SNP record has at least 10 columns") };
    if !format.contains("PL") {
        return Ok(row.to_owned());
    }
    let mut format_list: Vec<&str> = format.split(':').collect();
    let index = format_list
        .iter()
        .position(|&f| f == "PL")
        .ok_or_else(|| Error::PlNotInFormat { line, format: (*format).to_owned() })?;
    format_list.remove(index);
    let mut sample_list: Vec<&str> = sample.split(':').collect();
    if index >= sample_list.len() {
        return Err(Error::PlMissingInSample { line, index });
    }
    sample_list.remove(index);
    let (new_format, new_sample) = (format_list.join(":"), sample_list.join(":"));
    let mut out: Vec<&str> = columns[..columns.len() - 2].to_vec();
    out.push(&new_format);
    out.push(&new_sample);
    Ok(out.join("\t"))
}

/// `FiterHeteSnpPhasing` on a decompressed VCF: single-base REF/ALT records of `contig` with
/// genotype `0/1` or `1/0` (`|` read as `/`), minus those `filter` drops, with PL removed.
///
/// A later record at the same position replaces an earlier one. Header lines (`#...`) are
/// collected from anywhere in the file.
///
/// # Errors
///
/// Fails on unreadable input, an empty line, a record of `contig` with fewer than 10 columns,
/// an invalid POS or QUAL, a FORMAT with `PL` as a substring only, a sample column too short
/// for the PL index, or where the script raises on a non-finite proportion. With
/// [`HetSnpFilter::LowestProportion`], a NaN QUAL is refused, as Python's sort order for it
/// cannot be reproduced.
pub fn select_het_snps(vcf: impl BufRead, contig: Option<&str>, filter: HetSnpFilter) -> Result<HetSnps, Error> {
    // `if contig_name and contig_name != ctg_name`: an empty name filters nothing.
    let wanted = contig.filter(|c| !c.is_empty());
    let mut header = String::new();
    let mut snps: Vec<Snp> = Vec::new();
    let mut index: HashMap<i64, usize> = HashMap::new();
    let mut lines = Lines::new(vcf, Newlines::Universal);
    while let Some((line, raw)) = lines.next_line()? {
        let row = py::rstrip(&raw);
        if row.is_empty() {
            return Err(Error::EmptyLine { line });
        }
        if row.starts_with('#') {
            header.push_str(row);
            header.push('\n');
            continue;
        }
        let columns = py::split(row);
        if wanted.is_some_and(|c| c != columns[0]) {
            continue;
        }
        let pos_text = column(&columns, line, 1)?;
        let pos = py::parse_int(pos_text).ok_or_else(|| Error::InvalidNumber { line, value: pos_text.to_owned() })?;
        let (ref_base, alt_base) = (column(&columns, line, 3)?, column(&columns, line, 4)?);
        let genotype = column(&columns, line, 9)?.split(':').next().unwrap_or_default().replace('|', "/");
        if ref_base.chars().count() != 1 || alt_base.chars().count() != 1 || !(genotype == "0/1" || genotype == "1/0") {
            continue;
        }
        let row = without_pl(row, &columns, line)?;
        let qual_text = columns[5];
        let qual =
            py::parse_float(qual_text).ok_or_else(|| Error::InvalidNumber { line, value: qual_text.to_owned() })?;
        // `variant_dict[pos] = row; qual_set[pos] = qual`: a repeated key keeps its first slot.
        match index.get(&pos) {
            Some(&i) => snps[i] = Snp { pos, row, qual, line },
            None => {
                index.insert(pos, snps.len());
                snps.push(Snp { pos, row, qual, line });
            }
        }
    }

    let mut keep = vec![true; snps.len()];
    match filter {
        HetSnpFilter::BelowCutoff(cutoff) => {
            for (k, snp) in keep.iter_mut().zip(&snps) {
                // `v < cutoff`, not `v >= cutoff`: a NaN QUAL is kept.
                let dropped = snp.qual < cutoff;
                *k = !dropped;
            }
        }
        HetSnpFilter::LowestProportion(pct) => {
            if let Some(snp) = snps.iter().find(|s| s.qual.is_nan()) {
                return Err(Error::NanQuality { line: snp.line });
            }
            // sorted(qual_set.items(), key=lambda x: x[1]) is stable: ties stay in dict order.
            let mut order: Vec<usize> = (0..snps.len()).collect();
            order.sort_by(|&a, &b| snps[a].qual.partial_cmp(&snps[b].qual).expect("NaN qualities were refused"));
            let count = py::prefix_len(pct, order.len()).ok_or(Error::InvalidProportion { pct, count: order.len() })?;
            for &i in &order[..count] {
                keep[i] = false;
            }
        }
    }
    let mut kept: Vec<Snp> = snps.into_iter().zip(keep).filter_map(|(snp, k)| k.then_some(snp)).collect();
    kept.sort_unstable_by_key(|snp| snp.pos);
    Ok(HetSnps { contig: contig.map(str::to_owned), header, records: kept.into_iter().map(|snp| snp.row).collect() })
}

/// `SelectHetSnp --vcf_fn <vcf_fn> --split_folder <split_folder> --ctgName <contig>`: reads
/// the cut-off from `split_folder` (falling back to [`HetSnpFilter::LowestProportion`] when the
/// file is absent), selects the SNPs and writes `split_folder/{contig}.vcf`.
///
/// # Errors
///
/// Fails as [`read_phase_qual_cutoff`], [`open_vcf`] and [`select_het_snps`] do, or if the
/// output cannot be written. (The script, reading through a `pigz` pipe, would instead treat
/// an unreadable VCF as empty.)
pub fn select_het_snps_for_phasing(
    vcf_fn: &Path,
    split_folder: &Path,
    options: &HetSnpOptions,
) -> Result<HetSnps, Error> {
    select_het_snps_for_phasing_from(open_vcf(vcf_fn)?, split_folder, options)
}

/// As [`select_het_snps_for_phasing`], reading the pileup records from `vcf` (all of them, or
/// only the contig's).
///
/// # Errors
///
/// As [`select_het_snps_for_phasing`].
pub fn select_het_snps_for_phasing_from(
    vcf: impl BufRead,
    split_folder: &Path,
    options: &HetSnpOptions,
) -> Result<HetSnps, Error> {
    let filter = match read_phase_qual_cutoff(split_folder, &options.qual_fn)? {
        Some(cutoff) => HetSnpFilter::BelowCutoff(cutoff),
        None => HetSnpFilter::LowestProportion(options.var_pct_full),
    };
    let selection = select_het_snps(vcf, options.contig.as_deref(), filter)?;
    selection.write_in(split_folder)?;
    Ok(selection)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf_string(selection: &HetSnps) -> String {
        let mut out = Vec::new();
        selection.write_vcf(&mut out).expect("writing to a Vec cannot fail");
        String::from_utf8(out).expect("UTF-8")
    }

    #[test]
    fn selects_sorted_het_snps_above_cutoff() {
        let input = "##h  \n#CHROM\tPOS\n\
            chr1\t300\t.\tA\tC\t20.0\tPASS\tP\tGT:GQ\t0/1:20\n\
            chr1\t100\t.\tA\tC\t10.0\tPASS\tP\tGT:GQ\t1|0:10\n\
            chr2\t50\t.\tA\tC\t30.0\tPASS\tP\tGT:GQ\t0/1:30\n\
            chr1\t200\t.\tA\tC\t5.0\tPASS\tP\tGT:GQ\t0/1:5\n\
            chr1\t250\t.\tAT\tC\t50.0\tPASS\tP\tGT:GQ\t0/1:50\n\
            chr1\t260\t.\tA\tC\t50.0\tPASS\tP\tGT:GQ\t1/1:50\n";
        let selection =
            select_het_snps(input.as_bytes(), Some("chr1"), HetSnpFilter::BelowCutoff(10.0)).expect("valid");
        assert_eq!(selection.info_message(), "[INFO] Total heterozygous SNP positions selected: chr1: 2");
        assert_eq!(selection.file_name(), "chr1.vcf");
        assert_eq!(
            vcf_string(&selection),
            "##h\n#CHROM\tPOS\n\
             chr1\t100\t.\tA\tC\t10.0\tPASS\tP\tGT:GQ\t1|0:10\n\
             chr1\t300\t.\tA\tC\t20.0\tPASS\tP\tGT:GQ\t0/1:20\n"
        );
    }

    #[test]
    fn removes_pl_and_rejoins_with_tabs() {
        let input = "chr1  100 . A C 10.0 PASS P GT:PL:GQ 0/1:0,10,20:10  \n\
                     chr1\t200\t.\tA\tC\t10.0\tPASS\tP\tGT:GQ\t0/1:10\n";
        let selection = select_het_snps(input.as_bytes(), Some("chr1"), HetSnpFilter::BelowCutoff(0.0)).expect("valid");
        assert_eq!(
            vcf_string(&selection),
            "chr1\t100\t.\tA\tC\t10.0\tPASS\tP\tGT:GQ\t0/1:10\nchr1\t200\t.\tA\tC\t10.0\tPASS\tP\tGT:GQ\t0/1:10\n"
        );
        let bad = "chr1\t100\t.\tA\tC\t10.0\tPASS\tP\tGT:XPL\t0/1:10\n";
        assert!(matches!(
            select_het_snps(bad.as_bytes(), None, HetSnpFilter::BelowCutoff(0.0)),
            Err(Error::PlNotInFormat { line: 1, .. })
        ));
    }

    #[test]
    fn lowest_proportion_breaks_ties_by_first_appearance() {
        let input = "c\t30\t.\tA\tC\t5.0\t.\t.\tGT\t0/1\n\
                     c\t10\t.\tA\tC\t5.0\t.\t.\tGT\t0/1\n\
                     c\t20\t.\tA\tC\t5.0\t.\t.\tGT\t0/1\n\
                     c\t30\t.\tA\tC\t9.0\t.\t.\tGT\t0/1\n\
                     c\t40\t.\tA\tC\t1.0\t.\t.\tGT\t0/1\n";
        // Qualities in dict order: 30 -> 9.0, 10 -> 5.0, 20 -> 5.0, 40 -> 1.0; int(0.5 * 4) == 2
        // drops 40 and then 10, the first 5.0 in dict order.
        let selection = select_het_snps(input.as_bytes(), None, HetSnpFilter::LowestProportion(0.5)).expect("valid");
        assert_eq!(selection.info_message(), "[INFO] Total heterozygous SNP positions selected: None: 2");
        assert_eq!(selection.file_name(), "None.vcf");
        assert_eq!(vcf_string(&selection), "c\t20\t.\tA\tC\t5.0\t.\t.\tGT\t0/1\nc\t30\t.\tA\tC\t9.0\t.\t.\tGT\t0/1\n");
    }

    #[test]
    fn empty_line_fails_like_python() {
        let input = "#h\n   \nc\t1\t.\tA\tC\t5.0\t.\t.\tGT\t0/1\n";
        assert!(matches!(
            select_het_snps(input.as_bytes(), None, HetSnpFilter::BelowCutoff(0.0)),
            Err(Error::EmptyLine { line: 2 })
        ));
    }

    #[test]
    fn mode_dispatch() {
        assert_eq!(HetSnpMode::from_args(PhaseFlag::Absent, 1), HetSnpMode::Phasing);
        assert!(HetSnpMode::from_args(PhaseFlag::Absent, 0).ensure_ported().is_ok());
        assert!(matches!(
            HetSnpMode::from_args(PhaseFlag::Given, 1).ensure_ported(),
            Err(Error::UnsupportedMode { function: "FiterHeteSnp_FP" })
        ));
        assert!(matches!(
            HetSnpMode::from_args(PhaseFlag::Given, 0).ensure_ported(),
            Err(Error::UnsupportedMode { function: "FiterHeteSnp" })
        ));
    }
}

//! `preprocess/SelectCandidates.py`: the low-quality pileup calls re-called with the
//! full-alignment model, written as candidate BED chunks.
//!
//! Only the path the calling pipeline takes is ported: `--phasing_info_in_bam` not passed
//! (the script then reads phasing from the haplotagged BAM) and `--call_low_seq_entropy`
//! off, so the sequence-entropy code never runs.

use std::fs;
use std::io::BufRead;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use crate::lines::{Lines, Newlines, open_vcf};
use crate::{Error, FULL_QUAL_FILE, FullQual, log_warning, py};

/// `param.flankingBaseNum` (`shared/param_f.py`): bases on each side of a candidate in its
/// full-alignment window.
const FLANKING_BASE_NUM: i128 = 16;

/// `param.no_of_positions`: the full-alignment window width.
const NO_OF_POSITIONS: i128 = 2 * FLANKING_BASE_NUM + 1;

/// `param.split_region_size`: the block size Illumina candidates are rounded to.
const SPLIT_REGION_SIZE: i128 = 1000;

/// `--var_pct_full` and `--ref_pct_full` default in `SelectCandidates.py`, used only without a
/// `qual` file.
const DEFAULT_PCT_FULL: f64 = 0.3;

/// `--split_bed_size` default: the most candidates in one BED chunk.
pub const DEFAULT_SPLIT_BED_SIZE: NonZeroUsize = NonZeroUsize::new(10_000).expect("10000 is not zero");

/// The prefix of the file listing one contig's BED chunks: `FULL_ALN_FILE_{contig}`.
pub const FULL_ALN_FILE_PREFIX: &str = "FULL_ALN_FILE_";

/// The QUAL cut-offs `SelectQual` writes to `qual`: calls strictly below them are candidates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QualCutoffs {
    /// The cut-off for variant calls.
    pub var_cutoff: f64,
    /// The cut-off for reference calls.
    pub ref_cutoff: f64,
}

impl From<&FullQual> for QualCutoffs {
    /// The cut-offs as the `qual` file holds them (its float repr reads back exactly).
    fn from(full: &FullQual) -> Self {
        Self { var_cutoff: full.var_cutoff(), ref_cutoff: full.ref_cutoff() }
    }
}

/// Reads the cut-offs `SelectQual` left in `split_folder/qual_fn`, or `None` if there is no
/// such file.
///
/// # Errors
///
/// Fails if the file exists but cannot be read, or does not start with two numbers separated
/// by a single space (where the script's `split(' ')`, `float()` or indexing raises).
pub fn read_qual_cutoffs(split_folder: &Path, qual_fn: &str) -> Result<Option<QualCutoffs>, Error> {
    let path = split_folder.join(qual_fn);
    if !path.exists() {
        return Ok(None);
    }
    let text = fs::read_to_string(&path).map_err(|source| Error::Io { path: path.clone(), source })?;
    // line = f.read().rstrip().split(' '); float(line[0]), float(line[1])
    let mut fields = py::rstrip(&text).split(' ');
    let mut next_float = || fields.next().and_then(py::parse_float);
    match (next_float(), next_float()) {
        (Some(var_cutoff), Some(ref_cutoff)) => Ok(Some(QualCutoffs { var_cutoff, ref_cutoff })),
        _ => Err(Error::InvalidCutoffFile { path, value: text }),
    }
}

/// Which pileup calls become full-alignment candidates.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum CandidateFilter {
    /// Calls with QUAL strictly below the cut-off of their kind (a NaN QUAL is never below).
    BelowCutoffs(QualCutoffs),
    /// Without a `qual` file: the lowest `int(pct * n)` calls of each kind by QUAL, ties in
    /// file order.
    LowestProportions {
        /// The proportion of variant calls.
        var_pct_full: f64,
        /// The proportion of reference calls.
        ref_pct_full: f64,
    },
}

/// The window written for each candidate, chosen by `--platform`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum CandidateWindow {
    /// Any platform but `ilmn`: the 33-base window centred on the candidate.
    #[default]
    Flanking,
    /// `ilmn`: the 1000-base block holding the candidate, widened by 33 bases on each side.
    IlluminaBlock,
}

impl CandidateWindow {
    /// The window for a `--platform` value: `platform == 'ilmn'`, compared exactly.
    #[must_use]
    pub fn from_platform(platform: &str) -> Self {
        if platform == "ilmn" { Self::IlluminaBlock } else { Self::Flanking }
    }

    /// The BED `(start, end)` the script writes for a candidate at 1-based `pos`.
    fn bed_region(self, pos: i64) -> (i128, i128) {
        let pos = i128::from(pos);
        // The script builds a 1-based (start, end) pair and writes `max(start - 1, 0)` and
        // `end - 1`.
        let (start, end) = match self {
            Self::Flanking => (pos - FLANKING_BASE_NUM, pos + FLANKING_BASE_NUM + 2),
            Self::IlluminaBlock => {
                // `item // region_size * region_size`: floor division.
                let block = pos.div_euclid(SPLIT_REGION_SIZE) * SPLIT_REGION_SIZE;
                ((block - NO_OF_POSITIONS).max(1), block + SPLIT_REGION_SIZE + NO_OF_POSITIONS)
            }
        };
        ((start - 1).max(0), end - 1)
    }
}

/// Whether `--call_low_seq_entropy` is on.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum SeqEntropyCalling {
    /// Off (the default and the pipeline's choice).
    #[default]
    Off,
    /// On: also re-call the candidates with the lowest reference sequence entropy. Not ported.
    On,
}

/// Whether `--phasing_info_in_bam` was passed. It is a `store_false` flag: passing it turns
/// off the script's `phasing_info_in_bam` mode.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum PhasingInfoFlag {
    /// Not passed (as in the pipeline): candidates are written as plain windows.
    #[default]
    Absent,
    /// Passed: windows named by their first and last position, plus nearby phased het SNPs.
    /// Not ported.
    Given,
}

/// The arguments of `SelectCandidates`. (`--ref_fn`, `--samtools`, `--seq_entropy_pro` and
/// `--chr_prefix` do not affect the ported path.)
#[derive(Debug, Clone, PartialEq)]
pub struct CandidateOptions {
    /// `--ctgName`. `None` or empty selects every contig; the outputs are then named after it
    /// anyway (`None.0_1` or `.0_1`).
    pub contig: Option<String>,
    /// `--var_pct_full`: the proportion of variant calls used when there is no `qual` file.
    pub var_pct_full: f64,
    /// `--ref_pct_full`: the same for reference calls.
    pub ref_pct_full: f64,
    /// `--platform`, through [`CandidateWindow::from_platform`].
    pub window: CandidateWindow,
    /// `--split_bed_size`. (The script raises on zero and writes no chunk for a negative size.)
    pub split_bed_size: NonZeroUsize,
    /// `--qual_fn`: the cut-off file's name in the split folder.
    pub qual_fn: String,
    /// `--call_low_seq_entropy`; only [`SeqEntropyCalling::Off`] is ported.
    pub call_low_seq_entropy: SeqEntropyCalling,
    /// `--phasing_info_in_bam`; only [`PhasingInfoFlag::Absent`] is ported.
    pub phasing_info_in_bam: PhasingInfoFlag,
    /// `--phased_vcf_fn`; not ported, refused when given.
    pub phased_vcf_fn: Option<PathBuf>,
    /// `--all_alt_fn`; not ported, refused when given.
    pub all_alt_fn: Option<PathBuf>,
}

impl Default for CandidateOptions {
    /// The script's argparse defaults: every contig, 0.3, 0.3, `ont`, 10000 and `qual`.
    fn default() -> Self {
        Self {
            contig: None,
            var_pct_full: DEFAULT_PCT_FULL,
            ref_pct_full: DEFAULT_PCT_FULL,
            window: CandidateWindow::Flanking,
            split_bed_size: DEFAULT_SPLIT_BED_SIZE,
            qual_fn: FULL_QUAL_FILE.to_owned(),
            call_low_seq_entropy: SeqEntropyCalling::Off,
            phasing_info_in_bam: PhasingInfoFlag::Absent,
            phased_vcf_fn: None,
            all_alt_fn: None,
        }
    }
}

impl CandidateOptions {
    /// Succeeds when every option is on the ported path.
    ///
    /// # Errors
    ///
    /// [`Error::UnsupportedOption`] naming the first option that is not.
    pub fn ensure_ported(&self) -> Result<(), Error> {
        let unsupported = |option| Err(Error::UnsupportedOption { option });
        if self.call_low_seq_entropy == SeqEntropyCalling::On {
            return unsupported("--call_low_seq_entropy True");
        }
        if self.phasing_info_in_bam == PhasingInfoFlag::Given {
            return unsupported("--phasing_info_in_bam");
        }
        if self.phased_vcf_fn.is_some() {
            return unsupported("--phased_vcf_fn");
        }
        if self.all_alt_fn.is_some() {
            return unsupported("--all_alt_fn");
        }
        Ok(())
    }
}

/// One pileup call: its position, QUAL and line (for error messages).
struct Call {
    pos: i64,
    qual: f64,
    line: u64,
}

/// The positions of the calls with QUAL strictly below `cutoff`, in file order.
fn below_cutoff(calls: Vec<Call>, cutoff: f64) -> Vec<i64> {
    calls.into_iter().filter(|call| call.qual < cutoff).map(|call| call.pos).collect()
}

/// The positions of the lowest `int(pct * len(calls))` calls by QUAL.
fn lowest_proportion(mut calls: Vec<Call>, pct: f64) -> Result<Vec<i64>, Error> {
    if let Some(call) = calls.iter().find(|call| call.qual.is_nan()) {
        return Err(Error::NanQuality { line: call.line });
    }
    let count = py::prefix_len(pct, calls.len()).ok_or(Error::InvalidProportion { pct, count: calls.len() })?;
    // `sorted(..., key=lambda x: x[1])` is stable.
    calls.sort_by(|a, b| a.qual.partial_cmp(&b.qual).expect("NaN qualities were refused"));
    Ok(calls[..count].iter().map(|call| call.pos).collect())
}

fn column<'a>(columns: &[&'a str], line: u64, column: usize) -> Result<&'a str, Error> {
    columns.get(column).copied().ok_or(Error::MissingColumn { line, column })
}

/// The selected candidates of one contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Candidates {
    contig: Option<String>,
    low_qual_reference_calls: usize,
    low_qual_variants: usize,
    /// Sorted, without duplicates.
    positions: Vec<i64>,
}

impl Candidates {
    /// The number of candidate positions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Whether there is no candidate.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// The candidate positions (1-based), sorted and without duplicates.
    #[must_use]
    pub fn positions(&self) -> &[i64] {
        &self.positions
    }

    /// The number of reference calls selected, counting repeated positions.
    #[must_use]
    pub fn low_qual_reference_calls(&self) -> usize {
        self.low_qual_reference_calls
    }

    /// The number of variant calls selected, counting repeated positions.
    #[must_use]
    pub fn low_qual_variants(&self) -> usize {
        self.low_qual_variants
    }

    /// `str(contig_name)`: the contig, or `None` when none was given.
    fn contig_label(&self) -> &str {
        self.contig.as_deref().unwrap_or("None")
    }

    /// The lines (without newlines) the script prints to stderr before writing: two `[INFO]`
    /// counts, then an ANSI-coloured warning if there is no candidate.
    #[must_use]
    pub fn messages(&self) -> Vec<String> {
        let contig = self.contig_label();
        let mut messages = vec![
            format!(
                "[INFO] Low quality reference calls to be processed in {contig}: {}",
                self.low_qual_reference_calls
            ),
            format!("[INFO] Low quality variants to be processed in {contig}: {}", self.low_qual_variants),
        ];
        if self.is_empty() {
            messages.push(log_warning(&format!(
                "[WARNING] Cannot find any low-quality 0/0, 0/1 or 1/1 variant in pileup output in contig {contig}"
            )));
        }
        messages
    }

    /// The BED chunks: `ceil(len / split_bed_size)` files named `{contig}.{index}_{count}`,
    /// each holding the windows of up to `split_bed_size` consecutive candidates.
    ///
    /// # Errors
    ///
    /// [`Error::MissingContigName`] if there are candidates but no contig name, where the
    /// script's `'\t'.join` raises on `None` (after creating the first chunk, empty).
    pub fn bed_chunks(&self, window: CandidateWindow, split_bed_size: NonZeroUsize) -> Result<Vec<BedChunk>, Error> {
        if self.is_empty() {
            return Ok(Vec::new());
        }
        let Some(contig) = self.contig.as_deref() else {
            return Err(Error::MissingContigName);
        };
        let region_num = self.positions.len().div_ceil(split_bed_size.get());
        let chunks = self.positions.chunks(split_bed_size.get()).enumerate().map(|(idx, positions)| {
            let mut regions: Vec<(i128, i128)> = positions.iter().map(|&pos| window.bed_region(pos)).collect();
            // `sorted(split_output, key=lambda x: x[0])`, stable; the starts already ascend.
            regions.sort_by_key(|&(start, _)| start);
            let contents: String = regions.iter().map(|(start, end)| format!("{contig}\t{start}\t{end}\n")).collect();
            BedChunk { file_name: format!("{contig}.{idx}_{region_num}"), contents }
        });
        Ok(chunks.collect())
    }

    /// Writes the [`bed_chunks`](Self::bed_chunks) to `split_folder` and, if there is at least
    /// one, `split_folder/FULL_ALN_FILE_{contig}` listing their paths (`split_folder` joined
    /// with each name, one per line). Returns the path of that list, or `None` when there is no
    /// candidate and nothing is written.
    ///
    /// # Errors
    ///
    /// Fails as [`bed_chunks`](Self::bed_chunks) does, or if a file cannot be written.
    pub fn write_beds(
        &self,
        split_folder: &Path,
        window: CandidateWindow,
        split_bed_size: NonZeroUsize,
    ) -> Result<Option<PathBuf>, Error> {
        let chunks = self.bed_chunks(window, split_bed_size)?;
        if chunks.is_empty() {
            return Ok(None);
        }
        let mut list = Vec::new();
        for chunk in &chunks {
            let path = split_folder.join(&chunk.file_name);
            write_file(&path, chunk.contents.as_bytes())?;
            list.extend_from_slice(path.as_os_str().as_encoded_bytes());
            list.push(b'\n');
        }
        let list_path = split_folder.join(format!("{FULL_ALN_FILE_PREFIX}{}", self.contig_label()));
        write_file(&list_path, &list)?;
        Ok(Some(list_path))
    }
}

fn write_file(path: &Path, contents: &[u8]) -> Result<(), Error> {
    fs::write(path, contents).map_err(|source| Error::Io { path: path.to_owned(), source })
}

/// One candidate BED file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BedChunk {
    file_name: String,
    contents: String,
}

impl BedChunk {
    /// The file name, `{contig}.{index}_{count}`.
    #[must_use]
    pub fn file_name(&self) -> &str {
        &self.file_name
    }

    /// The BED lines (`contig`, `start`, `end`), each ending in `\n`.
    #[must_use]
    pub fn contents(&self) -> &str {
        &self.contents
    }
}

/// Selects one contig's candidates from a decompressed pileup VCF.
///
/// A record whose ALT is `.` or equal to REF is a reference call and any other a variant
/// call. A position selected more than once (say, as both kinds) is one candidate.
///
/// # Errors
///
/// Fails on unreadable input, a record of `contig` with fewer than 6 tab-separated columns,
/// an invalid POS or QUAL, or where the script raises on a non-finite proportion. With
/// [`CandidateFilter::LowestProportions`], a NaN QUAL is refused, as Python's sort order for it
/// cannot be reproduced.
pub fn select_candidates(
    vcf: impl BufRead,
    contig: Option<&str>,
    filter: CandidateFilter,
) -> Result<Candidates, Error> {
    // `if contig_name and contig_name != ctg_name`: an empty name filters nothing.
    let wanted = contig.filter(|c| !c.is_empty());
    let mut ref_calls = Vec::new();
    let mut variant_calls = Vec::new();
    let mut lines = Lines::new(vcf, Newlines::Universal);
    while let Some((line, row)) = lines.next_line()? {
        // `row[0] == '#'` on the line with its newline, so an empty line is a record.
        if row.starts_with('#') {
            continue;
        }
        // Tab-separated, unlike `SelectQual` and `SelectHetSnp`.
        let columns: Vec<&str> = py::rstrip(&row).split('\t').collect();
        if wanted.is_some_and(|c| c != columns[0]) {
            continue;
        }
        let pos_text = column(&columns, line, 1)?;
        let pos = py::parse_int(pos_text).ok_or_else(|| Error::InvalidNumber { line, value: pos_text.to_owned() })?;
        let (ref_base, alt_base) = (column(&columns, line, 3)?, column(&columns, line, 4)?);
        let qual_text = column(&columns, line, 5)?;
        let qual =
            py::parse_float(qual_text).ok_or_else(|| Error::InvalidNumber { line, value: qual_text.to_owned() })?;
        let call = Call { pos, qual, line };
        if alt_base == "." || ref_base == alt_base {
            ref_calls.push(call);
        } else {
            variant_calls.push(call);
        }
    }

    let (ref_positions, variant_positions) = match filter {
        CandidateFilter::BelowCutoffs(QualCutoffs { var_cutoff, ref_cutoff }) => {
            (below_cutoff(ref_calls, ref_cutoff), below_cutoff(variant_calls, var_cutoff))
        }
        CandidateFilter::LowestProportions { var_pct_full, ref_pct_full } => {
            (lowest_proportion(ref_calls, ref_pct_full)?, lowest_proportion(variant_calls, var_pct_full)?)
        }
    };
    let (low_qual_reference_calls, low_qual_variants) = (ref_positions.len(), variant_positions.len());
    let mut positions = ref_positions;
    positions.extend(variant_positions);
    positions.sort_unstable();
    positions.dedup();
    Ok(Candidates { contig: contig.map(str::to_owned), low_qual_reference_calls, low_qual_variants, positions })
}

/// What [`select_candidates_for_full_alignment`] selected and wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateBeds {
    candidates: Candidates,
    full_aln_file: Option<PathBuf>,
    merge_output: Option<PathBuf>,
    messages: Vec<String>,
}

impl CandidateBeds {
    /// The selected candidates.
    #[must_use]
    pub fn candidates(&self) -> &Candidates {
        &self.candidates
    }

    /// The `FULL_ALN_FILE_{contig}` written, or `None` when there was no candidate.
    #[must_use]
    pub fn full_aln_file(&self) -> Option<&Path> {
        self.full_aln_file.as_deref()
    }

    /// Where the pileup VCF was copied as the final output because there was no candidate,
    /// if it was attempted (see [`select_candidates_for_full_alignment`]).
    #[must_use]
    pub fn merge_output(&self) -> Option<&Path> {
        self.merge_output.as_deref()
    }

    /// Every line (without newlines) the script prints to stderr, in order. It prints nothing
    /// to stdout.
    #[must_use]
    pub fn messages(&self) -> &[String] {
        &self.messages
    }
}

/// `print_calling_step(output_fn=pileup_vcf_fn)`: when the pileup VCF is named `pileup.vcf`
/// or `pileup.vcf.gz`, copies `pileup.vcf.gz` from its directory to `merge_output.vcf.gz`
/// there. Returns that path and the warning printed.
fn copy_pileup_as_merge_output(pileup_vcf_fn: &Path) -> Option<(PathBuf, String)> {
    // `os.path.basename` of a path ending in `/` is empty.
    if pileup_vcf_fn.as_os_str().as_encoded_bytes().ends_with(b"/") {
        return None;
    }
    let name = pileup_vcf_fn.file_name()?;
    if name != "pileup.vcf" && name != "pileup.vcf.gz" {
        return None;
    }
    let dir = pileup_vcf_fn.parent().unwrap_or(Path::new(""));
    let merge_output = dir.join("merge_output.vcf.gz");
    let warning = log_warning(&format!("[WARNING] Copying pileup.vcf.gz to {}", merge_output.display()));
    // Reproduced: the script runs `cp` and ignores its failure (say, a plain `pileup.vcf`
    // with no `pileup.vcf.gz` beside it).
    let _ = fs::copy(dir.join("pileup.vcf.gz"), &merge_output);
    Some((merge_output, warning))
}

/// `SelectCandidates --pileup_vcf_fn <pileup_vcf_fn> --split_folder <split_folder> ...` on the
/// pipeline's path: reads the cut-offs from `split_folder/qual_fn` (falling back to
/// [`CandidateFilter::LowestProportions`] when the file is absent), selects the candidates
/// and writes the BED chunks and `FULL_ALN_FILE_{contig}` with [`Candidates::write_beds`].
///
/// With no candidate nothing is written in `split_folder`; instead, if the pileup VCF is
/// named `pileup.vcf` or `pileup.vcf.gz`, the `pileup.vcf.gz` beside it is copied to
/// `merge_output.vcf.gz` in the same directory, as the script does (copy errors ignored).
///
/// # Errors
///
/// [`Error::UnsupportedOption`] for options off the ported path, or fails as
/// [`read_qual_cutoffs`], [`open_vcf`], [`select_candidates`] and [`Candidates::write_beds`]
/// do. (The script skips a pileup VCF that does not exist and writes nothing; here that is an
/// error.)
pub fn select_candidates_for_full_alignment(
    pileup_vcf_fn: &Path,
    split_folder: &Path,
    options: &CandidateOptions,
) -> Result<CandidateBeds, Error> {
    options.ensure_ported()?;
    select_candidates_for_full_alignment_from(open_vcf(pileup_vcf_fn)?, pileup_vcf_fn, split_folder, options)
}

/// As [`select_candidates_for_full_alignment`], reading the records of `pileup_vcf_fn` from
/// `vcf` (all of them, or only the contig's).
///
/// # Errors
///
/// As [`select_candidates_for_full_alignment`].
pub fn select_candidates_for_full_alignment_from(
    vcf: impl BufRead,
    pileup_vcf_fn: &Path,
    split_folder: &Path,
    options: &CandidateOptions,
) -> Result<CandidateBeds, Error> {
    options.ensure_ported()?;
    let filter = match read_qual_cutoffs(split_folder, &options.qual_fn)? {
        Some(cutoffs) => CandidateFilter::BelowCutoffs(cutoffs),
        None => CandidateFilter::LowestProportions {
            var_pct_full: options.var_pct_full,
            ref_pct_full: options.ref_pct_full,
        },
    };
    let candidates = select_candidates(vcf, options.contig.as_deref(), filter)?;
    let mut messages = candidates.messages();
    let full_aln_file = candidates.write_beds(split_folder, options.window, options.split_bed_size)?;
    let mut merge_output = None;
    if full_aln_file.is_none()
        && let Some((path, warning)) = copy_pileup_as_merge_output(pileup_vcf_fn)
    {
        messages.push(warning);
        merge_output = Some(path);
    }
    Ok(CandidateBeds { candidates, full_aln_file, merge_output, messages })
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

    fn vcf(records: &[(&str, i64, &str, &str, &str)]) -> String {
        let mut text = HEADER.to_owned();
        for (ctg, pos, ref_base, alt, qual) in records {
            text += &format!("{ctg}\t{pos}\t.\t{ref_base}\t{alt}\t{qual}\tPASS\tP\tGT:GQ\t0/1:5\n");
        }
        text
    }

    const CUTOFFS: CandidateFilter = CandidateFilter::BelowCutoffs(QualCutoffs { var_cutoff: 10.0, ref_cutoff: 5.0 });

    #[test]
    fn selects_strictly_below_cutoffs_and_merges_positions() {
        let input = vcf(&[
            ("chr1", 300, "A", ".", "4.0"),
            ("chr1", 100, "A", "A", "1.0"),
            ("chr1", 200, "A", ".", "5.0"),
            ("chr1", 300, "A", "C", "9.0"),
            ("chr1", 400, "A", "C", "10.0"),
            ("chr1", 50, "AT", "A", "nan"),
            ("chr2", 10, "A", "C", "1.0"),
        ]);
        let candidates = select_candidates(input.as_bytes(), Some("chr1"), CUTOFFS).expect("valid");
        assert_eq!(candidates.positions(), [100, 300]);
        assert_eq!((candidates.low_qual_reference_calls(), candidates.low_qual_variants()), (2, 1));
        assert_eq!(
            candidates.messages(),
            [
                "[INFO] Low quality reference calls to be processed in chr1: 2",
                "[INFO] Low quality variants to be processed in chr1: 1"
            ]
        );
    }

    #[test]
    fn proportions_keep_ties_in_file_order() {
        let input = vcf(&[
            ("c", 30, "A", "C", "2.0"),
            ("c", 10, "A", "C", "1.0"),
            ("c", 20, "A", "C", "2.0"),
            ("c", 40, "A", "C", "3.0"),
        ]);
        let filter = CandidateFilter::LowestProportions { var_pct_full: 0.5, ref_pct_full: 0.3 };
        let candidates = select_candidates(input.as_bytes(), None, filter).expect("valid");
        assert_eq!(candidates.positions(), [10, 30]);
        let nan = vcf(&[("c", 30, "A", "C", "nan")]);
        assert!(matches!(select_candidates(nan.as_bytes(), None, filter), Err(Error::NanQuality { line: 3 })));
    }

    #[test]
    fn chunks_split_and_clamp_windows() {
        let candidates = Candidates {
            contig: Some("chr1".to_owned()),
            low_qual_reference_calls: 0,
            low_qual_variants: 3,
            positions: vec![1, 17, 2000],
        };
        let size = NonZeroUsize::new(2).expect("non-zero");
        let chunks = candidates.bed_chunks(CandidateWindow::Flanking, size).expect("named contig");
        let names: Vec<&str> = chunks.iter().map(BedChunk::file_name).collect();
        assert_eq!(names, ["chr1.0_2", "chr1.1_2"]);
        assert_eq!(chunks[0].contents(), "chr1\t0\t18\nchr1\t0\t34\n");
        assert_eq!(chunks[1].contents(), "chr1\t1983\t2017\n");

        let chunks = candidates.bed_chunks(CandidateWindow::IlluminaBlock, DEFAULT_SPLIT_BED_SIZE).expect("named");
        assert_eq!(chunks[0].contents(), "chr1\t0\t1032\nchr1\t0\t1032\nchr1\t1966\t3032\n");
    }

    #[test]
    fn candidates_without_contig_name_cannot_be_written() {
        let candidates =
            Candidates { contig: None, low_qual_reference_calls: 1, low_qual_variants: 0, positions: vec![5] };
        assert!(matches!(
            candidates.bed_chunks(CandidateWindow::Flanking, DEFAULT_SPLIT_BED_SIZE),
            Err(Error::MissingContigName)
        ));
        let empty = Candidates { positions: Vec::new(), ..candidates };
        assert_eq!(empty.bed_chunks(CandidateWindow::Flanking, DEFAULT_SPLIT_BED_SIZE).expect("empty"), []);
        assert!(empty.messages()[2].contains("variant in pileup output in contig None"));
    }

    #[test]
    fn platform_and_unported_options() {
        assert_eq!(CandidateWindow::from_platform("ilmn"), CandidateWindow::IlluminaBlock);
        assert_eq!(CandidateWindow::from_platform("ILMN"), CandidateWindow::Flanking);
        assert!(CandidateOptions::default().ensure_ported().is_ok());
        let entropy = CandidateOptions { call_low_seq_entropy: SeqEntropyCalling::On, ..CandidateOptions::default() };
        assert!(matches!(
            entropy.ensure_ported(),
            Err(Error::UnsupportedOption { option: "--call_low_seq_entropy True" })
        ));
        let phasing = CandidateOptions { phasing_info_in_bam: PhasingInfoFlag::Given, ..CandidateOptions::default() };
        assert!(matches!(phasing.ensure_ported(), Err(Error::UnsupportedOption { option: "--phasing_info_in_bam" })));
        let alt = CandidateOptions { all_alt_fn: Some(PathBuf::from("x")), ..CandidateOptions::default() };
        assert!(matches!(alt.ensure_ported(), Err(Error::UnsupportedOption { option: "--all_alt_fn" })));
    }

    #[test]
    fn short_and_empty_records() {
        // An empty line is a record of contig "", skipped when a contig is named.
        let input = format!("{HEADER}\nchr1\t100\t.\tA\tC\t5.0\n");
        assert_eq!(select_candidates(input.as_bytes(), Some("chr1"), CUTOFFS).expect("valid").positions(), [100]);
        assert!(matches!(
            select_candidates(input.as_bytes(), None, CUTOFFS),
            Err(Error::MissingColumn { line: 3, column: 1 })
        ));
        // Split on tabs only.
        let input = format!("{HEADER}chr1 100 . A C 5.0\n");
        assert_eq!(select_candidates(input.as_bytes(), Some("chr1"), CUTOFFS).expect("other contig").len(), 0);
    }
}

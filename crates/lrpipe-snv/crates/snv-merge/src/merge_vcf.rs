//! `MergeVcf.py`: merges the pileup and full-alignment calls of a contig into one VCF.

use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::bed::{BedError, BedTree, open_pigz};
use crate::gvcf::{GvcfMerge, merge_calls};
use crate::postprocess::{HaploidDecider, ParRegions, haploid_precise, haploid_sensitive, mark_low_qual};
use crate::pytext::{TextError, TextLines, parse_float, parse_int, rstrip, split_whitespace, strip};
use crate::sort::read_text;
use crate::{Error, Gender, Haploid, RefCalls, RowError};

/// `param.no_of_positions` (`2 * flankingBaseNum + 1`): how much `MergeVcf_illumina` trims off
/// both ends of every candidate BED interval.
const CANDIDATE_PADDING: i64 = 33;

/// How pileup and full-alignment calls are combined (`--platform`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum MergeMode {
    /// Every platform but `ilmn` (`MergeVcf`, the default): a full-alignment record replaces
    /// any pileup record at its contig and position, RefCall or not.
    #[default]
    Candidates,
    /// `--platform ilmn` (`MergeVcf_illumina`): full-alignment records inside the candidate
    /// regions, pileup records outside them.
    Regions {
        /// Directory of per-chunk candidate BEDs, `<contig>.<chunk>` (`--bed_fn_prefix`).
        /// The original concatenates the contig's files into `full_aln_regions_<contig>` there,
        /// and so does [`merge_vcf`].
        candidate_beds: PathBuf,
    },
}

/// Whether non-variant blocks are merged into a gVCF as well (`--gvcf`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum GvcfOutput {
    /// VCF only (the default).
    #[default]
    Off,
    /// Also write the contig's calls and non-variant blocks to an lz4-compressed gVCF
    /// (`mergeNonVariant`), after the VCF.
    On {
        /// The gVCF to write (`--gvcf_fn`).
        gvcf: PathBuf,
        /// The sorted, lz4-compressed non-variant blocks of all chunks (`--non_var_gvcf_fn`).
        non_variant: PathBuf,
        /// Reference FASTA, for the first base of a block a variant cuts short (`--ref_fn`).
        reference: PathBuf,
    },
}

/// What [`merge_vcf`] reads and writes, mirroring `MergeVcf.py`'s arguments.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct MergeOptions {
    /// Pileup calls, plain or gzipped (`--pileup_vcf_fn`).
    pub pileup_vcf: PathBuf,
    /// Full-alignment calls, plain or gzipped (`--full_alignment_vcf_fn`).
    pub full_alignment_vcf: PathBuf,
    /// The merged VCF (`--output_fn`).
    pub output: PathBuf,
    /// Only records of this contig are merged (`--ctgName`); all of them when `None`.
    pub contig: Option<String>,
    /// How the two call sets are combined (`--platform`, with `--bed_fn_prefix`).
    pub mode: MergeMode,
    /// RefCall records (`--print_ref_calls`, default hidden).
    pub ref_calls: RefCalls,
    /// Global haploid mode (`--haploid_precise` takes precedence over `--haploid_sensitive`).
    pub haploid: Haploid,
    /// Sample sex (`--gender`).
    pub gender: Gender,
    /// BED of pseudo-autosomal regions kept diploid for [`Gender::Male`] (`--par_regions_bed`;
    /// `EMPTY` is `None`). It is read even when unused, as the original does.
    pub par_regions: Option<PathBuf>,
    /// Variants with QUAL at or below this are marked `LowQual` (`--qual`, default 2); 0
    /// disables marking, as Python's truthiness does.
    pub qual: i32,
    /// gVCF output (`--gvcf`).
    pub gvcf: GvcfOutput,
}

impl MergeOptions {
    /// Options with `MergeVcf.py`'s defaults, merging `pileup_vcf` and `full_alignment_vcf`
    /// into `output`.
    pub fn new(
        pileup_vcf: impl Into<PathBuf>,
        full_alignment_vcf: impl Into<PathBuf>,
        output: impl Into<PathBuf>,
    ) -> Self {
        Self {
            pileup_vcf: pileup_vcf.into(),
            full_alignment_vcf: full_alignment_vcf.into(),
            output: output.into(),
            contig: None,
            mode: MergeMode::Candidates,
            ref_calls: RefCalls::Hide,
            haploid: Haploid::Off,
            gender: Gender::Unknown,
            par_regions: None,
            qual: 2,
            gvcf: GvcfOutput::Off,
        }
    }
}

/// What [`merge_vcf`] wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct MergeSummary {
    /// Non-empty lines written after the header.
    pub records: u64,
    /// Pileup records the original reports as processed: kept, including rows a haploid
    /// mode emptied.
    pub pileup: u64,
    /// Full-alignment records the original reports as processed, counted the same way. On
    /// [`MergeMode::Regions`] a record replacing an earlier one at its position counts too.
    pub full_alignment: u64,
    /// The concatenated candidate BED written on [`MergeMode::Regions`], if the contig had any
    /// candidate files.
    pub candidate_bed: Option<PathBuf>,
    /// Lines written to the gVCF after its header, with [`GvcfOutput::On`].
    pub gvcf_records: Option<u64>,
}

/// How a VCF line is split into columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Columns {
    /// `row.strip().split()`.
    Whitespace,
    /// `row.rstrip().split('\t')`.
    Tabs,
}

impl Columns {
    fn split(self, row: &str) -> Vec<&str> {
        match self {
            Self::Whitespace => split_whitespace(strip(row), None),
            Self::Tabs => rstrip(row).split('\t').collect(),
        }
    }
}

/// How `--ctgName` selects records. The two loops of `MergeVcf` test it differently.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ContigTest {
    /// `contig_name != None and ctg_name != contig_name` skips the record.
    IsNotNone,
    /// `contig_name and contig_name != ctg_name` skips the record: an empty name keeps all.
    Truthy,
}

/// A record kept for output: its position and the line to write (empty when a haploid mode
/// rejected its genotype).
type Kept = (i64, String);

/// The fields `MergeVcf` parses from every selected record.
#[derive(Debug)]
struct Site<'a> {
    contig: &'a str,
    pos: i64,
    qual: f64,
    is_reference: bool,
}

/// The per-record steps both loops share: parsing, contig selection, haploid genotypes,
/// `LowQual` and the RefCall filter.
#[derive(Debug)]
struct Calls<'a> {
    options: &'a MergeOptions,
    decider: HaploidDecider,
}

impl Calls<'_> {
    fn skips(&self, test: ContigTest, ctg_name: &str) -> bool {
        match (test, self.options.contig.as_deref()) {
            (_, None) => false,
            (ContigTest::IsNotNone, Some(contig)) => contig != ctg_name,
            (ContigTest::Truthy, Some(contig)) => !contig.is_empty() && contig != ctg_name,
        }
    }

    /// Parses a record's columns in the original's order; `Ok(None)` when it is for another contig.
    fn site<'r>(&self, columns: &[&'r str], test: ContigTest) -> Result<Option<Site<'r>>, RowError> {
        let column = |index: usize| columns.get(index).copied().ok_or(RowError::MissingColumn { index });
        let ctg_name = column(0)?;
        if self.skips(test, ctg_name) {
            return Ok(None);
        }
        let pos = column(1)?;
        let pos = parse_int(pos).ok_or_else(|| RowError::InvalidInteger { value: pos.to_owned() })?;
        let qual = column(5)?;
        let qual = parse_float(qual).ok_or_else(|| RowError::InvalidFloat { value: qual.to_owned() })?;
        let (ref_base, alt_base) = (column(3)?, column(4)?);
        let is_reference = alt_base == "." || ref_base == alt_base;
        Ok(Some(Site { contig: ctg_name, pos, qual, is_reference }))
    }

    /// The line to write for a record, if it is kept: haploid genotypes, then `LowQual` for
    /// variants; RefCalls only when shown. The haploid rewrite runs (and can fail) even for a
    /// RefCall that is then dropped.
    fn keep(&self, row: &str, columns: &[&str], site: &Site<'_>) -> Result<Option<Kept>, RowError> {
        // The haploid rewrite rejoins the split columns; otherwise the line is kept as read.
        let row = match self.decider.decide(site.contig, site.pos) {
            Haploid::Precise => haploid_precise(columns)?,
            Haploid::Sensitive => haploid_sensitive(columns)?,
            Haploid::Off => row.to_owned(),
        };
        let row = if !site.is_reference {
            Some(mark_low_qual(row, self.options.qual, site.qual)?)
        } else {
            (self.options.ref_calls == RefCalls::Show).then_some(row)
        };
        Ok(row.map(|row| (site.pos, row)))
    }
}

/// The lines of a VCF read through `pigz -fdc`, with their 1-based line numbers. A file that
/// cannot be opened has no lines, as the original ignores `pigz` failing.
fn vcf_lines(path: &Path) -> impl Iterator<Item = Result<(u64, String), Error>> {
    let reader = open_pigz(path).unwrap_or_else(|e| {
        log::warn!("cannot read {}: {e}; reading it as empty, as the original does", path.display());
        Box::new(io::empty())
    });
    TextLines::new(reader).zip(1..).map(move |(line, n)| match line {
        Ok(line) => Ok((n, line)),
        Err(TextError::Io(source)) => Err(Error::Io { path: path.to_owned(), source }),
        Err(TextError::Utf8) => Err(Error::Line { path: path.to_owned(), line: n, reason: RowError::InvalidUtf8 }),
    })
}

/// The lines of `path` a merge of `contig` reads: where `path` has a tabix index (`.tbi`) and
/// `contig` names one, its header lines and then only the contig's records, read from the
/// contig's first indexed block. The file is sorted by contig, so these are the lines of the
/// whole file the merge keeps, in order, without reading the other contigs' (each of the contigs
/// merged at once read the whole pileup VCF). Otherwise every line, as [`vcf_lines`]. Records are
/// numbered from the contig's first.
fn contig_lines<'a>(
    path: &'a Path,
    contig: Option<&'a str>,
) -> Box<dyn Iterator<Item = Result<(u64, String), Error>> + 'a> {
    let Some((contig, start)) = contig.filter(|c| !c.is_empty()).and_then(|c| Some((c, contig_start(path, c)?))) else {
        return Box::new(vcf_lines(path));
    };
    let header = vcf_lines(path).take_while(|line| line.as_ref().map_or(true, |(_, row)| row.starts_with('#')));
    let Some(start) = start else { return Box::new(header) };
    let records = (|| -> io::Result<_> {
        let mut reader = pasta::bgzf::io::Reader::new(File::open(path)?);
        reader.seek(start)?;
        Ok(TextLines::new(Box::new(reader) as Box<dyn io::BufRead>))
    })();
    let records = match records {
        Ok(lines) => lines,
        Err(source) => {
            return Box::new(header.chain(std::iter::once(Err(Error::Io { path: path.to_owned(), source }))));
        }
    };
    let is_contig = move |row: &String| row.split('\t').next() == Some(contig);
    let records = records
        .map(move |line| match line {
            Ok(row) => Ok(row),
            Err(TextError::Io(source)) => Err(Error::Io { path: path.to_owned(), source }),
            Err(TextError::Utf8) => Err(Error::Line { path: path.to_owned(), line: 0, reason: RowError::InvalidUtf8 }),
        })
        .skip_while(move |row| row.as_ref().is_ok_and(|row| !is_contig(row)))
        .take_while(move |row| row.as_ref().map_or(true, &is_contig))
        .zip(1..)
        .map(|(row, n)| row.map(|row| (n, row)));
    Box::new(header.chain(records))
}

/// Where `contig`'s records start in `path`, by its tabix index: `None` without a readable index,
/// `Some(None)` if the index has no records for the contig.
fn contig_start(path: &Path, contig: &str) -> Option<Option<pasta::bgzf::VirtualPosition>> {
    use pasta::csi::BinningIndex as _;
    let mut index_path = path.as_os_str().to_owned();
    index_path.push(".tbi");
    let index = pasta::tabix::fs::read(&index_path).ok()?;
    let id = index.header()?.reference_sequence_names().get_index_of(contig.as_bytes());
    let Some(id) = id else { return Some(None) };
    let sequence = index.reference_sequences().get(id)?;
    Some(sequence.bins().values().flat_map(|bin| bin.chunks()).map(|chunk| chunk.start()).min())
}

/// Adds a header line to the output's header, ending it with a newline if the file did not: a
/// header-only VCF of SortVcf's (`get_header`, as the per-contig full-alignment VCF of a contig
/// without full-alignment calls is) has none after its `#CHROM` line, and the records written
/// after it would join that line. `MergeVcf.py` copies the line as it is, but its inputs always
/// end with a newline.
fn push_header_line(header: &mut String, row: &str) {
    header.push_str(row);
    if !row.ends_with('\n') {
        header.push('\n');
    }
}

fn line_error(path: &Path, line: u64) -> impl Fn(RowError) -> Error + '_ {
    move |reason| Error::Line { path: path.to_owned(), line, reason }
}

/// Writes `rows` to `out` in the order `heapq.merge(full_alignment, pileup)` yields them:
/// repeatedly the smaller head by `(pos, row)`, the full-alignment one on a tie. Neither input
/// needs to be sorted.
fn heap_merge(
    full_alignment: Vec<Kept>,
    mut pileup: impl Iterator<Item = Result<Kept, Error>>,
    mut write: impl FnMut(&str) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut full_alignment = full_alignment.into_iter();
    let mut fa_head = full_alignment.next();
    let mut pileup_head = pileup.next().transpose()?;
    loop {
        let take_pileup = match (&fa_head, &pileup_head) {
            (Some((fa_pos, fa_row)), Some((pos, row))) => (*pos, row.as_str()) < (*fa_pos, fa_row.as_str()),
            (None, Some(_)) => true,
            (Some(_), None) => false,
            (None, None) => return Ok(()),
        };
        if take_pileup {
            let (_, row) = pileup_head.take().expect("the pileup head is set");
            write(&row)?;
            pileup_head = pileup.next().transpose()?;
        } else {
            let (_, row) = fa_head.take().expect("the full-alignment head is set");
            write(&row)?;
            fa_head = full_alignment.next();
        }
    }
}

/// The output file, counting the records written.
#[derive(Debug)]
struct Output<'a> {
    path: &'a Path,
    out: BufWriter<File>,
    records: u64,
}

impl<'a> Output<'a> {
    fn create(path: &'a Path) -> Result<Self, Error> {
        let file = File::create(path).map_err(|source| Error::Io { path: path.to_owned(), source })?;
        Ok(Self { path, out: BufWriter::with_capacity(1 << 20, file), records: 0 })
    }

    fn write(&mut self, text: &str) -> Result<(), Error> {
        self.out.write_all(text.as_bytes()).map_err(|source| Error::Io { path: self.path.to_owned(), source })
    }

    fn write_row(&mut self, row: &str) -> Result<(), Error> {
        if !row.is_empty() {
            self.records += 1;
        }
        self.write(row)
    }

    fn finish(self) -> Result<u64, Error> {
        let Self { path, out, records } = self;
        out.into_inner().map_err(|e| Error::Io { path: path.to_owned(), source: e.into_error() })?;
        Ok(records)
    }
}

/// `MergeVcf`: full-alignment records first, kept in file order; then the pileup records at
/// positions the full-alignment VCF has no record for, merged in by position.
fn merge_candidates(options: &MergeOptions, calls: &Calls<'_>) -> Result<MergeSummary, Error> {
    let fa_path = options.full_alignment_vcf.as_path();
    let mut header = String::new();
    let mut full_alignment: Vec<Kept> = Vec::new();
    let mut called: HashSet<(String, i64)> = HashSet::new();
    for line in contig_lines(fa_path, options.contig.as_deref()) {
        let (n, row) = line?;
        if row.starts_with('#') {
            push_header_line(&mut header, &row);
            continue;
        }
        let columns = Columns::Whitespace.split(&row);
        let line_err = line_error(fa_path, n);
        let Some(site) = calls.site(&columns, ContigTest::IsNotNone).map_err(&line_err)? else { continue };
        called.insert((site.contig.to_owned(), site.pos));
        full_alignment.extend(calls.keep(&row, &columns, &site).map_err(&line_err)?);
    }
    let from_full_alignment = full_alignment.len() as u64;

    let pileup_path = options.pileup_vcf.as_path();
    let mut out = Output::create(&options.output)?;
    out.write(&header)?;
    let mut from_pileup = 0;
    let pileup = contig_lines(pileup_path, options.contig.as_deref()).filter_map(|line| {
        let (n, row) = match line {
            Ok(line) => line,
            Err(e) => return Some(Err(e)),
        };
        if row.starts_with('#') {
            return None;
        }
        let columns = Columns::Tabs.split(&row);
        let line_err = line_error(pileup_path, n);
        let site = match calls.site(&columns, ContigTest::Truthy) {
            Ok(Some(site)) => site,
            Ok(None) => return None,
            Err(reason) => return Some(Err(line_err(reason))),
        };
        if called.contains(&(site.contig.to_owned(), site.pos)) {
            return None;
        }
        match calls.keep(&row, &columns, &site) {
            Ok(Some(kept)) => {
                from_pileup += 1;
                Some(Ok(kept))
            }
            Ok(None) => None,
            Err(reason) => Some(Err(line_err(reason))),
        }
    });
    heap_merge(full_alignment, pileup, |row| out.write_row(row))?;
    let records = out.finish()?;
    let contig = options.contig.as_deref().unwrap_or("None");
    log::info!("[INFO] Pileup variants processed in {contig}: {from_pileup}");
    log::info!("[INFO] Full-alignment variants processed in {contig}: {from_full_alignment}");
    Ok(MergeSummary {
        records,
        pileup: from_pileup,
        full_alignment: from_full_alignment,
        candidate_bed: None,
        gvcf_records: None,
    })
}

/// The candidate regions of `contig`: its files in `dir` (named `<contig>.*`, in directory
/// order) concatenated into `full_aln_regions_<contig>` and read with the padding.
fn candidate_regions(dir: &Path, contig: Option<&str>) -> Result<(BedTree, Option<PathBuf>), Error> {
    if !dir.exists() {
        return Err(Error::InputDir { path: dir.to_owned() });
    }
    let io_err = |source| Error::Io { path: dir.to_owned(), source };
    let mut names = Vec::new();
    for entry in fs::read_dir(dir).map_err(io_err)? {
        names.push(entry.map_err(io_err)?.file_name());
    }
    let contig = contig.ok_or(Error::ContigUnset)?;
    let prefix = format!("{contig}.");
    let files: Vec<PathBuf> =
        names.iter().filter(|name| name.to_string_lossy().starts_with(&prefix)).map(|name| dir.join(name)).collect();
    if files.is_empty() {
        return Ok((BedTree::default(), None));
    }
    let bed = dir.join(format!("full_aln_regions_{contig}"));
    let text = files.iter().map(|file| read_text(file)).collect::<Result<String, Error>>()?;
    fs::write(&bed, &text).map_err(|source| Error::Io { path: bed.clone(), source })?;
    let tree = BedTree::read(text.as_bytes(), Some(contig), CANDIDATE_PADDING).map_err(|e| match e {
        BedError::Io(source) => Error::Io { path: bed.clone(), source },
        BedError::Row { line, reason } => Error::CandidateBed { path: bed.clone(), line, reason },
    })?;
    Ok((tree, Some(bed)))
}

/// `MergeVcf_illumina`: pileup records outside the candidate regions and full-alignment
/// records inside them, by position (a later record at a position replaces an earlier one),
/// under the pileup VCF's header.
fn merge_regions(options: &MergeOptions, calls: &Calls<'_>, candidate_beds: &Path) -> Result<MergeSummary, Error> {
    let (tree, candidate_bed) = candidate_regions(candidate_beds, options.contig.as_deref())?;
    let in_regions = |site: &Site<'_>| !tree.is_empty() && tree.contains(site.contig, site.pos);
    let mut header = String::new();
    let mut rows: BTreeMap<i64, String> = BTreeMap::new();
    let mut counts = [0u64; 2];
    let sources = [(options.pileup_vcf.as_path(), false), (options.full_alignment_vcf.as_path(), true)];
    for (count, (path, inside)) in counts.iter_mut().zip(sources) {
        for line in contig_lines(path, options.contig.as_deref()) {
            let (n, row) = line?;
            if row.starts_with('#') {
                if !inside {
                    push_header_line(&mut header, &row);
                }
                continue;
            }
            let columns = Columns::Whitespace.split(&row);
            let line_err = line_error(path, n);
            let Some(site) = calls.site(&columns, ContigTest::IsNotNone).map_err(&line_err)? else { continue };
            let Some((pos, row)) = calls.keep(&row, &columns, &site).map_err(&line_err)? else { continue };
            if in_regions(&site) == inside {
                rows.insert(pos, row);
                *count += 1;
            }
        }
    }
    let [from_pileup, from_full_alignment] = counts;
    let contig = options.contig.as_deref().unwrap_or("None");
    log::info!("[INFO] Pileup positions variants proceeded in {contig}: {from_pileup}");
    log::info!("[INFO] Realigned positions variants proceeded in {contig}: {from_full_alignment}");
    let mut out = Output::create(&options.output)?;
    out.write(&header)?;
    for row in rows.values() {
        out.write_row(row)?;
    }
    let records = out.finish()?;
    Ok(MergeSummary {
        records,
        pileup: from_pileup,
        full_alignment: from_full_alignment,
        candidate_bed,
        gvcf_records: None,
    })
}

/// Merges the pileup and full-alignment calls of [`MergeOptions::contig`] into
/// [`MergeOptions::output`], as `MergeVcf.py` does.
///
/// With [`MergeMode::Candidates`] (every platform but `ilmn`) the header is the full-alignment
/// VCF's. Its records come first in file order, merged by position with the pileup records at
/// sites it has no record for. With [`MergeMode::Regions`] the header is the pileup VCF's and
/// each region of the candidate BEDs takes its calls from full alignment. Either way RefCalls
/// are dropped unless [`MergeOptions::ref_calls`] shows them, haploid genotypes are rewritten
/// (a genotype the mode rejects leaves nothing to write) and variants at or below
/// [`MergeOptions::qual`] are marked `LowQual`.
///
/// An input VCF that cannot be opened is read as empty, as the original ignores `pigz` failing.
///
/// With [`GvcfOutput::On`] the merged VCF is then interleaved with the non-variant blocks into
/// the gVCF (`mergeNonVariant`): variant records gain a `<NON_REF>` allele, RefCalls are left
/// out, and blocks are cut around the variants and at ploidy boundaries.
///
/// # Errors
///
/// Fails where the original raises or exits: a line it cannot parse, a malformed PAR or
/// candidate BED, a missing candidate BED directory or contig name on [`MergeMode::Regions`]. I/O
/// failures are errors too. After an error the output files may be missing or partly written.
pub fn merge_vcf(options: &MergeOptions) -> Result<MergeSummary, Error> {
    let par = match &options.par_regions {
        Some(path) => ParRegions::read(path)?,
        None => ParRegions::default(),
    };
    let calls = Calls { options, decider: HaploidDecider::new(options.haploid, options.gender, par) };
    let mut summary = match &options.mode {
        MergeMode::Candidates => merge_candidates(options, &calls)?,
        MergeMode::Regions { candidate_beds } => merge_regions(options, &calls, candidate_beds)?,
    };
    if let GvcfOutput::On { gvcf, non_variant, reference } = &options.gvcf {
        log::info!("[INFO] Merge variants and non-variants to GVCF");
        summary.gvcf_records = Some(merge_calls(&GvcfMerge {
            calls: &options.output,
            non_variant,
            output: gvcf,
            reference,
            contig: options.contig.as_deref(),
            decider: &calls.decider,
        })?);
    }
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kept(rows: &[(i64, &str)]) -> Vec<Kept> {
        rows.iter().map(|&(pos, row)| (pos, row.to_owned())).collect()
    }

    fn merged(full_alignment: &[(i64, &str)], pileup: &[(i64, &str)]) -> Vec<String> {
        let mut out = Vec::new();
        let pileup = kept(pileup).into_iter().map(Ok);
        heap_merge(kept(full_alignment), pileup, |row| {
            out.push(row.to_owned());
            Ok(())
        })
        .expect("no errors");
        out
    }

    #[test]
    fn heap_merge_matches_heapq() {
        assert_eq!(merged(&[(1, "a"), (5, "e")], &[(2, "b"), (5, "d"), (9, "z")]), ["a", "b", "d", "e", "z"]);
        // Equal (pos, row): full alignment first.
        assert_eq!(merged(&[(3, "x")], &[(3, "x")]), ["x", "x"]);
        // Unsorted input is not re-sorted: each head is compared as it comes.
        assert_eq!(merged(&[(5, "e"), (1, "a")], &[(3, "c"), (2, "b")]), ["c", "b", "e", "a"]);
        // An emptied (haploid) row sorts before other rows at its position.
        assert_eq!(merged(&[(4, "f")], &[(4, "")]), ["", "f"]);
        assert_eq!(merged(&[], &[(1, "p")]), ["p"]);
        assert_eq!(merged(&[(1, "f")], &[]), ["f"]);
    }

    #[test]
    fn heap_merge_stops_at_pileup_error() {
        let mut out = Vec::new();
        let pileup = vec![Ok((1, "p".to_owned())), Err(Error::ContigUnset)].into_iter();
        let result = heap_merge(kept(&[(2, "f")]), pileup, |row| {
            out.push(row.to_owned());
            Ok(())
        });
        assert!(matches!(result, Err(Error::ContigUnset)));
        assert_eq!(out, ["p"]);
    }

    /// `site` then `keep`, as both loops run them.
    fn process<'r>(
        calls: &Calls<'_>,
        row: &'r str,
        split: Columns,
        test: ContigTest,
    ) -> Result<Option<(Site<'r>, Option<Kept>)>, RowError> {
        let columns = split.split(row);
        let Some(site) = calls.site(&columns, test)? else { return Ok(None) };
        let kept = calls.keep(row, &columns, &site)?;
        Ok(Some((site, kept)))
    }

    fn calls(options: &MergeOptions) -> Calls<'_> {
        Calls { options, decider: HaploidDecider::new(options.haploid, options.gender, ParRegions::default()) }
    }

    #[test]
    fn records_are_filtered_and_marked() {
        let mut options = MergeOptions::new("p", "f", "o");
        options.contig = Some("chr1".to_owned());
        options.qual = 10;
        let calls = calls(&options);
        let row = "chr1\t7\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n";
        let (site, kept) = process(&calls, row, Columns::Tabs, ContigTest::Truthy).expect("valid").expect("selected");
        assert_eq!((site.contig, site.pos), ("chr1", 7));
        assert_eq!(kept, Some((7, row.replace("PASS", "LowQual"))));
        let refcall = "chr1\t7\t.\tA\t.\t5\tRefCall\tP\tGT\t0/0\n";
        let (_, kept) = process(&calls, refcall, Columns::Tabs, ContigTest::Truthy).expect("valid").expect("selected");
        assert_eq!(kept, None);
        let other = "chr2\t7\t.\tA\tC\t5\tPASS\n";
        assert!(process(&calls, other, Columns::Tabs, ContigTest::Truthy).expect("valid").is_none());
        let short = "chr1\t7\t.\tA\n";
        assert_eq!(
            process(&calls, short, Columns::Tabs, ContigTest::Truthy).err(),
            Some(RowError::MissingColumn { index: 5 })
        );
    }

    #[test]
    fn empty_contig_name_selects_differently() {
        let mut options = MergeOptions::new("p", "f", "o");
        options.contig = Some(String::new());
        let empty = calls(&options);
        assert!(empty.skips(ContigTest::IsNotNone, "chr1"));
        assert!(!empty.skips(ContigTest::Truthy, "chr1"));
        // A blank pileup line has contig "" and is skipped only by a non-empty --ctgName.
        let mut options = MergeOptions::new("p", "f", "o");
        options.contig = Some("chr1".to_owned());
        assert!(calls(&options).skips(ContigTest::Truthy, ""));
    }

    #[test]
    fn haploid_rewrite_rejoins_columns() {
        let mut options = MergeOptions::new("p", "f", "o");
        options.haploid = Haploid::Precise;
        let calls = calls(&options);
        let row = "chr1 7 . A C 30 PASS F GT:GQ 1/1:9  \n";
        let (_, kept) =
            process(&calls, row, Columns::Whitespace, ContigTest::IsNotNone).expect("valid").expect("selected");
        assert_eq!(kept, Some((7, "chr1\t7\t.\tA\tC\t30\tPASS\tF\tGT:GQ\t1:9\n".to_owned())));
    }
}

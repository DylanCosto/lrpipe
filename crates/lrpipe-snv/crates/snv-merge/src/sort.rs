//! `sort_vcf_from` (`SortVcf.py`): reading the chunk VCFs, ordering and writing the merge.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use snv_call::{Lz4Writer, open_lz4};

use crate::malformed::{Anchor, Checked, check_header, check_record};
use crate::postprocess::{HaploidDecider, ParRegions, RowFilter};
use crate::pytext::{TextError, TextLines, parse_int, rstrip, split_whitespace, strip, translate_newlines};
use crate::{
    ChunkLayout, EmptyReason, Error, GvcfContigs, IndexStatus, MalformedRecords, Postprocess, RowError, SortOptions,
    Summary, compress_and_index,
};

/// The rank of `contig` in `major_contigs_order`: `chr1`..`chr22`, `chrX`, `chrY`, then the
/// same names without `chr`.
fn major_contig_rank(contig: &str) -> Option<usize> {
    let (bare, offset) = match contig.strip_prefix("chr") {
        Some(bare) => (bare, 0),
        None => (contig, 24),
    };
    let rank = match bare {
        "X" => 22,
        "Y" => 23,
        // `str(n)` for 1..=22: no sign, no leading zeros.
        _ if !bare.starts_with(['0', '+']) => bare.parse::<usize>().ok().filter(|n| (1..=22).contains(n))? - 1,
        _ => return None,
    };
    Some(offset + rank)
}

/// `sorted(all_contigs_list, key=contigs_order.index)`, where `contigs_order` is
/// `major_contigs_order` followed by `all_contigs_list`: major contigs first in their fixed
/// order, then the others in file order. The sort is stable, so a contig listed twice is
/// merged (and written) twice.
fn contig_order(contigs: &[String]) -> Vec<String> {
    let first_index: HashMap<&str, usize> =
        contigs.iter().enumerate().rev().map(|(i, contig)| (contig.as_str(), i)).collect();
    let mut order = contigs.to_vec();
    order.sort_by_key(|contig| major_contig_rank(contig).unwrap_or_else(|| 48 + first_index[contig.as_str()]));
    order
}

/// One entry of `os.listdir(input_dir)`.
#[derive(Debug)]
struct InputFile {
    /// The file name as Python matches it (non-UTF-8 bytes replaced).
    name: String,
    path: PathBuf,
}

/// The directory entries in `readdir` order, which is the order `os.listdir` returns. That
/// order decides which of two rows at the same position wins, and the header line order.
/// macOS can write AppleDouble metadata beside chunks on filesystems such as exFAT.
fn list_inputs(dir: &Path) -> Result<Vec<InputFile>, Error> {
    let io_err = |source| Error::Io { path: dir.to_owned(), source };
    let mut files = Vec::new();
    for entry in fs::read_dir(dir).map_err(io_err)? {
        let entry = entry.map_err(io_err)?;
        let name = entry.file_name().to_string_lossy().into_owned();
        let path = entry.path();
        if name.starts_with("._") {
            let mut magic = [0; 4];
            if File::open(&path).and_then(|mut file| file.read_exact(&mut magic)).is_ok()
                && magic == [0x00, 0x05, 0x16, 0x07]
            {
                continue;
            }
        }
        files.push(InputFile { name, path });
    }
    Ok(files)
}

/// Reads a text file as `open(path).read()` does.
pub(crate) fn read_text(path: &Path) -> Result<String, Error> {
    let bytes = fs::read(path).map_err(|source| Error::Io { path: path.to_owned(), source })?;
    let text = String::from_utf8(bytes).map_err(|_| Error::Line {
        path: path.to_owned(),
        line: 1,
        reason: RowError::InvalidUtf8,
    })?;
    Ok(translate_newlines(&text))
}

/// `[item.rstrip() for item in open(contigs_fn).readlines()]`.
fn read_contigs(path: Option<&Path>) -> Result<Vec<String>, Error> {
    let Some(path) = path.filter(|p| !p.as_os_str().is_empty() && p.exists()) else {
        return Err(Error::ContigsFile { path: path.map(Path::to_owned) });
    };
    let text = read_text(path)?;
    Ok(text.split_inclusive('\n').map(|line| rstrip(line).to_owned()).collect())
}

/// The first two whitespace-separated columns of `row` (`row.split(maxsplit=3)[:2]`).
fn contig_and_pos(row: &str) -> Result<(&str, &str), RowError> {
    let columns = split_whitespace(row, Some(3));
    let contig = columns.first().ok_or(RowError::MissingColumn { index: 0 })?;
    let pos = columns.get(1).ok_or(RowError::MissingColumn { index: 1 })?;
    Ok((contig, pos))
}

/// `int(pos)`.
fn position(pos: &str) -> Result<i64, RowError> {
    parse_int(pos).ok_or_else(|| RowError::InvalidInteger { value: pos.to_owned() })
}

/// The header lines seen so far, in first-seen order (`if row not in header`).
#[derive(Debug, Default)]
struct Header {
    rows: Vec<String>,
    seen: HashSet<String>,
    written: bool,
}

impl Header {
    fn push(&mut self, row: String) {
        if self.seen.insert(row.clone()) {
            self.rows.push(row);
        }
    }
}

/// `check_header_in_gvcf`: drops the `##contig` lines of contigs not in `contigs`, reading the
/// contig as `row.split(',')[0].split('=')[2]`.
fn gvcf_header_rows(rows: Vec<String>, contigs: &[String]) -> Result<Vec<String>, Error> {
    let mut kept = Vec::with_capacity(rows.len());
    for row in rows {
        if row.starts_with("##contig=") {
            let Some(contig) = row.split(',').next().and_then(|s| s.split('=').nth(2)) else {
                return Err(Error::Record { row, reason: RowError::MissingColumn { index: 2 } });
            };
            if !contigs.iter().any(|c| c == contig) {
                continue;
            }
        }
        kept.push(row);
    }
    Ok(kept)
}

/// How the inputs and the output are stored, by `--vcf_fn_suffix`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    /// Plain VCFs in, a bgzipped VCF out.
    Vcf,
    /// lz4 gVCFs in (`.gvcf`), a bgzipped gVCF out.
    Gvcf,
    /// lz4 gVCFs in and out (`.tmp.gvcf` or any other suffix containing `gvcf`).
    Lz4Gvcf,
}

impl Format {
    fn of(suffix: &str) -> Self {
        if suffix == ".gvcf" {
            Self::Gvcf
        } else if suffix.contains("gvcf") {
            Self::Lz4Gvcf
        } else {
            Self::Vcf
        }
    }
}

/// The output file, plain or lz4-compressed (`lz4 -c`).
enum Sink {
    Plain(BufWriter<File>),
    Lz4(Lz4Writer<BufWriter<File>>),
}

impl Sink {
    fn writer(&mut self) -> &mut dyn Write {
        match self {
            Self::Plain(out) => out,
            Self::Lz4(out) => out,
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(out) => out.into_inner().map(drop).map_err(io::IntoInnerError::into_error),
            Self::Lz4(out) => out.finish().map(drop),
        }
    }
}

/// The state `sort_vcf_from` keeps while merging.
struct Merge<'a> {
    options: &'a SortOptions,
    filter: Option<RowFilter>,
    header: Header,
    /// Lines read, headers included (`row_count`).
    lines: u64,
    /// Whether any record was kept (`not no_vcf_output`).
    any_record: bool,
    /// Rows written after the header.
    written: u64,
    format: Format,
    /// The listed contigs, whose `##contig` lines a final gVCF keeps.
    gvcf_contigs: Option<Vec<String>>,
    out: Sink,
}

impl Merge<'_> {
    fn write(&mut self, row: &str) -> Result<(), Error> {
        self.out
            .writer()
            .write_all(row.as_bytes())
            .map_err(|source| Error::Io { path: self.options.output.clone(), source })
    }

    fn write_row(&mut self, row: &str) -> Result<(), Error> {
        if !row.is_empty() {
            self.written += 1;
        }
        self.write(row)
    }

    /// Writes the header the first time there is one, adding `##cmdline` from the command
    /// file if no input header had it.
    fn write_header_once(&mut self) -> Result<(), Error> {
        if self.header.written || self.header.rows.is_empty() {
            return Ok(());
        }
        let rows = &mut self.header.rows;
        if !rows.iter().any(|row| row.contains("##cmdline")) {
            // The original calls os.path.exists(cmd_fn) here, which raises for a missing --cmd_fn.
            let cmd_file = self.options.cmd_file.as_deref().ok_or(Error::CmdFileUnset)?;
            if cmd_file.exists() {
                let text = read_text(cmd_file)?;
                let cmd_line = rstrip(&text);
                if !cmd_line.is_empty() {
                    let at = if rows.len() >= 3 { 3 } else { rows.len() - 1 };
                    rows.insert(at, format!("##cmdline={cmd_line}\n"));
                }
            }
        }
        if let Some(contigs) = &self.gvcf_contigs {
            self.header.rows = gvcf_header_rows(std::mem::take(&mut self.header.rows), contigs)?;
        }
        let header = self.header.rows.concat();
        self.write(&header)?;
        self.header.written = true;
        Ok(())
    }

    /// The lines of `file` (through `lz4 -fdc` for gVCFs), with their 1-based line numbers.
    fn lines_of(&self, file: &InputFile) -> Result<impl Iterator<Item = Result<(u64, String), Error>> + use<>, Error> {
        let path = file.path.clone();
        let open_err = |source| Error::Io { path: path.clone(), source };
        let reader: Box<dyn BufRead> = match self.format {
            Format::Vcf => Box::new(BufReader::new(File::open(&path).map_err(open_err)?)),
            Format::Gvcf | Format::Lz4Gvcf => open_lz4(&path).map_err(open_err)?,
        };
        let lines = TextLines::new(reader).zip(1..).map(move |(line, n)| match line {
            Ok(line) => Ok((n, line)),
            Err(TextError::Io(source)) => Err(Error::Io { path: path.clone(), source }),
            Err(TextError::Utf8) => Err(Error::Line { path: path.clone(), line: n, reason: RowError::InvalidUtf8 }),
        });
        Ok(lines)
    }

    /// The non-GPU path: per contig, the files whose name contains it, each read up to its
    /// first record of another contig.
    fn by_contig(&mut self, files: &[InputFile], order: &[String]) -> Result<(), Error> {
        let repair = self.options.malformed == MalformedRecords::Repair;
        for contig in order {
            let mut rows: BTreeMap<i64, String> = BTreeMap::new();
            // A plain substring test, so `chr1` also opens `chr11` files (stopped below).
            for file in files.iter().filter(|f| f.name.contains(contig.as_str())) {
                for line in self.lines_of(file)? {
                    let (n, row) = line?;
                    let line_err = |reason| Error::Line { path: file.path.clone(), line: n, reason };
                    self.lines += 1;
                    if row.starts_with('#') {
                        if row.starts_with("##commandline=") {
                            continue;
                        }
                        let row = if repair {
                            let Checked { row, extra } =
                                check_header(&row, Anchor::Contig(contig), &file.path).map_err(line_err)?;
                            let Some(row) = row else { continue };
                            if let Some(extra) = extra {
                                let pos =
                                    contig_and_pos(&extra).and_then(|(_, pos)| position(pos)).map_err(line_err)?;
                                rows.insert(pos, extra);
                            }
                            row
                        } else {
                            row
                        };
                        self.header.push(row);
                        continue;
                    }
                    let (ctg_name, pos) = contig_and_pos(strip(&row)).map_err(line_err)?;
                    if ctg_name != contig {
                        break;
                    }
                    let pos = pos.to_owned();
                    let row = if repair {
                        let Checked { row, extra } = check_record(&row, contig, &file.path);
                        if let Some(extra) = extra {
                            let extra_pos =
                                contig_and_pos(&extra).and_then(|(_, pos)| position(pos)).map_err(line_err)?;
                            rows.insert(extra_pos, extra);
                        }
                        row
                    } else {
                        Some(row)
                    };
                    let Some(row) = row.filter(|row| !row.is_empty()) else { continue };
                    rows.insert(position(&pos).map_err(line_err)?, row);
                    self.any_record = true;
                }
            }
            self.write_header_once()?;
            for row in rows.values() {
                // The original writes the input row here, not the post-processed one, so only
                // its RefCall filter took effect on this path (LRP-012).
                self.write_processed(row)?;
            }
        }
        Ok(())
    }

    /// The GPU path: every file once, records grouped by their own contig.
    fn unordered(&mut self, files: &[InputFile], order: &[String]) -> Result<(), Error> {
        let repair = self.options.malformed == MalformedRecords::Repair;
        let sample = self.options.sample_name.as_str();
        let mut contigs: HashMap<String, BTreeMap<i64, String>> = HashMap::new();
        for file in files {
            for line in self.lines_of(file)? {
                let (n, row) = line?;
                let line_err = |reason| Error::Line { path: file.path.clone(), line: n, reason };
                self.lines += 1;
                if row.starts_with('#') {
                    if row.starts_with("##commandline=") {
                        continue;
                    }
                    let row = if repair {
                        let Checked { row, extra } =
                            check_header(&row, Anchor::Sample(sample), &file.path).map_err(line_err)?;
                        let Some(row) = row else { continue };
                        if let Some(extra) = extra {
                            let (ctg_name, pos) = contig_and_pos(&extra).map_err(line_err)?;
                            let pos = position(pos).map_err(line_err)?;
                            contigs.entry(ctg_name.to_owned()).or_default().insert(pos, extra);
                        }
                        row
                    } else {
                        row
                    };
                    self.header.push(row);
                    continue;
                }
                let (ctg_name, pos) = contig_and_pos(strip(&row)).map_err(line_err)?;
                let (ctg_name, pos) = (ctg_name.to_owned(), pos.to_owned());
                let row = if repair {
                    let Checked { row, extra } = check_record(&row, &ctg_name, &file.path);
                    if let Some(extra) = extra {
                        let extra_pos = contig_and_pos(&extra).and_then(|(_, pos)| position(pos)).map_err(line_err)?;
                        contigs.entry(ctg_name.clone()).or_default().insert(extra_pos, extra);
                    }
                    row
                } else {
                    Some(row)
                };
                let Some(row) = row.filter(|row| !row.is_empty()) else { continue };
                let pos = position(&pos).map_err(line_err)?;
                contigs.entry(ctg_name).or_default().insert(pos, row);
                self.any_record = true;
            }
        }
        self.write_header_once()?;
        for contig in order {
            let Some(rows) = contigs.get(contig) else { continue };
            for row in rows.values() {
                self.write_processed(row)?;
            }
        }
        Ok(())
    }

    /// Writes `row` as post-processing leaves it, if it leaves anything.
    fn write_processed(&mut self, row: &str) -> Result<(), Error> {
        match &self.filter {
            Some(filter) => {
                let processed = filter.apply(row).map_err(|reason| Error::Record { row: row.to_owned(), reason })?;
                match processed {
                    Some(processed) if !processed.trim().is_empty() => self.write_row(&processed),
                    _ => Ok(()),
                }
            }
            None => self.write_row(row),
        }
    }

    fn finish(self) -> Result<(u64, bool, u64), Error> {
        let Self { lines, any_record, written, out, options, .. } = self;
        out.finish().map_err(|source| Error::Io { path: options.output.clone(), source })?;
        Ok((lines, any_record, written))
    }
}

/// `print_calling_step`: an empty `pileup.vcf(.gz)` is also copied to `merge_output.vcf.gz` in
/// the same directory, so the pipeline has a final VCF. A failed copy is only logged, as
/// the original ignores `cp` failing.
fn copy_empty_pileup(output: &Path) -> Option<PathBuf> {
    let name = output.file_name()?;
    if name != "pileup.vcf" && name != "pileup.vcf.gz" {
        return None;
    }
    let dir = output.parent().unwrap_or(Path::new(""));
    let merge_output = dir.join("merge_output.vcf.gz");
    log::warn!("[WARNING] Copying pileup.vcf.gz to {}", merge_output.display());
    match fs::copy(dir.join("pileup.vcf.gz"), &merge_output) {
        Ok(_) => Some(merge_output),
        Err(e) => {
            log::warn!("cannot copy pileup.vcf.gz to {}: {e}", merge_output.display());
            None
        }
    }
}

/// Writes a complete header-only VCF or gVCF, then compresses and indexes it.
fn write_empty(options: &SortOptions, reason: EmptyReason) -> Result<Summary, Error> {
    let header = if options.output.extension().is_some_and(|extension| extension == "gvcf")
        && let Some(reference) = &options.reference
    {
        snv_call::gvcf_header(reference, &options.sample_name)?.0
    } else {
        snv_call::vcf_header(options.reference.as_deref(), options.cmd_file.as_deref(), &options.sample_name)?
    };
    fs::write(&options.output, header).map_err(|source| Error::Io { path: options.output.clone(), source })?;
    let (vcf, index) = compress_and_index(&options.output)?;
    let copied_to = copy_empty_pileup(&options.output);
    Ok(Summary { vcf, records: 0, empty: Some(reason), index, copied_to })
}

/// Merges the chunk VCFs in [`SortOptions::input_dir`] into `<output>.gz` and its `.tbi`, as
/// `SortVcf.py` (`sort_vcf_from`) does.
///
/// Records are ordered by contig (`chr1`..`chr22`, `chrX`, `chrY`, the same without `chr`, then
/// the other listed contigs in file order) and position; a later row at the same position
/// replaces an earlier one. Contigs missing from [`SortOptions::contigs`] are dropped. Header
/// lines are deduplicated in first-seen order and written once. When nothing matches, the
/// output is the header `get_header` builds from the reference.
///
/// # Errors
///
/// Fails where the original exits or raises: a missing input directory or contig file, an input
/// line it cannot parse, a malformed PAR BED or a missing `--cmd_fn` when the header needs one. I/O
/// and compression failures are errors too, but an indexing failure is reported in
/// [`Summary::index`], as the original ignores it.
pub fn sort_vcf(options: &SortOptions) -> Result<Summary, Error> {
    // The original builds the haploid decider, reading the PAR BED, before anything else.
    let par = match &options.par_regions {
        Some(path) => ParRegions::read(path)?,
        None => ParRegions::default(),
    };
    let decider = HaploidDecider::new(options.haploid, options.gender, par);
    let filter = match options.postprocess {
        Postprocess::Off => None,
        Postprocess::PileupOnly => Some(RowFilter { ref_calls: options.ref_calls, qual: options.qual, decider }),
    };

    let dir = &options.input_dir;
    if !dir.exists() {
        return Err(Error::InputDir { path: dir.clone() });
    }
    let mut files = list_inputs(dir)?;
    let prefix = options.prefix.as_deref();
    if let Some(prefix) = prefix {
        files.retain(|f| f.name.starts_with(prefix));
        if files.is_empty() {
            log::warn!("[WARNING] No vcf file found with prefix:{}/{prefix}, output empty vcf file", dir.display());
            return write_empty(options, EmptyReason::NoPrefixMatch);
        }
    }
    files.retain(|f| f.name.ends_with(options.suffix.as_str()));
    if files.is_empty() {
        // The original's message names the prefix here, not the suffix.
        log::warn!(
            "[WARNING] No vcf file found with suffix:{}/{}, output empty vcf file",
            dir.display(),
            prefix.unwrap_or("None")
        );
        return write_empty(options, EmptyReason::NoSuffixMatch);
    }

    let contigs = read_contigs(options.contigs.as_deref())?;
    let order = contig_order(&contigs);
    let format = Format::of(&options.suffix);
    let out = File::create(&options.output).map_err(|source| Error::Io { path: options.output.clone(), source })?;
    let out = BufWriter::with_capacity(1 << 20, out);
    let mut merge = Merge {
        options,
        filter,
        header: Header::default(),
        lines: 0,
        any_record: false,
        written: 0,
        format,
        gvcf_contigs: (format == Format::Gvcf && options.gvcf_contigs == GvcfContigs::Processed).then_some(contigs),
        out: if format == Format::Lz4Gvcf { Sink::Lz4(Lz4Writer::new(out)) } else { Sink::Plain(out) },
    };
    match options.layout {
        ChunkLayout::ByContig => merge.by_contig(&files, &order)?,
        ChunkLayout::Unordered => merge.unordered(&files, &order)?,
    }
    let (lines, any_record, written) = merge.finish()?;
    if format == Format::Lz4Gvcf {
        // The original returns as soon as `lz4` has written the output.
        let vcf = options.output.clone();
        return Ok(Summary { vcf, records: written, empty: None, index: IndexStatus::NotIndexed, copied_to: None });
    }

    if lines == 0 {
        log::warn!("[WARNING] No vcf file found, output empty vcf file");
        return write_empty(options, EmptyReason::NoLines);
    }
    if !any_record {
        log::warn!("[WARNING] No variant found, output empty vcf file");
        return write_empty(options, EmptyReason::NoRecords);
    }
    if format == Format::Gvcf {
        log::info!("[INFO] Need some time to compress and index GVCF file...");
    }
    let (vcf, index) = compress_and_index(&options.output)?;
    Ok(Summary { vcf, records: written, empty: None, index, copied_to: None })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| (*s).to_owned()).collect()
    }

    #[test]
    fn major_contig_ranks() {
        assert_eq!(major_contig_rank("chr1"), Some(0));
        assert_eq!(major_contig_rank("chr22"), Some(21));
        assert_eq!(major_contig_rank("chrX"), Some(22));
        assert_eq!(major_contig_rank("chrY"), Some(23));
        assert_eq!(major_contig_rank("1"), Some(24));
        assert_eq!(major_contig_rank("Y"), Some(47));
        for other in ["chr23", "chr0", "chr01", "chrM", "MT", "+1", "chr", ""] {
            assert_eq!(major_contig_rank(other), None, "{other:?}");
        }
    }

    #[test]
    fn gvcf_header_keeps_listed_contigs() {
        let rows =
            strings(&["##fileformat=VCFv4.2\n", "##contig=<ID=chr1,length=9>\n", "##contig=<ID=chr2,length=9>\n"]);
        let kept = gvcf_header_rows(rows, &strings(&["chr2"])).expect("valid rows");
        assert_eq!(kept, strings(&["##fileformat=VCFv4.2\n", "##contig=<ID=chr2,length=9>\n"]));
        assert!(gvcf_header_rows(strings(&["##contig=<chr1>\n"]), &[]).is_err());
    }

    #[test]
    fn contigs_sort_like_the_original() {
        let listed = strings(&["chrUn", "chrX", "chr2", "chrM", "chr20", "chr1", "2", "chrUn"]);
        let want = strings(&["chr1", "chr2", "chr20", "chrX", "2", "chrUn", "chrUn", "chrM"]);
        assert_eq!(contig_order(&listed), want);
    }
}

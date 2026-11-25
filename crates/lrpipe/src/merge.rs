//! `lrpipe merge`: the run's variant VCFs (SNVs and indels, SVs, CNVs, tandem repeats) as one
//! sorted, bgzipped and indexed VCF.
//!
//! Records keep their own columns, FORMAT included, and gain `INFO/SOURCE`, the analysis that
//! called them. The header holds every input's ALT, FILTER, INFO and FORMAT definitions once;
//! where two inputs define the same INFO or FORMAT ID with another Number or Type (lrpipe-cnv's
//! `GQ` is a Float, the others' an Integer), the later input's is renamed `SOURCE_ID` (`CNV_GQ`)
//! in its header line and records. The inputs' other header lines are kept prefixed with their
//! source (`##snv_cmdline=...`), and the contigs are the reference's, in its order.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;

/// The analysis an input VCF comes from, in the order its records come at one position.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Source {
    /// lrpipe snv (Clair3): SNVs and indels.
    Snv,
    /// lrpipe sv (Sniffles2): structural variants.
    Sv,
    /// lrpipe cnv (Spectre): copy number variants and loss of heterozygosity.
    Cnv,
    /// lrpipe tr: tandem repeat genotypes.
    Tr,
}

impl Source {
    /// The name in `INFO/SOURCE` and in prefixed header lines.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Snv => "snv",
            Self::Sv => "sv",
            Self::Cnv => "cnv",
            Self::Tr => "tr",
        }
    }
}

/// One input VCF.
#[derive(Clone, Debug)]
pub struct Input {
    /// The analysis it comes from.
    pub source: Source,
    /// The VCF, bgzipped or plain.
    pub path: PathBuf,
}

/// The `INFO/SOURCE` definition.
const SOURCE_INFO: &str = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"lrpipe analysis that called the \
                           record: snv (lrpipe snv, Clair3), sv (lrpipe sv, Sniffles2), cnv (lrpipe cnv, Spectre) \
                           or tr (lrpipe tr)\">";

/// Writes `inputs` as one VCF to `output` (`.vcf.gz`, with a `.tbi` beside it), for `sample`,
/// with the contigs of `reference`'s `.fai`.
///
/// # Errors
///
/// If an input cannot be read, has another sample, is not sorted in the reference's contig
/// order or names a contig the reference lacks, is the output itself, or the output cannot be
/// written.
pub fn merge(inputs: &[Input], reference: &Path, sample: &str, output: &Path) -> Result<()> {
    if !output.to_str().is_some_and(|o| o.ends_with(".gz")) {
        bail!("the merged VCF's name must end in .gz: {}", output.display());
    }
    refuse_an_input_as_output(inputs, output)?;
    let contigs = read_fai(reference)?;
    let order: HashMap<&str, usize> = contigs.iter().enumerate().map(|(i, (name, _))| (name.as_str(), i)).collect();
    let mut inputs: Vec<&Input> = inputs.iter().collect();
    inputs.sort_by_key(|i| i.source);

    let mut readers = Vec::with_capacity(inputs.len());
    let mut header = Header::default();
    for input in &inputs {
        let mut lines = open(&input.path)?.lines();
        let meta = header.read(input, &mut lines, sample)?;
        readers.push(Reader { input, lines, renames: meta, next: None, last: (0, 0), line: 0 });
    }

    // The VCF is written, compressed and indexed under names of this run's own, and only then
    // renamed over `output`: nothing beside it, an input among them, is created or removed
    // (LRUN-01), and a failed merge leaves an earlier `output` as it was.
    let staged = Staged::new(output)?;
    let plain = &staged.plain;
    let file = File::options()
        .write(true)
        .create_new(true)
        .open(plain)
        .with_context(|| format!("creating {}", plain.display()))?;
    let mut out = BufWriter::with_capacity(1 << 20, file);
    header.write(&mut out, reference, &contigs, sample)?;
    for reader in &mut readers {
        reader.advance(&order)?;
    }
    // The next record overall: the lowest (contig, position), the earlier source first.
    while let Some(i) = (0..readers.len())
        .filter(|&i| readers[i].next.is_some())
        .min_by_key(|&i| readers[i].next.as_ref().map(|(key, _)| (*key, i)))
    {
        let (_, line) = readers[i].next.take().expect("a pending record");
        readers[i].write(&mut out, &line)?;
        readers[i].advance(&order)?;
    }
    out.flush().with_context(|| format!("writing {}", plain.display()))?;
    drop(out);
    drop(readers);
    let (gz, index) = snv_merge::compress_and_index(plain)
        .with_context(|| format!("compressing and indexing {}", plain.display()))?;
    if index != snv_merge::IndexStatus::Built {
        bail!("tabix could not index the merged VCF {}", gz.display());
    }
    staged.commit(output)
}

/// Refuses an `output` that is one of the inputs, or its index, once both are made absolute
/// with links resolved.
fn refuse_an_input_as_output(inputs: &[Input], output: &Path) -> Result<()> {
    let written = [canonical(output), canonical(&with_suffix(output, ".tbi"))];
    for input in inputs {
        if written.contains(&canonical(&input.path)) {
            bail!(
                "--output {} would overwrite the --{} input {}; write the merged VCF to another file",
                output.display(),
                input.source.as_str(),
                input.path.display()
            );
        }
    }
    Ok(())
}

/// `path` absolute with links resolved, as far as it exists: a file not yet written is its
/// directory's canonical path and its name.
fn canonical(path: &Path) -> PathBuf {
    if let Ok(path) = path.canonicalize() {
        return path;
    }
    let parent = match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    };
    match (parent.canonicalize(), path.file_name()) {
        (Ok(dir), Some(name)) => dir.join(name),
        _ => path.to_path_buf(),
    }
}

/// `path` with `suffix` appended to its name.
fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// The merged VCF's files while it is written: hidden, beside the output (so renaming them is
/// atomic), and named after the output and this process, so no input's name can be one of
/// them. Those still there when it is dropped, after a failure, are removed.
struct Staged {
    /// The uncompressed VCF; compressing it removes it.
    plain: PathBuf,
    /// The bgzipped VCF.
    gz: PathBuf,
    /// Its tabix index.
    tbi: PathBuf,
}

/// The start of the names [`Staged`] gives the files of a merge into `output`, for removing a
/// killed merge's.
pub(crate) fn staging_prefix(output: &Path) -> std::ffi::OsString {
    let mut hidden = std::ffi::OsString::from(".");
    hidden.push(output.file_name().unwrap_or_default());
    hidden.push(".tmp-");
    hidden
}

impl Staged {
    fn new(output: &Path) -> Result<Self> {
        if output.file_name().is_none() {
            bail!("{} names no file", output.display());
        }
        // The process and the time: a killed merge's files, with a PID a container reuses, are
        // no obstacle.
        let since = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap_or_default();
        let mut hidden = staging_prefix(output);
        hidden.push(format!("{}-{}", std::process::id(), since.as_nanos()));
        let plain = output.with_file_name(hidden);
        let gz = with_suffix(&plain, ".gz");
        let tbi = with_suffix(&gz, ".tbi");
        for path in [&plain, &gz, &tbi] {
            if path.symlink_metadata().is_ok() {
                bail!("{} exists; remove it if no lrpipe merge is running", path.display());
            }
        }
        Ok(Self { plain, gz, tbi })
    }

    /// Renames the compressed VCF and its index to `output` and `output.tbi`. The earlier index
    /// goes first, so a merge stopped between the renames leaves no index rather than one of
    /// another VCF.
    fn commit(self, output: &Path) -> Result<()> {
        let index = with_suffix(output, ".tbi");
        match std::fs::remove_file(&index) {
            Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
                return Err(e).with_context(|| format!("removing {}", index.display()));
            }
            _ => {}
        }
        std::fs::rename(&self.gz, output)
            .with_context(|| format!("renaming {} to {}", self.gz.display(), output.display()))?;
        std::fs::rename(&self.tbi, &index)
            .with_context(|| format!("renaming {} to {}", self.tbi.display(), index.display()))?;
        Ok(())
    }
}

impl Drop for Staged {
    fn drop(&mut self) {
        for path in [&self.plain, &self.gz, &self.tbi] {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// The contigs of `reference.fai`, with their lengths, in order.
fn read_fai(reference: &Path) -> Result<Vec<(String, u64)>> {
    let fai = PathBuf::from(format!("{}.fai", reference.display()));
    let text = std::fs::read_to_string(&fai).with_context(|| format!("reading {}", fai.display()))?;
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut fields = l.split('\t');
            let name = fields.next().unwrap_or_default().to_owned();
            let length = fields.next().and_then(|f| f.parse().ok());
            length.map(|len| (name, len)).with_context(|| format!("{}: a line without a length: {l}", fai.display()))
        })
        .collect()
}

/// A VCF's lines, bgzipped or plain.
fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut magic = [0u8; 2];
    let gzipped =
        file.read(&mut magic).with_context(|| format!("reading {}", path.display()))? == 2 && magic == [0x1f, 0x8b];
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(if gzipped {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    })
}

/// An INFO or FORMAT ID renamed for one input.
type Renames = HashMap<(Kind, String), String>;

/// The kinds of header definitions.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Kind {
    Filter,
    Alt,
    Info,
    Format,
}

impl Kind {
    const ALL: [Kind; 4] = [Kind::Filter, Kind::Alt, Kind::Info, Kind::Format];

    fn of(line: &str) -> Option<Self> {
        Some(match line.split_once('=')?.0 {
            "##FILTER" => Kind::Filter,
            "##ALT" => Kind::Alt,
            "##INFO" => Kind::Info,
            "##FORMAT" => Kind::Format,
            _ => return None,
        })
    }
}

/// A structured header line's field: `ID=X`'s `X`.
fn field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let body = line.split_once("=<")?.1;
    let start = body.find(&format!("{key}=")).filter(|&p| p == 0 || body[..p].ends_with(','))?;
    let value = &body[start + key.len() + 1..];
    Some(value.split([',', '>']).next().unwrap_or(value))
}

/// The merged header, built input by input.
#[derive(Default)]
struct Header {
    /// Definitions by kind, in the order first seen.
    definitions: Vec<(Kind, String)>,
    /// (kind, ID) to (Number, Type), for spotting conflicts.
    seen: HashMap<(Kind, String), (String, String)>,
    /// The inputs' other header lines, prefixed with their source.
    other: Vec<String>,
}

impl Header {
    /// Reads one input's header up to its `#CHROM` line, adding its definitions, and returns
    /// the IDs renamed for it.
    fn read(&mut self, input: &Input, lines: &mut io::Lines<Box<dyn BufRead>>, sample: &str) -> Result<Renames> {
        let path = input.path.display();
        let source = input.source.as_str();
        let mut renames = Renames::new();
        for line in lines.by_ref() {
            let line = line.with_context(|| format!("reading {path}"))?;
            if let Some(columns) = line.strip_prefix("#CHROM") {
                let samples: Vec<&str> = columns.split('\t').skip(9).collect();
                if samples != [sample] {
                    bail!("{path} has samples {samples:?}, not {sample}");
                }
                return Ok(renames);
            }
            let Some(meta) = line.strip_prefix("##") else {
                bail!("{path}: a header line that does not start with ##: {line}");
            };
            if let Some(kind) = Kind::of(&line) {
                let id = field(&line, "ID").with_context(|| format!("{path}: no ID in {line}"))?.to_owned();
                let shape =
                    (field(&line, "Number").unwrap_or("").to_owned(), field(&line, "Type").unwrap_or("").to_owned());
                match self.seen.get(&(kind, id.clone())) {
                    None => {
                        self.seen.insert((kind, id.clone()), shape);
                        self.definitions.push((kind, line));
                    }
                    // FILTER and ALT have no shape; INFO and FORMAT keep the first of equal ones.
                    Some(first) if *first == shape || matches!(kind, Kind::Filter | Kind::Alt) => {}
                    Some(_) => {
                        let renamed = format!("{}_{id}", source.to_uppercase());
                        if self.seen.contains_key(&(kind, renamed.clone())) {
                            bail!("{path}: {id} conflicts with an earlier definition, and so does {renamed}");
                        }
                        self.seen.insert((kind, renamed.clone()), shape);
                        self.definitions.push((kind, line.replacen(&format!("ID={id}"), &format!("ID={renamed}"), 1)));
                        renames.insert((kind, id), renamed);
                    }
                }
            } else if !(meta.starts_with("fileformat=")
                || meta.starts_with("contig=")
                || meta.starts_with("reference="))
            {
                self.other.push(format!("##{source}_{meta}"));
            }
        }
        bail!("{path} ends before its #CHROM line")
    }

    fn write(&self, out: &mut impl Write, reference: &Path, contigs: &[(String, u64)], sample: &str) -> Result<()> {
        writeln!(out, "##fileformat=VCFv4.2")?;
        writeln!(out, "##source=lrpipe merge {}", env!("CARGO_PKG_VERSION"))?;
        writeln!(out, "##reference={}", reference.display())?;
        for (name, length) in contigs {
            writeln!(out, "##contig=<ID={name},length={length}>")?;
        }
        for kind in Kind::ALL {
            for (_, line) in self.definitions.iter().filter(|(k, _)| *k == kind) {
                writeln!(out, "{line}")?;
            }
            if kind == Kind::Info {
                writeln!(out, "{SOURCE_INFO}")?;
            }
        }
        for line in &self.other {
            writeln!(out, "{line}")?;
        }
        writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")?;
        Ok(())
    }
}

/// One input's records, one ahead.
struct Reader<'a> {
    input: &'a Input,
    lines: io::Lines<Box<dyn BufRead>>,
    renames: Renames,
    /// The next record, with its (contig index, position).
    next: Option<((usize, u64), String)>,
    /// The last record's (contig index, position), which the next may not precede.
    last: (usize, u64),
    /// Records read, for messages.
    line: usize,
}

impl Reader<'_> {
    fn advance(&mut self, order: &HashMap<&str, usize>) -> Result<()> {
        let path = self.input.path.display();
        let Some(line) = self.lines.next() else {
            return Ok(());
        };
        let line = line.with_context(|| format!("reading {path}"))?;
        self.line += 1;
        let mut fields = line.splitn(3, '\t');
        let (contig, position) = (fields.next().unwrap_or_default(), fields.next().unwrap_or_default());
        let Some(&contig) = order.get(contig) else {
            bail!("{path}: record {} is on {contig}, which the reference lacks", self.line);
        };
        let position: u64 =
            position.parse().with_context(|| format!("{path}: record {} has no position", self.line))?;
        let key = (contig, position);
        if key < self.last {
            bail!("{path} is not sorted in the reference's contig order (record {})", self.line);
        }
        self.last = key;
        self.next = Some((key, line));
        Ok(())
    }

    /// Writes `line` with `INFO/SOURCE` and its renamed INFO and FORMAT IDs.
    fn write(&self, out: &mut impl Write, line: &str) -> Result<()> {
        let source = self.input.source.as_str();
        let mut fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            bail!("{}: a record with {} columns", self.input.path.display(), fields.len());
        }
        let info = match fields[7] {
            "." | "" => format!("SOURCE={source}"),
            info if self.renames.is_empty() => format!("{info};SOURCE={source}"),
            info => {
                let renamed: Vec<String> = info.split(';').map(|kv| self.rename(Kind::Info, kv, '=')).collect();
                format!("{};SOURCE={source}", renamed.join(";"))
            }
        };
        let format = (fields.len() > 8 && !self.renames.is_empty())
            .then(|| fields[8].split(':').map(|key| self.rename(Kind::Format, key, ':')).collect::<Vec<_>>().join(":"));
        fields[7] = &info;
        if let Some(format) = &format {
            fields[8] = format;
        }
        writeln!(out, "{}", fields.join("\t"))?;
        Ok(())
    }

    /// `item` (an INFO `KEY=VALUE` or a FORMAT key) with its key renamed if it was.
    fn rename(&self, kind: Kind, item: &str, separator: char) -> String {
        let (key, rest) = match item.split_once(separator) {
            Some((key, _)) if kind == Kind::Info => (key, &item[key.len()..]),
            _ => (item, ""),
        };
        match self.renames.get(&(kind, key.to_owned())) {
            Some(renamed) => format!("{renamed}{rest}"),
            None => item.to_owned(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, text).expect("write a test file");
        path
    }

    #[test]
    fn merges_sorted_inputs_renaming_conflicting_definitions() {
        let dir = std::env::temp_dir().join(format!("lrpipe-merge-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("a test directory");
        let reference = write(&dir, "ref.fa", ">chr1\nACGT\n>chr2\nACGT\n");
        write(&dir, "ref.fa.fai", "chr1\t1000\t6\t4\t5\nchr2\t2000\t17\t4\t5\n");
        let snv = write(
            &dir,
            "snv.vcf",
            "##fileformat=VCFv4.2\n##source=Clair3\n##contig=<ID=chr1,length=1000>\n\
             ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
             ##INFO=<ID=P,Number=0,Type=Flag,Description=\"Pileup\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr1\t100\t.\tA\tG\t30\tPASS\tP\tGT:GQ\t0|1:30\n\
             chr2\t50\t.\tC\tT\t20\tPASS\t.\tGT:GQ\t1/1:20\n",
        );
        let cnv = write(
            &dir,
            "cnv.vcf",
            "##fileformat=VCFv4.2\n##source=lrpipe-cnv\n\
             ##FILTER=<ID=PASS,Description=\"All passed\">\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr1\t100\tc1\tN\t<DEL>\t10\tPASS\tEND=500\tGT:GQ\t0/1:5.5\n\
             chr1\t600\tc2\tN\t<DUP>\t10\tPASS\tEND=900\tGT:GQ\t0/1:7.5\n",
        );
        let out = dir.join("S.vcf.gz");
        let inputs = [Input { source: Source::Cnv, path: cnv }, Input { source: Source::Snv, path: snv }];
        merge(&inputs, &reference, "S", &out).expect("merged");
        let mut text = String::new();
        MultiGzDecoder::new(File::open(&out).expect("output")).read_to_string(&mut text).expect("gzip");
        let header: Vec<&str> = text.lines().filter(|l| l.starts_with("##")).collect();
        assert!(header.contains(&"##contig=<ID=chr2,length=2000>"));
        assert!(header.contains(&"##FORMAT=<ID=CNV_GQ,Number=1,Type=Float,Description=\"Genotype quality\">"));
        assert!(header.contains(&"##snv_source=Clair3"));
        assert_eq!(header.iter().filter(|l| l.starts_with("##FILTER=<ID=PASS")).count(), 1);
        let records: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            records,
            [
                "chr1\t100\t.\tA\tG\t30\tPASS\tP;SOURCE=snv\tGT:GQ\t0|1:30",
                "chr1\t100\tc1\tN\t<DEL>\t10\tPASS\tEND=500;SOURCE=cnv\tGT:CNV_GQ\t0/1:5.5",
                "chr1\t600\tc2\tN\t<DUP>\t10\tPASS\tEND=900;SOURCE=cnv\tGT:CNV_GQ\t0/1:7.5",
                "chr2\t50\t.\tC\tT\t20\tPASS\tSOURCE=snv\tGT:GQ\t1/1:20",
            ]
        );
        assert!(dir.join("S.vcf.gz.tbi").exists());

        // Another sample, or records out of order, are refused.
        let other = write(
            &dir,
            "other.vcf",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n",
        );
        let err = merge(&[Input { source: Source::Sv, path: other }], &reference, "S", &out).expect_err("other sample");
        assert!(err.to_string().contains("not S"), "{err}");
        let unsorted = write(
            &dir,
            "unsorted.vcf",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr2\t5\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\nchr1\t5\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\n",
        );
        let err = merge(&[Input { source: Source::Sv, path: unsorted }], &reference, "S", &out).expect_err("unsorted");
        assert!(err.to_string().contains("not sorted"), "{err}");
        std::fs::remove_dir_all(&dir).expect("remove the test directory");
    }

    /// A VCF of `records` SNVs of sample S on chr1, more than the 1 MiB an input's reader
    /// buffers when there are many.
    fn calls(records: usize) -> String {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        );
        for i in 0..records {
            text.push_str(&format!("chr1\t{}\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n", 1 + i));
        }
        text
    }

    fn records(path: &Path) -> usize {
        let mut text = String::new();
        MultiGzDecoder::new(File::open(path).expect("output")).read_to_string(&mut text).expect("gzip");
        text.lines().filter(|l| !l.starts_with('#')).count()
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("lrpipe-merge-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("a test directory");
        dir
    }

    fn listing(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .expect("listing")
            .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn compressing_calls_beside_them_keeps_every_record_and_the_input() {
        // LRUN-01: `--snv x.vcf --output x.vcf.gz` wrote the uncompressed merge to x.vcf, its own
        // input, then bgzip removed it: exit 0, input gone, output truncated.
        let dir = scratch("inplace");
        let reference = write(&dir, "ref.fa", ">chr1\nACGT\n");
        write(&dir, "ref.fa.fai", "chr1\t1000000\t6\t4\t5\n");
        let text = calls(100_000);
        let input = write(&dir, "x.vcf", &text);
        let out = dir.join("x.vcf.gz");
        merge(&[Input { source: Source::Snv, path: input.clone() }], &reference, "S", &out).expect("merged");
        assert_eq!(std::fs::read_to_string(&input).expect("the input is still there"), text);
        assert_eq!(records(&out), 100_000);
        assert_eq!(listing(&dir), ["ref.fa", "ref.fa.fai", "x.vcf", "x.vcf.gz", "x.vcf.gz.tbi"], "nothing else left");
        std::fs::remove_dir_all(&dir).expect("remove the test directory");
    }

    #[test]
    fn a_file_named_like_the_uncompressed_output_survives() {
        // LRUN-01 in `lrpipe run`: OUTPUT/S.vcf, a user's earlier calls, was overwritten and
        // removed by the merge stage writing OUTPUT/S.vcf.gz.
        let dir = scratch("sibling");
        let reference = write(&dir, "ref.fa", ">chr1\nACGT\n");
        write(&dir, "ref.fa.fai", "chr1\t1000\t6\t4\t5\n");
        let snv = write(&dir, "snv.vcf", &calls(3));
        write(&dir, "S.vcf", "the user's own calls");
        merge(&[Input { source: Source::Snv, path: snv }], &reference, "S", &dir.join("S.vcf.gz")).expect("merged");
        assert_eq!(std::fs::read_to_string(dir.join("S.vcf")).expect("still there"), "the user's own calls");
        assert_eq!(records(&dir.join("S.vcf.gz")), 3);
        std::fs::remove_dir_all(&dir).expect("remove the test directory");
    }

    #[test]
    fn an_input_given_as_the_output_is_refused_untouched() {
        let dir = scratch("same");
        let reference = write(&dir, "ref.fa", ">chr1\nACGT\n");
        write(&dir, "ref.fa.fai", "chr1\t1000\t6\t4\t5\n");
        let gz = dir.join("x.vcf.gz");
        std::fs::write(&gz, b"\x1f\x8b not really").expect("input");
        // The same file by another spelling.
        let spelled = dir.join(".").join("x.vcf.gz");
        let err = merge(&[Input { source: Source::Snv, path: spelled }], &reference, "S", &gz).expect_err("refused");
        assert!(err.to_string().contains("would overwrite the --snv input"), "{err}");
        assert_eq!(std::fs::read(&gz).expect("input"), b"\x1f\x8b not really");
        assert_eq!(listing(&dir), ["ref.fa", "ref.fa.fai", "x.vcf.gz"]);
        std::fs::remove_dir_all(&dir).expect("remove the test directory");
    }

    #[test]
    fn a_failed_merge_keeps_the_earlier_output_and_leaves_nothing() {
        let dir = scratch("failed");
        let reference = write(&dir, "ref.fa", ">chr1\nACGT\n");
        write(&dir, "ref.fa.fai", "chr1\t1000\t6\t4\t5\n");
        write(&dir, "S.vcf.gz", "the earlier merge");
        let unsorted = write(
            &dir,
            "unsorted.vcf",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr1\t9\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\nchr1\t5\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\n",
        );
        let out = dir.join("S.vcf.gz");
        merge(&[Input { source: Source::Sv, path: unsorted }], &reference, "S", &out).expect_err("unsorted");
        assert_eq!(std::fs::read_to_string(&out).expect("kept"), "the earlier merge");
        assert_eq!(listing(&dir), ["S.vcf.gz", "ref.fa", "ref.fa.fai", "unsorted.vcf"]);
        std::fs::remove_dir_all(&dir).expect("remove the test directory");
    }
}

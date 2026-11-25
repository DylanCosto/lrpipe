//! Compares two runs' outputs record by record, and their tools, times and peak memory: for
//! checking that a change to a tool (or a library under all of them) leaves every result as
//! it was, and what it did to speed.
//!
//! Lines starting with `#` are headers; the rest are records, compared in order and, where
//! they differ, also as multisets (how many records only one run has). Header lines that
//! change with every run (dates, command lines, and the same in the merged VCF, prefixed with
//! their source: `##snv_cmdline=`) are left out; the others are compared as sets and reported
//! apart from the records.
//!
//! A stage a run skipped (already done) keeps the measurements of the run that did it, which
//! are shown but not compared as the run's: no change in per cent is given for it, nor for the
//! run when it ran no stage.
//!
//! An `.snf` (lrpipe-sv's SV candidates) is a line of JSON and compressed blocks: its header
//! lines are the JSON's settings but for the command line, and its records its blocks, each
//! compared by its contig, start and contents once decompressed.
//!
//! A BAM (the aligned reads, the merged reads, the haplotagged reads) is compared by its records,
//! each as the bytes BAM stores it (every field and tag), those at the same place as a multiset
//! (their order there is not the sort's to keep), and its header lines, but for the command lines
//! of its programs (`@PG CL:`). Both BAMs are read once, side by side. A BAM without its BGZF
//! end-of-file block (cut short, perhaps at a record's end) is an error. A file named as text
//! (`.vcf`, `.bed.gz`, ...) that is not text is an error; any other file that is not text is
//! compared byte for byte, as one record.

use std::collections::BTreeMap;
use std::collections::hash_map::DefaultHasher;
use std::fmt;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::thread;

use flate2::read::MultiGzDecoder;
use serde::Deserialize;

use crate::provenance::ToolRecord;

/// Header lines that differ between runs of the same command: these keys, alone or after a
/// source's prefix as the merged VCF has them (`##sv_fileDate=`).
pub const VOLATILE_HEADERS: [&str; 3] = ["fileDate", "command", "cmdline"];

/// Whether the header `line` is one of [`VOLATILE_HEADERS`].
fn volatile(line: &str) -> bool {
    let Some((key, _)) = line.strip_prefix("##").and_then(|l| l.split_once('=')) else { return false };
    VOLATILE_HEADERS.iter().any(|v| key == *v || key.strip_suffix(v).is_some_and(|source| source.ends_with('_')))
}

/// Why a run could not be compared.
#[derive(Debug, thiserror::Error)]
pub enum CompareError {
    /// A file could not be read.
    #[error("reading {}", path.display())]
    Read {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A manifest is not lrpipe's.
    #[error("parsing {}", path.display())]
    Parse {
        /// The manifest.
        path: PathBuf,
        /// The cause.
        #[source]
        source: serde_json::Error,
    },
}

/// The parts of a manifest a comparison needs.
#[derive(Debug, Deserialize)]
struct RunManifest {
    #[serde(default)]
    output: Option<PathBuf>,
    #[serde(default)]
    wall_seconds: Option<f64>,
    stages: Vec<StageMeasurements>,
    tools: BTreeMap<String, ToolRecord>,
    outputs: BTreeMap<String, PathBuf>,
}

/// A stage's measurements as a manifest records them, by name, so that runs of lrpipe
/// versions with different stages compare.
#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct StageMeasurements {
    /// The stage's name.
    pub name: String,
    /// Wall-clock seconds.
    #[serde(default)]
    pub wall_seconds: Option<f64>,
    /// Peak resident memory in GiB.
    #[serde(default)]
    pub max_rss_gib: Option<f64>,
    /// How it ended (`succeeded`, `skipped`, ...).
    #[serde(default)]
    pub status: Option<String>,
}

impl StageMeasurements {
    /// Whether the run skipped it, keeping an earlier run's measurements.
    fn skipped(&self) -> bool {
        self.status.as_deref() == Some("skipped")
    }
}

/// A finished (or stopped) run, read from its output directory.
#[derive(Debug)]
pub struct Run {
    dir: PathBuf,
    manifest: RunManifest,
}

impl Run {
    /// Reads the manifest in `dir`.
    ///
    /// # Errors
    ///
    /// If `dir` has no readable lrpipe manifest.
    pub fn open(dir: impl Into<PathBuf>) -> Result<Self, CompareError> {
        let dir = dir.into();
        let path = dir.join("manifest.json");
        let text =
            std::fs::read_to_string(&path).map_err(|source| CompareError::Read { path: path.clone(), source })?;
        let manifest = serde_json::from_str(&text).map_err(|source| CompareError::Parse { path, source })?;
        Ok(Self { dir, manifest })
    }

    /// Where the run's output of `kind` is now: under `dir` even if the run was moved there.
    fn output(&self, kind: &str) -> Option<PathBuf> {
        let recorded = self.manifest.outputs.get(kind)?;
        let moved = self.manifest.output.as_ref().and_then(|root| recorded.strip_prefix(root).ok());
        Some(moved.map_or_else(|| recorded.clone(), |rel| self.dir.join(rel)))
    }
}

/// How two versions of one output compare.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// One run lacks the output (listed or on disk).
    Missing(Side),
    /// Both have it.
    Compared {
        /// Their records.
        records: Records,
        /// Header lines, other than [`VOLATILE_HEADERS`], only one has.
        headers: HeaderChange,
    },
}

/// Which run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    /// The first.
    Old,
    /// The second.
    New,
    /// Neither.
    Both,
}

/// How two outputs' records compare.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Records {
    /// The same records in the same order.
    Identical(u64),
    /// Not.
    Different {
        /// Records in the old output.
        old: u64,
        /// Records in the new output.
        new: u64,
        /// Records only the old output has (counting repeats).
        only_old: u64,
        /// Records only the new output has.
        only_new: u64,
        /// The first record where they part.
        first: Box<FirstDifference>,
    },
}

/// The first record where two outputs part.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FirstDifference {
    /// Its number, from 1.
    pub record: u64,
    /// The old output's record there, if it has one.
    pub old: Option<String>,
    /// The new output's.
    pub new: Option<String>,
}

/// Header lines only one output has.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HeaderChange {
    /// Only in the old.
    pub only_old: Vec<String>,
    /// Only in the new.
    pub only_new: Vec<String>,
}

impl HeaderChange {
    /// Whether the headers match.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.only_old.is_empty() && self.only_new.is_empty()
    }
}

/// One output of both runs.
#[derive(Clone, Debug)]
pub struct FileComparison {
    /// The output's kind, as the manifest names it.
    pub kind: String,
    /// The outcome.
    pub outcome: Outcome,
}

/// A tool that changed between the runs.
#[derive(Clone, Debug)]
pub struct ToolChange {
    /// Its name.
    pub name: String,
    /// The old run's.
    pub old: Option<ToolRecord>,
    /// The new run's.
    pub new: Option<ToolRecord>,
}

/// One stage of both runs.
#[derive(Clone, Debug)]
pub struct StageComparison {
    /// The stage's name.
    pub name: String,
    /// The old run's record.
    pub old: Option<StageMeasurements>,
    /// The new run's.
    pub new: Option<StageMeasurements>,
}

/// Everything two runs differ in.
#[derive(Clone, Debug)]
pub struct Comparison {
    /// Each output kind either run has.
    pub files: Vec<FileComparison>,
    /// Tools whose binary differs.
    pub tools: Vec<ToolChange>,
    /// Each stage either run has.
    pub stages: Vec<StageComparison>,
    /// The runs' wall-clock seconds.
    pub wall_seconds: (Option<f64>, Option<f64>),
}

impl Comparison {
    /// Whether every output has the same records in both runs.
    #[must_use]
    pub fn outputs_match(&self) -> bool {
        self.files.iter().all(|f| matches!(f.outcome, Outcome::Compared { records: Records::Identical(_), .. }))
    }
}

/// Compares `old` and `new`, reading each pair of outputs on its own thread.
///
/// # Errors
///
/// If an output both runs list and have cannot be read.
pub fn compare(old: &Run, new: &Run) -> Result<Comparison, CompareError> {
    let mut kinds: Vec<&String> = old.manifest.outputs.keys().chain(new.manifest.outputs.keys()).collect();
    kinds.sort();
    kinds.dedup();
    let files = thread::scope(|scope| {
        let jobs: Vec<_> = kinds
            .iter()
            .map(|&kind| {
                let pair = (old.output(kind).filter(|p| p.exists()), new.output(kind).filter(|p| p.exists()));
                scope.spawn(move || {
                    let outcome = match pair {
                        (Some(o), Some(n)) => compare_files(&o, &n)?,
                        (None, Some(_)) => Outcome::Missing(Side::Old),
                        (Some(_), None) => Outcome::Missing(Side::New),
                        (None, None) => Outcome::Missing(Side::Both),
                    };
                    Ok(FileComparison { kind: kind.clone(), outcome })
                })
            })
            .collect();
        jobs.into_iter().map(|j| j.join().expect("comparison threads do not panic")).collect::<Result<Vec<_>, _>>()
    })?;

    let mut names: Vec<&String> = old.manifest.tools.keys().chain(new.manifest.tools.keys()).collect();
    names.sort();
    names.dedup();
    let tools = names
        .into_iter()
        .filter_map(|name| {
            let (o, n) = (old.manifest.tools.get(name), new.manifest.tools.get(name));
            let same = matches!((o, n), (Some(o), Some(n)) if o.sha256 == n.sha256);
            (!same).then(|| ToolChange { name: name.clone(), old: o.cloned(), new: n.cloned() })
        })
        .collect();

    let mut stages: Vec<StageComparison> = Vec::new();
    for s in old.manifest.stages.iter().chain(&new.manifest.stages) {
        let name = s.name.as_str();
        if stages.iter().any(|c| c.name == name) {
            continue;
        }
        let find = |run: &Run| run.manifest.stages.iter().find(|r| r.name == s.name).cloned();
        stages.push(StageComparison { name: name.to_owned(), old: find(old), new: find(new) });
    }
    Ok(Comparison { files, tools, stages, wall_seconds: (old.manifest.wall_seconds, new.manifest.wall_seconds) })
}

/// Compares two versions of one output: plain text or gzip (including BGZF), an `.snf`, a BAM, or
/// any other file byte for byte.
///
/// # Errors
///
/// If either cannot be read.
pub fn compare_files(old: &Path, new: &Path) -> Result<Outcome, CompareError> {
    if is_bam(old) && is_bam(new) {
        return compare_bams(old, new);
    }
    let (a, b) = thread::scope(|scope| {
        let a = scope.spawn(|| Summary::read(old));
        let b = Summary::read(new);
        (a.join().expect("reading threads do not panic"), b)
    });
    let (a, b) = (a?, b?);
    let headers = HeaderChange {
        only_old: a.headers.iter().filter(|h| !b.headers.contains(h)).cloned().collect(),
        only_new: b.headers.iter().filter(|h| !a.headers.contains(h)).cloned().collect(),
    };
    let count = |s: &Summary| s.records.len() as u64;
    let records = match a.records.iter().zip(&b.records).position(|(x, y)| x != y) {
        None if a.records.len() == b.records.len() => Records::Identical(count(&a)),
        at => {
            let at = at.unwrap_or(a.records.len().min(b.records.len()));
            let (only_old, only_new) = multiset_difference(&a.records, &b.records);
            let first = Box::new(FirstDifference {
                record: at as u64 + 1,
                old: nth_record(old, at).map_err(|source| CompareError::Read { path: old.into(), source })?,
                new: nth_record(new, at).map_err(|source| CompareError::Read { path: new.into(), source })?,
            });
            Records::Different { old: count(&a), new: count(&b), only_old, only_new, first }
        }
    };
    Ok(Outcome::Compared { records, headers })
}

/// An output's headers, and a hash of each record in order.
#[derive(Debug)]
struct Summary {
    headers: Vec<String>,
    records: Vec<u64>,
}

impl Summary {
    fn read(path: &Path) -> Result<Self, CompareError> {
        let read = |source| CompareError::Read { path: path.into(), source };
        if is_snf(path) {
            let (headers, blocks) = snf(path).map_err(read)?;
            let records = blocks
                .iter()
                .map(|(label, data)| {
                    let mut h = DefaultHasher::new();
                    (label, data).hash(&mut h);
                    h.finish()
                })
                .collect();
            return Ok(Self { headers, records });
        }
        let mut summary = Self { headers: Vec::new(), records: Vec::new() };
        for line in lines(path).map_err(read)? {
            let line = match line {
                Ok(line) => line,
                // Not text, though its name says it is: said, not compared otherwise.
                Err(e) if e.kind() == io::ErrorKind::InvalidData && named_as_text(path) => {
                    return Err(read(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("it is not text, though its name says it is ({e})"),
                    )));
                }
                // Not text: the whole file, as one record.
                Err(e) if e.kind() == io::ErrorKind::InvalidData => {
                    return Ok(Self { headers: Vec::new(), records: vec![file_hash(path).map_err(read)?] });
                }
                Err(e) => return Err(read(e)),
            };
            if line.starts_with('#') {
                if !volatile(&line) {
                    summary.headers.push(line);
                }
            } else {
                let mut h = DefaultHasher::new();
                line.hash(&mut h);
                summary.records.push(h.finish());
            }
        }
        Ok(summary)
    }
}

fn lines(path: &Path) -> io::Result<io::Lines<BufReader<Box<dyn Read>>>> {
    let file = File::open(path)?;
    let gzip = path.extension().is_some_and(|e| e == "gz");
    let reader: Box<dyn Read> = if gzip { Box::new(MultiGzDecoder::new(file)) } else { Box::new(file) };
    Ok(BufReader::with_capacity(1 << 20, reader).lines())
}

/// Record `n` (from 0) of the output at `path`.
fn nth_record(path: &Path, n: usize) -> io::Result<Option<String>> {
    if is_snf(path) {
        let (_, mut blocks) = snf(path)?;
        return Ok((n < blocks.len()).then(|| {
            let (label, data) = blocks.swap_remove(n);
            format!("{label}: {} bytes decompressed", data.len())
        }));
    }
    let mut records = lines(path)?.filter(|l| l.as_ref().map_or(true, |l| !l.starts_with('#')));
    match records.nth(n).transpose() {
        Err(e) if e.kind() == io::ErrorKind::InvalidData => {
            let bytes = std::fs::metadata(path)?.len();
            Ok((n == 0).then(|| format!("not text: {bytes} bytes, compared byte for byte")))
        }
        other => other,
    }
}

/// A hash of the whole file at `path`.
fn file_hash(path: &Path) -> io::Result<u64> {
    let mut file = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut h = DefaultHasher::new();
    loop {
        let buf = file.fill_buf()?;
        if buf.is_empty() {
            return Ok(h.finish());
        }
        buf.hash(&mut h);
        let n = buf.len();
        file.consume(n);
    }
}

/// Whether `path` is a BAM (by its name, as lrpipe names them).
fn is_bam(path: &Path) -> bool {
    path.extension().is_some_and(|e| e == "bam")
}

/// Whether `path` is named as text, gzipped or not: a VCF, BED, bedMethyl, TSV or text file.
fn named_as_text(path: &Path) -> bool {
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    Path::new(name)
        .extension()
        .is_some_and(|e| ["vcf", "gvcf", "bed", "bedmethyl", "tsv", "txt", "csv"].iter().any(|t| e == *t))
}

/// Compares two BAMs, reading them side by side, once: their header lines, and their records
/// place by place (reference sequence and start), those at one place as a multiset. The first
/// difference is the first record where the two parts, the records at its place in hash order.
fn compare_bams(old: &Path, new: &Path) -> Result<Outcome, CompareError> {
    let err = |path: &Path| {
        let path = path.to_path_buf();
        move |source| CompareError::Read { path, source }
    };
    let (a, b) = thread::scope(|scope| {
        let a = scope.spawn(|| Bam::open(old));
        let b = Bam::open(new);
        (a.join().expect("reading threads do not panic"), b)
    });
    let (mut a, mut b) = (a.map_err(err(old))?, b.map_err(err(new))?);
    let headers = HeaderChange {
        only_old: a.headers.iter().filter(|h| !b.headers.contains(h)).cloned().collect(),
        only_new: b.headers.iter().filter(|h| !a.headers.contains(h)).cloned().collect(),
    };
    let (mut all_old, mut all_new) = (Vec::new(), Vec::new());
    let mut first: Option<FirstDifference> = None;
    loop {
        let describe = first.is_none();
        let (x, y) = thread::scope(|scope| {
            let x = scope.spawn(|| a.place(describe));
            let y = b.place(describe);
            (x.join().expect("reading threads do not panic"), y)
        });
        let (x, y) = (x.map_err(err(old))?, y.map_err(err(new))?);
        if x.is_empty() && y.is_empty() {
            break;
        }
        if first.is_none() {
            let at = x.iter().zip(&y).position(|(p, q)| p.0 != q.0);
            if let Some(at) = at.or((x.len() != y.len()).then(|| x.len().min(y.len()))) {
                first = Some(FirstDifference {
                    record: all_old.len() as u64 + at as u64 + 1,
                    old: x.get(at).and_then(|r| r.1.clone()),
                    new: y.get(at).and_then(|r| r.1.clone()),
                });
            }
        }
        all_old.extend(x.iter().map(|r| r.0));
        all_new.extend(y.iter().map(|r| r.0));
    }
    let records = match first {
        None => Records::Identical(all_old.len() as u64),
        Some(first) => {
            let (only_old, only_new) = multiset_difference(&all_old, &all_new);
            let (old, new) = (all_old.len() as u64, all_new.len() as u64);
            Records::Different { old, new, only_old, only_new, first: Box::new(first) }
        }
    };
    Ok(Outcome::Compared { records, headers })
}

/// BGZF's end-of-file block, which ends every whole BAM.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// Threads decompressing each BAM compared.
const BAM_WORKERS: usize = 4;

/// A BAM being read: its header lines, but for the `CL:` fields of its `@PG` lines, which name the
/// run's paths; its reference sequences' names; and its records, each the bytes BAM stores it as.
struct Bam {
    reader: Box<dyn Read + Send>,
    headers: Vec<String>,
    contigs: Vec<String>,
    /// A record read past the end of a place, which begins the next.
    pending: Option<Vec<u8>>,
}

impl Bam {
    fn open(path: &Path) -> io::Result<Self> {
        use std::io::{Seek, SeekFrom};
        let invalid = |what: &str| io::Error::new(io::ErrorKind::InvalidData, format!("not a BAM: {what}"));
        let mut file = File::open(path)?;
        let mut end = [0; BGZF_EOF.len()];
        let whole = file.metadata()?.len() >= end.len() as u64
            && file.seek(SeekFrom::End(-(end.len() as i64))).and_then(|_| file.read_exact(&mut end)).is_ok()
            && end == BGZF_EOF;
        if !whole {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "it is cut short: it has no BGZF end-of-file block, which every whole BAM ends with",
            ));
        }
        let workers = std::num::NonZero::new(BAM_WORKERS).expect("not zero");
        let mut reader: Box<dyn Read + Send> =
            Box::new(pasta::bgzf::io::MultithreadedReader::with_worker_count(workers, File::open(path)?));
        let mut magic = [0; 4];
        reader.read_exact(&mut magic)?;
        if magic != *b"BAM\x01" {
            return Err(invalid("no BAM magic"));
        }
        let text = read_block(&mut reader, "its header text")?;
        let headers = String::from_utf8_lossy(&text)
            .lines()
            .filter(|l| !l.is_empty())
            .map(|l| {
                if l.starts_with("@PG\t") {
                    l.split('\t').filter(|f| !f.starts_with("CL:")).collect::<Vec<_>>().join("\t")
                } else {
                    l.trim_end_matches('\0').to_owned()
                }
            })
            .filter(|l| !l.is_empty())
            .collect();
        let n = read_i32(&mut reader)?;
        let mut contigs = Vec::with_capacity(usize::try_from(n).unwrap_or(0));
        for _ in 0..n {
            let name = read_block(&mut reader, "a reference sequence's name")?;
            contigs.push(String::from_utf8_lossy(name.strip_suffix(b"\0").unwrap_or(&name)).into_owned());
            read_i32(&mut reader)?;
        }
        Ok(Self { reader, headers, contigs, pending: None })
    }

    /// The records at the next place (reference sequence and start, the first 8 bytes of each):
    /// each one's hash, and with `describe` how it is shown, in hash order; none at the end.
    fn place(&mut self, describe: bool) -> io::Result<Vec<(u64, Option<String>)>> {
        let mut out = Vec::new();
        let mut key: Option<[u8; 8]> = None;
        loop {
            let record = match self.pending.take() {
                Some(record) => record,
                None => {
                    let mut record = Vec::new();
                    if !self.next(&mut record)? {
                        break;
                    }
                    record
                }
            };
            let here: Option<[u8; 8]> = record.get(..8).and_then(|k| k.try_into().ok());
            if !out.is_empty() && (here.is_none() || here != key) {
                self.pending = Some(record);
                break;
            }
            key = here;
            let mut h = DefaultHasher::new();
            record.hash(&mut h);
            out.push((h.finish(), describe.then(|| self.describe(&record))));
        }
        out.sort_unstable_by_key(|r| r.0);
        Ok(out)
    }

    /// Reads the next record into `record`; false at the end.
    fn next(&mut self, record: &mut Vec<u8>) -> io::Result<bool> {
        let mut size = [0; 4];
        match self.reader.read_exact(&mut size) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
            Err(e) => return Err(e),
        }
        let size = usize::try_from(i32::from_le_bytes(size))
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "a BAM record of negative size"))?;
        record.resize(size, 0);
        self.reader.read_exact(record)?;
        Ok(true)
    }

    /// A record's name, place, flags and mapping quality, as shown where two BAMs part.
    fn describe(&self, record: &[u8]) -> String {
        if record.len() < 32 {
            return format!("a record of {} bytes", record.len());
        }
        let i32_at = |at: usize| i32::from_le_bytes(record[at..at + 4].try_into().expect("4 bytes"));
        let name_len = usize::from(record[8]);
        let name = record.get(32..32 + name_len.saturating_sub(1)).map(String::from_utf8_lossy).unwrap_or_default();
        let flags = u16::from_le_bytes([record[14], record[15]]);
        let place = usize::try_from(i32_at(0))
            .ok()
            .and_then(|i| self.contigs.get(i))
            .map_or_else(|| "unplaced".to_owned(), |contig| format!("{contig}:{}", i64::from(i32_at(4)) + 1));
        format!("{name} at {place}, flags {flags}, MAPQ {}, {} bytes", record[9], record.len())
    }
}

/// A little-endian `i32`.
fn read_i32(reader: &mut impl Read) -> io::Result<i32> {
    let mut n = [0; 4];
    reader.read_exact(&mut n)?;
    Ok(i32::from_le_bytes(n))
}

/// A BAM header's field of its length (an `i32`) and that many bytes.
fn read_block(reader: &mut impl Read, what: &str) -> io::Result<Vec<u8>> {
    let len = usize::try_from(read_i32(reader)?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("not a BAM: {what} of negative length")))?;
    let mut bytes = Vec::new();
    reader.take(len as u64).read_to_end(&mut bytes)?;
    if bytes.len() != len {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, format!("a BAM cut off in {what}")));
    }
    Ok(bytes)
}

/// Whether `path` is an `.snf` of SV candidates.
fn is_snf(path: &Path) -> bool {
    path.extension().is_some_and(|e| e == "snf")
}

/// A block of an `.snf`: its contig and start (`chr1 block 5000000`), and its contents.
type SnfBlock = (String, Vec<u8>);

/// An `.snf`'s settings as header lines (`key=value`), but for the command line, and its blocks in
/// file order, decompressed.
fn snf(path: &Path) -> io::Result<(Vec<String>, Vec<SnfBlock>)> {
    let invalid = |what: &str| io::Error::new(io::ErrorKind::InvalidData, format!("not an .snf: {what}"));
    let mut file = BufReader::new(File::open(path)?);
    let mut line = String::new();
    file.read_line(&mut line)?;
    // Sniffles2 writes Python's NaN and infinities, which JSON does not allow (as lrpipe-sv reads them).
    let text = line.replace(": NaN", ": null").replace(": Infinity", ": null").replace(": -Infinity", ": null");
    let header: serde_json::Value = serde_json::from_str(&text).map_err(|_| invalid("its header is not JSON"))?;
    let config = header["config"].as_object().ok_or_else(|| invalid("no config"))?;
    let mut headers: Vec<String> =
        config.iter().filter(|(k, _)| *k != "command").map(|(k, v)| format!("{k}={v}")).collect();
    headers.push(format!("snf_candidate_count={}", header["snf_candidate_count"]));
    let mut parts: Vec<(u64, u64, String)> = Vec::new();
    for (contig, blocks) in header["index"].as_object().ok_or_else(|| invalid("no index"))? {
        for (start, block) in blocks.as_object().ok_or_else(|| invalid("a bad index"))? {
            for part in block.as_array().ok_or_else(|| invalid("a bad index"))? {
                let (Some(offset), Some(len)) = (part[0].as_u64(), part[1].as_u64()) else {
                    return Err(invalid("a bad index"));
                };
                parts.push((offset, len, format!("{contig} block {start}")));
            }
        }
    }
    parts.sort_unstable();
    let mut blocks = Vec::with_capacity(parts.len());
    let mut at = 0;
    for (offset, len, label) in parts {
        if offset != at {
            return Err(invalid("its blocks overlap or leave gaps"));
        }
        let mut compressed = vec![0; usize::try_from(len).map_err(|_| invalid("a block too large"))?];
        file.read_exact(&mut compressed)?;
        let mut data = Vec::new();
        MultiGzDecoder::new(&compressed[..]).read_to_end(&mut data)?;
        blocks.push((label, data));
        at = offset + len;
    }
    Ok((headers, blocks))
}

/// How many of `a`'s hashes `b` lacks, and of `b`'s that `a` lacks, counting repeats.
fn multiset_difference(a: &[u64], b: &[u64]) -> (u64, u64) {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    a.sort_unstable();
    b.sort_unstable();
    let (mut i, mut j, mut only_a, mut only_b) = (0, 0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => (only_a, i) = (only_a + 1, i + 1),
            std::cmp::Ordering::Greater => (only_b, j) = (only_b + 1, j + 1),
            std::cmp::Ordering::Equal => (i, j) = (i + 1, j + 1),
        }
    }
    (only_a + (a.len() - i) as u64, only_b + (b.len() - j) as u64)
}

impl fmt::Display for Comparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "outputs")?;
        for file in &self.files {
            write!(f, "  {:<24} ", file.kind)?;
            match &file.outcome {
                Outcome::Missing(Side::Old) => writeln!(f, "only in the new run")?,
                Outcome::Missing(Side::New) => writeln!(f, "only in the old run")?,
                Outcome::Missing(Side::Both) => writeln!(f, "in neither run")?,
                Outcome::Compared { records, headers } => {
                    match records {
                        Records::Identical(n) => {
                            write!(f, "identical ({n} record{})", if *n == 1 { "" } else { "s" })?;
                        }
                        Records::Different { old, new, only_old, only_new, first } => {
                            write!(f, "DIFFERENT: {old} -> {new} records, {only_old} only old, {only_new} only new")?;
                            writeln!(f)?;
                            writeln!(f, "      first at record {}:", first.record)?;
                            let (o, n) =
                                (first.old.as_deref().unwrap_or("(none)"), first.new.as_deref().unwrap_or("(none)"));
                            let from = context_start(o, n);
                            writeln!(f, "      - {}", clip_from(o, from))?;
                            write!(f, "      + {}", clip_from(n, from))?;
                        }
                    }
                    if !headers.is_empty() {
                        write!(f, "; headers differ")?;
                        for h in &headers.only_old {
                            write!(f, "\n      - {}", clip(h))?;
                        }
                        for h in &headers.only_new {
                            write!(f, "\n      + {}", clip(h))?;
                        }
                    }
                    writeln!(f)?;
                }
            }
        }
        if !self.tools.is_empty() {
            writeln!(f, "tools changed")?;
            for t in &self.tools {
                let describe = |r: &Option<ToolRecord>| {
                    r.as_ref()
                        .map_or("-".to_owned(), |r| format!("{} ({})", r.version, &r.sha256[..12.min(r.sha256.len())]))
                };
                writeln!(f, "  {:<24} {} -> {}", t.name, describe(&t.old), describe(&t.new))?;
            }
        }
        writeln!(f, "stages{:16}{:^25}  {:^25}", "", "wall seconds", "peak GiB")?;
        let skipped = |r: &Option<StageMeasurements>| r.as_ref().is_some_and(StageMeasurements::skipped);
        for s in &self.stages {
            let wall = |r: &Option<StageMeasurements>| r.as_ref().and_then(|r| r.wall_seconds);
            let rss = |r: &Option<StageMeasurements>| r.as_ref().and_then(|r| r.max_rss_gib);
            let ran = !skipped(&s.old) && !skipped(&s.new);
            writeln!(
                f,
                "  {:<20}{}  {}{}",
                s.name,
                change(wall(&s.old), wall(&s.new), 0, ran),
                change(rss(&s.old), rss(&s.new), 1, ran),
                match (skipped(&s.old), skipped(&s.new)) {
                    (true, true) => "  (skipped in both)",
                    (true, false) => "  (skipped in the old run)",
                    (false, true) => "  (skipped in the new run)",
                    (false, false) => "",
                }
            )?;
        }
        // A run that ran no stage (a resume with nothing to do) took no time worth comparing.
        let idle = |side: fn(&StageComparison) -> &Option<StageMeasurements>| {
            self.stages.iter().all(|s| side(s).is_none() || skipped(side(s)))
        };
        let ran = !idle(|s| &s.old) && !idle(|s| &s.new);
        writeln!(f, "  {:<20}{}", "run", change(self.wall_seconds.0, self.wall_seconds.1, 0, ran))
    }
}

/// `old -> new (±x%)`, aligned; the change in per cent only when `compared` and the old value
/// shows as more than 0 at `decimals` (`0 -> 7` is no +2404%).
fn change(old: Option<f64>, new: Option<f64>, decimals: usize, compared: bool) -> String {
    let show = |x: Option<f64>| x.map_or("-".to_owned(), |x| format!("{x:.decimals$}"));
    let shown = 0.5 * 10f64.powi(-i32::try_from(decimals).unwrap_or(0));
    let pct = match (old, new) {
        (Some(o), Some(n)) if compared && o >= shown => format!("{:+.0}%", (n / o - 1.0) * 100.0),
        _ => String::new(),
    };
    format!("{:>7} -> {:<7}{pct:>7}", show(old), show(new))
}

/// Characters of a line shown.
const SHOWN: usize = 160;

fn clip(s: &str) -> String {
    clip_from(s, 0)
}

/// Up to [`SHOWN`] characters of `s` from character `from`, marking what is left out.
fn clip_from(s: &str, from: usize) -> String {
    let start = s.char_indices().nth(from).map_or(s.len(), |(i, _)| i);
    let rest = &s[start..];
    let head = if from > 0 { "..." } else { "" };
    match rest.char_indices().nth(SHOWN) {
        Some((i, _)) => format!("{head}{}...", &rest[..i]),
        None => format!("{head}{rest}"),
    }
}

/// Where to start showing two lines so that their first difference is in view, a quarter of
/// the way in.
fn context_start(a: &str, b: &str) -> usize {
    let differs_at = a.chars().zip(b.chars()).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len()));
    differs_at.saturating_sub(SHOWN / 4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("lrpipe-compare-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("temp dir");
        dir
    }

    fn gz(path: &Path, text: &str) {
        let mut e = flate2::write::GzEncoder::new(File::create(path).expect("created"), flate2::Compression::fast());
        e.write_all(text.as_bytes()).expect("written");
        e.finish().expect("finished");
    }

    #[test]
    fn volatile_headers_are_ignored_and_others_reported() {
        let dir = temp_dir("headers");
        let (a, b) = (dir.join("a.vcf.gz"), dir.join("b.vcf"));
        gz(&a, "##fileformat=VCFv4.2\n##fileDate=1\n##source=x_1\n#CHROM\nchr1\t1\n");
        std::fs::write(&b, "##fileformat=VCFv4.2\n##fileDate=2\n##command=y\n##source=x_2\n#CHROM\nchr1\t1\n")
            .expect("written");
        let outcome = compare_files(&a, &b).expect("compared");
        let Outcome::Compared { records, headers } = outcome else { panic!("both exist") };
        assert_eq!(records, Records::Identical(1));
        assert_eq!(headers.only_old, ["##source=x_1"]);
        assert_eq!(headers.only_new, ["##source=x_2"]);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn differences_are_located_and_counted() {
        let dir = temp_dir("records");
        let (a, b) = (dir.join("a.tsv"), dir.join("b.tsv"));
        std::fs::write(&a, "#h\nr1\nr2\nr3\nr4\n").expect("written");
        std::fs::write(&b, "#h\nr1\nr3\nr2\nr5\nr6\n").expect("written");
        let Outcome::Compared { records, .. } = compare_files(&a, &b).expect("compared") else { panic!("both exist") };
        let Records::Different { old, new, only_old, only_new, first } = records else { panic!("they differ") };
        assert_eq!((old, new, only_old, only_new), (4, 5, 1, 2), "r4 only old; r5 and r6 only new");
        assert_eq!(first.record, 2);
        assert_eq!((first.old.as_deref(), first.new.as_deref()), (Some("r2"), Some("r3")));

        // A prefix: the first difference is past the shorter one's end.
        std::fs::write(&b, "#h\nr1\nr2\n").expect("written");
        let Outcome::Compared { records: Records::Different { first, .. }, .. } =
            compare_files(&a, &b).expect("compared")
        else {
            panic!("they differ")
        };
        assert_eq!((first.record, first.old.as_deref(), first.new.as_deref()), (3, Some("r3"), None));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn long_lines_are_shown_around_their_difference() {
        let a = format!("{}A{}", "x".repeat(300), "y".repeat(300));
        let b = format!("{}B{}", "x".repeat(300), "y".repeat(300));
        let from = context_start(&a, &b);
        assert!(clip_from(&a, from).contains('A') && clip_from(&b, from).contains('B'));
        assert!(clip_from(&a, from).starts_with("...") && clip_from(&a, from).ends_with("..."));
        assert_eq!(clip_from("short", 0), "short");
    }

    /// An `.snf` with `blocks` (contig, start, contents) and `command` in its header.
    fn snf(path: &Path, command: &str, blocks: &[(&str, i64, &str)]) {
        let mut compressed = Vec::new();
        let mut index = serde_json::Map::new();
        for (contig, start, text) in blocks {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            e.write_all(text.as_bytes()).expect("written");
            let gz = e.finish().expect("finished");
            let entry = index.entry(contig.to_string()).or_insert_with(|| serde_json::json!({}));
            entry[start.to_string()] = serde_json::json!([[compressed.len(), gz.len()]]);
            compressed.extend(gz);
        }
        let header = serde_json::json!({
            "config": {"command": command, "sample_id": "S", "contig_coverages": {"chrN": "NAN", "chrI": "INF"}},
            "index": index,
            "snf_candidate_count": blocks.len(),
        });
        // As Sniffles2's json module writes a contig's coverage that is not a number.
        let text = serde_json::to_string(&header).expect("JSON");
        let text = text.replace(r#":"NAN""#, ": NaN").replace(r#":"INF""#, ": -Infinity");
        let mut out = text.into_bytes();
        out.push(b'\n');
        out.extend(compressed);
        std::fs::write(path, out).expect("written");
    }

    #[test]
    fn snf_blocks_are_compared_decompressed_and_the_command_is_ignored() {
        // F4: the read pass writes SAMPLE.snf, which is not text.
        let dir = temp_dir("snf");
        let (a, b) = (dir.join("a.snf"), dir.join("b.snf"));
        let blocks = [("chr1", 0, "candidates 1"), ("chr1", 5_000_000, "candidates 2"), ("chr2", 0, "3")];
        snf(&a, "lrpipe sv --snf a.snf", &blocks);
        snf(&b, "lrpipe sv --snf b.snf", &blocks);
        let Outcome::Compared { records, headers } = compare_files(&a, &b).expect("compared") else {
            panic!("both exist")
        };
        assert_eq!(records, Records::Identical(3));
        assert!(headers.is_empty(), "{headers:?}");

        snf(&b, "lrpipe sv --snf b.snf", &[blocks[0], ("chr1", 5_000_000, "other candidates"), blocks[2]]);
        let Outcome::Compared { records: Records::Different { only_old, only_new, first, .. }, .. } =
            compare_files(&a, &b).expect("compared")
        else {
            panic!("they differ")
        };
        assert_eq!((only_old, only_new, first.record), (1, 1, 2));
        assert_eq!(first.old.as_deref(), Some("chr1 block 5000000: 12 bytes decompressed"));
        assert_eq!(first.new.as_deref(), Some("chr1 block 5000000: 16 bytes decompressed"));
        let _ = std::fs::remove_dir_all(dir);
    }

    /// The ONT fixture's reads with `@PG CL:command` and every record but the `skip`th (from 0).
    fn fixture_bam(path: &Path, command: &str, skip: Option<usize>) {
        use pasta::bam;
        use pasta::sam::alignment::io::Write as _;
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture/ont.bam");
        let mut reader = bam::io::Reader::new(File::open(fixture).expect("the fixture"));
        let header = reader.read_header().expect("a header");
        let mut text = Vec::new();
        pasta::sam::io::Writer::new(&mut text).write_header(&header).expect("written");
        let text = format!("{}@PG\tID:test\tPN:test\tCL:{command}\n", String::from_utf8_lossy(&text));
        let header: pasta::sam::Header = text.parse().expect("a header");
        let mut writer = bam::io::Writer::new(File::create(path).expect("a BAM"));
        writer.write_header(&header).expect("written");
        for (i, record) in reader.record_bufs(&header).enumerate() {
            if Some(i) != skip {
                writer.write_alignment_record(&header, &record.expect("a record")).expect("written");
            }
        }
        writer.try_finish().expect("finished");
    }

    #[test]
    fn records_at_one_place_compare_as_a_multiset_and_a_bam_cut_short_is_refused() {
        use pasta::bam;
        use pasta::sam::alignment::io::Write as _;
        let dir = temp_dir("bam-places");
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture/ont.bam");
        let mut reader = bam::io::Reader::new(File::open(fixture).expect("the fixture"));
        let header = reader.read_header().expect("a header");
        let records: Vec<_> = reader.record_bufs(&header).take(4).map(|r| r.expect("a record")).collect();
        // Each record and a copy of it at the same place, the copy first or second.
        let write = |path: &Path, copy_first: bool, cut: usize| {
            let mut writer = bam::io::Writer::new(File::create(path).expect("a BAM"));
            writer.write_header(&header).expect("written");
            for record in &records {
                let mut copy = record.clone();
                *copy.name_mut() = Some(format!("copy-{}", record.name().expect("a name")).into());
                let pair = if copy_first { [&copy, record] } else { [record, &copy] };
                for r in pair {
                    writer.write_alignment_record(&header, r).expect("written");
                }
            }
            writer.try_finish().expect("finished");
            drop(writer);
            let bytes = std::fs::read(path).expect("written");
            std::fs::write(path, &bytes[..bytes.len() - cut]).expect("cut");
        };
        let (a, b, cut) = (dir.join("a.bam"), dir.join("b.bam"), dir.join("cut.bam"));
        write(&a, false, 0);
        write(&b, true, 0);
        let Outcome::Compared { records: found, .. } = compare_files(&a, &b).expect("compared") else {
            panic!("both are there")
        };
        assert_eq!(found, Records::Identical(8), "the same records at each place, in another order");
        // Cut at a record's end: its end-of-file block gone, it still reads as a BAM.
        write(&cut, false, BGZF_EOF.len());
        let e = format!("{:#}", anyhow::Error::from(compare_files(&a, &cut).expect_err("cut short")));
        assert!(e.contains("cut.bam") && e.contains("no BGZF end-of-file block"), "{e}");
        // A file named as text that is not: an error, not a comparison of its bytes.
        std::fs::write(dir.join("y.vcf"), [0xff, 0xfe, 0]).expect("written");
        let e = format!(
            "{:#}",
            anyhow::Error::from(compare_files(&dir.join("y.vcf"), &dir.join("y.vcf")).expect_err("not text"))
        );
        assert!(e.contains("not text, though its name says it is"), "{e}");
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn every_output_is_compared_bams_by_their_records() {
        // The aligned, merged and haplotagged reads are BAMs, which made `compare` fail ("stream
        // did not contain valid UTF-8"); another file that is not text is compared byte for byte.
        let dir = temp_dir("bams");
        let kinds = [
            ("aligned_bam", "align/S.bam"),
            ("merged_reads", "reads/S.bam"),
            ("haplotagged_reads", "phase/S.haplotagged.bam"),
        ];
        for (run, skip) in [("a", None), ("b", Some(3))] {
            let root = dir.join(run);
            let mut outputs = serde_json::Map::new();
            for (kind, file) in kinds {
                let path = root.join(file);
                std::fs::create_dir_all(path.parent().expect("a directory")).expect("directory");
                // Only the haplotagged reads differ; every @PG command line does.
                fixture_bam(&path, &format!("lrpipe {}", root.display()), skip.filter(|_| kind == "haplotagged_reads"));
                outputs.insert(kind.into(), path.into_os_string().into_string().expect("UTF-8").into());
            }
            std::fs::write(root.join("blob.bin"), [0xff, 0xfe, 0, if skip.is_some() { 1 } else { 0 }])
                .expect("written");
            outputs.insert("blob".into(), root.join("blob.bin").to_string_lossy().into_owned().into());
            let manifest = serde_json::json!({"output": root, "stages": [], "tools": {}, "outputs": outputs});
            std::fs::write(root.join("manifest.json"), manifest.to_string()).expect("written");
        }
        let comparison = compare(&Run::open(dir.join("a")).expect("a run"), &Run::open(dir.join("b")).expect("a run"))
            .expect("compared");
        let outcome = |kind: &str| comparison.files.iter().find(|f| f.kind == kind).map(|f| f.outcome.clone());
        let records = 152;
        for kind in ["aligned_bam", "merged_reads"] {
            assert_eq!(
                outcome(kind),
                Some(Outcome::Compared { records: Records::Identical(records), headers: HeaderChange::default() }),
                "{kind}"
            );
        }
        let Some(Outcome::Compared { records: Records::Different { old, new, only_old, only_new, first }, headers }) =
            outcome("haplotagged_reads")
        else {
            panic!("the haplotagged reads differ: {:?}", outcome("haplotagged_reads"))
        };
        assert_eq!((old, new, only_old, only_new, first.record), (records, records - 1, 1, 0, 4));
        assert!(first.old.as_deref().is_some_and(|r| r.contains(" at chr20:") && r.contains("MAPQ")), "{first:?}");
        assert!(headers.is_empty(), "the @PG command lines are left out: {headers:?}");
        let Some(Outcome::Compared { records: Records::Different { first, .. }, .. }) = outcome("blob") else {
            panic!("the blobs differ")
        };
        assert_eq!(first.old.as_deref(), Some("not text: 4 bytes, compared byte for byte"));
        assert!(!comparison.outputs_match());
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn the_merged_vcfs_prefixed_volatile_headers_are_ignored() {
        // LUX-12: two identical runs' merged VCFs differed in these.
        for line in ["##fileDate=1", "##snv_cmdline=x", "##sv_command=y", "##sv_fileDate=2", "##tr_command=z"] {
            assert!(volatile(line), "{line}");
        }
        for line in ["##fileformat=VCFv4.2", "##source=x", "##commandline_options=a", "##subcommand=b", "#CHROM"] {
            assert!(!volatile(line), "{line}");
        }
    }

    #[test]
    fn skipped_stages_and_idle_runs_get_no_change_in_per_cent() {
        let stage = |status: &str, wall: f64| StageMeasurements {
            name: "snv".into(),
            wall_seconds: Some(wall),
            max_rss_gib: Some(2.0),
            status: Some(status.into()),
        };
        let comparison = |old: &str, new: &str| Comparison {
            files: Vec::new(),
            tools: Vec::new(),
            stages: vec![StageComparison {
                name: "snv".into(),
                old: Some(stage(old, 100.0)),
                new: Some(stage(new, 50.0)),
            }],
            wall_seconds: (Some(102.0), Some(0.3)),
        };
        let text = comparison("succeeded", "skipped").to_string();
        assert!(!text.contains('%'), "{text}");
        assert!(text.contains("(skipped in the new run)"), "{text}");
        let text = comparison("succeeded", "succeeded").to_string();
        assert!(text.contains("-50%"), "{text}");
        // Under a second shows as 0: no per cent of it.
        assert!(!change(Some(0.3), Some(7.0), 0, true).contains('%'));
        assert!(change(Some(1.0), Some(7.0), 0, true).contains("+600%"));
    }

    #[test]
    fn multiset_difference_counts_repeats() {
        assert_eq!(multiset_difference(&[1, 1, 2], &[1, 2, 2, 3]), (1, 2));
        assert_eq!(multiset_difference(&[], &[5]), (0, 1));
    }
}

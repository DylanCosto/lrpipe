//! `--bedmethyl`: CpG methylation counted from the read pass's reads (lrpipe-methyl), written
//! in the contigs' order as lrpipe-methyl writes it from its own reading.

use std::collections::{HashMap, VecDeque};
use std::ffi::OsString;
use std::io::Write;
use std::thread::JoinHandle;

use anyhow::{Context, Result, bail};
use lrpipe_methyl::input::Fetched;
use lrpipe_methyl::output::Output;
use lrpipe_methyl::stream::{Contig, Counted, Pass};
use rayon::prelude::*;

use crate::fasta::Fasta;

/// Windows whose rows are made in parallel before they are written.
const WRITE_BATCH: usize = 32;

/// The rows of the contigs the pass does not count, by contig.
type Rest = HashMap<String, String>;

/// Counting and writing the table.
#[derive(Debug)]
pub struct Methylation {
    pub pass: Pass,
    writer: Writer,
}

/// The table, written a contig at a time.
pub struct Writer {
    path: String,
    out: Output,
    /// The contigs the pass does not count, in the header's order, with their number there.
    rest: VecDeque<(usize, String)>,
    /// Their rows, counted on other threads while the pass runs.
    counting: Option<JoinHandle<Result<(Rest, Fetched)>>>,
    rest_rows: Rest,
    fetched: Fetched,
}

impl std::fmt::Debug for Writer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Writer").field("path", &self.path).finish_non_exhaustive()
    }
}

impl Methylation {
    /// Estimates the pass thresholds (reading a sample of the reads of `regions`, every contig
    /// when empty), and starts counting the `wanted` contigs (header numbers, names and the
    /// regions of each to count, in order) that are not `counted` by the pass, on threads of
    /// their own.
    ///
    /// # Errors
    ///
    /// If the reads cannot be sampled or the table cannot be created.
    pub fn start(
        input: &str,
        reference: &str,
        path: &str,
        threads: usize,
        regions: &[String],
        wanted: &[(usize, String, Vec<String>)],
        counted: &[usize],
    ) -> Result<Self> {
        let args: Vec<OsString> = [
            input,
            path,
            "--ref",
            reference,
            "--threads",
            &threads.max(1).to_string(),
            "--modified-bases",
            "5mC",
            "5hmC",
            "--cpg",
            "--combine-strands",
        ]
        .map(OsString::from)
        .into_iter()
        .chain(regions.iter().flat_map(|r| [OsString::from("--region"), OsString::from(r)]))
        .collect();
        let job = lrpipe_methyl::cli::pileup_job_from(&args).context("--bedmethyl")?;
        let pass = Pass::new(&job)?;
        let rest: Vec<&(usize, String, Vec<String>)> = wanted.iter().filter(|(t, ..)| !counted.contains(t)).collect();
        let counting = (!rest.is_empty()).then(|| {
            let mut job = job.clone();
            job.regions = rest.iter().flat_map(|(.., regions)| regions.iter().cloned()).collect();
            std::thread::spawn(move || count_rest(&job))
        });
        let rest: VecDeque<(usize, String)> = rest.into_iter().map(|(t, name, _)| (*t, name.clone())).collect();
        let out = Output::create(path.as_ref()).with_context(|| format!("creating {path}"))?;
        let writer =
            Writer { path: path.to_owned(), out, rest, counting, rest_rows: Rest::new(), fetched: Fetched::default() };
        Ok(Methylation { pass, writer })
    }

    /// The pass's counts and the table's writer, apart.
    pub fn parts(&mut self) -> (&Pass, &mut Writer) {
        (&self.pass, &mut self.writer)
    }

    /// Writes the contigs left, and the table's end.
    ///
    /// # Errors
    ///
    /// If the table cannot be written.
    pub fn finish(self) -> Result<()> {
        let mut w = self.writer;
        w.write_rest_before(usize::MAX)?;
        let path = w.path;
        w.out.finish().with_context(|| format!("writing {path}"))?;
        lrpipe_methyl::cli::report_fetched(&w.fetched);
        Ok(())
    }
}

/// The rows of `job`'s regions (whole contigs, or the regions of `--bedmethyl-regions`), by contig, and what their fetches took.
fn count_rest(job: &lrpipe_methyl::run::Job) -> Result<(Rest, Fetched)> {
    let path = job.outputs.first().map(|(_, p)| p.clone()).unwrap_or_default();
    let (fetched, mut text) = lrpipe_methyl::run::run_into(job, || Ok(vec![(path, Vec::<u8>::new())]))?;
    let text = String::from_utf8(text.pop().map(|(_, t)| t).unwrap_or_default()).context("bedMethyl rows")?;
    let mut rows = Rest::new();
    for line in text.lines() {
        let contig = line.split('\t').next().unwrap_or_default();
        let t: &mut String = rows.entry(contig.to_owned()).or_default();
        t.push_str(line);
        t.push('\n');
    }
    Ok((rows, fetched))
}

/// The contig's bases `[start, end)`, clipped to it.
fn fetch(fasta: &Fasta, name: &str, start: i64, end: i64) -> std::io::Result<Vec<u8>> {
    fasta.contig(name).map_or(Ok(Vec::new()), |c| c.fetch(start, end))
}

/// The contig `name` (`len` bases) for the pass to count, when the reference has it.
///
/// # Errors
///
/// If the reference cannot be read.
pub fn contig(pass: &Pass, fasta: &Fasta, name: &str, len: i64) -> Result<Option<Contig>> {
    if fasta.contig(name).is_none() {
        return Ok(None);
    }
    let contig = pass.contig(name, len, |s, e| fetch(fasta, name, s, e)).with_context(|| format!("reading {name}"))?;
    Ok(Some(contig))
}

impl Writer {
    /// Writes contig `tid`'s rows from the pass's counts (none when the reference lacks it),
    /// after those of the contigs before it that the pass does not count.
    ///
    /// # Errors
    ///
    /// If the table cannot be written, or the reference read.
    pub fn write(&mut self, tid: usize, name: &str, counted: Option<Counted>, fasta: &Fasta) -> Result<()> {
        self.write_rest_before(tid)?;
        let Some(counted) = counted else { return Ok(()) };
        self.fetched += counted.fetched();
        let n = counted.windows();
        for batch in (0..n).step_by(WRITE_BATCH) {
            let texts: Vec<String> = (batch..(batch + WRITE_BATCH).min(n))
                .into_par_iter()
                .map(|k| {
                    let mut text = String::new();
                    counted.write_window(k, |s, e| fetch(fasta, name, s, e), &mut text)?;
                    Ok(text)
                })
                .collect::<std::io::Result<_>>()
                .with_context(|| format!("reading {name}"))?;
            for text in texts {
                self.out.write_all(text.as_bytes()).with_context(|| format!("writing {}", self.path))?;
            }
        }
        Ok(())
    }

    /// Writes the rows of the contigs the pass does not count that come before contig `tid`.
    fn write_rest_before(&mut self, tid: usize) -> Result<()> {
        while self.rest.front().is_some_and(|&(t, _)| t < tid) {
            if let Some(counting) = self.counting.take() {
                let Ok(result) = counting.join() else { bail!("counting methylation of the short contigs panicked") };
                let (rows, fetched) = result?;
                self.rest_rows = rows;
                self.fetched += fetched;
            }
            let (_, name) = self.rest.pop_front().expect("a contig");
            if let Some(text) = self.rest_rows.remove(&name) {
                self.out.write_all(text.as_bytes()).with_context(|| format!("writing {}", self.path))?;
            }
        }
        Ok(())
    }
}

//! Pileup over regions in windows, in parallel, written in order.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, mpsc};

use crate::cpg::{CpgCounts, Sites, Values};
use crate::input::{Fetched, Haplotags, InputError, ReadFilter, Reads, SharedReference};
use crate::mods::Mode;
use crate::output::Output;
use crate::pileup::{Calling, Counts, Partition, RefWindow, Rows, Strands, Thresholds, most_likely};
use crate::read::{Convention, ModRead, Strand};

/// Why a pileup could not be made.
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    /// Reads could not be read.
    #[error(transparent)]
    Input(#[from] InputError),
    /// An output could not be written.
    #[error("writing {}", path.display())]
    Write {
        /// The output.
        path: PathBuf,
        /// The cause.
        #[source]
        source: std::io::Error,
    },
    /// A region names a contig the reads do not have, or is malformed.
    #[error("bad region {0:?}")]
    Region(String),
    /// The reference lacks a contig of the reads, or has it shorter: it is not the one the
    /// reads were aligned to.
    #[error("--ref {}: {problem}; is it the reference the reads were aligned to?", reference.display())]
    Reference {
        /// The reference.
        reference: PathBuf,
        /// What does not match.
        problem: String,
    },
}

/// A stretch of one contig.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Window {
    /// Contig name.
    pub contig: String,
    /// 0-based start.
    pub start: i64,
    /// Exclusive end.
    pub end: i64,
}

/// How a region's `start-end` is read.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Coordinates {
    /// 1-based, both ends included (samtools).
    #[default]
    OneBased,
    /// A 0-based start, and the end excluded (modkit): one base later than
    /// [`Coordinates::OneBased`] at the start, the same at the end.
    ZeroBased,
}

/// Splits regions (`contig`, `contig:start-end` as `coordinates` has them, or all contigs of
/// `header` when none) into windows of at most `size` bases.
///
/// # Errors
///
/// If a region is malformed or names a contig the header does not have.
pub fn windows(
    header: &pasta::sam::Header,
    regions: &[String],
    coordinates: Coordinates,
    size: i64,
) -> Result<Vec<Window>, RunError> {
    let first = match coordinates {
        Coordinates::OneBased => 1,
        Coordinates::ZeroBased => 0,
    };
    let contigs: Vec<(String, i64)> = header
        .reference_sequences()
        .iter()
        .map(|(name, rs)| (name.to_string(), usize::from(rs.length()) as i64))
        .collect();
    let spans: Vec<(String, i64, i64)> = if regions.is_empty() {
        contigs.iter().map(|(n, l)| (n.clone(), 0, *l)).collect()
    } else {
        regions
            .iter()
            .map(|r| {
                let bad = || RunError::Region(r.clone());
                let (name, range) = match r.rsplit_once(':') {
                    Some((n, range)) if contigs.iter().any(|(c, _)| c == n) => (n, Some(range)),
                    _ => (r.as_str(), None),
                };
                let len = contigs.iter().find(|(c, _)| c == name).ok_or_else(bad)?.1;
                let (s, e) = match range {
                    None => (0, len),
                    Some(range) => {
                        let (s, e) = range.split_once('-').ok_or_else(bad)?;
                        let num = |x: &str| x.replace(',', "").parse::<i64>().map_err(|_| bad());
                        (num(s)? - first, num(e)?.min(len))
                    }
                };
                Ok((name.to_owned(), s.max(0), e))
            })
            .collect::<Result<_, RunError>>()?
    };
    let size = size.max(1);
    Ok(spans
        .into_iter()
        .flat_map(|(contig, s, e)| {
            (s..e).step_by(size as usize).map(move |w| Window {
                contig: contig.clone(),
                start: w,
                end: (w + size).min(e),
            })
        })
        .collect())
}

/// What a pileup run reads and writes.
#[derive(Clone, Debug)]
pub struct Job {
    /// Indexed BAM or CRAM.
    pub reads: PathBuf,
    /// Reference FASTA (needed for CRAM).
    pub reference: Option<PathBuf>,
    /// Regions; every contig when empty.
    pub regions: Vec<String>,
    /// How `regions` give their coordinates.
    pub coordinates: Coordinates,
    /// Window size in bases.
    pub window: i64,
    /// Worker threads.
    pub threads: usize,
    /// Which alignments count.
    pub filter: ReadFilter,
    /// What the tables report.
    pub report: Report,
    /// Tables to write: each partition of the reads, with its output.
    pub outputs: Vec<(Partition, PathBuf)>,
    /// The reads' haplotypes from a haplotag list, in place of the alignments' `HP` tags.
    pub haplotags: Option<Arc<Haplotags>>,
}

/// What a table reports.
#[derive(Clone, Debug)]
pub enum Report {
    /// modkit's bedMethyl: counts per position, strand and modification.
    BedMethyl {
        /// Thresholds and conventions for classifying calls.
        calling: Calling,
        /// Which rows the tables have.
        rows: Rows,
    },
    /// pb-CpG-tools' layouts: CpG sites from 5mC calls.
    CpgCount {
        /// Which sites.
        sites: Sites,
        /// How a site's level is estimated.
        estimate: Estimate,
        /// Sites with fewer calls are left out.
        min_coverage: u32,
        /// Recorded in the header.
        min_mapq: u8,
    },
}

/// How a CpG site's level is estimated.
#[derive(Clone, Debug)]
pub enum Estimate {
    /// The share of reads with a value over one half (pb-CpG-tools' count mode).
    Count,
    /// The expected level under a fixed model of methylated and unmethylated reads' values
    /// ([`crate::model`]).
    Model,
}

impl Report {
    /// Bases the window's counts need past its end: a CpG's G.
    pub(crate) fn overhang(&self) -> i64 {
        match self {
            Report::BedMethyl { rows, .. } if rows.strands == Strands::Combined => 1,
            Report::BedMethyl { .. } => 0,
            Report::CpgCount { .. } => 1,
        }
    }

    fn needs_reference(&self) -> bool {
        match self {
            Report::BedMethyl { rows, .. } => rows.modified_bases.is_some() || rows.motif != crate::pileup::Motif::Any,
            Report::CpgCount { .. } => true,
        }
    }

    /// The header lines of a table, if it has any.
    fn header(&self, partition: Partition) -> String {
        let _ = partition;
        match self {
            Report::BedMethyl { .. } => String::new(),
            Report::CpgCount { sites, estimate, min_coverage, min_mapq } => format!(
                "##lrpipe-methyl-version={}\n##pileup-mode={}\n##modsites-mode={}\n##min-coverage={min_coverage}\n##min-mapq={min_mapq}\n#chrom\tbegin\tend\tmod_score\ttype\tcov\t{}\n",
                env!("CARGO_PKG_VERSION"),
                match estimate {
                    Estimate::Count => "count",
                    Estimate::Model => "model",
                },
                match sites {
                    Sites::Reference => "reference",
                    Sites::Denovo => "denovo",
                },
                match estimate {
                    Estimate::Count => "mod_count\tunmod_count\tavg_mod_score\tavg_unmod_score",
                    Estimate::Model => "est_mod_count\test_unmod_count\tdiscretized_mod_score",
                },
            ),
        }
    }

    /// One window's text for each output, from its reads, which are counted as `reads`
    /// fetches them (only denovo CpG sites, which depend on every read, gather them in `buf`
    /// first). Returns what the fetch took and left out, and the texts.
    fn window(
        &self,
        job: &Job,
        reads: &mut Reads,
        w: &Window,
        reference: Option<&RefWindow>,
        buf: &mut Vec<ModRead>,
    ) -> Result<(Fetched, Vec<String>), InputError> {
        let end = w.end + self.overhang();
        match self {
            Report::BedMethyl { calling, rows } => {
                // With CpG rows only CpGs' C and G are counted.
                let mask: Option<Vec<bool>> = (rows.motif == crate::pileup::Motif::CpG).then(|| {
                    let base = |p: i64| reference.and_then(|r| r.base(p));
                    (w.start..end)
                        .map(|p| {
                            (base(p) == Some(b'C') && base(p + 1) == Some(b'G'))
                                || (base(p - 1) == Some(b'C') && base(p) == Some(b'G'))
                        })
                        .collect()
                });
                let mut counts: Vec<Counts> = job.outputs.iter().map(|_| Counts::empty(w.start, end)).collect();
                let took = reads.fetch_each(&w.contig, w.start, end, job.filter, |read| {
                    for (c, &(partition, _)) in counts.iter_mut().zip(&job.outputs) {
                        c.add(&read, *calling, partition, mask.as_deref());
                    }
                })?;
                let texts = counts
                    .into_iter()
                    .map(|c| {
                        let mut text = String::new();
                        c.finish().write_bedmethyl(&w.contig, w.start, w.end, rows, reference, &mut text);
                        text
                    })
                    .collect();
                Ok((took, texts))
            }
            Report::CpgCount { min_coverage, sites, estimate, .. } => {
                let empty = RefWindow { start: w.start, seq: Vec::new() };
                let reference = reference.unwrap_or(&empty);
                let partitions: Vec<Partition> = job.outputs.iter().map(|&(p, _)| p).collect();
                let keep = match estimate {
                    Estimate::Count => Values::Dropped,
                    Estimate::Model => Values::Kept,
                };
                let (took, counts) = match sites {
                    Sites::Reference => {
                        let mut counts = CpgCounts::at_reference_sites(w.start, w.end, reference, &partitions, keep);
                        let took = reads.fetch_each(&w.contig, w.start, end, job.filter, |read| counts.add(&read))?;
                        (took, counts)
                    }
                    Sites::Denovo => {
                        buf.clear();
                        let took = reads.fetch(&w.contig, w.start, end, job.filter, buf)?;
                        let counts = CpgCounts::new(buf, w.start, w.end, reference, &partitions, *sites, keep);
                        buf.clear();
                        (took, counts)
                    }
                };
                let texts = partitions
                    .iter()
                    .enumerate()
                    .map(|(k, p)| {
                        let label = match p {
                            Partition::All => "Total".to_owned(),
                            Partition::Haplotype(h) => format!("hap{h}"),
                        };
                        let mut text = String::new();
                        match estimate {
                            Estimate::Count => counts.write(&w.contig, k, &label, *min_coverage, &mut text),
                            Estimate::Model => counts.write_model(&w.contig, k, &label, *min_coverage, &mut text),
                        }
                        text
                    })
                    .collect();
                Ok((took, texts))
            }
        }
    }
}

fn open_reads(job: &Job, reference: Option<&SharedReference>) -> Result<Reads, RunError> {
    let mut reads = Reads::open(&job.reads, reference)?;
    if let Some(haplotags) = &job.haplotags {
        reads.set_haplotags(Arc::clone(haplotags));
    }
    Ok(reads)
}

/// Drops the contigs the reference holds, when there is one.
fn clear(reference: Option<&SharedReference>) {
    if let Some(r) = reference {
        r.clear();
    }
}

/// Fails unless the reference has every contig of `windows` that holds reads, as long as the
/// reads' header has it: a contig it lacks, or the part past its end, would get no rows, and a
/// contig of another length is another assembly's (CHM13's chr2 against GRCh38's), whose bases
/// would give the reads' positions the wrong CpGs. Contigs without reads, which get no rows
/// anyway, are not checked, so a reference of some contigs serves reads of those whose header
/// names a whole genome.
fn check_reference(
    job: &Job,
    reference: Option<&SharedReference>,
    reads: &Reads,
    windows: &[Window],
) -> Result<(), RunError> {
    let Some(reference) = reference else { return Ok(()) };
    let mismatch =
        |problem: String| RunError::Reference { reference: job.reference.clone().unwrap_or_default(), problem };
    let mut checked: Option<&str> = None;
    for w in windows {
        if checked == Some(w.contig.as_str()) {
            continue;
        }
        checked = Some(&w.contig);
        let Some(in_reads) =
            reads.header().reference_sequences().get(w.contig.as_bytes()).map(|rs| usize::from(rs.length()) as u64)
        else {
            continue;
        };
        let in_reference = reference.length(&w.contig);
        if in_reference == Some(in_reads) || reads.indexed_bytes(&w.contig, 0, in_reads as i64)? == 0 {
            continue;
        }
        return Err(mismatch(match in_reference {
            None => format!("no contig {} (the reads have it)", w.contig),
            Some(len) => format!("contig {} has {len} bases, but {in_reads} in the reads' header", w.contig),
        }));
    }
    Ok(())
}

fn shared_reference(job: &Job) -> Result<Option<SharedReference>, RunError> {
    match &job.reference {
        Some(r) => Ok(Some(SharedReference::open(r)?)),
        None => Ok(None),
    }
}

/// Estimates each canonical base's pass threshold as the `percentile` (0 to 1) of the
/// confidence of calls in reads sampled evenly across the job's regions: up to `per_window`
/// reads from each of up to `windows` windows.
///
/// # Errors
///
/// If the reads cannot be read.
pub fn estimate_thresholds(
    job: &Job,
    percentile: f64,
    windows_wanted: usize,
    per_window: usize,
) -> Result<(Thresholds, [u64; 5]), RunError> {
    let reference = shared_reference(job)?;
    let reads = open_reads(job, reference.as_ref())?;
    // Sampled over windows of a size of their own, not the job's, which is a speed setting that
    // depends on the format (larger for CRAM) and on --window: with the job's, a CRAM sampled
    // fewer windows than the same reads as BAM and estimated other thresholds.
    let all = windows(reads.header(), &job.regions, job.coordinates, SAMPLE_WINDOW)?;
    // Before the estimate, which can take a while, rather than after it.
    if job.report.needs_reference() {
        check_reference(job, reference.as_ref(), &reads, &all)?;
    }
    drop(reads);
    let step = (all.len() / windows_wanted.max(1)).max(1);
    let chosen: Vec<(usize, &Window)> = all.iter().enumerate().step_by(step).collect();
    let next = AtomicUsize::new(0);
    // Confidence histograms per base slot (A, C, G, T, other), on the ML scale, summed over
    // windows sampled on every thread.
    let hist = std::thread::scope(|s| -> Result<Vec<[u64; HIST_BINS]>, RunError> {
        let workers: Vec<_> = (0..job.threads.max(1))
            .map(|_| {
                s.spawn(|| -> Result<Vec<[u64; HIST_BINS]>, RunError> {
                    let mut reads = open_reads(job, reference.as_ref())?;
                    let mut hist = vec![[0u64; HIST_BINS]; 5];
                    let mut buf = Vec::new();
                    loop {
                        let k = next.fetch_add(1, Ordering::Relaxed);
                        let Some(&(i, w)) = chosen.get(k) else { break };
                        if k > 0 && chosen[k - 1].1.contig != w.contig {
                            clear(reference.as_ref());
                        }
                        // Windows sampled back to back share reads: each is sampled by the
                        // window it starts in (the first window takes those starting before it
                        // too).
                        let contiguous = step == 1 && i > 0 && all[i - 1].contig == w.contig;
                        let wanted = |r: &&ModRead| !contiguous || r.start >= w.start;
                        // The window's first reads in file order are among those of its start,
                        // unless that holds too few.
                        buf.clear();
                        reads.fetch(&w.contig, w.start, (w.start + SAMPLE_SPAN).min(w.end), job.filter, &mut buf)?;
                        if buf.iter().filter(wanted).count() < per_window && w.start + SAMPLE_SPAN < w.end {
                            buf.clear();
                            reads.fetch(&w.contig, w.start, w.end, job.filter, &mut buf)?;
                        }
                        for r in buf.iter().filter(wanted).take(per_window) {
                            sample(r, &mut hist);
                        }
                    }
                    Ok(hist)
                })
            })
            .collect();
        let mut total = vec![[0u64; HIST_BINS]; 5];
        for w in workers {
            for (t, h) in total.iter_mut().zip(w.join().expect("a sampling thread panicked")?) {
                for (a, b) in t.iter_mut().zip(h) {
                    *a += b;
                }
            }
        }
        Ok(total)
    })?;
    let mut t = Thresholds::NONE;
    let mut totals = [0u64; 5];
    for (slot, base) in b"ACGTN".iter().copied().enumerate() {
        let total: u64 = hist[slot].iter().sum();
        totals[slot] = total;
        if total == 0 {
            continue;
        }
        let want = (percentile * total as f64).floor() as u64;
        let mut seen = 0;
        for (q, &n) in hist[slot].iter().enumerate() {
            seen += n;
            if seen > want {
                t.set(base, crate::pileup::probability(q as u16));
                break;
            }
        }
    }
    Ok((t, totals))
}

/// Bases per window the threshold estimate samples reads from, whatever the format or
/// `--window` (the BAM default, which the estimates validated against modkit used).
const SAMPLE_WINDOW: i64 = 250_000;

/// The span at a sampled window's start whose reads are tried first.
const SAMPLE_SPAN: i64 = 10_000;

/// Histogram bins: the ML scale, 0 to 255.
const HIST_BINS: usize = 256;

/// Adds every call of `read`, aligned or not, to the confidence histograms: listed bases with
/// their most likely state's `ML`-scale value (the canonical base's being 255 less the
/// modifications', as modkit estimates thresholds), unlisted ones of implicit tracks as
/// certain.
fn sample(read: &ModRead, hist: &mut [[u64; HIST_BINS]]) {
    for t in &read.tracks {
        let slot = match t.base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        let mut listed = 0u64;
        for (_, values) in t.listed() {
            if let Some((_, q)) = most_likely(values, Convention::Reference) {
                // Back from the bin middle to the bin.
                hist[slot][((q * 256.0 - 0.5).round() as usize).min(HIST_BINS - 1)] += 1;
            }
            listed += 1;
        }
        if t.unlisted == Mode::Implicit {
            let in_seq = match t.strand {
                Strand::Forward => t.base,
                Strand::Reverse => complement(t.base),
            };
            let present = read.seq.iter().filter(|&&b| t.base == b'N' || b == in_seq).count() as u64;
            hist[slot][HIST_BINS - 1] += present.saturating_sub(listed);
        }
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// How much data, in compressed bytes by the index, the windows in progress may hold at
/// once: [`IN_FLIGHT`] times a median window's on every thread. A thread's memory peaks with
/// the reads and decoded CRAM slices of its window, which in deep regions (centromeres and
/// segmental duplications, hundreds to thousands of reads deep) are many times a normal
/// window's whatever its width; slow to process, those windows would otherwise gather every
/// thread at once.
fn data_budget(costs: &[u64], threads: usize) -> u64 {
    let mut sorted: Vec<u64> = costs.iter().copied().filter(|&c| c > 0).collect();
    sorted.sort_unstable();
    let median = sorted.get(sorted.len() / 2).copied().unwrap_or(0);
    (median as f64 * IN_FLIGHT * threads as f64) as u64
}

/// Median windows' worth of data in progress per thread, at most.
const IN_FLIGHT: f64 = 1.5;

/// Windows a worker may run ahead of the writer, per thread: each finished window's text
/// waits in memory until the windows before it are written.
const AHEAD_PER_THREAD: usize = 2;

/// Runs the pileup and writes every output. Returns what the fetches took and left out.
///
/// # Errors
///
/// If reads cannot be read or an output written.
///
/// # Panics
///
/// If a worker thread panics.
pub fn run(job: &Job) -> Result<Fetched, RunError> {
    let open = || {
        job.outputs
            .iter()
            .map(|&(partition, ref path)| {
                let write = |source| RunError::Write { path: path.clone(), source };
                let mut w = Output::create(path).map_err(write)?;
                w.write_all(job.report.header(partition).as_bytes()).map_err(write)?;
                Ok((path.clone(), w))
            })
            .collect()
    };
    let (fetched, writers) = run_into(job, open)?;
    for (path, w) in writers {
        w.finish().map_err(|source| RunError::Write { path, source })?;
    }
    Ok(fetched)
}

/// Runs the pileup, writing each output (in the order of the job's) to the writers `open`
/// gives once the regions are known to be valid, and returns what the fetches took and left
/// out and the writers. Nothing is written for a table's header.
///
/// # Errors
///
/// If reads cannot be read, `open` fails or an output cannot be written.
///
/// # Panics
///
/// If a worker thread panics.
pub fn run_into<W: Write>(
    job: &Job,
    open: impl FnOnce() -> Result<Vec<(PathBuf, W)>, RunError>,
) -> Result<(Fetched, Vec<(PathBuf, W)>), RunError> {
    let reference = shared_reference(job)?;
    let header_reads = open_reads(job, reference.as_ref())?;
    let all = windows(header_reads.header(), &job.regions, job.coordinates, job.window)?;
    if job.report.needs_reference() {
        check_reference(job, reference.as_ref(), &header_reads, &all)?;
    }
    let costs: Vec<u64> =
        all.iter().map(|w| header_reads.indexed_bytes(&w.contig, w.start, w.end)).collect::<Result<_, _>>()?;
    drop(header_reads);
    let mut writers = open()?;

    let next = AtomicUsize::new(0);
    let written = Mutex::new(0usize);
    let progress = Condvar::new();
    let fetched = Mutex::new(Fetched::default());
    let stop = AtomicBool::new(false);
    let threads = job.threads.max(1);
    let ahead = threads * AHEAD_PER_THREAD;
    let budget = data_budget(&costs, threads);
    let in_flight = Mutex::new(0u64);
    let freed = Condvar::new();
    let (tx, rx) = mpsc::channel::<Result<(usize, Vec<String>), RunError>>();
    std::thread::scope(|s| -> Result<(), RunError> {
        for _ in 0..threads {
            let tx = tx.clone();
            let (next, written, progress, fetched, stop, all, reference) =
                (&next, &written, &progress, &fetched, &stop, &all, &reference);
            let (costs, in_flight, freed) = (&costs, &in_flight, &freed);
            s.spawn(move || {
                let mut reads = match open_reads(job, reference.as_ref()) {
                    Ok(r) => r,
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                };
                let mut buf = Vec::new();
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(w) = all.get(i) else { break };
                    {
                        let mut done = written.lock().expect("lock poisoned");
                        while i >= *done + ahead && !stop.load(Ordering::Relaxed) {
                            done = progress.wait(done).expect("lock poisoned");
                        }
                    }
                    if stop.load(Ordering::Relaxed) {
                        break;
                    }
                    // Wait until the windows in progress leave room for this one's data (one
                    // window always goes ahead, however much it holds).
                    let cost = costs[i];
                    {
                        let mut held = in_flight.lock().expect("lock poisoned");
                        while *held > 0 && *held + cost > budget && !stop.load(Ordering::Relaxed) {
                            held = freed.wait(held).expect("lock poisoned");
                        }
                        *held += cost;
                    }
                    buf.clear();
                    let end = w.end + job.report.overhang();
                    let result = (|| -> Result<(usize, Vec<String>), RunError> {
                        let refwin = match (reference.as_ref(), job.report.needs_reference()) {
                            (Some(r), true) => r.window(&w.contig, w.start - 1, end + 1).map_err(|source| {
                                RunError::Input(InputError::Io {
                                    path: job.reference.clone().unwrap_or_default(),
                                    source,
                                })
                            })?,
                            _ => None,
                        };
                        let (took, texts) = job.report.window(job, &mut reads, w, refwin.as_ref(), &mut buf)?;
                        *fetched.lock().expect("lock poisoned") += took;
                        Ok((i, texts))
                    })();
                    buf.clear();
                    *in_flight.lock().expect("lock poisoned") -= cost;
                    freed.notify_all();
                    if tx.send(result).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx);
        let mut pending: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        let mut done = 0usize;
        let drain = || -> Result<(), RunError> {
            for result in rx {
                let (i, texts) = result?;
                pending.insert(i, texts);
                while let Some(texts) = pending.remove(&done) {
                    for ((path, w), text) in writers.iter_mut().zip(texts) {
                        w.write_all(text.as_bytes())
                            .map_err(|source| RunError::Write { path: path.clone(), source })?;
                    }
                    done += 1;
                    // Every window of a written contig is done: no worker reads it again.
                    if all.get(done).is_some_and(|next| next.contig != all[done - 1].contig) {
                        clear(reference.as_ref());
                    }
                    *written.lock().expect("lock poisoned") = done;
                    progress.notify_all();
                }
            }
            Ok(())
        };
        let result = drain();
        if result.is_err() {
            // Wake workers waiting for the writer so they see the stop and the scope can end.
            stop.store(true, Ordering::Relaxed);
            let _guard = written.lock().expect("lock poisoned");
            progress.notify_all();
        }
        result
    })?;
    Ok((fetched.into_inner().expect("lock poisoned"), writers))
}

/// The path with `suffix` added before the extension `.bed` (or `.bed.gz`).
#[must_use]
pub fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let (stem, ext) = match name.find(".bed") {
        Some(i) => (&name[..i], &name[i..]),
        None => (name.as_str(), ""),
    };
    path.with_file_name(format!("{stem}{suffix}{ext}"))
}

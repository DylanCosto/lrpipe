//! The whole-file driver: mosdepth's main loop, with contigs counted in parallel.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc;

use anyhow::{Result, anyhow};

use crate::counter::{ContigCounter, Filter, Mode};
use crate::reader::Alignments;
use crate::regions::By;
use crate::report::{ContigReport, ReportOptions, Scratch, Writer};

#[derive(Debug, Default)]
pub struct RunOptions {
    pub mode: Mode,
    pub filter: Filter,
    pub report: ReportOptions,
    /// `-c`: only this contig.
    pub chrom: Option<String>,
    /// Threads: contigs with reads are counted in parallel (each holds a 4-byte counter per
    /// base of its contig), and threads left over decompress BAM blocks in parallel.
    pub threads: usize,
}

/// Computes depth for the alignment file `path` and writes mosdepth's outputs under `prefix`.
pub fn run(path: &Path, reference: Option<&Path>, prefix: &str, opts: &RunOptions) -> Result<()> {
    let aln = Alignments::open(path, reference)?;
    let contigs = aln.contigs();
    let targets: Vec<usize> = match &opts.chrom {
        Some(c) => vec![contigs.iter().position(|(n, _)| n == c).ok_or_else(|| anyhow!("chromosome {c} not found"))?],
        None => (0..contigs.len()).collect(),
    };
    let max_len = contigs.iter().map(|c| c.1).max().unwrap_or(0);
    let mut writer = Writer::new(prefix, &opts.report, max_len)?;

    let work: Vec<usize> = targets.into_iter().filter(|&t| !opts.report.skips(&contigs[t].0)).collect();
    count_contigs(&aln, &work, opts, |_, report| writer.write_contig(&report))?;
    writer.finish()?;

    if let (Some(By::Bed(bed)), None) = (&opts.report.by, &opts.chrom) {
        let names: Vec<String> = contigs.iter().map(|c| c.0.clone()).collect();
        for (chrom, n) in bed.missing(&names) {
            eprintln!("[lrpipe-depth] warning chromosome:{chrom} from bed with {n} regions not found");
        }
    }
    Ok(())
}

/// Counts contigs `work` (indexes into `aln.contigs()`) as [`run`] does, in parallel, and calls
/// `each` with every contig's report, in the order of `work`.
pub fn count_contigs(
    aln: &Alignments,
    work: &[usize],
    opts: &RunOptions,
    mut each: impl FnMut(usize, ContigReport) -> Result<()>,
) -> Result<()> {
    let contigs = aln.contigs();
    // Contigs with reads first, largest first, so the last contig to finish is a small one.
    let with_reads: Vec<bool> = work.iter().map(|&t| aln.has_records(t)).collect();
    let mut order: Vec<usize> = (0..work.len()).collect();
    order.sort_by_key(|&k| (!with_reads[k], std::cmp::Reverse(contigs[work[k]].1)));

    // One worker per contig with reads, up to `threads`; the rest of the threads decompress.
    let total = opts.threads.max(1);
    let threads = total.min(with_reads.iter().filter(|&&r| r).count()).max(1);
    let decompression_threads = total / threads;
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let (tx, rx) = mpsc::channel::<(usize, Result<ContigReport>)>();

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..threads {
            let tx = tx.clone();
            let (order, next, failed) = (&order, &next, &failed);
            // pasta's CRAM decoding needs more than the default 2 MiB of stack.
            let spawned = std::thread::Builder::new().stack_size(64 << 20).spawn_scoped(s, move || {
                let mut reader = match aln.reader(decompression_threads) {
                    Ok(r) => r,
                    Err(e) => {
                        failed.store(true, Ordering::Relaxed);
                        let _ = tx.send((usize::MAX, Err(e)));
                        return;
                    }
                };
                let mut scratch = Scratch::new(&opts.report);
                let mut buf = Vec::new();
                let need_rg = !opts.filter.read_groups.is_empty();
                while !failed.load(Ordering::Relaxed) {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&k) = order.get(i) else { break };
                    let tid = work[k];
                    let (name, len) = (&aln.contigs()[tid].0, aln.contigs()[tid].1);
                    let mut counter = ContigCounter::with_buffer(
                        std::mem::take(&mut buf),
                        len as usize,
                        opts.mode,
                        opts.filter.clone(),
                    );
                    let result = reader.for_each(tid, need_rg, |a| counter.push(a)).map(|()| {
                        let cov = counter.finish();
                        let report = ContigReport::new(name, len as u32, &cov, &opts.report, &mut scratch);
                        buf = cov.depth;
                        report
                    });
                    reader.release_reference();
                    if result.is_err() {
                        failed.store(true, Ordering::Relaxed);
                    }
                    if tx.send((k, result)).is_err() {
                        break;
                    }
                }
            });
            spawned?;
        }
        drop(tx);

        // Hand out reports in the order of `work` as they arrive.
        let mut pending = BTreeMap::new();
        let mut written = 0;
        for (k, result) in rx {
            pending.insert(k, result?);
            while let Some(r) = pending.remove(&written) {
                each(work[written], r)?;
                written += 1;
            }
        }
        if written != work.len() {
            return Err(anyhow!("internal error: {} of {} contigs written", written, work.len()));
        }
        Ok(())
    })
}

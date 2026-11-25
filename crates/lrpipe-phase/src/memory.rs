//! Memory-aware scheduling of chromosomes that are phased in parallel.
//!
//! Most of the memory of a chromosome (besides the DP tables, which have
//! their own budget) holds its reads with detected alleles: the reference
//! sequence plus about 1.7 (simulated 35x paired-end short reads, HG003
//! chr20) to 1.9 (ONT, HG003 chr20) times the size of these reads at the end
//! of reading, when grouping the reads by name no longer copies them. Later
//! steps need about as much (up to 2.9 times measured on HG002 chr1 ONT).
//!
//! Before reading, the size of a chromosome's reads is unknown, but it is
//! about proportional to its number of variants, so the peak memory per
//! variant of chromosomes that have been (partly) read is used to estimate
//! the peak of the next one. A chromosome is only started while the
//! estimates of all chromosomes in progress fit into the budget.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use crate::variants::ReadProgress;

/// Peak memory of a chromosome relative to the size of its detected reads
/// (1.7-2.9 measured, see the module documentation)
const READ_FACTOR: f64 = 3.0;

/// Memory reserved per thread for everything but the chromosomes in
/// progress and the DP tables
pub const THREAD_MEMORY: u64 = 96 << 20;

/// A chromosome's own reading progress is used for its estimate once this
/// fraction of its variants has been passed.
const MIN_PROGRESS: f64 = 0.05;

/// Physical memory, or the memory limit of the cgroup this process runs in
/// (e.g. a container or a cluster job) if that is lower.
pub fn system_memory() -> Option<u64> {
    // SAFETY: sysconf has no preconditions.
    let (pages, page_size) = unsafe { (libc::sysconf(libc::_SC_PHYS_PAGES), libc::sysconf(libc::_SC_PAGESIZE)) };
    let physical = (pages > 0 && page_size > 0).then(|| pages as u64 * page_size as u64);
    match (physical, cgroup_memory_limit()) {
        (Some(p), Some(c)) => Some(p.min(c)),
        (p, c) => p.or(c),
    }
}

/// The lowest memory limit of this process's cgroup and its ancestors
/// (cgroup v2 `memory.max` or v1 `memory.limit_in_bytes`).
fn cgroup_memory_limit() -> Option<u64> {
    let cgroups = std::fs::read_to_string("/proc/self/cgroup").ok()?;
    let mut limit: Option<u64> = None;
    for line in cgroups.lines() {
        // hierarchy-id:controllers:path
        let mut fields = line.splitn(3, ':');
        let (Some(id), Some(controllers), Some(path)) = (fields.next(), fields.next(), fields.next()) else {
            continue;
        };
        let (root, file) = if id == "0" && controllers.is_empty() {
            ("/sys/fs/cgroup", "memory.max")
        } else if controllers.split(',').any(|c| c == "memory") {
            ("/sys/fs/cgroup/memory", "memory.limit_in_bytes")
        } else {
            continue;
        };
        let mut dir = std::path::PathBuf::from(format!("{root}{path}"));
        loop {
            let value = std::fs::read_to_string(dir.join(file)).ok().and_then(|s| s.trim().parse::<u64>().ok());
            if let Some(v) = value {
                limit = Some(limit.map_or(v, |l| l.min(v)));
            }
            if dir.as_os_str().len() <= root.len() || !dir.pop() {
                break;
            }
        }
    }
    limit
}

/// Memory use of one chromosome, updated by the worker phasing it.
#[derive(Debug)]
pub struct ChromosomeUsage {
    /// Variants of the chromosome in the VCF
    variants: u64,
    pub progress: Arc<ReadProgress>,
    reference_bytes: AtomicU64,
    finished: AtomicBool,
}

impl ChromosomeUsage {
    pub fn new(variants: usize) -> Self {
        ChromosomeUsage {
            variants: variants as u64,
            progress: Arc::new(ReadProgress::default()),
            reference_bytes: AtomicU64::new(0),
            finished: AtomicBool::new(false),
        }
    }

    /// Memory of the reference sequence (all copies)
    pub fn set_reference_bytes(&self, bytes: usize) {
        self.reference_bytes.store(bytes as u64, Ordering::Relaxed);
    }

    /// Called when the chromosome's memory has been freed.
    fn finish(&self) {
        self.finished.store(true, Ordering::Release);
    }

    fn is_finished(&self) -> bool {
        self.finished.load(Ordering::Acquire)
    }

    fn reference(&self) -> f64 {
        self.reference_bytes.load(Ordering::Relaxed) as f64
    }

    /// Memory needed at least, from the reads detected so far.
    fn current(&self) -> f64 {
        self.reference() + READ_FACTOR * self.progress.bytes() as f64
    }

    /// Estimated peak memory from the chromosome's own progress, once it is
    /// far enough to tell (None before reading has started).
    fn peak(&self) -> Option<f64> {
        if self.progress.loaded() {
            return Some(self.current());
        }
        let (done, total) = self.progress.variants();
        if total == 0 || (done as f64) < total as f64 * MIN_PROGRESS {
            return None;
        }
        // Reads of earlier samples (pedigrees) are complete.
        let (before, bytes) = self.progress.bytes_split();
        let reads = before as f64 + bytes as f64 * total as f64 / done as f64;
        Some(self.reference() + READ_FACTOR * reads)
    }
}

/// Marks a chromosome as finished when dropped, also if phasing it panics.
pub struct FinishOnDrop(pub Arc<ChromosomeUsage>);

impl Drop for FinishOnDrop {
    fn drop(&mut self) {
        self.0.finish();
    }
}

/// Decides when the next chromosome may start. Used by the thread that hands
/// out chromosomes, in VCF order.
#[derive(Debug)]
pub struct MemoryScheduler {
    budget: f64,
    running: Vec<Arc<ChromosomeUsage>>,
    /// Peak memory and variants of the finished chromosomes
    finished_bytes: f64,
    finished_variants: f64,
}

impl MemoryScheduler {
    pub fn new(budget_bytes: u64) -> Self {
        MemoryScheduler {
            budget: budget_bytes as f64,
            running: Vec::new(),
            finished_bytes: 0.0,
            finished_variants: 0.0,
        }
    }

    /// Whether a chromosome with `variants` variants fits next to the
    /// chromosomes in progress. Always true if none is in progress.
    pub fn fits(&mut self, variants: usize) -> bool {
        for c in self.running.extract_if(.., |c| c.is_finished()) {
            self.finished_bytes += c.peak().unwrap_or_else(|| c.current());
            self.finished_variants += c.variants as f64;
        }
        if self.running.is_empty() || variants == 0 {
            return true;
        }
        // Peak memory per variant of the chromosomes whose peak can be estimated
        let (mut bytes, mut n) = (self.finished_bytes, self.finished_variants);
        for c in &self.running {
            if let Some(estimate) = c.peak() {
                bytes += estimate;
                n += c.variants as f64;
            }
        }
        if n == 0.0 {
            return false;
        }
        let per_variant = bytes / n;
        let in_progress: f64 = self
            .running
            .iter()
            .map(|c| c.peak().unwrap_or_else(|| (per_variant * c.variants as f64).max(c.current())))
            .sum();
        in_progress + per_variant * variants as f64 <= self.budget
    }

    pub fn start(&mut self, usage: Arc<ChromosomeUsage>) {
        self.running.push(usage);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_chromosome_always_fits() {
        let mut s = MemoryScheduler::new(0);
        assert!(s.fits(1000));
    }

    #[test]
    fn waits_for_an_estimate_then_admits_what_fits() {
        let mut s = MemoryScheduler::new(6_000_000);
        let a = Arc::new(ChromosomeUsage::new(1000));
        s.start(Arc::clone(&a));
        // nothing known about memory per variant yet
        assert!(!s.fits(1000));
        assert!(s.fits(0));
        // 10% of the variants passed with 100 kB of reads: 1 MB of reads in
        // the end, so 3 MB estimated for each chromosome of 1000 variants
        a.progress.start(1000);
        a.progress.set_for_test(100, 100_000);
        assert!(s.fits(1000));
        assert!(!s.fits(1001));
        // reading finishes with 1 MB of reads
        a.progress.set_for_test(1000, 1_000_000);
        a.finish();
        assert!(s.fits(1_000_000));
        // the finished chromosome's estimate is kept for later ones
        let b = Arc::new(ChromosomeUsage::new(1000));
        s.start(Arc::clone(&b));
        assert!(s.fits(1000));
        assert!(!s.fits(1001));
    }

    #[test]
    fn loaded_reads_are_not_extrapolated() {
        let mut s = MemoryScheduler::new(6_000_000);
        let a = Arc::new(ChromosomeUsage::new(1000));
        s.start(Arc::clone(&a));
        a.progress.start(1000);
        a.progress.set_for_test(500, 1_000_000);
        // reading: 2 MB of reads projected, 6 MB peak; nothing else fits
        assert!(!s.fits(1));
        // loaded although only half the variants were passed: 3 MB, and
        // 3 MB per 1000 variants for the next
        a.progress.set_loaded_for_test();
        assert!(s.fits(1000));
        assert!(!s.fits(1001));
    }

    #[test]
    fn earlier_samples_are_not_extrapolated() {
        let a = ChromosomeUsage::new(1000);
        a.progress.start(1000);
        a.progress.set_for_test(1000, 1_000_000);
        // second sample: 10% of it read with 100 kB more
        a.progress.start(1000);
        a.progress.set_for_test(100, 1_100_000);
        assert_eq!(a.peak(), Some(READ_FACTOR * 2_000_000.0));
    }
}

//! Batch-level scheduling of CPU inference.
//!
//! The original runs one process per chunk, so a stage takes as long as its biggest chunk (on
//! HG003 chr20 one pileup chunk holds 124k of 333k candidates). Here each job's inputs are built on
//! the worker pool and every batch then becomes its own task. A batch is still the same
//! candidates the original puts in it, scored on one thread, so the probabilities do not
//! change; only which thread computes them does.

use std::cell::RefCell;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, mpsc};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use anyhow::{Result, anyhow};

/// CPU time spent building inputs and running the network, summed over threads.
#[derive(Debug, Default)]
pub(crate) struct StageTimes {
    build_ns: AtomicU64,
    score_ns: AtomicU64,
}

impl StageTimes {
    fn add(counter: &AtomicU64, since: Instant) {
        counter.fetch_add(since.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }

    /// `(building inputs, inference)` CPU time.
    pub(crate) fn totals(&self) -> (Duration, Duration) {
        (
            Duration::from_nanos(self.build_ns.load(Ordering::Relaxed)),
            Duration::from_nanos(self.score_ns.load(Ordering::Relaxed)),
        )
    }
}

/// How a job's flat input splits into batches.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Batching {
    /// Values per candidate tensor.
    pub(crate) tensor_len: usize,
    /// Candidates per batch.
    pub(crate) batch: usize,
}

/// A job whose batches are still being scored.
struct Pending<M, T> {
    meta: Mutex<Option<M>>,
    inputs: Vec<T>,
    probs: Mutex<Vec<Vec<f32>>>,
    remaining: AtomicUsize,
}

/// Runs `build` for every job, scores each job's inputs in batches (`shape`) with `score`, and
/// hands each job's probabilities, in candidate order, to `finish` as soon as its last batch is
/// done. `jobs` is drawn on the calling thread, so it may block until more jobs are ready.
///
/// # Errors
///
/// Returns the first error any step reports; remaining work is skipped.
pub(crate) fn score_in_batches<'a, J, M, T>(
    pool: &rayon::ThreadPool,
    jobs: impl IntoIterator<Item = &'a J>,
    shape: Batching,
    times: &StageTimes,
    build: impl Fn(&J) -> Result<(M, Vec<T>)> + Sync,
    score: impl Fn(&[T]) -> Result<Vec<f32>> + Sync,
    finish: impl Fn(&J, M, Vec<f32>) -> Result<()> + Sync,
) -> Result<()>
where
    J: Sync + 'a,
    M: Send,
    T: Send + Sync,
{
    let failure: Mutex<Option<anyhow::Error>> = Mutex::new(None);
    let failed = || failure.lock().expect("lock poisoned").is_some();
    let fail = |e: anyhow::Error| {
        failure.lock().expect("lock poisoned").get_or_insert(e);
    };
    let (build, score, finish) = (&build, &score, &finish);
    let Batching { tensor_len, batch } = shape;
    pool.in_place_scope(|s| {
        for job in jobs {
            s.spawn(move |s| {
                if failed() {
                    return;
                }
                let started = Instant::now();
                let built = build(job);
                StageTimes::add(&times.build_ns, started);
                let (meta, inputs) = match built {
                    Ok(b) => b,
                    Err(e) => return fail(e),
                };
                let n_batches = inputs.len().div_ceil(batch * tensor_len);
                if n_batches == 0 {
                    if let Err(e) = finish(job, meta, Vec::new()) {
                        fail(e);
                    }
                    return;
                }
                let pending = Arc::new(Pending {
                    meta: Mutex::new(Some(meta)),
                    inputs,
                    probs: Mutex::new(vec![Vec::new(); n_batches]),
                    remaining: AtomicUsize::new(n_batches),
                });
                for b in 0..n_batches {
                    let pending = Arc::clone(&pending);
                    s.spawn(move |_| {
                        if failed() {
                            return;
                        }
                        let start = b * batch * tensor_len;
                        let end = (start + batch * tensor_len).min(pending.inputs.len());
                        let started = Instant::now();
                        let scored = score(&pending.inputs[start..end]);
                        StageTimes::add(&times.score_ns, started);
                        match scored {
                            Ok(p) => pending.probs.lock().expect("lock poisoned")[b] = p,
                            Err(e) => return fail(e),
                        }
                        if pending.remaining.fetch_sub(1, Ordering::AcqRel) == 1 {
                            let probs = std::mem::take(&mut *pending.probs.lock().expect("lock poisoned")).concat();
                            let meta = pending.meta.lock().expect("lock poisoned").take();
                            let meta = meta.expect("the last batch finishes its job once");
                            if let Err(e) = finish(job, meta, probs) {
                                fail(e);
                            }
                        }
                    });
                }
            });
        }
    });
    match failure.into_inner().map_err(|_| anyhow!("lock poisoned"))? {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

/// One model per worker thread, loaded on first use and kept for later batches and stages.
pub(crate) struct ThreadModel<M: 'static> {
    cell: &'static std::thread::LocalKey<RefCell<Option<(PathBuf, M)>>>,
}

impl<M: 'static> ThreadModel<M> {
    pub(crate) const fn new(cell: &'static std::thread::LocalKey<RefCell<Option<(PathBuf, M)>>>) -> Self {
        Self { cell }
    }

    /// Runs `f` with this thread's model for `path`, loading it with `load` if needed.
    pub(crate) fn with<R>(
        &self,
        path: &Path,
        load: impl FnOnce(&Path) -> Result<M>,
        f: impl FnOnce(&mut M) -> Result<R>,
    ) -> Result<R> {
        self.cell.with(|cell| {
            let mut slot = cell.borrow_mut();
            if slot.as_ref().is_none_or(|(p, _)| p != path) {
                *slot = Some((path.to_owned(), load(path)?));
            }
            let (_, model) = slot.as_mut().expect("loaded above");
            f(model)
        })
    }
}

/// Gives the heap memory the allocator keeps after it is freed (glibc's arenas, macOS's malloc
/// zones) back to the system every few seconds, while it lives.
///
/// glibc keeps what a thread frees in that thread's arena for later allocations. In a
/// whole-genome run's phasing and full alignment that came to 10-17 GB, next to 2-5 GB in use,
/// and on a 31 GB machine pushed the BAM out of the page cache and into swap; with
/// `malloc_trim` every few seconds the peak fell from 16.7 to 11.3-11.7 GB and the swapping by
/// three quarters, with run times within their run-to-run spread. Pileup calling, which
/// allocates the most and is not where the peak is, is left out.
#[derive(Debug)]
pub(crate) struct MemoryTrimmer {
    stop: Option<mpsc::Sender<()>>,
    thread: Option<JoinHandle<()>>,
}

impl MemoryTrimmer {
    /// Trims every `every` until dropped.
    pub(crate) fn start(every: Duration) -> Self {
        let (stop, stopped) = mpsc::channel::<()>();
        let thread = std::thread::spawn(move || {
            while let Err(mpsc::RecvTimeoutError::Timeout) = stopped.recv_timeout(every) {
                trim_heap();
            }
        });
        Self { stop: Some(stop), thread: Some(thread) }
    }
}

impl Drop for MemoryTrimmer {
    fn drop(&mut self) {
        // Closing the channel stops the thread at its next wait.
        drop(self.stop.take());
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn trim_heap() {
    // SAFETY: `malloc_trim` has no preconditions and is thread-safe.
    unsafe { libc::malloc_trim(0) };
}

/// macOS's counterpart: every malloc zone (a null zone) gives back what it can.
#[cfg(target_os = "macos")]
fn trim_heap() {
    unsafe extern "C" {
        fn malloc_zone_pressure_relief(zone: *mut libc::c_void, goal: libc::size_t) -> libc::size_t;
    }
    // SAFETY: a null zone means every zone, and a goal of 0 as much as possible; the call is
    // thread-safe.
    unsafe { malloc_zone_pressure_relief(std::ptr::null_mut(), 0) };
}

#[cfg(not(any(all(target_os = "linux", target_env = "gnu"), target_os = "macos")))]
fn trim_heap() {}

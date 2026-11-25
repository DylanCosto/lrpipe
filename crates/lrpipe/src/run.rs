//! Runs a plan's stages: each as soon as the stages it follows have finished, so independent
//! stages run side by side; a failure kills the stages still running and starts no others.
//!
//! A stage that succeeds leaves a stamp of its command and inputs (each one's path and
//! contents, or for large files their size and time). Run again, a stage whose stamp matches
//! and whose outputs and recorded indexes are unchanged (an output's contents, or for one over
//! [`HASHED_OUTPUT_BYTES`] its size, time and pieces of it) is skipped, so a stopped run resumes where it
//! stopped, and a stage reruns, with every stage after it, when its command or an input changes.
//!
//! The stamp leaves out what does not change a stage's results: its thread count (and the align
//! stage's batch sizes, which the memory budget sets), how its paths were spelt (`./site.toml` or
//! `/data/site.toml`), and where the output directory is (a path in it is stamped as one under
//! the output directory), so a run resumed on a node with other cores or memory, from another
//! directory, or in an output directory moved or copied with its files' times, does not redo what
//! is done (LRUN-07). No stage reads a path out of another's outputs (only the aligned BAM's
//! `@PG` line names the directory, as a record), so a moved directory's skipped stages leave
//! nothing pointing back to where it was for those that run. A stamp of lrpipe 0.1.5 or before,
//! which held the output directory's path, is taken as long as the directory is where it was,
//! and is written anew when its stage is skipped.
//!
//! A stage's scratch directory is emptied before it starts and removed once it succeeds.
//!
//! Stages that run together stay within a memory budget, by each stage's expected peak: its
//! plan's estimate, or the peak measured when it last succeeded in this output directory, with
//! [`MEASURED_HEADROOM`] for run-to-run variation, if that is higher. A stage that was killed or
//! failed has not reached its peak, so what it measured is not used (LRUN-02).

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::MetadataExt;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, mpsc};
use std::thread;
use std::time::{Duration, Instant, SystemTime};

use sha2::{Digest, Sha256};

use crate::manifest::{self, Manifest, RunStatus, StageRecord, StageStatus};
use crate::plan::{BACKGROUND_NICE, Layout, Priority, Stage, StageName};
use crate::process::{Exit, Running};
use crate::provenance::{self, hex, utc};

/// Which stages run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Rerun {
    /// Those not already done with the same command and inputs.
    #[default]
    Changed,
    /// All of them.
    All,
}

/// How a run goes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Settings {
    /// Which stages run.
    pub rerun: Rerun,
    /// GiB the stages running together may use: a stage waits while the expected peaks of
    /// the running ones and its own would exceed it, unless nothing else is running.
    pub memory_gib: f64,
}

/// Progress, reported as it happens.
#[derive(Clone, Debug, PartialEq)]
pub enum Event {
    /// A stage was already done.
    Skipped(StageName),
    /// A stage is ready but waits for memory: it expects this many GiB, and the stages
    /// running expect this many.
    Waiting(StageName, f64, f64),
    /// A stage started, as this process.
    Started(StageName, u32),
    /// A stage succeeded after this many seconds.
    Succeeded(StageName, f64),
    /// A stage failed, for this reason.
    Failed(StageName, String),
    /// A stage was killed because another failed or the run was interrupted.
    Cancelled(StageName),
    /// The run was interrupted: the running stages are being stopped.
    Interrupted,
    /// A stage killed [`KILL_WAIT`] ago still runs, as this process: the run ends without it, and
    /// must not be run again until it has ended.
    Abandoned(StageName, u32),
}

/// Why a run did not complete.
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    /// A file of the run's own could not be written or read.
    #[error("{context}")]
    Io {
        /// What was being done.
        context: String,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A stage failed.
    #[error("stage {stage} failed: {reason} (log: {})", log.display())]
    Stage {
        /// Which.
        stage: StageName,
        /// How.
        reason: String,
        /// Its log.
        log: PathBuf,
    },
    /// Another run is using the output directory.
    #[error("another lrpipe run is using {}", .0.display())]
    Busy(PathBuf),
    /// The output directory holds another sample's run.
    #[error(
        "{} holds the run of sample {owner}, not {sample}: give each sample its own --output (if it holds no \
         run of {owner}'s, remove {})",
        dir.display(),
        record.display()
    )]
    OtherSample {
        /// The output directory.
        dir: PathBuf,
        /// The sample whose run it holds.
        owner: String,
        /// This run's sample.
        sample: String,
        /// The file that names `owner`: `.lrpipe/sample`, or the manifest of a run from before it.
        record: PathBuf,
    },
    /// The output directory's `snv/` holds files no lrpipe run wrote, which the snv stage would
    /// replace.
    #[error(
        "{} holds files no lrpipe run wrote (an earlier Clair3 run's, say), which the snv stage would replace: \
         move them away or give another --output",
        .0.display()
    )]
    ForeignSnv(PathBuf),
    /// The output directory holds files no lrpipe run wrote where this run writes its own.
    #[error(
        "{} is in the way: no lrpipe run wrote it, and this run would replace it. Move it away or give another --output",
        .0.display()
    )]
    Foreign(PathBuf),
    /// The run was interrupted (by SIGINT or SIGTERM).
    #[error("interrupted")]
    Interrupted,
    /// A stage follows one the plan does not have.
    #[error("stage {stage} follows {missing}, which is not in the plan")]
    UnknownDependency {
        /// The stage.
        stage: StageName,
        /// What it follows.
        missing: StageName,
    },
}

fn io_error(context: impl Into<String>) -> impl FnOnce(io::Error) -> RunError {
    let context = context.into();
    move |source| RunError::Io { context, source }
}

/// A stage started and not yet reaped.
#[derive(Debug)]
struct Active {
    process: Arc<Running>,
    started: Instant,
    fingerprint: Fingerprint,
    cancelled: bool,
}

/// The output directory's lock, held by one run at a time while it lives.
#[derive(Debug)]
pub struct OutputLock(File);

impl Drop for OutputLock {
    fn drop(&mut self) {
        // A concurrent fork can hold a duplicate descriptor until exec. Release the lock
        // when its owner finishes, even if that duplicate has not been closed yet.
        let _ = self.0.unlock();
    }
}

/// Takes the output directory for this run: two runs in one directory would write the same
/// files. Taken before the run writes anything there, the PAR BED its stages read among them,
/// so a run turned away changes nothing of the one that holds it (LRP-019).
///
/// The directory is one sample's: `snv/`, `logs/`, the manifest and the stage records are not
/// named after the sample, so another sample's run would replace them (LRUN-04). A directory
/// whose `.lrpipe/sample`, or else manifest, names another sample is refused with nothing of
/// its run changed; any other is claimed for this run's sample.
///
/// # Errors
///
/// [`RunError::Busy`] if another run holds it, [`RunError::OtherSample`] if it holds another
/// sample's run, or if `.lrpipe/lock` or `.lrpipe/sample` cannot be made.
pub fn lock(layout: &Layout) -> Result<OutputLock, RunError> {
    let dir = layout.root().join(".lrpipe");
    fs::create_dir_all(&dir).map_err(io_error(format!("creating {}", dir.display())))?;
    let lock_path = dir.join("lock");
    let lock = File::create(&lock_path).map_err(io_error(format!("creating {}", lock_path.display())))?;
    lock.try_lock().map_err(|e| match e {
        fs::TryLockError::WouldBlock => RunError::Busy(layout.root().to_path_buf()),
        fs::TryLockError::Error(source) => RunError::Io { context: format!("locking {}", lock_path.display()), source },
    })?;
    // Asked before `claim` records this run's sample, so a run turned away below leaves nothing
    // that would make the next one take the directory for lrpipe's.
    let earlier_run = earlier_run(layout);
    if !earlier_run {
        refuse_foreign(layout)?;
    }
    claim(layout)?;
    tidy(layout, earlier_run)?;
    Ok(OutputLock(lock))
}

/// What [`lock`] would refuse the output directory for, checked without writing anything, as a
/// dry run checks it: another run using it, another sample's run in it, or files in its `snv/`
/// that no lrpipe run wrote.
///
/// # Errors
///
/// [`RunError::Busy`], [`RunError::OtherSample`], [`RunError::ForeignSnv`] or [`RunError::Foreign`], as [`lock`] would
/// return them, or if `.lrpipe/lock` or `.lrpipe/sample` cannot be read.
pub fn check(layout: &Layout) -> Result<(), RunError> {
    let lock_path = layout.root().join(".lrpipe").join("lock");
    match File::open(&lock_path) {
        // Shared, and released when the file is dropped, at once: a run's exclusive lock refuses
        // it, and it refuses no run or other check (an exclusive lock on a file open read-only
        // also fails on NFS, where flock is fcntl's).
        Ok(file) => match file.try_lock_shared() {
            Ok(()) => {}
            Err(fs::TryLockError::WouldBlock) => return Err(RunError::Busy(layout.root().to_path_buf())),
            Err(fs::TryLockError::Error(source)) => {
                return Err(RunError::Io { context: format!("locking {}", lock_path.display()), source });
            }
        },
        Err(e) if e.kind() == io::ErrorKind::NotFound => {}
        Err(source) => return Err(RunError::Io { context: format!("opening {}", lock_path.display()), source }),
    }
    if !earlier_run(layout) {
        refuse_foreign(layout)?;
    }
    owner(layout).map(drop)
}

/// Whether the directory held an lrpipe run before this one: its `.lrpipe/sample`, or a manifest
/// as lrpipe writes them. Any other `manifest.json` is someone else's.
fn earlier_run(layout: &Layout) -> bool {
    layout.owner().exists() || manifest::is_lrpipe(&layout.manifest())
}

/// Refuses a directory that never held an lrpipe run but holds files where the run writes its
/// own, which it would replace: in `snv/` (or `snv/tmp/`) without lrpipe-snv's marker (LRUN-09),
/// in the other stages' directories, `logs/` and `qc/`, or the merged VCF, its index or
/// `manifest.json`. A directory with no file in it, only empty directories, is no one's, and
/// other files beside the run's (a project's own) are left alone.
fn refuse_foreign(layout: &Layout) -> Result<(), RunError> {
    let snv = layout.snv_dir();
    for dir in [snv.join("tmp"), snv] {
        if !dir.join(SNV_MARKER).exists() && holds_files(&dir) {
            return Err(RunError::ForeignSnv(dir));
        }
    }
    let root = layout.root();
    let merged = layout.merged_vcf();
    let mut index = merged.clone().into_os_string();
    index.push(".tbi");
    for file in [merged, PathBuf::from(index), layout.manifest()] {
        if file.symlink_metadata().is_ok() {
            return Err(RunError::Foreign(file));
        }
    }
    let mut dirs: Vec<PathBuf> = [
        layout.log(StageName::Snv),
        layout.qc_json(),
        layout.sv_vcf(),
        layout.depth_windows(),
        layout.repeats_vcf(),
        layout.bedmethyl(),
        layout.cnv_vcf(),
        layout.haplotags(),
        layout.aligned_bam(),
        layout.merged_reads(),
    ]
    .iter()
    .filter_map(|path| path.parent().map(Path::to_path_buf))
    .filter(|dir| dir != root)
    .collect();
    dirs.sort();
    dirs.dedup();
    for dir in dirs {
        if holds_files(&dir) {
            return Err(RunError::Foreign(dir));
        }
    }
    Ok(())
}

/// Whether `dir` holds a file (or link, or anything unreadable) at any depth.
fn holds_files(dir: &Path) -> bool {
    let Ok(entries) = fs::read_dir(dir) else { return false };
    entries
        .into_iter()
        .any(|entry| entry.map_or(true, |e| e.file_type().map_or(true, |t| !t.is_dir() || holds_files(&e.path()))))
}

/// lrpipe-snv's file saying it made what is beside it (its `OWNER_MARKER`).
const SNV_MARKER: &str = ".lrpipe-snv";

/// Readies the directory, which is this run's sample's, for this run:
/// - if it held an earlier lrpipe run (`earlier_run`: a `.lrpipe/sample` or manifest was there
///   before this run claimed it), marks that run's `snv/` and `snv/tmp/` as lrpipe-snv's, which
///   it then replaces when its stage reruns. lrpipe made them, but lrpipe-snv before 0.1.6 did
///   not mark them, and refuses a `tmp/` or results it did not mark (LRUN-09). A `snv/` in a
///   directory that held no lrpipe run is never marked: [`refuse_foreign`] turned the run
///   away if it held anything;
/// - removes the hidden files of a merge that was killed (`lrpipe merge` removes its own when it
///   fails, but not when it is killed, as stopping a run does).
fn tidy(layout: &Layout, earlier_run: bool) -> Result<(), RunError> {
    let snv = layout.snv_dir();
    for dir in [snv.clone(), snv.join("tmp")] {
        let marker = dir.join(SNV_MARKER);
        if earlier_run && dir.is_dir() && !marker.exists() {
            fs::write(&marker, "lrpipe-snv wrote what is here, and a run into this directory replaces it.\n")
                .map_err(io_error(format!("writing {}", marker.display())))?;
        }
    }
    let merged = layout.merged_vcf();
    let staged = crate::merge::staging_prefix(&merged);
    let entries = match fs::read_dir(layout.root()) {
        Ok(entries) => entries,
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(source) => return Err(RunError::Io { context: format!("listing {}", layout.root().display()), source }),
    };
    for entry in entries.flatten() {
        if entry.file_name().as_bytes().starts_with(staged.as_bytes()) {
            let path = entry.path();
            fs::remove_file(&path).map_err(io_error(format!("removing {}", path.display())))?;
        }
    }
    Ok(())
}

/// Refuses a directory that holds another sample's run; whether `.lrpipe/sample` records the
/// sample already.
fn owner(layout: &Layout) -> Result<bool, RunError> {
    let path = layout.owner();
    let (owner, recorded) = match fs::read_to_string(&path) {
        Ok(text) => (Some(text.trim_end_matches('\n').to_owned()), true),
        Err(e) if e.kind() == io::ErrorKind::NotFound => (manifest::earlier_sample(&layout.manifest()), false),
        Err(source) => return Err(RunError::Io { context: format!("reading {}", path.display()), source }),
    };
    if let Some(owner) = owner.filter(|o| o != layout.sample()) {
        return Err(RunError::OtherSample {
            dir: layout.root().to_path_buf(),
            owner,
            sample: layout.sample().to_owned(),
            record: if recorded { path } else { layout.manifest() },
        });
    }
    Ok(recorded)
}

/// Refuses a directory that holds another sample's run, and records this run's sample in it.
fn claim(layout: &Layout) -> Result<(), RunError> {
    if owner(layout)? {
        return Ok(());
    }
    let path = layout.owner();
    let tmp = path.with_extension("tmp");
    fs::write(&tmp, format!("{}\n", layout.sample()))
        .and_then(|()| fs::rename(&tmp, &path))
        .map_err(io_error(format!("writing {}", path.display())))
}

/// Runs `stages`, keeping `manifest`'s stage records and status current and writing it to the
/// layout's manifest after every change.
///
/// # Errors
///
/// If a stage fails or cannot start, `stop` is set (by a signal handler, say), or the run's
/// own files cannot be written. Stages still running are killed first, and the manifest
/// records why.
pub fn execute(
    stages: &[Stage],
    layout: &Layout,
    lock: &OutputLock,
    settings: Settings,
    stop: &AtomicBool,
    manifest: &mut Manifest,
    progress: impl FnMut(Event),
) -> Result<(), RunError> {
    let earlier = manifest::earlier_stages(&layout.manifest());
    execute_after(stages, layout, lock, settings, stop, manifest, &earlier, progress)
}

/// [`execute`] with `earlier` as the stage records of the runs before, in place of those of the
/// layout's manifest: a run that wrote the manifest in part before (the align stage alone, before
/// the others were planned) keeps the measurements of the run before that.
///
/// A stage `manifest` already records as succeeded, which this run ran before planning the others
/// (the align stage), keeps its record as it is, and is not run or reported again, if it is still
/// done: it is recorded as what this run did, not as skipped.
///
/// # Errors
///
/// As [`execute`].
#[expect(clippy::too_many_arguments, reason = "execute's, and the records it would read")]
pub fn execute_after(
    stages: &[Stage],
    layout: &Layout,
    _lock: &OutputLock,
    settings: Settings,
    stop: &AtomicBool,
    manifest: &mut Manifest,
    earlier: &[StageRecord],
    mut progress: impl FnMut(Event),
) -> Result<(), RunError> {
    let index: HashMap<StageName, usize> = stages.iter().enumerate().map(|(i, s)| (s.name, i)).collect();
    for stage in stages {
        if let Some(&missing) = stage.after.iter().find(|a| !index.contains_key(a)) {
            return Err(RunError::UnknownDependency { stage: stage.name, missing });
        }
    }
    let logs = layout.root().join("logs");
    fs::create_dir_all(&logs).map_err(io_error(format!("creating {}", logs.display())))?;
    let expected: Vec<f64> = stages.iter().map(|s| expected_peak(s, earlier)).collect();
    let mut waiting = vec![false; stages.len()];
    let ran: Vec<StageRecord> = std::mem::take(&mut manifest.stages);
    manifest.stages = stages
        .iter()
        .map(|s| {
            let this_run = ran.iter().find(|r| r.name == s.name && r.status == StageStatus::Succeeded);
            match this_run {
                Some(record) if is_done(s, layout) => record.clone(),
                _ => StageRecord::pending(s, layout.log(s.name)),
            }
        })
        .collect();
    manifest.status = RunStatus::Running;
    let write = |m: &Manifest| m.write(&layout.manifest()).map_err(io_error("writing the manifest"));
    let run_started = Instant::now();

    let (tx, rx) = mpsc::channel::<(usize, io::Result<Exit>)>();
    let mut active: Vec<Option<Active>> = stages.iter().map(|_| None).collect();
    let mut failure: Option<RunError> = None;
    // When the running stages were last killed, for giving up on one that does not end.
    let mut killed: Option<Instant> = None;
    'run: loop {
        // Start (or skip) every stage whose predecessors are done, until none is left.
        let mut changed = true;
        while changed && failure.is_none() {
            changed = false;
            for (i, stage) in stages.iter().enumerate() {
                // No stage starts once the run is asked to stop, not even one whose last
                // predecessor has just finished.
                if stop.load(Ordering::Relaxed) {
                    interrupt(&mut active, &mut failure, &mut progress);
                    killed = Some(Instant::now());
                    break;
                }
                let status = |n: &StageName| manifest.stages[index[n]].status;
                let ready = manifest.stages[i].status == StageStatus::Pending
                    && stage.after.iter().all(|n| matches!(status(n), StageStatus::Skipped | StageStatus::Succeeded));
                if !ready {
                    continue;
                }
                let fingerprint = fingerprint(stage, layout);
                if settings.rerun == Rerun::Changed && up_to_date(stage, layout, &fingerprint) {
                    changed = true;
                    let record = &mut manifest.stages[i];
                    record.status = StageStatus::Skipped;
                    let done = |r: &&StageRecord| {
                        r.name == stage.name && matches!(r.status, StageStatus::Succeeded | StageStatus::Skipped)
                    };
                    if let Some(e) = earlier.iter().find(done) {
                        record.keep_measurements(e);
                    }
                    progress(Event::Skipped(stage.name));
                    continue;
                }
                let in_use: f64 = active.iter().zip(&expected).filter(|(a, _)| a.is_some()).map(|(_, &m)| m).sum();
                if in_use > 0.0 && in_use + expected[i] > settings.memory_gib {
                    if !waiting[i] {
                        waiting[i] = true;
                        progress(Event::Waiting(stage.name, expected[i], in_use));
                    }
                    continue;
                }
                changed = true;
                match start(stage, layout) {
                    Ok(process) => {
                        let process = Arc::new(process);
                        let waiter = {
                            let (process, tx) = (Arc::clone(&process), tx.clone());
                            move || {
                                // The receiver outlives every stage, so sending cannot fail.
                                let _ = tx.send((i, process.wait()));
                            }
                        };
                        thread::spawn(waiter);
                        let pid = process.id();
                        active[i] = Some(Active { process, started: Instant::now(), fingerprint, cancelled: false });
                        let record = &mut manifest.stages[i];
                        record.status = StageStatus::Running;
                        record.started = Some(utc(SystemTime::now()));
                        progress(Event::Started(stage.name, pid));
                    }
                    Err(e) => {
                        let reason = format!("could not start {}: {e}", stage.program.display());
                        fail(manifest, &mut active, &mut failure, stage, layout, reason, &mut progress);
                        if stage.name.fails_the_run() {
                            killed = Some(Instant::now());
                            break;
                        }
                    }
                }
            }
        }
        if let Err(e) = write(manifest) {
            return Err(abandon(manifest, &active, layout, e));
        }
        if active.iter().all(Option::is_none) {
            break;
        }

        let (i, result) = loop {
            match rx.recv_timeout(STOP_POLL) {
                Ok(done) => break done,
                Err(mpsc::RecvTimeoutError::Timeout) => {
                    if stop.load(Ordering::Relaxed) && failure.is_none() {
                        interrupt(&mut active, &mut failure, &mut progress);
                        killed = Some(Instant::now());
                    }
                    if killed.is_some_and(|t| t.elapsed() > KILL_WAIT) {
                        // Killed and still not reaped (blocked in the kernel on a hung file
                        // system, say): the run ends without it rather than never.
                        for (i, a) in active.iter_mut().enumerate() {
                            if let Some(a) = a.take() {
                                let record = &mut manifest.stages[i];
                                record.status = StageStatus::Cancelled;
                                record.error = Some(format!(
                                    "killed, but still running as process {} {} s later: do not run again until it \
                                     has ended",
                                    a.process.id(),
                                    KILL_WAIT.as_secs()
                                ));
                                progress(Event::Abandoned(record.name, a.process.id()));
                            }
                        }
                        break 'run;
                    }
                }
                Err(mpsc::RecvTimeoutError::Disconnected) => {
                    unreachable!("a stage is running, so its waiter still holds a sender")
                }
            }
        };
        let stage = &stages[i];
        let done = active[i].take().expect("only running stages report");
        // A stage that died of SIGINT or SIGTERM was stopped, as the run was (SLURM's scancel
        // signals every process of the job), not failed; as is one that ends once the run is
        // asked to stop, before the stop is seen here.
        let stopped = stop.load(Ordering::Relaxed)
            || result.as_ref().is_ok_and(|e| matches!(e.status.signal(), Some(libc::SIGINT | libc::SIGTERM)));
        let record = &mut manifest.stages[i];
        record.wall_seconds = Some(round(done.started.elapsed().as_secs_f64()));
        let outcome = result.map_err(|e| format!("waiting for it failed: {e}")).and_then(|exit| {
            record.user_seconds = Some(round(exit.usage.user.as_secs_f64()));
            record.system_seconds = Some(round(exit.usage.system.as_secs_f64()));
            record.max_rss_gib = Some(round(exit.usage.max_rss_kib as f64 / 1_048_576.0));
            if exit.status.success() { Ok(()) } else { Err(describe(exit.status)) }
        });
        let outcome = outcome.and_then(|()| completed_stamp(stage, &done.fingerprint.current));
        if outcome.is_err() {
            clear_unfinished(stage);
        }
        match outcome {
            Ok(completed) => {
                if let Some(scratch) = &stage.scratch
                    && let Err(e) = remove_dir(scratch)
                {
                    let e = RunError::Io { context: format!("removing {}", scratch.display()), source: e };
                    return Err(abandon(manifest, &active, layout, e));
                }
                let stamp = layout.stamp(stage.name);
                if let Err(e) = fs::write(&stamp, &completed) {
                    record.status = StageStatus::Failed;
                    record.error = Some(format!("its stamp {} could not be written: {e}", stamp.display()));
                    let e = RunError::Io { context: format!("writing {}", stamp.display()), source: e };
                    return Err(abandon(manifest, &active, layout, e));
                }
                record.status = StageStatus::Succeeded;
                progress(Event::Succeeded(stage.name, record.wall_seconds.unwrap_or_default()));
            }
            Err(_) if done.cancelled => {
                record.status = StageStatus::Cancelled;
                record.error = Some(match failure {
                    Some(RunError::Interrupted) => "killed when the run was interrupted".into(),
                    _ => "killed after another stage failed".into(),
                });
                progress(Event::Cancelled(stage.name));
            }
            Err(reason) if stopped => {
                record.status = StageStatus::Cancelled;
                record.error = Some(format!("stopped when the run was interrupted ({reason})"));
                progress(Event::Cancelled(stage.name));
                if failure.is_none() {
                    interrupt(&mut active, &mut failure, &mut progress);
                    killed = Some(Instant::now());
                }
            }
            Err(reason) => {
                if failure.is_none() && stage.name.fails_the_run() {
                    killed = Some(Instant::now());
                }
                fail(manifest, &mut active, &mut failure, stage, layout, reason, &mut progress);
            }
        }
    }

    for record in manifest.stages.iter_mut().filter(|r| r.status == StageStatus::Pending) {
        record.status = StageStatus::NotRun;
    }
    manifest.status = match failure {
        None => RunStatus::Complete,
        Some(RunError::Interrupted) => RunStatus::Interrupted,
        Some(_) => RunStatus::Failed,
    };
    manifest.finished = Some(utc(SystemTime::now()));
    manifest.wall_seconds = Some(round(run_started.elapsed().as_secs_f64()));
    write(manifest)?;
    failure.map_or(Ok(()), Err)
}

/// How often a run waiting for its stages checks whether it was asked to stop.
const STOP_POLL: Duration = Duration::from_millis(200);

/// How long a run waits for a stage it killed to end before it records the stage as stopped and
/// ends without it.
const KILL_WAIT: Duration = Duration::from_secs(10);

/// Ends a run that cannot write its own files (a full disk, a directory made read-only): kills
/// the stages still running, and writes what the manifest can say of the run if it can, rather
/// than leave it `running` with stages that outlive the run. Returns `error`.
fn abandon(manifest: &mut Manifest, active: &[Option<Active>], layout: &Layout, error: RunError) -> RunError {
    for a in active.iter().flatten() {
        a.process.kill();
    }
    for record in &mut manifest.stages {
        match record.status {
            StageStatus::Running => {
                record.status = StageStatus::Cancelled;
                record.error = Some(format!("killed: the run could not go on ({error})"));
            }
            StageStatus::Pending => record.status = StageStatus::NotRun,
            _ => {}
        }
    }
    manifest.status = RunStatus::Failed;
    manifest.finished = Some(utc(SystemTime::now()));
    let _ = manifest.write(&layout.manifest());
    error
}

/// Records the run as interrupted and kills every stage still running.
fn interrupt(active: &mut [Option<Active>], failure: &mut Option<RunError>, progress: &mut impl FnMut(Event)) {
    progress(Event::Interrupted);
    *failure = Some(RunError::Interrupted);
    for a in active.iter_mut().flatten() {
        a.cancelled = true;
        a.process.kill();
    }
}

/// What a measured peak is scaled by for a stage's next run: lrpipe-snv's peak moves by 2 GiB
/// in 13 between runs of the same sample.
pub const MEASURED_HEADROOM: f64 = 1.15;

/// GiB `stage` is expected to use at its peak: its plan's estimate, or more if it measured more
/// when it last ran here. A peak of a run that succeeded (or of the one a skipped record carries)
/// is scaled by [`MEASURED_HEADROOM`]. One of a stage that was killed or failed, or is still
/// marked running, is only how far it got, which its next run passes: it raises the estimate
/// but never lowers it (LRUN-02).
fn expected_peak(stage: &Stage, earlier: &[StageRecord]) -> f64 {
    let Some(record) = earlier.iter().find(|r| r.name == stage.name) else { return stage.memory_gib };
    let Some(measured) = record.max_rss_gib else { return stage.memory_gib };
    let peak = match record.status {
        StageStatus::Succeeded | StageStatus::Skipped => measured * MEASURED_HEADROOM,
        _ => measured,
    };
    stage.memory_gib.max(peak)
}

/// Records `stage`'s failure and kills every stage still running, unless a failure was
/// already recorded (then only the record changes) or the stage's failure does not fail the run
/// ([`StageName::fails_the_run`]).
fn fail(
    manifest: &mut Manifest,
    active: &mut [Option<Active>],
    failure: &mut Option<RunError>,
    stage: &Stage,
    layout: &Layout,
    reason: String,
    progress: &mut impl FnMut(Event),
) {
    let record = manifest.stage_mut(stage.name);
    record.status = StageStatus::Failed;
    record.error = Some(reason.clone());
    progress(Event::Failed(stage.name, reason.clone()));
    if failure.is_none() && stage.name.fails_the_run() {
        *failure = Some(RunError::Stage { stage: stage.name, reason, log: layout.log(stage.name) });
        for a in active.iter_mut().flatten() {
            a.cancelled = true;
            a.process.kill();
        }
    }
}

/// Opens the stage's log, makes its output directories and spawns it.
fn start(stage: &Stage, layout: &Layout) -> io::Result<Running> {
    let _ = fs::remove_file(layout.stamp(stage.name));
    for dir in stage.outputs.iter().filter_map(|o| o.parent()) {
        fs::create_dir_all(dir)?;
    }
    // A killed run's scratch files are of no use to this one.
    if let Some(scratch) = &stage.scratch {
        remove_dir(scratch)?;
        fs::create_dir_all(scratch)?;
    }
    let log = File::create(layout.log(stage.name))?;
    // An earlier run's progress is not this one's.
    let progress = layout.progress(stage.name);
    let _ = fs::remove_file(&progress);
    if let Some(dir) = progress.parent() {
        fs::create_dir_all(dir)?;
    }
    let nice = match stage.priority {
        Priority::Critical => 0,
        Priority::Background => BACKGROUND_NICE,
    };
    Running::spawn(
        Command::new(&stage.program)
            .args(&stage.args)
            .env(lrpipe_reads::progress::ENV, progress)
            .stdin(Stdio::null())
            .stdout(log.try_clone()?)
            .stderr(log),
        nice,
    )
}

/// Removes what a stage that failed or was stopped leaves that is of no use: its outputs it had
/// opened but not written to (0 bytes, which no output of a finished stage is), and its scratch
/// directory, which its next run empties anyway. Its other outputs are left, as a failed tool
/// leaves them: its next run replaces them, and a stage is done only once its stamp is written.
fn clear_unfinished(stage: &Stage) {
    for output in &stage.outputs {
        if fs::symlink_metadata(output).is_ok_and(|m| m.is_file() && m.len() == 0) {
            let _ = fs::remove_file(output);
        }
    }
    if let Some(scratch) = &stage.scratch {
        let _ = remove_dir(scratch);
    }
}

/// Removes `dir` and everything in it, if it is there.
fn remove_dir(dir: &Path) -> io::Result<()> {
    match fs::remove_dir_all(dir) {
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        other => other,
    }
}

fn describe(status: ExitStatus) -> String {
    match (status.code(), status.signal()) {
        (Some(code), _) => format!("exit code {code}"),
        (None, Some(signal)) => format!("killed by signal {signal}"),
        (None, None) => status.to_string(),
    }
}

fn round(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

/// Inputs up to this size (models, catalogues, indexes) are fingerprinted by their contents,
/// so one replaced by a file of the same size and time still reruns its stages; larger ones (the
/// reads, the reference) by their size and modification time.
const HASHED_INPUT_BYTES: u64 = 64 << 20;

/// Options that change how a stage runs but not what it writes, left out of its fingerprint
/// with their values: the thread counts of lrpipe's tools, which a run takes from the node's
/// cores by default.
const EXECUTION_OPTIONS: [&str; 1] = ["--threads"];

/// The same for the align stage (bonefish): its threads, from the node's cores, and its batch
/// sizes, from the memory, none of which changes its output ([`crate::plan::align_stage`]).
const ALIGN_EXECUTION_OPTIONS: [&str; 3] = ["-t", "-K", "--sort-mem"];

/// A stage's fingerprint, and as lrpipe 0.1.5 and before took it.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Fingerprint {
    /// With the paths in the output directory as paths under it: what a stamp holds.
    current: String,
    /// With every path as [`provenance::canonical`] spells it, the output directory's included.
    legacy: String,
}

/// What a path under the output directory is stamped as, before its path under it: no absolute
/// path starts so.
const IN_OUTPUT: &[u8] = b"\0output/";

/// The two hashes of a [`Fingerprint`], fed together.
struct Hashes {
    current: Sha256,
    legacy: Sha256,
    /// The output directory, as [`provenance::canonical`] spells paths in it.
    root: PathBuf,
}

impl Hashes {
    /// `path` as [`provenance::canonical`] spells it, and for the current hash under the output
    /// directory when it is in it.
    fn path(&mut self, path: &Path) {
        let spelt = provenance::canonical(path);
        self.legacy.update(spelt.as_os_str().as_bytes());
        match spelt.strip_prefix(&self.root) {
            Ok(under) if spelt.is_absolute() => {
                self.current.update(IN_OUTPUT);
                self.current.update(under.as_os_str().as_bytes());
            }
            _ => self.current.update(spelt.as_os_str().as_bytes()),
        }
    }

    fn bytes(&mut self, bytes: &[u8]) {
        self.current.update(bytes);
        self.legacy.update(bytes);
    }
}

/// The stage's command, but for [`EXECUTION_OPTIONS`] and with its paths [`provenance::canonical`]
/// and those in `layout`'s output directory taken under it, and the state of its program and
/// inputs.
fn fingerprint(stage: &Stage, layout: &Layout) -> Fingerprint {
    let root = fs::canonicalize(layout.root()).unwrap_or_else(|_| provenance::canonical(layout.root()));
    let mut h = Hashes { current: Sha256::new(), legacy: Sha256::new(), root };
    let program = provenance::resolve(&stage.program).unwrap_or_else(|| stage.program.clone());
    for path in std::iter::once(&program).chain(&stage.inputs) {
        h.path(path);
        h.bytes(
            match fs::metadata(path) {
                Ok(m) if m.is_file() && m.len() <= HASHED_INPUT_BYTES => match fs::read(path) {
                    Ok(bytes) => format!("\0{}\0{}\n", m.len(), hex(&Sha256::digest(&bytes))),
                    Err(_) => "\0unreadable\n".into(),
                },
                Ok(m) => format!("\0{}\0{}.{}\n", m.len(), m.mtime(), m.mtime_nsec()),
                Err(_) => "\0missing\n".into(),
            }
            .as_bytes(),
        );
    }
    let execution: &[&str] = if stage.name == StageName::Align { &ALIGN_EXECUTION_OPTIONS } else { &EXECUTION_OPTIONS };
    let mut args = stage.args.iter();
    while let Some(arg) = args.next() {
        if execution.iter().any(|o| arg == *o) {
            args.next();
            continue;
        }
        if execution.iter().any(|o| arg.as_bytes().strip_prefix(o.as_bytes()).is_some_and(|v| v.starts_with(b"="))) {
            continue;
        }
        h.path(Path::new(arg));
        h.bytes(b"\0");
    }
    Fingerprint { current: hex(&h.current.finalize()), legacy: hex(&h.legacy.finalize()) }
}

/// Whether `stage` is done: it succeeded in `layout`'s directory with the same command and
/// inputs, and its outputs are all there, so a run with [`Rerun::Changed`] skips it. Nothing is
/// written (a dry run asks).
#[must_use]
pub fn is_done(stage: &Stage, layout: &Layout) -> bool {
    stamped(stage, layout, &fingerprint(stage, layout)).is_some()
}

/// Whether `stage` is done, by its `fingerprint`, as a run skipping it: a stamp of the legacy
/// fingerprint is written anew as the current one, so the directory can be moved from then on.
fn up_to_date(stage: &Stage, layout: &Layout, fingerprint: &Fingerprint) -> bool {
    match stamped(stage, layout, fingerprint) {
        Some(Stamp::Legacy) => {
            // Not written, it is taken again next time.
            if let Ok(completed) = completed_stamp(stage, &fingerprint.current) {
                let _ = fs::write(layout.stamp(stage.name), completed);
            }
            true
        }
        Some(Stamp::Current) => true,
        None => false,
    }
}

/// Which fingerprint a done stage's stamp holds.
enum Stamp {
    Current,
    Legacy,
}

/// Checks output indexes and adds their contents and paired data lengths to the completion stamp,
/// then the state of every output ([`output_state`]) as its last line.
fn completed_stamp(stage: &Stage, fingerprint: &str) -> Result<String, String> {
    let mut digest = None;
    for output in &stage.outputs {
        if !output.is_file() {
            return Err(format!("exited successfully without writing {}", output.display()));
        }
        let extension = output.extension().and_then(|part| part.to_str());
        let checked = match extension {
            Some("tbi") => pasta::tabix::fs::read(output).map(|_| ()),
            Some("bai") => pasta::bam::bai::fs::read(output).map(|_| ()),
            Some("csi") => pasta::csi::fs::read(output).map(|_| ()),
            _ => continue,
        };
        let error = |e| format!("invalid output index {}: {e}", output.display());
        checked.map_err(error)?;
        let hash = digest.get_or_insert_with(Sha256::new);
        let data = output.with_extension("");
        let index = File::open(output).map_err(error)?;
        hash.update(fs::metadata(data).map_err(error)?.len().to_le_bytes());
        hash.update(index.metadata().map_err(error)?.len().to_le_bytes());
        let mut reader = io::BufReader::new(index);
        let mut bytes = [0u8; 8192];
        loop {
            let n = reader.read(&mut bytes).map_err(error)?;
            if n == 0 {
                break;
            }
            hash.update(&bytes[..n]);
        }
    }
    let indexed = digest
        .map_or_else(|| fingerprint.to_owned(), |hash| format!("{fingerprint}\nindexes {}", hex(&hash.finalize())));
    let mut outputs = Sha256::new();
    for output in &stage.outputs {
        output_state(output, &mut outputs).map_err(|e| format!("reading {}: {e}", output.display()))?;
    }
    Ok(format!("{indexed}{OUTPUTS}{}", hex(&outputs.finalize())))
}

/// What starts a completion stamp's line of its outputs' state, which lrpipe 0.2.1 and before
/// did not write.
const OUTPUTS: &str = "\noutputs ";

/// Outputs up to this size (VCFs, gVCFs, BEDs, tables, the QC JSON: under 2 GB in all for a
/// whole genome) are stamped by their contents, which a resume reads again (about two seconds for a
/// whole genome on an NVMe disk); larger ones (the BAMs, 40 to 100 GB) by [`output_samples`].
const HASHED_OUTPUT_BYTES: u64 = 4 << 30;

/// How many evenly spaced pieces of a larger output, its first and last among them, are stamped.
const OUTPUT_SAMPLES: u64 = 33;

/// The size of each of those pieces.
const OUTPUT_SAMPLE_BYTES: u64 = 1 << 20;

/// The empty BGZF block that ends a BAM (SAM specification, 4.1.2).
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Feeds `hash` the state of `output`: its size and contents, or for one over
/// [`HASHED_OUTPUT_BYTES`] its size, modification time and [`OUTPUT_SAMPLES`] pieces, so a damaged
/// output (a flipped bit, a truncation) reruns its stage. A BAM over that size that does not end
/// with [`BGZF_EOF`] is an error: truncated.
fn output_state(output: &Path, hash: &mut Sha256) -> io::Result<()> {
    use std::os::unix::fs::FileExt;
    let file = File::open(output)?;
    let meta = file.metadata()?;
    let len = meta.len();
    hash.update(len.to_le_bytes());
    if len <= HASHED_OUTPUT_BYTES {
        io::copy(&mut io::BufReader::with_capacity(1 << 20, file), hash)?;
        return Ok(());
    }
    if output.extension().is_some_and(|e| e == "bam") {
        let mut tail = [0u8; BGZF_EOF.len()];
        file.read_exact_at(&mut tail, len - BGZF_EOF.len() as u64)?;
        if tail != BGZF_EOF {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "no BGZF end-of-file block: truncated"));
        }
    }
    hash.update(meta.mtime().to_le_bytes());
    hash.update(meta.mtime_nsec().to_le_bytes());
    let mut piece = vec![0u8; OUTPUT_SAMPLE_BYTES as usize];
    let last = len - OUTPUT_SAMPLE_BYTES;
    for i in 0..OUTPUT_SAMPLES {
        let offset = if i + 1 == OUTPUT_SAMPLES { last } else { last / (OUTPUT_SAMPLES - 1) * i };
        file.read_exact_at(&mut piece, offset)?;
        hash.update(&piece);
    }
    Ok(())
}

/// The stamp of `stage`, if its inputs and outputs still match. A stamp of lrpipe 0.2.1 or
/// before, without its outputs' state, is taken as [`Stamp::Legacy`] (and written anew) when the
/// rest matches: it holds nothing to check their contents against.
fn stamped(stage: &Stage, layout: &Layout, fingerprint: &Fingerprint) -> Option<Stamp> {
    let saved = fs::read_to_string(layout.stamp(stage.name)).ok()?;
    let inputs = saved.lines().next()?;
    if inputs != fingerprint.current && inputs != fingerprint.legacy {
        return None;
    }
    let completed = completed_stamp(stage, &fingerprint.current).ok()?;
    let without_outputs = completed.rsplit_once(OUTPUTS).map_or(completed.as_str(), |(indexed, _)| indexed);
    if saved == completed {
        Some(Stamp::Current)
    } else if saved == without_outputs || saved == fingerprint.current || saved == fingerprint.legacy {
        Some(Stamp::Legacy)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plan::Priority;

    #[test]
    fn a_dry_runs_check_takes_the_lock_shared() {
        // It took it exclusively: two dry runs, or a dry run and a run starting, refused each
        // other as busy.
        let root = std::env::temp_dir().join(format!("lrpipe-check-shared-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        let layout = Layout::new(&root, "S");
        drop(lock(&layout).expect("claimed"));
        let other = File::open(root.join(".lrpipe/lock")).expect("the lock file");
        other.try_lock_shared().expect("another check holds it");
        check(&layout).expect("checked beside another check");
        drop(other);
        let running = lock(&layout).expect("a run holds it");
        assert!(matches!(check(&layout), Err(RunError::Busy(_))), "a run's lock refuses it");
        drop(running);
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn a_large_output_is_stamped_by_its_size_time_and_pieces_and_a_bam_by_its_end() {
        use std::os::unix::fs::FileExt;
        let root = std::env::temp_dir().join(format!("lrpipe-large-output-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).expect("directory");
        // Sparse: it takes no space.
        let bam = root.join("S.bam");
        let file = File::create(&bam).expect("BAM");
        file.set_len(HASHED_OUTPUT_BYTES + 1).expect("sized");
        let state = || {
            let mut hash = Sha256::new();
            output_state(&bam, &mut hash).map(|()| hex(&hash.finalize()))
        };
        assert!(state().is_err(), "truncated: no end-of-file block");
        file.write_all_at(&BGZF_EOF, HASHED_OUTPUT_BYTES + 1).expect("ended");
        let time = file.metadata().and_then(|m| m.modified()).expect("a time");
        let intact = state().expect("ends with the block");
        assert_eq!(state().expect("ends with the block"), intact);
        file.write_all_at(&[1], 0).expect("first block damaged");
        file.set_modified(time).expect("time kept");
        assert_ne!(state().expect("ends with the block"), intact, "its first piece is stamped");
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn a_stamp_of_an_earlier_lrpipe_is_taken_and_written_anew() {
        let root = std::env::temp_dir().join(format!("lrpipe-legacy-stamp-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join(".lrpipe")).expect("scratch directory");
        let layout = Layout::new(&root, "S");
        let output = root.join("snv.txt");
        fs::write(&output, "done").expect("written");
        let stage = Stage {
            name: StageName::Snv,
            program: "sh".into(),
            args: vec!["--output".into(), output.clone().into()],
            inputs: Vec::new(),
            outputs: vec![output],
            after: Vec::new(),
            memory_gib: 1.0,
            priority: Priority::Critical,
            scratch: None,
        };
        let fingerprint = fingerprint(&stage, &layout);
        assert_ne!(fingerprint.current, fingerprint.legacy, "the legacy one holds the directory's path");
        fs::write(layout.stamp(StageName::Snv), &fingerprint.legacy).expect("written");
        assert!(is_done(&stage, &layout));
        assert_eq!(
            fs::read_to_string(layout.stamp(StageName::Snv)).expect("a stamp"),
            fingerprint.legacy,
            "asked only"
        );
        assert!(up_to_date(&stage, &layout, &fingerprint));
        assert_eq!(
            fs::read_to_string(layout.stamp(StageName::Snv)).expect("a stamp"),
            completed_stamp(&stage, &fingerprint.current).expect("complete")
        );
        fs::write(layout.stamp(StageName::Snv), "another").expect("written");
        assert!(!is_done(&stage, &layout));
        let _ = fs::remove_dir_all(root);
    }
}

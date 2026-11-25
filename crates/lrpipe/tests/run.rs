//! The scheduler, with shell commands standing in for the tools.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, MutexGuard, PoisonError};
use std::time::{Duration, Instant};

use lrpipe::manifest::{Manifest, RunStatus, SampleRecord, StageStatus};
use lrpipe::plan::{Layout, Priority, Stage, StageName};
use lrpipe::provenance::FileRecord;
use lrpipe::run::{Event, Rerun, RunError, Settings, execute, lock};

const CHANGED: Settings = Settings { rerun: Rerun::Changed, memory_gib: f64::INFINITY };
const ALL: Settings = Settings { rerun: Rerun::All, memory_gib: f64::INFINITY };

/// Never asked to stop.
static GO: AtomicBool = AtomicBool::new(false);

/// One test at a time. A run locks its directory with `flock`, and a stage another test forks
/// keeps a copy of that lock until it execs; on a loaded machine that is long enough for the
/// next run in this test's directory to find it busy.
static SERIAL: Mutex<()> = Mutex::new(());

fn serial() -> MutexGuard<'static, ()> {
    SERIAL.lock().unwrap_or_else(PoisonError::into_inner)
}

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir.canonicalize().expect("canonical temp dir")
}

fn manifest(input: &Path) -> Manifest {
    Manifest {
        lrpipe: "test".into(),
        command: Vec::new(),
        output: PathBuf::new(),
        memory_budget_gib: f64::INFINITY,
        started: String::new(),
        finished: None,
        status: RunStatus::Running,
        wall_seconds: None,
        error: None,
        sample: SampleRecord {
            name: "S".into(),
            platform: "hifi".into(),
            sex: "female".into(),
            sex_from_index: None,
            ploidy: None,
            reads: vec![FileRecord::of(input).expect("input exists")],
            unaligned_reads: Vec::new(),
            haplotypes: "none (input not haplotagged)".into(),
        },
        snv_model: None,
        resources: BTreeMap::new(),
        tools: BTreeMap::new(),
        stages: Vec::new(),
        outputs: BTreeMap::new(),
        warnings: Vec::new(),
    }
}

/// A stage running `script` in `sh`, which sees its outputs as `$1`, `$2`, ...
fn stage(name: StageName, script: &str, inputs: &[&Path], outputs: &[PathBuf], after: &[StageName]) -> Stage {
    let mut args = vec!["-c".into(), script.into(), "sh".into()];
    args.extend(outputs.iter().map(|o| o.clone().into_os_string()));
    Stage {
        name,
        program: "sh".into(),
        args,
        inputs: inputs.iter().map(|p| p.to_path_buf()).collect(),
        outputs: outputs.to_vec(),
        after: after.to_vec(),
        memory_gib: 1.0,
        priority: Priority::Critical,
        scratch: None,
    }
}

/// The pipeline's shape: two independent stages that sleep, then one after both. The first
/// writes something new each time, as the tools do (their headers name the date).
fn diamond(root: &Path, input: &Path) -> Vec<Stage> {
    let a = root.join("snv/a.txt");
    let b = root.join("reads/b.txt");
    vec![
        stage(StageName::Snv, "sleep 0.5; date +%s%N > $1", &[input], std::slice::from_ref(&a), &[]),
        stage(StageName::Reads, "sleep 0.5; echo b > $1", &[input], std::slice::from_ref(&b), &[]),
        stage(StageName::Cnv, "echo c > $1", &[&a, &b], &[root.join("cnv/c.txt")], &[StageName::Snv, StageName::Reads]),
    ]
}

fn statuses(m: &Manifest) -> Vec<StageStatus> {
    m.stages.iter().map(|s| s.status).collect()
}

#[test]
fn independent_stages_run_together_and_resume_skips_done_ones() {
    let _serial = serial();
    let root = temp_dir("diamond");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = diamond(&root, &input);

    let mut m = manifest(&input);
    let mut events = Vec::new();
    let t = Instant::now();
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |e| events.push(e)).expect("runs");
    assert!(
        t.elapsed() < Duration::from_millis(950),
        "the two 0.5 s stages ran one after the other: {:?}",
        t.elapsed()
    );
    assert_eq!(m.status, RunStatus::Complete);
    assert_eq!(statuses(&m), [StageStatus::Succeeded; 3]);
    assert!(m.stages.iter().all(|s| s.wall_seconds.is_some() && s.max_rss_gib.is_some()));
    let started_cnv = events.iter().position(|e| matches!(e, Event::Started(StageName::Cnv, _))).expect("cnv started");
    for first in [StageName::Snv, StageName::Reads] {
        let done = events.iter().position(|e| matches!(e, Event::Succeeded(n, _) if *n == first)).expect("done");
        assert!(done < started_cnv, "cnv started before {first} finished: {events:?}");
    }
    let written: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(layout.manifest()).expect("manifest written")).expect("json");
    assert_eq!(written["status"], "complete");
    assert_eq!(written["stages"][2]["name"], "cnv");

    // Nothing changed: everything is skipped.
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
    assert_eq!(statuses(&m), [StageStatus::Skipped; 3]);
    assert!(m.stages.iter().all(|s| s.wall_seconds.is_some()), "skipped stages keep the earlier run's times");

    // An output gone: that stage reruns, and so does the one reading it.
    fs::remove_file(root.join("snv/a.txt")).expect("removed");
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
    use StageStatus::{Skipped, Succeeded};
    assert_eq!(statuses(&m), [Succeeded, Skipped, Succeeded]);

    // The reads changed: everything reruns.
    std::thread::sleep(Duration::from_millis(10));
    fs::write(&input, "other reads").expect("input");
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("reruns");
    assert_eq!(statuses(&m), [Succeeded; 3]);

    // Rerun::All ignores the stamps.
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), ALL, &GO, &mut m, |_| {}).expect("reruns");
    assert_eq!(statuses(&m), [Succeeded; 3]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_failed_qc_stage_fails_nothing_else() {
    let _serial = serial();
    // F9: the QC summary's stage failing kills no other stage and leaves the run complete; its
    // record says why.
    let root = temp_dir("qc-fails");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = [
        stage(StageName::Snv, "sleep 0.5; echo a > $1", &[&input], &[root.join("snv/a.txt")], &[]),
        stage(StageName::Qc, "echo broken >&2; exit 3", &[&input], &[root.join("qc/S.qc.json")], &[]),
    ];
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("completes");
    assert_eq!(m.status, RunStatus::Complete);
    assert_eq!(statuses(&m), [StageStatus::Succeeded, StageStatus::Failed]);
    assert_eq!(m.stages[1].error.as_deref(), Some("exit code 3"));
    // Not done: run again, it runs again.
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("completes");
    assert_eq!(statuses(&m), [StageStatus::Skipped, StageStatus::Failed]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_model_replaced_with_one_of_the_same_size_and_time_reruns_its_stages() {
    let _serial = serial();
    // LRP-001: the stage reads `model` beside the reads; a new model written over the old one,
    // even with the same size and modification time, must rerun the stage and the one after it.
    let root = temp_dir("model");
    let (input, model) = (root.join("reads.bam"), root.join("pileup.pt"));
    fs::write(&input, "reads").expect("input");
    fs::write(&model, "MODEL-A").expect("model");
    let layout = Layout::new(&root, "S");
    let (calls, merged) = (root.join("snv/calls.txt"), root.join("merged.txt"));
    let stages = [
        stage(StageName::Snv, "date +%s%N > $1", &[&input, &model], std::slice::from_ref(&calls), &[]),
        stage(StageName::Merge, "echo merged > $1", &[&calls], std::slice::from_ref(&merged), &[StageName::Snv]),
    ];
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");

    let time = fs::metadata(&model).and_then(|m| m.modified()).expect("model time");
    fs::write(&model, "MODEL-B").expect("model");
    fs::File::options().write(true).open(&model).and_then(|f| f.set_modified(time)).expect("same time");
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("reruns");
    assert_eq!(statuses(&m), [StageStatus::Succeeded; 2]);

    // The same model again: nothing reruns.
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
    assert_eq!(statuses(&m), [StageStatus::Skipped; 2]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_missing_corrupt_or_replaced_index_reruns_its_stage() {
    let _serial = serial();
    // LRP-013: the index is one of the stage's outputs, so resume does not skip a stage whose
    // VCF is there without it.
    let root = temp_dir("index");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let vcf = root.join("S.vcf.gz");
    let outputs = [vcf.clone(), Layout::index(&vcf, "tbi")];
    let index = root.join("expected.tbi");
    write_empty_tabix(&index, pasta::csi::binning_index::index::header::Builder::vcf().build());
    let mut merge = stage(StageName::Merge, "echo v > \"$1\"; cp \"$3\" \"$2\"", &[&input], &outputs, &[]);
    merge.args.push(index.into());
    let stages = [merge];
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
    let resume = || {
        let mut m = manifest(&input);
        execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
        statuses(&m)
    };
    assert_eq!(resume(), [StageStatus::Skipped]);
    fs::remove_file(&outputs[1]).expect("index removed");
    assert_eq!(resume(), [StageStatus::Succeeded]);
    fs::write(&outputs[1], []).expect("index truncated");
    assert_eq!(resume(), [StageStatus::Succeeded]);
    write_empty_tabix(&outputs[1], pasta::csi::binning_index::index::header::Builder::bed().build());
    assert_eq!(resume(), [StageStatus::Succeeded]);
    assert_eq!(resume(), [StageStatus::Skipped]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_damaged_output_of_the_same_size_and_time_reruns_its_stage() {
    let _serial = serial();
    // lrpipe#46: a stamp checked only that outputs were there (and indexes' contents), so a
    // flipped bit in a VCF's BGZF block or a truncated QC JSON was taken as done.
    let root = temp_dir("damaged");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let expected = root.join("expected.vcf.gz");
    let mut writer = pasta::bgzf::io::Writer::new(fs::File::create(&expected).expect("VCF"));
    std::io::Write::write_all(&mut writer, b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        .expect("written");
    writer.finish().expect("finished");
    let vcf = root.join("S.vcf.gz");
    let json = root.join("S.qc.json");
    let mut snv = stage(StageName::Snv, "cp \"$2\" \"$1\"", &[&input], std::slice::from_ref(&vcf), &[]);
    snv.args.push(expected.into());
    let qc = stage(StageName::Qc, "echo '{\"reads\": 1}' > \"$1\"", &[&input], std::slice::from_ref(&json), &[]);
    let stages = [snv, qc];
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
    let resume = || {
        let mut m = manifest(&input);
        execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
        statuses(&m)
    };
    // Changes `path` in place, keeping its size and modification time.
    let damage = |path: &Path, change: &dyn Fn(&mut Vec<u8>)| {
        let time = fs::metadata(path).and_then(|m| m.modified()).expect("a time");
        let mut bytes = fs::read(path).expect("read");
        change(&mut bytes);
        fs::write(path, &bytes).expect("written");
        fs::File::options().write(true).open(path).and_then(|f| f.set_modified(time)).expect("time kept");
    };
    assert_eq!(resume(), [StageStatus::Skipped, StageStatus::Skipped]);
    // One bit of the first BGZF block's CRC32, the 8 bytes before its last 4.
    damage(&vcf, &|bytes| {
        let block = usize::from(u16::from_le_bytes([bytes[16], bytes[17]])) + 1;
        bytes[block - 8] ^= 1;
    });
    assert_eq!(resume(), [StageStatus::Succeeded, StageStatus::Skipped]);
    damage(&json, &|bytes| bytes[0] = b'[');
    assert_eq!(resume(), [StageStatus::Skipped, StageStatus::Succeeded]);
    let full = fs::read(&json).expect("read");
    fs::write(&json, &full[..full.len() / 2]).expect("truncated");
    assert_eq!(resume(), [StageStatus::Skipped, StageStatus::Succeeded]);
    assert_eq!(resume(), [StageStatus::Skipped, StageStatus::Skipped]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_stamp_without_its_outputs_state_is_taken_and_written_anew() {
    let _serial = serial();
    // lrpipe 0.2.1 and before stamped no outputs' state: a directory they wrote still resumes.
    let root = temp_dir("stamp-021");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = [stage(StageName::Qc, "echo '{}' > \"$1\"", &[&input], &[root.join("S.qc.json")], &[])];
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
    let current = fs::read_to_string(layout.stamp(StageName::Qc)).expect("a stamp");
    let (earlier, _) = current.rsplit_once("\noutputs ").expect("the outputs' state");
    fs::write(layout.stamp(StageName::Qc), earlier).expect("written");
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
    assert_eq!(statuses(&m), [StageStatus::Skipped]);
    assert_eq!(fs::read_to_string(layout.stamp(StageName::Qc)).expect("a stamp"), current, "written anew");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn stages_that_would_exceed_the_memory_budget_run_one_after_the_other() {
    let _serial = serial();
    let root = temp_dir("memory");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = diamond(&root, &input);
    let mut m = manifest(&input);
    let mut events = Vec::new();
    let t = Instant::now();
    let tight = Settings { rerun: Rerun::Changed, memory_gib: 1.5 };
    execute(&stages, &layout, &lock(&layout).expect("locked"), tight, &GO, &mut m, |e| events.push(e)).expect("runs");
    assert!(t.elapsed() >= Duration::from_millis(1000), "the two 0.5 s stages ran together: {:?}", t.elapsed());
    assert!(events.contains(&Event::Waiting(StageName::Reads, 1.0, 1.0)), "{events:?}");
    assert_eq!(statuses(&m), [StageStatus::Succeeded; 3]);

    // A stage bigger than the whole budget still runs, alone.
    let mut m = manifest(&input);
    let tiny = Settings { rerun: Rerun::All, memory_gib: 0.5 };
    execute(&stages, &layout, &lock(&layout).expect("locked"), tiny, &GO, &mut m, |_| {}).expect("runs");
    assert_eq!(statuses(&m), [StageStatus::Succeeded; 3]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_failure_stops_the_running_stage_and_starts_no_others() {
    let _serial = serial();
    let root = temp_dir("failure");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let mut stages = diamond(&root, &input);
    stages[0] = stage(StageName::Snv, "echo boom >&2; exit 7", &[&input], &[root.join("snv/a.txt")], &[]);
    stages[1] = stage(StageName::Reads, "sleep 30; echo b > $1", &[&input], &[root.join("reads/b.txt")], &[]);

    let mut m = manifest(&input);
    let t = Instant::now();
    let err =
        execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect_err("fails");
    assert!(t.elapsed() < Duration::from_secs(10), "the sleeping stage was not killed");
    assert!(matches!(&err, RunError::Stage { stage: StageName::Snv, reason, .. } if reason == "exit code 7"), "{err}");
    use StageStatus::{Cancelled, Failed, NotRun};
    assert_eq!(statuses(&m), [Failed, Cancelled, NotRun]);
    assert_eq!(m.status, RunStatus::Failed);
    assert_eq!(fs::read_to_string(layout.log(StageName::Snv)).expect("log"), "boom\n");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_failed_or_stopped_stage_leaves_no_empty_outputs_or_scratch() {
    let _serial = serial();
    // Ctrl-C while bonefish or the haplotagging had opened their outputs left 0-byte BAMs, and
    // the align stage's scratch directory, behind.
    let root = temp_dir("unfinished");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(root.join("out"), "S");
    let (empty, written) = (layout.root().join("align/S.bam"), layout.root().join("align/S.bam.bai"));
    let mut s = stage(
        StageName::Align,
        "touch $1; echo x > $2; mkdir -p $3; touch $3/run; exit 3",
        &[&input],
        &[empty.clone(), written.clone()],
        &[],
    );
    let scratch = layout.root().join("align/tmp");
    s.args.push(scratch.clone().into());
    s.scratch = Some(scratch.clone());
    let mut m = manifest(&input);
    let e = execute(&[s.clone()], &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {});
    assert!(matches!(e, Err(RunError::Stage { .. })), "{e:?}");
    assert!(!empty.exists() && !scratch.exists(), "the empty output and the scratch are removed");
    assert!(written.exists(), "what it wrote is left");
    // And it resumes: the stage runs again, and succeeds.
    let valid = root.join("fixture.bam");
    write_bam(&valid, &[]);
    s.args[1] = "cp \"$4\" \"$1\" && cp \"$5\" \"$2\" && mkdir -p \"$3\"".into();
    s.args.extend([valid.clone().into(), Layout::index(&valid, "bai").into()]);
    let mut m = manifest(&input);
    execute(&[s], &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
    assert_eq!(statuses(&m), [StageStatus::Succeeded]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn success_without_the_outputs_is_a_failure() {
    let _serial = serial();
    let root = temp_dir("no-output");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = vec![stage(StageName::Snv, "true", &[&input], &[root.join("snv/a.txt")], &[])];
    let mut m = manifest(&input);
    let err =
        execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect_err("fails");
    assert!(err.to_string().contains("without writing"), "{err}");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_program_that_does_not_exist_fails_the_run() {
    let _serial = serial();
    let root = temp_dir("no-program");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let mut s = stage(StageName::Snv, "true", &[&input], &[], &[]);
    s.program = root.join("no-such-tool");
    let mut m = manifest(&input);
    let err = execute(&[s], &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect_err("fails");
    assert!(err.to_string().contains("could not start"), "{err}");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn an_interrupted_run_stops_its_stages_and_resumes_later() {
    let _serial = serial();
    let root = temp_dir("interrupt");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let mut stages = diamond(&root, &input);
    stages[1] = stage(StageName::Reads, "sleep 30; echo b > $1", &[&input], &[root.join("reads/b.txt")], &[]);

    let stop = AtomicBool::new(false);
    let mut m = manifest(&input);
    let t = Instant::now();
    let err = std::thread::scope(|s| {
        s.spawn(|| {
            std::thread::sleep(Duration::from_millis(800));
            stop.store(true, Ordering::Relaxed);
        });
        execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &stop, &mut m, |_| {})
            .expect_err("interrupted")
    });
    assert!(t.elapsed() < Duration::from_secs(10), "the sleeping stage was not killed");
    assert!(matches!(err, RunError::Interrupted), "{err}");
    use StageStatus::{Cancelled, NotRun, Succeeded};
    assert_eq!(statuses(&m), [Succeeded, Cancelled, NotRun]);
    assert_eq!(m.status, RunStatus::Interrupted);
    assert_eq!(m.stages[1].error.as_deref(), Some("killed when the run was interrupted"));

    // Run again, the finished stage is skipped and the stopped one runs.
    stages[1] = stage(StageName::Reads, "echo b > $1", &[&input], &[root.join("reads/b.txt")], &[]);
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("resumes");
    use StageStatus::Skipped;
    assert_eq!(statuses(&m), [Skipped, Succeeded, Succeeded]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_second_run_in_the_same_directory_is_refused() {
    let _serial = serial();
    let root = temp_dir("busy");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = vec![stage(StageName::Snv, "sleep 1; echo a > $1", &[&input], &[root.join("snv/a.txt")], &[])];
    std::thread::scope(|s| {
        let first = s.spawn(|| {
            execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut manifest(&input), |_| {})
        });
        std::thread::sleep(Duration::from_millis(300));
        let second = lock(&layout);
        assert!(matches!(second, Err(RunError::Busy(_))), "{second:?}");
        first.join().expect("first run").expect("the first run finishes");
    });
    let _ = fs::remove_dir_all(root);
}

/// Writes the shell's niceness to `$1`: from `/proc` on Linux (its 19th field; `ps` can fail in
/// containers), else with `ps` (`nice` alone prints it only with Linux's coreutils).
const NICENESS: &str =
    "if [ -r /proc/$$/stat ]; then cut -d' ' -f19 /proc/$$/stat; else ps -o nice= -p $$; fi > \"$1\"";

#[test]
fn background_stages_run_at_a_lower_cpu_priority() {
    const CHILD: &str = "LRPIPE_TEST_INHERITED_PRIORITY";
    if std::env::var_os(CHILD).is_none() {
        let inherited = Command::new("nice")
            .args(["-n", "5"])
            .arg(std::env::current_exe().expect("test binary"))
            .args(["--exact", "background_stages_run_at_a_lower_cpu_priority"])
            .env(CHILD, "1")
            .output()
            .expect("runs with increased inherited niceness");
        assert!(inherited.status.success(), "{}", String::from_utf8_lossy(&inherited.stdout));
    }
    let _serial = serial();
    let dir = temp_dir("nice");
    let layout = Layout::new(dir.join("out"), "S");
    let input = dir.join("in.txt");
    fs::write(&input, "x").expect("input");
    let critical = dir.join("critical.txt");
    let background = dir.join("background.txt");
    let mut stages = vec![
        stage(StageName::Snv, NICENESS, &[&input], std::slice::from_ref(&critical), &[]),
        stage(StageName::Reads, NICENESS, &[&input], std::slice::from_ref(&background), &[]),
    ];
    stages[1].priority = Priority::Background;
    let mut manifest = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut manifest, |_| {}).expect("runs");
    let base: i32 = fs::read_to_string(&critical).expect("written").trim().parse().expect("a number");
    let low: i32 = fs::read_to_string(&background).expect("written").trim().parse().expect("a number");
    assert_eq!(low, (base + 10).min(19));
}

/// Writes an earlier run's manifest holding `stages` (name, status, peak GiB).
fn earlier_manifest(layout: &Layout, stages: &[(&str, &str, f64)]) {
    let stages: Vec<serde_json::Value> = stages
        .iter()
        .map(|(name, status, gib)| {
            serde_json::json!({"name": name, "status": status, "command": [], "log": "x.log", "max_rss_gib": gib})
        })
        .collect();
    fs::write(layout.manifest(), serde_json::json!({ "stages": stages }).to_string()).expect("manifest");
}

#[test]
fn a_stage_killed_early_keeps_its_estimated_memory() {
    let _serial = serial();
    // LRUN-02: a stage interrupted or failed soon after it started records a tiny peak; the
    // resumed run must not take that for its peak and start the other stages beside it.
    let root = temp_dir("killed-peak");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let stages = diamond(&root, &input);
    let waited = |events: &[Event]| events.iter().find(|e| matches!(e, Event::Waiting(StageName::Reads, ..))).cloned();
    let run = |budget: f64| {
        let mut events = Vec::new();
        let settings = Settings { rerun: Rerun::All, memory_gib: budget };
        execute(&stages, &layout, &lock(&layout).expect("locked"), settings, &GO, &mut manifest(&input), |e| {
            events.push(e)
        })
        .expect("runs");
        events
    };
    for status in ["cancelled", "failed", "running"] {
        earlier_manifest(&layout, &[("snv", status, 0.001)]);
        let events = run(1.5);
        assert_eq!(waited(&events), Some(Event::Waiting(StageName::Reads, 1.0, 1.0)), "{status}: {events:?}");
    }

    // A peak measured by a run that succeeded, lower than the estimate, does not lower it.
    earlier_manifest(&layout, &[("snv", "succeeded", 0.001), ("reads", "skipped", 0.001)]);
    assert_eq!(waited(&run(1.5)), Some(Event::Waiting(StageName::Reads, 1.0, 1.0)));
    // A higher one raises it, with room for the next run's to be higher still.
    earlier_manifest(&layout, &[("snv", "succeeded", 2.0)]);
    match waited(&run(3.0)) {
        Some(Event::Waiting(StageName::Reads, needs, in_use)) => {
            assert_eq!(needs, 1.0);
            assert!((in_use - 2.0 * lrpipe::run::MEASURED_HEADROOM).abs() < 1e-9, "{in_use}");
        }
        other => panic!("reads did not wait for snv's measured peak: {other:?}"),
    }
    // A stage that failed after passing its estimate needs at least what it reached.
    earlier_manifest(&layout, &[("snv", "failed", 2.0)]);
    assert_eq!(waited(&run(2.5)), Some(Event::Waiting(StageName::Reads, 1.0, 2.0)));
    let _ = fs::remove_dir_all(root);
}

#[test]
fn thread_counts_and_path_spellings_do_not_rerun_stages() {
    let _serial = serial();
    // LRUN-07: a run resumed on a node with other cores, or with its paths spelt otherwise,
    // skips what is done; a change to another option still reruns.
    let root = temp_dir("spelling");
    fs::create_dir_all(root.join("sub")).expect("scratch directory");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let snv = |threads: &str, input: &Path, output: PathBuf, extra: &[&str]| {
        let mut s = stage(StageName::Snv, "date +%s%N > $1", &[input], &[output], &[]);
        s.args.extend(["--threads".into(), threads.into(), "--bam".into(), input.into()]);
        s.args.extend(extra.iter().map(Into::into));
        vec![s]
    };
    let run = |stages: &[Stage]| {
        let mut m = manifest(&input);
        execute(stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
        statuses(&m)
    };
    assert_eq!(run(&snv("4", &input, root.join("snv/a.txt"), &[])), [StageStatus::Succeeded]);
    assert_eq!(run(&snv("8", &input, root.join("snv/a.txt"), &[])), [StageStatus::Skipped], "--threads 8");
    let (other_input, other_output) = (root.join("sub/../reads.bam"), root.join("sub/../snv/./a.txt"));
    assert_eq!(run(&snv("8", &other_input, other_output, &[])), [StageStatus::Skipped], "paths spelt otherwise");
    assert_eq!(run(&snv("8", &input, root.join("snv/a.txt"), &["--min-qual", "5"])), [StageStatus::Succeeded]);
    // bonefish's threads and batch sizes are the align stage's alone to leave out: another
    // stage's -t or -K is an option like any other.
    let out = root.join("snv/a.txt");
    assert_eq!(run(&snv("8", &input, out.clone(), &["-t", "1", "-K", "1"])), [StageStatus::Succeeded]);
    assert_eq!(run(&snv("8", &input, out.clone(), &["-t", "2", "-K", "1"])), [StageStatus::Succeeded], "-t");
    assert_eq!(run(&snv("8", &input, out, &["-t", "2", "-K", "2"])), [StageStatus::Succeeded], "-K");
    let align = |threads: &str, batch: &str, sort: &str| {
        let mut s = stage(StageName::Align, "date +%s%N > $1", &[&input], &[root.join("align/S.bam")], &[]);
        s.args.extend(["-t", threads, "-K", batch, "--sort-mem", sort, "-x", "dna"].map(Into::into));
        vec![s]
    };
    assert_eq!(run(&align("16", "100M", "4G")), [StageStatus::Succeeded]);
    assert_eq!(run(&align("4", "10M", "1G")), [StageStatus::Skipped], "bonefish writes the same");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_moved_or_copied_output_directory_resumes() {
    let _serial = serial();
    // The stamps held the output directory's path: a directory moved, or copied with its files'
    // times, redid every stage, the alignment included. Its paths, in the stages' arguments and
    // inputs (an earlier stage's output), are now stamped as paths under it.
    let base = temp_dir("moved");
    let input = base.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let stages = |root: &Path| {
        let (a, b) = (root.join("snv/a.txt"), root.join("cnv/b.txt"));
        let mut first = stage(StageName::Snv, "date +%s%N > $1", &[&input], std::slice::from_ref(&a), &[]);
        first.args.extend(["--output".into(), root.join("snv").into()]);
        let second = stage(StageName::Cnv, "echo b > $1", &[&a], &[b], &[StageName::Snv]);
        vec![first, second]
    };
    let run = |root: &Path| {
        let layout = Layout::new(root, "S");
        let mut m = manifest(&input);
        execute(&stages(root), &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
        statuses(&m)
    };
    let (first, moved, copied) = (base.join("out"), base.join("moved"), base.join("copied"));
    assert_eq!(run(&first), [StageStatus::Succeeded; 2]);
    assert_eq!(run(&first), [StageStatus::Skipped; 2]);
    fs::rename(&first, &moved).expect("moved");
    assert_eq!(run(&moved), [StageStatus::Skipped; 2], "moved");
    let status = std::process::Command::new("cp").arg("-a").arg(&moved).arg(&copied).status().expect("cp runs");
    assert!(status.success());
    assert_eq!(run(&copied), [StageStatus::Skipped; 2], "copied with its times");
    // An input changed still reruns.
    fs::write(&input, "other reads").expect("input");
    assert_eq!(run(&copied), [StageStatus::Succeeded; 2]);
    let _ = fs::remove_dir_all(base);
}

/// A stand-in for the lrpipe binary's tools, which writes each stage's outputs.
const FAKE_LRPIPE: &str = r#"#!/bin/sh
if [ "$1" = "--version" ]; then echo "fake-lrpipe 0"; exit 0; fi
fixtures=$(dirname "$0")
cmd=$1; shift
if [ "$cmd" = sv ] && [ "$1" = phase ]; then cmd=sv-phase; shift; fi
while [ $# -gt 0 ]; do
  case "$1" in
    --output|--vcf|--snf|--depth|--output-dir|--sample-id|--metadata|--bedmethyl|--repeats-vcf|--output-prefix|--phasing-evidence|--output-haplotag-list|--read-stats|--json|--tsv) eval "opt_$(echo "${1#--}" | tr - _)=\$2"; shift;;
    --gvcf=true) opt_gvcf=1;;
    --combined-only) opt_combined_only=1;;
  esac
  shift
done
w() { mkdir -p "$(dirname "$1")" && echo "${2:-x}" > "$1"; }
i() { cp "$fixtures/$2" "$1"; }
case $cmd in
  snv) # As lrpipe-snv's, its VCFs gain PL with --gvcf.
       for f in merge_output.vcf.gz phased_merge_output.vcf.gz ${opt_gvcf:+merge_output.gvcf.gz}; do
         w "$opt_output/$f" "x${opt_gvcf:+ PL}"; i "$opt_output/$f.tbi" empty.vcf.tbi; done;;
  haplotag) # A new list each time, as a new phasing gives; the BAM with its index. With
       # `slow-haplotag` in the run's directory, it waits there first.
       if [ -e slow-haplotag ]; then echo $$ > haplotag.pid; sleep 30; fi
       w "$opt_output_haplotag_list" "$(date +%s%N)"; if [ -n "$opt_output" ]; then w "$opt_output"; i "$opt_output.bai" empty.bam.bai; fi;;
  sv-phase) w "$opt_vcf"; w "$opt_snf";;
  tr) w "$opt_output"; i "$opt_output.tbi" empty.vcf.tbi;;
  sv) for f in "$opt_vcf" "$opt_snf" "$opt_bedmethyl" "$opt_repeats_vcf" "$opt_phasing_evidence" "$opt_read_stats"; do if [ -n "$f" ]; then w "$f"; fi; done
      if [ -n "$opt_repeats_vcf" ]; then i "$opt_repeats_vcf.tbi" empty.vcf.tbi; fi
      if [ -n "$opt_depth" ]; then
        w "$opt_depth.regions.bed.gz"; i "$opt_depth.regions.bed.gz.csi" empty.bed.csi; w "$opt_depth.mosdepth.summary.txt"
      fi;;
  cnv) w "$opt_output_dir/$opt_sample_id.vcf.gz"; i "$opt_output_dir/$opt_sample_id.vcf.gz.tbi" empty.vcf.tbi
       w "$opt_output_dir/${opt_sample_id}_cnv.bed.gz"; i "$opt_output_dir/${opt_sample_id}_cnv.bed.gz.tbi" empty.bed.tbi
       if [ ! -e "$opt_metadata" ]; then w "$opt_metadata"; fi;;
  merge) w "$opt_output"; i "$opt_output.tbi" empty.vcf.tbi;;
  merge-reads) w "$opt_output"; i "$opt_output.bai" empty.bam.bai;;
  qc) if [ -n "$FAKE_QC_FAIL" ]; then echo "$FAKE_QC_FAIL" >&2; exit 1; fi; w "$opt_json" '{}'; w "$opt_tsv";;
  methyl) tables="combined hap1 hap2"; if [ -n "$opt_combined_only" ]; then tables=combined; fi
          for t in $tables; do w "$opt_output_prefix.$t.bed.gz"; done;;
esac
exit 0
"#;

/// The read group dorado writes, as the ONT fixture's.
const ONT_READ_GROUP: &str = "@RG\tID:run1\tPU:PAY87954\tPM:PC48B260\tPL:ONT\t\
     DS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 modbase_models=dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v2.0.1 \
     runid=70eb\tLB:S\tSM:S";

/// A BAM of the fake reference's chr1 with `read_groups` in its header and no reads, and its
/// index (which `lrpipe run` reads for HP tags).
fn write_bam(path: &Path, read_groups: &[&str]) {
    use std::io::Write;
    let mut text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:26\n".to_owned();
    for rg in read_groups {
        text.push_str(rg);
        text.push('\n');
    }
    let mut bam = pasta::bgzf::io::Writer::new(fs::File::create(path).expect("BAM"));
    bam.write_all(b"BAM\x01").expect("written");
    bam.write_all(&(text.len() as i32).to_le_bytes()).expect("written");
    bam.write_all(text.as_bytes()).expect("written");
    bam.write_all(&1i32.to_le_bytes()).expect("written");
    bam.write_all(&5i32.to_le_bytes()).expect("written");
    bam.write_all(b"chr1\0").expect("written");
    bam.write_all(&26i32.to_le_bytes()).expect("written");
    bam.finish().expect("BAM written");
    match pasta::bam::fs::index(path).expect("indexed") {
        pasta::bam::Index::Bai(index) => {
            pasta::bam::bai::fs::write(format!("{}.bai", path.display()), &index).expect("index")
        }
        pasta::bam::Index::Csi(_) => unreachable!("a 26-base contig is indexed by a BAI"),
    }
}

/// The ONT fixture's reads (lrpipe-snv's), which have no HP tags.
fn fixture_reads() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture/ont.bam")
}

fn write_empty_tabix(path: &Path, header: pasta::csi::binning_index::index::Header) {
    let mut tabix = pasta::tabix::index::Indexer::default();
    tabix.set_header(header);
    pasta::tabix::fs::write(path, &tabix.build()).expect("empty tabix index");
}

/// A site for `lrpipe run` with [`FAKE_LRPIPE`] as its tools, ONT reads (the fixture's, as
/// `reads.bam`) and a reference with a run of Ns: `site.toml`, with `extra` added to it.
fn fake_site(name: &str, extra: &str) -> PathBuf {
    use pasta::csi::{self, binning_index::index::header::Builder};
    use std::os::unix::fs::PermissionsExt;
    let dir = temp_dir(name);
    fs::copy(fixture_reads(), dir.join("reads.bam")).expect("reads");
    fs::copy(Layout::index(&fixture_reads(), "bai"), dir.join("reads.bam.bai")).expect("index");
    write_bam(&dir.join("empty.bam"), &[]);
    write_empty_tabix(&dir.join("empty.vcf.tbi"), Builder::vcf().build());
    write_empty_tabix(&dir.join("empty.bed.tbi"), Builder::bed().build());
    let csi = csi::binning_index::Indexer::new(14, 5).expect("CSI indexer").set_header(Builder::bed().build()).build(0);
    csi::fs::write(dir.join("empty.bed.csi"), &csi).expect("empty CSI index");
    for (file, text) in [
        ("ref.fa", ">chr1\nACGTACGTNNNNNNNNNNACGTACGT\n"),
        // The fixture reads' chr20 too, which the FASTA need not hold: no stage reads it.
        ("ref.fa.fai", "chr1\t26\t6\t26\t27\nchr20\t64444167\t33\t60\t61\n"),
        ("model/pileup.pt", "pileup"),
        ("model/full_alignment.pt", "full alignment"),
        ("fake-lrpipe", FAKE_LRPIPE),
    ] {
        let path = dir.join(file);
        fs::create_dir_all(path.parent().expect("in the directory")).expect("directory");
        fs::write(path, text).expect("written");
    }
    fs::set_permissions(dir.join("fake-lrpipe"), fs::Permissions::from_mode(0o755)).expect("executable");
    let site =
        format!("reference = \"ref.fa\"\n[tools]\nlrpipe = \"fake-lrpipe\"\n[snv]\nont_model = \"model\"\n{extra}");
    fs::write(dir.join("site.toml"), site).expect("site configuration");
    dir
}

/// `lrpipe run` in `cwd` on the fake site's sample, with `args` after the sample's.
fn lrpipe_run(cwd: &Path, args: &[&str]) -> std::process::Command {
    lrpipe_run_on(cwd, "ont", args)
}

/// [`lrpipe_run`] with reads of `platform`.
fn lrpipe_run_on(cwd: &Path, platform: &str, args: &[&str]) -> std::process::Command {
    let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"));
    command.current_dir(cwd).args(["run", "--sample", "S", "--platform", platform, "--sex", "female"]);
    // 64 GiB, unless `args` give the memory.
    if !args.contains(&"--memory") {
        command.args(["--memory", "64"]);
    }
    command.args(args);
    command
}

fn written_manifest(path: &Path) -> serde_json::Value {
    serde_json::from_str(&fs::read_to_string(path).expect("a manifest")).expect("json")
}

#[test]
fn bacterial_run_calls_haploid_variants_on_short_nonhuman_contigs() {
    use pasta::bam;

    let _serial = serial();
    let dir = fake_site("bacterial", "");
    fs::write(dir.join("ref.fa"), ">chromosome\nACGTACGT\n>replicon\nACGT\n").expect("reference");
    fs::write(dir.join("ref.fa.fai"), "chromosome\t8\t12\t8\t9\nreplicon\t4\t31\t4\t5\n").expect("reference index");
    let bam = dir.join("bacterial.bam");
    let header =
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chromosome\tLN:8\n@SQ\tSN:replicon\tLN:4\n".parse().expect("BAM header");
    let mut writer = bam::io::Writer::new(fs::File::create(&bam).expect("BAM"));
    writer.write_header(&header).expect("written header");
    writer.try_finish().expect("finished BAM");
    let bam::Index::Bai(index) = bam::fs::index(&bam).expect("BAM indexed") else {
        panic!("short contigs use a BAI");
    };
    bam::bai::fs::write(Layout::index(&bam, "bai"), &index).expect("written index");
    fs::write(dir.join("targets.bed"), "replicon\t1\t3\n").expect("targets");

    let run = |extra: &[&str]| {
        std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
            .current_dir(&dir)
            .args([
                "run",
                "--bacterial",
                "--config",
                "site.toml",
                "--reads",
                "bacterial.bam",
                "--sample",
                "S",
                "--platform",
                "ont",
                "--output",
                "out",
                "--memory",
                "64",
                "--gvcf",
            ])
            .args(extra)
            .output()
            .expect("lrpipe runs")
    };
    let out = run(&[]);
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "complete");
    assert_eq!(manifest["sample"]["ploidy"], 1);
    assert_eq!(manifest["sample"]["sex"], "not-applicable");
    assert!(manifest["sample"].get("sex_from_index").is_none());
    let stages = manifest["stages"].as_array().expect("stages");
    assert_eq!(
        stages.iter().map(|s| s["name"].as_str().expect("stage name")).collect::<Vec<_>>(),
        ["snv", "merge", "qc"]
    );
    let snv = stages[0]["command"].as_array().expect("SNV command");
    for flag in ["--include-all-ctgs", "--haploid-precise", "--no-phasing-for-fa", "--gvcf=true"] {
        assert!(snv.iter().any(|arg| arg == flag), "{snv:?}");
    }
    assert!(!snv.iter().any(|arg| arg == "--gender" || arg == "--enable-phasing=true"), "{snv:?}");
    assert!(stages[2]["command"].as_array().expect("QC command").iter().any(|arg| arg == "not-applicable"));

    let out = run(&["--analyses", "all", "--regions", "targets.bed", "--dry-run"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(String::from_utf8_lossy(&out.stdout).contains("--bed-fn"));
    let out = run(&["--analyses", "all,sv"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("--bacterial supports only the snv analysis"));
    let out = run(&["--sex", "female"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("cannot be used with"));
    let _ = fs::remove_dir_all(dir);
}

/// Drops the state of its stage's outputs from `stamp`, as lrpipe 0.2.1 wrote it, so the stage
/// stays done with an output changed by hand: what reads that output still reruns. (With the
/// state there, the stage itself reruns and writes the stand-in's output again.)
fn stamp_without_outputs_state(stamp: &Path) {
    let current = fs::read_to_string(stamp).expect("a stamp");
    let (earlier, _) = current.rsplit_once("\noutputs ").expect("the outputs' state");
    fs::write(stamp, earlier).expect("written");
}

fn stage_statuses(manifest: &serde_json::Value) -> Vec<String> {
    manifest["stages"]
        .as_array()
        .expect("stages")
        .iter()
        .map(|s| s["status"].as_str().expect("a status").into())
        .collect()
}

#[test]
fn a_run_whose_stderr_is_gone_still_finishes() {
    let _serial = serial();
    // LRUN-05: `lrpipe run ... 2>&1 | tee run.log` with tee gone, or a closed terminal: the
    // progress lines fail to write, and the run must go on and write its manifest.
    let dir = fake_site("stderr", "");
    let (reader, writer) = std::io::pipe().expect("a pipe");
    drop(reader);
    let mut command = lrpipe_run(&dir, &["--config", "site.toml", "--reads", "reads.bam", "--output", "out"]);
    let status = command.stderr(writer).status().expect("lrpipe runs");
    assert_eq!(status.code(), Some(0), "{status:?}");
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "complete");
    assert!(stage_statuses(&manifest).iter().all(|s| s == "succeeded"), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn only_the_files_of_the_chosen_analyses_are_needed() {
    let _serial = serial();
    // LRUN-12: a repeat catalogue, tandem repeat BED and HiFi model that are not there do not
    // stop an ONT run of SNVs alone, which reads none of them.
    let extra = "hifi_model = \"gone/hifi\"\n[sv]\ntandem_repeats = \"gone/trf.bed\"\n[repeats]\ncatalog = \"gone/catalog.bed\"\n";
    let dir = fake_site("unused-files", extra);
    let out =
        lrpipe_run(&dir, &["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", "snv"])
            .output()
            .expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let resources = manifest["resources"].as_object().expect("resources");
    assert!(resources.contains_key("ont_model"), "{resources:?}");
    for unused in ["repeat_catalog", "tandem_repeats", "hifi_model"] {
        assert!(!resources.contains_key(unused), "{unused}: {resources:?}");
    }

    // A file the chosen analyses read is still checked.
    let out =
        lrpipe_run(&dir, &["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", "tr"])
            .output()
            .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    let e = String::from_utf8_lossy(&out.stderr);
    assert!(e.contains("the site configuration's [repeats] catalog ") && e.contains("does not exist"), "{e}");

    // So is it in a dry run of aligned reads, which checked only the plan, and the others the
    // configuration names by their keys.
    for (analyses, key) in [("tr", "[repeats] catalog"), ("sv", "[sv] tandem_repeats")] {
        let out = lrpipe_run(
            &dir,
            &["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", analyses, "--dry-run"],
        )
        .output()
        .expect("lrpipe runs");
        assert_eq!(out.status.code(), Some(2), "{analyses}");
        assert!(stderr(&out).contains(&format!("the site configuration's {key} ")), "{}", stderr(&out));
    }
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn memory_below_a_stages_peak_is_warned_of() {
    let _serial = serial();
    // --memory 4 ran each stage alone without a word, though lrpipe-snv expects 14.5 GiB.
    let dir = fake_site("low-memory", "");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--dry-run"];
    let out =
        lrpipe_run(&dir, &[&args[..], &["--memory", "8", "--analyses", "snv,sv"]].concat()).output().expect("runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains(
            "--memory 8.0 GiB is less than a stage is expected to use at its peak (snv 14.5 GiB): it runs with no other \
             stage beside it, and may run out of memory"
        ),
        "{}",
        stderr(&out)
    );
    let out = lrpipe_run(&dir, &args).output().expect("runs");
    assert!(!stderr(&out).contains("may run out of memory"), "64 GiB: {}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_dry_run_refuses_what_the_run_refuses_and_marks_what_is_done() {
    let _serial = serial();
    let dir = fake_site("dry-checks", "");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let dry = |extra: &[&str]| lrpipe_run(&dir, &[&args[..], extra, &["--dry-run"]].concat()).output().expect("runs");
    let plan = |out: &std::process::Output| String::from_utf8_lossy(&out.stdout).into_owned();
    let out = dry(&[]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!plan(&out).contains("done"), "nothing is done yet: {}", plan(&out));
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));

    // Every stage is done; a new gVCF reruns lrpipe-snv and every stage after it.
    let out = dry(&[]);
    let stages: Vec<String> =
        plan(&out).lines().filter(|l| l.starts_with("# ") && l.contains("(after:")).map(str::to_owned).collect();
    assert!(stages.len() > 3 && stages.iter().all(|l| l.ends_with(": done, would be skipped")), "{stages:?}");
    let out = dry(&["--gvcf"]);
    let text = plan(&out);
    assert!(text.contains("# snv (after: -)\n") && text.contains("# cnv (after: snv, reads, sv-phase)\n"), "{text}");
    assert!(text.contains("# reads (after: -): done, would be skipped"), "{text}");
    let out = dry(&["--rerun", "all"]);
    assert!(!plan(&out).contains("would be skipped"), "{}", plan(&out));

    // Another sample into this run's directory: refused, as the run refuses it.
    let other: Vec<&str> = ["run", "--sample", "T", "--platform", "ont", "--sex", "female", "--memory", "64"]
        .into_iter()
        .chain(args)
        .chain(["--dry-run"])
        .collect();
    let out =
        std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe")).current_dir(&dir).args(&other).output().expect("runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("holds the run of sample S, not T"), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_configured_cnv_metadata_is_made_if_absent() {
    let _serial = serial();
    // LRUN-12: as the README says, a configured `[cnv] metadata` that does not exist yet is
    // made from the reference, once: a resumed run skips every stage.
    let dir = fake_site("mdr", "[cnv]\nmetadata = \"shared/site.mdr\"\n");
    fs::create_dir_all(dir.join("shared")).expect("directory");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", "cnv"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let mdr = fs::read_to_string(dir.join("shared/site.mdr")).expect("the .mdr is made");
    assert!(mdr.contains("NSPOS\tchr1\t8\t17\n"), "{mdr}");
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert!(manifest["resources"]["cnv_metadata"]["sha256"].is_string(), "{manifest}");
    let cnv = manifest["stages"].as_array().expect("stages").iter().find(|s| s["name"] == "cnv").expect("cnv");
    assert!(cnv["command"].to_string().contains("shared/site.mdr"), "{cnv}");

    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let statuses = stage_statuses(&manifest);
    assert!(!statuses.is_empty() && statuses.iter().all(|s| s == "skipped"), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_run_resumed_with_other_threads_or_paths_skips_what_is_done() {
    let _serial = serial();
    // LRUN-07: the same run resumed with another --threads, the configuration named by another
    // path, and from another directory.
    let dir = fake_site("resume-spelling", "");
    fs::create_dir_all(dir.join("sub")).expect("directory");
    let (config, reads, output) = (dir.join("site.toml"), dir.join("reads.bam"), dir.join("out"));
    let runs: [(&Path, Vec<&str>); 3] = [
        (&dir, vec!["--config", "./site.toml", "--reads", "reads.bam", "--output", "out", "--threads", "2"]),
        (
            &dir.join("sub"),
            vec!["--config", "../site.toml", "--reads", "../reads.bam", "--output", "../out", "--threads", "1"],
        ),
        (
            Path::new("/"),
            vec![
                "--config",
                config.to_str().expect("UTF-8"),
                "--reads",
                reads.to_str().expect("UTF-8"),
                "--output",
                output.to_str().expect("UTF-8"),
            ],
        ),
    ];
    for (i, (cwd, args)) in runs.iter().enumerate() {
        let out = lrpipe_run(cwd, args).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        let statuses = stage_statuses(&written_manifest(&output.join("manifest.json")));
        let expect = if i == 0 { "succeeded" } else { "skipped" };
        assert!(statuses.iter().all(|s| s == expect), "run {i}: {statuses:?}");
    }
    let _ = fs::remove_dir_all(dir);
}

/// Clair3's models `names` under `models/` in `dir`, as `[snv] models` holds them.
fn install_models(dir: &Path, names: &[&str]) {
    for name in names {
        let model = dir.join("models").join(name);
        fs::create_dir_all(&model).expect("model directory");
        fs::write(model.join("pileup.pt"), format!("{name} pileup")).expect("model");
        fs::write(model.join("full_alignment.pt"), format!("{name} full alignment")).expect("model");
    }
}

/// A site whose `[snv]` is `snv`, in place of [`fake_site`]'s.
fn site_with_snv(dir: &Path, snv: &str) {
    let site = format!("reference = \"ref.fa\"\n[tools]\nlrpipe = \"fake-lrpipe\"\n[snv]\n{snv}");
    fs::write(dir.join("site.toml"), site).expect("site configuration");
}

fn stderr(out: &std::process::Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

#[test]
fn the_fixtures_model_is_chosen_from_its_header_and_recorded() {
    let _serial = serial();
    // F3: the ONT fixture's reads, with no model given, get Clair3's model for their basecall
    // model from `[snv] models`; the manifest says which and why, and the choice is the snv
    // stage's, so another model reruns it.
    let dir = fake_site("model-choice", "");
    install_models(&dir, &["r1041_e82_400bps_sup_v500", "r1041_e82_400bps_hac_v500", "hifi_revio"]);
    site_with_snv(&dir, "models = \"models\"\n");
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture/ont.bam");
    let fixture = fixture.to_str().expect("UTF-8");
    let args = ["--config", "site.toml", "--reads", fixture, "--output", "out", "--analyses", "snv"];
    let chosen = fs::canonicalize(&dir).expect("canonical").join("models/r1041_e82_400bps_sup_v500");

    // The BAM and the CRAM alike.
    let cram = fixture.replace("ont.bam", "ont.cram");
    for reads in [fixture, &cram] {
        let dry = ["--config", "site.toml", "--reads", reads, "--output", "out", "--analyses", "snv", "--dry-run"];
        let out = lrpipe_run(&dir, &dry).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", stderr(&out));
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            stdout.starts_with(&format!("# snv model: {} (snv.models: chosen from the header", chosen.display())),
            "{reads}: {stdout}"
        );
        assert!(stdout.contains("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), "{stdout}");
        // None but that the fixture's read group names its sample hg003, not S.
        let warnings = stderr(&out).lines().filter(|l| l.contains("warning") && !l.contains("SM:hg003")).count();
        assert_eq!(warnings, 0, "{}", stderr(&out));
    }

    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let model = &manifest["snv_model"];
    assert_eq!(model["path"], chosen.to_str().expect("UTF-8"), "{manifest}");
    assert_eq!(model["by"], "snv.models");
    assert!(model["reason"].as_str().expect("a reason").contains("sup@v5.0.0"), "{model}");
    assert!(model.get("warnings").is_none(), "{model}");
    assert!(manifest["resources"]["ont_model"]["sha256"].is_string(), "{manifest}");
    let snv = &manifest["stages"][0];
    assert!(snv["command"].to_string().contains(&format!("\"{}\"", chosen.display())), "{snv}");

    // Resumed: nothing to do. Another model given: snv reruns, with a warning naming the one
    // for these reads.
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert_eq!(stage_statuses(&written_manifest(&dir.join("out/manifest.json"))), ["skipped", "skipped", "skipped"]);
    let out_hac = lrpipe_run(&dir, &[&args[..], &["--model", "r1041_e82_400bps_hac_v500"]].concat())
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success() && out_hac.status.success(), "{}", stderr(&out_hac));
    assert!(
        stderr(&out_hac).contains(
            "warning: r1041_e82_400bps_hac_v500 is Clair3's model for Oxford Nanopore R10.4.1 E8.2 400 bps 5 kHz \
             reads basecalled with dorado v5.0.0 HAC, but these are Oxford Nanopore reads basecalled with \
             dna_r10.4.1_e8.2_400bps_sup@v5.0.0 (Clair3's model for them: r1041_e82_400bps_sup_v500)"
        ),
        "{}",
        stderr(&out_hac)
    );
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(stage_statuses(&manifest)[0], "succeeded", "{manifest}");
    assert_eq!(manifest["snv_model"]["by"], "--model");
    assert_eq!(manifest["snv_model"]["warnings"].as_array().map(Vec::len), Some(1), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_hifi_model_for_ont_reads_is_refused() {
    let _serial = serial();
    // F3: a HiFi model given for ONT reads, on the command line or in the site configuration,
    // stops the run before it writes anything; --model-mismatch warn runs it, warned.
    let dir = fake_site("model-mismatch", "");
    install_models(&dir, &["hifi_revio", "r1041_e82_400bps_sup_v500"]);
    site_with_snv(&dir, "models = \"models\"\nhifi_model = \"models/hifi_revio\"\n");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", "snv"];
    let out = lrpipe_run(&dir, &[&args[..], &["--model", "hifi_revio"]].concat()).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(
        stderr(&out).contains(
            "lrpipe: error: hifi_revio is Clair3's model for PacBio HiFi reads from REVIO, not ont reads; give the \
             model and --platform the reads were made for, or --model-mismatch warn to run anyway"
        ),
        "{}",
        stderr(&out)
    );
    assert!(!dir.join("out").exists(), "refused before writing anything");

    // The site's HiFi model, for ONT reads said to be HiFi.
    let out = lrpipe_run_on(&dir, "hifi", &args).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    for says in [
        "the header says these are Oxford Nanopore reads basecalled with dna_r10.4.1_e8.2_400bps_sup@v5.0.0, not hifi \
         reads as --platform says",
        "hifi_revio is Clair3's model for PacBio HiFi reads from REVIO, but these are Oxford Nanopore reads",
    ] {
        assert!(stderr(&out).contains(says), "{says}: {}", stderr(&out));
    }

    let out = lrpipe_run(&dir, &[&args[..], &["--model", "hifi_revio", "--model-mismatch", "warn"]].concat())
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("warning: hifi_revio is Clair3's model for PacBio HiFi"), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert!(manifest["snv_model"]["path"].as_str().expect("a path").ends_with("models/hifi_revio"), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn pacbio_headers_and_several_read_groups_choose_or_refuse() {
    let _serial = serial();
    let dir = fake_site("model-headers", "");
    install_models(&dir, &["hifi_revio", "hifi_sequel2", "r1041_e82_400bps_sup_v500", "r1041_e82_400bps_sup_v520"]);
    site_with_snv(&dir, "models = \"models\"\n");
    let revio =
        "@RG\tID:m84\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=102-739-100;BASECALLERVERSION=5.0\tPU:m84\tPM:REVIO\tSM:S";
    let sequel = "@RG\tID:m64\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=101-894-200\tPU:m64\tPM:SEQUELII\tSM:S";
    let dry = |reads: &str, platform: &str, more: &[&str]| {
        let mut command =
            lrpipe_run_on(&dir, platform, &["--config", "site.toml", "--reads", reads, "--output", "out"]);
        command.args(["--analyses", "snv", "--dry-run"]).args(more);
        command.output().expect("lrpipe runs")
    };
    let model_line =
        |out: &std::process::Output| String::from_utf8_lossy(&out.stdout).lines().next().unwrap_or_default().to_owned();
    for (rg, model) in [(revio, "models/hifi_revio "), (sequel, "models/hifi_sequel2 ")] {
        write_bam(&dir.join("pacbio.bam"), &[rg]);
        let out = dry("pacbio.bam", "hifi", &[]);
        assert!(out.status.success(), "{}", stderr(&out));
        assert!(model_line(&out).contains(model), "{}", model_line(&out));
    }

    // Two flow cells basecalled with different models: no one model fits, unless given.
    let newer = ONT_READ_GROUP.replace("ID:run1", "ID:run2").replace("sup@v5.0.0", "sup@v5.2.0");
    write_bam(&dir.join("two.bam"), &[ONT_READ_GROUP, &newer]);
    let out = dry("two.bam", "ont", &[]);
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("so no one model fits them all: give the model with --model"), "{}", stderr(&out));
    let out = dry("two.bam", "ont", &["--model", "r1041_e82_400bps_sup_v520"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(model_line(&out).contains("models/r1041_e82_400bps_sup_v520 (--model: "), "{}", model_line(&out));
    assert!(stderr(&out).contains("(Clair3's model for them: r1041_e82_400bps_sup_v500)"), "{}", stderr(&out));

    // A model the directory does not hold: the nearest it holds, warned; none of the chemistry,
    // an error naming the one to add.
    write_bam(&dir.join("hac.bam"), &[&ONT_READ_GROUP.replace("sup@", "hac@")]);
    let out = dry("hac.bam", "ont", &[]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(model_line(&out).contains("models/r1041_e82_400bps_sup_v500 "), "{}", model_line(&out));
    assert!(stderr(&out).contains("does not hold r1041_e82_400bps_hac_v500"), "{}", stderr(&out));
    write_bam(&dir.join("slow.bam"), &[&ONT_READ_GROUP.replace("400bps_sup@v5.0.0", "260bps_sup@v4.1.0")]);
    let out = dry("slow.bam", "ont", &[]);
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("add r1041_e82_260bps_sup_v410 (one of ONT's Rerio models"), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_model_given_for_a_run_without_snvs_is_warned_of() {
    let _serial = serial();
    let dir = fake_site("model-unused", "");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", "sv"];
    let out = lrpipe_run(&dir, &[&args[..], &["--model", "/m/r1041_e82_400bps_sup_v500"]].concat())
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("warning: --model /m/r1041_e82_400bps_sup_v500 is not used"), "{}", stderr(&out));
    assert!(written_manifest(&dir.join("out/manifest.json")).get("snv_model").is_none());
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_move_table_model_is_refused_before_any_stage_runs() {
    let _serial = serial();
    // lrpipe snv refuses a _with_mv model without --enable_dwell_time, which lrpipe run does
    // not pass: found when planning, not once the other stages have started.
    let dir = fake_site("model-move-table", "");
    install_models(&dir, &["r1041_e82_400bps_sup_v520_with_mv"]);
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let out = lrpipe_run(&dir, &[&args[..], &["--model", "models/r1041_e82_400bps_sup_v520_with_mv"]].concat())
        .output()
        .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("is a move-table (signal-aware, _with_mv) model"), "{}", stderr(&out));
    assert!(!dir.join("out").exists(), "refused before writing anything");
}

#[test]
fn a_gvcf_reruns_the_snv_stage_and_what_reads_its_calls_and_the_snf_is_always_written() {
    let _serial = serial();
    // F4: `--gvcf` is lrpipe-snv's option, so turning it on or off reruns that stage, and the
    // stages reading its calls, which it writes again with PL (as the stand-in does too): the
    // CNVs and the merge, but not the read pass. The read pass writes the SV candidates for
    // multi-sample calling with the SVs.
    let dir = fake_site("cohort", "");
    let run = |gvcf: bool| {
        let mut args = vec!["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
        if gvcf {
            args.push("--gvcf");
        }
        let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        written_manifest(&dir.join("out/manifest.json"))
    };
    let status = |manifest: &serde_json::Value, stage: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == stage).expect("a stage")["status"].as_str().expect("a status").to_owned()
    };
    let gvcf = dir.join("out/snv/merge_output.gvcf.gz");
    let manifest = run(false);
    assert!(stage_statuses(&manifest).iter().all(|s| s == "succeeded"), "{manifest}");
    assert_eq!(manifest["outputs"]["snf"], dir.join("out/sv/S.snf").to_str().expect("UTF-8"));
    assert!(dir.join("out/sv/S.snf").exists());
    assert!(manifest["outputs"].get("gvcf").is_none() && !gvcf.exists(), "{manifest}");

    let manifest = run(true);
    for stage in ["snv", "cnv", "merge"] {
        assert_eq!(status(&manifest, stage), "succeeded", "{stage}: {manifest}");
    }
    assert_eq!(status(&manifest, "reads"), "skipped", "{manifest}");
    assert_eq!(manifest["outputs"]["gvcf"], gvcf.to_str().expect("UTF-8"));
    assert!(gvcf.exists() && Layout::index(&gvcf, "tbi").exists());
    let snv = &manifest["stages"][0];
    assert!(snv["command"].to_string().contains("--gvcf=true"), "{snv}");

    // Resumed as it was, nothing reruns; without --gvcf, lrpipe-snv does again.
    assert!(stage_statuses(&run(true)).iter().all(|s| s == "skipped"));
    let manifest = run(false);
    assert_eq!(status(&manifest, "snv"), "succeeded", "{manifest}");
    assert_eq!(status(&manifest, "reads"), "skipped", "{manifest}");
    assert!(manifest["outputs"].get("gvcf").is_none(), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn regions_are_checked_sorted_and_rerun_the_stages_that_read_them() {
    let _serial = serial();
    // F7: `--regions` is checked against the reference, sorted in its order and merged into the
    // BED every restricted stage reads, which is part of their fingerprints: other regions rerun
    // them, the same regions (a BED touched, or written in another order) do not. The CNVs are
    // left out, saying why, and an error when asked for.
    let dir = fake_site("regions", "");
    fs::write(dir.join("targets.bed"), "chr20\t30190000\t30200000\nchr1\t0\t10\nchr20\t30180000\t30195000\n")
        .expect("BED");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--regions", "targets.bed"];
    let run = |args: &[&str]| {
        let out = lrpipe_run(&dir, args).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", stderr(&out));
        (written_manifest(&dir.join("out/manifest.json")), stderr(&out))
    };
    let status = |manifest: &serde_json::Value, stage: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == stage).expect("a stage")["status"].as_str().expect("a status").to_owned()
    };
    let (manifest, said) = run(&args);
    assert!(said.contains("cnv is not called in --regions"), "{said}");
    assert!(said.contains("--regions: 2 regions on 2 contigs, 20010 bases"), "{said}");
    // chr1 has no reads: accepted, as Clair3 takes it, with a warning.
    assert!(said.contains("warning: --regions on chr1 get no calls"), "{said}");
    let names: Vec<&str> =
        manifest["stages"].as_array().expect("stages").iter().map(|s| s["name"].as_str().expect("a name")).collect();
    assert!(!names.contains(&"cnv") && names.contains(&"snv") && names.contains(&"reads"), "{names:?}");
    let bed = dir.join("out/.lrpipe/regions.bed");
    assert_eq!(fs::read_to_string(&bed).expect("written"), "chr1\t0\t10\nchr20\t30180000\t30200000\n");
    let windows = fs::read_to_string(dir.join("out/.lrpipe/regions.windows.bed")).expect("written");
    assert_eq!(windows.lines().count(), 21, "{windows}");
    assert!(windows.starts_with("chr1\t0\t10\nchr20\t30180000\t30181000\n"), "{windows}");
    let command = |stage: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == stage).expect("a stage")["command"].to_string()
    };
    assert!(command("snv").contains(&format!("--bed-fn\",\"{}", bed.display())), "{}", command("snv"));
    for wanted in ["--regions", "--bedmethyl-regions"] {
        assert!(command("reads").contains(&format!("{wanted}\",\"{}", bed.display())), "{}", command("reads"));
    }
    assert!(command("haplotag").contains("--regions-bed"), "{}", command("haplotag"));
    let regions = &manifest["resources"]["regions"];
    assert!(regions["path"].as_str().is_some_and(|p| p.ends_with("targets.bed")), "{regions}");

    // The same regions, in another order and touched: nothing reruns.
    fs::write(dir.join("targets.bed"), "chr1\t0\t10\nchr20\t30180000\t30200000\n").expect("BED");
    let (manifest, _) = run(&args);
    assert!(stage_statuses(&manifest).iter().all(|s| s == "skipped"), "{manifest}");
    // Other regions rerun what reads them.
    fs::write(dir.join("targets.bed"), "chr20\t30180000\t30190000\n").expect("BED");
    let (manifest, _) = run(&args);
    for stage in ["snv", "reads", "haplotag"] {
        assert_eq!(status(&manifest, stage), "succeeded", "{stage}: {manifest}");
    }
    assert_eq!(fs::read_to_string(&bed).expect("written"), "chr20\t30180000\t30190000\n");

    // The CNVs asked for, and regions the reference or the reads do not have, stop the run
    // before it starts.
    let fai = fs::read_to_string(dir.join("ref.fa.fai")).expect("the index");
    // The fake tools never read these additional contigs. Their offsets lie beyond the stub FASTA.
    fs::write(dir.join("ref.fa.fai"), format!("{fai}chrNew\t2000000\t100\t60\t61\nchrM\t16569\t110\t60\t61\n"))
        .expect("the index");
    fs::write(dir.join("other.bed"), "chr21\t0\t100\n").expect("BED");
    fs::write(dir.join("unread.bed"), "chrNew\t0\t100\n").expect("BED");
    fs::write(dir.join("mito.bed"), "chrM\t0\t100\n").expect("BED");
    for (bed, analyses, wanted) in [
        ("targets.bed", "snv,cnv", "cnv cannot be called in --regions"),
        ("other.bed", "snv", "contig chr21 is not in the reference"),
        ("unread.bed", "snv", "has no contig chrNew (its header does not name it)"),
        ("mito.bed", "snv", "--regions on chrM are left out"),
        ("gone.bed", "snv", "gone.bed"),
    ] {
        let mut command = lrpipe_run(
            &dir,
            &[
                "--config",
                "site.toml",
                "--reads",
                "reads.bam",
                "--output",
                "refused",
                "--regions",
                bed,
                "--analyses",
                analyses,
            ],
        );
        let out = command.output().expect("lrpipe runs");
        assert!(!out.status.success() && stderr(&out).contains(wanted), "{bed}: {}", stderr(&out));
        assert!(!dir.join("refused/manifest.json").exists(), "{bed}");
    }
    // Every contig a run of the whole genome calls, whole, and the mitochondria, which it does
    // not call: the whole genome's run, with the CNVs, and no list of contigs.
    fs::write(dir.join("all.bed"), "chrM\t0\t16569\nchr1\t0\t26\nchr20\t0\t64444167\nchrNew\t0\t2000000\n")
        .expect("BED");
    let mut command = lrpipe_run(
        &dir,
        &["--config", "site.toml", "--reads", "reads.bam", "--output", "all", "--regions", "all.bed", "--dry-run"],
    );
    let out = command.output().expect("lrpipe runs");
    let (printed, said) = (String::from_utf8_lossy(&out.stdout), stderr(&out));
    assert!(out.status.success() && said.contains("the run is of the whole genome"), "{said}");
    assert!(printed.contains(" cnv ") && !printed.contains("--ctg-name") && !printed.contains("--contig"), "{printed}");
    let _ = fs::remove_dir_all(dir);
}

/// `reads.bam` in `dir`, indexed: reads on chr1 with base modification tags (`MM`/`ML`), the
/// second of them haplotagged (`HP:i:1`, or with `text` `HP:Z:1`).
fn haplotagged_bam(dir: &Path, text: bool) {
    reads_bam(dir, text, true);
}

/// [`haplotagged_bam`]'s reads, with base modification tags only with `modifications`.
fn reads_bam(dir: &Path, text: bool, modifications: bool) {
    use pasta::bam;
    use pasta::core::Position;
    use pasta::sam::alignment::io::Write as _;
    use pasta::sam::alignment::record::cigar::{Op, op::Kind};
    use pasta::sam::alignment::record::{Flags, MappingQuality, data::field::Tag};
    use pasta::sam::alignment::record_buf::{Cigar, RecordBuf, Sequence, data::field::Value};
    use pasta::sam::header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    };

    let length = std::num::NonZeroUsize::new(26).expect("non-empty");
    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build().expect("@HD");
    let header = pasta::sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(length))
        .build();
    let path = dir.join("reads.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("BAM"));
    writer.write_header(&header).expect("header");
    for (i, haplotype) in [None, Some(1u8), None].into_iter().enumerate() {
        let mut record = RecordBuf::builder()
            .set_name(format!("r{i}"))
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(i + 1).expect("non-zero"))
            .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect::<Cigar>())
            .set_sequence(Sequence::from(b"ACGT".to_vec()));
        let mut data = Vec::new();
        if modifications {
            // 5mC at the read's one C.
            data.push((Tag::new(b'M', b'M'), Value::String("C+m?,0;".into())));
            data.push((Tag::new(b'M', b'L'), Value::from(vec![200u8])));
        }
        if let Some(h) = haplotype {
            let value = if text { Value::String(h.to_string().into()) } else { Value::from(h) };
            data.push((Tag::new(b'H', b'P'), value));
        }
        record = record.set_data(data.into_iter().collect());
        writer.write_alignment_record(&header, &record.build()).expect("record");
    }
    writer.try_finish().expect("finished");
    drop(writer);
    let bam::Index::Bai(index) = bam::fs::index(&path).expect("indexed") else { panic!("a BAI") };
    fs::remove_file(dir.join("reads.bam.bai")).expect("the fixture's index removed");
    bam::bai::fs::write(dir.join("reads.bam.bai"), &index).expect("index written");
}

/// Every analysis of a run of the fake site but the phasing (which haplotags the reads).
const UNPHASED: [&str; 2] = ["--analyses", "snv,sv,cnv,depth,methyl"];

#[test]
fn reads_without_hp_tags_get_no_haplotype_tables_and_a_warning() {
    let _serial = serial();
    // F1a: the ONT fixture has no HP tags. Without phasing, which haplotags the reads, and with
    // pb-CpG-tools' tables asked for, only the combined one is planned and written (and an
    // earlier run's per-haplotype tables removed), the run warns at its start and in its
    // summary, and the manifest says so.
    let dir = fake_site("untagged", "[methylation]\nreport = \"cpg\"\n");
    let methylation = dir.join("out/methylation");
    fs::create_dir_all(&methylation).expect("directory");
    // An earlier run's directory, as lrpipe marks it.
    fs::create_dir_all(dir.join("out/.lrpipe")).expect("directory");
    fs::write(dir.join("out/.lrpipe/sample"), "S\n").expect("the earlier run's sample");
    for stale in ["S.cpg.hap1.bed.gz", "S.cpg.hap2.bed.gz"] {
        fs::write(methylation.join(stale), "an earlier run's").expect("written");
    }
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", UNPHASED[0], UNPHASED[1]];
    let dry = lrpipe_run(&dir, &[&args[..], &["--dry-run"]].concat()).output().expect("lrpipe runs");
    assert!(dry.status.success(), "{}", String::from_utf8_lossy(&dry.stderr));
    let plan = String::from_utf8_lossy(&dry.stdout);
    assert!(plan.contains("--combined-only"), "{plan}");

    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let stderr = String::from_utf8_lossy(&out.stderr);
    let warning = "is not haplotagged (no HP tag in any of 152 mapped reads sampled), and `lrpipe run` haplotags reads only \
                   when it phases the SNVs";
    assert_eq!(stderr.matches(warning).count(), 2, "at the start and in the summary: {stderr}");
    assert!(stderr.contains("the SVs are unphased, and no CpG tables per haplotype are written"), "{stderr}");
    assert!(stderr.contains("for them, add phase to --analyses"), "{stderr}");
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["sample"]["haplotypes"], "none (input not haplotagged)");
    // One shape, one file or several: a list.
    let reads = manifest["sample"]["reads"].as_array().expect("a list of the reads");
    assert_eq!(reads.len(), 1);
    assert!(reads[0]["path"].as_str().expect("a path").ends_with("/reads.bam"), "{reads:?}");
    let recorded = manifest["warnings"].as_array().expect("warnings");
    assert!(recorded.iter().any(|w| w.as_str().is_some_and(|w| w.contains(warning))), "{manifest}");
    let outputs = manifest["outputs"].as_object().expect("outputs");
    assert!(outputs.contains_key("cpg_combined"), "{outputs:?}");
    assert!(!outputs.contains_key("cpg_hap1") && !outputs.contains_key("cpg_hap2"), "{outputs:?}");
    assert!(methylation.join("S.cpg.combined.bed.gz").exists());
    assert!(!methylation.join("S.cpg.hap1.bed.gz").exists() && !methylation.join("S.cpg.hap2.bed.gz").exists());

    // Resumed, nothing reruns.
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(stage_statuses(&written_manifest(&dir.join("out/manifest.json"))).iter().all(|s| s == "skipped"));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn reads_haplotagged_elsewhere_get_the_runs_own_haplotypes_when_it_phases() {
    let _serial = serial();
    // F1b: reads haplotagged before the run get the per-haplotype tables and phased SVs from the
    // run's own phasing, which wins over their tags, with no warning (F1a warned that the
    // input's tags were used); without phasing, from their tags, as F1a did.
    let dir = fake_site("tagged", "[methylation]\nreport = \"cpg\"\n");
    haplotagged_bam(&dir, false);
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!stderr.contains("haplotagged"), "{stderr}");
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["sample"]["haplotypes"], "lrpipe run's phasing of the SNVs (its haplotag list)");
    assert!(manifest["warnings"].as_array().is_none_or(Vec::is_empty), "{manifest}");
    let command = |manifest: &serde_json::Value, name: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == name).expect("a stage")["command"].to_string()
    };
    assert!(command(&manifest, "methylation").contains("out/phase/S.haplotags.tsv.gz"), "{manifest}");
    for table in ["combined", "hap1", "hap2"] {
        assert!(dir.join(format!("out/methylation/S.cpg.{table}.bed.gz")).exists(), "{table}");
    }

    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", "snv,sv,methyl"]].concat()).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    // (The site's model is not one of Clair3's by name, which is warned of too.)
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!stderr.contains("haplotagged"), "{stderr}");
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["sample"]["haplotypes"], "the input's HP tags (not lrpipe run's phasing)");
    let methylation = command(&manifest, "methylation");
    assert!(!methylation.contains("--haplotags") && !methylation.contains("--combined-only"), "{methylation}");
    for table in ["combined", "hap1", "hap2"] {
        assert!(manifest["outputs"][format!("cpg_{table}")].is_string(), "{table}: {manifest}");
    }
    assert!(manifest["stages"].as_array().expect("stages").iter().all(|s| s["name"] != "haplotag"));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn the_repeats_follow_the_runs_haplotags_when_it_phases() {
    let _serial = serial();
    // With phasing, the tandem repeats are genotyped with the haplotag list after it, not in
    // the read pass from the reads' own tags; without, in the read pass, and reads without HP
    // tags are warned of.
    let dir = fake_site("repeats", "[repeats]\ncatalog = \"cat.bed\"\n");
    fs::write(dir.join("cat.bed"), "chr1\t1\t5\tID=a;MOTIFS=A\n").expect("written");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let stage = |manifest: &serde_json::Value, name: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == name).cloned()
    };
    let repeats = stage(&manifest, "repeats").expect("a repeats stage");
    assert!(repeats["command"].to_string().contains("out/phase/S.haplotags.tsv.gz"), "{repeats}");
    assert!(!stage(&manifest, "reads").expect("a read pass")["command"].to_string().contains("--repeats"));
    assert!(stage(&manifest, "merge").expect("a merge")["command"].to_string().contains("--tr"));
    assert!(dir.join("out/repeats/S.repeats.vcf.gz").exists());

    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", "snv,tr"]].concat()).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert!(stage(&manifest, "repeats").is_none() && stage(&manifest, "haplotag").is_none());
    assert!(stage(&manifest, "reads").expect("a read pass")["command"].to_string().contains("--repeats"));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("the tandem repeats' alleles are told apart by length alone"), "{stderr}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn phasing_haplotags_the_reads_for_the_svs_and_cpg_tables_and_resumes() {
    let _serial = serial();
    // F1b's stage graph: the read pass beside lrpipe-snv calls the SVs unphased; the haplotag
    // list follows the phased SNVs, the SVs are phased from it after the read pass, and the CpG
    // tables and the CNVs follow. A new phasing reruns the haplotagging and all after it.
    let dir = fake_site("haplotag", "[methylation]\nreport = \"cpg\"\n");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let dry = lrpipe_run(&dir, &[&args[..], &["--dry-run"]].concat()).output().expect("lrpipe runs");
    assert!(dry.status.success(), "{}", String::from_utf8_lossy(&dry.stderr));
    let plan = String::from_utf8_lossy(&dry.stdout);
    for line in [
        "# reads (after: -)",
        "# haplotag (after: snv)",
        "# sv-phase (after: reads, haplotag)",
        "# cnv (after: snv, reads, sv-phase)",
        "# methylation (after: haplotag)",
        "# merge (after: snv, sv-phase, cnv)",
    ] {
        assert!(plan.contains(line), "{line}: {plan}");
    }
    assert!(plan.contains("--list-only"), "no BAM by default: {plan}");

    let names = |manifest: &serde_json::Value| -> Vec<String> {
        manifest["stages"]
            .as_array()
            .expect("stages")
            .iter()
            .map(|s| s["name"].as_str().expect("a name").into())
            .collect()
    };
    let status = |manifest: &serde_json::Value, stage: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == stage).expect("a stage")["status"].as_str().expect("a status").to_owned()
    };
    let run = |extra: &[&str]| {
        let out = lrpipe_run(&dir, &[&args[..], extra].concat()).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        written_manifest(&dir.join("out/manifest.json"))
    };
    let manifest = run(&[]);
    assert_eq!(names(&manifest), ["snv", "reads", "haplotag", "sv-phase", "cnv", "methylation", "merge", "qc"]);
    assert!(stage_statuses(&manifest).iter().all(|s| s == "succeeded"), "{manifest}");
    let outputs = &manifest["outputs"];
    for (kind, file) in [("haplotags", "phase/S.haplotags.tsv.gz"), ("sv_vcf", "sv/S.sv.vcf"), ("snf", "sv/S.snf")] {
        assert_eq!(outputs[kind], dir.join("out").join(file).to_str().expect("UTF-8"), "{kind}: {outputs}");
    }
    assert!(outputs.get("haplotagged_reads").is_none(), "{outputs}");
    for file in ["sv/S.unphased.sv.vcf", "sv/S.phasing-evidence.gz", "methylation/S.cpg.hap2.bed.gz"] {
        assert!(dir.join("out").join(file).exists(), "{file}");
    }

    // Resumed, nothing reruns; a new phased VCF reruns the haplotagging and what reads its list
    // (and the merge, which reads the phased VCF), but not lrpipe-snv or the read pass.
    assert!(stage_statuses(&run(&[])).iter().all(|s| s == "skipped"));
    // Keep the stand-in's size and empty index so this exercises downstream content tracking.
    fs::write(dir.join("out/snv/phased_merge_output.vcf.gz"), "y\n").expect("written");
    stamp_without_outputs_state(&dir.join("out/.lrpipe/snv.done"));
    let manifest = run(&[]);
    for stage in ["snv", "reads"] {
        assert_eq!(status(&manifest, stage), "skipped", "{stage}: {manifest}");
    }
    for stage in ["haplotag", "sv-phase", "methylation", "merge"] {
        assert_eq!(status(&manifest, stage), "succeeded", "{stage}: {manifest}");
    }
    // The list gone: it is made again.
    fs::remove_file(dir.join("out/phase/S.haplotags.tsv.gz")).expect("removed");
    let manifest = run(&[]);
    for stage in ["haplotag", "sv-phase", "methylation"] {
        assert_eq!(status(&manifest, stage), "succeeded", "{stage}: {manifest}");
    }

    // The haplotagged reads: haplotag again, with the BAM.
    let manifest = run(&["--haplotagged-reads"]);
    assert_eq!(status(&manifest, "haplotag"), "succeeded", "{manifest}");
    assert_eq!(status(&manifest, "reads"), "skipped", "{manifest}");
    let bam = dir.join("out/phase/S.haplotagged.bam");
    assert_eq!(manifest["outputs"]["haplotagged_reads"], bam.to_str().expect("UTF-8"));
    assert!(bam.exists());
    let haplotag = &manifest["stages"][2]["command"].to_string();
    assert!(haplotag.contains("--write-index") && !haplotag.contains("--list-only"), "{haplotag}");

    // Without phasing there is nothing to haplotag from.
    let out = lrpipe_run(&dir, &[&args[..], &["--haplotagged-reads", UNPHASED[0], UNPHASED[1]]].concat())
        .output()
        .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("give phase among --analyses"));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn the_qc_summary_follows_what_it_reads_and_reruns_when_it_changes() {
    let _serial = serial();
    // F9: the qc stage comes last, after every other stage, and reads the files it summarises:
    // the read pass's statistics and depth, the SNVs and their phasing, the phased SVs and the
    // CpG tables; a change to one of them reruns it alone. Its files are in the manifest. Its
    // failure does not fail the run.
    let dir = fake_site("qc", "[methylation]\nreport = \"cpg\"\n");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let dry = lrpipe_run(&dir, &[&args[..], &["--dry-run"]].concat()).output().expect("lrpipe runs");
    assert!(dry.status.success(), "{}", stderr(&dry));
    let plan = String::from_utf8_lossy(&dry.stdout);
    assert!(plan.contains("# qc (after: snv, reads, haplotag, sv-phase, cnv, methylation, merge)"), "{plan}");
    let qc = plan.lines().skip_while(|l| !l.starts_with("# qc")).nth(1).expect("the qc command");
    for part in [
        " qc --sample S --sex female --reads ",
        " --read-stats out/qc/S.read-stats.json --depth-summary out/depth/S.mosdepth.summary.txt --snv-vcf \
         out/snv/merge_output.vcf.gz --phased-vcf out/snv/phased_merge_output.vcf.gz --sv-vcf out/sv/S.sv.vcf \
         --cpg-combined out/methylation/S.cpg.combined.bed.gz --cpg-hap1 out/methylation/S.cpg.hap1.bed.gz \
         --cpg-hap2 out/methylation/S.cpg.hap2.bed.gz --json out/qc/S.qc.json --tsv out/qc/S.qc.tsv",
    ] {
        let part = part.replace(" out/", &format!(" {}/out/", dir.display()));
        assert!(qc.contains(&part), "{part}\n{qc}");
    }
    let reads = plan.lines().skip_while(|l| !l.starts_with("# reads")).nth(1).expect("the read pass");
    assert!(reads.contains("--read-stats"), "{reads}");

    let run = |extra: &[&str]| {
        let out = lrpipe_run(&dir, &[&args[..], extra].concat()).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", stderr(&out));
        written_manifest(&dir.join("out/manifest.json"))
    };
    let status = |manifest: &serde_json::Value, stage: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == stage).expect("a stage")["status"].as_str().expect("a status").to_owned()
    };
    let manifest = run(&[]);
    assert_eq!(status(&manifest, "qc"), "succeeded", "{manifest}");
    for (kind, file) in [("qc_json", "qc/S.qc.json"), ("qc_tsv", "qc/S.qc.tsv"), ("read_stats", "qc/S.read-stats.json")]
    {
        assert_eq!(manifest["outputs"][kind], dir.join("out").join(file).to_str().expect("UTF-8"), "{kind}");
    }
    // A new CpG table: the QC alone reruns.
    fs::write(dir.join("out/methylation/S.cpg.hap1.bed.gz"), "other").expect("written");
    stamp_without_outputs_state(&dir.join("out/.lrpipe/methylation.done"));
    let manifest = run(&[]);
    let reran: Vec<&str> = manifest["stages"]
        .as_array()
        .expect("stages")
        .iter()
        .filter(|s| s["status"] == "succeeded")
        .map(|s| s["name"].as_str().expect("a name"))
        .collect();
    assert_eq!(reran, ["qc"]);
    // So does another sex, which only it and the stages calling chrX and chrY read.
    let mut male = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"));
    male.current_dir(&dir).args(["run", "--sample", "S", "--platform", "ont", "--sex", "male", "--memory", "64"]);
    let out = male.args(args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(status(&manifest, "qc"), "succeeded", "{manifest}");
    assert_eq!(status(&manifest, "reads"), "skipped", "{manifest}");

    // SNVs alone: no read pass, so the reads are not counted, and the earlier run's counts are
    // said to be left.
    let dry = lrpipe_run(&dir, &[&args[..], &["--dry-run", "--analyses", "snv"]].concat()).output().expect("runs");
    let plan = String::from_utf8_lossy(&dry.stdout);
    assert!(plan.contains("# qc (after: snv, merge)"), "{plan}");
    assert!(!plan.contains("--read-stats") && !plan.contains("--phased-vcf"), "{plan}");
    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", "snv"]].concat()).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let left = stderr(&out).lines().find(|l| l.contains("of an earlier run")).unwrap_or_default().to_owned();
    assert!(left.contains("qc/S.read-stats.json"), "{}", stderr(&out));

    // A QC that fails: recorded and warned of, but the run is complete and exits with 0.
    let out = lrpipe_run(&dir, &[&args[..], &["--rerun", "all"]].concat())
        .env("FAKE_QC_FAIL", "the summary broke")
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "complete", "{manifest}");
    let statuses = stage_statuses(&manifest);
    assert_eq!(statuses.last().map(String::as_str), Some("failed"), "{statuses:?}");
    assert!(statuses[..statuses.len() - 1].iter().all(|s| s == "succeeded"), "{statuses:?}");
    let warning = "no QC summary: its stage failed (exit code 1); see ";
    assert!(stderr(&out).contains(&format!("warning: {warning}")), "{}", stderr(&out));
    let warnings = manifest["warnings"].as_array().expect("warnings");
    assert!(warnings.iter().any(|w| w.as_str().is_some_and(|w| w.starts_with(warning))), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn hp_tags_that_are_not_integers_are_no_haplotypes() {
    let _serial = serial();
    // lrpipe-sv and lrpipe-methyl read an integer HP alone, so reads tagged HP:Z are not
    // haplotagged for them: no per-haplotype tables, and the warning says why.
    let dir = fake_site("text-hp", "[methylation]\nreport = \"cpg\"\n");
    haplotagged_bam(&dir, true);
    let out = lrpipe_run(
        &dir,
        &["--config", "site.toml", "--reads", "reads.bam", "--output", "out", UNPHASED[0], UNPHASED[1]],
    )
    .output()
    .expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let stderr = String::from_utf8_lossy(&out.stderr);
    let warning = "no integer HP tag in any of 3 mapped reads sampled, but HP tags of another type, read as none, on 1";
    assert!(stderr.contains(warning), "{stderr}");
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["sample"]["haplotypes"], "none (input not haplotagged)");
    assert!(manifest["outputs"].get("cpg_hap1").is_none(), "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_stage_killed_by_sigterm_is_stopped_with_the_run_not_failed() {
    let _serial = serial();
    // SLURM's scancel signals every process of the job, the stages with lrpipe: a stage that
    // dies of SIGTERM (or SIGINT) was stopped, and the run is interrupted, not failed.
    let root = temp_dir("sigterm");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let mut stages = diamond(&root, &input);
    stages[0] = stage(StageName::Snv, "sleep 0.3; kill -TERM $$", &[&input], &[root.join("snv/a.txt")], &[]);
    stages[1] = stage(StageName::Reads, "sleep 30; echo b > $1", &[&input], &[root.join("reads/b.txt")], &[]);
    let mut m = manifest(&input);
    let mut events = Vec::new();
    let t = Instant::now();
    let err = execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |e| events.push(e))
        .expect_err("interrupted");
    assert!(t.elapsed() < Duration::from_secs(10), "the sleeping stage was not killed");
    assert!(matches!(err, RunError::Interrupted), "{err}");
    use StageStatus::{Cancelled, NotRun};
    assert_eq!(statuses(&m), [Cancelled, Cancelled, NotRun]);
    assert_eq!(m.status, RunStatus::Interrupted);
    assert!(m.stages[0].error.as_deref().is_some_and(|e| e.contains("signal 15")), "{:?}", m.stages[0].error);
    assert!(!events.iter().any(|e| matches!(e, Event::Failed(..))), "{events:?}");
    let _ = fs::remove_dir_all(root);
}

/// Whether process `pid` is gone or a zombie, waiting a little for it.
fn ended(pid: &str) -> bool {
    let started = Instant::now();
    while started.elapsed() < Duration::from_secs(5) {
        let stat = fs::read_to_string(format!("/proc/{pid}/stat"));
        if stat.map_or(true, |s| s.rsplit_once(") ").is_some_and(|(_, rest)| rest.starts_with('Z'))) {
            return true;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    false
}

#[test]
fn what_a_stage_leaves_running_ends_with_it() {
    let _serial = serial();
    // A stage that starts a background process and exits, and one killed while its own child
    // runs: neither child outlives its stage's process group.
    let root = temp_dir("leftovers");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let (left, killed) = (root.join("left.pid"), root.join("killed.pid"));
    let script = format!("sleep 30 & echo $! > {}; echo a > $1", left.display());
    let stages = [stage(StageName::Snv, &script, &[&input], &[root.join("snv/a.txt")], &[])];
    let mut m = manifest(&input);
    execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {}).expect("runs");
    assert_eq!(statuses(&m), [StageStatus::Succeeded]);
    let pid = fs::read_to_string(&left).expect("its child's ID");
    assert!(ended(pid.trim()), "a finished stage's background process outlived it");

    let script = format!("sleep 30 & echo $! > {}; wait", killed.display());
    let stages = [stage(StageName::Reads, &script, &[&input], &[root.join("reads/b.txt")], &[])];
    let stop = AtomicBool::new(false);
    let mut m = manifest(&input);
    let err = std::thread::scope(|s| {
        s.spawn(|| {
            std::thread::sleep(Duration::from_millis(500));
            stop.store(true, Ordering::Relaxed);
        });
        execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &stop, &mut m, |_| {})
            .expect_err("interrupted")
    });
    assert!(matches!(err, RunError::Interrupted), "{err}");
    let pid = fs::read_to_string(&killed).expect("its child's ID");
    assert!(ended(pid.trim()), "a killed stage's child outlived it");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_stage_that_crashes_is_failed_not_stopped() {
    let _serial = serial();
    // Only SIGINT and SIGTERM are taken for the run being stopped.
    let root = temp_dir("crash");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    for (signal, number) in [("SEGV", 11), ("ABRT", 6)] {
        let stages = [stage(StageName::Snv, &format!("kill -{signal} $$"), &[&input], &[root.join("snv/a.txt")], &[])];
        let mut m = manifest(&input);
        let err = execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {})
            .expect_err("fails");
        let wanted = format!("killed by signal {number}");
        assert!(matches!(&err, RunError::Stage { reason, .. } if *reason == wanted), "{signal}: {err}");
        assert_eq!((statuses(&m), m.status), (vec![StageStatus::Failed], RunStatus::Failed), "{signal}");
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn no_stage_starts_once_the_run_is_asked_to_stop() {
    let _serial = serial();
    // Stopped as its first stage finished: the one after it must not start.
    let root = temp_dir("no-start");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let a = root.join("snv/a.txt");
    let stages = [
        stage(StageName::Snv, "echo a > $1", &[&input], std::slice::from_ref(&a), &[]),
        stage(StageName::Cnv, "sleep 30; echo c > $1", &[&a], &[root.join("cnv/c.txt")], &[StageName::Snv]),
    ];
    let stop = AtomicBool::new(false);
    let mut events = Vec::new();
    let mut m = manifest(&input);
    let err = execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &stop, &mut m, |e| {
        if matches!(e, Event::Succeeded(..)) {
            stop.store(true, Ordering::Relaxed);
        }
        events.push(e);
    })
    .expect_err("interrupted");
    assert!(matches!(err, RunError::Interrupted), "{err}");
    assert!(!events.iter().any(|e| matches!(e, Event::Started(StageName::Cnv, _))), "{events:?}");
    assert_eq!(statuses(&m), [StageStatus::Succeeded, StageStatus::NotRun]);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_run_that_cannot_write_its_files_kills_its_stages_and_says_so() {
    let _serial = serial();
    // The first stage puts a directory where its stamp goes, so the stamp cannot be written (by
    // root too, as in the release build's container, where a read-only directory would not stop
    // it): the other stage must be killed and the manifest finished, not left running.
    let root = temp_dir("unwritable");
    let input = root.join("reads.bam");
    fs::write(&input, "reads").expect("input");
    let layout = Layout::new(&root, "S");
    let dot = root.join(".lrpipe");
    let sleeper = root.join("sleeper.pid");
    let stages = vec![
        stage(
            StageName::Snv,
            &format!("sleep 0.5; mkdir -p {}; echo a > $1", dot.join("snv.done/x").display()),
            &[&input],
            &[root.join("snv/a.txt")],
            &[],
        ),
        stage(
            StageName::Reads,
            &format!("sleep 30 & echo $! > {}; wait", sleeper.display()),
            &[&input],
            &[root.join("reads/b.txt")],
            &[],
        ),
    ];
    let mut m = manifest(&input);
    let t = Instant::now();
    let err = execute(&stages, &layout, &lock(&layout).expect("locked"), CHANGED, &GO, &mut m, |_| {})
        .expect_err("cannot write the stamp");
    assert!(matches!(err, RunError::Io { .. }), "{err}");
    let written: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(layout.manifest()).expect("manifest")).expect("json");
    assert_eq!(written["status"], "failed", "{written}");
    assert_eq!(written["stages"][1]["status"], "cancelled", "{written}");
    assert!(t.elapsed() < Duration::from_secs(10));
    let pid = fs::read_to_string(&sleeper).expect("the stage's child's ID");
    assert!(ended(pid.trim()), "the other stage outlived the run");
    let _ = fs::remove_dir_all(root);
}

/// A stand-in for bonefish. Asked for an index's contigs (`-a INDEX -`), it prints `index.sq`
/// beside it, else the fake reference's. Otherwise it records its arguments in `bonefish.calls`,
/// writes a run file into its `--sort-tmp` prefix (whose directory must be there), "aligns" by
/// copying `aligned.bam` (of the fake reference's chr1, or as a test made it) and its index to
/// `-o`, and with `interrupt` beside it sends lrpipe SIGINT as it ends.
const FAKE_BONEFISH: &str = r#"#!/bin/sh
if [ "$1" = "--version" ]; then echo "bonefish 0.4.0"; exit 0; fi
here=$(dirname "$0")
if [ "$1" = "-a" ]; then cat "$here/index.sq" 2>/dev/null || printf '@SQ\tSN:chr1\tLN:26\n'; exit 0; fi
echo "$*" >> "$here/bonefish.calls"
if [ -e "$here/slow" ]; then
  trap 'echo INT > "$here/bonefish-saw-sigint"' INT
  sleep 30 & echo $! > "$here/bonefish-child.pid"; wait; exit 0
fi
while [ $# -gt 0 ]; do
  case "$1" in -o) out=$2; shift;; --sort-tmp) tmp=$2; shift;; esac
  shift
done
echo run > "$tmp.0000.run" || exit 1
cp "$here/aligned.bam" "$out" && cp "$here/aligned.bam.bai" "$out.bai" || exit 1
if [ -e "$here/interrupt" ]; then kill -INT $PPID; fi
"#;

/// An unaligned BAM, as dorado writes one: `read_groups` in its header, no `@SQ` lines, and one
/// unmapped read with base modification tags.
fn unaligned_bam(path: &Path, read_groups: &[&str]) {
    use pasta::sam::alignment::io::Write as _;
    use pasta::sam::alignment::record::{Flags, data::field::Tag};
    use pasta::sam::alignment::record_buf::data::field::{Value, value::Array};
    use pasta::sam::alignment::record_buf::{RecordBuf, Sequence};
    let mut text = "@HD\tVN:1.6\tSO:unknown\n".to_owned();
    for rg in read_groups {
        text.push_str(rg);
        text.push('\n');
    }
    let header: pasta::sam::Header = text.parse().expect("header");
    let mut writer = pasta::bam::io::Writer::new(fs::File::create(path).expect("BAM"));
    writer.write_header(&header).expect("header");
    let tags = [
        (Tag::new(b'M', b'M'), Value::String("C+m?,0;".into())),
        (Tag::new(b'M', b'L'), Value::Array(Array::UInt8(vec![230]))),
    ];
    let record = RecordBuf::builder()
        .set_name("r1")
        .set_flags(Flags::UNMAPPED)
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_data(tags.into_iter().collect())
        .build();
    writer.write_alignment_record(&header, &record).expect("record");
    writer.try_finish().expect("finished");
}

/// [`fake_site`] with [`FAKE_BONEFISH`] as its bonefish, two unaligned BAMs of the fixture's read
/// group (`flowcell1.bam`, `flowcell2.bam`), what the stand-in "aligns" them into
/// (`aligned.bam`), and Clair3's models in `[snv] models`.
fn fake_aligning_site(name: &str, extra: &str) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let dir = fake_site(name, "");
    // The reference of what the stand-in aligns to: chr1 alone, without the fixture reads' chr20.
    fs::write(dir.join("ref.fa.fai"), "chr1\t26\t6\t26\t27\n").expect("written");
    fs::write(dir.join("fake-bonefish"), FAKE_BONEFISH).expect("written");
    fs::set_permissions(dir.join("fake-bonefish"), fs::Permissions::from_mode(0o755)).expect("executable");
    install_models(&dir, &["r1041_e82_400bps_sup_v500", "hifi_revio"]);
    let site = format!(
        "reference = \"ref.fa\"\n[tools]\nlrpipe = \"fake-lrpipe\"\nbonefish = \"fake-bonefish\"\n[snv]\nmodels = \
         \"models\"\n{extra}"
    );
    fs::write(dir.join("site.toml"), site).expect("site configuration");
    for flowcell in ["flowcell1.bam", "flowcell2.bam"] {
        unaligned_bam(&dir.join(flowcell), &[ONT_READ_GROUP]);
    }
    write_bam(&dir.join("aligned.bam"), &[ONT_READ_GROUP]);
    dir
}

fn bonefish_calls(dir: &Path) -> Vec<String> {
    fs::read_to_string(dir.join("bonefish.calls")).unwrap_or_default().lines().map(str::to_owned).collect()
}

#[test]
fn unaligned_reads_are_aligned_once_and_every_stage_reads_the_alignments() {
    let _serial = serial();
    // F6: two unaligned BAMs of a sample are aligned by bonefish into one BAM in the output
    // directory, first and once; every other stage reads it, the model is chosen from the
    // unaligned BAMs' header, and a resumed run, with other threads or memory, does not realign.
    let dir = fake_aligning_site("align", "");
    let args = ["--config", "site.toml", "--reads", "flowcell1.bam", "--reads", "flowcell2.bam", "--output", "out"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let calls = bonefish_calls(&dir);
    assert_eq!(calls.len(), 1, "{calls:?}");
    let (bam, prefix) = (dir.join("out/align/S.bam"), dir.join("out/align/tmp/S"));
    let want = format!(
        "-x dna -y --sort-mem 1G -o {} --sort-tmp {} {} {} {}",
        bam.display(),
        prefix.display(),
        dir.join("ref.fa").display(),
        dir.join("flowcell1.bam").display(),
        dir.join("flowcell2.bam").display()
    );
    // 64 GiB: the defaults, but for the sorter's batches; the threads are the machine's.
    let want = want.replace("--sort-mem 1G ", "");
    assert!(calls[0].starts_with("-t ") && calls[0].ends_with(&want), "{}\n{want}", calls[0]);
    assert!(bam.exists() && dir.join("out/align/S.bam.bai").exists());
    assert!(!dir.join("out/align/tmp").exists(), "the scratch directory is removed once the stage succeeds");

    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "complete");
    let stages = manifest["stages"].as_array().expect("stages");
    let names: Vec<&str> = stages.iter().map(|s| s["name"].as_str().expect("a name")).collect();
    assert_eq!(names, ["align", "snv", "reads", "haplotag", "sv-phase", "cnv", "merge", "qc"]);
    assert!(stage_statuses(&manifest).iter().all(|s| s == "succeeded"), "{manifest}");
    assert_eq!(stages[0]["status"], "succeeded", "aligned before the others were planned, by this run");
    assert!(stages[0]["wall_seconds"].is_number(), "with the alignment's measurements");
    assert!(!stderr(&out).contains("align: already done"), "{}", stderr(&out));
    for (stage, flag) in [(1, "--bam"), (2, "--input")] {
        let command = stages[stage]["command"].as_array().expect("a command");
        let at = command.iter().position(|w| w == flag).expect("the reads");
        assert_eq!(command[at + 1], bam.to_str().expect("UTF-8"), "{}", stages[stage]);
    }
    // The haplotagging, planned after the alignment, tags the aligned BAM.
    let haplotag = stages[3]["command"].as_array().expect("a command");
    assert_eq!(haplotag.last().and_then(|w| w.as_str()), bam.to_str(), "{}", stages[3]);
    assert_eq!(manifest["tools"]["bonefish"]["version"], "bonefish 0.4.0");
    assert!(manifest["tools"]["bonefish"]["sha256"].is_string() && manifest["tools"]["lrpipe"].is_object());
    assert_eq!(manifest["sample"]["reads"][0]["path"], bam.to_str().expect("UTF-8"), "the aligned reads");
    let unaligned = manifest["sample"]["unaligned_reads"].as_array().expect("the unaligned reads");
    assert_eq!(unaligned.len(), 2);
    assert_eq!(manifest["outputs"]["aligned_bam"], bam.to_str().expect("UTF-8"));
    assert_eq!(manifest["snv_model"]["by"], "snv.models", "chosen from the unaligned BAMs' @RG");
    assert!(manifest["snv_model"]["path"].as_str().expect("a model").ends_with("r1041_e82_400bps_sup_v500"));
    assert!(manifest["wall_seconds"].is_number());

    // Resumed, with other threads and less memory (which change bonefish's settings but not what
    // it writes): nothing runs again.
    for extra in [&[][..], &["--threads", "1", "--memory", "12"]] {
        let out = lrpipe_run(&dir, &[&args[..], extra].concat()).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", stderr(&out));
        let statuses = stage_statuses(&written_manifest(&dir.join("out/manifest.json")));
        assert!(statuses.iter().all(|s| s == "skipped"), "{extra:?}: {statuses:?}");
    }
    assert_eq!(bonefish_calls(&dir).len(), 1, "resuming never realigns");
    // The dry run of a run aligned before lists every stage.
    let dry = lrpipe_run(&dir, &[&args[..], &["--dry-run"]].concat()).output().expect("lrpipe runs");
    let plan = String::from_utf8_lossy(&dry.stdout);
    assert!(plan.contains("# align (after: -)") && plan.contains("# snv (after: align)"), "{plan}");

    // New reads realign; --rerun all realigns once, and reruns the rest.
    std::thread::sleep(Duration::from_millis(10));
    unaligned_bam(&dir.join("flowcell2.bam"), &[ONT_READ_GROUP, &ONT_READ_GROUP.replace("ID:run1", "ID:run2")]);
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert_eq!(bonefish_calls(&dir).len(), 2, "new reads are aligned");
    let out = lrpipe_run(&dir, &[&args[..], &["--rerun", "all"]].concat()).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert_eq!(bonefish_calls(&dir).len(), 3, "--rerun all aligns once");
    let statuses = stage_statuses(&written_manifest(&dir.join("out/manifest.json")));
    assert_eq!(statuses[1..], ["succeeded"; 7], "{statuses:?}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_dry_run_before_the_alignment_lists_the_align_stage_alone() {
    let _serial = serial();
    let dir = fake_aligning_site("align-dry", "");
    let args = ["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "out", "--dry-run", "--memory", "16"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let plan = String::from_utf8_lossy(&out.stdout);
    let (model, plan) = plan.split_once('\n').expect("lines");
    assert!(model.starts_with("# snv model: ") && model.contains("r1041_e82_400bps_sup_v500"), "{model}");
    assert!(plan.starts_with("# align (after: -)\n"), "{plan}");
    // 16 GiB: bonefish's smallest settings.
    assert!(plan.contains(" -x dna -y --sort-mem 1G -K 10M -o "), "{plan}");
    assert!(plan.contains("# then the other stages, planned from") && !plan.contains("# snv"), "{plan}");
    assert!(bonefish_calls(&dir).is_empty() && !dir.join("out").exists(), "a dry run runs nothing");

    // A configured `.mdr` not made yet is made by the run, not refused before aligning.
    let site = fs::read_to_string(dir.join("site.toml")).expect("site");
    fs::write(dir.join("site.toml"), format!("{site}[cnv]\nmetadata = \"site.mdr\"\n")).expect("site");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    // A missing index of bonefish's is named by its key before bonefish runs.
    fs::write(dir.join("site.toml"), format!("{site}[align]\nindex = \"gone.raw\"\n")).expect("site");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(
        stderr(&out).contains("the site configuration's [align] index ")
            && stderr(&out).contains("lrpipe setup --bonefish"),
        "{}",
        stderr(&out)
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn aligned_and_unaligned_reads_together_are_refused() {
    let _serial = serial();
    let dir = fake_aligning_site("align-mixed", "");
    // The aligned reads here are the fixture's, on chr20; nothing is aligned.
    fs::write(dir.join("ref.fa.fai"), "chr1\t26\t6\t26\t27\nchr20\t64444167\t33\t60\t61\n").expect("written");
    for reads in [["reads.bam", "flowcell1.bam"], ["flowcell1.bam", "reads.bam"]] {
        let args = ["--config", "site.toml", "--reads", reads[0], "--reads", reads[1], "--output", "out"];
        let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
        assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
        assert!(
            stderr(&out).contains("reads.bam is aligned and") && stderr(&out).contains("is not"),
            "{}",
            stderr(&out)
        );
    }
    // Several aligned files are merged, not aligned, and one file is not given twice.
    split_fixture(&dir, None);
    let args = [
        "--config",
        "site.toml",
        "--reads",
        "reads.bam",
        "--reads",
        "flowcell1.bam",
        "--reads",
        "./flowcell1.bam",
        "--output",
        "out",
    ];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("flowcell1.bam is given twice"), "{}", stderr(&out));
    let args = ["--config", "site.toml", "--reads", "a.bam", "--reads", "b.bam", "--output", "out", "--dry-run"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(String::from_utf8_lossy(&out.stdout).contains("# merge-reads (after: -)\n"), "{}", stderr(&out));
    assert!(bonefish_calls(&dir).is_empty() && !dir.join("out").exists(), "refused before anything runs");

    // Aligned reads never go through bonefish, which need not even be there.
    fs::remove_file(dir.join("fake-bonefish")).expect("removed");
    let out = lrpipe_run(&dir, &["--config", "site.toml", "--reads", "reads.bam", "--output", "out"])
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert!(manifest["stages"][0]["name"] == "snv" && manifest["tools"].get("bonefish").is_none(), "{manifest}");
    assert!(manifest["sample"].get("unaligned_reads").is_none());
    // Unaligned reads without bonefish: refused, saying where it was looked for.
    let out = lrpipe_run(&dir, &["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "out2"])
        .output()
        .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(
        stderr(&out).contains("fake-bonefish --version") && stderr(&out).contains("[tools] bonefish"),
        "{}",
        stderr(&out)
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn fastq_without_modification_tags_gives_no_methylation() {
    let _serial = serial();
    // F6: FASTQ carries methylation only as MM/ML tags in its comments. Without them, methyl asked
    // for is refused; the default analyses leave it out, saying so; with them, it is called.
    let dir = fake_aligning_site("align-fastq", "ont_model = \"model\"\n");
    fs::write(dir.join("plain.fastq"), "@r1 qs:f:20\nACGT\n+\n!!!!\n").expect("FASTQ");
    fs::write(dir.join("tagged.fastq"), "@r1\tqs:f:20\tMM:Z:C+m?,0;\tML:B:C,200\nACGT\n+\n!!!!\n").expect("FASTQ");
    let run = |reads: &str, extra: &[&str]| {
        let args = ["--config", "site.toml", "--reads", reads, "--output", "out", "--dry-run"];
        lrpipe_run(&dir, &[&args[..], extra].concat()).output().expect("lrpipe runs")
    };
    let out = run("plain.fastq", &["--analyses", "snv,methyl"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(
        stderr(&out).contains(
            "plain.fastq has no base modification (MM/ML) tags on its first 100 reads, so methyl would find no methylation"
        ),
        "{}",
        stderr(&out)
    );
    assert!(stderr(&out).contains("samtools fastq -T MM,ML"), "{}", stderr(&out));

    // The default analyses: no methylation, with a warning, and a run that goes on.
    fs::remove_dir_all(dir.join("out")).ok();
    let args = ["--config", "site.toml", "--reads", "plain.fastq", "--output", "out"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let warning = "plain.fastq has no base modification (MM/ML) tags on its first 100 reads: no methylation is called";
    assert_eq!(stderr(&out).matches(warning).count(), 2, "at the start and in the summary: {}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert!(manifest["outputs"].get("bedmethyl").is_none(), "{manifest}");
    assert!(manifest["warnings"].to_string().contains(warning), "{manifest}");
    // The site's model, as FASTQ has no header to choose one from.
    assert_eq!(manifest["snv_model"]["by"], "snv.ont_model");

    // FASTQ with the tags: bedMethyl, from the read pass.
    let out = lrpipe_run(&dir, &["--config", "site.toml", "--reads", "tagged.fastq", "--output", "tagged"])
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("tagged/manifest.json"));
    assert!(manifest["outputs"]["bedmethyl"].is_string(), "{manifest}");
    // Some files with them and some without: methylation from the ones with, and a warning.
    let out = run("tagged.fastq", &["--reads", "plain.fastq", "--analyses", "methyl"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains(
            "plain.fastq has no base modification (MM/ML) tags on its first 100 reads: their reads give no methylation"
        ),
        "{}",
        stderr(&out)
    );
    // The command as a shell takes it: the read group's `\t`s quoted.
    let plan = String::from_utf8_lossy(&out.stdout);
    assert!(plan.contains(r" -R '@RG\tID:S\tSM:S' "), "{plan}");

    // FASTQ alone has no header to choose an ONT model from.
    site_with_snv(&dir, "models = \"models\"\n");
    let site = fs::read_to_string(dir.join("site.toml"))
        .expect("site")
        .replace("[tools]\n", "[tools]\nbonefish = \"fake-bonefish\"\n");
    fs::write(dir.join("site.toml"), site).expect("site");
    let out = run("tagged.fastq", &[]);
    assert_eq!(out.status.code(), Some(2));
    assert!(
        stderr(&out).contains("FASTQ has no header saying how its reads were basecalled")
            && stderr(&out).contains("snv.ont_model"),
        "{}",
        stderr(&out)
    );
    // A model given for it is said to be unchecked for want of a header, not of `@RG` fields.
    let out = run("tagged.fastq", &["--model", "r1041_e82_400bps_sup_v500"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains(
            "FASTQ has no header saying what made the reads, so r1041_e82_400bps_sup_v500 (given by --model) is not \
             checked against them"
        ),
        "{}",
        stderr(&out)
    );
    // dorado's FASTQ names the basecall model in each read's RG:Z tag: the model is chosen from it.
    fs::write(
        dir.join("dorado.fastq"),
        "@r1\tqs:f:20\tRG:Z:70eb6e24d3fa466002d6394f3ad558a46ab22393_dna_r10.4.1_e8.2_400bps_sup@v5.0.0\nACGT\n+\n!!!!\n",
    )
    .expect("FASTQ");
    let out = run("dorado.fastq", &[]);
    assert!(out.status.success(), "{}", stderr(&out));
    let plan = String::from_utf8_lossy(&out.stdout);
    let model = plan.lines().next().unwrap_or_default();
    assert!(
        model.contains("models/r1041_e82_400bps_sup_v500 (snv.models: chosen from the header")
            && model.contains("RG:Z tags"),
        "{model}"
    );
    // An empty FASTQ is said to be empty.
    fs::write(dir.join("empty.fastq"), "").expect("FASTQ");
    let out = run("empty.fastq", &[]);
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("empty.fastq: empty, with no reads"), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_resumed_run_from_unaligned_reads_keeps_every_stages_measurements() {
    let _serial = serial();
    // A run that aligns again (--rerun all) writes a manifest of the align stage alone before
    // planning the rest; the other stages' measured peaks of the run before must still be theirs
    // (LRUN-02), and a resume keeps every stage's times and peaks.
    let dir = fake_aligning_site("align-measured", "");
    let args = ["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "out", "--analyses", "snv,sv"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    for stage in manifest["stages"].as_array().expect("stages") {
        assert_eq!(stage["status"], "skipped", "{stage}");
        assert!(stage["wall_seconds"].is_number() && stage["max_rss_gib"].is_number(), "kept: {stage}");
    }
    assert_eq!(bonefish_calls(&dir).len(), 1);

    // The read pass measured 60 GiB: beside lrpipe-snv's 14.5 it does not fit in 64, even after
    // a realignment.
    let mut manifest = manifest;
    let reads = manifest["stages"].as_array_mut().expect("stages").iter_mut().find(|s| s["name"] == "reads");
    reads.expect("the read pass")["max_rss_gib"] = 60.0.into();
    fs::write(dir.join("out/manifest.json"), manifest.to_string()).expect("written");
    let out = lrpipe_run(&dir, &[&args[..], &["--rerun", "all"]].concat()).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert_eq!(bonefish_calls(&dir).len(), 2, "realigned");
    assert!(stderr(&out).contains("reads: waiting for memory (expects 69.0 GiB"), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn what_would_fail_after_aligning_fails_before() {
    let _serial = serial();
    // Plan errors, an index of other contigs than the reference's, and files the plan needs stop
    // the run before bonefish runs, in a dry run too.
    let dir = fake_aligning_site("align-early", "");
    let args = ["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "out"];
    for dry in [false, true] {
        let extra: &[&str] = if dry { &["--analyses", "tr", "--dry-run"] } else { &["--analyses", "tr"] };
        let out = lrpipe_run(&dir, &[&args[..], extra].concat()).output().expect("lrpipe runs");
        assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
        assert!(stderr(&out).contains("need a [repeats] catalog"), "{}", stderr(&out));
        assert!(bonefish_calls(&dir).is_empty() && !dir.join("out").exists(), "before aligning");
    }
    // --sex auto with a reference of one contig, whose reads cannot give the sex.
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(&dir)
        .args(["run", "--sample", "S", "--platform", "ont", "--sex", "auto", "--memory", "64"])
        .args(args)
        .output()
        .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(
        stderr(&out)
            .contains("ref.fa has 1 of the 22 autosomes (chr1) and no chrX: the sex is estimated from a whole genome"),
        "{}",
        stderr(&out)
    );
    assert!(bonefish_calls(&dir).is_empty() && !dir.join("out").exists(), "before aligning");
    let site = fs::read_to_string(dir.join("site.toml")).expect("site");
    fs::write(dir.join("site.toml"), format!("{site}[repeats]\ncatalog = \"gone.bed\"\n")).expect("site");
    // A file missing is said to be before the sex is: the files are checked first.
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(&dir)
        .args(["run", "--sample", "S", "--platform", "ont", "--sex", "auto", "--memory", "64", "--analyses", "tr"])
        .args(args)
        .output()
        .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("[repeats] catalog"), "{}", stderr(&out));
    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", "tr"]].concat()).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("[repeats] catalog") && stderr(&out).contains("gone.bed"), "{}", stderr(&out));

    // An index of another reference.
    fs::write(dir.join("ref.raw"), "an index").expect("index");
    fs::write(dir.join("site.toml"), format!("{site}[align]\nindex = \"ref.raw\"\n")).expect("site");
    fs::write(dir.join("index.sq"), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n").expect("contigs");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    let e = stderr(&out);
    assert!(
        e.contains("ref.raw has other contigs than") && e.contains("chr1 of 248956422 bases against chr1 of 26"),
        "{e}"
    );
    assert!(e.contains("lrpipe setup --bonefish"), "{e}");
    assert!(bonefish_calls(&dir).is_empty() && !dir.join("out").exists(), "nothing ran");
    // The index's own contigs: it is aligned to, with the reference an input of the stage too.
    fs::remove_file(dir.join("index.sq")).expect("removed");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let calls = bonefish_calls(&dir);
    assert!(
        calls[0].contains(&format!("{} {}", dir.join("ref.raw").display(), dir.join("flowcell1.bam").display())),
        "{calls:?}"
    );
    assert!(manifest["resources"]["align_index"]["sha256"].is_string(), "{manifest}");
    // Another reference (here the same file, changed) realigns, index or not.
    std::thread::sleep(Duration::from_millis(10));
    fs::write(dir.join("ref.fa"), ">chr1\nACGTACGTNNNNNNNNNNACGTACGA\n").expect("reference");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert_eq!(bonefish_calls(&dir).len(), 2, "another reference realigns");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_failure_or_signal_after_aligning_finishes_the_manifest() {
    let _serial = serial();
    // The aligned reads are of another contig than the reference's: the run fails after the align
    // stage, before the others, and its manifest says so rather than `running`.
    let dir = fake_aligning_site("align-between", "");
    let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:26\n";
    {
        use std::io::Write;
        let mut bam = pasta::bgzf::io::Writer::new(fs::File::create(dir.join("aligned.bam")).expect("BAM"));
        bam.write_all(b"BAM\x01").expect("written");
        bam.write_all(&(header.len() as i32).to_le_bytes()).expect("written");
        bam.write_all(header.as_bytes()).expect("written");
        bam.write_all(&1i32.to_le_bytes()).expect("written");
        bam.write_all(&5i32.to_le_bytes()).expect("written");
        bam.write_all(b"chr2\0").expect("written");
        bam.write_all(&26i32.to_le_bytes()).expect("written");
        bam.finish().expect("BAM written");
    }
    let args = ["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "out"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("S.bam has other contigs than"), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "failed", "{manifest}");
    assert!(manifest["error"].as_str().expect("why").contains("other contigs"), "{manifest}");
    assert!(manifest["finished"].is_string() && manifest["wall_seconds"].is_number());
    assert_eq!(stage_statuses(&manifest), ["succeeded"]);

    // SIGINT as the align stage ends: the run ends interrupted, and no other stage starts.
    write_bam(&dir.join("aligned.bam"), &[ONT_READ_GROUP]);
    fs::write(dir.join("interrupt"), "").expect("written");
    let args = ["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "signalled"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(130), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("signalled/manifest.json"));
    assert_eq!(manifest["status"], "interrupted", "{manifest}");
    assert_eq!(stage_statuses(&manifest), ["succeeded"], "{manifest}");
    assert!(!dir.join("signalled/logs/snv.log").exists(), "no stage after the align stage started");
    let _ = fs::remove_dir_all(dir);
}

/// Whether process `pid` is gone or a zombie, waiting a little for it.
fn gone(pid: &str) -> bool {
    let started = Instant::now();
    while started.elapsed() < Duration::from_secs(5) {
        let stat = fs::read_to_string(format!("/proc/{pid}/stat"));
        if stat.map_or(true, |s| s.rsplit_once(") ").is_some_and(|(_, rest)| rest.starts_with('Z'))) {
            return true;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    false
}

#[test]
fn a_terminals_signals_while_aligning_interrupt_the_run() {
    let _serial = serial();
    use std::os::unix::process::CommandExt;
    // Ctrl-C (SIGINT to the job's process group) and a closed terminal (SIGHUP) while bonefish
    // runs: the run ends interrupted with 130, the align stage stopped, bonefish never seeing the
    // terminal's signal, and what it started killed with its group.
    for signal in ["INT", "HUP"] {
        let dir = fake_aligning_site(&format!("align-{signal}"), "");
        fs::write(dir.join("slow"), "").expect("written");
        let mut command = lrpipe_run(&dir, &["--config", "site.toml", "--reads", "flowcell1.bam", "--output", "out"]);
        command.stdout(std::process::Stdio::null()).stderr(std::process::Stdio::piped()).process_group(0);
        let child = command.spawn().expect("lrpipe runs");
        let started = Instant::now();
        while !dir.join("bonefish-child.pid").exists() {
            assert!(started.elapsed() < Duration::from_secs(30), "bonefish never started");
            std::thread::sleep(Duration::from_millis(50));
        }
        std::thread::sleep(Duration::from_millis(200));
        let group = format!("-{}", child.id());
        let status = std::process::Command::new("kill").args([&format!("-{signal}"), "--", &group]).status();
        assert!(status.expect("kill runs").success());
        let out = child.wait_with_output().expect("lrpipe ends");
        assert_eq!(out.status.code(), Some(130), "{signal}: {}", stderr(&out));
        let manifest = written_manifest(&dir.join("out/manifest.json"));
        assert_eq!(manifest["status"], "interrupted", "{signal}: {manifest}");
        assert_eq!(stage_statuses(&manifest), ["cancelled"], "{signal}: {manifest}");
        assert!(!dir.join("bonefish-saw-sigint").exists(), "the terminal's SIGINT reached bonefish");
        let pid = fs::read_to_string(dir.join("bonefish-child.pid")).expect("its child's ID");
        assert!(gone(pid.trim()), "{signal}: bonefish's child outlived the run");
        let _ = fs::remove_dir_all(dir);
    }
}

#[test]
fn a_sigint_while_haplotagging_interrupts_the_run() {
    let _serial = serial();
    use std::os::unix::process::CommandExt;
    // Ctrl-C while the haplotag stage runs: it is stopped with its process group, the run ends
    // interrupted with 130, and the stages after it do not start; resumed, it runs again.
    let dir = fake_site("haplotag-sigint", "");
    fs::write(dir.join("slow-haplotag"), "").expect("written");
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out"];
    let mut command = lrpipe_run(&dir, &args);
    command.stdout(std::process::Stdio::null()).stderr(std::process::Stdio::piped()).process_group(0);
    let child = command.spawn().expect("lrpipe runs");
    let started = Instant::now();
    while !dir.join("haplotag.pid").exists() {
        assert!(started.elapsed() < Duration::from_secs(30), "haplotag never started");
        std::thread::sleep(Duration::from_millis(50));
    }
    std::thread::sleep(Duration::from_millis(200));
    let status = std::process::Command::new("kill").args(["-INT", "--", &format!("-{}", child.id())]).status();
    assert!(status.expect("kill runs").success());
    let out = child.wait_with_output().expect("lrpipe ends");
    assert_eq!(out.status.code(), Some(130), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "interrupted", "{manifest}");
    let status = |name: &str| {
        let stages = manifest["stages"].as_array().expect("stages");
        stages.iter().find(|s| s["name"] == name).expect("a stage")["status"].as_str().expect("a status").to_owned()
    };
    assert_eq!(status("haplotag"), "cancelled", "{manifest}");
    assert_eq!(status("sv-phase"), "not_run", "{manifest}");
    let pid = fs::read_to_string(dir.join("haplotag.pid")).expect("its ID");
    assert!(gone(pid.trim()), "the haplotag stage outlived the run");

    fs::remove_file(dir.join("slow-haplotag")).expect("removed");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    assert_eq!(manifest["status"], "complete");
    let _ = fs::remove_dir_all(dir);
}

/// The ONT fixture's reads split in two by read name (every alignment of a read in the same
/// file), as `a.bam` and `b.bam` (or with `cram`, `b.cram`, against `reference`) in `dir`, each
/// indexed, with the fixture's header.
fn split_fixture(dir: &Path, cram: Option<&Path>) -> [PathBuf; 2] {
    split_bam(&fixture_reads(), dir, cram)
}

/// The BAM `source` split as [`split_fixture`] splits the ONT fixture.
fn split_bam(source: &Path, dir: &Path, cram: Option<&Path>) -> [PathBuf; 2] {
    use pasta::{bam, cram as cram_io};
    let mut reader = bam::io::Reader::new(fs::File::open(source).expect("the BAM"));
    let header = reader.read_header().expect("a header");
    let (a, b) = (dir.join("a.bam"), dir.join(if cram.is_some() { "b.cram" } else { "b.bam" }));
    let mut first = bam::io::Writer::new(fs::File::create(&a).expect("a.bam"));
    first.write_header(&header).expect("a header");
    let reference = cram.map(|r| lrpipe_reads::Reference::open(r).expect("the reference"));
    let mut second: Box<dyn pasta::sam::alignment::io::Write> = match &reference {
        Some(reference) => Box::new(
            cram_io::io::writer::Builder::default()
                .set_reference_sequence_repository(reference.repository().clone())
                .build_from_path(&b)
                .expect("b.cram"),
        ),
        None => Box::new(bam::io::Writer::new(fs::File::create(&b).expect("b.bam"))),
    };
    second.write_alignment_header(&header).expect("a header");
    let mut record = bam::Record::default();
    let mut counts = [0, 0];
    while reader.read_record(&mut record).expect("a record") > 0 {
        let name = record.name().expect("a name");
        let half = name.iter().map(|&c| usize::from(c)).sum::<usize>() % 2;
        counts[half] += 1;
        if half == 0 {
            first.write_record(&header, &record).expect("written");
        } else {
            second.write_alignment_record(&header, &record).expect("written");
        }
    }
    assert!(counts[0] > 20 && counts[1] > 20, "both halves have reads: {counts:?}");
    first.try_finish().expect("finished");
    second.finish(&header).expect("finished");
    drop((first, second));
    for path in [&a, &b] {
        if path.extension().is_some_and(|e| e == "cram") {
            let index = cram_io::fs::index(path).expect("indexed");
            cram_io::crai::fs::write(Layout::index(path, "crai"), &index).expect("index written");
        } else {
            let bam::Index::Bai(index) = bam::fs::index(path).expect("indexed") else { panic!("a BAI") };
            bam::bai::fs::write(Layout::index(path, "bai"), &index).expect("index written");
        }
    }
    [a, b]
}

#[test]
fn several_reads_are_merged_by_a_first_stage_that_resume_skips() {
    let _serial = serial();
    // F5: two BAMs of the sample: a merge-reads stage first, which every other stage follows
    // and reads the output of; resumed, nothing reruns, and an input's index replaced reruns
    // the merge and so the rest.
    let dir = fake_site("several-reads", "");
    split_fixture(&dir, None);
    let args = ["--config", "site.toml", "--reads", "a.bam", "--reads", "b.bam", "--output", "out"];
    let dry = lrpipe_run(&dir, &[&args[..], &["--dry-run"]].concat()).output().expect("lrpipe runs");
    assert!(dry.status.success(), "{}", stderr(&dry));
    let plan = String::from_utf8_lossy(&dry.stdout);
    let merged = dir.join("out/reads/S.bam");
    assert!(
        plan.contains(&format!(
            "# merge-reads (after: -)\n{} merge-reads --output {} --threads",
            dir.join("fake-lrpipe").display(),
            merged.display()
        )),
        "{plan}"
    );
    assert!(plan.contains(&format!("{} {}\n", dir.join("a.bam").display(), dir.join("b.bam").display())), "{plan}");
    assert!(plan.contains("# snv (after: merge-reads)") && plan.contains("# reads (after: merge-reads)"), "{plan}");
    assert!(plan.contains(&format!("--bam {} ", merged.display())), "{plan}");
    assert!(plan.contains(&format!("--input {} ", merged.display())), "{plan}");

    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let names: Vec<&str> =
        manifest["stages"].as_array().expect("stages").iter().map(|s| s["name"].as_str().expect("a name")).collect();
    assert_eq!(names[0], "merge-reads", "{names:?}");
    assert!(stage_statuses(&manifest).iter().all(|s| s == "succeeded"), "{manifest}");
    assert_eq!(manifest["outputs"]["merged_reads"], merged.to_str().expect("UTF-8"));
    let reads = manifest["sample"]["reads"].as_array().expect("both files recorded");
    assert!(reads[1]["path"].as_str().expect("a path").ends_with("/b.bam"), "{reads:?}");
    // With phasing the run haplotags the merged BAM itself: no haplotype warning; without, the
    // warning names both files.
    assert!(!stderr(&out).contains("haplotagged"), "{}", stderr(&out));
    let haplotag = manifest["stages"].as_array().expect("stages").iter().find(|s| s["name"] == "haplotag");
    let haplotag = haplotag.expect("a haplotag stage")["command"].as_array().expect("a command").clone();
    assert_eq!(haplotag.last().and_then(|w| w.as_str()), merged.to_str(), "it tags the merged BAM: {haplotag:?}");
    let unphased = lrpipe_run(&dir, &[&args[..], &UNPHASED, &["--dry-run"]].concat()).output().expect("lrpipe runs");
    assert!(stderr(&unphased).contains("a.bam, "), "the haplotype warning names both: {}", stderr(&unphased));

    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let resumed = written_manifest(&dir.join("out/manifest.json"));
    assert!(stage_statuses(&resumed).iter().all(|s| s == "skipped"));
    // The earlier manifest, merge-reads and all, is read back for the stages' times and peaks.
    for stage in resumed["stages"].as_array().expect("stages") {
        assert!(stage["wall_seconds"].is_number() && stage["max_rss_gib"].is_number(), "kept: {stage}");
    }

    // An index of other contents (BAI's optional count of unplaced reads, appended).
    let mut bai = fs::read(dir.join("b.bam.bai")).expect("the index");
    bai.extend(0u64.to_le_bytes());
    fs::write(dir.join("b.bam.bai"), bai).expect("index replaced");
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let statuses = stage_statuses(&written_manifest(&dir.join("out/manifest.json")));
    assert_eq!(statuses[0], "succeeded", "the merge reruns: {statuses:?}");
    // The stand-in writes the same merged file, so the stages reading it have nothing new; but
    // the QC, which reads the inputs' indexes for the sex.
    let last = statuses.len() - 1;
    assert!(statuses[1..last].iter().all(|s| s == "skipped") && statuses[last] == "succeeded", "{statuses:?}");

    // One file: no merge, as before.
    let one = lrpipe_run(&dir, &["--config", "site.toml", "--reads", "a.bam", "--output", "one", "--dry-run"])
        .output()
        .expect("lrpipe runs");
    assert!(one.status.success(), "{}", stderr(&one));
    assert!(!String::from_utf8_lossy(&one.stdout).contains("merge-reads"));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn an_aligned_bam_requires_explicit_aliases_for_each_other_sample_name() {
    let _serial = serial();
    let dir = fake_site("sample-aliases", "");
    let first_alias = ONT_READ_GROUP.replace("ID:run1", "ID:run2").replace("SM:S", "SM:S-old");
    let second_alias = ONT_READ_GROUP.replace("ID:run1", "ID:run3").replace("SM:S", "SM:S-cell3");
    write_bam(&dir.join("aliased.bam"), &[ONT_READ_GROUP, &first_alias, &second_alias]);
    let run = |extra: &[&str]| {
        let args = ["--config", "site.toml", "--reads", "aliased.bam", "--analyses", "sv", "--output", "out"];
        lrpipe_run(&dir, &args).args(extra).output().expect("lrpipe runs")
    };

    for extra in [&[][..], &["--dry-run"], &["--sample-alias", "S-old"]] {
        let out = run(extra);
        assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
        let error = stderr(&out);
        assert!(error.contains("several samples (@RG SM: S, S-old, S-cell3)"), "{error}");
        assert!(error.contains("give --sample-alias"), "{error}");
        assert!(!dir.join("out").exists(), "refused before writing anything");
    }

    let out = run(&["--sample-alias", "S-old", "--sample-alias", "S-cell3", "--reads", "empty.bam"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("--sample-alias requires one aligned BAM or CRAM"), "{}", stderr(&out));

    let out = run(&["--sample-alias", "S-old", "--sample-alias", "S-cell3"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("as aliases for --sample S"), "{}", stderr(&out));
    let manifest_path = dir.join("out/manifest.json");
    let manifest = written_manifest(&manifest_path);
    assert_eq!(manifest["status"], "complete");
    let out = run(&[]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert_eq!(written_manifest(&manifest_path), manifest, "resume still requires the alias opt-in");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn reads_that_cannot_be_merged_are_refused() {
    let _serial = serial();
    // F5: other @SQ lines, or a read group of one ID with other fields, are refused before
    // anything is written, naming the first difference; the same read group twice is one, and
    // the same reads twice are refused, and a few of them warned of.
    let dir = fake_site("unmergeable", "");
    split_fixture(&dir, None);
    write_bam(&dir.join("other.bam"), &[ONT_READ_GROUP]);
    let run = |reads: &[&str]| {
        let mut args = vec!["--config", "site.toml", "--output", "out", "--analyses", "sv", "--dry-run"];
        for r in reads {
            args.extend(["--reads", r]);
        }
        lrpipe_run(&dir, &args).output().expect("lrpipe runs")
    };
    let out = run(&["a.bam", "other.bam"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(
        stderr(&out).contains(
            "are aligned to different references (@SQ line 1 is chr1 of 248956422 bp in the first and chr1 of 26 bp in \
             the other): every --reads must be aligned to the same one"
        ),
        "{}",
        stderr(&out)
    );

    write_bam(&dir.join("same.bam"), &[ONT_READ_GROUP]);
    let out = run(&["other.bam", "same.bam"]);
    assert!(out.status.success(), "the same read group twice: {}", stderr(&out));
    write_bam(&dir.join("clash.bam"), &[&ONT_READ_GROUP.replace("PU:PAY87954", "PU:PAY00001")]);
    let out = run(&["other.bam", "same.bam", "clash.bam"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(
        stderr(&out).contains(&format!(
            "{} and {} both have read group run1 (@RG ID:run1), with different fields: give the read groups \
             distinct IDs (`samtools addreplacerg`), or the same fields (`samtools reheader`)",
            dir.join("other.bam").display(),
            dir.join("clash.bam").display()
        )),
        "{}",
        stderr(&out)
    );
    assert!(!dir.join("out").exists(), "refused before writing anything");

    // Another sample's flow cell.
    write_bam(&dir.join("other-sample.bam"), &[&ONT_READ_GROUP.replace("ID:run1", "ID:run9").replace("SM:S", "SM:T")]);
    let out = run(&["other.bam", "other-sample.bam"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(
        stderr(&out).contains(&format!(
            "the --reads are of several samples (@RG SM: S in {}, T in {}): `lrpipe run` calls one",
            dir.join("other.bam").display(),
            dir.join("other-sample.bam").display()
        )),
        "{}",
        stderr(&out)
    );

    // A copy of a file holds the same reads: refused.
    fs::copy(dir.join("a.bam"), dir.join("again.bam")).expect("copied");
    fs::copy(dir.join("a.bam.bai"), dir.join("again.bam.bai")).expect("copied");
    let out = run(&["a.bam", "b.bam", "again.bam"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    let error = format!(
        "reads sampled from {} are in {} too: they hold the same reads",
        dir.join("again.bam").display(),
        dir.join("a.bam").display()
    );
    assert!(stderr(&out).contains(&error), "{}", stderr(&out));

    // A few of a.bam's reads in another file: a warning, and the run goes on.
    {
        use pasta::bam;
        let mut reader = bam::io::Reader::new(fs::File::open(dir.join("b.bam")).expect("b.bam"));
        let header = reader.read_header().expect("a header");
        let mut writer = bam::io::Writer::new(fs::File::create(dir.join("some.bam")).expect("some.bam"));
        writer.write_header(&header).expect("a header");
        let mut source = bam::io::Reader::new(fs::File::open(fixture_reads()).expect("the fixture"));
        source.read_header().expect("a header");
        let mut record = bam::Record::default();
        let (mut i, mut from_a) = (0, 0);
        // The fixture's reads in order: b.bam's, and every tenth of a.bam's.
        while source.read_record(&mut record).expect("a record") > 0 {
            let half = record.name().expect("a name").iter().map(|&c| usize::from(c)).sum::<usize>() % 2;
            if half == 1 || i % 10 == 0 {
                from_a += usize::from(half == 0);
                writer.write_record(&header, &record).expect("written");
            }
            i += 1;
        }
        writer.try_finish().expect("finished");
        drop(writer);
        assert!(from_a > 0);
        let bam::Index::Bai(index) = bam::fs::index(dir.join("some.bam")).expect("indexed") else { panic!("a BAI") };
        bam::bai::fs::write(dir.join("some.bam.bai"), &index).expect("index written");
    }
    let out = run(&["a.bam", "some.bam"]);
    assert!(out.status.success(), "{}", stderr(&out));
    // Of a.bam's sample, the smaller.
    let warning = format!(
        "reads sampled from {} are in {} too: a read given twice counts twice",
        dir.join("a.bam").display(),
        dir.join("some.bam").display()
    );
    assert!(stderr(&out).contains(&warning), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

/// The homopolymers of 12 bases or more of `reference`'s chr20 from `start` (the first 30 of the
/// next 150 kb, where the ONT fixture's reads are from chr20:30.15 Mb), as a repeat catalogue.
fn homopolymer_catalog(reference: &Path, start: usize) -> String {
    let sequence = lrpipe_reads::Reference::open(reference)
        .expect("the reference")
        .sequence("chr20")
        .expect("read")
        .expect("chr20");
    let window = &AsRef::<[u8]>::as_ref(sequence.as_ref())[start..start + 150_000];
    let mut catalog = String::new();
    let mut i = 0;
    while i < window.len() {
        let run = window[i..].iter().take_while(|&&b| b == window[i]).count();
        if run >= 12 && window[i] != b'N' && catalog.lines().count() < 30 {
            let at = start + i;
            catalog += &format!("chr20\t{at}\t{}\tID=h{at};MOTIFS={};STRUC=(x)n\n", at + run, window[i] as char);
        }
        i += run;
    }
    assert!(catalog.lines().count() > 5, "repeat loci: {catalog}");
    catalog
}

/// A VCF's text (bgzipped or not).
fn vcf_text(path: &Path) -> String {
    use std::io::Read;
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    if path.extension().is_some_and(|e| e == "gz") {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut text).expect("gzip");
        text
    } else {
        String::from_utf8(bytes).expect("UTF-8")
    }
}

/// A VCF's records, and its header lines but those naming the command, the run's files or the
/// date.
fn vcf_comparable(text: &str) -> Vec<&str> {
    let named = ["command", "cmdline", "fileDate", "/tmp/"];
    text.lines().filter(|l| !(l.starts_with("##") && named.iter().any(|n| l.contains(n)))).collect()
}

/// A VCF's calls: where, which alleles, the filter and the genotype (not the quality).
fn vcf_calls(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let gt = f.get(9).and_then(|s| s.split(':').next()).unwrap_or_default();
            format!("{} {} {} {} {} {gt}", f[0], f[1], f[3], f[4], f[6])
        })
        .collect()
}

#[test]
fn split_reads_give_the_outputs_of_the_whole() {
    let _serial = serial();
    // F5 acceptance, with the real tools: the ONT fixture split in two by read name, as a BAM
    // and a CRAM, gives what the fixture's BAM gives through `lrpipe run`, byte for byte but for
    // the headers naming the command, files and date: the read pass's SVs, depth, repeats and
    // bedMethyl, and the SNVs and CNVs. The fixture uses a bgzipped reference. The SNVs need
    // SNV_MODEL, the Clair3 model r1041_e82_400bps_sup_v500 (tools/check.sh sets it);
    // without it they are skipped, saying so.
    let dir = temp_dir("split-whole");
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let reference = fixture.join("chr20_window.fa.gz");
    split_fixture(&dir, Some(&reference));
    fs::write(dir.join("catalog.bed"), homopolymer_catalog(&reference, 30_150_000)).expect("the catalog");
    let model = std::env::var_os("SNV_MODEL").map(PathBuf::from).filter(|m| m.join("pileup.pt").exists());
    let snv = match &model {
        Some(model) => format!("ont_model = \"{}\"\ndevice = \"cpu\"\n", model.display()),
        None => {
            eprintln!(
                "SNV_MODEL not set (or not a Clair3 model): the SNVs and CNVs of the split reads are not compared"
            );
            String::new()
        }
    };
    let site = format!("reference = \"{}\"\n[repeats]\ncatalog = \"catalog.bed\"\n[snv]\n{snv}", reference.display());
    fs::write(dir.join("site.toml"), site).expect("site configuration");
    let analyses = if model.is_some() { "snv,phase,sv,cnv,tr,depth,methyl" } else { "sv,tr,depth,methyl" };
    let whole = fixture.join("ont.bam");
    for (out, reads) in [("whole", vec![whole.to_str().expect("UTF-8")]), ("split", vec!["a.bam", "b.cram"])] {
        let mut args = vec!["--config", "site.toml", "--output", out, "--analyses", analyses, "--threads", "4"];
        for r in reads {
            args.extend(["--reads", r]);
        }
        let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"));
        command.current_dir(&dir).args(["run", "--sample", "S", "--platform", "ont", "--sex", "female"]);
        let run = command.args(["--memory", "64"]).args(&args).output().expect("lrpipe runs");
        assert!(run.status.success(), "{out}: {}", stderr(&run));
    }
    assert!(dir.join("split/reads/S.bam.bai").exists(), "the merged reads are indexed");
    let (whole, split) = (dir.join("whole"), dir.join("split"));
    let manifest = written_manifest(&whole.join("manifest.json"));
    let mut compared = 0;
    for (kind, path) in manifest["outputs"].as_object().expect("outputs") {
        let path = Path::new(path.as_str().expect("a path"));
        let relative = path.strip_prefix(&whole).expect("in the run's directory");
        let other = split.join(relative);
        compared += 1;
        match kind.as_str() {
            // Their headers name the command and the date.
            "sv_vcf" | "repeats_vcf" | "cnv_vcf" => {
                assert_eq!(vcf_comparable(&vcf_text(&other)), vcf_comparable(&vcf_text(path)), "{kind}");
                // The fixture's reads hold no SVs (the SV test below has some).
                let records = vcf_text(path).lines().any(|l| !l.starts_with('#'));
                assert!(records || kind == "sv_vcf", "{kind} has records");
            }
            "vcf" | "snv_vcf" | "phased_snv_vcf" => {
                let (text, split_text) = (vcf_text(path), vcf_text(&other));
                let calls = vcf_calls(&text);
                assert!(calls.len() > 10, "{kind}: {calls:?}");
                assert_eq!(vcf_calls(&split_text), calls, "{kind}: the calls");
                let differs = vcf_comparable(&split_text).into_iter().zip(vcf_comparable(&text)).find(|(a, b)| a != b);
                assert_eq!(differs, None, "{kind}: the records");
            }
            // The candidates' header names the command too; their calls are compared below.
            "snf" => {}
            _ => {
                let bytes = fs::read(path).expect("an output");
                assert!(!bytes.is_empty(), "{kind}");
                assert!(fs::read(&other).expect("the split's") == bytes, "{kind}: {} differs", relative.display());
            }
        }
    }
    assert!(compared >= if model.is_some() { 11 } else { 7 }, "{compared} outputs: {manifest}");
    // The SV candidates give the same calls.
    for run in [&whole, &split] {
        let vcf = run.join("sv/from-snf.vcf");
        let out = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
            .args(["sv", "--input"])
            .arg(run.join("sv/S.snf"))
            .arg("--vcf")
            .arg(&vcf)
            .output()
            .expect("lrpipe sv runs");
        assert!(out.status.success(), "{}", stderr(&out));
    }
    let from_snf = |run: &Path| vcf_comparable(&vcf_text(&run.join("sv/from-snf.vcf"))).join("\n");
    assert_eq!(from_snf(&split), from_snf(&whole), "calls from the .snf");
    let bedmethyl = vcf_text(&whole.join("methylation/S.bedmethyl.gz"));
    assert!(bedmethyl.lines().count() > 100, "bedMethyl rows from the bgzipped reference: {}", bedmethyl.len());
    let _ = fs::remove_dir_all(dir);
}

/// Where the ONT fixture's chr20 starts in [`two_contig_fixture`]'s, and how long that is.
const TWO_OFFSET: usize = 29_500_000;
const TWO_LENGTH: usize = 1_200_000;

/// The ONT fixture as a genome of two contigs, in `dir`: `two.fa` (with its `.fai`), chr20 the
/// fixture's 1.2 Mb from chr20:29.5 Mb, where its reads are, and chr21 the same bases; and
/// `two.bam` (indexed), the reads on chr20 at their positions less [`TWO_OFFSET`], and every
/// other one again on chr21 (without `SA` tags, which name the fixture's positions). chr21 has
/// half the depth, so its calls' qualities are lower: a run of both has other quality cut-offs
/// than a run of chr20 alone.
fn two_contig_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    use pasta::bam;
    use pasta::sam::alignment::RecordBuf;
    use pasta::sam::alignment::io::Write as _;
    use pasta::sam::alignment::record::data::field::Tag;
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let sequence = lrpipe_reads::Reference::open(&fixture.join("chr20_window.fa.gz"))
        .expect("the reference")
        .sequence("chr20")
        .expect("read")
        .expect("chr20");
    let bases = &AsRef::<[u8]>::as_ref(sequence.as_ref())[TWO_OFFSET..TWO_OFFSET + TWO_LENGTH];
    let (mut fasta, mut fai) = (Vec::new(), String::new());
    for name in ["chr20", "chr21"] {
        fasta.extend_from_slice(format!(">{name}\n").as_bytes());
        fai += &format!("{name}\t{TWO_LENGTH}\t{}\t60\t61\n", fasta.len());
        for line in bases.chunks(60) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
    }
    let reference = dir.join("two.fa");
    fs::write(&reference, fasta).expect("the reference");
    fs::write(dir.join("two.fa.fai"), fai).expect("its index");

    let mut reader = bam::io::Reader::new(fs::File::open(fixture.join("ont.bam")).expect("the fixture"));
    let header = reader.read_header().expect("a header");
    let mut text = Vec::new();
    pasta::sam::io::Writer::new(&mut text).write_header(&header).expect("the header's text");
    let read_groups: String = String::from_utf8(text)
        .expect("UTF-8")
        .lines()
        .filter(|l| l.starts_with("@RG"))
        .map(|l| format!("{l}\n"))
        .collect();
    let header: pasta::sam::Header = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20\tLN:{TWO_LENGTH}\n@SQ\tSN:chr21\tLN:{TWO_LENGTH}\n{read_groups}"
    )
    .parse()
    .expect("a header");
    let records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.expect("a record")).collect();
    let reads = dir.join("two.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&reads).expect("two.bam"));
    writer.write_header(&header).expect("a header");
    for (contig, every) in [(0, 1), (1, 2)] {
        for record in records.iter().step_by(every) {
            let mut moved = record.clone();
            *moved.reference_sequence_id_mut() = Some(contig);
            let start = usize::from(record.alignment_start().expect("aligned")) - TWO_OFFSET;
            *moved.alignment_start_mut() = Some(pasta::core::Position::try_from(start).expect("a position"));
            *moved.mate_reference_sequence_id_mut() = None;
            *moved.mate_alignment_start_mut() = None;
            moved.data_mut().remove(&Tag::OTHER_ALIGNMENTS);
            writer.write_alignment_record(&header, &moved).expect("written");
        }
    }
    writer.try_finish().expect("finished");
    drop(writer);
    let bam::Index::Bai(index) = bam::fs::index(&reads).expect("indexed") else { panic!("a BAI") };
    bam::bai::fs::write(Layout::index(&reads, "bai"), &index).expect("index written");
    (reference, reads)
}

/// The lines of `text` that are records (not `#` headers) on chr20 at a position (column 2)
/// within `[start, end)`, or every chr20 record without a window.
fn chr20_records(text: &str, window: Option<(u64, u64)>) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#') && l.starts_with("chr20\t"))
        .filter(|l| {
            let at = l.split('\t').nth(1).and_then(|p| p.parse::<u64>().ok()).expect("a position");
            window.is_none_or(|(s, e)| s <= at && at < e)
        })
        .map(str::to_owned)
        .collect()
}

#[test]
fn regions_give_the_whole_genomes_records_in_them() {
    let _serial = serial();
    // F7 acceptance, with the real tools on the ONT fixture as a genome of two contigs
    // ([`two_contig_fixture`]): a run in chr20 gives a run of both's chr20 records, byte for byte
    // but for the headers, for the depth, the SVs and the tandem repeats' loci, and as near as
    // lrpipe-snv's calling allows for the SNVs: it takes its quality cut-offs (the share of
    // pileup calls to call again with full alignment, and to phase from) from every contig it
    // calls, as Clair3 does, so some calls near those cut-offs differ, and what follows from them
    // (the phased VCF, the haplotag list, the haplotypes of the SVs and repeats). A run in a
    // 20 kb window gives the depth's windows and the repeats' loci in it, and the SNVs Clair3
    // 2.0.3 calls with that BED (`lrpipe-snv/ci/fixture/expected/window`). The bedMethyl's pass
    // thresholds are estimated from the regions' reads, so its rows are those of `lrpipe methyl
    // pileup` of the same regions. The SNVs need SNV_MODEL (tools/check.sh sets it); without it
    // they are skipped, saying so.
    let dir = temp_dir("regions-whole");
    let (reference, reads) = two_contig_fixture(&dir);
    let start = 30_150_000 - TWO_OFFSET;
    let window = (30_180_000 - TWO_OFFSET as u64, 30_200_000 - TWO_OFFSET as u64);
    // And a locus across the window's start, which its run genotypes whole.
    let catalog = homopolymer_catalog(&reference, start)
        + &format!("chr20\t{}\t{}\tID=edge;MOTIFS=AT;STRUC=(AT)n\n", window.0 - 10, window.0 + 10);
    fs::write(dir.join("catalog.bed"), &catalog).expect("the catalog");
    let model = std::env::var_os("SNV_MODEL").map(PathBuf::from).filter(|m| m.join("pileup.pt").exists());
    let snv = match &model {
        Some(model) => format!("ont_model = \"{}\"\ndevice = \"cpu\"\n", model.display()),
        None => {
            eprintln!("SNV_MODEL not set (or not a Clair3 model): the SNVs in regions are not compared");
            String::new()
        }
    };
    let site = format!("reference = \"{}\"\n[repeats]\ncatalog = \"catalog.bed\"\n[snv]\n{snv}", reference.display());
    fs::write(dir.join("site.toml"), site).expect("site configuration");
    fs::write(dir.join("chr20.bed"), format!("chr20\t0\t{TWO_LENGTH}\n")).expect("BED");
    fs::write(dir.join("window.bed"), format!("chr20\t{}\t{}\n", window.0, window.1)).expect("BED");
    let analyses = if model.is_some() { "snv,phase,sv,tr,depth,methyl" } else { "sv,tr,depth,methyl" };
    for (out, bed) in [("whole", None), ("chr20", Some("chr20.bed")), ("window", Some("window.bed"))] {
        let mut args = vec!["--config", "site.toml", "--output", out, "--analyses", analyses, "--threads", "4"];
        args.extend(["--reads", reads.to_str().expect("UTF-8")]);
        if let Some(bed) = bed {
            args.extend(["--regions", bed]);
        }
        let run = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
        assert!(run.status.success(), "{out}: {}", stderr(&run));
    }
    let manifest = written_manifest(&dir.join("whole/manifest.json"));
    let text = |run: &str, relative: &Path| vcf_text(&dir.join(run).join(relative));
    // The repeats' loci a window reaches, by the catalogue: one may start before it.
    let reached: Vec<String> = catalog
        .lines()
        .filter(|l| {
            let f: Vec<u64> = l.split('\t').skip(1).take(2).map(|n| n.parse().expect("a coordinate")).collect();
            f[0] < window.1 && window.0 < f[1]
        })
        .map(|l| {
            l.split('\t')
                .nth(3)
                .and_then(|i| i.split(';').next())
                .and_then(|i| i.strip_prefix("ID="))
                .expect("an ID")
                .to_owned()
        })
        .collect();
    assert!(reached.len() >= 2, "repeat loci in the window: {reached:?}");
    let pass_calls = |text: &str, w| -> Vec<String> {
        vcf_calls(&chr20_records(text, w).join("\n")).into_iter().filter(|c| c.contains(" PASS ")).collect()
    };
    let mut compared = 0;
    for (kind, path) in manifest["outputs"].as_object().expect("outputs") {
        // Its header names the command; its candidates are of the SVs compared.
        if kind == "snf" {
            continue;
        }
        let relative = Path::new(path.as_str().expect("a path")).strip_prefix(dir.join("whole")).expect("in the run");
        let whole = text("whole", relative);
        let (chr20, in_window) = (text("chr20", relative), text("window", relative));
        compared += 1;
        match kind.as_str() {
            "sv_vcf" => {
                assert_eq!(chr20_records(&chr20, None), chr20_records(&whole, None), "{kind}: chr20");
                assert_eq!(chr20_records(&in_window, None), chr20_records(&whole, Some(window)), "{kind}");
            }
            "repeats_vcf" => {
                let loci = |t: &str| -> Vec<String> {
                    chr20_records(t, None)
                        .iter()
                        .map(|l| l.split('\t').take(5).collect::<Vec<_>>().join("\t"))
                        .collect()
                };
                assert_eq!(loci(&chr20), loci(&whole), "{kind}: chr20's loci");
                let ids: Vec<String> = chr20_records(&in_window, None)
                    .iter()
                    .map(|l| {
                        let info = l.split('\t').nth(7).expect("INFO");
                        info.split(';').find_map(|f| f.strip_prefix("TRID=")).expect("a TRID").to_owned()
                    })
                    .collect();
                let mut reached = reached.clone();
                reached.sort_by_key(|id| ids.iter().position(|i| i == id));
                assert_eq!(ids, reached, "{kind}: the loci the window reaches");
            }
            // Without SNVs, the merged VCF has the SVs and repeats, as each.
            "vcf" if model.is_none() => {
                assert_eq!(chr20_records(&chr20, None), chr20_records(&whole, None), "{kind}: chr20");
            }
            "snv_vcf" | "phased_snv_vcf" | "vcf" => {
                let (theirs, ours) = (pass_calls(&whole, None), pass_calls(&chr20, None));
                let shared = ours.iter().filter(|c| theirs.contains(c)).count();
                assert!(
                    theirs.len() > 400 && shared * 100 >= theirs.len().max(ours.len()) * 98,
                    "{kind}: {shared} PASS calls of chr20 shared, of {} and {}",
                    theirs.len(),
                    ours.len()
                );
                if kind == "snv_vcf" {
                    let expected = fs::read_to_string(
                        Path::new(env!("CARGO_MANIFEST_DIR"))
                            .join("../lrpipe-snv/ci/fixture/expected/window/merge_output.vcf"),
                    )
                    .expect("Clair3's calls");
                    let shifted: Vec<String> = pass_calls(&expected, None)
                        .into_iter()
                        .map(|c| {
                            let mut f: Vec<String> = c.split(' ').map(str::to_owned).collect();
                            f[1] = (f[1].parse::<usize>().expect("a position") - TWO_OFFSET).to_string();
                            f.join(" ")
                        })
                        .collect();
                    assert!(shifted.len() > 50, "{shifted:?}");
                    assert_eq!(pass_calls(&in_window, None), shifted, "{kind}: Clair3 2.0.3's calls with the window");
                }
            }
            "depth_windows" => {
                assert_eq!(chr20_records(&chr20, None), chr20_records(&whole, None), "{kind}: chr20");
                let cut: Vec<String> = chr20_records(&whole, Some(window));
                assert!(cut.len() == 20 && chr20_records(&in_window, None) == cut, "{kind}: {cut:?}");
            }
            "depth_summary" => {
                let chr20_rows =
                    |t: &str| t.lines().filter(|l| l.starts_with("chr20")).map(str::to_owned).collect::<Vec<_>>();
                assert_eq!(chr20_rows(&chr20), chr20_rows(&whole), "{kind}");
                // The window's own row alone: chr20's depth is not counted.
                let rows: Vec<&str> = in_window.lines().map(|l| l.split('\t').next().expect("a row")).collect();
                assert_eq!(rows, ["chrom", "chr20_region", "total_region"], "{in_window}");
            }
            "haplotags" => {
                let names: Vec<&str> =
                    in_window.lines().skip(1).map(|l| l.split('\t').next().expect("a name")).collect();
                let unique: std::collections::BTreeSet<&str> = names.iter().copied().collect();
                assert!(!names.is_empty() && unique.len() == names.len(), "each read listed once: {names:?}");
            }
            "bedmethyl" => {
                let pileup = |args: &[&str]| {
                    let out = dir.join("pileup.bed");
                    let run = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
                        .args(["methyl", "pileup"])
                        .arg(&reads)
                        .arg(&out)
                        .arg("--ref")
                        .arg(&reference)
                        .args(["--modified-bases", "5mC", "5hmC", "--cpg", "--combine-strands"])
                        .args(args)
                        .current_dir(&dir)
                        .output()
                        .expect("lrpipe methyl runs");
                    assert!(run.status.success(), "{}", stderr(&run));
                    fs::read_to_string(out).expect("a bedMethyl")
                };
                assert!(chr20.lines().count() > 100, "{kind}");
                assert_eq!(chr20, pileup(&["--region", "chr20"]), "{kind}: chr20");
                assert_eq!(in_window, pileup(&["--include-bed", "window.bed"]), "{kind}: the window");
            }
            // F9: the read pass's counts and the QC summary are of the regions, and say so.
            "read_stats" => {
                let json = |t: &str| serde_json::from_str::<serde_json::Value>(t).expect("JSON");
                let (whole, chr20, in_window) = (json(&whole), json(&chr20), json(&in_window));
                let restricted = [&whole, &chr20, &in_window].map(|v| v["restricted"].as_bool());
                assert_eq!(restricted, [Some(false), Some(true), Some(true)], "{kind}");
                let reads = [&whole, &chr20, &in_window].map(|v| v["reads"].as_u64().expect("reads"));
                assert!(reads[0] >= reads[1] && reads[1] > reads[2] && reads[2] > 0, "{kind}: {reads:?}");
            }
            "qc_json" => {
                let json = |t: &str| serde_json::from_str::<serde_json::Value>(t).expect("JSON");
                let (whole, chr20, in_window) = (json(&whole), json(&chr20), json(&in_window));
                assert!(whole["regions"].is_null(), "{kind}: {whole}");
                assert_eq!(chr20["regions"], serde_json::json!({"regions": 1, "bases": TWO_LENGTH}), "{kind}");
                let bases = window.1 - window.0;
                assert_eq!(in_window["regions"], serde_json::json!({"regions": 1, "bases": bases}), "{kind}");
                assert!(whole["reads"]["mapped_percent"].is_number(), "{kind}: {whole}");
                assert!(in_window["reads"]["mapped_percent"].is_null(), "{kind}: {in_window}");
                let depth = |v: &serde_json::Value| v["coverage"]["per_autosome"]["chr20"].clone();
                assert_eq!(depth(&chr20), depth(&whole), "{kind}: chr20");
                assert!(chr20["coverage"]["per_autosome"]["chr21"].is_null(), "{kind}: {chr20}");
            }
            "qc_tsv" => {
                assert!(
                    whole.contains("regions.regions\tNA\n") && in_window.contains("regions.regions\t1\n"),
                    "{kind}"
                );
            }
            other => panic!("{other} is not compared"),
        }
    }
    assert!(compared >= if model.is_some() { 9 } else { 6 }, "{compared} outputs: {manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[path = "../../lrpipe-sv/tests/common/mod.rs"]
mod sv_common;

#[test]
fn svs_in_regions_are_sniffles2s() {
    let _serial = serial();
    // F7: the SVs of a run in regions are Sniffles 2.8.1's with `--regions`
    // (lrpipe-sv/tests/data/expected/regions_run.vcf), on lrpipe-sv's fixture (HG002 ONT reads
    // around chr17:46.2-46.4 Mb): an SV is called from the signals that start in a region, so of
    // the duplication at 46,205,411-46,209,199, which reaches the first region from before it,
    // nothing, of the insertion in it and the deletion in the second their calls, and of the
    // deletion at 46,237,501, between them, nothing.
    let dir = temp_dir("regions-svs");
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-sv/tests/data");
    for (from, to) in
        [("sample1.bam", "reads.bam"), ("sample1.bam.bai", "reads.bam.bai"), ("run_regions.bed", "regions.bed")]
    {
        fs::copy(data.join(from), dir.join(to)).expect("copied");
    }
    sv_common::write_reference(&data.join("reference.window.fa"), &dir.join("ref.fa"), &[]);
    fs::write(dir.join("site.toml"), "reference = \"ref.fa\"\n[snv]\n").expect("site configuration");
    let args =
        ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--analyses", "sv", "--threads", "4"];
    let out = lrpipe_run(&dir, &args).args(["--regions", "regions.bed"]).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let records =
        |text: &str| -> Vec<String> { text.lines().filter(|l| !l.starts_with('#')).map(str::to_owned).collect() };
    let expected = records(&fs::read_to_string(data.join("expected/regions_run.vcf")).expect("Sniffles2's calls"));
    let called: Vec<&str> = expected.iter().map(|l| l.split('\t').nth(1).expect("a position")).collect();
    assert_eq!(called, ["46210499", "46277765"], "{expected:?}");
    assert_eq!(records(&vcf_text(&dir.join("out/sv/S.sv.vcf"))), expected);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn split_reads_give_the_svs_of_the_whole() {
    let _serial = serial();
    // F5: lrpipe-sv's fixture (HG002 ONT reads around chr17:46.2-46.4 Mb, haplotagged, with SVs)
    // split in two by read name and merged by `lrpipe merge-reads` gives lrpipe-sv's phased
    // calls from the whole BAM, record for record.
    let dir = temp_dir("split-svs");
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-sv/tests/data/sample1.bam");
    let [a, b] = split_bam(&source, &dir, None);
    let merged = dir.join("merged.bam");
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .args(["merge-reads", "--threads", "3", "--output"])
        .args([&merged, &a, &b])
        .output()
        .expect("lrpipe merge-reads runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(Layout::index(&merged, "bai").exists(), "indexed");
    let calls = |reads: &Path, name: &str| {
        let vcf = dir.join(name);
        let out = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
            .args(["sv", "--phase", "--threads", "2", "--input"])
            .arg(reads)
            .arg("--vcf")
            .arg(&vcf)
            .output()
            .expect("lrpipe sv runs");
        assert!(out.status.success(), "{}", stderr(&out));
        vcf_text(&vcf).lines().filter(|l| !l.starts_with('#')).map(str::to_owned).collect::<Vec<_>>()
    };
    let whole = calls(&source, "whole.vcf");
    assert!(whole.len() > 5 && whole.iter().any(|l| l.contains('|')), "phased SV calls: {whole:?}");
    assert_eq!(calls(&merged, "merged.vcf"), whole);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn every_inputs_read_groups_and_haplotags_are_read() {
    let _serial = serial();
    // F5 with F3: the model is chosen from every input's read groups: two flow cells of one
    // basecall model take its model, of two models no one model. With F1a: each input is
    // sampled for HP tags, and one haplotagged beside one not is warned of.
    let dir = fake_site("several-read-groups", "");
    install_models(&dir, &["r1041_e82_400bps_sup_v500", "r1041_e82_400bps_sup_v520"]);
    site_with_snv(&dir, "models = \"models\"\n");
    let second = ONT_READ_GROUP.replace("ID:run1", "ID:run2");
    write_bam(&dir.join("one.bam"), &[ONT_READ_GROUP]);
    write_bam(&dir.join("two.bam"), &[&second]);
    write_bam(&dir.join("newer.bam"), &[&second.replace("sup@v5.0.0", "sup@v5.2.0")]);
    let dry = |reads: &[&str], analyses: &str| {
        let mut args = vec!["--config", "site.toml", "--output", "out", "--analyses", analyses, "--dry-run"];
        for r in reads {
            args.extend(["--reads", r]);
        }
        lrpipe_run(&dir, &args).output().expect("lrpipe runs")
    };
    let out = dry(&["one.bam", "two.bam"], "snv");
    assert!(out.status.success(), "{}", stderr(&out));
    let plan = String::from_utf8_lossy(&out.stdout);
    assert!(plan.lines().next().is_some_and(|l| l.contains("models/r1041_e82_400bps_sup_v500 ")), "{plan}");
    let out = dry(&["one.bam", "newer.bam"], "snv");
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr(&out).contains("so no one model fits them all: give the model with --model"), "{}", stderr(&out));

    haplotagged_bam(&dir, false);
    let out = dry(&["reads.bam", "one.bam"], "sv");
    assert!(out.status.success(), "{}", stderr(&out));
    let warning = format!(
        "{} is haplotagged and {} not (no HP tag in their sampled reads): the SVs' phases come from the haplotagged \
         reads alone",
        dir.join("reads.bam").display(),
        dir.join("one.bam").display()
    );
    assert!(stderr(&out).contains(&warning), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_merge_never_writes_over_an_input_nor_leaves_a_partial_bam() {
    let _serial = serial();
    // F5: reads merged where one of them is (`--reads out/reads/S.bam` with `--output out`, or
    // `lrpipe merge-reads --output a.bam a.bam b.bam`) are refused before anything is written;
    // a merge that fails leaves nothing that looks like a merged BAM.
    let dir = fake_site("merge-into-input", "");
    split_fixture(&dir, None);
    fs::create_dir_all(dir.join("out/reads")).expect("directory");
    fs::copy(dir.join("a.bam"), dir.join("out/reads/S.bam")).expect("copied");
    fs::copy(dir.join("a.bam.bai"), dir.join("out/reads/S.bam.bai")).expect("copied");
    let before = fs::read(dir.join("out/reads/S.bam")).expect("the input");
    let out = lrpipe_run(
        &dir,
        &[
            "--config",
            "site.toml",
            "--reads",
            "out/reads/S.bam",
            "--reads",
            "b.bam",
            "--output",
            "./out",
            "--analyses",
            "sv",
        ],
    )
    .output()
    .expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("would overwrite the input"), "{}", stderr(&out));
    assert!(!dir.join("out/manifest.json").exists(), "refused before running");

    let merge = |args: &[&str]| {
        std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
            .current_dir(&dir)
            .arg("merge-reads")
            .args(args)
            .output()
            .expect("lrpipe merge-reads runs")
    };
    let out = merge(&["--output", "out/reads/S.bam", "b.bam", "out/reads/S.bam"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("would overwrite the input"), "{}", stderr(&out));
    assert_eq!(fs::read(dir.join("out/reads/S.bam")).expect("the input"), before, "untouched");

    // A CRAM input without a reference fails the merge after it has begun writing.
    fs::copy(fixture_reads().with_extension("cram"), dir.join("c.cram")).expect("copied");
    fs::copy(Layout::index(&fixture_reads().with_extension("cram"), "crai"), dir.join("c.cram.crai")).expect("copied");
    let out = merge(&["--output", "m.bam", "a.bam", "c.cram"]);
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("needs a reference"), "{}", stderr(&out));
    let left: Vec<String> = fs::read_dir(&dir)
        .expect("listed")
        .map(|e| e.expect("an entry").file_name().to_string_lossy().into_owned())
        .filter(|n| n.contains("m.bam"))
        .collect();
    assert!(left.is_empty(), "nothing of the failed merge is left: {left:?}");

    let out = merge(&["--output", "m.bam", "a.bam", "b.bam"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(dir.join("m.bam").exists() && dir.join("m.bam.bai").exists() && !dir.join(".m.bam.merging").exists());
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn reads_in_two_files_are_found_whatever_else_one_holds() {
    let _serial = serial();
    // F5: the files are sampled on the same contigs: a file whose first contigs hold other reads
    // (here copies of the fixture's, renamed, on chr1 to chr4) is still found to share chr20's
    // reads with another.
    use pasta::bam;
    use pasta::sam::alignment::RecordBuf;
    use pasta::sam::alignment::io::Write as _;
    let dir = fake_site("shared-reads", "");
    // The reads are on chr1 to chr4 and chr20 too, of the fixture header's lengths.
    let fai: String = [("chr1", 248_956_422), ("chr2", 242_193_529), ("chr3", 198_295_559), ("chr4", 190_214_555)]
        .iter()
        .chain(&[("chr20", 64_444_167)])
        .map(|(name, len)| format!("{name}\t{len}\t6\t60\t61\n"))
        .collect();
    fs::write(dir.join("ref.fa.fai"), fai).expect("written");
    split_fixture(&dir, None);
    let mut reader = bam::io::Reader::new(fs::File::open(fixture_reads()).expect("the fixture"));
    let header = reader.read_header().expect("a header");
    let records: Vec<RecordBuf> = reader.record_bufs(&header).map(|r| r.expect("a record")).collect();
    let spread = dir.join("spread.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&spread).expect("spread.bam"));
    writer.write_header(&header).expect("a header");
    for contig in 0..4 {
        for record in &records {
            let mut copy = record.clone();
            *copy.reference_sequence_id_mut() = Some(contig);
            *copy.name_mut() = Some(format!("copy{contig}-{}", record.name().expect("a name")).into());
            writer.write_alignment_record(&header, &copy).expect("written");
        }
    }
    for record in &records {
        writer.write_alignment_record(&header, record).expect("written");
    }
    writer.try_finish().expect("finished");
    drop(writer);
    let bam::Index::Bai(index) = bam::fs::index(&spread).expect("indexed") else { panic!("a BAI") };
    bam::bai::fs::write(Layout::index(&spread, "bai"), &index).expect("index written");

    let out = lrpipe_run(
        &dir,
        &["--config", "site.toml", "--reads", "spread.bam", "--reads", "a.bam", "--output", "out", "--dry-run"],
    )
    .output()
    .expect("lrpipe runs");
    // a.bam's reads are all in spread.bam: the same reads, refused.
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    let error = format!(
        "reads sampled from {} are in {} too: they hold the same reads",
        dir.join("a.bam").display(),
        dir.join("spread.bam").display()
    );
    assert!(stderr(&out).contains(&error), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn reads_without_modification_tags_get_no_methylation_stage() {
    let _serial = serial();
    // Reads aligned without MM/ML tags (basecalled without modified bases, or aligned from a
    // FASTQ that dropped them) gave empty CpG tables and an exit status of 0; the default
    // analyses now leave methylation out with a warning, and methyl asked for is an error.
    let dir = fake_site("no-modifications", "[methylation]\nreport = \"cpg\"\n");
    reads_bam(&dir, false, false);
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--dry-run"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let plan = String::from_utf8_lossy(&out.stdout);
    assert!(!plan.contains("# methylation"), "{plan}");
    assert!(
        stderr(&out).contains("has no base modification (MM/ML) tags on its 3 mapped reads sampled: no methylation"),
        "{}",
        stderr(&out)
    );

    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", "snv,methyl"]].concat()).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("so methyl would find no methylation"), "{}", stderr(&out));

    // With the tags, methylation is planned.
    reads_bam(&dir, false, true);
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(String::from_utf8_lossy(&out.stdout).contains("# methylation"), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn haplotag_off_phases_the_snvs_but_leaves_the_reads_tags_as_they_are() {
    let _serial = serial();
    // --haplotag off: the SNVs are phased, but the reads are not haplotagged from them, so the
    // SVs, tandem repeats and CpG tables per haplotype follow the reads' own tags and run beside
    // lrpipe-snv, as in 0.1.5; and tags of another pipeline beside the run's phasing are warned
    // of.
    let dir = fake_site("haplotag-off", "[methylation]\nreport = \"cpg\"\n");
    haplotagged_bam(&dir, false);
    let args = ["--config", "site.toml", "--reads", "reads.bam", "--output", "out", "--haplotag", "off"];
    let out = lrpipe_run(&dir, &args).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains("follow the reads' own HP tags (--haplotag off), not the run's phasing"),
        "{}",
        stderr(&out)
    );
    let manifest = written_manifest(&dir.join("out/manifest.json"));
    let stages = manifest["stages"].as_array().expect("stages");
    for absent in ["haplotag", "sv-phase", "repeats"] {
        assert!(stages.iter().all(|s| s["name"] != absent), "no {absent} stage: {manifest}");
    }
    let snv = stages.iter().find(|s| s["name"] == "snv").expect("snv")["command"].to_string();
    assert!(snv.contains("--enable-phasing=true"), "the SNVs are phased: {snv}");
    assert_eq!(manifest["sample"]["haplotypes"], "the input's HP tags (not lrpipe run's phasing)");
    for table in ["combined", "hap1", "hap2"] {
        assert!(dir.join(format!("out/methylation/S.cpg.{table}.bed.gz")).exists(), "{table}");
    }

    // Reads without HP tags: the remedy is to leave --haplotag on.
    let untagged = fake_site("haplotag-off-untagged", "");
    let out = lrpipe_run(&untagged, &[&args[..], &["--dry-run"]].concat()).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("for them, leave out --haplotag off"), "{}", stderr(&out));

    // The haplotagged BAM needs the haplotagging.
    let out =
        lrpipe_run(&dir, &[&args[..], &["--haplotagged-reads", "--dry-run"]].concat()).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("which --haplotag off leaves out"), "{}", stderr(&out));
    for d in [dir, untagged] {
        let _ = fs::remove_dir_all(d);
    }
}

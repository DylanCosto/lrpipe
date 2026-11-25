//! `lrpipe run` as its user meets it: what it says when it is stopped, when its inputs are
//! wrong, and when a stage fails. The tools are a shell script standing in for the lrpipe
//! binary's, on the ONT fixture's reads.

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};

/// A stand-in for the lrpipe binary's tools: lrpipe-snv's outputs (after `$FAKE_SNV_SLEEP`
/// seconds, as the stage's own process, when it is set; with `$FAKE_MARK` too, as a shell that
/// waits for a `sleep` of its own, whose ID it writes to `$FAKE_MARK.pid`, and writes `$FAKE_MARK`
/// if it is sent SIGINT; failing with `$FAKE_SNV_FAIL` as its last words when that is set), and
/// the merge's and the QC's.
const FAKE_LRPIPE: &str = r#"#!/bin/sh
if [ "$1" = "--version" ]; then echo "fake-lrpipe 0"; exit 0; fi
index="$(dirname "$0")/empty.vcf.tbi"
cmd=$1; shift
if [ "$cmd" = snv ] && [ -n "$FAKE_MARK" ]; then
  trap 'echo INT > "$FAKE_MARK"' INT
  sleep "$FAKE_SNV_SLEEP" & echo $! > "$FAKE_MARK.pid"; wait; exit 0
fi
if [ "$cmd" = snv ] && [ -n "$FAKE_SNV_SLEEP" ]; then exec sleep "$FAKE_SNV_SLEEP"; fi
if [ "$cmd" = snv ] && [ -n "$FAKE_SNV_FAIL" ]; then
  for i in 1 2 3 4 5 6 7 8 9 10 11 12; do echo "progress line $i"; done
  echo; echo "Error: $FAKE_SNV_FAIL" >&2; exit 1
fi
while [ $# -gt 0 ]; do
  case "$1" in --output) out=$2; shift;; --json) json=$2; shift;; --tsv) tsv=$2; shift;; esac
  shift
done
case $cmd in
  snv) mkdir -p "$out"; echo x > "$out/merge_output.vcf.gz"; cp "$index" "$out/merge_output.vcf.gz.tbi";;
  merge) echo x > "$out"; cp "$index" "$out.tbi";;
  qc) mkdir -p "$(dirname "$json")"; echo '{}' > "$json"; echo x > "$tsv";;
esac
exit 0
"#;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture").join(name)
}

/// A scratch site: the fixture's reads as `reads.bam`, a reference whose `.fai` has the reads'
/// chr20 (the FASTA itself is never read), a model, [`FAKE_LRPIPE`] as the tools and
/// `site.toml` naming them, with `snv` added to its `[snv]`.
fn site(name: &str, snv: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-ux-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("model")).expect("scratch directory");
    fs::copy(fixture("ont.bam"), dir.join("reads.bam")).expect("reads");
    fs::copy(fixture("ont.bam.bai"), dir.join("reads.bam.bai")).expect("index");
    let mut index = pasta::tabix::index::Indexer::default();
    index.set_header(pasta::csi::binning_index::index::header::Builder::vcf().build());
    pasta::tabix::fs::write(dir.join("empty.vcf.tbi"), &index.build()).expect("empty tabix index");
    for (file, text) in [
        ("ref.fa", ">chr20\nACGT\n"),
        ("ref.fa.fai", "chr20\t64444167\t7\t60\t61\n"),
        ("model/pileup.pt", "pileup"),
        ("model/full_alignment.pt", "full alignment"),
        ("fake-lrpipe", FAKE_LRPIPE),
    ] {
        fs::write(dir.join(file), text).expect("written");
    }
    fs::set_permissions(dir.join("fake-lrpipe"), fs::Permissions::from_mode(0o755)).expect("executable");
    let site =
        format!("reference = \"ref.fa\"\n[tools]\nlrpipe = \"fake-lrpipe\"\n[snv]\nont_model = \"model\"\n{snv}");
    fs::write(dir.join("site.toml"), site).expect("site configuration");
    dir
}

/// `lrpipe run` in `dir` on the site's sample, SNVs alone, with `args` after the sample's.
fn lrpipe_run(dir: &Path, args: &[&str]) -> Command {
    run_on(dir, "reads.bam", "S", &[&["--analyses", "snv"], args].concat())
}

/// `lrpipe run` in `dir` on `reads` as `sample`, into `out`, with `args` after them.
fn run_on(dir: &Path, reads: &str, sample: &str, args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_lrpipe"));
    command.current_dir(dir).args(["run", "--config", "site.toml", "--reads", reads, "--sample", sample]);
    command.args(["--platform", "ont", "--sex", "female", "--memory", "64", "--output", "out"]);
    command.args(args);
    command
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

#[test]
fn a_foreign_fai_does_not_replace_the_completed_manifest() {
    let dir = site("foreign-fai", "");
    let out = lrpipe_run(&dir, &[]).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    let path = dir.join("out/manifest.json");
    let completed = fs::read(&path).expect("completed manifest");

    // The offset from a FASTA with a longer defline points into a sequence here.
    fs::write(dir.join("ref.fa.fai"), "chr20\t64444167\t8\t60\t61\n").expect("write a foreign index");
    let out = lrpipe_run(&dir, &["--rerun", "all"]).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2), "{}", stderr(&out));
    assert!(stderr(&out).contains("its .fai is not this FASTA's index"), "{}", stderr(&out));
    assert_eq!(fs::read(&path).expect("completed manifest"), completed);
    fs::remove_dir_all(dir).expect("remove scratch directory");
}

fn manifest(dir: &Path) -> serde_json::Value {
    serde_json::from_str(&fs::read_to_string(dir.join("out/manifest.json")).expect("a manifest")).expect("json")
}

/// Waits until the manifest says the snv stage is running.
fn wait_for_snv(dir: &Path) {
    let started = Instant::now();
    while started.elapsed() < Duration::from_secs(30) {
        let running = fs::read_to_string(dir.join("out/manifest.json"))
            .ok()
            .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
            .is_some_and(|m| m["stages"][0]["status"] == "running");
        if running {
            // The stage's process has started; give it the moment it needs to exec.
            std::thread::sleep(Duration::from_millis(300));
            return;
        }
        std::thread::sleep(Duration::from_millis(50));
    }
    panic!("the snv stage never started");
}

/// Sends SIGsignal to `target`: a process, or with a leading `-` a process group.
fn signal(signal: &str, target: &str) {
    let status = Command::new("kill").args([&format!("-{signal}"), "--", target]).status().expect("kill runs");
    assert!(status.success(), "kill -{signal} {target}");
}

/// Starts `command` as a terminal starts a job: in a process group of its own.
fn spawn_job(command: &mut Command) -> std::process::Child {
    command.env("FAKE_SNV_SLEEP", "30").stdout(Stdio::null()).stderr(Stdio::piped()).process_group(0);
    command.spawn().expect("lrpipe runs")
}

/// A terminal's Ctrl-C sends SIGINT to the job's process group, the stages with lrpipe when they
/// share it: the run must end as interrupted, its stages stopped, never as a failed stage.
/// The stage itself never sees the terminal's SIGINT (it runs in a process group of its own),
/// and lrpipe's kill of its group reaches what it started.
#[test]
fn ctrl_c_to_the_process_group_interrupts_the_run() {
    let dir = site("ctrl-c", "");
    let mark = dir.join("stage-saw-sigint");
    let child = spawn_job(lrpipe_run(&dir, &[]).env("FAKE_MARK", &mark));
    wait_for_snv(&dir);
    let grandchild = fs::read_to_string(mark.with_extension("pid")).expect("the stage's sleep").trim().to_owned();
    signal("INT", &format!("-{}", child.id()));
    let out = child.wait_with_output().expect("lrpipe ends");
    let text = stderr(&out);
    assert_eq!(out.status.code(), Some(130), "{text}");
    assert!(!text.contains("FAILED"), "{text}");
    assert!(text.contains("interrupted; running it again resumes"), "{text}");
    let m = manifest(&dir);
    assert_eq!(m["status"], "interrupted", "{m}");
    assert_eq!(m["stages"][0]["status"], "cancelled", "{m}");
    assert!(!mark.exists(), "the terminal's SIGINT reached the stage");
    assert!(ended(&grandchild), "the stage's own child outlived it");
    let _ = fs::remove_dir_all(dir);
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

/// A closed terminal or a dropped ssh session sends SIGHUP to lrpipe's group, which no longer
/// holds the stages: lrpipe must stop them and finish the manifest, not die and leave them.
#[test]
fn a_hangup_interrupts_the_run() {
    let dir = site("hangup", "");
    let child = spawn_job(&mut lrpipe_run(&dir, &[]));
    wait_for_snv(&dir);
    signal("HUP", &format!("-{}", child.id()));
    let out = child.wait_with_output().expect("lrpipe ends");
    assert_eq!(out.status.code(), Some(130), "{}", stderr(&out));
    let m = manifest(&dir);
    assert_eq!(m["status"], "interrupted", "{m}");
    assert_eq!(m["stages"][0]["status"], "cancelled", "{m}");
    let _ = fs::remove_dir_all(dir);
}

/// Under `nohup` (SIGHUP ignored), a hangup leaves the run going (LRUN-05).
#[test]
fn under_nohup_a_hangup_changes_nothing() {
    let dir = site("nohup", "");
    let run = lrpipe_run(&dir, &[]);
    let mut command = Command::new("nohup");
    command.current_dir(&dir).arg(run.get_program()).args(run.get_args());
    let child = spawn_job(&mut command);
    wait_for_snv(&dir);
    signal("HUP", &format!("-{}", child.id()));
    std::thread::sleep(Duration::from_millis(800));
    assert_eq!(manifest(&dir)["stages"][0]["status"], "running", "a hangup stopped a nohup run");
    signal("INT", &child.id().to_string());
    let out = child.wait_with_output().expect("lrpipe ends");
    assert_eq!(out.status.code(), Some(130), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

/// Two signals in quick succession (`timeout -s INT` sends one to lrpipe and one to its group):
/// the run still finishes its manifest.
#[test]
fn a_second_signal_still_finishes_the_manifest() {
    let dir = site("two-signals", "");
    let child = spawn_job(&mut lrpipe_run(&dir, &[]));
    wait_for_snv(&dir);
    signal("INT", &child.id().to_string());
    signal("INT", &child.id().to_string());
    let out = child.wait_with_output().expect("lrpipe ends");
    assert_eq!(out.status.code(), Some(130), "{}", stderr(&out));
    let m = manifest(&dir);
    assert_eq!(m["status"], "interrupted", "{m}");
    assert_eq!(m["stages"][0]["status"], "cancelled", "{m}");
    assert!(m["finished"].is_string(), "{m}");
    let _ = fs::remove_dir_all(dir);
}

/// Runs `command`, which must fail before any stage with exit code 2 and say `wanted`; it
/// writes nothing into the output directory.
fn refused(dir: &Path, command: &mut Command, wanted: &str) {
    let out = command.output().expect("lrpipe runs");
    let text = stderr(&out);
    assert_eq!(out.status.code(), Some(2), "{text}");
    assert!(text.contains(wanted), "wanted {wanted:?} in: {text}");
    assert!(!dir.join("out/logs").exists(), "a stage ran: {text}");
}

/// LUX-02: a reference whose contigs are not the reads' is an error before any stage, which
/// gave an empty VCF and exit code 0 for one naming them `20` against reads on `chr20`. The
/// reads' header names a whole genome and the site's reference only chr20, where the reads are:
/// that one is right.
#[test]
fn a_reference_without_the_reads_contigs_is_refused() {
    let dir = site("contigs", "");
    let out = lrpipe_run(&dir, &["--dry-run"]).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    fs::write(dir.join("ref.fa"), ">20\nACGT\n").expect("reference for the unprefixed contig");
    fs::write(dir.join("ref.fa.fai"), "20\t64444167\t4\t60\t61\n").expect("written");
    refused(
        &dir,
        &mut lrpipe_run(&dir, &["--dry-run"]),
        "reads.bam has reads on: it has 20, so it names its contigs without `chr`",
    );
    refused(&dir, &mut lrpipe_run(&dir, &[]), "has no contig chr20");
    fs::write(dir.join("ref.fa"), ">chr20\nACGT\n").expect("reference for the prefixed contig");
    fs::write(dir.join("ref.fa.fai"), "chr20\t64444000\t7\t60\t61\n").expect("written");
    refused(&dir, &mut lrpipe_run(&dir, &["--dry-run"]), "contig chr20 has 64444000 bases in the reference");
    let _ = fs::remove_dir_all(dir);
}

/// LUX-09: a missing reference, or one without its `.fai`, is said to be so, with the command
/// that indexes it.
#[test]
fn a_missing_reference_or_index_is_named() {
    let dir = site("reference", "");
    fs::remove_file(dir.join("ref.fa.fai")).expect("removed");
    refused(&dir, &mut lrpipe_run(&dir, &["--dry-run"]), "has no index");
    refused(&dir, &mut lrpipe_run(&dir, &["--dry-run"]), "samtools faidx");
    fs::remove_file(dir.join("ref.fa")).expect("removed");
    refused(
        &dir,
        &mut lrpipe_run(&dir, &["--dry-run"]),
        "ref.fa (the site configuration's `reference`) does not exist",
    );
    let _ = fs::remove_dir_all(dir);
}

/// LUX-04: a sample name makes file names and fills the sample column, so it may not leave
/// the output directory or break either.
#[test]
fn unsafe_sample_names_are_refused() {
    let dir = site("sample", "");
    for bad in ["../evil", "a/b", "", "HG 003", "HG\t003", ".hidden"] {
        refused(&dir, &mut run_on(&dir, "reads.bam", bad, &["--analyses", "snv"]), "a sample name is letters, digits");
    }
    assert!(!dir.join("evil.vcf.gz").exists());
    let _ = fs::remove_dir_all(dir);
}

/// A BAM of the header `text` (its `@SQ` lines' contigs, which are the BAM's) and no reads.
fn write_bam(path: &Path, text: &str) {
    use std::io::Write;
    let contigs: Vec<(&str, i32)> = text
        .lines()
        .filter_map(|l| l.strip_prefix("@SQ\tSN:"))
        .filter_map(|l| l.split_once("\tLN:"))
        .map(|(name, len)| (name, len.parse().expect("a length")))
        .collect();
    let mut bam = pasta::bgzf::io::Writer::new(fs::File::create(path).expect("BAM"));
    bam.write_all(b"BAM\x01").expect("written");
    bam.write_all(&i32::try_from(text.len()).expect("short").to_le_bytes()).expect("written");
    bam.write_all(text.as_bytes()).expect("written");
    bam.write_all(&i32::try_from(contigs.len()).expect("few").to_le_bytes()).expect("written");
    for (name, len) in contigs {
        bam.write_all(&i32::try_from(name.len() + 1).expect("short").to_le_bytes()).expect("written");
        bam.write_all(name.as_bytes()).expect("written");
        bam.write_all(&[0]).expect("written");
        bam.write_all(&len.to_le_bytes()).expect("written");
    }
    bam.finish().expect("BAM written");
}

/// A BAM index naming one contig, which belongs to no file with the fixture's header.
const ONE_CONTIG_BAI: &[u8] = b"BAI\x01\x01\0\0\0\0\0\0\0\0\0\0\0";

/// LUX-08: the reads are judged by their contents before their index: a missing file, SAM, BAM
/// named as CRAM, and an index that is not theirs each get their own message.
#[test]
fn reads_files_are_checked_before_their_index() {
    let dir = site("reads", "");
    let run = |reads: &str| run_on(&dir, reads, "S", &["--analyses", "snv", "--dry-run"]);
    refused(&dir, &mut run("gone.bam"), "gone.bam does not exist");
    fs::write(dir.join("reads.sam"), "@HD\tVN:1.6\n").expect("written");
    refused(&dir, &mut run("reads.sam"), "SAM text, which lrpipe does not read: convert, sort and index it");
    fs::copy(dir.join("reads.bam"), dir.join("bam.cram")).expect("copied");
    fs::copy(dir.join("reads.bam.bai"), dir.join("bam.cram.bai")).expect("copied");
    refused(&dir, &mut run("bam.cram"), "is BAM, but named as CRAM");
    fs::copy(dir.join("reads.bam"), dir.join("other.bam")).expect("copied");
    refused(&dir, &mut run("other.bam"), "has no index");
    fs::write(dir.join("other.bam.bai"), ONE_CONTIG_BAI).expect("written");
    refused(&dir, &mut run("other.bam"), "other.bam (it indexes 1 contig, the reads' header has 195)");
    write_bam(&dir.join("by-name.bam"), "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr20\tLN:64444167\n");
    refused(&dir, &mut run("by-name.bam"), "is sorted by read name (@HD SO:queryname)");
    fs::write(dir.join("other.bam.bai"), "not an index").expect("written");
    refused(&dir, &mut run("other.bam"), "other.bam (reading it:");
    let _ = fs::remove_dir_all(dir);
}

/// LUX-03, LUX-13: configuration mistakes that reached a stage are errors at planning, with
/// --dry-run too: an `[snv] device` lrpipe-snv does not take, and a model that is not there.
#[test]
fn a_bad_device_or_a_missing_model_fails_the_dry_run() {
    let dir = site("device", "device = \"gpu\"\n");
    refused(&dir, &mut lrpipe_run(&dir, &["--dry-run"]), "[snv] device = \"gpu\": lrpipe snv takes \"cpu\"");
    let dir = site("model", "models = \".\"\n");
    let mut command = lrpipe_run(&dir, &["--dry-run", "--model", "r1041_bogus"]);
    refused(&dir, &mut command, "r1041_bogus (given by --model) does not exist; ");
    refused(&dir, &mut lrpipe_run(&dir, &["--dry-run", "--model", "r1041_bogus"]), "holds model");
    fs::remove_file(dir.join("model/full_alignment.pt")).expect("removed");
    refused(&dir, &mut lrpipe_run(&dir, &["--dry-run"]), "(given by snv.ont_model) has no full_alignment.pt");
    let _ = fs::remove_dir_all(dir);
}

/// LUX-14: --memory takes a unit and refuses what no run fits in; --threads 0 is an error and
/// more threads than cores a warning.
#[test]
fn memory_and_threads_are_checked() {
    let dir = site("resources", "");
    let run = |args: &[&str]| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_lrpipe"));
        command.current_dir(&dir).args(["run", "--config", "site.toml", "--reads", "reads.bam", "--sample", "S"]);
        command.args(["--platform", "ont", "--sex", "female", "--output", "out", "--analyses", "snv", "--dry-run"]);
        command.args(args).output().expect("lrpipe runs")
    };
    for (args, wanted) in [
        (&["--memory", "0"][..], "less than the 1 GiB a run needs"),
        (&["--memory", "8 parsecs"], "give GiB, as 32, 32G or 32GiB"),
        (&["--threads", "0"], "at least 1 thread"),
    ] {
        let out = run(args);
        assert_eq!(out.status.code(), Some(2), "{args:?}: {}", stderr(&out));
        assert!(stderr(&out).contains(wanted), "{args:?}: {}", stderr(&out));
    }
    assert!(run(&["--memory", "8G"]).status.success());
    let out = run(&["--threads", "100000"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("warning: --threads 100000 is more than the"), "{}", stderr(&out));
    let _ = fs::remove_dir_all(dir);
}

/// LUX-07: `--analyses all` is every analysis, as the help says the default is.
#[test]
fn analyses_all_is_the_default() {
    let dir = site("all", "");
    let stages = |args: &[&str]| {
        let out = run_on(&dir, "reads.bam", "S", &[args, &["--dry-run"]].concat()).output().expect("lrpipe runs");
        assert!(out.status.success(), "{}", stderr(&out));
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .filter(|l| l.starts_with("# "))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    let default = stages(&[]);
    assert!(default.iter().any(|l| l.starts_with("# cnv")), "{default:?}");
    assert_eq!(stages(&["--analyses", "all"]), default);
    let _ = fs::remove_dir_all(dir);
}

/// LUX-03: a failed stage's last words are on the terminal, not only in its log.
#[test]
fn a_failed_stage_shows_the_end_of_its_log() {
    let dir = site("failed", "");
    let out = lrpipe_run(&dir, &[]).env("FAKE_SNV_FAIL", "the reads are on fire").output().expect("lrpipe runs");
    let text = stderr(&out);
    assert_eq!(out.status.code(), Some(2), "{text}");
    assert!(text.contains("stage snv failed: exit code 1"), "{text}");
    assert!(text.contains("snv.log:\n  progress line 4\n"), "the last 10 lines that are not blank: {text}");
    assert!(text.ends_with("  progress line 12\n  Error: the reads are on fire\n"), "{text}");
    assert!(!text.contains("progress line 3\n"), "{text}");
    let _ = fs::remove_dir_all(dir);
}

/// LUX-15: a narrower run into a directory with another run's outputs warns of those it leaves;
/// LUX-11: a sample other than the read groups' is warned of.
#[test]
fn outputs_left_by_another_run_and_another_sample_name_are_warned_of() {
    let dir = site("leftovers", "");
    let out = lrpipe_run(&dir, &[]).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains("of an earlier run"), "{}", stderr(&out));
    assert!(stderr(&out).contains("--sample S is not the sample the reads' read groups name (SM:hg003)"));
    fs::create_dir_all(dir.join("out/sv")).expect("directory");
    fs::write(dir.join("out/sv/S.sv.vcf"), "").expect("written");
    let out = lrpipe_run(&dir, &[]).output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains("holds an output of an earlier run that this one does not write (sv/S.sv.vcf)"),
        "{}",
        stderr(&out)
    );
    assert!(dir.join("out/sv/S.sv.vcf").exists(), "left, not removed");
    let _ = fs::remove_dir_all(dir);
}

/// An output directory that is a file is said to be one, not "File exists (os error 17)".
#[test]
fn an_output_that_is_a_file_is_said_to_be_one() {
    let dir = site("output-file", "");
    fs::write(dir.join("out"), "").expect("written");
    let out = lrpipe_run(&dir, &[]).output().expect("lrpipe runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(
        stderr(&out).contains("--output") && stderr(&out).contains("is a file, not a directory"),
        "{}",
        stderr(&out)
    );
    let _ = fs::remove_dir_all(dir);
}

/// Several aligned --reads: each is checked, its contigs against the reference too, and the
/// sample name against all their read groups.
#[test]
fn several_aligned_reads_are_each_checked() {
    let dir = site("several", "");
    // The fixture's reads under other names (a copy would hold the same reads, which is refused).
    {
        use pasta::bam;
        use pasta::sam::alignment::io::Write as _;
        let mut reader = bam::io::Reader::new(fs::File::open(dir.join("reads.bam")).expect("the reads"));
        let header = reader.read_header().expect("a header");
        let mut writer = bam::io::Writer::new(fs::File::create(dir.join("second.bam")).expect("second.bam"));
        writer.write_header(&header).expect("a header");
        for record in reader.record_bufs(&header) {
            let mut record = record.expect("a record");
            *record.name_mut() = Some(format!("second-{}", record.name().expect("a name")).into());
            writer.write_alignment_record(&header, &record).expect("written");
        }
        writer.try_finish().expect("finished");
    }
    let run = || run_on(&dir, "reads.bam", "S", &["--reads", "second.bam", "--analyses", "snv", "--dry-run"]);
    refused(&dir, &mut run(), "second.bam has no index");
    let pasta::bam::Index::Bai(index) = pasta::bam::fs::index(dir.join("second.bam")).expect("indexed") else {
        panic!("a BAI")
    };
    pasta::bam::bai::fs::write(dir.join("second.bam.bai"), &index).expect("index written");
    let out = run().output().expect("lrpipe runs");
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("--sample S is not the sample the reads' read groups name (SM:hg003)"));
    fs::write(dir.join("ref.fa"), ">20\nACGT\n").expect("reference for the unprefixed contig");
    fs::write(dir.join("ref.fa.fai"), "20\t64444167\t4\t60\t61\n").expect("written");
    refused(&dir, &mut run(), "reads.bam has reads on");
    let _ = fs::remove_dir_all(dir);
}

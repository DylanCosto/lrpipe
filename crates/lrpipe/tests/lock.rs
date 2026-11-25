//! A run turned away from an output directory another run holds changes nothing there
//! (LRP-019): not the PAR BED the running one's stages read.

use std::fs;
use std::path::Path;
use std::process::Command;

use lrpipe::plan::Layout;
use lrpipe::run::lock;

#[test]
fn a_busy_output_directory_keeps_the_running_runs_pars() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("busy-output");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("model")).expect("scratch directory");
    // A GRCh37 reference (by chrX's length) for the second run; the first wrote GRCh38's PARs.
    fs::write(dir.join("ref.fa"), ">chrX\nACGT\n").expect("reference");
    // With the reads' chr20, as the reference they were aligned to has it.
    fs::write(dir.join("ref.fa.fai"), "chrX\t155270560\t6\t4\t5\nchr20\t64444167\t6\t4\t5\n").expect("reference index");
    reads(&dir);
    for name in ["model/pileup.pt", "model/full_alignment.pt"] {
        fs::write(dir.join(name), "").expect("input");
    }
    fs::write(dir.join("site.toml"), "reference = \"ref.fa\"\n[snv]\nont_model = \"model\"\n").expect("site");

    let layout = Layout::new(dir.join("out"), "S");
    let running = lock(&layout).expect("the first run holds the directory");
    let par = dir.join("out/.lrpipe/par.bed");
    let grch38 = "chrX\t10000\t2781479\nchrX\t155701382\t156030895\n";
    fs::write(&par, grch38).expect("the first run's PARs");

    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(&dir)
        .args(["run", "--config", "site.toml", "--reads", "reads.bam", "--sample", "S", "--platform", "ont"])
        .args(["--sex", "male", "--analyses", "snv", "--output", "out"])
        .output()
        .expect("running lrpipe");
    assert_eq!(out.status.code(), Some(2), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stderr).contains("another lrpipe run is using"));
    assert_eq!(fs::read_to_string(&par).expect("PARs"), grch38, "the running run's PARs are untouched");
    drop(running);
}

/// `reads.bam` and its index in `dir`: the ONT fixture's, which a run reads for its header (the
/// Clair3 model) and HP tags before it takes the output directory.
fn reads(dir: &Path) {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture/ont.bam");
    fs::copy(&fixture, dir.join("reads.bam")).expect("reads");
    fs::copy(Layout::index(&fixture, "bai"), dir.join("reads.bam.bai")).expect("index");
}

/// A scratch directory with a site whose inputs exist (but for the reads they are never read:
/// the run is turned away first), and the command line of a run of `sample` into `out`.
fn site(name: &str) -> std::path::PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("model")).expect("scratch directory");
    fs::write(dir.join("ref.fa"), ">chrX\nACGT\n").expect("reference");
    fs::write(dir.join("ref.fa.fai"), "chrX\t156040895\t6\t4\t5\nchr20\t64444167\t6\t4\t5\n").expect("reference index");
    reads(&dir);
    for name in ["model/pileup.pt", "model/full_alignment.pt"] {
        fs::write(dir.join(name), "").expect("input");
    }
    fs::write(dir.join("site.toml"), "reference = \"ref.fa\"\n[snv]\nont_model = \"model\"\n").expect("site");
    dir
}

fn run_sample(dir: &Path, sample: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(dir)
        .args(["run", "--config", "site.toml", "--reads", "reads.bam", "--sample", sample, "--platform", "ont"])
        .args(["--sex", "female", "--analyses", "snv", "--output", "out"])
        .output()
        .expect("running lrpipe")
}

/// Every file under `dir` but the lock, with its contents.
fn contents(dir: &Path) -> Vec<(std::path::PathBuf, Vec<u8>)> {
    let mut files = Vec::new();
    let mut todo = vec![dir.to_path_buf()];
    while let Some(d) = todo.pop() {
        for entry in fs::read_dir(&d).expect("listing") {
            let path = entry.expect("entry").path();
            if path.is_dir() {
                todo.push(path);
            } else if path.file_name() != Some("lock".as_ref()) {
                files.push((path.clone(), fs::read(&path).expect("file")));
            }
        }
    }
    files.sort();
    files
}

#[test]
fn another_samples_run_in_the_same_directory_is_refused_before_writing() {
    // LRUN-04: sample B run into sample A's directory replaced A's snv/, manifest and logs.
    let dir = site("other-sample");
    let out = dir.join("out");
    drop(lock(&Layout::new(&out, "A")).expect("A's run claims the directory"));
    fs::create_dir_all(out.join("snv")).expect("snv");
    fs::create_dir_all(out.join("logs")).expect("logs");
    fs::write(out.join("snv/merge_output.vcf.gz"), "calls of A").expect("A's calls");
    fs::write(out.join("logs/snv.log"), "A's log").expect("A's log");
    fs::write(out.join("manifest.json"), "{\"sample\": {\"name\": \"A\"}}").expect("A's manifest");
    let before = contents(&out);

    let run = run_sample(&dir, "B");
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("holds the run of sample A, not B"), "{stderr}");
    assert_eq!(contents(&out), before, "nothing of A's run changed");

    // A directory from before `.lrpipe/sample` is known by its manifest, which the error names
    // as the file to remove (there is no .lrpipe/sample to remove).
    fs::remove_file(out.join(".lrpipe/sample")).expect("an older run's directory");
    let before = contents(&out);
    let run = run_sample(&dir, "B");
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("holds the run of sample A, not B"), "{stderr}");
    assert!(stderr.contains("manifest.json") && !stderr.contains(".lrpipe/sample"), "{stderr}");
    assert_eq!(contents(&out), before, "nothing of A's run changed");

    // A's own directory is still A's to resume.
    drop(lock(&Layout::new(&out, "A")).expect("A resumes"));
    assert_eq!(fs::read_to_string(out.join(".lrpipe/sample")).expect("claimed"), "A\n");
}

#[test]
fn an_snv_directory_no_lrpipe_run_wrote_is_refused_and_never_marked() {
    // LRUN-09 through `lrpipe run`: OUTPUT/snv/ holds an earlier Clair3 run, and the directory
    // never held an lrpipe run. Marking it for lrpipe-snv let the snv stage delete it.
    let out = Path::new(env!("CARGO_TARGET_TMPDIR")).join("foreign-snv");
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(out.join("snv/tmp/old")).expect("snv/tmp");
    fs::write(out.join("snv/tmp/old/keep.txt"), "Clair3's chunk").expect("chunk");
    fs::write(out.join("snv/merge_output.vcf.gz"), "Clair3's calls").expect("calls");
    let snv_before = contents(&out.join("snv"));

    for attempt in ["first", "second"] {
        let err = lock(&Layout::new(&out, "HG002")).expect_err("refused");
        assert!(matches!(err, lrpipe::run::RunError::ForeignSnv(_)), "{attempt}: {err}");
        assert!(err.to_string().contains("no lrpipe run wrote"), "{attempt}: {err}");
        assert_eq!(contents(&out.join("snv")), snv_before, "{attempt}: snv/ untouched, not marked");
        assert!(!out.join(".lrpipe/sample").exists(), "{attempt}: the directory is not claimed");
    }

    // Without tmp/, the results alone are refused too.
    fs::remove_dir_all(out.join("snv/tmp")).expect("tmp");
    assert!(matches!(lock(&Layout::new(&out, "HG002")), Err(lrpipe::run::RunError::ForeignSnv(_))));

    // An snv/ with no files, only empty directories, is no one's: the run goes ahead and does
    // not mark it.
    fs::remove_file(out.join("snv/merge_output.vcf.gz")).expect("calls");
    fs::create_dir_all(out.join("snv/tmp/pileup_output")).expect("empty directories");
    drop(lock(&Layout::new(&out, "HG002")).expect("an empty snv/ is taken"));
    assert!(!out.join("snv/.lrpipe-snv").exists(), "lrpipe-snv marks what it writes itself");
}

#[test]
fn an_earlier_runs_snv_directory_is_marked_for_lrpipe_snv_and_a_killed_merge_is_cleared() {
    // A directory of 0.1.5's: snv/ and snv/tmp/ without lrpipe-snv's markers, which lrpipe-snv
    // now refuses to replace; and the hidden files of a merge killed with the run.
    let out = Path::new(env!("CARGO_TARGET_TMPDIR")).join("earlier-snv");
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(out.join("snv/tmp/pileup_output")).expect("snv/tmp");
    fs::write(out.join("snv/tmp/pileup_output/pileup_chr1.1_1.vcf"), "chunk").expect("chunk");
    fs::write(out.join("snv/merge_output.vcf.gz"), "calls").expect("calls");
    fs::write(out.join("manifest.json"), "{\"sample\": {\"name\": \"A\"}}").expect("A's manifest");
    for name in [".A.vcf.gz.tmp-7-1", ".A.vcf.gz.tmp-7-1.gz", ".A.vcf.gz.tmp-7-1.gz.tbi"] {
        fs::write(out.join(name), "killed merge").expect("staged file");
    }
    fs::write(out.join("A.vcf.gz"), "merged").expect("merged");
    fs::write(out.join(".B.vcf.gz.tmp-7-1"), "not this merge's").expect("other file");

    drop(lock(&Layout::new(&out, "A")).expect("A resumes"));
    assert!(out.join("snv/.lrpipe-snv").is_file() && out.join("snv/tmp/.lrpipe-snv").is_file(), "marked");
    assert_eq!(fs::read_to_string(out.join("snv/merge_output.vcf.gz")).expect("kept"), "calls");
    assert!(out.join("snv/tmp/pileup_output/pileup_chr1.1_1.vcf").exists(), "lrpipe-snv clears it, not lrpipe");
    assert!(!out.join(".A.vcf.gz.tmp-7-1.gz").exists(), "the killed merge's files are gone");
    assert!(out.join("A.vcf.gz").exists() && out.join(".B.vcf.gz.tmp-7-1").exists(), "nothing else is");
}

#[test]
fn files_no_lrpipe_run_wrote_where_the_run_writes_are_refused() {
    // Edge-case testing: only snv/ was guarded, and any manifest.json made the directory count
    // as an earlier lrpipe run's, so a researcher's SAMPLE.vcf.gz, manifest.json, qc/ and logs/
    // were replaced without a word.
    let out = Path::new(env!("CARGO_TARGET_TMPDIR")).join("foreign-files");
    let layout = || Layout::new(&out, "HG002");
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(out.join("qc")).expect("qc/");
    fs::write(out.join("notes.txt"), "the project's own").expect("notes");

    let refused = |what: &str| {
        let err = lock(&layout()).expect_err(what);
        assert!(matches!(err, lrpipe::run::RunError::Foreign(_)), "{what}: {err}");
        assert!(err.to_string().contains("no lrpipe run wrote"), "{what}: {err}");
        assert!(!out.join(".lrpipe/sample").exists(), "{what}: the directory is not claimed");
    };
    fs::write(out.join("manifest.json"), r#"{"tool":"other"}"#).expect("manifest");
    refused("someone else's manifest.json");
    fs::remove_file(out.join("manifest.json")).expect("manifest");
    fs::write(out.join("HG002.vcf.gz"), "earlier calls").expect("vcf");
    refused("a merged VCF");
    fs::remove_file(out.join("HG002.vcf.gz")).expect("vcf");
    fs::write(out.join("qc/HG002.qc.tsv"), "mine").expect("qc");
    refused("a file in qc/");
    fs::remove_file(out.join("qc/HG002.qc.tsv")).expect("qc");
    fs::create_dir_all(out.join("logs")).expect("logs/");
    fs::write(out.join("logs/snv.log"), "mine").expect("log");
    refused("a file in logs/");
    fs::remove_file(out.join("logs/snv.log")).expect("log");

    // Files beside the run's, and directories with nothing in them, are no one's in the way.
    drop(lock(&layout()).expect("the project's own files are left beside the run"));
    assert_eq!(fs::read_to_string(out.join("notes.txt")).expect("notes"), "the project's own");
}

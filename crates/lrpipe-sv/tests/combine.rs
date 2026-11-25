//! Multi-sample calling from more `.snf` files than the open-file limit allows (LSV-003): each
//! file is opened again for every block read, with `--combine-close-handles` and also without it
//! when the inputs would not fit under the limit, and the calls are those of a run that keeps
//! every file open.

use std::{
    fs,
    os::unix::process::CommandExt,
    path::{Path, PathBuf},
    process::{Command, Output},
};

/// More inputs than [`OPEN_FILES`] allows open at once.
const SAMPLES: usize = 300;
/// The soft and hard limit on open files of the limited runs.
const OPEN_FILES: libc::rlim_t = 128;

fn lrpipe_sv(dir: &Path, args: &[&str], open_files: Option<libc::rlim_t>) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"));
    cmd.current_dir(dir).args(args).arg("--allow-overwrite");
    if let Some(n) = open_files {
        let limit = libc::rlimit { rlim_cur: n, rlim_max: n };
        // SAFETY: between fork and exec the closure only calls setrlimit, a direct wrapper of the
        // system call that neither allocates nor locks (POSIX does not list it as async-signal-safe,
        // but glibc's and musl's are), and reads nothing but its own copy of `limit`.
        unsafe {
            cmd.pre_exec(move || {
                if libc::setrlimit(libc::RLIMIT_NOFILE, &raw const limit) == 0 {
                    Ok(())
                } else {
                    Err(std::io::Error::last_os_error())
                }
            });
        }
    }
    cmd.output().expect("running lrpipe-sv")
}

fn assert_success(what: &str, output: &Output) {
    assert!(output.status.success(), "{what} failed:\n{}", String::from_utf8_lossy(&output.stderr));
}

/// The VCF without the header lines that name the command or date.
fn records(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    let kept: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with("##command=") && !l.starts_with("##fileDate=") && !l.starts_with("##source="))
        .collect();
    kept.join("\n")
}

#[test]
fn more_snf_files_than_the_open_file_limit() {
    let dir: PathBuf = Path::new(env!("CARGO_TARGET_TMPDIR")).join("combine-many");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("many")).unwrap();
    let bam = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sample1.bam");
    let output = lrpipe_sv(&dir, &["-i", bam.to_str().unwrap(), "--snf", "one.snf", "-t", "4"], None);
    assert_success("writing the .snf", &output);
    let mut list = String::new();
    for i in 1..=SAMPLES {
        fs::copy(dir.join("one.snf"), dir.join(format!("many/s{i}.snf"))).unwrap();
        list.push_str(&format!("many/s{i}.snf\tS{i}\n"));
    }
    fs::write(dir.join("many.tsv"), list).unwrap();

    let all_open = lrpipe_sv(&dir, &["-i", "many.tsv", "-v", "open.vcf", "-t", "4"], None);
    assert_success("combining with every file open", &all_open);
    let text = records(&dir.join("open.vcf"));
    let chrom = text.lines().find(|l| l.starts_with("#CHROM")).expect("a #CHROM line");
    assert_eq!(chrom.split('\t').count(), 9 + SAMPLES, "one column per sample");
    assert!(text.lines().any(|l| !l.starts_with('#')), "combined calls");

    for (name, flag) in [("flag", Some("--combine-close-handles")), ("auto", None)] {
        let vcf = format!("{name}.vcf");
        let mut args = vec!["-i", "many.tsv", "-v", &vcf, "-t", "4"];
        args.extend(flag);
        let output = lrpipe_sv(&dir, &args, Some(OPEN_FILES));
        assert_success(&format!("combining {SAMPLES} files under ulimit -n {OPEN_FILES} ({name})"), &output);
        assert_eq!(records(&dir.join(&vcf)), text, "{name}: the calls of the run with every file open");
        // Without the flag the run saw that the files would not all fit under the limit; with it
        // the files were never kept open, so there was nothing to see.
        let stderr = String::from_utf8_lossy(&output.stderr);
        let noticed = stderr.contains(&format!("limit of {OPEN_FILES} open files"));
        assert_eq!(noticed, flag.is_none(), "{name}: {stderr}");
    }
}

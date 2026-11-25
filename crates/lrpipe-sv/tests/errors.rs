//! Bad input fails with an error that says what is wrong, not a panic.

use std::{fs, path::Path, process::Command};

fn data(name: &str) -> String {
    format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))
}

/// Runs lrpipe-sv in `dir` and returns its stderr, asserting that it failed cleanly.
fn fails(dir: &Path, args: &[&str]) -> String {
    let output =
        Command::new(env!("CARGO_BIN_EXE_lrpipe-sv")).current_dir(dir).args(args).output().expect("running lrpipe-sv");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(output.status.code(), Some(1), "lrpipe-sv {} should exit with 1:\n{stderr}", args.join(" "));
    assert!(!stderr.contains("panicked"), "lrpipe-sv {} panicked:\n{stderr}", args.join(" "));
    stderr
}

fn scratch(name: &str) -> std::path::PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("errors").join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("creating scratch directory");
    dir
}

#[test]
fn missing_input() {
    let stderr = fails(&scratch("missing_input"), &["--input", "nope.bam", "--vcf", "out.vcf"]);
    assert!(stderr.contains("nope.bam"), "{stderr}");
}

#[test]
fn cram_without_reference() {
    let stderr = fails(&scratch("cram_without_reference"), &["--input", &data("sample1.cram"), "--vcf", "out.vcf"]);
    assert!(stderr.contains("requires --reference"), "{stderr}");
}

#[test]
fn existing_output_without_overwrite() {
    let dir = scratch("existing_output");
    fs::write(dir.join("out.vcf"), "").unwrap();
    let stderr = fails(&dir, &["--input", &data("sample1.bam"), "--vcf", "out.vcf"]);
    assert!(stderr.contains("--allow-overwrite"), "{stderr}");
}

#[test]
fn truncated_bam_names_file_and_region() {
    let dir = scratch("truncated_bam");
    let bam = fs::read(data("sample1.bam")).unwrap();
    fs::write(dir.join("trunc.bam"), &bam[..bam.len() / 4]).unwrap();
    fs::copy(data("sample1.bam.bai"), dir.join("trunc.bam.bai")).unwrap();
    let stderr = fails(&dir, &["--input", "trunc.bam", "--vcf", "out.vcf"]);
    assert!(stderr.contains("trunc.bam") && stderr.contains("truncated or corrupt"), "{stderr}");
}

#[test]
fn phasing_options_with_multi_sample_calling() {
    // They phase the SVs of a sample's reads; with .snf inputs they would be ignored.
    let dir = scratch("phasing_combine");
    fs::write(dir.join("a.snf"), "").expect("written");
    fs::write(dir.join("x"), "").expect("written");
    for option in ["--phasing-evidence", "--haplotags"] {
        let stderr = fails(&dir, &["--input", "a.snf", "--vcf", "out.vcf", option, "x"]);
        assert!(stderr.contains("not multi-sample calling"), "{option}: {stderr}");
    }
}

/// A reference compressed with gzip but not bgzip cannot be read by its `.fai`: an error that
/// says so, as lrpipe-reads gives (a bgzipped one is read, `parity.rs`).
#[test]
fn gzip_reference_is_not_bgzip() {
    let dir = scratch("gzip_reference");
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(&mut gz, b">chr17\nACGT\n").unwrap();
    fs::write(dir.join("ref.fa.gz"), gz.finish().unwrap()).unwrap();
    fs::write(dir.join("ref.fa.gz.fai"), "chr17\t4\t7\t4\t5\n").unwrap();
    let stderr = fails(&dir, &["--input", &data("sample1.bam"), "--vcf", "out.vcf", "--reference", "ref.fa.gz"]);
    assert!(stderr.contains("ref.fa.gz") && stderr.contains("not bgzip"), "{stderr}");
}

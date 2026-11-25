//! A sample whose calls are haploid leaves another sample's phasing as it is (LRP-017): on
//! lrpipe-snv's ONT fixture, SAMPLE phased alone and beside a sample M haploid at every site.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture")
}

fn phase(dir: &Path, vcf: &Path, out: &str, extra: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-phase"))
        .args(["phase", "--ignore-read-groups", "--reference"])
        .arg(fixture().join("chr20_window.fa.gz"))
        .args(extra)
        .arg("-o")
        .arg(dir.join(out))
        .arg(vcf)
        .arg(fixture().join("ont.bam"))
        .output()
        .expect("running lrpipe-phase");
    assert!(output.status.success(), "lrpipe-phase failed:\n{}", String::from_utf8_lossy(&output.stderr));
    fs::read_to_string(dir.join(out)).expect("reading the phased VCF")
}

/// Each record's columns up to `last` (0-based, inclusive).
fn columns(text: &str, last: usize) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').take(last + 1).collect::<Vec<_>>().join("\t"))
        .collect()
}

#[test]
fn an_unselected_haploid_sample_leaves_the_phasing_as_it_is() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("haploid-sample");
    fs::create_dir_all(&dir).expect("scratch directory");
    let calls = fixture().join("expected/merge_output.vcf");
    let alone = phase(&dir, &calls, "alone.vcf", &[]);

    let mut two = String::new();
    for line in fs::read_to_string(&calls).expect("reading the calls").lines() {
        match line {
            l if l.starts_with("##") => two.push_str(&format!("{l}\n")),
            l if l.starts_with('#') => two.push_str(&format!("{l}\tM\n")),
            l => two.push_str(&format!("{l}\t1\n")),
        }
    }
    let two_path = dir.join("two.vcf");
    fs::write(&two_path, two).expect("writing the two-sample VCF");
    let beside = phase(&dir, &two_path, "beside.vcf", &["--sample", "SAMPLE"]);

    let phased = columns(&alone, 9).iter().filter(|r| r.rsplit('\t').next().is_some_and(|s| s.contains('|'))).count();
    assert!(phased > 400, "SAMPLE is phased alone: {phased} calls");
    assert_eq!(columns(&beside, 9), columns(&alone, 9), "SAMPLE's records are the same beside M");
    for record in columns(&beside, 10) {
        let m = record.rsplit('\t').next().expect("M's column");
        assert_eq!(m.split(':').next(), Some("1"), "M's haploid call stays: {record}");
    }
}

//! A run with nothing to call still writes every file it was asked for (LRP-020): with
//! `--enable_phasing`, the phased VCF beside the merged one, both indexed and without records.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../ci/fixture")
}

/// A BAM of the fixture's chr20 without reads, and its BAI, which then has no bins.
fn write_empty_bam(path: &Path) {
    let name = b"chr20\0";
    let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20\tLN:64444167\n";
    let mut bam = pasta::bgzf::io::Writer::new(fs::File::create(path).expect("BAM"));
    bam.write_all(b"BAM\x01").unwrap();
    bam.write_all(&(text.len() as i32).to_le_bytes()).unwrap();
    bam.write_all(text).unwrap();
    bam.write_all(&1i32.to_le_bytes()).unwrap();
    bam.write_all(&(name.len() as i32).to_le_bytes()).unwrap();
    bam.write_all(name).unwrap();
    bam.write_all(&64_444_167i32.to_le_bytes()).unwrap();
    bam.finish().expect("BAM written");
    let mut bai = b"BAI\x01".to_vec();
    for n in [1i32, 0, 0] {
        bai.extend(n.to_le_bytes());
    }
    fs::write(path.with_extension("bam.bai"), bai).expect("BAI written");
}

#[test]
fn a_run_without_reads_writes_the_phased_vcf_it_was_asked_for() {
    // The run checks its model before anything else (a broken one would otherwise be found only
    // after the earlier results were removed), so it needs a real one: SNV_MODEL, as
    // tools/check.sh sets it.
    let Some(model) = std::env::var_os("SNV_MODEL").map(PathBuf::from) else {
        eprintln!("skipped: no Clair3 model (set SNV_MODEL)");
        return;
    };
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("empty-input");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    let bam = dir.join("empty.bam");
    write_empty_bam(&bam);
    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe-snv"))
        .args(["run", "--bam"])
        .arg(&bam)
        .arg("--ref")
        .arg(fixture().join("chr20_window.fa.gz"))
        .args(["--platform", "ont", "--enable_phasing", "--model-path"])
        .arg(&model)
        .arg("--output")
        .arg(dir.join("out"))
        .output()
        .expect("running lrpipe-snv");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    for name in ["merge_output", "phased_merge_output"] {
        let vcf = dir.join("out").join(format!("{name}.vcf.gz"));
        assert!(vcf.with_extension("gz.tbi").exists(), "{name}.vcf.gz.tbi");
        let mut text = String::new();
        std::io::Read::read_to_string(
            &mut flate2::read::MultiGzDecoder::new(fs::File::open(&vcf).expect(name)),
            &mut text,
        )
        .expect("gzip");
        assert!(text.lines().any(|l| l.starts_with("#CHROM")), "{name} has a header");
        assert!(text.lines().all(|l| l.starts_with('#')), "{name} has no records");
    }
}

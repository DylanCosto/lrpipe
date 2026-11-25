//! A `.fai` that is another FASTA's index is refused before anything is read or written, by every
//! subcommand that reads the reference (as `lrpipe run` and `lrpipe cnv` already do).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture")
}

/// A tiny FASTA whose `.fai` was written for a FASTA with other line lengths.
fn foreign_fai_reference(dir: &Path) -> PathBuf {
    let fasta = dir.join("ref.fa");
    fs::write(&fasta, ">chr20\nACGTACGTAC\nGTACGTACGT\nACGTACGTAC\nGT\n").expect("writing the FASTA");
    fs::write(dir.join("ref.fa.fai"), "chr20\t32\t7\t20\t21\n").expect("writing the .fai");
    fasta
}

#[test]
fn foreign_fai_is_refused_before_any_output() {
    let dir = std::env::temp_dir().join(format!("lrpipe-phase-foreign-fai-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("creating the scratch directory");
    let reference = foreign_fai_reference(&dir);
    let vcf = fixture().join("expected/merge_output.vcf");
    let bam = fixture().join("ont.bam");
    let out = dir.join("out.vcf");
    let tsv = dir.join("tags.tsv");
    let cases: [(&str, Vec<&Path>); 3] = [
        ("phase", vec![Path::new("-o"), &out, &vcf, &bam]),
        ("haplotag", vec![Path::new("--list-only"), Path::new("--output-haplotag-list"), &tsv, &vcf, &bam]),
        ("genotype", vec![Path::new("-o"), &out, &vcf, &bam]),
    ];
    for (command, args) in cases {
        let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-phase"))
            .args([command, "--ignore-read-groups", "--reference"])
            .arg(&reference)
            .args(args)
            .output()
            .expect("running lrpipe-phase");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(!output.status.success(), "{command} accepted a foreign .fai");
        assert!(stderr.contains("its .fai is not this FASTA's index"), "{command}: {stderr}");
        assert!(!out.exists() && !tsv.exists(), "{command} wrote an output");
    }
    fs::remove_dir_all(&dir).ok();
}

//! F1b end to end with the real tools, as `lrpipe run`'s stages run them: on lrpipe-snv's ONT
//! fixture (152 reads without HP tags) and its calls phased (`lrpipe-phase/tests/ont_fixture`),
//! `lrpipe haplotag --list-only` writes WhatsHap 2.8's haplotag list, the CpG tables per
//! haplotype from it have rows and are those of the reads haplotagged into a BAM, and so are the
//! SV candidates phased from it after the read pass. The fixture's reads give no SV calls: the
//! calls are checked on lrpipe-sv's split-read fixture, whose list is of its reads' own tags.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn snv_fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture").join(name)
}

fn phase_fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-phase/tests/ont_fixture").join(name)
}

fn lrpipe(dir: &Path, args: &[&str]) {
    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe")).current_dir(dir).args(args).output().expect("lrpipe runs");
    assert!(out.status.success(), "lrpipe {} failed:\n{}", args.join(" "), String::from_utf8_lossy(&out.stderr));
}

/// A gzipped (or BGZF) file's text.
fn gunzip(path: &Path) -> String {
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display())))
        .read_to_string(&mut text)
        .expect("gzipped text");
    text
}

/// A table's or VCF's lines but its header.
fn rows(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

#[test]
fn the_runs_haplotag_list_gives_the_haplotypes_of_the_tagged_reads() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("f1b-haplotag");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    let (bam, reference, vcf) =
        (snv_fixture("ont.bam"), snv_fixture("chr20_window.fa.gz"), phase_fixture("phased.vcf.gz"));
    let path = |p: &Path| p.to_str().expect("UTF-8").to_owned();
    let (bam, reference, vcf) = (path(&bam), path(&reference), path(&vcf));
    // As the haplotag stage runs it, and with the reads written haplotagged.
    let haplotag = ["haplotag", "--reference", &reference, "--ignore-read-groups", "--untag-missing-contigs"];
    lrpipe(&dir, &[&haplotag[..], &["--output-haplotag-list", "list.tsv.gz", "--list-only", &vcf, &bam]].concat());
    lrpipe(&dir, &[&haplotag[..], &["--output", "tagged.bam", "--write-index", &vcf, &bam]].concat());
    let list = gunzip(&dir.join("list.tsv.gz"));
    let expected = fs::read_to_string(phase_fixture("haplotags.whatshap.tsv")).expect("WhatsHap's list");
    assert_eq!(list, expected, "WhatsHap 2.8's list");
    assert_eq!(rows(&list).iter().filter(|l| l.contains("\tH1\t")).count(), 74);
    assert_eq!(rows(&list).iter().filter(|l| l.contains("\tH2\t")).count(), 58);

    // The CpG tables per haplotype: from the list on the untagged reads, and from the tags.
    // (Within the reference's bases: the rest of chr20 is N, which takes the tables long to pass.)
    let cpg = ["methyl", "cpg", "--ref", &reference, "--threads", "2", "--region", "chr20:30140001-30310000", "--bam"];
    lrpipe(&dir, &[&cpg[..], &[&bam, "--output-prefix", "listed", "--haplotags", "list.tsv.gz"]].concat());
    lrpipe(&dir, &[&cpg[..], &["tagged.bam", "--output-prefix", "tagged"]].concat());
    for table in ["combined", "hap1", "hap2"] {
        let listed = gunzip(&dir.join(format!("listed.{table}.bed.gz")));
        assert!(rows(&listed).len() > 300, "{table}: {} rows", rows(&listed).len());
        assert_eq!(listed, gunzip(&dir.join(format!("tagged.{table}.bed.gz"))), "{table}");
    }
    let hap1 = gunzip(&dir.join("listed.hap1.bed.gz"));
    assert!(rows(&hap1).iter().all(|r| r.contains("\thap1\t")), "{hap1}");

    // The SVs: called unphased beside lrpipe-snv, then phased from the list, as from the tags.
    // (lrpipe sv reads no bgzipped FASTA, and a BAM needs none.)
    let sv = ["sv", "--sample-id", "S", "--input"];
    lrpipe(&dir, &[&sv[..], &[&bam, "--vcf", "u.vcf", "--snf", "u.snf", "--phasing-evidence", "e.gz"]].concat());
    let phase = ["sv", "phase", "--evidence", "e.gz", "--haplotags", "list.tsv.gz", "--input-vcf", "u.vcf"];
    lrpipe(&dir, &[&phase[..], &["--vcf", "p.vcf", "--input-snf", "u.snf", "--snf", "p.snf"]].concat());
    lrpipe(&dir, &[&sv[..], &["tagged.bam", "--vcf", "t.vcf", "--snf", "t.snf"]].concat());
    // The fixture's few reads give no SV calls, but candidates, phased in the SNF.
    let snf = |name: &str| {
        let data = fs::read(dir.join(name)).expect("an SNF");
        let end = data.iter().position(|&b| b == b'\n').expect("a header") + 1;
        let mut header: serde_json::Value = serde_json::from_slice(&data[..end]).expect("JSON");
        header["config"].as_object_mut().expect("a config").remove("command");
        (header, data[end..].to_vec())
    };
    assert!(snf("t.snf").0["snf_candidate_count"].as_u64() > Some(0));
    assert_eq!(snf("p.snf"), snf("t.snf"));
    assert_ne!(snf("u.snf").1, snf("t.snf").1, "the unphased candidates are not");

    // SV calls: lrpipe-sv's reads split into primary and supplementary alignments, tagged HP 1
    // throughout by another pipeline, phased from the list of `split.bam`'s tags (as `haplotag`
    // tags them), as `lrpipe run` phases them.
    let sv_data = |name: &str| path(&Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-sv/tests/data").join(name));
    let (foreign, list, tagged) = (sv_data("split.foreign.bam"), sv_data("split.haplotags.tsv"), sv_data("split.bam"));
    lrpipe(&dir, &[&sv[..], &[&foreign, "--vcf", "su.vcf", "--snf", "su.snf", "--phasing-evidence", "se.gz"]].concat());
    let phase = ["sv", "phase", "--evidence", "se.gz", "--haplotags", &list, "--input-vcf", "su.vcf"];
    lrpipe(&dir, &[&phase[..], &["--vcf", "sp.vcf", "--input-snf", "su.snf", "--snf", "sp.snf"]].concat());
    lrpipe(&dir, &[&sv[..], &[&tagged, "--vcf", "st.vcf", "--snf", "st.snf"]].concat());
    let records = |name: &str| fs::read_to_string(dir.join(name)).expect("a VCF");
    let phased = records("sp.vcf");
    assert_eq!(rows(&phased).len(), 5);
    assert_eq!(rows(&phased), rows(&records("st.vcf")));
    assert_ne!(rows(&records("su.vcf")), rows(&phased), "the read pass's are unphased");
    assert_eq!(snf("sp.snf"), snf("st.snf"));
    let _ = fs::remove_dir_all(dir);
}

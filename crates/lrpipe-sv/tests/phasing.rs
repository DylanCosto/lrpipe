//! SVs phased after calling (`--phasing-evidence`, then `lrpipe-sv phase`) and from a haplotag
//! list in the read pass (`--haplotags`), on the parity fixture (`tests/data`), whose reads are
//! haplotagged: from their own tags as a list, both must give what Sniffles 2.8.1 writes from the
//! tags (`expected/phase.vcf`) and what lrpipe-sv writes from them (VCF and SNF), as from any
//! other list.

mod common;

use std::collections::BTreeSet;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::write_reference;

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// A scratch directory with the fixture's BAM and full-length reference.
fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    for name in [
        "sample1.bam",
        "sample1.bam.bai",
        "split.bam",
        "split.bam.bai",
        "split.foreign.bam",
        "split.foreign.bam.bai",
        "split.haplotags.tsv",
    ] {
        fs::copy(data_dir().join(name), dir.join(name)).expect("fixture");
    }
    write_reference(&data_dir().join("reference.window.fa"), &dir.join("reference.fa"), &[]);
    dir
}

fn lrpipe_sv(dir: &Path, args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv")).current_dir(dir).args(args).output().expect("lrpipe-sv");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "lrpipe-sv {} failed:\n{stderr}", args.join(" "));
    stderr
}

/// The fixture's reads' own HP and PS tags as a haplotag list (`haplotag` writes primary
/// alignments, and the fixture has no other).
fn list_of_tags(dir: &Path) -> PathBuf {
    use pasta::sam::alignment::record::data::field::Tag;
    let mut reader = pasta::bam::io::reader::Builder.build_from_path(dir.join("sample1.bam")).expect("the BAM");
    let header = reader.read_header().expect("its header");
    let mut list = String::from("#readname\thaplotype\tphaseset\tchromosome\n");
    for record in reader.record_bufs(&header) {
        let record = record.expect("a record");
        assert!(!record.flags().is_supplementary() && !record.flags().is_secondary());
        let int = |tag: &[u8; 2]| record.data().get(&Tag::new(tag[0], tag[1])).and_then(|v| v.as_int());
        let name = String::from_utf8_lossy(record.name().expect("a name").as_ref()).into_owned();
        let assigned = match (int(b"HP"), int(b"PS")) {
            (Some(hp), Some(ps)) => format!("H{hp}\t{ps}"),
            _ => "none\tnone".into(),
        };
        list.push_str(&format!("{name}\t{assigned}\tchr17\n"));
    }
    let path = dir.join("tags.tsv");
    fs::write(&path, list).expect("the list");
    path
}

/// A VCF's lines but those naming the program, command or date.
fn comparable(path: &Path) -> Vec<String> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    text.lines()
        .filter(|l| !["##source=", "##command=", "##fileDate="].iter().any(|p| l.starts_with(p)))
        .map(str::to_owned)
        .collect()
}

fn records(path: &Path) -> Vec<String> {
    comparable(path).into_iter().filter(|l| !l.starts_with('#')).collect()
}

/// An SNF's header without the command, and its blocks.
fn snf(path: &Path) -> (serde_json::Value, Vec<u8>) {
    let data = fs::read(path).expect("an SNF");
    let end = data.iter().position(|&b| b == b'\n').expect("a header line") + 1;
    let mut header: serde_json::Value = serde_json::from_slice(&data[..end]).expect("JSON");
    header["config"].as_object_mut().expect("a config").remove("command");
    (header, data[end..].to_vec())
}

/// [`phase_later`] on `input`.
fn phase_later_on(dir: &Path, input: &str, list: &str) -> String {
    let pass = ["--input", input, "--vcf", "later.vcf", "--snf", "later.snf", "--phasing-evidence", "evidence.gz"];
    lrpipe_sv(dir, &[&pass[..], &["--allow-overwrite"]].concat());
    let phase = ["phase", "--evidence", "evidence.gz", "--haplotags", list, "--input-vcf", "later.vcf"];
    lrpipe_sv(dir, &[&phase[..], &["--vcf", "phased.vcf", "--input-snf", "later.snf", "--snf", "phased.snf"]].concat())
}

#[test]
fn split_reads_and_secondary_alignments_are_phased_as_sniffles2_phases_them() {
    // `split.bam` (tools/make_split_fixture.py): reads split at a 724 bp deletion on haplotype 2
    // into a primary and a supplementary alignment, and secondary alignments, tagged as
    // `haplotag` tags them (the primary alignments only); `split.foreign.bam`, the same with every
    // alignment tagged HP 1 by another pipeline. From `split.bam`'s list, the foreign reads give
    // what Sniffles 2.8.1 writes from `split.bam`: the split reads' leads take their primary
    // alignment's haplotype, the supplementary and secondary alignments none.
    let dir = scratch("phasing-split");
    let expected = comparable(&data_dir().join("expected/split.vcf"));
    lrpipe_sv(&dir, &["--input", "split.bam", "--vcf", "tagged.vcf", "--snf", "tagged.snf"]);
    assert_eq!(comparable(&dir.join("tagged.vcf")), expected, "lrpipe-sv on the tagged reads is Sniffles2");
    assert!(expected.iter().any(|l| l.contains("SUPPORT_SA=11") && l.contains("\t0|1:")), "{expected:?}");

    lrpipe_sv(&dir, &["--input", "split.foreign.bam", "--vcf", "foreign.vcf"]);
    assert_ne!(records(&dir.join("foreign.vcf")), records(&dir.join("tagged.vcf")), "the foreign tags differ");
    lrpipe_sv(&dir, &["--input", "split.foreign.bam", "--vcf", "listed.vcf", "--haplotags", "split.haplotags.tsv"]);
    assert_eq!(comparable(&dir.join("listed.vcf")), expected);
    phase_later_on(&dir, "split.foreign.bam", "split.haplotags.tsv");
    assert_eq!(comparable(&dir.join("phased.vcf")), expected);
    assert_eq!(snf(&dir.join("phased.snf")), snf(&dir.join("tagged.snf")));
    let _ = fs::remove_dir_all(dir);
}

/// The read pass with `--phasing-evidence`, then `lrpipe-sv phase` with `list`: the phased VCF
/// and SNF (`phased.*`), beside the read pass's (`later.*`); the phasing's stderr.
fn phase_later(dir: &Path, list: &Path, reference: bool) -> String {
    let mut pass = vec!["--input", "sample1.bam", "--vcf", "later.vcf", "--snf", "later.snf"];
    pass.extend(["--phasing-evidence", "evidence.gz", "--allow-overwrite"]);
    if reference {
        pass.extend(["--reference", "reference.fa"]);
    }
    lrpipe_sv(dir, &pass);
    let list = list.to_str().expect("UTF-8");
    lrpipe_sv(
        dir,
        &[
            "phase",
            "--evidence",
            "evidence.gz",
            "--haplotags",
            list,
            "--input-vcf",
            "later.vcf",
            "--vcf",
            "phased.vcf",
            "--input-snf",
            "later.snf",
            "--snf",
            "phased.snf",
        ],
    )
}

#[test]
fn phasing_later_from_the_reads_tags_is_sniffles2s_phasing() {
    let dir = scratch("phasing-later");
    let list = list_of_tags(&dir);
    phase_later(&dir, &list, false);
    let expected: Vec<String> = comparable(&data_dir().join("expected/phase.vcf"));
    assert_eq!(comparable(&dir.join("phased.vcf")), expected, "Sniffles2's phased VCF, header and all");
    assert_ne!(records(&dir.join("later.vcf")), records(&dir.join("phased.vcf")), "the read pass's are unphased");
    assert!(
        records(&dir.join("phased.vcf")).iter().filter(|r| r.contains("\t1|0:") || r.contains("\t0|1:")).count() > 3
    );

    // The same calls from the list in the read pass, and the SNF lrpipe-sv writes from the tags.
    lrpipe_sv(&dir, &["--input", "sample1.bam", "--vcf", "listed.vcf", "--haplotags", "tags.tsv"]);
    assert_eq!(comparable(&dir.join("listed.vcf")), expected);
    lrpipe_sv(&dir, &["--input", "sample1.bam", "--vcf", "tagged.vcf", "--snf", "tagged.snf"]);
    assert_eq!(snf(&dir.join("phased.snf")), snf(&dir.join("tagged.snf")));
    assert_ne!(snf(&dir.join("later.snf")).1, snf(&dir.join("tagged.snf")).1);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn an_empty_list_leaves_the_read_passs_calls_as_they_are() {
    let dir = scratch("phasing-empty");
    let empty = dir.join("empty.tsv");
    fs::write(&empty, "#readname\thaplotype\tphaseset\tchromosome\n").expect("a list");
    let stderr = phase_later(&dir, &empty, true);
    assert_eq!(comparable(&dir.join("phased.vcf")), comparable(&dir.join("later.vcf")));
    assert_eq!(fs::read(dir.join("phased.snf")).expect("SNF"), fs::read(dir.join("later.snf")).expect("SNF"));
    // A list that names none of the reads is warned of, in either mode.
    assert!(stderr.contains("warning: ") && stderr.contains("empty.tsv names 0 of "), "{stderr}");
    let listed = lrpipe_sv(&dir, &["--input", "sample1.bam", "--vcf", "l.vcf", "--haplotags", "empty.tsv"]);
    assert!(listed.contains("warning: empty.tsv names 0 of "), "{listed}");
    let _ = fs::remove_dir_all(dir);
}

/// The evidence's reads of `contig` and the reads of call `id`.
fn evidence_reads(path: &Path, contig: &str, id: &str) -> (Vec<String>, BTreeSet<usize>) {
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).expect("evidence"))
        .read_to_string(&mut text)
        .expect("gzipped text");
    let (mut names, mut supporting, mut current) = (Vec::new(), BTreeSet::new(), String::new());
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "C" => current = f[1].to_owned(),
            "R" if current == contig => names.push(f[1].to_owned()),
            "V" if f[1] == id => supporting = f[7].split(',').filter_map(|t| t.parse().ok()).collect(),
            _ => {}
        }
    }
    (names, supporting)
}

#[test]
fn the_phasing_rescue_follows_the_list() {
    // A MOSAIC_VAF insertion at chr17:46,207,509 whose reads are all on haplotype 1, and every
    // other read on haplotype 2: Sniffles2's `rescue_phasing` passes it as 1|0.
    let dir = scratch("phasing-rescue");
    let empty = dir.join("empty.tsv");
    fs::write(&empty, "").expect("a list");
    phase_later(&dir, &empty, true);
    let (names, supporting) = evidence_reads(&dir.join("evidence.gz"), "chr17", "INS.4S10");
    assert!(supporting.len() >= 4, "{supporting:?}");
    let mut list = String::new();
    for (k, name) in names.iter().enumerate() {
        list.push_str(&format!("{name}\t{}\t7\tchr17\n", if supporting.contains(&k) { "H1" } else { "H2" }));
    }
    let list_path = dir.join("rescue.tsv");
    fs::write(&list_path, list).expect("a list");
    phase_later(&dir, &list_path, true);
    let phased = records(&dir.join("phased.vcf"));
    let rescued: Vec<&String> = phased.iter().filter(|r| r.contains("Sniffles2.INS.4S10\t")).collect();
    assert_eq!(rescued.len(), 1, "re-passed: {phased:?}");
    assert!(rescued[0].contains("\tPASS\t") && rescued[0].ends_with("1|0:60:82:11:7"), "{}", rescued[0]);
    assert_eq!(phased.len(), records(&dir.join("later.vcf")).len() + 1);

    // As the read pass phases from the same list.
    lrpipe_sv(
        &dir,
        &["--input", "sample1.bam", "--vcf", "listed.vcf", "--snf", "listed.snf", "--reference", "reference.fa"]
            .into_iter()
            .chain(["--haplotags", "rescue.tsv", "--allow-overwrite"])
            .collect::<Vec<_>>(),
    );
    assert_eq!(records(&dir.join("phased.vcf")), records(&dir.join("listed.vcf")));
    assert_eq!(snf(&dir.join("phased.snf")), snf(&dir.join("listed.snf")));
    let _ = fs::remove_dir_all(dir);
}

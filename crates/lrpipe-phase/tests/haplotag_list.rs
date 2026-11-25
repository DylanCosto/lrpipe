//! `haplotag`'s own options for `lrpipe run`: the haplotag list alone (`--list-only`), the
//! output BAM's index (`--write-index`), and reads on contigs the VCF lacks written untagged
//! (`--untag-missing-contigs`), on lrpipe-snv's ONT fixture with its calls phased.
//!
//! `ont_fixture/`: `phased.vcf.gz`, the fixture's Clair3 calls (`expected/merge_output.vcf`)
//! phased by `lrpipe-phase phase --ignore-read-groups`; `haplotags.whatshap.tsv`, the list
//! WhatsHap 2.8's `haplotag --ignore-read-groups --reference chr20_window.fa.gz` wrote for it;
//! `no_chr20.vcf.gz`, a VCF without the reads' contig.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture").join(name)
}

fn ont(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/ont_fixture").join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

fn haplotag(vcf: &Path, args: &[&str]) -> Output {
    haplotag_reads(vcf, "ont.bam", args)
}

fn haplotag_reads(vcf: &Path, reads: &str, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_lrpipe-phase"))
        .args(["haplotag", "--ignore-read-groups", "--reference"])
        .arg(fixture("chr20_window.fa.gz"))
        .args(args)
        .arg(vcf)
        .arg(fixture(reads))
        .output()
        .expect("running lrpipe-phase")
}

fn ok(out: Output) {
    assert!(out.status.success(), "lrpipe-phase failed:\n{}", String::from_utf8_lossy(&out.stderr));
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
}

#[test]
fn the_list_alone_is_the_list_written_beside_the_alignments_and_whatshaps() {
    let dir = scratch("haplotag-list-only");
    let (only, beside, bam) = (dir.join("only.tsv"), dir.join("beside.tsv"), dir.join("tagged.bam"));
    let path = |p: &Path| p.to_str().expect("UTF-8").to_owned();
    ok(haplotag(&ont("phased.vcf.gz"), &["--list-only", "--output-haplotag-list", &path(&only)]));
    ok(haplotag(
        &ont("phased.vcf.gz"),
        &["-o", &path(&bam), "--write-index", "--output-haplotag-list", &path(&beside)],
    ));
    let expected = read(&ont("haplotags.whatshap.tsv"));
    assert_eq!(read(&only), expected);
    assert_eq!(read(&beside), expected);
    // From the same reads as CRAM, whose records the list alone decodes in part, on threads.
    for threads in ["1", "4"] {
        let cram = dir.join(format!("cram{threads}.tsv"));
        let args = ["--list-only", "--threads", threads, "--output-haplotag-list", &path(&cram)];
        ok(haplotag_reads(&ont("phased.vcf.gz"), "ont.cram", &args));
        assert_eq!(read(&cram), expected, "{threads} threads");
        fs::remove_file(cram).expect("removed");
    }
    assert_eq!(expected.lines().filter(|l| l.contains("\tH1\t") || l.contains("\tH2\t")).count(), 132);
    assert!(!dir.join("only.bam").exists() && fs::read_dir(&dir).expect("listed").count() == 4);

    // The BAM's index, as `samtools index` writes it: it finds every record.
    let index = pasta::bam::bai::fs::read(bam.with_added_extension("bai")).expect("a BAI");
    let mut reader = pasta::bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(&bam)
        .expect("the tagged BAM opens with it");
    let header = reader.read_header().expect("a header");
    let region = "chr20:30000000-30400000".parse().expect("a region");
    let mut query = reader.query(&header, &region).expect("queried");
    let (mut record, mut found) = (pasta::bam::Record::default(), 0);
    while query.read_record(&mut record).expect("a record") != 0 {
        found += 1;
    }
    assert_eq!(found, 152);

    // The list alone needs the list, and takes no alignment output.
    assert!(!haplotag(&ont("phased.vcf.gz"), &["--list-only"]).status.success());
    let out = haplotag(&ont("phased.vcf.gz"), &["--list-only", "--output-haplotag-list", &path(&only), "-o", "x.bam"]);
    assert!(!out.status.success());
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn reads_on_contigs_the_vcf_lacks_can_be_written_untagged() {
    let dir = scratch("haplotag-missing-contig");
    let path = |name: &str| dir.join(name).to_str().expect("UTF-8").to_owned();
    // As WhatsHap: an error, or with --skip-missing-contigs the reads left out.
    assert!(
        !haplotag(&ont("no_chr20.vcf.gz"), &["--list-only", "--output-haplotag-list", &path("a.tsv")]).status.success()
    );
    ok(haplotag(
        &ont("no_chr20.vcf.gz"),
        &["--skip-missing-contigs", "--list-only", "--output-haplotag-list", &path("skip.tsv")],
    ));
    assert_eq!(read(&dir.join("skip.tsv")).lines().count(), 1, "the header alone");

    ok(haplotag(
        &ont("no_chr20.vcf.gz"),
        &["--untag-missing-contigs", "-o", &path("untagged.bam"), "--output-haplotag-list", &path("untagged.tsv")],
    ));
    let list = read(&dir.join("untagged.tsv"));
    assert_eq!(list.lines().skip(1).filter(|l| l.ends_with("\tnone\tnone\tchr20")).count(), 152, "{list}");
    let mut reader = pasta::bam::io::reader::Builder.build_from_path(dir.join("untagged.bam")).expect("a BAM");
    let header = reader.read_header().expect("a header");
    let records: Vec<_> = reader.record_bufs(&header).collect::<Result<_, _>>().expect("records");
    assert_eq!(records.len(), 152, "every read is written");
    let hp = pasta::sam::alignment::record::data::field::Tag::new(b'H', b'P');
    assert!(records.iter().all(|r| r.data().get(&hp).is_none()), "and none is tagged");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn regions_from_a_bed_write_each_read_once_in_order() {
    // lrpipe run's regions: `--regions-bed`, sorted and merged, writes and lists a read that
    // reaches two of them once, where WhatsHap's `--regions` writes it again for each, after reads
    // that start later. Its haplotypes are those WhatsHap gives in the same regions (from the
    // variants there, so not all the whole contig's).
    let dir = scratch("haplotag-regions-bed");
    let path = |name: &str| dir.join(name).to_str().expect("UTF-8").to_owned();
    fs::write(dir.join("regions.bed"), "chr20\t30186000\t30190000\nchr20\t30180000\t30185000\n").expect("BED");
    let records = |bam: &str| -> Vec<(String, i64)> {
        let mut reader = pasta::bam::io::reader::Builder.build_from_path(dir.join(bam)).expect("a BAM");
        let header = reader.read_header().expect("a header");
        let records: Vec<_> = reader.record_bufs(&header).collect::<Result<_, _>>().expect("records");
        let start = |r: &pasta::sam::alignment::RecordBuf| r.alignment_start().map_or(0, |p| usize::from(p) as i64);
        records.iter().map(|r| (String::from_utf8_lossy(r.name().expect("a name")).into_owned(), start(r))).collect()
    };
    ok(haplotag(
        &ont("phased.vcf.gz"),
        &["--regions-bed", &path("regions.bed"), "-o", &path("bed.bam"), "--output-haplotag-list", &path("bed.tsv")],
    ));
    ok(haplotag(
        &ont("phased.vcf.gz"),
        &[
            "--regions",
            "chr20:30180001-30185000",
            "--regions",
            "chr20:30186001-30190000",
            "-o",
            &path("whatshap.bam"),
            "--output-haplotag-list",
            &path("whatshap.tsv"),
        ],
    ));
    let (once, again) = (records("bed.bam"), records("whatshap.bam"));
    let mut sorted = once.clone();
    sorted.sort_by_key(|(_, start)| *start);
    assert_eq!(once, sorted, "in order");
    let mut unique = once.clone();
    unique.dedup();
    assert_eq!(unique.len(), once.len(), "each read once");
    assert!(again.len() > once.len() && once.iter().all(|r| again.contains(r)), "{} and {}", again.len(), once.len());
    let mut first: Vec<&str> = Vec::new();
    let whatshap = read(&dir.join("whatshap.tsv"));
    for line in whatshap.lines() {
        if !first.contains(&line) {
            first.push(line);
        }
    }
    let list = read(&dir.join("bed.tsv"));
    assert_eq!(list.lines().count(), once.len() + 1, "{list}");
    assert_eq!(list.lines().collect::<Vec<_>>(), first, "WhatsHap's, each read once");
    // Not with --regions as well.
    assert!(
        !haplotag(&ont("phased.vcf.gz"), &["--regions-bed", &path("regions.bed"), "--regions", "chr20"])
            .status
            .success()
    );
    let _ = fs::remove_dir_all(dir);
}

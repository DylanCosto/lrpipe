//! `--depth` writes the same files as lrpipe-depth (and so mosdepth) for the same input and options,
//! and does not change the SV calls.

mod common;

use std::{
    fs,
    io::Read,
    path::{Path, PathBuf},
    process::Command,
};

const OUTPUTS: [&str; 5] = [
    "mosdepth.summary.txt",
    "mosdepth.global.dist.txt",
    "mosdepth.region.dist.txt",
    "regions.bed.gz",
    "thresholds.bed.gz",
];

fn bam() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sample1.bam")
}

fn text(path: &Path) -> Option<Vec<u8>> {
    let raw = fs::read(path).ok()?;
    if path.extension().is_some_and(|e| e == "gz") {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut out).expect("decompressing");
        Some(out)
    } else {
        Some(raw)
    }
}

/// Runs lrpipe-sv with `--depth` and lrpipe-depth with the matching options, and compares.
fn check(name: &str, sv_depth_args: &[&str], opts: lrpipe_depth::run::RunOptions) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("depth").join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let sv_prefix = dir.join("sv");
    let vcf = dir.join("calls.vcf");
    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
        .args(["-i", bam().to_str().unwrap(), "-t", "4", "-v", vcf.to_str().unwrap()])
        .args(["--depth", sv_prefix.to_str().unwrap()])
        .args(sv_depth_args)
        .output()
        .expect("running lrpipe-sv");
    assert!(output.status.success(), "lrpipe-sv failed: {}", String::from_utf8_lossy(&output.stderr));

    let db_prefix = dir.join("lrpipe-depth");
    lrpipe_depth::run::run(&bam(), None, db_prefix.to_str().unwrap(), &opts).expect("running lrpipe-depth");
    for o in OUTPUTS {
        let (sv, db) = (text(&dir.join(format!("sv.{o}"))), text(&dir.join(format!("lrpipe-depth.{o}"))));
        assert_eq!(sv.is_some(), db.is_some(), "{name}: {o} written by only one of them");
        assert!(sv == db, "{name}: {o} differs from lrpipe-depth's");
    }

    // The same calls as without --depth.
    let plain = dir.join("plain.vcf");
    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
        .args(["-i", bam().to_str().unwrap(), "-t", "4", "-v", plain.to_str().unwrap()])
        .output()
        .expect("running lrpipe-sv");
    assert!(output.status.success());
    let records = |p: &Path| {
        fs::read_to_string(p).unwrap().lines().filter(|l| !l.starts_with("##")).collect::<Vec<_>>().join("\n")
    };
    assert_eq!(records(&vcf), records(&plain), "{name}: --depth changed the SV calls");
}

#[test]
fn windows_in_cigar_mode() {
    let by = lrpipe_depth::By::parse("500").unwrap();
    let opts = lrpipe_depth::run::RunOptions {
        report: lrpipe_depth::ReportOptions { by: Some(by), thresholds: vec![1, 10, 30], ..Default::default() },
        threads: 2,
        ..Default::default()
    };
    check("cigar", &["--depth-by", "500", "--depth-thresholds", "30,1,10"], opts);
}

#[test]
fn windows_in_fast_mode_with_mapq() {
    let opts = lrpipe_depth::run::RunOptions {
        mode: lrpipe_depth::Mode::Fast,
        filter: lrpipe_depth::Filter { min_mapq: 20, ..Default::default() },
        report: lrpipe_depth::ReportOptions {
            by: Some(lrpipe_depth::By::parse("1000").unwrap()),
            ..Default::default()
        },
        threads: 2,
        ..Default::default()
    };
    check("fast", &["--depth-by", "1000", "--depth-fast", "--depth-mapq", "20"], opts);
}

#[test]
fn bed_regions_with_medians() {
    let bed = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/regions.bed");
    let opts = lrpipe_depth::run::RunOptions {
        filter: lrpipe_depth::Filter { exclude_flags: 0, ..Default::default() },
        report: lrpipe_depth::ReportOptions {
            by: Some(lrpipe_depth::By::parse(bed.to_str().unwrap()).unwrap()),
            use_median: true,
            ..Default::default()
        },
        threads: 2,
        ..Default::default()
    };
    check("bed", &["--depth-by", bed.to_str().unwrap(), "--depth-median", "--depth-flag", "0"], opts);
}

#[test]
fn summary_only() {
    check("summary", &[], lrpipe_depth::run::RunOptions { threads: 1, ..Default::default() });
}

/// `--contig` restricts `--depth` and so `--cnv` to its contigs: with another reference contig
/// over 1 Mbp, CNV calling failed after the whole pass for that contig's missing coverage
/// (LSV-002). A `--contig` set with no contig lrpipe-cnv calls fails before the pass.
#[test]
fn copy_numbers_of_the_contig_option() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("depth").join("cnv-contig");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let reference = dir.join("reference.fa");
    let window = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reference.window.fa");
    common::write_reference(&window, &reference, &[("extra", 1_500_000), ("short", 50_000)]);
    let run = |name: &str, contig: &str| {
        let vcf = dir.join(format!("{name}.vcf"));
        let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
            .args(["-i", bam().to_str().unwrap(), "-t", "4", "-v", vcf.to_str().unwrap()])
            .args(["--reference", reference.to_str().unwrap(), "--contig", contig])
            .args(["--depth", dir.join(name).to_str().unwrap()])
            .args(["--depth-by", "1000", "--depth-fast", "--depth-mapq", "20"])
            .args(["--cnv", dir.join(format!("{name}.cnv")).to_str().unwrap()])
            .output()
            .expect("running lrpipe-sv");
        (output, vcf)
    };

    let (output, _) = run("chr17", "chr17");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "lrpipe-sv --contig chr17 --cnv failed: {stderr}");
    let calls = text(&dir.join("chr17.cnv/SAMPLE_cnv.bed.gz")).expect("the CNV calls");
    let calls = String::from_utf8(calls).unwrap();
    assert!(calls.lines().all(|l| l.starts_with('#') || l.starts_with("chr17\t")), "calls on chr17 only: {calls}");

    let (output, vcf) = run("short", "short");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!output.status.success(), "no contig to call");
    assert!(stderr.contains("no --contig is one (short)"), "{stderr}");
    assert!(!vcf.exists(), "the error comes before the pass");
}

/// `--depth-in-regions` (`lrpipe run --regions`): counted from the reads of `--regions` alone,
/// in the pass, the `--depth-by` regions inside them get lrpipe-depth's values from the whole
/// contig, a read reaching two regions counted once; the summary keeps the regions' rows alone,
/// and no global distribution is written.
#[test]
fn depth_in_regions_counts_their_reads_once() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("depth").join("in-regions");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    // Two regions 1 kb apart, which reads of several kb both reach, and their 1 kb windows.
    let regions = dir.join("regions.bed");
    fs::write(&regions, "chr17\t46230000\t46240000\nchr17\t46241000\t46250500\n").unwrap();
    let mut windows = String::new();
    for (start, end) in [(46_230_000, 46_240_000), (46_241_000, 46_250_500)] {
        for s in (start..end).step_by(1000) {
            windows += &format!("chr17\t{s}\t{}\n", (s + 1000).min(end));
        }
    }
    let by = dir.join("windows.bed");
    fs::write(&by, windows).unwrap();
    let path = |p: &Path| p.to_str().unwrap().to_owned();
    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
        .args(["-i", &path(&bam()), "-t", "4", "--regions", &path(&regions), "--depth", &path(&dir.join("sv"))])
        .args(["--depth-by", &path(&by), "--depth-fast", "--depth-mapq", "20", "--depth-in-regions"])
        .output()
        .expect("running lrpipe-sv");
    assert!(output.status.success(), "lrpipe-sv failed: {}", String::from_utf8_lossy(&output.stderr));
    let opts = lrpipe_depth::run::RunOptions {
        mode: lrpipe_depth::Mode::Fast,
        filter: lrpipe_depth::Filter { min_mapq: 20, ..Default::default() },
        report: lrpipe_depth::ReportOptions {
            by: Some(lrpipe_depth::By::parse(&path(&by)).unwrap()),
            ..Default::default()
        },
        threads: 2,
        ..Default::default()
    };
    let whole = dir.join("whole");
    lrpipe_depth::run::run(&bam(), None, &path(&whole), &opts).expect("running lrpipe-depth");
    let file = |prefix: &str, o: &str| text(&dir.join(format!("{prefix}.{o}")));
    for o in ["regions.bed.gz", "mosdepth.region.dist.txt"] {
        assert!(file("sv", o).is_some() && file("sv", o) == file("whole", o), "{o} differs from lrpipe-depth's");
    }
    let rows = |prefix: &str, regions_only: bool| -> Vec<String> {
        String::from_utf8(file(prefix, "mosdepth.summary.txt").unwrap())
            .unwrap()
            .lines()
            .filter(|l| !regions_only || l.starts_with("chrom") || l.split('\t').next().unwrap().ends_with("_region"))
            .map(str::to_owned)
            .collect()
    };
    let summary = rows("sv", false);
    assert_eq!(summary.len(), 3, "{summary:?}");
    assert_eq!(summary, rows("whole", true));
    assert!(file("sv", "mosdepth.global.dist.txt").is_none(), "no global distribution");

    // The regions must be sorted and apart, and hold the depth's.
    for (regions, wanted) in [
        ("chr17\t46241000\t46250500\nchr17\t46230000\t46242000\n", "sorted and apart"),
        ("chr17\t46230000\t46240000\n", "is not inside one of --regions"),
    ] {
        fs::write(dir.join("bad.bed"), regions).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
            .args(["-i", &path(&bam()), "--regions", &path(&dir.join("bad.bed")), "--depth", &path(&dir.join("bad"))])
            .args(["--depth-by", &path(&by), "--depth-in-regions"])
            .output()
            .expect("running lrpipe-sv");
        let said = String::from_utf8_lossy(&output.stderr);
        assert!(!output.status.success() && said.contains(wanted), "{said}");
    }
}

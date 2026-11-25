//! `--read-stats` counts every record once, as reading the whole file does, from BAM and CRAM
//! alike, and does not change the SV calls.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("read-stats").join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("creating the scratch directory");
    dir
}

/// Runs lrpipe-sv on `input` with `args`, and returns its `--read-stats` JSON.
fn read_stats(dir: &Path, input: &Path, args: &[&str]) -> serde_json::Value {
    let json = dir.join("stats.json");
    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
        .current_dir(dir)
        .arg("--input")
        .arg(input)
        .args(["--threads", "4", "--allow-overwrite", "--read-stats"])
        .arg(&json)
        .args(args)
        .output()
        .expect("running lrpipe-sv");
    assert!(output.status.success(), "lrpipe-sv failed: {}", String::from_utf8_lossy(&output.stderr));
    serde_json::from_str(&fs::read_to_string(&json).expect("the statistics")).expect("JSON")
}

/// The statistics of the BAM at `path`, from every one of its records in turn.
fn expected(path: &Path) -> serde_json::Value {
    expected_in(path, None)
}

/// [`expected`] of the records that reach `regions` (contig, 0-based start, end), if given.
fn expected_in(path: &Path, regions: Option<&[(&str, i64, i64)]>) -> serde_json::Value {
    use pasta::sam::alignment::Record as _;
    let mut reader = pasta::bam::io::Reader::new(fs::File::open(path).expect("the BAM"));
    let header = reader.read_header().expect("a header");
    let mut record = pasta::bam::Record::default();
    let (mut lengths, mut secondary, mut supplementary) = (Vec::new(), 0u64, 0u64);
    // The bases with qualities (the reads of a BAM have them, or none of their bases has one)
    // and their sum.
    let (mut quality_bases, mut quality) = (0u64, 0u64);
    while reader.read_record(&mut record).expect("a record") != 0 {
        let flags = record.flags();
        if let Some(regions) = regions {
            let tid = record.reference_sequence_id().transpose().expect("a reference");
            let contig = tid.map(|t| header.reference_sequences().get_index(t).expect("a contig").0.to_string());
            let start = record.alignment_start().transpose().expect("a start").map_or(0, |p| usize::from(p) as i64 - 1);
            let end = start + record.alignment_span().transpose().expect("a span").unwrap_or(0) as i64;
            let reaches = |&&(c, s, e): &&(&str, i64, i64)| contig.as_deref() == Some(c) && start < e && end > s;
            if !regions.iter().any(|r| reaches(&r)) {
                continue;
            }
        }
        if flags.is_unmapped() {
            continue;
        } else if flags.is_secondary() {
            secondary += 1;
        } else if flags.is_supplementary() {
            supplementary += 1;
        } else {
            lengths.push(record.sequence().len() as u64);
            let qualities = record.quality_scores();
            if !qualities.is_empty() && qualities.as_ref()[0] != 0xff {
                quality_bases += qualities.len() as u64;
                quality += qualities.as_ref().iter().map(|&q| u64::from(q)).sum::<u64>();
            }
        }
    }
    let bases: u64 = lengths.iter().sum();
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let mut held = 0;
    let n50 = lengths.iter().find(|&&l| {
        held += l;
        2 * held >= bases
    });
    let round = |x: f64| (x * 1000.0).round() / 1000.0;
    serde_json::json!({
        "reads": lengths.len(),
        "bases": bases,
        "n50": n50,
        "mean_length": round(bases as f64 / lengths.len() as f64),
        "mean_quality": (quality_bases > 0).then(|| round(quality as f64 / quality_bases as f64)),
        "quality_bases": quality_bases,
        "secondary": secondary,
        "supplementary": supplementary,
        "unmapped_placed": 0,
        "restricted": regions.is_some(),
    })
}

#[test]
fn every_record_is_counted_once_and_the_calls_stay_the_same() {
    // HG002 ONT reads with secondary and supplementary alignments, read in 2 Mb chunks.
    let dir = scratch("split");
    let bam = data("split.bam");
    let stats = read_stats(&dir, &bam, &["--vcf", "with.vcf"]);
    assert_eq!(stats, expected(&bam));
    assert!(stats["secondary"].as_u64() > Some(0) && stats["supplementary"].as_u64() > Some(0), "{stats}");

    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
        .current_dir(&dir)
        .arg("--input")
        .arg(&bam)
        .args(["--threads", "4", "--vcf", "without.vcf"])
        .output()
        .expect("running lrpipe-sv");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let records = |name: &str| {
        let text = fs::read_to_string(dir.join(name)).expect("a VCF");
        text.lines().filter(|l| !l.starts_with("##")).map(str::to_owned).collect::<Vec<_>>()
    };
    assert!(records("with.vcf").len() > 1, "SV calls");
    assert_eq!(records("with.vcf"), records("without.vcf"), "--read-stats changed the SV calls");
}

#[test]
fn a_cram_gives_its_bams_statistics_but_the_mean_quality() {
    // lrpipe-snv's fixture, whose reads have qualities, as BAM and CRAM.
    let dir = scratch("cram");
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let bam = read_stats(&dir, &fixture.join("ont.bam"), &[]);
    assert_eq!(bam, expected(&fixture.join("ont.bam")));
    assert!(bam["mean_quality"].as_f64().is_some_and(|q| q > 10.0), "{bam}");
    // But for the qualities, which the pass leaves undecoded.
    let reference = fixture.join("chr20_window.fa.gz");
    let mut cram = read_stats(&dir, &fixture.join("ont.cram"), &["--reference", reference.to_str().expect("UTF-8")]);
    assert_eq!((&cram["mean_quality"], &cram["quality_bases"]), (&serde_json::Value::Null, &serde_json::json!(0)));
    for key in ["mean_quality", "quality_bases"] {
        cram[key] = bam[key].clone();
    }
    assert_eq!(cram, bam);
    // Reads without qualities have no mean quality.
    let none = read_stats(&dir, &data("sample1.bam"), &[]);
    assert_eq!(none, expected(&data("sample1.bam")));
    assert!(none["mean_quality"].is_null(), "{none}");
}

/// A BAM of 40 contigs of 3,000 bases with three reads each, with qualities, and the same reads
/// as CRAM (whose small contigs share slices), against `reference.fa` of `dir`.
fn short_contigs(dir: &Path) -> (PathBuf, PathBuf) {
    use pasta::core::Position;
    use pasta::sam::alignment::io::Write as _;
    use pasta::sam::alignment::record::cigar::{Op, op::Kind};
    use pasta::sam::alignment::record::{Flags, MappingQuality};
    use pasta::sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};
    use pasta::sam::header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    };

    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build().expect("@HD");
    let (mut fasta, mut fai) = (String::new(), String::new());
    let mut builder = pasta::sam::Header::builder().set_header(hd);
    for i in 0..40 {
        let name = format!("c{i}");
        let offset = fasta.len() + name.len() + 2;
        fasta += &format!(">{name}\n{}\n", "ACGGT".repeat(600));
        fai += &format!("{name}\t3000\t{offset}\t3000\t3001\n");
        let length = std::num::NonZeroUsize::new(3000).expect("non-empty");
        builder = builder.add_reference_sequence(name, Map::<ReferenceSequence>::new(length));
    }
    fs::write(dir.join("reference.fa"), fasta).expect("written");
    fs::write(dir.join("reference.fa.fai"), fai).expect("written");
    let header = builder.build();
    let reference = lrpipe_reads::Reference::open(&dir.join("reference.fa")).expect("the reference");
    let (bam, cram) = (dir.join("short.bam"), dir.join("short.cram"));
    let mut bam_writer = pasta::bam::io::Writer::new(fs::File::create(&bam).expect("BAM"));
    bam_writer.write_header(&header).expect("header");
    let mut cram_writer = pasta::cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference.repository().clone())
        .build_from_path(&cram)
        .expect("CRAM");
    cram_writer.write_alignment_header(&header).expect("header");
    for tid in 0..40 {
        for k in 0..3 {
            let len = 50 + 30 * k + tid;
            let flags = if k == 2 && tid % 5 == 0 { Flags::SUPPLEMENTARY } else { Flags::empty() };
            let record = RecordBuf::builder()
                .set_name(format!("r{tid}-{k}"))
                .set_flags(flags)
                .set_reference_sequence_id(tid)
                .set_alignment_start(Position::new(1 + 500 * k).expect("non-zero"))
                .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
                .set_cigar([Op::new(Kind::SoftClip, 5), Op::new(Kind::Match, len - 5)].into_iter().collect::<Cigar>())
                .set_sequence(Sequence::from(b"ACGGT".repeat(len.div_ceil(5))[..len].to_vec()))
                .set_quality_scores(QualityScores::from(vec![20; len]))
                .build();
            bam_writer.write_alignment_record(&header, &record).expect("written");
            cram_writer.write_alignment_record(&header, &record).expect("written");
        }
    }
    bam_writer.try_finish().expect("finished");
    cram_writer.finish(&header).expect("finished");
    drop((bam_writer, cram_writer));
    let pasta::bam::Index::Bai(index) = pasta::bam::fs::index(&bam).expect("indexed") else { panic!("a BAI") };
    pasta::bam::bai::fs::write(dir.join("short.bam.bai"), &index).expect("index written");
    let index = pasta::cram::fs::index(&cram).expect("indexed");
    pasta::cram::crai::fs::write(dir.join("short.cram.crai"), &index).expect("index written");
    (bam, cram)
}

#[test]
fn short_contigs_of_a_cram_are_counted_without_their_reference() {
    // Contigs under 1 Mb, which the pass reads for the statistics alone: from CRAM their reads'
    // lengths come from their CIGARs, so a reference without them (as here, after writing) does
    // not stop the counting.
    let dir = scratch("short");
    let (bam, cram) = short_contigs(&dir);
    let from_bam = read_stats(&dir, &bam, &[]);
    assert_eq!(from_bam, expected(&bam));
    assert_eq!((&from_bam["reads"], &from_bam["supplementary"]), (&serde_json::json!(112), &serde_json::json!(8)));
    fs::write(dir.join("other.fa"), ">other\nACGT\n").expect("written");
    fs::write(dir.join("other.fa.fai"), "other\t4\t7\t4\t5\n").expect("written");
    let mut from_cram = read_stats(&dir, &cram, &["--reference", "other.fa"]);
    assert!(from_cram["mean_quality"].is_null(), "{from_cram}");
    for key in ["mean_quality", "quality_bases"] {
        from_cram[key] = from_bam[key].clone();
    }
    assert_eq!(from_cram, from_bam);
}

#[test]
fn in_regions_each_read_reaching_them_is_counted_once() {
    // The reads of --regions (sorted and apart, as `lrpipe run` gives them), as the depth in
    // them counts them: those starting before a region that reach it too, each once.
    let dir = scratch("regions");
    let regions = [("chr17", 46_250_000, 46_300_000), ("chr17", 46_300_100, 46_310_000)];
    let bed: String = regions.iter().map(|(c, s, e)| format!("{c}\t{s}\t{e}\n")).collect();
    fs::write(dir.join("regions.bed"), bed).expect("written");
    let bam = data("split.bam");
    let stats = read_stats(&dir, &bam, &["--regions", "regions.bed"]);
    assert_eq!(stats, expected_in(&bam, Some(&regions)));
    let whole = expected(&bam);
    assert!(stats["reads"].as_u64() < whole["reads"].as_u64() && stats["reads"].as_u64() > Some(10), "{stats}");
    assert_eq!(read_stats(&dir, &bam, &["--contig", "chr17"])["restricted"], true);
}

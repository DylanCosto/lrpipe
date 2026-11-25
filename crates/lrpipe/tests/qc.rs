//! F9: the QC summary of real runs, each value against one computed here from the run's inputs
//! and outputs: a hand-built BAM of a male sample run as female, and the ONT fixture.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-qc-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir.canonicalize().expect("canonical temp dir")
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&fs::read_to_string(path).expect("written")).expect("JSON")
}

/// A file's text, gunzipped when it is gzipped.
fn text(path: &Path) -> String {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut text).expect("gzip");
        text
    } else {
        String::from_utf8(bytes).expect("UTF-8")
    }
}

/// `lrpipe run` of sample S in `dir`, with `args`.
fn lrpipe_run(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(dir)
        .args(["run", "--sample", "S", "--memory", "16", "--threads", "4"])
        .args(args)
        .output()
        .expect("lrpipe runs")
}

/// What each primary alignment of a BAM gives: its length and its qualities' sum, if it has
/// qualities; and the secondary and supplementary alignments, and the unmapped reads.
#[derive(Debug, Default)]
struct Records {
    reads: Vec<(u64, Option<u64>)>,
    secondary: u64,
    supplementary: u64,
    unmapped: u64,
}

fn records(bam: &Path) -> Records {
    let mut reader = pasta::bam::io::Reader::new(fs::File::open(bam).expect("the BAM"));
    reader.read_header().expect("a header");
    let mut record = pasta::bam::Record::default();
    let mut out = Records::default();
    while reader.read_record(&mut record).expect("a record") != 0 {
        let flags = record.flags();
        if flags.is_unmapped() {
            out.unmapped += 1;
        } else if flags.is_secondary() {
            out.secondary += 1;
        } else if flags.is_supplementary() {
            out.supplementary += 1;
        } else {
            let q = record.quality_scores();
            let sum = (!q.is_empty()).then(|| q.as_ref().iter().map(|&b| u64::from(b)).sum());
            out.reads.push((record.sequence().len() as u64, sum));
        }
    }
    out
}

/// The length N50 of `lengths`.
fn n50(lengths: &[u64]) -> u64 {
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = sorted.iter().sum();
    let mut held = 0;
    *sorted
        .iter()
        .find(|&&l| {
            held += l;
            2 * held >= total
        })
        .expect("reads")
}

/// The QC summary's reads against those of `bam`.
fn check_reads(qc: &serde_json::Value, bam: &Path, unplaced: u64) {
    let r = records(bam);
    let lengths: Vec<u64> = r.reads.iter().map(|&(l, _)| l).collect();
    let bases: u64 = lengths.iter().sum();
    let quality: u64 = r.reads.iter().filter_map(|&(_, q)| q).sum();
    let reads = &qc["reads"];
    assert_eq!(reads["mapped_reads"], lengths.len() as u64, "{reads}");
    assert_eq!(reads["mapped_bases"], bases, "{reads}");
    assert_eq!(reads["n50"], n50(&lengths), "{reads}");
    let close = |v: &serde_json::Value, want: f64| (v.as_f64().expect("a number") - want).abs() < 1e-3;
    assert!(close(&reads["mean_length"], bases as f64 / lengths.len() as f64), "{reads}");
    assert!(close(&reads["mean_quality"], quality as f64 / bases as f64), "{reads}");
    assert_eq!(reads["unmapped"], r.unmapped, "{reads}");
    assert_eq!(r.unmapped, unplaced, "the BAM's unmapped reads");
    let all = lengths.len() as f64 + r.unmapped as f64;
    assert!(close(&reads["mapped_percent"], 100.0 * lengths.len() as f64 / all), "{reads}");
    assert_eq!((&reads["secondary"], &reads["supplementary"]), (&r.secondary.into(), &r.supplementary.into()));
    let alignments = (lengths.len() as u64 + r.secondary + r.supplementary) as f64;
    assert!(close(&reads["supplementary_percent"], 100.0 * r.supplementary as f64 / alignments), "{reads}");
}

/// The contigs of the hand-built sample, each [`LEN`] bases.
fn contigs() -> Vec<String> {
    (1..=22).map(|i| format!("chr{i}")).chain(["chrX".into(), "chrY".into()]).collect()
}

const LEN: usize = 2000;

/// A BAM of a male sample (reads on chrX and chrY at half the autosomes' density) with its BAI:
/// reads of 100 to 190 bases with qualities 10 to 19, a secondary and a supplementary alignment,
/// and two unplaced unmapped reads; and the reference.
fn male_bam(dir: &Path) -> PathBuf {
    use pasta::bam;
    use pasta::core::Position;
    use pasta::sam::alignment::io::Write as _;
    use pasta::sam::alignment::record::cigar::{Op, op::Kind};
    use pasta::sam::alignment::record::{Flags, MappingQuality};
    use pasta::sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};
    use pasta::sam::header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    };

    let mut fasta = String::new();
    let mut fai = String::new();
    let mut builder = pasta::sam::Header::builder()
        .set_header(Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build().unwrap());
    for name in contigs() {
        let offset = fasta.len() + name.len() + 2;
        fasta += &format!(">{name}\n{}\n", "ACGT".repeat(LEN / 4));
        fai += &format!("{name}\t{LEN}\t{offset}\t{LEN}\t{}\n", LEN + 1);
        let length = std::num::NonZeroUsize::new(LEN).expect("non-empty");
        builder = builder.add_reference_sequence(name, Map::<ReferenceSequence>::new(length));
    }
    fs::write(dir.join("ref.fa"), fasta).expect("written");
    fs::write(dir.join("ref.fa.fai"), fai).expect("written");
    let header = builder.build();
    let path = dir.join("male.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("BAM"));
    writer.write_header(&header).expect("header");
    let record = |name: String, contig: Option<(usize, usize)>, flags: Flags, len: usize, quality: u8| {
        let mut r = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_sequence(Sequence::from(b"ACGT".repeat(len.div_ceil(4))[..len].to_vec()))
            .set_quality_scores(QualityScores::from(vec![quality; len]));
        if let Some((tid, start)) = contig {
            r = r
                .set_reference_sequence_id(tid)
                .set_alignment_start(Position::new(start).expect("non-zero"))
                .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
                .set_cigar([Op::new(Kind::Match, len)].into_iter().collect::<Cigar>());
        }
        r.build()
    };
    for (tid, name) in contigs().iter().enumerate() {
        let n = if name.starts_with("chrX") || name.starts_with("chrY") { 5 } else { 10 };
        for k in 0..n {
            let r = record(format!("{name}-{k}"), Some((tid, 1 + 150 * k)), Flags::empty(), 100 + 10 * k, 10 + k as u8);
            writer.write_alignment_record(&header, &r).expect("record");
            // Another alignment of chr1's and chr2's last reads.
            if k == 9 && tid < 2 {
                let flags = if tid == 0 { Flags::SECONDARY } else { Flags::SUPPLEMENTARY };
                let other = record(format!("{name}-{k}"), Some((tid, 1500)), flags, 50, 30);
                writer.write_alignment_record(&header, &other).expect("record");
            }
        }
    }
    for i in 0..2 {
        writer.write_alignment_record(&header, &record(format!("u{i}"), None, Flags::UNMAPPED, 80, 5)).expect("record");
    }
    writer.try_finish().expect("finished");
    drop(writer);
    let bam::Index::Bai(index) = bam::fs::index(&path).expect("indexed") else { panic!("a BAI") };
    bam::bai::fs::write(dir.join("male.bam.bai"), &index).expect("index written");
    path
}

#[test]
fn a_wrong_sex_is_warned_of_and_the_reads_and_coverage_are_counted() {
    let dir = temp_dir("sex");
    let bam = male_bam(&dir);
    fs::write(dir.join("site.toml"), "reference = \"ref.fa\"\n[snv]\n").expect("site configuration");
    let args = ["--platform", "hifi", "--config", "site.toml", "--reads", "male.bam", "--output", "out"];
    let out = lrpipe_run(&dir, &[&args[..], &["--sex", "female", "--analyses", "depth"]].concat());
    assert!(out.status.success(), "{}", stderr(&out));
    let err = stderr(&out);
    let warning = "warning: the reads' index says the sample is male (chrX at 0.500 and chrY at 0.500 of the \
                   autosomes' reads per base), but it is called female: if the index is right, run again with --sex \
                   male (or auto)";
    // As the run starts, and again at its end.
    assert_eq!(err.matches(warning).count(), 2, "{err}");
    let qc = json(&dir.join("out/qc/S.qc.json"));
    let sex = &qc["sex"];
    assert_eq!((&sex["called"], &sex["index"], &sex["agrees"]), (&"female".into(), &"male".into(), &false.into()));
    assert_eq!(qc["warnings"].as_array().map(Vec::len), Some(1), "{qc}");
    let manifest = json(&dir.join("out/manifest.json"));
    assert_eq!(manifest["warnings"].as_array().map(Vec::len), Some(1), "the warning once: {manifest}");
    assert!(manifest["warnings"][0].as_str().is_some_and(|w| warning.ends_with(w)), "{manifest}");

    check_reads(&qc, &bam, 2);
    let line = "QC: 230 mapped reads (33.1 kb), read N50 150 b, mean length 143 b, mean base quality 15.0 (Phred, per \
                base), 99.14% mapped, of the alignments 0.43% supplementary, 0.43% secondary\n";
    assert!(err.contains(line), "{err}");
    // Each read's span, MAPQ 20 or more, the supplementary alignment too (mosdepth's -F 1796).
    let coverage = &qc["coverage"];
    // An autosome's 1450 bases of reads; chr2's supplementary alignment 50 more. Each contig's
    // mean is mosdepth's, to two decimals; the autosomes' is their bases over their length.
    let near = |v: &serde_json::Value, want: f64| v.as_f64().is_some_and(|x| (x - want).abs() < 0.0051);
    let autosome = 1450.0 / LEN as f64;
    assert!(near(&coverage["per_autosome"]["chr1"], autosome), "{coverage}");
    assert!(near(&coverage["per_autosome"]["chr2"], autosome + 50.0 / LEN as f64), "{coverage}");
    assert_eq!(coverage["per_autosome"].as_object().map(|m| m.len()), Some(22), "{coverage}");
    assert_eq!(coverage["autosomes"], 0.726, "31,950 bases over 44,000: {coverage}");
    assert_eq!((&coverage["chrx"], &coverage["chry"]), (&0.3.into(), &0.3.into()), "{coverage}");
    for missing in ["snvs", "indels", "svs", "phasing", "methylation"] {
        assert!(qc[missing].is_null(), "{missing}: {qc}");
    }
    assert_eq!(qc["not_measured"].as_array().map(Vec::len), Some(4), "{qc}");
    let tsv = text(&dir.join("out/qc/S.qc.tsv"));
    for row in
        ["sample\tS\n", "reads.mapped_reads\t230\n", "sex.agrees\tfalse\n", "snvs.pass\tNA\n", "coverage.chr22\t"]
    {
        assert!(tsv.contains(row), "{row}: {tsv}");
    }

    // Run as male: no warning, and the QC alone runs again (the depth does not read the sex).
    let out = lrpipe_run(&dir, &[&args[..], &["--sex", "male", "--analyses", "depth"]].concat());
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(!stderr(&out).contains("warning: the reads' index"), "{}", stderr(&out));
    assert!(stderr(&out).contains("; sex male, as the index says\n"), "{}", stderr(&out));
    let qc = json(&dir.join("out/qc/S.qc.json"));
    assert_eq!(qc["sex"]["agrees"], true);
    assert_eq!(qc["warnings"], serde_json::json!([]));
    let manifest = json(&dir.join("out/manifest.json"));
    let statuses: Vec<&str> =
        manifest["stages"].as_array().expect("stages").iter().map(|s| s["status"].as_str().unwrap()).collect();
    assert_eq!(statuses, ["skipped", "succeeded"], "{manifest}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn the_unplaced_unmapped_reads_are_counted_from_a_bai_and_a_crams_shared_container() {
    use pasta::cram;
    use pasta::sam::alignment::io::Write as _;
    let dir = temp_dir("unmapped");
    let bam = male_bam(&dir);
    assert_eq!(lrpipe::sex::unplaced_unmapped(&bam).expect("read"), Some(2));
    // The same reads as CRAM, whose index counts no records. Their containers' headers do, but
    // this small file's one container holds the mapped reads too: its records are decoded
    // (tests/sex.rs counts containers of unmapped reads alone).
    let reference = lrpipe_reads::Reference::open(&dir.join("ref.fa")).expect("the reference");
    let path = dir.join("male.cram");
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(reference.repository().clone())
        .build_from_path(&path)
        .expect("a CRAM");
    let mut reader = pasta::bam::io::Reader::new(fs::File::open(&bam).expect("the BAM"));
    let header = reader.read_header().expect("a header");
    writer.write_alignment_header(&header).expect("a header");
    let mut record = pasta::bam::Record::default();
    while reader.read_record(&mut record).expect("a record") != 0 {
        writer.write_alignment_record(&header, &record).expect("written");
    }
    writer.finish(&header).expect("finished");
    drop(writer);
    let index = cram::fs::index(&path).expect("indexed");
    cram::crai::fs::write(dir.join("male.cram.crai"), &index).expect("index written");
    assert_eq!(lrpipe::sex::unplaced_unmapped(&path).expect("read"), Some(2));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn what_cannot_be_measured_is_said_and_a_broken_index_is_warned_of() {
    let dir = temp_dir("unmeasured");
    let bam = male_bam(&dir);
    // The read pass's statistics of reads without qualities.
    let stats = r#"{"reads": 3, "bases": 30, "n50": 10, "mean_length": 10.0, "mean_quality": null,
        "quality_bases": 0, "secondary": 0, "supplementary": 0, "unmapped_placed": 0}"#;
    fs::write(dir.join("stats.json"), stats).expect("written");
    let qc = |dir: &Path| {
        let out = Command::new(env!("CARGO_BIN_EXE_lrpipe"))
            .current_dir(dir)
            .args(["qc", "--sample", "S", "--sex", "male", "--reads"])
            .arg(&bam)
            .args(["--read-stats", "stats.json", "--json", "qc.json", "--tsv", "qc.tsv"])
            .output()
            .expect("lrpipe qc runs");
        assert!(out.status.success(), "{}", stderr(&out));
        json(&dir.join("qc.json"))
    };
    let summary = qc(&dir);
    let notes: Vec<&str> =
        summary["not_measured"].as_array().expect("notes").iter().filter_map(|n| n.as_str()).collect();
    assert!(notes.contains(&"reads.mean_quality: the reads have no base qualities"), "{notes:?}");
    assert_eq!(summary["reads"]["unmapped"], 2);
    assert_eq!(summary["warnings"], serde_json::json!([]));
    // A BAI cut short: the unmapped reads are not measured, and why is a warning, not hidden.
    let bai = fs::read(dir.join("male.bam.bai")).expect("the index");
    fs::write(dir.join("male.bam.bai"), &bai[..bai.len() / 2]).expect("written");
    let summary = qc(&dir);
    assert!(summary["reads"]["unmapped"].is_null(), "{summary}");
    let warnings = summary["warnings"].as_array().expect("warnings");
    assert!(
        warnings.iter().any(|w| w.as_str().is_some_and(|w| w.starts_with("the reads' unplaced unmapped reads cannot"))),
        "{summary}"
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_run_in_regions_is_summarised_within_them() {
    // F9 with F7: the read pass counts the reads reaching the regions, the depth is of the
    // regions alone, and the summary says it is of them; the unmapped reads, of the whole file,
    // are not compared with them.
    let dir = temp_dir("regions");
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let reference = fixture.join("chr20_window.fa.gz");
    fs::write(dir.join("site.toml"), format!("reference = \"{}\"\n[snv]\n", reference.display())).expect("site");
    let (start, end) = (30_180_000_i64, 30_200_000_i64);
    fs::write(dir.join("window.bed"), format!("chr20\t{start}\t{end}\n")).expect("BED");
    let bam = fixture.join("ont.bam");
    let bam_arg = bam.to_str().expect("UTF-8");
    let args = ["--platform", "ont", "--sex", "female", "--config", "site.toml", "--reads", bam_arg, "--output", "out"];
    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", "sv,depth", "--regions", "window.bed"]].concat());
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stderr(&out).contains("QC: within the region (20.0 kb), every value but the sex"), "{}", stderr(&out));
    let qc = json(&dir.join("out/qc/S.qc.json"));
    assert_eq!(qc["regions"], serde_json::json!({"regions": 1, "bases": 20_000}));
    // The primary alignments reaching the window.
    use pasta::sam::alignment::Record as _;
    let mut reader = pasta::bam::io::Reader::new(fs::File::open(&bam).expect("the BAM"));
    reader.read_header().expect("a header");
    let mut record = pasta::bam::Record::default();
    let mut reaching = 0;
    while reader.read_record(&mut record).expect("a record") != 0 {
        let flags = record.flags();
        let from = record.alignment_start().transpose().expect("a start").map_or(0, |p| usize::from(p) as i64 - 1);
        let to = from + record.alignment_span().transpose().expect("a span").unwrap_or(0) as i64;
        if !flags.is_unmapped() && !flags.is_secondary() && !flags.is_supplementary() && from < end && to > start {
            reaching += 1;
        }
    }
    assert!(reaching > 10 && reaching < 152, "{reaching}");
    assert_eq!(qc["reads"]["mapped_reads"], reaching, "{}", qc["reads"]);
    assert!(qc["reads"]["unmapped"].is_null() && qc["reads"]["mapped_percent"].is_null(), "{}", qc["reads"]);
    let notes = qc["not_measured"].to_string();
    assert!(notes.contains("reads.unmapped: the reads are counted within the regions"), "{notes}");
    // The depth of the window's reads over the window.
    let depth = text(&dir.join("out/depth/S.mosdepth.summary.txt"));
    assert!(depth.lines().all(|l| !l.starts_with("chr20\t")), "regions' rows alone: {depth}");
    let row = depth.lines().find(|l| l.starts_with("chr20_region\t")).expect("the window's row");
    let mean: f64 = row.split('\t').nth(3).expect("a mean").parse().expect("a number");
    assert_eq!(qc["coverage"]["per_autosome"]["chr20"], mean);
    assert!(mean > 5.0, "{row}");
    let _ = fs::remove_dir_all(dir);
}

/// A VCF's PASS records: REF, ALTs and genotype.
fn pass_records(path: &Path) -> Vec<(String, Vec<String>, String)> {
    text(path)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[6] == "PASS")
        .map(|f| {
            let gt = f[9].split(':').next().expect("GT first").to_owned();
            (f[3].to_owned(), f[4].split(',').map(str::to_owned).collect(), gt)
        })
        .collect()
}

#[test]
fn the_fixtures_summary_agrees_with_its_files() {
    // F9 acceptance on the ONT fixture: every value of the summary against its own count here.
    // The SNVs and their phasing need SNV_MODEL (tools/check.sh sets it); without it they are
    // left out, saying so.
    let dir = temp_dir("fixture");
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let model = std::env::var_os("SNV_MODEL").map(PathBuf::from).filter(|m| m.join("pileup.pt").exists());
    let snv = match &model {
        Some(model) => format!("ont_model = \"{}\"\ndevice = \"cpu\"\n", model.display()),
        None => {
            eprintln!("SNV_MODEL not set (or not a Clair3 model): the SNVs and their phasing are not checked");
            String::new()
        }
    };
    let reference = fixture.join("chr20_window.fa.gz");
    fs::write(dir.join("site.toml"), format!("reference = \"{}\"\n[snv]\n{snv}", reference.display())).expect("site");
    let analyses = if model.is_some() { "snv,phase,sv,depth,methyl" } else { "sv,depth,methyl" };
    let bam = fixture.join("ont.bam");
    let bam_arg = bam.to_str().expect("UTF-8");
    let args = ["--platform", "ont", "--sex", "female", "--config", "site.toml", "--reads", bam_arg, "--output", "out"];
    let out = lrpipe_run(&dir, &[&args[..], &["--analyses", analyses]].concat());
    assert!(out.status.success(), "{}", stderr(&out));
    let err = stderr(&out);
    assert!(err.contains("QC: 152 mapped reads (3.35 Mb), read N50 40.2 kb"), "{err}");
    let out_dir = dir.join("out");
    let qc = json(&out_dir.join("qc/S.qc.json"));
    let manifest = json(&out_dir.join("manifest.json"));
    assert_eq!(manifest["outputs"]["qc_json"], out_dir.join("qc/S.qc.json").to_str().expect("UTF-8"));

    check_reads(&qc, &bam, 0);
    // chr20's mean depth: the reads' spans (MAPQ 20 or more) over its length.
    use pasta::sam::alignment::Record as _;
    let mut reader = pasta::bam::io::Reader::new(fs::File::open(&bam).expect("the BAM"));
    reader.read_header().expect("a header");
    let mut record = pasta::bam::Record::default();
    let mut span = 0usize;
    while reader.read_record(&mut record).expect("a record") != 0 {
        let flags = record.flags();
        let mapq = record.mapping_quality().map_or(255, u8::from);
        if mapq >= 20 && !flags.is_secondary() && !flags.is_unmapped() {
            span += record.alignment_span().expect("placed").expect("a span");
        }
    }
    let coverage = &qc["coverage"];
    let mean = coverage["autosomes"].as_f64().expect("the autosomes' depth");
    assert!((mean - span as f64 / 64_444_167.0).abs() < 5e-4, "{mean} against {span} bases");
    assert_eq!(coverage["per_autosome"]["chr20"], 0.05, "mosdepth's two decimals: {coverage}");
    // The depth is of the reference's contigs, which are chr20 alone.
    assert!(coverage["chrx"].is_null() && coverage["chry"].is_null(), "{coverage}");
    // A chr20 window cannot tell the sex.
    let sex = &qc["sex"];
    assert!(sex["index"].is_null() && sex["agrees"].is_null(), "{sex}");
    assert!(sex["note"].as_str().is_some_and(|n| n.contains("1 of the 22 autosomes (chr20)")), "{sex}");
    assert_eq!(qc["warnings"], serde_json::json!([]));

    // The fixture has no SVs.
    assert_eq!(qc["svs"]["pass"], 0, "{}", qc["svs"]);
    // bedMethyl: a row of 5mC and one of 5hmC at each CpG.
    let rows = text(&out_dir.join("methylation/S.bedmethyl.gz"));
    let mut cpgs: Vec<(&str, &str)> = rows
        .lines()
        .map(|l| {
            let mut f = l.split('\t');
            (f.next().unwrap(), f.next().unwrap())
        })
        .collect();
    let rows_n = cpgs.len();
    cpgs.dedup();
    assert!(cpgs.len() > 100 && cpgs.len() < rows_n, "{} CpGs of {rows_n} rows", cpgs.len());
    assert_eq!(qc["methylation"]["combined"], cpgs.len(), "{}", qc["methylation"]);
    assert!(qc["methylation"]["hap1"].is_null());

    if model.is_none() {
        assert!(qc["snvs"].is_null() && qc["phasing"].is_null());
        let _ = fs::remove_dir_all(dir);
        return;
    }
    let records = pass_records(&out_dir.join("snv/merge_output.vcf.gz"));
    let snv = |r: &(String, Vec<String>, String)| r.0.len() == 1 && r.1.iter().all(|a| a.len() == 1);
    let snvs: Vec<_> = records.iter().filter(|r| snv(r)).collect();
    let transitions = snvs
        .iter()
        .flat_map(|r| r.1.iter().map(move |a| (r.0.as_str(), a.as_str())))
        .filter(|p| matches!(p, ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C")))
        .count();
    let alleles: usize = snvs.iter().map(|r| r.1.len()).sum();
    let het = |r: &&&(String, Vec<String>, String)| {
        let a: Vec<&str> = r.2.split(['/', '|']).collect();
        a.len() == 2 && a[0] != a[1]
    };
    let hom = |r: &&&(String, Vec<String>, String)| matches!(r.2.as_str(), "1/1" | "1|1" | "2/2" | "2|2");
    let s = &qc["snvs"];
    assert!(snvs.len() > 10, "{s}");
    assert_eq!(s["pass"], snvs.len(), "{s}");
    assert_eq!((&s["transitions"], &s["transversions"]), (&transitions.into(), &(alleles - transitions).into()));
    assert_eq!(s["het"], snvs.iter().filter(het).count(), "{s}");
    assert_eq!(s["hom_alt"], snvs.iter().filter(hom).count(), "{s}");
    let ti_tv = s["ti_tv"].as_f64().expect("Ti/Tv");
    assert!((ti_tv - transitions as f64 / (alleles - transitions) as f64).abs() < 1e-3, "{s}");
    assert_eq!(qc["indels"]["pass"], records.len() - snvs.len(), "{}", qc["indels"]);

    // Phasing: the heterozygous biallelic SNVs, and those in a phase set of two or more
    // variants.
    let phased = text(&out_dir.join("snv/phased_merge_output.vcf.gz"));
    let mut sets: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut het_snvs = Vec::new();
    for l in phased.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split('\t').collect();
        let (keys, values) = (f[8].split(':').collect::<Vec<_>>(), f[9].split(':').collect::<Vec<_>>());
        let gt = values[0];
        let ps = keys.iter().position(|k| *k == "PS").and_then(|i| values.get(i)).filter(|p| **p != ".");
        let biallelic = !f[4].contains(',');
        let alleles: Vec<&str> = gt.split(['/', '|']).collect();
        if !biallelic || alleles.len() != 2 || alleles[0] == alleles[1] {
            continue;
        }
        let set = ps.filter(|_| gt.contains('|')).map(|p| (*p).to_owned());
        if let Some(set) = &set {
            *sets.entry(set.clone()).or_default() += 1;
        }
        if f[3].len() == 1 && f[4].len() == 1 {
            het_snvs.push(set);
        }
    }
    let in_blocks = het_snvs.iter().filter(|s| s.as_ref().is_some_and(|s| sets[s] > 1)).count();
    let p = &qc["phasing"];
    assert_eq!(p["het_snvs_any_filter"], het_snvs.len(), "{p}");
    assert_eq!(p["phased_snvs"], in_blocks, "{p}");
    assert!(in_blocks > 0 && p["blocks"].as_u64() > Some(0), "{p}");
    // Blocks of the window span far less than half of chr20: its NG50 is 0, as WhatsHap's.
    assert_eq!(p["block_n50"], 0, "{p}");
    let _ = fs::remove_dir_all(dir);
}

#[cfg(unix)]
#[test]
fn an_output_linked_to_an_input_is_refused() {
    let dir = temp_dir("linked");
    let stats = r#"{"reads": 3}"#;
    fs::write(dir.join("stats.json"), stats).expect("written");
    std::os::unix::fs::symlink("stats.json", dir.join("qc.json")).expect("linked");
    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(&dir)
        .args(["qc", "--sample", "S", "--sex", "male", "--reads", "absent.bam"])
        .args(["--read-stats", "stats.json", "--json", "qc.json", "--tsv", "qc.tsv"])
        .output()
        .expect("lrpipe qc runs");
    assert!(!out.status.success());
    assert!(stderr(&out).contains("--json qc.json would overwrite the --read-stats input"), "{}", stderr(&out));
    assert_eq!(fs::read_to_string(dir.join("stats.json")).expect("kept"), stats);
    let _ = fs::remove_dir_all(dir);
}

//! `lrpipe setup --bonefish` and `lrpipe run` from unaligned reads with the real bonefish: a tiny
//! reference, and two unaligned BAMs of reads drawn from it with base modification (`MM`/`ML`/`MN`)
//! tags and dorado's read groups, as dorado writes them.
//!
//! bonefish is `$BONEFISH`, else `bonefish` on `PATH`; without it the test is skipped, saying so,
//! unless `LRPIPE_REQUIRE_BONEFISH` is set (as `tools/check.sh` sets it where bonefish is), when it
//! fails.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use pasta::bam;
use pasta::sam::alignment::io::Write as _;
use pasta::sam::alignment::record::Flags;
use pasta::sam::alignment::record::data::field::Tag;
use pasta::sam::alignment::record_buf::data::field::Value;
use pasta::sam::alignment::record_buf::data::field::value::Array;
use pasta::sam::alignment::record_buf::{QualityScores, RecordBuf, Sequence};

use lrpipe::plan::Platform;
use lrpipe::setup::{self, Options};
use lrpipe::snv_model::{self, Mismatch, Request};

/// The read groups of the two unaligned BAMs: two flowcells, basecalled alike.
const READ_GROUPS: [&str; 2] = [
    "@RG\tID:fc1_sup\tPU:FC1\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
     modbase_models=dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v2.0.1 runid=fc1\tLB:S\tSM:S",
    "@RG\tID:fc2_sup\tPU:FC2\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
     modbase_models=dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v2.0.1 runid=fc2\tLB:S\tSM:S",
];

/// The reference's length: long enough for bonefish to place each read once.
const REFERENCE_LENGTH: usize = 40_000;

/// The bonefish to test with, if any.
fn bonefish() -> Option<PathBuf> {
    match std::env::var_os("BONEFISH") {
        Some(path) => Some(PathBuf::from(path)),
        None => lrpipe::provenance::resolve(Path::new("bonefish")),
    }
    .filter(|p| p.is_file())
}

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-align-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// A random sequence of `n` bases (a fixed linear congruential generator).
fn random_bases(n: usize, seed: &mut u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[(*seed >> 62) as usize]
        })
        .collect()
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// The read's 5mC calls as dorado writes them: every C of a CpG (in the read's own orientation)
/// called methylated, with `MM`'s counts of the Cs skipped between them, one `ML` probability
/// each.
fn modifications(read: &[u8]) -> (String, Vec<u8>) {
    let (mut mm, mut ml, mut skipped) = (String::from("C+m?"), Vec::new(), 0);
    for (i, &b) in read.iter().enumerate() {
        if b != b'C' {
            continue;
        }
        if read.get(i + 1) == Some(&b'G') {
            mm.push_str(&format!(",{skipped}"));
            ml.push(230);
            skipped = 0;
        } else {
            skipped += 1;
        }
    }
    mm.push(';');
    (mm, ml)
}

/// A read's tags as they were written, to compare with the alignments'.
#[derive(Debug, PartialEq, Eq)]
struct Tags {
    mm: String,
    ml: Vec<u8>,
    mn: i64,
    rg: String,
}

fn tags(record: &RecordBuf) -> Tags {
    let data = record.data();
    let get = |tag: &[u8; 2]| data.get(&Tag::new(tag[0], tag[1]));
    let text = |tag: &[u8; 2]| match get(tag) {
        Some(Value::String(s)) => s.to_string(),
        other => panic!("{}: {other:?}", String::from_utf8_lossy(tag)),
    };
    let ml = match get(b"ML") {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("ML: {other:?}"),
    };
    let mn = get(b"MN").and_then(Value::as_int).expect("MN");
    Tags { mm: text(b"MM"), ml, mn, rg: text(b"RG") }
}

/// Writes the reference and its `.fai`, and two unaligned BAMs of reads from it, half of them
/// from its reverse strand, and the same reads as dorado's FASTQ (`reads.fastq`: the tags in each
/// record's comment, with `RG` naming a read group no header holds). Returns the reads' tags by
/// name.
fn sample(dir: &Path) -> BTreeMap<String, Tags> {
    let mut seed = 2026;
    let reference = random_bases(REFERENCE_LENGTH, &mut seed);
    let mut fasta = String::from(">chr1\n");
    for line in reference.chunks(60) {
        fasta.push_str(std::str::from_utf8(line).expect("ASCII"));
        fasta.push('\n');
    }
    fs::write(dir.join("ref.fa"), fasta).expect("reference");
    fs::write(dir.join("ref.fa.fai"), format!("chr1\t{REFERENCE_LENGTH}\t6\t60\t61\n")).expect("index");

    let mut written = BTreeMap::new();
    let mut fastq = String::new();
    for (file, read_group) in READ_GROUPS.iter().enumerate() {
        let id = read_group.split('\t').nth(1).and_then(|f| f.strip_prefix("ID:")).expect("an ID");
        let header: pasta::sam::Header = format!("@HD\tVN:1.6\tSO:unknown\n{read_group}\n").parse().expect("header");
        let path = dir.join(format!("flowcell{}.bam", file + 1));
        let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("BAM"));
        writer.write_header(&header).expect("header");
        for i in 0..20 {
            let start = (seed >> 33) as usize % (REFERENCE_LENGTH - 6_000);
            let length = 3_000 + (seed >> 40) as usize % 3_000;
            random_bases(1, &mut seed);
            let forward = &reference[start..start + length];
            let read = if i % 2 == 0 { forward.to_vec() } else { reverse_complement(forward) };
            let (mm, ml) = modifications(&read);
            let name = format!("fc{}_read{i}", file + 1);
            let data = [
                (Tag::new(b'M', b'M'), Value::String(mm.clone().into())),
                (Tag::new(b'M', b'L'), Value::Array(Array::UInt8(ml.clone()))),
                (Tag::new(b'M', b'N'), Value::from(i32::try_from(length).expect("short"))),
                (Tag::new(b'R', b'G'), Value::String(id.into())),
                (Tag::new(b'q', b's'), Value::from(20)),
            ];
            let record = RecordBuf::builder()
                .set_name(name.clone())
                .set_flags(Flags::UNMAPPED)
                .set_sequence(Sequence::from(read))
                .set_quality_scores(QualityScores::from(vec![30; length]))
                .set_data(data.into_iter().collect())
                .build();
            writer.write_alignment_record(&header, &record).expect("record");
            let probabilities: Vec<String> = ml.iter().map(u8::to_string).collect();
            fastq += &format!(
                "@{name}\tqs:i:20\tMN:i:{length}\tRG:Z:{id}\tMM:Z:{mm}\tML:B:C,{}\n{}\n+\n{}\n",
                probabilities.join(","),
                String::from_utf8_lossy(record.sequence().as_ref()),
                "?".repeat(length)
            );
            let mn = i64::try_from(length).expect("short");
            written.insert(name, Tags { mm, ml, mn, rg: id.to_owned() });
        }
        writer.try_finish().expect("finished");
    }
    fs::write(dir.join("reads.fastq"), fastq).expect("FASTQ");
    written
}

/// A manifest of the tiny genome: the reference and an ONT model, as `file://` URLs in `dir`.
fn manifest(dir: &Path) -> String {
    let mut text = String::from("format = 1\ngenome = \"tiny\"\n\n[site]\nreference = \"reference\"\n\n");
    text += "[models]\nont = \"r1041_e82_400bps_sup_v500\"\n\n";
    let mut entry = |id: &str, file: &str, path: &str, extra: &str| {
        let data = fs::read(dir.join(file)).expect("published");
        let sha256 = lrpipe::provenance::sha256(&dir.join(file)).expect("hashed");
        text += &format!(
            "[[file]]\nid = \"{id}\"\nurl = \"file://{}\"\npath = \"{path}\"\nsize = {}\nsha256 = \"{sha256}\"\nlicence = \
             \"test\"\n{extra}\n",
            dir.join(file).display(),
            data.len()
        );
    };
    entry("reference", "ref.fa", "tiny.fa", "");
    for weights in ["pileup", "full_alignment"] {
        fs::write(dir.join(format!("{weights}.pt")), weights).expect("weights");
        let model = "r1041_e82_400bps_sup_v500";
        entry(
            &format!("{model}/{weights}"),
            &format!("{weights}.pt"),
            &format!("models/{model}/{weights}.pt"),
            &format!("model = \"{model}\"\nplatform = \"ont\"\n"),
        );
    }
    text
}

#[test]
fn unaligned_reads_are_aligned_by_bonefish_keeping_their_tags_and_read_groups() {
    let Some(bonefish) = bonefish() else {
        assert!(
            std::env::var_os("LRPIPE_REQUIRE_BONEFISH").is_none(),
            "LRPIPE_REQUIRE_BONEFISH is set, but bonefish is not found ($BONEFISH or on PATH)"
        );
        eprintln!("skipped: no bonefish ($BONEFISH, or bonefish on PATH)");
        return;
    };
    let dir = temp_dir("bonefish");
    let written = sample(&dir);

    // lrpipe setup builds bonefish's index of the reference and names both in site.toml.
    let options = Options {
        dir: dir.join("site"),
        platforms: vec![Platform::Ont],
        accept_licences: false,
        bonefish: Some(bonefish.clone()),
    };
    let site = setup::setup(&manifest(&dir), &options, &setup::fetch::Fetcher::new()).expect("set up");
    let index = fs::canonicalize(dir.join("site")).expect("set up").join("tiny.fa.bonefish.raw");
    assert!(index.is_file(), "bonefish's index");

    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(&dir)
        .args(["run", "--config", site.to_str().expect("UTF-8")])
        .args(["--reads", "flowcell1.bam", "--reads", "flowcell2.bam", "--sample", "S", "--platform", "ont"])
        .args(["--sex", "female", "--analyses", "sv,depth,methyl", "--output", "out", "--memory", "64"])
        .args(["--threads", "4"])
        .output()
        .expect("lrpipe runs");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "{stderr}");

    // One sorted, indexed BAM of both files, with both read groups, which name the basecall model.
    let aligned = dir.join("out/align/S.bam");
    assert!(dir.join("out/align/S.bam.bai").is_file(), "indexed");
    assert!(!dir.join("out/align/tmp").exists(), "the scratch directory is removed");
    let mut reader = bam::io::Reader::new(fs::File::open(&aligned).expect("aligned BAM"));
    let header = reader.read_header().expect("header");
    assert_eq!(header.reference_sequences().len(), 1);
    let text = snv_model::header_text(&aligned).expect("header");
    for rg in READ_GROUPS {
        assert!(text.contains(rg), "{rg} kept:\n{text}");
    }
    // F3 chooses the model from the aligned BAM's header as from the unaligned ones'.
    let models = dir.join("models");
    fs::create_dir_all(models.join("r1041_e82_400bps_sup_v500")).expect("models");
    let request = Request {
        platform: Platform::Ont,
        given: None,
        models: Some(&models),
        mismatch: Mismatch::Error,
        fastq: false,
    };
    let choice = snv_model::choose(&text, &request).expect("a model");
    assert!(choice.path.ends_with("r1041_e82_400bps_sup_v500"), "{choice:?}");

    // Every read mapped, its base modifications and read group as written.
    let mut seen = BTreeMap::new();
    for record in reader.record_bufs(&header) {
        let record = record.expect("a record");
        if record.flags().is_secondary() || record.flags().is_supplementary() {
            continue;
        }
        assert!(!record.flags().is_unmapped(), "{:?} unmapped", record.name());
        let name = record.name().expect("a name").to_string();
        seen.insert(name, tags(&record));
    }
    assert_eq!(seen, written, "every read's MM, ML, MN and RG kept");

    // The later stages read the aligned BAM: depth and bedMethyl of its reads.
    let manifest: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("out/manifest.json")).expect("manifest")).expect("JSON");
    assert_eq!(manifest["status"], "complete", "{manifest}");
    let stages = manifest["stages"].as_array().expect("stages");
    assert_eq!(stages[0]["name"], "align");
    let command = stages[0]["command"].to_string();
    assert!(command.contains("\"-x\",\"dna\",\"-y\"") && command.contains(index.to_str().expect("UTF-8")), "{command}");
    assert!(stages[1]["command"].to_string().contains(aligned.to_str().expect("UTF-8")));
    assert!(manifest["tools"]["bonefish"]["version"].as_str().expect("a version").starts_with("bonefish "));
    let bedmethyl = dir.join("out/methylation/S.bedmethyl.gz");
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut flate2::read::MultiGzDecoder::new(fs::File::open(&bedmethyl).expect("bedMethyl")),
        &mut text,
    )
    .expect("read");
    let methylated = text.lines().filter(|l| l.split('\t').nth(3) == Some("m")).count();
    assert!(methylated > 100, "bedMethyl of the reads' CpGs: {methylated} rows");
    let summary = fs::read_to_string(dir.join("out/depth/S.mosdepth.summary.txt")).expect("depth");
    let mean: f64 = summary.lines().nth(1).and_then(|l| l.split('\t').nth(3)?.parse().ok()).expect("chr1's depth");
    assert!(mean > 2.0, "{summary}");

    // dorado's FASTQ: its reads' RG tags name read groups no header holds, so bonefish gives them
    // one named after the sample, and keeps their base modifications.
    let out = Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .current_dir(&dir)
        .args(["run", "--config", site.to_str().expect("UTF-8"), "--reads", "reads.fastq", "--sample", "S"])
        .args(["--platform", "ont", "--sex", "female", "--analyses", "depth,methyl", "--output", "fastq"])
        .args(["--memory", "64", "--threads", "4"])
        .output()
        .expect("lrpipe runs");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let aligned = dir.join("fastq/align/S.bam");
    let text = snv_model::header_text(&aligned).expect("header");
    let groups: Vec<&str> = text.lines().filter(|l| l.starts_with("@RG")).collect();
    assert_eq!(groups, ["@RG\tID:S\tSM:S"], "{text}");
    let mut reader = bam::io::Reader::new(fs::File::open(&aligned).expect("aligned BAM"));
    let header = reader.read_header().expect("header");
    let mut reads = 0;
    for record in reader.record_bufs(&header) {
        let record = record.expect("a record");
        if record.flags().is_secondary() || record.flags().is_supplementary() {
            continue;
        }
        let name = record.name().expect("a name").to_string();
        let want = &written[&name];
        assert_eq!(
            tags(&record),
            Tags { rg: "S".into(), mm: want.mm.clone(), ml: want.ml.clone(), mn: want.mn },
            "{name}"
        );
        reads += 1;
    }
    assert_eq!(reads, written.len());
    assert!(dir.join("fastq/methylation/S.bedmethyl.gz").is_file(), "the FASTQ's tags give methylation");
    let _ = fs::remove_dir_all(dir);
}

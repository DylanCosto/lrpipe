//! Runs lrpipe-sv on the fixture in `tests/data` (a slice of HG002 ONT reads around
//! chr17:46.2-46.4 Mb) and compares its VCFs with the ones Sniffles 2.8.1 wrote for the same input
//! and options (`tests/data/expected`, made by `tools/make_fixture.sh`). Header lines that name the
//! program, command or date are ignored; everything else must match exactly.

mod common;

use std::{
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

use common::{bgzip_reference, write_reference};

const IGNORED_HEADER_PREFIXES: &[&str] =
    &["##source=", "##command=", "##fileDate=", "##genotypeSource=", "##genotypeCommand=", "##genotypeFileDate="];

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// A scratch directory with the fixture inputs and the full-length reference, built once.
fn work_dir() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("parity");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).expect("creating the scratch directory");
        for name in [
            "sample1.bam",
            "sample1.bam.bai",
            "sample1.cram",
            "sample1.cram.crai",
            "sample2.cram",
            "sample2.cram.crai",
            "regions.bed",
        ] {
            fs::copy(data_dir().join(name), dir.join(name)).unwrap_or_else(|e| panic!("copying {name}: {e}"));
        }
        write_reference(&data_dir().join("reference.window.fa"), &dir.join("reference.fa"), &[]);
        bgzip_reference(&dir.join("reference.fa"), &dir.join("reference.fa.gz"));
        // Bgzipped, under a name that does not say so.
        for suffix in ["", ".fai", ".gzi"] {
            fs::copy(dir.join(format!("reference.fa.gz{suffix}")), dir.join(format!("reference.bgzf.fa{suffix}")))
                .expect("copying the bgzipped reference");
        }
        dir
    })
}

/// Runs lrpipe-sv in the scratch directory and panics with its output if it fails.
fn sniffles(args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_lrpipe-sv"))
        .current_dir(work_dir())
        .args(args)
        .arg("--allow-overwrite")
        .output()
        .expect("running lrpipe-sv");
    assert!(
        output.status.success(),
        "lrpipe-sv {} failed:\n{}{}",
        args.join(" "),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn comparable(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !IGNORED_HEADER_PREFIXES.iter().any(|p| l.starts_with(p))).collect()
}

/// Asserts that `actual` (in the scratch directory) matches the expected VCF `expected`.
fn assert_vcf_eq(expected: &str, actual: &str) {
    assert_fixed_vcf_eq(expected, actual, str::to_owned, |_| true);
}

/// As [`assert_vcf_eq`] where lrpipe-sv fixes a Sniffles2 bug: each of Sniffles2's expected lines
/// becomes `fixed(line)`, and only lrpipe-sv's lines for which `compared` holds are compared.
fn assert_fixed_vcf_eq(expected: &str, actual: &str, fixed: impl Fn(&str) -> String, compared: impl Fn(&str) -> bool) {
    let expected_text = fs::read_to_string(data_dir().join("expected").join(expected)).expect("reading expected VCF");
    let actual_text = read_vcf(&work_dir().join(actual));
    let e: Vec<String> = comparable(&expected_text).into_iter().map(fixed).collect();
    let a: Vec<&str> = comparable(&actual_text).into_iter().filter(|l| compared(l)).collect();
    if e != a {
        let first = e.iter().zip(&a).position(|(x, y)| x != y).unwrap_or(e.len().min(a.len()));
        panic!(
            "{actual} differs from {expected} ({} vs {} lines), first at line {}:\n  expected: {}\n  actual:   {}",
            e.len(),
            a.len(),
            first + 1,
            e.get(first).map_or("<end>", String::as_str),
            a.get(first).unwrap_or(&"<end>"),
        );
    }
}

fn read_vcf(path: &Path) -> String {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    if path.extension().is_some_and(|e| e == "gz") {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..]).read_to_string(&mut text).expect("decompressing VCF");
        text
    } else {
        String::from_utf8(bytes).expect("VCF is UTF-8")
    }
}

#[test]
fn default_options_any_thread_count() {
    for threads in ["1", "2", "4"] {
        let out = format!("default.t{threads}.vcf");
        sniffles(&["--input", "sample1.bam", "--vcf", &out, "--threads", threads]);
        assert_vcf_eq("default.vcf", &out);
    }
}

#[test]
fn reference() {
    sniffles(&["--input", "sample1.bam", "--vcf", "reference.vcf", "--reference", "reference.fa", "--threads", "4"]);
    assert_vcf_eq("reference.vcf", "reference.vcf");
}

#[test]
fn cram() {
    sniffles(&["--input", "sample1.cram", "--vcf", "cram.vcf", "--reference", "reference.fa", "--threads", "4"]);
    assert_vcf_eq("cram.vcf", "cram.vcf");
}

#[test]
fn phase() {
    sniffles(&["--input", "sample1.bam", "--vcf", "phase.vcf", "--phase"]);
    assert_vcf_eq("phase.vcf", "phase.vcf");
}

#[test]
fn mosaic() {
    sniffles(&["--input", "sample1.bam", "--vcf", "mosaic.vcf", "--mosaic"]);
    assert_vcf_eq("mosaic.vcf", "mosaic.vcf");
}

#[test]
fn no_qc() {
    sniffles(&["--input", "sample1.bam", "--vcf", "no_qc.vcf", "--no-qc"]);
    assert_vcf_eq("no_qc.vcf", "no_qc.vcf");
}

#[test]
fn symbolic() {
    sniffles(&["--input", "sample1.bam", "--vcf", "symbolic.vcf", "--symbolic"]);
    assert_vcf_eq("symbolic.vcf", "symbolic.vcf");
}

#[test]
fn regions() {
    sniffles(&["--input", "sample1.bam", "--vcf", "regions.vcf", "--regions", "regions.bed"]);
    assert_vcf_eq("regions.vcf", "regions.vcf");
}

#[test]
fn genotype_vcf() {
    let known = data_dir().join("expected/default.vcf");
    let known = known.to_str().expect("UTF-8 path");
    sniffles(&[
        "--input",
        "sample2.cram",
        "--vcf",
        "genotype.vcf",
        "--reference",
        "reference.fa",
        "--genotype-vcf",
        known,
    ]);
    // Sniffles2 leaves PS out of FORMAT though every genotype ends with it (LRP-008). It gives
    // the breakend the support of another candidate with a mate on chr17, 1 read; the breakend
    // itself, as sample2's own calls have it at this position with this ALT, has 8 (LRP-015).
    assert_fixed_vcf_eq(
        "genotype.vcf",
        "genotype.vcf",
        |l| {
            l.replace("\tGT:GQ:DR:DV\t", "\tGT:GQ:DR:DV:PS\t").replace(
                "\t0/0:26:13:1:.",
                if l.contains("Sniffles2.BND.20S10") { "\t0/1:37:6:8:." } else { "\t0/0:26:13:1:." },
            )
        },
        |_| true,
    );
    assert_genotyped_vcf_is_whole(&read_vcf(&work_dir().join("genotype.vcf")));
}

/// The genotyped VCF's header describes its records: one sample column on the `#CHROM` line as in
/// every record, and every FORMAT key declared (LRP-008).
fn assert_genotyped_vcf_is_whole(text: &str) {
    let declared: Vec<&str> =
        text.lines().filter_map(|l| l.strip_prefix("##FORMAT=<ID=")).filter_map(|l| l.split(',').next()).collect();
    let chrom = text.lines().find(|l| l.starts_with("#CHROM")).expect("a #CHROM line");
    assert_eq!(chrom.split('\t').count(), 10, "one sample: {chrom}");
    for record in text.lines().filter(|l| !l.starts_with('#')) {
        let fields: Vec<&str> = record.split('\t').collect();
        assert_eq!(fields.len(), 10, "as many columns as the header: {record}");
        assert_eq!(fields[8].split(':').count(), fields[9].split(':').count(), "FORMAT and sample agree: {record}");
        for key in fields[8].split(':') {
            assert!(declared.contains(&key), "FORMAT {key} is declared: {record}");
        }
    }
}

/// A sites-only input (8 columns, no FORMAT definitions) and a two-sample one give the same genotypes as the single-sample
/// input, each in a VCF whose header has the one sample its records have (LRP-008).
#[test]
fn genotype_vcf_sites_only_and_multi_sample_inputs() {
    let known = fs::read_to_string(data_dir().join("expected/default.vcf")).expect("reading default.vcf");
    let (mut sites, mut two) = (String::new(), String::new());
    for line in known.lines() {
        if line.starts_with("##") {
            // A sites-only VCF declares no FORMAT fields, PS among them.
            if !line.starts_with("##FORMAT=") {
                sites.push_str(&format!("{line}\n"));
            }
            two.push_str(&format!("{line}\n"));
        } else {
            let fields: Vec<&str> = line.split('\t').collect();
            sites.push_str(&format!("{}\n", fields[..8].join("\t")));
            let other = if line.starts_with('#') { "OTHER" } else { fields[9] };
            two.push_str(&format!("{line}\t{other}\n"));
        }
    }
    let mut genotypes = Vec::new();
    for (name, text) in [("default", known.as_str()), ("sites", &sites), ("two", &two)] {
        let input = format!("known.{name}.vcf");
        let output = format!("known.{name}.genotyped.vcf");
        fs::write(work_dir().join(&input), text).expect("writing the input");
        sniffles(&[
            "--input",
            "sample2.cram",
            "--vcf",
            &output,
            "--reference",
            "reference.fa",
            "--genotype-vcf",
            &input,
        ]);
        let out = read_vcf(&work_dir().join(&output));
        assert_genotyped_vcf_is_whole(&out);
        let records: Vec<String> = out
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                format!("{}\t{}\t{}", f[1], f[2], f[9])
            })
            .collect();
        assert!(!records.is_empty(), "{name}: genotyped records");
        genotypes.push((name, records));
    }
    for (name, records) in &genotypes[1..] {
        assert_eq!(records, &genotypes[0].1, "{name} input: the same genotypes");
    }
}

/// A sites-only input whose INFO ends with SVTYPE, as catalogs often write it, gets the genotypes
/// of the same records with a sample column, under a header with FORMAT and one sample. Sniffles2
/// read the last INFO value with the line's newline (`INS\n`), skipped the record as an
/// unsupported type and wrote it 0/0, under its input's 8-column `#CHROM` line (LSV-004).
#[test]
fn genotype_vcf_sites_only_svtype_last() {
    let known = fs::read_to_string(data_dir().join("expected/default.vcf")).expect("reading default.vcf");
    let (mut sites, mut full) = (String::new(), String::new());
    for line in known.lines() {
        if line.starts_with('#') {
            if !line.starts_with("##FORMAT=") {
                let chrom: Vec<&str> = line.split('\t').collect();
                sites.push_str(&format!("{}\n", chrom[..chrom.len().min(8)].join("\t")));
            }
            full.push_str(&format!("{line}\n"));
            continue;
        }
        let mut fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        let (svtype, rest): (Vec<&str>, Vec<&str>) = fields[7].split(';').partition(|f| f.starts_with("SVTYPE="));
        fields[7] = [rest, svtype].concat().join(";");
        assert!(fields[7].contains(";SVTYPE="), "SVTYPE last: {}", fields[7]);
        sites.push_str(&format!("{}\n", fields[..8].join("\t")));
        full.push_str(&format!("{}\n", fields.join("\t")));
    }
    let mut genotypes = Vec::new();
    for (name, text) in [("full", &full), ("sites", &sites)] {
        let (input, output) = (format!("svtype_last.{name}.vcf"), format!("svtype_last.{name}.genotyped.vcf"));
        fs::write(work_dir().join(&input), text).expect("writing the input");
        sniffles(&[
            "--input",
            "sample2.cram",
            "--vcf",
            &output,
            "--reference",
            "reference.fa",
            "--genotype-vcf",
            &input,
        ]);
        let out = read_vcf(&work_dir().join(&output));
        assert_genotyped_vcf_is_whole(&out);
        let records: Vec<String> = out
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                format!("{}\t{}\t{}", f[1], f[2], f[9])
            })
            .collect();
        genotypes.push(records);
    }
    assert_eq!(genotypes[1], genotypes[0], "the sites-only input: the same genotypes");
    // The insertions and the breakend are called in sample2, so not all of them can be 0/0.
    for svtype in ["INS", "BND"] {
        assert!(
            genotypes[1].iter().any(|r| r.contains(&format!(".{svtype}.")) && !r.contains("\t0/0:")),
            "a {svtype} is genotyped: {:?}",
            genotypes[1]
        );
    }
}

#[test]
fn multi_sample_from_snf() {
    sniffles(&["--input", "sample1.bam", "--snf", "s1.snf", "--sample-id", "s1"]);
    sniffles(&["--input", "sample2.cram", "--snf", "s2.snf", "--sample-id", "s2", "--reference", "reference.fa"]);
    sniffles(&["--input", "s1.snf", "s2.snf", "--vcf", "combine.vcf"]);
    // Sniffles2 drops every breakend when combining (LRP-007); lrpipe-sv keeps them, and is
    // otherwise Sniffles2's.
    let bnd = |l: &str| !l.starts_with('#') && l.contains("SVTYPE=BND");
    assert_fixed_vcf_eq("combine.vcf", "combine.vcf", str::to_owned, |l| !bnd(l));
    let text = read_vcf(&work_dir().join("combine.vcf"));
    assert_eq!(text.lines().filter(|l| bnd(l)).count(), 1, "the breakend both samples share is kept");
}

#[test]
fn bgzipped_output_with_index() {
    sniffles(&["--input", "sample1.bam", "--vcf", "default.vcf.gz"]);
    assert_vcf_eq("default.vcf", "default.vcf.gz");
    assert!(work_dir().join("default.vcf.gz.tbi").exists(), "missing tabix index");
}

/// Not a Sniffles2 comparison: `--accuracy` changes the output on purpose, and this snapshot makes
/// any further change to it deliberate (regenerate with `tools/make_fixture.sh`).
#[test]
fn accuracy_snapshot() {
    sniffles(&["--input", "sample1.bam", "--vcf", "accuracy.vcf", "--reference", "reference.fa", "--accuracy"]);
    assert_vcf_eq("accuracy.vcf", "accuracy.vcf");
}

/// Repeat loci: the window's homopolymers of 10 bases or more, in reference coordinates.
fn homopolymer_catalog() -> String {
    let window = fs::read_to_string(data_dir().join("reference.window.fa")).expect("reading the window");
    let mut lines = window.lines();
    let header = lines.next().expect("FASTA header");
    let start: usize =
        header.split(':').nth(1).and_then(|s| s.split('-').next()).and_then(|s| s.parse().ok()).expect("window start");
    let seq: Vec<u8> = lines.flat_map(str::bytes).collect();
    let mut catalog = String::new();
    let mut i = 0;
    while i < seq.len() {
        let run = seq[i..].iter().take_while(|&&b| b == seq[i]).count();
        if run >= 10 && seq[i] != b'N' {
            let (s, e) = (start - 1 + i, start - 1 + i + run);
            catalog += &format!("chr17\t{s}\t{e}\tID=h{s};MOTIFS={};STRUC=(x)n\n", seq[i] as char);
        }
        i += run;
    }
    assert!(catalog.lines().count() > 20, "repeat loci in the window: {catalog}");
    catalog
}

/// `lrpipe run --regions`'s options for the pass: `--repeats-regions` genotypes the catalog's
/// loci that a region reaches, as the pass without it does; `--bedmethyl-regions` writes the
/// rows of the regions, with the pass thresholds of their reads, as `lrpipe methyl pileup
/// --include-bed` (modkit's option) writes them, here from lrpipe-snv's ONT fixture, whose reads
/// have modification calls.
#[test]
fn regions_restrict_the_repeats_and_the_bedmethyl() {
    let dir = work_dir();
    let catalog = homopolymer_catalog();
    fs::write(dir.join("regions-catalog.bed"), &catalog).expect("writing the catalog");
    let loci: Vec<(i64, i64)> = catalog
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[1].parse().expect("a start"), f[2].parse().expect("an end"))
        })
        .collect();
    // From inside the 5th locus to inside the 15th, and the last base of the 20th.
    let (from, to) = ((loci[4].0 + loci[4].1) / 2, (loci[14].0 + loci[14].1) / 2);
    let bed = format!("chr17\t{from}\t{to}\nchr17\t{}\t{}\n", loci[19].1 - 1, loci[19].1);
    fs::write(dir.join("regions-pass.bed"), bed).expect("writing the regions");
    let pass = |out: &str, regions: Option<&str>| {
        let mut args = vec!["--input", "sample1.bam", "--reference", "reference.fa", "--threads", "4"];
        args.extend(["--repeats", "regions-catalog.bed", "--repeats-platform", "ont"]);
        let repeats = format!("{out}.repeats.vcf");
        args.extend(["--repeats-vcf", &repeats]);
        args.extend(regions.iter().flat_map(|r| ["--repeats-regions", r]));
        sniffles(&args);
    };
    pass("regions-all", None);
    pass("regions-in", Some("regions-pass.bed"));
    let read = |path: &Path| fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    let records =
        |text: String| -> Vec<String> { text.lines().filter(|l| !l.starts_with('#')).map(str::to_owned).collect() };
    let start = |l: &str| l.split('\t').nth(1).and_then(|p| p.parse::<i64>().ok()).expect("a position");
    // A locus's record is at its start, with the base before it.
    let reached: Vec<String> = records(read(&dir.join("regions-all.repeats.vcf")))
        .into_iter()
        .filter(|l| {
            let at = start(l);
            loci.iter().any(|&(s, e)| s == at && ((s < to && from < e) || (s < loci[19].1 && loci[19].1 - 1 < e)))
        })
        .collect();
    assert_eq!(reached.len(), 12, "{reached:?}");
    assert_eq!(records(read(&dir.join("regions-in.repeats.vcf"))), reached, "the loci the regions reach");

    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let (reads, reference) = (fixture.join("ont.bam"), fixture.join("chr20_window.fa.gz"));
    let regions = dir.join("regions-chr20.bed");
    fs::write(&regions, "chr20\t30195000\t30196000\nchr20\t30180000\t30190000\n").expect("writing the regions");
    let (pass, pileup) = (dir.join("regions-pass.bedmethyl"), dir.join("regions-pileup.bedmethyl"));
    let path = |p: &Path| p.to_str().expect("UTF-8").to_owned();
    sniffles(&[
        "--input",
        &path(&reads),
        "--reference",
        &path(&reference),
        "--bedmethyl",
        &path(&pass),
        "--bedmethyl-regions",
        &path(&regions),
    ]);
    let args = ["lrpipe-methyl", "pileup", &path(&reads), &path(&pileup), "--ref", &path(&reference)];
    let flags = ["--modified-bases", "5mC", "5hmC", "--cpg", "--combine-strands", "--include-bed", &path(&regions)];
    lrpipe_methyl::cli::main(args.iter().chain(&flags).map(Into::into).collect()).expect("lrpipe methyl pileup runs");
    let rows = read(&pass);
    assert!(rows.lines().count() > 100, "{rows}");
    let inside = |at: i64| (30_180_000..30_190_000).contains(&at) || (30_195_000..30_196_000).contains(&at);
    assert!(rows.lines().all(|l| inside(start(l))), "{rows}");
    assert_eq!(rows, read(&pileup), "the pileup of the regions");
}

/// A tabix index of the bgzipped VCF at `path`, as `tabix -p vcf` writes it, saved to `index`.
fn tabix_index(path: &Path, index: &Path) {
    use pasta::variant::index::tabix::{TabixBuilder, TabixConfig};
    use std::io::BufRead;
    let mut reader = pasta::bgzf::io::Reader::new(fs::File::open(path).expect("opening the VCF"));
    let mut builder = TabixBuilder::new(TabixConfig::VCF, 0);
    let mut line = Vec::new();
    while reader.read_until(b'\n', &mut line).expect("reading the VCF") > 0 {
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        builder.push_line(&line, u64::from(reader.virtual_position())).expect("indexing a line");
        line.clear();
    }
    let format = builder.format();
    builder.finish(u64::from(reader.virtual_position())).index.save(index, format).expect("saving the index");
}

/// The pass's bgzipped repeat VCF gets an index of its own: after `lrpipe run` genotyped a panel's
/// repeats in its separate stage (which indexes them), a whole-genome run's pass rewrote the VCF
/// in the same directory and left the panel's index, which found 1 of its records (lrpipe#42).
#[test]
fn the_repeat_vcf_gets_a_fresh_index_over_a_stale_one() {
    let dir = work_dir();
    let catalog = homopolymer_catalog();
    fs::write(dir.join("stale-catalog.bed"), &catalog).expect("writing the catalog");
    let first: Vec<&str> = catalog.lines().next().expect("a locus").split('\t').collect();
    fs::write(dir.join("stale-panel.bed"), format!("chr17\t{}\t{}\n", first[1], first[2])).expect("the panel");
    let vcf = dir.join("stale.repeats.vcf.gz");
    let index = dir.join("stale.repeats.vcf.gz.tbi");
    let pass = |regions: Option<&str>| {
        let mut args = vec!["--input", "sample1.bam", "--reference", "reference.fa", "--threads", "4"];
        args.extend(["--repeats", "stale-catalog.bed", "--repeats-platform", "ont"]);
        args.extend(["--repeats-vcf", "stale.repeats.vcf.gz"]);
        args.extend(regions.iter().flat_map(|r| ["--repeats-regions", r]));
        sniffles(&args);
    };
    // The panel's VCF and its index, as the separate repeat stage leaves them.
    pass(Some("stale-panel.bed"));
    tabix_index(&vcf, &index);
    let stale = fs::read(&index).expect("the panel's index");
    pass(None);
    let records = read_vcf(&vcf).lines().filter(|l| !l.starts_with('#')).count();
    assert!(records > 20, "repeat records: {records}");
    let fresh = dir.join("stale.fresh.tbi");
    tabix_index(&vcf, &fresh);
    let written = fs::read(&index).expect("the index");
    assert!(written != stale, "the panel's index was left beside the whole VCF");
    assert!(written == fs::read(&fresh).expect("the fresh index"), "the index is not the whole VCF's");
}

/// Reading a CRAM decodes only what the read pass reads (not the base qualities, and of the
/// tags `NM`, `SA`, `HP` and `PS`, with methylation's for `--bedmethyl`): the whole pass, SVs
/// phased with `HP` and `PS`, depth, repeats and bedMethyl, must write from `sample1.cram`
/// what it writes from the same reads as BAM.
#[test]
fn cram_gives_the_pass_bam_gives() {
    let dir = work_dir();
    fs::write(dir.join("homopolymers.bed"), homopolymer_catalog()).expect("writing the catalog");

    for (reads, out) in [("sample1.bam", "pass-bam"), ("sample1.cram", "pass-cram")] {
        sniffles(&[
            "--input",
            reads,
            "--reference",
            "reference.fa",
            "--threads",
            "4",
            "--phase",
            "--vcf",
            &format!("{out}.vcf"),
            "--depth",
            out,
            "--depth-by",
            "500",
            "--repeats",
            "homopolymers.bed",
            "--repeats-vcf",
            &format!("{out}.repeats.vcf"),
            "--repeats-platform",
            "ont",
            "--bedmethyl",
            &format!("{out}.bedmethyl"),
        ]);
    }
    let read = |name: &str| fs::read_to_string(dir.join(name)).unwrap_or_else(|e| panic!("reading {name}: {e}"));
    let (bam, cram) = (read("pass-bam.vcf"), read("pass-cram.vcf"));
    assert!(bam.lines().any(|l| !l.starts_with('#') && l.contains("|")), "phased SV calls: {bam}");
    let (bam_calls, cram_calls) = (comparable(&bam), comparable(&cram));
    if let Some((b, c)) = bam_calls.iter().zip(&cram_calls).find(|(b, c)| b != c) {
        panic!("the SV calls differ:\nBAM:  {b}\nCRAM: {c}");
    }
    assert_eq!(bam_calls.len(), cram_calls.len(), "SV calls");
    let repeats = read("pass-bam.repeats.vcf");
    assert!(repeats.lines().filter(|l| !l.starts_with('#')).count() > 20, "repeat records: {repeats}");
    assert_eq!(read("pass-cram.repeats.vcf"), repeats, "repeat genotypes");
    for suffix in ["mosdepth.summary.txt", "mosdepth.region.dist.txt"] {
        assert_eq!(read(&format!("pass-cram.{suffix}")), read(&format!("pass-bam.{suffix}")), "{suffix}");
    }
    assert_eq!(read("pass-cram.bedmethyl"), read("pass-bam.bedmethyl"), "bedMethyl");
}

/// A bgzipped reference (`.fa.gz` with its `.fai` and `.gzi`, as references are distributed, or
/// under a name that does not say so) gives what the plain one gives: Sniffles2's calls from the
/// BAM and the CRAM, forced genotypes,
/// and the whole read pass from the CRAM (phased SVs, depth, copy numbers, repeats, bedMethyl:
/// empty, as these reads have no methylation tags; `lrpipe run`'s tests read the ONT fixture's).
/// lrpipe-sv refused it before, and named `.fa` its repeats and copy numbers read it as text.
#[test]
fn a_bgzipped_reference_gives_what_the_plain_one_gives() {
    let dir = work_dir();
    sniffles(&["--input", "sample1.bam", "--vcf", "bgzf.reference.vcf", "--reference", "reference.fa.gz"]);
    assert_vcf_eq("reference.vcf", "bgzf.reference.vcf");
    sniffles(&["--input", "sample1.cram", "--vcf", "bgzf.cram.vcf", "--reference", "reference.fa.gz"]);
    assert_vcf_eq("cram.vcf", "bgzf.cram.vcf");

    let known = data_dir().join("expected/default.vcf");
    let catalog = dir.join("bgzf.catalog.bed");
    fs::write(
        &catalog,
        "chr17\t46250000\t46250030\tID=a;MOTIFS=CA;STRUC=(CA)n\nchr17\t46300000\t46300012\tID=b;MOTIFS=T;STRUC=(T)n\n",
    )
    .expect("writing the catalog");
    let runs = [("reference.fa", "plain"), ("reference.fa.gz", "bgzf"), ("reference.bgzf.fa", "named")];
    for (reference, out) in runs {
        sniffles(&[
            "--input",
            "sample2.cram",
            "--vcf",
            &format!("{out}.genotype.vcf"),
            "--reference",
            reference,
            "--genotype-vcf",
            known.to_str().expect("UTF-8 path"),
        ]);
        sniffles(&[
            "--input",
            "sample1.cram",
            "--reference",
            reference,
            "--threads",
            "4",
            "--phase",
            "--vcf",
            &format!("{out}.pass.vcf"),
            "--depth",
            &format!("{out}.pass"),
            "--depth-by",
            "1000",
            "--depth-fast",
            "--depth-mapq",
            "20",
            "--cnv",
            &format!("{out}.pass.cnv"),
            "--repeats",
            "bgzf.catalog.bed",
            "--repeats-vcf",
            &format!("{out}.pass.repeats.vcf"),
            "--repeats-platform",
            "ont",
            "--bedmethyl",
            &format!("{out}.pass.bedmethyl"),
        ]);
    }
    let read = |name: &str| fs::read(dir.join(name)).unwrap_or_else(|e| panic!("reading {name}: {e}"));
    for name in ["genotype.vcf", "pass.vcf"] {
        let plain = read_vcf(&dir.join(format!("plain.{name}")));
        assert!(plain.lines().any(|l| !l.starts_with('#')), "calls in {name}");
        for bgzf in ["bgzf", "named"] {
            let other = read_vcf(&dir.join(format!("{bgzf}.{name}")));
            assert_eq!(comparable(&other), comparable(&plain), "{bgzf}.{name}");
        }
    }
    let repeats = String::from_utf8(read("plain.pass.repeats.vcf")).expect("UTF-8");
    assert_eq!(repeats.lines().filter(|l| !l.starts_with('#')).count(), 2, "repeat records: {repeats}");
    let mut compared = vec![
        "pass.repeats.vcf".to_owned(),
        "pass.bedmethyl".to_owned(),
        "pass.mosdepth.summary.txt".to_owned(),
        "pass.regions.bed.gz".to_owned(),
        "pass.mosdepth.region.dist.txt".to_owned(),
    ];
    let mut cnv: Vec<String> = fs::read_dir(dir.join("plain.pass.cnv"))
        .expect("the CNV directory")
        .map(|e| format!("pass.cnv/{}", e.expect("an entry").file_name().to_string_lossy()))
        // lrpipe-cnv's note of the reference its metadata was checked against names that file.
        .filter(|f| !f.ends_with(".mdr.checked"))
        .collect();
    cnv.sort();
    assert!(cnv.iter().any(|f| f.ends_with("_cnv.bed.gz")), "CNV calls: {cnv:?}");
    compared.extend(cnv);
    for name in compared {
        for bgzf in ["bgzf", "named"] {
            assert!(read(&format!("{bgzf}.{name}")) == read(&format!("plain.{name}")), "{bgzf}.{name} differs");
        }
    }
}

/// `bai` as a CSI (BGZF-compressed, `min_shift` 14, depth 5: the BAI's bins), with every bin's
/// linear offset 0, which only makes queries read from each bin's first chunk.
fn csi_from_bai(bai: &[u8]) -> Vec<u8> {
    let int = |o: usize| u32::from_le_bytes(bai[o..o + 4].try_into().expect("four bytes")) as usize;
    let mut csi = b"CSI\x01".to_vec();
    for v in [14u32, 5, 0] {
        csi.extend(v.to_le_bytes());
    }
    csi.extend_from_slice(&bai[4..8]);
    let mut o = 8;
    for _ in 0..int(4) {
        csi.extend_from_slice(&bai[o..o + 4]);
        let bins = int(o);
        o += 4;
        for _ in 0..bins {
            let chunks = int(o + 4);
            csi.extend_from_slice(&bai[o..o + 4]);
            csi.extend(0u64.to_le_bytes());
            csi.extend_from_slice(&bai[o + 4..o + 8 + 16 * chunks]);
            o += 8 + 16 * chunks;
        }
        o += 4 + 8 * int(o);
    }
    csi.extend_from_slice(&bai[o..]);
    bgzf_block(&csi).into_iter().chain(bgzf_block(&[])).collect()
}

/// One BGZF block holding `data` (under 64 KiB).
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut deflate = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
    deflate.write_all(data).expect("compressing");
    let body = deflate.finish().expect("compressed");
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let size = u16::try_from(18 + body.len() + 8 - 1).expect("one block");
    let mut out = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0];
    out.extend(size.to_le_bytes());
    out.extend(body);
    out.extend(crc.sum().to_le_bytes());
    out.extend(u32::try_from(data.len()).expect("small").to_le_bytes());
    out
}

#[test]
fn a_csi_index_reads_as_its_bai_does() {
    // LRP-009: the BAM with only a CSI gives the calls it gives with its BAI.
    let dir = work_dir().join("csi");
    fs::create_dir_all(&dir).expect("csi directory");
    fs::copy(data_dir().join("sample1.bam"), dir.join("sample1.bam")).expect("BAM copied");
    let bai = fs::read(data_dir().join("sample1.bam.bai")).expect("BAI read");
    fs::write(dir.join("sample1.bam.csi"), csi_from_bai(&bai)).expect("CSI written");
    sniffles(&["--input", "csi/sample1.bam", "--vcf", "csi.vcf"]);
    assert_vcf_eq("default.vcf", "csi.vcf");
}

/// Forced genotyping of a breakend takes the support of the called breakend only: the same mate
/// contig, the mate within `--cluster-merge-bnd` (1000) and the same orientation on both sides
/// (LRP-015). sample2 calls `N]chr17:45585676]` at chr17:46292138 with 8 reads.
#[test]
fn genotype_vcf_breakend_identity() {
    let header: String = fs::read_to_string(data_dir().join("expected/default.vcf"))
        .expect("reading default.vcf")
        .lines()
        .filter(|l| l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    let targets = [
        ("same", "N]chr17:45585676]", true),
        ("reference_base", "C]chr17:45585676]", true),
        ("mate_at_tolerance", "N]chr17:45586676]", true),
        ("mate_past_tolerance", "N]chr17:45586677]", false),
        ("mate_far", "N]chr17:46485676]", false),
        ("mate_other_contig", "N]chr1:45585676]", false),
        ("remote_side_opposite", "N[chr17:45585676[", false),
        ("local_side_opposite", "]chr17:45585676]N", false),
        ("both_opposite", "[chr17:45585676[N", false),
    ];
    let mut vcf = header;
    for (id, alt, _) in targets {
        vcf.push_str(&format!("chr17\t46292138\t{id}\tN\t{alt}\t.\tPASS\tSVTYPE=BND\n"));
    }
    fs::write(work_dir().join("breakends.vcf"), vcf).expect("writing the targets");
    sniffles(&[
        "--input",
        "sample2.cram",
        "--vcf",
        "breakends.genotyped.vcf",
        "--reference",
        "reference.fa",
        "--genotype-vcf",
        "breakends.vcf",
    ]);
    let text = read_vcf(&work_dir().join("breakends.genotyped.vcf"));
    let records: Vec<Vec<&str>> =
        text.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect()).collect();
    assert_eq!(records.len(), targets.len(), "one record per target");
    for ((id, _, matched), record) in targets.iter().zip(&records) {
        assert_eq!(record[2], *id, "targets in input order");
        if *matched {
            assert_eq!(record[9], "0/1:37:6:8:.", "{id} takes the breakend's 8 reads");
        } else {
            let sample: Vec<&str> = record[9].split(':').collect();
            assert_eq!((sample[0], sample[1], sample[3]), ("0/0", "0", "0"), "{id} takes no support: {}", record[9]);
        }
    }
}

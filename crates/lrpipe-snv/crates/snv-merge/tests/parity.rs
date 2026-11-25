//! Compares `sort_vcf` with the original's own `SortVcf.py` outputs.
//!
//! Opt-in, as the oracle files are large:
//!
//! - `SNV_M3_ORACLE`: a directory with `pileup_output/pileup_<contig>_<chunk>.vcf`,
//!   `CONTIGS`, `CMD`, and the `pileup.vcf.gz` (+ `.tbi`) SortVcf made from them with the
//!   pipeline's default flags. `SNV_M3_REF` optionally names the `--ref_fn` it was run with.
//! - `SNV_MERGE_CASES`: a directory with `cases.tsv`, `inputs/` and `py/<case>/` (the original's
//!   outputs and exit `status`), as `tests/oracle/make_inputs.sh` and `run_python.sh` build
//!   them. Each line of `cases.tsv` is `name, input dir, contigs file, cmd file or None, REF or
//!   -, extra SortVcf flags`.
//!
//! Outputs go to `SNV_MERGE_OUT` (default: cargo's test temp dir). Decompressed VCFs must be
//! identical; whether the `.gz` and `.tbi` bytes are identical too is printed.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use pasta::core::{Position, Region};
use pasta::csi::{self, BinningIndex as _};
use pasta::tabix;
use snv_merge::{
    ChunkLayout, Gender, Haploid, IndexStatus, MalformedRecords, Postprocess, RefCalls, SortOptions, sort_vcf,
};

fn gunzip(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    MultiGzDecoder::new(file).read_to_end(&mut out).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    out
}

fn same_bytes(a: &Path, b: &Path) -> bool {
    matches!((fs::read(a), fs::read(b)), (Ok(x), Ok(y)) if x == y)
}

fn out_dir(name: &str) -> PathBuf {
    let base =
        std::env::var_os("SNV_MERGE_OUT").map_or_else(|| PathBuf::from(env!("CARGO_TARGET_TMPDIR")), PathBuf::from);
    let dir = base.join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("output directory can be created");
    dir
}

/// Every record tabix returns for `contig:start-end` (0-based, half-open).
fn query(gz: &Path, contig: &str, start: u64, end: u64) -> Vec<Vec<u8>> {
    let index = tabix::fs::read(gz.with_added_extension("tbi")).unwrap_or_else(|e| panic!("{}.tbi: {e}", gz.display()));
    let known = index.header().is_some_and(|h| h.reference_sequence_names().contains(contig.as_bytes()));
    if !known {
        return Vec::new();
    }
    let file = fs::File::open(gz).unwrap_or_else(|e| panic!("{}: {e}", gz.display()));
    let mut reader = csi::io::IndexedReader::new(file, index);
    let position = |p: u64| Position::new(usize::try_from(p).expect("position fits")).expect("positive");
    let region = Region::new(contig, position(start + 1)..=position(end));
    let records = reader.query(&region).expect("query succeeds");
    records.map(|r| AsRef::<str>::as_ref(&r.expect("record reads")).as_bytes().to_vec()).collect()
}

/// Checks that tabix queries over both indexed files return the same records, over windows
/// spanning each contig's records. Returns the number of records compared.
fn compare_queries(text: &[u8], got: &Path, want: &Path) -> usize {
    let mut spans: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    for line in text.split(|&b| b == b'\n').filter(|l| !l.is_empty() && l[0] != b'#') {
        let line = String::from_utf8_lossy(line);
        let mut cols = line.split('\t');
        let (Some(contig), Some(Ok(pos))) = (cols.next(), cols.next().map(str::parse::<u64>)) else { continue };
        let span = spans.entry(contig.to_owned()).or_insert((pos, pos));
        *span = (span.0.min(pos), span.1.max(pos));
    }
    let mut compared = 0;
    for (contig, (first, last)) in &spans {
        let step = ((last - first) / 16).max(1);
        let mut start = first.saturating_sub(1);
        while start <= *last {
            let end = start + step / 2 + 1;
            let (a, b) = (query(got, contig, start, end), query(want, contig, start, end));
            assert_eq!(a, b, "tabix {contig}:{start}-{end}");
            compared += a.len();
            start += step;
        }
        let whole = query(got, contig, 0, last + 1);
        assert_eq!(whole, query(want, contig, 0, last + 1), "tabix {contig}");
        compared += whole.len();
    }
    compared
}

fn str2bool(value: &str) -> bool {
    match value.to_ascii_lowercase().as_str() {
        "yes" | "true" | "t" | "y" | "1" => true,
        "no" | "false" | "f" | "n" | "0" => false,
        _ => panic!("not a boolean: {value:?}"),
    }
}

/// `SortOptions` for the pipeline's pileup SortVcf call plus `extra` flags.
fn pipeline_options(input: &Path, output: &Path, contigs: &Path, extra: &str, inputs: &Path) -> SortOptions {
    let mut options = SortOptions::new(input, output);
    options.prefix = Some("pileup".to_owned());
    options.contigs = Some(contigs.to_owned());
    let (mut precise, mut sensitive) = (false, false);
    let flags: Vec<&str> = extra.split_whitespace().collect();
    for pair in flags.chunks(2) {
        let [flag, value] = pair else { panic!("flag without value in {extra:?}") };
        match *flag {
            "--pileup_only" => {
                options.postprocess = if str2bool(value) { Postprocess::PileupOnly } else { Postprocess::Off };
            }
            "--print_ref_calls" => options.ref_calls = if str2bool(value) { RefCalls::Show } else { RefCalls::Hide },
            "--haploid_precise" => precise = str2bool(value),
            "--haploid_sensitive" => sensitive = str2bool(value),
            "--gender" => {
                options.gender = match *value {
                    "male" => Gender::Male,
                    "female" => Gender::Female,
                    _ => Gender::Unknown,
                };
            }
            "--par_regions_bed" => {
                options.par_regions = (*value != "EMPTY").then(|| inputs.join(value));
            }
            "--qual" => options.qual = value.parse().expect("--qual is an integer"),
            "--use_gpu" => {
                options.layout = if str2bool(value) { ChunkLayout::Unordered } else { ChunkLayout::ByContig };
            }
            "--check_malformed_records" => {
                options.malformed = if str2bool(value) { MalformedRecords::Repair } else { MalformedRecords::Keep };
            }
            _ => panic!("unknown flag {flag}"),
        }
    }
    options.haploid = match (precise, sensitive) {
        (true, _) => Haploid::Precise,
        (false, true) => Haploid::Sensitive,
        (false, false) => Haploid::Off,
    };
    options
}

#[test]
fn m3_oracle_matches_the_original() {
    let Some(oracle) = std::env::var_os("SNV_M3_ORACLE").map(PathBuf::from) else {
        eprintln!("skipped: set SNV_M3_ORACLE to the m3 oracle directory");
        return;
    };
    let out = out_dir("m3");
    let mut options =
        pipeline_options(&oracle.join("pileup_output"), &out.join("pileup.vcf"), &oracle.join("CONTIGS"), "", &oracle);
    options.cmd_file = Some(oracle.join("CMD"));
    options.reference = std::env::var_os("SNV_M3_REF").map(PathBuf::from);
    let summary = sort_vcf(&options).expect("sort_vcf succeeds");
    assert_eq!(summary.index, IndexStatus::Built);
    assert!(!out.join("pileup.vcf").exists(), "the uncompressed output is removed");

    let want_gz = oracle.join("pileup.vcf.gz");
    let text = gunzip(&summary.vcf);
    assert!(text == gunzip(&want_gz), "decompressed VCF differs from the oracle");
    let compared = compare_queries(&text, &summary.vcf, &want_gz);
    assert!(compared > 0, "no records compared through tabix");
    eprintln!(
        "m3 oracle: {} records, VCF identical; tabix queries agree ({compared} records); .gz bytes identical: {}; .tbi bytes identical: {}",
        summary.records,
        same_bytes(&summary.vcf, &want_gz),
        same_bytes(&PathBuf::from(format!("{}.tbi", summary.vcf.display())), &oracle.join("pileup.vcf.gz.tbi")),
    );
}

#[test]
fn extra_cases_match_the_original() {
    let Some(cases) = std::env::var_os("SNV_MERGE_CASES").map(PathBuf::from) else {
        eprintln!("skipped: set SNV_MERGE_CASES to a directory with cases.tsv, inputs/ and py/");
        return;
    };
    let inputs = cases.join("inputs");
    let reference = fs::read_to_string(inputs.join("REF")).expect("inputs/REF names the reference");
    let table = fs::read_to_string(cases.join("cases.tsv")).expect("cases.tsv is readable");
    let mut checked = 0;
    for line in table.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
        let cols: Vec<&str> = line.split('\t').collect();
        let [name, input, contigs, cmd, reference_token, extra] = cols[..] else { panic!("bad case line {line:?}") };
        let py = cases.join("py").join(name);
        let out = out_dir(&format!("cases/{name}"));
        let mut options =
            pipeline_options(&inputs.join(input), &out.join("pileup.vcf"), &inputs.join(contigs), extra, &inputs);
        options.cmd_file = (cmd != "None").then(|| inputs.join(cmd));
        options.reference = (reference_token == "REF").then(|| PathBuf::from(reference.trim()));
        let status = fs::read_to_string(py.join("status")).expect("py/<case>/status exists");
        let result = sort_vcf(&options);
        if status.trim() != "0" {
            assert!(result.is_err(), "{name}: the original failed but sort_vcf returned {result:?}");
            eprintln!("{name}: both fail ({})", result.expect_err("checked above"));
            checked += 1;
            continue;
        }
        let summary = result.unwrap_or_else(|e| panic!("{name}: {e}"));
        let want_gz = py.join("pileup.vcf.gz");
        let text = gunzip(&summary.vcf);
        assert!(text == gunzip(&want_gz), "{name}: decompressed VCF differs");
        let want_tbi = py.join("pileup.vcf.gz.tbi");
        let got_tbi = PathBuf::from(format!("{}.tbi", summary.vcf.display()));
        assert_eq!(got_tbi.exists(), want_tbi.exists(), "{name}: index presence differs");
        assert_eq!(summary.index == IndexStatus::Built, want_tbi.exists(), "{name}: {:?}", summary.index);
        let compared = if want_tbi.exists() { compare_queries(&text, &summary.vcf, &want_gz) } else { 0 };
        let want_copy = py.join("merge_output.vcf.gz");
        assert_eq!(summary.copied_to.is_some(), want_copy.exists(), "{name}: merge_output copy differs");
        if let Some(copy) = &summary.copied_to {
            assert!(gunzip(copy) == gunzip(&want_copy), "{name}: merge_output.vcf.gz differs");
        }
        eprintln!(
            "{name}: {} rows, VCF identical, tabix agrees ({compared} records); .gz identical: {}; .tbi identical: {}",
            summary.records,
            same_bytes(&summary.vcf, &want_gz),
            if want_tbi.exists() {
                same_bytes(&got_tbi, &want_tbi).to_string()
            } else {
                "n/a (neither indexed)".to_owned()
            },
        );
        checked += 1;
    }
    assert!(checked > 0, "no cases in cases.tsv");
}

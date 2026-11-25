//! Compares `merge_vcf` with the original's own `MergeVcf.py` outputs.
//!
//! Opt-in, as the oracle files are large:
//!
//! - `SNV_M4_ORACLE`: the original's output directory (`--remove_intermediate_dir=False`) with
//!   `pileup.vcf.gz`, `full_alignment.vcf.gz`, `tmp/merge_output/merge_<contig>.vcf`,
//!   `tmp/CONTIGS`, `tmp/CMD` and `merge_output.vcf.gz` (+ `.tbi`). `SNV_M4_CONTIG` names
//!   the contig (default `chr20`) and `SNV_M4_REF` the `--ref_fn` of the final SortVcf.
//!   The merged VCF must be byte-identical, then `sort_vcf` with the pipeline's final-sort
//!   options must give the same decompressed `merge_output.vcf.gz`.
//! - `SNV_MERGEVCF_CASES`: a directory with `cases.tsv`, `inputs/` and `py/<case>/` (the original's
//!   `merge.vcf`, exit `status` and copied BED directories), as
//!   `tests/oracle/merge_vcf_inputs.sh` and `merge_vcf_python.sh` build them.
//!
//! Outputs go to `SNV_MERGE_OUT` (default: cargo's test temp dir).

use std::fs;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use snv_merge::{Gender, Haploid, IndexStatus, MergeMode, MergeOptions, RefCalls, SortOptions, merge_vcf, sort_vcf};

fn gunzip(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    MultiGzDecoder::new(file).read_to_end(&mut out).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    out
}

fn read(path: &Path) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
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

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).expect("directory can be created");
    for entry in fs::read_dir(from).expect("directory is readable") {
        let entry = entry.expect("entry is readable");
        fs::copy(entry.path(), to.join(entry.file_name())).expect("file can be copied");
    }
}

fn str2bool(value: &str) -> bool {
    match value.to_ascii_lowercase().as_str() {
        "yes" | "true" | "t" | "y" | "1" => true,
        "no" | "false" | "f" | "n" | "0" => false,
        _ => panic!("not a boolean: {value:?}"),
    }
}

fn first_difference(got: &[u8], want: &[u8]) -> String {
    let line = got.iter().zip(want).take_while(|(a, b)| a == b).filter(|(a, _)| **a == b'\n').count() + 1;
    format!("first difference on line {line} ({} vs {} bytes)", got.len(), want.len())
}

/// `MergeOptions` from MergeVcf flags. `--bed_fn_prefix` and `--par_regions_bed` are relative
/// to `beds` and `inputs`; `''` is an empty value.
fn options_from(flags: &str, inputs: &Path, beds: &Path, options: &mut MergeOptions) {
    let (mut precise, mut sensitive, mut platform, mut bed_dir) = (false, false, "ont".to_owned(), None);
    let flags: Vec<&str> = flags.split_whitespace().map(|f| if f == "''" { "" } else { f }).collect();
    for pair in flags.chunks(2) {
        let [flag, value] = *pair else { panic!("flag without value in {flags:?}") };
        match flag {
            "--ctgName" => options.contig = Some(value.to_owned()),
            "--platform" => value.clone_into(&mut platform),
            "--bed_fn_prefix" => bed_dir = Some(beds.join(value)),
            "--print_ref_calls" => options.ref_calls = if str2bool(value) { RefCalls::Show } else { RefCalls::Hide },
            "--haploid_precise" => precise = str2bool(value),
            "--haploid_sensitive" => sensitive = str2bool(value),
            "--gender" => {
                options.gender = match value {
                    "male" => Gender::Male,
                    "female" => Gender::Female,
                    _ => Gender::Unknown,
                };
            }
            "--par_regions_bed" => options.par_regions = (value != "EMPTY").then(|| inputs.join(value)),
            "--qual" => options.qual = value.parse().expect("--qual is an integer"),
            _ => panic!("unknown flag {flag}"),
        }
    }
    options.haploid = match (precise, sensitive) {
        (true, _) => Haploid::Precise,
        (false, true) => Haploid::Sensitive,
        (false, false) => Haploid::Off,
    };
    if platform == "ilmn" {
        options.mode = MergeMode::Regions { candidate_beds: bed_dir.expect("ilmn cases name --bed_fn_prefix") };
    }
}

#[test]
fn m4_oracle_matches_the_original() {
    let Some(oracle) = std::env::var_os("SNV_M4_ORACLE").map(PathBuf::from) else {
        eprintln!("skipped: set SNV_M4_ORACLE to the original's output directory with its tmp/");
        return;
    };
    let contig = std::env::var("SNV_M4_CONTIG").unwrap_or_else(|_| "chr20".to_owned());
    let out = out_dir("m4");
    let merge_dir = out.join("merge_output");
    fs::create_dir_all(&merge_dir).expect("merge_output can be created");
    let name = format!("merge_{contig}.vcf");
    // The pipeline's MergeVcf call with its default flags.
    let mut options =
        MergeOptions::new(oracle.join("pileup.vcf.gz"), oracle.join("full_alignment.vcf.gz"), merge_dir.join(&name));
    options.contig = Some(contig.clone());
    let summary = merge_vcf(&options).expect("merge_vcf succeeds");
    let (got, want) = (read(&merge_dir.join(&name)), read(&oracle.join("tmp/merge_output").join(&name)));
    assert!(got == want, "{name} differs from the oracle: {}", first_difference(&got, &want));
    eprintln!(
        "m4 oracle: {name} byte-identical ({} records; {} from pileup, {} from full alignment)",
        summary.records, summary.pileup, summary.full_alignment
    );

    // The pipeline's final SortVcf call.
    let mut sort = SortOptions::new(&merge_dir, out.join("merge_output.vcf"));
    sort.prefix = Some("merge".to_owned());
    sort.reference = std::env::var_os("SNV_M4_REF").map(PathBuf::from);
    sort.cmd_file = Some(oracle.join("tmp/CMD"));
    sort.contigs = Some(oracle.join("tmp/CONTIGS"));
    let sorted = sort_vcf(&sort).expect("sort_vcf succeeds");
    assert_eq!(sorted.index, IndexStatus::Built);
    let want_gz = oracle.join("merge_output.vcf.gz");
    let (got, want) = (gunzip(&sorted.vcf), gunzip(&want_gz));
    assert!(got == want, "merge_output.vcf.gz differs: {}", first_difference(&got, &want));
    eprintln!(
        "m4 oracle: merge_output.vcf.gz decompressed identical ({} records); .gz bytes identical: {}; .tbi bytes identical: {}",
        sorted.records,
        same_bytes(&sorted.vcf, &want_gz),
        same_bytes(&PathBuf::from(format!("{}.tbi", sorted.vcf.display())), &oracle.join("merge_output.vcf.gz.tbi")),
    );
}

#[test]
fn merge_vcf_cases_match_the_original() {
    let Some(cases) = std::env::var_os("SNV_MERGEVCF_CASES").map(PathBuf::from) else {
        eprintln!("skipped: set SNV_MERGEVCF_CASES to a directory with cases.tsv, inputs/ and py/");
        return;
    };
    let inputs = cases.join("inputs");
    let table = fs::read_to_string(cases.join("cases.tsv")).expect("cases.tsv is readable");
    let mut checked = 0;
    for line in table.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
        let cols: Vec<&str> = line.split('\t').collect();
        let [name, pileup, full_alignment, flags] = cols[..] else { panic!("bad case line {line:?}") };
        let py = cases.join("py").join(name);
        let out = out_dir(&format!("merge_vcf_cases/{name}"));
        let mut options = MergeOptions::new(inputs.join(pileup), inputs.join(full_alignment), out.join("merge.vcf"));
        options_from(flags, &inputs, &out, &mut options);
        if let MergeMode::Regions { candidate_beds } = &options.mode {
            let source = inputs.join(candidate_beds.file_name().expect("a directory name"));
            if source.is_dir() {
                copy_dir(&source, candidate_beds);
            }
        }
        let status = fs::read_to_string(py.join("status")).expect("py/<case>/status exists");
        let result = merge_vcf(&options);
        if status.trim() != "0" {
            let error = result.expect_err(&format!("{name}: the original failed but merge_vcf succeeded"));
            eprintln!("{name}: both fail ({error})");
            checked += 1;
            continue;
        }
        let summary = result.unwrap_or_else(|e| panic!("{name}: {e}"));
        let (got, want) = (read(&out.join("merge.vcf")), read(&py.join("merge.vcf")));
        assert!(got == want, "{name}: merge.vcf differs: {}", first_difference(&got, &want));
        let mut regions = String::new();
        if let (MergeMode::Regions { candidate_beds }, Some(contig)) = (&options.mode, &options.contig) {
            let file = format!("full_aln_regions_{contig}");
            let want = py.join(candidate_beds.file_name().expect("a directory name")).join(&file);
            assert_eq!(summary.candidate_bed.is_some(), want.exists(), "{name}: {file} presence differs");
            if let Some(got) = &summary.candidate_bed {
                assert!(read(got) == read(&want), "{name}: {file} differs");
                regions = format!(", {file} identical");
            }
        }
        eprintln!("{name}: {} records, merge.vcf identical{regions}", summary.records);
        checked += 1;
    }
    assert!(checked > 0, "no cases in cases.tsv");
}

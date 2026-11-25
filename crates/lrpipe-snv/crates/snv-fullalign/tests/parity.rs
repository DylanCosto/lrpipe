//! Compares `FullAlignmentJob::run` with the original's own `CreateTensorFullAlignment`.
//!
//! Opt-in, as the oracle files are large: set `SNV_FA_CASES` to a `cases.tsv` whose lines
//! are `name, bam, reference, contig, full_aln_regions, phased VCF or -, min MQ, platform
//! (ont/hifi/ilmn), long indels (0/1)`. Next to it, `NAME.tensors.npy`, `NAME.info.tsv` and
//! optionally `NAME.seconds` are what `tools/oracle_fullalign.py` wrote for that job.
//! `SNV_FA_ONLY` restricts the run to a comma-separated list of names. Tensors must be
//! bit-identical and the strings identical; the time of each side is printed (run with
//! `--release --nocapture` for meaningful numbers).

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use ndarray::Array4;
use snv_fullalign::{FeatureOptions, FullAlignment, FullAlignmentJob, IndelLimit, Platform};

struct Case {
    name: String,
    job: FullAlignmentJob,
}

fn parse_case(line: &str) -> Case {
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!(f.len(), 9, "case line needs 9 fields: {line}");
    let mut options = FeatureOptions::new();
    options.min_mapq = f[6].parse().expect("min MQ");
    options.platform = match f[7] {
        "ont" => Platform::Ont,
        "hifi" => Platform::Hifi,
        "ilmn" => Platform::Ilmn,
        other => panic!("unknown platform {other}"),
    };
    options.indel_limit = if f[8] == "1" { IndelLimit::Long } else { IndelLimit::Standard };
    let mut job = FullAlignmentJob::new(f[1], f[2], f[3], f[4]).options(options);
    if f[5] != "-" {
        job = job.phased_vcf(f[5]);
    }
    Case { name: f[0].to_owned(), job }
}

/// The first difference between two equal-shape tensors, as (candidate, row, position, channel).
fn first_diff(got: &[i8], want: &[i8], shape: [usize; 4]) -> Option<String> {
    let i = got.iter().zip(want).position(|(a, b)| a != b)?;
    let [_, d, p, c] = shape;
    let (cand, rest) = (i / (d * p * c), i % (d * p * c));
    let (row, rest) = (rest / (p * c), rest % (p * c));
    let n = got.iter().zip(want).filter(|(a, b)| a != b).count();
    Some(format!(
        "{n} values differ; first at candidate {cand} row {row} position {} channel {}: {} vs {}",
        rest / c,
        rest % c,
        got[i],
        want[i]
    ))
}

fn compare(name: &str, dir: &Path, got: &FullAlignment) -> Vec<String> {
    let mut problems = Vec::new();
    let want: Array4<i8> =
        ndarray_npy::read_npy(dir.join(format!("{name}.tensors.npy"))).expect("oracle tensors are readable");
    let shape = got.shape();
    if want.shape() != shape {
        problems.push(format!("shape {:?} vs oracle {:?}", shape, want.shape()));
    } else {
        let want = want.as_standard_layout();
        let want = want.as_slice().expect("standard layout");
        if let Some(diff) = first_diff(got.tensors(), want, shape) {
            problems.push(diff);
        }
    }
    let info = fs::read_to_string(dir.join(format!("{name}.info.tsv"))).expect("oracle info is readable");
    let want_lines: Vec<&str> = info.lines().collect();
    let got_lines: Vec<String> = got.positions().iter().zip(got.alt_info()).map(|(p, a)| format!("{p}\t{a}")).collect();
    if got_lines.len() != want_lines.len() {
        problems.push(format!("{} strings vs oracle {}", got_lines.len(), want_lines.len()));
    }
    if let Some(i) = got_lines.iter().zip(&want_lines).position(|(g, w)| g != w) {
        problems.push(format!("string {i}: {:?} vs oracle {:?}", got_lines[i], want_lines[i]));
    }
    problems
}

#[test]
fn matches_the_original() {
    let Some(cases) = std::env::var_os("SNV_FA_CASES").map(PathBuf::from) else {
        eprintln!("skipping: set SNV_FA_CASES to a cases.tsv with oracle outputs to compare");
        return;
    };
    let only: Option<Vec<String>> =
        std::env::var("SNV_FA_ONLY").ok().map(|s| s.split(',').map(str::to_owned).collect());
    let dir = cases.parent().expect("cases.tsv has a directory").to_owned();
    let text = fs::read_to_string(&cases).expect("cases.tsv is readable");
    let mut failures = Vec::new();
    let mut compared = 0;
    for line in text.lines().filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let case = parse_case(line);
        if only.as_ref().is_some_and(|names| !names.contains(&case.name)) {
            continue;
        }
        let start = Instant::now();
        let got = case.job.run().unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        let seconds = start.elapsed().as_secs_f64();
        let oracle_seconds = fs::read_to_string(dir.join(format!("{}.seconds", case.name)))
            .ok()
            .and_then(|s| s.trim().parse::<f64>().ok());
        let problems = compare(&case.name, &dir, &got);
        let status = if problems.is_empty() { "identical" } else { "DIFFERS" };
        let oracle = oracle_seconds.map_or_else(|| "-".to_owned(), |s| format!("{s:.3} s"));
        eprintln!("{}: {} candidates, {status}; rust {seconds:.3} s, original {oracle}", case.name, got.len());
        for p in &problems {
            eprintln!("    {p}");
        }
        if !problems.is_empty() {
            failures.push(case.name.clone());
        }
        compared += 1;
    }
    assert!(compared > 0, "no case was compared");
    assert!(failures.is_empty(), "cases differ from the original: {failures:?}");
}

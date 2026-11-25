//! Small hand-made VCFs against outputs of the original's own scripts (`tools/oracle_select.sh`).
//!
//! Each `tests/data/<case>` has `input.vcf`, `runs.txt` and `expected/<n>-<kind>/` with the
//! files, stdout and stderr of run `n`. Every case is run from the plain file and from a
//! two-member gzip copy, as `pigz -fdc` accepts both. In `FULL_ALN_FILE_*` the split folder
//! is written as `@SPLIT_FOLDER@`.

use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;
use snv_select::{
    CandidateOptions, CandidateVcf, CandidateWindow, FULL_ALN_FILE_PREFIX, FULL_QUAL_FILE, FullQualOptions,
    HetSnpOptions, PHASE_QUAL_FILE, open_vcf, select_candidates_for_full_alignment, select_full_qual,
    select_het_snps_for_phasing, select_phase_qual,
};

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// The lines as Python's `print` writes them.
fn printed<S: AsRef<str>>(lines: &[S]) -> String {
    lines.iter().map(|l| format!("{}\n", l.as_ref())).collect()
}

/// `input` as two concatenated gzip members, split mid-file like BGZF blocks.
fn gzip_copy(input: &Path, out: &Path) {
    let text = fs::read(input).expect("read input");
    let (a, b) = text.split_at(text.len() / 2);
    let mut file = Vec::new();
    for part in [a, b] {
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(part).expect("compress");
        file.extend(gz.finish().expect("compress"));
    }
    fs::write(out, file).expect("write gzip copy");
}

/// Replaces the split folder in the `FULL_ALN_FILE_*` lists in `dir` with `@SPLIT_FOLDER@`, as
/// `tools/oracle_select.sh` does.
fn normalize_full_aln_lists(dir: &Path) {
    let prefix = format!("{}/", dir.to_str().expect("UTF-8 path"));
    for entry in fs::read_dir(dir).expect("output dir") {
        let path = entry.expect("entry").path();
        let is_list = path.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.starts_with(FULL_ALN_FILE_PREFIX));
        if is_list {
            let text = fs::read_to_string(&path).expect("read list");
            let normalized: String = text
                .lines()
                .map(|line| {
                    format!(
                        "{}\n",
                        line.strip_prefix(&prefix).map_or(line.to_owned(), |n| format!("@SPLIT_FOLDER@/{n}"))
                    )
                })
                .collect();
            fs::write(&path, normalized).expect("write list");
        }
    }
}

/// Runs every line of `runs.txt` on `input`, writing run `n` to `work/<n>-<kind>`.
fn run_case(case: &Path, input: &Path, work: &Path) {
    let runs = fs::read_to_string(case.join("runs.txt")).expect("runs.txt");
    let mut phase_qual: Option<PathBuf> = None;
    let mut full_qual: Option<PathBuf> = None;
    for (i, run) in runs.lines().enumerate() {
        let args: Vec<&str> = run.split_whitespace().collect();
        let out = work.join(format!("{}-{}", i + 1, args[0]));
        fs::create_dir_all(&out).expect("create run dir");
        let pct = |k: usize| -> f64 { args[k].parse().expect("proportion in runs.txt") };
        let (stdout, stderr) = match args[0] {
            "phase" => {
                let result = select_phase_qual(open_vcf(input).expect("open"), pct(1)).expect("phase");
                phase_qual = Some(result.write_in(&out).expect("write phase_qual"));
                (printed(&[result.info_message()]), printed(result.warnings()))
            }
            "full" => {
                let vcf_fn = case.join(args[3]);
                let options = FullQualOptions {
                    var_pct_full: pct(1),
                    ref_pct_full: pct(2),
                    candidate_vcf: CandidateVcf::from_arg(vcf_fn.to_str().expect("UTF-8 path")),
                };
                let result = select_full_qual(open_vcf(input).expect("open"), &options).expect("full");
                full_qual = Some(result.write_in(&out).expect("write qual"));
                (printed(&result.info_messages()), printed(result.warnings()))
            }
            "het" => {
                if let Some(phase_qual) = &phase_qual {
                    fs::copy(phase_qual, out.join(PHASE_QUAL_FILE)).expect("copy phase_qual");
                }
                let contig = (args[1] != "-").then(|| args[1].to_owned());
                let options = HetSnpOptions { contig, ..HetSnpOptions::default() };
                let result = select_het_snps_for_phasing(input, &out, &options).expect("het");
                let _ = fs::remove_file(out.join(PHASE_QUAL_FILE));
                (printed(&[result.info_message()]), String::new())
            }
            "cand" => {
                let qual = out.join(FULL_QUAL_FILE);
                match args[6] {
                    "latest" => {
                        fs::copy(full_qual.as_ref().expect("a full run before `latest`"), &qual).expect("copy qual");
                    }
                    "none" => {}
                    file => {
                        fs::copy(case.join(file), &qual).expect("copy qual file");
                    }
                }
                let options = CandidateOptions {
                    contig: (args[1] != "-").then(|| args[1].to_owned()),
                    var_pct_full: pct(2),
                    ref_pct_full: pct(3),
                    window: CandidateWindow::from_platform(args[4]),
                    split_bed_size: args[5].parse().expect("split_bed_size in runs.txt"),
                    ..CandidateOptions::default()
                };
                let result = select_candidates_for_full_alignment(input, &out, &options).expect("cand");
                assert_eq!(result.merge_output(), None, "the fixtures' pileup VCF is not named pileup.vcf");
                let _ = fs::remove_file(&qual);
                normalize_full_aln_lists(&out);
                (String::new(), printed(result.messages()))
            }
            other => panic!("unknown run kind {other:?} in {}", case.display()),
        };
        fs::write(out.join("stdout"), stdout).expect("write stdout");
        fs::write(out.join("stderr"), stderr).expect("write stderr");
    }
}

/// Every file under `expected` must exist in `got` with the same bytes, and nothing more.
fn assert_same_tree(expected: &Path, got: &Path) {
    let mut want: Vec<_> =
        fs::read_dir(expected).expect("expected dir").map(|e| e.expect("entry").file_name()).collect();
    let mut have: Vec<_> = fs::read_dir(got).expect("output dir").map(|e| e.expect("entry").file_name()).collect();
    want.sort();
    have.sort();
    assert_eq!(have, want, "files in {}", got.display());
    for name in want {
        let (e, g) = (expected.join(&name), got.join(&name));
        if e.is_dir() {
            assert_same_tree(&e, &g);
        } else {
            let (e_bytes, g_bytes) = (fs::read(&e).expect("read"), fs::read(&g).expect("read"));
            assert!(
                e_bytes == g_bytes,
                "{} differs from {}:\n--- expected\n{}\n--- got\n{}",
                g.display(),
                e.display(),
                String::from_utf8_lossy(&e_bytes),
                String::from_utf8_lossy(&g_bytes)
            );
        }
    }
}

#[test]
fn fixtures_match_the_original_scripts() {
    let root = Path::new(env!("CARGO_TARGET_TMPDIR")).join("snv-select-fixtures");
    let _ = fs::remove_dir_all(&root);
    let mut cases: Vec<PathBuf> =
        fs::read_dir(data_dir()).expect("data dir").map(|e| e.expect("entry").path()).collect();
    cases.sort();
    assert!(cases.len() >= 6, "expected the checked-in cases, found {cases:?}");
    for case in cases {
        let name = case.file_name().expect("case name").to_owned();
        let plain = root.join(&name).join("plain");
        run_case(&case, &case.join("input.vcf"), &plain);
        assert_same_tree(&case.join("expected"), &plain);

        let gz_input = root.join(&name).join("input.vcf.gz");
        gzip_copy(&case.join("input.vcf"), &gz_input);
        let gz = root.join(&name).join("gz");
        run_case(&case, &gz_input, &gz);
        assert_same_tree(&case.join("expected"), &gz);
    }
}

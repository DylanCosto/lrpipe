//! `SelectCandidates` through the public API: the no-candidate fallback, and byte parity with
//! a run of the original on a real pileup VCF.
//!
//! The parity test needs `SNV_CANDIDATES_ORACLE`, the original's output directory (`--output`)
//! from a run with `--platform ont` (`--var_pct_full 0.7 --ref_pct_full 0.1`): its
//! `pileup.vcf.gz`, `tmp/full_alignment_output/candidate_bed/` (`qual`, the BED chunks and
//! `FULL_ALN_FILE_*`) and `log/5_select_candidate.log`. `SNV_CANDIDATES_OUT` (optional) is
//! where to write the Rust outputs. It is skipped, with a message, when the oracle is unset.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use snv_select::{
    CandidateOptions, CandidateVcf, Error, FULL_ALN_FILE_PREFIX, FULL_QUAL_FILE, FullQualOptions, open_vcf,
    select_candidates_for_full_alignment, select_full_qual,
};

/// Every record of a small pileup VCF has QUAL 40 or more.
const HIGH_QUALITY_VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
    chr1\t100\t.\tA\tC\t40.0\tPASS\tP\tGT:GQ\t0/1:40\n\
    chr1\t200\t.\tA\t.\t45.0\tRefCall\tP\tGT:GQ\t0/0:45\n";

fn fresh_dir(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("bed")).expect("create test dir");
    fs::write(dir.join("bed").join(FULL_QUAL_FILE), "0 0\n").expect("write qual");
    dir
}

fn chr1() -> CandidateOptions {
    CandidateOptions { contig: Some("chr1".to_owned()), ..CandidateOptions::default() }
}

/// As the original's `print_calling_step`: checked against the script with a `pileup.vcf.gz` and a
/// plain `pileup.vcf` (`cp` then fails and nothing is copied, but the warning is printed).
#[test]
fn no_candidate_copies_pileup_as_merge_output() {
    let dir = fresh_dir("candidates-fallback-gz");
    let pileup = dir.join("pileup.vcf.gz");
    fs::write(&pileup, HIGH_QUALITY_VCF).expect("write pileup");
    let result = select_candidates_for_full_alignment(&pileup, &dir.join("bed"), &chr1()).expect("select");
    let merge_output = dir.join("merge_output.vcf.gz");
    assert!(result.candidates().is_empty());
    assert_eq!(result.full_aln_file(), None);
    assert_eq!(result.merge_output(), Some(merge_output.as_path()));
    assert_eq!(fs::read(&merge_output).expect("copied"), HIGH_QUALITY_VCF.as_bytes());
    assert_eq!(
        result.messages(),
        [
            "[INFO] Low quality reference calls to be processed in chr1: 0".to_owned(),
            "[INFO] Low quality variants to be processed in chr1: 0".to_owned(),
            "\x1b[93m[WARNING] Cannot find any low-quality 0/0, 0/1 or 1/1 variant in pileup output in contig chr1\x1b[0m"
                .to_owned(),
            format!("\x1b[93m[WARNING] Copying pileup.vcf.gz to {}\x1b[0m", merge_output.display()),
        ]
    );
    let bed: Vec<_> = fs::read_dir(dir.join("bed")).expect("bed dir").map(|e| e.expect("entry").file_name()).collect();
    assert_eq!(bed, [FULL_QUAL_FILE], "nothing but qual in the split folder");

    let dir = fresh_dir("candidates-fallback-plain");
    let pileup = dir.join("pileup.vcf");
    fs::write(&pileup, HIGH_QUALITY_VCF).expect("write pileup");
    let result = select_candidates_for_full_alignment(&pileup, &dir.join("bed"), &chr1()).expect("select");
    assert_eq!(result.messages().len(), 4);
    assert!(!dir.join("merge_output.vcf.gz").exists());

    let dir = fresh_dir("candidates-fallback-other-name");
    let pileup = dir.join("calls.vcf");
    fs::write(&pileup, HIGH_QUALITY_VCF).expect("write pileup");
    let result = select_candidates_for_full_alignment(&pileup, &dir.join("bed"), &chr1()).expect("select");
    assert_eq!((result.merge_output(), result.messages().len()), (None, 3));
}

#[test]
fn unported_paths_fail_clearly() {
    let dir = fresh_dir("candidates-unported");
    let pileup = dir.join("calls.vcf");
    fs::write(&pileup, HIGH_QUALITY_VCF).expect("write pileup");
    let options = CandidateOptions { phased_vcf_fn: Some(pileup.clone()), ..chr1() };
    let err = select_candidates_for_full_alignment(&pileup, &dir.join("bed"), &options).expect_err("unported");
    assert_eq!(err.to_string(), "SelectCandidates option --phased_vcf_fn is not ported");

    // Candidates but no --ctgName: the script raises in `'\t'.join`.
    fs::write(dir.join("bed").join(FULL_QUAL_FILE), "50 50").expect("write qual");
    let options = CandidateOptions { contig: None, ..CandidateOptions::default() };
    let err = select_candidates_for_full_alignment(&pileup, &dir.join("bed"), &options).expect_err("no contig");
    assert!(matches!(err, Error::MissingContigName), "{err:?}");

    fs::write(dir.join("bed").join(FULL_QUAL_FILE), "50  50").expect("write qual");
    let err = select_candidates_for_full_alignment(&pileup, &dir.join("bed"), &chr1()).expect_err("bad qual");
    assert!(matches!(err, Error::InvalidCutoffFile { .. }), "{err:?}");
}

fn assert_same_bytes(expected: &[u8], got: &[u8], what: &str) {
    assert!(expected == got, "{what} differs ({} vs {} bytes)", got.len(), expected.len());
}

#[test]
fn candidate_beds_match_the_original() {
    let Ok(oracle) = env::var("SNV_CANDIDATES_ORACLE") else {
        eprintln!("skipping: set SNV_CANDIDATES_ORACLE to the original's output directory (see tests/candidates.rs)");
        return;
    };
    let oracle = PathBuf::from(oracle);
    let oracle_bed = oracle.join("tmp/full_alignment_output/candidate_bed");
    let out = env::var("SNV_CANDIDATES_OUT")
        .map_or_else(|_| Path::new(env!("CARGO_TARGET_TMPDIR")).join("snv-candidates-parity"), PathBuf::from);
    let _ = fs::remove_dir_all(&out);
    let bed = out.join("candidate_bed");
    fs::create_dir_all(&bed).expect("create output dir");
    // A link named pileup.vcf.gz in the output, so a no-candidate fallback copies there and
    // never into the oracle.
    let pileup = out.join("pileup.vcf.gz");
    std::os::unix::fs::symlink(oracle.join("pileup.vcf.gz"), &pileup).expect("link pileup VCF");

    // SelectQual --var_pct_full 0.7 --ref_pct_full 0.1 --platform ont --vcf_fn EMPTY
    let options =
        FullQualOptions { var_pct_full: 0.7, ref_pct_full: 0.1, candidate_vcf: CandidateVcf::from_arg("EMPTY") };
    let full = select_full_qual(open_vcf(&pileup).expect("open pileup VCF"), &options).expect("SelectQual");
    let qual = full.write_in(&bed).expect("write qual");
    assert_same_bytes(
        &fs::read(oracle_bed.join(FULL_QUAL_FILE)).expect("oracle qual"),
        &fs::read(qual).expect("qual"),
        "qual",
    );
    let mut log: Vec<String> = full.info_messages().to_vec();

    let mut contigs: Vec<String> = fs::read_dir(&oracle_bed)
        .expect("oracle candidate_bed")
        .filter_map(|e| e.expect("entry").file_name().to_str()?.strip_prefix(FULL_ALN_FILE_PREFIX).map(str::to_owned))
        .collect();
    contigs.sort();
    assert!(!contigs.is_empty(), "no FULL_ALN_FILE_* in {}", oracle_bed.display());

    let (oracle_prefix, out_prefix) = (format!("{}/", oracle_bed.display()), format!("{}/", bed.display()));
    let mut regions = 0;
    for contig in &contigs {
        // SelectCandidates --var_pct_full 0.7 --ref_pct_full 0.1 --platform ont --ctgName <contig>
        let options = CandidateOptions {
            contig: Some(contig.clone()),
            var_pct_full: 0.7,
            ref_pct_full: 0.1,
            ..Default::default()
        };
        let result = select_candidates_for_full_alignment(&pileup, &bed, &options).expect("SelectCandidates");
        log.extend(result.messages().iter().cloned());
        let list_name = format!("{FULL_ALN_FILE_PREFIX}{contig}");
        let list = result.full_aln_file().expect("the oracle has candidates here");
        let want = fs::read_to_string(oracle_bed.join(&list_name)).expect("oracle list");
        // The lists name files by absolute path: compare with the folder swapped.
        assert_eq!(fs::read_to_string(list).expect("list"), want.replace(&oracle_prefix, &out_prefix), "{list_name}");
        for line in want.lines() {
            let name = line.strip_prefix(&oracle_prefix).expect("an oracle chunk path");
            let expected = fs::read(oracle_bed.join(name)).expect("oracle chunk");
            assert_same_bytes(&expected, &fs::read(bed.join(name)).expect("chunk"), name);
            regions += expected.iter().filter(|&&b| b == b'\n').count();
        }
        println!("{list_name}: {} chunks identical", want.lines().count());
    }
    println!("{regions} regions identical");

    // The same files, apart from the pipeline's concatenated FULL_ALN_FILES.
    let names = |dir: &Path| {
        let mut names: Vec<_> = fs::read_dir(dir)
            .expect("dir")
            .map(|e| e.expect("entry").file_name())
            .filter(|n| n != "FULL_ALN_FILES")
            .collect();
        names.sort();
        names
    };
    assert_eq!(names(&bed), names(&oracle_bed));

    // Stage 5's log: SelectQual's stdout, then every SelectCandidates job's stderr (in job
    // completion order when there are several).
    let want = fs::read_to_string(oracle.join("log/5_select_candidate.log")).expect("oracle log");
    let mut want: Vec<&str> = want.lines().collect();
    let mut got: Vec<&str> = log.iter().map(String::as_str).collect();
    if contigs.len() > 1 {
        want.sort_unstable();
        got.sort_unstable();
    }
    assert_eq!(got, want, "5_select_candidate.log");
}

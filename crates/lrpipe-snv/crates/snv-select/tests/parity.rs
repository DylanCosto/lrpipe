//! Byte parity with the original's own `SelectQual` and `SelectHetSnp` on a real pileup VCF.
//!
//! - `SNV_SELECT_ORACLE`: a directory with `pileup.vcf.gz`, `CONTIGS` (one name per line)
//!   and `phase_vcf/{phase_qual,<contig>.vcf}` from the pipeline's
//!   `SelectQual --phase --var_pct_phasing 0.7` and `SelectHetSnp --ctgName <contig>`.
//! - `SNV_SELECT_FULL_ORACLE` (optional): a directory with `qual` (and optionally `stdout`)
//!   from `SelectQual --var_pct_full 0.7 --ref_pct_full 0.1 --platform ont --vcf_fn EMPTY`,
//!   the ONT defaults the original's pipeline script (`run_*.py`) passes.
//! - `SNV_SELECT_OUT` (optional): where to write the Rust outputs.
//!
//! Skipped, with a message, when `SNV_SELECT_ORACLE` is unset.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use snv_select::{
    CandidateVcf, FullQualOptions, HetSnpOptions, PHASE_QUAL_FILE, open_vcf, select_full_qual,
    select_het_snps_for_phasing, select_phase_qual,
};

fn assert_same_file(expected: &Path, got: &Path) {
    let (e, g) = (fs::read(expected).expect("read oracle"), fs::read(got).expect("read output"));
    assert!(e == g, "{} differs from {} ({} vs {} bytes)", got.display(), expected.display(), g.len(), e.len());
    println!("identical: {} ({} bytes)", got.display(), g.len());
}

#[test]
fn pileup_vcf_selection_matches_the_original() {
    let Ok(oracle) = env::var("SNV_SELECT_ORACLE") else {
        eprintln!("skipping: set SNV_SELECT_ORACLE to an oracle directory (see tests/parity.rs)");
        return;
    };
    let oracle = PathBuf::from(oracle);
    let out = env::var("SNV_SELECT_OUT")
        .map_or_else(|_| Path::new(env!("CARGO_TARGET_TMPDIR")).join("snv-select-parity"), PathBuf::from);
    let pileup = oracle.join("pileup.vcf.gz");

    // Step 2/7: SelectQual --phase, then SelectHetSnp per contig.
    let phase_dir = out.join("phase_vcf");
    fs::create_dir_all(&phase_dir).expect("create output dir");
    let phase = select_phase_qual(open_vcf(&pileup).expect("open pileup VCF"), 0.7).expect("SelectQual --phase");
    println!("{}", phase.info_message());
    let phase_qual = phase.write_in(&phase_dir).expect("write phase_qual");
    assert_same_file(&oracle.join("phase_vcf").join(PHASE_QUAL_FILE), &phase_qual);

    let contigs = fs::read_to_string(oracle.join("CONTIGS")).expect("read CONTIGS");
    for contig in contigs.split_whitespace() {
        let options = HetSnpOptions { contig: Some(contig.to_owned()), ..HetSnpOptions::default() };
        let snps = select_het_snps_for_phasing(&pileup, &phase_dir, &options).expect("SelectHetSnp");
        println!("{}", snps.info_message());
        assert_same_file(&oracle.join("phase_vcf").join(snps.file_name()), &phase_dir.join(snps.file_name()));
    }

    // Step 5/7: SelectQual for full-alignment calling.
    let Ok(full_oracle) = env::var("SNV_SELECT_FULL_ORACLE") else {
        eprintln!("skipping full-alignment SelectQual: SNV_SELECT_FULL_ORACLE is unset");
        return;
    };
    let full_oracle = PathBuf::from(full_oracle);
    let full_dir = out.join("candidate_bed");
    fs::create_dir_all(&full_dir).expect("create output dir");
    let options =
        FullQualOptions { var_pct_full: 0.7, ref_pct_full: 0.1, candidate_vcf: CandidateVcf::from_arg("EMPTY") };
    let full = select_full_qual(open_vcf(&pileup).expect("open pileup VCF"), &options).expect("SelectQual");
    let qual = full.write_in(&full_dir).expect("write qual");
    assert_same_file(&full_oracle.join("qual"), &qual);
    if let Ok(stdout) = fs::read_to_string(full_oracle.join("stdout")) {
        let got: String = full.info_messages().iter().map(|l| format!("{l}\n")).collect();
        assert_eq!(got, stdout, "SelectQual stdout");
    }
}

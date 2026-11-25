//! Synthetic candidates against the original's `output_with`, with and without indel-length heads.
//!
//! `data/call_synthetic.tsv` comes from `tools/oracle_call_synthetic.py`: random reference bases,
//! alt infos and float32 rows that reach every outcome type, retries after zeroed-out outcomes,
//! exact ties and inputs the original rejects, with the line it writes for each.

use snv_call::{CallSource, IupacBases, OutputConfig, vcf_record, vcf_record_with_indel_lengths};

const CASES: &str = include_str!("data/call_synthetic.tsv");

#[test]
fn synthetic_cases_match_the_original() {
    let mut checked = [0usize; 2];
    let original = OutputConfig { iupac: IupacBases::ToN, ..OutputConfig::default() };
    for (n, case) in CASES.lines().enumerate() {
        let [mode, position, alt_info, probabilities, expected] = case.splitn(5, '\t').collect::<Vec<_>>()[..] else {
            panic!("line {}: malformed case", n + 1);
        };
        let probabilities: Vec<f32> = probabilities.split(',').map(|p| p.parse().expect("float32 repr")).collect();
        let got = match mode {
            "fa" => {
                checked[0] += 1;
                let config = OutputConfig { source: CallSource::FullAlignment, ..original.clone() };
                vcf_record_with_indel_lengths(position, alt_info, &probabilities, &config)
            }
            "p" => {
                checked[1] += 1;
                vcf_record(position, alt_info, &probabilities, &original)
            }
            other => panic!("line {}: unknown mode {other:?}", n + 1),
        };
        let got = match got {
            Ok(Some(line)) => line.strip_suffix('\n').expect("records end with a newline").to_owned(),
            Ok(None) => "None".to_owned(),
            Err(_) => "ERROR".to_owned(),
        };
        assert_eq!(got, expected, "line {}: {mode} {position} {alt_info}", n + 1);
    }
    assert!(checked.iter().all(|&c| c > 300), "too few cases: {checked:?}");
}

#[test]
fn short_rows_are_rejected() {
    let config = OutputConfig { source: CallSource::FullAlignment, ..OutputConfig::default() };
    assert!(vcf_record_with_indel_lengths("chr1:10:A", "10-RA 10", &[0.0; 89], &config).is_err());
    assert!(vcf_record("chr1:10:A", "10-RA 10", &[0.0; 23], &config).is_err());
}

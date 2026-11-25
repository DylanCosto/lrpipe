//! VCF float text against htslib 1.18 (`tools/htslib_float_vectors.py`).

use lrpipe_cnv::output::htsfloat::vcf_float;

#[test]
fn floats_print_like_htslib() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/htslib_floats.json");
    let pairs: Vec<(String, String)> = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let failures: Vec<String> = pairs
        .iter()
        .filter(|(text, want)| vcf_float(text) != *want)
        .map(|(text, want)| format!("{text}: got {}, want {want}", vcf_float(text)))
        .collect();
    assert!(
        failures.is_empty(),
        "{} of {} differ:\n{}",
        failures.len(),
        pairs.len(),
        failures[..failures.len().min(20)].join("\n")
    );
}

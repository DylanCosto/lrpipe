//! Full-alignment chunk VCFs rebuilt from the original's own inputs, compared byte for byte.
//!
//! `SNV_FA_ORACLE` names a directory written by `tools/oracle_fa_output.py`: for each job,
//! `PREFIX.info.tsv` (position and alt-info strings), `PREFIX.probs.npy` (float32 rows of 90
//! probabilities) and `PREFIX.job.tsv` (the VCF the original wrote, its reference, `--cmd_fn` and
//! sample). The test is skipped when the variable is unset.

use std::fs;
use std::path::{Path, PathBuf};

use snv_call::{CallSource, OutputConfig, vcf_header, vcf_record_with_indel_lengths};

/// Columns of a full-alignment probability row.
const ROW: usize = 90;

/// Reads a C-order little-endian float32 2-D `.npy` file (format 1.0, as `np.save` writes it).
fn read_npy_f32(path: &Path) -> (Vec<f32>, usize) {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    assert_eq!(&bytes[..6], b"\x93NUMPY", "{} is not a .npy file", path.display());
    assert_eq!(bytes[6], 1, "unsupported .npy version in {}", path.display());
    let header_len = usize::from(u16::from_le_bytes([bytes[8], bytes[9]]));
    let header = std::str::from_utf8(&bytes[10..10 + header_len]).expect("ASCII .npy header");
    assert!(header.contains("'descr': '<f4'") && header.contains("'fortran_order': False"), "{header}");
    let shape = &header[header.find("'shape': (").expect("shape in header") + 10..];
    let rows: usize = shape[..shape.find(',').expect("2-D shape")].trim().parse().expect("row count");
    let data: Vec<f32> = bytes[10 + header_len..].as_chunks::<4>().0.iter().map(|&b| f32::from_le_bytes(b)).collect();
    assert_eq!(data.len(), rows * ROW, "{} is not {rows} x {ROW}", path.display());
    (data, rows)
}

/// Rebuilds one job's VCF and compares it with the one the original wrote; returns the record
/// count.
fn check_job(prefix: &str) -> usize {
    let job = fs::read_to_string(format!("{prefix}.job.tsv")).expect("job.tsv");
    let [vcf, reference, cmd_fn, sample] = job.trim_end().split('\t').collect::<Vec<_>>()[..] else {
        panic!("malformed {prefix}.job.tsv");
    };
    let info = fs::read_to_string(format!("{prefix}.info.tsv")).expect("info.tsv");
    let (probs, rows) = read_npy_f32(Path::new(&format!("{prefix}.probs.npy")));
    let candidates: Vec<(&str, &str)> =
        info.lines().map(|l| l.split_once('\t').expect("position and alt info")).collect();
    assert_eq!(candidates.len(), rows, "{prefix}: candidates vs probability rows");

    let config = OutputConfig { source: CallSource::FullAlignment, ..OutputConfig::default() };
    let mut out = vcf_header(Some(Path::new(reference)), Some(Path::new(cmd_fn)), sample).expect("header");
    let mut records = 0;
    for ((position, alt_info), row) in candidates.iter().zip(probs.as_chunks::<ROW>().0) {
        if let Some(line) = vcf_record_with_indel_lengths(position, alt_info, row, &config).expect("record") {
            out.push_str(&line);
            records += 1;
        }
    }
    let expected = fs::read_to_string(vcf).expect("the original's VCF");
    if out != expected {
        let first = out.lines().zip(expected.lines()).position(|(a, b)| a != b);
        let (got, want) = match first {
            Some(i) => (out.lines().nth(i), expected.lines().nth(i)),
            None => (Some("<length differs>"), None),
        };
        panic!("{prefix}: differs from {vcf}\n  got:  {got:?}\n  want: {want:?}");
    }
    records
}

#[test]
fn full_alignment_vcfs_match_the_original() {
    let Some(dir) = std::env::var_os("SNV_FA_ORACLE").map(PathBuf::from) else {
        eprintln!("skipped: set SNV_FA_ORACLE to a tools/oracle_fa_output.py output directory");
        return;
    };
    let mut prefixes: Vec<String> = fs::read_dir(&dir)
        .expect("oracle directory")
        .filter_map(|e| e.ok()?.path().to_str()?.strip_suffix(".job.tsv").map(str::to_owned))
        .collect();
    prefixes.sort();
    assert!(!prefixes.is_empty(), "no *.job.tsv in {}", dir.display());
    let mut total = 0;
    for prefix in &prefixes {
        let n = check_job(prefix);
        eprintln!("{prefix}: {n} records identical");
        total += n;
    }
    eprintln!("{} jobs, {total} records identical", prefixes.len());
}

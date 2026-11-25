//! Checks lrpipe-depth's outputs against mosdepth 0.3.14's on every case in `tools/cases.tsv`.
//!
//! The expected hashes in `tests/data/expected.sha256` come from real mosdepth runs
//! (`tools/parity.sh --record`), so this test needs no mosdepth.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::Command;

use sha2::{Digest, Sha256};

const OUTPUTS: [&str; 5] = [
    "mosdepth.summary.txt",
    "mosdepth.global.dist.txt",
    "mosdepth.region.dist.txt",
    "regions.bed.gz",
    "thresholds.bed.gz",
];

fn text_sha256(path: &Path) -> String {
    let raw = fs::read(path).unwrap();
    let mut text = Vec::new();
    if path.extension().is_some_and(|e| e == "gz") {
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut text).unwrap();
    } else {
        text = raw;
    }
    Sha256::digest(&text).iter().map(|b| format!("{b:02x}")).collect()
}

#[test]
fn outputs_match_mosdepth() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut expected: BTreeMap<(String, String), String> = BTreeMap::new();
    for line in fs::read_to_string(root.join("tests/data/expected.sha256")).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        expected.insert((f[0].to_owned(), f[1].to_owned()), f[2].to_owned());
    }
    let work = std::env::temp_dir().join(format!("lrpipe-depth-parity-{}", std::process::id()));
    fs::create_dir_all(&work).unwrap();

    let mut failures = Vec::new();
    let mut cases = 0;
    for line in fs::read_to_string(root.join("tools/cases.tsv")).unwrap().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (name, input, prec, opts) = (f[0], f[1], f[2], f[3]);
        cases += 1;
        let prefix = work.join(name);
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_lrpipe-depth"));
        cmd.current_dir(root).args(["-t", "2"]).args(opts.split_whitespace());
        cmd.arg(&prefix).arg(root.join("tests/data").join(input));
        cmd.env_remove("MOSDEPTH_PRECISION");
        if prec != "-" {
            cmd.env("MOSDEPTH_PRECISION", prec);
        }
        let out = cmd.output().unwrap();
        if !out.status.success() {
            failures.push(format!("{name}: lrpipe-depth failed: {}", String::from_utf8_lossy(&out.stderr)));
            continue;
        }
        for o in OUTPUTS {
            let path = Path::new(&format!("{}.{o}", prefix.display())).to_owned();
            match (expected.get(&(name.to_owned(), o.to_owned())), path.exists()) {
                (Some(want), true) => {
                    let got = text_sha256(&path);
                    if &got != want {
                        failures.push(format!("{name}: {o} differs from mosdepth's"));
                    }
                }
                (Some(_), false) => failures.push(format!("{name}: {o} not written")),
                (None, true) => failures.push(format!("{name}: {o} written, but mosdepth writes none")),
                (None, false) => {}
            }
        }
        if !expected.keys().any(|(n, _)| n == name) {
            failures.push(format!("{name}: no recorded mosdepth output (run tools/parity.sh --record)"));
        }
    }
    let _ = fs::remove_dir_all(&work);
    assert!(failures.is_empty(), "{} of {cases} cases:\n{}", failures.len(), failures.join("\n"));
}

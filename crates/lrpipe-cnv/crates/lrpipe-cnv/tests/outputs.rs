//! Runs the whole pipeline on each reference run and compares the output
//! files with Spectre's: VCF and BED byte for byte (decompressed), the `.spc` with each
//! `merged_sample_references` list sorted (Spectre writes those Python sets in
//! hash-randomized order). Also checks `.mdr` creation against Spectre's.
//!
//! Needs `SPECTRE_BENCH` and `SPECTRE_REFERENCE` (see tests/reference.rs); skipped without, unless
//! `LRPIPE_REQUIRE_SPECTRE` is set (as `tools/check.sh` sets it where the data is), when it fails.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use lrpipe_cnv::metadata::write_mdr;
use lrpipe_cnv::pipeline::{Options, OutputPaths, call_cnvs, write_outputs};

fn env() -> Option<(PathBuf, PathBuf)> {
    let found = std::env::var_os("SPECTRE_BENCH").zip(std::env::var_os("SPECTRE_REFERENCE"));
    assert!(
        found.is_some() || std::env::var_os("LRPIPE_REQUIRE_SPECTRE").is_none(),
        "LRPIPE_REQUIRE_SPECTRE is set, but SPECTRE_BENCH or SPECTRE_REFERENCE is not"
    );
    found.map(|(bench, reference)| (PathBuf::from(bench), PathBuf::from(reference)))
}

fn gunzip(path: &Path) -> String {
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap()).read_to_string(&mut s).unwrap();
    s
}

/// The `.spc` text with the lines of each `merged_sample_references` list sorted.
fn canonical_spc(text: &str) -> String {
    let lines: Vec<&str> = text.split('\n').collect();
    let mut out: Vec<String> = Vec::with_capacity(lines.len());
    let mut i = 0;
    while i < lines.len() {
        out.push(lines[i].to_owned());
        if lines[i].ends_with("\"merged_sample_references\": [") {
            let mut j = i + 1;
            while !lines[j].trim_start().starts_with(']') {
                j += 1;
            }
            let mut items: Vec<&str> = lines[i + 1..j].iter().map(|l| l.trim_end_matches(',')).collect();
            items.sort_unstable();
            let n = items.len();
            out.extend(
                items.into_iter().enumerate().map(|(k, it)| if k + 1 < n { format!("{it},") } else { it.to_owned() }),
            );
            i = j;
            continue;
        }
        i += 1;
    }
    out.join("\n")
}

/// The first differing line of two texts, for the failure message.
fn first_difference(a: &str, b: &str) -> String {
    match a.lines().zip(b.lines()).enumerate().find(|(_, (x, y))| x != y) {
        Some((n, (x, y))) => format!("line {}: got {x:?}, want {y:?}", n + 1),
        None => format!("lengths differ: {} versus {} lines", a.lines().count(), b.lines().count()),
    }
}

#[test]
fn output_files_match_spectre() {
    let Some((bench, reference)) = env() else {
        eprintln!("SPECTRE_BENCH / SPECTRE_REFERENCE not set: skipping the output comparison");
        return;
    };
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/reference_runs.tsv");
    let mut failures = Vec::new();
    for line in fs::read_to_string(manifest).unwrap().lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let (name, sample, coverage, snv, only_chr) = (f[0], f[1], f[2], f[3], f[4]);
        let out = std::env::temp_dir().join(format!("lrpipe-cnv-outputs-{}-{name}", std::process::id()));
        let mut opts = Options::new(bench.join(coverage), sample.to_owned(), out.clone(), reference.clone());
        opts.metadata = Some(bench.join("inputs/grch38_reference.mdr"));
        opts.blacklist = Some(bench.join("inputs/grch38_blacklist_spectre.bed"));
        opts.only_chr = (only_chr != "-").then(|| only_chr.split(',').map(str::to_owned).collect());
        opts.snv = (snv != "-").then(|| bench.join(snv));
        opts.dip_seed = if name.ends_with("seed7") { 7 } else { 0 };
        let run = call_cnvs(&opts).unwrap();
        let got = write_outputs(&run, &opts).unwrap();
        let want = OutputPaths::new(&bench.join("ref").join(name).join("out"), sample);
        for (what, g, w) in [("VCF", &got.vcf, &want.vcf), ("BED", &got.bed, &want.bed)] {
            let (g, w) = (gunzip(g), gunzip(w));
            if g != w {
                failures.push(format!("{name} {what}: {}", first_difference(&g, &w)));
            }
        }
        let (g, w) = (canonical_spc(&gunzip(&got.spc)), canonical_spc(&gunzip(&want.spc)));
        if g != w {
            failures.push(format!("{name} .spc: {}", first_difference(&g, &w)));
        }
        let n: usize = run.calls.iter().map(|(_, c)| c.len()).sum();
        eprintln!("{name}: {n} calls; VCF, BED and .spc compared");
        let _ = fs::remove_dir_all(&out);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn mdr_matches_spectre() {
    let Some((bench, reference)) = env() else {
        eprintln!("SPECTRE_BENCH / SPECTRE_REFERENCE not set: skipping the .mdr comparison");
        return;
    };
    let out = std::env::temp_dir().join(format!("lrpipe-cnv-{}.mdr", std::process::id()));
    write_mdr(&reference, 5, &out).unwrap();
    let (got, want) = (fs::read(&out).unwrap(), fs::read(bench.join("inputs/grch38_reference.mdr")).unwrap());
    let _ = fs::remove_file(&out);
    assert!(got == want, "{}", first_difference(&String::from_utf8_lossy(&got), &String::from_utf8_lossy(&want)));
}

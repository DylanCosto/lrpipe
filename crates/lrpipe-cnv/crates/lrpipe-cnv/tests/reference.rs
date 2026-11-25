//! Compares each ported stage with Spectre 0.2.1's state in the reference runs listed in
//! `tools/reference_runs.tsv` (dumps from `tools/run_spectre.py`).
//!
//! The reference data is too large to commit; set `SPECTRE_BENCH` to the directory holding
//! `inputs/` and `ref/`, and `SPECTRE_REFERENCE` to the reference FASTA (with `.fai`).
//! Without them the test prints a note and passes, unless `LRPIPE_REQUIRE_SPECTRE` is set (as
//! `tools/check.sh` sets it where the data is), when it fails.

use std::fs;
use std::path::{Path, PathBuf};

use lrpipe_cnv::call::call;
use lrpipe_cnv::candidate::{Arena, Candidate, CandidateId};
use lrpipe_cnv::coverage::read_regions;
use lrpipe_cnv::genome::Genome;
use lrpipe_cnv::metadata::Metadata;
use lrpipe_cnv::normalize::{NormalizeOptions, Normalized, Stats, normalize};
use lrpipe_cnv::number::Number;
use lrpipe_cnv::refine::{ContigContext, RefineOptions, refine};
use lrpipe_cnv::snv::{SnvVcf, cn_neutral};
use lrpipe_cnv::thresholds::{Clamping, Thresholds, WindowSample, seed_rng};
use pynum::numpy_random::RandomState;
use serde_json::Value;
use sha2::{Digest, Sha256};

struct Run {
    name: String,
    sample: String,
    coverage: PathBuf,
    snv: Option<PathBuf>,
    only_chr: Option<String>,
}

fn runs(bench: &Path) -> Vec<Run> {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/reference_runs.tsv");
    fs::read_to_string(manifest)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Run {
                name: f[0].into(),
                sample: f[1].into(),
                coverage: bench.join(f[2]),
                snv: (f[3] != "-").then(|| bench.join(f[3])),
                only_chr: (f[4] != "-").then(|| f[4].to_owned()),
            }
        })
        .collect()
}

/// The dump written after stage `stage` (file names are numbered by order of stages run).
fn dump(bench: &Path, run: &str, stage: &str) -> Value {
    let dir = bench.join("ref").join(run).join("dump");
    let path = fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .find(|p| p.file_stem().unwrap().to_str().unwrap().ends_with(&format!("_{stage}")))
        .unwrap_or_else(|| panic!("no {stage} dump in {}", dir.display()));
    serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap()
}

/// `tools/run_spectre.py`'s array summary: SHA-256 of float64 little-endian bytes, NaNs normalized.
fn sha(values: impl IntoIterator<Item = f64>) -> String {
    let mut h = Sha256::new();
    for v in values {
        let v = if v.is_nan() { f64::from_bits(0x7ff8_0000_0000_0000) } else { v };
        h.update(v.to_le_bytes());
    }
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}

struct Check {
    run: String,
    failures: Vec<String>,
}

impl Check {
    fn float(&mut self, what: &str, got: f64, want: &Value) {
        let ok = match (want.as_f64(), want.as_str()) {
            (Some(w), _) => got.to_bits() == w.to_bits(),
            (None, Some("inf")) => got == f64::INFINITY,
            (None, Some("-inf")) => got == f64::NEG_INFINITY,
            _ => want.is_null() && got.is_nan(),
        };
        if !ok {
            self.failures.push(format!("{}: {what}: got {got:?}, want {want}", self.run));
        }
    }

    fn number(&mut self, what: &str, got: Number, want: &Value) {
        let same_kind = matches!(got, Number::Int(_)) == want.is_i64();
        if !same_kind {
            self.failures.push(format!("{}: {what}: got {got:?}, want {want} (int versus float)", self.run));
        }
        self.float(what, got.value(), want);
    }

    fn array(&mut self, what: &str, got: impl IntoIterator<Item = f64>, want: &Value) {
        // Arrays are {"array": summary}; pandas columns hold the summary under "series".
        let want = if want.get("array").is_some() { &want["array"] } else { want };
        let got: Vec<f64> = got.into_iter().collect();
        let (len, nan) = (got.len() as u64, got.iter().filter(|v| v.is_nan()).count() as u64);
        if want["len"] != len || want["nan"] != nan || want["sha256"] != sha(got.iter().copied()) {
            self.failures.push(format!(
                "{}: {what}: got len {len}, {nan} NaN, head {:?}; want len {}, {} NaN, head {}",
                self.run,
                &got[..got.len().min(3)],
                want["len"],
                want["nan"],
                want["head"]
            ));
        }
    }

    fn stats(&mut self, what: &str, got: &Stats, want: &Value) {
        self.float(&format!("{what}.average"), got.average, &want["average"]);
        self.float(&format!("{what}.std_dev"), got.std_dev, &want["std_dev"]);
        self.float(&format!("{what}.median"), got.median, &want["median"]);
        self.float(&format!("{what}.min"), got.min, &want["min"]);
        self.float(&format!("{what}.max"), got.max, &want["max"]);
    }
}

/// Every object in a dump by its `@id`, so `{"@ref": ...}` entries can be followed.
fn objects(v: &Value, out: &mut std::collections::HashMap<String, Value>) {
    match v {
        Value::Object(m) => {
            if let Some(Value::String(id)) = m.get("@id") {
                out.insert(id.clone(), v.clone());
            }
            m.values().for_each(|x| objects(x, out));
        }
        Value::Array(a) => a.iter().for_each(|x| objects(x, out)),
        _ => {}
    }
}

/// The candidates of `list` (a dump's `{contig: [candidate | ref]}`), references resolved.
fn candidates(dump: &Value, list: &str) -> Vec<(String, Vec<Value>)> {
    let mut by_id = std::collections::HashMap::new();
    objects(dump, &mut by_id);
    let resolve = |v: &Value| match v.get("@ref").and_then(Value::as_str) {
        Some(r) => by_id[r].clone(),
        None => v.clone(),
    };
    dump[list]
        .as_object()
        .map(|m| m.iter().map(|(k, v)| (k.clone(), v.as_array().unwrap().iter().map(resolve).collect())).collect())
        .unwrap_or_default()
}

impl Check {
    fn int(&mut self, what: &str, got: i64, want: &Value) {
        if want.as_i64() != Some(got) {
            self.failures.push(format!("{}: {what}: got {got}, want {want}", self.run));
        }
    }

    fn text(&mut self, what: &str, got: &str, want: &Value) {
        if want.as_str() != Some(got) {
            self.failures.push(format!("{}: {what}: got {got:?}, want {want}", self.run));
        }
    }

    fn candidate(&mut self, what: &str, c: &Candidate, w: &Value) {
        let what = format!("{what} {}", c.id);
        self.text(&format!("{what} id"), &c.id, &w["id"]);
        self.text(&format!("{what} type"), c.cnv_type.as_str(), &w["type"]);
        self.int(&format!("{what} start"), c.start, &w["start"]);
        self.int(&format!("{what} end"), c.end, &w["end"]);
        self.int(&format!("{what} size"), c.size, &w["size"]);
        self.int(&format!("{what} size_kb"), c.size_kb, &w["size_kb"]);
        self.int(&format!("{what} cn_status"), c.cn_status, &w["cn_status"]);
        self.text(&format!("{what} gt"), c.gt, &w["gt"]);
        self.float(&format!("{what} median_cov_norm"), c.median_cov_norm, &w["median_cov_norm"]);
        self.float(&format!("{what} median_raw_cov"), c.median_raw_cov, &w["median_raw_cov"]);
        self.array(&format!("{what} pos"), c.pos.iter().map(|&p| p as f64), &w["pos"]);
        self.array(&format!("{what} cov"), c.cov.iter().copied(), &w["cov"]);
        let mut got: Vec<&str> = c.merged_sample_references.iter().map(String::as_str).collect();
        got.sort_unstable();
        let mut want: Vec<&str> = match &w["merged_sample_references"] {
            Value::Array(a) => a.iter().filter_map(Value::as_str).collect(),
            _ => Vec::new(),
        };
        want.sort_unstable();
        if got != want {
            self.failures.push(format!("{}: {what} merged_sample_references: got {got:?}, want {want:?}", self.run));
        }
    }

    /// Our lists (contig, candidates) against a dump's.
    fn lists(&mut self, what: &str, arena: &Arena, got: &[(String, Vec<CandidateId>)], want: &[(String, Vec<Value>)]) {
        for (contig, wants) in want {
            let Some((_, gots)) = got.iter().find(|(c, _)| c == contig) else {
                self.failures.push(format!("{}: {what}: no candidates for {contig}", self.run));
                continue;
            };
            if gots.len() != wants.len() {
                self.failures.push(format!(
                    "{}: {what} {contig}: {} candidates, want {}",
                    self.run,
                    gots.len(),
                    wants.len()
                ));
            }
            for (&id, w) in gots.iter().zip(wants) {
                self.candidate(&format!("{what} {contig}"), arena.get(id), w);
            }
        }
    }
}

fn normalize_run(bench: &Path, reference: &Path, run: &Run) -> Normalized {
    let genome = Genome::from_fasta_index(reference, 1_000_000, 2, &[]).unwrap();
    let only_chr: Vec<String> = match &run.only_chr {
        Some(list) => list.split(',').map(str::to_owned).collect(),
        None => genome.contigs.iter().map(|c| c.name.clone()).collect(),
    };
    let snv_cn_neutral = run.snv.as_ref().map(|path| {
        let vcf = SnvVcf::read(path).unwrap();
        let seed = if run.name.ends_with("seed7") { 7 } else { 0 };
        let cn = cn_neutral(&vcf, &genome, &mut RandomState::new(seed)).unwrap();
        let want = &dump(bench, &run.name, "snv_cn_neutral");
        for (what, got) in [
            ("cn_neutral_coverage_med", cn.median),
            ("cn_neutral_coverage_mean", cn.mean),
            ("cn_neutral_coverage_stdev", cn.stdev),
        ] {
            assert_eq!(Some(got.to_bits()), want[what].as_f64().map(f64::to_bits), "{}: {what}", run.name);
        }
        cn.median
    });
    let opts = NormalizeOptions { snv_cn_neutral, ..Default::default() };
    let depths = read_regions(&run.coverage).unwrap();
    let mdr = bench.join("inputs/grch38_reference.mdr");
    let blacklist = bench.join("inputs/grch38_blacklist_spectre.bed");
    let (norm, _) =
        normalize(&genome, depths, &only_chr, &opts, |bin| Metadata::load(&mdr, Some(&blacklist), bin)).unwrap();
    norm
}

#[test]
fn stages_match_spectre() {
    let (Some(bench), Some(reference)) = (std::env::var_os("SPECTRE_BENCH"), std::env::var_os("SPECTRE_REFERENCE"))
    else {
        assert!(
            std::env::var_os("LRPIPE_REQUIRE_SPECTRE").is_none(),
            "LRPIPE_REQUIRE_SPECTRE is set, but SPECTRE_BENCH or SPECTRE_REFERENCE is not"
        );
        eprintln!("SPECTRE_BENCH / SPECTRE_REFERENCE not set: skipping the reference comparison");
        return;
    };
    let (bench, reference) = (PathBuf::from(bench), PathBuf::from(reference));
    let mut all_failures = Vec::new();
    for run in runs(&bench) {
        let mut check = Check { run: run.name.clone(), failures: Vec::new() };
        let norm = normalize_run(&bench, &reference, &run);
        let want = dump(&bench, &run.name, "normalization");
        check.float("bin_size", norm.bin_size as f64, &want["bin_size"]);
        check.float("genome_wide_median", norm.genome_wide_median, &want["genome_wide_median"]);
        check.float("genome_wide_mean", norm.genome_wide_mean, &want["genome_wide_mean"]);
        check.float("genome_std", norm.genome_std, &want["genome_std"]);
        check.float("normalization_value", norm.normalization_value, &want["normalization_value"]);
        let analysis = want["coverage_analysis"].as_object().unwrap();
        // The dump's keys are sorted, so only the set of contigs can be compared here.
        let mut got_names: Vec<&str> = norm.contigs.iter().map(|c| c.name.as_str()).collect();
        got_names.sort_unstable();
        let mut want_names: Vec<&str> = analysis.keys().map(String::as_str).collect();
        want_names.sort_unstable();
        if got_names != want_names {
            check.failures.push(format!("{}: analysed contigs {got_names:?}, want {want_names:?}", run.name));
        }
        for c in &norm.contigs {
            let Some(w) = analysis.get(&c.name) else { continue };
            let data = &w["cov_data"];
            check.array(&format!("{} positions", c.name), c.positions.iter().map(|&p| p as f64), &data["positions"]);
            check.array(&format!("{} coverage_raw", c.name), c.raw.iter().copied(), &data["coverage_raw"]);
            check.array(&format!("{} normalized_cov", c.name), c.normalized.iter().copied(), &data["normalized_cov"]);
            check.array(
                &format!("{} normalized_cov_ploidy", c.name),
                c.normalized_ploidy.iter().copied(),
                &data["normalized_cov_ploidy"],
            );
            check.stats(&format!("{} statistics", c.name), &c.stats, &w["statistics"]);
            check.stats(&format!("{} norm_statistics", c.name), &c.norm_stats, &w["norm_statistics"]);
        }

        let mut rng = seed_rng(&run.sample);
        let sample = WindowSample::draw(&norm, &run.sample, &mut rng);
        let metrics = &want["cnv_metrics"];
        let frame = &metrics["df_coverage_candidate_no_blacklist_random_samples"]["frame"];
        check.array("sample coverage", sample.windows.iter().map(|w| w.coverage), &frame["coverage"]["series"]);
        check.array("sample positions", sample.windows.iter().map(|w| w.position as f64), &frame["position"]["series"]);
        let t = Thresholds::from_quantiles(&sample, 2, Clamping::PloidyBounds);
        check.number("del_threshold", t.del, &metrics["del_threshold"]);
        check.number("dup_threshold", t.dup, &metrics["dup_threshold"]);
        check.number("lower_2n_threshold", t.del, &want["lower_2n_threshold"]);
        check.number("upper_2n_threshold", t.dup, &want["upper_2n_threshold"]);

        // Raw candidates, then refinement.
        let mut arena = Arena::default();
        let raw = call(&norm, &t, &mut rng, &mut arena);
        check.lists("raw", &arena, &raw.contigs, &candidates(&dump(&bench, &run.name, "raw_calls"), "cnv_calls_list"));
        let metadata = Metadata::load(
            &bench.join("inputs/grch38_reference.mdr"),
            Some(&bench.join("inputs/grch38_blacklist_spectre.bed")),
            norm.bin_size,
        )
        .unwrap();
        let strict = Thresholds::from_sd(&sample, 1.85, 2, Clamping::PloidyBounds);
        let mut refined = Vec::new();
        for (contig, ids) in &raw.contigs {
            let coverage = norm.contigs.iter().find(|c| &c.name == contig).unwrap();
            let ctx = ContigContext {
                coverage,
                metadata: &metadata,
                bin_size: norm.bin_size,
                thresholds: &t,
                strict: &strict,
            };
            refined.push((contig.clone(), refine(ids, &mut arena, &ctx, &RefineOptions::default()).unwrap()));
        }
        let after = dump(&bench, &run.name, "refined");
        check.lists("refined", &arena, &refined, &candidates(&after, "cnv_calls_list"));
        check.lists("raw after refinement", &arena, &raw.contigs, &candidates(&after, "raw_cnv_calls_list"));
        let count = |l: &[(String, Vec<CandidateId>)]| l.iter().map(|(_, v)| v.len()).sum::<usize>();
        eprintln!(
            "{}: {} raw candidates, {} calls; {} mismatches",
            run.name,
            count(&raw.contigs),
            count(&refined),
            check.failures.len()
        );
        for f in &check.failures {
            eprintln!("  {f}");
        }
        all_failures.extend(check.failures);
    }
    assert!(all_failures.is_empty(), "{} mismatches:\n{}", all_failures.len(), all_failures.join("\n"));
}

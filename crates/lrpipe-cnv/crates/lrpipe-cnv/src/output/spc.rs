//! `{sample}.spc.gz`: Spectre's intermediate file, for population calling (analysis.py:784-840).

use std::io::Write;
use std::path::Path;

use pynum::json::Json;

use crate::Error;
use crate::candidate::{Arena, Candidate, CandidateId};
use crate::genome::Genome;
use crate::loh::{LohCalls, LohCandidate};
use crate::number::Number;
use crate::thresholds::{Thresholds, WindowSample};

/// Settings recorded in the file's `analysis_metrics`.
#[derive(Clone, Debug)]
pub struct Metrics<'a> {
    pub bin_size: i64,
    pub min_cnv_len: i64,
    pub normalization_value: f64,
    pub max_std_outlier_rm: i64,
    pub mosdepth_cov_genome_chr_diff: f64,
    pub cov_diff_threshold: f64,
    pub dist_proportion: f64,
    pub dist_min_overwrite: i64,
    pub thresholds: &'a Thresholds,
    pub strict: &'a Thresholds,
    pub sample: &'a WindowSample,
}

fn number(n: Number) -> Json {
    match n {
        Number::Int(i) => Json::Int(i),
        Number::Float(f) => Json::Float(f),
    }
}

/// `dict(candidate.__dict__)`, keys in `CNVCandidate.__init__` order.
fn candidate(c: &Candidate, sample_id: &str) -> Json {
    let statistics = match &c.statistics {
        Some(z) => Json::Dict(vec![(
            "z-score".into(),
            Json::Dict(vec![
                ("statistics".into(), Json::Null),
                ("score".into(), number(z.score)),
                ("pvalue".into(), number(z.pvalue)),
                ("sample_score".into(), Json::Int(z.sample_score)),
            ]),
        )]),
        None => Json::Dict(vec![]),
    };
    Json::Dict(vec![
        ("chromosome".into(), Json::Str(c.chromosome.clone())),
        ("start".into(), Json::Int(c.start)),
        ("end".into(), Json::Int(c.end)),
        ("size".into(), Json::Int(c.size)),
        ("bin_size".into(), Json::Int(c.bin_size)),
        ("pos".into(), Json::List(c.pos.iter().map(|&p| Json::Int(p)).collect())),
        ("cov".into(), Json::List(c.cov.iter().map(|&v| Json::Float(v)).collect())),
        ("normalization_value".into(), Json::Float(1.0)),
        ("type".into(), Json::Str(c.cnv_type.as_str().into())),
        ("cn_status".into(), Json::Int(c.cn_status)),
        ("gt".into(), Json::Str(c.gt.into())),
        ("median_cov_norm".into(), Json::Float(c.median_cov_norm)),
        ("median_raw_cov".into(), Json::Float(c.median_raw_cov)),
        ("size_kb".into(), Json::Int(c.size_kb)),
        ("is_init".into(), Json::Bool(true)),
        ("filter".into(), Json::Str(".".into())),
        ("quality".into(), Json::Str(".".into())),
        ("id".into(), Json::Str(c.id.clone())),
        ("het_score".into(), Json::Float(0.0)),
        ("statistics".into(), statistics),
        ("support_cnv_calls".into(), Json::Dict(vec![])),
        ("sample_origin".into(), Json::Str(sample_id.into())),
        (
            "merged_sample_references".into(),
            Json::List(c.merged_sample_references.iter().map(|r| Json::Str(r.clone())).collect()),
        ),
        ("sv_support".into(), Json::Bool(false)),
        ("merged_sample".into(), Json::Bool(false)),
    ])
}

/// A `LoHCandidate.__dict__`: the `CNVCandidate` keys, then LoH's own.
fn loh_candidate(c: &LohCandidate) -> Json {
    let nan = || Json::Str("nan".into());
    let statistics = Json::Dict(vec![(
        "z-score".into(),
        Json::Dict(vec![
            ("statistics".into(), nan()),
            ("score".into(), c.score.map_or_else(nan, Json::Float)),
            ("pvalue".into(), nan()),
            ("sample_score".into(), c.score.map_or_else(nan, Json::Float)),
        ]),
    )]);
    Json::Dict(vec![
        ("chromosome".into(), Json::Str(c.chromosome.clone())),
        ("start".into(), Json::Int(c.start)),
        ("end".into(), Json::Int(c.end)),
        ("size".into(), Json::Int(c.size)),
        ("bin_size".into(), Json::Int(1000)),
        ("pos".into(), Json::List(c.pos.iter().map(|&p| Json::Int(p)).collect())),
        ("cov".into(), Json::List(c.cov.iter().map(|&d| Json::Int(d)).collect())),
        ("normalization_value".into(), Json::Float(1.0)),
        ("type".into(), Json::Str("LOH".into())),
        ("cn_status".into(), Json::Int(c.cn_status)),
        ("gt".into(), Json::Str(if c.loh_pass { "1/1" } else { "0/0" }.into())),
        ("median_cov_norm".into(), Json::Float(c.median_cov_norm())),
        ("median_raw_cov".into(), Json::Float(f64::NAN)),
        ("size_kb".into(), Json::Int(c.size_kb)),
        ("is_init".into(), Json::Bool(true)),
        ("filter".into(), Json::Str(c.filter.clone())),
        ("quality".into(), c.quality.map_or_else(|| Json::Str(".".into()), Json::Float)),
        ("id".into(), Json::Str(c.id.clone())),
        ("het_score".into(), Json::Float(c.het_score)),
        ("statistics".into(), statistics),
        ("support_cnv_calls".into(), Json::Dict(vec![])),
        ("sample_origin".into(), Json::Str(String::new())),
        ("merged_sample_references".into(), Json::List(vec![])),
        ("sv_support".into(), Json::Bool(false)),
        ("merged_sample".into(), Json::Bool(false)),
        ("aflist".into(), Json::List(c.aflist.iter().map(|&v| Json::Float(v)).collect())),
        ("quallist".into(), Json::List(c.quallist.iter().map(|&v| Json::Float(v)).collect())),
        ("median_cov".into(), Json::Float(c.median_cov)),
        ("copy_number".into(), Json::Float(c.copy_number)),
        ("pretty_size".into(), Json::Str(String::new())),
        ("number_snvs".into(), Json::Int(c.number_snvs)),
        ("gt_alt_prop".into(), Json::Float(c.gt_alt_prop)),
        ("gt_alt_count".into(), Json::Int(c.gt_alt_count)),
        ("gt_all_count".into(), Json::Int(c.gt_all_count)),
        ("loh_pass".into(), Json::Bool(c.loh_pass)),
        ("af".into(), Json::Float(f64::NAN)),
    ])
}

fn loh_candidates(lists: &[(String, Vec<LohCandidate>)]) -> Json {
    Json::Dict(lists.iter().map(|(c, v)| (c.clone(), Json::List(v.iter().map(loh_candidate).collect()))).collect())
}

fn candidates(arena: &Arena, lists: &[(String, Vec<CandidateId>)], sample_id: &str) -> Json {
    Json::Dict(
        lists
            .iter()
            .map(|(contig, ids)| {
                (contig.clone(), Json::List(ids.iter().map(|&id| candidate(arena.get(id), sample_id)).collect()))
            })
            .collect(),
    )
}

/// The `.spc` content. `merged_sample_references` lists are in merge order; Spectre writes
/// them in Python set order, which depends on per-process hash randomization.
pub fn document(
    sample_id: &str,
    genome: &Genome,
    arena: &Arena,
    raw: &[(String, Vec<CandidateId>)],
    refined: &[(String, Vec<CandidateId>)],
    loh: Option<&LohCalls>,
    m: &Metrics<'_>,
) -> Json {
    let names = genome.contigs.iter().map(|c| Json::Str(c.name.clone())).collect();
    let lengths = genome.contigs.iter().map(|c| (c.name.clone(), Json::Int(c.len as i64))).collect();
    let ploidy = genome.contigs.iter().map(|c| (c.name.clone(), Json::Int(c.ploidy_multiplier))).collect();
    let cov: Vec<f64> = m.sample.windows.iter().map(|w| w.coverage).collect();
    Json::Dict(vec![
        (
            "metadata".into(),
            Json::Dict(vec![
                ("source".into(), Json::Str("spectre".into())),
                ("spectre_version".into(), Json::Str("0.2.alpha".into())),
                ("bin_size".into(), Json::Int(m.bin_size)),
                ("min_cn_len".into(), Json::Int(m.min_cnv_len)),
                ("normalization_value".into(), Json::Float(m.normalization_value)),
                ("sample_id".into(), Json::Str(sample_id.into())),
            ]),
        ),
        (
            "genome_info".into(),
            Json::Dict(vec![
                ("chromosomes".into(), Json::List(names)),
                ("chr_lengths".into(), Json::Dict(lengths)),
                ("chr_ploidy".into(), Json::Dict(ploidy)),
            ]),
        ),
        ("raw_cnvs".into(), candidates(arena, raw, sample_id)),
        ("refined_cnvs".into(), candidates(arena, refined, sample_id)),
        ("refined_loh_cnvs".into(), loh.map_or_else(|| Json::Dict(vec![]), |l| loh_candidates(&l.passing))),
        ("raw_loh_cnvs".into(), loh.map_or_else(|| Json::Dict(vec![]), |l| loh_candidates(&l.failing))),
        (
            "analysis_metrics".into(),
            Json::Dict(vec![
                ("min_chr_length".into(), Json::Int(100_000)),
                ("max_std_outlier_rm".into(), Json::Int(m.max_std_outlier_rm)),
                ("mosdepth_cov_genome_chr_diff".into(), Json::Float(m.mosdepth_cov_genome_chr_diff)),
                ("lower_2n_threshold".into(), Json::Float(m.thresholds.del.value())),
                ("upper_2n_threshold".into(), Json::Float(m.thresholds.dup.value())),
                ("strict_del_threshold".into(), Json::Float(m.strict.del.value())),
                ("strict_dup_threshold".into(), Json::Float(m.strict.dup.value())),
                ("cov_diff_threshold".into(), Json::Float(m.cov_diff_threshold)),
                ("dist_proportion".into(), Json::Float(m.dist_proportion)),
                ("dist_min_overwrite".into(), Json::Int(m.dist_min_overwrite)),
                ("candidate_final_threshold".into(), Json::Int(m.min_cnv_len)),
                ("genome_mean".into(), Json::Float(m.sample.mean)),
                ("genome_sd".into(), Json::Float(m.sample.sd)),
                ("genome_var".into(), Json::Float(pynum::sum::series_var(&cov, 0))),
            ]),
        ),
    ])
}

/// Writes `doc` gzipped, as `json.dump(doc, gzip.open(path, "wt"), indent="\t")`.
///
/// # Errors
///
/// If the file cannot be written.
pub fn write(path: &Path, doc: &Json) -> Result<(), Error> {
    let mut text = String::new();
    doc.write_indented("\t", &mut text);
    let file = std::fs::File::create(path).map_err(Error::io(path))?;
    let mut gz = flate2::write::GzEncoder::new(std::io::BufWriter::new(file), flate2::Compression::default());
    gz.write_all(text.as_bytes()).map_err(Error::io(path))?;
    gz.finish().and_then(|mut w| w.flush()).map_err(Error::io(path))?;
    Ok(())
}

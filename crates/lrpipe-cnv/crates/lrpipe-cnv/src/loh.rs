//! Loss of heterozygosity: runs of homozygous ALT SNVs (spec §3.13; snv_analysis.py:255-543,
//! classes/loh_candidate.py).

use pynum::{random::Random, round::np_round, sum};

use crate::call::new_ids;
use crate::snv::{Af, SnvVcf};

/// LoH filter settings (`--loh-min-snv-perkb`, `--loh-min-snv-total`, `--loh-min-region-size`).
#[derive(Clone, Debug)]
pub struct LohOptions {
    pub min_snv_perkb: i64,
    pub min_snv_total: i64,
    pub min_region_size: f64,
}

impl Default for LohOptions {
    fn default() -> Self {
        LohOptions { min_snv_perkb: 10, min_snv_total: 1000, min_region_size: 1_000_000.0 }
    }
}

impl LohOptions {
    /// Spectre's floors on these settings (snv_analysis.py:50-58).
    fn floored(&self) -> Self {
        LohOptions {
            min_snv_perkb: self.min_snv_perkb.max(1),
            min_snv_total: self.min_snv_total.max(10),
            min_region_size: self.min_region_size.max(10_000.0),
        }
    }
}

/// A LoH candidate (`LoHCandidate`). Coordinates are 1-based SNV positions.
#[derive(Clone, Debug, PartialEq)]
pub struct LohCandidate {
    pub chromosome: String,
    pub id: String,
    pub start: i64,
    pub end: i64,
    pub size: i64,
    pub size_kb: i64,
    pub pos: Vec<i64>,
    /// SNV depths.
    pub cov: Vec<i64>,
    pub aflist: Vec<f64>,
    pub quallist: Vec<f64>,
    pub number_snvs: i64,
    pub gt_alt_count: i64,
    pub gt_all_count: i64,
    pub median_cov: f64,
    pub copy_number: f64,
    pub cn_status: i64,
    pub gt_alt_prop: f64,
    pub filter: String,
    pub loh_pass: bool,
    /// For passing candidates: mean AF and mean QUAL (rounded), and the score used as GQ.
    pub het_score: f64,
    pub quality: Option<f64>,
    pub score: Option<f64>,
}

impl LohCandidate {
    fn start(chromosome: &str, pos: i64, dp: i64, af: f64, qual: f64) -> Self {
        LohCandidate {
            chromosome: chromosome.to_owned(),
            id: String::new(),
            start: pos,
            end: pos,
            size: 1,
            size_kb: 0,
            pos: vec![pos],
            cov: vec![dp],
            aflist: vec![af],
            quallist: vec![qual],
            number_snvs: 1,
            gt_alt_count: 1,
            gt_all_count: 1,
            median_cov: f64::NAN,
            copy_number: f64::NAN,
            cn_status: 0,
            gt_alt_prop: 0.0,
            filter: ".".into(),
            loh_pass: false,
            het_score: 0.0,
            quality: None,
            score: None,
        }
    }

    fn push(&mut self, pos: i64, dp: i64, af: f64, qual: f64) {
        self.pos.push(pos);
        self.end = pos;
        self.size = self.end - self.start;
        self.size_kb = (self.size as f64 / 1000.0) as i64;
        self.number_snvs += 1;
        self.aflist.push(af);
        self.quallist.push(qual);
        self.gt_alt_count += 1;
        self.gt_all_count += 1;
        self.cov.push(dp);
    }

    /// `set_copy_number` with Spectre's SNV ploidy, always 2.
    fn set_copy_number(&mut self, cn_neutral: f64) {
        let cov: Vec<f64> = self.cov.iter().map(|&d| d as f64).collect();
        self.median_cov = sum::nanmedian(&cov);
        let cn = self.median_cov / cn_neutral * 2.0;
        self.copy_number = np_round(cn, 2);
        self.cn_status = np_round(cn, 0) as i64;
    }

    /// `median_cov_norm`, which LoH candidates set to their copy number.
    pub fn median_cov_norm(&self) -> f64 {
        self.copy_number
    }
}

const MIN_SNV_QUAL: f64 = 5.0;
const MIN_GT_ALT_PROP: f64 = 0.85;
const MIN_CHR_LENGTH: f64 = 1e6;
const MAX_SKIP: i64 = 2;
const MIN_SAVE: usize = 10;
const MIN_ALT_SITES_CAND: f64 = 0.70;

/// `compute_gt_alt_prop`: the share of records overlapping the region that are homozygous ALT
/// with QUAL > 5 (all records count in the denominator).
fn gt_alt_prop(vcf: &SnvVcf, contig: &str, start: i64, end: i64) -> f64 {
    let (mut alt, mut all) = (0u64, 0u64);
    for r in vcf.region(contig, start, end) {
        all += 1;
        if r.is_hom_alt() && r.qual > MIN_SNV_QUAL {
            alt += 1;
        }
    }
    alt as f64 / all as f64
}

/// LoH calls: passing and failing candidates per `only_chr` contig (`loh_pass_and_non_pass`).
#[derive(Clone, Debug, Default)]
pub struct LohCalls {
    pub passing: Vec<(String, Vec<LohCandidate>)>,
    pub failing: Vec<(String, Vec<LohCandidate>)>,
}

/// `SNVAnalysis.loh`: scan, merge neighbours, filter. IDs are drawn from `rng` (Python's
/// `random` as the window sample left it) avoiding `existing_ids`.
pub fn call_loh(
    vcf: &SnvVcf,
    only_chr: &[String],
    cn_neutral: f64,
    opts: &LohOptions,
    rng: &mut Random,
    existing_ids: &mut std::collections::HashSet<String>,
) -> LohCalls {
    let opts = opts.floored();
    // candidate_loh, in VCF header contig order.
    let mut by_contig: Vec<(String, Vec<LohCandidate>)> = Vec::new();
    for (contig, len) in &vcf.header_contigs {
        if !(len.is_some_and(|l| l as f64 > MIN_CHR_LENGTH) && only_chr.contains(contig)) {
            continue;
        }
        let mut found = Vec::new();
        let mut cand: Option<LohCandidate> = None;
        let (mut not_loh, mut loh_per, mut sites) = (0i64, 0i64, 0i64);
        let mut save = |mut c: LohCandidate, loh_per: i64, sites: i64, found: &mut Vec<LohCandidate>| {
            c.set_copy_number(cn_neutral);
            let id = new_ids(rng, existing_ids, 1).pop().expect("one ID");
            c.id = format!("Spectre.LOH.{id}");
            c.gt_alt_prop = loh_per as f64 / sites as f64;
            found.push(c);
        };
        for r in vcf.records(contig) {
            if r.qual <= MIN_SNV_QUAL {
                continue;
            }
            let Af::Value(af) = r.af else { continue };
            if r.is_hom_alt() {
                loh_per += 1;
                sites += 1;
                match &mut cand {
                    None => cand = Some(LohCandidate::start(contig, r.pos, r.dp, af, r.qual)),
                    Some(c) => {
                        not_loh = 0;
                        c.push(r.pos, r.dp, af, r.qual);
                    }
                }
            } else {
                not_loh += 1;
                if not_loh >= MAX_SKIP {
                    if let Some(c) = cand.take()
                        && c.pos.len() >= MIN_SAVE
                        && loh_per as f64 >= sites as f64 * MIN_ALT_SITES_CAND
                    {
                        save(c, loh_per, sites, &mut found);
                    }
                    sites = 0;
                    loh_per = 0;
                    not_loh = 0;
                } else {
                    sites += 1;
                }
            }
        }
        if let Some(c) = cand.take()
            && c.pos.len() > MIN_SAVE
        {
            save(c, loh_per, sites, &mut found);
        }
        by_contig.push((contig.clone(), found));
    }

    // candidate_loh_neighboring: per contig, merge passes until one merges nothing.
    let mut merged: Vec<LohCandidate> = Vec::new();
    for (_, list) in by_contig {
        let mut list = list;
        for _ in 0..500 {
            let (next, merged_any) = merge_pass(&list, vcf, cn_neutral);
            list = next;
            if !merged_any {
                break;
            }
        }
        merged.extend(list);
    }

    // candidate_loh_filtering.
    for c in &mut merged {
        let perkb = c.size_kb * opts.min_snv_perkb;
        if c.size as f64 > opts.min_region_size {
            if c.number_snvs > perkb || c.number_snvs > opts.min_snv_total {
                if c.gt_alt_prop >= MIN_GT_ALT_PROP {
                    c.loh_pass = true;
                    c.filter = "PASS".into();
                    let mean_af = sum::nanmean(&c.aflist);
                    c.het_score = np_round(mean_af, 2);
                    c.quality = Some(np_round(sum::nanmean(&c.quallist), 2));
                    // get_score_sigmoid: Python's min keeps 60 (an int) when it is smaller.
                    let s = 60.0 / (1.0 + (-(c.size as f64 / 1e6 - 1.0)).exp());
                    c.score = Some(if 60.0 < s { 60.0 } else { s } * mean_af);
                } else {
                    c.filter = format!("ALT GT propotion: {}|0.85", pynum::repr::py_repr(c.gt_alt_prop));
                }
            } else {
                c.filter = format!("Number of SVs: {}|{perkb}/kb|{}", c.number_snvs, c.size_kb);
            }
        } else {
            c.filter = format!("LoH region size kb: {}|{}|{perkb}/kb", c.size_kb, c.number_snvs);
        }
    }

    let mut out = LohCalls {
        passing: only_chr.iter().map(|c| (c.clone(), Vec::new())).collect(),
        failing: only_chr.iter().map(|c| (c.clone(), Vec::new())).collect(),
    };
    for c in merged {
        let side = if c.loh_pass { &mut out.passing } else { &mut out.failing };
        if let Some((_, list)) = side.iter_mut().find(|(name, _)| *name == c.chromosome) {
            list.push(c);
        }
    }
    out
}

/// `merge_candidates`: one pass joining neighbours that are close (for their size) and whose
/// span is mostly homozygous ALT.
fn merge_pass(list: &[LohCandidate], vcf: &SnvVcf, cn_neutral: f64) -> (Vec<LohCandidate>, bool) {
    let mut out = Vec::new();
    let mut combine: Vec<&LohCandidate> = Vec::new();
    let mut merged_any = false;
    for (k, cand) in list.iter().enumerate() {
        let in_combine = |c: &LohCandidate, combine: &[&LohCandidate]| combine.iter().any(|x| x.id == c.id);
        let Some(next) = list.get(k + 1) else {
            if combine.is_empty() {
                out.push(cand.clone());
            } else if in_combine(cand, &combine) {
                out.push(combined(&combine, vcf, cn_neutral));
            }
            break;
        };
        let near = (next.start - cand.end) as f64 <= merge_dist(cand.size, next.size);
        let prop = if near { gt_alt_prop(vcf, &cand.chromosome, cand.start, next.end) } else { 0.0 };
        if near && prop >= MIN_GT_ALT_PROP {
            merged_any = true;
            if combine.is_empty() {
                combine = vec![cand, next];
            } else {
                combine.push(next);
            }
        } else {
            if combine.is_empty() {
                out.push(cand.clone());
            } else if in_combine(cand, &combine) {
                out.push(combined(&combine, vcf, cn_neutral));
            }
            combine.clear();
        }
    }
    (out, merged_any)
}

/// `merge_dist`: how far apart two candidates may be, given their sizes.
fn merge_dist(size1: i64, size2: i64) -> f64 {
    let m = sum::mean(&[size1 as f64, size2 as f64]);
    let l = m.log10();
    m * (1.0 / l.powf(1.75) + (6.0 / l) / 100.0)
}

/// `LoHCandidate.combine`: the first candidate's ID, all their SNVs, the span of the first to
/// the last, and the homozygous ALT share recomputed over that span.
fn combined(parts: &[&LohCandidate], vcf: &SnvVcf, cn_neutral: f64) -> LohCandidate {
    let (first, last) = (parts[0], parts[parts.len() - 1]);
    let mut c = LohCandidate::start(&first.chromosome, first.start, 0, 0.0, 0.0);
    c.id = first.id.clone();
    c.pos = parts.iter().flat_map(|p| p.pos.iter().copied()).collect();
    c.cov = parts.iter().flat_map(|p| p.cov.iter().copied()).collect();
    c.aflist = parts.iter().flat_map(|p| p.aflist.iter().copied()).collect();
    c.quallist = parts.iter().flat_map(|p| p.quallist.iter().copied()).collect();
    c.number_snvs = parts.iter().map(|p| p.number_snvs).sum();
    c.end = last.end;
    c.size = c.end - c.start;
    c.size_kb = (c.size as f64 / 1000.0) as i64;
    // A new LoHCandidate: its counters start at 0 and combine does not set them.
    c.gt_alt_count = 0;
    c.gt_all_count = 0;
    c.set_copy_number(cn_neutral);
    c.gt_alt_prop = gt_alt_prop(vcf, &c.chromosome, c.start, c.end);
    c
}

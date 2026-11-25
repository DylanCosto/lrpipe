//! Raw candidates: runs of windows beyond the thresholds (spec §3.5; call_cnv_coverage.py,
//! analysis.py:311-331).

use pynum::random::Random;
use std::collections::HashSet;

use crate::candidate::{Arena, Candidate, CandidateId, CnvType};
use crate::normalize::{ContigCoverage, Normalized};
use crate::thresholds::Thresholds;

/// Runs shorter than this many windows are not candidates (`min_run`).
const MIN_RUN: usize = 10;

/// The characters of candidate IDs (`string.ascii_uppercase + string.digits`).
const ID_CHARS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

/// Raw candidates per analysed contig, in contig order.
#[derive(Clone, Debug, Default)]
pub struct RawCalls {
    pub contigs: Vec<(String, Vec<CandidateId>)>,
    /// Every ID drawn (`existing_cnv_ids`); later IDs must differ from these.
    pub existing_ids: HashSet<String>,
}

/// `call_cnv_coverage`: segments each contig and gives every candidate an ID drawn from
/// `rng`, which continues from the window sample's draws.
pub fn call(norm: &Normalized, thresholds: &Thresholds, rng: &mut Random, arena: &mut Arena) -> RawCalls {
    let mut out = RawCalls::default();
    for c in &norm.contigs {
        let candidates = segment(c, norm.bin_size, thresholds);
        let ids = new_ids(rng, &mut out.existing_ids, candidates.len());
        let ids: Vec<CandidateId> = candidates
            .into_iter()
            .zip(ids)
            .map(|(mut cand, id)| {
                cand.id = format!("Spectre.{}.{id}", cand.cnv_type.as_str());
                arena.push(cand)
            })
            .collect();
        out.contigs.push((c.name.clone(), ids));
    }
    out
}

/// `CNV_ID.n_id_generator`: `n` new 8-character IDs, each redrawn until unseen.
pub fn new_ids(rng: &mut Random, existing: &mut HashSet<String>, n: usize) -> Vec<String> {
    (0..n)
        .map(|_| {
            loop {
                let id: String = (0..8).map(|_| char::from(*rng.choice(ID_CHARS))).collect();
                if existing.insert(id.clone()) {
                    break id;
                }
            }
        })
        .collect()
}

/// `CNVCall.cnv_coverage` for one contig.
fn segment(c: &ContigCoverage, bin_size: i64, t: &Thresholds) -> Vec<Candidate> {
    let (del, dup) = (t.del.value(), t.dup.value());
    let mut out = Vec::new();
    let (mut pos_list, mut cov_list): (Vec<i64>, Vec<f64>) = (Vec::new(), Vec::new());
    let (mut run_start, mut run_current_pos) = (0i64, 0i64);
    let mut current: Option<CnvType> = None;
    let mut emit = |pos: &mut Vec<i64>, cov: &mut Vec<f64>, current: Option<CnvType>| {
        if pos.len() >= MIN_RUN {
            let t = current.expect("a run has a type");
            out.push(Candidate::new(&c.name, t, bin_size, std::mem::take(pos), std::mem::take(cov)));
        }
    };
    for (&v, &p) in c.normalized_ploidy.iter().zip(&c.positions) {
        if v.is_nan() {
            continue;
        }
        if v > dup || v < del {
            let t = if v > dup { CnvType::Dup } else { CnvType::Del };
            if run_start == 0 {
                // Quirk Q-RUNSTART: the window that starts the first run is not kept.
                current = Some(t);
                run_start = p;
            } else if p - run_current_pos < bin_size + 1 && current == Some(t) {
                pos_list.push(p);
                cov_list.push(v);
            } else {
                emit(&mut pos_list, &mut cov_list, current);
                pos_list = vec![p];
                cov_list = vec![v];
                run_start = p;
                current = Some(t);
            }
        } else {
            emit(&mut pos_list, &mut cov_list, current);
            pos_list.clear();
            cov_list.clear();
            run_start = p;
            current = None;
        }
        run_current_pos = p;
    }
    emit(&mut pos_list, &mut cov_list, current);
    out
}

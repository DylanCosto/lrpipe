//! CNV candidates (`CNVCandidate`, classes/cnv_candidate.py).
//!
//! Spectre's candidate lists share candidate objects: the raw list written to the `.spc`
//! holds the same objects that merging and refinement change in place (quirk Q-ALIAS). An
//! [`Arena`] owns every candidate and lists refer to them by [`CandidateId`], which keeps
//! that sharing.

use pynum::{round::np_round, sum};

use crate::breakpoints::Breakpoints;
use crate::score::ZScore;

/// Deletion or duplication.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CnvType {
    Del,
    Dup,
}

impl CnvType {
    pub fn as_str(self) -> &'static str {
        match self {
            CnvType::Del => "DEL",
            CnvType::Dup => "DUP",
        }
    }
}

/// A candidate's index in its [`Arena`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CandidateId(pub usize);

/// A CNV candidate: a run of windows, with the fields Spectre keeps on it.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub chromosome: String,
    pub cnv_type: CnvType,
    /// `Spectre.<TYPE>.<8 characters>`.
    pub id: String,
    pub bin_size: i64,
    /// Window starts, in order.
    pub pos: Vec<i64>,
    /// Normalized, ploidy-scaled depth of each window in `pos` (NaN after refinement where
    /// a window was masked).
    pub cov: Vec<f64>,
    /// First and last window start (not window end).
    pub start: i64,
    pub end: i64,
    pub size: i64,
    pub size_kb: i64,
    pub median_cov_norm: f64,
    pub median_raw_cov: f64,
    /// Copy number from the raw run's median; never updated (bug B3).
    pub cn_status: i64,
    pub gt: &'static str,
    /// IDs of the candidates merged into this one (a Python set: order is not meaningful).
    pub merged_sample_references: Vec<String>,
    /// `statistics["z-score"]`, set by the quality scoring.
    pub statistics: Option<ZScore>,
    /// Base-pair breakpoints from an SV call (an extension; never set in Spectre's method).
    pub breakpoints: Option<Breakpoints>,
}

impl Candidate {
    /// `CNVCandidate.push_candidates`.
    pub fn new(chromosome: &str, cnv_type: CnvType, bin_size: i64, pos: Vec<i64>, cov: Vec<f64>) -> Self {
        let median = sum::nanmedian(&cov);
        // `set_copy_number_status`: an infinite median gives CN 2 and median 2.0.
        let (cn_status, median_cov_norm) =
            if median == f64::INFINITY { (2, 2.0) } else { (np_round(median, 0) as i64, median) };
        let mut c = Candidate {
            chromosome: chromosome.to_owned(),
            cnv_type,
            id: String::new(),
            bin_size,
            start: 0,
            end: 0,
            size: 0,
            size_kb: 0,
            pos,
            cov,
            median_cov_norm,
            median_raw_cov: median_cov_norm,
            cn_status,
            gt: gt(cn_status),
            merged_sample_references: Vec::new(),
            statistics: None,
            breakpoints: None,
        };
        c.update_span();
        c
    }

    fn update_span(&mut self) {
        self.start = self.pos[0];
        self.end = *self.pos.last().expect("a candidate has windows");
        self.size = self.end - self.start;
        self.size_kb = (self.size as f64 / self.bin_size as f64) as i64;
    }

    /// `add_candidates`: appends (or prepends) another candidate's windows and records its ID
    /// and the IDs merged into it.
    pub fn absorb(&mut self, pos: &[i64], cov: &[f64], id: &str, refs: &[String], side: Side) {
        match side {
            Side::Front => {
                self.pos.splice(0..0, pos.iter().copied());
                self.cov.splice(0..0, cov.iter().copied());
            }
            Side::Back => {
                self.pos.extend_from_slice(pos);
                self.cov.extend_from_slice(cov);
            }
        }
        self.update_span();
        self.median_cov_norm = sum::nanmedian(&self.cov);
        self.median_raw_cov = self.median_cov_norm;
        for r in std::iter::once(id).chain(refs.iter().map(String::as_str)) {
            if !self.merged_sample_references.iter().any(|x| x == r) {
                self.merged_sample_references.push(r.to_owned());
            }
        }
    }

    /// `median_coverage_candidates_merged`: the median of `cov`, 0.0 when empty.
    pub fn update_median(&mut self) {
        self.median_cov_norm = if self.cov.is_empty() { 0.0 } else { sum::nanmedian(&self.cov) };
    }
}

/// Where [`Candidate::absorb`] puts the other candidate's windows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Front,
    Back,
}

/// `set_gt`.
/// Spectre's genotype for a copy number on a diploid contig.
pub fn gt(cn_status: i64) -> &'static str {
    match cn_status {
        0 | 4 => "1/1",
        1 | 3 => "0/1",
        2 => "0/0",
        _ => "./.",
    }
}

/// Owns every candidate of a run.
#[derive(Clone, Debug, Default)]
pub struct Arena {
    candidates: Vec<Candidate>,
}

impl Arena {
    pub fn push(&mut self, c: Candidate) -> CandidateId {
        self.candidates.push(c);
        CandidateId(self.candidates.len() - 1)
    }

    pub fn get(&self, id: CandidateId) -> &Candidate {
        &self.candidates[id.0]
    }

    pub fn get_mut(&mut self, id: CandidateId) -> &mut Candidate {
        &mut self.candidates[id.0]
    }

    /// `a` mutably and `b` shared, for merging `b` into `a`.
    ///
    /// # Panics
    ///
    /// If `a == b`.
    pub fn pair_mut(&mut self, a: CandidateId, b: CandidateId) -> (&mut Candidate, &Candidate) {
        assert_ne!(a, b, "merging a candidate into itself");
        if a.0 < b.0 {
            let (lo, hi) = self.candidates.split_at_mut(b.0);
            (&mut lo[a.0], &hi[0])
        } else {
            let (lo, hi) = self.candidates.split_at_mut(a.0);
            (&mut hi[0], &lo[b.0])
        }
    }
}

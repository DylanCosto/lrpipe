//! Port of `whatshap/merge.py` (`phase --merge-reads`).

use std::collections::{BTreeMap, HashMap};

use anyhow::{Result, bail};

use crate::core::read::{Read, ReadSet};

/// `ReadMerger` with the parameters `phase` uses.
pub struct ReadMerger {
    /// Probability that a nucleotide is wrong
    pub error_rate: f64,
    /// Maximum error rate of an edge of the merge graph
    pub max_error_rate: f64,
    /// Threshold on the ratio of the probabilities that two reads come from
    /// the same or from different haplotypes
    pub positive_threshold: f64,
    /// The same, for evidence that two reads come from different haplotypes
    pub negative_threshold: f64,
}

impl Default for ReadMerger {
    fn default() -> Self {
        ReadMerger {
            error_rate: 0.15,
            max_error_rate: 0.25,
            positive_threshold: 1_000_000.0,
            negative_threshold: 1000.0,
        }
    }
}

/// A read in the overlap queue: genomic position of its first variant,
/// that position plus its number of variants (as in WhatsHap, which mixes
/// the two), and its alleles.
struct Queued {
    begin: i64,
    end: i64,
    alleles: Vec<i32>,
}

/// `eval_overlap`: matches and mismatches of `n2`'s alleles against those of
/// `n1` from index `n2.begin - n1.begin` on (a Python slice).
fn eval_overlap(n1: &Queued, n2: &Queued) -> (usize, usize) {
    let hang = n2.begin - n1.begin;
    let len = n1.alleles.len() as i64;
    let from = if hang < 0 { (len + hang).max(0) } else { hang.min(len) } as usize;
    let (mut matches, mut mismatches) = (0, 0);
    for (a, b) in n1.alleles[from..].iter().zip(&n2.alleles) {
        if a == b {
            matches += 1;
        } else {
            mismatches += 1;
        }
    }
    (matches, mismatches)
}

/// Union-find over read indices.
fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

impl ReadMerger {
    /// `ReadMerger.merge`: merges reads that are likely to come from the
    /// same haplotype into super reads.
    ///
    /// WhatsHap also builds a graph of evidence that two reads come from
    /// different haplotypes and cuts the merge graph between such reads. An
    /// edge of that graph needs `mismatches - matches >= thr_neg_diff`, but it
    /// is only considered for pairs with `matches - mismatches >= thr_diff`,
    /// so with positive thresholds it is always empty and the super reads are
    /// the connected components of the merge graph.
    pub fn merge(&self, readset: &ReadSet) -> Result<ReadSet> {
        log::info!(
            "Merging {} reads with error rate {:.2}, maximum error rate {:.2}, positive threshold {} and negative threshold {} ...",
            readset.len(),
            self.error_rate,
            self.max_error_rate,
            self.positive_threshold,
            self.negative_threshold
        );
        let e = self.error_rate;
        let base = ((1.0 - e) / (e / 3.0)).ln();
        let thr_diff = 1 + (self.positive_threshold.ln() / base) as i64;
        let thr_neg_diff = 1 + (self.negative_threshold.ln() / base) as i64;
        assert!(thr_diff > 0 && thr_neg_diff > 0, "the evidence against merging is not modelled");
        log::debug!("Thr. Diff.: {thr_diff} - Thr. Neg. Diff.: {thr_neg_diff}");

        let n = readset.len();
        let mut parent: Vec<usize> = (0..n).collect();
        // dict of queued reads in insertion order
        let mut queue: Vec<(usize, Queued)> = Vec::new();
        for (i, read) in readset.iter().enumerate() {
            let mut alleles = Vec::with_capacity(read.len());
            for v in &read.variants {
                if v.allele != 0 && v.allele != 1 {
                    bail!("Read merging needs biallelic variants (read {} has allele {})", read.name, v.allele);
                }
                alleles.push(v.allele);
            }
            let begin = read.first_position() as i64;
            let end = begin + alleles.len() as i64;
            queue.push((i, Queued { begin, end, alleles }));
            queue.retain(|(_, q)| q.end > begin);
            let current = &queue.last().expect("read i was just queued").1;
            for (j, q) in &queue[..queue.len() - 1] {
                let (matches, mismatches) = eval_overlap(q, current);
                let total = matches + mismatches;
                if total as i64 >= thr_neg_diff
                    && matches.min(mismatches) as f64 / total as f64 <= self.max_error_rate
                    && matches as i64 - mismatches as i64 >= thr_diff
                {
                    let (a, b) = (find(&mut parent, *j), find(&mut parent, i));
                    parent[a.max(b)] = a.min(b);
                }
            }
        }

        // Representative (smallest index) of each component with more than
        // one read, and the summed qualities per position and allele
        let roots: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
        let mut size: HashMap<usize, usize> = HashMap::new();
        for &r in &roots {
            *size.entry(r).or_default() += 1;
        }
        let mut superreads: HashMap<usize, BTreeMap<i32, [i64; 2]>> = HashMap::new();
        for (i, read) in readset.iter().enumerate() {
            if size[&roots[i]] > 1 {
                let sums = superreads.entry(roots[i]).or_default();
                for v in &read.variants {
                    sums.entry(v.position).or_insert([0, 0])[v.allele as usize] += v.quality as i64;
                }
            }
        }
        let mut merged = ReadSet::new();
        for (i, read) in readset.iter().enumerate() {
            let mut out = Read::new(&format!("read{i}"), 0, 0, 0);
            if size[&roots[i]] > 1 {
                if roots[i] != i {
                    continue;
                }
                for (&position, z) in &superreads[&i] {
                    let allele = if z[0] >= z[1] { 0 } else { 1 };
                    out.add_variant(position, allele, (z[1] - z[0]).abs() as i32);
                }
            } else {
                for v in &read.variants {
                    out.add_variant(v.position, v.allele, v.quality as i32);
                }
            }
            merged.add(out)?;
        }
        log::info!("... after merging: merged {} reads into {} reads", readset.len(), merged.len());
        Ok(merged)
    }
}

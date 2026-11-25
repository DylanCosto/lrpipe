//! Insertion consensus sequence (`consensus.novel_from_reads`).

use std::collections::{HashMap, HashSet};

use crate::pymath;

const CONSENSUS_MIN: usize = 2;
const MINSPAN: f64 = 0.2;
const MINALNS: f64 = 0.25;
const MINIDENT: f64 = 0.5;
const MINIDENT_ABS: usize = 5;
const MINBESTDIFF: usize = 3;

/// k-mer start positions as `iter_kmers` yields them (the last k-mer is never included).
fn kmer_starts(len: usize, klen: usize, skip: usize) -> impl Iterator<Item = usize> {
    (0..len.saturating_sub(klen)).step_by(skip.max(1))
}

/// Corrects `best` by pseudo-aligning the other insertion sequences to it through unique k-mers.
pub fn novel_from_reads(best: &[u8], others: &[&[u8]], klen: usize, skip: usize, skip_repetitive: usize) -> Vec<u8> {
    let maxshift = klen as i64;
    let mut anchors: HashMap<&[u8], usize> = HashMap::new();
    let mut taboo: HashSet<&[u8]> = HashSet::new();
    for i in kmer_starts(best.len(), klen, skip_repetitive) {
        let kmer = &best[i..i + klen];
        if taboo.contains(kmer) {
            continue;
        }
        if anchors.remove(kmer).is_some() {
            taboo.insert(kmer);
            continue;
        }
        anchors.insert(kmer, i);
    }

    let n = best.len();
    let mut alignments: Vec<Vec<u8>> = Vec::new();
    for seq in others {
        let mut last: Option<(usize, usize)> = None;
        let mut conseq: Vec<u8> = Vec::new();
        let mut span = 0usize;
        for j in kmer_starts(seq.len(), klen, skip) {
            let Some(&i) = anchors.get(&seq[j..j + klen]) else { continue };
            if (i as i64 - j as i64).abs() > maxshift {
                continue;
            }
            match last {
                Some((li, _)) if i <= li => continue,
                None => {
                    if j > 0 {
                        conseq = vec![b'-'; i];
                    }
                }
                Some((li, lj)) => {
                    let fwd_i = (i - li) as i64;
                    let mut fwd_j = (j - lj) as i64;
                    if conseq.len() as i64 + fwd_j > n as i64 {
                        fwd_j = n as i64 - conseq.len() as i64;
                    }
                    if fwd_i == fwd_j && fwd_j > 0 {
                        span += j - lj;
                        let m = (1..=(j - lj)).filter(|&l| seq[lj + l] == best[li + l]).count();
                        let ident = m as f64 / (j - lj) as f64;
                        if ident >= MINIDENT {
                            conseq.extend_from_slice(&seq[lj..j][..fwd_j as usize]);
                        } else {
                            conseq.extend(std::iter::repeat_n(b'-', fwd_j as usize));
                        }
                    } else {
                        conseq.extend(std::iter::repeat_n(b'-', fwd_j.max(0) as usize));
                    }
                }
            }
            last = Some((i, j));
        }
        if conseq.len() < n {
            conseq.resize(n, b'-');
        }
        // Keep only runs that mostly agree with the best sequence.
        let mut out: Vec<u8> = Vec::with_capacity(n);
        let mut h = 0;
        while h < n {
            if conseq[h] == b'-' {
                out.push(b'-');
                h += 1;
            } else {
                let begin = h;
                let mut ident = 0;
                while h < n && conseq[h] != b'-' {
                    ident += usize::from(best[h] == conseq[h]);
                    h += 1;
                }
                let run = &conseq[begin..h];
                if ident as f64 / run.len() as f64 > MINIDENT && ident > MINIDENT_ABS {
                    out.extend_from_slice(run);
                } else {
                    out.extend(std::iter::repeat_n(b'-', run.len()));
                }
            }
        }
        if span as f64 / n as f64 > MINSPAN {
            alignments.push(out);
        }
    }

    // Every alignment column is a base or '-', never '^' or '_'.
    let maxal = if n > 0 { (1 + alignments.len()) as f64 } else { 1.0 };
    (0..n)
        .map(|i| {
            let al: Vec<u8> = alignments.iter().map(|a| a[i]).filter(|&c| c != b'-').collect();
            if al.len() < CONSENSUS_MIN || al.len() as f64 / maxal < MINALNS {
                return best[i];
            }
            let top = pymath::most_common(std::iter::once(best[i]).chain(al));
            if top.len() > 1 && top[0].0 - top[1].0 >= MINBESTDIFF { top[0].1 } else { best[i] }
        })
        .collect()
}

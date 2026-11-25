//! Port of `preprocess/realign/debruijn_graph.cpp` (`get_consensus`): candidate haplotypes of a
//! window from a de Bruijn graph of the reference and the reads.
//!
//! The haplotypes are sorted before they are returned, so the order in which the Boost graph
//! visits vertices does not show in the result, with one exception: the 256-path cap of
//! `CandidatePaths` is checked while the breadth-first search runs, and here adjacent vertices
//! are visited in creation order (Boost orders them by `listS` node address, which is
//! allocation order in practice).

use std::collections::{BTreeSet, HashMap, VecDeque};

/// Path count above which `CandidatePaths` gives up.
const MAX_PATHS: usize = 256;
/// Smallest k tried (`KMinMaxFromReference`).
const MIN_K: usize = 10;
/// Largest k tried.
const MAX_K: usize = 101;

/// A de Bruijn graph over k-mers; vertices are indices in creation order.
struct Graph {
    k: usize,
    kmers: Vec<Vec<u8>>,
    /// Out-edges by target, with `(weight, is_ref)`.
    out: Vec<Vec<(usize, i32, bool)>>,
    alive: Vec<bool>,
    index: HashMap<Vec<u8>, usize>,
    source: usize,
    sink: usize,
}

impl Graph {
    /// `EnsureVertex`.
    fn ensure_vertex(&mut self, kmer: &[u8]) -> usize {
        if let Some(&v) = self.index.get(kmer) {
            return v;
        }
        let v = self.kmers.len();
        self.kmers.push(kmer.to_vec());
        self.out.push(Vec::new());
        self.alive.push(true);
        self.index.insert(kmer.to_vec(), v);
        v
    }

    /// `AddEdge`: one more observation of `from -> to`.
    fn add_edge(&mut self, from: usize, to: usize, is_ref: bool) {
        let edges = &mut self.out[from];
        match edges.iter_mut().find(|(t, _, _)| *t == to) {
            Some((_, weight, r)) => {
                *weight += 1;
                *r |= is_ref;
            }
            None => {
                // `setS` keeps out-edges ordered by target.
                let at = edges.partition_point(|(t, _, _)| *t < to);
                edges.insert(at, (to, 1, is_ref));
            }
        }
    }

    /// `bases.substr(i, k)`, shorter at the end of `bases`.
    fn substr(bases: &[u8], i: usize, k: usize) -> &[u8] {
        let start = i.min(bases.len());
        &bases[start..(start + k).min(bases.len())]
    }

    /// `AddKmersAndEdges(bases, start, end, is_ref)`.
    fn add_kmers_and_edges(&mut self, bases: &[u8], start: i64, end: i64, is_ref: bool) {
        if end <= 0 {
            return;
        }
        let k = self.k;
        let mut prev = self.ensure_vertex(Self::substr(bases, start as usize, k));
        for i in start + 1..=end {
            let cur = self.ensure_vertex(Self::substr(bases, i as usize, k));
            self.add_edge(prev, cur, is_ref);
            prev = cur;
        }
    }

    /// The constructor: the reference path, then every read's runs of good bases.
    fn build(reference: &[u8], reads: &[&[u8]], low_quality: &[BTreeSet<usize>], k: usize) -> Option<Self> {
        let mut g = Self {
            k,
            kmers: Vec::new(),
            out: Vec::new(),
            alive: Vec::new(),
            index: HashMap::new(),
            source: 0,
            sink: 0,
        };
        g.add_kmers_and_edges(reference, 0, reference.len() as i64 - k as i64, true);
        // `VertexForKmer` throws when the reference is a single k-mer.
        g.source = *g.index.get(&reference[..k])?;
        g.sink = *g.index.get(&reference[reference.len() - k..])?;
        let empty = BTreeSet::new();
        for (r, read) in reads.iter().enumerate() {
            let bad = low_quality.get(r).unwrap_or(&empty);
            let next_bad = |start: usize| {
                (start..read.len()).find(|&i| !b"ACGT".contains(&read[i]) || bad.contains(&i)).unwrap_or(read.len())
            };
            let stop = read.len() as i64 - k as i64;
            let mut i = 0i64;
            while i < stop {
                let bad_at = next_bad(i as usize) as i64;
                g.add_kmers_and_edges(read, i, bad_at - k as i64, false);
                i = bad_at + 1;
            }
        }
        Some(g)
    }

    /// `HasCycle`: any back edge in a depth-first search over every vertex.
    fn has_cycle(&self) -> bool {
        #[derive(Clone, Copy, PartialEq)]
        enum Color {
            White,
            Gray,
            Black,
        }
        let mut color = vec![Color::White; self.kmers.len()];
        for root in 0..self.kmers.len() {
            if !self.alive[root] || color[root] != Color::White {
                continue;
            }
            let mut stack: Vec<(usize, usize)> = vec![(root, 0)];
            color[root] = Color::Gray;
            while let Some(&mut (v, ref mut next)) = stack.last_mut() {
                if let Some(&(to, _, _)) = self.out[v].get(*next) {
                    *next += 1;
                    match color[to] {
                        Color::Gray => return true,
                        Color::White => {
                            color[to] = Color::Gray;
                            stack.push((to, 0));
                        }
                        Color::Black => {}
                    }
                } else {
                    color[v] = Color::Black;
                    stack.pop();
                }
            }
        }
        false
    }

    /// Vertices reachable from `start` along `edges`.
    fn reachable(&self, start: usize, edges: &[Vec<usize>]) -> Vec<bool> {
        let mut seen = vec![false; self.kmers.len()];
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(v) = stack.pop() {
            for &to in &edges[v] {
                if !seen[to] {
                    seen[to] = true;
                    stack.push(to);
                }
            }
        }
        seen
    }

    /// `Prune`: drop non-reference edges seen once, then vertices not on a source-to-sink path.
    fn prune(&mut self) {
        for edges in &mut self.out {
            edges.retain(|&(_, weight, is_ref)| is_ref || weight >= 2);
        }
        let forward: Vec<Vec<usize>> = self.out.iter().map(|e| e.iter().map(|&(t, _, _)| t).collect()).collect();
        let mut backward: Vec<Vec<usize>> = vec![Vec::new(); self.kmers.len()];
        for (from, targets) in forward.iter().enumerate() {
            for &to in targets {
                backward[to].push(from);
            }
        }
        let from_source = self.reachable(self.source, &forward);
        let to_sink = self.reachable(self.sink, &backward);
        for v in 0..self.kmers.len() {
            if !(from_source[v] && to_sink[v]) {
                self.alive[v] = false;
                self.out[v].clear();
            }
        }
        let alive = self.alive.clone();
        for edges in &mut self.out {
            edges.retain(|&(t, _, _)| alive[t]);
        }
    }

    /// `CandidatePaths`: breadth-first, a path ends at the sink or at a vertex without
    /// out-edges; more than 256 paths give none.
    fn candidate_paths(&self) -> Vec<Vec<usize>> {
        let mut done: Vec<Vec<usize>> = Vec::new();
        let mut queue: VecDeque<Vec<usize>> = VecDeque::from([vec![self.source]]);
        while let Some(path) = queue.pop_front() {
            if done.len() + queue.len() + 1 > MAX_PATHS {
                return Vec::new();
            }
            let last = *path.last().expect("paths are never empty");
            for &(to, _, _) in &self.out[last] {
                let mut extended = path.clone();
                extended.push(to);
                if to == self.sink || self.out[to].is_empty() {
                    done.push(extended);
                } else {
                    queue.push_back(extended);
                }
            }
        }
        done
    }

    /// `HaplotypeForPath`.
    fn haplotype(&self, path: &[usize]) -> Vec<u8> {
        let mut out: Vec<u8> = path.iter().map(|&v| self.kmers[v][0]).collect();
        if let Some(&last) = path.last() {
            out.extend_from_slice(&self.kmers[last][1..]);
        }
        out
    }
}

/// `KMinMaxFromReference`: the smallest k from 10 at which the reference has no repeated k-mer.
fn k_bounds(reference: &[u8]) -> Option<(usize, usize)> {
    let max_k = (reference.len() as i64 - 1).min(MAX_K as i64);
    if max_k < MIN_K as i64 {
        return None;
    }
    let max_k = max_k as usize;
    (MIN_K..=max_k)
        .find(|&k| {
            let mut seen = std::collections::HashSet::new();
            reference.windows(k).all(|kmer| seen.insert(kmer))
        })
        .map(|min_k| (min_k, max_k))
}

/// `DeBruijnGraph::Build(ref, reads, base_quality)`: the sorted haplotypes of the first k whose
/// graph has no cycle, or none. `low_quality[i]` holds read `i`'s positions to skip.
pub fn candidate_haplotypes(reference: &[u8], reads: &[&[u8]], low_quality: &[BTreeSet<usize>]) -> Vec<Vec<u8>> {
    let Some((min_k, max_k)) = k_bounds(reference) else { return Vec::new() };
    for k in min_k..=max_k {
        let Some(mut graph) = Graph::build(reference, reads, low_quality, k) else { return Vec::new() };
        if graph.has_cycle() {
            continue;
        }
        graph.prune();
        let mut haplotypes: Vec<Vec<u8>> = graph.candidate_paths().iter().map(|p| graph.haplotype(p)).collect();
        haplotypes.sort();
        return haplotypes;
    }
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Opt-in: `SNV_DBG_CASES` holds `get_consensus` inputs as Python passes them
    /// (`reference<TAB>reads joined by ,<TAB>low-quality positions`) and `SNV_DBG_EXPECTED`
    /// the haplotypes the original's `debruijn_graph.cpp` returns for each, joined by `,`.
    #[test]
    fn matches_the_cpp_graph() {
        let (Ok(cases), Ok(expected)) = (std::env::var("SNV_DBG_CASES"), std::env::var("SNV_DBG_EXPECTED")) else {
            return;
        };
        let cases = std::fs::read_to_string(cases).expect("cases");
        let expected = std::fs::read_to_string(expected).expect("expected");
        let mut checked = 0;
        for (case, want) in cases.lines().zip(expected.lines()) {
            let cols: Vec<&str> = case.split('\t').collect();
            let reads: Vec<&[u8]> = cols[1].split(',').map(str::as_bytes).collect();
            let low: Vec<BTreeSet<usize>> = cols[2]
                .split(',')
                .map(|s| s.split_whitespace().map(|v| v.parse().expect("position")).collect())
                .collect();
            let got: Vec<String> = candidate_haplotypes(cols[0].as_bytes(), &reads, &low)
                .into_iter()
                .map(|h| String::from_utf8(h).expect("ASCII"))
                .collect();
            assert_eq!(got.join(","), want, "case {case}");
            checked += 1;
        }
        assert!(checked > 0);
    }

    #[test]
    fn a_supported_snp_gives_two_haplotypes() {
        let reference = b"ACGTTGCAAGGCTTACCGATAGGCTAACGTTAGCCA";
        let mut alt = reference.to_vec();
        alt[18] = b'T';
        let reads: Vec<&[u8]> = vec![&alt, &alt, reference];
        let none = vec![BTreeSet::new(); reads.len()];
        let haplotypes = candidate_haplotypes(reference, &reads, &none);
        assert_eq!(haplotypes.len(), 2);
        assert!(haplotypes.contains(&reference.to_vec()));
        assert!(haplotypes.contains(&alt));
    }

    #[test]
    fn a_single_read_is_pruned() {
        let reference = b"ACGTTGCAAGGCTTACCGATAGGCTAACGTTAGCCA";
        let mut alt = reference.to_vec();
        alt[18] = b'T';
        let reads: Vec<&[u8]> = vec![&alt];
        let haplotypes = candidate_haplotypes(reference, &reads, &[BTreeSet::new()]);
        assert_eq!(haplotypes, [reference.to_vec()]);
    }
}

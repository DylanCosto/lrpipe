//! Port of `readselect.pyx`, `priorityqueue.pyx` and `coverage.py`.
//!
//! Read selection is sensitive to the iteration order of Python sets and of a
//! C++ `unordered_set` (they decide how ties in the priority queue are
//! broken), which is why the emulations from `compat` are used here.

use std::collections::{BTreeSet, HashSet};

use anyhow::{Result, bail};

use crate::compat::cxx::CxxUnorderedSetI32;
use crate::compat::pyset::PySet;
use crate::core::read::ReadSet;
use crate::graph::ComponentFinder;

type Score = [i32; 3];

/// Max-heap with "change score" support (`priorityqueue.pyx`).
struct PriorityQueue {
    heap: Vec<(Score, usize)>,
    /// Each item's index in `heap`, or [`PriorityQueue::ABSENT`]; items are read indices.
    positions: Vec<usize>,
}

impl PriorityQueue {
    const ABSENT: usize = usize::MAX;

    /// A queue for items below `items`.
    fn new(items: usize) -> Self {
        PriorityQueue { heap: Vec::new(), positions: vec![Self::ABSENT; items] }
    }

    fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }

    fn push(&mut self, score: Score, item: usize) {
        let newindex = self.heap.len();
        self.heap.push((score, item));
        self.positions[item] = newindex;
        self.sift_up(newindex);
    }

    fn swap(&mut self, i1: usize, i2: usize) {
        let (e1, e2) = (self.heap[i1], self.heap[i2]);
        self.positions.swap(e1.1, e2.1);
        self.heap[i1] = e2;
        self.heap[i2] = e1;
    }

    fn score_lower(&self, i1: usize, i2: usize) -> bool {
        self.heap[i1].0 < self.heap[i2].0
    }

    fn sift_up(&mut self, index: usize) {
        if index == 0 {
            return;
        }
        let parent = (index - 1) / 2;
        if self.score_lower(parent, index) {
            self.swap(parent, index);
            self.sift_up(parent);
        }
    }

    fn sift_down(&mut self, index: usize) {
        let r = 2 * index + 2;
        let l = 2 * index + 1;
        let n = self.heap.len();
        if r < n {
            if self.score_lower(l, r) {
                if self.score_lower(index, r) {
                    self.swap(r, index);
                    self.sift_down(r);
                }
            } else if self.score_lower(index, l) {
                self.swap(l, index);
                self.sift_down(l);
            }
        } else if l < n && self.score_lower(index, l) {
            self.swap(l, index);
            self.sift_down(l);
        }
    }

    fn pop(&mut self) -> (Score, usize) {
        assert!(!self.heap.is_empty(), "PriorityQueue empty.");
        let last = *self.heap.last().unwrap();
        let first = self.heap[0];
        if self.heap.len() == 1 {
            self.positions[first.1] = Self::ABSENT;
            self.heap.pop();
        } else {
            self.heap[0] = last;
            self.heap.pop();
            self.positions[last.1] = 0;
            self.positions[first.1] = Self::ABSENT;
            self.sift_down(0);
        }
        first
    }

    fn get_score_by_item(&self, item: usize) -> Option<Score> {
        let position = self.positions[item];
        if position == Self::ABSENT { None } else { Some(self.heap[position].0) }
    }

    fn change_score(&mut self, item: usize, new_score: Score) {
        let position = self.positions[item];
        let old = self.heap[position].0;
        self.heap[position].0 = new_score;
        if old < new_score {
            self.sift_up(position);
        } else {
            self.sift_down(position);
        }
    }
}

/// `coverage.CovMonitor`
struct CovMonitor {
    coverage: Vec<i32>,
}

impl CovMonitor {
    fn max_coverage_in_range(&self, begin: usize, end: usize) -> i32 {
        *self.coverage[begin..end].iter().max().expect("max() arg is an empty sequence")
    }

    fn add_read(&mut self, begin: usize, end: usize) {
        for c in &mut self.coverage[begin..end] {
            *c += 1;
        }
    }
}

struct Indexes {
    /// The variants' positions, ascending.
    positions: Vec<i32>,
    /// Each read's variants' indices in `positions`.
    read_variants: Vec<Vec<usize>>,
    variant_to_reads_map: Vec<Vec<usize>>,
}

impl Indexes {
    /// `vcf_indices[position]`: the index of a variant's position.
    fn of(&self, position: i32) -> usize {
        self.positions.binary_search(&position).expect("a variant of the read set")
    }

    /// The range of variant indices a read spans.
    fn span(&self, read: usize) -> (usize, usize) {
        let variants = &self.read_variants[read];
        (variants[0], variants[variants.len() - 1] + 1)
    }
}

fn compute_score_for_read(readset: &ReadSet, index: usize, idx: &Indexes) -> Score {
    let read = readset.get(index);
    let mut min_quality: i32 = -1;
    for (i, v) in read.variants.iter().enumerate() {
        let quality = v.quality as i32;
        if i == 0 {
            min_quality = quality;
        } else {
            min_quality = min_quality.min(quality);
        }
    }
    // Every variant of the read set has an index.
    let covered = &idx.read_variants[index];
    let good = covered.len() as i32;
    let span = (covered[covered.len() - 1] as i64 - covered[0] as i64 + 1) as i32;
    let bad = if good != span { span - good } else { 0 };
    [good - bad, good - bad, min_quality]
}

fn update_score_for_reads(former: Score, readset: &ReadSet, index: usize, covered: &CxxUnorderedSetI32) -> Score {
    let mut first = former[0];
    for v in &readset.get(index).variants {
        // Faithful to WhatsHap: decrements for variants *not* in the set.
        if !covered.contains(v.position) {
            first -= 1;
        }
    }
    [first, former[1], former[2]]
}

fn construct_priorityqueue(readset: &ReadSet, read_indices: &PySet, idx: &Indexes) -> PriorityQueue {
    let mut pq = PriorityQueue::new(readset.len());
    for index in read_indices.iter() {
        let score = compute_score_for_read(readset, index, idx);
        pq.push(score, index);
    }
    pq
}

fn slice_read_selection(
    pq: &mut PriorityQueue,
    coverages: &mut CovMonitor,
    max_cov: i32,
    readset: &ReadSet,
    idx: &Indexes,
) -> (PySet, PySet) {
    let mut already_covered = vec![false; idx.positions.len()];
    let mut reads_in_slice = PySet::new();
    let mut reads_violating_coverage = PySet::new();
    let mut variants_covered_by_this_read = CxxUnorderedSetI32::new();
    while !pq.is_empty() {
        variants_covered_by_this_read.clear();
        let (_max_score, max_item) = pq.pop();
        let read = readset.get(max_item);
        let mut covers_new_variant = false;
        for (v, &vi) in read.variants.iter().zip(&idx.read_variants[max_item]) {
            if already_covered[vi] {
                continue;
            }
            covers_new_variant = true;
            variants_covered_by_this_read.insert(v.position);
        }
        let (begin, end) = idx.span(max_item);
        if coverages.max_coverage_in_range(begin, end) >= max_cov {
            reads_violating_coverage.add(max_item);
        } else if covers_new_variant {
            coverages.add_read(begin, end);
            reads_in_slice.add(max_item);
            let mut to_update = PySet::new();
            for pos in variants_covered_by_this_read.iter() {
                let vi = idx.of(pos);
                already_covered[vi] = true;
                to_update.update_iter(idx.variant_to_reads_map[vi].iter().copied());
            }
            // WhatsHap subtracts a copy of `reads_in_slice`; a set difference reads the other
            // set only through its length and membership, so the set itself serves.
            let d_set = to_update.difference(&reads_in_slice);
            for element in d_set.iter() {
                if let Some(old) = pq.get_score_by_item(element) {
                    let new = update_score_for_reads(old, readset, element, &variants_covered_by_this_read);
                    pq.change_score(element, new);
                }
            }
        }
    }
    (reads_in_slice, reads_violating_coverage)
}

#[allow(clippy::too_many_arguments)]
fn readselection_helper(
    coverages: &mut CovMonitor,
    max_cov: i32,
    readset: &ReadSet,
    idx: &Indexes,
    selected_reads: &mut BTreeSet<usize>,
    undecided_reads: &mut PySet,
    bridging: bool,
) {
    let mut loop_count = 0;
    while !undecided_reads.is_empty() {
        let mut pq = construct_priorityqueue(readset, undecided_reads, idx);
        let (reads_in_slice, reads_violating_coverage) =
            slice_read_selection(&mut pq, coverages, max_cov, readset, idx);
        selected_reads.extend(reads_in_slice.iter());
        undecided_reads.difference_update(&reads_in_slice);
        undecided_reads.difference_update(&reads_violating_coverage);

        let mut component_finder = ComponentFinder::new(idx.positions.iter().copied());
        for read_index in reads_in_slice.iter() {
            let read = readset.get(read_index);
            for v in &read.variants[1..] {
                component_finder.merge(read.variants[0].position, v.position);
            }
        }

        let mut n_bridging = 0;
        if bridging {
            let mut pq = construct_priorityqueue(readset, undecided_reads, idx);
            while !pq.is_empty() {
                let (_score, read_index) = pq.pop();
                let read = readset.get(read_index);
                let mut covered_blocks = HashSet::new();
                for v in &read.variants {
                    covered_blocks.insert(component_finder.find(v.position));
                }
                let (begin, end) = idx.span(read_index);
                if coverages.max_coverage_in_range(begin, end) >= max_cov {
                    undecided_reads.remove(read_index);
                    continue;
                }
                if covered_blocks.len() < 2 {
                    continue;
                }
                n_bridging += 1;
                selected_reads.insert(read_index);
                coverages.add_read(begin, end);
                undecided_reads.remove(read_index);
                for v in &read.variants[1..] {
                    component_finder.merge(read.variants[0].position, v.position);
                }
            }
        }
        loop_count += 1;
        log::debug!(
            "... iteration {}: selected {} reads to cover positions and {} reads for bridging; {} reads left undecided",
            loop_count,
            reads_in_slice.len(),
            n_bridging,
            undecided_reads.len()
        );
    }
}

/// Return the indices of the selected reads (sorted).
pub fn readselection(
    readset: &ReadSet,
    max_cov: usize,
    preferred_source_ids: Option<&HashSet<i32>>,
    bridging: bool,
) -> Result<Vec<usize>> {
    let positions = readset.get_positions();
    let read_variants: Vec<Vec<usize>> = readset
        .iter()
        .map(|read| {
            read.variants.iter().map(|v| positions.binary_search(&v.position).expect("a read set position")).collect()
        })
        .collect();
    let mut variant_to_reads_map = vec![Vec::new(); positions.len()];
    let mut preferred_reads = PySet::new();
    for (index, read) in readset.iter().enumerate() {
        if let Some(pref) = preferred_source_ids
            && pref.contains(&read.source_id)
        {
            preferred_reads.add(index);
        }
        for &vi in &read_variants[index] {
            variant_to_reads_map[vi].push(index);
        }
    }
    let idx = Indexes { positions, read_variants, variant_to_reads_map };
    log::debug!(
        "Running read selection for {} reads covering {} variants (bridging {})",
        readset.len(),
        idx.positions.len(),
        if bridging { "ON" } else { "OFF" }
    );
    let mut coverages = CovMonitor { coverage: vec![0; idx.positions.len()] };
    let mut selected_reads = BTreeSet::new();
    for r in readset.iter() {
        if r.len() < 2 {
            bail!("readselection expects reads that cover at least two variants");
        }
    }
    let mut undecided_reads = PySet::from_iter_adds(0..readset.len());
    let max_cov = max_cov as i32;
    if !preferred_reads.is_empty() {
        // `preferred_reads` is consumed (emptied) by the helper, exactly as in
        // the Python code, so the subsequent difference_update removes nothing.
        readselection_helper(
            &mut coverages,
            max_cov,
            readset,
            &idx,
            &mut selected_reads,
            &mut preferred_reads,
            bridging,
        );
        undecided_reads.difference_update(&preferred_reads);
    }
    readselection_helper(&mut coverages, max_cov, readset, &idx, &mut selected_reads, &mut undecided_reads, bridging);
    Ok(selected_reads.into_iter().collect())
}

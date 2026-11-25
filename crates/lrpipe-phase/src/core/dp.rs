//! Port of the exact (PedMEC) dynamic programming algorithm:
//! `graycodes`, `columnindexingscheme`, `columnindexingiterator`,
//! `columniterator`, `pedigreecolumncostcomputer` and `pedigreedptable`.

use std::cmp::Reverse;
use std::collections::HashMap;
use std::ops::Range;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Result, bail};
use rayon::prelude::*;

use super::genotype::Genotype;
use super::pedigree::{Pedigree, PedigreePartitions};
use super::read::{ALT_ALLELE, BLANK, EQUAL_SCORES, Entry, REF_ALLELE, Read, ReadSet};

const UMAX: u32 = u32::MAX;

/// Gray code enumeration as in `graycodes.cpp` (Mossige 1977), which yields
/// the reflected binary Gray code: step k is `k ^ (k >> 1)` and changes bit
/// `trailing_zeros(k)`. (The test checks this against the original algorithm.)
struct GrayCodes {
    count: u64,
    k: u64,
}

impl GrayCodes {
    fn new(length: usize) -> Self {
        assert!(length <= 32);
        GrayCodes { count: 1 << length, k: 0 }
    }

    fn has_next(&self) -> bool {
        self.k < self.count
    }

    /// The next code and the bit that changed (-1 for the first code).
    fn get_next(&mut self, changed_bit: &mut i32) -> u32 {
        let k = self.k;
        self.k += 1;
        *changed_bit = if k == 0 { -1 } else { k.trailing_zeros() as i32 };
        (k ^ (k >> 1)) as u32
    }
}

pub(crate) struct ColumnIndexingScheme {
    read_ids: Vec<u32>,
    backward_projection_width: u32,
    /// For each read of this column, the bit index in the forward projection,
    /// or -1 if the read does not continue into the next column. `None` for the
    /// last column.
    forward_projection_mask: Option<Vec<i32>>,
    forward_projection_bits: u32,
}

impl ColumnIndexingScheme {
    pub(crate) fn new(previous: Option<&ColumnIndexingScheme>, read_ids: Vec<u32>) -> Self {
        let mut backward_projection_width = 0;
        if let Some(prev) = previous {
            let (mut i, mut j) = (0, 0);
            while i < prev.read_ids.len() && j < read_ids.len() {
                if prev.read_ids[i] == read_ids[j] {
                    backward_projection_width += 1;
                    i += 1;
                    j += 1;
                } else if prev.read_ids[i] < read_ids[j] {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
        ColumnIndexingScheme {
            read_ids,
            backward_projection_width,
            forward_projection_mask: None,
            forward_projection_bits: 0,
        }
    }

    pub(crate) fn set_next_column(&mut self, next_read_ids: &[u32]) {
        let mut mask = vec![-1i32; self.read_ids.len()];
        let (mut i, mut j, mut n) = (0, 0, 0);
        while i < next_read_ids.len() && j < self.read_ids.len() {
            if next_read_ids[i] == self.read_ids[j] {
                mask[j] = n;
                n += 1;
                i += 1;
                j += 1;
            } else if next_read_ids[i] < self.read_ids[j] {
                i += 1;
            } else {
                j += 1;
            }
        }
        self.forward_projection_mask = Some(mask);
        self.forward_projection_bits = n as u32;
    }

    /// Number of rows needed in the forward projection column. (C++ allocates
    /// 2^(#reads in column) rows, but only the first 2^(#shared reads) can be
    /// addressed, so this allocates just those.)
    pub(crate) fn forward_projection_size(&self) -> usize {
        1usize << self.forward_projection_bits
    }

    pub(crate) fn iter(&self) -> ColumnIndexingIterator<'_> {
        ColumnIndexingIterator {
            parent: self,
            graycodes: GrayCodes::new(self.read_ids.len()),
            index: UMAX,
            forward_projection: UMAX,
        }
    }

    fn index_backward_projection(&self, i: u32) -> u32 {
        i & ((1u64 << self.backward_projection_width) - 1) as u32
    }
}

pub(crate) struct ColumnIndexingIterator<'a> {
    parent: &'a ColumnIndexingScheme,
    graycodes: GrayCodes,
    index: u32,
    forward_projection: u32,
}

impl ColumnIndexingIterator<'_> {
    pub(crate) fn has_next(&self) -> bool {
        self.graycodes.has_next()
    }

    /// The current partition (bit set: read is in partition 1)
    pub(crate) fn partition(&self) -> u32 {
        self.index
    }

    pub(crate) fn forward_projection(&self) -> u32 {
        self.forward_projection
    }

    pub(crate) fn advance(&mut self) -> i32 {
        let mut bit = -1;
        self.index = self.graycodes.get_next(&mut bit);
        if bit == -1 {
            assert_eq!(self.index, 0);
            if self.parent.forward_projection_mask.is_some() {
                self.forward_projection = 0;
            }
        } else if let Some(mask) = &self.parent.forward_projection_mask {
            let bit_index = mask[bit as usize];
            if bit_index >= 0 {
                self.forward_projection ^= 1u32 << bit_index;
            }
        }
        bit
    }

    pub(crate) fn backward_projection(&self) -> u32 {
        self.parent.index_backward_projection(self.index)
    }
}

/// Iterates over the columns of the read/variant matrix (`columniterator.cpp`).
pub struct ColumnIterator<'a> {
    set: &'a ReadSet,
    /// Reads with indices from `read_end` on are not visited.
    read_end: usize,
    n: usize,
    next_read_index: usize,
    /// (read index, active entry)
    active_reads: Vec<(usize, usize)>,
    positions: Vec<i32>,
    first_reads: Vec<usize>,
}

impl<'a> ColumnIterator<'a> {
    pub fn new(set: &'a ReadSet, positions: Option<&[i32]>) -> Result<Self> {
        let positions: Vec<i32> = match positions {
            Some(p) => p.to_vec(),
            None => set.get_positions(),
        };
        Self::with_reads(set, 0..set.len(), positions)
    }

    /// Iterate over the given columns using only the reads in `reads`, which
    /// must all lie within those columns.
    fn with_reads(set: &'a ReadSet, reads: Range<usize>, positions: Vec<i32>) -> Result<Self> {
        let position_map: HashMap<i32, usize> = positions.iter().enumerate().map(|(i, &p)| (p, i)).collect();
        let mut first_reads = vec![usize::MAX; positions.len()];
        let mut pos = 0;
        for i in reads.clone() {
            let read = set.get(i);
            if read.first_position() < pos {
                bail!("ColumnIterator: reads in ReadSet are not sorted.");
            }
            if !read.is_sorted() {
                bail!("ColumnIterator: encountered read with unsorted variants.");
            }
            let first = position_map[&read.first_position()];
            let last = position_map[&read.last_position()];
            for fr in first_reads.iter_mut().take(last + 1).skip(first) {
                if *fr == usize::MAX {
                    *fr = i;
                }
            }
            pos = read.first_position();
        }
        if first_reads.len() >= 2 {
            let mut next_index = first_reads[first_reads.len() - 1];
            for i in (0..first_reads.len() - 1).rev() {
                if first_reads[i] == usize::MAX {
                    first_reads[i] = next_index;
                } else {
                    next_index = first_reads[i];
                }
            }
        }
        Ok(ColumnIterator {
            set,
            read_end: reads.end,
            n: 0,
            next_read_index: reads.start,
            active_reads: Vec::new(),
            positions,
            first_reads,
        })
    }

    pub fn column_count(&self) -> usize {
        self.positions.len()
    }

    pub fn positions(&self) -> &[i32] {
        &self.positions
    }

    pub fn has_next(&self) -> bool {
        self.n < self.positions.len()
    }

    pub fn get_next(&mut self) -> Vec<Entry> {
        let next_pos = self.positions[self.n];
        let set = self.set;
        self.active_reads.retain_mut(|(ri, ae)| {
            let read: &Read = set.get(*ri);
            if read.last_position() < next_pos {
                return false;
            }
            while read.variants[*ae].position < next_pos {
                *ae += 1;
            }
            true
        });
        while self.next_read_index < self.read_end {
            let start = set.get(self.next_read_index).first_position();
            if start == next_pos {
                self.active_reads.push((self.next_read_index, 0));
                self.next_read_index += 1;
            } else {
                assert!(start > next_pos);
                break;
            }
        }
        let mut result = Vec::with_capacity(self.active_reads.len());
        for &(ri, ae) in &self.active_reads {
            let read = set.get(ri);
            if read.variants[ae].position == next_pos {
                result.push(read.entry(ae));
            } else {
                result.push(Entry { read_id: read.id as u32, allele: BLANK, quality: 0 });
            }
        }
        self.n += 1;
        result
    }

    pub fn jump_to_column(&mut self, k: usize) {
        if k == self.n {
            return;
        }
        assert!(k < self.positions.len());
        self.active_reads.clear();
        self.n = k;
        self.next_read_index = self.first_reads[k];
        let pos = self.positions[k];
        while self.next_read_index < self.read_end {
            let read = self.set.get(self.next_read_index);
            if read.last_position() < pos {
                self.next_read_index += 1;
                continue;
            }
            if read.first_position() <= pos {
                let mut ae = 0;
                while read.variants[ae].position < pos {
                    ae += 1;
                }
                self.active_reads.push((self.next_read_index, ae));
                self.next_read_index += 1;
            } else {
                break;
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhasedVariant {
    pub allele0: i32,
    pub allele1: i32,
    pub quality: u32,
}

/// `PedigreeColumnCostComputer`
struct ColumnCostComputer<'a> {
    column: &'a [Entry],
    read_marks: &'a [u32],
    /// Member `partitioning` of the C++ class. Note that `set_partitioning`
    /// in C++ assigns its parameter to itself, so this member is only ever
    /// modified by `update_partitioning`. We keep that behaviour.
    partitioning: u32,
    pedigree: &'a Pedigree,
    cost_partition: Vec<[u32; 2]>,
    partitions: &'a PedigreePartitions,
    /// (assignment, cost)
    allele_assignments: Vec<(u32, u32)>,
}

impl<'a> ColumnCostComputer<'a> {
    fn new(
        column: &'a [Entry],
        column_index: usize,
        read_marks: &'a [u32],
        pedigree: &'a Pedigree,
        partitions: &'a PedigreePartitions,
        distrust_genotypes: bool,
    ) -> Self {
        let mut allele_assignments = Vec::new();
        for i in 0u32..(1u32 << partitions.count()) {
            let mut compatible = true;
            let mut cost: u32 = 0;
            for ind in 0..pedigree.size() {
                let p0 = partitions.haplotype_to_partition(ind, 0);
                let p1 = partitions.haplotype_to_partition(ind, 1);
                let a0 = (i >> p0) & 1;
                let a1 = (i >> p1) & 1;
                let genotype = Genotype::new(&[a0, a1]).unwrap();
                if distrust_genotypes {
                    let gls = pedigree.genotype_likelihoods(ind, column_index).expect("genotype likelihoods missing");
                    // `unsigned int += double` in C++
                    cost = (cost as f64 + gls.get(&genotype)) as u32;
                } else if genotype != *pedigree.genotype(ind, column_index) {
                    compatible = false;
                    break;
                }
            }
            if compatible {
                allele_assignments.push((i, cost));
            }
        }
        ColumnCostComputer {
            column,
            read_marks,
            partitioning: 0,
            pedigree,
            cost_partition: vec![[0, 0]; partitions.count()],
            partitions,
            allele_assignments,
        }
    }

    fn set_partitioning(&mut self, partitioning: u32) {
        for c in self.cost_partition.iter_mut() {
            *c = [0, 0];
        }
        let mut partitioning = partitioning;
        for entry in self.column {
            let in_p1 = (partitioning & 1) == 0;
            let ind = self.read_marks[entry.read_id as usize] as usize;
            let hap = if in_p1 { 0 } else { 1 };
            let p = self.partitions.haplotype_to_partition(ind, hap);
            match entry.allele {
                REF_ALLELE => self.cost_partition[p][1] = self.cost_partition[p][1].wrapping_add(entry.quality),
                ALT_ALLELE => self.cost_partition[p][0] = self.cost_partition[p][0].wrapping_add(entry.quality),
                BLANK => {}
                _ => panic!("unexpected allele in column"),
            }
            partitioning >>= 1;
        }
    }

    fn update_partitioning(&mut self, bit_to_flip: usize) {
        let entry = self.column[bit_to_flip];
        self.partitioning ^= 1u32 << bit_to_flip;
        let in_p1 = (self.partitioning & (1u32 << bit_to_flip)) == 0;
        let ind = self.read_marks[entry.read_id as usize] as usize;
        let (from, to) = if in_p1 { (1, 0) } else { (0, 1) };
        let pf = self.partitions.haplotype_to_partition(ind, from);
        let pt = self.partitions.haplotype_to_partition(ind, to);
        let slot = match entry.allele {
            REF_ALLELE => 1,
            ALT_ALLELE => 0,
            BLANK => return,
            _ => panic!("unexpected allele in column"),
        };
        self.cost_partition[pf][slot] = self.cost_partition[pf][slot].wrapping_sub(entry.quality);
        self.cost_partition[pt][slot] = self.cost_partition[pt][slot].wrapping_add(entry.quality);
    }

    fn assignment_cost(&self, assignment: u32, base: u32) -> u32 {
        let mut cost = base;
        for (p, cp) in self.cost_partition.iter().enumerate() {
            let allele = ((assignment >> p) & 1) as usize;
            cost = cost.wrapping_add(cp[allele]);
        }
        cost
    }

    fn get_cost(&self) -> u32 {
        let mut best = UMAX;
        for &(a, c) in &self.allele_assignments {
            let cost = self.assignment_cost(a, c);
            if cost < best {
                best = cost;
            }
        }
        best
    }

    fn get_alleles(&self) -> Result<Vec<PhasedVariant>> {
        let n = self.pedigree.size();
        let mut best_cost = UMAX;
        let mut pop_haps = vec![PhasedVariant { allele0: BLANK, allele1: BLANK, quality: 0 }; n];
        // best_cost_for_allele[individual][haplotype][to_allele]
        let mut best_for_allele = vec![[[UMAX; 2]; 2]; n];
        for &(a, c) in &self.allele_assignments {
            let cost = self.assignment_cost(a, c);
            let mut new_best = false;
            if cost <= best_cost {
                best_cost = cost;
                new_best = true;
            }
            for ind in 0..n {
                let p0 = self.partitions.haplotype_to_partition(ind, 0);
                let p1 = self.partitions.haplotype_to_partition(ind, 1);
                let a0 = ((a >> p0) & 1) as usize;
                let a1 = ((a >> p1) & 1) as usize;
                if new_best {
                    pop_haps[ind] = PhasedVariant {
                        allele0: if a0 == 0 { REF_ALLELE } else { ALT_ALLELE },
                        allele1: if a1 == 0 { REF_ALLELE } else { ALT_ALLELE },
                        quality: 0,
                    };
                }
                if cost < best_for_allele[ind][0][a0] {
                    best_for_allele[ind][0][a0] = cost;
                }
                if cost < best_for_allele[ind][1][a1] {
                    best_for_allele[ind][1][a1] = cost;
                }
            }
        }
        if best_cost == UMAX {
            bail!("Error: Mendelian conflict");
        }
        for (ind, best) in best_for_allele.iter().enumerate() {
            for (hap, b) in best.iter().enumerate() {
                let q = ((b[0] as i32).wrapping_sub(b[1] as i32)).wrapping_abs();
                pop_haps[ind].quality = q as u32;
                if q == 0 {
                    if hap == 0 {
                        pop_haps[ind].allele0 = EQUAL_SCORES;
                    } else {
                        pop_haps[ind].allele1 = EQUAL_SCORES;
                    }
                }
            }
        }
        Ok(pop_haps)
    }
}

#[cfg_attr(test, derive(Debug, PartialEq))]
struct Vec2D {
    size1: usize,
    v: Vec<u32>,
}

impl Vec2D {
    fn new(size0: usize, size1: usize, init: u32) -> Self {
        Vec2D { size1, v: vec![init; size0 * size1] }
    }
    #[inline]
    fn at(&self, i: usize, j: usize) -> u32 {
        self.v[i * self.size1 + j]
    }
    #[inline]
    fn set(&mut self, i: usize, j: usize, x: u32) {
        self.v[i * self.size1 + j] = x;
    }
}

#[derive(Clone, Copy, Default)]
struct IndexAndInheritance {
    index: u32,
    inheritance_value: u32,
}

/// Bytes of backtrace tables that may be held at once (over all threads) by
/// segments that keep every column instead of recomputing columns from
/// checkpoints during the backtrace. Either way gives the same result.
static TABLE_MEMORY_BUDGET: AtomicU64 = AtomicU64::new(1 << 30);

/// Change the memory budget for backtrace tables (see [`solve`]).
pub fn set_table_memory_budget(bytes: u64) {
    TABLE_MEMORY_BUDGET.store(bytes, Ordering::Relaxed);
}

/// Part of [`TABLE_MEMORY_BUDGET`], returned when dropped.
struct TableMemory(u64);

impl TableMemory {
    fn reserve(bytes: u64) -> Option<TableMemory> {
        TABLE_MEMORY_BUDGET
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |free| free.checked_sub(bytes))
            .ok()
            .map(|_| TableMemory(bytes))
    }
}

impl Drop for TableMemory {
    fn drop(&mut self) {
        TABLE_MEMORY_BUDGET.fetch_add(self.0, Ordering::Relaxed);
    }
}

/// Inputs shared by all segments of one (Ped)MEC instance.
struct Problem<'a> {
    read_set: &'a ReadSet,
    read_sources: Vec<u32>,
    recombcost: &'a [u32],
    pedigree: &'a Pedigree,
    distrust_genotypes: bool,
    pedigree_partitions: Vec<PedigreePartitions>,
    positions: Vec<i32>,
}

impl Problem<'_> {
    fn transmission_configurations(&self) -> usize {
        self.pedigree_partitions.len()
    }
}

/// Consecutive columns together with the reads starting in them, such that no
/// read crosses the segment boundaries.
struct Segment {
    columns: Range<usize>,
    reads: Range<usize>,
    /// Size of the backtrace tables if all columns are kept
    table_bytes: u64,
    /// Rough amount of work (number of DP rows)
    work: u64,
}

/// Split the matrix at every column boundary that no read crosses.
///
/// With a single transmission configuration (no trios), the DP of one segment
/// only sees the previous segments through a single projection value that is
/// added to every entry, and the backtrace enters the previous segment at the
/// first optimal row of its last column, which is also what solving that
/// segment alone picks. So solving the segments separately gives the same
/// result as solving everything at once (see [`offsets_fit`] for the overflow
/// caveat). With trios, recombination costs couple the segments.
fn find_segments(problem: &Problem, split: bool) -> Vec<Segment> {
    let read_set = problem.read_set;
    let column_count = problem.positions.len();
    let column_of: HashMap<i32, usize> = problem.positions.iter().enumerate().map(|(i, &p)| (p, i)).collect();
    // crossing[c]: reads active in both column c and c + 1; active[c]: reads
    // active in column c (both as difference arrays first)
    let mut crossing = vec![0i64; column_count + 1];
    let mut active = vec![0i64; column_count + 1];
    let mut first_columns = Vec::with_capacity(read_set.len());
    for read in read_set.iter() {
        let first = column_of[&read.first_position()];
        let last = column_of[&read.last_position()];
        crossing[first] += 1;
        crossing[last] -= 1;
        active[first] += 1;
        active[last + 1] -= 1;
        first_columns.push(first);
    }
    for c in 1..=column_count {
        crossing[c] += crossing[c - 1];
        active[c] += active[c - 1];
    }
    let row_bytes = 3 * size_of::<u32>() as u64 * problem.transmission_configurations() as u64;
    let mut segments = Vec::new();
    let (mut start, mut read_start, mut table_bytes, mut work) = (0, 0, 0u64, 0u64);
    for c in 0..column_count {
        work = work.saturating_add(1u64 << active[c].min(63));
        if c + 1 < column_count && !(split && crossing[c] == 0) {
            table_bytes = table_bytes.saturating_add(row_bytes.saturating_mul(1u64 << crossing[c].min(63)));
            continue;
        }
        let mut read_end = read_start;
        while read_end < first_columns.len() && first_columns[read_end] <= c {
            read_end += 1;
        }
        segments.push(Segment { columns: start..c + 1, reads: read_start..read_end, table_bytes, work });
        (start, read_start, table_bytes, work) = (c + 1, read_end, 0, 0);
    }
    segments
}

struct SegmentResult {
    index_path: Vec<IndexAndInheritance>,
    optimal_score: u32,
    /// Largest finite value in the DP table
    max_value: u32,
}

/// The exact PedMEC DP (`PedigreeDPTable`) for one segment.
struct SegmentTable<'p, 'a> {
    problem: &'p Problem<'a>,
    /// Index of the first column of this segment
    column_offset: usize,
    /// Keep all columns instead of recomputing them during the backtrace
    keep_all_columns: bool,
    indexers: Vec<ColumnIndexingScheme>,
    optimal_score: u32,
    optimal_score_index: u32,
    optimal_transmission_value: u32,
    previous_transmission_value: u32,
    max_value: u32,
    projection_column_table: Vec<Option<Vec2D>>,
    index_backtrace_table: Vec<Option<Vec2D>>,
    transmission_backtrace_table: Vec<Option<Vec2D>>,
    input_column_iterator: ColumnIterator<'a>,
    /// Buffers of [`SingleColumn::rows_binary`], kept from column to column
    scratch: SingleScratch,
}

impl<'p, 'a> SegmentTable<'p, 'a> {
    fn solve(problem: &'p Problem<'a>, segment: &Segment) -> Result<SegmentResult> {
        let reservation = TableMemory::reserve(segment.table_bytes);
        let positions = problem.positions[segment.columns.clone()].to_vec();
        let mut table = SegmentTable {
            problem,
            column_offset: segment.columns.start,
            keep_all_columns: reservation.is_some(),
            indexers: Vec::new(),
            optimal_score: 0,
            optimal_score_index: 0,
            optimal_transmission_value: 0,
            previous_transmission_value: 0,
            max_value: 0,
            projection_column_table: Vec::new(),
            index_backtrace_table: Vec::new(),
            transmission_backtrace_table: Vec::new(),
            input_column_iterator: ColumnIterator::with_reads(problem.read_set, segment.reads.clone(), positions)?,
            scratch: SingleScratch::default(),
        };
        let index_path = table.compute_table()?;
        Ok(SegmentResult { index_path, optimal_score: table.optimal_score, max_value: table.max_value })
    }

    fn compute_table(&mut self) -> Result<Vec<IndexAndInheritance>> {
        let column_count = self.input_column_iterator.column_count();
        self.projection_column_table = (0..column_count).map(|_| None).collect();
        self.index_backtrace_table = (0..column_count).map(|_| None).collect();
        self.transmission_backtrace_table = (0..column_count).map(|_| None).collect();
        self.indexers = Vec::with_capacity(column_count);
        self.optimal_score = UMAX;
        self.optimal_score_index = 0;
        self.optimal_transmission_value = 0;
        self.previous_transmission_value = 0;

        if column_count == 0 {
            self.optimal_score = 0;
            self.optimal_score_index = 0;
            return Ok(Vec::new());
        }

        // Build all indexing schemes up front (equivalent to the C++ code,
        // which builds them one column ahead of the computation).
        self.input_column_iterator.jump_to_column(0);
        let first = self.input_column_iterator.get_next();
        let first_ids: Vec<u32> = first.iter().map(|e| e.read_id).collect();
        self.indexers.push(ColumnIndexingScheme::new(None, first_ids));
        let mut next_column = Some(first);

        // Checkpoint distance; with k = 1 every column is kept.
        let k = if self.keep_all_columns { 1 } else { (column_count as f64).sqrt() as usize };
        for column_index in 0..column_count {
            let current = next_column.take().unwrap();
            if self.input_column_iterator.has_next() {
                let next = self.input_column_iterator.get_next();
                let next_ids: Vec<u32> = next.iter().map(|e| e.read_id).collect();
                let next_indexer = ColumnIndexingScheme::new(Some(&self.indexers[column_index]), next_ids.clone());
                self.indexers[column_index].set_next_column(&next_ids);
                self.indexers.push(next_indexer);
                next_column = Some(next);
            }
            self.compute_column(column_index, Some(current))?;
            if k > 1 && column_index > 0 && (column_index - 1) % k != 0 {
                self.index_backtrace_table[column_index - 1] = None;
                self.transmission_backtrace_table[column_index - 1] = None;
                self.projection_column_table[column_index - 1] = None;
            }
        }

        // Backtrace
        let n = self.indexers.len();
        let mut index_path = vec![IndexAndInheritance::default(); n];
        let mut prev_inheritance_value = self.previous_transmission_value;
        let mut v =
            IndexAndInheritance { index: self.optimal_score_index, inheritance_value: self.optimal_transmission_value };
        index_path[n - 1] = v;
        let mut i = n - 1;
        while i > 0 {
            if self.projection_column_table[i - 1].is_none() {
                let j0 = (i - 1) / k * k;
                assert!(self.projection_column_table[j0].is_some());
                for j in (j0 + 1)..i {
                    self.compute_column(j, None)?;
                }
            }
            let backtrace_index = self.indexers[i].index_backward_projection(v.index) as usize;
            v.index = self.index_backtrace_table[i - 1]
                .as_ref()
                .unwrap()
                .at(backtrace_index, prev_inheritance_value as usize);
            v.inheritance_value = prev_inheritance_value;
            prev_inheritance_value = self.transmission_backtrace_table[i - 1]
                .as_ref()
                .unwrap()
                .at(backtrace_index, v.inheritance_value as usize);
            index_path[i - 1] = v;
            if i.is_multiple_of(k) {
                let mut j = i;
                while j < i + k && j < column_count - 1 {
                    assert!(self.projection_column_table[j].is_some());
                    self.index_backtrace_table[j] = None;
                    self.transmission_backtrace_table[j] = None;
                    self.projection_column_table[j] = None;
                    j += 1;
                }
            }
            i -= 1;
        }
        Ok(index_path)
    }

    fn compute_column(&mut self, column_index: usize, current_input_column: Option<Vec<Entry>>) -> Result<()> {
        let column_count = self.input_column_iterator.column_count();
        if self.projection_column_table[column_index].is_some() {
            return Ok(());
        }
        let problem = self.problem;
        let configs = problem.transmission_configurations();
        let current_input_column = match current_input_column {
            Some(c) => c,
            None => {
                self.input_column_iterator.jump_to_column(column_index);
                self.input_column_iterator.get_next()
            }
        };
        if configs == 1 && problem.pedigree.size() == 1 {
            return self.compute_column_single(column_index, current_input_column);
        }
        let indexer = &self.indexers[column_index];
        // The C++ code stores a full DP column, but each row is only read back
        // right after it was computed, so one row suffices.
        let mut dp_row = vec![0u32; configs];
        let previous_projection_column =
            if column_index > 0 { self.projection_column_table[column_index - 1].as_ref() } else { None };

        let last_column = column_index + 1 >= column_count;
        let (mut current_projection_column, mut transmission_backtrace_column, mut index_backtrace_column) =
            if !last_column {
                let fps = indexer.forward_projection_size();
                (
                    Some(Vec2D::new(fps, configs, UMAX)),
                    Some(Vec2D::new(fps, configs, UMAX)),
                    Some(Vec2D::new(fps, configs, UMAX)),
                )
            } else {
                (None, None, None)
            };

        let global_column_index = self.column_offset + column_index;
        let mut cost_computers: Vec<ColumnCostComputer> = problem
            .pedigree_partitions
            .iter()
            .map(|partitions| {
                ColumnCostComputer::new(
                    &current_input_column,
                    global_column_index,
                    &problem.read_sources,
                    problem.pedigree,
                    partitions,
                    problem.distrust_genotypes,
                )
            })
            .collect();

        let recomb = problem.recombcost[global_column_index] as u64;
        let mut min_recomb_index = vec![0u32; configs];
        let mut max_value = self.max_value;
        let mut it = indexer.iter();
        while it.has_next() {
            let bit_changed = it.advance();
            if bit_changed >= 0 {
                for cc in cost_computers.iter_mut() {
                    cc.update_partitioning(bit_changed as usize);
                }
            } else {
                for cc in cost_computers.iter_mut() {
                    cc.set_partitioning(it.index);
                }
            }
            let backward_projection_index = if column_index > 0 { it.backward_projection() as usize } else { 0 };

            let mut found_valid = false;
            for i in 0..configs {
                let current_cost = cost_computers[i].get_cost();
                let mut min = UMAX;
                let mut min_index = 0u32;
                if current_cost < UMAX {
                    found_valid = true;
                }
                for j in 0..configs {
                    let previous_cost = match previous_projection_column {
                        Some(p) => p.at(backward_projection_index, j),
                        None => 0,
                    };
                    let mut val = if current_cost < UMAX && previous_cost < UMAX {
                        current_cost.wrapping_add(previous_cost)
                    } else {
                        UMAX
                    };
                    let x = (i ^ j) as u64;
                    if val < UMAX {
                        val = (val as u64).wrapping_add(x.count_ones() as u64 * recomb) as u32;
                    }
                    if val < min {
                        min = val;
                        min_index = j as u32;
                    }
                }
                dp_row[i] = min;
                min_recomb_index[i] = min_index;
                if min < UMAX && min > max_value {
                    max_value = min;
                }
            }
            if !found_valid {
                bail!("Error: Mendelian conflict");
            }

            match current_projection_column.as_mut() {
                None => {
                    for i in 0..configs {
                        if dp_row[i] < self.optimal_score {
                            self.optimal_score = dp_row[i];
                            self.optimal_score_index = it.index;
                            self.optimal_transmission_value = i as u32;
                            self.previous_transmission_value = min_recomb_index[i];
                        }
                    }
                }
                Some(proj) => {
                    let forward_index = it.forward_projection as usize;
                    let it_idx = it.index;
                    let tb = transmission_backtrace_column.as_mut().unwrap();
                    let ib = index_backtrace_column.as_mut().unwrap();
                    for i in 0..configs {
                        let d = dp_row[i];
                        if d < proj.at(forward_index, i) {
                            proj.set(forward_index, i, d);
                            ib.set(forward_index, i, it_idx);
                            tb.set(forward_index, i, min_recomb_index[i]);
                        }
                    }
                }
            }
        }
        drop(cost_computers);
        self.max_value = max_value;
        if let Some(proj) = current_projection_column {
            self.index_backtrace_table[column_index] = index_backtrace_column;
            self.transmission_backtrace_table[column_index] = transmission_backtrace_column;
            self.projection_column_table[column_index] = Some(proj);
        }
        Ok(())
    }
}

impl SegmentTable<'_, '_> {
    /// `compute_column` for a single individual without trios. There is one
    /// transmission configuration and the two partitions are the haplotypes,
    /// so the column cost of a partitioning is four sums ([`SingleColumn`]).
    /// The arithmetic (including u32 wrapping) and the outcome of the
    /// comparisons are those of the general code.
    fn compute_column_single(&mut self, column_index: usize, column: Vec<Entry>) -> Result<()> {
        let problem = self.problem;
        let column_count = self.input_column_iterator.column_count();
        let global_column_index = self.column_offset + column_index;
        let assignments = ColumnCostComputer::new(
            &column,
            global_column_index,
            &problem.read_sources,
            problem.pedigree,
            &problem.pedigree_partitions[0],
            problem.distrust_genotypes,
        )
        .allele_assignments;
        let entries: Vec<Option<(usize, u32)>> = column
            .iter()
            .map(|e| match e.allele {
                REF_ALLELE => Some((1, e.quality)),
                ALT_ALLELE => Some((0, e.quality)),
                BLANK => None,
                _ => panic!("unexpected allele in column"),
            })
            .collect();
        let single = SingleColumn {
            entries: &entries,
            assignments: &assignments,
            indexer: &self.indexers[column_index],
            previous: if column_index > 0 { self.projection_column_table[column_index - 1].as_ref() } else { None },
            last: column_index + 1 >= column_count,
        };
        let rows = if entries.len() <= SingleColumn::MAX_BINARY_READS {
            single.rows_binary(self.max_value, &mut self.scratch)?
        } else {
            single.rows_gray(self.max_value)?
        };
        self.max_value = rows.max_value;
        match rows.tables {
            None => {
                if let Some((value, index)) = rows.best
                    && value < self.optimal_score
                {
                    self.optimal_score = value;
                    self.optimal_score_index = index;
                    self.optimal_transmission_value = 0;
                    self.previous_transmission_value = 0;
                }
            }
            Some((projection, index_backtrace, transmission_backtrace)) => {
                self.index_backtrace_table[column_index] = Some(index_backtrace);
                self.transmission_backtrace_table[column_index] = Some(transmission_backtrace);
                self.projection_column_table[column_index] = Some(projection);
            }
        }
        Ok(())
    }
}

/// One column of the single-individual DP: the cost of a partitioning of the
/// column's reads (bit set: haplotype 1) is the cheapest allele assignment's,
/// with `cost[haplotype][slot]` the summed qualities of the reads on that
/// haplotype whose observation costs if the haplotype gets allele `slot`.
struct SingleColumn<'c> {
    /// Per read (bit): the slot its quality counts against, and the quality
    entries: &'c [Option<(usize, u32)>],
    /// (assignment, base cost); bit 0 of an assignment is haplotype 0's allele
    assignments: &'c [(u32, u32)],
    indexer: &'c ColumnIndexingScheme,
    previous: Option<&'c Vec2D>,
    last: bool,
}

/// Per partitioning of a column's reads, in binary order.
#[derive(Default)]
struct SingleScratch {
    /// Summed qualities of the reads on haplotype 1, per slot
    sums: [Vec<u32>; 2],
    /// Rank in Gray code order
    rank: Vec<u32>,
    forward_index: Vec<u32>,
    values: Vec<u32>,
}

impl SingleScratch {
    /// Sizes every buffer for `n` partitionings; their contents are left over.
    fn resize(&mut self, n: usize) {
        let [sums0, sums1] = &mut self.sums;
        for buffer in [sums0, sums1, &mut self.rank, &mut self.forward_index, &mut self.values] {
            if buffer.len() < n {
                buffer.resize(n, 0);
            }
        }
    }
}

/// What a [`SingleColumn`] yields.
#[cfg_attr(test, derive(Debug, PartialEq))]
struct SingleRows {
    /// Projection, index backtrace and transmission backtrace, unless last
    tables: Option<(Vec2D, Vec2D, Vec2D)>,
    /// In the last column, the first partitioning in Gray code order with the
    /// lowest finite value, and that value
    best: Option<(u32, u32)>,
    max_value: u32,
}

impl SingleColumn<'_> {
    /// Columns of up to this many reads are computed in binary order, which
    /// keeps a few words per partitioning.
    const MAX_BINARY_READS: usize = 20;

    fn tables(&self) -> Option<(Vec2D, Vec2D, Vec2D)> {
        (!self.last).then(|| {
            let fps = self.indexer.forward_projection_size();
            (Vec2D::new(fps, 1, UMAX), Vec2D::new(fps, 1, UMAX), Vec2D::new(fps, 1, UMAX))
        })
    }

    /// The partitionings in Gray code order, as `PedigreeDPTable` visits them:
    /// a row's projection entry is replaced only by a lower value, so it keeps
    /// the first lowest one in that order.
    fn rows_gray(&self, mut max_value: u32) -> Result<SingleRows> {
        let mut tables = self.tables();
        let mut best = None;
        let mut optimal = UMAX;
        let mut cost = [[0u32; 2]; 2];
        let mut it = self.indexer.iter();
        while it.has_next() {
            let bit_changed = it.advance();
            if bit_changed >= 0 {
                let bit = bit_changed as usize;
                if let Some((allele, quality)) = self.entries[bit] {
                    let (from, to) = if (it.index >> bit) & 1 == 0 { (1, 0) } else { (0, 1) };
                    cost[from][allele] = cost[from][allele].wrapping_sub(quality);
                    cost[to][allele] = cost[to][allele].wrapping_add(quality);
                }
            } else {
                cost = [[0; 2]; 2];
                for (k, entry) in self.entries.iter().enumerate() {
                    if let Some((allele, quality)) = *entry {
                        let haplotype = ((it.index >> k) & 1) as usize;
                        cost[haplotype][allele] = cost[haplotype][allele].wrapping_add(quality);
                    }
                }
            }
            let mut current_cost = UMAX;
            for &(a, base) in self.assignments {
                let c = base.wrapping_add(cost[0][(a & 1) as usize]).wrapping_add(cost[1][((a >> 1) & 1) as usize]);
                if c < current_cost {
                    current_cost = c;
                }
            }
            if current_cost == UMAX {
                bail!("Error: Mendelian conflict");
            }
            let previous_cost = match self.previous {
                Some(p) => p.at(it.backward_projection() as usize, 0),
                None => 0,
            };
            let value = if previous_cost < UMAX { current_cost.wrapping_add(previous_cost) } else { UMAX };
            if value < UMAX && value > max_value {
                max_value = value;
            }
            match tables.as_mut() {
                None => {
                    if value < optimal {
                        optimal = value;
                        best = Some((value, it.index));
                    }
                }
                Some((projection, index_backtrace, transmission_backtrace)) => {
                    let forward_index = it.forward_projection as usize;
                    if value < projection.at(forward_index, 0) {
                        projection.set(forward_index, 0, value);
                        index_backtrace.set(forward_index, 0, it.index);
                        transmission_backtrace.set(forward_index, 0, 0);
                    }
                }
            }
        }
        Ok(SingleRows { tables, best, max_value })
    }

    /// The same as [`Self::rows_gray`], computed with the partitionings in
    /// binary order: each partitioning's sums, Gray code rank and forward
    /// projection index come from those of the partitioning without its top
    /// bit, the costs in passes the compiler vectorizes, and ties are broken by
    /// the Gray code rank.
    fn rows_binary(&self, mut max_value: u32, scratch: &mut SingleScratch) -> Result<SingleRows> {
        let k = self.entries.len();
        let n = 1usize << k;
        let forward = if self.last { None } else { self.indexer.forward_projection_mask.as_deref() };
        scratch.resize(n);
        let SingleScratch { sums: [sums0, sums1], rank, forward_index, values } = scratch;
        let (sums0, sums1, rank) = (&mut sums0[..n], &mut sums1[..n], &mut rank[..n]);
        let (forward_index, values) = (&mut forward_index[..n], &mut values[..n]);
        (sums0[0], sums1[0], rank[0], forward_index[0]) = (0, 0, 0, 0);
        let mut totals = [0u32; 2];
        for (bit, entry) in self.entries.iter().enumerate() {
            let half = 1usize << bit;
            let (slot, quality) = entry.unwrap_or((0, 0));
            totals[slot] = totals[slot].wrapping_add(quality);
            for (s, sum) in [&mut *sums0, &mut *sums1].into_iter().enumerate() {
                let (low, high) = sum[..2 * half].split_at_mut(half);
                let add = if s == slot { quality } else { 0 };
                for (h, &l) in high.iter_mut().zip(low.iter()) {
                    *h = l.wrapping_add(add);
                }
            }
            // The Gray code at rank r is r ^ (r >> 1), so the rank of a
            // partitioning is the XOR of the ranks of its bits, 2^(bit+1) - 1.
            let (low, high) = rank[..2 * half].split_at_mut(half);
            let flip = (2 * half - 1) as u32;
            for (h, &l) in high.iter_mut().zip(low.iter()) {
                *h = l ^ flip;
            }
            if let Some(mask) = forward {
                let (low, high) = forward_index[..2 * half].split_at_mut(half);
                let set = if mask[bit] >= 0 { 1u32 << mask[bit] } else { 0 };
                for (h, &l) in high.iter_mut().zip(low.iter()) {
                    *h = l | set;
                }
            }
        }
        values.fill(UMAX);
        for &(a, base) in self.assignments {
            // haplotype 0 holds the reads not on haplotype 1
            let (slot0, slot1) = ((a & 1) as usize, ((a >> 1) & 1) as usize);
            let (on0, total0) = if slot0 == 0 { (&*sums0, totals[0]) } else { (&*sums1, totals[1]) };
            let on1 = if slot1 == 0 { &*sums0 } else { &*sums1 };
            for ((value, &other0), &cost1) in values.iter_mut().zip(on0.iter()).zip(on1.iter()) {
                let c = base.wrapping_add(total0.wrapping_sub(other0)).wrapping_add(cost1);
                *value = (*value).min(c);
            }
        }
        if values.contains(&UMAX) {
            bail!("Error: Mendelian conflict");
        }
        if let Some(previous) = self.previous {
            let backward = &previous.v[..1usize << self.indexer.backward_projection_width];
            // The reads shared with the previous column are the low bits.
            for rows in values.chunks_mut(backward.len()) {
                for (value, &previous_cost) in rows.iter_mut().zip(backward) {
                    *value = if previous_cost < UMAX { value.wrapping_add(previous_cost) } else { UMAX };
                }
            }
        }
        for &value in values.iter() {
            if value < UMAX && value > max_value {
                max_value = value;
            }
        }
        let mut tables = self.tables();
        let mut best = None;
        match tables.as_mut() {
            None => {
                let (mut lowest, mut lowest_rank) = (UMAX, UMAX);
                for (i, &value) in values.iter().enumerate() {
                    if value < lowest || (value == lowest && value < UMAX && rank[i] < lowest_rank) {
                        (lowest, lowest_rank) = (value, rank[i]);
                        best = Some((value, i as u32));
                    }
                }
            }
            Some((projection, index_backtrace, transmission_backtrace)) if projection.v.len() == n => {
                // Every read continues, so each partitioning has an entry of its own.
                debug_assert!(forward_index.iter().enumerate().all(|(i, &f)| f as usize == i));
                for (i, &value) in values.iter().enumerate() {
                    if value < UMAX {
                        projection.v[i] = value;
                        index_backtrace.v[i] = i as u32;
                        transmission_backtrace.v[i] = 0;
                    }
                }
            }
            Some((projection, index_backtrace, transmission_backtrace)) => {
                let mut lowest_rank = vec![UMAX; projection.v.len()];
                for (i, &value) in values.iter().enumerate() {
                    if value == UMAX {
                        continue;
                    }
                    let f = forward_index[i] as usize;
                    let held = projection.v[f];
                    if value < held || (value == held && rank[i] < lowest_rank[f]) {
                        projection.v[f] = value;
                        index_backtrace.v[f] = i as u32;
                        transmission_backtrace.v[f] = 0;
                        lowest_rank[f] = rank[i];
                    }
                }
            }
        }
        Ok(SingleRows { tables, best, max_value })
    }
}

/// Optimal solution of a (Ped)MEC instance.
pub struct Solution {
    pub optimal_score: u32,
    /// One pair of super reads per individual (in pedigree order)
    pub superreads: Vec<ReadSet>,
    pub transmission_vector: Vec<u32>,
    /// Haplotype (0 or 1) for each read, as returned by the Python wrapper's
    /// `get_optimal_partitioning`
    pub partitioning: Vec<u8>,
}

/// Solve the (Ped)MEC problem exactly, as WhatsHap's `PedigreeDPTable` does.
///
/// `read_set` must be sorted and have had `reassign_read_ids()` called (the
/// C++ constructor does that itself; in Rust the caller does it because we only
/// borrow the set).
///
/// Without trios, independent segments of the matrix are solved in parallel on
/// the current rayon thread pool.
pub fn solve(
    read_set: &ReadSet,
    recombcost: &[u32],
    pedigree: &Pedigree,
    distrust_genotypes: bool,
    positions: Option<&[i32]>,
) -> Result<Solution> {
    for (i, r) in read_set.iter().enumerate() {
        assert_eq!(r.id, i as i32, "call reassign_read_ids() before solving");
    }
    let n_configs = 4usize.pow(pedigree.triple_count() as u32);
    let pedigree_partitions = (0..n_configs).map(|i| PedigreePartitions::new(pedigree, i as u32)).collect();
    let read_sources =
        read_set.iter().map(|r| Ok(pedigree.id_to_index(r.sample_id as u32)? as u32)).collect::<Result<Vec<_>>>()?;
    let mut column_iterator = ColumnIterator::new(read_set, positions)?;
    let problem = Problem {
        read_set,
        read_sources,
        recombcost,
        pedigree,
        distrust_genotypes,
        pedigree_partitions,
        positions: column_iterator.positions().to_vec(),
    };
    let segments = find_segments(&problem, n_configs == 1);
    log::debug!(
        "Solving {} independent segment(s); largest has {:.1}% of the work and {} MB of backtrace tables",
        segments.len(),
        100.0 * segments.iter().map(|s| s.work).max().unwrap_or(0) as f64
            / segments.iter().map(|s| s.work).sum::<u64>().max(1) as f64,
        segments.iter().map(|s| s.table_bytes).max().unwrap_or(0) >> 20
    );
    let mut results = solve_segments(&problem, &segments)?;
    if segments.len() > 1 && !offsets_fit(&results) {
        log::debug!("Segment costs too large to be added up exactly, solving the MEC problem in one piece");
        results = solve_segments(&problem, &find_segments(&problem, false))?;
    }
    let optimal_score = results.iter().fold(0u32, |sum, r| sum.wrapping_add(r.optimal_score));
    let index_path: Vec<IndexAndInheritance> = results.into_iter().flat_map(|r| r.index_path).collect();

    // `get_super_reads` and `get_optimal_partitioning`
    let mut transmission_vector = Vec::new();
    let mut superreads: Vec<(Read, Read)> = (0..pedigree.size())
        .map(|i| {
            let id = pedigree.index_to_id(i) as i32;
            (Read::new(&format!("superread_0_{i}"), -1, -1, id), Read::new(&format!("superread_1_{i}"), -1, -1, id))
        })
        .collect();
    let mut in_first_partition = vec![false; read_set.len()];
    let mut i = 0;
    while column_iterator.has_next() {
        let v = index_path[i];
        let column = column_iterator.get_next();
        let mut cc = ColumnCostComputer::new(
            &column,
            i,
            &problem.read_sources,
            pedigree,
            &problem.pedigree_partitions[v.inheritance_value as usize],
            distrust_genotypes,
        );
        cc.set_partitioning(v.index);
        let alleles = cc.get_alleles()?;
        for (k, a) in alleles.iter().enumerate() {
            superreads[k].0.add_variant(problem.positions[i], a.allele0, a.quality as i32);
            superreads[k].1.add_variant(problem.positions[i], a.allele1, a.quality as i32);
        }
        transmission_vector.push(v.inheritance_value);
        let mut mask: u32 = 1;
        for entry in &column {
            if v.index & mask == 0 {
                in_first_partition[entry.read_id as usize] = true;
            }
            mask = mask.wrapping_shl(1);
        }
        i += 1;
    }
    let mut superread_sets = Vec::new();
    for (r0, r1) in superreads {
        let mut rs = ReadSet::new();
        rs.add(r0)?;
        rs.add(r1)?;
        superread_sets.push(rs);
    }
    let partitioning = in_first_partition.iter().map(|&x| if x { 0 } else { 1 }).collect();
    Ok(Solution { optimal_score, superreads: superread_sets, transmission_vector, partitioning })
}

/// Solve segments in parallel, largest first; the first error (in column
/// order) is returned.
fn solve_segments(problem: &Problem, segments: &[Segment]) -> Result<Vec<SegmentResult>> {
    let mut order: Vec<usize> = (0..segments.len()).collect();
    order.sort_by_key(|&i| Reverse(segments[i].work));
    let mut results: Vec<Option<Result<SegmentResult>>> = (0..segments.len()).map(|_| None).collect();
    let solved: Vec<(usize, Result<SegmentResult>)> =
        order.par_iter().with_max_len(1).map(|&i| (i, SegmentTable::solve(problem, &segments[i]))).collect();
    for (i, result) in solved {
        results[i] = Some(result);
    }
    results.into_iter().map(|r| r.expect("every segment was solved")).collect()
}

/// Whether the DP values of all segments, each increased by the optimal scores
/// of the segments before it (which is what a single DP over all columns would
/// compute), stay below `UMAX`. Only then are the separate solutions guaranteed
/// to be identical to the joint one.
fn offsets_fit(results: &[SegmentResult]) -> bool {
    let mut offset = 0u64;
    for r in results {
        if offset + r.max_value as u64 >= UMAX as u64 {
            return false;
        }
        offset += r.optimal_score as u64;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `graycodes.cpp`
    struct Mossige {
        length: i32,
        i: i32,
        s: u32,
        c: u32,
        changed_bit: i32,
    }

    impl Mossige {
        fn next(&mut self) -> Option<(u32, i32)> {
            if self.i >= self.length {
                return None;
            }
            let (result, changed) = (self.c, self.changed_bit);
            self.i = 0;
            while self.i < self.length {
                let mask = 1u32 << self.i;
                if ((self.c & mask) ^ (self.s & mask)) != 0 {
                    self.c ^= mask;
                    self.changed_bit = self.i;
                    break;
                }
                self.s ^= mask;
                self.i += 1;
            }
            Some((result, changed))
        }
    }

    #[test]
    fn gray_codes_match_graycodes_cpp() {
        for length in 0..=18usize {
            let mut reference = Mossige { length: length as i32, i: -1, s: !0, c: 0, changed_bit: -1 };
            let mut g = GrayCodes::new(length);
            while let Some((code, changed)) = reference.next() {
                assert!(g.has_next(), "length {length}");
                let mut bit = 0;
                assert_eq!(g.get_next(&mut bit), code, "length {length}");
                assert_eq!(bit, changed, "length {length}");
            }
            assert!(!g.has_next(), "length {length}");
        }
    }

    /// A small deterministic generator for the tests below.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: u64) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) % n
        }
    }

    #[test]
    fn single_column_rows_in_binary_order_match_gray_code_order() {
        let mut rng = Lcg(7);
        let mut scratch = SingleScratch::default();
        for case in 0..3000 {
            // Reads 0..k of this column; the previous column had reads from
            // `first` on and the next one keeps some of these.
            let k = rng.below(11) as u32;
            let first = rng.below(4) as u32;
            let previous_ids: Vec<u32> = (first..first + rng.below(8) as u32).collect();
            let read_ids: Vec<u32> = (2..2 + k).collect();
            let has_previous = rng.below(4) != 0;
            let mut indexer = ColumnIndexingScheme::new(
                has_previous.then_some(&ColumnIndexingScheme::new(None, previous_ids)),
                read_ids.clone(),
            );
            let last = rng.below(5) == 0;
            if !last {
                let next: Vec<u32> = read_ids.iter().copied().filter(|_| rng.below(3) != 0).chain(100..102).collect();
                indexer.set_next_column(&next);
            }
            // Low qualities and few reads make ties common.
            let entries: Vec<Option<(usize, u32)>> = (0..k)
                .map(|_| match rng.below(5) {
                    0 => None,
                    _ => Some((rng.below(2) as usize, rng.below(4) as u32)),
                })
                .collect();
            let mut assignments = Vec::new();
            for a in 0..4u32 {
                if rng.below(5) != 0 {
                    let base = if rng.below(10) == 0 { UMAX - rng.below(3) as u32 } else { rng.below(5) as u32 };
                    assignments.push((a, base));
                }
            }
            let previous = has_previous.then(|| {
                let size = 1usize << indexer.backward_projection_width;
                let mut table = Vec2D::new(size, 1, UMAX);
                for i in 0..size {
                    if rng.below(6) != 0 {
                        table.set(i, 0, rng.below(6) as u32);
                    }
                }
                table
            });
            let column = SingleColumn {
                entries: &entries,
                assignments: &assignments,
                indexer: &indexer,
                previous: previous.as_ref(),
                last,
            };
            let max_value = rng.below(4) as u32;
            match (column.rows_gray(max_value), column.rows_binary(max_value, &mut scratch)) {
                (Ok(gray), Ok(binary)) => assert_eq!(gray, binary, "case {case}"),
                (Err(_), Err(_)) => {}
                (gray, binary) => panic!("case {case}: {:?} against {:?}", gray.is_ok(), binary.is_ok()),
            }
        }
    }
}

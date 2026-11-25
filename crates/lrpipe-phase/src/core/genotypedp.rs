//! Port of WhatsHap's genotyping: `genotyper.cpp` (`compute_genotypes`, the
//! prior genotypes from all reads) and `genotypedptable.cpp` with
//! `transitionprobabilitycomputer.cpp` and `genotypecolumncostcomputer.cpp`
//! (forward-backward algorithm over the read/variant matrix).
//!
//! The DP computes in 80-bit `long double` (`compat::f80`). Its powers of
//! ten come from glibc's `powl`, which is not correctly rounded; the values
//! for phred scores and recombination costs below 1024 are taken from glibc
//! (`powl_table.rs`), and integer powers below 4 are computed like glibc
//! does. Other powers (costs of 1024 and more, or pedigrees with two or more
//! trios in one family) are approximations that can differ from WhatsHap in
//! the last bit.

use std::collections::HashMap;

use anyhow::{Result, bail};

use super::dp::{ColumnIndexingScheme, ColumnIterator};
use super::genotype::Genotype;
use super::pedigree::{Pedigree, PedigreePartitions};
use super::powl_table::POW10_NEG_TENTHS;
use super::read::{BLANK, Entry, REF_ALLELE, ReadSet};
use crate::compat::f80::{Divisor, F80};

/// glibc `powl(10, -c/10.0L)`
fn pow10_neg_tenths(c: u32) -> F80 {
    match POW10_NEG_TENTHS.get(c as usize) {
        Some(&(e, m)) => F80::from_bits(e, m),
        None => F80::from_f64(10f64.powf(-(c as f64) / 10.0)),
    }
}

/// glibc `powl(x, n)` for a small natural exponent: repeated squaring in
/// extended precision (exact for n < 4, as in glibc; an approximation of
/// glibc's general algorithm for larger n)
fn powl_int(x: F80, mut n: u32) -> F80 {
    if n == 0 {
        return F80::ONE;
    }
    let (mut r, mut b) = (F80::ONE, x);
    while n > 0 {
        if n & 1 == 1 {
            r *= b;
        }
        b *= b;
        n >>= 1;
    }
    r
}

/// `get_phred_probability`: `pow(10, -(int)phred_score/10.0L)`
fn phred_probability(phred_score: u32) -> F80 {
    if phred_score == 0 {
        return F80::from_f64(0.9999);
    }
    // scores of 2^31 and more are negative ints
    match phred_score as i32 {
        c if c >= 0 => pow10_neg_tenths(c as u32),
        c => F80::from_f64(10f64.powf(-(c as f64) / 10.0)),
    }
}

/// A phred probability p, 1 - p and both as divisors
#[derive(Clone, Copy)]
struct PhredProbability {
    p: F80,
    one_minus_p: F80,
    div_p: Divisor,
    div_one_minus_p: Divisor,
}

impl PhredProbability {
    fn new(phred_score: u32) -> Self {
        let p = phred_probability(phred_score);
        let one_minus_p = F80::ONE - p;
        PhredProbability { p, one_minus_p, div_p: Divisor::new(p), div_one_minus_p: Divisor::new(one_minus_p) }
    }
}

/// Phred probabilities of scores below 1024, computed on first use
fn phred_probabilities() -> &'static [PhredProbability] {
    static TABLE: std::sync::OnceLock<Vec<PhredProbability>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| (0..1024).map(PhredProbability::new).collect())
}

fn phred(phred_score: u32) -> PhredProbability {
    match phred_probabilities().get(phred_score as usize) {
        Some(p) => *p,
        None => PhredProbability::new(phred_score),
    }
}

/// `Vector2D<long double>`
#[derive(Clone)]
struct Table {
    width: usize,
    data: Vec<F80>,
}

impl Table {
    fn new(rows: usize, cols: usize) -> Self {
        Table { width: cols, data: vec![F80::ZERO; rows * cols] }
    }

    fn at(&self, i: usize, j: usize) -> F80 {
        self.data[i * self.width + j]
    }

    fn at_mut(&mut self, i: usize, j: usize) -> &mut F80 {
        &mut self.data[i * self.width + j]
    }

    fn divide_entries_by(&mut self, v: F80) {
        let d = Divisor::new(v);
        for x in &mut self.data {
            *x = x.div_by(&d);
        }
    }
}

/// `TransitionProbabilityComputer`
struct TransitionProbabilities {
    transmissions: Table,
    allele_assignments: Table,
    /// transmission(j, i) * allele_assignment(i, a) at [(i * A + a) * T + j],
    /// as the backward pass multiplies them
    products: Vec<F80>,
}

impl TransitionProbabilities {
    fn new(
        column_index: usize,
        recombcost: u32,
        pedigree: &Pedigree,
        partitions: &[PedigreePartitions],
    ) -> Result<Self> {
        let trio_count = pedigree.triple_count() as u32;
        let transmission_configurations = 4usize.pow(trio_count);
        let allele_assignments = 1usize << partitions[0].count();
        let mut transmissions = Table::new(transmission_configurations, transmission_configurations);
        let recomb_prob = pow10_neg_tenths(recombcost);
        let bernoulli: Vec<F80> = (0..=2 * trio_count)
            .map(|i| powl_int(recomb_prob, i) * powl_int(F80::ONE - recomb_prob, 2 * trio_count - i))
            .collect();
        for i in 0..transmission_configurations {
            let mut normalization_sum = F80::ZERO;
            for j in 0..transmission_configurations {
                let prob = bernoulli[(i ^ j).count_ones() as usize];
                *transmissions.at_mut(i, j) = prob;
                normalization_sum += prob;
            }
            for j in 0..transmission_configurations {
                *transmissions.at_mut(i, j) /= normalization_sum;
            }
        }
        let mut assignments = Table::new(transmission_configurations, allele_assignments);
        for (i, partition) in partitions.iter().enumerate().take(transmission_configurations) {
            let mut genotypes_to_haplotype_counts: HashMap<Vec<Genotype>, u64> = HashMap::new();
            let mut haplotypes_to_genotypes: Vec<Vec<Genotype>> = Vec::with_capacity(allele_assignments);
            for a in 0..allele_assignments {
                let mut prob = F80::ONE;
                let mut genotype_vector = Vec::with_capacity(pedigree.size());
                for individual in 0..pedigree.size() {
                    let allele0 = (a >> partition.haplotype_to_partition(individual, 0)) & 1;
                    let allele1 = (a >> partition.haplotype_to_partition(individual, 1)) & 1;
                    let genotype = Genotype::new(&[allele0 as u32, allele1 as u32])?;
                    let Some(gls) = pedigree.genotype_likelihoods(individual, column_index) else {
                        bail!("genotype likelihoods missing for the genotyping DP");
                    };
                    prob *= F80::from_f64(gls.get(&genotype));
                    genotype_vector.push(genotype);
                }
                *genotypes_to_haplotype_counts.entry(genotype_vector.clone()).or_default() += 1;
                *assignments.at_mut(i, a) = prob;
                haplotypes_to_genotypes.push(genotype_vector);
            }
            let mut normalization_sum = F80::ZERO;
            for a in 0..allele_assignments {
                *assignments.at_mut(i, a) /= F80::from_u64(genotypes_to_haplotype_counts[&haplotypes_to_genotypes[a]]);
                normalization_sum += assignments.at(i, a);
            }
            for a in 0..allele_assignments {
                *assignments.at_mut(i, a) /= normalization_sum;
            }
        }
        let mut products =
            Vec::with_capacity(transmission_configurations * allele_assignments * transmission_configurations);
        for i in 0..transmission_configurations {
            for a in 0..allele_assignments {
                for j in 0..transmission_configurations {
                    products.push(transmissions.at(j, i) * assignments.at(i, a));
                }
            }
        }
        Ok(TransitionProbabilities { transmissions, allele_assignments: assignments, products })
    }

    /// `get_prob_transmission(j, i) * get_prob_allele_assignment(i, a)` for all j
    fn products(&self, i: usize, a: usize) -> &[F80] {
        let t = self.transmissions.width;
        let start = (i * self.allele_assignments.width + a) * t;
        &self.products[start..start + t]
    }

    fn transmission(&self, t1: usize, t2: usize) -> F80 {
        self.transmissions.at(t1, t2)
    }

    fn allele_assignment(&self, t: usize, a: usize) -> F80 {
        self.allele_assignments.at(t, a)
    }
}

/// `GenotypeColumnCostComputer`
struct ColumnCosts<'a> {
    column: &'a [Entry],
    read_marks: &'a [usize],
    partitioning: u32,
    cost_partition: Vec<[F80; 2]>,
    partitions: &'a PedigreePartitions,
}

impl<'a> ColumnCosts<'a> {
    fn new(column: &'a [Entry], read_marks: &'a [usize], partitions: &'a PedigreePartitions) -> Self {
        ColumnCosts {
            column,
            read_marks,
            partitioning: 0,
            cost_partition: vec![[F80::ONE, F80::ONE]; partitions.count()],
            partitions,
        }
    }

    fn set_partitioning(&mut self, mut p: u32) {
        self.cost_partition = vec![[F80::ONE, F80::ONE]; self.partitions.count()];
        self.partitioning = p;
        for entry in self.column {
            // (blank entries do not consume a bit, as in WhatsHap)
            if entry.allele == BLANK {
                continue;
            }
            let in_partition1 = p & 1 == 0;
            let individual = self.read_marks[entry.read_id as usize];
            let is_ref = entry.allele == REF_ALLELE;
            let proba = phred(entry.quality);
            let c =
                &mut self.cost_partition[self.partitions.haplotype_to_partition(individual, in_partition1 as usize)];
            c[usize::from(!is_ref)] *= proba.one_minus_p;
            c[usize::from(is_ref)] *= proba.p;
            p >>= 1;
        }
    }

    fn update_partitioning(&mut self, bit_to_flip: i32) {
        let entry = &self.column[bit_to_flip as usize];
        if entry.allele == BLANK {
            return;
        }
        self.partitioning ^= 1 << bit_to_flip;
        let in_partition1 = self.partitioning & (1 << bit_to_flip) == 0;
        let individual = self.read_marks[entry.read_id as usize];
        let is_ref = entry.allele == REF_ALLELE;
        let proba = phred(entry.quality);
        let new = self.partitions.haplotype_to_partition(individual, in_partition1 as usize);
        let old = self.partitions.haplotype_to_partition(individual, (!in_partition1) as usize);
        self.cost_partition[new][usize::from(!is_ref)] *= proba.one_minus_p;
        self.cost_partition[new][usize::from(is_ref)] *= proba.p;
        let c = &mut self.cost_partition[old];
        c[usize::from(!is_ref)] = c[usize::from(!is_ref)].div_by(&proba.div_one_minus_p);
        c[usize::from(is_ref)] = c[usize::from(is_ref)].div_by(&proba.div_p);
    }

    fn cost(&self, allele_assignment: usize) -> F80 {
        // (WhatsHap starts from 1.0, and multiplying by one is exact)
        let mut parts = self.cost_partition.iter().enumerate();
        let Some((_, first)) = parts.next() else { return F80::ONE };
        let mut cost = first[allele_assignment & 1];
        for (p, c) in parts {
            cost *= c[(allele_assignment >> p) & 1];
        }
        cost
    }
}

/// `GenotypeDPTable`: genotype likelihoods per individual and position
pub struct GenotypeDPTable {
    /// likelihoods[individual index][column] for 0/0, 0/1, 1/1
    likelihoods: Vec<Vec<[F80; 3]>>,
    individual_index: HashMap<u32, usize>,
}

struct Dp<'a> {
    columns: Vec<Vec<Entry>>,
    read_sources: Vec<usize>,
    pedigree: &'a Pedigree,
    partitions: Vec<PedigreePartitions>,
    indexers: Vec<ColumnIndexingScheme>,
    transitions: Vec<TransitionProbabilities>,
    forward_projection: Option<Table>,
    backward_projection: Vec<Option<Table>>,
    scaling: Vec<F80>,
    likelihoods: Vec<Vec<[F80; 3]>>,
}

impl GenotypeDPTable {
    /// `read_set` must be sorted; `positions` are the columns.
    pub fn new(read_set: &mut ReadSet, recombcost: &[u32], pedigree: &Pedigree, positions: &[i32]) -> Result<Self> {
        read_set.reassign_read_ids();
        let mut iterator = ColumnIterator::new(read_set, Some(positions))?;
        let mut columns = Vec::with_capacity(iterator.column_count());
        while iterator.has_next() {
            columns.push(iterator.get_next());
        }
        let column_count = columns.len();
        let trio_count = pedigree.triple_count() as u32;
        let partitions: Vec<PedigreePartitions> =
            (0..4u32.pow(trio_count)).map(|i| PedigreePartitions::new(pedigree, i)).collect();
        let read_sources = (0..read_set.len())
            .map(|i| pedigree.id_to_index(read_set.get(i).sample_id as u32))
            .collect::<Result<Vec<_>>>()?;
        let mut dp = Dp {
            columns,
            read_sources,
            pedigree,
            partitions,
            indexers: Vec::with_capacity(column_count),
            transitions: Vec::with_capacity(column_count),
            forward_projection: None,
            backward_projection: Vec::new(),
            scaling: vec![F80::from_f64(-1.0); column_count],
            likelihoods: vec![vec![[F80::ZERO; 3]; column_count]; pedigree.size()],
        };
        dp.compute_index(recombcost)?;
        dp.compute_backward_prob();
        dp.compute_forward_prob();
        let individual_index = (0..pedigree.size()).map(|i| (pedigree.index_to_id(i), i)).collect();
        Ok(GenotypeDPTable { likelihoods: dp.likelihoods, individual_index })
    }

    /// `get_genotype_likelihoods` (converted to double)
    pub fn genotype_likelihoods(&self, individual_id: u32, position: usize) -> [f64; 3] {
        self.likelihoods[self.individual_index[&individual_id]][position].map(|l| l.to_f64())
    }
}

impl Dp<'_> {
    fn transmission_configurations(&self) -> usize {
        self.partitions.len()
    }

    fn read_ids(column: &[Entry]) -> Vec<u32> {
        column.iter().map(|e| e.read_id).collect()
    }

    fn compute_index(&mut self, recombcost: &[u32]) -> Result<()> {
        let column_count = self.columns.len();
        if column_count == 0 {
            return Ok(());
        }
        self.indexers.push(ColumnIndexingScheme::new(None, Self::read_ids(&self.columns[0])));
        self.transitions.push(TransitionProbabilities::new(0, recombcost[0], self.pedigree, &self.partitions)?);
        for column_index in 0..column_count {
            if column_index + 1 < column_count {
                let next_read_ids = Self::read_ids(&self.columns[column_index + 1]);
                let next = ColumnIndexingScheme::new(Some(&self.indexers[column_index]), next_read_ids.clone());
                self.indexers[column_index].set_next_column(&next_read_ids);
                self.indexers.push(next);
                self.transitions.push(TransitionProbabilities::new(
                    column_index + 1,
                    recombcost[column_index + 1],
                    self.pedigree,
                    &self.partitions,
                )?);
            }
        }
        Ok(())
    }

    /// floor(sqrt(column_count)): every k-th backward column is kept
    fn checkpoint_distance(&self) -> usize {
        (self.columns.len() as f64).sqrt() as usize
    }

    fn compute_backward_prob(&mut self) {
        let column_count = self.columns.len();
        self.backward_projection = vec![None; column_count];
        if column_count == 0 {
            return;
        }
        let k = self.checkpoint_distance();
        for column_index in (0..column_count).rev() {
            self.compute_backward_column(column_index);
            if k > 1 && column_index < column_count - 1 && (column_index + 1) % k != 0 {
                self.backward_projection[column_index + 1] = None;
            }
        }
    }

    fn compute_forward_prob(&mut self) {
        self.forward_projection = None;
        for column_index in 0..self.columns.len() {
            self.compute_forward_column(column_index);
        }
    }

    fn compute_backward_column(&mut self, column_index: usize) {
        let column_count = self.columns.len();
        if column_index > 0 && self.backward_projection[column_index - 1].is_some() {
            return;
        }
        let t_count = self.transmission_configurations();
        let mut previous =
            if column_index < column_count - 1 { self.backward_projection[column_index].take() } else { None };
        let mut current = if column_index > 0 {
            Some(Table::new(self.indexers[column_index - 1].forward_projection_size(), t_count))
        } else {
            None
        };
        let column = &self.columns[column_index];
        let mut costs: Vec<ColumnCosts> =
            self.partitions.iter().map(|p| ColumnCosts::new(column, &self.read_sources, p)).collect();
        let transitions = &self.transitions[column_index];
        let mut scaling_sum = F80::ZERO;
        let mut iterator = self.indexers[column_index].iter();
        while iterator.has_next() {
            let bit_changed = iterator.advance();
            if bit_changed >= 0 {
                for c in &mut costs {
                    c.update_partitioning(bit_changed);
                }
            } else {
                for c in &mut costs {
                    c.set_partitioning(iterator.partition());
                }
            }
            let mut backward_prob = F80::ONE;
            for (i, partition) in self.partitions.iter().enumerate() {
                let number_of_allele_assignments = 1usize << partition.count();
                if column_index + 1 < column_count {
                    let forward_projection_index = iterator.forward_projection() as usize;
                    backward_prob = previous.as_ref().expect("backward projection").at(forward_projection_index, i);
                }
                let backward_projection_index = iterator.backward_projection() as usize;
                for a in 0..number_of_allele_assignments {
                    if let Some(current) = current.as_mut() {
                        // (backward_prob * local_cost) * transition_prob, as in WhatsHap
                        let weighted = backward_prob * costs[i].cost(a);
                        for (j, &transition_prob) in transitions.products(i, a).iter().enumerate() {
                            *current.at_mut(backward_projection_index, j) += weighted * transition_prob;
                        }
                    }
                    scaling_sum += backward_prob;
                }
            }
        }
        if let Some(p) = previous.as_mut() {
            p.divide_entries_by(scaling_sum);
        }
        if column_index < column_count - 1 {
            self.backward_projection[column_index] = previous;
        }
        if let Some(mut c) = current {
            c.divide_entries_by(scaling_sum);
            self.backward_projection[column_index - 1] = Some(c);
        }
        self.scaling[column_index] = scaling_sum;
    }

    fn compute_forward_column(&mut self, column_index: usize) {
        let column_count = self.columns.len();
        let t_count = self.transmission_configurations();
        let k = self.checkpoint_distance();
        if column_index + 1 < column_count && self.backward_projection[column_index].is_none() {
            let next = ((column_index + k) / k * k).min(column_count - 1);
            for i in ((column_index + 1)..=next).rev() {
                self.compute_backward_column(i);
            }
            let scaling = self.scaling[column_index];
            self.backward_projection[column_index]
                .as_mut()
                .expect("recomputed backward column")
                .divide_entries_by(scaling);
        }
        let backward =
            if column_index + 1 < column_count { self.backward_projection[column_index].take() } else { None };
        let previous = if column_index > 0 { self.forward_projection.take() } else { None };
        let mut current = if column_index + 1 < column_count {
            Some(Table::new(self.indexers[column_index].forward_projection_size(), t_count))
        } else {
            None
        };
        let column = &self.columns[column_index];
        let mut costs: Vec<ColumnCosts> =
            self.partitions.iter().map(|p| ColumnCosts::new(column, &self.read_sources, p)).collect();
        let transitions = &self.transitions[column_index];
        let scaling = Divisor::new(self.scaling[column_index]);
        let mut normalization = F80::ZERO;
        let mut iterator = self.indexers[column_index].iter();
        while iterator.has_next() {
            let bit_changed = iterator.advance();
            if bit_changed >= 0 {
                for c in &mut costs {
                    c.update_partitioning(bit_changed);
                }
            } else {
                for c in &mut costs {
                    c.set_partitioning(iterator.partition());
                }
            }
            let backward_projection_index = if column_index > 0 { iterator.backward_projection() as usize } else { 0 };
            for (i, partition) in self.partitions.iter().enumerate() {
                let mut sum_prev_values = F80::ZERO;
                let number_of_allele_assignments = 1usize << partition.count();
                if let Some(previous) = &previous {
                    for j in 0..t_count {
                        sum_prev_values += previous.at(backward_projection_index, j) * transitions.transmission(j, i);
                    }
                } else {
                    sum_prev_values = F80::ONE;
                }
                for a in 0..number_of_allele_assignments {
                    let mut backward_probability = F80::ONE;
                    if let Some(backward) = &backward {
                        backward_probability = backward.at(iterator.forward_projection() as usize, i);
                    }
                    let forward_probability =
                        (sum_prev_values * costs[i].cost(a) * transitions.allele_assignment(i, a)).div_by(&scaling);
                    let forward_backward = forward_probability * backward_probability;
                    normalization += forward_backward;
                    for (individual, likelihoods) in self.likelihoods.iter_mut().enumerate() {
                        let allele0 = (a >> partition.haplotype_to_partition(individual, 0)) & 1;
                        let allele1 = (a >> partition.haplotype_to_partition(individual, 1)) & 1;
                        likelihoods[column_index][allele0 + allele1] += forward_backward;
                    }
                    if let Some(current) = current.as_mut() {
                        *current.at_mut(iterator.forward_projection() as usize, i) += forward_probability;
                    }
                }
            }
        }
        if current.is_some() {
            self.forward_projection = current;
        }
        // the backward column of this column is no longer needed
        let normalization = Divisor::new(normalization);
        for likelihoods in &mut self.likelihoods {
            for l in &mut likelihoods[column_index] {
                *l = l.div_by(&normalization);
            }
        }
    }
}

/// `GenotypeDistribution`
#[derive(Clone, Copy)]
struct GenotypeDistribution([f64; 3]);

impl GenotypeDistribution {
    fn uniform() -> Self {
        GenotypeDistribution([1.0 / 3.0; 3])
    }

    fn mul(self, other: GenotypeDistribution) -> Self {
        let mut d = self.0;
        let mut sum = 0.0;
        for (x, o) in d.iter_mut().zip(other.0) {
            *x *= o;
            sum += *x;
        }
        for x in &mut d {
            *x /= sum;
        }
        GenotypeDistribution(d)
    }

    fn normalize(&mut self) {
        let p_sum: f64 = self.0[0] + self.0[1] + self.0[2];
        if p_sum <= 0.0 {
            *self = GenotypeDistribution::uniform();
        } else {
            for x in &mut self.0 {
                *x /= p_sum;
            }
        }
    }

    fn best(&self) -> usize {
        let (mut best_index, mut best) = (0, 0.0);
        for (i, &p) in self.0.iter().enumerate() {
            if p > best {
                best = p;
                best_index = i;
            }
        }
        best_index
    }

    fn error_probability(&self) -> f64 {
        let best_index = self.best();
        let mut p = 0.0;
        for (i, &x) in self.0.iter().enumerate() {
            if i != best_index {
                p += x;
            }
        }
        p
    }
}

/// `compute_genotypes`: genotype (none if uncertain) and the probabilities of
/// 0/0, 0/1 and 1/1 per position, from the alleles the reads show
pub fn compute_genotypes(readset: &ReadSet, positions: &[i32]) -> Result<(Vec<Genotype>, Vec<[f64; 3]>)> {
    let mut iterator = ColumnIterator::new(readset, Some(positions))?;
    let mut genotypes = Vec::with_capacity(positions.len());
    let mut likelihoods = Vec::with_capacity(positions.len());
    while iterator.has_next() {
        let column = iterator.get_next();
        let mut distribution = GenotypeDistribution::uniform();
        for e in &column {
            let p_wrong = f64::max(0.05, 10f64.powf(-(e.quality as f64) / 10.0));
            match e.allele {
                REF_ALLELE => {
                    distribution = distribution.mul(GenotypeDistribution([
                        2.0 / 3.0 - 1.0 / 3.0 * p_wrong,
                        1.0 / 3.0,
                        1.0 / 3.0 * p_wrong,
                    ]))
                }
                super::read::ALT_ALLELE => {
                    distribution = distribution.mul(GenotypeDistribution([
                        1.0 / 3.0 * p_wrong,
                        1.0 / 3.0,
                        2.0 / 3.0 - 1.0 / 3.0 * p_wrong,
                    ]))
                }
                _ => {}
            }
        }
        distribution.normalize();
        let genotype = if distribution.error_probability() < 0.1 {
            Genotype::from_index(distribution.best() as u64, 2)?
        } else {
            Genotype::none()
        };
        genotypes.push(genotype);
        likelihoods.push(distribution.0);
    }
    Ok((genotypes, likelihoods))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glibc_integer_powers() {
        let x = F80::from_f64(0.3);
        assert_eq!(powl_int(x, 0), F80::ONE);
        assert_eq!(powl_int(x, 1), x);
        assert_eq!(powl_int(x, 2), x * x);
        assert_eq!(powl_int(x, 3), x * (x * x));
    }
}

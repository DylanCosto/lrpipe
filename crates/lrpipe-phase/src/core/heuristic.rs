//! Port of `pedmecheuristic.{h,cpp}` (`phase --algorithm heuristic`).
//!
//! Scores are `float` (`MecScore`) in C++, so all score arithmetic is done in
//! f32 in the same order; comparisons use the C++ semantics of `std::min`,
//! `std::max` and `std::max_element`.

use std::collections::{BTreeSet, HashMap};

use anyhow::{Result, bail};

use super::dp::Solution;
use super::pedigree::Pedigree;
use super::read::{Read, ReadSet};

/// `MAX_ROW_LIMIT`
pub const MAX_ROW_LIMIT: u32 = 65535;

type Balance = Vec<f32>;

#[derive(Clone)]
struct PedSolution {
    bp: Vec<bool>,
    bp_new: Vec<bool>,
    trans: u32,
    score: f32,
    mutation_score: f32,
    bt_row: u16,
    balances: Vec<Balance>,
}

impl PedSolution {
    fn new(bp: Vec<bool>, trans: u32, score: f32, num_samples: usize) -> Self {
        PedSolution {
            bp,
            bp_new: Vec::new(),
            trans,
            score,
            mutation_score: 0.0,
            bt_row: 0,
            balances: vec![vec![0.0]; num_samples * 2],
        }
    }

    /// The constructor used by `extendSolutions`, which appends `bp_new` to
    /// `bp` (and keeps `bp_new`, so `finalize` appends it a second time).
    fn extended(
        bp: &[bool],
        bp_new: &[bool],
        trans: u32,
        score: f32,
        mutation_score: f32,
        bt_row: u16,
        balances: &[Balance],
    ) -> Self {
        let mut full = bp.to_vec();
        full.extend_from_slice(bp_new);
        PedSolution {
            bp: full,
            bp_new: bp_new.to_vec(),
            trans,
            score,
            mutation_score,
            bt_row,
            balances: balances.to_vec(),
        }
    }

    fn finalize(&mut self) {
        let new = self.bp_new.clone();
        self.bp.extend(new);
    }
}

/// `std::min(a, b)`
fn min(a: f32, b: f32) -> f32 {
    if b < a { b } else { a }
}

/// `std::max(a, b)`
fn max(a: f32, b: f32) -> f32 {
    if a < b { b } else { a }
}

fn b2f(b: bool) -> f32 {
    b as u8 as f32
}

struct Heuristic<'a> {
    rs: &'a ReadSet,
    recomb_cost: Vec<f32>,
    distrust_genotypes: bool,
    positions: Vec<u32>,
    row_limit: u32,
    tm_bits: u32,
    num_samples: usize,
    global_sample_ids: Vec<u32>,
    sample_map: HashMap<u32, usize>,
    trios: Vec<[usize; 3]>,
    /// per sample per position: 0 = 0/0, 1 = 0/1, 2 = 1/1
    genotypes: Vec<Vec<i32>>,
    pos_map: HashMap<u32, usize>,
    mutation_cost: Vec<f32>,
}

impl<'a> Heuristic<'a> {
    fn new(
        rs: &'a ReadSet,
        recombcost: &[u32],
        pedigree: &Pedigree,
        distrust_genotypes: bool,
        positions: &[i32],
        row_limit: u32,
    ) -> Result<Self> {
        let n = recombcost.len();
        if n == 0 {
            bail!("The heuristic needs at least one variant to phase");
        }
        // allowMutations is always true in `phase`
        let mut recomb_cost = vec![0f32; n];
        let mut mutation_cost = vec![f32::INFINITY; n];
        for i in 1..n {
            recomb_cost[i] = recombcost[i] as f32;
            mutation_cost[i - 1] = (0.75 * (recomb_cost[i - 1] + recomb_cost[i]) as f64) as f32;
        }
        mutation_cost[n - 1] = (recomb_cost[n - 1] as f64 * 1.5) as f32;

        let positions: Vec<u32> = positions.iter().map(|&p| p as u32).collect();
        let mut pos_map = HashMap::new();
        for (i, &p) in positions.iter().enumerate() {
            pos_map.insert(p, i);
        }
        // Sample ids of the reads and the *indices* of the trio members, as
        // in C++ (they coincide for families whose ids start at 0).
        let mut sample_set: BTreeSet<u32> = rs.iter().map(|r| r.sample_id as u32).collect();
        for trio in pedigree.triples() {
            sample_set.extend(trio.iter().map(|&s| s as u32));
        }
        let global_sample_ids: Vec<u32> = sample_set.into_iter().collect();
        let num_samples = global_sample_ids.len();
        let sample_map: HashMap<u32, usize> = global_sample_ids.iter().enumerate().map(|(i, &s)| (s, i)).collect();
        let trios: Vec<[usize; 3]> = pedigree
            .triples()
            .iter()
            .map(|t| [sample_map[&(t[0] as u32)], sample_map[&(t[1] as u32)], sample_map[&(t[2] as u32)]])
            .collect();

        let n_pos = positions.len();
        let mut genotypes = Vec::with_capacity(num_samples);
        for s in 0..num_samples {
            if s >= pedigree.size() {
                bail!("The heuristic needs consecutive numeric sample ids starting at 0");
            }
            let mut g = Vec::with_capacity(n_pos);
            for p in 0..n_pos {
                let a = pedigree.genotype(s, p).as_vector();
                if a.len() < 2 {
                    bail!("The heuristic needs a diploid genotype for every sample at every variant");
                }
                g.push((a[0] + a[1]) as i32);
            }
            genotypes.push(g);
        }
        Ok(Heuristic {
            rs,
            recomb_cost,
            distrust_genotypes,
            positions,
            row_limit: row_limit.min(MAX_ROW_LIMIT),
            tm_bits: 2 * pedigree.triple_count() as u32,
            num_samples,
            global_sample_ids,
            sample_map,
            trios,
            genotypes,
            pos_map,
            mutation_cost,
        })
    }

    /// `posMap[position]` (operator[] inserts 0 for unknown positions)
    fn pos(&self, position: i32) -> usize {
        self.pos_map.get(&(position as u32)).copied().unwrap_or(0)
    }

    fn sample(&self, read: &Read) -> usize {
        self.sample_map.get(&(read.sample_id as u32)).copied().unwrap_or(0)
    }

    fn update_solution(new_sol: &mut PedSolution, old_sol: &PedSolution, old_idx: u16) {
        if new_sol.score > old_sol.score {
            new_sol.score = old_sol.score;
            new_sol.bt_row = old_idx;
            new_sol.balances =
                old_sol.balances.iter().map(|a| if a.is_empty() { Vec::new() } else { a[1..].to_vec() }).collect();
        }
    }

    fn recombination_cost(&self, t1: u32, t2: u32, p: usize) -> f32 {
        self.recomb_cost[p] * (t1 ^ t2).count_ones() as f32
    }

    #[expect(clippy::needless_range_loop, reason = "mirrors the C++ loop over columns of several balance vectors")]
    fn mutation_cost_of(&self, balances: &[Balance], t: u32, p: usize, allow_flips: bool, ahead: usize) -> f32 {
        let mut cost = 0f32;
        let last = ahead.min(balances[0].len().wrapping_sub(1));
        let mc = self.mutation_cost[p];
        for i in 0..=last {
            for (k, trio) in self.trios.iter().enumerate() {
                let m2c = ((t >> (2 * k)) & 1) as usize;
                let f2c = ((t >> (2 * k + 1)) & 1) as usize;
                let cm = balances[2 * trio[2]][i];
                let cf = balances[2 * trio[2] + 1][i];
                let m = balances[2 * trio[0] + m2c][i];
                let f = balances[2 * trio[1] + f2c][i];
                if allow_flips {
                    if cm * m < 0.0 {
                        cost += min(mc, min(cm.abs(), m.abs()));
                    }
                    if cf * f < 0.0 {
                        cost += min(mc, min(cf.abs(), f.abs()));
                    }
                } else {
                    cost += b2f(cm * m < 0.0) * mc;
                    cost += b2f(cf * f < 0.0) * mc;
                }
            }
        }
        cost
    }

    /// `getOptPhasing`: the cost of the best phasing of position `p` and, if
    /// requested, the phasing (two alleles per sample) and mutated alleles.
    fn opt_phasing(&self, balances: &[f32], t: u32, p: usize, mut out: Option<(&mut [i8], &mut [bool])>) -> f32 {
        let ns = self.num_samples;
        let mc = self.mutation_cost[p];
        let mut phase_cost = vec![[0f32; 5]; ns];
        for (s, pc) in phase_cost.iter_mut().enumerate() {
            let (a0, a1) = (balances[2 * s], balances[2 * s + 1]);
            pc[0] = a0 * b2f(a0 > 0.0) + a1 * b2f(a1 > 0.0);
            pc[1] = -a0 * b2f(a0 < 0.0) + a1 * b2f(a1 > 0.0);
            pc[2] = a0 * b2f(a0 > 0.0) - a1 * b2f(a1 < 0.0);
            pc[3] = -a0 * b2f(a0 < 0.0) - a1 * b2f(a1 < 0.0);
            let mut largest = 0;
            for i in 1..4 {
                if pc[largest] < pc[i] {
                    largest = i;
                }
            }
            pc[4] = pc[largest];
        }
        // allowed phasings per sample: 0 = 0|0, 1 = 0|1, 2 = 1|0, 3 = 1|1
        let mut phases: Vec<Vec<usize>> = vec![Vec::new(); ns];
        for s in 0..ns {
            if self.distrust_genotypes {
                for i in 0..4 {
                    if phase_cost[s][i] < phase_cost[s][4] + 2.0 * mc {
                        phases[s].push(i);
                    }
                }
            } else {
                match self.genotypes[s][p] {
                    0 => phases[s].push(0),
                    2 => phases[s].push(3),
                    _ => phases[s].extend([1, 2]),
                }
            }
        }
        let mut min_cost = f32::INFINITY;
        let mut v = vec![0usize; ns];
        let mut mutated = vec![false; 2 * ns];
        while v[ns - 1] < phases[ns - 1].len() {
            let mut cost = 0f32;
            mutated.fill(false);
            for (k, trio) in self.trios.iter().enumerate() {
                let m2c = (t >> (2 * k)) & 1;
                let f2c = (t >> (2 * k + 1)) & 1;
                let child = phases[trio[2]][v[trio[2]]];
                let acm = child & 1;
                let acf = (child & 2) >> 1;
                let am = (phases[trio[0]][v[trio[0]]] & (1 + m2c as usize)) >> m2c;
                let af = (phases[trio[1]][v[trio[1]]] & (1 + f2c as usize)) >> f2c;
                cost += b2f(am != acm) * mc;
                cost += b2f(af != acf) * mc;
                mutated[2 * trio[2]] = am != acm;
                mutated[2 * trio[2] + 1] = af != acf;
            }
            for s in 0..ns {
                cost += phase_cost[s][phases[s][v[s]]];
            }
            if cost < min_cost {
                min_cost = cost;
                if let Some((phasing, muts)) = out.as_mut() {
                    for s in 0..ns {
                        phasing[2 * s] = (phases[s][v[s]] & 1) as i8;
                        phasing[2 * s + 1] = ((phases[s][v[s]] & 2) >> 1) as i8;
                    }
                    muts.copy_from_slice(&mutated);
                }
            }
            v[0] += 1;
            for j in 0..ns - 1 {
                if v[j] >= phases[j].len() {
                    v[j] = 0;
                    v[j + 1] += 1;
                }
            }
        }
        min_cost
    }

    /// `addBalance` for `balances[basis]`, with `balances[co_basis]` the other
    /// partition of the same sample.
    fn add_balance(&self, balances: &mut [Balance], basis: usize, co_basis: usize, add: &[f32], target: &[i32]) -> f32 {
        let mut penalty = 0f32;
        for i in 0..add.len() {
            let b = balances[basis][i];
            if self.distrust_genotypes {
                if b * add[i] < 0.0 {
                    penalty += min(b.abs(), add[i].abs());
                }
            } else if target[i] == 1 {
                let co = balances[co_basis][i];
                if add[i] <= 0.0 {
                    penalty += min(-add[i], max(b - co, 0.0));
                } else {
                    penalty += min(add[i], max(co - b, 0.0));
                }
            } else {
                penalty += add[i].abs() * b2f(add[i] * ((target[i] - 1) as f32) < 0.0);
            }
            balances[basis][i] += add[i];
        }
        penalty
    }

    fn extend_solutions(&self, sols: &mut Vec<PedSolution>, to_ext: usize, p: usize) {
        let ms = self.mutation_cost_of(&sols[to_ext].balances, sols[to_ext].trans, p, false, 0);
        sols[to_ext].mutation_score = ms;
        if ms > 0.0 {
            let limit = 2f64.powi(self.tm_bits as i32);
            let mut t: u32 = 0;
            while (t as f64) < limit {
                let base = &sols[to_ext];
                if t != base.trans {
                    let rc = self.recombination_cost(base.trans, t, p);
                    #[allow(clippy::neg_cmp_op_on_partial_ord, reason = "C++ `continue`s on >= (NaN falls through)")]
                    if !(rc >= ms) {
                        let m = self.mutation_cost_of(&base.balances, t, p, false, 0);
                        #[allow(clippy::neg_cmp_op_on_partial_ord, reason = "as above")]
                        if !(m + rc >= ms) {
                            let new = PedSolution::extended(
                                &base.bp,
                                &base.bp_new,
                                t,
                                base.score + rc,
                                m,
                                base.bt_row,
                                &base.balances,
                            );
                            sols.push(new);
                        }
                    }
                }
                t += 1;
            }
        }
    }

    fn filter_solutions(&self, sols: &mut Vec<PedSolution>) {
        let mut scores: Vec<f32> = sols.iter().map(|s| s.score + s.mutation_score).collect();
        scores.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let too_high =
            if scores.len() > self.row_limit as usize { scores[self.row_limit as usize] } else { f32::INFINITY };
        let mut kept = Vec::new();
        for (i, s) in sols.iter().enumerate() {
            let score = s.score + s.mutation_score;
            if (score < too_high || score == scores[0]) && kept.len() < MAX_ROW_LIMIT as usize {
                kept.push(i);
            }
        }
        for (k, &i) in kept.iter().enumerate() {
            if k != i {
                sols[k] = sols[i].clone();
            }
        }
        sols.truncate(kept.len());
    }

    fn solve(&self) -> Result<Solution> {
        let rs = self.rs;
        let m = rs.len();
        let n = self.positions.len();
        let ns = self.num_samples;

        let mut start_index = vec![0usize];
        let mut q = 0;
        for p in 0..n {
            while q < m && self.pos(rs.get(q).first_position()) <= p {
                q += 1;
            }
            start_index.push(q);
        }
        // the first read of a not yet seen non-child sample always goes into
        // partition 0
        let mut seen = vec![false; ns];
        for trio in &self.trios {
            seen[trio[2]] = true;
        }

        let mut last_col: Vec<PedSolution> = vec![PedSolution::new(Vec::new(), 0, 0.0, ns)];
        let mut active: Vec<usize> = Vec::new();
        let mut m_bt: Vec<Vec<u16>> = vec![Vec::new(); n];
        let mut m_bp: Vec<(usize, Vec<bool>)> = Vec::new();
        let mut m_tm: Vec<Vec<bool>> = Vec::new();

        let mut right = 0;
        for p in 0..n {
            let active_last = std::mem::take(&mut active);
            let mut kept = Vec::new();
            for (i, &r) in active_last.iter().enumerate() {
                if self.pos(rs.get(r).last_position()) >= p {
                    active.push(r);
                    kept.push(i);
                }
            }

            // copy bipartitions without the reads that ended
            let mut sols: Vec<PedSolution> = Vec::new();
            let mut buckets: HashMap<u64, Vec<u16>> = HashMap::new();
            for (i, last) in last_col.iter().enumerate() {
                let mut b = vec![false; kept.len()];
                let mut bucket = last.trans as u64;
                for (a, &k) in kept.iter().enumerate() {
                    let bo = last.bp[k];
                    b[a] = bo;
                    bucket = (bucket << 1).wrapping_add(bo as u64);
                }
                let list = buckets.entry(bucket).or_default();
                let mut duplicate = false;
                for &j in list.iter() {
                    let j = j as usize;
                    if last.trans == sols[j].trans && b == sols[j].bp {
                        duplicate = true;
                        Self::update_solution(&mut sols[j], last, i as u16);
                        break;
                    }
                }
                if !duplicate {
                    list.push(sols.len() as u16);
                    let mut sol = PedSolution::new(b, last.trans, f32::INFINITY, ns);
                    Self::update_solution(&mut sol, last, i as u16);
                    sols.push(sol);
                }
            }

            right = right.max(p);
            for r in start_index[p]..start_index[p + 1] {
                right = right.max(self.pos(rs.get(r).last_position()));
            }
            let width = right + 1 - p;
            for s in sols.iter_mut() {
                for b in s.balances.iter_mut() {
                    b.resize(width, 0.0);
                }
            }

            // summarize equal reads
            let num_new = start_index[p + 1] - start_index[p];
            let mut equal_to: Vec<i64> = vec![-1; num_new];
            let mut balances: Vec<Balance> = Vec::with_capacity(num_new);
            let mut sample_ids: Vec<usize> = Vec::with_capacity(num_new);
            for i in 0..num_new {
                let r = start_index[p] + i;
                active.push(r);
                let read = rs.get(r);
                let mut b = vec![0f32; width];
                sample_ids.push(self.sample(read));
                for v in &read.variants {
                    let o = self.pos(v.position).wrapping_sub(p);
                    let a = v.allele as i8 as i32;
                    let q = v.quality as f32;
                    b[o] += q * a as f32 - q * (1 - a) as f32;
                }
                for j in 0..i {
                    if equal_to[j] != -1 || sample_ids[j] != sample_ids[i] {
                        continue;
                    }
                    let equal =
                        (0..width).all(|k| !(balances[j][k] * b[k] < 0.0 || (balances[j][k] != 0.0) != (b[k] != 0.0)));
                    if equal {
                        equal_to[i] = j as i64;
                        for k in 0..width {
                            balances[j][k] += b[k];
                        }
                        break;
                    }
                }
                balances.push(b);
            }

            // two extensions of every solution for every new read
            for i in 0..num_new {
                let balance = &balances[i];
                let sample = sample_ids[i];
                let target = &self.genotypes[sample][p..=right];
                let sol_end = sols.len();
                for sol in 0..sol_end {
                    if equal_to[i] >= 0 {
                        let same = sols[sol].bp_new[equal_to[i] as usize];
                        sols[sol].bp_new.push(same);
                        continue;
                    }
                    let mut useful = false;
                    if self.distrust_genotypes {
                        let mut j = 0;
                        while j < balance.len() && !useful {
                            let s0 = sols[sol].balances[sample * 2][j];
                            let s1 = sols[sol].balances[sample * 2 + 1][j];
                            useful |= (balance[j] != 0.0 && s0 * s1 < 0.0)
                                || ((balance[j] + s0) * s0 <= 0.0 && (balance[j] + s1) * s1 <= 0.0);
                            j += 1;
                        }
                    } else {
                        let mut j = 0;
                        while j < balance.len() && !useful {
                            useful |= self.genotypes[sample][p + j] == 1 && balance[j] != 0.0;
                            j += 1;
                        }
                    }
                    let mut sol1 = 0;
                    if seen[sample] {
                        let mut s1 = sols[sol].clone();
                        s1.score += self.add_balance(&mut s1.balances, 2 * sample + 1, 2 * sample, balance, target);
                        s1.mutation_score = self.mutation_cost_of(&s1.balances, s1.trans, p, true, 5);
                        s1.bp_new.push(true);
                        sols.push(s1);
                        sol1 = sols.len() - 1;
                    }
                    let s0 = &mut sols[sol];
                    let add = self.add_balance(&mut s0.balances, 2 * sample, 2 * sample + 1, balance, target);
                    s0.score += add;
                    s0.mutation_score = self.mutation_cost_of(&s0.balances, s0.trans, p, true, 5);
                    s0.bp_new.push(false);
                    if sol1 != 0 && !useful {
                        if sols[sol].score + sols[sol].mutation_score > sols[sol1].score + sols[sol1].mutation_score {
                            sols[sol] = sols[sol1].clone();
                        }
                        sols.pop();
                    }
                }
                seen[sample] = true;
                if sols.len() > self.row_limit as usize {
                    self.filter_solutions(&mut sols);
                }
            }
            // extend transmissions
            let sol_end = sols.len();
            for i in 0..sol_end {
                self.extend_solutions(&mut sols, i, p);
            }
            if sols.len() > self.row_limit as usize {
                self.filter_solutions(&mut sols);
            }
            // add costs for mutations
            for s in sols.iter_mut() {
                let firsts: Vec<f32> = s.balances.iter().map(|b| b[0]).collect();
                s.score += self.opt_phasing(&firsts, s.trans, p, None);
            }

            last_col.clear();
            let mut bt_vector = Vec::new();
            let mut tm_vector = Vec::new();
            for mut s in sols {
                s.finalize();
                m_bt[p].push(s.bt_row);
                bt_vector.extend_from_slice(&s.bp_new);
                for j in 0..self.tm_bits {
                    tm_vector.push((s.trans >> j) & 1 == 1);
                }
                last_col.push(s);
            }
            m_bp.push((num_new, bt_vector));
            m_tm.push(tm_vector);
        }

        // best score in the last column
        let mut opt_trans = vec![0u32; n];
        let mut best = f32::INFINITY;
        let mut ri: usize = 0;
        for (i, s) in last_col.iter().enumerate() {
            if s.score < best {
                best = s.score;
                ri = i;
            }
        }

        // backtracking
        let mut opt_bipart = vec![false; m];
        let tm_bits = self.tm_bits as usize;
        for p in (0..n).rev() {
            let offset = start_index[p];
            let new_count = m_bp[p].0;
            let current = &m_bp[p].1[new_count * ri..new_count * (ri + 1)];
            opt_bipart[offset..offset + current.len()].copy_from_slice(current);
            let tm = &m_tm[p][tm_bits * ri..tm_bits * (ri + 1)];
            for t in (0..tm_bits).rev() {
                opt_trans[p] = (opt_trans[p] << 1) + tm[t] as u32;
            }
            ri = m_bt[p][ri] as usize;
        }

        // allele votes
        let mut balances = vec![vec![0f32; 2 * ns]; n];
        for (r, read) in rs.iter().enumerate() {
            let sid = self.sample(read);
            for v in &read.variants {
                let a = v.allele as i8 as i32;
                let q = v.quality as f32;
                if a >= 0 {
                    balances[self.pos(v.position)][2 * sid + opt_bipart[r] as usize] += (2 * a - 1) as f32 * q;
                }
            }
        }

        // optimal phasing per position
        let mut opt_haps = vec![[vec![-1i8; n], vec![-1i8; n]]; ns];
        let mut phasing = vec![0i8; 2 * ns];
        let mut mutated = vec![false; 2 * ns];
        for p in 0..n {
            phasing.fill(0);
            mutated.fill(false);
            self.opt_phasing(&balances[p], opt_trans[p], p, Some((&mut phasing, &mut mutated)));
            for (sid, haps) in opt_haps.iter_mut().enumerate() {
                haps[0][p] = phasing[2 * sid];
                haps[1][p] = phasing[2 * sid + 1];
            }
        }

        let mut superreads = Vec::with_capacity(ns);
        for (sid, haps) in opt_haps.iter().enumerate() {
            let sample_id = self.global_sample_ids[sid] as i32;
            let mut read0 = Read::new("superread_0", -1, -1, sample_id);
            let mut read1 = Read::new("superread_1", -1, -1, sample_id);
            for (p, &position) in self.positions.iter().enumerate() {
                read0.add_variant(position as i32, haps[0][p] as i32, 30);
                read1.add_variant(position as i32, haps[1][p] as i32, 30);
            }
            let mut set = ReadSet::new();
            set.add(read0)?;
            set.add(read1)?;
            superreads.push(set);
        }
        Ok(Solution {
            // getOptScore() returns optScore, which solve() never sets
            optimal_score: 0,
            superreads,
            transmission_vector: opt_trans,
            partitioning: opt_bipart.iter().map(|&x| if x { 0 } else { 1 }).collect(),
        })
    }
}

/// Phase with `PedMecHeuristic` as `phase --algorithm heuristic` does.
/// `read_set` must be sorted.
pub fn solve(
    read_set: &ReadSet,
    recombcost: &[u32],
    pedigree: &Pedigree,
    distrust_genotypes: bool,
    positions: &[i32],
    row_limit: u32,
) -> Result<Solution> {
    Heuristic::new(read_set, recombcost, pedigree, distrust_genotypes, positions, row_limit)?.solve()
}

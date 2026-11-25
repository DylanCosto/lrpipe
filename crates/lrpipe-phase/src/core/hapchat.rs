//! Port of `hapchat/*` (`phase --algorithm hapchat`, HapCHAT by Beretta,
//! Dell'Acqua, Patterson, Pirola and Zaccaria).
//!
//! As in WhatsHap, only the first block of overlapping reads is phased (the
//! C++ code calls `has_next_block()` once), and the super reads cover the
//! positions of all reads: for positions after the first block, the C++ code
//! reads past the end of the haplotype buffers (undefined behaviour); this
//! port uses allele 0 there.
//!
//! `std::bitset<64>` is a `u64`, `Cost` saturates at `u32::MAX`, binomial
//! coefficients wrap around like the C++ `unsigned int`s, and `<<` on
//! `unsigned int` takes the shift count modulo 32 as on x86.

#![expect(clippy::too_many_arguments, reason = "mirrors the C++ functions")]

use std::sync::OnceLock;

use anyhow::{Result, bail};

use super::dp::{ColumnIterator, Solution};
use super::read::{ALT_ALLELE, BLANK, REF_ALLELE, Read, ReadSet};

const MAX_COVERAGE: usize = 64;
const MAX_COUNTER: u32 = u32::MAX;
const INFTY: u32 = u32::MAX;
/// HapChatCore's fixed parameters
const ERROR_RATE: f64 = 0.05;
const ALPHA: f64 = 0.01;

type BitColumn = u64;

#[derive(Clone, Copy, Debug)]
struct HEntry {
    read_id: i32,
    allele: i32,
    phred: u32,
}

const BLANK_ENTRY: HEntry = HEntry { read_id: -1, allele: BLANK, phred: 0 };

type Column = Vec<HEntry>;

#[derive(Clone, Copy)]
struct Backtrace1 {
    jump: i32,
    index: u32,
}

impl Default for Backtrace1 {
    fn default() -> Self {
        Backtrace1 { jump: -1, index: 0 }
    }
}

/// `Cost::operator+`
fn add(a: u32, b: u32) -> u32 {
    a.saturating_add(b)
}

fn test(b: BitColumn, i: usize) -> bool {
    i < 64 && (b >> i) & 1 == 1
}

/// `bitset::set(pos, value)`, which throws for positions past the end.
fn set(b: &mut BitColumn, i: usize, value: bool) -> Result<()> {
    if i >= 64 {
        bail!("bitset::set: __position (which is {i}) >= _Nb (which is 64)");
    }
    if value {
        *b |= 1 << i;
    } else {
        *b &= !(1 << i);
    }
    Ok(())
}

fn shr(b: BitColumn, s: u32) -> BitColumn {
    b.checked_shr(s).unwrap_or(0)
}

/// `ffsl`: 1 + index of the lowest set bit (`b` is not 0).
fn ffsl(b: BitColumn) -> u32 {
    b.trailing_zeros() + 1
}

/// The low `len` bits (`(ones << len).flip()`).
fn low_bits(len: u32) -> BitColumn {
    if len >= 64 { !0 } else { (1u64 << len) - 1 }
}

/// `BinomialCoefficient` initialized for n = k = MAX_COVERAGE.
struct Binomial {
    b: Vec<Vec<u32>>,
    c: Vec<Vec<u32>>,
}

fn binomials() -> &'static Binomial {
    static TABLE: OnceLock<Binomial> = OnceLock::new();
    TABLE.get_or_init(|| {
        let (n, k) = (MAX_COVERAGE, MAX_COVERAGE);
        let mut b = vec![vec![0u32; n + 1]; n + 1];
        for i in 0..=n {
            for j in 0..=i {
                b[i][j] = if i == 0 || j == 0 || j == i { 1 } else { b[i - 1][j - 1].wrapping_add(b[i - 1][j]) };
            }
        }
        let c: Vec<Vec<u32>> = (0..=n)
            .map(|i| (0..=k).map(|j| b[i][..=j].iter().fold(0u32, |sum, &x| sum.wrapping_add(x))).collect())
            .collect();
        Binomial { b, c }
    })
}

impl Binomial {
    fn binomial(&self, n: u32, k: u32) -> u32 {
        self.b[n as usize][k as usize]
    }

    fn cumulative(&self, n: u32, k: u32) -> u32 {
        self.c[n as usize][k as usize]
    }

    fn indexof(&self, mut comb: BitColumn) -> u32 {
        let (mut k, mut c_k, mut result) = (0i32, 0i32, 0u32);
        while comb != 0 {
            let temp = ffsl(comb);
            c_k += temp as i32;
            k += 1;
            result = result.wrapping_add(self.binomial((c_k - 1) as u32, k as u32));
            comb = shr(comb, temp);
        }
        result
    }

    fn cumulative_indexof(&self, comb: BitColumn, n_elements: u32) -> u32 {
        let k = comb.count_ones();
        let mut result = self.indexof(comb);
        for i in 0..k {
            result = result.wrapping_add(self.binomial(n_elements, i));
        }
        result
    }
}

/// `computeK`: the number of corrections allowed for a coverage, computed
/// once for coverages up to MAX_COVERAGE with alpha 0.01 and error rate 0.05.
fn compute_k(cov: u32) -> u32 {
    static KS: OnceLock<Vec<u32>> = OnceLock::new();
    let ks = KS.get_or_init(|| {
        let bin = binomials();
        let mut ks = vec![0u32; MAX_COVERAGE + 1];
        for (i, slot) in ks.iter_mut().enumerate().skip(1) {
            let i = i as u32;
            let mut k = 0;
            let mut cumulative = (1.0 - ERROR_RATE).powf(i as f64);
            while 1.0 - cumulative > ALPHA && k < i {
                k += 1;
                cumulative +=
                    bin.binomial(i, k) as f64 * ERROR_RATE.powf(k as f64) * (1.0 - ERROR_RATE).powf((i - k) as f64);
            }
            *slot = k;
        }
        ks
    });
    ks[cov as usize]
}

/// The `Combinations` generator, as used by the DP (`initialize_cumulative`).
struct Combinations {
    n: u32,
    t: u32,
    count: u32,
    end: bool,
    combination: Vec<i32>,
    result: BitColumn,
    ending_mask: BitColumn,
    j: i32,
    x: i32,
    cumulative: bool,
    t_cumulative: u32,
    first: bool,
}

impl Combinations {
    fn new() -> Self {
        Combinations {
            n: 0,
            t: 0,
            count: 0,
            end: false,
            combination: vec![-1; 2],
            result: 0,
            ending_mask: 0,
            j: 0,
            x: 0,
            cumulative: false,
            t_cumulative: 0,
            first: false,
        }
    }

    fn initialize_cumulative(&mut self, n: u32, t: u32) {
        self.n = n;
        self.t = 0;
        self.count = u32::MAX;
        self.end = false;
        self.combination.clear();
        self.combination.resize(t as usize + 2, -1);
        self.result = 0;
        self.ending_mask = 0;
        self.j = 0;
        self.x = 0;
        self.cumulative = true;
        self.t_cumulative = t;
        self.first = true;
    }

    fn has_next(&self) -> bool {
        !self.end
    }

    fn next(&mut self) -> Result<()> {
        if !self.end {
            if self.first {
                self.first = false;
                self.initializing_step()?;
            } else {
                self.basic_step()?;
            }
            if self.cumulative && self.end && self.t != self.t_cumulative {
                self.end = false;
                self.t += 1;
                self.combination.clear();
                self.combination.resize(self.t as usize + 2, -1);
                self.ending_mask = 0;
                self.j = 0;
                self.x = 0;
                self.first = true;
            }
        }
        Ok(())
    }

    fn combination(&self) -> BitColumn {
        self.result
    }

    fn initializing_step(&mut self) -> Result<()> {
        let t = self.t as usize;
        for j in 1..=t {
            self.combination[j - 1] = j as i32 - 1;
            set(&mut self.ending_mask, (self.n as usize).wrapping_sub(j), true)?;
        }
        self.combination[t] = self.n as i32;
        self.combination[t + 1] = 0;
        self.j = t as i32;
        self.make_combination()?;
        self.check_end();
        Ok(())
    }

    fn basic_step(&mut self) -> Result<()> {
        let c = &mut self.combination;
        if self.j > 0 {
            self.x = self.j;
            c[self.j as usize - 1] = self.x;
            self.j -= 1;
        } else if c[0] + 1 < c[1] {
            c[0] += 1;
        } else {
            self.j = 2;
            loop {
                let j = self.j as usize;
                c[j - 2] = j as i32 - 2;
                self.x = c[j - 1] + 1;
                if self.x == c[j] {
                    self.j += 1;
                } else {
                    break;
                }
            }
            if self.j > self.t as i32 {
                self.end = true;
            } else {
                c[self.j as usize - 1] = self.x;
                self.j -= 1;
            }
        }
        self.make_combination()?;
        self.check_end();
        Ok(())
    }

    fn make_combination(&mut self) -> Result<()> {
        self.result = 0;
        for k in 0..self.t as usize {
            set(&mut self.result, self.combination[k] as usize, true)?;
        }
        self.count = self.count.wrapping_add(1);
        Ok(())
    }

    fn check_end(&mut self) {
        self.end = self.result == self.ending_mask;
    }
}

/// `HapChatColumnIterator`, with the read blocks computed by `set_block`.
fn read_blocks(read_set: &ReadSet) -> Vec<Vec<usize>> {
    let mut min: Vec<u32> = Vec::new();
    let mut max: Vec<u32> = Vec::new();
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for (i, read) in read_set.iter().enumerate() {
        let mut minn = read.first_position() as u32;
        let mut maxx = read.last_position() as u32;
        if min.is_empty() {
            min.push(minn);
            max.push(maxx);
            blocks.push(vec![i]);
            continue;
        }
        for j in 0..min.len() {
            let over = minn < min[j] && maxx > max[j];
            if (minn >= min[j] && minn <= max[j]) || (maxx >= min[j] && maxx <= max[j]) || over {
                minn = min[j].min(minn);
                maxx = max[j].max(maxx);
                min[j] = minn;
                max[j] = maxx;
                blocks[j].push(i);
                // (a read starting at position 0 is also never a new block)
                minn = 0;
                break;
            }
        }
        if minn != 0 {
            min.push(minn);
            max.push(maxx);
            blocks.push(vec![i]);
        }
    }
    blocks
}

/// The column source of the DP: the block's `ColumnIterator`, shared by all
/// copies of `HapChatColumnIterator` in C++ (each copy has its own `end`).
struct Columns<'a> {
    it: ColumnIterator<'a>,
}

impl Columns<'_> {
    fn has_next(&self) -> bool {
        self.it.has_next()
    }

    /// `get_column`; an empty column (and `end` set) once exhausted.
    fn get_column(&mut self, end: &mut bool) -> Column {
        if self.it.has_next() {
            self.it
                .get_next()
                .into_iter()
                .map(|e| HEntry { read_id: e.read_id as i32, allele: e.allele, phred: e.quality })
                .collect()
        } else {
            *end = true;
            Vec::new()
        }
    }

    fn reset(&mut self) {
        self.it.jump_to_column(0);
    }

    fn column_count(&self) -> usize {
        self.it.column_count()
    }
}

fn next(pointer: i32, total: usize, shift: i64) -> i32 {
    ((pointer as i64 + shift) % total as i64) as i32
}

fn prev(pointer: i32, total: usize, shift: i64) -> i32 {
    ((pointer as i64 + total as i64 - shift) % total as i64) as i32
}

fn u(p: i32) -> usize {
    p as usize
}

/// `compute_active_common`: common reads of the two columns and the number
/// of them that are gaps in `col_j`.
fn compute_active_common(col_j: &Column, col_q: &Column) -> (u32, u32) {
    let (mut active, mut gaps) = (0u32, 0u32);
    let (mut i, mut j) = (0, 0);
    while i < col_q.len() && j < col_j.len() && col_j[j].read_id != -1 && col_q[i].read_id != -1 {
        if col_q[i].read_id == col_j[j].read_id {
            if col_j[j].allele == BLANK {
                gaps += 1;
            }
            i += 1;
            j += 1;
            active += 1;
        } else if col_q[i].read_id < col_j[j].read_id {
            i += 1;
        } else {
            j += 1;
        }
    }
    (active, gaps)
}

struct Dp<'a> {
    bin: &'static Binomial,
    columns: Columns<'a>,
}

/// State of the DP's input window, updated by `insert_col_and_update`.
struct Window {
    input: Vec<Column>,
    k_j: Vec<u32>,
    homo_cost: Vec<u32>,
    homo_weight: Vec<u32>,
}

impl Window {
    fn insert_col_and_update(&mut self, pointer: usize, column: &Column, kind_homozygous: &mut [bool], step: usize) {
        let (mut count_major, mut weight_major, mut count_minor, mut weight_minor) = (0u32, 0u32, 0u32, 0u32);
        let slot = &mut self.input[pointer];
        let mut i = 0;
        while i < column.len() {
            let e = column[i];
            let gap = e.allele == BLANK;
            slot[i] = HEntry { read_id: e.read_id, allele: if gap { REF_ALLELE } else { e.allele }, phred: e.phred };
            if !gap {
                if e.allele == ALT_ALLELE {
                    count_minor += 1;
                    weight_minor = add(weight_minor, e.phred);
                } else if e.allele == REF_ALLELE {
                    count_major += 1;
                    weight_major = add(weight_major, e.phred);
                }
            }
            i += 1;
        }
        if i < slot.len() && slot[i].read_id != -1 {
            slot[i] = BLANK_ENTRY;
        }
        let k = compute_k(count_minor + count_major);
        self.k_j[pointer] = k;
        self.homo_cost[pointer] = MAX_COUNTER;
        self.homo_weight[pointer] = INFTY;
        if count_minor <= k && weight_minor < self.homo_weight[pointer] {
            self.homo_cost[pointer] = count_minor;
            self.homo_weight[pointer] = weight_minor;
            if step < kind_homozygous.len() {
                kind_homozygous[step] = true;
            }
        }
        if count_major <= k && weight_major < self.homo_weight[pointer] {
            self.homo_cost[pointer] = count_major;
            self.homo_weight[pointer] = weight_major;
            if step < kind_homozygous.len() {
                kind_homozygous[step] = false;
            }
        }
        self.homo_cost[pointer] = MAX_COVERAGE as u32 + 1;
    }
}

/// Inputs of the DP computed by `computeInputParams`.
struct InputParams {
    num_cols: usize,
    max_cov: u32,
    max_l: u32,
    max_k: u32,
    sum_successive_l: Vec<u32>,
    scheme_backtrace: Vec<Vec<u32>>,
}

impl Dp<'_> {
    fn compute_input_params(&mut self) -> Result<InputParams> {
        let bin = self.bin;
        let mut end = false;
        self.columns.reset();
        let num_cols = self.columns.column_count() + 1;
        let mut input: Vec<Column> = vec![Vec::new(); num_cols];
        let mut homo_cost = vec![0u32; num_cols];
        let mut k_j = vec![0u32; num_cols];
        let mut rows = vec![0u32; num_cols * MAX_COVERAGE];
        let (mut max_cov, mut max_k, mut _max_gaps) = (0u32, 0u32, 0u32);
        let mut idx = 0;
        loop {
            let read_column = if idx == 0 { Vec::new() } else { self.columns.get_column(&mut end) };
            if read_column.len() > MAX_COVERAGE {
                // the C++ code fails once such a column is put into a bitset
                bail!("bitset::set: __position (which is {MAX_COVERAGE}) >= _Nb (which is {MAX_COVERAGE})");
            }
            let (mut major, mut minor, mut gaps) = (0u32, 0u32, 0u32);
            for e in &read_column {
                if e.allele != BLANK {
                    if e.allele == REF_ALLELE {
                        major += 1;
                    } else if e.allele == ALT_ALLELE {
                        minor += 1;
                    } else {
                        // the C++ code prints an error and exits
                        bail!("ERROR: read invalid entry of type: {}", e.allele);
                    }
                } else {
                    gaps += 1;
                }
                rows[e.read_id as usize] += 1;
            }
            input[idx] = read_column;
            homo_cost[idx] = MAX_COVERAGE as u32 + 1;
            k_j[idx] = compute_k(minor + major);
            max_cov = max_cov.max(input[idx].len() as u32);
            max_k = max_k.max(k_j[idx]);
            _max_gaps = _max_gaps.max(gaps);
            idx += 1;
            if !(idx != num_cols && self.columns.has_next()) {
                break;
            }
        }
        let mut max_l = rows.iter().copied().max().unwrap_or(0).max(2);
        let mut max_cons_homo = 0;
        let mut sum_successive_l = vec![0u32; max_l as usize];
        let mut scheme_backtrace = vec![Vec::new(); num_cols];
        for i in 0..input.len() {
            let mut y = 1;
            let k_temp = k_j[i];
            let mut current_cons_homo = 0;
            let mut flag = true;
            scheme_backtrace[i].push(0);
            while y < max_l as usize && i + y < input.len() {
                let (active, common_gaps) = compute_active_common(&input[i], &input[i + y]);
                let result = bin.cumulative(active.wrapping_sub(common_gaps), k_temp).wrapping_shl(common_gaps);
                sum_successive_l[y] = sum_successive_l[y].max(result);
                if flag {
                    if homo_cost[i + y] <= k_j[i + y] && active != 0 {
                        current_cons_homo += 1;
                    } else {
                        flag = false;
                    }
                    scheme_backtrace[i].push(result);
                }
                y += 1;
            }
            max_cons_homo = max_cons_homo.max(current_cons_homo);
        }
        max_l = max_l.min(max_cons_homo + 2);
        Ok(InputParams { num_cols, max_cov, max_l, max_k, sum_successive_l, scheme_backtrace })
    }

    /// `dp`: the optimal cost and the two haplotypes of the block's columns.
    fn run(&mut self) -> Result<(u32, Vec<bool>, Vec<bool>)> {
        let bin = self.bin;
        let params = self.compute_input_params()?;
        let InputParams { num_cols: num_col, max_cov, max_l, max_k, sum_successive_l, scheme_backtrace } = params;
        let ml = max_l as usize;
        let size = 2 * (ml - 1) + 1;

        let mut w = Window {
            input: vec![vec![BLANK_ENTRY; max_cov as usize]; size],
            k_j: vec![max_k; size],
            homo_cost: vec![MAX_COUNTER; size],
            homo_weight: vec![INFTY; size],
        };
        let mut back_indexer = vec![vec![-1i32; max_cov as usize]; size];
        let mut forw_indexer = vec![vec![-1i32; max_cov as usize]; size];
        let indexer_pointer = max_l as i32 - 1;
        let mut pos_gaps = vec![0 as BitColumn; size];
        let mut num_pos_gaps = vec![0u32; size];
        let mut prevision: Vec<Vec<Vec<u32>>> =
            (0..ml).map(|_| (0..ml).map(|q| vec![INFTY; sum_successive_l[q] as usize]).collect()).collect();
        let mut prevision_pointer: i32 = 0;
        let mut opt = vec![INFTY; ml + 1];
        let mut opt_pointer: i32 = 0;

        let mut bt1: Vec<Vec<Vec<Backtrace1>>> = Vec::with_capacity(num_col);
        let mut bt2_hap: Vec<Vec<Vec<bool>>> = Vec::with_capacity(num_col);
        let mut bt2_nb: Vec<Vec<Vec<bool>>> = Vec::with_capacity(num_col);
        for scheme in scheme_backtrace.iter().take(num_col) {
            bt1.push(scheme.iter().map(|&n| vec![Backtrace1::default(); n as usize]).collect());
            bt2_hap.push(scheme.iter().map(|&n| vec![false; n as usize]).collect());
            bt2_nb.push(scheme.iter().map(|&n| vec![false; n as usize]).collect());
        }
        let mut is_homozygous = vec![false; num_col];
        let mut homo_haplotypes = vec![false; num_col];
        let mut best1 = vec![Backtrace1::default(); num_col];
        let mut best2_hap = vec![false; num_col];
        let mut best2_nb = vec![false; num_col];

        let mut generator = Combinations::new();
        self.columns.reset();
        let mut end = false;
        let mut step: usize = 0;

        // the first MAX_L columns
        let mut input_pointer: i32 = 0;
        let mut l = 0;
        let mut flag = true;
        loop {
            let new_l_pointer = next(input_pointer, size, l as i64);
            let column = if l == 0 {
                Vec::new()
            } else {
                flag = self.columns.has_next();
                self.columns.get_column(&mut end)
            };
            w.insert_col_and_update(u(new_l_pointer), &column, &mut homo_haplotypes, step + l);
            l += 1;
            if !(l < ml && flag) {
                break;
            }
        }

        // base case
        opt[u(opt_pointer)] = 0;
        w.k_j[u(input_pointer)] = 0;
        w.homo_weight[u(input_pointer)] = 0;
        w.homo_cost[u(input_pointer)] = 0;
        let mut current_cost: u32 = 0;
        let mut p: usize = 1;
        loop {
            let np = u(next(input_pointer, size, p as i64 - 1));
            let feasibility = p - 1 == 0 || w.homo_cost[np] <= w.k_j[np];
            if p >= ml || forw_indexer[u(indexer_pointer) + p][0] == -1 || !feasibility {
                break;
            }
            let npp = u(next(prevision_pointer, prevision.len(), p as i64));
            prevision[npp][p].resize(1, 0);
            prevision[npp][p][0] = current_cost;
            p += 1;
        }

        let mut k_j_inc = w.k_j[u(input_pointer)] as f64;
        let mut re_run_k = false;
        let mut solution_existence = true;
        let mut cov_j = 0u32;

        while !(!re_run_k && end && w.input[u(next(input_pointer, size, 1))][0].read_id == -1) {
            let mut current_best = INFTY;
            solution_existence = false;
            let mut temp_jump: u32 = u32::MAX;
            let mut temp_index: u32 = 0;
            let mut temp_haplotypes = false;
            let mut temp_new_block = false;
            step += 1;

            input_pointer = next(input_pointer, size, 1);
            let new_input_pointer = next(input_pointer, size, ml as i64 - 1);
            if !re_run_k {
                let column = self.columns.get_column(&mut end);
                w.insert_col_and_update(u(new_input_pointer), &column, &mut homo_haplotypes, step + ml - 1);
                k_j_inc = w.k_j[u(input_pointer)] as f64;
            }

            // indexers
            let mut q = 1;
            while q < ml {
                let qp = u(indexer_pointer) + q;
                let (cq, cj) = (u(next(input_pointer, size, q as i64)), u(input_pointer));
                intersect(
                    &w.input[cq],
                    &w.input[cj],
                    qp,
                    &mut forw_indexer,
                    &mut back_indexer,
                    &mut pos_gaps,
                    &mut num_pos_gaps,
                    false,
                )?;
                if forw_indexer[qp][0] == -1 {
                    for p in q + 1..ml {
                        forw_indexer[u(indexer_pointer) + p][0] = -1;
                        back_indexer[u(indexer_pointer) + p][0] = -1;
                    }
                    q = ml;
                }
                q += 1;
            }
            let mut q = 1;
            while q < ml {
                let qp = u(indexer_pointer) - q;
                let (cq, cj) = (u(prev(input_pointer, size, q as i64)), u(input_pointer));
                intersect(
                    &w.input[cq],
                    &w.input[cj],
                    qp,
                    &mut forw_indexer,
                    &mut back_indexer,
                    &mut pos_gaps,
                    &mut num_pos_gaps,
                    true,
                )?;
                if forw_indexer[qp][0] == -1 {
                    for p in q + 1..ml {
                        forw_indexer[u(indexer_pointer) - p][0] = -1;
                        back_indexer[u(indexer_pointer) - p][0] = -1;
                    }
                    q = ml;
                }
                q += 1;
            }

            // prevision
            prevision_pointer = next(prevision_pointer, prevision.len(), 1);
            let npp = u(next(prevision_pointer, prevision.len(), ml as i64 - 1));
            let last_input_pointer = next(input_pointer, size, ml as i64 - 1);
            for i in 1..prevision[npp].len() {
                let prec = u(prev(last_input_pointer, size, i as i64));
                let (active, common_gaps) = compute_active_common(&w.input[prec], &w.input[u(last_input_pointer)]);
                let combinations =
                    bin.cumulative(active.wrapping_sub(common_gaps), w.k_j[prec]).wrapping_shl(common_gaps) as usize;
                if re_run_k && prevision[npp][i].len() < combinations {
                    prevision[npp][i].resize(combinations, 0);
                }
                if re_run_k && i < bt1[step].len() && bt1[step][i].len() < combinations {
                    bt1[step][i].resize(combinations, Backtrace1::default());
                    bt2_hap[step][i].resize(combinations, false);
                    bt2_nb[step][i].resize(combinations, false);
                }
                for c in prevision[npp][i].iter_mut().take(combinations) {
                    *c = INFTY;
                }
            }

            opt_pointer = next(opt_pointer, opt.len(), 1);
            opt[u(opt_pointer)] = INFTY;

            let ip = u(input_pointer);
            let (colj, cov, gaps_mask, num_gaps) = represent_column(&w.input[ip])?;
            cov_j = cov;
            let proj = project(colj, gaps_mask, cov_j)?;
            let _ = proj;

            if w.homo_cost[ip] <= w.k_j[ip] {
                let temp = add(w.homo_weight[ip], opt[u(prev(opt_pointer, opt.len(), 1))]);
                if temp < opt[u(opt_pointer)] {
                    opt[u(opt_pointer)] = temp;
                    solution_existence = true;
                    is_homozygous[step] = true;
                }
            }

            generator.initialize_cumulative(cov_j - num_gaps, w.k_j[ip]);
            let mut loop_var = generator.has_next();
            while loop_var {
                generator.next()?;
                let comb_no_gaps = generator.combination();
                let mut comb_gaps_int: u32 = 0;
                loop {
                    let comb_gaps = comb_gaps_int as BitColumn;
                    let mask = make_mask(gaps_mask, cov_j, comb_gaps, comb_no_gaps)?;
                    current_cost = INFTY;
                    let weight_mask = compute_weight_mask(mask, &w.input[ip]);

                    let mut q: usize = 1;
                    let mut new_homo_pointer = u(prev(input_pointer, size, q as i64 - 1));
                    let mut cumulative_homo: u32 = 0;
                    let mut feasibility;
                    loop {
                        feasibility = q - 1 == 0 || w.homo_cost[new_homo_pointer] <= w.k_j[new_homo_pointer];
                        if q >= ml || forw_indexer[u(indexer_pointer) - q][0] == -1 || !feasibility {
                            break;
                        }
                        let iq = u(indexer_pointer) - q;
                        let new_q_pointer = u(prev(input_pointer, size, q as i64));
                        let (mut mask_qj, active_qj) = extract_common_mask(
                            &w.input[new_q_pointer],
                            &w.input[ip],
                            mask,
                            &back_indexer[iq],
                            &forw_indexer[iq],
                        )?;
                        let ungaps = (mask_qj & !pos_gaps[iq]).count_ones();
                        if ungaps <= w.k_j[new_q_pointer] {
                            let index = compute_index_of(bin, mask_qj, active_qj, num_pos_gaps[iq], pos_gaps[iq])?;
                            let temp = add(
                                add(prevision[u(prevision_pointer)][q][index as usize], weight_mask),
                                cumulative_homo,
                            );
                            if temp < current_cost {
                                current_cost = temp;
                                solution_existence = true;
                                temp_jump = q as u32;
                                temp_index = index;
                                temp_haplotypes = bt2_hap[step - q][q][index as usize];
                                temp_new_block = false;
                            }
                        }
                        mask_qj ^= low_bits(active_qj);
                        if active_qj.wrapping_sub(num_pos_gaps[iq]).wrapping_sub(ungaps) <= w.k_j[new_q_pointer] {
                            let index = compute_index_of(bin, mask_qj, active_qj, num_pos_gaps[iq], pos_gaps[iq])?;
                            let temp = add(
                                add(prevision[u(prevision_pointer)][q][index as usize], weight_mask),
                                cumulative_homo,
                            );
                            if temp < current_cost {
                                current_cost = temp;
                                solution_existence = true;
                                temp_jump = q as u32;
                                temp_index = index;
                                temp_haplotypes = !bt2_hap[step - q][q][index as usize];
                                temp_new_block = false;
                            }
                        }
                        q += 1;
                        new_homo_pointer = u(prev(input_pointer, size, q as i64 - 1));
                        cumulative_homo = add(cumulative_homo, w.homo_weight[new_homo_pointer]);
                    }

                    if q <= ml && feasibility {
                        let temp =
                            add(add(opt[u(prev(opt_pointer, opt.len(), q as i64))], weight_mask), cumulative_homo);
                        if temp < current_cost {
                            current_cost = temp;
                            solution_existence = true;
                            temp_jump = q as u32;
                            temp_index = 0;
                            temp_haplotypes = false;
                            temp_new_block = true;
                        }
                    }

                    // predictions for the following columns
                    let mut p: usize = 1;
                    loop {
                        let nhp = u(next(input_pointer, size, p as i64 - 1));
                        let feasibility = p - 1 == 0 || w.homo_cost[nhp] <= w.k_j[nhp];
                        if p >= ml || forw_indexer[u(indexer_pointer) + p][0] == -1 || !feasibility {
                            break;
                        }
                        let ipp = u(indexer_pointer) + p;
                        let (cut_mask, active_pj) = cut(mask, &back_indexer[ipp])?;
                        let index =
                            compute_index_of(bin, cut_mask, active_pj, num_pos_gaps[ipp], pos_gaps[ipp])? as usize;
                        let npp = u(next(prevision_pointer, prevision.len(), p as i64));
                        let temp = &mut prevision[npp][p][index];
                        if current_cost < *temp {
                            *temp = current_cost;
                            bt1[step][p][index] = Backtrace1 { jump: temp_jump as i32, index: temp_index };
                            bt2_hap[step][p][index] = temp_haplotypes;
                            bt2_nb[step][p][index] = temp_new_block;
                        }
                        p += 1;
                    }

                    if current_cost < current_best {
                        current_best = current_cost;
                        best1[step] = Backtrace1 { jump: temp_jump as i32, index: temp_index };
                        best2_hap[step] = temp_haplotypes;
                        best2_nb[step] = temp_new_block;
                    }
                    if current_cost < opt[u(opt_pointer)] {
                        opt[u(opt_pointer)] = current_cost;
                        is_homozygous[step] = false;
                    }
                    comb_gaps_int = comb_gaps_int.wrapping_add(1);
                    if comb_gaps_int >= (1i32.wrapping_shl(num_gaps)) as u32 {
                        break;
                    }
                }
                loop_var = generator.has_next();
            }

            // INC-K: allow more corrections if the column has no solution
            if solution_existence {
                re_run_k = false;
            } else {
                k_j_inc = if k_j_inc <= 0.0 { 1.0 } else { k_j_inc + k_j_inc.log2() + 1.0 };
                w.k_j[ip] = k_j_inc.floor() as u32;
                input_pointer = prev(input_pointer, size, 1);
                prevision_pointer = prev(prevision_pointer, prevision.len(), 1);
                opt_pointer = prev(opt_pointer, opt.len(), 1);
                step -= 1;
                re_run_k = true;
            }
        }

        if !solution_existence {
            // the C++ code calls exit(EXIT_FAILURE)
            bail!(
                "No feasible solution exists with these parameters -- alpha = {ALPHA} and error rate = {ERROR_RATE} (column {step}, coverage {cov_j})"
            );
        }
        let optimum = opt[u(opt_pointer)];
        let (h1, h2) = reconstruct_haplotypes(
            &bt1,
            &bt2_hap,
            &bt2_nb,
            &is_homozygous,
            &homo_haplotypes,
            &best1,
            &best2_hap,
            &best2_nb,
        );
        Ok((optimum, h1, h2))
    }
}

fn intersect(
    col_q: &Column,
    col_j: &Column,
    q: usize,
    forw_indexer: &mut [Vec<i32>],
    back_indexer: &mut [Vec<i32>],
    pos_gaps: &mut [BitColumn],
    num_pos_gaps: &mut [u32],
    q_is_back: bool,
) -> Result<()> {
    let (mut i, mut j, mut count) = (0, 0, 0);
    pos_gaps[q] = 0;
    num_pos_gaps[q] = 0;
    while i < col_q.len() && j < col_j.len() && col_j[j].read_id != -1 && col_q[i].read_id != -1 {
        if col_q[i].read_id == col_j[j].read_id {
            forw_indexer[q][count] = i as i32;
            back_indexer[q][count] = j as i32;
            let gap = if q_is_back { col_q[i].allele == BLANK } else { col_j[j].allele == BLANK };
            if gap {
                set(&mut pos_gaps[q], count, true)?;
                num_pos_gaps[q] += 1;
            }
            i += 1;
            j += 1;
            count += 1;
        } else if col_q[i].read_id < col_j[j].read_id {
            i += 1;
        } else {
            j += 1;
        }
    }
    if count < forw_indexer[q].len() {
        forw_indexer[q][count] = -1;
        back_indexer[q][count] = -1;
    }
    Ok(())
}

/// (column as bits (1 = ALT), coverage, gaps, number of gaps)
fn represent_column(column: &Column) -> Result<(BitColumn, u32, BitColumn, u32)> {
    let (mut result, mut gaps_mask, mut cov, mut num_gaps) = (0, 0, 0usize, 0u32);
    while cov < column.len() && column[cov].read_id != -1 {
        if column[cov].allele == BLANK {
            set(&mut gaps_mask, cov, true)?;
            num_gaps += 1;
        }
        set(&mut result, cov, column[cov].allele == ALT_ALLELE)?;
        cov += 1;
    }
    Ok((result, cov as u32, gaps_mask, num_gaps))
}

fn make_mask(mask_gaps: BitColumn, cov: u32, comb_gaps: BitColumn, comb_no_gaps: BitColumn) -> Result<BitColumn> {
    let mut mask = 0;
    let (mut i_no_gaps, mut i_gaps) = (0, 0);
    for i in 0..cov as usize {
        if test(mask_gaps, i) {
            set(&mut mask, i, test(comb_gaps, i_gaps))?;
            i_gaps += 1;
        } else {
            set(&mut mask, i, test(comb_no_gaps, i_no_gaps))?;
            i_no_gaps += 1;
        }
    }
    Ok(mask)
}

fn project(col: BitColumn, mask_out: BitColumn, cov: u32) -> Result<BitColumn> {
    let mut projection = 0;
    let mut j = 0;
    for i in 0..cov as usize {
        if !test(mask_out, i) {
            set(&mut projection, j, test(col, i))?;
            j += 1;
        }
    }
    Ok(projection)
}

fn compute_index_of(bin: &Binomial, mask: BitColumn, cov: u32, num_gaps: u32, pos_gaps: BitColumn) -> Result<u32> {
    let (mut comb_gaps, mut comb_no_gaps) = (0, 0);
    let (mut i_gaps, mut i_no_gaps) = (0, 0);
    for i in 0..cov as usize {
        if test(pos_gaps, i) {
            set(&mut comb_gaps, i_gaps, test(mask, i))?;
            i_gaps += 1;
        } else {
            set(&mut comb_no_gaps, i_no_gaps, test(mask, i))?;
            i_no_gaps += 1;
        }
    }
    Ok(bin.cumulative_indexof(comb_no_gaps, cov.wrapping_sub(num_gaps)).wrapping_shl(num_gaps) | comb_gaps as u32)
}

fn cut(in_col: BitColumn, indexer: &[i32]) -> Result<(BitColumn, u32)> {
    let mut cut_mask = 0;
    let mut active = 0;
    while active < indexer.len() && indexer[active] != -1 {
        set(&mut cut_mask, active, test(in_col, indexer[active] as usize))?;
        active += 1;
    }
    Ok((cut_mask, active as u32))
}

fn extract_common_mask(
    column_q: &Column,
    column_j: &Column,
    mask_colj: BitColumn,
    back_indexer_q: &[i32],
    forw_indexer_q: &[i32],
) -> Result<(BitColumn, u32)> {
    let mut mask_qj = 0;
    let mut active = 0;
    while active < back_indexer_q.len() && back_indexer_q[active] != -1 {
        let (fq, bj) = (forw_indexer_q[active] as usize, back_indexer_q[active] as usize);
        if (column_q[fq].allele != column_j[bj].allele) != test(mask_colj, bj) {
            set(&mut mask_qj, active, true)?;
        }
        active += 1;
    }
    Ok((mask_qj, active as u32))
}

fn compute_weight_mask(mask: BitColumn, column: &Column) -> u32 {
    let mut comb = mask;
    let mut weight = 0;
    let mut pos = 0usize;
    while comb != 0 {
        let temp = ffsl(comb);
        pos += temp as usize;
        weight = add(weight, column[pos - 1].phred);
        comb = shr(comb, temp);
    }
    weight
}

fn reconstruct_haplotypes(
    bt1: &[Vec<Vec<Backtrace1>>],
    bt2_hap: &[Vec<Vec<bool>>],
    bt2_nb: &[Vec<Vec<bool>>],
    is_homozygous: &[bool],
    homo_haplotypes: &[bool],
    best1: &[Backtrace1],
    best2_hap: &[bool],
    best2_nb: &[bool],
) -> (Vec<bool>, Vec<bool>) {
    let mut col = bt1.len() - 1;
    let mut h1 = vec![false; col];
    let mut h2 = vec![false; col];
    while col > 0 {
        while is_homozygous[col] {
            let v = !homo_haplotypes[col];
            h1[col - 1] = v;
            h2[col - 1] = v;
            col -= 1;
        }
        let mut back1 = best1[col];
        let mut back2_hap = best2_hap[col];
        let mut back2_nb = best2_nb[col];
        let mut flag = col > 0;
        while flag {
            h1[col - 1] = !back2_hap;
            h2[col - 1] = back2_hap;
            for _ in 0..back1.jump - 1 {
                col -= 1;
                let v = !homo_haplotypes[col];
                h1[col - 1] = v;
                h2[col - 1] = v;
            }
            col -= 1;
            if back2_nb || col == 0 {
                flag = false;
            } else {
                let (j, i) = (back1.jump as usize, back1.index as usize);
                back2_hap = bt2_hap[col][j][i];
                back2_nb = bt2_nb[col][j][i];
                back1 = bt1[col][j][i];
            }
        }
    }
    (h1, h2)
}

/// Phase with `HapChatCore` as `phase --algorithm hapchat` does. The
/// partitioning is not implemented by HapChat (the read list shows `*`).
pub fn solve(read_set: &ReadSet) -> Result<Solution> {
    let sample_id = if read_set.is_empty() { 0 } else { read_set.get(0).sample_id };
    let mut read0 = Read::new("superread_0_0", -1, -1, sample_id);
    let mut read1 = Read::new("superread_1_0", -1, -1, sample_id);
    let mut optimum = 0;
    let blocks = read_blocks(read_set);
    if let Some(first) = blocks.first() {
        let mut block = ReadSet::with_capacity(first.len());
        for &i in first {
            block.add(read_set.get(i).clone())?;
        }
        block.reassign_read_ids();
        block.sort();
        let it = ColumnIterator::new(&block, None)?;
        let mut dp = Dp { bin: binomials(), columns: Columns { it } };
        let column_count = dp.columns.column_count();
        let (mut h1, mut h2) = (vec![false; column_count], vec![false; column_count]);
        if column_count > 0 {
            let (opt, a, b) = dp.run()?;
            optimum = opt;
            h1 = a;
            h2 = b;
        }
        for (i, position) in read_set.get_positions().into_iter().enumerate() {
            let allele = |h: &[bool]| h.get(i).copied().unwrap_or(false) as i32;
            read0.add_variant(position, allele(&h1), 30);
            read1.add_variant(position, allele(&h2), 30);
        }
    }
    let mut superreads = ReadSet::new();
    superreads.add(read0)?;
    superreads.add(read1)?;
    Ok(Solution {
        optimal_score: optimum,
        superreads: vec![superreads],
        transmission_vector: Vec::new(),
        partitioning: Vec::new(),
    })
}

//! Port of `switchflipcalculator.cpp` (polyphase): the minimum number of
//! switches and flips that turn one polyploid phasing into another.
//!
//! Ties between equally good permutations are broken by the iteration order
//! of libstdc++'s `std::unordered_map`, which is emulated so that the split
//! into switches and flips is the same as in WhatsHap.

use crate::compat::cxx::CxxHashMap;

type Score = f64;

/// A permutation of haplotypes: 4 bits per position (`Permutation`).
type Permutation = u64;

const BITS_PER_BLOCK: u32 = 4;
const INVALID: Permutation = 0xf000_0000_0000_0000;

fn get(p: Permutation, i: u32) -> u32 {
    ((p >> (i * BITS_PER_BLOCK)) & 15) as u32
}

#[derive(Clone, Copy)]
struct Entry {
    score: Score,
    pred: Permutation,
}

impl Default for Entry {
    fn default() -> Self {
        Entry { score: f64::INFINITY, pred: INVALID }
    }
}

pub struct SwitchFlipCalculator {
    ploidy: u32,
    switch_cost: f64,
    flip_cost: f64,
}

/// Result of `compare`: total switches and flips (both infinite if no path
/// was found), and per position the switches, the flipped haplotypes and
/// the permutation of the optimal path.
pub struct SwitchFlips {
    pub switches: Score,
    pub flips: Score,
    pub switches_in_column: Vec<u32>,
    pub flipped_haps_in_column: Vec<Vec<u32>>,
    pub perm_in_column: Vec<Vec<u32>>,
}

impl SwitchFlipCalculator {
    pub fn new(ploidy: u32, switch_cost: f64, flip_cost: f64) -> Self {
        SwitchFlipCalculator { ploidy, switch_cost, flip_cost }
    }

    /// `compare`: phasings are given position-wise (one allele per haplotype).
    pub fn compare(&self, phasing0: &[Vec<u32>], phasing1: &[Vec<u32>]) -> SwitchFlips {
        let permutations = self.permutations();
        let num_vars = phasing0.len();
        let mut m: Vec<CxxHashMap<Permutation, Entry>> = Vec::with_capacity(num_vars);

        // initialize first column
        let mut column: CxxHashMap<Permutation, Entry> = CxxHashMap::new();
        for &p in &permutations {
            *column.entry_or_default(p) =
                Entry { score: self.flip_cost * self.num_flips(p, &phasing0[0], &phasing1[0]), pred: INVALID };
        }
        m.push(column.copied());

        for pos in 1..num_vars {
            column.clear();
            let mut min_exists = false;
            let mut minimum_in_column = f64::INFINITY;
            for &row in &permutations {
                let mut minimum = f64::INFINITY;
                let mut minimum_pred = INVALID;
                for (pred, entry) in m[pos - 1].iter() {
                    let s = entry.score + self.switch_cost * self.num_switches(row, pred);
                    if s < minimum {
                        min_exists = true;
                        minimum = s;
                        minimum_pred = pred;
                    }
                }
                let flip_cost = self.flip_cost * self.num_flips(row, &phasing0[pos], &phasing1[pos]);
                let entry = if min_exists {
                    Entry { score: minimum + flip_cost, pred: minimum_pred }
                } else {
                    Entry { score: flip_cost, pred: INVALID }
                };
                *column.entry_or_default(row) = entry;
                if entry.score < minimum_in_column {
                    minimum_in_column = entry.score;
                }
            }

            // remove non-profitable entries
            let score = |column: &CxxHashMap<Permutation, Entry>, p| column.get(p).map_or(f64::INFINITY, |e| e.score);
            let mut profitable = Vec::new();
            let mut open = Vec::new();
            for &t in &permutations {
                if score(&column, t) <= minimum_in_column {
                    profitable.push(t);
                } else {
                    open.push(t);
                }
            }
            for t in open {
                let st = score(&column, t);
                let is_profitable =
                    !profitable.iter().any(|&p| st >= score(&column, p) + self.switch_cost * self.num_switches(t, p));
                if is_profitable {
                    if profitable.len() < self.ploidy as usize {
                        profitable.push(t);
                    }
                } else {
                    column.remove(t);
                }
            }
            m.push(column.copied());
        }

        // backtracking start
        let mut result = SwitchFlips {
            switches: 0.0,
            flips: 0.0,
            switches_in_column: Vec::new(),
            flipped_haps_in_column: Vec::new(),
            perm_in_column: Vec::new(),
        };
        let last = num_vars - 1;
        let mut current = INVALID;
        let mut minimum = f64::INFINITY;
        for (p, entry) in m[last].iter() {
            if entry.score < minimum {
                minimum = entry.score;
                current = p;
            }
        }
        if current == INVALID {
            println!("No minimum at end of compared block!");
            result.switches = f64::INFINITY;
            result.flips = f64::INFINITY;
            return result;
        }
        let pred = m[last].get(current).map_or(INVALID, |e| e.pred);
        self.record(&mut result, current, self.num_switches(current, pred), &phasing0[last], &phasing1[last]);

        // backtracking iteration
        for pos in (0..last).rev() {
            current = m[pos + 1].get(current).map_or(INVALID, |e| e.pred);
            let local_switches = if pos == 0 {
                0.0
            } else {
                self.num_switches(current, m[pos].get(current).map_or(INVALID, |e| e.pred))
            };
            self.record(&mut result, current, local_switches, &phasing0[pos], &phasing1[pos]);
        }

        // reverse as we constructed it back to front
        result.perm_in_column.reverse();
        result.flipped_haps_in_column.reverse();
        result.switches_in_column.reverse();
        result
    }

    fn record(&self, result: &mut SwitchFlips, row: Permutation, switches: Score, phase0: &[u32], phase1: &[u32]) {
        result.perm_in_column.push((0..self.ploidy).map(|i| get(row, i)).collect());
        result
            .flipped_haps_in_column
            .push((0..self.ploidy).filter(|&i| phase0[get(row, i) as usize] != phase1[i as usize]).collect());
        result.switches_in_column.push(switches as u32);
        result.flips += self.num_flips(row, phase0, phase1);
        result.switches += switches;
    }

    fn num_flips(&self, p: Permutation, phase0: &[u32], phase1: &[u32]) -> Score {
        (0..self.ploidy).filter(|&i| phase0[get(p, i) as usize] != phase1[i as usize]).count() as Score
    }

    /// Number of positions at which two permutations differ
    fn num_switches(&self, p1: Permutation, p2: Permutation) -> Score {
        const MASK: u64 = 0x0111_1111_1111_1111;
        let x = p1 ^ p2;
        let y = (x & MASK) | ((x >> 3) & MASK) | ((x >> 2) & MASK) | ((x >> 1) & MASK);
        y.count_ones() as Score
    }

    /// All permutations in lexicographic order (`std::next_permutation`)
    fn permutations(&self) -> Vec<Permutation> {
        let mut cur: Vec<u32> = (0..self.ploidy).collect();
        let mut perms = Vec::new();
        loop {
            perms.push(cur.iter().enumerate().map(|(i, &h)| u64::from(h) << (i as u32 * BITS_PER_BLOCK)).sum());
            // next_permutation
            let Some(i) = (1..cur.len()).rev().find(|&i| cur[i - 1] < cur[i]) else { break };
            let j = (i..cur.len()).rev().find(|&j| cur[j] > cur[i - 1]).expect("a larger element exists");
            cur.swap(i - 1, j);
            cur[i..].reverse();
        }
        perms
    }
}

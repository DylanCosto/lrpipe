//! glibc's `rand()`, which the original uses to pick reads when a candidate is deeper than the
//! matrix.
//!
//! The original never calls `srand`, so every full-alignment job (one process each) starts from
//! glibc's default state, `srandom(1)`, and draws from it in candidate order.

/// `RAND_MAX` in glibc.
pub(crate) const RAND_MAX: u32 = 0x7fff_ffff;

/// Degree of glibc's default `TYPE_3` additive feedback generator.
const DEGREE: usize = 31;
/// Separation between the front and rear pointers (`SEP_3`).
const SEPARATION: usize = 3;

/// glibc `random_r` with a `TYPE_3` state, as `rand()` uses it.
#[derive(Debug, Clone)]
pub(crate) struct GlibcRand {
    state: [i32; DEGREE],
    front: usize,
    rear: usize,
}

impl Default for GlibcRand {
    fn default() -> Self {
        Self::new()
    }
}

impl GlibcRand {
    /// The state of a process that has not called `srand` (`srandom_r(1)`).
    pub(crate) fn new() -> Self {
        Self::with_seed(1)
    }

    /// `srandom_r(seed)`.
    pub(crate) fn with_seed(seed: u32) -> Self {
        let mut state = [0i32; DEGREE];
        // A seed of 0 would give an all-zero state, so glibc uses 1 instead.
        let mut word = if seed == 0 { 1 } else { i64::from(seed as i32) };
        state[0] = word as i32;
        for slot in state.iter_mut().skip(1) {
            // Park-Miller "minimal standard" step via Schrage's method, in `long int`.
            let hi = word / 127_773;
            let lo = word % 127_773;
            word = 16_807 * lo - 2_836 * hi;
            if word < 0 {
                word += 2_147_483_647;
            }
            *slot = word as i32;
        }
        let mut rng = Self { state, front: SEPARATION, rear: 0 };
        for _ in 0..DEGREE * 10 {
            rng.next_u31();
        }
        rng
    }

    /// The next `rand()` value, in `0..=RAND_MAX`.
    pub(crate) fn next_u31(&mut self) -> u32 {
        let val = (self.state[self.front] as u32).wrapping_add(self.state[self.rear] as u32);
        self.state[self.front] = val as i32;
        self.front += 1;
        if self.front >= DEGREE {
            self.front = 0;
            self.rear += 1;
        } else {
            self.rear += 1;
            if self.rear >= DEGREE {
                self.rear = 0;
            }
        }
        val >> 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values from glibc (`tools/khash_oracle_fullalign.c rand N`).
    #[test]
    fn matches_glibc_default_sequence() {
        let mut rng = GlibcRand::new();
        let first: Vec<u32> = (0..5).map(|_| rng.next_u31()).collect();
        assert_eq!(first, [1_804_289_383, 846_930_886, 1_681_692_777, 1_714_636_915, 1_957_747_793]);
        let mut rng = GlibcRand::new();
        let mut nth = |n: usize| (0..n).map(|_| rng.next_u31()).last().expect("n > 0");
        assert_eq!(nth(1000), 1_143_565_421);
        assert_eq!(nth(99_000), 660_658_275);
        assert_eq!(nth(900_000), 429_357_853);
    }

    #[test]
    fn stays_in_range() {
        let mut rng = GlibcRand::with_seed(0);
        assert!((0..10_000).all(|_| rng.next_u31() <= RAND_MAX));
    }
}

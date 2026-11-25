//! Python's `random.seed(0); random.sample(range(n), k)` (CPython/PyPy 3.11: MT19937 and
//! `_randbelow_with_getrandbits`), which `sorted_by_hap_read_name` uses to pick reads.

/// MT19937 as `_random.Random` seeds it from an integer.
struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Mt19937 {
    /// `random_seed` with a non-negative integer: `init_by_array` over its 32-bit words.
    fn from_seed(seed: u64) -> Self {
        let key: Vec<u32> = if seed >> 32 == 0 { vec![seed as u32] } else { vec![seed as u32, (seed >> 32) as u32] };
        let mut mt = Self::init_genrand(19_650_218);
        let n = 624;
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..n.max(key.len()) {
            let prev = mt.state[i - 1] ^ (mt.state[i - 1] >> 30);
            mt.state[i] = (mt.state[i] ^ prev.wrapping_mul(1_664_525)).wrapping_add(key[j]).wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= n {
                mt.state[0] = mt.state[n - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..n - 1 {
            let prev = mt.state[i - 1] ^ (mt.state[i - 1] >> 30);
            mt.state[i] = (mt.state[i] ^ prev.wrapping_mul(1_566_083_941)).wrapping_sub(i as u32);
            i += 1;
            if i >= n {
                mt.state[0] = mt.state[n - 1];
                i = 1;
            }
        }
        mt.state[0] = 0x8000_0000;
        mt.index = n;
        mt
    }

    fn init_genrand(seed: u32) -> Self {
        let mut state = [0u32; 624];
        state[0] = seed;
        for i in 1..624 {
            state[i] = 1_812_433_253u32.wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30)).wrapping_add(i as u32);
        }
        Self { state, index: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        const UPPER: u32 = 0x8000_0000;
        const LOWER: u32 = 0x7fff_ffff;
        const MATRIX_A: u32 = 0x9908_b0df;
        if self.index >= 624 {
            for k in 0..624 {
                let y = (self.state[k] & UPPER) | (self.state[(k + 1) % 624] & LOWER);
                let mag = if y & 1 == 0 { 0 } else { MATRIX_A };
                self.state[k] = self.state[(k + 397) % 624] ^ (y >> 1) ^ mag;
            }
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// `getrandbits(k)` for `0 < k <= 32`.
    fn getrandbits(&mut self, k: u32) -> u64 {
        u64::from(self.next_u32() >> (32 - k))
    }

    /// `_randbelow_with_getrandbits(n)`.
    fn randbelow(&mut self, n: u64) -> u64 {
        let k = 64 - n.leading_zeros();
        assert!(k <= 32, "randbelow({n}) needs more than 32 bits");
        loop {
            let r = self.getrandbits(k);
            if r < n {
                return r;
            }
        }
    }
}

/// `random.seed(0); random.sample(range(n), k)`, in selection order.
///
/// # Panics
///
/// Panics if `k > n` (Python raises `ValueError`).
pub(crate) fn sample_seed0(n: usize, k: usize) -> Vec<usize> {
    assert!(k <= n, "sample of {k} from {n}");
    let mut rng = Mt19937::from_seed(0);
    let mut setsize = 21usize;
    if k > 5 {
        // `4 ** _ceil(_log(k * 3, 4))`.
        let exp = ((k * 3) as f64).ln() / 4f64.ln();
        setsize += 4usize.pow(exp.ceil() as u32);
    }
    if n <= setsize {
        let mut pool: Vec<usize> = (0..n).collect();
        (0..k)
            .map(|i| {
                let j = rng.randbelow((n - i) as u64) as usize;
                let picked = pool[j];
                pool[j] = pool[n - i - 1];
                picked
            })
            .collect()
    } else {
        let mut selected = std::collections::HashSet::new();
        (0..k)
            .map(|_| {
                loop {
                    let j = rng.randbelow(n as u64) as usize;
                    if selected.insert(j) {
                        break j;
                    }
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_python() {
        // `random.seed(0); random.random()` is 0.8444218515250481; its first word is
        // 3626764237 (from `getrandbits(32)` after `random.seed(0)`).
        let mut rng = Mt19937::from_seed(0);
        assert_eq!(rng.next_u32(), 3_626_764_237);
        // `random.seed(0); random.sample(range(10), 3)` is [6, 9, 0].
        assert_eq!(sample_seed0(10, 3), [6, 9, 0]);
        // `random.seed(0); random.sample(range(n), 55)[:10]` with a pool (n = 60) and with a
        // set (n = 300 and 1000).
        assert_eq!(sample_seed0(60, 55)[..10], [54, 24, 48, 56, 26, 2, 16, 32, 31, 25]);
        assert_eq!(sample_seed0(300, 55)[..10], [197, 215, 20, 132, 261, 248, 207, 155, 244, 183]);
        assert_eq!(sample_seed0(1000, 55)[..10], [864, 394, 776, 911, 430, 41, 265, 988, 523, 497]);
    }
}

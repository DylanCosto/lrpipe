//! CPython's `random.Random` (`Lib/random.py` and `Modules/_randommodule.c`, Python 3.10).

use sha2::{Digest, Sha512};

use crate::mt19937::Mt19937;

/// A CPython `random.Random` generator.
#[derive(Clone, Debug)]
pub struct Random {
    mt: Mt19937,
}

impl Random {
    /// `random.Random(); seed(s)` for a `str` seed (version 2): the seed is the big-endian
    /// integer of the UTF-8 bytes followed by their SHA-512 digest.
    pub fn from_str_seed(seed: &str) -> Self {
        let mut bytes = seed.as_bytes().to_vec();
        bytes.extend_from_slice(&Sha512::digest(seed.as_bytes()));
        Random { mt: Mt19937::from_key(&int_key(&bytes)) }
    }

    /// `getrandbits(k)` for `k <= 32`.
    ///
    /// # Panics
    ///
    /// If `k > 32`.
    pub fn getrandbits(&mut self, k: u32) -> u32 {
        assert!(k <= 32, "getrandbits({k}) needs more than one 32-bit word");
        if k == 0 { 0 } else { self.mt.next_u32() >> (32 - k) }
    }

    /// `_randbelow(n)`: uniform in `0..n`, by rejection from `n.bit_length()` random bits.
    ///
    /// # Panics
    ///
    /// If `n` is 0 or needs more than 32 bits.
    pub fn randbelow(&mut self, n: u64) -> u64 {
        assert!(n > 0 && n <= 1 << 32, "randbelow({n}) out of range");
        let k = 64 - n.leading_zeros();
        loop {
            let r = u64::from(self.getrandbits(k));
            if r < n {
                return r;
            }
        }
    }

    /// `choice(seq)`.
    ///
    /// # Panics
    ///
    /// If `seq` is empty.
    pub fn choice<'a, T>(&mut self, seq: &'a [T]) -> &'a T {
        &seq[self.randbelow(seq.len() as u64) as usize]
    }

    /// `sample(range(n), k)`: `k` distinct values of `0..n`, in draw order.
    ///
    /// # Panics
    ///
    /// If `k > n`.
    pub fn sample_range(&mut self, n: u64, k: usize) -> Vec<u64> {
        assert!(k as u64 <= n, "sample larger than population ({k} > {n})");
        let mut setsize = 21usize;
        if k > 5 {
            // `4 ** _ceil(_log(k * 3, 4))`, with math.log(x, base) = log(x) / log(base).
            let e = (((k * 3) as f64).ln() / 4f64.ln()).ceil() as u32;
            setsize += 4usize.pow(e);
        }
        if n <= setsize as u64 {
            let mut pool: Vec<u64> = (0..n).collect();
            (0..k)
                .map(|i| {
                    let j = self.randbelow(n - i as u64) as usize;
                    let v = pool[j];
                    pool[j] = pool[(n - i as u64 - 1) as usize];
                    v
                })
                .collect()
        } else {
            let mut selected = std::collections::HashSet::with_capacity(k);
            (0..k)
                .map(|_| {
                    loop {
                        let j = self.randbelow(n);
                        if selected.insert(j) {
                            break j;
                        }
                    }
                })
                .collect()
        }
    }
}

/// The 32-bit little-endian words of a non-negative big-endian integer (CPython's
/// `random_seed`: `_PyLong_AsByteArray(..., little_endian=1)`), at least one word.
fn int_key(big_endian: &[u8]) -> Vec<u32> {
    let mut le: Vec<u8> = big_endian.iter().rev().copied().collect();
    while le.last() == Some(&0) {
        le.pop();
    }
    le.resize(le.len().div_ceil(4).max(1) * 4, 0);
    le.as_chunks::<4>().0.iter().map(|&c| u32::from_le_bytes(c)).collect()
}

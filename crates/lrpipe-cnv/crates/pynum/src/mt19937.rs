//! The 32-bit Mersenne Twister (Matsumoto and Nishimura's MT19937, 2002 seeding), which both
//! CPython's `random` and numpy's legacy `RandomState` use.

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

#[derive(Clone)]
pub(crate) struct Mt19937 {
    state: [u32; N],
    index: usize,
}

impl std::fmt::Debug for Mt19937 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mt19937").field("index", &self.index).finish_non_exhaustive()
    }
}

impl Mt19937 {
    /// `init_genrand(s)`.
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = [0u32; N];
        state[0] = seed;
        for i in 1..N {
            state[i] = 1_812_433_253u32.wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30)).wrapping_add(i as u32);
        }
        Mt19937 { state, index: N }
    }

    /// `init_by_array(key)`.
    pub(crate) fn from_key(key: &[u32]) -> Self {
        let mut mt = Mt19937::new(19_650_218);
        let s = &mut mt.state;
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..N.max(key.len()) {
            s[i] = (s[i] ^ (s[i - 1] ^ (s[i - 1] >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                s[0] = s[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..N - 1 {
            s[i] = (s[i] ^ (s[i - 1] ^ (s[i - 1] >> 30)).wrapping_mul(1_566_083_941)).wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                s[0] = s[N - 1];
                i = 1;
            }
        }
        s[0] = 0x8000_0000;
        mt
    }

    /// `genrand_uint32()`.
    pub(crate) fn next_u32(&mut self) -> u32 {
        if self.index >= N {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    fn twist(&mut self) {
        let s = &mut self.state;
        for k in 0..N {
            let y = (s[k] & UPPER_MASK) | (s[(k + 1) % N] & LOWER_MASK);
            let mag = if y & 1 == 0 { 0 } else { MATRIX_A };
            s[k] = s[(k + M) % N] ^ (y >> 1) ^ mag;
        }
        self.index = 0;
    }
}

//! numpy's legacy global generator (`np.random.seed`, `np.random.uniform`), numpy 1.24.

use crate::mt19937::Mt19937;

/// `np.random.RandomState` seeded with an integer.
#[derive(Clone, Debug)]
pub struct RandomState {
    mt: Mt19937,
}

impl RandomState {
    /// `np.random.seed(seed)`: `init_genrand(seed)`.
    pub fn new(seed: u32) -> Self {
        RandomState { mt: Mt19937::new(seed) }
    }

    /// One `np.random.uniform()` draw in `[0, 1)`: 53 random bits from two outputs.
    pub fn uniform(&mut self) -> f64 {
        let a = self.mt.next_u32() >> 5;
        let b = self.mt.next_u32() >> 6;
        (f64::from(a) * 67_108_864.0 + f64::from(b)) / 9_007_199_254_740_992.0
    }
}

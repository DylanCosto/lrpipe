//! Elementwise functions the compiler can vectorize: `exp` without a library call, and the LSTM
//! cell update built on it, compiled for AVX2 and FMA when the CPU has them.

/// `exp(x)`, within 2 units in the last place of `f32`'s for `x` in `[-87, 88]`, and `x` clamped
/// to that range (Cephes' `expf`: `2^n * exp(r)` with `|r| <= ln 2 / 2`).
#[inline(always)]
fn exp(x: f32) -> f32 {
    /// Adding 1.5 * 2^23 rounds a float of magnitude below 2^22 to an integer, which then sits
    /// in the low bits of the sum.
    const ROUND: f32 = 12_582_912.0;
    const LN2_HIGH: f32 = 0.693_359_4;
    const LN2_LOW: f32 = -2.121_944_4e-4;
    let x = x.clamp(-87.0, 88.0);
    let shifted = x * std::f32::consts::LOG2_E + ROUND;
    let n = shifted - ROUND;
    let r = x - n * LN2_HIGH - n * LN2_LOW;
    let mut p = 1.987_569_2e-4;
    p = p * r + 1.398_199_9e-3;
    p = p * r + 8.333_452e-3;
    p = p * r + 4.166_579_6e-2;
    p = p * r + 1.666_666_5e-1;
    p = p * r + 0.5;
    p = p * r * r + r + 1.0;
    // `2^n`: `n` is in [-126, 127], so `n + 127` is a normal float's exponent.
    let exponent = shifted.to_bits().wrapping_sub(ROUND.to_bits()).wrapping_add(127) << 23;
    p * f32::from_bits(exponent)
}

#[inline(always)]
fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + exp(-x))
}

/// `tanh(x)` to within about `1e-7`.
#[inline(always)]
fn tanh(x: f32) -> f32 {
    let e = exp(2.0 * x);
    (e - 1.0) / (e + 1.0)
}

/// One step of an LSTM cell for one sequence: from the pre-activation `gates` (input, forget,
/// cell and output gates, `hidden` values each, PyTorch's order), updates `cell` and writes the
/// new hidden state to `hidden_state`.
#[inline(always)]
fn cell_step(gates: &[f32], cell: &mut [f32], hidden_state: &mut [f32]) {
    let hidden = cell.len();
    let (input, rest) = gates.split_at(hidden);
    let (forget, rest) = rest.split_at(hidden);
    let (candidate, output) = rest.split_at(hidden);
    let gates = input.iter().zip(forget).zip(candidate).zip(&output[..hidden]);
    for ((((i, f), g), o), (c, h)) in gates.zip(cell.iter_mut().zip(hidden_state.iter_mut())) {
        *c = sigmoid(*f) * *c + sigmoid(*i) * tanh(*g);
        *h = sigmoid(*o) * tanh(*c);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
fn cell_step_avx2(gates: &[f32], cell: &mut [f32], hidden_state: &mut [f32]) {
    cell_step(gates, cell, hidden_state);
}

/// [`cell_step`], vectorized for the CPU it runs on.
///
/// # Panics
///
/// If `hidden_state` is not as long as `cell`, or `gates` not four times as long.
pub(super) fn lstm_cell_step(gates: &[f32], cell: &mut [f32], hidden_state: &mut [f32]) {
    assert!(gates.len() == 4 * cell.len() && hidden_state.len() == cell.len(), "LSTM cell shapes");
    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma") {
        // SAFETY: the CPU has AVX2 and FMA.
        unsafe { cell_step_avx2(gates, cell, hidden_state) };
        return;
    }
    cell_step(gates, cell, hidden_state);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exp_is_within_two_ulp() {
        for i in -8700..=8800 {
            let x = i as f32 / 100.0;
            let (got, want) = (exp(x), f64::from(x).exp());
            let nearest = want as f32;
            let ulp = f64::from(f32::from_bits(nearest.to_bits() + 1) - nearest);
            assert!((f64::from(got) - want).abs() <= 2.0 * ulp, "exp({x}) = {got}, want {want}");
        }
        assert_eq!(exp(0.0), 1.0);
        assert!(exp(-1000.0) > 0.0 && exp(1000.0).is_finite());
    }

    #[test]
    fn sigmoid_and_tanh_are_close() {
        for i in -2000..=2000 {
            let x = i as f32 / 100.0;
            let sigmoid_x = 1.0 / (1.0 + (-f64::from(x)).exp());
            assert!((f64::from(sigmoid(x)) - sigmoid_x).abs() < 2e-7, "sigmoid({x})");
            assert!((f64::from(tanh(x)) - f64::from(x).tanh()).abs() < 2e-7, "tanh({x})");
        }
    }

    #[test]
    fn the_cell_step_follows_pytorch() {
        let hidden = 19;
        let gates: Vec<f32> = (0..4 * hidden).map(|i| (i as f32 * 0.7).sin() * 3.0).collect();
        let mut cell: Vec<f32> = (0..hidden).map(|i| (i as f32 * 0.3).cos()).collect();
        let mut want_cell = cell.clone();
        let mut state = vec![0.0; hidden];
        lstm_cell_step(&gates, &mut cell, &mut state);
        let s = |x: f32| 1.0 / (1.0 + (-x).exp());
        for j in 0..hidden {
            let g = |k: usize| gates[k * hidden + j];
            want_cell[j] = s(g(1)) * want_cell[j] + s(g(0)) * g(2).tanh();
            assert!((cell[j] - want_cell[j]).abs() < 1e-6, "cell {j}");
            assert!((state[j] - s(g(3)) * want_cell[j].tanh()).abs() < 1e-6, "hidden {j}");
        }
    }
}

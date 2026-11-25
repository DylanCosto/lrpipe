//! The Rust networks against libtorch's probabilities on the same synthetic weights, which
//! libtorch 2.13.0 (one CPU thread) wrote to `tests/data/libtorch_*.f32` before lrpipe-snv
//! dropped it (see `tests/data/README.md`).

use std::path::Path;

use super::synthetic::{full_alignment_state, largest_difference, pileup_state};
use super::{CpuFullAlignmentModel, CpuPileupModel, Weights};
use crate::{
    FULL_ALIGNMENT_CHANNELS, FULL_ALIGNMENT_OUTPUTS, FULL_ALIGNMENT_TENSOR_LEN, FullAlignmentModel, PILEUP_OUTPUTS,
    PILEUP_TENSOR_LEN, PileupModel, Precision,
};

/// How far the Rust networks may be from libtorch's: the same `float32` arithmetic, summed in
/// another order.
const LIBTORCH: f32 = 1e-4;

/// libtorch's probabilities in `tests/data/name`.
fn libtorch(name: &str) -> Vec<f32> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name);
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    bytes.as_chunks::<4>().0.iter().map(|&b| f32::from_le_bytes(b)).collect()
}

#[test]
fn pileup_matches_libtorch() {
    let mut rust = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), pileup_state())).expect("rust");
    // Counts like the pileup tensors': small integers, some negative, over three sub-batches.
    let n = 120;
    let inputs: Vec<i32> = (0..n * PILEUP_TENSOR_LEN).map(|i| ((i * 37) % 61) as i32 - 20).collect();
    let got = rust.predict(&inputs).expect("rust");
    assert_eq!(got.len(), n * PILEUP_OUTPUTS);
    let difference = largest_difference(&got, &libtorch("libtorch_pileup.f32"));
    assert!(difference < LIBTORCH, "largest difference {difference}");
    assert!(rust.predict(&[]).expect("empty").is_empty());
    assert!(rust.predict(&[0; PILEUP_TENSOR_LEN + 1]).is_err());
}

#[test]
fn eight_bit_pileup_stays_close() {
    if !crate::int8::supported() {
        return;
    }
    let state = || pileup_state();
    let mut exact = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), state())).expect("rust");
    let mut mixed = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), state()))
        .expect("rust")
        .with_precision(Precision::Mixed);
    let n = 50;
    let inputs: Vec<i32> = (0..n * PILEUP_TENSOR_LEN).map(|i| ((i * 37) % 61) as i32 - 20).collect();
    let (want, got) = (exact.predict(&inputs).expect("rust"), mixed.predict(&inputs).expect("rust"));
    assert_eq!(got.len(), want.len());
    let difference = largest_difference(&got, &want);
    assert!(difference < 0.06, "largest difference in mixed precision: {difference}");
}

#[test]
fn full_alignment_matches_libtorch() {
    let state = full_alignment_state();
    let mut rust = CpuFullAlignmentModel::from_state(&mut Weights::new(Path::new("synthetic"), state)).expect("rust");
    let n = 3;
    let inputs: Vec<i8> = (0..n * FULL_ALIGNMENT_TENSOR_LEN).map(|i| ((i * 37) % 201) as i8).collect();
    let got = rust.predict(&inputs).expect("rust");
    assert_eq!(got.len(), n * FULL_ALIGNMENT_OUTPUTS);
    let difference = largest_difference(&got, &libtorch("libtorch_full_alignment.f32"));
    assert!(difference < LIBTORCH, "largest difference {difference}");
    assert_eq!(rust.channels(), FULL_ALIGNMENT_CHANNELS);

    // HiFi's shallower matrices.
    let mut shallow = rust.with_matrix_depth(55);
    let shallow_inputs = &inputs[..2 * 55 * 33 * FULL_ALIGNMENT_CHANNELS];
    let got = shallow.predict(shallow_inputs).expect("rust");
    let difference = largest_difference(&got, &libtorch("libtorch_full_alignment_55.f32"));
    assert!(difference < LIBTORCH, "largest difference at depth 55: {difference}");

    // 8-bit convolutions, with close calls scored again in float32: near the exact network.
    let mut exact = shallow.with_matrix_depth(89);
    let inputs: Vec<i8> = (0..n * FULL_ALIGNMENT_TENSOR_LEN).map(|i| ((i * 13) % 101) as i8).collect();
    let want = exact.predict(&inputs).expect("rust");
    let mut mixed = exact.with_precision(Precision::Mixed);
    let difference = largest_difference(&mixed.predict(&inputs).expect("rust"), &want);
    assert!(difference < 1e-2, "largest difference in mixed precision: {difference}");
}

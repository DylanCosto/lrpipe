//! The Metal networks against the CPU ones on the same synthetic weights. They run only where a
//! Metal GPU is found, and say so otherwise. Both are `float32`, summing in other orders, so the
//! GPU is held to 1e-5 of the CPU (CUDA's TF32 operands get 1e-2); it has come within 6e-7.

use std::path::Path;

use super::{MetalFullAlignmentModel, MetalPileupModel};
use crate::cpu::synthetic::{full_alignment_state, full_alignment_state_with, largest_difference, pileup_state};
use crate::cpu::{CpuFullAlignmentModel, CpuPileupModel, Weights};
use crate::{
    Error, FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME, FULL_ALIGNMENT_DEPTH, FULL_ALIGNMENT_OUTPUTS,
    FULL_ALIGNMENT_TENSOR_LEN, FullAlignmentModel, PILEUP_OUTPUTS, PILEUP_TENSOR_LEN, PileupModel,
};

/// How far the GPU's `float32` may be from the CPU's.
const FLOAT32: f32 = 1e-5;

/// The model, or `None` (with a note) where there is no Metal GPU. Kernels that do not compile
/// fail the test, and so does a missing GPU with `SNV_REQUIRE_GPU` set.
fn on_gpu<M>(model: Result<M, Error>) -> Option<M> {
    match model {
        Ok(model) => Some(model),
        Err(ref err @ Error::Backend { ref message }) if message.contains("compiling the Metal kernels") => {
            panic!("loading: {err}")
        }
        Err(err @ (Error::Backend { .. } | Error::Device { .. })) if std::env::var_os("SNV_REQUIRE_GPU").is_none() => {
            eprintln!("skipped: no usable Metal GPU ({err})");
            None
        }
        Err(err) => panic!("loading: {err}"),
    }
}

#[test]
fn pileup_matches_the_cpu() {
    let cpu = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), pileup_state())).expect("cpu");
    let Some(mut gpu) = on_gpu(MetalPileupModel::from_cpu(&cpu)) else { return };
    let mut cpu = cpu;
    // More candidates than one 64-row tile, and not a multiple of it.
    let n = 150;
    let inputs: Vec<i32> = (0..n * PILEUP_TENSOR_LEN).map(|i| ((i * 37) % 61) as i32 - 20).collect();
    let want = cpu.predict(&inputs).expect("cpu");
    let got = gpu.predict(&inputs).expect("gpu");
    assert_eq!(got.len(), n * PILEUP_OUTPUTS);
    let difference = largest_difference(&got, &want);
    assert!(difference < FLOAT32, "largest difference {difference}");
    // A smaller batch after a larger one reuses the buffers.
    let few = &inputs[..3 * PILEUP_TENSOR_LEN];
    let difference = largest_difference(&gpu.predict(few).expect("gpu"), &want[..3 * PILEUP_OUTPUTS]);
    assert!(difference < FLOAT32, "smaller batch: largest difference {difference}");
    assert!(gpu.predict(&[]).expect("empty").is_empty());
    assert!(gpu.predict(&[0; PILEUP_TENSOR_LEN + 1]).is_err());
}

#[test]
fn full_alignment_matches_the_cpu() {
    let state = full_alignment_state();
    let cpu = CpuFullAlignmentModel::from_state(&mut Weights::new(Path::new("synthetic"), state)).expect("cpu");
    let Some(gpu) = on_gpu(MetalFullAlignmentModel::from_cpu(&cpu)) else { return };
    let mut cpu = cpu;
    let n = 5;
    let inputs: Vec<i8> = (0..n * FULL_ALIGNMENT_TENSOR_LEN).map(|i| ((i * 37) % 201) as i8).collect();
    let want = cpu.predict(&inputs).expect("cpu");
    let mut gpu = gpu;
    let got = gpu.predict(&inputs).expect("gpu");
    assert_eq!(got.len(), n * FULL_ALIGNMENT_OUTPUTS);
    let difference = largest_difference(&got, &want);
    assert!(difference < FLOAT32, "largest difference {difference}");
    assert_eq!(gpu.channels(), cpu.channels());

    // HiFi's shallower matrices.
    let (mut cpu, mut shallow) = (cpu.with_matrix_depth(55), gpu.with_matrix_depth(55));
    let inputs = &inputs[..2 * 55 * 33 * cpu.channels()];
    let difference = largest_difference(&shallow.predict(inputs).expect("gpu"), &cpu.predict(inputs).expect("cpu"));
    assert!(difference < FLOAT32, "at depth 55: largest difference {difference}");
}

#[test]
fn dwell_time_inputs_match_the_cpu() {
    // Nine channels: the convolution's channel index wraps mid-tile.
    let channels = FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME;
    let state = full_alignment_state_with(channels);
    let cpu = CpuFullAlignmentModel::from_state(&mut Weights::new(Path::new("synthetic"), state)).expect("cpu");
    let Some(mut gpu) = on_gpu(MetalFullAlignmentModel::from_cpu(&cpu)) else { return };
    let mut cpu = cpu;
    let n = 3;
    let inputs: Vec<i8> = (0..n * crate::full_alignment_tensor_len(FULL_ALIGNMENT_DEPTH, channels))
        .map(|i| ((i * 29) % 151) as i8)
        .collect();
    let difference = largest_difference(&gpu.predict(&inputs).expect("gpu"), &cpu.predict(&inputs).expect("cpu"));
    assert!(difference < FLOAT32, "largest difference {difference}");
}

/// Split-K sums (a long product over few rows) give what one threadgroup's sum gives.
#[test]
fn split_sums_match_the_cpu() {
    let cpu = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), pileup_state())).expect("cpu");
    let Some(mut gpu) = on_gpu(MetalPileupModel::from_cpu(&cpu)) else { return };
    let mut cpu = cpu;
    // One candidate: L4's 33 * 2 * hidden inputs over a single row tile, so the sum is split.
    let inputs: Vec<i32> = (0..PILEUP_TENSOR_LEN).map(|i| ((i * 13) % 41) as i32 - 10).collect();
    let difference = largest_difference(&gpu.predict(&inputs).expect("gpu"), &cpu.predict(&inputs).expect("cpu"));
    assert!(difference < FLOAT32, "largest difference {difference}");
}

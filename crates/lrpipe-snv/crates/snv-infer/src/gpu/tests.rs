//! The GPU networks against the CPU ones on the same synthetic weights. They run only where a
//! CUDA driver, a device and NVRTC are found, and say so otherwise. The GPU's float32 rounds some
//! operands to TF32 as libtorch's does, so it is held to about 1e-3 of the CPU's, not 1e-5.

use std::path::Path;

use super::{GpuFullAlignmentModel, GpuPileupModel};
use crate::cpu::synthetic::{full_alignment_state, full_alignment_state_with, largest_difference, pileup_state};
use crate::cpu::{CpuFullAlignmentModel, CpuPileupModel, Weights};
use crate::{
    Error, FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME, FULL_ALIGNMENT_DEPTH, FULL_ALIGNMENT_OUTPUTS,
    FULL_ALIGNMENT_TENSOR_LEN, FullAlignmentModel, PILEUP_OUTPUTS, PILEUP_TENSOR_LEN, PileupModel,
};

/// How far the GPU's float32 (TF32 operands where libtorch has them) may be from the CPU's.
const FLOAT32: f32 = 0.01;

/// The model, or `None` (with a note) where there is no usable GPU. Kernels that do not compile
/// fail the test, and so does a missing GPU with `SNV_REQUIRE_GPU` set (as `ci/check.sh` sets
/// it where it finds one).
fn on_gpu<M>(model: Result<M, Error>) -> Option<M> {
    match model {
        Ok(model) => Some(model),
        Err(ref err @ Error::Backend { ref message }) if message.contains("compiling the GPU kernels") => {
            panic!("loading: {err}")
        }
        Err(err @ (Error::Backend { .. } | Error::Device { .. })) if std::env::var_os("SNV_REQUIRE_GPU").is_none() => {
            eprintln!("skipped: no usable GPU ({err})");
            None
        }
        Err(err) => panic!("loading: {err}"),
    }
}

#[test]
fn pileup_matches_the_cpu() {
    let cpu = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), pileup_state())).expect("cpu");
    let Some(mut gpu) = on_gpu(GpuPileupModel::from_cpu(&cpu, 0)) else { return };
    let mut cpu = cpu;
    // More candidates than one LSTM block's 64 rows, and not a multiple of it.
    let n = 150;
    let inputs: Vec<i32> = (0..n * PILEUP_TENSOR_LEN).map(|i| ((i * 37) % 61) as i32 - 20).collect();
    let want = cpu.predict(&inputs).expect("cpu");
    let got = gpu.predict(&inputs).expect("gpu");
    assert_eq!(got.len(), n * PILEUP_OUTPUTS);
    let difference = largest_difference(&got, &want);
    assert!(difference < FLOAT32, "largest difference {difference}");
    assert!(gpu.predict(&[]).expect("empty").is_empty());
    assert!(gpu.predict(&[0; PILEUP_TENSOR_LEN + 1]).is_err());
}

#[test]
fn full_alignment_matches_the_cpu() {
    let state = full_alignment_state();
    let cpu = CpuFullAlignmentModel::from_state(&mut Weights::new(Path::new("synthetic"), state)).expect("cpu");
    let Some(gpu) = on_gpu(GpuFullAlignmentModel::from_cpu(&cpu, 0)) else { return };
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
    let Some(mut gpu) = on_gpu(GpuFullAlignmentModel::from_cpu(&cpu, 0)) else { return };
    let mut cpu = cpu;
    let n = 3;
    let inputs: Vec<i8> = (0..n * crate::full_alignment_tensor_len(FULL_ALIGNMENT_DEPTH, channels))
        .map(|i| ((i * 29) % 151) as i8)
        .collect();
    let difference = largest_difference(&gpu.predict(&inputs).expect("gpu"), &cpu.predict(&inputs).expect("cpu"));
    assert!(difference < FLOAT32, "largest difference {difference}");
}

/// On older GPUs (here: the kernels compiled for Pascal 6.0 and Turing 7.5, which this GPU's
/// driver runs too), the kernels round to TF32 in integer arithmetic. The compiler's code for
/// each architecture moves the last bits (compiled for 8.0, with the hardware rounding, it does
/// too), so the probabilities are held to 1e-4 of the native ones.
#[test]
fn older_gpus_give_the_same_probabilities() {
    const ARCHITECTURES: f32 = 1e-4;
    let pileup = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), pileup_state())).expect("cpu");
    let state = full_alignment_state();
    let full = CpuFullAlignmentModel::from_state(&mut Weights::new(Path::new("synthetic"), state)).expect("cpu");
    let Some(mut native) = on_gpu(GpuPileupModel::from_cpu(&pileup, 0)) else { return };
    let mut native_full = GpuFullAlignmentModel::from_cpu(&full, 0).expect("the same GPU");
    let pileup_inputs: Vec<i32> = (0..70 * PILEUP_TENSOR_LEN).map(|i| ((i * 37) % 61) as i32 - 20).collect();
    let full_inputs: Vec<i8> = (0..3 * FULL_ALIGNMENT_TENSOR_LEN).map(|i| ((i * 37) % 201) as i8).collect();
    let want_pileup = native.predict(&pileup_inputs).expect("gpu");
    let want_full = native_full.predict(&full_inputs).expect("gpu");
    for capability in [(6, 0), (7, 5)] {
        let mut older = GpuPileupModel::compiled_for(&pileup, 0, capability).expect("compiled for an older GPU");
        let difference = largest_difference(&older.predict(&pileup_inputs).expect("gpu"), &want_pileup);
        assert!(difference < ARCHITECTURES, "pileup on {capability:?}: largest difference {difference}");

        let mut older = GpuFullAlignmentModel::compiled_for(&full, 0, capability).expect("compiled for an older GPU");
        let difference = largest_difference(&older.predict(&full_inputs).expect("gpu"), &want_full);
        assert!(difference < ARCHITECTURES, "full alignment on {capability:?}: largest difference {difference}");
    }
    let too_old = GpuPileupModel::compiled_for(&pileup, 0, (5, 2)).expect_err("Maxwell is not supported");
    assert!(too_old.to_string().contains("6.0"), "{too_old}");
}

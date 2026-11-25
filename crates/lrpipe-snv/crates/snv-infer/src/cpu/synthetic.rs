//! Synthetic checkpoints with the shapes of the original's, for tests of the backends.

use crate::FULL_ALIGNMENT_CHANNELS;
use crate::checkpoint::{CheckpointTensor, TensorValues};

/// A tensor of `shape` holding `base + scale * sin(i + phase)` over its elements `i`.
fn wave(shape: &[usize], base: f32, scale: f32, phase: f32) -> CheckpointTensor {
    let numel = shape.iter().product::<usize>();
    let values = (0..numel).map(|i| base + scale * (i as f32 + phase).sin()).collect();
    CheckpointTensor { shape: shape.to_vec(), values: TensorValues::F32(values) }
}

/// An `nn.Linear` of `inputs -> outputs` as `{name}.weight` and `{name}.bias`.
fn linear(state: &mut Vec<(String, CheckpointTensor)>, name: &str, outputs: usize, inputs: usize) {
    let scale = 2.0 / (inputs as f32).sqrt();
    state.push((format!("{name}.weight"), wave(&[outputs, inputs], 0.0, scale, 6.0)));
    state.push((format!("{name}.bias"), wave(&[outputs], 0.0, 0.1, 7.0)));
}

/// Deterministic weights with the shapes of the original's pileup checkpoint.
pub(crate) fn pileup_state() -> Vec<(String, CheckpointTensor)> {
    let mut state = Vec::new();
    for (name, inputs, hidden) in [("LSTM1", 18, 128), ("LSTM2", 256, 160)] {
        let scale = 1.0 / (hidden as f32).sqrt();
        for (d, suffix) in ["", "_reverse"].iter().enumerate() {
            let phase = d as f32;
            state.push((format!("{name}.weight_ih_l0{suffix}"), wave(&[4 * hidden, inputs], 0.0, scale, phase)));
            state.push((format!("{name}.weight_hh_l0{suffix}"), wave(&[4 * hidden, hidden], 0.0, scale, phase + 1.0)));
            state.push((format!("{name}.bias_ih_l0{suffix}"), wave(&[4 * hidden], 0.0, 0.1, phase + 2.0)));
            state.push((format!("{name}.bias_hh_l0{suffix}"), wave(&[4 * hidden], 0.0, 0.1, phase + 3.0)));
        }
    }
    for (name, outputs, inputs) in [
        ("L4", 128, 33 * 320),
        ("L5_1", 128, 128),
        ("L5_2", 128, 128),
        ("Y_gt21_logits", 21, 128),
        ("Y_genotype_logits", 3, 128),
    ] {
        linear(&mut state, name, outputs, inputs);
    }
    state
}

/// Deterministic weights with the shapes of the original's full-alignment checkpoint.
pub(crate) fn full_alignment_state() -> Vec<(String, CheckpointTensor)> {
    full_alignment_state_with(FULL_ALIGNMENT_CHANNELS)
}

/// As [`full_alignment_state`], for inputs of `channels` channels (9 with dwell time).
pub(crate) fn full_alignment_state_with(channels: usize) -> Vec<(String, CheckpointTensor)> {
    let mut state = Vec::new();
    let mut conv_bn = |conv: String, bn: String, outputs: usize, inputs: usize| {
        let scale = 2.0 / ((inputs * 9) as f32).sqrt();
        state.push((format!("{conv}.weight"), wave(&[outputs, inputs, 3, 3], 0.0, scale, 0.0)));
        state.push((format!("{conv}.bias"), wave(&[outputs], 0.0, 0.1, 1.0)));
        state.push((format!("{bn}.weight"), wave(&[outputs], 1.0, 0.2, 2.0)));
        state.push((format!("{bn}.bias"), wave(&[outputs], 0.0, 0.1, 3.0)));
        state.push((format!("{bn}.running_mean"), wave(&[outputs], 0.0, 0.1, 4.0)));
        state.push((format!("{bn}.running_var"), wave(&[outputs], 1.0, 0.5, 5.0)));
    };
    for (conv, block, inputs, outputs) in [
        ("conv1", "res_block1.0", channels, 64),
        ("conv3", "res_block2.0", 64, 128),
        ("conv5", "res_block3.0", 128, 256),
    ] {
        conv_bn(format!("{conv}.conv"), format!("{conv}.bn"), outputs, inputs);
        for i in 1..=2 {
            conv_bn(format!("{block}.conv{i}"), format!("{block}.bn{i}"), outputs, outputs);
        }
    }
    for (name, outputs, inputs) in [
        ("L4", 256, 14 * 256),
        ("L5_1", 128, 256),
        ("L5_2", 128, 256),
        ("L5_3", 128, 256),
        ("L5_4", 128, 256),
        ("Y_gt21_logits", 21, 128),
        ("Y_genotype_logits", 3, 128),
        ("Y_indel_length_logits_1", 33, 128),
        ("Y_indel_length_logits_2", 33, 128),
    ] {
        linear(&mut state, name, outputs, inputs);
    }
    state
}

/// The largest difference between two probability vectors.
pub(crate) fn largest_difference(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).fold(0f32, |m, (x, y)| m.max((x - y).abs()))
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;
    use crate::cpu::{CpuFullAlignmentModel, CpuPileupModel, Weights};
    use crate::{
        FULL_ALIGNMENT_OUTPUTS, FULL_ALIGNMENT_TENSOR_LEN, FullAlignmentModel, PILEUP_OUTPUTS, PILEUP_TENSOR_LEN,
        PileupModel,
    };

    /// Each head of each row is a probability distribution.
    fn assert_distributions(probabilities: &[f32], heads: &[usize]) {
        let width: usize = heads.iter().sum();
        for (candidate, row) in probabilities.chunks_exact(width).enumerate() {
            let mut rest = row;
            for &classes in heads {
                let (head, tail) = rest.split_at(classes);
                let sum: f32 = head.iter().sum();
                assert!((sum - 1.0).abs() < 1e-5, "candidate {candidate}: head of {classes} sums to {sum}");
                assert!(head.iter().all(|p| (0.0..=1.0).contains(p)), "candidate {candidate}: {head:?}");
                rest = tail;
            }
        }
    }

    #[test]
    fn the_networks_give_deterministic_distributions() {
        let mut pileup = CpuPileupModel::from_state(&mut Weights::new(Path::new("synthetic"), pileup_state()));
        let pileup = pileup.as_mut().expect("pileup weights cover the network");
        let inputs: Vec<i32> = (0..3 * PILEUP_TENSOR_LEN).map(|i| ((i * 37) % 61) as i32 - 20).collect();
        let first = pileup.predict(&inputs).expect("predict");
        assert_eq!(first.len(), 3 * PILEUP_OUTPUTS);
        assert_distributions(&first, &[21, 3]);
        assert_eq!(largest_difference(&first, &pileup.predict(&inputs).expect("predict")), 0.0);

        let state = full_alignment_state();
        let mut full = CpuFullAlignmentModel::from_state(&mut Weights::new(Path::new("synthetic"), state));
        let full = full.as_mut().expect("full-alignment weights cover the network");
        let inputs: Vec<i8> = (0..2 * FULL_ALIGNMENT_TENSOR_LEN).map(|i| ((i * 37) % 201) as i8).collect();
        let first = full.predict(&inputs).expect("predict");
        assert_eq!(first.len(), 2 * FULL_ALIGNMENT_OUTPUTS);
        assert_distributions(&first, &[21, 3, 33, 33]);
        assert_ne!(first[..FULL_ALIGNMENT_OUTPUTS], first[FULL_ALIGNMENT_OUTPUTS..]);
        assert_eq!(largest_difference(&first, &full.predict(&inputs).expect("predict")), 0.0);
    }
}

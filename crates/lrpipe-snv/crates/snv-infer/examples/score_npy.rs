//! Scores saved model inputs on a device and writes the probabilities, for comparing them with
//! the original's (`tools/score_npy_torch.py`).
//!
//! usage: score_npy pileup|full_alignment MODEL.pt INPUTS.npy BATCH DEVICE OUT.npy
//!
//! `SNV_PRECISION=mixed` runs at [`Precision::Mixed`](snv_infer::Precision::Mixed).

use std::path::Path;
use std::time::Instant;

use ndarray::{Array2, ArrayD};
use ndarray_npy::{read_npy, write_npy};
use snv_infer::{InferenceDevice, ModelConfig, Precision};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let [_, kind, model, inputs, batch, device, out] = &args[..] else {
        return Err("usage: score_npy pileup|full_alignment MODEL INPUTS BATCH DEVICE OUT".into());
    };
    let batch: usize = batch.parse()?;
    let device: InferenceDevice = device.parse()?;
    let precision = match std::env::var("SNV_PRECISION").as_deref() {
        Ok("mixed") => Precision::Mixed,
        _ => Precision::Float32,
    };
    let config = ModelConfig { device, precision };
    let model = Path::new(model);
    let (n, probs, elapsed) = if kind == "pileup" {
        let x: ArrayD<i32> = read_npy(inputs)?;
        let n = x.shape()[0];
        let x = x.as_standard_layout();
        let flat = x.as_slice().ok_or("inputs not contiguous")?;
        let mut m = snv_infer::load_pileup_model(model, config)?;
        let chunk = batch * flat.len() / n.max(1);
        // Warm-up, untimed: the GPU kernels load on the first batch.
        m.predict(&flat[..chunk.min(flat.len())])?;
        let started = Instant::now();
        let mut p = Vec::new();
        for b in flat.chunks(chunk) {
            p.extend(m.predict(b)?);
        }
        (n, p, started.elapsed())
    } else {
        let x: ArrayD<i8> = read_npy(inputs)?;
        let n = x.shape()[0];
        let x = x.as_standard_layout();
        let flat = x.as_slice().ok_or("inputs not contiguous")?;
        let mut m = snv_infer::load_full_alignment_model(model, config, snv_infer::FULL_ALIGNMENT_DEPTH)?;
        let chunk = batch * flat.len() / n.max(1);
        // Warm-up, untimed: the GPU kernels load on the first batch.
        m.predict(&flat[..chunk.min(flat.len())])?;
        let started = Instant::now();
        let mut p = Vec::new();
        for b in flat.chunks(chunk) {
            p.extend(m.predict(b)?);
        }
        (n, p, started.elapsed())
    };
    eprintln!("{kind} on {device:?}: {n} candidates in {elapsed:.2?}");
    let cols = probs.len() / n.max(1);
    write_npy(out, &Array2::from_shape_vec((n, cols), probs)?)?;
    Ok(())
}

//! The bi-LSTM pileup network on the GPU.

use std::path::Path;

use cudarc::driver::{CudaSlice, PushKernelArg};

use super::{Activation, Address, Buffer, DeviceDense, Direction, Gpu, Rounding, Staging, backend};
use crate::cpu::CpuPileupModel;
use crate::cpu::pileup::BiLstm;
use crate::{Error, PILEUP_OUTPUTS, PILEUP_POSITIONS, PILEUP_TENSOR_LEN, PileupModel};

/// A bidirectional LSTM layer on the device.
#[derive(Debug)]
struct DeviceLstm {
    /// Both directions' `W_ih` side by side, with both biases: `[inputs][2 * 4 * hidden]`.
    input: DeviceDense,
    /// `W_hh` transposed, per direction: `[2][hidden][4 * hidden]`.
    hidden_weight: CudaSlice<f32>,
    hidden: usize,
}

impl DeviceLstm {
    fn new(gpu: &Gpu, lstm: &BiLstm) -> Result<Self, Error> {
        let hidden_weight: Vec<f32> = lstm.hidden_weight.concat();
        let mut input = DeviceDense::from_parts(gpu, &lstm.input_weight, &lstm.bias)?;
        // cuDNN computes the input projections with TF32 tensor cores, the recurrence in float32.
        input.rounding = Rounding::Tf32;
        Ok(Self { input, hidden_weight: gpu.upload(&hidden_weight)?, hidden: lstm.hidden })
    }
}

/// Activations, reused across calls.
#[derive(Debug)]
struct Buffers {
    upload: Staging<i32>,
    download: Staging<f32>,
    counts: Buffer<i32>,
    input: Buffer<f32>,
    /// Input projections of the layer being run: `[n][steps][2 * 4 * hidden]`.
    pre: Buffer<f32>,
    cell: Buffer<f32>,
    lstm1: Buffer<f32>,
    lstm2: Buffer<f32>,
    l4: Buffer<f32>,
    l5: Buffer<f32>,
    logits: Buffer<f32>,
    probabilities: Buffer<f32>,
}

/// The original's bi-LSTM pileup network (its `*_P` class), in Rust kernels on a CUDA GPU.
#[derive(Debug)]
pub struct GpuPileupModel {
    gpu: Gpu,
    lstm: [DeviceLstm; 2],
    l4: DeviceDense,
    /// `L5_1` and `L5_2` side by side.
    l5: DeviceDense,
    gt21: DeviceDense,
    genotype: DeviceDense,
    buffers: Buffers,
}

impl GpuPileupModel {
    /// Loads a pileup checkpoint (`pileup.pt`) onto CUDA device `device`.
    ///
    /// # Errors
    ///
    /// As [`CpuPileupModel::load`], and fails if the CUDA driver or NVRTC is missing, the device
    /// does not exist, or the kernels cannot be compiled or the weights copied.
    pub fn load(path: impl AsRef<Path>, device: usize) -> Result<Self, Error> {
        Self::from_cpu(&CpuPileupModel::load(path)?, device)
    }

    /// The network of `cpu` on CUDA device `device`.
    pub(crate) fn from_cpu(cpu: &CpuPileupModel, device: usize) -> Result<Self, Error> {
        Self::on(cpu, Gpu::new(device)?)
    }

    /// The network of `cpu` with its kernels compiled for compute capability `capability`, older
    /// than device `device`'s, as on an older GPU.
    #[cfg(test)]
    pub(crate) fn compiled_for(cpu: &CpuPileupModel, device: usize, capability: (i32, i32)) -> Result<Self, Error> {
        Self::on(cpu, Gpu::compiled_for(device, capability)?)
    }

    /// The network of `cpu` on `gpu`.
    fn on(cpu: &CpuPileupModel, gpu: Gpu) -> Result<Self, Error> {
        Ok(Self {
            lstm: [DeviceLstm::new(&gpu, &cpu.lstm1)?, DeviceLstm::new(&gpu, &cpu.lstm2)?],
            l4: DeviceDense::new(&gpu, &cpu.l4)?,
            l5: DeviceDense::side_by_side(&gpu, &cpu.gt21.hidden, &cpu.genotype.hidden)?,
            gt21: DeviceDense::new(&gpu, &cpu.gt21.logits)?,
            genotype: DeviceDense::new(&gpu, &cpu.genotype.logits)?,
            gpu,
            buffers: Buffers {
                upload: Staging::new(Direction::ToDevice),
                download: Staging::new(Direction::ToHost),
                counts: Buffer::default(),
                input: Buffer::default(),
                pre: Buffer::default(),
                cell: Buffer::default(),
                lstm1: Buffer::default(),
                lstm2: Buffer::default(),
                l4: Buffer::default(),
                l5: Buffer::default(),
                logits: Buffer::default(),
                probabilities: Buffer::default(),
            },
        })
    }

    /// One layer over `x` (`[n][steps][inputs]`) into `out` (`[n][steps][2 * hidden]`).
    fn lstm(&mut self, layer: usize, x: u64, n: usize, out: u64) -> Result<(), Error> {
        let (gpu, lstm, b) = (&self.gpu, &self.lstm[layer], &mut self.buffers);
        let steps = PILEUP_POSITIONS;
        let gates = 8 * lstm.hidden;
        let pre = b.pre.ptr(&gpu.stream, n * steps * gates)?;
        gpu.dense(x, lstm.input.inputs, &lstm.input, (pre, gates), n * steps, Activation::None)?;
        let cell = b.cell.ptr(&gpu.stream, 2 * n * lstm.hidden)?;
        let whh = lstm.hidden_weight.ptr(&gpu.stream);
        let int = |v: usize| i32::try_from(v).map_err(backend);
        let (n_, steps_, hidden) = (int(n)?, int(steps)?, int(lstm.hidden)?);
        let round = Rounding::Float32 as i32;
        let config = Gpu::config((n.div_ceil(64), lstm.hidden.div_ceil(32), 2))?;
        for s in 0..steps_ {
            let mut launch = gpu.stream.launch_builder(&gpu.kernels.lstm_step);
            launch.arg(&pre).arg(&whh).arg(&out).arg(&cell).arg(&n_).arg(&steps_).arg(&hidden).arg(&s).arg(&round);
            // SAFETY: the arguments match `lstm_step`'s parameters; `pre`, `out` and `cell` hold
            // `n` sequences of `steps` steps of this layer.
            unsafe { launch.launch(config) }.map_err(backend)?;
        }
        Ok(())
    }

    fn forward(&mut self, inputs: &[i32], n: usize) -> Result<Vec<f32>, Error> {
        let stream = self.gpu.stream.clone();
        let b = &mut self.buffers;
        let counts = b.counts.ptr(&stream, inputs.len())?;
        b.upload.upload(&stream, inputs, counts)?;
        let x = b.input.ptr(&stream, inputs.len())?;
        self.gpu.to_float(counts, x, inputs.len(), None)?;
        let width1 = PILEUP_POSITIONS * 2 * self.lstm[0].hidden;
        let width2 = PILEUP_POSITIONS * 2 * self.lstm[1].hidden;
        let out1 = self.buffers.lstm1.ptr(&stream, n * width1)?;
        let out2 = self.buffers.lstm2.ptr(&stream, n * width2)?;
        self.lstm(0, x, n, out1)?;
        self.lstm(1, out1, n, out2)?;

        let (gpu, b) = (&self.gpu, &mut self.buffers);
        let l4 = b.l4.ptr(&stream, n * self.l4.outputs)?;
        gpu.dense(out2, width2, &self.l4, (l4, self.l4.outputs), n, Activation::Selu)?;
        let l5 = b.l5.ptr(&stream, n * self.l5.outputs)?;
        gpu.dense(l4, self.l4.outputs, &self.l5, (l5, self.l5.outputs), n, Activation::Selu)?;
        let probabilities = b.probabilities.ptr(&stream, n * PILEUP_OUTPUTS)?;
        let logits = b.logits.ptr(&stream, n * self.gt21.outputs.max(self.genotype.outputs))?;
        let half = self.gt21.inputs;
        let mut offset = 0;
        for (head, at) in [(&self.gt21, 0), (&self.genotype, half)] {
            let input = l5 + (at * size_of::<f32>()) as u64;
            gpu.dense(input, self.l5.outputs, head, (logits, head.outputs), n, Activation::Selu)?;
            gpu.softmax(logits, probabilities, n, head.outputs, (PILEUP_OUTPUTS, offset))?;
            offset += head.outputs;
        }
        let out = b.download.download(&stream, probabilities, n * PILEUP_OUTPUTS)?;
        Ok(out)
    }
}

impl PileupModel for GpuPileupModel {
    fn predict(&mut self, inputs: &[i32]) -> Result<Vec<f32>, Error> {
        if !inputs.len().is_multiple_of(PILEUP_TENSOR_LEN) {
            return Err(Error::InputLength { len: inputs.len(), tensor_len: PILEUP_TENSOR_LEN });
        }
        let n = inputs.len() / PILEUP_TENSOR_LEN;
        if n == 0 {
            return Ok(Vec::new());
        }
        self.forward(inputs, n)
    }
}

//! The bi-LSTM pileup network on the GPU.

use std::path::Path;

use super::{Activation, Buffer, DeviceDense, Metal, MetalBuffer, Pass, View, int};
use crate::cpu::CpuPileupModel;
use crate::cpu::pileup::BiLstm;
use crate::{Error, PILEUP_OUTPUTS, PILEUP_POSITIONS, PILEUP_TENSOR_LEN, PileupModel};

/// A bidirectional LSTM layer on the GPU.
#[derive(Debug)]
struct DeviceLstm {
    /// Both directions' `W_ih` side by side, with both biases: `[inputs][2 * 4 * hidden]`.
    input: DeviceDense,
    /// `W_hh` transposed, per direction: `[hidden][4 * hidden]`.
    hidden_weight: [MetalBuffer; 2],
    hidden: usize,
}

impl DeviceLstm {
    fn new(metal: &Metal, lstm: &BiLstm) -> Result<Self, Error> {
        let [forward, reverse] = &lstm.hidden_weight;
        Ok(Self {
            input: DeviceDense::from_parts(metal, &lstm.input_weight, &lstm.bias)?,
            hidden_weight: [metal.upload(forward)?, metal.upload(reverse)?],
            hidden: lstm.hidden,
        })
    }

    /// The layer over `x` (`[n][steps][inputs]`) into `out` (`[n][steps][2 * hidden]`), with the
    /// gates in `pre` and the cell states in `cell`.
    fn forward(
        &self,
        pass: &Pass<'_>,
        x: &View,
        n: usize,
        out: &View,
        (pre, cell): (&View, &View),
    ) -> Result<(), Error> {
        let steps = PILEUP_POSITIONS;
        let (hidden, gates) = (self.hidden, 4 * self.hidden);
        pass.dense((x, self.input.inputs), &self.input, (pre, 2 * gates), n * steps, Activation::None)?;
        for s in 0..steps {
            if s > 0 {
                for (d, w_hh) in self.hidden_weight.iter().enumerate() {
                    // The gates of step s, plus the previous step's hidden state times `W_hh`.
                    let (t, prev) = if d == 0 { (s, s - 1) } else { (steps - 1 - s, steps - s) };
                    let h_prev = out.at(prev * 2 * hidden + d * hidden);
                    let g = pre.at(t * 2 * gates + d * gates);
                    let shape = (n, gates, hidden);
                    pass.product(
                        (&h_prev, steps * 2 * hidden),
                        &View::of(w_hh),
                        None,
                        (&g, steps * 2 * gates),
                        Some(&g),
                        shape,
                        Activation::None,
                    )?;
                }
            }
            let params = [int(n)?, int(steps)?, int(hidden)?, int(s)?];
            pass.elementwise(&pass.metal.kernels.lstm_cell, &[pre, out, cell], &params, n * hidden, 2);
        }
        Ok(())
    }
}

/// Activations, reused across calls.
#[derive(Debug, Default)]
struct Buffers {
    counts: Buffer,
    input: Buffer,
    /// Gates of the layer being run: `[n][steps][2 * 4 * hidden]`.
    pre: Buffer,
    cell: Buffer,
    lstm1: Buffer,
    lstm2: Buffer,
    l4: Buffer,
    l5: Buffer,
    logits: Buffer,
    probabilities: Buffer,
}

/// The original's bi-LSTM pileup network (its `*_P` class), in Rust and Metal kernels on an
/// Apple GPU.
#[derive(Debug)]
pub struct MetalPileupModel {
    metal: Metal,
    lstm: [DeviceLstm; 2],
    l4: DeviceDense,
    /// `L5_1` and `L5_2` side by side.
    l5: DeviceDense,
    gt21: DeviceDense,
    genotype: DeviceDense,
    buffers: Buffers,
}

// SAFETY: Metal's devices, command queues, pipeline states and buffers may be used from any
// thread (Apple's "Metal objects are thread-safe" apart from command buffers and encoders, which
// a call creates and finishes on its own thread). The model is used by one thread at a time, as
// `&mut`, so its `RefCell` is never shared.
unsafe impl Send for MetalPileupModel {}

impl MetalPileupModel {
    /// Loads a pileup checkpoint (`pileup.pt`) onto the system's GPU.
    ///
    /// # Errors
    ///
    /// As [`CpuPileupModel::load`], and fails if there is no Metal GPU or the kernels cannot be
    /// compiled or the weights copied.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, Error> {
        Self::from_cpu(&CpuPileupModel::load(path)?)
    }

    /// The network of `cpu` on the system's GPU.
    pub(crate) fn from_cpu(cpu: &CpuPileupModel) -> Result<Self, Error> {
        let metal = Metal::new()?;
        Ok(Self {
            lstm: [DeviceLstm::new(&metal, &cpu.lstm1)?, DeviceLstm::new(&metal, &cpu.lstm2)?],
            l4: DeviceDense::new(&metal, &cpu.l4)?,
            l5: DeviceDense::side_by_side(&metal, &cpu.gt21.hidden, &cpu.genotype.hidden)?,
            gt21: DeviceDense::new(&metal, &cpu.gt21.logits)?,
            genotype: DeviceDense::new(&metal, &cpu.genotype.logits)?,
            metal,
            buffers: Buffers::default(),
        })
    }

    fn forward(&mut self, inputs: &[i32], n: usize) -> Result<Vec<f32>, Error> {
        let (metal, b) = (&self.metal, &mut self.buffers);
        let steps = PILEUP_POSITIONS;
        let counts = b.counts.view::<i32>(metal, inputs.len())?;
        counts.write(inputs);
        let x = b.input.view::<f32>(metal, inputs.len())?;
        let widest = self.lstm.iter().map(|l| l.hidden).max().unwrap_or(0);
        let pre = b.pre.view::<f32>(metal, n * steps * 8 * widest)?;
        let cell = b.cell.view::<f32>(metal, 2 * n * widest)?;
        let width1 = steps * 2 * self.lstm[0].hidden;
        let width2 = steps * 2 * self.lstm[1].hidden;
        let out1 = b.lstm1.view::<f32>(metal, n * width1)?;
        let out2 = b.lstm2.view::<f32>(metal, n * width2)?;
        let l4 = b.l4.view::<f32>(metal, n * self.l4.outputs)?;
        let l5 = b.l5.view::<f32>(metal, n * self.l5.outputs)?;
        let logits = b.logits.view::<f32>(metal, n * self.gt21.outputs.max(self.genotype.outputs))?;
        let probabilities = b.probabilities.view::<f32>(metal, n * PILEUP_OUTPUTS)?;

        let pass = metal.pass()?;
        let len = u32::try_from(inputs.len()).map_err(|_| super::backend(format!("{} inputs", inputs.len())))?;
        pass.elementwise(&metal.kernels.int32_to_float, &[&counts, &x], &len, inputs.len(), 1);
        self.lstm[0].forward(&pass, &x, n, &out1, (&pre, &cell))?;
        self.lstm[1].forward(&pass, &out1, n, &out2, (&pre, &cell))?;
        pass.dense((&out2, width2), &self.l4, (&l4, self.l4.outputs), n, Activation::Selu)?;
        pass.dense((&l4, self.l4.outputs), &self.l5, (&l5, self.l5.outputs), n, Activation::Selu)?;
        let mut offset = 0;
        for (head, at) in [(&self.gt21, 0), (&self.genotype, self.gt21.inputs)] {
            pass.dense((&l5.at(at), self.l5.outputs), head, (&logits, head.outputs), n, Activation::Selu)?;
            pass.softmax(&logits, &probabilities, n, head.outputs, (PILEUP_OUTPUTS, offset))?;
            offset += head.outputs;
        }
        pass.finish()?;
        Ok(probabilities.read(n * PILEUP_OUTPUTS))
    }
}

impl PileupModel for MetalPileupModel {
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

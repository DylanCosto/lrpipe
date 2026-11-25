//! The residual CNN full-alignment network on the GPU.

use std::path::Path;

use super::{Activation, Buffer, DeviceConv, DeviceDense, Metal, backend, int};
use crate::cpu::CpuFullAlignmentModel;
use crate::{
    Error, FULL_ALIGNMENT_DEPTH, FULL_ALIGNMENT_OUTPUTS, FULL_ALIGNMENT_POSITIONS, FullAlignmentModel,
    full_alignment_tensor_len,
};

/// `param_f.NORMALIZE_NUM`: the inputs are divided by it after the cast to `float32`.
const NORMALIZE_NUM: f32 = 100.0;

/// Activations, reused across calls.
#[derive(Debug, Default)]
struct Buffers {
    tensors: Buffer,
    input: Buffer,
    maps: [Buffer; 3],
    pooled: Buffer,
    l4: Buffer,
    l5: Buffer,
    logits: Buffer,
    probabilities: Buffer,
}

/// The original's residual CNN full-alignment network (its `*_F` class), in Rust and Metal
/// kernels on an Apple GPU.
#[derive(Debug)]
pub struct MetalFullAlignmentModel {
    metal: Metal,
    conv1: DeviceConv,
    /// Three residual blocks' convolutions.
    blocks: [(DeviceConv, DeviceConv); 3],
    /// The strided convolutions after the first two blocks.
    strided: [DeviceConv; 2],
    l4: DeviceDense,
    /// `L5_1` to `L5_4` side by side.
    l5: DeviceDense,
    /// GT21, genotype, indel length 1 and 2 logits.
    heads: [DeviceDense; 4],
    matrix_depth: usize,
    channels: usize,
    buffers: Buffers,
}

// SAFETY: Metal's devices, command queues, pipeline states and buffers may be used from any
// thread (Apple's "Metal objects are thread-safe" apart from command buffers and encoders, which
// a call creates and finishes on its own thread). The model is used by one thread at a time, as
// `&mut`, so its `RefCell` is never shared.
unsafe impl Send for MetalFullAlignmentModel {}

impl MetalFullAlignmentModel {
    /// Loads a full-alignment checkpoint (`full_alignment.pt`) onto the system's GPU.
    ///
    /// # Errors
    ///
    /// As [`CpuFullAlignmentModel::load`], and fails if there is no Metal GPU or the kernels
    /// cannot be compiled or the weights copied.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, Error> {
        Self::from_cpu(&CpuFullAlignmentModel::load(path)?)
    }

    /// The network of `cpu` on the system's GPU.
    pub(crate) fn from_cpu(cpu: &CpuFullAlignmentModel) -> Result<Self, Error> {
        let metal = Metal::new()?;
        let conv = |c| DeviceConv::new(&metal, c);
        let block = |(a, b)| -> Result<_, Error> { Ok((conv(a)?, conv(b)?)) };
        let [b1, b2, b3] = &cpu.trunk.blocks;
        let [s1, s2] = &cpu.trunk.strided;
        let l5_parts: Vec<f32> = {
            let width: usize = cpu.heads.iter().map(|h| h.hidden.bias.len()).sum();
            let inputs = cpu.l4.bias.len();
            let mut weight = Vec::with_capacity(inputs * width);
            for row in 0..inputs {
                for head in &cpu.heads {
                    let n = head.hidden.bias.len();
                    weight.extend_from_slice(&head.hidden.weight[row * n..(row + 1) * n]);
                }
            }
            weight
        };
        let l5_bias: Vec<f32> = cpu.heads.iter().flat_map(|h| h.hidden.bias.iter().copied()).collect();
        let [h1, h2, h3, h4] = &cpu.heads;
        Ok(Self {
            conv1: conv(&cpu.conv1)?,
            blocks: [block((&b1.0, &b1.1))?, block((&b2.0, &b2.1))?, block((&b3.0, &b3.1))?],
            strided: [conv(s1)?, conv(s2)?],
            l4: DeviceDense::new(&metal, &cpu.l4)?,
            l5: DeviceDense::from_parts(&metal, &l5_parts, &l5_bias)?,
            heads: [
                DeviceDense::new(&metal, &h1.logits)?,
                DeviceDense::new(&metal, &h2.logits)?,
                DeviceDense::new(&metal, &h3.logits)?,
                DeviceDense::new(&metal, &h4.logits)?,
            ],
            matrix_depth: FULL_ALIGNMENT_DEPTH,
            channels: cpu.conv1.in_channels,
            metal,
            buffers: Buffers::default(),
        })
    }

    /// Scores matrices of `depth` rows instead of [`FULL_ALIGNMENT_DEPTH`].
    #[must_use]
    pub fn with_matrix_depth(mut self, depth: usize) -> Self {
        self.matrix_depth = depth;
        self
    }

    fn forward(&mut self, inputs: &[i8], n: usize) -> Result<Vec<f32>, Error> {
        let (metal, b) = (&self.metal, &mut self.buffers);
        let tensors = b.tensors.view::<i8>(metal, inputs.len())?;
        tensors.write(inputs);
        let x = b.input.view::<f32>(metal, inputs.len())?;
        // The largest map is conv1's output.
        let (h1, w1) = ((self.matrix_depth - 1) / 2 + 1, (FULL_ALIGNMENT_POSITIONS - 1) / 2 + 1);
        let largest = n * h1 * w1 * self.conv1.outputs.max(self.blocks[0].0.outputs);
        let [m0, m1, m2] = &mut b.maps;
        let mut maps =
            [m0.view::<f32>(metal, largest)?, m1.view::<f32>(metal, largest)?, m2.view::<f32>(metal, largest)?];
        let channels = self.blocks[2].1.outputs;
        let features = 14 * channels;
        let pooled = b.pooled.view::<f32>(metal, n * features)?;
        let l4 = b.l4.view::<f32>(metal, n * self.l4.outputs)?;
        let l5 = b.l5.view::<f32>(metal, n * self.l5.outputs)?;
        let classes = self.heads.iter().map(|h| h.outputs).max().unwrap_or(0);
        let logits = b.logits.view::<f32>(metal, n * classes)?;
        let probabilities = b.probabilities.view::<f32>(metal, n * FULL_ALIGNMENT_OUTPUTS)?;

        let pass = metal.pass()?;
        let len = u32::try_from(inputs.len()).map_err(|_| backend(format!("{} inputs", inputs.len())))?;
        // `x.float() / NORMALIZE_NUM`.
        let normalize = NormalizeParams { len, divisor: NORMALIZE_NUM };
        pass.elementwise(&metal.kernels.int8_to_float, &[&tensors, &x], &normalize, inputs.len(), 1);
        let mut dims = pass.conv(&self.conv1, &x, (n, self.matrix_depth, FULL_ALIGNMENT_POSITIONS), &maps[0], None)?;
        for (i, (conv1, conv2)) in self.blocks.iter().enumerate() {
            let [x, y, z] = &maps;
            pass.conv(conv1, x, (n, dims.0, dims.1), y, None)?;
            pass.conv(conv2, y, (n, dims.0, dims.1), z, Some(x))?;
            match self.strided.get(i) {
                Some(strided) => dims = pass.conv(strided, z, (n, dims.0, dims.1), x, None)?,
                None => maps.swap(0, 2),
            }
        }
        let params = [int(n)?, int(dims.0)?, int(dims.1)?, int(channels)?];
        pass.elementwise(&metal.kernels.pyramid_pool, &[&maps[0], &pooled], &params, n * features, 1);
        pass.dense((&pooled, features), &self.l4, (&l4, self.l4.outputs), n, Activation::Selu)?;
        pass.dense((&l4, self.l4.outputs), &self.l5, (&l5, self.l5.outputs), n, Activation::Selu)?;
        let (mut hidden_at, mut offset) = (0, 0);
        for head in &self.heads {
            pass.dense((&l5.at(hidden_at), self.l5.outputs), head, (&logits, head.outputs), n, Activation::Selu)?;
            pass.softmax(&logits, &probabilities, n, head.outputs, (FULL_ALIGNMENT_OUTPUTS, offset))?;
            hidden_at += head.inputs;
            offset += head.outputs;
        }
        pass.finish()?;
        Ok(probabilities.read(n * FULL_ALIGNMENT_OUTPUTS))
    }
}

/// `NormalizeParams` of the kernels: `int8_to_float`'s length and divisor.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct NormalizeParams {
    len: u32,
    divisor: f32,
}

impl FullAlignmentModel for MetalFullAlignmentModel {
    fn predict(&mut self, inputs: &[i8]) -> Result<Vec<f32>, Error> {
        let tensor_len = full_alignment_tensor_len(self.matrix_depth, self.channels);
        if !inputs.len().is_multiple_of(tensor_len) {
            return Err(Error::InputLength { len: inputs.len(), tensor_len });
        }
        let n = inputs.len() / tensor_len;
        if n == 0 {
            return Ok(Vec::new());
        }
        self.forward(inputs, n)
    }

    fn channels(&self) -> usize {
        self.channels
    }
}

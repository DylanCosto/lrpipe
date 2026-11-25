//! The residual CNN full-alignment network on the GPU.

use std::path::Path;

use cudarc::driver::PushKernelArg;

use super::{Activation, Buffer, DeviceConv, DeviceDense, Direction, Gpu, NULL, Staging, backend};
use crate::cpu::CpuFullAlignmentModel;
use crate::{
    Error, FULL_ALIGNMENT_DEPTH, FULL_ALIGNMENT_OUTPUTS, FULL_ALIGNMENT_POSITIONS, FullAlignmentModel,
    full_alignment_tensor_len,
};

/// `param_f.NORMALIZE_NUM`: the inputs are divided by it after the cast to `float32`.
const NORMALIZE_NUM: f32 = 100.0;

/// Activations, reused across calls.
#[derive(Debug)]
struct Buffers {
    upload: Staging<i8>,
    download: Staging<f32>,
    tensors: Buffer<i8>,
    input: Buffer<f32>,
    maps: [Buffer<f32>; 3],
    pooled: Buffer<f32>,
    l4: Buffer<f32>,
    l5: Buffer<f32>,
    logits: Buffer<f32>,
    probabilities: Buffer<f32>,
}

/// The original's residual CNN full-alignment network (its `*_F` class), in Rust kernels on a
/// CUDA GPU.
#[derive(Debug)]
pub struct GpuFullAlignmentModel {
    gpu: Gpu,
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

impl GpuFullAlignmentModel {
    /// Loads a full-alignment checkpoint (`full_alignment.pt`) onto CUDA device `device`.
    ///
    /// # Errors
    ///
    /// As [`CpuFullAlignmentModel::load`], and fails if the CUDA driver or NVRTC is missing, the
    /// device does not exist, or the kernels cannot be compiled or the weights copied.
    pub fn load(path: impl AsRef<Path>, device: usize) -> Result<Self, Error> {
        Self::from_cpu(&CpuFullAlignmentModel::load(path)?, device)
    }

    /// The network of `cpu` on CUDA device `device`.
    pub(crate) fn from_cpu(cpu: &CpuFullAlignmentModel, device: usize) -> Result<Self, Error> {
        Self::on(cpu, Gpu::new(device)?)
    }

    /// The network of `cpu` with its kernels compiled for compute capability `capability`, older
    /// than device `device`'s, as on an older GPU.
    #[cfg(test)]
    pub(crate) fn compiled_for(
        cpu: &CpuFullAlignmentModel,
        device: usize,
        capability: (i32, i32),
    ) -> Result<Self, Error> {
        Self::on(cpu, Gpu::compiled_for(device, capability)?)
    }

    /// The network of `cpu` on `gpu`.
    fn on(cpu: &CpuFullAlignmentModel, gpu: Gpu) -> Result<Self, Error> {
        let conv = |c| DeviceConv::new(&gpu, c);
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
            l4: DeviceDense::new(&gpu, &cpu.l4)?,
            l5: DeviceDense::from_parts(&gpu, &l5_parts, &l5_bias)?,
            heads: [
                DeviceDense::new(&gpu, &h1.logits)?,
                DeviceDense::new(&gpu, &h2.logits)?,
                DeviceDense::new(&gpu, &h3.logits)?,
                DeviceDense::new(&gpu, &h4.logits)?,
            ],
            matrix_depth: FULL_ALIGNMENT_DEPTH,
            channels: cpu.conv1.in_channels,
            gpu,
            buffers: Buffers {
                upload: Staging::new(Direction::ToDevice),
                download: Staging::new(Direction::ToHost),
                tensors: Buffer::default(),
                input: Buffer::default(),
                maps: Default::default(),
                pooled: Buffer::default(),
                l4: Buffer::default(),
                l5: Buffer::default(),
                logits: Buffer::default(),
                probabilities: Buffer::default(),
            },
        })
    }

    /// Scores matrices of `depth` rows instead of [`FULL_ALIGNMENT_DEPTH`].
    #[must_use]
    pub fn with_matrix_depth(mut self, depth: usize) -> Self {
        self.matrix_depth = depth;
        self
    }

    fn forward(&mut self, inputs: &[i8], n: usize) -> Result<Vec<f32>, Error> {
        let (gpu, b) = (&self.gpu, &mut self.buffers);
        let stream = gpu.stream.clone();
        let tensors = b.tensors.ptr(&stream, inputs.len())?;
        b.upload.upload(&stream, inputs, tensors)?;
        let x = b.input.ptr(&stream, inputs.len())?;
        gpu.to_float(tensors, x, inputs.len(), Some(NORMALIZE_NUM))?;

        // The largest map is conv1's output.
        let (h1, w1) = ((self.matrix_depth - 1) / 2 + 1, (FULL_ALIGNMENT_POSITIONS - 1) / 2 + 1);
        let largest = n * h1 * w1 * self.conv1.outputs.max(self.blocks[0].0.outputs);
        let [m0, m1, m2] = &mut b.maps;
        let mut maps = [m0.ptr(&stream, largest)?, m1.ptr(&stream, largest)?, m2.ptr(&stream, largest)?];
        let mut dims = gpu.conv(&self.conv1, x, (n, self.matrix_depth, FULL_ALIGNMENT_POSITIONS), maps[0], NULL)?;
        for (i, (conv1, conv2)) in self.blocks.iter().enumerate() {
            let [x, y, z] = maps;
            gpu.conv(conv1, x, (n, dims.0, dims.1), y, NULL)?;
            gpu.conv(conv2, y, (n, dims.0, dims.1), z, x)?;
            match self.strided.get(i) {
                Some(strided) => dims = gpu.conv(strided, z, (n, dims.0, dims.1), x, NULL)?,
                None => maps = [z, y, x],
            }
        }
        let channels = self.blocks[2].1.outputs;
        let features = 14 * channels;
        let pooled = b.pooled.ptr(&stream, n * features)?;
        let int = |v: usize| i32::try_from(v).map_err(backend);
        let (n_, h, w, c) = (int(n)?, int(dims.0)?, int(dims.1)?, int(channels)?);
        let mut launch = stream.launch_builder(&gpu.kernels.pyramid_pool);
        launch.arg(&maps[0]).arg(&pooled).arg(&n_).arg(&h).arg(&w).arg(&c);
        // SAFETY: the arguments match `pyramid_pool`'s parameters; the map holds `n` maps of
        // `h x w x c` and `pooled` `n x 14 c` values.
        unsafe { launch.launch(Gpu::config((Gpu::blocks(n * features), 1, 1))?) }.map_err(backend)?;

        let l4 = b.l4.ptr(&stream, n * self.l4.outputs)?;
        gpu.dense(pooled, features, &self.l4, (l4, self.l4.outputs), n, Activation::Selu)?;
        let l5 = b.l5.ptr(&stream, n * self.l5.outputs)?;
        gpu.dense(l4, self.l4.outputs, &self.l5, (l5, self.l5.outputs), n, Activation::Selu)?;
        let probabilities = b.probabilities.ptr(&stream, n * FULL_ALIGNMENT_OUTPUTS)?;
        let classes = self.heads.iter().map(|h| h.outputs).max().unwrap_or(0);
        let logits = b.logits.ptr(&stream, n * classes)?;
        let (mut hidden_at, mut offset) = (0, 0);
        for head in &self.heads {
            let input = l5 + (hidden_at * size_of::<f32>()) as u64;
            gpu.dense(input, self.l5.outputs, head, (logits, head.outputs), n, Activation::Selu)?;
            gpu.softmax(logits, probabilities, n, head.outputs, (FULL_ALIGNMENT_OUTPUTS, offset))?;
            hidden_at += head.inputs;
            offset += head.outputs;
        }
        let out = b.download.download(&stream, probabilities, n * FULL_ALIGNMENT_OUTPUTS)?;
        Ok(out)
    }
}

impl FullAlignmentModel for GpuFullAlignmentModel {
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

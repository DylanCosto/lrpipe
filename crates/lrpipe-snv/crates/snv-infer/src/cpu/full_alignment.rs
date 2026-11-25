//! The residual CNN full-alignment network (the `*_F` class in the original's `model.py`).
//!
//! Each candidate runs through the convolutions on its own, in `[height][width][channel]`
//! layout; pyramid pooling, `L4` and the heads then run on the whole batch.

use std::path::Path;

use super::conv::ConvScratch;
use super::{Conv, Dense, Epilogue, Head, Weights, selu};
use crate::int8::{self, QuantConv, Scratch};
use crate::{
    Error, FULL_ALIGNMENT_CHANNELS, FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME, FULL_ALIGNMENT_DEPTH,
    FULL_ALIGNMENT_OUTPUTS, FULL_ALIGNMENT_POSITIONS, FullAlignmentModel, Precision, full_alignment_tensor_len,
    top_two_gap,
};

/// `param_f.NORMALIZE_NUM`: the inputs are divided by it after the cast to `float32`.
const NORMALIZE_NUM: f32 = 100.0;

/// With 8-bit convolutions, candidates whose two most likely GT21 or genotype classes are
/// closer than this are scored again in `float32`.
const RESCORE_GAP: f32 = 0.1;

/// Candidates that go through each convolution together, so that the deeper layers, whose
/// weights outweigh one candidate's activations, read their weights once for several.
const TRUNK_BATCH: usize = 4;

/// `PyramidPolling.spatial_pool_size`: output grid sizes of the three max pools.
const PYRAMID_POOL_SIZES: [usize; 3] = [3, 2, 1];

/// A residual block's two convolutions (`BasicBlock` with stride 1 and an identity shortcut).
pub(crate) type Block<C> = (C, C);

/// The convolutions after `conv1`: three residual blocks with the two strided convolutions
/// between them.
#[derive(Debug)]
pub(crate) struct Trunk<C> {
    pub(crate) blocks: [Block<C>; 3],
    pub(crate) strided: [C; 2],
}

impl<C> Trunk<C> {
    fn map<D>(&self, mut f: impl FnMut(&C) -> D) -> Trunk<D> {
        Trunk { blocks: self.blocks.each_ref().map(|(a, b)| (f(a), f(b))), strided: self.strided.each_ref().map(f) }
    }
}

/// One convolution kind's forward pass, for [`Trunk::forward`].
pub(crate) trait TrunkConv {
    type Scratch: Default + std::fmt::Debug;

    fn output_size(&self, height: usize, width: usize) -> (usize, usize);

    /// As [`Conv::forward`]: `images` maps of `height` x `width` in `input`.
    fn forward(
        &self,
        input: &[f32],
        dims: (usize, usize, usize),
        epilogue: Epilogue<'_>,
        output: &mut Vec<f32>,
        scratch: &mut Self::Scratch,
    );
}

impl TrunkConv for Conv {
    type Scratch = ConvScratch;

    fn output_size(&self, height: usize, width: usize) -> (usize, usize) {
        Conv::output_size(self, height, width)
    }

    fn forward(
        &self,
        x: &[f32],
        dims: (usize, usize, usize),
        e: Epilogue<'_>,
        out: &mut Vec<f32>,
        s: &mut ConvScratch,
    ) {
        Conv::forward(self, x, dims, e, out, s);
    }
}

impl TrunkConv for QuantConv {
    type Scratch = Scratch;

    fn output_size(&self, height: usize, width: usize) -> (usize, usize) {
        QuantConv::output_size(self, height, width)
    }

    /// Each image's activations are quantized on their own scale.
    fn forward(&self, x: &[f32], dims: (usize, usize, usize), e: Epilogue<'_>, out: &mut Vec<f32>, s: &mut Scratch) {
        QuantConv::forward(self, x, dims, e, out, s);
    }
}

/// Buffers reused across candidates.
#[derive(Debug, Default)]
struct Buffers<S> {
    scratch: S,
    activations: [Vec<f32>; 3],
}

impl<C: TrunkConv> Trunk<C> {
    /// The `conv1` output of `images` candidates (`[image][height][width][64]`, in
    /// `buffers.activations[0]`) through the trunk; the result is left in
    /// `buffers.activations[0]`. Returns its height and width.
    fn forward(&self, buffers: &mut Buffers<C::Scratch>, images: usize, mut dims: (usize, usize)) -> (usize, usize) {
        let [x, y, z] = &mut buffers.activations;
        for (i, (conv1, conv2)) in self.blocks.iter().enumerate() {
            let at = (images, dims.0, dims.1);
            conv1.forward(x, at, Epilogue::Relu, y, &mut buffers.scratch);
            conv2.forward(y, at, Epilogue::ResidualRelu(x), z, &mut buffers.scratch);
            match self.strided.get(i) {
                Some(strided) => {
                    strided.forward(z, at, Epilogue::Relu, x, &mut buffers.scratch);
                    dims = strided.output_size(dims.0, dims.1);
                }
                None => std::mem::swap(x, z),
            }
        }
        dims
    }
}

/// `PyramidPolling.forward` on one `[height][width][channels]` feature map: TF-style "same" max
/// pools to 3x3, 2x2 and 1x1 grids (zero padding split with the extra row or column at the end),
/// each flattened in `[row][column][channel]` order, appended to `out`.
fn pyramid_pooling(x: &[f32], (height, width): (usize, usize), channels: usize, out: &mut Vec<f32>) {
    for pool_size in PYRAMID_POOL_SIZES {
        // `int(np.ceil(a / b))` is exact for these small positive integers.
        let (window_h, window_w) = (height.div_ceil(pool_size), width.div_ceil(pool_size));
        let (out_h, out_w) = (height.div_ceil(window_h), width.div_ceil(window_w));
        let pad_top = ((out_h - 1) * window_h + window_h).saturating_sub(height) / 2;
        let pad_left = ((out_w - 1) * window_w + window_w).saturating_sub(width) / 2;
        for oy in 0..out_h {
            let rows =
                (oy * window_h).saturating_sub(pad_top)..((oy + 1) * window_h).saturating_sub(pad_top).min(height);
            let padded_rows = rows.len() < window_h;
            for ox in 0..out_w {
                let columns =
                    (ox * window_w).saturating_sub(pad_left)..((ox + 1) * window_w).saturating_sub(pad_left).min(width);
                // A window over padding also takes the padding's zero.
                let start = if padded_rows || columns.len() < window_w { 0.0 } else { f32::NEG_INFINITY };
                let at = out.len();
                out.resize(at + channels, start);
                let pooled = &mut out[at..];
                for iy in rows.clone() {
                    for ix in columns.clone() {
                        for (m, &v) in pooled.iter_mut().zip(&x[(iy * width + ix) * channels..][..channels]) {
                            *m = m.max(v);
                        }
                    }
                }
            }
        }
    }
}

/// The original's residual CNN full-alignment network (its `*_F` class with
/// `add_indel_length=True`, `predict=True`), with or without the dwell-time channel, in Rust on
/// the CPU.
#[derive(Debug)]
pub struct CpuFullAlignmentModel {
    pub(crate) conv1: Conv,
    pub(crate) trunk: Trunk<Conv>,
    pub(crate) l4: Dense,
    /// GT21, genotype, indel length 1 and indel length 2, in output order.
    pub(crate) heads: [Head; 4],
    /// Rows of each input matrix: [`FULL_ALIGNMENT_DEPTH`] for ONT, 55 for HiFi and Illumina.
    matrix_depth: usize,
    /// The 8-bit trunk, for [`Precision::Mixed`] on a CPU with AVX-VNNI or Arm's `sdot`.
    quantized: Option<Trunk<QuantConv>>,
    buffers: Buffers<ConvScratch>,
    quantized_buffers: Buffers<Scratch>,
}

impl CpuFullAlignmentModel {
    /// Loads a full-alignment checkpoint (`full_alignment.pt`).
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read as a PyTorch state dict, lacks a required tensor, has
    /// tensors of shapes the network cannot use, or was trained for a number of input channels
    /// other than [`FULL_ALIGNMENT_CHANNELS`] or [`FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME`].
    pub fn load(path: impl AsRef<Path>) -> Result<Self, Error> {
        Self::from_state(&mut Weights::load(path.as_ref())?)
    }

    pub(crate) fn from_state(state: &mut Weights) -> Result<Self, Error> {
        let conv1 = Conv::from_state(state, "conv1.conv", "conv1.bn", 2)?;
        if ![FULL_ALIGNMENT_CHANNELS, FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME].contains(&conv1.in_channels) {
            return Err(state.bad(format!(
                "conv1 takes {} input channels, expected {FULL_ALIGNMENT_CHANNELS} or \
                 {FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME}",
                conv1.in_channels
            )));
        }
        let mut block = |name: &str| -> Result<Block<Conv>, Error> {
            Ok((
                Conv::from_state(state, &format!("{name}.conv1"), &format!("{name}.bn1"), 1)?,
                Conv::from_state(state, &format!("{name}.conv2"), &format!("{name}.bn2"), 1)?,
            ))
        };
        let blocks = [block("res_block1.0")?, block("res_block2.0")?, block("res_block3.0")?];
        let strided = [
            Conv::from_state(state, "conv3.conv", "conv3.bn", 2)?,
            Conv::from_state(state, "conv5.conv", "conv5.bn", 2)?,
        ];
        let model = Self {
            conv1,
            trunk: Trunk { blocks, strided },
            l4: Dense::from_state(state, "L4")?,
            heads: [
                Head::from_state(state, "L5_1", "Y_gt21_logits")?,
                Head::from_state(state, "L5_2", "Y_genotype_logits")?,
                Head::from_state(state, "L5_3", "Y_indel_length_logits_1")?,
                Head::from_state(state, "L5_4", "Y_indel_length_logits_2")?,
            ],
            matrix_depth: FULL_ALIGNMENT_DEPTH,
            quantized: None,
            buffers: Buffers::default(),
            quantized_buffers: Buffers::default(),
        };
        let mut channels = model.conv1.out_channels;
        for (i, (a, b)) in model.trunk.blocks.iter().enumerate() {
            let chained = [a.in_channels, a.out_channels, b.in_channels, b.out_channels].iter().all(|&c| c == channels);
            if !chained {
                return Err(state.bad(format!("residual block {} does not keep {channels} channels", i + 1)));
            }
            if let Some(strided) = model.trunk.strided.get(i) {
                if strided.in_channels != channels {
                    return Err(state.bad(format!("strided convolution {} does not take {channels} channels", i + 1)));
                }
                channels = strided.out_channels;
            }
        }
        let pooled = PYRAMID_POOL_SIZES.iter().map(|s| s * s).sum::<usize>() * channels;
        let classes: usize = model.heads.iter().map(Head::classes).sum();
        if model.l4.inputs() != pooled || classes != FULL_ALIGNMENT_OUTPUTS {
            return Err(state.bad("layer sizes do not chain into the full-alignment network".to_owned()));
        }
        Ok(model)
    }

    /// Input channels per position: [`FULL_ALIGNMENT_CHANNELS`], or
    /// [`FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME`] for a `_with_mv` model.
    pub fn channels(&self) -> usize {
        self.conv1.in_channels
    }

    /// Scores matrices of `depth` rows instead of [`FULL_ALIGNMENT_DEPTH`] (`param_f.matrix_depth_dict`:
    /// 55 for HiFi and Illumina). The network pools to a fixed size, so any depth fits it.
    #[must_use]
    pub fn with_matrix_depth(mut self, depth: usize) -> Self {
        self.matrix_depth = depth;
        self
    }

    /// Runs the network at `precision` instead of [`Precision::Float32`]: for
    /// [`Precision::Mixed`] on a CPU with AVX-VNNI or Arm's `sdot`, the convolutions after the
    /// first in 8 bits, with candidates whose call is close scored again in `float32`.
    #[must_use]
    pub fn with_precision(mut self, precision: Precision) -> Self {
        self.quantized = (precision == Precision::Mixed && int8::supported()).then(|| self.trunk.map(QuantConv::new));
        self
    }

    /// Probabilities of the `n` candidates of `inputs`, with the 8-bit trunk if `quantized`.
    fn probabilities(&mut self, inputs: &[i8], n: usize, quantized: Option<&Trunk<QuantConv>>) -> Vec<f32> {
        let tensor_len = full_alignment_tensor_len(self.matrix_depth, self.channels());
        let dims = (self.matrix_depth, FULL_ALIGNMENT_POSITIONS);
        let conv1_dims = self.conv1.output_size(dims.0, dims.1);
        let mut pooled = Vec::new();
        let mut input = Vec::new();
        for group in inputs.chunks(TRUNK_BATCH * tensor_len) {
            let images = group.len() / tensor_len;
            input.clear();
            // `x.float() / NORMALIZE_NUM`; `[depth][position][channel]` is already the layout.
            input.extend(group.iter().map(|&v| f32::from(v) / NORMALIZE_NUM));
            let at = (images, dims.0, dims.1);
            let columns = &mut self.buffers.scratch;
            let (size, maps) = match quantized {
                Some(trunk) => {
                    let buffers = &mut self.quantized_buffers;
                    self.conv1.forward(&input, at, Epilogue::Relu, &mut buffers.activations[0], columns);
                    (trunk.forward(buffers, images, conv1_dims), &buffers.activations[0])
                }
                None => {
                    self.conv1.forward(&input, at, Epilogue::Relu, &mut self.buffers.activations[0], columns);
                    (self.trunk.forward(&mut self.buffers, images, conv1_dims), &self.buffers.activations[0])
                }
            };
            let map_len = maps.len() / images;
            for map in maps.chunks_exact(map_len) {
                pyramid_pooling(map, size, map_len / (size.0 * size.1), &mut pooled);
            }
        }
        let x = self.l4.forward(&pooled, n, selu);
        let mut out = vec![0.0; n * FULL_ALIGNMENT_OUTPUTS];
        let mut at = 0;
        for head in &self.heads {
            head.forward_into(&x, n, &mut out, FULL_ALIGNMENT_OUTPUTS, at);
            at += head.classes();
        }
        out
    }
}

impl FullAlignmentModel for CpuFullAlignmentModel {
    fn predict(&mut self, inputs: &[i8]) -> Result<Vec<f32>, Error> {
        let tensor_len = full_alignment_tensor_len(self.matrix_depth, self.channels());
        if !inputs.len().is_multiple_of(tensor_len) {
            return Err(Error::InputLength { len: inputs.len(), tensor_len });
        }
        let n = inputs.len() / tensor_len;
        let quantized = self.quantized.take();
        let mut out = self.probabilities(inputs, n, quantized.as_ref());
        if quantized.is_some() {
            let close: Vec<usize> = (0..n)
                .filter(|&i| {
                    let p = &out[i * FULL_ALIGNMENT_OUTPUTS..];
                    top_two_gap(&p[..21]) < RESCORE_GAP || top_two_gap(&p[21..24]) < RESCORE_GAP
                })
                .collect();
            let subset: Vec<i8> =
                close.iter().flat_map(|&i| &inputs[i * tensor_len..(i + 1) * tensor_len]).copied().collect();
            let exact = self.probabilities(&subset, close.len(), None);
            for (&i, row) in close.iter().zip(exact.as_chunks::<FULL_ALIGNMENT_OUTPUTS>().0) {
                out[i * FULL_ALIGNMENT_OUTPUTS..(i + 1) * FULL_ALIGNMENT_OUTPUTS].copy_from_slice(row);
            }
        }
        self.quantized = quantized;
        Ok(out)
    }

    fn channels(&self) -> usize {
        self.conv1.in_channels
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pyramid_pooling_pads_right_and_flattens_nhwc() {
        // Two channels over a 12x5 grid, the network's final feature map size: channel 0 holds
        // `5 * h + w`, channel 1 the same plus 100.
        let x: Vec<f32> = (0..60).flat_map(|i| [i as f32, i as f32 + 100.0]).collect();
        let mut got = Vec::new();
        pyramid_pooling(&x, (12, 5), 2, &mut got);

        // Pool to 3x3: 4x2 windows, one zero column padded on the right. Pool to 2x2: 6x3
        // windows, also padded on the right. Pool to 1x1: one 12x5 window.
        let mut want = Vec::new();
        for (rows, cols) in [(&[3, 7, 11][..], &[1, 3, 4][..]), (&[5, 11], &[2, 4]), (&[11], &[4])] {
            for &h in rows {
                for &w in cols {
                    let value = (5 * h + w) as f32;
                    want.extend([value, value + 100.0]);
                }
            }
        }
        assert_eq!(got, want);
    }

    #[test]
    fn padding_counts_as_zero() {
        // A 1x5 map of negative values: pooled to 1x3 with 1x2 windows and to 1x2 with 1x3
        // windows, the last window of each over one column of padding.
        let mut got = Vec::new();
        pyramid_pooling(&[-1.0, -2.0, -3.0, -4.0, -5.0], (1, 5), 1, &mut got);
        assert_eq!(got, [-1.0, -3.0, 0.0, -1.0, 0.0, -1.0]);
    }
}

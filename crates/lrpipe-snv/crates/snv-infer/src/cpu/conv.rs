//! 3x3 convolutions with padding 1 and a folded batch norm, as im2col and one matrix product.

use super::{Product, Strided, Weights, matmul};
use crate::Error;

/// `eps` of every `nn.BatchNorm2d` in the full-alignment network.
pub(crate) const BATCH_NORM_EPS: f32 = 1e-3;

/// Narrower inputs go through im2col even at stride 1: the unit-stride products waste two
/// columns of every `width + 2`, and a narrow input's im2col rows are few.
const MIN_UNIT_STRIDE_WIDTH: usize = 8;

/// A batch norm in eval mode: `(x - mean) / sqrt(var + eps) * weight + bias`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct BatchNorm<'a> {
    pub(crate) weight: &'a [f32],
    pub(crate) bias: &'a [f32],
    pub(crate) mean: &'a [f32],
    pub(crate) var: &'a [f32],
    pub(crate) eps: f32,
}

/// What follows a convolution's bias.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Epilogue<'a> {
    /// ReLU.
    Relu,
    /// Adds this (the block input, laid out like the output), then ReLU.
    ResidualRelu(&'a [f32]),
}

/// An `nn.Conv2d` (3x3, padding 1, bias) followed by an eval-mode `nn.BatchNorm2d`, folded into
/// one convolution over `[height][width][channel]` inputs.
#[derive(Debug, Clone)]
pub(crate) struct Conv {
    pub(crate) in_channels: usize,
    pub(crate) out_channels: usize,
    pub(crate) stride: usize,
    /// `[k][out]` with `k = 9 * in` ordered as (row, column, channel), batch norm folded in.
    pub(crate) weight: Vec<f32>,
    /// Per output channel, the bias with the batch norm folded in.
    pub(crate) bias: Vec<f32>,
    /// Per output channel, the batch norm's scale folded into `weight`.
    #[cfg_attr(not(feature = "gpu"), expect(dead_code, reason = "only the GPU backend unfolds the batch norm"))]
    pub(crate) scale: Vec<f32>,
}

impl Conv {
    /// From PyTorch's `[out][in][3][3]` weights and `[out]` bias, followed by `bn`.
    ///
    /// # Panics
    ///
    /// If the slices do not match the channel counts.
    pub(crate) fn new(
        weight: &[f32],
        bias: &[f32],
        bn: BatchNorm<'_>,
        in_channels: usize,
        out_channels: usize,
        stride: usize,
    ) -> Self {
        assert_eq!(weight.len(), out_channels * in_channels * 9, "conv weight shape");
        assert!(
            [bias, bn.weight, bn.bias, bn.mean, bn.var].iter().all(|v| v.len() == out_channels),
            "conv bias and batch norm shapes"
        );
        let k = 9 * in_channels;
        let mut folded = vec![0.0; k * out_channels];
        let mut biases = Vec::with_capacity(out_channels);
        let mut scales = Vec::with_capacity(out_channels);
        for n in 0..out_channels {
            let bn_scale = bn.weight[n] / (bn.var[n] + bn.eps).sqrt();
            scales.push(bn_scale);
            biases.push((bias[n] - bn.mean[n]) * bn_scale + bn.bias[n]);
            let filter = &weight[n * k..(n + 1) * k];
            for c in 0..in_channels {
                for tap in 0..9 {
                    folded[(tap * in_channels + c) * out_channels + n] = filter[c * 9 + tap] * bn_scale;
                }
            }
        }
        Self { in_channels, out_channels, stride, weight: folded, bias: biases, scale: scales }
    }

    /// The convolution `{prefix}.weight`/`.bias` of a checkpoint with the batch norm `{bn}.*`.
    pub(crate) fn from_state(state: &mut Weights, conv: &str, bn: &str, stride: usize) -> Result<Self, Error> {
        let name = format!("{conv}.weight");
        let (shape, weight) = state.take(&name)?;
        let [out_channels, in_channels, kh, kw] = state.dims(&name, &shape)?;
        if (kh, kw) != (3, 3) {
            return Err(state.bad(format!("{name} has shape {shape:?}, expected a 3x3 kernel")));
        }
        let bias = state.take_vector(&format!("{conv}.bias"), out_channels)?;
        let mut take = |kind: &str| state.take_vector(&format!("{bn}.{kind}"), out_channels);
        let (gamma, beta, mean, var) = (take("weight")?, take("bias")?, take("running_mean")?, take("running_var")?);
        let norm = BatchNorm { weight: &gamma, bias: &beta, mean: &mean, var: &var, eps: BATCH_NORM_EPS };
        Ok(Self::new(&weight, &bias, norm, in_channels, out_channels, stride))
    }

    /// Output height and width for an input of `height` x `width`.
    pub(crate) fn output_size(&self, height: usize, width: usize) -> (usize, usize) {
        ((height - 1) / self.stride + 1, (width - 1) / self.stride + 1)
    }

    /// Convolves `images` inputs (`[image][height][width][in]`) into `output`
    /// (`[image][out_height][out_width][out]`), then applies `epilogue`.
    ///
    /// # Panics
    ///
    /// If `input` or a residual does not have the expected length.
    pub(crate) fn forward(
        &self,
        input: &[f32],
        (images, height, width): (usize, usize, usize),
        epilogue: Epilogue<'_>,
        output: &mut Vec<f32>,
        scratch: &mut ConvScratch,
    ) {
        assert_eq!(input.len(), images * height * width * self.in_channels, "conv input shape");
        let (out_height, out_width) = self.output_size(height, width);
        output.clear();
        output.resize(images * out_height * out_width * self.out_channels, 0.0);
        if self.stride == 1 && width >= MIN_UNIT_STRIDE_WIDTH {
            self.unit_stride(input, (height, width), output, scratch);
        } else {
            self.im2col(input, (images, height, width), output, &mut scratch.input);
        }
        let cout = self.out_channels;
        match epilogue {
            Epilogue::Relu => {
                for row in output.chunks_exact_mut(cout) {
                    for (y, b) in row.iter_mut().zip(&self.bias) {
                        *y = (*y + b).max(0.0);
                    }
                }
            }
            Epilogue::ResidualRelu(residual) => {
                assert_eq!(residual.len(), output.len(), "residual shape");
                for (row, residual) in output.chunks_exact_mut(cout).zip(residual.chunks_exact(cout)) {
                    for ((y, b), r) in row.iter_mut().zip(&self.bias).zip(residual) {
                        *y = (*y + b + r).max(0.0);
                    }
                }
            }
        }
    }

    /// Stride 1 without im2col: in a zero-padded `[height + 2][width + 2][in]` copy of an image,
    /// the three taps of a kernel row are `3 * in` adjacent values, so each kernel row is one
    /// product over an overlapping view of the padded image. The products cover the padded
    /// width; the two extra columns of each output row are dropped.
    fn unit_stride(
        &self,
        input: &[f32],
        (height, width): (usize, usize),
        output: &mut [f32],
        scratch: &mut ConvScratch,
    ) {
        let (cin, cout) = (self.in_channels, self.out_channels);
        let padded_width = width + 2;
        let rows = height * padded_width;
        // The last view (kernel row 2) reads two positions past the padded image's end.
        let positions = (height + 2) * padded_width + 2;
        let ConvScratch { input: padded, wide } = scratch;
        for (image, out) in input.chunks_exact(height * width * cin).zip(output.chunks_exact_mut(height * width * cout))
        {
            padded.clear();
            padded.resize(positions * cin, 0.0);
            for (y, row) in image.chunks_exact(width * cin).enumerate() {
                padded[((y + 1) * padded_width + 1) * cin..][..width * cin].copy_from_slice(row);
            }
            wide.clear();
            wide.resize(rows * cout, 0.0);
            for ky in 0..3 {
                let view = Strided { rows, cols: 3 * cin, stride: cin };
                let weights = &self.weight[ky * 3 * cin * cout..][..3 * cin * cout];
                let product = if ky == 0 { Product::Overwrite } else { Product::Accumulate };
                matmul(wide, Strided::dense(rows, cout), &padded[ky * padded_width * cin..], view, weights, product);
            }
            for (out_row, wide_row) in out.chunks_exact_mut(width * cout).zip(wide.chunks_exact(padded_width * cout)) {
                out_row.copy_from_slice(&wide_row[..width * cout]);
            }
        }
    }

    /// Any stride: one row of `9 * in` values per output position (im2col), then one product.
    fn im2col(
        &self,
        input: &[f32],
        (images, height, width): (usize, usize, usize),
        output: &mut [f32],
        columns: &mut Vec<f32>,
    ) {
        let (cin, cout) = (self.in_channels, self.out_channels);
        let (out_height, out_width) = self.output_size(height, width);
        let rows = images * out_height * out_width;
        let k = 9 * cin;
        columns.clear();
        columns.resize(rows * k, 0.0);
        let image_columns = columns.chunks_exact_mut(out_height * out_width * k);
        for (image, columns) in input.chunks_exact(height * width * cin).zip(image_columns) {
            for oy in 0..out_height {
                for ox in 0..out_width {
                    let row = &mut columns[(oy * out_width + ox) * k..][..k];
                    for ky in 0..3 {
                        let Some(iy) = (oy * self.stride + ky).checked_sub(1).filter(|&y| y < height) else { continue };
                        for kx in 0..3 {
                            let Some(ix) = (ox * self.stride + kx).checked_sub(1).filter(|&x| x < width) else {
                                continue;
                            };
                            row[(ky * 3 + kx) * cin..][..cin].copy_from_slice(&image[(iy * width + ix) * cin..][..cin]);
                        }
                    }
                }
            }
        }
        matmul(output, Strided::dense(rows, cout), columns, Strided::dense(rows, k), &self.weight, Product::Overwrite);
    }
}

/// Buffers a convolution reuses across calls.
#[derive(Debug, Default)]
pub(crate) struct ConvScratch {
    /// The im2col rows (strided convolutions) or the zero-padded image (stride 1).
    input: Vec<f32>,
    /// Stride 1: the output over the padded width.
    wide: Vec<f32>,
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A float reference: 3x3 convolution with padding 1, batch norm, ReLU.
    pub(crate) fn reference(conv: (&[f32], &[f32]), bn: BatchNorm<'_>, input: &[f32], dims: [usize; 5]) -> Vec<f32> {
        let [height, width, cin, cout, stride] = dims;
        let (weight, bias) = conv;
        let (oh, ow) = ((height - 1) / stride + 1, (width - 1) / stride + 1);
        let mut out = vec![0f32; oh * ow * cout];
        for oy in 0..oh {
            for ox in 0..ow {
                for n in 0..cout {
                    let mut sum = bias[n];
                    for c in 0..cin {
                        for ky in 0..3 {
                            for kx in 0..3 {
                                let (iy, ix) = ((oy * stride + ky) as isize - 1, (ox * stride + kx) as isize - 1);
                                if iy < 0 || ix < 0 || iy >= height as isize || ix >= width as isize {
                                    continue;
                                }
                                let x = input[(iy as usize * width + ix as usize) * cin + c];
                                sum += x * weight[((n * cin + c) * 3 + ky) * 3 + kx];
                            }
                        }
                    }
                    let y = (sum - bn.mean[n]) / (bn.var[n] + bn.eps).sqrt() * bn.weight[n] + bn.bias[n];
                    out[(oy * ow + ox) * cout + n] = y.max(0.0);
                }
            }
        }
        out
    }

    /// `scale * sin(0.37 i + phase)` for `i` in `0..n`.
    pub(crate) fn wave(n: usize, scale: f32, phase: f32) -> Vec<f32> {
        (0..n).map(|i| ((i as f32 * 0.37 + phase).sin()) * scale).collect()
    }

    /// Weights, bias and batch norm vectors for a `cin -> cout` convolution.
    pub(crate) fn conv_params(cin: usize, cout: usize) -> [Vec<f32>; 6] {
        let gamma = wave(cout, 0.2, 2.0).iter().map(|g| g + 1.0).collect();
        let var = wave(cout, 0.3, 5.0).iter().map(|v| v.abs() + 0.5).collect();
        [wave(cout * cin * 9, 0.3, 0.0), wave(cout, 0.1, 1.0), gamma, wave(cout, 0.1, 3.0), wave(cout, 0.1, 4.0), var]
    }

    #[test]
    fn matches_a_direct_convolution() {
        let (height, width, cin, cout) = (7, 5, 8, 32);
        let [weight, bias, gamma, beta, mean, var] = conv_params(cin, cout);
        let bn = BatchNorm { weight: &gamma, bias: &beta, mean: &mean, var: &var, eps: BATCH_NORM_EPS };
        let input: Vec<f32> = wave(height * width * cin, 1.0, 6.0).iter().map(|x| x.max(0.0)).collect();
        for stride in [1, 2] {
            let conv = Conv::new(&weight, &bias, bn, cin, cout, stride);
            let (mut out, mut columns) = (Vec::new(), ConvScratch::default());
            conv.forward(&input, (1, height, width), Epilogue::Relu, &mut out, &mut columns);
            let want = reference((&weight, &bias), bn, &input, [height, width, cin, cout, stride]);
            assert_eq!(out.len(), want.len());
            let error = out.iter().zip(&want).fold(0f32, |m, (a, b)| m.max((a - b).abs()));
            assert!(error < 1e-5, "stride {stride}: error {error}");
        }
    }

    #[test]
    fn convolves_several_images_as_each_alone() {
        let (images, height, width, cin, cout) = (3, 6, 5, 4, 8);
        let [weight, bias, gamma, beta, mean, var] = conv_params(cin, cout);
        let bn = BatchNorm { weight: &gamma, bias: &beta, mean: &mean, var: &var, eps: BATCH_NORM_EPS };
        let input = wave(images * height * width * cin, 1.0, 0.5);
        for stride in [1, 2] {
            let conv = Conv::new(&weight, &bias, bn, cin, cout, stride);
            let (mut together, mut columns) = (Vec::new(), ConvScratch::default());
            let residual =
                wave(images * conv.output_size(height, width).0 * conv.output_size(height, width).1 * cout, 1.0, 2.5);
            let epilogue = Epilogue::ResidualRelu(&residual);
            conv.forward(&input, (images, height, width), epilogue, &mut together, &mut columns);
            let mut alone = Vec::new();
            let residuals = residual.chunks_exact(residual.len() / images);
            for (image, residual) in input.chunks_exact(height * width * cin).zip(residuals) {
                let mut out = Vec::new();
                conv.forward(image, (1, height, width), Epilogue::ResidualRelu(residual), &mut out, &mut columns);
                alone.extend(out);
            }
            assert_eq!(together, alone, "stride {stride}");
        }
    }

    #[test]
    fn adds_the_residual_before_the_relu() {
        let (height, width, c) = (3, 4, 4);
        let [weight, bias, gamma, beta, mean, var] = conv_params(c, c);
        let bn = BatchNorm { weight: &gamma, bias: &beta, mean: &mean, var: &var, eps: BATCH_NORM_EPS };
        let conv = Conv::new(&weight, &bias, bn, c, c, 1);
        let input = wave(height * width * c, 1.0, 0.5);
        let residual = wave(height * width * c, 2.0, 1.5);
        let (mut plain, mut summed, mut columns) = (Vec::new(), Vec::new(), ConvScratch::default());
        // A residual of zero is a plain ReLU; otherwise relu(conv + residual).
        conv.forward(
            &input,
            (1, height, width),
            Epilogue::ResidualRelu(&vec![0.0; residual.len()]),
            &mut plain,
            &mut columns,
        );
        let mut relu = Vec::new();
        conv.forward(&input, (1, height, width), Epilogue::Relu, &mut relu, &mut columns);
        assert_eq!(plain, relu);
        conv.forward(&input, (1, height, width), Epilogue::ResidualRelu(&residual), &mut summed, &mut columns);
        // The pre-activation values, from a convolution whose bias is shifted so its ReLU never
        // clips.
        let mut linear = Vec::new();
        let shifted = Conv { bias: conv.bias.iter().map(|b| b + 1e3).collect(), ..conv.clone() };
        shifted.forward(&input, (1, height, width), Epilogue::Relu, &mut linear, &mut columns);
        for ((s, l), r) in summed.iter().zip(&linear).zip(&residual) {
            assert!((s - (l - 1e3 + r).max(0.0)).abs() < 1e-3, "{s} vs {l} + {r}");
        }
    }
}

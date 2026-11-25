//! 8-bit 3x3 convolutions and matrix products for CPUs with AVX-VNNI or AVX-512 VNNI (x86-64) or
//! the dot product instructions (aarch64), used by the full-alignment and pileup networks in
//! [`Precision::Mixed`](crate::Precision::Mixed).
//!
//! A convolution takes non-negative `f32` activations (every one in the network follows a ReLU),
//! scales each image's to `u8` by its largest value, multiplies them by `i8` weights scaled per
//! output channel (batch norm folded in) and sums the products exactly in `i32` with
//! `vpdpbusd`, or on aarch64 with `sdot`, which multiplies signed bytes: on the inputs less 128,
//! with 128 times each weight column's sum added back, the same sums. The sums are scaled back
//! to `f32` for the bias, the residual and the ReLU.
//!
//! The products read the quantized images in place (an implicit im2col): with the images stored
//! `[height][width][channel]` and a zero border around each, the three pixels a 3x3 window takes
//! from one input row are adjacent, so each output position's 9 * `in` inputs are three runs of
//! 3 * `in` bytes. The images of a batch go through each block of weights together.
//!
//! [`QuantMatrix`] multiplies activations in (-1, 1) (the LSTMs' outputs) by a weight matrix:
//! they are scaled to `-127..=127` and offset by 128 into `u8`, and the offset is taken back out
//! with each weight column's sum.

use crate::cpu::{Conv, Epilogue};

/// Output rows each microkernel call computes.
const MR: usize = 6;

/// Output channels each microkernel call computes: two registers of eight `i32`.
const NR: usize = 16;

/// The `vpdpbusd` the products use. Both give the same sums: the integer arithmetic is exact.
#[cfg(target_arch = "x86_64")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Vnni {
    /// AVX-512 VNNI (Cascade Lake, Ice Lake, Sapphire Rapids and later Xeons; AMD Zen 4 and
    /// later): a block's 16 output channels in one 512-bit register.
    Avx512,
    /// AVX-VNNI (Alder Lake and later Intel desktop CPUs): the 16 channels in two 256-bit
    /// registers.
    Avx,
}

/// The `vpdpbusd` this CPU has, with AVX2 for the rest of the kernels, if any.
#[cfg(target_arch = "x86_64")]
fn vnni() -> Option<Vnni> {
    if is_x86_feature_detected!("avx2") {
        if is_x86_feature_detected!("avx512vnni") {
            return Some(Vnni::Avx512);
        }
        if is_x86_feature_detected!("avxvnni") {
            return Some(Vnni::Avx);
        }
    }
    None
}

/// Whether this CPU has the instructions the kernels use.
pub(crate) fn supported() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        vnni().is_some()
    }
    #[cfg(target_arch = "aarch64")]
    {
        std::arch::is_aarch64_feature_detected!("dotprod")
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        false
    }
}

/// Buffers reused across calls.
#[derive(Debug, Default)]
pub(crate) struct Scratch {
    /// The quantized images, each `[height + 2][width + 2][in]` with a zero border.
    padded: Vec<u8>,
    /// Per output row (padded to a multiple of [`MR`]), where its window starts in `padded`.
    rows: Vec<usize>,
    /// Per image and output channel, the factor that turns its sums back into `f32`.
    factors: Vec<f32>,
}

/// A 3x3 convolution with padding 1 and a folded batch norm, with weights quantized to `i8`.
#[derive(Debug, Clone)]
pub(crate) struct QuantConv {
    in_channels: usize,
    pub(crate) out_channels: usize,
    stride: usize,
    /// `[out / NR][k / 4][NR][4]`, `k = 9 * in` ordered as (row, column, channel).
    packed: Vec<i8>,
    /// Per output channel, `128 * sum(weights)`, which `sdot` on the inputs less 128 leaves out.
    offset: Vec<i32>,
    /// Per output channel, the weight scale.
    scale: Vec<f32>,
    /// Per output channel, the bias with the batch norm folded in.
    bias: Vec<f32>,
}

impl QuantConv {
    /// `conv` with its weights quantized per output channel.
    ///
    /// # Panics
    ///
    /// If `conv`'s output channels are not a multiple of 16 or its input channels of 4.
    pub(crate) fn new(conv: &Conv) -> Self {
        let (in_channels, out_channels) = (conv.in_channels, conv.out_channels);
        assert!(out_channels.is_multiple_of(NR) && in_channels.is_multiple_of(4), "{in_channels} -> {out_channels}");
        let k = 9 * in_channels;
        let mut packed = vec![0i8; out_channels * k];
        let (mut offset, mut scales) = (Vec::with_capacity(out_channels), Vec::with_capacity(out_channels));
        for n in 0..out_channels {
            let filter = || conv.weight.iter().skip(n).step_by(out_channels);
            let largest = filter().fold(0f32, |m, &w| m.max(w.abs()));
            let scale = if largest > 0.0 { largest / 127.0 } else { 1.0 };
            scales.push(scale);
            let mut sum = 0i32;
            for (kk, &w) in filter().enumerate() {
                let at = ((n / NR * (k / 4) + kk / 4) * NR + n % NR) * 4 + kk % 4;
                packed[at] = (w / scale).round().clamp(-127.0, 127.0) as i8;
                sum += i32::from(packed[at]);
            }
            offset.push(128 * sum);
        }
        let bias = conv.bias.clone();
        Self { in_channels, out_channels, stride: conv.stride, packed, offset, scale: scales, bias }
    }

    /// Output height and width for an input of `height` x `width`.
    pub(crate) fn output_size(&self, height: usize, width: usize) -> (usize, usize) {
        ((height - 1) / self.stride + 1, (width - 1) / self.stride + 1)
    }

    /// Convolves `images` maps of `height` x `width` in `input` (`[image][height][width][in]`,
    /// non-negative, each quantized on its own scale) into `output`
    /// (`[image][out_height][out_width][out]`), then applies `epilogue`, whose residual has the
    /// output's layout.
    pub(crate) fn forward(
        &self,
        input: &[f32],
        (images, height, width): (usize, usize, usize),
        epilogue: Epilogue<'_>,
        output: &mut Vec<f32>,
        scratch: &mut Scratch,
    ) {
        let (cin, cout) = (self.in_channels, self.out_channels);
        let image_len = height * width * cin;
        assert_eq!(input.len(), images * image_len, "conv input shape");
        let (out_height, out_width) = self.output_size(height, width);
        let rows = images * out_height * out_width;
        let padded_width = width + 2;
        let padded_len = (height + 2) * padded_width * cin;

        scratch.padded.clear();
        scratch.padded.resize(images * padded_len, 0);
        scratch.factors.clear();
        for (i, image) in input.chunks_exact(image_len.max(1)).take(images).enumerate() {
            let largest = largest(image);
            let input_scale = if largest > 0.0 { largest / 255.0 } else { 1.0 };
            let padded = &mut scratch.padded[i * padded_len..][..padded_len];
            for (y, row) in image.chunks_exact(width * cin).enumerate() {
                quantize(row, input_scale, &mut padded[((y + 1) * padded_width + 1) * cin..][..width * cin]);
            }
            scratch.factors.extend(self.scale.iter().map(|&s| input_scale * s));
        }

        // Window starts, padded to a whole number of microkernel calls with the last row again.
        scratch.rows.clear();
        for i in 0..images {
            for oy in 0..out_height {
                for ox in 0..out_width {
                    let (y, x) = (oy * self.stride, ox * self.stride);
                    scratch.rows.push(i * padded_len + (y * padded_width + x) * cin);
                }
            }
        }
        let last = scratch.rows.last().copied().unwrap_or(0);
        scratch.rows.resize(rows.div_ceil(MR) * MR, last);

        output.clear();
        output.resize(rows * cout, 0.0);
        let windows = Windows { data: &scratch.padded, row_stride: padded_width * cin, run: 3 * cin, runs: 3 };
        let k = 9 * cin;
        let rows_per_image = out_height * out_width;
        let mut sums = [[0i32; NR]; MR];
        for block in 0..cout / NR {
            let weights = &self.packed[block * k * NR..][..k * NR];
            let offset = &self.offset[block * NR..][..NR];
            for (call, starts) in scratch.rows.as_chunks::<MR>().0.iter().enumerate() {
                multiply(&windows, starts, (weights, offset), &mut sums);
                let first = call * MR;
                for (r, row_sums) in sums.iter().enumerate().take(rows.saturating_sub(first).min(MR)) {
                    let row = first + r;
                    let at = row * cout + block * NR;
                    let factors = &scratch.factors[(row / rows_per_image) * cout + block * NR..][..NR];
                    let residual = match epilogue {
                        Epilogue::Relu => None,
                        Epilogue::ResidualRelu(residual) => Some(&residual[at..at + NR]),
                    };
                    finish(row_sums, factors, &self.bias[block * NR..][..NR], residual, &mut output[at..at + NR]);
                }
            }
        }
    }
}

/// A weight matrix (`[k][n]`, as [`crate::cpu::Dense`] stores it) quantized to `i8` per column,
/// for activations in (-1, 1).
#[derive(Debug, Clone)]
pub(crate) struct QuantMatrix {
    k: usize,
    n: usize,
    /// `[n / NR][k / 4][NR][4]`.
    packed: Vec<i8>,
    /// Per column, `128 * sum(column)`: what the inputs' offset adds to each sum.
    offset: Vec<i32>,
    /// Per column, the weight scale over the activations' (127).
    factor: Vec<f32>,
}

/// Activations in (-1, 1) are scaled by this to `-127..=127`.
const ACTIVATION_SCALE: f32 = 127.0;

impl QuantMatrix {
    /// `weight` (`k` rows of `n`) with each column quantized on its own scale.
    ///
    /// # Panics
    ///
    /// If `n` is not a multiple of 16, `k` of 4, or `weight` does not hold `k * n` values.
    pub(crate) fn new(weight: &[f32], k: usize, n: usize) -> Self {
        assert!(n.is_multiple_of(NR) && k.is_multiple_of(4) && weight.len() == k * n, "{k} x {n}");
        let mut packed = vec![0i8; k * n];
        let (mut offset, mut factor) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for j in 0..n {
            let column = || weight.iter().skip(j).step_by(n);
            let largest = column().fold(0f32, |m, &w| m.max(w.abs()));
            let scale = if largest > 0.0 { largest / 127.0 } else { 1.0 };
            let mut sum = 0i32;
            for (kk, &w) in column().enumerate() {
                let q = (w / scale).round().clamp(-127.0, 127.0) as i8;
                packed[((j / NR * (k / 4) + kk / 4) * NR + j % NR) * 4 + kk % 4] = q;
                sum += i32::from(q);
            }
            offset.push(128 * sum);
            factor.push(scale / ACTIVATION_SCALE);
        }
        Self { k, n, packed, offset, factor }
    }

    /// Rows of the matrix: the length of each input row.
    pub(crate) fn inputs(&self) -> usize {
        self.k
    }

    /// `out[r][j] += x[r] . column j` for `rows` rows of `k` values in `x` (`x_stride` apart,
    /// each in (-1, 1)), into rows of `out` `out_stride` apart.
    pub(crate) fn accumulate(
        &self,
        x: &[f32],
        (rows, x_stride): (usize, usize),
        out: &mut [f32],
        out_stride: usize,
        scratch: &mut Vec<u8>,
    ) {
        let k = self.k;
        if rows == 0 {
            return;
        }
        assert!(x.len() >= (rows - 1) * x_stride + k, "quantized matrix input");
        prepare_rows(scratch, rows, k);
        for r in 0..rows {
            quantize_signed(&x[r * x_stride..][..k], &mut scratch[r * k..][..k]);
        }
        self.multiply_rows(scratch, rows, out, out_stride);
    }

    /// As [`QuantMatrix::accumulate`], with the rows already quantized (by [`quantize_signed`]
    /// into a buffer from [`prepare_rows`]).
    pub(crate) fn multiply_rows(&self, quantized: &[u8], rows: usize, out: &mut [f32], out_stride: usize) {
        let (k, n) = (self.k, self.n);
        if rows == 0 {
            return;
        }
        assert!(
            quantized.len() >= rows.div_ceil(MR) * MR * k
                && out.len() >= (rows - 1) * out_stride + n
                && out_stride >= n,
            "quantized matrix operands"
        );
        let windows = Windows { data: quantized, row_stride: 0, run: k, runs: 1 };
        let mut sums = [[0i32; NR]; MR];
        for block in 0..n / NR {
            let weights = &self.packed[block * k * NR..][..k * NR];
            let (offset, factor) = (&self.offset[block * NR..][..NR], &self.factor[block * NR..][..NR]);
            for first in (0..rows).step_by(MR) {
                let starts = std::array::from_fn(|r| (first + r) * k);
                multiply(&windows, &starts, (weights, offset), &mut sums);
                for (r, row_sums) in sums.iter().enumerate().take((rows - first).min(MR)) {
                    let out = &mut out[(first + r) * out_stride + block * NR..][..NR];
                    for j in 0..NR {
                        out[j] += (row_sums[j] - offset[j]) as f32 * factor[j];
                    }
                }
            }
        }
    }
}

/// Makes `buffer` hold `rows` rows of `k` quantized values, padded to whole microkernel calls
/// with zeros' code.
pub(crate) fn prepare_rows(buffer: &mut Vec<u8>, rows: usize, k: usize) {
    buffer.clear();
    buffer.resize(rows.div_ceil(MR) * MR * k, 128);
}

/// `out = round(x * 127) + 128`, clamped to `1..=255`, for `x` in (-1, 1).
pub(crate) fn quantize_signed(x: &[f32], out: &mut [u8]) {
    #[cfg(target_arch = "x86_64")]
    if supported() {
        // SAFETY: the CPU has AVX2.
        unsafe { quantize_signed_avx2(x, out) };
        return;
    }
    for (o, &v) in out.iter_mut().zip(x) {
        *o = (v * ACTIVATION_SCALE + 128.0).round().clamp(1.0, 255.0) as u8;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn quantize_signed_avx2(x: &[f32], out: &mut [u8]) {
    use std::arch::x86_64::{
        _mm_packus_epi16, _mm_packus_epi32, _mm_storel_epi64, _mm256_add_ps, _mm256_castsi256_si128,
        _mm256_cvtps_epi32, _mm256_extracti128_si256, _mm256_loadu_ps, _mm256_max_ps, _mm256_min_ps, _mm256_mul_ps,
        _mm256_set1_ps,
    };
    let (chunks, rest) = x.as_chunks::<8>();
    let (scale, offset) = (_mm256_set1_ps(ACTIVATION_SCALE), _mm256_set1_ps(128.0));
    let (low, high) = (_mm256_set1_ps(1.0), _mm256_set1_ps(255.0));
    for (c, o) in chunks.iter().zip(out.as_chunks_mut::<8>().0) {
        // SAFETY: `c` holds eight values and `o` eight bytes.
        unsafe {
            let v = _mm256_add_ps(_mm256_mul_ps(_mm256_loadu_ps(c.as_ptr()), scale), offset);
            let i = _mm256_cvtps_epi32(_mm256_min_ps(_mm256_max_ps(v, low), high));
            let words = _mm_packus_epi32(_mm256_castsi256_si128(i), _mm256_extracti128_si256::<1>(i));
            _mm_storel_epi64(o.as_mut_ptr().cast(), _mm_packus_epi16(words, words));
        }
    }
    let done = chunks.len() * 8;
    for (o, &v) in out[done..].iter_mut().zip(rest) {
        *o = (v * ACTIVATION_SCALE + 128.0).round().clamp(1.0, 255.0) as u8;
    }
}

/// Quantized inputs and how each output row's inputs lie in them: `runs` runs of `run` bytes,
/// `row_stride` apart (three for a convolution's window, one for a matrix row).
struct Windows<'a> {
    data: &'a [u8],
    row_stride: usize,
    run: usize,
    runs: usize,
}

/// The largest value of `x`, or 0.
fn largest(x: &[f32]) -> f32 {
    #[cfg(target_arch = "x86_64")]
    if supported() {
        // SAFETY: the CPU has AVX2.
        return unsafe { largest_avx2(x) };
    }
    x.iter().fold(0f32, |m, &v| m.max(v))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn largest_avx2(x: &[f32]) -> f32 {
    use std::arch::x86_64::{_mm256_loadu_ps, _mm256_max_ps, _mm256_setzero_ps, _mm256_storeu_ps};
    let (chunks, rest) = x.as_chunks::<8>();
    let mut m = _mm256_setzero_ps();
    for c in chunks {
        // SAFETY: `c` holds eight values.
        m = _mm256_max_ps(m, unsafe { _mm256_loadu_ps(c.as_ptr()) });
    }
    let mut lanes = [0f32; 8];
    // SAFETY: `lanes` holds eight values.
    unsafe { _mm256_storeu_ps(lanes.as_mut_ptr(), m) };
    lanes.iter().chain(rest).fold(0f32, |m, &v| m.max(v))
}

/// `out = round(x / scale)` clamped to `0..=255`, rounding halves away from zero as `f32::round`.
fn quantize(x: &[f32], scale: f32, out: &mut [u8]) {
    assert_eq!(x.len(), out.len(), "quantize lengths");
    #[cfg(target_arch = "x86_64")]
    if supported() {
        // SAFETY: the CPU has AVX2.
        unsafe { quantize_avx2(x, scale, out) };
        return;
    }
    for (o, &v) in out.iter_mut().zip(x) {
        *o = (v / scale).round().clamp(0.0, 255.0) as u8;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn quantize_avx2(x: &[f32], scale: f32, out: &mut [u8]) {
    use std::arch::x86_64::{
        _CMP_GE_OQ, _MM_FROUND_TO_ZERO, _mm_packus_epi16, _mm_packus_epi32, _mm_storel_epi64, _mm256_add_ps,
        _mm256_and_ps, _mm256_castsi256_si128, _mm256_cmp_ps, _mm256_cvttps_epi32, _mm256_div_ps,
        _mm256_extracti128_si256, _mm256_loadu_ps, _mm256_max_ps, _mm256_min_ps, _mm256_round_ps, _mm256_set1_ps,
        _mm256_setzero_ps, _mm256_sub_ps,
    };
    let (chunks, rest) = x.as_chunks::<8>();
    let (scale_v, half, one) = (_mm256_set1_ps(scale), _mm256_set1_ps(0.5), _mm256_set1_ps(1.0));
    let (low, high) = (_mm256_setzero_ps(), _mm256_set1_ps(255.0));
    for (c, o) in chunks.iter().zip(out.as_chunks_mut::<8>().0) {
        // SAFETY: `c` holds eight values and `o` eight bytes.
        unsafe {
            let v = _mm256_div_ps(_mm256_loadu_ps(c.as_ptr()), scale_v);
            // Half away from zero, as `f32::round` for the non-negative values here: the integer
            // part, plus one when the fraction (exact in `f32`) reaches a half.
            let whole = _mm256_round_ps::<_MM_FROUND_TO_ZERO>(v);
            let up = _mm256_and_ps(_mm256_cmp_ps::<_CMP_GE_OQ>(_mm256_sub_ps(v, whole), half), one);
            let r = _mm256_min_ps(_mm256_max_ps(_mm256_add_ps(whole, up), low), high);
            let i = _mm256_cvttps_epi32(r);
            let words = _mm_packus_epi32(_mm256_castsi256_si128(i), _mm256_extracti128_si256::<1>(i));
            _mm_storel_epi64(o.as_mut_ptr().cast(), _mm_packus_epi16(words, words));
        }
    }
    let done = chunks.len() * 8;
    for (o, &v) in out[done..].iter_mut().zip(rest) {
        *o = (v / scale).round().clamp(0.0, 255.0) as u8;
    }
}

/// `out = max(sums * factors + bias (+ residual), 0)` for one row of [`NR`] channels, each step
/// rounded as the scalar expression would be.
fn finish(sums: &[i32; NR], factors: &[f32], bias: &[f32], residual: Option<&[f32]>, out: &mut [f32]) {
    #[cfg(target_arch = "x86_64")]
    if supported() {
        // SAFETY: the CPU has AVX2.
        unsafe { finish_avx2(sums, factors, bias, residual, out) };
        return;
    }
    for j in 0..NR {
        let mut y = sums[j] as f32 * factors[j] + bias[j];
        if let Some(residual) = residual {
            y += residual[j];
        }
        out[j] = y.max(0.0);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn finish_avx2(sums: &[i32; NR], factors: &[f32], bias: &[f32], residual: Option<&[f32]>, out: &mut [f32]) {
    use std::arch::x86_64::{
        __m256i, _mm256_add_ps, _mm256_cvtepi32_ps, _mm256_loadu_ps, _mm256_loadu_si256, _mm256_max_ps, _mm256_mul_ps,
        _mm256_setzero_ps, _mm256_storeu_ps,
    };
    assert!(factors.len() >= NR && bias.len() >= NR && out.len() >= NR && residual.is_none_or(|r| r.len() >= NR));
    for half in [0, 8] {
        // SAFETY: the assert above bounds every eight-value access.
        unsafe {
            let s = _mm256_cvtepi32_ps(_mm256_loadu_si256(sums.as_ptr().add(half).cast::<__m256i>()));
            // A product, then a sum (no fused multiply-add), as the scalar code rounds.
            let mut y = _mm256_add_ps(
                _mm256_mul_ps(s, _mm256_loadu_ps(factors.as_ptr().add(half))),
                _mm256_loadu_ps(bias.as_ptr().add(half)),
            );
            if let Some(residual) = residual {
                y = _mm256_add_ps(y, _mm256_loadu_ps(residual.as_ptr().add(half)));
            }
            _mm256_storeu_ps(out.as_mut_ptr().add(half), _mm256_max_ps(y, _mm256_setzero_ps()));
        }
    }
}

/// `sums[r][j]` = row `starts[r]` . weight column `j`, over the row's runs. `offset[j]` is 128
/// times column `j`'s sum.
fn multiply(
    windows: &Windows<'_>,
    starts: &[usize; MR],
    (weights, offset): (&[i8], &[i32]),
    sums: &mut [[i32; NR]; MR],
) {
    let Windows { data, row_stride, run, runs } = *windows;
    let span = (runs - 1) * row_stride + run;
    assert!(
        runs >= 1
            && run.is_multiple_of(4)
            && weights.len() >= runs * run * NR
            && offset.len() >= NR
            && starts.iter().all(|&s| s + span <= data.len()),
        "kernel operands"
    );
    #[cfg(target_arch = "x86_64")]
    match vnni() {
        // SAFETY: the CPU has AVX-512 VNNI, and the assert above bounds every read.
        Some(Vnni::Avx512) => return unsafe { multiply_vnni512(windows, starts, weights, sums) },
        // SAFETY: the CPU has AVX2 and AVX-VNNI, and the assert above bounds every read.
        Some(Vnni::Avx) => return unsafe { multiply_vnni(windows, starts, weights, sums) },
        None => {}
    }
    #[cfg(target_arch = "aarch64")]
    if supported() {
        // SAFETY: the CPU has the dot product instructions, and the assert above bounds every read.
        unsafe { multiply_sdot(windows, starts, (weights, offset), sums) };
        return;
    }
    let _ = offset;
    for (r, row) in sums.iter_mut().enumerate() {
        for (j, sum) in row.iter_mut().enumerate() {
            *sum = (0..runs * run)
                .map(|kk| {
                    let input = data[starts[r] + (kk / run) * row_stride + kk % run];
                    i32::from(input) * i32::from(weights[((kk / 4) * NR + j) * 4 + kk % 4])
                })
                .sum();
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,avxvnni")]
unsafe fn multiply_vnni(windows: &Windows<'_>, starts: &[usize; MR], weights: &[i8], sums: &mut [[i32; NR]; MR]) {
    use std::arch::x86_64::{
        __m256i, _mm256_dpbusd_avx_epi32, _mm256_loadu_si256, _mm256_set1_epi32, _mm256_setzero_si256,
        _mm256_storeu_si256,
    };
    let Windows { data, row_stride, run, runs } = *windows;
    let mut acc = [[_mm256_setzero_si256(); 2]; MR];
    for ky in 0..runs {
        // Pointers to each row's run `ky` and its weights; `multiply` checked that every row's
        // runs lie in `data` and that `weights` holds `runs * run * NR` bytes.
        let a = starts.map(|s| data[s + ky * row_stride..].as_ptr());
        let b = weights[ky * run * NR..].as_ptr();
        for step in 0..run / 4 {
            // SAFETY: as above, every read is inside run `ky` of a row and its weights.
            unsafe {
                let b0 = _mm256_loadu_si256(b.add(step * 4 * NR).cast::<__m256i>());
                let b1 = _mm256_loadu_si256(b.add(step * 4 * NR + 32).cast::<__m256i>());
                for (row, a) in acc.iter_mut().zip(&a) {
                    let av = _mm256_set1_epi32(a.add(step * 4).cast::<i32>().read_unaligned());
                    row[0] = _mm256_dpbusd_avx_epi32(row[0], av, b0);
                    row[1] = _mm256_dpbusd_avx_epi32(row[1], av, b1);
                }
            }
        }
    }
    for (row, out) in acc.iter().zip(sums.iter_mut()) {
        // SAFETY: each half of `out` holds eight `i32`.
        unsafe {
            _mm256_storeu_si256(out.as_mut_ptr().cast::<__m256i>(), row[0]);
            _mm256_storeu_si256(out.as_mut_ptr().add(8).cast::<__m256i>(), row[1]);
        }
    }
}

/// [`multiply_vnni`] with AVX-512: a step's weights for the block's [`NR`] channels are one
/// 64-byte register, so each row takes one `vpdpbusd` per step instead of two.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512vnni")]
unsafe fn multiply_vnni512(windows: &Windows<'_>, starts: &[usize; MR], weights: &[i8], sums: &mut [[i32; NR]; MR]) {
    use std::arch::x86_64::{
        _mm512_dpbusd_epi32, _mm512_loadu_si512, _mm512_set1_epi32, _mm512_setzero_si512, _mm512_storeu_si512,
    };
    let Windows { data, row_stride, run, runs } = *windows;
    let mut acc = [_mm512_setzero_si512(); MR];
    for ky in 0..runs {
        // As in `multiply_vnni`: `multiply` checked every read.
        let a = starts.map(|s| data[s + ky * row_stride..].as_ptr());
        let b = weights[ky * run * NR..].as_ptr();
        for step in 0..run / 4 {
            // SAFETY: as above, every read is inside run `ky` of a row and its weights.
            unsafe {
                let w = _mm512_loadu_si512(b.add(step * 4 * NR).cast());
                for (row, a) in acc.iter_mut().zip(&a) {
                    let av = _mm512_set1_epi32(a.add(step * 4).cast::<i32>().read_unaligned());
                    *row = _mm512_dpbusd_epi32(*row, av, w);
                }
            }
        }
    }
    for (row, out) in acc.iter().zip(sums.iter_mut()) {
        // SAFETY: `out` holds sixteen `i32`.
        unsafe { _mm512_storeu_si512(out.as_mut_ptr().cast(), *row) };
    }
}

/// [`multiply`] with `sdot`, which multiplies signed bytes: on the inputs less 128 (each byte's
/// top bit flipped), plus `offset`, 128 times each column's sum, the same sums.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon,dotprod")]
unsafe fn multiply_sdot(
    windows: &Windows<'_>,
    starts: &[usize; MR],
    (weights, offset): (&[i8], &[i32]),
    sums: &mut [[i32; NR]; MR],
) {
    use std::arch::aarch64::{
        int8x16_t, int32x4_t, vaddq_s32, vdupq_n_s32, vdupq_n_u32, vld1q_s8, vld1q_s32, vreinterpretq_s8_u32, vst1q_s32,
    };

    /// `acc` plus the dot products of each four bytes of `a` and `b` (`sdot`, here as assembly:
    /// its intrinsic needs Rust 1.98).
    #[inline(always)]
    fn sdot(acc: int32x4_t, a: int8x16_t, b: int8x16_t) -> int32x4_t {
        let mut acc = acc;
        // SAFETY: `sdot` only reads and writes these registers; the caller has the dot product
        // instructions.
        unsafe {
            std::arch::asm!(
                "sdot {acc:v}.4s, {a:v}.16b, {b:v}.16b",
                acc = inout(vreg) acc, a = in(vreg) a, b = in(vreg) b,
                options(pure, nomem, nostack, preserves_flags),
            );
        }
        acc
    }

    let Windows { data, row_stride, run, runs } = *windows;
    let mut acc = [[vdupq_n_s32(0); NR / 4]; MR];
    for ky in 0..runs {
        // Pointers to each row's run `ky` and its weights; `multiply` checked that every row's
        // runs lie in `data` and that `weights` holds `runs * run * NR` bytes.
        let a = starts.map(|s| data[s + ky * row_stride..].as_ptr());
        let b = weights[ky * run * NR..].as_ptr();
        for step in 0..run / 4 {
            // SAFETY: as above, every read is inside run `ky` of a row and its weights.
            unsafe {
                let bv = [0, 1, 2, 3].map(|q| vld1q_s8(b.add(step * 4 * NR + q * 16)));
                for (row, a) in acc.iter_mut().zip(&a) {
                    let four = a.add(step * 4).cast::<u32>().read_unaligned() ^ 0x8080_8080;
                    let av = vreinterpretq_s8_u32(vdupq_n_u32(four));
                    for (lane, &w) in row.iter_mut().zip(&bv) {
                        *lane = sdot(*lane, av, w);
                    }
                }
            }
        }
    }
    // SAFETY: `multiply` checked that `offset` holds NR `i32`.
    let offset = [0, 1, 2, 3].map(|q| unsafe { vld1q_s32(offset.as_ptr().add(q * 4)) });
    for (row, out) in acc.iter().zip(sums.iter_mut()) {
        for (q, (&lane, &o)) in row.iter().zip(&offset).enumerate() {
            // SAFETY: each quarter of `out` holds four `i32`.
            unsafe { vst1q_s32(out.as_mut_ptr().add(q * 4), vaddq_s32(lane, o)) };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpu::BatchNorm;
    use crate::cpu::conv::tests::{conv_params, reference, wave};

    /// The convolution as it was computed before the implicit im2col: one image at a time, an
    /// explicit im2col and scalar arithmetic. The new code must give the same bits.
    fn explicit_im2col(
        conv: &QuantConv,
        input: &[f32],
        (height, width): (usize, usize),
        residual: Option<&[f32]>,
    ) -> Vec<f32> {
        let (cin, cout, k) = (conv.in_channels, conv.out_channels, 9 * conv.in_channels);
        let (out_height, out_width) = conv.output_size(height, width);
        let largest = input.iter().fold(0f32, |m, &x| m.max(x));
        let input_scale = if largest > 0.0 { largest / 255.0 } else { 1.0 };
        let quantized: Vec<u8> = input.iter().map(|&x| (x / input_scale).round().clamp(0.0, 255.0) as u8).collect();
        let mut out = vec![0f32; out_height * out_width * cout];
        for oy in 0..out_height {
            for ox in 0..out_width {
                let mut column = vec![0u8; k];
                for ky in 0..3 {
                    let Some(iy) = (oy * conv.stride + ky).checked_sub(1).filter(|&y| y < height) else { continue };
                    for kx in 0..3 {
                        let Some(ix) = (ox * conv.stride + kx).checked_sub(1).filter(|&x| x < width) else { continue };
                        column[(ky * 3 + kx) * cin..][..cin]
                            .copy_from_slice(&quantized[(iy * width + ix) * cin..][..cin]);
                    }
                }
                for n in 0..cout {
                    let sum: i32 = (0..k)
                        .map(|kk| {
                            let w = conv.packed[((n / NR * (k / 4) + kk / 4) * NR + n % NR) * 4 + kk % 4];
                            i32::from(column[kk]) * i32::from(w)
                        })
                        .sum();
                    let at = (oy * out_width + ox) * cout + n;
                    let mut y = sum as f32 * (input_scale * conv.scale[n]) + conv.bias[n];
                    if let Some(residual) = residual {
                        y += residual[at];
                    }
                    out[at] = y.max(0.0);
                }
            }
        }
        out
    }

    #[test]
    fn matches_a_float_convolution_closely() {
        let (height, width, cin, cout) = (7, 5, 8, 32);
        let [weight, bias, gamma, beta, mean, var] = conv_params(cin, cout);
        let bn = BatchNorm { weight: &gamma, bias: &beta, mean: &mean, var: &var, eps: 1e-3 };
        let input: Vec<f32> = wave(height * width * cin, 1.0, 6.0).iter().map(|x| x.max(0.0)).collect();
        for stride in [1, 2] {
            let conv = QuantConv::new(&Conv::new(&weight, &bias, bn, cin, cout, stride));
            let (mut out, mut scratch) = (Vec::new(), Scratch::default());
            conv.forward(&input, (1, height, width), Epilogue::Relu, &mut out, &mut scratch);
            let want = reference((&weight, &bias), bn, &input, [height, width, cin, cout, stride]);
            assert_eq!(out.len(), want.len());
            let largest = want.iter().fold(0f32, |m, &y| m.max(y));
            let error = out.iter().zip(&want).fold(0f32, |m, (a, b)| m.max((a - b).abs()));
            assert!(error < 0.02 * largest, "stride {stride}: error {error} of {largest}");
        }
    }

    #[test]
    fn batches_give_each_images_bits() {
        // Three images on different scales (one all zero), odd sizes, with and without a
        // residual: each image's output must be what it was alone, with an explicit im2col.
        let (images, height, width, cin, cout) = (3, 9, 7, 8, 32);
        let [weight, bias, gamma, beta, mean, var] = conv_params(cin, cout);
        let bn = BatchNorm { weight: &gamma, bias: &beta, mean: &mean, var: &var, eps: 1e-3 };
        let image_len = height * width * cin;
        let mut input: Vec<f32> = wave(images * image_len, 1.0, 6.0).iter().map(|x| x.max(0.0)).collect();
        input[image_len..2 * image_len].fill(0.0);
        for (i, v) in input[2 * image_len..].iter_mut().enumerate() {
            *v *= 37.5 + (i % 3) as f32;
        }
        for stride in [1, 2] {
            let conv = QuantConv::new(&Conv::new(&weight, &bias, bn, cin, cout, stride));
            let (out_height, out_width) = conv.output_size(height, width);
            let out_len = out_height * out_width * cout;
            let residual: Vec<f32> = wave(images * out_len, 0.5, 2.0);
            for with_residual in [false, true] {
                let epilogue = if with_residual { Epilogue::ResidualRelu(&residual) } else { Epilogue::Relu };
                let (mut out, mut scratch) = (Vec::new(), Scratch::default());
                conv.forward(&input, (images, height, width), epilogue, &mut out, &mut scratch);
                for i in 0..images {
                    let image = &input[i * image_len..][..image_len];
                    let residual = with_residual.then(|| &residual[i * out_len..][..out_len]);
                    let want = explicit_im2col(&conv, image, (height, width), residual);
                    let got = &out[i * out_len..][..out_len];
                    let same = got.iter().zip(&want).all(|(a, b)| a.to_bits() == b.to_bits());
                    assert!(same, "stride {stride}, image {i}, residual {with_residual}");
                }
            }
        }
    }

    #[test]
    fn quantizing_rounds_halves_away_from_zero() {
        // Exact halves, values just below a half, the largest value and zeros.
        let x: Vec<f32> = [0.5, 1.5, 2.5, 0.499_999_97, 254.5, 255.0, 0.0, 3.25, 7.75, 128.5, 1e-7]
            .iter()
            .cycle()
            .take(43)
            .copied()
            .collect();
        let mut got = vec![0u8; x.len()];
        quantize(&x, 1.0, &mut got);
        let want: Vec<u8> = x.iter().map(|&v| v.round().clamp(0.0, 255.0) as u8).collect();
        assert_eq!(got, want);
    }

    #[test]
    fn quantized_matrices_match_float_products_closely() {
        // 13 rows (two microkernel calls and a part) of 40 inputs in (-1, 1), strided, by a
        // 40 x 32 matrix, added to what `out` holds.
        let (rows, k, n, stride, out_stride) = (13, 40, 32, 45, 35);
        let x: Vec<f32> = wave(rows * stride, 0.9, 1.0);
        let weight: Vec<f32> = wave(k * n, 0.3, 3.0);
        let mut out: Vec<f32> = wave(rows * out_stride, 1.0, 5.0);
        let mut want = out.clone();
        for r in 0..rows {
            for j in 0..n {
                want[r * out_stride + j] += (0..k).map(|kk| x[r * stride + kk] * weight[kk * n + j]).sum::<f32>();
            }
        }
        let matrix = QuantMatrix::new(&weight, k, n);
        matrix.accumulate(&x, (rows, stride), &mut out, out_stride, &mut Vec::new());
        let error = out.iter().zip(&want).fold(0f32, |m, (a, b)| m.max((a - b).abs()));
        assert!(error < 0.05, "largest error {error}");
        // Columns past `n` are left alone.
        assert_eq!(out[n], want[n]);
    }

    /// Each kernel this CPU can run, the chosen one (through `multiply`) first. On a CPU without
    /// AVX-512 VNNI its kernel is not tested; Intel's SDE (`sde64 -icx -- cargo test`) runs it.
    fn kernels() -> Vec<(&'static str, Kernel)> {
        let mut kernels: Vec<(&'static str, Kernel)> = vec![("chosen", multiply)];
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("avxvnni") {
                // SAFETY: the CPU has AVX2 and AVX-VNNI; the test's operands are in bounds.
                kernels.push(("avx-vnni", |w, s, (b, _), out| unsafe { multiply_vnni(w, s, b, out) }));
            }
            if is_x86_feature_detected!("avx512vnni") {
                // SAFETY: the CPU has AVX-512 VNNI; the test's operands are in bounds.
                kernels.push(("avx512-vnni", |w, s, (b, _), out| unsafe { multiply_vnni512(w, s, b, out) }));
            }
        }
        #[cfg(target_arch = "aarch64")]
        if std::arch::is_aarch64_feature_detected!("dotprod") {
            // SAFETY: the CPU has the dot product instructions; the test's operands are in bounds.
            kernels.push(("sdot", |w, s, b, out| unsafe { multiply_sdot(w, s, b, out) }));
        }
        kernels
    }

    type Kernel = fn(&Windows<'_>, &[usize; MR], (&[i8], &[i32]), &mut [[i32; NR]; MR]);

    #[test]
    fn the_kernels_match_scalar_sums() {
        let (row_stride, run) = (40, 24);
        let data: Vec<u8> = (0..2 * row_stride + run + 5 * 8).map(|i| (i * 7 % 256) as u8).collect();
        let weights: Vec<i8> = (0..3 * run * NR).map(|i| ((i * 13 % 255) as i32 - 127) as i8).collect();
        let starts = [0, 8, 16, 24, 32, 40];
        let windows = Windows { data: &data, row_stride, run, runs: 3 };
        let offset: [i32; NR] = std::array::from_fn(|j| {
            128 * (0..3 * run).map(|kk| i32::from(weights[((kk / 4) * NR + j) * 4 + kk % 4])).sum::<i32>()
        });
        for (name, kernel) in kernels() {
            let mut sums = [[0i32; NR]; MR];
            kernel(&windows, &starts, (&weights, &offset), &mut sums);
            for (r, row) in sums.iter().enumerate() {
                for (j, &sum) in row.iter().enumerate() {
                    let want: i32 = (0..3 * run)
                        .map(|kk| {
                            let input = data[starts[r] + (kk / run) * row_stride + kk % run];
                            i32::from(input) * i32::from(weights[((kk / 4) * NR + j) * 4 + kk % 4])
                        })
                        .sum();
                    assert_eq!(sum, want, "{name}: row {r} column {j}");
                }
            }
        }
    }
}

//! Rust Metal backend: the networks of [`crate::cpu`] as Metal kernels, for Apple GPUs.
//!
//! The kernels (`metal/kernels.metal`) are compiled when the first model loads, by the Metal
//! framework, which `objc2-metal` calls from Rust; building needs no Xcode. Weights come from the
//! CPU models, already transposed, with the batch norms folded into the convolutions as on the
//! CPU. Products are register-tiled on the GPU's cores as CUDA's are (convolutions as implicit
//! im2col, long sums split across threadgroups), not run on the simdgroup matrix units: on an M4
//! Pro those ran `float32` at half the rate of plain multiply-adds. Each LSTM step is the
//! recurrent product of each direction, added to the gates, then one kernel for both directions'
//! cell updates. A batch is one command buffer (which keeps the buffers it binds alive), and its
//! inputs and outputs go through memory the CPU and GPU share, without copies.
//!
//! The arithmetic is `float32` throughout, with precise math: Apple GPUs have no TF32, so the
//! probabilities follow the CPU's [`Precision::Float32`](crate::Precision::Float32) path, the
//! sums in another order. There is no faster arithmetic on the GPU:
//! [`Precision::Mixed`](crate::Precision::Mixed) is the CPU's.

mod full_alignment;
mod pileup;
#[cfg(test)]
mod tests;

use std::ffi::c_void;
use std::ptr::NonNull;

use objc2::available;
use objc2::rc::Retained;
use objc2::runtime::ProtocolObject;
use objc2_foundation::NSString;
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder, MTLCommandQueue, MTLCompileOptions,
    MTLComputeCommandEncoder, MTLComputePipelineState, MTLCreateSystemDefaultDevice, MTLDevice, MTLLibrary,
    MTLMathMode, MTLResourceOptions, MTLSize,
};

use crate::Error;
use crate::cpu::{Conv, Dense};

pub use full_alignment::MetalFullAlignmentModel;
pub use pileup::MetalPileupModel;

const SOURCE: &str = include_str!("metal/kernels.metal");

// `MTLCreateSystemDefaultDevice` finds no GPU in a command-line program unless CoreGraphics is
// loaded.
#[link(name = "CoreGraphics", kind = "framework")]
unsafe extern "C" {}

/// Sums at least this long are split across threadgroups when a product has too few tiles.
const SPLIT_K_MIN: usize = 1024;

/// Threadgroups that keep the GPU busy: four per core of a 16-core M4 Pro.
const BUSY_GROUPS: usize = 64;

/// Threads per threadgroup of every kernel.
const THREADS: usize = 256;

/// What follows a product's bias (and residual): the `ACT_*` codes of the kernels.
#[derive(Debug, Clone, Copy)]
enum Activation {
    None = 0,
    Relu = 1,
    Selu = 2,
}

type Device = ProtocolObject<dyn MTLDevice>;
type Pipeline = Retained<ProtocolObject<dyn MTLComputePipelineState>>;
type MetalBuffer = Retained<ProtocolObject<dyn MTLBuffer>>;

fn backend(err: impl std::fmt::Display) -> Error {
    Error::Backend { message: format!("Metal: {err}") }
}

/// The kernels of `kernels.metal`.
#[derive(Debug)]
struct Kernels {
    gemm_128x64: Pipeline,
    gemm_64x32: Pipeline,
    conv_128x64: Pipeline,
    splitk_reduce: Pipeline,
    lstm_cell: Pipeline,
    int8_to_float: Pipeline,
    int32_to_float: Pipeline,
    pyramid_pool: Pipeline,
    softmax_rows: Pipeline,
}

/// The system's GPU, its command queue and the compiled kernels.
#[derive(Debug)]
pub(crate) struct Metal {
    device: Retained<Device>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    kernels: Kernels,
    /// Bound where a kernel takes no bias or residual.
    unused: MetalBuffer,
    /// Partial sums of split-K products.
    workspace: std::cell::RefCell<Buffer>,
}

impl Metal {
    /// Opens the system's GPU and compiles the kernels for it.
    fn new() -> Result<Self, Error> {
        let device = MTLCreateSystemDefaultDevice().ok_or_else(|| Error::Device { device: "metal".to_owned() })?;
        let library = library(&device)?;
        let f = |name: &str| -> Result<Pipeline, Error> {
            let function = library
                .newFunctionWithName(&NSString::from_str(name))
                .ok_or_else(|| backend(format!("no kernel {name}")))?;
            device.newComputePipelineStateWithFunction_error(&function).map_err(|e| backend(e.localizedDescription()))
        };
        let kernels = Kernels {
            gemm_128x64: f("gemm_128x64")?,
            gemm_64x32: f("gemm_64x32")?,
            conv_128x64: f("conv_128x64")?,
            splitk_reduce: f("splitk_reduce")?,
            lstm_cell: f("lstm_cell")?,
            int8_to_float: f("int8_to_float")?,
            int32_to_float: f("int32_to_float")?,
            pyramid_pool: f("pyramid_pool")?,
            softmax_rows: f("softmax_rows")?,
        };
        let products = [&kernels.gemm_128x64, &kernels.gemm_64x32, &kernels.conv_128x64];
        if products.iter().any(|p| p.maxTotalThreadsPerThreadgroup() < THREADS) {
            return Err(backend(format!("the products cannot run {THREADS} threads per threadgroup")));
        }
        let queue = device.newCommandQueue().ok_or_else(|| backend("no command queue"))?;
        let unused = new_buffer(&device, 16)?;
        Ok(Self { device, queue, kernels, unused, workspace: Default::default() })
    }

    /// A buffer holding `values`.
    fn upload(&self, values: &[f32]) -> Result<MetalBuffer, Error> {
        let buffer = new_buffer(&self.device, size_of_val(values))?;
        // SAFETY: the buffer is shared memory of at least `values`' size, and nothing else uses
        // it yet.
        unsafe {
            std::ptr::copy_nonoverlapping(values.as_ptr(), buffer.contents().cast::<f32>().as_ptr(), values.len())
        };
        Ok(buffer)
    }

    /// A command buffer and encoder for one batch.
    fn pass(&self) -> Result<Pass<'_>, Error> {
        let commands = self.queue.commandBuffer().ok_or_else(|| backend("no command buffer"))?;
        let encoder = commands.computeCommandEncoder().ok_or_else(|| backend("no compute encoder"))?;
        Ok(Pass { metal: self, commands, encoder, ended: false })
    }
}

/// The kernels compiled for `device` (the framework caches the compiled library across loads).
fn library(device: &Device) -> Result<Retained<ProtocolObject<dyn MTLLibrary>>, Error> {
    let options = MTLCompileOptions::new();
    // Precise math: the networks' exp, log and tanh as on the CPU, and no reassociation.
    // mathMode is macOS 15's; calling it on an earlier macOS raises an Objective-C exception,
    // so there fast math is turned off with the property it replaced.
    if available!(macos = 15.0) {
        options.setMathMode(MTLMathMode::Safe);
    } else {
        #[allow(deprecated, reason = "its replacement, mathMode, is what macOS 14 and earlier lack")]
        options.setFastMathEnabled(false);
    }
    device
        .newLibraryWithSource_options_error(&NSString::from_str(SOURCE), Some(&options))
        .map_err(|e| Error::Backend { message: format!("compiling the Metal kernels: {}", e.localizedDescription()) })
}

/// A buffer of `bytes` in memory the CPU and GPU share.
fn new_buffer(device: &Device, bytes: usize) -> Result<MetalBuffer, Error> {
    device
        .newBufferWithLength_options(bytes.max(16), MTLResourceOptions::StorageModeShared)
        .ok_or_else(|| backend(format!("allocating {bytes} bytes")))
}

/// A buffer that grows to the largest size asked for.
#[derive(Debug, Default)]
struct Buffer {
    buffer: Option<MetalBuffer>,
}

impl Buffer {
    /// The buffer, with room for at least `len` values of `T`.
    fn view<T>(&mut self, metal: &Metal, len: usize) -> Result<View, Error> {
        let bytes = len * size_of::<T>();
        if self.buffer.as_ref().is_none_or(|b| b.length() < bytes) {
            self.buffer = None;
            self.buffer = Some(new_buffer(&metal.device, bytes)?);
        }
        Ok(View { buffer: self.buffer.clone().expect("allocated above"), offset: 0 })
    }
}

/// A place in a buffer, as the kernels take their pointers.
#[derive(Debug, Clone)]
struct View {
    buffer: MetalBuffer,
    /// In bytes.
    offset: usize,
}

impl View {
    fn of(buffer: &MetalBuffer) -> Self {
        Self { buffer: buffer.clone(), offset: 0 }
    }

    /// `floats` values further on.
    fn at(&self, floats: usize) -> Self {
        Self { buffer: self.buffer.clone(), offset: self.offset + floats * size_of::<f32>() }
    }

    /// Copies `values` in. The GPU must not be using the buffer.
    fn write<T: Copy>(&self, values: &[T]) {
        assert!(self.offset + size_of_val(values) <= self.buffer.length(), "write past the buffer");
        // SAFETY: shared memory, in bounds (asserted), not in use by the GPU (the caller's
        // promise: each batch waits for its command buffer).
        unsafe {
            let dst = self.buffer.contents().as_ptr().cast::<u8>().add(self.offset).cast::<T>();
            std::ptr::copy_nonoverlapping(values.as_ptr(), dst, values.len());
        }
    }

    /// Copies `len` floats out. The GPU must be done with the buffer.
    fn read(&self, len: usize) -> Vec<f32> {
        assert!(self.offset + len * size_of::<f32>() <= self.buffer.length(), "read past the buffer");
        // SAFETY: as in `write`; the command buffer that wrote the values has completed.
        unsafe {
            let src = self.buffer.contents().as_ptr().cast::<u8>().add(self.offset).cast::<f32>();
            std::slice::from_raw_parts(src, len).to_vec()
        }
    }
}

/// `ProductParams` of the kernels.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
struct ProductParams {
    m: i32,
    n: i32,
    k: i32,
    lda: i32,
    ldc: i32,
    ldr: i32,
    act: i32,
    has_bias: i32,
    has_residual: i32,
    h: i32,
    w: i32,
    c: i32,
    oh: i32,
    ow: i32,
    stride: i32,
}

/// One batch's kernels, run in order (a serial compute encoder orders the dispatches).
struct Pass<'a> {
    metal: &'a Metal,
    commands: Retained<ProtocolObject<dyn MTLCommandBuffer>>,
    encoder: Retained<ProtocolObject<dyn MTLComputeCommandEncoder>>,
    /// Whether the encoder has ended: Metal asserts if one is released before it ends, as it
    /// would be when an error returns from the middle of a pass.
    ended: bool,
}

impl Drop for Pass<'_> {
    fn drop(&mut self) {
        if !self.ended {
            self.encoder.endEncoding();
        }
    }
}

fn int(v: usize) -> Result<i32, Error> {
    i32::try_from(v).map_err(|_| backend(format!("size {v}")))
}

impl Pass<'_> {
    /// Runs `pipeline` on `groups` threadgroups of `threads`, with `buffers` at indices 0.. and
    /// `params` after them.
    fn dispatch<P: Copy>(
        &self,
        pipeline: &Pipeline,
        buffers: &[&View],
        params: &P,
        groups: (usize, usize, usize),
        threads: usize,
    ) {
        let e = &self.encoder;
        e.setComputePipelineState(pipeline);
        for (i, view) in buffers.iter().enumerate() {
            // SAFETY: the kernels' bounds keep within the sizes the callers allocated, and the
            // pass holds the buffers (their views) until it completes.
            unsafe { e.setBuffer_offset_atIndex(Some(&view.buffer), view.offset, i) };
        }
        // SAFETY: `params` is a `#[repr(C)]` struct of `i32`s (or a scalar) laid out as the
        // kernel's parameter at this index, and Metal copies the bytes.
        unsafe {
            e.setBytes_length_atIndex(NonNull::from(params).cast::<c_void>(), size_of::<P>(), buffers.len());
        }
        let size = |width, height, depth| MTLSize { width, height, depth };
        e.dispatchThreadgroups_threadsPerThreadgroup(
            size(groups.0.max(1), groups.1.max(1), groups.2.max(1)),
            size(threads, 1, 1),
        );
    }

    /// An elementwise kernel over `len` threads (and `depth` rows of them).
    fn elementwise<P: Copy>(&self, pipeline: &Pipeline, buffers: &[&View], params: &P, len: usize, depth: usize) {
        self.dispatch(pipeline, buffers, params, (len.div_ceil(THREADS), depth, 1), THREADS);
    }

    /// `out = act(a * b + bias (+ residual))`: `a` is `m x k` with rows `lda` apart, `b` the
    /// dense `k x n` weight, `out` (and `residual`) have rows `ldc` apart.
    #[expect(clippy::too_many_arguments, reason = "a product's operands")]
    fn product(
        &self,
        (a, lda): (&View, usize),
        b: &View,
        bias: Option<&View>,
        (out, ldc): (&View, usize),
        residual: Option<&View>,
        (m, n, k): (usize, usize, usize),
        act: Activation,
    ) -> Result<(), Error> {
        // On an M4 Pro, 128 x 64 tiles ran fastest at every size (2.3 to 2.5 TFLOPS on large
        // products, the multiply-add rate); 128 x 128 ones held too many sums per thread (0.6).
        let (function, bm, bn) =
            if n <= 32 { (&self.metal.kernels.gemm_64x32, 64, 32) } else { (&self.metal.kernels.gemm_128x64, 128, 64) };
        let tiles = n.div_ceil(bn) * m.div_ceil(bm);
        // Too few tiles to fill the GPU, over a long sum: split the sum across threadgroups.
        let splits = if k >= SPLIT_K_MIN && tiles < BUSY_GROUPS && residual.is_none() && bias.is_some() {
            BUSY_GROUPS.div_ceil(tiles).min(k / (SPLIT_K_MIN / 2))
        } else {
            1
        };
        let unused = View::of(&self.metal.unused);
        let params = ProductParams {
            m: int(m)?,
            n: int(n)?,
            k: int(k)?,
            lda: int(lda)?,
            ldc: int(ldc)?,
            ldr: int(ldc)?,
            act: act as i32,
            has_bias: i32::from(bias.is_some()),
            has_residual: i32::from(residual.is_some()),
            ..ProductParams::default()
        };
        let groups = (n.div_ceil(bn), m.div_ceil(bm), splits);
        if splits == 1 {
            let buffers = [a, b, out, bias.unwrap_or(&unused), residual.unwrap_or(&unused)];
            self.dispatch(function, &buffers, &params, groups, THREADS);
            return Ok(());
        }
        let workspace = self.metal.workspace.borrow_mut().view::<f32>(self.metal, splits * m * n)?;
        let buffers = [a, b, &workspace, &unused, &unused];
        self.dispatch(function, &buffers, &params, groups, THREADS);
        let reduce = [int(m)?, int(n)?, int(ldc)?, int(splits)?, act as i32];
        let bias = bias.expect("split only with a bias");
        self.elementwise(&self.metal.kernels.splitk_reduce, &[&workspace, out, bias], &reduce, m * n, 1);
        Ok(())
    }

    /// `layer` on `m` rows of `a` (`lda` apart) into `out` (rows `ldc` apart).
    fn dense(
        &self,
        (a, lda): (&View, usize),
        layer: &DeviceDense,
        (out, ldc): (&View, usize),
        m: usize,
        act: Activation,
    ) -> Result<(), Error> {
        let (weight, bias) = (View::of(&layer.weight), View::of(&layer.bias));
        self.product((a, lda), &weight, Some(&bias), (out, ldc), None, (m, layer.outputs, layer.inputs), act)
    }

    /// `conv` over `images` maps of `height x width` at `input`, into `out`, then `residual` and
    /// ReLU. Returns the output size.
    fn conv(
        &self,
        conv: &DeviceConv,
        input: &View,
        (images, height, width): (usize, usize, usize),
        out: &View,
        residual: Option<&View>,
    ) -> Result<(usize, usize), Error> {
        let (oh, ow) = ((height - 1) / conv.stride + 1, (width - 1) / conv.stride + 1);
        let (m, n) = (images * oh * ow, conv.outputs);
        let unused = View::of(&self.metal.unused);
        let params = ProductParams {
            m: int(m)?,
            n: int(n)?,
            k: int(9 * conv.inputs)?,
            ldc: int(n)?,
            ldr: int(n)?,
            act: Activation::Relu as i32,
            has_bias: 1,
            has_residual: i32::from(residual.is_some()),
            h: int(height)?,
            w: int(width)?,
            c: int(conv.inputs)?,
            oh: int(oh)?,
            ow: int(ow)?,
            stride: int(conv.stride)?,
            ..ProductParams::default()
        };
        let (weight, bias) = (View::of(&conv.weight), View::of(&conv.bias));
        let buffers = [input, &weight, out, &bias, residual.unwrap_or(&unused)];
        self.dispatch(
            &self.metal.kernels.conv_128x64,
            &buffers,
            &params,
            (n.div_ceil(64), m.div_ceil(128), 1),
            THREADS,
        );
        Ok((oh, ow))
    }

    /// Softmax of each row of `logits` (`n x classes`) into columns `offset..` of `out`'s rows
    /// of `width`.
    fn softmax(
        &self,
        logits: &View,
        out: &View,
        n: usize,
        classes: usize,
        (width, offset): (usize, usize),
    ) -> Result<(), Error> {
        let params = [int(n)?, int(classes)?, int(width)?, int(offset)?];
        self.elementwise(&self.metal.kernels.softmax_rows, &[logits, out], &params, n, 1);
        Ok(())
    }

    /// Ends the pass, runs it and waits for it.
    fn finish(mut self) -> Result<(), Error> {
        self.encoder.endEncoding();
        self.ended = true;
        self.commands.commit();
        self.commands.waitUntilCompleted();
        if self.commands.status() == MTLCommandBufferStatus::Completed {
            return Ok(());
        }
        let reason =
            self.commands.error().map_or_else(|| "unknown error".to_owned(), |e| e.localizedDescription().to_string());
        Err(backend(format!("a batch failed: {reason}")))
    }
}

/// A dense layer on the GPU: weight `[inputs][outputs]`, bias `[outputs]`.
#[derive(Debug)]
struct DeviceDense {
    weight: MetalBuffer,
    bias: MetalBuffer,
    inputs: usize,
    outputs: usize,
}

impl DeviceDense {
    fn new(metal: &Metal, layer: &Dense) -> Result<Self, Error> {
        Self::from_parts(metal, &layer.weight, &layer.bias)
    }

    fn from_parts(metal: &Metal, weight: &[f32], bias: &[f32]) -> Result<Self, Error> {
        let (inputs, outputs) = (weight.len() / bias.len(), bias.len());
        Ok(Self { weight: metal.upload(weight)?, bias: metal.upload(bias)?, inputs, outputs })
    }

    /// Two layers on the same input side by side: `[inputs][a outputs | b outputs]`.
    fn side_by_side(metal: &Metal, a: &Dense, b: &Dense) -> Result<Self, Error> {
        let (na, nb) = (a.bias.len(), b.bias.len());
        let weight: Vec<f32> = a
            .weight
            .chunks_exact(na)
            .zip(b.weight.chunks_exact(nb))
            .flat_map(|(x, y)| x.iter().chain(y))
            .copied()
            .collect();
        let bias: Vec<f32> = a.bias.iter().chain(&b.bias).copied().collect();
        Self::from_parts(metal, &weight, &bias)
    }
}

/// A 3x3 convolution with its batch norm folded in, as on the CPU: weight `[9 * inputs][outputs]`.
#[derive(Debug)]
struct DeviceConv {
    weight: MetalBuffer,
    bias: MetalBuffer,
    inputs: usize,
    outputs: usize,
    stride: usize,
}

impl DeviceConv {
    fn new(metal: &Metal, conv: &Conv) -> Result<Self, Error> {
        Ok(Self {
            weight: metal.upload(&conv.weight)?,
            bias: metal.upload(&conv.bias)?,
            inputs: conv.in_channels,
            outputs: conv.out_channels,
            stride: conv.stride,
        })
    }
}

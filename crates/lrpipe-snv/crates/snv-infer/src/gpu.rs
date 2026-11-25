//! Rust GPU backend: the networks of [`crate::cpu`] as CUDA kernels.
//!
//! The kernels (`gpu/kernels.cu`) are compiled when the first model loads, by NVRTC for the
//! device's compute capability, and launched through the CUDA driver API; both libraries are
//! loaded at run time (`cudarc`), so building needs no CUDA. Weights come from the CPU models,
//! already transposed. Products are register-tiled on the CUDA cores (convolutions as implicit
//! im2col, long sums split across blocks), and each LSTM step is one kernel that computes the
//! recurrent product and the cell update of both directions.
//!
//! [`Precision::Float32`](crate::Precision::Float32) follows libtorch's own GPU arithmetic, so
//! that the calls agree with its: cuDNN, which libtorch lets use TF32 tensor cores, multiplies
//! TF32-rounded operands in the LSTMs' input projections and in the convolutions, the latter by
//! the convolution's own weights with the batch norm after. Those operands are rounded the same
//! way here (the products of TF32 values are exact in `float32`, so only the order of the sums
//! differs), and the recurrent products and dense layers stay `float32`, as in libtorch. There is
//! no faster arithmetic on the GPU: [`Precision::Mixed`](crate::Precision::Mixed) is the CPU's.

mod full_alignment;
mod pileup;
#[cfg(test)]
mod tests;

use std::sync::{Arc, OnceLock};

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, DevicePtr, DeviceRepr, LaunchConfig, PinnedHostSlice,
    PushKernelArg, ValidAsZeroBits,
};
use cudarc::nvrtc::{CompileOptions, Ptx};

use crate::Error;
use crate::cpu::{Conv, Dense};

pub use full_alignment::GpuFullAlignmentModel;
pub use pileup::GpuPileupModel;

const SOURCE: &str = include_str!("gpu/kernels.cu");

/// Sums at least this long are split across blocks when a product has too few tiles.
const SPLIT_K_MIN: usize = 1024;

/// Threads per block of every kernel.
const THREADS: u32 = 256;

/// The oldest compute capability the kernels run on: Pascal (6.0). Before Ampere the TF32
/// rounding is done in integer arithmetic, with the same results.
const MIN_CAPABILITY: (i32, i32) = (6, 0);

/// The oldest compute capability with TF32 tensor cores (`mma.sync.m16n8k8.tf32`): Ampere (8.0).
/// There the TF32 products (the convolutions and the LSTMs' input products) run on them, as
/// libtorch's cuDNN runs them; before it, on the CUDA cores with the operands rounded alike.
const TENSOR_CORES: (i32, i32) = (8, 0);

/// What follows a product's bias (and residual): the `ACT_*` codes of the kernels.
#[derive(Debug, Clone, Copy)]
enum Activation {
    None = 0,
    Relu = 1,
    Selu = 2,
}

/// How a float32 product rounds its operands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Rounding {
    /// Full float32.
    Float32 = 0,
    /// To TF32 (10 mantissa bits) first, as cuDNN's tensor cores do in libtorch's float32
    /// convolutions and LSTMs (`allow_tf32`), so the results agree closely with libtorch's.
    Tf32 = 1,
}

/// A raw device address, as the kernels take their pointers.
type DevicePointer = u64;

/// No residual.
const NULL: DevicePointer = 0;

fn backend(err: impl std::fmt::Display) -> Error {
    Error::Backend { message: format!("CUDA: {err}") }
}

/// The kernels compiled for one compute capability, kept for the life of the process.
fn compiled(capability: (i32, i32)) -> Result<String, Error> {
    static PTX: OnceLock<((i32, i32), Result<String, String>)> = OnceLock::new();
    let compile = || {
        let options = CompileOptions {
            options: vec![format!("--gpu-architecture=compute_{}{}", capability.0, capability.1)],
            ..CompileOptions::default()
        };
        cudarc::nvrtc::compile_ptx_with_opts(SOURCE, options).map(|ptx| ptx.to_src()).map_err(|e| format!("{e:?}"))
    };
    let (compiled_for, ptx) = PTX.get_or_init(|| (capability, compile()));
    if *compiled_for == capability {
        ptx.clone().map_err(|e| Error::Backend { message: format!("compiling the GPU kernels: {e}") })
    } else {
        compile().map_err(|e| Error::Backend { message: format!("compiling the GPU kernels: {e}") })
    }
}

/// NVRTC's library, by the name it is loaded under.
const NVRTC: &str = "libnvrtc.so.12";

/// Loads NVRTC from `path`, so that the kernels compile with it: for an NVRTC outside the
/// library path (a Python environment's, say).
///
/// # Errors
///
/// Fails if `path` cannot be loaded as a shared library.
pub fn use_nvrtc(path: &std::path::Path) -> Result<(), Error> {
    let name = std::ffi::CString::new(path.as_os_str().as_encoded_bytes())
        .map_err(|_| Error::Backend { message: format!("{}: a path with a NUL byte", path.display()) })?;
    // SAFETY: dlopen with a NUL-terminated path. The library stays loaded for the life of the
    // process (the handle is never closed), so `cudarc`, which opens it again by its soname,
    // gets this copy.
    let handle = unsafe { libc::dlopen(name.as_ptr(), libc::RTLD_NOW | libc::RTLD_GLOBAL) };
    if handle.is_null() {
        // SAFETY: dlerror returns NULL or a NUL-terminated message valid until the next dl call.
        let reason = unsafe { libc::dlerror() };
        let reason = if reason.is_null() {
            "unknown error".to_owned()
        } else {
            // SAFETY: as above.
            unsafe { std::ffi::CStr::from_ptr(reason) }.to_string_lossy().into_owned()
        };
        return Err(Error::Backend { message: format!("loading NVRTC from {}: {reason}", path.display()) });
    }
    Ok(())
}

/// Loads NVRTC from where CUDA usually installs it, when it is not on the library path:
/// `$CUDA_HOME`, `$CUDA_PATH` and `/usr/local/cuda`. Returns whether it is loaded.
fn find_nvrtc() -> bool {
    let loaded = |name: &std::ffi::CStr| {
        // SAFETY: dlopen by name, with a NUL-terminated name; a handle is kept for the process.
        let handle = unsafe { libc::dlopen(name.as_ptr(), libc::RTLD_NOW | libc::RTLD_GLOBAL) };
        !handle.is_null()
    };
    if loaded(c"libnvrtc.so.12") {
        return true;
    }
    let roots = ["CUDA_HOME", "CUDA_PATH"].into_iter().filter_map(std::env::var_os).map(std::path::PathBuf::from);
    roots
        .chain([std::path::PathBuf::from("/usr/local/cuda")])
        .flat_map(|root| [root.join("lib64").join(NVRTC), root.join("targets/x86_64-linux/lib").join(NVRTC)])
        .any(|path| path.exists() && use_nvrtc(&path).is_ok())
}

/// Sets the primary context of CUDA device `device` to blocking sync, so threads waiting for it
/// sleep. Returns whether it was set.
pub(crate) fn block_while_waiting(device: usize) -> bool {
    use cudarc::driver::sys::{self, CUctx_flags, CUresult};
    // SAFETY: this only looks the shared library up by name.
    if !unsafe { sys::is_culib_present() } {
        return false;
    }
    let Ok(ordinal) = i32::try_from(device) else { return false };
    let mut handle = 0;
    // SAFETY: plain driver calls; `handle` is a valid place for the device handle.
    unsafe {
        sys::cuInit(0) == CUresult::CUDA_SUCCESS
            && sys::cuDeviceGet(&raw mut handle, ordinal) == CUresult::CUDA_SUCCESS
            && sys::cuDevicePrimaryCtxSetFlags_v2(handle, CUctx_flags::CU_CTX_SCHED_BLOCKING_SYNC as u32)
                == CUresult::CUDA_SUCCESS
    }
}

/// The kernels of `kernels.cu`.
#[derive(Debug)]
struct Kernels {
    gemm_128x128: CudaFunction,
    gemm_128x64: CudaFunction,
    gemm_64x32: CudaFunction,
    conv_128x128: CudaFunction,
    conv_128x64: CudaFunction,
    gemm_tf32_128x128: CudaFunction,
    gemm_tf32_128x64: CudaFunction,
    conv_tf32_128x128: CudaFunction,
    conv_tf32_128x64: CudaFunction,
    splitk_reduce: CudaFunction,
    lstm_step: CudaFunction,
    int8_to_float: CudaFunction,
    int32_to_float: CudaFunction,
    pyramid_pool: CudaFunction,
    softmax_rows: CudaFunction,
}

/// Checks that CUDA device `device` can run the networks without compiling anything: the CUDA
/// driver and NVRTC are found, the device exists and its compute capability (or `target`'s) is
/// supported. Returns its context and that capability.
fn usable_device(device: usize, target: Option<(i32, i32)>) -> Result<(Arc<CudaContext>, (i32, i32)), Error> {
    // SAFETY: these only look the shared libraries up by name.
    if !unsafe { cudarc::driver::sys::is_culib_present() } {
        return Err(Error::Backend { message: "no CUDA driver (libcuda.so) found".to_owned() });
    }
    // SAFETY: as above.
    if !find_nvrtc() || !unsafe { cudarc::nvrtc::sys::is_culib_present() } {
        return Err(Error::Backend {
            message: "NVRTC (libnvrtc.so.12), which compiles the GPU kernels, not found: pass --nvrtc \
                      PATH, set CUDA_HOME, or add its directory to LD_LIBRARY_PATH"
                .to_owned(),
        });
    }
    let count = CudaContext::device_count().map_err(backend)?;
    if usize::try_from(count).map_or(true, |count| device >= count) {
        let found =
            if count > 0 { format!("cuda:0 to cuda:{} found", count - 1) } else { "no CUDA device found".to_owned() };
        return Err(Error::Device { device: format!("cuda:{device} ({found})") });
    }
    let ctx = CudaContext::new(device).map_err(backend)?;
    let capability = target.unwrap_or(ctx.compute_capability().map_err(backend)?);
    if capability < MIN_CAPABILITY {
        return Err(Error::Backend {
            message: format!(
                "cuda:{device} has compute capability {}.{}; the Rust GPU backend needs {}.{} (Pascal) or later",
                capability.0, capability.1, MIN_CAPABILITY.0, MIN_CAPABILITY.1
            ),
        });
    }
    Ok((ctx, capability))
}

/// Checks that CUDA device `device` can run the networks, as loading a model there would, without
/// compiling the kernels: a run checks its device before it changes anything.
///
/// # Errors
///
/// Fails without the CUDA driver or NVRTC, for a device index the machine lacks, or for a GPU
/// older than Pascal.
pub fn check_cuda_device(device: usize) -> Result<(), Error> {
    usable_device(device, None).map(|_| ())
}

/// A device's stream and the compiled kernels.
#[derive(Debug)]
pub(crate) struct Gpu {
    stream: Arc<CudaStream>,
    kernels: Kernels,
    /// Streaming multiprocessors of the device.
    sms: usize,
    /// Whether TF32 products run on the tensor cores ([`TENSOR_CORES`]).
    tensor_cores: bool,
    /// Partial sums of split-K products.
    workspace: std::cell::RefCell<Buffer<f32>>,
}

impl Gpu {
    /// Opens CUDA device `device` and compiles the kernels for it.
    fn new(device: usize) -> Result<Self, Error> {
        Self::open(device, None)
    }

    /// Opens CUDA device `device` and compiles the kernels for an older compute capability than
    /// its own, as on an older GPU (the driver runs code compiled for older ones).
    #[cfg(test)]
    pub(crate) fn compiled_for(device: usize, capability: (i32, i32)) -> Result<Self, Error> {
        Self::open(device, Some(capability))
    }

    /// Opens CUDA device `device`, with the kernels compiled for `target` (its own compute
    /// capability when `None`).
    fn open(device: usize, target: Option<(i32, i32)>) -> Result<Self, Error> {
        let (ctx, capability) = usable_device(device, target)?;
        let ptx = compiled(capability)?;
        let module = ctx.load_module(Ptx::from_src(ptx)).map_err(backend)?;
        let f = |name: &str| module.load_function(name).map_err(backend);
        let kernels = Kernels {
            gemm_128x128: f("gemm_128x128")?,
            gemm_128x64: f("gemm_128x64")?,
            gemm_64x32: f("gemm_64x32")?,
            conv_128x128: f("conv_128x128")?,
            conv_128x64: f("conv_128x64")?,
            gemm_tf32_128x128: f("gemm_tf32_128x128")?,
            gemm_tf32_128x64: f("gemm_tf32_128x64")?,
            conv_tf32_128x128: f("conv_tf32_128x128")?,
            conv_tf32_128x64: f("conv_tf32_128x64")?,
            splitk_reduce: f("splitk_reduce")?,
            lstm_step: f("lstm_step")?,
            int8_to_float: f("int8_to_float")?,
            int32_to_float: f("int32_to_float")?,
            pyramid_pool: f("pyramid_pool")?,
            softmax_rows: f("softmax_rows")?,
        };
        let sms = ctx
            .attribute(cudarc::driver::sys::CUdevice_attribute::CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
            .map_err(backend)?;
        Ok(Self {
            stream: ctx.default_stream(),
            kernels,
            sms: usize::try_from(sms).unwrap_or(1),
            tensor_cores: capability >= TENSOR_CORES,
            workspace: Default::default(),
        })
    }

    fn config(grid: (usize, usize, usize)) -> Result<LaunchConfig, Error> {
        let dim = |d: usize| u32::try_from(d.max(1)).map_err(|_| backend(format!("grid of {d} blocks")));
        Ok(LaunchConfig {
            grid_dim: (dim(grid.0)?, dim(grid.1)?, dim(grid.2)?),
            block_dim: (THREADS, 1, 1),
            shared_mem_bytes: 0,
        })
    }

    /// Blocks to cover `len` elements in a grid-stride loop.
    fn blocks(len: usize) -> usize {
        len.div_ceil(THREADS as usize).min(4096)
    }

    /// Copies `values` to the device.
    fn upload<T: DeviceRepr>(&self, values: &[T]) -> Result<CudaSlice<T>, Error> {
        self.stream.clone_htod(values).map_err(backend)
    }

    /// `out = act(a * b + bias (+ residual))`: `a` is `m x k` with rows `lda` apart, `b` the
    /// dense `k x n` weight of `layer`, `out` has rows `ldc` apart.
    fn dense(
        &self,
        a: DevicePointer,
        lda: usize,
        layer: &DeviceDense,
        (out, ldc): (DevicePointer, usize),
        m: usize,
        act: Activation,
    ) -> Result<(), Error> {
        let (n, k) = (layer.outputs, layer.inputs);
        let tf32 = layer.rounding == Rounding::Tf32 && self.tensor_cores;
        let (function, bm, bn) = if tf32 {
            if n <= 64 {
                (&self.kernels.gemm_tf32_128x64, 128, 64)
            } else {
                (&self.kernels.gemm_tf32_128x128, 128, 128)
            }
        } else if n <= 32 || m <= 2048 {
            (&self.kernels.gemm_64x32, 64, 32)
        } else if n <= 64 {
            (&self.kernels.gemm_128x64, 128, 64)
        } else {
            (&self.kernels.gemm_128x128, 128, 128)
        };
        let b = layer.weight.ptr(&self.stream);
        let bias = layer.bias.ptr(&self.stream);
        // Too few tiles to fill the device, over a long sum: split the sum across blocks.
        let tiles = n.div_ceil(bn) * m.div_ceil(bm);
        let splits = if k >= SPLIT_K_MIN && tiles < 4 * self.sms {
            (4 * self.sms).div_ceil(tiles).min(k / (SPLIT_K_MIN / 2))
        } else {
            1
        };
        let (target, target_ld) =
            if splits > 1 { (self.workspace.borrow_mut().ptr(&self.stream, splits * m * n)?, n) } else { (out, ldc) };
        let int = |v: usize| i32::try_from(v).map_err(|_| backend(format!("size {v}")));
        let (lda, target_ld, m_, n_, k) = (int(lda)?, int(target_ld)?, int(m)?, int(n)?, int(k)?);
        let (act_, residual, round) = (act as i32, NULL, layer.rounding as i32);
        let mut launch = self.stream.launch_builder(function);
        launch.arg(&a).arg(&lda).arg(&b).arg(&target).arg(&target_ld).arg(&bias).arg(&residual);
        launch.arg(&m_).arg(&n_).arg(&k).arg(&act_).arg(&round);
        let config = Self::config((n.div_ceil(bn), m.div_ceil(bm), splits))?;
        // SAFETY: the arguments match `DENSE_GEMM`'s parameters, and the caller's buffers hold
        // `m` rows of `lda` and `ldc` values.
        unsafe { launch.launch(config) }.map_err(backend)?;
        if splits > 1 {
            let (ldc, splits_) = (int(ldc)?, int(splits)?);
            let mut launch = self.stream.launch_builder(&self.kernels.splitk_reduce);
            launch.arg(&target).arg(&out).arg(&m_).arg(&n_).arg(&ldc).arg(&splits_).arg(&bias).arg(&act_);
            // SAFETY: the arguments match `splitk_reduce`'s parameters; the workspace holds
            // `splits` sums of `m x n` and `out` `m` rows of `ldc`.
            unsafe { launch.launch(Self::config((Self::blocks(m * n), 1, 1))?) }.map_err(backend)?;
        }
        Ok(())
    }

    /// `conv` over `images` maps of `height x width` at `input`, into `out`, then `residual` and
    /// ReLU. Returns the output size.
    fn conv(
        &self,
        conv: &DeviceConv,
        input: DevicePointer,
        (images, height, width): (usize, usize, usize),
        out: DevicePointer,
        residual: DevicePointer,
    ) -> Result<(usize, usize), Error> {
        let (oh, ow) = ((height - 1) / conv.stride + 1, (width - 1) / conv.stride + 1);
        let (m, n) = (images * oh * ow, conv.outputs);
        let (function, bn) = match (conv.rounding == Rounding::Tf32 && self.tensor_cores, n <= 64) {
            (true, true) => (&self.kernels.conv_tf32_128x64, 64),
            (true, false) => (&self.kernels.conv_tf32_128x128, 128),
            (false, true) => (&self.kernels.conv_128x64, 64),
            (false, false) => (&self.kernels.conv_128x128, 128),
        };
        let int = |v: usize| i32::try_from(v).map_err(|_| backend(format!("size {v}")));
        let (h, w, c, oh_, ow_, stride) =
            (int(height)?, int(width)?, int(conv.inputs)?, int(oh)?, int(ow)?, int(conv.stride)?);
        let (m_, n_, act) = (int(m)?, int(n)?, Activation::Relu as i32);
        let b = conv.weight.ptr(&self.stream);
        let bias = conv.bias.ptr(&self.stream);
        let scale = conv.scale.ptr(&self.stream);
        let mut launch = self.stream.launch_builder(function);
        launch.arg(&input).arg(&h).arg(&w).arg(&c).arg(&oh_).arg(&ow_).arg(&stride).arg(&b).arg(&out);
        let round = conv.rounding as i32;
        launch.arg(&scale).arg(&bias).arg(&residual).arg(&m_).arg(&n_).arg(&act).arg(&round);
        let config = Self::config((n.div_ceil(bn), m.div_ceil(128), 1))?;
        // SAFETY: the arguments match `CONV_GEMM`'s parameters; `input` holds `images` maps of
        // `height x width x inputs` and `out` (and `residual`) `m x n` values.
        unsafe { launch.launch(config) }.map(|_| ()).map_err(backend)?;
        Ok((oh, ow))
    }

    /// Softmax of each row of `logits` (`n x classes`) into columns `offset..` of `out`'s rows
    /// of `width`.
    fn softmax(
        &self,
        logits: DevicePointer,
        out: DevicePointer,
        n: usize,
        classes: usize,
        (width, offset): (usize, usize),
    ) -> Result<(), Error> {
        let int = |v: usize| i32::try_from(v).map_err(|_| backend(format!("size {v}")));
        let (n_, classes, width, offset) = (int(n)?, int(classes)?, int(width)?, int(offset)?);
        let mut launch = self.stream.launch_builder(&self.kernels.softmax_rows);
        launch.arg(&logits).arg(&out).arg(&n_).arg(&classes).arg(&width).arg(&offset);
        // SAFETY: the arguments match `softmax_rows`'s parameters and the buffers' sizes.
        unsafe { launch.launch(Self::config((n.div_ceil(THREADS as usize), 1, 1))?) }.map(|_| ()).map_err(backend)
    }

    /// Converts `len` integers at `input` to floats at `out`, divided by `divisor` for `i8`.
    fn to_float(
        &self,
        input: DevicePointer,
        out: DevicePointer,
        len: usize,
        divisor: Option<f32>,
    ) -> Result<(), Error> {
        let len_ = i64::try_from(len).map_err(backend)?;
        let config = Self::config((Self::blocks(len), 1, 1))?;
        let mut launch;
        match divisor {
            Some(divisor) => {
                launch = self.stream.launch_builder(&self.kernels.int8_to_float);
                launch.arg(&input).arg(&out).arg(&len_).arg(&divisor);
                // SAFETY: the arguments match `int8_to_float`'s parameters and the buffers' sizes.
                unsafe { launch.launch(config) }.map(|_| ()).map_err(backend)
            }
            None => {
                launch = self.stream.launch_builder(&self.kernels.int32_to_float);
                launch.arg(&input).arg(&out).arg(&len_);
                // SAFETY: the arguments match `int32_to_float`'s parameters and the buffers' sizes.
                unsafe { launch.launch(config) }.map(|_| ()).map_err(backend)
            }
        }
    }
}

/// A device buffer that grows to the largest size asked for.
#[derive(Debug, Default)]
struct Buffer<T> {
    slice: Option<CudaSlice<T>>,
}

impl<T: DeviceRepr + ValidAsZeroBits> Buffer<T> {
    /// The buffer with room for at least `len` values.
    fn ensure(&mut self, stream: &Arc<CudaStream>, len: usize) -> Result<&mut CudaSlice<T>, Error> {
        if self.slice.as_ref().is_none_or(|s| s.len() < len) {
            self.slice = None;
            self.slice = Some(stream.alloc_zeros(len.max(1)).map_err(backend)?);
        }
        Ok(self.slice.as_mut().expect("allocated above"))
    }

    /// The buffer's address, after making room for `len` values.
    fn ptr(&mut self, stream: &Arc<CudaStream>, len: usize) -> Result<DevicePointer, Error> {
        Ok(self.ensure(stream, len)?.ptr(stream))
    }
}

/// Which way a [`Staging`] buffer copies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    ToDevice,
    ToHost,
}

/// Page-locked host memory that copies go through, grown as needed: the GPU reads and writes it
/// directly, at full bus speed.
#[derive(Debug)]
struct Staging<T> {
    pinned: Option<PinnedHostSlice<T>>,
    direction: Direction,
}

impl<T: DeviceRepr + ValidAsZeroBits + Copy> Staging<T> {
    fn new(direction: Direction) -> Self {
        Self { pinned: None, direction }
    }

    fn ensure(&mut self, stream: &Arc<CudaStream>, len: usize) -> Result<&mut PinnedHostSlice<T>, Error> {
        if self.pinned.as_ref().is_none_or(|p| p.len() < len) {
            self.pinned = None;
            // Write-combined memory is fast to write and to copy to the device, but slow to read.
            let flags = match self.direction {
                Direction::ToDevice => cudarc::driver::sys::CU_MEMHOSTALLOC_WRITECOMBINED,
                Direction::ToHost => 0,
            };
            // SAFETY: the memory is only read after it has been written (`upload` copies the
            // values in first; `download` reads what the device wrote, after synchronizing).
            let pinned = unsafe { stream.context().alloc_pinned_with_flags(len.max(1), flags) };
            self.pinned = Some(pinned.map_err(backend)?);
        }
        Ok(self.pinned.as_mut().expect("allocated above"))
    }

    /// Queues a copy of `values` to `dst` on the stream; the stream must be synchronized before
    /// the next call.
    fn upload(&mut self, stream: &Arc<CudaStream>, values: &[T], dst: DevicePointer) -> Result<(), Error> {
        let pinned = self.ensure(stream, values.len())?.as_mut_slice().map_err(backend)?;
        let pinned = &mut pinned[..values.len()];
        pinned.copy_from_slice(values);
        // SAFETY: `pinned` is page-locked and stays untouched until the stream is synchronized;
        // the caller's `dst` holds `values.len()` values.
        unsafe { cudarc::driver::result::memcpy_htod_async(dst, pinned, stream.cu_stream()) }.map_err(backend)
    }

    /// Copies `len` values from `src` once the stream has run up to here.
    fn download(&mut self, stream: &Arc<CudaStream>, src: DevicePointer, len: usize) -> Result<Vec<T>, Error> {
        let pinned = &mut self.ensure(stream, len)?.as_mut_slice().map_err(backend)?[..len];
        // SAFETY: `pinned` is page-locked and read only after the stream is synchronized below;
        // the caller's `src` holds `len` values.
        unsafe { cudarc::driver::result::memcpy_dtoh_async(pinned, src, stream.cu_stream()) }.map_err(backend)?;
        stream.synchronize().map_err(backend)?;
        Ok(pinned.to_vec())
    }
}

/// The device address of a slice.
trait Address {
    fn ptr(&self, stream: &CudaStream) -> DevicePointer;
}

impl<T> Address for CudaSlice<T> {
    fn ptr(&self, stream: &CudaStream) -> DevicePointer {
        self.device_ptr(stream).0
    }
}

/// A dense layer on the device: weight `[inputs][outputs]`, bias `[outputs]`.
#[derive(Debug)]
struct DeviceDense {
    weight: CudaSlice<f32>,
    bias: CudaSlice<f32>,
    inputs: usize,
    outputs: usize,
    /// How the float32 product rounds.
    rounding: Rounding,
}

impl DeviceDense {
    fn new(gpu: &Gpu, layer: &Dense) -> Result<Self, Error> {
        Self::from_parts(gpu, &layer.weight, &layer.bias)
    }

    fn from_parts(gpu: &Gpu, weight: &[f32], bias: &[f32]) -> Result<Self, Error> {
        let (inputs, outputs) = (weight.len() / bias.len(), bias.len());
        Ok(Self { weight: gpu.upload(weight)?, bias: gpu.upload(bias)?, inputs, outputs, rounding: Rounding::Float32 })
    }

    /// Two layers on the same input side by side: `[inputs][a outputs | b outputs]`.
    fn side_by_side(gpu: &Gpu, a: &Dense, b: &Dense) -> Result<Self, Error> {
        let (na, nb) = (a.bias.len(), b.bias.len());
        let weight: Vec<f32> = a
            .weight
            .chunks_exact(na)
            .zip(b.weight.chunks_exact(nb))
            .flat_map(|(x, y)| x.iter().chain(y))
            .copied()
            .collect();
        let bias: Vec<f32> = a.bias.iter().chain(&b.bias).copied().collect();
        Self::from_parts(gpu, &weight, &bias)
    }
}

/// A 3x3 convolution and its batch norm on the device: weight `[9 * inputs][outputs]`.
#[derive(Debug)]
struct DeviceConv {
    weight: CudaSlice<f32>,
    bias: CudaSlice<f32>,
    inputs: usize,
    outputs: usize,
    stride: usize,
    /// How the float32 product rounds.
    rounding: Rounding,
    /// The batch norm's scale, applied after the product: `weight` holds the convolution's own
    /// weights, as libtorch multiplies them.
    scale: CudaSlice<f32>,
}

impl DeviceConv {
    fn new(gpu: &Gpu, conv: &Conv) -> Result<Self, Error> {
        // libtorch rounds the convolution's own weights to TF32 and applies the batch norm after,
        // so the batch norm's scale comes back out of the folded weights. A scale of zero leaves
        // nothing to divide, and those channels keep the folded (zero) weights.
        let unfolded: Vec<f32> = conv
            .weight
            .chunks_exact(conv.out_channels)
            .flat_map(|row| row.iter().zip(&conv.scale).map(|(&w, &s)| if s == 0.0 { w } else { w / s }))
            .collect();
        Ok(Self {
            scale: gpu.upload(&conv.scale)?,
            weight: gpu.upload(&unfolded)?,
            bias: gpu.upload(&conv.bias)?,
            inputs: conv.in_channels,
            outputs: conv.out_channels,
            stride: conv.stride,
            rounding: Rounding::Tf32,
        })
    }
}

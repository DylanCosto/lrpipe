//! Network inference behind a backend-neutral trait.
//!
//! The callers hand over feature tensors exactly as the original builds them (row-major integers)
//! and get back the softmax probabilities the Python code would see, computed in Rust: on the CPU
//! ([`CpuPileupModel`], [`CpuFullAlignmentModel`]), and with the `gpu` feature on CUDA GPUs
//! (`GpuPileupModel`, `GpuFullAlignmentModel`: kernels compiled at run time by NVRTC and launched
//! through the CUDA driver, both loaded at run time), and with the `metal` feature on Apple GPUs
//! (`MetalPileupModel`, `MetalFullAlignmentModel`: kernels compiled at run time by the Metal
//! framework). No C or C++ library is linked. The probabilities agree with the original's libtorch
//! to about `1e-5` on the CPU and on Metal (`float32` throughout) and, on CUDA, round the same
//! operands to TF32 as libtorch's cuDNN does.

mod checkpoint;
mod cpu;
#[cfg(feature = "gpu")]
mod gpu;
mod int8;
#[cfg(all(feature = "metal", target_os = "macos"))]
mod metal;

use std::path::{Path, PathBuf};

pub use checkpoint::{CheckpointTensor, TensorValues, read_state_dict};
pub use cpu::{CpuFullAlignmentModel, CpuPileupModel};
#[cfg(feature = "gpu")]
pub use gpu::{GpuFullAlignmentModel, GpuPileupModel, check_cuda_device, use_nvrtc};
#[cfg(all(feature = "metal", target_os = "macos"))]
pub use metal::{MetalFullAlignmentModel, MetalPileupModel};

/// Whether this CPU has the 8-bit instructions `Precision::Mixed` runs the networks with
/// (AVX-VNNI or AVX-512 VNNI with AVX2 on x86-64, dotprod on aarch64). Without them the networks
/// run in float32 whatever the precision asked for.
#[must_use]
pub fn eight_bit_supported() -> bool {
    int8::supported()
}

/// Makes the threads that wait for GPU `device` (copying a batch's results back, or
/// synchronizing) sleep until it is done, instead of spinning on a core as CUDA does by default.
/// Only how the CPU waits changes, not what the GPU computes. Returns whether it was set: without
/// the `gpu` feature or a CUDA driver it does nothing.
#[must_use]
pub fn block_while_waiting(device: usize) -> bool {
    #[cfg(feature = "gpu")]
    return gpu::block_while_waiting(device);
    #[cfg(not(feature = "gpu"))]
    {
        let _ = device;
        false
    }
}

/// How [`load_pileup_model`] and [`load_full_alignment_model`] run a network.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ModelConfig {
    /// Where it runs.
    pub device: InferenceDevice,
    /// Its arithmetic.
    pub precision: Precision,
}

impl ModelConfig {
    /// Fails unless this build can run on the configured device.
    ///
    /// # Errors
    ///
    /// For CUDA without the `gpu` feature, or Metal without the `metal` feature on macOS.
    pub fn check(&self) -> Result<(), Error> {
        match self.device {
            InferenceDevice::Cuda(_) if !cfg!(feature = "gpu") => {
                Err(Error::Backend { message: "this build has no GPU backend (the `gpu` feature)".to_owned() })
            }
            InferenceDevice::Metal if !cfg!(all(feature = "metal", target_os = "macos")) => Err(Error::Backend {
                message: "this build has no Metal backend (the `metal` feature, on macOS)".to_owned(),
            }),
            _ => Ok(()),
        }
    }
}

/// Loads a pileup checkpoint (`pileup.pt`) for `config`.
///
/// # Errors
///
/// Fails if this build cannot use the device, the GPU cannot be used (see
/// `GpuPileupModel::load` and `MetalPileupModel::load`), or the checkpoint cannot be loaded (see [`CpuPileupModel::load`]).
pub fn load_pileup_model(path: &Path, config: ModelConfig) -> Result<Box<dyn PileupModel + Send>, Error> {
    config.check()?;
    match config.device {
        #[cfg(feature = "gpu")]
        InferenceDevice::Cuda(i) => Ok(Box::new(GpuPileupModel::load(path, i)?)),
        #[cfg(all(feature = "metal", target_os = "macos"))]
        InferenceDevice::Metal => Ok(Box::new(MetalPileupModel::load(path)?)),
        _ => Ok(Box::new(CpuPileupModel::load(path)?.with_precision(config.precision))),
    }
}

/// Loads a full-alignment checkpoint (`full_alignment.pt`) for `config`, for matrices of `depth`
/// rows ([`FULL_ALIGNMENT_DEPTH`] for ONT, 55 for HiFi and Illumina).
///
/// # Errors
///
/// Fails if this build cannot use the device, the GPU cannot be used (see
/// `GpuFullAlignmentModel::load`), or the checkpoint cannot be loaded (see
/// [`CpuFullAlignmentModel::load`]).
pub fn load_full_alignment_model(
    path: &Path,
    config: ModelConfig,
    depth: usize,
) -> Result<Box<dyn FullAlignmentModel + Send>, Error> {
    config.check()?;
    match config.device {
        #[cfg(feature = "gpu")]
        InferenceDevice::Cuda(i) => Ok(Box::new(GpuFullAlignmentModel::load(path, i)?.with_matrix_depth(depth))),
        #[cfg(all(feature = "metal", target_os = "macos"))]
        InferenceDevice::Metal => Ok(Box::new(MetalFullAlignmentModel::load(path)?.with_matrix_depth(depth))),
        _ => Ok(Box::new(CpuFullAlignmentModel::load(path)?.with_matrix_depth(depth).with_precision(config.precision))),
    }
}

/// The arithmetic of a network's forward pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum Precision {
    /// `float32` throughout, as the original runs its networks (the default).
    #[default]
    Float32,
    /// Faster arithmetic on the CPU that changes a few calls: with AVX-VNNI (x86-64) or `sdot`
    /// (aarch64), the full-alignment convolutions after the first and the pileup network's
    /// products (all but its first input product and its heads') in 8 bits, with candidates whose
    /// two most likely classes are close scored again in `float32`. A GPU runs `float32` either way: there, faster arithmetic saved
    /// little.
    Mixed,
}

/// The gap between the two largest probabilities of a head.
#[must_use]
pub fn top_two_gap(head: &[f32]) -> f32 {
    let (mut first, mut second) = (f32::MIN, f32::MIN);
    for &p in head {
        if p > first {
            (first, second) = (p, first);
        } else if p > second {
            second = p;
        }
    }
    first - second
}

/// Where a model runs.
///
/// CUDA results follow libtorch's GPU arithmetic (TF32 operands where cuDNN has them), so they
/// differ from CPU ones in the last bits; compare them with the original's own `--use_gpu` run.
/// Metal results are `float32`, as the CPU's, in another order of sums.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum InferenceDevice {
    /// The CPU (the default).
    #[default]
    Cpu,
    /// A CUDA GPU, by index (`cuda:N`).
    Cuda(usize),
    /// The system's Apple GPU, through Metal (`metal`).
    Metal,
}

impl std::str::FromStr for InferenceDevice {
    type Err = Error;

    /// Parses `cpu` or `cuda:N` (`cuda` alone means `cuda:0`), as the original's `--device`, or
    /// `metal`.
    fn from_str(s: &str) -> Result<Self, Error> {
        match s {
            "cpu" => Ok(Self::Cpu),
            "metal" => Ok(Self::Metal),
            "cuda" => Ok(Self::Cuda(0)),
            _ => s
                .strip_prefix("cuda:")
                .and_then(|i| i.parse().ok())
                .map(Self::Cuda)
                .ok_or_else(|| Error::Device { device: s.to_owned() }),
        }
    }
}

/// Genome positions in a pileup tensor: 16 flanking bases on each side of the candidate.
pub const PILEUP_POSITIONS: usize = 33;

/// Feature channels per position in a pileup tensor (`shared/param_p.py`, `channel`).
pub const PILEUP_CHANNELS: usize = 18;

/// Number of `i32` values in one pileup tensor.
pub const PILEUP_TENSOR_LEN: usize = PILEUP_POSITIONS * PILEUP_CHANNELS;

/// Probabilities per candidate from the pileup model: 21 GT21 classes then 3 genotype classes.
pub const PILEUP_OUTPUTS: usize = 24;

/// A network that scores pileup tensors.
pub trait PileupModel {
    /// Scores a batch of pileup tensors.
    ///
    /// `inputs` holds `n` tensors of [`PILEUP_TENSOR_LEN`] values back to back; the result holds
    /// `n` rows of [`PILEUP_OUTPUTS`] probabilities. Results can depend on how candidates are
    /// split into batches, so callers that need parity with the original must batch as it
    /// does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InputLength`] if `inputs` is not a whole number of tensors, and
    /// [`Error::Backend`] if the backend fails.
    fn predict(&mut self, inputs: &[i32]) -> Result<Vec<f32>, Error>;
}

impl<M: PileupModel + ?Sized> PileupModel for Box<M> {
    fn predict(&mut self, inputs: &[i32]) -> Result<Vec<f32>, Error> {
        (**self).predict(inputs)
    }
}

/// Rows (read slots) in an ONT full-alignment tensor (`shared/param_f.py`, `matrix_depth_dict`).
pub const FULL_ALIGNMENT_DEPTH: usize = 89;

/// Genome positions in a full-alignment tensor: 16 flanking bases on each side of the candidate.
pub const FULL_ALIGNMENT_POSITIONS: usize = 33;

/// Feature channels per read and position in a full-alignment tensor without dwell time
/// (`shared/param_f.py`, `channel`).
pub const FULL_ALIGNMENT_CHANNELS: usize = 8;

/// Feature channels with dwell time (`--enable_dwell_time`, the `_with_mv` models).
pub const FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME: usize = 9;

/// Number of `i8` values in one ONT full-alignment tensor.
pub const FULL_ALIGNMENT_TENSOR_LEN: usize = full_alignment_tensor_len(FULL_ALIGNMENT_DEPTH, FULL_ALIGNMENT_CHANNELS);

/// Number of `i8` values in a full-alignment tensor of `depth` rows and `channels` channels.
pub const fn full_alignment_tensor_len(depth: usize, channels: usize) -> usize {
    depth * FULL_ALIGNMENT_POSITIONS * channels
}

/// Probabilities per candidate from the full-alignment model: 21 GT21 classes, 3 genotype
/// classes, then 33 classes for each of the two indel lengths.
pub const FULL_ALIGNMENT_OUTPUTS: usize = 90;

/// A network that scores full-alignment tensors.
pub trait FullAlignmentModel {
    /// Scores a batch of full-alignment tensors.
    ///
    /// `inputs` holds `n` tensors of [`full_alignment_tensor_len`] values back to back
    /// ([`FULL_ALIGNMENT_TENSOR_LEN`] for ONT), each
    /// `[depth][position][channel]` row-major and `i8` as `CreateTensorFullAlignment` makes them;
    /// the result holds `n` rows of [`FULL_ALIGNMENT_OUTPUTS`] probabilities. Results can depend
    /// on how candidates are split into batches, so callers that need parity with the original
    /// must batch as it does (`param_f.predictBatchSize * threads` candidates per batch).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InputLength`] if `inputs` is not a whole number of tensors, and
    /// [`Error::Backend`] if the backend fails.
    fn predict(&mut self, inputs: &[i8]) -> Result<Vec<f32>, Error>;

    /// Input channels per position: [`FULL_ALIGNMENT_CHANNELS`], or
    /// [`FULL_ALIGNMENT_CHANNELS_WITH_DWELL_TIME`] for a `_with_mv` model.
    fn channels(&self) -> usize;
}

impl<M: FullAlignmentModel + ?Sized> FullAlignmentModel for Box<M> {
    fn predict(&mut self, inputs: &[i8]) -> Result<Vec<f32>, Error> {
        (**self).predict(inputs)
    }

    fn channels(&self) -> usize {
        (**self).channels()
    }
}

/// Errors from loading or running a model.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// The input is not a whole number of tensors.
    #[error("input of {len} values is not a multiple of the tensor size {tensor_len}")]
    InputLength {
        /// Number of values passed in.
        len: usize,
        /// Number of values in one tensor.
        tensor_len: usize,
    },
    /// A checkpoint file could not be read.
    #[error("cannot read checkpoint {path:?}: {what}")]
    Checkpoint {
        /// The checkpoint file.
        path: PathBuf,
        /// What went wrong.
        what: String,
    },
    /// The checkpoint lacks a tensor the network needs.
    #[error("checkpoint has no tensor named {name:?}")]
    MissingWeight {
        /// Name of the missing state-dict entry.
        name: String,
    },
    /// A device name could not be parsed, or the device is not available in this build.
    #[error("unusable inference device {device:?}: use cpu, cuda, cuda:N or metal")]
    Device {
        /// The device as given.
        device: String,
    },
    /// The inference backend reported an error.
    #[error("inference backend: {message}")]
    Backend {
        /// The backend's error message.
        message: String,
    },
}

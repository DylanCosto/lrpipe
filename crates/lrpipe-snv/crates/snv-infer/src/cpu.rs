//! Pure-Rust CPU backend.
//!
//! The networks of the original's `model.py`, computed with the `gemm` crate's SIMD matrix
//! products: the LSTMs as one product for every input step plus one per recurrent step, the
//! convolutions as im2col products with their batch norms folded in. The arithmetic is `float32`
//! as in PyTorch, but sums run in another order and batch norm is folded into the convolution
//! weights, so probabilities agree with libtorch's to a few units in the sixth decimal rather than
//! bit for bit.

mod activation;
pub(crate) mod conv;
pub(crate) mod full_alignment;
pub(crate) mod pileup;
#[cfg(test)]
pub(crate) mod synthetic;
#[cfg(test)]
mod tests;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use gemm::Parallelism;

use crate::Error;
use crate::checkpoint::{CheckpointTensor, TensorValues, read_state_dict};

#[cfg(test)]
pub(crate) use conv::BatchNorm;
pub(crate) use conv::{Conv, Epilogue};
pub use full_alignment::CpuFullAlignmentModel;
pub use pileup::CpuPileupModel;

/// Whether a matrix product replaces the output or is added to it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Product {
    /// `out = a * b`.
    Overwrite,
    /// `out += a * b`.
    Accumulate,
}

/// A row-major matrix of `rows` rows, `cols` columns and rows `stride` values apart.
#[derive(Debug, Clone, Copy)]
struct Strided {
    rows: usize,
    cols: usize,
    stride: usize,
}

impl Strided {
    const fn dense(rows: usize, cols: usize) -> Self {
        Self { rows, cols, stride: cols }
    }

    /// Values a slice must hold for the matrix to start at its beginning.
    fn span(self) -> usize {
        if self.rows == 0 { 0 } else { (self.rows - 1) * self.stride + self.cols }
    }
}

/// `out (m x n) = a (m x k) * b (k x n)`, or `+=` with [`Product::Accumulate`], on one thread;
/// `b` is dense and row-major.
///
/// # Panics
///
/// If the shapes disagree or a slice is too short for its shape.
fn matmul(out: &mut [f32], out_shape: Strided, a: &[f32], a_shape: Strided, b: &[f32], product: Product) {
    let (m, k, n) = (a_shape.rows, a_shape.cols, out_shape.cols);
    assert_eq!(out_shape.rows, m, "matmul rows");
    assert!(a.len() >= a_shape.span() && out.len() >= out_shape.span() && b.len() >= k * n, "matmul operands");
    // Output rows must not overlap; input rows may (the unit-stride convolution's views do).
    assert!(out_shape.stride >= n, "matmul output rows overlap");
    if m == 0 || n == 0 {
        return;
    }
    let isize = |x: usize| isize::try_from(x).expect("matrix strides fit in isize");
    // SAFETY: the asserts above keep every element `gemm` touches inside the three slices and
    // give each output element its own place; `out` cannot overlap `a` or `b` as it is borrowed
    // mutably, and `gemm` only reads `a` and `b`.
    unsafe {
        gemm::gemm(
            m,
            n,
            k,
            out.as_mut_ptr(),
            1,
            isize(out_shape.stride),
            product == Product::Accumulate,
            a.as_ptr(),
            1,
            isize(a_shape.stride),
            b.as_ptr(),
            1,
            isize(n),
            1.0,
            1.0,
            false,
            false,
            false,
            Parallelism::None,
        );
    }
}

/// `torch.selu`.
fn selu(x: f32) -> f32 {
    const ALPHA: f32 = 1.673_263_2;
    const SCALE: f32 = 1.050_701;
    if x > 0.0 { SCALE * x } else { SCALE * ALPHA * x.exp_m1() }
}

/// `torch.softmax` over one row, in place.
fn softmax(row: &mut [f32]) {
    let largest = row.iter().fold(f32::NEG_INFINITY, |m, &x| m.max(x));
    let mut sum = 0.0;
    for x in row.iter_mut() {
        *x = (*x - largest).exp();
        sum += *x;
    }
    for x in row.iter_mut() {
        *x /= sum;
    }
}

/// An `nn.Linear` layer, its weight transposed to `[in][out]`.
#[derive(Debug)]
pub(crate) struct Dense {
    pub(crate) weight: Vec<f32>,
    pub(crate) bias: Vec<f32>,
}

impl Dense {
    fn from_state(state: &mut Weights, prefix: &str) -> Result<Self, Error> {
        let (shape, weight) = state.take(&format!("{prefix}.weight"))?;
        let [outputs, inputs] = state.dims(&format!("{prefix}.weight"), &shape)?;
        let bias = state.take_vector(&format!("{prefix}.bias"), outputs)?;
        Ok(Self { weight: transpose(&weight, outputs, inputs), bias })
    }

    pub(crate) fn inputs(&self) -> usize {
        self.weight.len() / self.bias.len()
    }

    /// `F.linear` over `n` rows of `x`, followed by `activation` on each value.
    fn forward(&self, x: &[f32], n: usize, activation: fn(f32) -> f32) -> Vec<f32> {
        let outputs = self.bias.len();
        let mut y = vec![0.0; n * outputs];
        matmul(
            &mut y,
            Strided::dense(n, outputs),
            x,
            Strided::dense(n, self.inputs()),
            &self.weight,
            Product::Overwrite,
        );
        for row in y.chunks_exact_mut(outputs) {
            for (v, b) in row.iter_mut().zip(&self.bias) {
                *v = activation(*v + b);
            }
        }
        y
    }
}

/// A dense head: `Linear`, SELU, `Linear`, SELU, softmax.
#[derive(Debug)]
pub(crate) struct Head {
    pub(crate) hidden: Dense,
    pub(crate) logits: Dense,
}

impl Head {
    fn from_state(state: &mut Weights, hidden: &str, logits: &str) -> Result<Self, Error> {
        Ok(Self { hidden: Dense::from_state(state, hidden)?, logits: Dense::from_state(state, logits)? })
    }

    pub(crate) fn classes(&self) -> usize {
        self.logits.bias.len()
    }

    /// Probabilities of `n` rows of `x`, written to columns `at..at + classes` of `out`'s rows of
    /// `width` values.
    fn forward_into(&self, x: &[f32], n: usize, out: &mut [f32], width: usize, at: usize) {
        let hidden = self.hidden.forward(x, n, selu);
        let logits = self.logits.forward(&hidden, n, selu);
        let classes = self.classes();
        for (row, logits) in out.chunks_exact_mut(width).zip(logits.chunks_exact(classes)) {
            let row = &mut row[at..at + classes];
            row.copy_from_slice(logits);
            softmax(row);
        }
    }
}

/// `m` (`rows x cols`, row-major) transposed.
fn transpose(m: &[f32], rows: usize, cols: usize) -> Vec<f32> {
    let mut t = vec![0.0; m.len()];
    for (r, row) in m.chunks_exact(cols).enumerate() {
        for (c, &v) in row.iter().enumerate() {
            t[c * rows + r] = v;
        }
    }
    t
}

/// The `float32` tensors of a checkpoint, removed as layers claim them.
#[derive(Debug)]
pub(crate) struct Weights {
    path: PathBuf,
    tensors: HashMap<String, CheckpointTensor>,
}

impl Weights {
    pub(crate) fn load(path: &Path) -> Result<Self, Error> {
        Ok(Self::new(path, read_state_dict(path)?))
    }

    pub(crate) fn new(path: &Path, tensors: Vec<(String, CheckpointTensor)>) -> Self {
        Self { path: path.to_owned(), tensors: tensors.into_iter().collect() }
    }

    fn bad(&self, what: String) -> Error {
        Error::Checkpoint { path: self.path.clone(), what }
    }

    /// The shape and values of the `float32` tensor `name`.
    fn take(&mut self, name: &str) -> Result<(Vec<usize>, Vec<f32>), Error> {
        let tensor = self.tensors.remove(name).ok_or_else(|| Error::MissingWeight { name: name.to_owned() })?;
        match tensor.values {
            TensorValues::F32(values) => Ok((tensor.shape, values)),
            TensorValues::I64(_) => Err(self.bad(format!("{name} is not float32"))),
        }
    }

    /// The `float32` vector `name`, which must hold `len` values.
    fn take_vector(&mut self, name: &str, len: usize) -> Result<Vec<f32>, Error> {
        let (shape, values) = self.take(name)?;
        if shape != [len] {
            return Err(self.bad(format!("{name} has shape {shape:?}, expected [{len}]")));
        }
        Ok(values)
    }

    /// `shape` as `N` dimensions.
    fn dims<const N: usize>(&self, name: &str, shape: &[usize]) -> Result<[usize; N], Error> {
        shape.try_into().map_err(|_| self.bad(format!("{name} has shape {shape:?}, expected {N} dimensions")))
    }
}

#[cfg(test)]
mod unit_tests {
    use super::*;

    #[test]
    fn matmul_reads_and_writes_strided_rows() {
        // a: rows 0 and 2 of a 3 x 4 buffer, first two columns; b: 2 x 3.
        let a = [1.0, 2.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 3.0, 4.0, 9.0, 9.0];
        let b = [1.0, 0.0, 2.0, 0.0, 1.0, 3.0];
        let mut out = vec![1.0; 2 * 5];
        let (out_shape, a_shape) = (Strided { rows: 2, cols: 3, stride: 5 }, Strided { rows: 2, cols: 2, stride: 8 });
        matmul(&mut out, out_shape, &a, a_shape, &b, Product::Accumulate);
        assert_eq!(out, [2.0, 3.0, 9.0, 1.0, 1.0, 4.0, 5.0, 19.0, 1.0, 1.0]);
        matmul(&mut out, out_shape, &a, a_shape, &b, Product::Overwrite);
        assert_eq!(out, [1.0, 2.0, 8.0, 1.0, 1.0, 3.0, 4.0, 18.0, 1.0, 1.0]);
    }

    #[test]
    fn softmax_and_selu_match_torch() {
        let mut row = [1.0, 2.0, 3.0];
        softmax(&mut row);
        // torch.softmax(torch.tensor([1., 2., 3.]), 0)
        for (got, want) in row.iter().zip([0.090_030_57, 0.244_728_48, 0.665_240_94]) {
            assert!((got - want).abs() < 1e-7, "{row:?}");
        }
        // torch.selu(torch.tensor([-1., 0., 2.]))
        for (x, want) in [(-1.0, -1.111_330_7), (0.0, 0.0), (2.0, 2.101_402)] {
            assert!((selu(x) - want).abs() < 1e-6, "selu({x}) = {}", selu(x));
        }
    }

    #[test]
    fn transpose_swaps_rows_and_columns() {
        assert_eq!(transpose(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2, 3), [1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
    }
}

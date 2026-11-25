//! The bi-LSTM pileup network (the `*_P` class in the original's `model.py`).

use std::path::Path;

use super::activation::lstm_cell_step;
use super::{Dense, Head, Product, Strided, Weights, matmul, selu, transpose};
use crate::int8::{self, QuantMatrix, prepare_rows, quantize_signed};
use crate::{
    Error, PILEUP_CHANNELS, PILEUP_OUTPUTS, PILEUP_POSITIONS, PILEUP_TENSOR_LEN, PileupModel, Precision, top_two_gap,
};

/// With 8-bit products, candidates whose two most likely GT21 or genotype classes are closer than
/// this are scored again in `float32`. On ONT's model the 8-bit network moved no probability by
/// more than 0.053 on 13,285 chr20 candidates, so a call could only change where the two were
/// within 0.106; lrpipe-snv checks a run's own candidates against half this gap.
const RESCORE_GAP: f32 = 0.16;

/// Candidates scored together: enough rows for efficient matrix products, few enough that a
/// layer's gates stay in cache.
const SUB_BATCH: usize = 24;

/// A bidirectional, single-layer `nn.LSTM` with `batch_first=True`.
#[derive(Debug)]
pub(crate) struct BiLstm {
    pub(crate) inputs: usize,
    pub(crate) hidden: usize,
    /// Both directions' `weight_ih` transposed side by side: `[inputs][2 * 4 * hidden]`.
    pub(crate) input_weight: Vec<f32>,
    /// `bias_ih + bias_hh` of both directions: `[2 * 4 * hidden]`.
    pub(crate) bias: Vec<f32>,
    /// Each direction's `weight_hh` transposed: `[hidden][4 * hidden]`.
    pub(crate) hidden_weight: [Vec<f32>; 2],
}

impl BiLstm {
    fn from_state(state: &mut Weights, prefix: &str, hidden: usize) -> Result<Self, Error> {
        let gates = 4 * hidden;
        let mut input_weights = Vec::with_capacity(2);
        let mut hidden_weight = Vec::with_capacity(2);
        let mut bias = Vec::with_capacity(2 * gates);
        for suffix in ["", "_reverse"] {
            let name = format!("{prefix}.weight_ih_l0{suffix}");
            let (shape, w_ih) = state.take(&name)?;
            let [rows, inputs] = state.dims(&name, &shape)?;
            if rows != gates {
                return Err(state.bad(format!("{name} has shape {shape:?}, expected {gates} rows")));
            }
            let name = format!("{prefix}.weight_hh_l0{suffix}");
            let (shape, w_hh) = state.take(&name)?;
            if shape != [gates, hidden] {
                return Err(state.bad(format!("{name} has shape {shape:?}, expected [{gates}, {hidden}]")));
            }
            let b_ih = state.take_vector(&format!("{prefix}.bias_ih_l0{suffix}"), gates)?;
            let b_hh = state.take_vector(&format!("{prefix}.bias_hh_l0{suffix}"), gates)?;
            bias.extend(b_ih.iter().zip(&b_hh).map(|(a, b)| a + b));
            input_weights.push((inputs, transpose(&w_ih, gates, inputs)));
            hidden_weight.push(transpose(&w_hh, gates, hidden));
        }
        let inputs = input_weights[0].0;
        if input_weights[1].0 != inputs {
            return Err(state.bad(format!("{prefix}: the two directions take different input sizes")));
        }
        let mut input_weight = Vec::with_capacity(inputs * 2 * gates);
        for row in 0..inputs {
            for (_, w) in &input_weights {
                input_weight.extend_from_slice(&w[row * gates..(row + 1) * gates]);
            }
        }
        let hidden_weight = hidden_weight.try_into().expect("two directions");
        Ok(Self { inputs, hidden, input_weight, bias, hidden_weight })
    }

    /// `nn.LSTM.forward` with `hx=None` on `n` sequences of `steps` steps (`x` is
    /// `[n][steps][inputs]`); writes `[n][steps][2 * hidden]` to `out`, forward then reverse.
    /// With `quantized`, its products run in 8 bits.
    fn forward(
        &self,
        x: &[f32],
        (n, steps): (usize, usize),
        out: &mut Vec<f32>,
        scratch: &mut LstmScratch,
        quantized: Option<&QuantLstm>,
    ) {
        let recurrent = match quantized {
            Some(QuantLstm::Fused(products)) => return self.forward_fused(x, (n, steps), out, scratch, products),
            Some(QuantLstm::Recurrent(products)) => Some(products),
            None => None,
        };
        let (hidden, gates) = (self.hidden, 4 * self.hidden);
        let rows = n * steps;
        // Input contributions and both biases, for every step of both directions at once.
        let LstmScratch { pre, cell, bytes } = scratch;
        pre.clear();
        for _ in 0..rows {
            pre.extend_from_slice(&self.bias);
        }
        let x_shape = Strided::dense(rows, self.inputs);
        matmul(pre, Strided::dense(rows, 2 * gates), x, x_shape, &self.input_weight, Product::Accumulate);

        out.clear();
        out.resize(rows * 2 * hidden, 0.0);
        cell.resize(n * hidden, 0.0);
        for (direction, w_hh) in self.hidden_weight.iter().enumerate() {
            cell.fill(0.0);
            let step_at = |s: usize| if direction == 0 { s } else { steps - 1 - s };
            for s in 0..steps {
                let t = step_at(s);
                let gates_at = t * 2 * gates + direction * gates;
                if s > 0 {
                    // The previous step's hidden state, a strided view of `out`.
                    let h_at = step_at(s - 1) * 2 * hidden + direction * hidden;
                    match recurrent {
                        Some(products) => products[direction].accumulate(
                            &out[h_at..],
                            (n, steps * 2 * hidden),
                            &mut pre[gates_at..],
                            steps * 2 * gates,
                            bytes,
                        ),
                        None => {
                            let h_shape = Strided { rows: n, cols: hidden, stride: steps * 2 * hidden };
                            let g_shape = Strided { rows: n, cols: gates, stride: steps * 2 * gates };
                            matmul(&mut pre[gates_at..], g_shape, &out[h_at..], h_shape, w_hh, Product::Accumulate);
                        }
                    }
                }
                for i in 0..n {
                    let g = &pre[i * steps * 2 * gates + gates_at..][..gates];
                    let c = &mut cell[i * hidden..][..hidden];
                    let h = &mut out[(i * steps + t) * 2 * hidden + direction * hidden..][..hidden];
                    lstm_cell_step(g, c, h);
                }
            }
        }
    }
}

/// Buffers an LSTM layer reuses across calls.
#[derive(Debug, Default)]
struct LstmScratch {
    /// Gate pre-activations, `[n][steps][2 * 4 * hidden]`.
    pre: Vec<f32>,
    /// One direction's cell state, `[n][hidden]`.
    cell: Vec<f32>,
    /// Quantized inputs of an 8-bit product.
    bytes: Vec<u8>,
}

/// An LSTM layer's products in 8 bits, per direction.
#[derive(Debug)]
enum QuantLstm {
    /// The recurrent products only, for a layer whose inputs are not in (-1, 1) (the first
    /// layer's counts): its input product stays `float32`, for every step at once.
    Recurrent([QuantMatrix; 2]),
    /// For a layer whose inputs are in (-1, 1) (the first layer's outputs): each step's input
    /// and recurrent products as one, of `[x_t, h_prev]` by `[W_ih; W_hh]`, just before that
    /// step's cell update, so a step's gates stay in cache.
    Fused([QuantMatrix; 2]),
}

impl QuantLstm {
    fn recurrent(lstm: &BiLstm) -> Self {
        let gates = 4 * lstm.hidden;
        Self::Recurrent(lstm.hidden_weight.each_ref().map(|w| QuantMatrix::new(w, lstm.hidden, gates)))
    }

    fn fused(lstm: &BiLstm) -> Self {
        let gates = 4 * lstm.hidden;
        Self::Fused(std::array::from_fn(|direction| {
            let mut stacked = Vec::with_capacity((lstm.inputs + lstm.hidden) * gates);
            for row in lstm.input_weight.chunks_exact(2 * gates) {
                stacked.extend_from_slice(&row[direction * gates..][..gates]);
            }
            stacked.extend_from_slice(&lstm.hidden_weight[direction]);
            QuantMatrix::new(&stacked, lstm.inputs + lstm.hidden, gates)
        }))
    }
}

impl BiLstm {
    /// [`BiLstm::forward`] with [`QuantLstm::Fused`] products.
    fn forward_fused(
        &self,
        x: &[f32],
        (n, steps): (usize, usize),
        out: &mut Vec<f32>,
        scratch: &mut LstmScratch,
        products: &[QuantMatrix; 2],
    ) {
        let (inputs, hidden, gates) = (self.inputs, self.hidden, 4 * self.hidden);
        let k = inputs + hidden;
        let LstmScratch { pre, cell, bytes } = scratch;
        out.clear();
        out.resize(n * steps * 2 * hidden, 0.0);
        cell.resize(n * hidden, 0.0);
        pre.resize(n * gates, 0.0);
        for (direction, product) in products.iter().enumerate() {
            debug_assert_eq!(product.inputs(), k);
            cell.fill(0.0);
            let bias = &self.bias[direction * gates..][..gates];
            let step_at = |s: usize| if direction == 0 { s } else { steps - 1 - s };
            for s in 0..steps {
                let t = step_at(s);
                // Each row `[x_t, h_prev]`; the first step's `h_prev` is zero, which the
                // buffer's fill already encodes.
                prepare_rows(bytes, n, k);
                for i in 0..n {
                    let row = &mut bytes[i * k..][..k];
                    quantize_signed(&x[(i * steps + t) * inputs..][..inputs], &mut row[..inputs]);
                    if s > 0 {
                        let h_at = (i * steps + step_at(s - 1)) * 2 * hidden + direction * hidden;
                        quantize_signed(&out[h_at..][..hidden], &mut row[inputs..]);
                    }
                }
                for g in pre.chunks_exact_mut(gates) {
                    g.copy_from_slice(bias);
                }
                product.multiply_rows(bytes, n, pre, gates);
                for i in 0..n {
                    let h = &mut out[(i * steps + t) * 2 * hidden + direction * hidden..][..hidden];
                    lstm_cell_step(&pre[i * gates..][..gates], &mut cell[i * hidden..][..hidden], h);
                }
            }
        }
    }
}

/// The pileup network's products in 8 bits, for [`Precision::Mixed`]: every product but the first
/// layer's input product (whose inputs are counts) and the heads'.
#[derive(Debug)]
struct QuantPileup {
    lstm1: QuantLstm,
    lstm2: QuantLstm,
    l4: QuantMatrix,
}

/// Activations the network reuses across calls.
#[derive(Debug, Default)]
struct Buffers {
    input: Vec<f32>,
    lstm1: Vec<f32>,
    lstm2: Vec<f32>,
    scratch: LstmScratch,
}

/// The original's bi-LSTM pileup network (its `*_P` class with `predict=True`, no indel-length
/// heads), in Rust on the CPU.
#[derive(Debug)]
pub struct CpuPileupModel {
    pub(crate) lstm1: BiLstm,
    pub(crate) lstm2: BiLstm,
    pub(crate) l4: Dense,
    pub(crate) gt21: Head,
    pub(crate) genotype: Head,
    /// The 8-bit products, for [`Precision::Mixed`] on a CPU with AVX-VNNI or Arm's `sdot`.
    quantized: Option<QuantPileup>,
    buffers: Buffers,
}

impl CpuPileupModel {
    /// Loads a pileup checkpoint (`pileup.pt`).
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read as a PyTorch state dict, lacks a required tensor, or has
    /// tensors of shapes the network cannot use.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, Error> {
        Self::from_state(&mut Weights::load(path.as_ref())?)
    }

    pub(crate) fn from_state(state: &mut Weights) -> Result<Self, Error> {
        let model = Self {
            lstm1: BiLstm::from_state(state, "LSTM1", 128)?,
            lstm2: BiLstm::from_state(state, "LSTM2", 160)?,
            l4: Dense::from_state(state, "L4")?,
            gt21: Head::from_state(state, "L5_1", "Y_gt21_logits")?,
            genotype: Head::from_state(state, "L5_2", "Y_genotype_logits")?,
            quantized: None,
            buffers: Buffers::default(),
        };
        let flattened = PILEUP_POSITIONS * 2 * model.lstm2.hidden;
        if model.lstm1.inputs != PILEUP_CHANNELS
            || model.lstm2.inputs != 2 * model.lstm1.hidden
            || model.l4.inputs() != flattened
            || model.gt21.classes() + model.genotype.classes() != PILEUP_OUTPUTS
        {
            return Err(state.bad("layer sizes do not chain into the pileup network".to_owned()));
        }
        Ok(model)
    }

    /// Runs the network at `precision` instead of [`Precision::Float32`]: for
    /// [`Precision::Mixed`] on a CPU with AVX-VNNI or Arm's `sdot`, its products in 8 bits (all
    /// but the first layer's input product and the heads'), with candidates whose call is close
    /// scored again in `float32`.
    #[must_use]
    pub fn with_precision(mut self, precision: Precision) -> Self {
        self.quantized = (precision == Precision::Mixed && int8::supported()).then(|| QuantPileup {
            lstm1: QuantLstm::recurrent(&self.lstm1),
            lstm2: QuantLstm::fused(&self.lstm2),
            l4: QuantMatrix::new(&self.l4.weight, self.l4.inputs(), self.l4.bias.len()),
        });
        self
    }

    /// The network on `n` candidates of `inputs`, with the 8-bit products if `quantized`; `n` rows
    /// of [`PILEUP_OUTPUTS`] probabilities.
    fn forward(&mut self, inputs: &[i32], n: usize, quantized: Option<&QuantPileup>) -> Vec<f32> {
        let Buffers { input, lstm1, lstm2, scratch } = &mut self.buffers;
        input.clear();
        input.extend(inputs.iter().map(|&v| v as f32));
        self.lstm1.forward(input, (n, PILEUP_POSITIONS), lstm1, scratch, quantized.map(|q| &q.lstm1));
        self.lstm2.forward(lstm1, (n, PILEUP_POSITIONS), lstm2, scratch, quantized.map(|q| &q.lstm2));
        let x = match quantized {
            Some(q) => {
                let outputs = self.l4.bias.len();
                let mut x = vec![0.0; n * outputs];
                q.l4.accumulate(lstm2, (n, self.l4.inputs()), &mut x, outputs, &mut scratch.bytes);
                for row in x.chunks_exact_mut(outputs) {
                    for (v, b) in row.iter_mut().zip(&self.l4.bias) {
                        *v = selu(*v + b);
                    }
                }
                x
            }
            None => self.l4.forward(lstm2, n, selu),
        };
        let mut out = vec![0.0; n * PILEUP_OUTPUTS];
        self.gt21.forward_into(&x, n, &mut out, PILEUP_OUTPUTS, 0);
        self.genotype.forward_into(&x, n, &mut out, PILEUP_OUTPUTS, self.gt21.classes());
        out
    }
}

impl PileupModel for CpuPileupModel {
    fn predict(&mut self, inputs: &[i32]) -> Result<Vec<f32>, Error> {
        if !inputs.len().is_multiple_of(PILEUP_TENSOR_LEN) {
            return Err(Error::InputLength { len: inputs.len(), tensor_len: PILEUP_TENSOR_LEN });
        }
        let quantized = self.quantized.take();
        let mut out = Vec::with_capacity(inputs.len() / PILEUP_TENSOR_LEN * PILEUP_OUTPUTS);
        for batch in inputs.chunks(SUB_BATCH * PILEUP_TENSOR_LEN) {
            out.extend(self.forward(batch, batch.len() / PILEUP_TENSOR_LEN, quantized.as_ref()));
        }
        if quantized.is_some() {
            let close: Vec<usize> = (0..out.len() / PILEUP_OUTPUTS)
                .filter(|&i| {
                    let p = &out[i * PILEUP_OUTPUTS..][..PILEUP_OUTPUTS];
                    top_two_gap(&p[..21]) < RESCORE_GAP || top_two_gap(&p[21..]) < RESCORE_GAP
                })
                .collect();
            let subset: Vec<i32> = close
                .iter()
                .flat_map(|&i| &inputs[i * PILEUP_TENSOR_LEN..(i + 1) * PILEUP_TENSOR_LEN])
                .copied()
                .collect();
            let mut exact = Vec::with_capacity(close.len() * PILEUP_OUTPUTS);
            for batch in subset.chunks(SUB_BATCH * PILEUP_TENSOR_LEN) {
                exact.extend(self.forward(batch, batch.len() / PILEUP_TENSOR_LEN, None));
            }
            for (&i, row) in close.iter().zip(exact.as_chunks::<PILEUP_OUTPUTS>().0) {
                out[i * PILEUP_OUTPUTS..(i + 1) * PILEUP_OUTPUTS].copy_from_slice(row);
            }
        }
        self.quantized = quantized;
        Ok(out)
    }
}

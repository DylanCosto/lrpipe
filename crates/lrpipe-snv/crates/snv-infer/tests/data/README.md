# Test data

`libtorch_*.f32` are the probabilities libtorch 2.13.0 (conda-forge's `cpu_mkl` build, one CPU
thread) computed for the synthetic weights and inputs of `src/cpu/synthetic.rs` and
`src/cpu/tests.rs`, as little-endian `f32`, written on 2026-11-01 before lrpipe-snv dropped
libtorch:

- `libtorch_pileup.f32`: the pileup network, 120 candidates of 24 outputs.
- `libtorch_full_alignment.f32`: the full-alignment network, 3 candidates of 90 outputs (89 rows).
- `libtorch_full_alignment_55.f32`: the same at 55 rows (HiFi), 2 candidates.

The Rust networks are held to 1e-4 of them.

//! Long-read copy number variant calling, ported from Spectre 0.2.1.
//!
//! Each module is one stage of Spectre's `CNVCaller` and reproduces it exactly, quirks
//! included; `docs/SPECTRE_SPEC.md` describes the original, and the section numbers in the
//! docs here refer to it.
//!
//! So far: the inputs ([`genome`], [`metadata`], [`coverage`]), normalization
//! ([`normalize`]), the DEL/DUP thresholds ([`thresholds`]), raw candidates ([`call`],
//! [`candidate`]), their refinement into calls ([`refine`]), quality scores ([`score`]) and
//! the output files ([`output`]); [`pipeline`] runs them in Spectre's order.
//!
//! [`extensions`] holds changes to the method, each off by default.

pub mod breakpoints;
pub mod call;
pub mod candidate;
pub mod coverage;
pub mod dip;
pub mod error;
pub mod extensions;
pub mod gc;
pub mod genome;
pub mod hmm;
pub mod loh;
pub mod merged;
pub mod metadata;
pub mod normalize;
pub mod number;
pub mod output;
pub mod pipeline;
pub mod refine;
pub mod score;
pub mod sex;
pub mod snv;
pub mod thresholds;

pub use error::Error;

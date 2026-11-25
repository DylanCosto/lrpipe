//! lrpipe-snv: small variant calling from long reads (a Rust port of Clair3 2.0.3), as a
//! library so that `lrpipe snv` runs it; [`cli::main`] takes its command line.

mod chunk;
pub mod cli;
mod contigs;

pub use contigs::is_default_contig;
mod genotyping;
mod gvcf;
mod parallel;
mod pipeline;
mod python;
mod rescore;
mod settings;

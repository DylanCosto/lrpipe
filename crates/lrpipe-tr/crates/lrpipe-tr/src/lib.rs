//! Tandem repeat genotyping from long reads.
//!
//! For each locus of a catalog ([`catalog`]), the reads that span it with both flanks
//! ([`reads`]) give one repeat sequence each ([`extract`]): the bases between the flanks, cut
//! at the reference boundaries through the alignment and, where an indel lies near a boundary,
//! by aligning the flank itself ([`align`]). The sequences are split into alleles (by haplotag,
//! else by length), each allele's sequence voted from its reads ([`consensus`]) and its motifs
//! counted ([`genotype`]), and written as VCF ([`vcf`]). [`run`] drives it over a catalog in parallel;
//! [`stream`] lets another program's pass over the reads feed it instead.

pub mod align;
pub mod catalog;
pub mod cli;
pub mod consensus;
pub mod error;
pub mod extract;
pub mod genotype;
pub mod reads;
pub mod run;
pub mod stream;
pub mod vcf;

pub use error::Error;

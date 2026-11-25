//! A Rust port of WhatsHap (<https://github.com/whatshap/whatshap>), aiming
//! for output that is byte-identical to WhatsHap 2.8.
pub mod align;
pub mod alignment;
pub mod bam;
pub mod cli;
pub mod compare;
pub mod compat;
pub mod core;
pub mod fasta;
pub mod fastx;
pub mod genotype;
pub mod graph;
pub mod hapcut2vcf;
pub mod haplotag;
pub mod haplotagphase;
pub mod hts;
pub mod memory;
pub mod merge;
pub mod pedigree;
pub mod phase;
pub mod readselect;
pub mod split;
pub mod stats;
pub mod unphase;
pub mod variants;
pub mod vcf;
pub mod xopen;

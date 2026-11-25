//! lrpipe-sv: structural variant calling from long reads, ported from Sniffles2 2.8.1.

mod alleles;
mod call;
mod cluster;
mod combine;
mod config;
mod consensus;
mod fasta;
mod genotype;
mod methyl;
mod phasing;
mod pickle;
mod pymath;
mod read_stats;
mod reads;
mod snf;
mod vcf;

/// The command line.
pub mod cli;

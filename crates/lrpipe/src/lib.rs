//! lrpipe: one command from a sample's long reads, aligned or not (bonefish aligns them), to its
//! phased SNVs and indels (lrpipe-snv), SVs, depth, tandem repeats and methylation (`lrpipe sv`'s
//! single read pass) and CNVs (lrpipe-cnv), merged into one VCF ([`merge`]), with a manifest of
//! every input, tool, command, time and peak memory.
//!
//! Each stage runs its tool as a program (the lrpipe binary's subcommands, or lrpipe-snv). Two runs
//! compare record by record (see [`compare`]), to check that a change to a tool keeps its
//! results.

pub mod checks;
pub mod cli;
pub mod compare;
pub mod config;
pub mod cores;
pub mod input;
pub mod manifest;
pub mod merge;
pub mod merge_reads;
pub mod plan;
pub mod process;
pub mod provenance;
pub mod qc;
pub mod run;
pub mod setup;
pub mod sex;
pub mod snv_model;
pub mod tui;

//! Binned read depth from BAM/CRAM, with output identical to mosdepth 0.3.14.
//!
//! Two ways in:
//!
//! - [`run::run`] reads an indexed BAM/CRAM and writes mosdepth's files, counting several
//!   contigs at once.
//! - [`ContigCounter`] takes records one at a time, so a tool that already reads the
//!   alignments (an SV caller, say) can get depth from the same pass: fill an [`Alignment`]
//!   per record, [`push`](ContigCounter::push) it, then turn the finished [`Coverage`] into a
//!   [`ContigReport`] and hand it to a [`Writer`], or read the per-base depth directly.

pub mod alignment;
/// The `lrpipe-depth` command line.
pub mod cli;
pub mod counter;
pub mod reader;
pub mod regions;
pub mod report;
pub mod run;
pub mod stats;

pub use alignment::{Alignment, Cigar};
pub use counter::{ContigCounter, Coverage, Filter, Mode};
pub use regions::{BedRegions, By, Region};
pub use report::{ContigReport, ReportOptions, Scratch, Writer};

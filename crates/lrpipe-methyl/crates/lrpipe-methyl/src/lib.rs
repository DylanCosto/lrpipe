//! lrpipe-methyl: base modification pileups (5mC, 5hmC, 6mA and any other `MM` code) from aligned
//! long reads, ONT and PacBio alike.

pub mod cli;
pub mod cpg;
pub mod input;
pub mod model;
pub mod mods;
pub mod output;
pub mod pileup;
pub mod read;
pub mod run;
pub mod stream;

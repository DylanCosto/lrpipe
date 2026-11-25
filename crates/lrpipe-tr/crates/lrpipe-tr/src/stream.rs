//! Genotyping from reads another program streams past, so one pass over the alignments can
//! feed several callers.
//!
//! A program reading a contig in pieces gives each piece a [`Collector`] over the contig's
//! loci, offers it every record ([`Collector::spans_any`] first, to skip decoding reads no
//! locus needs, from its alignment span and soft clips), then merges the pieces' observations in
//! read order ([`Observations::merge`])
//! and genotypes the loci, in slices it may spread over threads ([`Observations::records`]).
//! The records are the ones [`crate::run::run`] writes from the same reads.

use std::ops::Range;

use crate::catalog::Locus;
use crate::extract::{RefWindow, may_reach, soft_clips};
use crate::genotype::{Observation, one_per_read};
use crate::reads::Read;
use crate::run::{RunOptions, evenly, genotype_locus, observe};

/// Observations of one piece of reads, by locus index, in read order.
#[derive(Debug)]
pub struct Collector<'c> {
    loci: &'c [Locus],
    reference: &'c RefWindow,
    opts: &'c RunOptions,
    obs: Vec<(u32, Observation)>,
    /// The longest locus, which bounds how far before a read a locus it reaches can start.
    longest: i64,
}

impl<'c> Collector<'c> {
    /// A collector for `loci` (one contig's, sorted by start) with `reference` covering them.
    pub fn new(loci: &'c [Locus], reference: &'c RefWindow, opts: &'c RunOptions) -> Self {
        let longest = loci.iter().map(|l| l.end - l.start).max().unwrap_or(0);
        Collector { loci, reference, opts, obs: Vec::new(), longest }
    }

    /// The loci a read over `[start, end)` could reach: those starting from `longest` before it
    /// to its end ([`may_reach`] decides).
    fn candidates(&self, start: i64, end: i64) -> Range<usize> {
        let lo = self.loci.partition_point(|l| l.start < start - self.longest);
        let hi = self.loci.partition_point(|l| l.start <= end);
        lo..hi.max(lo)
    }

    /// Whether a read aligned over `[start, end)`, with `clips` soft-clipped bases at its start
    /// and end, could give an allele of any locus: spanning it with both flanks, or one flank
    /// and the other in a clip.
    pub fn spans_any(&self, start: i64, end: i64, clips: [usize; 2]) -> bool {
        let flank = self.opts.extract.flank;
        self.loci[self.candidates(start, end)].iter().any(|l| may_reach(start, end, clips, l.start, l.end, flank))
    }

    /// Records the read's observation of every locus it reaches.
    pub fn push(&mut self, read: &Read) {
        let (flank, clips) = (self.opts.extract.flank, soft_clips(&read.cigar));
        for i in self.candidates(read.start, read.end) {
            let locus = &self.loci[i];
            if !may_reach(read.start, read.end, clips, locus.start, locus.end, flank) {
                continue;
            }
            if let Some(o) = observe(read, locus, self.reference, self.opts) {
                self.obs.push((i as u32, o));
            }
        }
    }

    /// The piece's observations, for [`Observations::merge`].
    pub fn into_observations(self) -> Vec<(u32, Observation)> {
        self.obs
    }
}

/// A contig's observations per locus, in read order.
#[derive(Debug, Default)]
pub struct Observations {
    by_locus: Vec<Vec<Observation>>,
}

impl Observations {
    /// Merges the pieces' observations, given in read order, for `n_loci` loci; each locus
    /// keeps one observation per read and `max_reads` of them spread evenly, as [`crate::run::run`] does.
    pub fn merge(pieces: impl IntoIterator<Item = Vec<(u32, Observation)>>, n_loci: usize, max_reads: usize) -> Self {
        let mut by_locus: Vec<Vec<Observation>> = vec![Vec::new(); n_loci];
        for piece in pieces {
            for (i, o) in piece {
                by_locus[i as usize].push(o);
            }
        }
        let by_locus = by_locus.into_iter().map(|v| evenly(one_per_read(v), max_reads)).collect();
        Observations { by_locus }
    }

    /// The VCF records of loci `range` (indexes into the contig's loci), and how many were
    /// called.
    pub fn records(
        &self,
        loci: &[Locus],
        range: Range<usize>,
        reference: &RefWindow,
        opts: &RunOptions,
    ) -> (String, usize) {
        let mut text = String::new();
        let mut called = 0;
        for i in range {
            let (record, was_called) = genotype_locus(&loci[i], &self.by_locus[i], reference, opts);
            called += usize::from(was_called);
            text.push_str(&record);
        }
        (text, called)
    }
}

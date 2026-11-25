//! The per-position counts `CreateTensorPileup` passes to gVCF output (`--gvcf`).
//!
//! With `--gvcf` the C pileup also records, for every position of the region, the reads
//! supporting the reference and the reads in all ([`SiteDepth`]), and `CreateTensorPileup`
//! walks the chunk's positions with the reference base to build the non-variant blocks.

use crate::Error;
use crate::counts::{PileupCounts, SiteDepth};
use crate::hts::Fasta;

/// Reference bases fetched on each side of the chunk (`param.expandReferenceRegion`).
const REFERENCE_PADDING: i64 = 1000;

/// One position of a chunk, as `CreateTensorPileup` passes it to `make_gvcf_online`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Site {
    /// 1-based position.
    pub pos: i64,
    /// Uppercase reference base.
    pub reference: u8,
    /// Reads with the reference base.
    pub n_ref: u64,
    /// Reads in all (see [`SiteDepth::total`]).
    pub n_total: u64,
}

/// A pileup chunk's positions for gVCF output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChunkSites {
    contig: String,
    start: i64,
    end: i64,
    /// 0-based position read for 1-based position `pos` is `pos + shift`.
    shift: i64,
    /// The C code's buffers: 0-based positions `bounds.0..bounds.1`.
    bounds: (i64, i64),
    /// First buffer index of the emptiness test (`start`, inclusive) and its end.
    empty_test: (i64, i64),
    positions: Vec<i64>,
    depths: Vec<SiteDepth>,
    /// `reference_sequence_from` of `reference_start..=end + 1000`, uppercased.
    reference: Vec<u8>,
    reference_start: i64,
}

impl ChunkSites {
    /// Takes the gVCF counts out of `counts`, the pileup of the chunk `start..end` (1-based,
    /// end exclusive, as `CreateTensorPileup` iterates it) extended to `extend_start`, and
    /// fetches the reference around the chunk.
    ///
    /// # Errors
    ///
    /// Fails if the reference cannot be read.
    pub(crate) fn new(
        counts: &mut PileupCounts,
        fasta: &Fasta,
        contig: &str,
        start: i64,
        end: i64,
        extend_start: i64,
    ) -> Result<Self, Error> {
        // `offset = 0 if ctg_start == 1 else 1`: right when the region starts at 2 or later,
        // one position off for chunks starting at 2..=34, whose region starts at 1.
        let offset = if start == 1 { 0 } else { 1 };
        let reference_start = (start - REFERENCE_PADDING).max(1);
        let reference_end = end + REFERENCE_PADDING;
        let mut reference = fasta.fetch(contig, (reference_start - 1) as i32, (reference_end - 1) as i32)?;
        reference.make_ascii_uppercase();
        let bounds = counts.bounds;
        Ok(Self {
            contig: contig.to_owned(),
            start,
            end,
            shift: bounds.0 - extend_start + offset,
            bounds,
            empty_test: (start - extend_start + offset, end + 1 - extend_start + offset),
            positions: std::mem::take(&mut counts.positions),
            depths: std::mem::take(&mut counts.depths),
            reference,
            reference_start,
        })
    }

    /// The contig.
    pub fn contig(&self) -> &str {
        &self.contig
    }

    /// The chunk's first position (`ctg_start`, at least 1).
    pub fn start(&self) -> i64 {
        self.start
    }

    /// The chunk's end (`ctg_end`); [`sites`](Self::sites) stops before it.
    pub fn end(&self) -> i64 {
        self.end
    }

    /// Whether no read covers the chunk (`empty_pileup_flag`), so the original writes one `N`
    /// block (`write_empty_pileup`) instead of the sites. The test spans one more position than the
    /// sites do.
    pub fn is_empty(&self) -> bool {
        let len = self.bounds.1 - self.bounds.0;
        let (from, to) = (self.empty_test.0.clamp(0, len), self.empty_test.1.clamp(0, len));
        let (from, to) = (self.bounds.0 + from, self.bounds.0 + to);
        let first = self.positions.partition_point(|&p| p < from);
        self.positions[first..]
            .iter()
            .zip(&self.depths[first..])
            .take_while(|&(&p, _)| p < to)
            .all(|(_, d)| d.total == 0)
    }

    /// The sites of `start..end`, skipping those past the end of the reference.
    pub fn sites(&self) -> impl Iterator<Item = Site> + '_ {
        let mut cursor = 0;
        (self.start..self.end).filter_map(move |pos| {
            let reference = *self.reference.get(usize::try_from(pos - self.reference_start).ok()?)?;
            let at = pos + self.shift;
            while cursor < self.positions.len() && self.positions[cursor] < at {
                cursor += 1;
            }
            let depth = match self.positions.get(cursor) {
                Some(&p) if p == at => self.depths[cursor],
                _ => SiteDepth::default(),
            };
            Some(Site { pos, reference, n_ref: depth.reference, n_total: depth.total })
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sites(start: i64, end: i64, extend_start: i64, bounds: (i64, i64), covered: &[(i64, u64)]) -> ChunkSites {
        let offset = if start == 1 { 0 } else { 1 };
        ChunkSites {
            contig: "chr1".to_owned(),
            start,
            end,
            shift: bounds.0 - extend_start + offset,
            bounds,
            empty_test: (start - extend_start + offset, end + 1 - extend_start + offset),
            positions: covered.iter().map(|&(p, _)| p).collect(),
            depths: covered.iter().map(|&(_, n)| SiteDepth { reference: n, total: n }).collect(),
            reference: vec![b'A'; 10_000],
            reference_start: (start - REFERENCE_PADDING).max(1),
        }
    }

    #[test]
    fn sites_read_the_previous_zero_based_position() {
        // Chunk 101..200: region "chr1:67-233", which htslib starts at 0-based 66.
        let s = sites(101, 200, 68, (66, 233), &[(99, 5), (100, 7)]);
        let got: Vec<_> = s.sites().take(2).map(|site| (site.pos, site.n_total)).collect();
        assert_eq!(got, [(101, 7), (102, 0)]);
        assert!(!s.is_empty());
    }

    #[test]
    fn early_chunks_read_one_position_late() {
        // Chunk 10..20 extends to 1, region "chr1:0-53" starts at 0-based 0, but the offset is 1.
        let s = sites(10, 20, 1, (0, 53), &[(9, 3), (10, 4)]);
        assert_eq!(s.sites().next().map(|site| site.n_total), Some(4));
    }

    #[test]
    fn emptiness_includes_the_end_position() {
        // Coverage only at 1-based 200 (0-based 199), the chunk's end.
        let s = sites(101, 200, 68, (66, 233), &[(199, 2)]);
        assert!(!s.is_empty());
        assert!(s.sites().all(|site| site.n_total == 0));
        assert!(sites(101, 200, 68, (66, 233), &[(200, 2)]).is_empty());
    }
}

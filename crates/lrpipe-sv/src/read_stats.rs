//! `--read-stats`: the reads' count, yield, length N50, mean length and mean base quality, and
//! the secondary, supplementary and unmapped records, counted from the records the pass reads
//! anyway (for `lrpipe run`'s QC summary).
//!
//! The reads are the primary alignments of mapped reads, as cramino counts them: an unplaced
//! unmapped read (past every contig, where the pass does not read) is not seen, and `lrpipe run`
//! takes their number from the index. CRAM's records come without their qualities, which the
//! pass does not decode (decoding them made it 65% slower on HG002 ONT chr20, 3.8 s against
//! 2.3 s): their mean quality is unknown.

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, Result};

const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// BAM's quality of a base whose quality is missing: a read with it first has none.
const MISSING_QUALITY: u8 = 0xff;

/// Counts of the records read.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadStats {
    /// Primary alignments of mapped reads: the reads.
    pub reads: u64,
    /// Their bases (sequence lengths).
    pub bases: u64,
    /// The bases of the reads that have qualities.
    pub quality_bases: u64,
    /// Those bases' qualities, added up.
    pub quality_sum: u64,
    /// Secondary alignments.
    pub secondary: u64,
    /// Supplementary alignments.
    pub supplementary: u64,
    /// Unmapped reads placed at a position (a mate's), which the pass reads with the contig's.
    pub unmapped: u64,
    /// Whether the reads are those of `--contig` or `--regions` alone (in regions, each read
    /// that reaches them, once when they are sorted and apart), not the whole file's.
    pub restricted: bool,
    /// How many reads have each length.
    lengths: BTreeMap<u64, u64>,
}

impl ReadStats {
    /// Counts one record, with `flags`, its sequence's length and its qualities.
    pub fn add(&mut self, flags: u16, length: usize, qualities: &[u8]) {
        if flags & FLAG_UNMAPPED != 0 {
            self.unmapped += 1;
        } else if flags & FLAG_SECONDARY != 0 {
            self.secondary += 1;
        } else if flags & FLAG_SUPPLEMENTARY != 0 {
            self.supplementary += 1;
        } else {
            let length = length as u64;
            self.reads += 1;
            self.bases += length;
            *self.lengths.entry(length).or_default() += 1;
            if qualities.len() as u64 == length && qualities.first().is_some_and(|&q| q != MISSING_QUALITY) {
                self.quality_bases += length;
                // In u32 within a read's pieces of 2^24 bases (which cannot overflow it), which
                // vectorises wider than u64.
                let piece_sum = |piece: &[u8]| u64::from(piece.iter().map(|&q| u32::from(q)).sum::<u32>());
                self.quality_sum += qualities.chunks(1 << 24).map(piece_sum).sum::<u64>();
            }
        }
    }

    /// Adds `other`'s counts to these.
    pub fn merge(&mut self, other: ReadStats) {
        self.reads += other.reads;
        self.bases += other.bases;
        self.quality_bases += other.quality_bases;
        self.quality_sum += other.quality_sum;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.unmapped += other.unmapped;
        for (length, n) in other.lengths {
            *self.lengths.entry(length).or_default() += n;
        }
    }

    /// The read length N50: the length of the shortest of the longest reads that together hold
    /// at least half the bases; 0 without reads.
    #[must_use]
    pub fn n50(&self) -> u64 {
        let mut held = 0;
        for (&length, &n) in self.lengths.iter().rev() {
            held += length * n;
            if 2 * held >= self.bases {
                return length;
            }
        }
        0
    }

    /// The counts as JSON: the reads (primary alignments of mapped reads), their bases, N50,
    /// mean length and mean base quality (null for none), and the secondary, supplementary and
    /// placed unmapped records, and whether they are of `--contig` or `--regions` alone.
    #[must_use]
    pub fn to_json(&self) -> serde_json::Value {
        let mean = |sum: u64, n: u64| (n > 0).then(|| round(sum as f64 / n as f64));
        serde_json::json!({
            "reads": self.reads,
            "bases": self.bases,
            "n50": self.n50(),
            "mean_length": mean(self.bases, self.reads),
            "mean_quality": mean(self.quality_sum, self.quality_bases),
            "quality_bases": self.quality_bases,
            "secondary": self.secondary,
            "supplementary": self.supplementary,
            "unmapped_placed": self.unmapped,
            "restricted": self.restricted,
        })
    }

    /// Writes [`Self::to_json`] to `path`.
    ///
    /// # Errors
    ///
    /// If `path` cannot be written.
    pub fn write(&self, path: &Path) -> Result<()> {
        let text = serde_json::to_string_pretty(&self.to_json()).expect("JSON of numbers") + "\n";
        std::fs::write(path, text).with_context(|| format!("writing {}", path.display()))
    }
}

/// `x` to three decimals.
fn round(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn of(lengths: &[usize]) -> ReadStats {
        let mut stats = ReadStats::default();
        for &n in lengths {
            stats.add(0, n, &vec![20; n]);
        }
        stats
    }

    #[test]
    fn n50_is_the_shortest_of_the_longest_reads_holding_half_the_bases() {
        // 100 bases: 40 + 30 = 70 hold half, 40 alone does not.
        assert_eq!(of(&[10, 40, 20, 30]).n50(), 30);
        // Exactly half counts: 50 of 100.
        assert_eq!(of(&[50, 25, 25]).n50(), 50);
        assert_eq!(of(&[5, 5, 5, 5]).n50(), 5);
        assert_eq!(of(&[7]).n50(), 7);
        assert_eq!(of(&[]).n50(), 0);
        // Merged pieces give the whole's.
        let mut merged = of(&[10, 40]);
        merged.merge(of(&[20, 30]));
        assert_eq!(merged, of(&[10, 40, 20, 30]));
    }

    #[test]
    fn counts_each_kind_of_record_and_the_qualities_there_are() {
        let mut stats = ReadStats::default();
        stats.add(0, 4, &[10, 20, 30, 40]);
        stats.add(0x10, 2, &[MISSING_QUALITY, MISSING_QUALITY]);
        stats.add(FLAG_SECONDARY, 3, &[]);
        stats.add(FLAG_SUPPLEMENTARY | 0x10, 3, &[1, 1, 1]);
        stats.add(FLAG_UNMAPPED, 9, &[1; 9]);
        let json = stats.to_json();
        assert_eq!(json["reads"], 2);
        assert_eq!(json["bases"], 6);
        assert_eq!(json["mean_length"], 3.0);
        // The second read has no qualities.
        assert_eq!(json["mean_quality"], 25.0);
        assert_eq!(json["quality_bases"], 4);
        assert_eq!((json["secondary"].as_u64(), json["supplementary"].as_u64()), (Some(1), Some(1)));
        assert_eq!(json["unmapped_placed"], 1);
        assert_eq!(ReadStats::default().to_json()["mean_quality"], serde_json::Value::Null);
    }
}

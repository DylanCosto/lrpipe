//! `bed_tree_from` and `is_region_in` (`shared/interval_tree.py`) for point queries, and the
//! `pigz -fdc` stream the original reads its inputs through.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::RowError;
use crate::pytext::{TextError, TextLines, parse_int, split_whitespace, strip};

/// Opens `path` as `pigz -fdc` streams it: gzip (any number of members) is decompressed and
/// anything else is passed through unchanged.
///
/// An error here is where `pigz` fails before writing anything; the original never checks its exit
/// status, so callers read such a file as empty.
pub(crate) fn open_pigz(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut reader = BufReader::with_capacity(1 << 20, File::open(path)?);
    let is_gzip = reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    Ok(if is_gzip {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(reader)))
    } else {
        Box::new(reader)
    })
}

/// The intervals of one contig, sorted by start, with the running maximum end so a point
/// query is a binary search.
#[derive(Debug, Clone, Default)]
pub(crate) struct Intervals {
    /// `(start, end)` sorted by start.
    rows: Vec<(i64, i64)>,
    /// `max_end[i]` is the largest end among `rows[..=i]`.
    max_end: Vec<i64>,
}

impl Intervals {
    pub(crate) fn new(mut rows: Vec<(i64, i64)>) -> Self {
        rows.sort_unstable();
        let max_end = rows
            .iter()
            .scan(i64::MIN, |max, &(_, end)| {
                *max = (*max).max(end);
                Some(*max)
            })
            .collect();
        Self { rows, max_end }
    }

    /// `IntervalTree.at(pos)` is non-empty: some interval has `start <= pos < end`.
    fn contains(&self, pos: i64) -> bool {
        let starts_at_or_before = self.rows.partition_point(|&(start, _)| start <= pos);
        starts_at_or_before > 0 && self.max_end[starts_at_or_before - 1] > pos
    }
}

/// Why [`BedTree::read`] stopped.
#[derive(Debug)]
pub(crate) enum BedError {
    /// Reading the stream failed.
    Io(io::Error),
    /// A row is rejected where the original raises or exits.
    Row {
        /// 1-based line number.
        line: usize,
        /// What is wrong with it.
        reason: RowError,
    },
}

/// The tree `bed_tree_from` returns: 0-based half-open intervals per contig.
#[derive(Debug, Clone, Default)]
pub(crate) struct BedTree {
    contigs: HashMap<String, Intervals>,
}

impl BedTree {
    /// `bed_tree_from(bed, contig_name=contig, padding=padding)` over the lines of `reader`.
    ///
    /// Rows of other contigs are skipped when `contig` is given. A non-zero `padding` (zero is
    /// falsy in Python) shrinks each interval by that much at both ends after validation; an
    /// interval that ends up empty-or-reversed is rejected as `IntervalTree` does.
    pub(crate) fn read(reader: impl BufRead, contig: Option<&str>, padding: i64) -> Result<Self, BedError> {
        let mut rows: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (row_id, line) in TextLines::new(reader).enumerate() {
            let bad = |reason| BedError::Row { line: row_id + 1, reason };
            let line = match line {
                Ok(line) => line,
                Err(TextError::Io(source)) => return Err(BedError::Io(source)),
                Err(TextError::Utf8) => return Err(bad(RowError::InvalidUtf8)),
            };
            if line.starts_with('#') {
                continue;
            }
            let columns = split_whitespace(strip(&line), None);
            let ctg_name = *columns.first().ok_or(bad(RowError::MissingColumn { index: 0 }))?;
            if contig.is_some_and(|contig| contig != ctg_name) {
                continue;
            }
            // The original creates the contig's tree before parsing the coordinates.
            let contig_rows = rows.entry(ctg_name.to_owned()).or_default();
            let coordinate = |index: usize| {
                let column = columns.get(index).ok_or(bad(RowError::MissingColumn { index }))?;
                parse_int(column).ok_or_else(|| bad(RowError::InvalidInteger { value: (*column).to_owned() }))
            };
            let (mut start, mut end) = (coordinate(1)?, coordinate(2)?);
            if end < start || start < 0 || end < 0 {
                return Err(bad(RowError::InvalidInterval { start, end }));
            }
            if padding != 0 {
                start = start.saturating_add(padding);
                end = end.saturating_sub(padding);
            }
            if start == end {
                end += 1;
            }
            if start >= end {
                // `IntervalTree.addi` raises ValueError on a null interval.
                return Err(bad(RowError::NullInterval { start, end }));
            }
            contig_rows.push((start, end));
        }
        let contigs = rows.into_iter().map(|(contig, rows)| (contig, Intervals::new(rows))).collect();
        Ok(Self { contigs })
    }

    /// `is_region_in(tree, contig, pos)`.
    pub(crate) fn contains(&self, contig: &str, pos: i64) -> bool {
        self.contigs.get(contig).is_some_and(|intervals| intervals.contains(pos))
    }

    /// The next position where interval membership changes, including across overlapping rows.
    pub(crate) fn next_boundary(&self, contig: &str, pos: i64) -> Option<i64> {
        let intervals = self.contigs.get(contig)?;
        let after = intervals.rows.partition_point(|&(start, _)| start <= pos);
        if let Some(mut end) = after.checked_sub(1).map(|i| intervals.max_end[i]).filter(|&end| end > pos) {
            for &(start, stop) in &intervals.rows[after..] {
                if start > end {
                    break;
                }
                end = end.max(stop);
            }
            Some(end)
        } else {
            intervals.rows.get(after).map(|&(start, _)| start)
        }
    }

    /// `not tree`: no row named a (selected) contig.
    pub(crate) fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_queries() {
        let iv = Intervals::new(vec![(10, 100), (20, 30), (150, 151)]);
        for (pos, want) in [(9, false), (10, true), (35, true), (99, true), (100, false), (150, true), (151, false)] {
            assert_eq!(iv.contains(pos), want, "{pos}");
        }
    }

    #[test]
    fn reads_like_bed_tree_from() {
        let bed = "#c\nchr1\t10\t20\nchr2 5 5\nchr1\t30\t100\n";
        let tree = BedTree::read(bed.as_bytes(), None, 0).expect("valid BED");
        assert!(tree.contains("chr1", 10) && !tree.contains("chr1", 20) && tree.contains("chr2", 5));
        let only = BedTree::read(bed.as_bytes(), Some("chr2"), 0).expect("valid BED");
        assert!(!only.contains("chr1", 10) && only.contains("chr2", 5));
        // Padding 33 leaves [63, 67) of [30, 100); [10, 20) would be reversed.
        let padded = BedTree::read("chr1\t30\t100\nchr1\t40\t106\n".as_bytes(), None, 33).expect("valid BED");
        assert!(!padded.contains("chr1", 62) && padded.contains("chr1", 63) && padded.contains("chr1", 66));
        assert!(!padded.contains("chr1", 67));
        assert!(matches!(
            BedTree::read(bed.as_bytes(), None, 33),
            Err(BedError::Row { line: 2, reason: RowError::NullInterval { start: 43, end: -13 } })
        ));
        assert!(matches!(
            BedTree::read("chr1\t9\t3\n".as_bytes(), None, 0),
            Err(BedError::Row { line: 1, reason: RowError::InvalidInterval { start: 9, end: 3 } })
        ));
        assert!(matches!(
            BedTree::read("\n".as_bytes(), None, 0),
            Err(BedError::Row { line: 1, reason: RowError::MissingColumn { index: 0 } })
        ));
        assert!(BedTree::read("chr9\t1\t2\n".as_bytes(), Some("chr1"), 0).expect("valid BED").is_empty());
    }
}

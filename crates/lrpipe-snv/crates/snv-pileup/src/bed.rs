//! The BED handling of `shared/interval_tree.py::bed_tree_from` for one contig.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

use crate::Error;

/// A half-open interval `[start, end)`, 0-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    /// First position.
    pub start: i64,
    /// One past the last position.
    pub end: i64,
}

impl Interval {
    /// Whether this interval overlaps `[start, end)`, as `IntervalTree.overlap` tests it.
    pub fn overlaps(&self, start: i64, end: i64) -> bool {
        self.start < end && start < self.end
    }
}

/// The rows of a BED file on one contig, with their raw coordinates.
#[derive(Debug, Clone, Default)]
pub struct Bed {
    rows: Vec<(i64, i64)>,
}

impl Bed {
    /// Reads the rows of `contig` from a plain or gzipped BED file (`pigz -fdc` reads both).
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read or a row is malformed.
    pub fn read_contig(path: &Path, contig: &str) -> Result<Self, Error> {
        let io_err = |source| Error::Io { path: path.to_owned(), source };
        let mut file = File::open(path).map_err(io_err)?;
        let mut magic = [0u8; 2];
        let n = file.read(&mut magic).map_err(io_err)?;
        let file = File::open(path).map_err(io_err)?;
        let reader: Box<dyn BufRead> = if n == 2 && magic == [0x1f, 0x8b] {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        let mut rows = Vec::new();
        for (row_id, line) in reader.lines().enumerate() {
            let line = line.map_err(io_err)?;
            if line.starts_with('#') {
                continue;
            }
            let mut cols = line.split_whitespace();
            if cols.next() != Some(contig) {
                continue;
            }
            let bad = || Error::Bed { path: PathBuf::from(path), line: row_id + 1 };
            let mut coord = || cols.next().and_then(|c| c.parse::<i64>().ok()).ok_or_else(bad);
            let (start, end) = (coord()?, coord()?);
            if end < start || start < 0 {
                return Err(bad());
            }
            rows.push((start, end));
        }
        Ok(Self { rows })
    }

    /// `(bed_start, bed_end)`: the smallest start and largest end.
    pub fn span(&self) -> Option<(i64, i64)> {
        let start = self.rows.iter().map(|r| r.0).min()?;
        let end = self.rows.iter().map(|r| r.1).max()?;
        Some((start, end))
    }

    /// The tree `bed_tree_from(..., bed_ctg_start, bed_ctg_end)` builds: rows touching
    /// `[range_start, range_end]`, with empty rows widened to one base.
    pub fn overlapping(&self, range_start: i64, range_end: i64) -> Vec<Interval> {
        let filter = range_start != 0 && range_end != 0;
        self.rows
            .iter()
            .filter(|(s, e)| !filter || !(*e < range_start || *s > range_end))
            .map(|&(start, end)| Interval { start, end: if start == end { end + 1 } else { end } })
            .collect()
    }
}

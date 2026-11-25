//! BED regions as samtools reads them (`bedidx.c`) for `mpileup -l` and `view -L`.

use std::collections::HashMap;
use std::path::Path;

use crate::Error;

/// Half-open 0-based intervals per contig, as `bed_read` keeps them.
#[derive(Debug, Clone, Default)]
pub struct Regions {
    by_contig: HashMap<String, Vec<(i64, i64)>>,
    /// Contig names by target id, for lookups by `tid`.
    names: Vec<String>,
}

impl Regions {
    /// Reads a plain or gzipped BED file. Lines with one coordinate are 1-based positions
    /// (`end = beg--`); `track` and `browser` lines are skipped.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read or a line cannot be parsed.
    pub fn read(path: &Path) -> Result<Self, Error> {
        let text = crate::read_text(path)?;
        let mut by_contig: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (index, line) in text.lines().enumerate() {
            let line = line.trim_start();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split_whitespace();
            let name = cols.next().unwrap_or_default();
            let nums: Vec<i64> = cols.take(2).map_while(|c| c.parse().ok()).collect();
            let (beg, end) = match nums[..] {
                [beg, end] => (beg, end),
                [pos] => (pos - 1, pos),
                _ => (0, -1),
            };
            if nums.is_empty() || end < beg {
                if name == "browser" || name == "track" {
                    continue;
                }
                return Err(Error::Bed { path: path.to_owned(), line: index + 1 });
            }
            by_contig.entry(name.to_owned()).or_default().push((beg, end));
        }
        for intervals in by_contig.values_mut() {
            intervals.sort_unstable();
        }
        Ok(Self { by_contig, names: Vec::new() })
    }

    /// Resolves target ids against `names` (the header's contig names, in order).
    pub fn with_target_names(mut self, names: Vec<String>) -> Self {
        self.names = names;
        self
    }

    /// `bed_overlap(bed, chr, beg, end)`: whether any interval overlaps `[beg, end)` on target
    /// `tid`.
    pub fn overlaps(&self, tid: i32, beg: i64, end: i64) -> bool {
        let Some(name) = usize::try_from(tid).ok().and_then(|t| self.names.get(t)) else { return false };
        self.overlaps_contig(name, beg, end)
    }

    /// The same by contig name.
    pub fn overlaps_contig(&self, contig: &str, beg: i64, end: i64) -> bool {
        let Some(intervals) = self.by_contig.get(contig) else { return false };
        intervals.iter().take_while(|&&(b, _)| b < end).any(|&(b, e)| e > beg && b < end)
    }
}

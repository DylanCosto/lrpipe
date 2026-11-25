//! Regions of the genome from a BED file: each line's contig, 0-based start and exclusive end,
//! plain or gzipped. The intervals of each contig are sorted, and those that overlap or touch
//! are merged (as `bedtools merge` merges them), so that each base is in at most one interval
//! and the intervals of a contig are in order; the contigs keep the order the file first names
//! them in. Tools restricted to regions read them with [`Regions::read`] and ask
//! [`Regions::overlaps`] of what they would report.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Why a BED file of regions could not be read.
#[derive(Debug, thiserror::Error)]
pub enum RegionsError {
    /// The file could not be opened or read.
    #[error("reading {}", path.display())]
    Io {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A line is not a contig, start and end.
    #[error("{}, line {line}: {reason}", path.display())]
    Format {
        /// The file.
        path: PathBuf,
        /// Its line, from 1.
        line: usize,
        /// What is wrong with it.
        reason: String,
    },
}

/// The regions of a BED file, by contig.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Regions {
    /// Each contig in the order first named, with its intervals sorted and merged.
    contigs: Vec<(String, Vec<(i64, i64)>)>,
}

impl Regions {
    /// Reads a BED file, gzipped (told apart by its first bytes) or not, with or without a UTF-8
    /// byte order mark. Empty lines and `#`, `track` and `browser` lines are skipped; columns
    /// after the third are ignored, and the columns may be separated by spaces as well as tabs.
    /// A zero-length row (`start == end`, an insertion point) is the base at its start, as
    /// Clair3 reads it.
    ///
    /// # Errors
    ///
    /// If the file cannot be read, or a line does not have a contig, a start and an end that are
    /// integers with `0 <= start <= end`.
    pub fn read(path: &Path) -> Result<Self, RegionsError> {
        let io_error = |source| RegionsError::Io { path: path.into(), source };
        let mut file = BufReader::new(File::open(path).map_err(io_error)?);
        let gzipped = file.fill_buf().map_err(io_error)?.starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> =
            if gzipped { Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))) } else { Box::new(file) };
        let mut rows = Vec::new();
        for (n, line) in reader.lines().enumerate() {
            let line = line.map_err(io_error)?;
            let line = if n == 0 { line.strip_prefix('\u{feff}').unwrap_or(&line) } else { &line };
            let fail = |reason: &str| RegionsError::Format { path: path.into(), line: n + 1, reason: reason.into() };
            let mut columns = line.split_whitespace();
            let contig = match columns.next() {
                None | Some("track" | "browser") => continue,
                Some(c) if c.starts_with('#') => continue,
                Some(c) => c,
            };
            let (Some(start), Some(end)) = (columns.next(), columns.next()) else {
                return Err(fail("expected a contig, a start and an end"));
            };
            let (Ok(start), Ok(end)) = (start.parse::<i64>(), end.parse::<i64>()) else {
                return Err(fail("the start and end must be integers"));
            };
            if start < 0 || end < start {
                return Err(fail(&format!("the interval {start}-{end} ends before it starts, or is negative")));
            }
            rows.push((contig.to_owned(), start, end.max(start + 1)));
        }
        Ok(Self::from_intervals(rows))
    }

    /// The regions of `intervals` (contig, 0-based start, exclusive end), sorted and merged.
    #[must_use]
    pub fn from_intervals(intervals: impl IntoIterator<Item = (String, i64, i64)>) -> Self {
        let mut contigs: Vec<(String, Vec<(i64, i64)>)> = Vec::new();
        for (contig, start, end) in intervals {
            match contigs.iter_mut().rev().find(|(c, _)| *c == contig) {
                Some((_, v)) => v.push((start, end)),
                None => contigs.push((contig, vec![(start, end)])),
            }
        }
        for (_, intervals) in &mut contigs {
            intervals.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::with_capacity(intervals.len());
            for &(start, end) in intervals.iter() {
                match merged.last_mut() {
                    Some(last) if start <= last.1 => last.1 = last.1.max(end),
                    _ => merged.push((start, end)),
                }
            }
            *intervals = merged;
        }
        Self { contigs }
    }

    /// The same regions with their contigs in the order of `names` (a header's or a reference's),
    /// those it does not name after them in the order first named, as tools write contigs.
    #[must_use]
    pub fn in_order<'a>(&self, names: impl IntoIterator<Item = &'a str>) -> Self {
        let mut contigs = Vec::with_capacity(self.contigs.len());
        for name in names {
            if let Some(c) = self.contigs.iter().find(|(c, _)| c == name)
                && !contigs.iter().any(|(n, _): &(String, Vec<(i64, i64)>)| n == name)
            {
                contigs.push(c.clone());
            }
        }
        let rest: Vec<_> = self.contigs.iter().filter(|(c, _)| !contigs.iter().any(|(n, _)| n == c)).cloned().collect();
        contigs.extend(rest);
        Self { contigs }
    }

    /// Each contig, in the order first named, with its intervals in order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &[(i64, i64)])> {
        self.contigs.iter().map(|(c, v)| (c.as_str(), v.as_slice()))
    }

    /// The intervals of `contig`, in order (none if the regions do not name it).
    #[must_use]
    pub fn contig(&self, contig: &str) -> &[(i64, i64)] {
        self.contigs.iter().find(|(c, _)| c == contig).map_or(&[], |(_, v)| v.as_slice())
    }

    /// Whether the regions name `contig`.
    #[must_use]
    pub fn has(&self, contig: &str) -> bool {
        self.contigs.iter().any(|(c, _)| c == contig)
    }

    /// Whether `[start, end)` of `contig` overlaps a region: shares a base with one. An empty
    /// span (`start == end`) is taken as the base at `start`.
    #[must_use]
    pub fn overlaps(&self, contig: &str, start: i64, end: i64) -> bool {
        let intervals = self.contig(contig);
        let end = end.max(start + 1);
        // The first interval ending after `start`.
        let i = intervals.partition_point(|&(_, e)| e <= start);
        intervals.get(i).is_some_and(|&(s, _)| s < end)
    }

    /// The number of intervals.
    #[must_use]
    pub fn len(&self) -> usize {
        self.contigs.iter().map(|(_, v)| v.len()).sum()
    }

    /// Whether there are none.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The bases the regions cover.
    #[must_use]
    pub fn bases(&self) -> u64 {
        self.contigs.iter().flat_map(|(_, v)| v).map(|&(s, e)| (e - s).unsigned_abs()).sum()
    }

    /// The regions as a BED file's text: `contig`, start and end, tab-separated, a line each.
    #[must_use]
    pub fn bed(&self) -> String {
        let mut out = String::new();
        for (contig, intervals) in &self.contigs {
            for (start, end) in intervals {
                out.push_str(&format!("{contig}\t{start}\t{end}\n"));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temp(name: &str, text: &[u8]) -> PathBuf {
        let path = std::env::temp_dir().join(format!("lrpipe-reads-regions-{}-{name}", std::process::id()));
        File::create(&path).and_then(|mut f| f.write_all(text)).expect("written");
        path
    }

    #[test]
    fn intervals_are_sorted_and_merged_per_contig_in_the_files_order() {
        let text = "\u{feff}track name=x\n# a comment\nchr2\t50\t60\tname\nchr1\t30\t40\nchr2 10 20\nchr1\t0\t10\n\
                    chr1\t10\t15\nchr1\t35\t45\nchr2\t55\t58\nchr3\t7\t7\n";
        let path = temp("merge.bed", text.as_bytes());
        let regions = Regions::read(&path).expect("read");
        // A byte order mark is not the first contig's name, and an insertion point is a base.
        assert_eq!(regions.bed(), "chr2\t10\t20\nchr2\t50\t60\nchr1\t0\t15\nchr1\t30\t45\nchr3\t7\t8\n");
        assert_eq!((regions.len(), regions.bases()), (5, 51));
        let bom = temp("bom.bed", "\u{feff}chr1\t0\t5\n".as_bytes());
        assert_eq!(Regions::read(&bom).expect("read").bed(), "chr1\t0\t5\n");
        let _ = std::fs::remove_file(bom);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn contigs_follow_a_headers_order() {
        let regions = Regions::from_intervals([("chr2".into(), 1, 2), ("chrZ".into(), 0, 1), ("chr1".into(), 3, 4)]);
        let ordered = regions.in_order(["chr1", "chr2", "chr3"]);
        assert_eq!(ordered.bed(), "chr1\t3\t4\nchr2\t1\t2\nchrZ\t0\t1\n");
    }

    #[test]
    fn overlaps_share_a_base() {
        let regions = Regions::from_intervals([("c".into(), 10, 20), ("c".into(), 30, 40)]);
        for (start, end, overlaps) in [
            (0, 10, false),
            (0, 11, true),
            (19, 25, true),
            (20, 30, false),
            (25, 35, true),
            (39, 39, true),
            (40, 41, false),
            (0, 100, true),
        ] {
            assert_eq!(regions.overlaps("c", start, end), overlaps, "{start}-{end}");
        }
        assert!(!regions.overlaps("d", 0, 100));
    }

    #[test]
    fn gzipped_files_are_read_and_bad_lines_named() {
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gz.write_all(b"chr1\t5\t9\n").expect("written");
        let path = temp("gz.bed.gz", &gz.finish().expect("finished"));
        assert_eq!(Regions::read(&path).expect("read").bed(), "chr1\t5\t9\n");
        for (text, reason) in [
            ("chr1\t5\n", "a contig, a start and an end"),
            ("chr1\t5\tx\n", "integers"),
            ("chr1\t5\t4\n", "ends before it starts"),
            ("chr1\t-1\t5\n", "negative"),
        ] {
            let path = temp("bad.bed", format!("chr1\t0\t1\n{text}").as_bytes());
            match Regions::read(&path) {
                Err(RegionsError::Format { line: 2, reason: r, .. }) => assert!(r.contains(reason), "{r}"),
                other => panic!("{text}: {other:?}"),
            }
        }
        assert!(matches!(Regions::read(Path::new("/nonexistent/x.bed")), Err(RegionsError::Io { .. })));
        let _ = std::fs::remove_file(path);
    }
}

//! Regions to report depth over: fixed windows or BED intervals, as mosdepth makes them.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};

/// A 0-based half-open interval on a contig, with an optional BED name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    pub start: u32,
    pub stop: u32,
    pub name: Option<String>,
}

/// mosdepth's `window_gen` (mosdepth.nim:382-388): `window`-sized intervals, the last one
/// shortened to end at the contig end.
pub fn windows(window: u32, contig_len: u32) -> impl Iterator<Item = Region> {
    let mut start = 0u64;
    let (window, len) = (u64::from(window), u64::from(contig_len));
    std::iter::from_fn(move || {
        if start + window < len {
            let r = Region { start: start as u32, stop: (start + window) as u32, name: None };
            start += window;
            Some(r)
        } else if start != len {
            let r = Region { start: start as u32, stop: len as u32, name: None };
            start = len;
            Some(r)
        } else {
            None
        }
    })
}

/// BED intervals grouped by contig, each group sorted by start (stable, like Nim's sort).
#[derive(Debug, Default)]
pub struct BedRegions {
    by_contig: HashMap<String, Vec<Region>>,
}

impl BedRegions {
    /// Reads a BED file, plain or gzip/BGZF-compressed (mosdepth.nim:362-380 and 195-209).
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening BED file {}", path.display()))?;
        let mut reader = BufReader::new(file);
        let gz = reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> =
            if gz { Box::new(BufReader::new(flate2::bufread::MultiGzDecoder::new(reader))) } else { Box::new(reader) };
        Self::from_reader(reader).with_context(|| format!("reading BED file {}", path.display()))
    }

    pub fn from_reader(mut reader: impl BufRead) -> Result<Self> {
        let mut by_contig: HashMap<String, Vec<Region>> = HashMap::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            // htslib's hts_getline returns the line length, and mosdepth's loop stops at 0:
            // reading ends at the first empty line.
            if line.is_empty() {
                break;
            }
            if line.starts_with(b"track ") || line[0] == b'#' {
                continue;
            }
            let line = String::from_utf8_lossy(&line);
            let cols: Vec<&str> = line.splitn(6, '\t').collect();
            if cols.len() < 3 {
                eprintln!("[lrpipe-depth] skipping bad bed line:{}", line.trim());
                continue;
            }
            let parse = |s: &str| -> Result<u32> {
                s.parse::<u32>().with_context(|| format!("invalid BED coordinate '{s}' in line: {line}"))
            };
            let (start, stop) = (parse(cols[1])?, parse(cols[2])?);
            if start > stop {
                bail!("start > end in bed line:{line}");
            }
            let name = cols.get(3).map(|n| n.trim()).filter(|n| !n.is_empty()).map(str::to_owned);
            by_contig.entry(cols[0].to_owned()).or_default().push(Region { start, stop, name });
        }
        for regions in by_contig.values_mut() {
            regions.sort_by_key(|r| r.start);
        }
        Ok(BedRegions { by_contig })
    }

    pub fn contains(&self, contig: &str) -> bool {
        self.by_contig.contains_key(contig)
    }

    pub fn get(&self, contig: &str) -> &[Region] {
        self.by_contig.get(contig).map_or(&[], Vec::as_slice)
    }

    /// Contigs and their region counts, for contigs not in `known`.
    pub fn missing<'a>(&'a self, known: &'a [String]) -> impl Iterator<Item = (&'a str, usize)> + 'a {
        self.by_contig.iter().filter(|(c, _)| !known.contains(c)).map(|(c, r)| (c.as_str(), r.len()))
    }
}

/// What to report depth over (mosdepth `--by`).
#[derive(Debug)]
pub enum By {
    Windows(u32),
    Bed(BedRegions),
}

impl By {
    /// An all-digit argument is a window size, anything else a BED path (mosdepth.nim:679).
    pub fn parse(arg: &str) -> Result<Self> {
        if !arg.is_empty() && arg.bytes().all(|b| b.is_ascii_digit()) {
            let w: u32 = arg.parse().context("window size out of range")?;
            if w == 0 {
                bail!("window size must be positive");
            }
            Ok(By::Windows(w))
        } else {
            Ok(By::Bed(BedRegions::from_path(Path::new(arg))?))
        }
    }

    pub fn is_windows(&self) -> bool {
        matches!(self, By::Windows(_))
    }

    /// The regions on one contig.
    pub fn regions(&self, contig: &str, contig_len: u32) -> Box<dyn Iterator<Item = Region> + '_> {
        match self {
            By::Windows(w) => Box::new(windows(*w, contig_len)),
            By::Bed(bed) => Box::new(bed.get(contig).iter().cloned()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(w: u32, len: u32) -> Vec<(u32, u32)> {
        windows(w, len).map(|r| (r.start, r.stop)).collect()
    }

    #[test]
    fn windows_like_mosdepth() {
        assert_eq!(spans(10, 25), [(0, 10), (10, 20), (20, 25)]);
        assert_eq!(spans(10, 20), [(0, 10), (10, 20)]);
        assert_eq!(spans(10, 5), [(0, 5)]);
        assert_eq!(spans(10, 0), []);
    }

    #[test]
    fn bed_parsing_quirks() {
        let bed = b"track name=x\n#c\nchr1\t20\t30\tb \nchr1\t5\t10\nchr1\t5\t8\t \tx\nchr2\t1\n\nchr3\t0\t1\n";
        let r = BedRegions::from_reader(&bed[..]).unwrap();
        let c1 = r.get("chr1");
        assert_eq!(c1.len(), 3);
        assert_eq!((c1[0].start, c1[0].stop, c1[0].name.as_deref()), (5, 10, None));
        assert_eq!((c1[1].start, c1[1].stop, c1[1].name.as_deref()), (5, 8, None));
        assert_eq!(c1[2].name.as_deref(), Some("b"));
        assert!(!r.contains("chr2"));
        // Reading stops at the empty line.
        assert!(!r.contains("chr3"));
    }
}

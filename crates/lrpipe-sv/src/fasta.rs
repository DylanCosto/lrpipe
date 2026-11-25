//! `.fai`-indexed FASTA access through positioned reads: plain, or compressed with bgzip and read
//! through its `.gzi` (told apart by contents, as lrpipe-reads and htslib do).

use std::collections::HashMap;
use std::io;
use std::path::Path;

use lrpipe_reads::FastaFile;

use anyhow::{Context, Result, bail};

/// Bytes read at a time when scanning a contig.
const SCAN_BLOCK: u64 = 1 << 20;

#[derive(Debug, Clone, Copy)]
struct FaiEntry {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

#[derive(Debug)]
pub struct Fasta {
    file: FastaFile,
    index: HashMap<String, FaiEntry>,
}

impl Fasta {
    pub fn open(path: &str) -> Result<Self> {
        let fai = std::fs::read_to_string(format!("{path}.fai"))
            .with_context(|| format!("reading {path}.fai (index the reference with samtools faidx)"))?;
        let mut index = HashMap::new();
        for line in fai.lines().filter(|l| !l.is_empty()) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                bail!("malformed .fai line: {line}");
            }
            index.insert(
                f[0].to_string(),
                FaiEntry {
                    len: f[1].parse()?,
                    offset: f[2].parse()?,
                    line_bases: f[3].parse()?,
                    line_width: f[4].parse()?,
                },
            );
        }
        let file = FastaFile::open(Path::new(path)).with_context(|| format!("opening the reference {path}"))?;
        lrpipe_reads::check_fai(Path::new(path))
            .map_err(lrpipe_reads::OpenError::into_io_error)
            .with_context(|| format!("checking the reference {path}"))?;
        Ok(Fasta { file, index })
    }

    /// One contig of the reference, or `None` if there is no such contig.
    pub fn contig(&self, name: &str) -> Option<Contig<'_>> {
        self.index.get(name).map(|&e| Contig { file: &self.file, e })
    }

    /// Checks the reference against the contigs SVs are called on (`source`'s, with their
    /// lengths): a contig of another length is another assembly, whose bases would make every
    /// REF and insertion anchor wrong (or empty past its end), so it is refused, as Sniffles2
    /// fails reading past it. Returns the contigs the reference lacks, whose SVs Sniffles2 writes
    /// with an `N` REF, for the caller to warn of.
    ///
    /// # Errors
    ///
    /// Fails on a contig whose length differs from the reference's.
    pub fn check_contigs<'c>(&self, contigs: &'c [(String, i64)], source: &str) -> Result<Vec<&'c str>> {
        let mut missing = Vec::new();
        for (name, len) in contigs {
            match self.contig(name) {
                Some(contig) if contig.len() != *len => bail!(
                    "contig {name} is {len} bp in {source} but {} bp in the reference: give the reference the reads \
                     were aligned to",
                    contig.len()
                ),
                Some(_) => {}
                None => missing.push(name.as_str()),
            }
        }
        Ok(missing)
    }
}

/// A warning for SVs called on `missing`, contigs the reference lacks, or `None`.
pub fn missing_contigs_warning(missing: &[&str]) -> Option<String> {
    let (first, rest) = missing.split_first()?;
    let more = if rest.is_empty() { String::new() } else { format!(" and {} more", rest.len()) };
    Some(format!(
        "the reference has no contig {first}{more}: SVs there are written with N as REF (does the reference \
         name contigs as the reads do, chr20 or 20?)"
    ))
}

/// One contig of an indexed FASTA.
#[derive(Debug, Clone, Copy)]
pub struct Contig<'a> {
    file: &'a FastaFile,
    e: FaiEntry,
}

impl Contig<'_> {
    pub fn len(&self) -> i64 {
        self.e.len as i64
    }

    /// Bases [start, end), clipped to the contig (`pysam.FastaFile.fetch`).
    pub fn fetch(&self, start: i64, end: i64) -> io::Result<Vec<u8>> {
        let e = self.e;
        let start = start.clamp(0, e.len as i64) as u64;
        let end = end.clamp(start as i64, e.len as i64) as u64;
        if start == end {
            return Ok(Vec::new());
        }
        let file_start = e.offset + start / e.line_bases * e.line_width + start % e.line_bases;
        let last = end - 1;
        let file_end = e.offset + last / e.line_bases * e.line_width + last % e.line_bases + 1;
        let mut raw = vec![0; (file_end - file_start) as usize];
        self.file.read_exact_at(&mut raw, file_start)?;
        raw.retain(|&b| b != b'\n' && b != b'\r');
        Ok(raw)
    }

    /// Base position of a file offset inside this contig's sequence lines.
    fn base_pos(&self, file_off: u64) -> i64 {
        let rel = file_off - self.e.offset;
        (rel / self.e.line_width * self.e.line_bases + rel % self.e.line_width) as i64
    }

    /// All runs of `N` bases, found by scanning the raw file in blocks.
    pub fn n_runs(&self) -> io::Result<NRuns> {
        let e = self.e;
        let file_end = if e.len == 0 {
            e.offset
        } else {
            let last = e.len - 1;
            e.offset + last / e.line_bases * e.line_width + last % e.line_bases + 1
        };
        let mut runs = Vec::new();
        let mut run_start: Option<i64> = None;
        let mut buf = vec![0; SCAN_BLOCK as usize];
        let mut off = e.offset;
        while off < file_end {
            let n = ((file_end - off) as usize).min(buf.len());
            let block = &mut buf[..n];
            self.file.read_exact_at(block, off)?;
            let mut k = 0;
            while k < n {
                match run_start {
                    // Inside a run: it continues over N bases and line breaks.
                    Some(s) => match block[k..].iter().position(|&b| !matches!(b, b'N' | b'\n' | b'\r')) {
                        Some(j) => {
                            runs.push((s, self.base_pos(off + (k + j) as u64)));
                            run_start = None;
                            k += j + 1;
                        }
                        None => k = n,
                    },
                    None => match memchr::memchr(b'N', &block[k..]) {
                        Some(j) => {
                            run_start = Some(self.base_pos(off + (k + j) as u64));
                            k += j + 1;
                        }
                        None => k = n,
                    },
                }
            }
            off += n as u64;
        }
        if let Some(s) = run_start {
            runs.push((s, self.len()));
        }
        Ok(NRuns { runs })
    }
}

/// Maximal runs of `N` bases on a contig, sorted.
#[derive(Debug, Default)]
pub struct NRuns {
    runs: Vec<(i64, i64)>,
}

impl NRuns {
    /// Whether the base at `pos` is `N`.
    pub fn is_n(&self, pos: i64) -> bool {
        let k = self.runs.partition_point(|&(s, _)| s <= pos);
        k > 0 && pos < self.runs[k - 1].1
    }

    /// Bases in [start, end) that lie in runs of at least `min_len` `N`s.
    pub fn long_overlap(&self, start: i64, end: i64, min_len: i64) -> i64 {
        let first = self.runs.partition_point(|&(_, e)| e <= start);
        self.runs[first..]
            .iter()
            .take_while(|&&(s, _)| s < end)
            .filter(|&&(s, e)| e - s >= min_len)
            .map(|&(s, e)| e.min(end) - s.max(start))
            .sum()
    }

    /// The parts of these runs inside the given regions (which may overlap).
    pub fn clipped_to(&self, regions: &[(i64, i64)]) -> NRuns {
        let mut runs: Vec<(i64, i64)> = Vec::new();
        for &(rs, re) in regions {
            let first = self.runs.partition_point(|&(_, e)| e <= rs);
            runs.extend(self.runs[first..].iter().take_while(|&&(s, _)| s < re).map(|&(s, e)| (s.max(rs), e.min(re))));
        }
        runs.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::with_capacity(runs.len());
        for (s, e) in runs {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        NRuns { runs: merged }
    }

    pub fn iter(&self) -> impl Iterator<Item = (i64, i64)> + '_ {
        self.runs.iter().copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_reference_of_other_contig_lengths_is_refused() {
        let path = std::env::temp_dir().join(format!("lrpipe-sv-contigs-{}.fa", std::process::id()));
        std::fs::write(&path, ">chrA\nACGTACGT\n").expect("FASTA");
        std::fs::write(path.with_extension("fa.fai"), "chrA\t8\t6\t8\t9\n").expect(".fai");
        let fasta = Fasta::open(path.to_str().expect("UTF-8")).expect("the reference");
        let (name, len) = ("chrA".to_owned(), 8);
        let same = [(name.clone(), len)];
        assert!(fasta.check_contigs(&same, "the BAM").expect("matching lengths").is_empty());
        let longer = [(name.clone(), len + 1)];
        let error = fasta.check_contigs(&longer, "the BAM").expect_err("another length").to_string();
        assert!(error.contains(&format!("contig {name} is {} bp in the BAM but {len} bp", len + 1)), "{error}");
        let absent = [("chrNone".to_owned(), 5)];
        assert_eq!(fasta.check_contigs(&absent, "the BAM").expect("absent is no error"), ["chrNone"]);
        assert!(missing_contigs_warning(&["chrNone"]).expect("a warning").contains("no contig chrNone"));
        assert!(missing_contigs_warning(&[]).is_none());
        let _ = std::fs::remove_file(path.with_extension("fa.fai"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn finds_n_runs_and_overlaps() {
        let runs = NRuns { runs: vec![(10, 20), (100, 250)] };
        assert!(!runs.is_n(9));
        assert!(runs.is_n(10));
        assert!(runs.is_n(19));
        assert!(!runs.is_n(20));
        assert!(runs.is_n(200));
        // Only the 150-bp run counts at a 100-bp minimum.
        assert_eq!(runs.long_overlap(0, 1000, 100), 150);
        assert_eq!(runs.long_overlap(0, 1000, 5), 160);
        assert_eq!(runs.long_overlap(15, 120, 5), 5 + 20);
    }
}

//! Methylation counted from the reads of another pass over the alignments (`lrpipe sv`'s read
//! pass), so the reads are decoded once for both.
//!
//! The pass reads a contig in pieces, each piece the reads that start in it, in parallel. A
//! piece counts its reads' calls at the contig's CpG sites alone ([`Contig::piece`]), the
//! pieces are merged ([`Contig::merge`]), and the contig's table is written window by window
//! ([`Counted::write_window`]) as [`run::run`](crate::run::run) writes it: each window with the
//! rows its reads bring, in the order they bring them.

use std::collections::BTreeMap;
use std::io;
use std::sync::Mutex;

use pasta::bam;

use crate::input::{Converted, Fetched, ReadFilter, convert};
use crate::pileup::{Calling, Counts, Kinds, Motif, Partition, RefWindow, Rows, add_into, tally};
use crate::read::ModRead;
use crate::run::{Job, Report};

/// Why a job's table cannot be counted within another pass.
#[derive(Debug, thiserror::Error)]
pub enum PassError {
    /// Only one bedMethyl table of every read, at CpG sites, is counted so.
    #[error("only one bedMethyl table of every read, at CpG sites (--cpg), can be counted in a pass")]
    Unsupported,
}

/// How a pass counts and writes a bedMethyl table: a [`Job`]'s.
#[derive(Clone, Debug)]
pub struct Pass {
    calling: Calling,
    rows: Rows,
    filter: ReadFilter,
    window: i64,
    overhang: i64,
}

impl Pass {
    /// The pass that counts `job`'s table, on the contigs the pass reads whole (the job's
    /// regions are left to it).
    ///
    /// # Errors
    ///
    /// Unless the job writes one bedMethyl table of every read with CpG rows.
    pub fn new(job: &Job) -> Result<Self, PassError> {
        let Report::BedMethyl { calling, rows } = &job.report else { return Err(PassError::Unsupported) };
        let one_table = matches!(job.outputs.as_slice(), [(Partition::All, _)]);
        if rows.motif != Motif::CpG || !one_table {
            return Err(PassError::Unsupported);
        }
        Ok(Pass {
            calling: *calling,
            rows: rows.clone(),
            filter: job.filter,
            window: job.window.max(1),
            overhang: job.report.overhang(),
        })
    }

    /// A contig `len` bases long to count, its CpG sites found in the bases `fetch` gives for
    /// `[start, end)` (clipped to the contig, in any case).
    ///
    /// # Errors
    ///
    /// If `fetch` fails.
    pub fn contig(
        &self,
        name: &str,
        len: i64,
        fetch: impl FnMut(i64, i64) -> io::Result<Vec<u8>>,
    ) -> io::Result<Contig> {
        let sites = Sites::scan(len, fetch)?;
        let windows = usize::try_from((len.max(0) + self.window - 1) / self.window).expect("a window count");
        let merged = Merged {
            kinds: Kinds::default(),
            counts: vec![0; sites.count() * 2 * Kinds::default().width()],
            windows: vec![Vec::new(); windows],
            fetched: Fetched::default(),
        };
        Ok(Contig { pass: self.clone(), name: name.to_owned(), len, windows, sites, merged: Mutex::new(merged) })
    }
}

/// A contig's CpG sites (each CpG's C and G), numbered in order.
#[derive(Debug)]
struct Sites {
    len: i64,
    /// One bit per position.
    bits: Vec<u64>,
    /// Sites before each word of `bits`.
    ranks: Vec<u32>,
}

/// Bases scanned at a time.
const SCAN_BLOCK: i64 = 1 << 20;

impl Sites {
    fn scan(len: i64, mut fetch: impl FnMut(i64, i64) -> io::Result<Vec<u8>>) -> io::Result<Self> {
        let len = len.max(0);
        let mut bits = vec![0u64; (len as usize).div_ceil(64)];
        let mut set = |p: i64| {
            if p < len {
                bits[(p >> 6) as usize] |= 1 << (p & 63);
            }
        };
        for start in (0..len).step_by(SCAN_BLOCK as usize) {
            // One base more, for a CpG whose C ends the block.
            let seq = fetch(start, (start + SCAN_BLOCK + 1).min(len))?;
            for (i, pair) in seq.windows(2).enumerate().take(SCAN_BLOCK as usize) {
                if pair[0] | 0x20 == b'c' && pair[1] | 0x20 == b'g' {
                    set(start + i as i64);
                    set(start + i as i64 + 1);
                }
            }
        }
        let mut ranks = Vec::with_capacity(bits.len());
        let mut sum = 0u32;
        for w in &bits {
            ranks.push(sum);
            sum += w.count_ones();
        }
        Ok(Sites { len, bits, ranks })
    }

    fn count(&self) -> usize {
        self.ranks.last().zip(self.bits.last()).map_or(0, |(&r, w)| r as usize + w.count_ones() as usize)
    }

    /// The number of `pos`, if it is a site.
    fn slot(&self, pos: i64) -> Option<usize> {
        if pos < 0 || pos >= self.len {
            return None;
        }
        let (w, b) = ((pos >> 6) as usize, pos & 63);
        let word = self.bits[w];
        (word >> b & 1 == 1).then(|| self.ranks[w] as usize + (word & ((1 << b) - 1)).count_ones() as usize)
    }

    /// The number of sites before `pos`.
    fn rank(&self, pos: i64) -> usize {
        let pos = pos.clamp(0, self.len);
        let (w, b) = ((pos >> 6) as usize, pos & 63);
        match self.bits.get(w) {
            Some(&word) => self.ranks[w] as usize + (word & ((1 << b) - 1)).count_ones() as usize,
            None => self.count(),
        }
    }

    /// The sites in `[start, end)`, in order.
    fn positions(&self, start: i64, end: i64) -> Vec<i64> {
        let (start, end) = (start.clamp(0, self.len), end.clamp(0, self.len));
        let mut out = Vec::new();
        if start >= end {
            return out;
        }
        for w in (start >> 6) as usize..=((end - 1) >> 6) as usize {
            let mut word = self.bits[w];
            while word != 0 {
                let p = ((w as i64) << 6) + i64::from(word.trailing_zeros());
                if p >= start && p < end {
                    out.push(p);
                }
                word &= word - 1;
            }
        }
        out
    }
}

/// One contig's counts, as its pieces are merged.
#[derive(Debug)]
pub struct Contig {
    pass: Pass,
    name: String,
    len: i64,
    /// Windows of the table.
    windows: usize,
    sites: Sites,
    merged: Mutex<Merged>,
}

#[derive(Debug)]
struct Merged {
    kinds: Kinds,
    /// Counters of each site and strand, laid out as `kinds` has it.
    counts: Vec<u32>,
    /// Per window, the rows each piece's reads brought it, with the piece's number.
    windows: Vec<Vec<(usize, Kinds)>>,
    fetched: Fetched,
}

impl Contig {
    /// Counts for the reads starting at or after `start`, up to the next piece's.
    #[must_use]
    pub fn piece(&self, start: i64) -> Piece<'_> {
        Piece {
            contig: self,
            start,
            first: self.sites.rank(start),
            kinds: Kinds::default(),
            counts: Vec::new(),
            windows: BTreeMap::new(),
            fetched: Fetched::default(),
        }
    }

    /// Adds a piece's counts: `index` is its place among the contig's pieces, in order of
    /// their reads.
    ///
    /// # Panics
    ///
    /// If another thread panicked while merging.
    pub fn merge(&self, index: usize, piece: Piece<'_>) {
        let mut merged = self.merged.lock().expect("lock poisoned");
        let merged = &mut *merged;
        let mut kinds = merged.kinds.clone();
        kinds.learn_kinds(&piece.kinds);
        if kinds != merged.kinds {
            let mut counts = vec![0; self.sites.count() * 2 * kinds.width()];
            add_into(&mut counts, &kinds, &merged.counts, &merged.kinds);
            merged.counts = counts;
            merged.kinds = kinds;
        }
        let at = piece.first * 2 * merged.kinds.width();
        add_into(&mut merged.counts[at..], &merged.kinds, &piece.counts, &piece.kinds);
        for (w, kinds) in piece.windows {
            merged.windows[w].push((index, kinds));
        }
        merged.fetched += piece.fetched;
    }

    /// The counts once every piece is merged.
    ///
    /// # Panics
    ///
    /// If another thread panicked while merging.
    #[must_use]
    pub fn finish(self) -> Counted {
        let merged = self.merged.into_inner().expect("lock poisoned");
        // Each window's rows in the order its reads brought them: the pieces' in order.
        let windows = merged
            .windows
            .into_iter()
            .map(|mut pieces| {
                pieces.sort_by_key(|&(i, _)| i);
                pieces.into_iter().fold(Kinds::default(), |mut all, (_, k)| {
                    all.learn_kinds(&k);
                    all
                })
            })
            .collect();
        Counted {
            pass: self.pass,
            name: self.name,
            len: self.len,
            sites: self.sites,
            kinds: merged.kinds,
            counts: merged.counts,
            windows,
            fetched: merged.fetched,
        }
    }
}

/// The counts of a contig's reads that start in one piece of it.
#[derive(Debug)]
pub struct Piece<'c> {
    contig: &'c Contig,
    start: i64,
    /// The first site the piece's reads can reach.
    first: usize,
    kinds: Kinds,
    /// Counters of each site from `first` on and strand, laid out as `kinds` has it.
    counts: Vec<u32>,
    /// Per window the reads reach, the rows they bring.
    windows: BTreeMap<usize, Kinds>,
    fetched: Fetched,
}

impl Piece<'_> {
    /// Counts the calls of `record`, if the pass's filter lets it through and it has any. The
    /// piece's records come in file order.
    ///
    /// # Errors
    ///
    /// If the record cannot be read.
    ///
    /// # Panics
    ///
    /// If the record starts before the piece.
    pub fn add(&mut self, record: &bam::Record) -> io::Result<()> {
        match convert(record, self.contig.pass.filter)? {
            Converted::Read(read) => {
                self.fetched.kept += 1;
                self.count(&read);
            }
            Converted::Untagged => self.fetched.untagged += 1,
            Converted::Unplaceable => self.fetched.unplaceable += 1,
            Converted::Filtered => {}
        }
        Ok(())
    }

    fn count(&mut self, read: &ModRead) {
        assert!(read.start >= self.start, "a read at {} in a piece from {}", read.start, self.start);
        let Pass { calling, window, overhang, .. } = self.contig.pass;
        // The windows a window's fetch would bring the read to: those it overlaps, each with
        // its overhang, an alignment without reference bases covering its start.
        let reach = read.start + (read.end - read.start).max(1);
        let last_window = self.contig.windows.saturating_sub(1);
        let first = usize::try_from((read.start - overhang).max(0) / window).expect("a window");
        let last = usize::try_from((reach - 1) / window).expect("a window").min(last_window);
        for w in first..=last {
            let kinds = self.windows.entry(w).or_default();
            if !kinds.know(read) {
                kinds.learn(read);
            }
        }
        if !self.kinds.know(read) {
            let mut kinds = self.kinds.clone();
            kinds.learn(read);
            let mut counts = vec![0; self.counts.len() / self.kinds.width() * kinds.width()];
            add_into(&mut counts, &kinds, &self.counts, &self.kinds);
            self.kinds = kinds;
            self.counts = counts;
        }
        let sites = &self.contig.sites;
        let needed = (sites.rank(read.end) - self.first) * 2 * self.kinds.width();
        if self.counts.len() < needed {
            self.counts.resize(needed, 0);
        }
        let first = self.first;
        tally(&self.kinds, &mut self.counts, read, calling, (read.start, read.end), |pos| {
            sites.slot(pos).map(|i| i - first)
        });
    }
}

/// A contig's merged counts, ready to write.
#[derive(Debug)]
pub struct Counted {
    pass: Pass,
    name: String,
    len: i64,
    sites: Sites,
    kinds: Kinds,
    counts: Vec<u32>,
    /// Each window's rows.
    windows: Vec<Kinds>,
    fetched: Fetched,
}

impl Counted {
    /// The contig's windows.
    #[must_use]
    pub fn windows(&self) -> usize {
        self.windows.len()
    }

    /// What the pieces took and left out (each alignment once).
    #[must_use]
    pub fn fetched(&self) -> Fetched {
        self.fetched
    }

    /// Appends window `index`'s bedMethyl rows to `out`, with the contig's bases `fetch` gives
    /// for `[start, end)` (clipped to the contig, in any case).
    ///
    /// # Errors
    ///
    /// If `fetch` fails.
    pub fn write_window(
        &self,
        index: usize,
        mut fetch: impl FnMut(i64, i64) -> io::Result<Vec<u8>>,
        out: &mut String,
    ) -> io::Result<()> {
        let kinds = &self.windows[index];
        if kinds.bases.is_empty() {
            // No reads with calls reach it.
            return Ok(());
        }
        let window = self.pass.window;
        let start = index as i64 * window;
        let end = (start + window).min(self.len);
        let stop = end + self.pass.overhang;
        let positions = self.sites.positions(start, stop);
        let first = self.sites.rank(start);
        let (from, width) = (self.kinds.width(), kinds.width());
        let mut counts = vec![0; positions.len() * 2 * width];
        let span = &self.counts[first * 2 * from..(first + positions.len()) * 2 * from];
        // The window's rows are among the contig's.
        let columns = self.kinds.columns_of(kinds);
        for (dst, src) in counts.chunks_exact_mut(width).zip(span.chunks_exact(from)) {
            for (d, &c) in dst.iter_mut().zip(&columns) {
                *d = src[c];
            }
        }
        let ref_start = (start - 1).max(0);
        let mut seq = fetch(ref_start, stop + 1)?;
        seq.make_ascii_uppercase();
        let reference = RefWindow { start: ref_start, seq };
        Counts::at_sites(start, stop, kinds.clone(), positions, counts).finish().write_bedmethyl(
            &self.name,
            start,
            end,
            &self.pass.rows,
            Some(&reference),
            out,
        );
        Ok(())
    }
}

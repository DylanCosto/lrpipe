//! CpG methylation from 5mC calls by counting, as pb-CpG-tools' count mode reports it.
//!
//! A read counts at a reference CpG when it has a base aligned to the C. Its value
//! is the 5mC call it places there, or 0 when it has none (another base, a deletion, or an
//! unlisted C). A forward read places a call where the called C aligns; a reverse read where
//! the base after it in the read as sequenced (the CpG's G, read on that strand) aligns, so a
//! call on an inserted base still lands on the site. Reads with fewer than [`EDGE`] aligned
//! bases (insertions, deletions and clips aside) between that base and either end of their
//! alignment are left out. A value over one half is modified. A site reports the share of modified reads, the
//! counts, and the mean value of the modified and of the unmodified reads.

use std::fmt::Write as _;

use crate::model::{Histogram, NO_CALL};
use crate::mods::Code;
use crate::pileup::{Partition, RefWindow, probability};
use crate::read::{Convention, Event, MISSING, ModRead, Op, Strand};

/// 5mC.
const METHYL: Code = Code::Letter(b'm');

/// Reads with fewer aligned bases than this between their base at a site and either end of
/// their alignment do not count there, as in pb-CpG-tools.
pub const EDGE: i64 = 20;

/// Which CpG sites are reported.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Sites {
    /// CpG in the reference.
    #[default]
    Reference,
    /// CpG in most reads covering it.
    Denovo,
}

/// Per-site sums for one partition of the reads.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Site {
    modified: u32,
    unmodified: u32,
    /// Sums in single precision, as pb-CpG-tools keeps them (it shows in the third decimal).
    modified_p: f32,
    unmodified_p: f32,
}

/// What one read shows at one site.
#[derive(Clone, Copy, Debug, Default)]
struct Seen {
    /// A base aligned to the C.
    aligned: bool,
    /// The call the read places at the site, and its SEQ index.
    value: Option<(u16, u32)>,
    /// SEQ index of the base aligned to the C.
    c_index: Option<u32>,
}

/// The bases counted reads have at a site's C and at the position after it.
#[derive(Clone, Copy, Debug, Default)]
struct Bases {
    here: [u32; 4],
    next: [u32; 4],
}

impl Bases {
    /// Whether the consensus is a CpG: C the most common base at the site and G after it,
    /// strictly (a tie is no consensus).
    fn is_cpg(&self) -> bool {
        let most = |counts: &[u32; 4], k: usize| {
            counts[k] > 0 && counts.iter().enumerate().all(|(j, &n)| j == k || n < counts[k])
        };
        most(&self.here, 1) && most(&self.next, 2)
    }
}

/// Denovo sites of a window of `len` positions from `start`: where the reads that would count
/// there (a base aligned to the C, not near an end) agree, in some partition, on a C and on a G
/// at the next position.
fn denovo_sites(reads: &[ModRead], start: i64, len: usize, partitions: &[Partition]) -> Vec<bool> {
    let mut bases = vec![vec![Bases::default(); len]; partitions.len()];
    let mut placed: Vec<(i64, u32)> = Vec::new();
    for read in reads {
        let wanted: Vec<usize> = (0..partitions.len()).filter(|&k| partitions[k].admits(read)).collect();
        let Some((track, k)) = methyl_track(read) else { continue };
        if wanted.is_empty() {
            continue;
        }
        let reference_of = seq_to_reference(read);
        let (rank, aligned) = aligned_ranks(&reference_of);
        // Where the read places calls: a site's counting index is the placed call's.
        placed.clear();
        placed
            .extend(read.tracks[track].listed().filter_map(|(s, values)| placement(read, &reference_of, s, values[k])));
        placed.sort_unstable();
        // The previous position, if the read counts there.
        let mut previous: Option<(i64, usize)> = None;
        read.walk(start, start + len as i64 + 1, Convention::Reference, |pos, index, ti, _| {
            if ti != track {
                return;
            }
            let base = index.and_then(|x| read.seq.get(x as usize)).and_then(|&b| base_slot(b));
            if let (Some((p, i)), Some(b)) = (previous, base)
                && p == pos - 1
            {
                for &w in &wanted {
                    bases[w][i].next[b] += 1;
                }
            }
            previous = None;
            let (Some(index), Some(i)) = (index, index_in(pos, start, len)) else { return };
            let used = placed.binary_search_by_key(&pos, |&(p, _)| p).map_or(index, |at| placed[at].1);
            let r = rank[used as usize];
            if r < EDGE || aligned - 1 - r < EDGE {
                return;
            }
            previous = Some((pos, i));
            if let Some(b) = base {
                for &w in &wanted {
                    bases[w][i].here[b] += 1;
                }
            }
        });
    }
    (0..len).map(|i| bases.iter().any(|p| p[i].is_cpg())).collect()
}

/// Each SEQ base's rank among the read's aligned bases, and how many there are.
fn aligned_ranks(reference_of: &[Option<i64>]) -> (Vec<i64>, i64) {
    let mut rank = Vec::with_capacity(reference_of.len());
    let mut aligned = 0i64;
    for r in reference_of {
        rank.push(aligned);
        aligned += i64::from(r.is_some());
    }
    (rank, aligned)
}

/// A, C, G, T as 0 to 3.
fn base_slot(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Counts of the CpG sites of one window `[start, end)` (sites are placed at their C).
/// Reads are [added](CpgCounts::add) one at a time, so a window's memory does not grow with
/// its depth.
#[derive(Debug)]
pub struct CpgCounts {
    start: i64,
    end: i64,
    partitions: Vec<Partition>,
    /// Which sites are reported.
    mode: Sites,
    /// Per position in the window: whether the reference has a CpG there.
    reference_cpg: Vec<bool>,
    /// Per position (denovo): whether some partition's counted reads' consensus is a CpG.
    consensus_cpg: Vec<bool>,
    /// Per partition, per position.
    sites: Vec<Vec<Site>>,
    /// With [`Values::Kept`]: per partition, per position, the counted reads' values.
    values: Option<Vec<Vec<Histogram>>>,
    /// The sites of the read being added (scratch).
    seen: Vec<(usize, Seen)>,
}

/// Whether the counted reads' values are kept (binned), for estimating levels by model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Values {
    /// Only counts and sums.
    #[default]
    Dropped,
    /// Each site's [`Histogram`] of values.
    Kept,
}

impl CpgCounts {
    /// Counts the 5mC calls of `reads` at the CpG sites `mode` selects in `[start, end)`, for
    /// each of `partitions`. `reference` must cover `[start, end + 1)`.
    #[must_use]
    pub fn new(
        reads: &[ModRead],
        start: i64,
        end: i64,
        reference: &RefWindow,
        partitions: &[Partition],
        mode: Sites,
        keep: Values,
    ) -> Self {
        let consensus_cpg = match mode {
            Sites::Reference => Vec::new(),
            Sites::Denovo => denovo_sites(reads, start, (end - start).max(0) as usize, partitions),
        };
        let mut counts = Self::with_sites(start, end, reference, partitions, mode, keep, consensus_cpg);
        for read in reads {
            counts.add(read);
        }
        counts
    }

    /// No reads yet at the reference's CpG sites in `[start, end)`, for each of
    /// `partitions`: [`add`](Self::add) them one at a time. (Denovo sites depend on every
    /// read, so they are counted by [`new`](Self::new).) `reference` must cover
    /// `[start, end + 1)`.
    #[must_use]
    pub fn at_reference_sites(
        start: i64,
        end: i64,
        reference: &RefWindow,
        partitions: &[Partition],
        keep: Values,
    ) -> Self {
        Self::with_sites(start, end, reference, partitions, Sites::Reference, keep, Vec::new())
    }

    fn with_sites(
        start: i64,
        end: i64,
        reference: &RefWindow,
        partitions: &[Partition],
        mode: Sites,
        keep: Values,
        consensus_cpg: Vec<bool>,
    ) -> Self {
        let len = (end - start).max(0) as usize;
        let reference_cpg: Vec<bool> =
            (start..end).map(|p| reference.base(p) == Some(b'C') && reference.base(p + 1) == Some(b'G')).collect();
        let sites = vec![vec![Site::default(); len]; partitions.len()];
        let values = (keep == Values::Kept).then(|| vec![vec![Histogram::default(); len]; partitions.len()]);
        CpgCounts {
            start,
            end,
            partitions: partitions.to_vec(),
            mode,
            reference_cpg,
            consensus_cpg,
            sites,
            values,
            seen: Vec::new(),
        }
    }

    /// Adds `read`'s 5mC calls at the window's sites, to each partition that admits it.
    pub fn add(&mut self, read: &ModRead) {
        let (start, end) = (self.start, self.end);
        let len = (end - start).max(0) as usize;
        let partitions = &self.partitions;
        let is_site = match self.mode {
            Sites::Reference => &self.reference_cpg,
            Sites::Denovo => &self.consensus_cpg,
        };
        let (sites, values, seen) = (&mut self.sites, &mut self.values, &mut self.seen);
        let wanted: Vec<usize> = (0..partitions.len()).filter(|&k| partitions[k].admits(read)).collect();
        if wanted.is_empty() {
            return;
        }
        let Some((track, k)) = methyl_track(read) else { return };
        seen.clear();
        // The reverse strand's G is one past the window's last C.
        read.walk(start, end + 1, Convention::Reference, |pos, index, ti, event| {
            if ti != track {
                return;
            }
            // This position is a site's G, a site's C, or neither (in site order, which
            // keeps `seen` sorted when denovo sites are adjacent).
            for site in [pos - 1, pos] {
                let Ok(i) = usize::try_from(site - start) else { continue };
                if i >= len || !is_site[i] {
                    continue;
                }
                let entry = match seen.binary_search_by_key(&i, |&(j, _)| j) {
                    Ok(at) => &mut seen[at].1,
                    Err(at) => {
                        seen.insert(at, (i, Seen::default()));
                        &mut seen[at].1
                    }
                };
                // Counted when a base aligns to the C, whichever strand the read is on.
                if site == pos && !matches!(event, Event::Deletion) {
                    entry.aligned = true;
                    entry.c_index = index;
                }
            }
        });
        // Place the read's calls.
        let reference_of = seq_to_reference(read);
        for (s, values) in read.tracks[track].listed() {
            let v = values[k];
            let Some((pos, anchor)) = placement(read, &reference_of, s, v) else { continue };
            let Some(i) = index_in(pos, start, len).filter(|&i| is_site[i]) else { continue };
            match seen.binary_search_by_key(&i, |&(j, _)| j) {
                Ok(at) => seen[at].1.value = Some((v, anchor)),
                Err(at) => seen.insert(at, (i, Seen { aligned: true, value: Some((v, anchor)), ..Seen::default() })),
            }
        }
        let (rank, aligned) = aligned_ranks(&reference_of);
        for &(i, e) in seen.iter() {
            // The SEQ index of the base the read places at the C.
            let Some(index) = e.value.map(|(_, s)| s).or(e.c_index) else { continue };
            let Some(&r) = rank.get(index as usize) else { continue };
            // Near an end: fewer than EDGE aligned bases (insertions and deletions aside)
            // between the base and either end of the alignment.
            let near_end = r < EDGE || aligned - 1 - r < EDGE;
            if !e.aligned || near_end {
                continue;
            }
            let p = e.value.map_or(0.0, |(v, _)| probability(v));
            if let Some(values) = values.as_mut() {
                let v = e.value.map_or(NO_CALL, |(v, _)| v);
                for &w in &wanted {
                    values[w][i].add(v);
                }
            }
            for &w in &wanted {
                let s = &mut sites[w][i];
                if p > 0.5 {
                    s.modified += 1;
                    s.modified_p += p;
                } else {
                    s.unmodified += 1;
                    s.unmodified_p += p;
                }
            }
        }
    }

    /// Whether position `i` is reported: a reference CpG, or (denovo) a CpG in most reads of
    /// some partition.
    fn reported(&self, i: usize) -> bool {
        match self.mode {
            Sites::Reference => self.reference_cpg[i],
            Sites::Denovo => self.consensus_cpg[i],
        }
    }

    /// Writes partition `k`'s sites with at least `min_coverage` reads as pb-CpG-tools' model
    /// mode lays them out, levels estimated by [`Histogram::level`] (needs [`Values::Kept`]): level,
    /// label, coverage, estimated modified and unmodified reads, and the level those give.
    ///
    /// # Panics
    ///
    /// If the values were not kept.
    pub fn write_model(&self, contig: &str, k: usize, label: &str, min_coverage: u32, out: &mut String) {
        let values = &self.values.as_ref().expect("values kept for model estimates")[k];
        for (i, s) in self.sites[k].iter().enumerate() {
            let cov = s.modified + s.unmodified;
            if cov < min_coverage || cov == 0 || !self.reported(i) {
                continue;
            }
            let pos = self.start + i as i64;
            let level = values[i].level();
            let est_mod = (level * f64::from(cov)).round() as u32;
            let _ = writeln!(
                out,
                "{contig}\t{pos}\t{}\t{:.1}\t{label}\t{cov}\t{est_mod}\t{}\t{:.1}",
                pos + 1,
                half_up(100.0 * level, 1),
                cov - est_mod,
                half_up(100.0 * f64::from(est_mod) / f64::from(cov), 1),
            );
        }
    }

    /// Writes partition `k`'s sites with at least `min_coverage` reads as pb-CpG-tools' count
    /// mode does, labelled `label` (`Total`, `hap1`, `hap2`).
    pub fn write(&self, contig: &str, k: usize, label: &str, min_coverage: u32, out: &mut String) {
        for (i, s) in self.sites[k].iter().enumerate() {
            let cov = s.modified + s.unmodified;
            if cov < min_coverage || cov == 0 || !self.reported(i) {
                continue;
            }
            let pos = self.start + i as i64;
            let mean = |sum: f32, n: u32| if n == 0 { 0.0 } else { sum / n as f32 };
            // The share is rounded half up, the means half to even, as pb-CpG-tools prints them.
            let _ = writeln!(
                out,
                "{contig}\t{pos}\t{}\t{:.1}\t{label}\t{cov}\t{}\t{}\t{:.3}\t{:.3}",
                pos + 1,
                half_up(100.0 * f64::from(s.modified) / f64::from(cov), 1),
                s.modified,
                s.unmodified,
                mean(s.modified_p, s.modified),
                mean(s.unmodified_p, s.unmodified),
            );
        }
    }
}

/// The read's own 5mC track (C on its own strand) and the index of the 5mC code in it.
fn methyl_track(read: &ModRead) -> Option<(usize, usize)> {
    let own = if read.reverse { Strand::Reverse } else { Strand::Forward };
    read.tracks.iter().enumerate().find_map(|(ti, t)| {
        let k = t.codes.iter().position(|&c| c == METHYL)?;
        (t.base == b'C' && t.strand == own).then_some((ti, k))
    })
}

/// Where a read places its call at SEQ index `s` with value `v`: the reference position of
/// its C (forward) or of the base after it as sequenced (reverse), and that base's SEQ index.
fn placement(read: &ModRead, reference_of: &[Option<i64>], s: u32, v: u16) -> Option<(i64, u32)> {
    if v == MISSING {
        return None;
    }
    let anchor = if read.reverse { s.checked_sub(1)? } else { s };
    Some((reference_of.get(anchor as usize).copied().flatten()?, anchor))
}

/// `pos` as an index into a window starting at `start` of `len` positions.
fn index_in(pos: i64, start: i64, len: usize) -> Option<usize> {
    usize::try_from(pos - start).ok().filter(|&i| i < len)
}

/// The reference position each SEQ base aligns to (`None` for inserted and clipped bases).
fn seq_to_reference(read: &ModRead) -> Vec<Option<i64>> {
    let mut out = Vec::with_capacity(read.seq.len());
    let mut r = read.start;
    for &(op, len) in &read.cigar {
        match op {
            Op::Match => {
                out.extend((0..i64::from(len)).map(|k| Some(r + k)));
                r += i64::from(len);
            }
            Op::Read => out.extend(std::iter::repeat_n(None, len as usize)),
            Op::Deletion | Op::Skip => r += i64::from(len),
        }
    }
    out
}

/// `x` rounded to `digits` decimals, halves away from zero.
fn half_up(x: f64, digits: i32) -> f64 {
    let scale = 10f64.powi(digits);
    (x * scale).round() / scale
}

//! Per-position counts of modification calls, written as bedMethyl.
//!
//! Each call is classified by its most likely state, the canonical base or one of its
//! modifications, where the canonical base's `ML`-scale value is 255 less the
//! modifications'. A call whose most likely state is less likely than the canonical base's
//! pass threshold fails. Each row (position, strand, modification) counts its own calls
//! (`Nmod`), unmodified calls (`Ncanonical`), calls of another modification of the same base
//! (`Nother_mod`), deletions, failed calls, reads with another base (and no valid call of that
//! other base there) and reads with the base but no call.

use std::fmt::Write as _;

use crate::mods::Code;
use crate::read::{Convention, Event, MISSING, ModRead, Strand};

/// Pass thresholds per canonical base, as probabilities: a call passes when its most likely
/// state's probability is at least its base's threshold.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Thresholds([f32; 5]);

impl Thresholds {
    /// Every call passes.
    pub const NONE: Thresholds = Thresholds([0.0; 5]);

    /// The threshold for `base`.
    #[must_use]
    pub fn get(&self, base: u8) -> f32 {
        self.0[base_slot(base)]
    }

    /// Sets the threshold for `base` (`N` for every base).
    pub fn set(&mut self, base: u8, value: f32) {
        if base == b'N' {
            self.0 = [value; 5];
        } else {
            self.0[base_slot(base)] = value;
        }
    }
}

impl Default for Thresholds {
    fn default() -> Self {
        Thresholds::NONE
    }
}

fn base_slot(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

/// How calls are classified: pass thresholds and modkit's convention.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Calling {
    /// Pass thresholds.
    pub thresholds: Thresholds,
    /// Conventions for tallying.
    pub convention: Convention,
}

/// How a call is counted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    /// Unmodified.
    Canonical,
    /// The track's code at this index.
    Modified(usize),
    /// Below the pass threshold.
    Fail,
    /// No value for any code.
    NoCall,
}

/// The probability an `ML` value stands for, in single precision: the middle of its bin.
#[must_use]
pub fn probability(ml: u16) -> f32 {
    (f32::from(ml) + 0.5) / 256.0
}

/// The most likely state of a call and its probability, when the call has any value. The
/// canonical base's probability is what the modifications' leave, reckoned as `convention`
/// has it.
#[must_use]
pub fn most_likely(values: &[u16], convention: Convention) -> Option<(Class, f32)> {
    let mut sum = 0f32;
    let mut ml_sum = 0u32;
    let mut best: Option<(usize, f32)> = None;
    for (k, &v) in values.iter().enumerate() {
        if v == MISSING {
            continue;
        }
        let p = probability(v);
        sum += p;
        ml_sum += u32::from(v);
        if best.is_none_or(|(_, b)| p > b) {
            best = Some((k, p));
        }
    }
    let (k, p) = best?;
    let canonical = match convention {
        Convention::Read => 1.0 - sum,
        Convention::Reference => probability(255u32.saturating_sub(ml_sum) as u16),
    };
    // A tie with the canonical base goes to the modification with --modified-bases (modkit's
    // Reference convention, checked against its calls with --no-filtering or a threshold below
    // one half). modkit's default mode (Read) gives it to the canonical base here: either rule
    // gave the same calls in 3 Mb of HG002 ONT, where no such tie arose, so it is unchecked.
    let canonical_wins = match convention {
        Convention::Reference => canonical > p,
        Convention::Read => canonical >= p,
    };
    Some(if canonical_wins { (Class::Canonical, canonical) } else { (Class::Modified(k), p) })
}

/// Classifies a call against `threshold`.
#[must_use]
pub fn classify(values: &[u16], threshold: f32, convention: Convention) -> Class {
    match most_likely(values, convention) {
        None => Class::NoCall,
        Some((_, q)) if q < threshold => Class::Fail,
        Some((class, _)) => class,
    }
}

/// The rows a window can have: canonical bases and modification codes seen in its reads.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Kinds {
    /// Canonical bases, in the order the reads brought them.
    pub(crate) bases: Vec<u8>,
    /// (base index, code): in the order the reads brought them until [`Counts::finish`] sorts
    /// them by code.
    pub(crate) codes: Vec<(usize, Code)>,
}

impl Kinds {
    /// Whether every base and code of `read` is already a row.
    pub(crate) fn know(&self, read: &ModRead) -> bool {
        read.tracks.iter().all(|t| {
            self.bases
                .iter()
                .position(|&b| b == t.base)
                .is_some_and(|b| t.codes.iter().all(|&c| self.codes.contains(&(b, c))))
        })
    }

    /// Adds `read`'s bases and codes that are new, after those already here.
    pub(crate) fn learn(&mut self, read: &ModRead) {
        for t in &read.tracks {
            let b = self.base_index(t.base);
            for &c in &t.codes {
                if !self.codes.contains(&(b, c)) {
                    self.codes.push((b, c));
                }
            }
        }
    }

    /// Adds `other`'s bases and codes that are new, in its order: as if the reads that brought
    /// `other`'s were added after those that brought these.
    pub(crate) fn learn_kinds(&mut self, other: &Kinds) {
        for &base in &other.bases {
            self.base_index(base);
        }
        for &(b, c) in &other.codes {
            let b = self.base_index(other.bases[b]);
            if !self.codes.contains(&(b, c)) {
                self.codes.push((b, c));
            }
        }
    }

    /// The index of `base`, added last when new.
    fn base_index(&mut self, base: u8) -> usize {
        match self.bases.iter().position(|&b| b == base) {
            Some(b) => b,
            None => {
                self.bases.push(base);
                self.bases.len() - 1
            }
        }
    }

    /// Where each of `other`'s counters is in these kinds' layout: `other`'s must be among
    /// these.
    pub(crate) fn columns_of(&self, other: &Kinds) -> Vec<usize> {
        let base = |b: usize| self.bases.iter().position(|&x| x == other.bases[b]).expect("a base these have");
        let mut columns = Vec::with_capacity(other.width());
        for b in 0..other.bases.len() {
            columns.extend((0..BASE_FIELDS).map(|f| base(b) * BASE_FIELDS + f));
        }
        for &(b, c) in &other.codes {
            let k = self.codes.iter().position(|&x| x == (base(b), c)).expect("a code these have");
            columns.push(self.bases.len() * BASE_FIELDS + k);
        }
        columns.push(self.width() - 1);
        columns
    }

    /// Counters per position and strand: per base canonical, fail, no call and other base;
    /// per code its calls; and deletions.
    pub(crate) fn width(&self) -> usize {
        self.bases.len() * BASE_FIELDS + self.codes.len() + 1
    }
}

const CANONICAL: usize = 0;
const FAIL: usize = 1;
const NO_CALL: usize = 2;
const OTHER_BASE: usize = 3;
const BASE_FIELDS: usize = 4;

/// Which reads a table counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Partition {
    /// All of them.
    All,
    /// Those of haplotype 1 or 2.
    Haplotype(u8),
}

impl Partition {
    /// Whether the partition has `read`.
    #[must_use]
    pub fn admits(self, read: &ModRead) -> bool {
        match self {
            Partition::All => true,
            Partition::Haplotype(h) => read.haplotype == Some(h),
        }
    }
}

/// Other-base events of one read at one position, held until its valid calls there are known.
#[derive(Debug, Default)]
struct Pending {
    pos: Option<i64>,
    /// Whether the read has a valid call on each strand.
    valid: [bool; 2],
    /// Other-base counters to add to, and their strands.
    other: Vec<(usize, usize)>,
}

impl Pending {
    fn flush(&mut self, counts: &mut [u32]) {
        for &(i, s) in &self.other {
            if !self.valid[s] {
                counts[i] += 1;
            }
        }
        self.other.clear();
        self.valid = [false; 2];
    }
}

/// Counts of one window `[start, end)` of a contig: reads are [added](Counts::add) one at a
/// time, so a window's memory does not grow with its depth.
#[derive(Debug)]
pub struct Counts {
    start: i64,
    len: usize,
    kinds: Kinds,
    /// Counters of each position (or of each of `sites`) and strand, [`Kinds::width`] each.
    counts: Vec<u32>,
    /// The positions counted, when only some are, in order; every position of the window
    /// otherwise.
    sites: Option<Vec<i64>>,
}

/// Adds the counters of `from`, laid out as `from_kinds` has it, to those of `to`, laid out as
/// `to_kinds` has it (which has every base and code of `from_kinds`), row by row (a row is a
/// slot's strand) for as many rows as `from` has.
pub(crate) fn add_into(to: &mut [u32], to_kinds: &Kinds, from: &[u32], from_kinds: &Kinds) {
    let columns = to_kinds.columns_of(from_kinds);
    for (src, dst) in from.chunks_exact(from_kinds.width()).zip(to.chunks_exact_mut(to_kinds.width())) {
        for (&n, &c) in src.iter().zip(&columns) {
            dst[c] += n;
        }
    }
}

/// Adds `read`'s calls over `[from, to)`, classified as `calling` says, to `counts`, laid out
/// as `kinds` has it (which knows the read's bases and codes): at the positions `slot` gives a
/// counter slot (each slot has [`Kinds::width`] counters per strand), the others left out.
pub(crate) fn tally(
    kinds: &Kinds,
    counts: &mut [u32],
    read: &ModRead,
    calling: Calling,
    (from, to): (i64, i64),
    slot: impl Fn(i64) -> Option<usize>,
) {
    let Calling { thresholds, convention } = calling;
    let width = kinds.width();
    let mut code_index: Vec<usize> = Vec::new();
    // Track index -> (base index, code indices).
    let maps: Vec<(usize, Vec<usize>)> = read
        .tracks
        .iter()
        .map(|t| {
            let b = kinds.bases.iter().position(|&x| x == t.base).expect("bases come from these reads");
            code_index.clear();
            code_index.extend(
                t.codes
                    .iter()
                    .map(|&c| kinds.codes.iter().position(|&k| k == (b, c)).expect("codes come from these reads")),
            );
            (b, code_index.clone())
        })
        .collect();
    // A read counts as another base for a row only where it has no valid call of
    // another base on that strand: other-base events wait until the position is done.
    let mut pending = Pending::default();
    let count = |pos: i64, _, ti: usize, event: Event<'_>| {
        if pending.pos != Some(pos) {
            pending.flush(counts);
            pending.pos = Some(pos);
        }
        let t = &read.tracks[ti];
        let (b, codes) = &maps[ti];
        let s = t.strand.index();
        let row = (slot(pos).expect("only counted positions are walked") * 2 + s) * width;
        let base = row + b * BASE_FIELDS;
        match event {
            Event::Call(values) => match classify(values, thresholds.get(t.base), convention) {
                Class::Canonical => {
                    counts[base + CANONICAL] += 1;
                    pending.valid[s] = true;
                }
                Class::Modified(k) => {
                    counts[row + kinds.bases.len() * BASE_FIELDS + codes[k]] += 1;
                    pending.valid[s] = true;
                }
                Class::Fail => counts[base + FAIL] += 1,
                Class::NoCall => counts[base + NO_CALL] += 1,
            },
            Event::Unmodified => {
                counts[base + CANONICAL] += 1;
                pending.valid[s] = true;
            }
            Event::NoCall => counts[base + NO_CALL] += 1,
            Event::Uncounted => pending.valid[s] = true,
            Event::OtherBase if convention == Convention::Reference => counts[base + OTHER_BASE] += 1,
            Event::OtherBase => pending.other.push((base + OTHER_BASE, s)),
            Event::Deletion => counts[row + width - 1] += 1,
        }
    };
    read.walk_where(from, to, |pos| slot(pos).is_some(), convention, count);
    pending.flush(counts);
}

impl Counts {
    /// Counts the calls of `reads` (those `partition` admits) over `[start, end)`, classified as
    /// `calling` says; with `only`, at the positions it marks (by offset from `start`) alone.
    #[must_use]
    pub fn new(
        reads: &[ModRead],
        start: i64,
        end: i64,
        calling: Calling,
        partition: Partition,
        only: Option<&[bool]>,
    ) -> Self {
        let mut counts = Counts::empty(start, end);
        for read in reads {
            counts.add(read, calling, partition, only);
        }
        counts.finish()
    }

    /// No reads yet over `[start, end)`: [`add`](Self::add) them, then [`finish`](Self::finish).
    #[must_use]
    pub fn empty(start: i64, end: i64) -> Self {
        let len = (end - start).max(0) as usize;
        let kinds = Kinds::default();
        let counts = vec![0; len * 2 * kinds.width()];
        Counts { start, len, kinds, counts, sites: None }
    }

    /// Counts of `[start, end)` at `sites` alone (in order, inside the window), already
    /// tallied: `counts` has [`Kinds::width`] counters per site and strand. Positions that are
    /// not sites have no calls.
    pub(crate) fn at_sites(start: i64, end: i64, kinds: Kinds, sites: Vec<i64>, counts: Vec<u32>) -> Self {
        debug_assert_eq!(counts.len(), sites.len() * 2 * kinds.width());
        Counts { start, len: (end - start).max(0) as usize, kinds, counts, sites: Some(sites) }
    }

    /// The counters' slot (position or site) of `pos`, if it is counted.
    fn slot(&self, pos: i64) -> Option<usize> {
        match &self.sites {
            None => usize::try_from(pos - self.start).ok().filter(|&i| i < self.len),
            Some(sites) => sites.binary_search(&pos).ok(),
        }
    }

    /// Counter slots: positions, or sites.
    fn slots(&self) -> usize {
        self.sites.as_ref().map_or(self.len, Vec::len)
    }

    /// Adds `read`'s calls if `partition` admits it, classified as `calling` says; with
    /// `only`, at the positions it marks (by offset from the window's start) alone. Its bases
    /// and codes become rows of the window whatever its partition, as every table of a window
    /// has the same rows.
    pub fn add(&mut self, read: &ModRead, calling: Calling, partition: Partition, only: Option<&[bool]>) {
        if !self.kinds.know(read) {
            let mut kinds = self.kinds.clone();
            kinds.learn(read);
            self.relayout(kinds);
        }
        if !partition.admits(read) {
            return;
        }
        let (start, end) = (self.start, self.start + self.len as i64);
        if let Some(keep) = only {
            assert!(keep.len() >= self.len, "{} marks for {} positions", keep.len(), self.len);
        }
        let slot = |pos: i64| {
            let i = (pos - start) as usize;
            only.is_none_or(|keep| keep[i]).then_some(i)
        };
        tally(&self.kinds, &mut self.counts, read, calling, (start, end), slot);
    }

    /// Moves the counters to `kinds`' layout, which has every base and code of the current
    /// one.
    fn relayout(&mut self, kinds: Kinds) {
        let mut counts = vec![0u32; self.slots() * 2 * kinds.width()];
        add_into(&mut counts, &kinds, &self.counts, &self.kinds);
        self.kinds = kinds;
        self.counts = counts;
    }

    /// Puts the codes in output order (by code, then base) once every read is added.
    #[must_use]
    pub fn finish(mut self) -> Self {
        let mut order: Vec<usize> = (0..self.kinds.codes.len()).collect();
        order.sort_by_key(|&k| {
            let (b, c) = self.kinds.codes[k];
            (c, b)
        });
        if order.iter().enumerate().any(|(i, &k)| i != k) {
            let width = self.kinds.width();
            let first = self.kinds.bases.len() * BASE_FIELDS;
            let mut sorted = vec![0u32; order.len()];
            for row in self.counts.chunks_exact_mut(width) {
                for (to, &from) in sorted.iter_mut().zip(&order) {
                    *to = row[first + from];
                }
                row[first..first + order.len()].copy_from_slice(&sorted);
            }
            self.kinds.codes = order.iter().map(|&k| self.kinds.codes[k]).collect();
        }
        self
    }

    /// The counts of one position and strand for each row kind: per code (or, with
    /// [`Mods::Combined`], per canonical base), in output order.
    fn tallies(&self, pos: i64, strand: Strand, mods: Mods, out: &mut Vec<(Label, Tally)>) {
        out.clear();
        let width = self.kinds.width();
        let nb = self.kinds.bases.len();
        // A position that is not counted has no calls.
        let none;
        let row = match self.slot(pos) {
            Some(i) => &self.counts[(i * 2 + strand.index()) * width..][..width],
            None => {
                none = vec![0; width];
                &none[..]
            }
        };
        let delete = row[width - 1];
        let base_tally = |b: usize| {
            let base = &row[b * BASE_FIELDS..][..BASE_FIELDS];
            Tally {
                canonical: base[CANONICAL],
                delete,
                fail: base[FAIL],
                diff: base[OTHER_BASE],
                nocall: base[NO_CALL],
                ..Tally::default()
            }
        };
        match mods {
            Mods::Separate => {
                for &(b, code) in &self.kinds.codes {
                    let mut t = base_tally(b);
                    for (k, &(x, c)) in self.kinds.codes.iter().enumerate() {
                        if x != b {
                            continue;
                        }
                        let n = row[nb * BASE_FIELDS + k];
                        if c == code { t.modified += n } else { t.other += n }
                    }
                    out.push((Label::Code(self.kinds.bases[b], code), t));
                }
            }
            Mods::Combined => {
                let mut bases: Vec<usize> = (0..nb).collect();
                bases.sort_by_key(|&b| self.kinds.bases[b]);
                for b in bases {
                    let mut t = base_tally(b);
                    t.modified = (0..self.kinds.codes.len())
                        .filter(|&k| self.kinds.codes[k].0 == b)
                        .map(|k| row[nb * BASE_FIELDS + k])
                        .sum();
                    out.push((Label::Base(self.kinds.bases[b]), t));
                }
            }
        }
    }

    /// Writes the bedMethyl rows (18 columns, as modkit's `pileup`) of positions `[from, to)`
    /// to `out`: one per position, strand and modification (or base) with at least one valid
    /// call, as `rows` selects and combines them. `reference` must cover `[from - 1, to + 1)`
    /// when `rows` filters by the reference; with combined strands the counts must cover
    /// `to`.
    pub fn write_bedmethyl(
        &self,
        contig: &str,
        from: i64,
        to: i64,
        rows: &Rows,
        reference: Option<&RefWindow>,
        out: &mut String,
    ) {
        let mut tallies = Vec::new();
        let mut partner = Vec::new();
        let base_at = |pos: i64, strand: Strand| -> Option<u8> {
            let b = reference?.base(pos)?;
            Some(match strand {
                Strand::Forward => b,
                Strand::Reverse => complement(b),
            })
        };
        let reference_base = |pos: i64| reference.and_then(|r| r.base(pos));
        let cpg = |pos: i64, strand: Strand| match strand {
            Strand::Forward => reference_base(pos) == Some(b'C') && reference_base(pos + 1) == Some(b'G'),
            Strand::Reverse => reference_base(pos) == Some(b'G') && reference_base(pos - 1) == Some(b'C'),
        };
        let wanted = |label: Label, pos: i64, strand: Strand| -> bool {
            let Some(selected) = &rows.modified_bases else { return true };
            let ok = match label {
                Label::Code(b, c) => selected.iter().any(|&(sb, sc)| sb == b && sc.is_none_or(|sc| sc == c)),
                Label::Base(b) => selected.iter().any(|&(sb, _)| sb == b),
            };
            ok && base_at(pos, strand) == Some(label.base())
        };
        for pos in from..to {
            match rows.strands {
                Strands::Separate => {
                    for strand in [Strand::Forward, Strand::Reverse] {
                        if rows.motif == Motif::CpG && !cpg(pos, strand) {
                            continue;
                        }
                        self.tallies(pos, strand, rows.mods, &mut tallies);
                        for &(label, t) in &tallies {
                            if t.valid() > 0 && wanted(label, pos, strand) {
                                t.write(contig, pos, strand.symbol(), label, out);
                            }
                        }
                    }
                }
                Strands::Combined => {
                    if !cpg(pos, Strand::Forward) {
                        continue;
                    }
                    self.tallies(pos, Strand::Forward, rows.mods, &mut tallies);
                    self.tallies(pos + 1, Strand::Reverse, rows.mods, &mut partner);
                    for &(label, t) in &tallies {
                        let t = partner.iter().find(|(l, _)| *l == label).map_or(t, |&(_, p)| t + p);
                        if t.valid() > 0 && wanted(label, pos, Strand::Forward) {
                            t.write(contig, pos, '.', label, out);
                        }
                    }
                    // Kinds the reverse strand has and the forward one does not.
                    for &(label, p) in &partner {
                        if !tallies.iter().any(|(l, _)| *l == label)
                            && p.valid() > 0
                            && wanted(label, pos, Strand::Forward)
                        {
                            p.write(contig, pos, '.', label, out);
                        }
                    }
                }
            }
        }
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// A stretch of reference sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RefWindow {
    /// 0-based position of the first base.
    pub start: i64,
    /// Uppercase bases.
    pub seq: Vec<u8>,
}

impl RefWindow {
    /// The base at `pos`, if the window has it.
    #[must_use]
    pub fn base(&self, pos: i64) -> Option<u8> {
        usize::try_from(pos - self.start).ok().and_then(|i| self.seq.get(i).copied())
    }
}

/// What a row counts: one modification code of a canonical base, or all of a base's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Label {
    /// One code.
    Code(u8, Code),
    /// Every modification of a base.
    Base(u8),
}

impl Label {
    fn base(self) -> u8 {
        match self {
            Label::Code(b, _) | Label::Base(b) => b,
        }
    }
}

impl std::fmt::Display for Label {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Label::Code(_, c) => write!(f, "{c}"),
            Label::Base(b) => write!(f, "{}", char::from(*b)),
        }
    }
}

/// One row's counts.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Tally {
    modified: u32,
    canonical: u32,
    other: u32,
    delete: u32,
    fail: u32,
    diff: u32,
    nocall: u32,
}

impl std::ops::Add for Tally {
    type Output = Tally;
    fn add(self, o: Tally) -> Tally {
        Tally {
            modified: self.modified + o.modified,
            canonical: self.canonical + o.canonical,
            other: self.other + o.other,
            delete: self.delete + o.delete,
            fail: self.fail + o.fail,
            diff: self.diff + o.diff,
            nocall: self.nocall + o.nocall,
        }
    }
}

impl Tally {
    fn valid(&self) -> u32 {
        self.modified + self.other + self.canonical
    }

    fn write(&self, contig: &str, pos: i64, strand: char, label: Label, out: &mut String) {
        let valid = self.valid();
        // As modkit reckons it: in single precision, the share first (rounding differs in the
        // last digit for some shares otherwise, such as 94 of 320).
        let pct = self.modified as f32 / valid as f32 * 100.0;
        let e = pos + 1;
        let _ = writeln!(
            out,
            "{contig}\t{pos}\t{e}\t{label}\t{valid}\t{strand}\t{pos}\t{e}\t255,0,0\t{valid}\t{pct:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.modified, self.canonical, self.other, self.delete, self.fail, self.diff, self.nocall
        );
    }
}

/// Positions a table reports.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Motif {
    /// Every position.
    #[default]
    Any,
    /// CpG sites (the C on the forward strand, the G on the reverse).
    CpG,
}

/// Whether strands are reported apart.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Strands {
    /// One row per strand.
    #[default]
    Separate,
    /// A CpG's two strands summed onto its C (strand `.`).
    Combined,
}

/// Whether modifications of a base are reported apart.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mods {
    /// One row per modification.
    #[default]
    Separate,
    /// One row per canonical base, every modification counting as modified.
    Combined,
}

/// Which rows a bedMethyl table has.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Rows {
    /// Only these canonical bases (and codes, `None` for all of a base's), and only where
    /// the reference has the base; every row when `None`.
    pub modified_bases: Option<Vec<(u8, Option<Code>)>>,
    /// Positions reported.
    pub motif: Motif,
    /// Strands apart or combined (combined needs [`Motif::CpG`]).
    pub strands: Strands,
    /// Modifications apart or combined.
    pub mods: Mods,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::{Bases, decode};
    use crate::read::Op;

    /// A forward read over `ACGACGACGA` at 0, with these `MM`/`ML` tags.
    fn read(mm: &str, ml: &[u8]) -> ModRead {
        let seq = b"ACGACGACGA".to_vec();
        let groups = decode(mm.as_bytes(), ml, &seq, false, Bases::ALL).expect("valid tags");
        ModRead::new(0, false, vec![(Op::Match, 10)], seq, None, groups)
    }

    fn table(reads: &[ModRead]) -> String {
        let calling = Calling { thresholds: Thresholds::NONE, convention: Convention::Read };
        let counts = Counts::new(reads, 0, 10, calling, Partition::All, None);
        let mut text = String::new();
        counts.write_bedmethyl("chr", 0, 10, &Rows::default(), None, &mut text);
        text
    }

    #[test]
    fn rows_do_not_depend_on_which_read_brings_them_first() {
        // The first read has only 5hmC calls, the second 5mC as well: the rows grow as reads
        // are added, and the table is the same whichever comes first.
        let only_h = read("C+h?,0,1;", &[220, 30]);
        let both = read("C+m?,0,0,0;C+h?,0,0,0;", &[200, 10, 90, 20, 150, 60]);
        let forward = table(&[only_h.clone(), both.clone()]);
        assert_eq!(forward, table(&[both, only_h]));
        assert!(forward.contains("\th\t") && forward.contains("\tm\t"), "{forward}");
    }

    #[test]
    fn canonical_takes_what_the_modifications_leave() {
        assert_eq!(most_likely(&[200], Convention::Read), Some((Class::Modified(0), 200.5 / 256.0)));
        assert_eq!(most_likely(&[20], Convention::Read), Some((Class::Canonical, 1.0 - 20.5 / 256.0)));
        assert_eq!(most_likely(&[100, 120], Convention::Read).map(|(c, _)| c), Some(Class::Modified(1)));
        assert_eq!(most_likely(&[10, MISSING], Convention::Read).map(|(c, _)| c), Some(Class::Canonical));
        assert_eq!(most_likely(&[MISSING, MISSING], Convention::Read), None);
    }

    #[test]
    fn calls_below_threshold_fail() {
        let t = 209.5 / 256.0;
        assert_eq!(classify(&[200], t, Convention::Read), Class::Fail);
        assert_eq!(classify(&[209], t, Convention::Read), Class::Modified(0));
        assert_eq!(classify(&[46], t, Convention::Read), Class::Canonical);
        // Two modifications leave the canonical base half a bin less than one would, unless
        // reckoned as a bin of its own.
        assert_eq!(classify(&[20, 26], t, Convention::Read), Class::Fail);
        assert_eq!(classify(&[20, 26], t, Convention::Reference), Class::Canonical);
        assert_eq!(classify(&[MISSING], 0.0, Convention::Read), Class::NoCall);
    }
}

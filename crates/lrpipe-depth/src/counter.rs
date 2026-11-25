//! Per-base depth for one contig, counted exactly as mosdepth 0.3.14 counts it.
//!
//! Line references are to `mosdepth.nim` at v0.3.14.

use std::collections::HashMap;

use crate::alignment::{Alignment, Cigar, flag, op};

/// How a record adds depth (mosdepth `-x` and `-a`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mode {
    /// Count bases the CIGAR aligns (M, =, X), not deletions or skips; correct mate overlaps.
    #[default]
    Cigar,
    /// `-x`: count the whole span from start to end; no mate-overlap correction.
    Fast,
    /// `-a`: count the full fragment of read 1 of proper pairs.
    Fragment,
}

/// Which records count (mosdepth `-Q -l -u -F -i -R`).
#[derive(Clone, Debug)]
pub struct Filter {
    /// Records with a lower mapping quality are skipped (`-Q`, default 0).
    pub min_mapq: u8,
    /// Records with |TLEN| below this are skipped (`-l`, default -1).
    pub min_frag_len: i64,
    /// Records with |TLEN| above this are skipped (`-u`, default unbounded).
    pub max_frag_len: i64,
    /// Records with any of these flag bits are skipped (`-F`, default 1796).
    pub exclude_flags: u16,
    /// When nonzero, records need at least one of these bits (`-i`, default 0).
    pub include_flags: u16,
    /// When non-empty, only records with one of these RG values count (`-R`).
    pub read_groups: Vec<Vec<u8>>,
}

impl Default for Filter {
    fn default() -> Self {
        Filter {
            min_mapq: 0,
            min_frag_len: -1,
            max_frag_len: i64::MAX,
            exclude_flags: 1796,
            include_flags: 0,
            read_groups: Vec::new(),
        }
    }
}

impl Filter {
    /// mosdepth.nim:276-285.
    #[inline]
    pub fn passes(&self, a: &Alignment<'_>) -> bool {
        if a.mapq < self.min_mapq {
            return false;
        }
        let isize = i64::from(a.template_len).abs();
        if isize < self.min_frag_len || isize > self.max_frag_len {
            return false;
        }
        if a.flags & self.exclude_flags != 0 {
            return false;
        }
        if self.include_flags != 0 && a.flags & self.include_flags == 0 {
            return false;
        }
        if !self.read_groups.is_empty() {
            match a.read_group {
                Some(rg) if self.read_groups.iter().any(|g| g == rg) => {}
                _ => return false,
            }
        }
        true
    }
}

/// The first mate of an overlapping proper pair, kept until its mate arrives.
struct Mate {
    pos: i64,
    end: i64,
    /// Operations as (op, len); one entry means a single-operation CIGAR.
    cigar: Vec<(u8, u32)>,
}

/// A record whose mate was not found in the same chunk (see [`ContigCounter::chunk`]).
struct Missed {
    name: Vec<u8>,
    pos: i64,
    cigar: Vec<(u8, u32)>,
}

/// +1/-1 at read starts and ends over `[offset, offset + diffs.len())` of a contig, summed into
/// depth when finished (mosdepth's approach). Positions past the contig end go to the element at
/// the contig length, which only some outputs read and where they cancel out; mosdepth would
/// index out of bounds instead.
struct Window {
    len: usize,
    offset: usize,
    diffs: Vec<i32>,
}

impl Window {
    #[inline]
    fn add(&mut self, pos: i64, v: i32) {
        let p = pos.clamp(0, self.len as i64) as usize;
        if self.diffs.is_empty() {
            self.offset = p;
        } else if p < self.offset {
            // Only a chunk counter's window grows to the left (fragment mode, mate overlaps).
            self.diffs.splice(0..0, std::iter::repeat_n(0, self.offset - p));
            self.offset = p;
        }
        let i = p - self.offset;
        if i >= self.diffs.len() {
            self.diffs.resize(i + 1, 0);
        }
        self.diffs[i] += v;
    }

    fn add_window(&mut self, o: &Window) {
        if o.diffs.is_empty() {
            return;
        }
        self.add(o.offset as i64, 0);
        self.add((o.offset + o.diffs.len() - 1) as i64, 0);
        let at = o.offset - self.offset;
        for (a, b) in self.diffs[at..].iter_mut().zip(&o.diffs) {
            *a += b;
        }
    }

    /// Depth at every position of the contig, plus the element at the contig length.
    fn into_depth(self) -> Vec<i32> {
        let Window { len, offset, mut diffs } = self;
        let n = diffs.len();
        diffs.resize(len + 1, 0);
        if offset > 0 {
            diffs.copy_within(0..n, offset);
            diffs[..offset.min(n)].fill(0);
        }
        cumsum(&mut diffs);
        diffs
    }
}

/// Accumulates depth over one contig from records pushed in coordinate order.
///
/// A counter from [`new`](Self::new) covers the whole contig (4 bytes per base, allocated at the
/// first record, as mosdepth does). To count a contig in parallel pieces, give each piece of
/// reads a counter from [`chunk`](Self::chunk), which only stores the bases its reads touch, and
/// combine them in order with [`append`](Self::append): the result is identical to pushing every
/// record into one counter, including mate-overlap correction for pairs split between pieces.
pub struct ContigCounter {
    w: Window,
    mode: Mode,
    filter: Filter,
    seen: HashMap<Vec<u8>, Mate>,
    found: bool,
    pairs: Vec<(i64, i32)>,
    /// Records reads whose mate lookup failed, for [`append`](Self::append).
    chunked: bool,
    missed: Vec<Missed>,
}

impl std::fmt::Debug for ContigCounter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ContigCounter")
            .field("len", &self.w.len)
            .field("offset", &self.w.offset)
            .field("window", &self.w.diffs.len())
            .field("mode", &self.mode)
            .field("found", &self.found)
            .finish_non_exhaustive()
    }
}

impl ContigCounter {
    /// A counter for a contig of `len` bases.
    pub fn new(len: usize, mode: Mode, filter: Filter) -> Self {
        Self::with_buffer(Vec::new(), len, mode, filter)
    }

    /// Like [`new`](Self::new), reusing the allocation of a finished [`Coverage`]'s `depth`.
    pub fn with_buffer(mut buf: Vec<i32>, len: usize, mode: Mode, filter: Filter) -> Self {
        buf.clear();
        ContigCounter {
            w: Window { len, offset: 0, diffs: buf },
            mode,
            filter,
            seen: HashMap::new(),
            found: false,
            pairs: Vec::new(),
            chunked: false,
            missed: Vec::new(),
        }
    }

    /// A counter for one piece of a contig's reads, to be combined with the other pieces'
    /// counters by [`append`](Self::append).
    pub fn chunk(len: usize, mode: Mode, filter: Filter) -> Self {
        ContigCounter { chunked: true, ..Self::new(len, mode, filter) }
    }

    /// Reuses this counter's allocation for another contig.
    pub fn reset(&mut self, len: usize) {
        self.w.len = len;
        self.w.offset = 0;
        self.w.diffs.clear();
        self.seen.clear();
        self.missed.clear();
        self.found = false;
    }

    pub fn contig_len(&self) -> usize {
        self.w.len
    }

    /// Whether any record was pushed, filtered or not. mosdepth reports a contig with no
    /// records at all differently from one whose records were all filtered out.
    pub fn found(&self) -> bool {
        self.found
    }

    /// Counts one record (mosdepth.nim:272-356). Pass every record of the contig, including
    /// ones the filter will reject.
    pub fn push(&mut self, a: &Alignment<'_>) {
        if !self.found {
            self.found = true;
            if !self.chunked {
                self.w.diffs.resize(self.w.len + 1, 0);
            }
        }
        if !self.filter.passes(a) {
            return;
        }
        let end = a.end();

        // mosdepth.nim:292-340: remove double counting where mates of a proper pair overlap.
        if self.mode == Mode::Cigar && a.flags & flag::PROPER_PAIR != 0 && a.flags & flag::SUPPLEMENTARY == 0 {
            if a.mate_tid == Some(a.tid)
                && end > a.mate_pos
                && (a.pos < a.mate_pos || (a.pos == a.mate_pos && !self.seen.contains_key(a.name)))
            {
                self.seen.insert(a.name.to_vec(), Mate { pos: a.pos, end, cigar: a.cigar.iter().collect() });
            } else if let Some(mate) = self.seen.remove(a.name) {
                self.remove_overlap(a.pos, a.cigar, &mate);
            } else if self.chunked {
                self.missed.push(Missed { name: a.name.to_vec(), pos: a.pos, cigar: a.cigar.iter().collect() });
            }
        }

        match self.mode {
            Mode::Fast => {
                self.w.add(a.pos, 1);
                self.w.add(end, -1);
            }
            Mode::Fragment => {
                // Only read 1 of proper pairs (mosdepth.nim:346-353).
                if a.flags & flag::READ2 != 0 || a.flags & flag::PROPER_PAIR == 0 || a.flags & flag::SUPPLEMENTARY != 0
                {
                    return;
                }
                let start = a.pos.min(a.mate_pos);
                self.w.add(start, 1);
                self.w.add(start + i64::from(a.template_len).abs(), -1);
            }
            Mode::Cigar => {
                let w = &mut self.w;
                start_ends(a.cigar, a.pos, |p, v| w.add(p, v));
            }
        }
    }

    /// mosdepth.nim:302-340, for a read at `pos` whose earlier, overlapping mate is `mate`.
    fn remove_overlap(&mut self, pos: i64, cigar: Cigar<'_>, mate: &Mate) {
        if cigar.len() == 1 && mate.cigar.len() == 1 {
            self.w.add(pos, -1);
            self.w.add(mate.end, 1);
            return;
        }
        let mut pairs = std::mem::take(&mut self.pairs);
        pairs.clear();
        start_ends(cigar, pos, |p, v| pairs.push((p, v)));
        start_ends(Cigar::Ops(&mate.cigar), mate.pos, |p, v| pairs.push((p, v)));
        // Nim's `algorithm.sort` is a stable merge sort, as is `sort_by_key`.
        pairs.sort_by_key(|&(p, _)| p);
        let mut depth = 0;
        let mut last = 0;
        for &(p, v) in &pairs {
            if v == -1 && depth == 2 {
                self.w.add(last, -1);
                self.w.add(p, 1);
            }
            depth += v;
            last = p;
        }
        self.pairs = pairs;
    }

    /// Adds another counter's records for the same contig, ignoring mates split between the
    /// two. Exact in `Mode::Fast`, or for reads that are not proper pairs (all long reads).
    pub fn merge(&mut self, other: &ContigCounter) {
        assert_eq!(self.w.len, other.w.len, "merging counters of different contigs");
        self.found |= other.found;
        self.w.add_window(&other.w);
    }

    /// Adds the counter of the next piece of the contig's reads: `next` holds the records that
    /// come after all of this counter's (and the counters already appended), in coordinate
    /// order. Mates split between the pieces are resolved as one counter would have.
    pub fn append(&mut self, next: ContigCounter) {
        for m in &next.missed {
            if let Some(mate) = self.seen.remove(&m.name) {
                self.remove_overlap(m.pos, Cigar::Ops(&m.cigar), &mate);
            }
        }
        // Reads still unmatched have no earlier mate anywhere.
        self.missed.clear();
        self.merge(&next);
        self.seen.extend(next.seen);
    }

    /// Per-base depth: `depth[i]` for 0-based position `i`, plus a trailing 0 at the contig
    /// length (mosdepth keeps this extra element and some of its outputs read it). Empty when
    /// no record was pushed.
    pub fn finish(self) -> Coverage {
        if !self.found {
            let mut depth = self.w.diffs;
            depth.clear();
            return Coverage { depth, found: false };
        }
        Coverage { depth: self.w.into_depth(), found: true }
    }
}

fn cumsum(v: &mut [i32]) {
    let mut s = 0i32;
    for x in v {
        s = s.wrapping_add(*x);
        *x = s;
    }
}

/// mosdepth's `gen_start_ends` (mosdepth.nim:146-168): +1 where an aligned block starts and
/// -1 where it ends. Adjacent blocks (insertions between them) are merged; deletions and
/// skips split blocks. The yield order matters for overlap correction, so it is kept.
#[inline]
fn start_ends(cigar: Cigar<'_>, pos: i64, mut f: impl FnMut(i64, i32)) {
    if cigar.len() == 1 {
        let (o, len) = cigar.iter().next().unwrap();
        if o == op::MATCH {
            f(pos, 1);
            f(pos + i64::from(len), -1);
            return;
        }
    }
    let mut pos = pos;
    let mut last_stop = -1;
    for (o, len) in cigar.iter() {
        let (query, reference) = op::consumes(o);
        if !reference {
            continue;
        }
        let len = i64::from(len);
        if query {
            if pos != last_stop {
                f(pos, 1);
                if last_stop != -1 {
                    f(last_stop, -1);
                }
            }
            last_stop = pos + len;
        }
        pos += len;
    }
    if last_stop != -1 {
        f(last_stop, -1);
    }
}

/// Per-base depth of one contig.
pub struct Coverage {
    /// Depth at each 0-based position, plus a trailing element at the contig length; empty
    /// when the contig had no records.
    pub depth: Vec<i32>,
    /// Whether the contig had any records (see [`ContigCounter::found`]).
    pub found: bool,
}

impl Coverage {
    pub fn empty() -> Self {
        Coverage { depth: Vec::new(), found: false }
    }

    /// Depth over the contig's bases, without the trailing element.
    pub fn bases(&self) -> &[i32] {
        &self.depth[..self.depth.len().saturating_sub(1)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::op::*;

    fn aln<'a>(pos: i64, cigar: &'a [(u8, u32)], flags: u16, name: &'a [u8], mate_pos: i64) -> Alignment<'a> {
        Alignment {
            tid: 0,
            pos,
            flags,
            mapq: 60,
            mate_tid: Some(0),
            mate_pos,
            template_len: 0,
            name,
            cigar: Cigar::Ops(cigar),
            read_group: None,
        }
    }

    fn depth(len: usize, mode: Mode, recs: &[Alignment<'_>]) -> Vec<i32> {
        let mut c = ContigCounter::new(len, mode, Filter::default());
        for r in recs {
            c.push(r);
        }
        c.finish().bases().to_vec()
    }

    #[test]
    fn deletions_are_not_counted_insertions_do_not_split() {
        let cig = [(SOFT_CLIP, 2), (MATCH, 3), (INS, 2), (MATCH, 2), (DEL, 2), (MATCH, 1)];
        let d = depth(12, Mode::Cigar, &[aln(1, &cig, 0, b"a", -1)]);
        assert_eq!(d, [0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0]);
        let d = depth(12, Mode::Fast, &[aln(1, &cig, 0, b"a", -1)]);
        assert_eq!(d, [0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]);
    }

    #[test]
    fn overlapping_mates_count_once() {
        let m = [(MATCH, 10)];
        let pp = flag::PAIRED | flag::PROPER_PAIR;
        let recs = [aln(0, &m, pp, b"p", 5), aln(5, &m, pp | flag::READ2, b"p", 0)];
        let d = depth(20, Mode::Cigar, &recs);
        assert_eq!(d[..15], [1; 15]);
        assert_eq!(d[15], 0);
        let d = depth(20, Mode::Fast, &recs);
        assert_eq!(d[5..10], [2; 5]);
    }

    #[test]
    fn overlapping_mates_with_split_cigars() {
        let a = [(MATCH, 4), (DEL, 2), (MATCH, 4)];
        let b = [(MATCH, 3), (INS, 1), (MATCH, 5)];
        let pp = flag::PAIRED | flag::PROPER_PAIR;
        let recs = [aln(0, &a, pp, b"p", 3), aln(3, &b, pp | flag::READ2, b"p", 0)];
        let d = depth(14, Mode::Cigar, &recs);
        assert_eq!(d, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]);
    }

    #[test]
    fn records_past_the_end_do_not_panic() {
        let m = [(MATCH, 10)];
        let d = depth(5, Mode::Cigar, &[aln(3, &m, 0, b"a", -1)]);
        assert_eq!(d, [0, 0, 0, 1, 1]);
    }

    #[test]
    fn filtered_records_still_mark_the_contig_found() {
        let m = [(MATCH, 10)];
        let mut c = ContigCounter::new(20, Mode::Cigar, Filter::default());
        c.push(&aln(0, &m, flag::UNMAPPED, b"a", -1));
        let cov = c.finish();
        assert!(cov.found);
        assert!(cov.bases().iter().all(|&d| d == 0));
    }
}

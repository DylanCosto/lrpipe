//! One aligned read's modification calls, walked along its alignment to give, at each
//! reference position, what the read says there: a call, no call, a different base, or a
//! deletion.

use crate::mods::{Code, Group, Mode, Side};

/// How calls are tallied: modkit's two conventions, which its `pileup` uses without and with
/// `--modified-bases`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Convention {
    /// A call's canonical probability is one less its modifications'; the SEQ base right
    /// after a read's last listed call is not counted; a read with a valid call of another
    /// base is not counted as having another base.
    #[default]
    Read,
    /// A call's canonical probability is the middle of the `ML` bin its modifications leave;
    /// every base is counted; every read with another base is counted as such.
    Reference,
}

/// A CIGAR operation, reduced to what the walk needs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    /// Aligned bases (M, =, X).
    Match,
    /// Read bases not on the reference (I, S).
    Read,
    /// Reference bases not in the read (D).
    Deletion,
    /// Reference bases skipped (N).
    Skip,
}

/// Reference strand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    /// Forward.
    Forward,
    /// Reverse.
    Reverse,
}

impl Strand {
    /// `+` or `-`.
    #[must_use]
    pub fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }

    /// 0 for forward, 1 for reverse.
    #[must_use]
    pub fn index(self) -> usize {
        match self {
            Strand::Forward => 0,
            Strand::Reverse => 1,
        }
    }
}

/// The modification calls of one canonical base on one strand of a read: the codes of every
/// `MM` group for that base and strand, and the calls at each listed base.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Track {
    /// The canonical base, as read on `strand` of the reference (`N` for any).
    pub base: u8,
    /// The reference strand the modifications are on.
    pub strand: Strand,
    /// Modification codes.
    pub codes: Vec<Code>,
    /// Whether an unlisted base is unmodified (every group implicit) or not called.
    pub unlisted: Mode,
    /// Listed bases' indices in `SEQ`, ascending.
    calls: Vec<u32>,
    /// `codes.len()` values per listed base, in the order of `calls`; [`MISSING`] where a
    /// group did not list it.
    values: Vec<u16>,
}

/// A value no group gave.
pub const MISSING: u16 = u16::MAX;

/// What a read says about one canonical base at one reference position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Event<'t> {
    /// Its `ML` values, one per code of the track (some may be [`MISSING`]).
    Call(&'t [u16]),
    /// Unlisted by implicit groups: unmodified with certainty.
    Unmodified,
    /// A base left out of every count: under [`Convention::Read`] the first aligned SEQ base
    /// after the read's last listed call (in SEQ order); in either convention a base whose track has
    /// only unspecified modifications (`C+C`).
    Uncounted,
    /// The base is there but no group says anything about it.
    NoCall,
    /// The read has another base there.
    OtherBase,
    /// The read has a deletion there (given once per strand, with the first track of that
    /// strand).
    Deletion,
}

/// One aligned read with modification calls.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModRead {
    /// 0-based alignment start.
    pub start: i64,
    /// Exclusive alignment end.
    pub end: i64,
    /// Whether `SEQ` is the reverse complement of the read as sequenced.
    pub reverse: bool,
    /// Alignment operations and lengths, hard clips and padding left out.
    pub cigar: Vec<(Op, u32)>,
    /// `SEQ`, uppercase.
    pub seq: Vec<u8>,
    /// Haplotype (the `HP` tag), when 1 or 2.
    pub haplotype: Option<u8>,
    /// SAM flags.
    pub flags: u16,
    /// Mapping quality (255 when unavailable).
    pub mapq: u8,
    /// Calls grouped by canonical base and strand.
    pub tracks: Vec<Track>,
    /// The first aligned SEQ index past the last listed call of any track.
    after_last_listed: Option<u32>,
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

impl Track {
    /// The listed bases: SEQ index and one value per code.
    pub fn listed(&self) -> impl Iterator<Item = (u32, &[u16])> + '_ {
        let n = self.codes.len();
        self.calls.iter().enumerate().map(move |(j, &s)| (s, &self.values[j * n..(j + 1) * n]))
    }
}

impl ModRead {
    /// A read from its alignment and decoded `MM` groups.
    #[must_use]
    pub fn new(
        start: i64,
        reverse: bool,
        cigar: Vec<(Op, u32)>,
        seq: Vec<u8>,
        haplotype: Option<u8>,
        groups: Vec<Group>,
    ) -> Self {
        let span: i64 = cigar
            .iter()
            .filter(|(op, _)| matches!(op, Op::Match | Op::Deletion | Op::Skip))
            .map(|&(_, l)| i64::from(l))
            .sum();
        let tracks = tracks(groups, reverse, seq.len());
        let after_last_listed =
            tracks.iter().filter_map(|t| t.calls.last().map(|&s| s + 1)).max().and_then(|s| next_aligned(&cigar, s));
        ModRead {
            start,
            end: start + span,
            reverse,
            cigar,
            seq,
            haplotype,
            flags: 0,
            mapq: 255,
            tracks,
            after_last_listed,
        }
    }

    /// Walks the alignment over reference positions `[from, to)`, calling `f` with each
    /// position, the read's SEQ index there (`None` in a deletion), the track (its index in
    /// [`ModRead::tracks`]) and the event, in position order (tracks in order at each
    /// position), as `convention` has them.
    pub fn walk(&self, from: i64, to: i64, convention: Convention, f: impl FnMut(i64, Option<u32>, usize, Event<'_>)) {
        self.walk_where(from, to, |_| true, convention, f);
    }

    /// [`ModRead::walk`] at only the positions `keep` marks (`keep[0]` is `from`).
    ///
    /// # Panics
    ///
    /// If `keep` is shorter than `[from, to)`.
    pub fn walk_at(
        &self,
        from: i64,
        to: i64,
        keep: &[bool],
        convention: Convention,
        f: impl FnMut(i64, Option<u32>, usize, Event<'_>),
    ) {
        assert!(keep.len() as i64 >= to - from, "{} marks for {} positions", keep.len(), to - from);
        self.walk_where(from, to, |pos| keep[(pos - from) as usize], convention, f);
    }

    /// [`ModRead::walk`] at only the positions `kept` accepts.
    pub(crate) fn walk_where(
        &self,
        from: i64,
        to: i64,
        kept: impl Fn(i64) -> bool,
        convention: Convention,
        mut f: impl FnMut(i64, Option<u32>, usize, Event<'_>),
    ) {
        let mut cursors = vec![0usize; self.tracks.len()];
        let (mut r, mut q) = (self.start, 0u32);
        for &(op, len) in &self.cigar {
            if r >= to {
                break;
            }
            match op {
                Op::Match => {
                    let skip = from.saturating_sub(r).clamp(0, i64::from(len)) as u32;
                    for k in skip..len {
                        let pos = r + i64::from(k);
                        if pos >= to {
                            break;
                        }
                        if !kept(pos) {
                            continue;
                        }
                        let s = q + k;
                        let b = self.seq[s as usize];
                        for (ti, (t, cursor)) in self.tracks.iter().zip(cursors.iter_mut()).enumerate() {
                            let base = if t.strand == Strand::Forward { b } else { complement(b) };
                            while t.calls.get(*cursor).is_some_and(|&i| i < s) {
                                *cursor += 1;
                            }
                            let event = if convention == Convention::Read && Some(s) == self.after_last_listed {
                                Event::Uncounted
                            } else if t.base != b'N' && base != t.base {
                                Event::OtherBase
                            } else if t.codes.is_empty() {
                                // Only unspecified modifications: the base, and no call.
                                Event::Uncounted
                            } else {
                                match t.calls.get(*cursor) {
                                    Some(&i) if i == s => {
                                        let n = t.codes.len();
                                        Event::Call(&t.values[*cursor * n..(*cursor + 1) * n])
                                    }
                                    _ if t.unlisted == Mode::Implicit => Event::Unmodified,
                                    _ => Event::NoCall,
                                }
                            };
                            f(pos, Some(s), ti, event);
                        }
                    }
                    r += i64::from(len);
                    q += len;
                }
                Op::Read => q += len,
                Op::Deletion => {
                    let lo = r.max(from);
                    let hi = (r + i64::from(len)).min(to);
                    // Once per strand: a read with tracks for two bases has one deletion.
                    let first_of_strand =
                        |ti: usize| !self.tracks[..ti].iter().any(|t| t.strand == self.tracks[ti].strand);
                    for pos in (lo..hi).filter(|&p| kept(p)) {
                        for ti in (0..self.tracks.len()).filter(|&ti| first_of_strand(ti)) {
                            f(pos, None, ti, Event::Deletion);
                        }
                    }
                    r += i64::from(len);
                }
                Op::Skip => r += i64::from(len),
            }
        }
    }
}

/// The first SEQ index at or after `s` that aligns to the reference (past an insertion).
fn next_aligned(cigar: &[(Op, u32)], s: u32) -> Option<u32> {
    let mut q = 0u32;
    for &(op, len) in cigar {
        match op {
            Op::Match if q + len > s => return Some(s.max(q)),
            Op::Match | Op::Read => q += len,
            Op::Deletion | Op::Skip => {}
        }
    }
    None
}

/// Merges groups into one track per canonical base and strand.
/// The `j`-th call of group `g` in SEQ order: its SEQ index and its values.
fn nth_call(g: &Group, j: usize, reverse: bool, len: usize) -> (u32, &[u8]) {
    let n = g.positions.len();
    let i = if reverse { n - 1 - j } else { j };
    let s = if reverse { (len - 1 - g.positions[i] as usize) as u32 } else { g.positions[i] };
    let w = g.codes.len();
    (s, &g.probs[i * w..(i + 1) * w])
}

/// Merges groups into one track per canonical base and strand. Each group's positions are in
/// order (ascending as sequenced, so descending in SEQ for a reverse read), so a track's calls
/// come from merging its groups, without sorting.
fn tracks(groups: Vec<Group>, reverse: bool, len: usize) -> Vec<Track> {
    // Each track, and the indices of its groups with each one's first code's index.
    let mut out: Vec<(Track, Vec<(usize, usize)>)> = Vec::new();
    // The read's own strand is the reference's reverse when SEQ was reversed.
    let read_strand = if reverse { Strand::Reverse } else { Strand::Forward };
    for (gi, g) in groups.iter().enumerate() {
        let strand = match (g.side, read_strand) {
            (Side::Same, s) => s,
            (Side::Opposite, Strand::Forward) => Strand::Reverse,
            (Side::Opposite, Strand::Reverse) => Strand::Forward,
        };
        // The base as read on that strand of the reference.
        let base = match g.side {
            Side::Same => g.base,
            Side::Opposite => complement(g.base),
        };
        let k = match out.iter().position(|(t, _)| t.base == base && t.strand == strand) {
            Some(k) => k,
            None => {
                let track = Track {
                    base,
                    strand,
                    codes: Vec::new(),
                    unlisted: Mode::Implicit,
                    calls: Vec::new(),
                    values: Vec::new(),
                };
                out.push((track, Vec::new()));
                out.len() - 1
            }
        };
        let (track, parts) = &mut out[k];
        parts.push((gi, track.codes.len()));
        track.codes.extend(&g.codes);
        if g.mode == Mode::Explicit {
            track.unlisted = Mode::Explicit;
        }
    }
    let nth = |g: &Group, j: usize| nth_call(g, j, reverse, len).0;
    out.into_iter()
        .map(|(mut t, parts)| {
            let n = t.codes.len();
            let capacity = parts.iter().map(|&(gi, _)| groups[gi].positions.len()).max().unwrap_or(0);
            t.calls.reserve(capacity);
            t.values.reserve(capacity * n);
            // Groups listing the same bases (ONT's `C+h` and `C+m`) need no merging.
            let first = &groups[parts[0].0];
            if parts.iter().all(|&(gi, _)| groups[gi].positions == first.positions) {
                for j in 0..first.positions.len() {
                    let s = nth(first, j);
                    t.calls.push(s);
                    for &(gi, _) in &parts {
                        t.values.extend(nth_call(&groups[gi], j, reverse, len).1.iter().map(|&v| u16::from(v)));
                    }
                }
                return t;
            }
            let mut cursor = vec![0usize; parts.len()];
            loop {
                let next = parts
                    .iter()
                    .zip(&cursor)
                    .filter(|&(&(gi, _), &c)| c < groups[gi].positions.len())
                    .map(|(&(gi, _), &c)| nth(&groups[gi], c))
                    .min();
                let Some(s) = next else { break };
                let at = t.values.len();
                t.calls.push(s);
                t.values.resize(at + n, MISSING);
                for (&(gi, offset), c) in parts.iter().zip(cursor.iter_mut()) {
                    let g = &groups[gi];
                    if *c < g.positions.len() {
                        let (gs, probs) = nth_call(g, *c, reverse, len);
                        if gs == s {
                            for (k, &v) in probs.iter().enumerate() {
                                t.values[at + offset + k] = u16::from(v);
                            }
                            *c += 1;
                        }
                    }
                }
            }
            t
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::decode;

    fn read(seq: &[u8], reverse: bool, cigar: Vec<(Op, u32)>, mm: &[u8], ml: &[u8]) -> ModRead {
        let groups = decode(mm, ml, seq, reverse, crate::mods::Bases::ALL).expect("valid tags");
        ModRead::new(100, reverse, cigar, seq.to_vec(), None, groups)
    }

    fn events(r: &ModRead) -> Vec<(i64, char, String)> {
        events_by(r, Convention::Read)
    }

    fn events_by(r: &ModRead, convention: Convention) -> Vec<(i64, char, String)> {
        let mut out = Vec::new();
        r.walk(i64::MIN, i64::MAX, convention, |pos, _, ti, e| {
            let t = &r.tracks[ti];
            let e = match e {
                Event::Call(v) => format!("{v:?}"),
                other => format!("{other:?}"),
            };
            out.push((pos, t.strand.symbol(), e));
        });
        out
    }

    #[test]
    fn forward_read_calls_land_on_the_forward_strand() {
        // CGACG with an insertion of T after the first base: C at reference 100 and 103.
        let r = read(b"CTGACG", false, vec![(Op::Match, 1), (Op::Read, 1), (Op::Match, 4)], b"C+m?,1;", &[200]);
        let e = events(&r);
        assert_eq!(e[0], (100, '+', "NoCall".into()), "the first C is not listed (explicit mode)");
        assert_eq!(e[1], (101, '+', "OtherBase".into()));
        assert_eq!(e[3], (103, '+', "[200]".into()));
    }

    #[test]
    fn reverse_read_calls_land_on_the_reverse_strand() {
        // SEQ CGG reversed and complemented is CCG: as sequenced, C's at 0 and 1, which are
        // SEQ indices 2 and 1 (G's), reference 102 and 101.
        let r = read(b"CGG", true, vec![(Op::Match, 3)], b"C+m.,1;", &[77]);
        let e = events_by(&r, Convention::Reference);
        assert_eq!(e, [(100, '-', "OtherBase".into()), (101, '-', "[77]".into()), (102, '-', "Unmodified".into())]);
        // modkit's read convention leaves out the SEQ base after the last listed call.
        let e = events(&r);
        assert_eq!(e[2], (102, '-', "Uncounted".into()));
    }

    #[test]
    fn the_uncounted_base_is_the_next_aligned_one() {
        // CTC with the T inserted: the base after the listed C is the insertion, so the next
        // aligned base, the second C (reference 101), is the one left out.
        let r = read(b"CTC", false, vec![(Op::Match, 1), (Op::Read, 1), (Op::Match, 1)], b"C+m.,0;", &[9]);
        assert_eq!(events(&r), [(100, '+', "[9]".into()), (101, '+', "Uncounted".into())]);
        assert_eq!(events_by(&r, Convention::Reference)[1], (101, '+', "Unmodified".into()));
    }

    #[test]
    fn unspecified_modifications_count_nowhere() {
        // `C+C` lists both C's but says nothing of what they are.
        let r = read(b"CAC", false, vec![(Op::Match, 3)], b"C+C.,0,0;", &[1, 2]);
        let e = events_by(&r, Convention::Reference);
        assert_eq!(e, [(100, '+', "Uncounted".into()), (101, '+', "OtherBase".into()), (102, '+', "Uncounted".into())]);
    }

    #[test]
    fn deletions_and_windows() {
        let r = read(b"CC", false, vec![(Op::Match, 1), (Op::Deletion, 2), (Op::Match, 1)], b"C+m.,0,0;", &[1, 2]);
        let e = events(&r);
        assert_eq!(e[1], (101, '+', "Deletion".into()));
        assert_eq!(e[3], (103, '+', "[2]".into()));
        let mut inside = Vec::new();
        r.walk(101, 103, Convention::Read, |pos, _, _, _| inside.push(pos));
        assert_eq!(inside, [101, 102]);
        assert_eq!(r.end, 104);
    }

    #[test]
    fn groups_of_one_base_share_a_track() {
        let r = read(b"CAC", false, vec![(Op::Match, 3)], b"C+h?,0,0;C+m?,1;A+a.,0;", &[1, 2, 3, 4]);
        assert_eq!(r.tracks.len(), 2);
        let c = &r.tracks[0];
        assert_eq!(c.codes, [Code::Letter(b'h'), Code::Letter(b'm')]);
        let e = events(&r);
        let at = |p: i64, s: &str| e.iter().filter(|(pos, _, _)| *pos == p).any(|(_, _, x)| x == s);
        assert!(at(100, &format!("{:?}", [1, MISSING])), "{e:?}");
        assert!(at(102, "[2, 3]"), "{e:?}");
        assert!(at(101, "[4]"), "{e:?}");
    }
}

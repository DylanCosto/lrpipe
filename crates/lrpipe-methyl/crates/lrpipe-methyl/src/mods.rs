//! A read's base modification calls, decoded from its `MM` and `ML` tags as the SAMtags
//! specification defines them.
//!
//! `MM` lists groups such as `C+m?,5,12,0;` and `C+hm.,1,0;`: a canonical base, the strand
//! the modification is on relative to the read (`+` the read's own, `-` the opposite one), one
//! or more modification codes (a letter or a ChEBI number), and an optional mode. Each number
//! skips that many occurrences of the canonical base in the read as sequenced (the reverse
//! complement of `SEQ` when the read is mapped to the reverse strand) to reach the next
//! modified one. `ML` holds one byte per listed base and code, in order: the probability the
//! base carries that modification, in 256 bins.
//!
//! The mode says what an unlisted occurrence of the base means: with `.` (the default) it is
//! unmodified; with `?` nothing is known about it.

use std::fmt;

/// A modification code: a single letter such as `m` (5mC), `h` (5hmC) or `a` (6mA), or a
/// ChEBI identifier.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Code {
    /// A letter code.
    Letter(u8),
    /// A ChEBI number.
    Chebi(u32),
}

impl fmt::Display for Code {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Code::Letter(c) => write!(f, "{}", char::from(*c)),
            Code::Chebi(n) => write!(f, "{n}"),
        }
    }
}

/// Which strand a group's modifications are on, relative to the read as sequenced.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Side {
    /// The read's own strand (`+`).
    Same,
    /// The opposite strand (`-`, duplex reads).
    Opposite,
}

/// What an occurrence of the canonical base that a group does not list means.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Mode {
    /// Unmodified (`.`, and the default).
    Implicit,
    /// Unknown (`?`).
    Explicit,
}

/// Why `MM`/`ML` tags could not be decoded.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum DecodeError {
    /// A group is not `<base><strand><codes>[mode]`.
    #[error("malformed MM group {0:?}")]
    Group(String),
    /// A position skips past the end of the read.
    #[error("MM group {0:?} runs past the read's {1} bases")]
    PastEnd(String, usize),
    /// `ML` has a different number of values than `MM` lists.
    #[error("ML has {ml} values; MM lists {mm}")]
    Count {
        /// Values in `ML`.
        ml: usize,
        /// Probabilities `MM` calls for.
        mm: usize,
    },
}

/// One `MM` group, its positions resolved.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Group {
    /// The canonical base, as sequenced (`N` matches any base).
    pub base: u8,
    /// The strand of the modifications.
    pub side: Side,
    /// The modification codes, in `ML` order.
    pub codes: Vec<Code>,
    /// What unlisted occurrences of the base mean.
    pub mode: Mode,
    /// Indices of the listed bases in the read as sequenced (0 is its first base).
    pub positions: Vec<u32>,
    /// `ML` values: `codes.len()` per position, in order.
    pub probs: Vec<u8>,
}

/// The probability an `ML` value stands for: the middle of its bin, `(v + 0.5) / 256`.
#[must_use]
pub fn probability(ml: u8) -> f64 {
    (f64::from(ml) + 0.5) / 256.0
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' | b'U' | b'u' => b'A',
        _ => b'N',
    }
}

/// Canonical bases whose groups are decoded (others are skipped, their `ML` values passed
/// over).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Bases(u8);

impl Bases {
    /// Every base.
    pub const ALL: Bases = Bases(0b11111);

    /// Only `bases` (of A, C, G, T/U, N).
    #[must_use]
    pub fn only(bases: &[u8]) -> Bases {
        Bases(bases.iter().fold(0, |m, &b| m | Self::bit(b)))
    }

    fn bit(b: u8) -> u8 {
        match b.to_ascii_uppercase() {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' | b'U' => 8,
            _ => 16,
        }
    }

    /// Whether `base` is decoded.
    #[must_use]
    pub fn contains(self, base: u8) -> bool {
        self.0 & Self::bit(base) != 0
    }
}

impl Default for Bases {
    fn default() -> Self {
        Bases::ALL
    }
}

/// Decodes the groups of `mm` with the probabilities of `ml`, for a read whose `SEQ` is `seq`
/// (reverse-complemented from the read as sequenced when `reverse`). Positions index the read
/// as sequenced. `MM` and `ML` given in their legacy `Mm`/`Ml` names decode the same way.
/// Groups whose modified base, on the strand the group names, is outside `bases` are left out:
/// `G-m` is cytosine methylation on the opposite strand, so it is kept for `C`.
///
/// # Errors
///
/// If a group is malformed, a position falls past the read's end, or `ML` does not hold one
/// value per listed base and code.
pub fn decode(mm: &[u8], ml: &[u8], seq: &[u8], reverse: bool, bases: Bases) -> Result<Vec<Group>, DecodeError> {
    let n = seq.len();
    // Where the read has each base, as sequenced, built once for the groups that need it.
    let mut found: [Option<Vec<u32>>; 4] = Default::default();
    let mut groups = Vec::new();
    let mut used = 0;
    for text in mm.split(|&b| b == b';').filter(|t| !t.is_empty()) {
        let bad = || DecodeError::Group(String::from_utf8_lossy(text).into_owned());
        let mut parts = text.split(|&b| b == b',');
        let head = parts.next().ok_or_else(bad)?;
        let (&base, rest) = head.split_first().ok_or_else(bad)?;
        let base = base.to_ascii_uppercase();
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
            return Err(bad());
        }
        let (&strand, rest) = rest.split_first().ok_or_else(bad)?;
        let side = match strand {
            b'+' => Side::Same,
            b'-' => Side::Opposite,
            _ => return Err(bad()),
        };
        let (rest, mode) = match rest.last() {
            Some(b'.') => (&rest[..rest.len() - 1], Mode::Implicit),
            Some(b'?') => (&rest[..rest.len() - 1], Mode::Explicit),
            _ => (rest, Mode::Implicit),
        };
        let codes = parse_codes(rest).ok_or_else(bad)?;
        let target = if base == b'U' { b'T' } else { base };
        // The modified base is the group's base on its own strand; the deltas still count the
        // read's `target`s. modkit compares `target` itself, so with `--modified-bases 5mC` it
        // drops `G-m` groups (LRP-006).
        let modified = match side {
            Side::Same => target,
            Side::Opposite => complement(target),
        };
        if !bases.contains(modified) {
            // Pass over its ML values without resolving its positions.
            used += parts.filter(|d| !d.is_empty()).count() * codes.len();
            continue;
        }
        // A code that is the base itself (`C+C`) means a modification of unknown kind: the
        // group says where the read has the base and which bases it lists, nothing about their
        // calls (as in modkit).
        let unspecified = codes.contains(&Code::Letter(base));
        let at: Option<&[u32]> = match target {
            b'N' => None,
            _ => {
                let slot = &mut found[usize::from(Bases::bit(target).trailing_zeros() as u8)];
                Some(slot.get_or_insert_with(|| occurrences(seq, reverse, target)))
            }
        };
        let mut positions = Vec::new();
        // The occurrence of the base each delta counts from.
        let mut k = 0usize;
        for delta in parts {
            // Skip `delta` occurrences of the base, then take the next one.
            k += parse_delta(delta).ok_or_else(bad)?;
            let i = match at {
                None => (k < n).then_some(k as u32),
                Some(at) => at.get(k).copied(),
            };
            let Some(i) = i else {
                return Err(DecodeError::PastEnd(String::from_utf8_lossy(text).into_owned(), n));
            };
            positions.push(i);
            k += 1;
        }
        let want = positions.len() * codes.len();
        if unspecified {
            used += want;
            groups.push(Group { base: target, side, codes: Vec::new(), mode, positions, probs: Vec::new() });
            continue;
        }
        let probs = ml.get(used..used + want).map(<[u8]>::to_vec).unwrap_or_default();
        used += want;
        groups.push(Group { base: target, side, codes, mode, positions, probs });
    }
    if used != ml.len() {
        return Err(DecodeError::Count { ml: ml.len(), mm: used });
    }
    Ok(groups)
}

/// The indices, as sequenced and in that order, where the read has `base` (A, C, G or T).
fn occurrences(seq: &[u8], reverse: bool, base: u8) -> Vec<u32> {
    // Written without branches on the base: every index is stored, and kept when it matches.
    let mut out = vec![0u32; seq.len() + 1];
    let mut k = 0;
    if reverse {
        // As sequenced, index i is the complement of SEQ index n - 1 - i.
        let want = complement(base);
        for (i, &b) in seq.iter().rev().enumerate() {
            out[k] = i as u32;
            k += usize::from(b.to_ascii_uppercase() == want);
        }
    } else {
        for (i, &b) in seq.iter().enumerate() {
            out[k] = i as u32;
            k += usize::from(b.to_ascii_uppercase() == base);
        }
    }
    out.truncate(k);
    out
}

/// A delta of `MM`: decimal digits.
fn parse_delta(text: &[u8]) -> Option<usize> {
    if text.is_empty() || text.len() > 9 {
        return None;
    }
    text.iter().try_fold(0usize, |v, &d| d.is_ascii_digit().then(|| v * 10 + usize::from(d - b'0')))
}

/// Letters (`m`, `hm`) or one ChEBI number.
fn parse_codes(text: &[u8]) -> Option<Vec<Code>> {
    if text.is_empty() {
        return None;
    }
    if text.iter().all(u8::is_ascii_digit) {
        let n = std::str::from_utf8(text).ok()?.parse().ok()?;
        return Some(vec![Code::Chebi(n)]);
    }
    text.iter().all(u8::is_ascii_alphabetic).then(|| text.iter().map(|&c| Code::Letter(c)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positions_skip_occurrences_of_the_base() {
        // C's at 1, 3, 5, 9: skip 1 (take 3), then 0 (take 5), then 0 (take 9).
        let seq = b"ACACACGTAC";
        let g = decode(b"C+m,1,0,0;", &[10, 20, 30], seq, false, Bases::ALL).expect("valid");
        assert_eq!(g[0].positions, [3, 5, 9]);
        assert_eq!(g[0].probs, [10, 20, 30]);
        assert_eq!((g[0].base, g[0].side, g[0].mode), (b'C', Side::Same, Mode::Implicit));
    }

    #[test]
    fn reverse_reads_count_on_the_read_as_sequenced() {
        // SEQ GTGTG is CACAC reversed and complemented: as sequenced, C's at 0, 2, 4.
        let g = decode(b"C+m?,1,0;", &[200, 100], b"GTGTG", true, Bases::ALL).expect("valid");
        assert_eq!(g[0].positions, [2, 4]);
        assert_eq!(g[0].mode, Mode::Explicit);
    }

    #[test]
    fn several_codes_interleave_their_probabilities() {
        let g = decode(b"C+hm?,0,1;", &[1, 2, 3, 4], b"CCC", false, Bases::ALL).expect("valid");
        assert_eq!(g[0].codes, [Code::Letter(b'h'), Code::Letter(b'm')]);
        assert_eq!(g[0].positions, [0, 2]);
        assert_eq!(g[0].probs, [1, 2, 3, 4]);
    }

    #[test]
    fn separate_groups_and_chebi_codes() {
        let g = decode(b"C+h?,0;C+m?,0;A+a.,1;N+76792,0;", &[5, 6, 7, 8], b"CAAT", false, Bases::ALL).expect("valid");
        assert_eq!(g.len(), 4);
        assert_eq!(g[2].positions, [2]);
        assert_eq!(g[3].codes, [Code::Chebi(76792)]);
        assert_eq!(g[3].positions, [0]);
    }

    #[test]
    fn empty_groups_and_bad_tags() {
        let g = decode(b"C+m?;", &[], b"CC", false, Bases::ALL).expect("an empty group is valid");
        assert!(g[0].positions.is_empty());
        assert!(matches!(decode(b"C+m,5;", &[1], b"CC", false, Bases::ALL), Err(DecodeError::PastEnd(..))));
        assert!(matches!(
            decode(b"C+m,0;", &[1, 2], b"CC", false, Bases::ALL),
            Err(DecodeError::Count { ml: 2, mm: 1 })
        ));
        assert!(matches!(decode(b"Cxm,0;", &[1], b"CC", false, Bases::ALL), Err(DecodeError::Group(_))));
        assert!(matches!(decode(b"C+,0;", &[1], b"CC", false, Bases::ALL), Err(DecodeError::Group(_))));
    }

    #[test]
    fn opposite_strand_groups_are_kept_for_their_modified_base() {
        // LRP-006: `G-m` is 5mC on the opposite strand, so choosing C keeps it, decoded on the
        // read's G's, and choosing G does not.
        let tag = b"C+m?,0;G-m?,1;";
        let g = decode(tag, &[1, 2], b"CGAG", false, Bases::only(b"C")).expect("valid");
        assert_eq!(g.len(), 2);
        assert_eq!(
            (g[1].base, g[1].side, g[1].positions.as_slice(), g[1].probs.as_slice()),
            (b'G', Side::Opposite, &[3][..], &[2][..])
        );
        let g = decode(tag, &[1, 2], b"CGAG", false, Bases::only(b"G")).expect("valid");
        assert!(g.is_empty(), "{g:?}");
    }

    #[test]
    fn unwanted_groups_are_passed_over() {
        let g = decode(b"A+a.,0,0;C+m?,0;", &[1, 2, 3], b"AAC", false, Bases::only(b"C")).expect("valid");
        assert_eq!(g.len(), 1);
        assert_eq!((g[0].base, g[0].probs.as_slice()), (b'C', &[3][..]));
    }

    #[test]
    fn unspecified_modifications_keep_only_their_base() {
        let g = decode(b"C+C.,0;C+m?,1;", &[9, 3], b"CC", false, Bases::ALL).expect("valid");
        assert_eq!(g.len(), 2);
        assert!(g[0].codes.is_empty() && g[0].probs.is_empty());
        assert_eq!(g[0].positions, [0]);
        assert_eq!((g[1].codes.as_slice(), g[1].positions.as_slice()), (&[Code::Letter(b'm')][..], &[1][..]));
    }

    #[test]
    fn probabilities_are_bin_midpoints() {
        assert_eq!(probability(0), 0.5 / 256.0);
        assert_eq!(probability(209), 209.5 / 256.0);
    }
}

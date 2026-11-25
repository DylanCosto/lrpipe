//! The fields of an alignment record that depth counting reads.
//!
//! [`Alignment`] borrows plain values rather than a pasta record type, so callers that
//! already decode records their own way (lrpipe-sv, for one) can fill it without converting.

use std::fmt;

/// BAM CIGAR operation codes.
pub mod op {
    pub const MATCH: u8 = 0;
    pub const INS: u8 = 1;
    pub const DEL: u8 = 2;
    pub const REF_SKIP: u8 = 3;
    pub const SOFT_CLIP: u8 = 4;
    pub const HARD_CLIP: u8 = 5;
    pub const PAD: u8 = 6;
    pub const EQUAL: u8 = 7;
    pub const DIFF: u8 = 8;
    pub const BACK: u8 = 9;

    /// htslib's `bam_cigar_type`: bit 0 = consumes query, bit 1 = consumes reference.
    #[inline]
    pub fn consumes(op: u8) -> (bool, bool) {
        // "\0\1\3\2\1\0\0\3\3\0" in htslib, indexed by operation (0x3C1A7 >> (op << 1)) & 3.
        let t = (0x3C1A7u32 >> (u32::from(op) << 1)) & 3;
        (t & 1 != 0, t & 2 != 0)
    }
}

/// A CIGAR: packed little-endian BAM words (`len << 4 | op`), or already-split operations.
#[derive(Clone, Copy)]
pub enum Cigar<'a> {
    /// The raw CIGAR bytes of a BAM record (4 bytes per operation).
    Packed(&'a [u8]),
    /// Operations as (op, len) pairs; `op` is a BAM code from [`op`].
    Ops(&'a [(u8, u32)]),
}

impl<'a> Cigar<'a> {
    pub fn len(&self) -> usize {
        match self {
            Cigar::Packed(b) => b.len() / 4,
            Cigar::Ops(ops) => ops.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Iterates (op, len) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (u8, u32)> + 'a {
        let (packed, ops): (&'a [u8], &'a [(u8, u32)]) = match *self {
            Cigar::Packed(b) => (b, &[]),
            Cigar::Ops(o) => (&[], o),
        };
        packed
            .as_chunks::<4>()
            .0
            .iter()
            .map(|&c| {
                let v = u32::from_le_bytes(c);
                ((v & 0xF) as u8, v >> 4)
            })
            .chain(ops.iter().copied())
    }

    /// The number of reference bases the alignment spans (htslib's `bam_cigar2rlen`).
    pub fn reference_len(&self) -> i64 {
        self.iter().filter(|&(op, _)| op::consumes(op).1).map(|(_, len)| i64::from(len)).sum()
    }
}

impl fmt::Debug for Cigar<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (op, len) in self.iter() {
            write!(f, "{len}{}", b"MIDNSHP=XB".get(usize::from(op)).map_or('?', |&c| char::from(c)))?;
        }
        Ok(())
    }
}

/// SAM flag bits used here.
pub mod flag {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const READ2: u16 = 0x80;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// One alignment record, as depth counting sees it.
#[derive(Clone, Copy, Debug)]
pub struct Alignment<'a> {
    /// 0-based reference sequence index.
    pub tid: usize,
    /// 0-based leftmost position.
    pub pos: i64,
    pub flags: u16,
    /// Mapping quality; 255 when unavailable.
    pub mapq: u8,
    /// Mate's reference sequence index, if any.
    pub mate_tid: Option<usize>,
    /// Mate's 0-based position; -1 when unavailable.
    pub mate_pos: i64,
    /// TLEN.
    pub template_len: i32,
    pub name: &'a [u8],
    pub cigar: Cigar<'a>,
    /// The RG tag value. Only needed when filtering by read group.
    pub read_group: Option<&'a [u8]>,
}

impl Alignment<'_> {
    /// htslib's `bam_endpos`: one past the last reference base, and at least `pos + 1`.
    #[inline]
    pub fn end(&self) -> i64 {
        let rlen = if self.flags & flag::UNMAPPED != 0 { 0 } else { self.cigar.reference_len() };
        self.pos + rlen.max(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consumes_matches_htslib() {
        // M I D N S H P = X B
        let want = [(true, true), (true, false), (false, true), (false, true), (true, false)];
        for (op, w) in want.iter().enumerate() {
            assert_eq!(op::consumes(op as u8), *w, "op {op}");
        }
        assert_eq!(op::consumes(op::HARD_CLIP), (false, false));
        assert_eq!(op::consumes(op::PAD), (false, false));
        assert_eq!(op::consumes(op::EQUAL), (true, true));
        assert_eq!(op::consumes(op::DIFF), (true, true));
        assert_eq!(op::consumes(op::BACK), (false, false));
    }

    #[test]
    fn decodes_packed_cigar() {
        let words = [(10u32 << 4) | 4, (5 << 4), (2 << 4) | 2, (7 << 4) | 1];
        let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
        let c = Cigar::Packed(&bytes);
        assert_eq!(format!("{c:?}"), "10S5M2D7I");
        assert_eq!(c.reference_len(), 7);
    }

    #[test]
    fn end_is_at_least_one_base() {
        let ops = [(op::SOFT_CLIP, 10)];
        let a = Alignment {
            tid: 0,
            pos: 100,
            flags: 0,
            mapq: 60,
            mate_tid: None,
            mate_pos: -1,
            template_len: 0,
            name: b"r",
            cigar: Cigar::Ops(&ops),
            read_group: None,
        };
        assert_eq!(a.end(), 101);
    }
}

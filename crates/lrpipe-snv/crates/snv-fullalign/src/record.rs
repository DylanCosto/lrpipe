//! The parts of an alignment record and of the reference that the feature code reads.

/// CIGAR operation codes (`BAM_C*`).
pub(crate) mod op {
    pub(crate) const MATCH: u32 = 0;
    pub(crate) const INS: u32 = 1;
    pub(crate) const DEL: u32 = 2;
    pub(crate) const REF_SKIP: u32 = 3;
    pub(crate) const SOFT_CLIP: u32 = 4;
    pub(crate) const EQUAL: u32 = 7;
    pub(crate) const DIFF: u32 = 8;

    /// `M`, `=` or `X`.
    pub(crate) fn is_match(op: u32) -> bool {
        matches!(op, MATCH | EQUAL | DIFF)
    }
}

/// `seq_nt16_str`: the character of each 4-bit base code.
pub(crate) const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// `(bam_cigar_op(c), bam_cigar_oplen(c))`.
pub(crate) fn cigar_op(c: u32) -> (u32, u64) {
    (c & 0xf, u64::from(c >> 4))
}

/// One alignment record as `struct Read` sees it while it is current.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadView<'rec> {
    /// 0-based leftmost reference position.
    pub(crate) start: u64,
    /// Raw CIGAR words.
    pub(crate) cigar: &'rec [u32],
    /// 4-bit packed sequence (`bam_get_seq`).
    pub(crate) seq: &'rec [u8],
    /// Base qualities (`bam_get_qual`).
    pub(crate) qual: &'rec [u8],
    /// `l_qseq`.
    pub(crate) len: usize,
}

impl ReadView<'_> {
    /// `bam_seqi(seq, i)`. The original reads past the sequence when the CIGAR is longer than
    /// it (undefined behaviour); here such bases read as code 0.
    pub(crate) fn base_code(&self, i: u64) -> u8 {
        let Ok(i) = usize::try_from(i) else { return 0 };
        if i >= self.len {
            return 0;
        }
        self.seq.get(i >> 1).map_or(0, |byte| (byte >> ((!i & 1) << 2)) & 0xf)
    }

    /// Quality of base `i`, 0 past the sequence (see [`Self::base_code`]).
    pub(crate) fn qual(&self, i: u64) -> u8 {
        usize::try_from(i).ok().and_then(|i| self.qual.get(i)).copied().unwrap_or(0)
    }

    /// `get_read_end`: the start plus every `M=XDN` length.
    pub(crate) fn end(&self) -> u64 {
        self.cigar
            .iter()
            .map(|&c| cigar_op(c))
            .filter(|&(o, _)| op::is_match(o) || o == op::DEL || o == op::REF_SKIP)
            .fold(self.start, |pos, (_, len)| pos.wrapping_add(len))
    }
}

/// The reference sequence the original loads for a job: `faidx_fetch_seq` from `ref_start`.
#[derive(Debug, Clone)]
pub(crate) struct RefWindow {
    seq: Vec<u8>,
    start: u64,
}

impl RefWindow {
    pub(crate) fn new(seq: Vec<u8>, start: u64) -> Self {
        Self { seq, start }
    }

    /// `ref_seq[pos - ref_start]`. The C string ends in a NUL, which is what positions past
    /// the end read here (the original reads further only for reads past the contig end).
    pub(crate) fn at(&self, pos: u64) -> u8 {
        usize::try_from(pos.wrapping_sub(self.start)).ok().and_then(|i| self.seq.get(i)).copied().unwrap_or(0)
    }

    /// The C string starting at `pos - ref_start`, at most `max` bytes (`strncpy` or
    /// `%.*s` semantics: it stops early at the terminating NUL).
    pub(crate) fn c_str(&self, pos: u64, max: u64) -> &[u8] {
        let Ok(i) = usize::try_from(pos.wrapping_sub(self.start)) else { return &[] };
        let Some(rest) = self.seq.get(i..) else { return &[] };
        let n = usize::try_from(max).unwrap_or(usize::MAX).min(rest.len());
        let s = &rest[..n];
        s.iter().position(|&b| b == 0).map_or(s, |nul| &s[..nul])
    }
}

/// `upper_base` from `medaka_common.c`.
pub(crate) fn upper_base(c: u8) -> u8 {
    if c.is_ascii_lowercase() { c - 32 } else { c }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_window() {
        let r = RefWindow::new(b"ACGTacgt".to_vec(), 100);
        assert_eq!(r.at(100), b'A');
        assert_eq!(r.at(107), b't');
        assert_eq!(r.at(108), 0);
        assert_eq!(r.at(99), 0);
        assert_eq!(r.c_str(102, 3), b"GTa");
        assert_eq!(r.c_str(106, 10), b"gt");
        assert_eq!(r.c_str(200, 10), b"");
        assert_eq!(upper_base(b'g'), b'G');
        assert_eq!(upper_base(b'N'), b'N');
    }

    #[test]
    fn packed_bases() {
        // "ACGTN" packed two bases per byte, high nibble first.
        let seq = [0x12, 0x48, 0xf0];
        let read = ReadView { start: 10, cigar: &[(5 << 4) | op::MATCH], seq: &seq, qual: &[1, 2, 3, 4, 5], len: 5 };
        let bases: Vec<u8> = (0..6).map(|i| NT16[usize::from(read.base_code(i))]).collect();
        assert_eq!(bases, b"ACGTN=");
        assert_eq!(read.end(), 15);
        assert_eq!(read.qual(4), 5);
        assert_eq!(read.qual(5), 0);
    }
}

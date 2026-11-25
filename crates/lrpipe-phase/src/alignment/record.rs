//! An alignment record in BAM layout (htslib's `bam1_t`).
//!
//! Records read from BAM files keep their bytes, so records written back
//! (haplotag) are identical apart from the tags that are changed. Reading
//! and writing apply the same transformations as htslib's `bam_read1` and
//! `bam_write1`: the `bin` field is recomputed, a missing NUL after the read
//! name is added, and CIGARs with more than 65535 operations are moved
//! between the CG tag and the CIGAR field.

use std::fmt;

use anyhow::{Result, bail};

pub const FPAIRED: u16 = 0x1;
pub const FUNMAP: u16 = 0x4;
pub const FREVERSE: u16 = 0x10;
pub const FSECONDARY: u16 = 0x100;
pub const FDUP: u16 = 0x400;
pub const FSUPPLEMENTARY: u16 = 0x800;

/// CIGAR operation codes (`BAM_C*`)
pub const CMATCH: u8 = 0;
pub const CINS: u8 = 1;
pub const CDEL: u8 = 2;
pub const CREF_SKIP: u8 = 3;
pub const CSOFT_CLIP: u8 = 4;
pub const CEQUAL: u8 = 7;
pub const CDIFF: u8 = 8;

const CIGAR_CHARS: &[u8; 16] = b"MIDNSHP=XBXXXXXX";
const SEQ_CHARS: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// The two bases of each byte of a packed sequence (high nibble first).
const SEQ_PAIRS: [[u8; 2]; 256] = {
    let mut pairs = [[0u8; 2]; 256];
    let mut b = 0;
    while b < 256 {
        pairs[b] = [SEQ_CHARS[b >> 4], SEQ_CHARS[b & 0xf]];
        b += 1;
    }
    pairs
};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Record {
    pub tid: i32,
    pub pos: i32,
    pub bin: u16,
    pub mapq: u8,
    pub flag: u16,
    pub l_seq: i32,
    pub mtid: i32,
    pub mpos: i32,
    pub isize: i32,
    n_cigar: u32,
    /// Length of the read name including its NUL
    l_qname: usize,
    /// read name, CIGAR (`n_cigar` u32), sequence (4 bits per base),
    /// qualities, aux data
    data: Vec<u8>,
}

/// An aux field value (types as in the BAM file).
#[derive(Clone, Debug, PartialEq)]
pub enum Aux<'a> {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64),
    String(&'a str),
    Hex(&'a [u8]),
    /// `B` arrays: subtype and raw little-endian values
    Array(u8, &'a [u8]),
}

fn le_u32(b: &[u8]) -> u32 {
    u32::from_le_bytes(b[..4].try_into().unwrap())
}

fn le_i32(b: &[u8]) -> i32 {
    le_u32(b) as i32
}

/// `hts_reg2bin(beg, end, 14, 5)`
pub fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let (min_shift, n_lvls) = (14, 5);
    let mut l = n_lvls;
    let mut s = min_shift;
    let mut t = ((1i64 << ((n_lvls << 1) + n_lvls)) - 1) / 7;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u16;
        }
        l -= 1;
        s += 3;
        t -= 1i64 << ((l << 1) + l);
    }
    0
}

/// Size of an aux value of `typ` starting at `p` (after the type byte), or
/// None if it is malformed or runs past the end of the record (htslib's
/// `skip_aux`).
fn aux_value_size(typ: u8, p: &[u8]) -> Option<usize> {
    let size = match typ {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        b'Z' | b'H' => p.iter().position(|&c| c == 0)? + 1,
        b'B' => {
            let sub = *p.first()?;
            let n = le_u32(p.get(1..5)?) as usize;
            let size = match sub {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            n.checked_mul(size)?.checked_add(5)?
        }
        _ => return None,
    };
    (size <= p.len()).then_some(size)
}

impl Record {
    /// Parses a record (without its block size) as `bam_read1` does.
    pub fn from_bam(block: &[u8]) -> Result<Record> {
        if block.len() < 32 {
            bail!("truncated BAM record");
        }
        let x2 = le_u32(&block[8..]);
        let x3 = le_u32(&block[12..]);
        let mut rec = Record {
            tid: le_i32(&block[0..]),
            pos: le_i32(&block[4..]),
            bin: (x2 >> 16) as u16,
            mapq: (x2 >> 8) as u8,
            flag: (x3 >> 16) as u16,
            l_seq: le_i32(&block[16..]),
            mtid: le_i32(&block[20..]),
            mpos: le_i32(&block[24..]),
            isize: le_i32(&block[28..]),
            n_cigar: x3 & 0xffff,
            l_qname: (x2 & 0xff) as usize,
            data: block[32..].to_vec(),
        };
        let l_seq = usize::try_from(rec.l_seq).ok();
        let fits = l_seq.is_some_and(|l| {
            rec.l_qname >= 1 && rec.n_cigar as usize * 4 + rec.l_qname + l.div_ceil(2) + l <= rec.data.len()
        });
        if !fits {
            bail!("invalid BAM record");
        }
        if rec.data[rec.l_qname - 1] != 0 {
            // fixup_missing_qname_nul
            rec.data.insert(rec.l_qname, 0);
            rec.l_qname += 1;
        }
        rec.tag2cigar();
        if rec.n_cigar > 0 {
            let (rlen, qlen) = rec.cigar_lengths();
            let rlen = if rec.flag & FUNMAP != 0 || rlen == 0 { 1 } else { rlen };
            rec.bin = reg2bin(rec.pos as i64, rec.pos as i64 + rlen);
            if rec.l_seq > 0 && rec.flag & FUNMAP == 0 && qlen != rec.l_seq as i64 {
                bail!("CIGAR and query sequence lengths differ for {}", rec.qname_str());
            }
        }
        Ok(rec)
    }

    /// Sets `bin` as `bam_set1` computes it (htslib's CRAM decoding).
    pub fn set_bin_from_alignment(&mut self) {
        let rlen = if self.flag & FUNMAP != 0 { 0 } else { self.cigar_lengths().0 };
        let rlen = if rlen == 0 { 1 } else { rlen };
        self.bin = reg2bin(self.pos as i64, self.pos as i64 + rlen);
    }

    /// `bam_tag2cigar`: a CIGAR stored in the CG tag (`<l_seq>S` CIGAR).
    fn tag2cigar(&mut self) {
        let test = (self.l_seq as u32) << 4 | CSOFT_CLIP as u32;
        if self.n_cigar == 0 || le_u32(&self.data[self.l_qname..]) != test || self.tid < 0 || self.pos < 0 {
            return;
        }
        let Some(cg) = self.aux_offset(b"CG") else { return };
        let d = &self.data;
        if d[cg + 2] != b'B' || !(d[cg + 3] == b'I' || d[cg + 3] == b'i') {
            return;
        }
        let cg_len = le_u32(&d[cg + 4..]);
        if cg_len < self.n_cigar || cg_len >= 1 << 29 {
            return;
        }
        let cigar_start = self.l_qname;
        let fake = self.n_cigar as usize * 4;
        let real: Vec<u8> = d[cg + 8..cg + 8 + cg_len as usize * 4].to_vec();
        let cg_end = cg + 8 + real.len();
        let mut data = Vec::with_capacity(d.len() - fake - 8);
        data.extend_from_slice(&d[..cigar_start]);
        data.extend_from_slice(&real);
        data.extend_from_slice(&d[cigar_start + fake..cg]);
        data.extend_from_slice(&d[cg_end..]);
        self.data = data;
        self.n_cigar = cg_len;
    }

    /// The record in BAM layout with its block size, as `bam_write1` writes it.
    pub fn write_bam(&self, out: &mut Vec<u8>) -> Result<()> {
        if self.l_qname > 255 {
            bail!("QNAME \"{}\" is longer than 254 characters", self.qname_str());
        }
        let long_cigar = self.n_cigar > 0xffff;
        let block_len = 32 + self.data.len() + if long_cigar { 16 } else { 0 };
        out.extend_from_slice(&(block_len as u32).to_le_bytes());
        out.extend_from_slice(&self.tid.to_le_bytes());
        out.extend_from_slice(&self.pos.to_le_bytes());
        let x2 = (self.bin as u32) << 16 | (self.mapq as u32) << 8 | self.l_qname as u32;
        out.extend_from_slice(&x2.to_le_bytes());
        let n = if long_cigar { 2 } else { self.n_cigar };
        out.extend_from_slice(&((self.flag as u32) << 16 | n).to_le_bytes());
        out.extend_from_slice(&self.l_seq.to_le_bytes());
        out.extend_from_slice(&self.mtid.to_le_bytes());
        out.extend_from_slice(&self.mpos.to_le_bytes());
        out.extend_from_slice(&self.isize.to_le_bytes());
        out.extend_from_slice(&self.data[..self.l_qname]);
        if !long_cigar {
            out.extend_from_slice(&self.data[self.l_qname..]);
        } else {
            let (rlen, _) = self.cigar_lengths();
            if rlen >= 1 << 28 {
                bail!(
                    "Record {} with {} CIGAR ops and ref length {rlen} cannot be written in BAM.  Try writing SAM or CRAM instead.",
                    self.qname_str(),
                    self.n_cigar
                );
            }
            let cigar_end = self.l_qname + self.n_cigar as usize * 4;
            out.extend_from_slice(&((self.l_seq as u32) << 4 | CSOFT_CLIP as u32).to_le_bytes());
            out.extend_from_slice(&((rlen as u32) << 4 | CREF_SKIP as u32).to_le_bytes());
            out.extend_from_slice(&self.data[cigar_end..]);
            out.extend_from_slice(b"CGBI");
            out.extend_from_slice(&self.n_cigar.to_le_bytes());
            out.extend_from_slice(&self.data[self.l_qname..cigar_end]);
        }
        Ok(())
    }

    pub fn qname(&self) -> &[u8] {
        let name = &self.data[..self.l_qname];
        let end = name.iter().position(|&c| c == 0).unwrap_or(name.len());
        &name[..end]
    }

    fn qname_str(&self) -> String {
        String::from_utf8_lossy(self.qname()).into_owned()
    }

    pub fn pos(&self) -> i64 {
        self.pos as i64
    }

    pub fn tid(&self) -> i32 {
        self.tid
    }

    pub fn mapq(&self) -> u8 {
        self.mapq
    }

    pub fn flags(&self) -> u16 {
        self.flag
    }

    pub fn is_paired(&self) -> bool {
        self.flag & FPAIRED != 0
    }

    pub fn is_unmapped(&self) -> bool {
        self.flag & FUNMAP != 0
    }

    pub fn is_reverse(&self) -> bool {
        self.flag & FREVERSE != 0
    }

    pub fn is_secondary(&self) -> bool {
        self.flag & FSECONDARY != 0
    }

    pub fn is_duplicate(&self) -> bool {
        self.flag & FDUP != 0
    }

    pub fn is_supplementary(&self) -> bool {
        self.flag & FSUPPLEMENTARY != 0
    }

    pub fn cigar(&self) -> CigarView<'_> {
        CigarView { raw: &self.data[self.l_qname..self.l_qname + self.n_cigar as usize * 4], pos: self.pos as i64 }
    }

    /// `bam_cigar2rqlens`: (reference length, query length)
    fn cigar_lengths(&self) -> (i64, i64) {
        let (mut rlen, mut qlen) = (0, 0);
        for (op, len) in self.cigar().iter() {
            let len = len as i64;
            match op {
                CMATCH | CEQUAL | CDIFF => {
                    rlen += len;
                    qlen += len;
                }
                CINS | CSOFT_CLIP => qlen += len,
                CDEL | CREF_SKIP => rlen += len,
                _ => {}
            }
        }
        (rlen, qlen)
    }

    /// `bam_endpos`: end of the alignment on the reference (pos + 1 for
    /// unmapped reads or reads without reference-consuming operations).
    pub fn end_pos(&self) -> i64 {
        let (rlen, _) = self.cigar_lengths();
        if self.flag & FUNMAP != 0 || rlen == 0 { self.pos as i64 + 1 } else { self.pos as i64 + rlen }
    }

    fn seq_offset(&self) -> usize {
        self.l_qname + self.n_cigar as usize * 4
    }

    fn qual_offset(&self) -> usize {
        self.seq_offset() + (self.l_seq.max(0) as usize).div_ceil(2)
    }

    fn aux_start(&self) -> usize {
        self.qual_offset() + self.l_seq.max(0) as usize
    }

    /// The bases as ASCII characters.
    pub fn seq(&self) -> Vec<u8> {
        let packed = &self.data[self.seq_offset()..self.qual_offset()];
        let length = self.l_seq.max(0) as usize;
        let mut seq = Vec::with_capacity(length + 1);
        for &b in &packed[..length / 2] {
            seq.extend_from_slice(&SEQ_PAIRS[b as usize]);
        }
        if length % 2 == 1 {
            seq.push(SEQ_CHARS[(packed[length / 2] >> 4) as usize]);
        }
        seq
    }

    /// Base qualities (255 if missing).
    pub fn qual(&self) -> &[u8] {
        &self.data[self.qual_offset()..self.aux_start()]
    }

    /// Offset of the aux field `tag` (its tag bytes) in `data`.
    fn aux_offset(&self, tag: &[u8; 2]) -> Option<usize> {
        let d = &self.data;
        let mut p = self.aux_start();
        while p + 3 <= d.len() {
            let size = aux_value_size(d[p + 2], &d[p + 3..])?;
            if &d[p..p + 2] == tag {
                return Some(p);
            }
            p += 3 + size;
        }
        None
    }

    /// `bam_aux_get` + conversion (Err if absent or malformed).
    pub fn aux(&self, tag: &[u8; 2]) -> Result<Aux<'_>> {
        let Some(p) = self.aux_offset(tag) else { bail!("tag {} not present", String::from_utf8_lossy(tag)) };
        let d = &self.data[p + 3..];
        Ok(match self.data[p + 2] {
            b'A' => Aux::Char(d[0]),
            b'c' => Aux::I8(d[0] as i8),
            b'C' => Aux::U8(d[0]),
            b's' => Aux::I16(i16::from_le_bytes([d[0], d[1]])),
            b'S' => Aux::U16(u16::from_le_bytes([d[0], d[1]])),
            b'i' => Aux::I32(le_i32(d)),
            b'I' => Aux::U32(le_u32(d)),
            b'f' => Aux::Float(f32::from_le_bytes(d[..4].try_into().unwrap())),
            b'd' => Aux::Double(f64::from_le_bytes(d[..8].try_into().unwrap())),
            b'Z' => {
                let end = d.iter().position(|&c| c == 0).unwrap_or(d.len());
                Aux::String(std::str::from_utf8(&d[..end])?)
            }
            b'H' => Aux::Hex(&d[..d.iter().position(|&c| c == 0).unwrap_or(d.len())]),
            b'B' => {
                let size = aux_value_size(b'B', d).expect("aux_offset checked the size");
                Aux::Array(d[0], &d[5..size])
            }
            t => bail!("unknown aux type {}", t as char),
        })
    }

    /// Tag and offset (in `data`) of the last aux field.
    fn last_aux(&self) -> Option<([u8; 2], usize)> {
        let d = &self.data;
        let mut p = self.aux_start();
        let mut last = None;
        while p + 3 <= d.len() {
            let size = aux_value_size(d[p + 2], &d[p + 3..])?;
            last = Some(([d[p], d[p + 1]], p));
            p += 3 + size;
        }
        last
    }

    /// Inserts encoded aux fields before a final RG field (where htslib puts
    /// tags generated while decoding CRAM), or at the end.
    pub fn insert_aux_before_rg(&mut self, fields: &[u8]) {
        let at = match self.last_aux() {
            Some((tag, p)) if &tag == b"RG" => p,
            _ => self.data.len(),
        };
        self.data.splice(at..at, fields.iter().copied());
    }

    /// `bam_aux_del` of `tag` if present.
    pub fn remove_aux(&mut self, tag: &[u8; 2]) {
        if let Some(p) = self.aux_offset(tag) {
            let size = aux_value_size(self.data[p + 2], &self.data[p + 3..]).unwrap_or(0);
            self.data.drain(p..p + 3 + size);
        }
    }

    /// `bam_aux_append` of an integer with the type pysam's `set_tag` picks
    /// (which refuses values outside the 32-bit BAM integer types).
    pub fn push_aux_int(&mut self, tag: &[u8; 2], v: i64) -> Result<()> {
        if v < i32::MIN as i64 || v > u32::MAX as i64 {
            bail!("ValueError: at least one integer out of range of BAM/SAM specification");
        }
        self.data.extend_from_slice(tag);
        if v < 0 {
            if v >= i8::MIN as i64 {
                self.data.push(b'c');
                self.data.push(v as i8 as u8);
            } else if v >= i16::MIN as i64 {
                self.data.push(b's');
                self.data.extend_from_slice(&(v as i16).to_le_bytes());
            } else {
                self.data.push(b'i');
                self.data.extend_from_slice(&(v as i32).to_le_bytes());
            }
        } else if v <= u8::MAX as i64 {
            self.data.push(b'C');
            self.data.push(v as u8);
        } else if v <= u16::MAX as i64 {
            self.data.push(b'S');
            self.data.extend_from_slice(&(v as u16).to_le_bytes());
        } else {
            self.data.push(b'I');
            self.data.extend_from_slice(&(v as u32).to_le_bytes());
        }
        Ok(())
    }
}

/// The CIGAR of a record.
#[derive(Clone, Copy)]
pub struct CigarView<'a> {
    raw: &'a [u8],
    pos: i64,
}

impl CigarView<'_> {
    /// (operation code, length)
    pub fn iter(&self) -> impl Iterator<Item = (u8, u32)> + '_ {
        self.raw.as_chunks::<4>().0.iter().map(|c| {
            let v = u32::from_le_bytes(*c);
            ((v & 0xf) as u8, v >> 4)
        })
    }

    pub fn len(&self) -> usize {
        self.raw.len() / 4
    }

    pub fn is_empty(&self) -> bool {
        self.raw.is_empty()
    }

    /// End position of the CIGAR on the reference (from the alignment start).
    pub fn end_pos(&self) -> i64 {
        self.pos
            + self
                .iter()
                .filter(|(op, _)| matches!(*op, CMATCH | CDEL | CREF_SKIP | CEQUAL | CDIFF))
                .map(|(_, l)| l as i64)
                .sum::<i64>()
    }
}

impl fmt::Display for CigarView<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (op, len) in self.iter() {
            write!(f, "{len}{}", CIGAR_CHARS[op as usize] as char)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A record with one base, CIGAR 1S and the given aux data.
    fn block(aux: &[u8]) -> Vec<u8> {
        let mut b = Vec::new();
        for v in [0i32, 0] {
            b.extend_from_slice(&v.to_le_bytes());
        }
        b.extend_from_slice(&2u32.to_le_bytes()); // bin, mapq 0, l_qname 2
        b.extend_from_slice(&1u32.to_le_bytes()); // flag 0, one CIGAR operation
        for v in [1i32, -1, -1, 0] {
            b.extend_from_slice(&v.to_le_bytes());
        }
        b.extend_from_slice(b"r\0");
        b.extend_from_slice(&(1u32 << 4 | CSOFT_CLIP as u32).to_le_bytes());
        b.extend_from_slice(&[0x10, 30]);
        b.extend_from_slice(aux);
        b
    }

    #[test]
    fn truncated_aux_fields_are_ignored() {
        // a CG array claiming 1000 values without any
        let mut aux = b"CGBI".to_vec();
        aux.extend_from_slice(&1000u32.to_le_bytes());
        let rec = Record::from_bam(&block(&aux)).unwrap();
        assert_eq!(rec.cigar().to_string(), "1S");
        assert!(rec.aux(b"CG").is_err());
        // an integer tag cut short
        let rec = Record::from_bam(&block(b"XXi\x01\x00")).unwrap();
        assert!(rec.aux(b"XX").is_err());
        let rec = Record::from_bam(&block(b"XXi\x01\x00\x00\x00")).unwrap();
        assert_eq!(rec.aux(b"XX").unwrap(), Aux::I32(1));
    }

    #[test]
    fn integer_tags_out_of_range_are_refused() {
        let mut rec = Record::from_bam(&block(b"")).unwrap();
        assert!(rec.push_aux_int(b"PS", u32::MAX as i64 + 1).is_err());
        assert!(rec.push_aux_int(b"PS", i32::MIN as i64 - 1).is_err());
        rec.push_aux_int(b"PS", u32::MAX as i64).unwrap();
        assert_eq!(rec.aux(b"PS").unwrap(), Aux::U32(u32::MAX));
    }
}

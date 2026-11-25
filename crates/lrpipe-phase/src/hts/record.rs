//! A VCF/BCF record as htslib 1.23 models it (`bcf1_t`), with its VCF
//! parser (`vcf_parse`), VCF formatter (`vcf_format`), BCF encoding and the
//! field updates (`bcf_update_info`, `bcf_update_format`).
//!
//! INFO and FORMAT fields keep their BCF encoding (type code and
//! little-endian values), which decides how values are printed and which
//! integer width pysam's setters write into.

use anyhow::{Result, bail};

use super::header::{
    BCF_DT_CTG, BCF_DT_ID, BCF_HL_FMT, BCF_HL_INFO, BCF_HT_FLAG, BCF_HT_INT, BCF_HT_REAL, BCF_HT_STR, HeaderCore, VCF44,
};
use super::num::{atof, put_double, put_int, str2dbl, str2int, str2uint};

pub const BCF_BT_NULL: u8 = 0;
pub const BCF_BT_INT8: u8 = 1;
pub const BCF_BT_INT16: u8 = 2;
pub const BCF_BT_INT32: u8 = 3;
pub const BCF_BT_INT64: u8 = 4;
pub const BCF_BT_FLOAT: u8 = 5;
pub const BCF_BT_CHAR: u8 = 7;

const TYPE_SHIFT: [u32; 16] = [0, 0, 1, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

pub const INT8_MISSING: i8 = i8::MIN;
pub const INT8_VECTOR_END: i8 = i8::MIN + 1;
pub const INT16_MISSING: i16 = i16::MIN;
pub const INT16_VECTOR_END: i16 = i16::MIN + 1;
pub const INT32_MISSING: i32 = i32::MIN;
pub const INT32_VECTOR_END: i32 = i32::MIN + 1;
const INT64_MISSING: i64 = i64::MIN;
pub const FLOAT_MISSING: u32 = 0x7F80_0001;
pub const FLOAT_VECTOR_END: u32 = 0x7F80_0002;
const STR_MISSING: u8 = 0x07;

const MIN_INT8: i32 = -120;
const MAX_INT8: i32 = 127;
const MIN_INT16: i32 = -32760;
const MAX_INT16: i32 = 32767;
const MIN_INT32: i64 = -2147483640;
const MAX_INT32: i64 = 2147483647;

pub const BCF_ERR_CTG_UNDEF: i32 = 1;
pub const BCF_ERR_TAG_UNDEF: i32 = 2;
pub const BCF_ERR_NCOLS: i32 = 4;
pub const BCF_ERR_LIMITS: i32 = 8;
pub const BCF_ERR_CHAR: i32 = 16;
pub const BCF_ERR_CTG_INVALID: i32 = 32;
pub const BCF_ERR_TAG_INVALID: i32 = 64;

const MAX_N_FMT: usize = 255;

// --- typed value encoding (bcf_enc_*) ---------------------------------------

/// `bcf_enc_size`
pub fn enc_size(s: &mut Vec<u8>, size: usize, typ: u8) {
    if size < 15 {
        s.push(((size as u8) << 4) | typ);
    } else {
        s.push(15 << 4 | typ);
        if size < 128 {
            s.push(1 << 4 | BCF_BT_INT8);
            s.push(size as u8);
        } else if size < 32768 {
            s.push(1 << 4 | BCF_BT_INT16);
            s.extend_from_slice(&(size as i16).to_le_bytes());
        } else {
            s.push(1 << 4 | BCF_BT_INT32);
            s.extend_from_slice(&(size as i32).to_le_bytes());
        }
    }
}

/// `bcf_enc_int1`
pub fn enc_int1(s: &mut Vec<u8>, x: i32) {
    if x == INT32_VECTOR_END {
        s.extend_from_slice(&[1 << 4 | BCF_BT_INT8, INT8_VECTOR_END as u8]);
    } else if x == INT32_MISSING {
        s.extend_from_slice(&[1 << 4 | BCF_BT_INT8, INT8_MISSING as u8]);
    } else if (MIN_INT8..=MAX_INT8).contains(&x) {
        s.extend_from_slice(&[1 << 4 | BCF_BT_INT8, x as u8]);
    } else if (MIN_INT16..=MAX_INT16).contains(&x) {
        s.push(1 << 4 | BCF_BT_INT16);
        s.extend_from_slice(&(x as i16).to_le_bytes());
    } else {
        s.push(1 << 4 | BCF_BT_INT32);
        s.extend_from_slice(&x.to_le_bytes());
    }
}

/// `bcf_enc_vint`
pub fn enc_vint(s: &mut Vec<u8>, a: &[i32], wsize: i32) {
    let n = a.len();
    if n == 0 {
        enc_size(s, 0, BCF_BT_NULL);
        return;
    }
    if n == 1 {
        enc_int1(s, a[0]);
        return;
    }
    let wsize = if wsize <= 0 { n } else { wsize as usize };
    let mut max = i32::MIN;
    let mut min = i32::MAX;
    for &x in a {
        max = max.max(x);
        if min > x && x > INT32_VECTOR_END {
            min = x;
        }
    }
    if max <= MAX_INT8 && min >= MIN_INT8 {
        enc_size(s, wsize, BCF_BT_INT8);
        s.extend(a.iter().map(|&x| match x {
            INT32_VECTOR_END => INT8_VECTOR_END as u8,
            INT32_MISSING => INT8_MISSING as u8,
            _ => x as u8,
        }));
    } else if max <= MAX_INT16 && min >= MIN_INT16 {
        enc_size(s, wsize, BCF_BT_INT16);
        for &x in a {
            let v = match x {
                INT32_VECTOR_END => INT16_VECTOR_END,
                INT32_MISSING => INT16_MISSING,
                _ => x as i16,
            };
            s.extend_from_slice(&v.to_le_bytes());
        }
    } else {
        enc_size(s, wsize, BCF_BT_INT32);
        for &x in a {
            s.extend_from_slice(&x.to_le_bytes());
        }
    }
}

/// `bcf_enc_vfloat` (values as bit patterns, to keep missing values)
pub fn enc_vfloat(s: &mut Vec<u8>, a: &[u32]) {
    enc_size(s, a.len(), BCF_BT_FLOAT);
    for &x in a {
        s.extend_from_slice(&x.to_le_bytes());
    }
}

/// `bcf_enc_vchar`
pub fn enc_vchar(s: &mut Vec<u8>, a: &[u8]) {
    enc_size(s, a.len(), BCF_BT_CHAR);
    s.extend_from_slice(a);
}

// --- decoding ---------------------------------------------------------------

/// A bounds-checked reader over a BCF data block.
struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl Cursor<'_> {
    fn byte(&mut self) -> Result<u8> {
        let Some(&b) = self.buf.get(self.pos) else { bail!("truncated BCF record") };
        self.pos += 1;
        Ok(b)
    }

    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let Some(s) = self.buf.get(self.pos..self.pos.saturating_add(n)) else { bail!("truncated BCF record") };
        self.pos += n;
        Ok(s)
    }

    /// `bcf_dec_int1`
    fn int1(&mut self, typ: u8) -> Result<i64> {
        Ok(match typ {
            BCF_BT_INT8 => self.byte()? as i8 as i64,
            BCF_BT_INT16 => i16::from_le_bytes(self.take(2)?.try_into()?) as i64,
            BCF_BT_INT32 => i32::from_le_bytes(self.take(4)?.try_into()?) as i64,
            BCF_BT_INT64 => i64::from_le_bytes(self.take(8)?.try_into()?),
            _ => bail!("invalid BCF integer type {typ}"),
        })
    }

    /// `bcf_dec_typed_int1`
    fn typed_int1(&mut self) -> Result<i64> {
        let t = self.byte()? & 0xf;
        self.int1(t)
    }

    /// `bcf_dec_size`: (size, type)
    fn size(&mut self) -> Result<(usize, u8)> {
        let b = self.byte()?;
        let typ = b & 0xf;
        if b >> 4 != 15 {
            return Ok(((b >> 4) as usize, typ));
        }
        let n = self.typed_int1()?;
        if n < 0 {
            bail!("invalid BCF vector size {n}");
        }
        Ok((n as usize, typ))
    }
}

fn le_i16(p: &[u8], i: usize) -> i16 {
    i16::from_le_bytes([p[2 * i], p[2 * i + 1]])
}

fn le_i32(p: &[u8], i: usize) -> i32 {
    i32::from_le_bytes([p[4 * i], p[4 * i + 1], p[4 * i + 2], p[4 * i + 3]])
}

fn le_u32(p: &[u8], i: usize) -> u32 {
    le_i32(p, i) as u32
}

fn le_i64(p: &[u8], i: usize) -> i64 {
    i64::from_le_bytes(p[8 * i..8 * i + 8].try_into().expect("8 bytes"))
}

/// Integer `i` of a vector of type `typ` (`None` for non-integer types).
fn int_at(typ: u8, p: &[u8], i: usize) -> Option<i64> {
    Some(match typ {
        BCF_BT_INT8 => p[i] as i8 as i64,
        BCF_BT_INT16 => le_i16(p, i) as i64,
        BCF_BT_INT32 => le_i32(p, i) as i64,
        BCF_BT_INT64 => le_i64(p, i),
        _ => return None,
    })
}

/// Like `int_at`, with missing values as 0 (as `get_rlen` reads them).
fn int_or_zero(typ: u8, p: &[u8], i: usize) -> Option<i64> {
    let v = int_at(typ, p, i)?;
    let missing = match typ {
        BCF_BT_INT8 => INT8_MISSING as i64,
        BCF_BT_INT16 => INT16_MISSING as i64,
        BCF_BT_INT32 => INT32_MISSING as i64,
        _ => INT64_MISSING,
    };
    Some(if v == missing { 0 } else { v })
}

/// `bcf_fmt_array`
fn fmt_array(s: &mut Vec<u8>, n: usize, typ: u8, data: &[u8]) {
    if n == 0 {
        s.push(b'.');
        return;
    }
    if typ == BCF_BT_CHAR {
        let d = &data[..n.min(data.len())];
        let end = d.iter().position(|&c| c == 0).unwrap_or(d.len());
        s.extend_from_slice(&d[..end]);
        return;
    }
    for j in 0..n {
        let (missing, vector_end) = match typ {
            BCF_BT_INT8 => (data[j] as i8 == INT8_MISSING, data[j] as i8 == INT8_VECTOR_END),
            BCF_BT_INT16 => (le_i16(data, j) == INT16_MISSING, le_i16(data, j) == INT16_VECTOR_END),
            BCF_BT_INT32 => (le_i32(data, j) == INT32_MISSING, le_i32(data, j) == INT32_VECTOR_END),
            BCF_BT_FLOAT => (le_u32(data, j) == FLOAT_MISSING, le_u32(data, j) == FLOAT_VECTOR_END),
            _ => panic!("Unexpected type {typ}"),
        };
        if vector_end {
            break;
        }
        if j > 0 {
            s.push(b',');
        }
        if missing {
            s.push(b'.');
        } else if typ == BCF_BT_FLOAT {
            put_double(f32::from_bits(le_u32(data, j)) as f64, s);
        } else {
            put_int(int_at(typ, data, j).expect("integer type"), s);
        }
    }
}

/// `bcf_fmt_array1` (a FORMAT field with one value per sample)
fn fmt_array1(s: &mut Vec<u8>, typ: u8, p: &[u8]) {
    match typ {
        BCF_BT_CHAR => s.push(if p[0] == STR_MISSING { b'.' } else { p[0] }),
        BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 => {
            let (v, missing, vector_end) = match typ {
                BCF_BT_INT8 => (p[0] as i8 as i64, INT8_MISSING as i64, INT8_VECTOR_END as i64),
                BCF_BT_INT16 => (le_i16(p, 0) as i64, INT16_MISSING as i64, INT16_VECTOR_END as i64),
                _ => (le_i32(p, 0) as i64, INT32_MISSING as i64, INT32_VECTOR_END as i64),
            };
            if v != vector_end {
                if v == missing {
                    s.push(b'.');
                } else {
                    put_int(v, s);
                }
            }
        }
        BCF_BT_FLOAT => {
            let v = le_u32(p, 0);
            if v != FLOAT_VECTOR_END {
                if v == FLOAT_MISSING {
                    s.push(b'.');
                } else {
                    put_double(f32::from_bits(v) as f64, s);
                }
            }
        }
        _ => log::error!("Unexpected type {typ}"),
    }
}

/// An INFO field (`bcf_info_t`): `raw` is its BCF encoding (key, type and
/// size, values), the values starting at `off`.
#[derive(Clone, Debug)]
pub struct Info {
    pub key: i32,
    pub typ: u8,
    pub len: i32,
    raw: Vec<u8>,
    off: usize,
    /// marked for removal (`vptr == NULL`)
    pub deleted: bool,
}

impl Info {
    /// `bcf_unpack_info_core1` on an encoded field.
    fn decode(c: &mut Cursor) -> Result<Info> {
        let start = c.pos;
        let key = c.typed_int1()? as i32;
        let (len, typ) = c.size()?;
        let off = c.pos - start;
        c.take(len << TYPE_SHIFT[typ as usize])?;
        Ok(Info { key, typ, len: len as i32, raw: c.buf[start..c.pos].to_vec(), off, deleted: false })
    }

    fn from_encoded(raw: Vec<u8>) -> Result<Info> {
        let mut c = Cursor { buf: &raw, pos: 0 };
        let info = Info::decode(&mut c)?;
        Ok(info)
    }

    pub fn data(&self) -> &[u8] {
        &self.raw[self.off..]
    }
}

/// A FORMAT field (`bcf_fmt_t`): `n` values of type `typ` per sample, `size`
/// bytes per sample, after the encoded key, type and size in `raw[..off]`.
#[derive(Clone, Debug)]
pub struct Fmt {
    pub id: i32,
    pub typ: u8,
    pub n: i32,
    pub size: usize,
    raw: Vec<u8>,
    off: usize,
    /// marked for removal (`p == NULL`)
    pub deleted: bool,
}

impl Fmt {
    /// `bcf_unpack_fmt_core1`
    fn decode(c: &mut Cursor, n_sample: usize) -> Result<Fmt> {
        let start = c.pos;
        let id = c.typed_int1()? as i32;
        let (n, typ) = c.size()?;
        let size = n << TYPE_SHIFT[typ as usize];
        let off = c.pos - start;
        c.take(n_sample * size)?;
        Ok(Fmt { id, typ, n: n as i32, size, raw: c.buf[start..c.pos].to_vec(), off, deleted: false })
    }

    pub fn p(&self) -> &[u8] {
        &self.raw[self.off..]
    }

    pub fn p_mut(&mut self) -> &mut [u8] {
        &mut self.raw[self.off..]
    }

    pub fn sample(&self, j: usize) -> &[u8] {
        &self.p()[j * self.size..(j + 1) * self.size]
    }

    pub fn sample_mut(&mut self, j: usize) -> &mut [u8] {
        let size = self.size;
        &mut self.p_mut()[j * size..(j + 1) * size]
    }
}

/// Values for `Record::update_format`.
pub enum FmtValues<'a> {
    Int(&'a [i32]),
    /// float bit patterns
    Float(&'a [u32]),
    Str(&'a [u8]),
}

/// Values for `Record::update_info`.
pub enum InfoValues<'a> {
    Int(&'a [i32]),
}

#[derive(Clone, Debug)]
pub struct Record {
    pub rid: i32,
    pub pos: i64,
    pub rlen: i64,
    /// QUAL as float bits
    pub qual: u32,
    pub n_sample: u32,
    /// ID as printed ("." if missing)
    pub id: Vec<u8>,
    pub alleles: Vec<Vec<u8>>,
    pub flt: Vec<i32>,
    pub info: Vec<Info>,
    pub fmt: Vec<Fmt>,
    pub errcode: i32,
    /// positions or values that do not fit BCF (`BCF_IS_64BIT`)
    pub is_64bit: bool,
    shared_dirty: bool,
}

impl Default for Record {
    fn default() -> Self {
        Record {
            rid: 0,
            pos: 0,
            rlen: 0,
            qual: FLOAT_MISSING,
            n_sample: 0,
            id: Vec::new(),
            alleles: Vec::new(),
            flt: Vec::new(),
            info: Vec::new(),
            fmt: Vec::new(),
            errcode: 0,
            is_64bit: false,
            shared_dirty: false,
        }
    }
}

/// htslib's `svlen_on_ref_for_vcf_alt`: whether SVLEN gives the reference
/// length of the ALT allele (<CNV>, <DEL>, <DUP>, <INV> and subtypes such as
/// <DEL:ME>).
pub fn svlen_on_ref_for_vcf_alt(alt: &[u8]) -> bool {
    alt.len() >= 5
        && (alt[4] == b'>' || alt[4] == b':')
        && [b"<CNV", b"<DEL", b"<DUP", b"<INV"].iter().any(|p| alt.starts_with(*p))
        && alt.ends_with(b">")
}

/// The fixed-size part of a BCF record (`bcf_read1_core`).
pub struct BcfCore {
    pub rid: i32,
    pub pos: i64,
    pub rlen: i64,
    pub qual: u32,
    pub n_info: usize,
    pub n_allele: usize,
    pub n_sample: u32,
    pub n_fmt: usize,
}

/// A C-string token tokenizer over a buffer that is NUL-terminated (htslib's
/// `kstrtok` with one separator; separators are overwritten with NULs).
struct Tok {
    pos: usize,
    finished: bool,
}

impl Tok {
    fn next(&mut self, s: &mut [u8], sep: u8) -> Option<(usize, usize)> {
        if self.finished {
            return None;
        }
        let start = self.pos;
        let mut p = start;
        while s[p] != 0 && s[p] != sep {
            p += 1;
        }
        if s[p] == 0 {
            self.finished = true;
        }
        s[p] = 0;
        self.pos = p + 1;
        Some((start, p))
    }
}

/// The C string starting at `p`.
fn cstr(s: &[u8], p: usize) -> &[u8] {
    let end = s[p..].iter().position(|&c| c == 0).map_or(s.len(), |e| p + e);
    &s[p..end]
}

fn lossy(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

/// Adds a header line for an undeclared contig/tag (the error recovery of
/// `vcf_parse`): whether the header now has it.
fn add_dummy(h: &mut HeaderCore, line: &str) -> bool {
    let (hrec, _) = HeaderCore::parse_line(line.as_bytes());
    let Some(hrec) = hrec else { return false };
    match h.add_hrec(hrec) {
        Ok(r) => {
            if r > 0 {
                h.sync();
            }
            true
        }
        Err(_) => false,
    }
}

/// Per-FORMAT-field state of `vcf_parse_format` (`fmt_aux_t`).
#[derive(Clone, Default)]
struct FmtAux {
    key: i32,
    max_m: usize,
    /// bytes per sample, `None` for fields dropped (duplicates)
    size: Option<usize>,
    is_gt: bool,
    max_g: usize,
    max_l: usize,
    y: u32,
    buf: Vec<u8>,
}

impl Record {
    pub fn new() -> Record {
        Record::default()
    }

    pub fn n_allele(&self) -> usize {
        self.alleles.len()
    }

    /// Decodes a BCF record (`bcf_unpack` of everything).
    pub fn from_bcf(core: &BcfCore, shared: &[u8], indiv: &[u8]) -> Result<Record> {
        let mut r = Record {
            rid: core.rid,
            pos: core.pos,
            rlen: core.rlen,
            qual: core.qual,
            n_sample: core.n_sample,
            ..Record::default()
        };
        let mut c = Cursor { buf: shared, pos: 0 };
        let sized_string = |c: &mut Cursor| -> Result<Vec<u8>> {
            let (n, typ) = c.size()?;
            let data = c.take(n << TYPE_SHIFT[typ as usize])?;
            let mut out = Vec::new();
            fmt_array(&mut out, n, typ, data);
            Ok(out)
        };
        r.id = sized_string(&mut c)?;
        for _ in 0..core.n_allele {
            r.alleles.push(sized_string(&mut c)?);
        }
        if c.buf.get(c.pos).is_some_and(|&b| b >> 4 != 0) {
            let (n, typ) = c.size()?;
            for _ in 0..n {
                r.flt.push(c.int1(typ)? as i32);
            }
        } else {
            c.byte()?;
        }
        for _ in 0..core.n_info {
            r.info.push(Info::decode(&mut c)?);
        }
        if core.n_sample > 0 {
            let mut c = Cursor { buf: indiv, pos: 0 };
            for _ in 0..core.n_fmt {
                r.fmt.push(Fmt::decode(&mut c, core.n_sample as usize)?);
            }
        }
        Ok(r)
    }

    /// `vcf_parse`: parses a VCF line (without the newline), adding header
    /// lines for undeclared contigs and tags like htslib. Returns the record
    /// and 0, or -2 on errors (with `errcode` set for most of them).
    pub fn parse_vcf(line: &[u8], h: &mut HeaderCore) -> (Record, i32) {
        let mut v = Record::new();
        let mut shared = Vec::new();
        let mut indiv = Vec::new();
        let mut core =
            BcfCore { rid: 0, pos: 0, rlen: 0, qual: FLOAT_MISSING, n_info: 0, n_allele: 0, n_sample: 0, n_fmt: 0 };
        let ok = v.parse_vcf_inner(line, h, &mut core, &mut shared, &mut indiv);
        if ok.is_none() {
            return (v, -2);
        }
        let (errcode, is_64bit) = (v.errcode, v.is_64bit);
        let len_ref = core.rlen;
        match Record::from_bcf(&core, &shared, &indiv) {
            Ok(mut r) => {
                r.errcode = errcode;
                r.is_64bit = is_64bit;
                r.rlen = r.get_rlen(h, len_ref);
                (r, 0)
            }
            Err(_) => (v, -2),
        }
    }

    fn parse_vcf_inner(
        &mut self,
        line: &[u8],
        h: &mut HeaderCore,
        core: &mut BcfCore,
        str: &mut Vec<u8>,
        indiv: &mut Vec<u8>,
    ) -> Option<()> {
        let l = line.len();
        let mut s = Vec::with_capacity(l + 4);
        s.extend_from_slice(line);
        s.extend_from_slice(&[0, 0, 0, 0]);
        let mut tok = Tok { pos: 0, finished: false };

        // CHROM
        let (p, _) = tok.next(&mut s, b'\t')?;
        let name = lossy(cstr(&s, p));
        let mut rid = h.id2int(BCF_DT_CTG, &name);
        if rid < 0 {
            log::warn!("Contig '{name}' is not defined in the header. (Quick workaround: index the file with tabix.)");
            self.errcode = BCF_ERR_CTG_UNDEF;
            add_dummy(h, &format!("##contig=<ID={name}>"));
            rid = h.id2int(BCF_DT_CTG, &name);
            if rid < 0 {
                log::error!("Could not add dummy header for contig '{name}'");
                self.errcode |= BCF_ERR_CTG_INVALID;
                return None;
            }
        }
        core.rid = rid;
        self.rid = rid;

        // POS
        let (p, _) = tok.next(&mut s, b'\t')?;
        let (pos, used, overflow) = str2uint(&s[p..], 62);
        if overflow {
            log::error!("Position value '{}' is too large", lossy(cstr(&s, p)));
            return None;
        } else if s[p + used] != 0 {
            log::error!("Could not parse the position '{}'", lossy(cstr(&s, p)));
            return None;
        }
        let pos = pos as i64 - 1;
        core.pos = pos;
        self.pos = pos;
        if pos >= i32::MAX as i64 {
            self.is_64bit = true;
        }

        let not_dot = |s: &[u8], p: usize| !(s[p] == b'.' && s[p + 1] == 0);

        // ID
        let (p, q) = tok.next(&mut s, b'\t')?;
        if not_dot(&s, p) {
            enc_vchar(str, &s[p..q]);
        } else {
            enc_size(str, 0, BCF_BT_CHAR);
        }

        // REF
        let (p, q) = tok.next(&mut s, b'\t')?;
        enc_vchar(str, &s[p..q]);
        core.n_allele = 1;
        core.rlen = (q - p) as i64;

        // ALT
        let (p, q) = tok.next(&mut s, b'\t')?;
        if not_dot(&s, p) {
            let mut t = p;
            let mut r = p;
            loop {
                if s[r] == b',' || s[r] == 0 {
                    if core.n_allele == u16::MAX as usize {
                        log::error!("Too many ALT alleles");
                        self.errcode |= BCF_ERR_LIMITS;
                        return None;
                    }
                    enc_vchar(str, &s[t..r]);
                    t = r + 1;
                    core.n_allele += 1;
                }
                if r == q {
                    break;
                }
                r += 1;
            }
        }

        // QUAL
        let (p, _) = tok.next(&mut s, b'\t')?;
        core.qual = if not_dot(&s, p) { (atof(cstr(&s, p)) as f32).to_bits() } else { FLOAT_MISSING };

        // FILTER
        let (p, q) = tok.next(&mut s, b'\t')?;
        if not_dot(&s, p) {
            self.parse_filter(&mut s, p, q, h, str)?;
        } else {
            enc_vint(str, &[], -1);
        }

        // INFO
        let (p, q) = tok.next(&mut s, b'\t')?;
        if not_dot(&s, p) {
            core.n_info = self.parse_info(&mut s, p, q, h, str)?;
        }

        // FORMAT; optional
        if let Some((p, q)) = tok.next(&mut s, b'\t') {
            self.parse_format(&mut s, l, p, q, h, core, indiv)?;
        }
        Some(())
    }

    fn parse_filter(&mut self, s: &mut [u8], p: usize, q: usize, h: &mut HeaderCore, str: &mut Vec<u8>) -> Option<()> {
        if q > 0 && s[q - 1] == b';' {
            s[q - 1] = 0;
        }
        let mut flt = Vec::new();
        let mut tok = Tok { pos: p, finished: false };
        while let Some((t, _)) = tok.next(s, b';') {
            let name = lossy(cstr(s, t));
            let mut id = h.id2int(BCF_DT_ID, &name);
            if id < 0 {
                log::warn!("FILTER '{name}' is not defined in the header");
                let added = add_dummy(h, &format!("##FILTER=<ID={name},Description=\"Dummy\">"));
                id = h.id2int(BCF_DT_ID, &name);
                self.errcode |= BCF_ERR_TAG_UNDEF;
                if !added || id < 0 {
                    log::error!("Could not add dummy header for FILTER '{name}'");
                    self.errcode |= BCF_ERR_TAG_INVALID;
                    return None;
                }
            }
            flt.push(id);
        }
        enc_vint(str, &flt, -1);
        Some(())
    }

    /// Returns the number of INFO fields.
    fn parse_info(&mut self, s: &mut [u8], p: usize, q: usize, h: &mut HeaderCore, str: &mut Vec<u8>) -> Option<usize> {
        let mut n_info = 0;
        if q > 0 && s[q - 1] == b';' {
            s[q - 1] = 0;
        }
        let mut r = p;
        let mut key = p;
        loop {
            while s[r] > b'=' || (s[r] != b';' && s[r] != b'=' && s[r] != 0) {
                r += 1;
            }
            if n_info == u16::MAX as usize {
                log::error!("Too many INFO entries");
                self.errcode |= BCF_ERR_LIMITS;
                return None;
            }
            let mut val = None;
            let end;
            let mut c = s[r];
            s[r] = 0;
            if c == b'=' {
                let vs = r + 1;
                let mut e = vs;
                while s[e] != b';' && s[e] != 0 {
                    e += 1;
                }
                c = s[e];
                s[e] = 0;
                val = Some(vs);
                end = e;
            } else {
                end = r;
            }
            if s[key] == 0 {
                // faulty VCF, ";;" in the INFO
                if c == 0 {
                    break;
                }
                r = end + 1;
                key = r;
                continue;
            }
            let name = lossy(cstr(s, key));
            let mut id = h.id2int(BCF_DT_ID, &name);
            if id < 0 || h.info(BCF_HL_INFO, id) == 15 {
                log::warn!("INFO '{name}' is not defined in the header, assuming Type=String");
                let added = add_dummy(h, &format!("##INFO=<ID={name},Number=1,Type=String,Description=\"Dummy\">"));
                id = h.id2int(BCF_DT_ID, &name);
                self.errcode |= BCF_ERR_TAG_UNDEF;
                if !added || id < 0 {
                    log::error!("Could not add dummy header for INFO '{name}'");
                    self.errcode |= BCF_ERR_TAG_INVALID;
                    return None;
                }
            }
            let y = h.info(BCF_HL_INFO, id);
            n_info += 1;
            enc_int1(str, id);
            let typ = (y >> 4) & 0xf;
            match val {
                None => enc_size(str, 0, BCF_BT_NULL),
                Some(val) if typ == BCF_HT_FLAG || typ == BCF_HT_STR => enc_vchar(str, &s[val..end]),
                Some(val) => {
                    let n_val = 1 + cstr(s, val).iter().filter(|&&c| c == b',').count();
                    if typ == BCF_HT_INT {
                        let mut a = Vec::with_capacity(n_val);
                        let mut t = val;
                        for _ in 0..n_val {
                            let (x, used, overflow) = str2int(&s[t..], 64);
                            let x = if used == 0 {
                                INT32_MISSING
                            } else if overflow || !(MIN_INT32..=MAX_INT32).contains(&x) {
                                log::warn!("Extreme INFO/{name} value encountered and set to missing");
                                INT32_MISSING
                            } else {
                                x as i32
                            };
                            a.push(x);
                            t += used;
                            while s[t] != 0 && s[t] != b',' {
                                t += 1;
                            }
                            t += 1;
                        }
                        if n_val == 1 {
                            enc_int1(str, a[0]);
                        } else {
                            enc_vint(str, &a, -1);
                        }
                    } else if typ == BCF_HT_REAL {
                        let mut a = Vec::with_capacity(n_val);
                        let mut t = val;
                        for _ in 0..n_val {
                            let (x, used, overflow) = str2dbl(&s[t..]);
                            a.push(if used == 0 || overflow { FLOAT_MISSING } else { (x as f32).to_bits() });
                            t += used;
                            while s[t] != 0 && s[t] != b',' {
                                t += 1;
                            }
                            t += 1;
                        }
                        enc_vfloat(str, &a);
                    }
                }
            }
            if c == 0 {
                break;
            }
            r = end + 1;
            key = r;
        }
        Some(n_info)
    }

    /// `vcf_parse_format`; `l` is the length of the line.
    #[allow(clippy::too_many_arguments)]
    fn parse_format(
        &mut self,
        s: &mut [u8],
        l: usize,
        p: usize,
        q: usize,
        h: &mut HeaderCore,
        core: &mut BcfCore,
        indiv: &mut Vec<u8>,
    ) -> Option<()> {
        let nsamples = h.samples.len();
        if nsamples == 0 {
            return Some(());
        }
        let end = l;
        // empty1
        if q >= end {
            log::error!("FORMAT column with no sample columns");
            self.errcode |= BCF_ERR_NCOLS;
            return Some(());
        }
        core.n_fmt = 0;
        if s[p] == b'.' && s[p + 1] == 0 {
            core.n_sample = nsamples as u32;
            return Some(());
        }

        // dict2
        let mut fmt: Vec<FmtAux> = Vec::new();
        let mut tok = Tok { pos: p, finished: false };
        while let Some((t, _)) = tok.next(s, b':') {
            if fmt.len() >= MAX_N_FMT {
                self.errcode |= BCF_ERR_LIMITS;
                log::error!("FORMAT column lists more identifiers than htslib can handle");
                return None;
            }
            let name = lossy(cstr(s, t));
            let mut id = h.id2int(BCF_DT_ID, &name);
            if id < 0 || h.info(BCF_HL_FMT, id) == 15 {
                if name == "." {
                    log::error!("Invalid FORMAT tag name '.'");
                    self.errcode |= BCF_ERR_TAG_INVALID;
                    return None;
                }
                log::warn!("FORMAT '{name}' is not defined in the header, assuming Type=String");
                let added = add_dummy(h, &format!("##FORMAT=<ID={name},Number=1,Type=String,Description=\"Dummy\">"));
                id = h.id2int(BCF_DT_ID, &name);
                self.errcode |= BCF_ERR_TAG_UNDEF;
                if !added || id < 0 {
                    log::error!("Could not add dummy header for FORMAT '{name}'");
                    self.errcode |= BCF_ERR_TAG_INVALID;
                    return None;
                }
            }
            fmt.push(FmtAux { key: id, is_gt: name == "GT", y: h.info(BCF_HL_FMT, id), ..FmtAux::default() });
        }
        let n_fmt = fmt.len();

        // max3
        let mut n_sample = 0usize;
        let mut r = q + 1;
        let (mut m, mut g) = (1usize, 1usize);
        let meta = |c: u8| matches!(c, 0 | b'\t' | b',' | b'/' | b':' | b'|');
        while r < end {
            let mut j = 0;
            let mut r_start = r;
            loop {
                while !meta(s[r]) {
                    r += 1;
                }
                match s[r] {
                    b',' => m += 1,
                    b'|' | b'/' => {
                        if fmt[j].is_gt {
                            g += 1;
                        }
                    }
                    c => {
                        if c == b'\t' {
                            s[r] = 0;
                        }
                        let len = r - r_start;
                        r_start = r;
                        let f = &mut fmt[j];
                        f.max_m = f.max_m.max(m);
                        f.max_l = f.max_l.max(len);
                        if f.is_gt {
                            f.max_g = f.max_g.max(g);
                        }
                        m = 1;
                        g = 1;
                        if c == b':' {
                            j += 1;
                            if j >= n_fmt {
                                log::error!("Incorrect number of FORMAT fields");
                                self.errcode |= BCF_ERR_NCOLS;
                                return None;
                            }
                        } else {
                            break;
                        }
                    }
                }
                if r >= end {
                    break;
                }
                r += 1;
            }
            n_sample += 1;
            if n_sample == nsamples {
                break;
            }
            r += 1;
        }

        // alloc4
        for f in fmt.iter_mut() {
            if f.max_m == 0 {
                f.max_m = 1;
            }
            let htype = (f.y >> 4) & 0xf;
            let size = if htype == BCF_HT_STR {
                if f.is_gt { f.max_g << 2 } else { f.max_l }
            } else if htype == BCF_HT_REAL || htype == BCF_HT_INT {
                f.max_m << 2
            } else {
                log::error!("The format type {htype} is currently not supported");
                self.errcode |= BCF_ERR_TAG_INVALID;
                return None;
            };
            f.size = Some(size);
            f.buf = vec![0; n_sample * size];
        }
        for i in 1..n_fmt {
            if fmt[i].size.is_none() {
                continue;
            }
            if (0..i).any(|j| fmt[j].size.is_some() && fmt[j].key == fmt[i].key) {
                log::warn!("Duplicate FORMAT tag {}", h.id_key(BCF_DT_ID, fmt[i].key).unwrap_or(""));
                self.errcode |= BCF_ERR_TAG_INVALID;
                fmt[i].size = None;
                fmt[i].buf = Vec::new();
            }
        }

        // fill5
        let ver = h.version();
        let mut t = q + 1;
        let mut ms = 0usize;
        while t < end {
            if ms == nsamples {
                break;
            }
            let mut j = 0;
            while t < end {
                let z = &mut fmt[j];
                j += 1;
                let htype = (z.y >> 4) & 0xf;
                match z.size {
                    None => {
                        while s[t] != b':' && s[t] != 0 {
                            t += 1;
                        }
                    }
                    Some(size) if htype == BCF_HT_STR && z.is_gt => {
                        let x = &mut z.buf[size * ms..size * (ms + 1)];
                        let mut is_phased = 0u32;
                        let mut unreadable = false;
                        let mut max = 0u32;
                        let mut overflow = false;
                        let (mut ploidy, mut anyunphased, mut phasingprfx, mut unknown1) = (0, false, false, false);
                        let put =
                            |x: &mut [u8], l: usize, v: u32| x[4 * l..4 * l + 4].copy_from_slice(&v.to_le_bytes());
                        let get = |x: &[u8], l: usize| u32::from_le_bytes(x[4 * l..4 * l + 4].try_into().unwrap());
                        if ver >= VCF44 && (s[t] == b'|' || s[t] == b'/') {
                            is_phased = (s[t] == b'|') as u32;
                            t += 1;
                            phasingprfx = true;
                        }
                        let mut l = 0;
                        loop {
                            ploidy += 1;
                            if s[t] == b'.' {
                                t += 1;
                                put(x, l, is_phased);
                                l += 1;
                                if l == 1 {
                                    unknown1 = true;
                                }
                            } else {
                                let val;
                                if s[t].is_ascii_digit() && !s[t + 1].is_ascii_digit() {
                                    val = (s[t] - b'0') as u32;
                                    t += 1;
                                } else {
                                    let (v, used, of) = str2uint(&s[t..], 4 * 127 - 2);
                                    overflow |= of;
                                    unreadable |= used == 0;
                                    t += used;
                                    val = v as u32;
                                }
                                max = max.max(val);
                                put(x, l, (val.wrapping_add(1) << 1) | is_phased);
                                l += 1;
                            }
                            anyunphased |= ploidy != 1 && is_phased == 0;
                            is_phased = (s[t] == b'|') as u32;
                            if s[t] != b'|' && s[t] != b'/' {
                                break;
                            }
                            t += 1;
                        }
                        if !phasingprfx {
                            if ploidy == 1 {
                                if !unknown1 {
                                    let v = get(x, 0) | 1;
                                    put(x, 0, v);
                                }
                            } else if !anyunphased {
                                let v = get(x, 0) | 1;
                                put(x, 0, v);
                            }
                        }
                        if overflow || max > (i32::MAX >> 1) as u32 - 1 {
                            log::error!("Couldn't read GT data: value too large");
                            return None;
                        }
                        if unreadable {
                            log::error!("Couldn't read GT data: value not a number or '.'");
                            return None;
                        }
                        while l < size >> 2 {
                            put(x, l, INT32_VECTOR_END as u32);
                            l += 1;
                        }
                    }
                    Some(size) if htype == BCF_HT_STR => {
                        let x = &mut z.buf[size * ms..size * (ms + 1)];
                        let mut l = 0;
                        while s[t] != b':' && s[t] != 0 {
                            x[l] = s[t];
                            l += 1;
                            t += 1;
                        }
                        x[l..].fill(0);
                    }
                    Some(size) if htype == BCF_HT_INT => {
                        let x = &mut z.buf[size * ms..size * (ms + 1)];
                        let mut l = 0;
                        loop {
                            let v = if s[t] == b'.' {
                                t += 1;
                                INT32_MISSING
                            } else {
                                let (v, used, overflow) = str2int(&s[t..], 64);
                                t += used;
                                if used == 0 || overflow || !(MIN_INT32..=MAX_INT32).contains(&v) {
                                    log::warn!("Extreme FORMAT value encountered and set to missing");
                                    INT32_MISSING
                                } else {
                                    v as i32
                                }
                            };
                            x[4 * l..4 * l + 4].copy_from_slice(&v.to_le_bytes());
                            l += 1;
                            if s[t] != b',' {
                                break;
                            }
                            t += 1;
                        }
                        while l < size >> 2 {
                            x[4 * l..4 * l + 4].copy_from_slice(&INT32_VECTOR_END.to_le_bytes());
                            l += 1;
                        }
                    }
                    Some(size) if htype == BCF_HT_REAL => {
                        let x = &mut z.buf[size * ms..size * (ms + 1)];
                        let mut l = 0;
                        loop {
                            let v = if s[t] == b'.' && !s[t + 1].is_ascii_digit() {
                                t += 1;
                                FLOAT_MISSING
                            } else {
                                let (v, used, overflow) = str2dbl(&s[t..]);
                                if used == 0 || overflow {
                                    log::warn!("Extreme FORMAT value encountered");
                                }
                                t += used;
                                (v as f32).to_bits()
                            };
                            x[4 * l..4 * l + 4].copy_from_slice(&v.to_le_bytes());
                            l += 1;
                            if s[t] != b',' {
                                break;
                            }
                            t += 1;
                        }
                        while l < size >> 2 {
                            x[4 * l..4 * l + 4].copy_from_slice(&FLOAT_VECTOR_END.to_le_bytes());
                            l += 1;
                        }
                    }
                    Some(_) => {
                        log::error!("Unknown FORMAT field type {htype}");
                        self.errcode |= BCF_ERR_TAG_INVALID;
                        return None;
                    }
                }
                if s[t] == 0 {
                    break;
                } else if s[t] == b':' {
                    t += 1;
                } else {
                    log::error!("Invalid character '{}' in FORMAT field", s[t] as char);
                    self.errcode |= BCF_ERR_CHAR;
                    return None;
                }
            }
            // fill end-of-vector values
            for z in fmt[j..].iter_mut() {
                let Some(size) = z.size else { continue };
                let x = &mut z.buf[size * ms..size * (ms + 1)];
                let htype = (z.y >> 4) & 0xf;
                let (first, rest): (u32, u32) = if htype == BCF_HT_STR && !z.is_gt {
                    if size > 0 {
                        x[0] = b'.';
                        x[1..].fill(0);
                    }
                    continue;
                } else if htype == BCF_HT_REAL {
                    (FLOAT_MISSING, FLOAT_VECTOR_END)
                } else {
                    (INT32_MISSING as u32, INT32_VECTOR_END as u32)
                };
                if size >= 4 {
                    x[..4].copy_from_slice(&first.to_le_bytes());
                }
                for l in 1..size >> 2 {
                    x[4 * l..4 * l + 4].copy_from_slice(&rest.to_le_bytes());
                }
            }
            ms += 1;
            t += 1;
        }

        // gt6
        if n_sample > 0 {
            for z in &fmt {
                let Some(size) = z.size else { continue };
                enc_int1(indiv, z.key);
                let htype = (z.y >> 4) & 0xf;
                if htype == BCF_HT_STR && !z.is_gt {
                    enc_size(indiv, size, BCF_BT_CHAR);
                    indiv.extend_from_slice(&z.buf);
                } else if htype == BCF_HT_INT || z.is_gt {
                    let a: Vec<i32> = z.buf.as_chunks::<4>().0.iter().map(|c| i32::from_le_bytes(*c)).collect();
                    enc_vint(indiv, &a, (size >> 2) as i32);
                } else {
                    enc_size(indiv, size >> 2, BCF_BT_FLOAT);
                    indiv.extend_from_slice(&z.buf);
                }
            }
        }
        core.n_fmt = fmt.iter().filter(|z| z.size.is_some()).count();
        core.n_sample = n_sample as u32;

        // check7
        if n_sample != nsamples {
            log::error!("Number of columns does not match the number of samples ({n_sample} vs {nsamples})");
            self.errcode |= BCF_ERR_NCOLS;
            return None;
        }
        Some(())
    }

    /// The INFO field with header ID `id` (`bcf_get_info_id`).
    fn info_index(&self, id: i32) -> Option<usize> {
        self.info.iter().position(|i| i.key == id)
    }

    /// `bcf_get_info`
    pub fn get_info(&self, h: &HeaderCore, key: &str) -> Option<&Info> {
        let id = h.id2int(BCF_DT_ID, key);
        if !h.idinfo_exists(BCF_HL_INFO, id) {
            return None;
        }
        self.info_index(id).map(|i| &self.info[i])
    }

    /// `bcf_get_fmt`: index of the FORMAT field (possibly marked removed).
    pub fn get_fmt(&self, h: &HeaderCore, key: &str) -> Option<usize> {
        let id = h.id2int(BCF_DT_ID, key);
        if !h.idinfo_exists(BCF_HL_FMT, id) {
            return None;
        }
        self.fmt.iter().position(|f| f.id == id)
    }

    /// `get_rlen` (htslib 1.23): the longest of REF (`len_ref`), INFO/END, the
    /// largest |SVLEN| of <DEL>/<DUP>/<CNV>/<INV> alleles and, for gVCF
    /// alleles, FORMAT/LEN.
    pub fn get_rlen(&self, h: &HeaderCore, len_ref: i64) -> i64 {
        let pos = self.pos;
        let n_allele = self.alleles.len();
        let mut svlenals = vec![0u8; 1 + n_allele / 8 + 1];
        let mut use_svlen = false;
        let mut gvcf = false;
        for (i, a) in self.alleles.iter().enumerate().skip(1) {
            if a.first() != Some(&b'<') {
                continue;
            }
            if svlen_on_ref_for_vcf_alt(a) {
                svlenals[i >> 3] |= 1 << (i & 7);
                use_svlen = true;
            } else if a == b"<*>" || a == b"<NON_REF>" {
                gvcf = true;
            }
        }
        let endid = h.id2int(BCF_DT_ID, "END");
        let svlenid = if use_svlen { h.id2int(BCF_DT_ID, "SVLEN") } else { -1 };
        let (mut endinfo, mut svleninfo) = (None, None);
        if svlenid >= 0 || endid >= 0 {
            endinfo = self.get_info(h, "END");
            svleninfo = self.get_info(h, "SVLEN");
        }
        let lenid = if gvcf { h.id2int(BCF_DT_ID, "LEN") } else { -1 };
        let lenfmt = if lenid >= 0 { self.get_fmt(h, "LEN").map(|i| &self.fmt[i]) } else { None };

        let mut end = 0;
        if let Some(e) = endinfo.filter(|e| !e.deleted) {
            end = if e.len == 1 { int_or_zero(e.typ, e.data(), 0).unwrap_or(0) } else { 0 };
        }
        let mut len: i64 = 0;
        if let Some(sv) = svleninfo.filter(|i| !i.deleted) {
            let mut i = 0;
            while (i as i32) < sv.len && i + 1 < n_allele {
                // (htslib checks byte i >> 3 for allele i + 1)
                if svlenals[i >> 3] & (1 << ((i + 1) & 7)) != 0 {
                    let Some(tmp) = int_or_zero(sv.typ, sv.data(), i) else {
                        len = 0;
                        break;
                    };
                    len = len.max(tmp.abs());
                }
                i += 1;
            }
        }
        if (svleninfo.is_none() || len == 0) && end != 0 {
            len = if end > pos { end - pos - 1 } else { 0 };
        }
        let end_svlen = pos + len + 1;

        len = 0;
        if let Some(f) = lenfmt.filter(|f| !f.deleted) {
            'samples: for j in 0..self.n_sample as usize {
                for k in 0..f.n as usize {
                    let Some(tmp) = int_or_zero(f.typ, f.sample(j), k) else {
                        len = 0;
                        break 'samples;
                    };
                    len = len.max(tmp);
                }
            }
        }
        if (lenfmt.is_none() || len == 0) && end != 0 {
            len = if end > pos { end - pos } else { 0 };
        }
        let end_fmtlen = pos + len;
        let hpos = end.max(end_svlen).max(end_fmtlen);
        (hpos - pos).max(len_ref)
    }

    fn len_ref(&self) -> i64 {
        self.alleles.first().map_or(0, |a| a.len() as i64)
    }

    /// `bcf_update_info` with `n` values (none: remove the field).
    pub fn update_info(&mut self, h: &HeaderCore, key: &str, values: Option<InfoValues>) -> Result<()> {
        let id = h.id2int(BCF_DT_ID, key);
        if !h.idinfo_exists(BCF_HL_INFO, id) {
            bail!("no INFO/{key} in the header");
        }
        let is_end = key == "END";
        let is_svlen = key == "SVLEN";
        let idx = self.info_index(id);
        let Some(InfoValues::Int(values)) = values.filter(|InfoValues::Int(v)| !v.is_empty()) else {
            if let Some(i) = idx {
                self.shared_dirty = true;
                self.info[i].deleted = true;
            }
            if is_end || is_svlen {
                self.rlen = self.get_rlen(h, self.len_ref());
            }
            return Ok(());
        };
        if is_end && values.len() != 1 {
            self.errcode |= BCF_ERR_TAG_INVALID;
            bail!("END info tag should only have one value");
        }
        let mut str = Vec::new();
        enc_int1(&mut str, id);
        enc_vint(&mut str, values, -1);
        let info = Info::from_encoded(str)?;
        match idx {
            Some(i) => self.info[i] = info,
            None => self.info.push(info),
        }
        self.shared_dirty = true;
        if is_end || is_svlen {
            self.rlen = self.get_rlen(h, self.len_ref());
        }
        Ok(())
    }

    /// `bcf_update_format` (`None`: remove the field).
    pub fn update_format(&mut self, h: &HeaderCore, key: &str, values: Option<FmtValues>) -> Result<()> {
        let id = h.id2int(BCF_DT_ID, key);
        if !h.idinfo_exists(BCF_HL_FMT, id) {
            if values.is_none() {
                return Ok(());
            }
            bail!("no FORMAT/{key} in the header");
        }
        let idx = self.fmt.iter().position(|f| f.id == id);
        let is_len = key == "LEN";
        let Some(values) = values else {
            if let Some(i) = idx {
                self.fmt[i].deleted = true;
            }
            if is_len {
                self.rlen = self.get_rlen(h, self.len_ref());
            }
            return Ok(());
        };
        self.n_sample = h.samples.len() as u32;
        let ns = self.n_sample as usize;
        if ns == 0 {
            bail!("no samples to update FORMAT/{key} for");
        }
        let mut str = Vec::new();
        enc_int1(&mut str, id);
        match values {
            FmtValues::Int(a) => enc_vint(&mut str, a, (a.len() / ns) as i32),
            FmtValues::Float(a) => {
                enc_size(&mut str, a.len() / ns, BCF_BT_FLOAT);
                for &x in a {
                    str.extend_from_slice(&x.to_le_bytes());
                }
            }
            FmtValues::Str(a) => {
                enc_size(&mut str, a.len() / ns, BCF_BT_CHAR);
                str.extend_from_slice(a);
            }
        }
        let fmt = Fmt::decode(&mut Cursor { buf: &str, pos: 0 }, ns)?;
        match idx {
            Some(i) => self.fmt[i] = fmt,
            None if !self.fmt.is_empty() && key == "GT" => self.fmt.insert(0, fmt),
            None => self.fmt.push(fmt),
        }
        if is_len {
            self.rlen = self.get_rlen(h, self.len_ref());
        }
        Ok(())
    }

    /// `vcf_format`: the VCF line of this record.
    pub fn format_vcf(&self, h: &HeaderCore, s: &mut Vec<u8>) -> Result<()> {
        let Some(chrom) = h.id_key(BCF_DT_CTG, self.rid) else {
            bail!("Invalid BCF, CONTIG id={} not present in the header", self.rid);
        };
        s.extend_from_slice(chrom.as_bytes());
        s.push(b'\t');
        put_int(self.pos + 1, s);
        s.push(b'\t');
        s.extend_from_slice(&self.id);
        s.push(b'\t');
        match self.alleles.first() {
            Some(a) => s.extend_from_slice(a),
            None => s.push(b'.'),
        }
        s.push(b'\t');
        if self.alleles.len() > 1 {
            for (i, a) in self.alleles[1..].iter().enumerate() {
                if i > 0 {
                    s.push(b',');
                }
                s.extend_from_slice(a);
            }
        } else {
            s.push(b'.');
        }
        s.push(b'\t');
        if self.qual == FLOAT_MISSING {
            s.push(b'.');
        } else {
            put_double(f32::from_bits(self.qual) as f64, s);
        }
        s.push(b'\t');
        if self.flt.is_empty() {
            s.push(b'.');
        }
        for (i, &idx) in self.flt.iter().enumerate() {
            let Some(key) = h.id_key(BCF_DT_ID, idx) else {
                bail!("Invalid BCF, the FILTER tag id={idx} not present in the header");
            };
            if i > 0 {
                s.push(b';');
            }
            s.extend_from_slice(key.as_bytes());
        }
        s.push(b'\t');
        let mut first = true;
        for z in self.info.iter().filter(|z| !z.deleted) {
            let Some(key) = h.id_key(BCF_DT_ID, z.key) else {
                bail!("Invalid BCF, the INFO tag id={} not present in the header", z.key);
            };
            if !first {
                s.push(b';');
            }
            first = false;
            s.extend_from_slice(key.as_bytes());
            if z.len <= 0 {
                continue;
            }
            s.push(b'=');
            let d = z.data();
            if z.len != 1 {
                fmt_array(s, z.len as usize, z.typ, d);
            } else {
                match z.typ {
                    BCF_BT_FLOAT => {
                        if le_u32(d, 0) == FLOAT_MISSING {
                            s.push(b'.');
                        } else {
                            put_double(f32::from_bits(le_u32(d, 0)) as f64, s);
                        }
                    }
                    BCF_BT_CHAR => s.push(d[0]),
                    BCF_BT_NULL | BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 | BCF_BT_INT64 => {
                        let (v, missing) = match z.typ {
                            BCF_BT_INT8 => (d[0] as i8 as i64, INT8_MISSING as i64),
                            BCF_BT_INT16 => (le_i16(d, 0) as i64, INT16_MISSING as i64),
                            BCF_BT_INT32 => (le_i32(d, 0) as i64, INT32_MISSING as i64),
                            BCF_BT_INT64 => (le_i64(d, 0), INT64_MISSING),
                            // v1 of a NULL value stays 0
                            _ => (0, 0),
                        };
                        if v == missing {
                            s.push(b'.');
                        } else {
                            put_int(v, s);
                        }
                    }
                    t => bail!("Unexpected type {t}"),
                }
            }
        }
        if first {
            s.push(b'.');
        }

        if self.n_sample > 0 {
            if self.fmt.is_empty() {
                for _ in 0..=self.n_sample {
                    s.extend_from_slice(b"\t.");
                }
            } else {
                let mut gt_i = None;
                let mut first = true;
                for (i, z) in self.fmt.iter().enumerate() {
                    if z.deleted {
                        continue;
                    }
                    s.push(if first { b'\t' } else { b':' });
                    first = false;
                    let Some(key) = h.id_key(BCF_DT_ID, z.id) else {
                        bail!("Invalid BCF, the FORMAT tag id={} not present in the header", z.id);
                    };
                    s.extend_from_slice(key.as_bytes());
                    if key == "GT" {
                        gt_i = Some(i);
                    }
                }
                if first {
                    s.extend_from_slice(b"\t.");
                }
                for j in 0..self.n_sample as usize {
                    s.push(b'\t');
                    let mut first = true;
                    for (i, f) in self.fmt.iter().enumerate() {
                        if f.deleted {
                            continue;
                        }
                        if !first {
                            s.push(b':');
                        }
                        first = false;
                        if gt_i == Some(i) {
                            format_gt(h, f, j, s)?;
                        } else if f.n == 1 {
                            fmt_array1(s, f.typ, f.sample(j));
                        } else {
                            fmt_array(s, f.n as usize, f.typ, f.sample(j));
                        }
                    }
                    if first {
                        s.push(b'.');
                    }
                }
            }
        }
        s.push(b'\n');
        Ok(())
    }

    /// The BCF encoding of the record (`bcf1_sync` + `bcf_write`).
    pub fn encode_bcf(&mut self, out: &mut Vec<u8>) -> Result<()> {
        if self.is_64bit {
            bail!("Data contains 64-bit values not representable in BCF. Please use VCF instead");
        }
        let mut shared = Vec::new();
        if self.id == b"." {
            enc_size(&mut shared, 0, BCF_BT_CHAR);
        } else {
            enc_vchar(&mut shared, &self.id);
        }
        for a in &self.alleles {
            enc_vchar(&mut shared, a);
        }
        if self.shared_dirty && self.rlen == 0 && !self.alleles.is_empty() {
            self.rlen = self.alleles[0].len() as i64;
        }
        enc_vint(&mut shared, &self.flt, -1);
        let mut n_info = 0u16;
        for i in self.info.iter().filter(|i| !i.deleted) {
            shared.extend_from_slice(&i.raw);
            n_info += 1;
        }
        let mut indiv = Vec::new();
        let mut n_fmt = 0u32;
        if self.n_sample > 0 {
            for f in self.fmt.iter().filter(|f| !f.deleted) {
                indiv.extend_from_slice(&f.raw);
                n_fmt += 1;
            }
        }
        out.extend_from_slice(&(shared.len() as u32 + 24).to_le_bytes());
        out.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
        out.extend_from_slice(&self.rid.to_le_bytes());
        out.extend_from_slice(&(self.pos as u32).to_le_bytes());
        out.extend_from_slice(&(self.rlen as u32).to_le_bytes());
        out.extend_from_slice(&self.qual.to_le_bytes());
        out.extend_from_slice(&n_info.to_le_bytes());
        out.extend_from_slice(&(self.alleles.len() as u16).to_le_bytes());
        out.extend_from_slice(&(n_fmt << 24 | (self.n_sample & 0xff_ffff)).to_le_bytes());
        out.extend_from_slice(&shared);
        out.extend_from_slice(&indiv);
        Ok(())
    }

    /// `bcf_record_check` on a decoded BCF record: IDs must be in the
    /// header, and before VCFv4.4 the first allele of GT gets its phasing
    /// bit as VCFv4.4 infers it.
    pub fn check_bcf(&mut self, h: &HeaderCore) -> Result<()> {
        let mut err = 0;
        if h.id_key(BCF_DT_CTG, self.rid).is_none() {
            log::warn!("Bad BCF record at {}: Invalid CONTIG id {}", self.pos + 1, self.rid);
            err |= BCF_ERR_CTG_INVALID;
        }
        if self.alleles.is_empty() {
            err |= BCF_ERR_TAG_UNDEF;
        }
        let valid_type = |t: u8, n: usize| {
            matches!(t, BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 | BCF_BT_FLOAT | BCF_BT_CHAR)
                || (t == BCF_BT_NULL && n == 0)
        };
        let known = |id: i32| h.id_key(BCF_DT_ID, id).is_some();
        if !self.flt.iter().all(|&id| known(id)) {
            err |= BCF_ERR_TAG_UNDEF;
        }
        for i in &self.info {
            if !known(i.key) {
                err |= BCF_ERR_TAG_UNDEF;
            }
            if !valid_type(i.typ, i.len as usize) {
                err |= BCF_ERR_TAG_INVALID;
            }
        }
        let idgt = if h.version() < VCF44 { h.id2int(BCF_DT_ID, "GT") } else { -1 };
        let n_sample = self.n_sample as usize;
        for f in self.fmt.iter_mut() {
            if !known(f.id) {
                err |= BCF_ERR_TAG_UNDEF;
            }
            if !valid_type(f.typ, f.n as usize) {
                err |= BCF_ERR_TAG_INVALID;
            }
            if idgt >= 0 && f.id == idgt {
                // updatephasing: on the least significant byte of each value
                let inc = 1usize << TYPE_SHIFT[f.typ as usize];
                let ploidy = f.n as usize;
                let p = f.p_mut();
                for j in 0..n_sample {
                    let b = j * ploidy * inc;
                    match ploidy {
                        0 => {}
                        1 => {
                            if p[b] != 0 {
                                p[b] |= 1;
                            }
                        }
                        _ => {
                            let mut allphased = 1;
                            for k in 1..ploidy {
                                allphased &= p[b + inc * k];
                            }
                            p[b] |= allphased & 1;
                        }
                    }
                }
            }
        }
        if err == 0 && self.rlen < 0 {
            log::warn!("BCF record at {} has invalid RLEN ({})", self.pos + 1, self.rlen);
            self.rlen = self.get_rlen(h, self.len_ref()).max(0);
        }
        self.errcode |= err;
        if err != 0 {
            bail!("Bad BCF record");
        }
        Ok(())
    }
}

/// `bcf_format_gt_v2`
fn format_gt(h: &HeaderCore, f: &Fmt, j: usize, s: &mut Vec<u8>) -> Result<()> {
    let pos = s.len();
    let mut anyunphased = false;
    let mut val0 = 0i64;
    let ploidy = match f.typ {
        BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32 => {
            let vector_end = match f.typ {
                BCF_BT_INT8 => INT8_VECTOR_END as i64,
                BCF_BT_INT16 => INT16_VECTOR_END as i64,
                _ => INT32_VECTOR_END as i64,
            };
            let p = f.sample(j);
            let mut i = 0;
            while i < f.n as usize {
                let val = int_at(f.typ, p, i).expect("integer type");
                if val == vector_end {
                    break;
                }
                if i == 0 {
                    val0 = val;
                } else {
                    s.push(if val & 1 != 0 { b'|' } else { b'/' });
                    anyunphased |= val & 1 == 0;
                }
                if val >> 1 == 0 {
                    s.push(b'.');
                } else {
                    put_int((val >> 1) - 1, s);
                }
                i += 1;
            }
            if i == 0 {
                s.push(b'.');
            }
            i
        }
        BCF_BT_NULL => {
            s.push(b'.');
            1
        }
        t => bail!("Unexpected type {t}"),
    };
    if h.version_raw() >= VCF44 {
        let prefix = if val0 & 1 != 0 {
            ((ploidy > 1 && anyunphased) || (ploidy <= 1 && val0 >> 1 == 0)).then_some(b'|')
        } else {
            ((ploidy <= 1 && val0 != 0) || (ploidy > 1 && !anyunphased)).then_some(b'/')
        };
        if let Some(c) = prefix {
            s.insert(pos, c);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(text: &str) -> HeaderCore {
        let mut h = HeaderCore::new();
        h.parse(text.as_bytes()).unwrap();
        h
    }

    const HDR: &str = "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
        ##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">\n##contig=<ID=chr1,length=1000>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

    fn roundtrip(h: &mut HeaderCore, line: &str) -> String {
        let (r, ret) = Record::parse_vcf(line.as_bytes(), h);
        assert_eq!(ret, 0, "{line}");
        let mut out = Vec::new();
        r.format_vcf(h, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn parse_and_format() {
        let mut h = header(HDR);
        let line = "chr1\t10\trs1\tA\tC,G\t50.5\tPASS\tDP=7;AF=0.5,1e-05\tGT:PL\t0|1:1,2,3,4,5,6\t./.:.";
        assert_eq!(roundtrip(&mut h, line), format!("{line}\n"));
        let line = "chr1\t10\t.\tA\t.\t.\t.\t.\tGT\t0/1\t1";
        assert_eq!(roundtrip(&mut h, line), format!("{line}\n"));
    }

    #[test]
    fn undeclared_tags_are_added() {
        let mut h = header(HDR);
        let out = roundtrip(&mut h, "chr2\t5\t.\tA\tT\t1\tq10\tXX=abc\tGT:YY\t0/1:3\t1/1");
        assert_eq!(out, "chr2\t5\t.\tA\tT\t1\tq10\tXX=abc\tGT:YY\t0/1:3\t1/1:.\n");
        assert!(h.id2int(BCF_DT_CTG, "chr2") >= 0);
    }

    #[test]
    fn bcf_roundtrip() {
        let mut h = header(HDR);
        let (mut r, _) = Record::parse_vcf(b"chr1\t10\trs1\tA\tC\t50\tPASS\tDP=300\tGT:PL\t0|1:1,2,3\t1/1:.", &mut h);
        let mut out = Vec::new();
        r.encode_bcf(&mut out).unwrap();
        let shared_len = u32::from_le_bytes(out[0..4].try_into().unwrap()) as usize - 24;
        let core = BcfCore { rid: 0, pos: 9, rlen: 1, qual: r.qual, n_info: 1, n_allele: 2, n_sample: 2, n_fmt: 2 };
        let r2 = Record::from_bcf(&core, &out[32..32 + shared_len], &out[32 + shared_len..]).unwrap();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        r.format_vcf(&h, &mut a).unwrap();
        r2.format_vcf(&h, &mut b).unwrap();
        assert_eq!(a, b);
    }
}

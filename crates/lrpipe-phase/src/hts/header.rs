//! The VCF/BCF header as htslib models it (`bcf_hdr_t`): header records in
//! order, the dictionaries of IDs (FILTER/INFO/FORMAT), contigs and samples
//! with their numeric IDs, and the dictionary of generic and structured
//! records used to detect duplicates.

use std::collections::HashMap;

use anyhow::{Result, bail};

pub const BCF_HL_FLT: i32 = 0;
pub const BCF_HL_INFO: i32 = 1;
pub const BCF_HL_FMT: i32 = 2;
pub const BCF_HL_CTG: i32 = 3;
pub const BCF_HL_STR: i32 = 4;
pub const BCF_HL_GEN: i32 = 5;

pub const BCF_DT_ID: usize = 0;
pub const BCF_DT_CTG: usize = 1;
pub const BCF_DT_SAMPLE: usize = 2;

pub const BCF_HT_FLAG: u32 = 0;
pub const BCF_HT_INT: u32 = 1;
pub const BCF_HT_REAL: u32 = 2;
pub const BCF_HT_STR: u32 = 3;

pub const BCF_VL_FIXED: u32 = 0;
pub const BCF_VL_VAR: u32 = 1;
pub const BCF_VL_A: u32 = 2;
pub const BCF_VL_G: u32 = 3;
pub const BCF_VL_R: u32 = 4;
const BCF_VL_P: u32 = 5;
const BCF_VL_LA: u32 = 6;
const BCF_VL_LR: u32 = 7;
const BCF_VL_LG: u32 = 8;
const BCF_VL_M: u32 = 9;

pub const VCF_DEF: i32 = 4_002_000;
pub const VCF44: i32 = 4_004_000;

/// A header record (`bcf_hrec_t`).
#[derive(Clone, Debug)]
pub struct Hrec {
    pub typ: i32,
    pub key: String,
    /// Value of a generic line (`##key=value`)
    pub value: Option<String>,
    pub keys: Vec<String>,
    pub vals: Vec<Option<String>>,
    uid: u64,
}

impl Hrec {
    pub fn new(key: &str) -> Hrec {
        Hrec { typ: -1, key: key.to_string(), value: None, keys: Vec::new(), vals: Vec::new(), uid: 0 }
    }

    /// `bcf_hrec_add_key`
    pub fn add_key(&mut self, key: &str) {
        self.keys.push(key.to_string());
        self.vals.push(None);
    }

    /// `bcf_hrec_set_val`
    pub fn set_val(&mut self, i: usize, val: &str, quoted: bool) {
        self.vals[i] = Some(if quoted { format!("\"{val}\"") } else { val.to_string() });
    }

    /// `bcf_hrec_find_key` (case-insensitive)
    pub fn find_key(&self, key: &str) -> Option<usize> {
        self.keys.iter().position(|k| k.eq_ignore_ascii_case(key))
    }

    fn val(&self, i: usize) -> &str {
        self.vals[i].as_deref().unwrap_or("")
    }

    fn set_type(&mut self) {
        self.typ = match self.key.as_str() {
            "contig" => BCF_HL_CTG,
            "INFO" => BCF_HL_INFO,
            "FILTER" => BCF_HL_FLT,
            "FORMAT" => BCF_HL_FMT,
            _ if !self.keys.is_empty() => BCF_HL_STR,
            _ => BCF_HL_GEN,
        };
    }

    /// `_bcf_hrec_format`
    fn format(&self, is_bcf: bool, out: &mut String) {
        match &self.value {
            None => {
                out.push_str("##");
                out.push_str(&self.key);
                out.push_str("=<");
                let mut n = 0;
                for (k, v) in self.keys.iter().zip(&self.vals) {
                    if !is_bcf && k == "IDX" {
                        continue;
                    }
                    if n > 0 {
                        out.push(',');
                    }
                    out.push_str(k);
                    out.push('=');
                    out.push_str(v.as_deref().unwrap_or("(null)"));
                    n += 1;
                }
                out.push_str(">\n");
            }
            Some(v) => {
                out.push_str("##");
                out.push_str(&self.key);
                out.push('=');
                out.push_str(v);
                out.push('\n');
            }
        }
    }
}

/// `bcf_idinfo_t`
#[derive(Clone, Copy, Debug)]
pub struct IdInfo {
    pub info: [u32; 3],
    hrec: [Option<u64>; 3],
    pub id: i32,
}

const IDINFO_DEF: IdInfo = IdInfo { info: [15, 15, 15], hrec: [None, None, None], id: -1 };

#[derive(Clone, Debug, Default)]
pub struct HeaderCore {
    pub hrecs: Vec<Hrec>,
    dicts: [HashMap<String, IdInfo>; 3],
    /// numeric ID -> key, per dictionary (`h->id[i]`, length `h->n[i]`)
    ids: [Vec<Option<String>>; 3],
    /// generic and ID-carrying structured records (`aux->gen`)
    gen_keys: HashMap<String, u64>,
    pub samples: Vec<String>,
    /// cached version (`aux->version`, 0 if unknown)
    version: i32,
    next_uid: u64,
    pub dirty: bool,
}

const VALID_CTG: &[u8] = b"!#$%&*+./:;=?@^_|~-";

fn valid_ctg(c: u8) -> bool {
    c.is_ascii_alphanumeric() || VALID_CTG.contains(&c)
}

fn valid_tag(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// `bcf_get_version` for a version string.
pub fn version_from_str(version: &str) -> i32 {
    let fail = || {
        log::warn!("Couldn't get VCF version, considering as 4.2");
        VCF_DEF
    };
    let Some(p) = version.find("VCFv") else { return fail() };
    let major = &version[p + 4..];
    let Some(dot) = major.find('.') else { return fail() };
    let strtol = |s: &str| -> i64 {
        let s = s.trim_start();
        let (neg, digits) = match s.as_bytes().first() {
            Some(b'-') => (true, &s[1..]),
            Some(b'+') => (false, &s[1..]),
            _ => (false, s),
        };
        let n: i64 = digits.bytes().take_while(u8::is_ascii_digit).fold(0, |a, d| a * 10 + (d - b'0') as i64);
        if neg { -n } else { n }
    };
    let ma = strtol(major);
    let mi = strtol(&major[dot + 1..]);
    (ma * 100 * 10000 + mi * 1000) as i32
}

impl HeaderCore {
    pub fn new() -> HeaderCore {
        HeaderCore::default()
    }

    pub fn n(&self, dict: usize) -> usize {
        self.ids[dict].len()
    }

    /// Keys by numeric ID in dictionary `dict` (`h->id[dict]`).
    pub fn ids(&self, dict: usize) -> &[Option<String>] {
        &self.ids[dict]
    }

    /// The version cached when `##fileformat` was added (`aux->version`, 0
    /// if none was).
    pub fn version_raw(&self) -> i32 {
        self.version
    }

    /// Key of numeric ID `id` in dictionary `dict`.
    pub fn id_key(&self, dict: usize, id: i32) -> Option<&str> {
        usize::try_from(id).ok().and_then(|i| self.ids[dict].get(i)).and_then(|k| k.as_deref())
    }

    /// `bcf_hdr_id2int`
    pub fn id2int(&self, dict: usize, key: &str) -> i32 {
        self.dicts[dict].get(key).map_or(-1, |v| v.id)
    }

    fn idinfo_of(&self, id: i32) -> Option<&IdInfo> {
        let key = self.id_key(BCF_DT_ID, id)?;
        self.dicts[BCF_DT_ID].get(key)
    }

    /// The info word of `id` for header line type `hl`.
    pub fn info(&self, hl: i32, id: i32) -> u32 {
        self.idinfo_of(id).map_or(15, |v| v.info[hl as usize])
    }

    /// `bcf_hdr_idinfo_exists`
    pub fn idinfo_exists(&self, hl: i32, id: i32) -> bool {
        id >= 0 && self.info(hl, id) & 0xf != 0xf
    }

    pub fn id2type(&self, hl: i32, id: i32) -> u32 {
        (self.info(hl, id) >> 4) & 0xf
    }

    pub fn id2length(&self, hl: i32, id: i32) -> u32 {
        (self.info(hl, id) >> 8) & 0xf
    }

    pub fn id2number(&self, hl: i32, id: i32) -> u32 {
        self.info(hl, id) >> 12
    }

    /// Contig length (`info[0]` of the contig dictionary).
    pub fn contig_length(&self, rid: i32) -> u64 {
        self.id_key(BCF_DT_CTG, rid).and_then(|k| self.dicts[BCF_DT_CTG].get(k)).map_or(0, |v| v.info[0] as u64)
    }

    pub fn contig_has_hrec(&self, name: &str) -> bool {
        self.dicts[BCF_DT_CTG].get(name).is_some_and(|v| v.hrec[0].is_some())
    }

    /// `bcf_hdr_get_version`-based `bcf_get_version(hdr, NULL)`
    pub fn version(&self) -> i32 {
        if self.version != 0 {
            return self.version;
        }
        match self.hrecs.iter().find(|h| h.typ == BCF_HL_GEN && h.key == "fileformat") {
            Some(h) => version_from_str(h.value.as_deref().unwrap_or("")),
            None => {
                log::warn!("No version string found, assuming VCFv4.2");
                version_from_str("VCFv4.2")
            }
        }
    }

    /// `bcf_hdr_parse_line`: the record and the number of bytes consumed (0
    /// if `line` does not start with "##").
    pub fn parse_line(line: &[u8]) -> (Option<Hrec>, usize) {
        let at = |i: usize| line.get(i).copied().unwrap_or(0);
        if at(0) != b'#' || at(1) != b'#' {
            return (None, 0);
        }
        let mut p = 2;
        let mut q = p;
        while at(q) != 0 && at(q) != b'=' && at(q) != b'\n' {
            q += 1;
        }
        let line_end = |mut q: usize| {
            while at(q) != 0 && at(q) != b'\n' {
                q += 1;
            }
            q + if at(q) != 0 { 1 } else { 0 }
        };
        let malformed = |q: usize| {
            log::error!(
                "Could not parse the header line: \"{}\"",
                String::from_utf8_lossy(&line[..line_end(q).min(line.len())]).trim_end()
            );
            (None, line_end(q))
        };
        if at(q) != b'=' || q == p {
            return malformed(q);
        }
        let s = |a: usize, b: usize| String::from_utf8_lossy(&line[a..b]).into_owned();
        let mut hrec = Hrec::new(&s(p, q));
        q += 1;
        p = q;
        if at(p) != b'<' {
            while at(q) != 0 && at(q) != b'\n' {
                q += 1;
            }
            hrec.value = Some(s(p, q));
            return (Some(hrec), q + if at(q) != 0 { 1 } else { 0 });
        }
        let mut nopen = 1;
        while at(q) != 0 && at(q) != b'\n' && nopen > 0 {
            q += 1;
            p = q;
            while at(q) == b' ' {
                p += 1;
                q += 1;
            }
            if p == q && at(q) != 0 && (at(q).is_ascii_alphabetic() || at(q) == b'_') {
                q += 1;
                while at(q).is_ascii_alphanumeric() || at(q) == b'_' || at(q) == b'.' {
                    q += 1;
                }
            }
            let n = q - p;
            let mut m = 0;
            while at(q) == b' ' {
                q += 1;
                m += 1;
            }
            if at(q) != b'=' || n == 0 {
                return malformed(q);
            }
            hrec.add_key(&s(p, q - m));
            q += 1;
            p = q;
            while at(q) == b' ' {
                p += 1;
                q += 1;
            }
            let (mut quoted, ending) = match at(p) {
                b'"' => {
                    p += 1;
                    (true, b'"')
                }
                b'[' => (true, b']'),
                _ => (false, 0),
            };
            if quoted {
                q += 1;
            }
            while at(q) != 0 && at(q) != b'\n' {
                if quoted {
                    if at(q) == ending {
                        // is_escaped
                        let mut nb = 0;
                        let mut k = q;
                        while k > p && line[k - 1] == b'\\' {
                            nb += 1;
                            k -= 1;
                        }
                        if nb % 2 == 0 {
                            break;
                        }
                    }
                } else {
                    if at(q) == b'<' {
                        nopen += 1;
                    }
                    if at(q) == b'>' {
                        nopen -= 1;
                    }
                    if nopen == 0 {
                        break;
                    }
                    if at(q) == b',' && nopen == 1 {
                        break;
                    }
                }
                q += 1;
            }
            let mut r = q;
            if quoted && ending == b']' {
                if at(q) == ending {
                    r += 1;
                    q += 1;
                    quoted = false;
                } else {
                    log::error!("Missing ']' in header line {}", String::from_utf8_lossy(&line[..q]));
                    return (None, usize::MAX);
                }
            }
            while r > p && line[r - 1] == b' ' {
                r -= 1;
            }
            let i = hrec.keys.len() - 1;
            hrec.set_val(i, &s(p, r), quoted);
            if quoted && at(q) == ending {
                q += 1;
            }
            if at(q) == b'>' {
                if nopen > 0 {
                    nopen -= 1;
                }
                q += 1;
            }
        }
        if nopen > 0 {
            log::warn!("Incomplete header line, trying to proceed anyway");
        }
        let mut nonspace = false;
        while at(q) != 0 && at(q) != b'\n' {
            nonspace |= !matches!(at(q), b' ' | b'\t' | b'\r' | 0x0b | 0x0c);
            q += 1;
        }
        if nonspace {
            log::warn!("Dropped trailing junk from header line");
        }
        (Some(hrec), q + if at(q) != 0 { 1 } else { 0 })
    }

    /// `bcf_hrec_check`: warnings only (and sets the type).
    fn hrec_check(hrec: &mut Hrec) {
        hrec.set_type();
        if hrec.typ == BCF_HL_CTG || hrec.typ == BCF_HL_INFO || hrec.typ == BCF_HL_FMT {
            let Some(i) = hrec.find_key("ID") else {
                log::warn!("Missing ID attribute in one or more header lines");
                return;
            };
            let v = hrec.val(i).as_bytes();
            let ok = if hrec.typ == BCF_HL_CTG {
                !v.is_empty() && v[0] != b'*' && v[0] != b'=' && v.iter().all(|&c| valid_ctg(c))
            } else if hrec.typ == BCF_HL_INFO && v == b"1000G" {
                true
            } else {
                !v.is_empty() && v[0] != b'.' && !v[0].is_ascii_digit() && v.iter().all(|&c| valid_tag(c))
            };
            if !ok {
                log::warn!(
                    "Invalid {} name: \"{}\"",
                    if hrec.typ == BCF_HL_CTG { "contig" } else { "tag" },
                    hrec.val(i)
                );
            }
        }
    }

    /// `bcf_hdr_set_idx`
    fn set_idx(&mut self, dict: usize, key: &str, idinfo: &mut IdInfo) -> Result<()> {
        let n = self.ids[dict].len() as i32;
        if idinfo.id == -1 {
            idinfo.id = n;
        } else if idinfo.id < n && self.ids[dict][idinfo.id as usize].is_some() {
            bail!("Conflicting IDX={} lines in the header dictionary, the new tag is {key}", idinfo.id);
        }
        let new_n = if idinfo.id >= n { idinfo.id + 1 } else { n };
        self.ids[dict].resize(new_n as usize, None);
        self.ids[dict][idinfo.id as usize] = Some(key.to_string());
        Ok(())
    }

    fn add_idx(hrec: &mut Hrec, idx: i32) {
        hrec.keys.push("IDX".into());
        hrec.vals.push(Some(idx.to_string()));
    }

    /// `bcf_hdr_register_hrec`: 1 if registered, 0 if not.
    fn register_hrec(&mut self, hrec: &mut Hrec) -> Result<i32> {
        hrec.set_type();
        if hrec.typ == BCF_HL_CTG {
            let len: u64 = match hrec.find_key("length") {
                None => 0,
                Some(i) => {
                    let v = hrec.val(i);
                    let digits = v.trim_start();
                    let (neg, d) = match digits.as_bytes().first() {
                        Some(b'-') => (true, &digits[1..]),
                        Some(b'+') => (false, &digits[1..]),
                        _ => (false, digits),
                    };
                    let nd = d.bytes().take_while(u8::is_ascii_digit).count();
                    if nd == 0 || neg && d[..nd].bytes().any(|c| c != b'0') {
                        return Ok(0);
                    }
                    d[..nd].parse::<u64>().unwrap_or(u64::MAX)
                }
            };
            let Some(i) = hrec.find_key("ID") else { return Ok(0) };
            let name = hrec.val(i).to_string();
            let replacing = match self.dicts[BCF_DT_CTG].get(&name) {
                Some(v) if v.hrec[0].is_some() => return Ok(0),
                Some(_) => true,
                None => false,
            };
            let mut idx = -1;
            if let Some(j) = hrec.find_key("IDX") {
                match hrec.val(j).parse::<i64>() {
                    Ok(v) if (0..(i32::MAX - 1) as i64).contains(&v) => idx = v as i32,
                    _ => {
                        log::warn!("Error parsing the IDX tag, skipping");
                        return Ok(0);
                    }
                }
            }
            let mut info = IDINFO_DEF;
            info.id = idx;
            info.info[0] = len as u32;
            info.hrec[0] = Some(hrec.uid);
            if let Err(e) = self.set_idx(BCF_DT_CTG, &name, &mut info) {
                if !replacing {
                    self.dicts[BCF_DT_CTG].remove(&name);
                }
                return Err(e);
            }
            self.dicts[BCF_DT_CTG].insert(name, info);
            if idx == -1 {
                Self::add_idx(hrec, info.id);
            }
            return Ok(1);
        }
        if hrec.typ == BCF_HL_STR {
            return Ok(1);
        }
        if hrec.typ != BCF_HL_INFO && hrec.typ != BCF_HL_FLT && hrec.typ != BCF_HL_FMT {
            return Ok(0);
        }
        let mut id: Option<String> = None;
        let mut typ = u32::MAX;
        let mut var = u32::MAX;
        let mut num: i32 = -1;
        let mut idx: i32 = -1;
        for (k, v) in hrec.keys.iter().zip(&hrec.vals) {
            let v = v.as_deref().unwrap_or("");
            match k.as_str() {
                "ID" => id = Some(v.to_string()),
                "IDX" => match v.parse::<i64>() {
                    Ok(x) if (0..(i32::MAX - 1) as i64).contains(&x) => idx = x as i32,
                    _ => {
                        log::warn!("Error parsing the IDX tag, skipping");
                        return Ok(0);
                    }
                },
                "Type" => {
                    typ = match v {
                        "Integer" => BCF_HT_INT,
                        "Float" => BCF_HT_REAL,
                        "String" | "Character" => BCF_HT_STR,
                        "Flag" => BCF_HT_FLAG,
                        _ => {
                            log::warn!("The type \"{v}\" is not supported, assuming \"String\"");
                            BCF_HT_STR
                        }
                    }
                }
                "Number" => {
                    let is_fmt = hrec.typ == BCF_HL_FMT;
                    var = match v {
                        "A" => BCF_VL_A,
                        "R" => BCF_VL_R,
                        "G" => BCF_VL_G,
                        "." => BCF_VL_VAR,
                        "P" if is_fmt => BCF_VL_P,
                        "LA" if is_fmt => BCF_VL_LA,
                        "LR" if is_fmt => BCF_VL_LR,
                        "LG" if is_fmt => BCF_VL_LG,
                        "M" if is_fmt => BCF_VL_M,
                        _ => {
                            // sscanf("%d")
                            let t = v.trim_start();
                            let (neg, d) = match t.as_bytes().first() {
                                Some(b'-') => (true, &t[1..]),
                                Some(b'+') => (false, &t[1..]),
                                _ => (false, t),
                            };
                            let nd = d.bytes().take_while(u8::is_ascii_digit).count();
                            if nd > 0 {
                                let x = d[..nd].parse::<i64>().unwrap_or(i64::MAX) as i32;
                                num = if neg { x.wrapping_neg() } else { x };
                                BCF_VL_FIXED
                            } else {
                                var
                            }
                        }
                    };
                    if var != BCF_VL_FIXED {
                        num = 0xfffff;
                    }
                }
                _ => {}
            }
        }
        if hrec.typ == BCF_HL_INFO || hrec.typ == BCF_HL_FMT {
            if typ == u32::MAX {
                log::warn!("A {} field has no Type defined. Assuming String", hrec.key);
                typ = BCF_HT_STR;
            }
            if var == u32::MAX {
                log::warn!("A {} field has no Number defined. Assuming '.'", hrec.key);
                var = BCF_VL_VAR;
            }
            if typ == BCF_HT_FLAG && (var != BCF_VL_FIXED || num != 0) {
                log::warn!(
                    "The definition of Flag \"{}/{}\" is invalid, forcing Number=0",
                    hrec.key,
                    id.as_deref().unwrap_or("")
                );
                var = BCF_VL_FIXED;
                num = 0;
            }
        }
        let info = ((num as u32) & 0xfffff) << 12 | (var & 0xf) << 8 | (typ & 0xf) << 4 | ((hrec.typ as u32) & 0xf);
        let Some(id) = id else { return Ok(0) };
        let t = (info & 0xf) as usize;
        if let Some(existing) = self.dicts[BCF_DT_ID].get_mut(&id) {
            if existing.hrec[t].is_some() {
                return Ok(0);
            }
            existing.info[t] = info;
            existing.hrec[t] = Some(hrec.uid);
            let eid = existing.id;
            if idx == -1 {
                Self::add_idx(hrec, eid);
            }
            return Ok(1);
        }
        let mut v = IDINFO_DEF;
        v.info[t] = info;
        v.hrec[t] = Some(hrec.uid);
        v.id = idx;
        self.set_idx(BCF_DT_ID, &id, &mut v)?;
        self.dicts[BCF_DT_ID].insert(id, v);
        if idx == -1 {
            Self::add_idx(hrec, v.id);
        }
        Ok(1)
    }

    fn gen_key(hrec: &Hrec) -> Option<String> {
        match hrec.typ {
            BCF_HL_GEN => Some(format!("##{}={}", hrec.key, hrec.value.as_deref().unwrap_or("(null)"))),
            BCF_HL_STR => hrec.find_key("ID").map(|i| format!("##{}=<ID={}>", hrec.key, hrec.val(i))),
            _ => None,
        }
    }

    /// `bcf_hdr_add_hrec`: Ok(1) if the header needs a sync, Ok(0) otherwise.
    pub fn add_hrec(&mut self, mut hrec: Hrec) -> Result<i32> {
        self.next_uid += 1;
        hrec.uid = self.next_uid;
        Self::hrec_check(&mut hrec);
        let res = self.register_hrec(&mut hrec)?;
        let mut key: Option<String> = None;
        if res == 0 {
            if hrec.typ != BCF_HL_GEN {
                return Ok(0);
            }
            let k = Self::gen_key(&hrec).expect("generic record");
            if self.gen_keys.contains_key(&k) {
                return Ok(0);
            }
            if hrec.key == "fileformat" {
                self.version = version_from_str(hrec.value.as_deref().unwrap_or(""));
            }
            key = Some(k);
        }
        if hrec.typ == BCF_HL_STR
            && let Some(k) = Self::gen_key(&hrec)
        {
            if self.gen_keys.contains_key(&k) {
                return Ok(0);
            }
            key = Some(k);
        }
        if let Some(k) = key {
            self.gen_keys.insert(k, hrec.uid);
        }
        let ret = if hrec.typ == BCF_HL_GEN { 0 } else { 1 };
        self.hrecs.push(hrec);
        self.dirty = true;
        Ok(ret)
    }

    /// `bcf_hdr_append`
    pub fn append(&mut self, line: &str) -> Result<()> {
        let (hrec, _) = Self::parse_line(line.as_bytes());
        let Some(hrec) = hrec else { bail!("Invalid header line") };
        self.add_hrec(hrec)?;
        Ok(())
    }

    fn hrec_index(&self, uid: u64) -> Option<usize> {
        self.hrecs.iter().position(|h| h.uid == uid)
    }

    /// `bcf_hdr_get_hrec` (index into `hrecs`)
    pub fn get_hrec(&self, typ: i32, key: &str, value: Option<&str>, str_class: Option<&str>) -> Option<usize> {
        match typ {
            BCF_HL_GEN => match value {
                Some(v) => self.gen_keys.get(&format!("##{key}={v}")).and_then(|&u| self.hrec_index(u)),
                None => self.hrecs.iter().position(|h| h.typ == typ && h.key == key),
            },
            BCF_HL_STR => {
                let class = str_class?;
                let value = value?;
                if key == "ID" {
                    self.gen_keys.get(&format!("##{class}=<{key}={value}>")).and_then(|&u| self.hrec_index(u))
                } else {
                    self.hrecs.iter().position(|h| {
                        h.typ == typ && h.key == class && h.find_key(key).is_some_and(|j| h.val(j) == value)
                    })
                }
            }
            _ => {
                let d = if typ == BCF_HL_CTG { BCF_DT_CTG } else { BCF_DT_ID };
                let slot = if typ == BCF_HL_CTG { 0 } else { typ as usize };
                self.dicts[d].get(value?).and_then(|v| v.hrec[slot]).and_then(|u| self.hrec_index(u))
            }
        }
    }

    /// `bcf_hdr_remove` with a NULL key: all records of a type
    pub fn remove_all(&mut self, typ: i32) {
        let mut i = 0;
        while i < self.hrecs.len() {
            if self.hrecs[i].typ != typ {
                i += 1;
                continue;
            }
            // bcf_hdr_unregister_hrec
            if matches!(typ, BCF_HL_FLT | BCF_HL_INFO | BCF_HL_FMT | BCF_HL_CTG)
                && let Some(j) = self.hrecs[i].find_key("ID")
            {
                let id = self.hrecs[i].val(j).to_string();
                let d = if typ == BCF_HL_CTG { BCF_DT_CTG } else { BCF_DT_ID };
                let slot = if typ == BCF_HL_CTG { 0 } else { typ as usize };
                if let Some(v) = self.dicts[d].get_mut(&id) {
                    v.hrec[slot] = None;
                }
            }
            // bcf_hdr_remove_from_hdict
            if let Some(k) = Self::gen_key(&self.hrecs[i])
                && self.gen_keys.get(&k) == Some(&self.hrecs[i].uid)
            {
                self.gen_keys.remove(&k);
            }
            self.hrecs.remove(i);
            self.dirty = true;
        }
    }

    /// `bcf_hdr_remove`
    pub fn remove(&mut self, typ: i32, key: &str) {
        loop {
            let i = if matches!(typ, BCF_HL_FLT | BCF_HL_INFO | BCF_HL_FMT | BCF_HL_CTG) {
                let Some(i) = self.get_hrec(typ, "ID", Some(key), None) else { return };
                let d = if typ == BCF_HL_CTG { BCF_DT_CTG } else { BCF_DT_ID };
                let slot = if typ == BCF_HL_CTG { 0 } else { typ as usize };
                if let Some(v) = self.dicts[d].get_mut(key) {
                    v.hrec[slot] = None;
                }
                i
            } else {
                let found = self.hrecs.iter().position(|h| {
                    h.typ == typ
                        && if typ == BCF_HL_GEN {
                            h.key == key
                        } else {
                            h.find_key("ID").is_some_and(|j| h.val(j) == key)
                        }
                });
                let Some(i) = found else { return };
                if let Some(k) = Self::gen_key(&self.hrecs[i])
                    && self.gen_keys.get(&k) == Some(&self.hrecs[i].uid)
                {
                    self.gen_keys.remove(&k);
                }
                i
            };
            self.hrecs.remove(i);
            self.dirty = true;
        }
    }

    /// `bcf_hdr_add_sample_len`
    fn add_sample(&mut self, name: &str) -> Result<()> {
        if name.bytes().all(|c| matches!(c, b' ' | b'\t' | b'\r' | 0x0b | 0x0c | b'\n')) {
            bail!("Empty sample name: trailing spaces/tabs in the header line?");
        }
        if self.dicts[BCF_DT_SAMPLE].contains_key(name) {
            bail!("Duplicated sample name '{name}'");
        }
        let mut v = IDINFO_DEF;
        v.id = self.samples.len() as i32;
        self.dicts[BCF_DT_SAMPLE].insert(name.to_string(), v);
        self.samples.push(name.to_string());
        self.dirty = true;
        Ok(())
    }

    /// `bcf_hdr_parse_sample_line`
    fn parse_sample_line(&mut self, line: &[u8]) -> Result<()> {
        const MANDATORY: &[u8] = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        if !line.starts_with(MANDATORY) {
            bail!(
                "Could not parse the \"#CHROM..\" line, either the fields are incorrect or spaces are present instead of tabs"
            );
        }
        let rest = &line[MANDATORY.len()..];
        if rest.is_empty() || rest[0] == b'\n' || rest[0] == 0 {
            return Ok(());
        }
        if !rest.starts_with(b"\tFORMAT\t") {
            bail!(
                "Could not parse the \"#CHROM..\" line, either FORMAT is missing or spaces are present instead of tabs"
            );
        }
        let mut beg = 8;
        while beg < rest.len() && rest[beg] != 0 {
            let mut end = beg;
            while end < rest.len() && rest[end] != 0 && rest[end] != b'\t' && rest[end] != b'\n' {
                end += 1;
            }
            self.add_sample(&String::from_utf8_lossy(&rest[beg..end]))?;
            if end >= rest.len() || rest[end] == 0 || rest[end] == b'\n' {
                break;
            }
            beg = end + 1;
        }
        Ok(())
    }

    /// `bcf_hdr_sync`
    pub fn sync(&mut self) {
        self.ids[BCF_DT_SAMPLE] = self.samples.iter().cloned().map(Some).collect();
        self.dirty = false;
    }

    /// `bcf_hdr_parse`
    pub fn parse(&mut self, text: &[u8]) -> Result<()> {
        let (first, _) = Self::parse_line(text);
        if first.as_ref().is_none_or(|h| !h.key.eq_ignore_ascii_case("fileformat")) {
            log::warn!("The first line should be ##fileformat; is the VCF/BCF header broken?");
        }
        if let Some(h) = first {
            self.add_hrec(h)?;
        }
        let (pass, _) = Self::parse_line(b"##FILTER=<ID=PASS,Description=\"All filters passed\">");
        self.add_hrec(pass.expect("valid line"))?;
        let mut p = 0;
        loop {
            let (hrec, len) = Self::parse_line(&text[p..]);
            if let Some(h) = hrec {
                self.add_hrec(h)?;
                p += len;
                continue;
            }
            if len == usize::MAX {
                bail!("Could not parse header line");
            }
            if len > 0 {
                p += len;
                continue;
            }
            let rest = &text[p..];
            if rest.starts_with(b"#CHROM\t") || rest.starts_with(b"#CHROM ") {
                break;
            }
            match rest.iter().position(|&c| c == b'\n') {
                Some(eol) => {
                    if !rest.is_empty() && rest[0] != 0 {
                        log::warn!("Could not parse header line: \"{}\"", String::from_utf8_lossy(&rest[..eol]));
                    }
                    p += eol + 1;
                }
                None => bail!("Could not parse the header, sample line not found"),
            }
        }
        self.parse_sample_line(&text[p..])?;
        self.sync();
        Ok(())
    }

    /// `bcf_hdr_format`
    pub fn format(&self, is_bcf: bool) -> String {
        let mut out = String::new();
        for h in &self.hrecs {
            h.format(is_bcf, &mut out);
        }
        out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            out.push_str("\tFORMAT");
            for s in &self.samples {
                out.push('\t');
                out.push_str(s);
            }
        }
        out.push('\n');
        out
    }

    /// `bcf_hdr_dup`
    pub fn dup(&self) -> Result<HeaderCore> {
        let mut h = HeaderCore::new();
        h.parse(self.format(true).as_bytes())?;
        Ok(h)
    }

    /// Adds `##contig=<ID=name>` (`add_missing_contig_hrec`).
    pub fn add_missing_contig(&mut self, name: &str) -> Result<()> {
        let mut h = Hrec::new("contig");
        h.add_key("ID");
        h.set_val(0, name, false);
        self.add_hrec(h)?;
        Ok(())
    }

    /// Contig names in ID order (`bcf_hdr_seqnames`: removed ones skipped).
    pub fn seqnames(&self) -> Vec<String> {
        let mut names: Vec<(i32, &str)> = self.dicts[BCF_DT_CTG]
            .iter()
            .filter(|(_, v)| v.hrec[0].is_some())
            .map(|(k, v)| (v.id, k.as_str()))
            .collect();
        names.sort();
        names.into_iter().map(|(_, k)| k.to_string()).collect()
    }
}

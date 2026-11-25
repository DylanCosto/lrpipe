//! Base-pair breakpoints for CNV calls from a structural variant caller's VCF (lrpipe-sv,
//! Sniffles2): a call takes the coordinates of the deletion or duplication that overlaps it
//! best, if they overlap reciprocally by at least half.

use std::io::BufRead;
use std::path::Path;

use crate::Error;
use crate::candidate::{Arena, CandidateId, CnvType};

/// A deletion or duplication from the SV VCF: `[start, end)` 0-based, and its ID.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sv {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub cnv_type: CnvType,
    pub id: String,
}

/// Where a call's breakpoints came from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Breakpoints {
    /// 0-based start and exclusive end.
    pub start: i64,
    pub end: i64,
    /// The SV's ID.
    pub sv_id: String,
}

/// Reads the PASS (or unfiltered) symbolic or sequence-resolved deletions and duplications of
/// at least `min_len` bases from a VCF (plain or gzip).
///
/// # Errors
///
/// If the file cannot be read or a record's POS or END is not a number.
pub fn read(path: &Path, min_len: i64) -> Result<Vec<Sv>, Error> {
    let file = std::fs::File::open(path).map_err(Error::io(path))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let mut out = Vec::new();
    for (n, line) in reader.lines().enumerate() {
        let line = line.map_err(Error::io(path))?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 || !matches!(f[6], "PASS" | ".") {
            continue;
        }
        let info = |key: &str| f[7].split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
        let cnv_type = match info("SVTYPE") {
            Some("DEL") => CnvType::Del,
            Some("DUP") => CnvType::Dup,
            _ => continue,
        };
        let bad = |m: &str| Error::format(path, n + 1, m);
        let pos: i64 = f[1].parse().map_err(|_| bad("bad POS"))?;
        let end: i64 = match (info("END"), info("SVLEN")) {
            (Some(e), _) => e.parse().map_err(|_| bad("bad END"))?,
            (None, Some(l)) => pos + l.parse::<i64>().map_err(|_| bad("bad SVLEN"))?.abs(),
            (None, None) => continue,
        };
        // POS is the base before the event (symbolic alleles); END its last base.
        let (start, end) = (pos, end);
        if end - start >= min_len {
            out.push(Sv { contig: f[0].to_owned(), start, end, cnv_type, id: f[2].to_owned() });
        }
    }
    Ok(out)
}

/// Sets the breakpoints of each call that a same-type SV overlaps reciprocally by at least
/// half (the best such SV). Returns how many calls got breakpoints.
pub fn snap(arena: &mut Arena, calls: &[(String, Vec<CandidateId>)], svs: &[Sv]) -> usize {
    let mut n = 0;
    for &id in calls.iter().flat_map(|(_, ids)| ids) {
        let c = arena.get_mut(id);
        let (cs, ce) = (c.start, c.end + c.bin_size);
        let best = svs
            .iter()
            .filter(|s| s.contig == c.chromosome && s.cnv_type == c.cnv_type)
            .map(|s| {
                let ov = (ce.min(s.end) - cs.max(s.start)).max(0) as f64;
                let r = (ov / (ce - cs) as f64).min(ov / (s.end - s.start).max(1) as f64);
                (r, s)
            })
            .filter(|(r, _)| *r >= 0.5)
            .max_by(|a, b| a.0.total_cmp(&b.0));
        if let Some((_, s)) = best {
            c.breakpoints = Some(Breakpoints { start: s.start, end: s.end, sv_id: s.id.clone() });
            n += 1;
        }
    }
    n
}

//! Scores CNV calls against a truth set.
//!
//! Truth is a VCF (sequence-resolved, as GIAB's: deletions and insertions by allele length, or
//! symbolic DEL/DUP with END/SVLEN) or a BED of `contig start end DEL|DUP`. Calls are a
//! Spectre or lrpipe-cnv VCF, whose POS is a 0-based window start and END the last window's
//! start.
//!
//! A deletion call matches a truth deletion by reciprocal overlap. A duplication call
//! matches a truth duplication by reciprocal overlap, or a truth insertion (how
//! sequence-resolved sets record duplications) of similar size placed within the call or
//! near it.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::Parser;

/// Score CNV calls against a truth set.
#[derive(Parser)]
#[command(version, about)]
struct Args {
    /// Calls: a Spectre / lrpipe-cnv VCF
    calls: PathBuf,
    /// Truth: a VCF, or a BED of `contig start end DEL|DUP`
    truth: PathBuf,
    /// Smallest event scored, in both sets
    #[arg(long, default_value_t = 50_000)]
    min_size: i64,
    /// Window size of the calls (their END is the last window's start)
    #[arg(long, default_value_t = 1000)]
    bin_size: i64,
    /// Reciprocal overlap a match needs
    #[arg(long, default_value_t = 0.5)]
    overlap: f64,
    /// Only these contigs (comma-separated); default all but chrX, chrY and chrM
    #[arg(long)]
    contigs: Option<String>,
    /// Only truth events inside these regions (BED), and calls overlapping them
    #[arg(long)]
    include: Option<PathBuf>,
    /// List unmatched truth events and calls
    #[arg(long)]
    list: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Kind {
    Del,
    Dup,
    Ins,
}

impl Kind {
    fn label(self) -> &'static str {
        match self {
            Kind::Del => "DEL",
            Kind::Dup => "DUP",
            Kind::Ins => "INS",
        }
    }
}

/// An event: `[start, end)` on `contig`; insertions have `start == end` and a `size`.
#[derive(Clone, Debug)]
struct Event {
    contig: String,
    start: i64,
    end: i64,
    size: i64,
    kind: Kind,
}

fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut magic = [0u8; 2];
    let gz = file.read(&mut magic)? == 2 && magic == [0x1f, 0x8b];
    let file = std::fs::File::open(path)?;
    Ok(if gz {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    })
}

fn info<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

/// Events from a VCF: symbolic DEL/DUP (END or SVLEN), or alleles differing in length.
fn read_vcf(path: &Path, calls: bool, bin_size: i64) -> Result<Vec<Event>> {
    let mut out = Vec::new();
    // lrpipe-cnv's standard-style VCF ends events at their last base, as truth sets do.
    let mut standard = false;
    for line in open(path)?.lines() {
        let line = line?;
        if line.starts_with('#') {
            standard |= line == "##source=lrpipe-cnv";
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 || (f.len() > 6 && f[6] != "PASS" && f[6] != "." && !calls) {
            continue;
        }
        let pos: i64 = f[1].parse().context("POS")?;
        let svtype = info(f[7], "SVTYPE").unwrap_or("");
        let kind = match (svtype, f[4]) {
            ("DEL", _) | (_, "<DEL>") => Some(Kind::Del),
            ("DUP", _) | (_, "<DUP>") => Some(Kind::Dup),
            _ => None,
        };
        if let Some(kind) = kind {
            let end: i64 = match info(f[7], "END") {
                Some(e) => e.parse()?,
                None => pos + info(f[7], "SVLEN").map_or(Ok(0), str::parse::<i64>)?.abs(),
            };
            // Calls: 0-based window starts, END the last window's start, unless they took SV
            // breakpoints (SVID). Truth: POS the base before the event, END its last base.
            let (start, end) =
                if calls && !standard && info(f[7], "SVID").is_none() { (pos, end + bin_size) } else { (pos, end) };
            out.push(Event { contig: f[0].into(), start, end, size: end - start, kind });
            continue;
        }
        if calls {
            continue;
        }
        for alt in f[4].split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }
            let (r, a) = (f[3].len() as i64, alt.len() as i64);
            if r > a {
                let start = pos - 1 + 1; // the anchor base precedes the deleted bases
                out.push(Event { contig: f[0].into(), start, end: start + (r - a), size: r - a, kind: Kind::Del });
            } else if a > r {
                out.push(Event { contig: f[0].into(), start: pos, end: pos, size: a - r, kind: Kind::Ins });
            }
        }
    }
    Ok(out)
}

/// Events from a BED of `contig start end DEL|DUP`.
fn read_bed(path: &Path) -> Result<Vec<Event>> {
    let mut out = Vec::new();
    for line in open(path)?.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            bail!("{}: BED lines need contig, start, end and DEL or DUP", path.display());
        }
        let (start, end): (i64, i64) = (f[1].parse()?, f[2].parse()?);
        let kind = match f[3] {
            "DEL" => Kind::Del,
            "DUP" => Kind::Dup,
            other => bail!("{}: unknown type {other}", path.display()),
        };
        out.push(Event { contig: f[0].into(), start, end, size: end - start, kind });
    }
    Ok(out)
}

fn read_regions(path: &Path) -> Result<HashMap<String, Vec<(i64, i64)>>> {
    let mut out: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for line in open(path)?.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 3 && !line.starts_with('#') {
            out.entry(f[0].into()).or_default().push((f[1].parse()?, f[2].parse()?));
        }
    }
    Ok(out)
}

fn reciprocal(a: &Event, b: &Event) -> f64 {
    let ov = (a.end.min(b.end) - a.start.max(b.start)).max(0) as f64;
    (ov / a.size.max(1) as f64).min(ov / b.size.max(1) as f64)
}

/// Whether call `c` and truth `t` describe the same event.
fn matches(c: &Event, t: &Event, overlap: f64) -> bool {
    if c.contig != t.contig {
        return false;
    }
    match (c.kind, t.kind) {
        (Kind::Del, Kind::Del) | (Kind::Dup, Kind::Dup) => reciprocal(c, t) >= overlap,
        (Kind::Dup, Kind::Ins) => {
            // An insertion of the duplicated sequence lies next to (or within) the copy.
            let slack = c.size.max(t.size);
            let near = t.start >= c.start - slack && t.start <= c.end + slack;
            let ratio = c.size.min(t.size) as f64 / c.size.max(t.size) as f64;
            near && ratio >= 0.5
        }
        _ => false,
    }
}

fn band(size: i64) -> &'static str {
    match size {
        s if s < 50_000 => "<50k",
        s if s < 100_000 => "50-100k",
        s if s < 1_000_000 => "100k-1M",
        _ => ">=1M",
    }
}

fn main() -> Result<()> {
    let a = Args::parse();
    let keep_contig = |c: &str| match &a.contigs {
        Some(list) => list.split(',').any(|x| x == c),
        None => !matches!(c, "chrX" | "chrY" | "chrM" | "X" | "Y" | "MT"),
    };
    let include = a.include.as_deref().map(read_regions).transpose()?;
    let inside =
        |e: &Event, whole: bool| {
            include.as_ref().is_none_or(|r| {
                r.get(&e.contig).is_some_and(|v| {
                    v.iter().any(|&(s, t)| {
                        if whole { e.start >= s && e.end.max(e.start + 1) <= t } else { e.start < t && e.end > s }
                    })
                })
            })
        };
    let truth_all = match a.truth.extension().and_then(|e| e.to_str()) {
        Some("bed") => read_bed(&a.truth)?,
        _ => read_vcf(&a.truth, false, a.bin_size)?,
    };
    let truth: Vec<Event> =
        truth_all.into_iter().filter(|e| e.size >= a.min_size && keep_contig(&e.contig) && inside(e, true)).collect();
    let calls: Vec<Event> = read_vcf(&a.calls, true, a.bin_size)?
        .into_iter()
        .filter(|e| e.size >= a.min_size && keep_contig(&e.contig) && inside(e, false))
        .collect();

    let truth_found: Vec<bool> = truth.iter().map(|t| calls.iter().any(|c| matches(c, t, a.overlap))).collect();
    let call_true: Vec<bool> = calls.iter().map(|c| truth.iter().any(|t| matches(c, t, a.overlap))).collect();

    let mut out = String::from("type\tsize\ttruth\tfound\trecall\tcalls\ttrue\tprecision\tmedian_breakpoint_error\n");
    for kind in [Kind::Del, Kind::Dup] {
        for b in ["50-100k", "100k-1M", ">=1M", "all"] {
            let in_band = |e: &Event| b == "all" || band(e.size) == b;
            let truth_kinds =
                |t: &Event| if kind == Kind::Dup { matches!(t.kind, Kind::Dup | Kind::Ins) } else { t.kind == kind };
            let t: Vec<usize> = (0..truth.len()).filter(|&i| truth_kinds(&truth[i]) && in_band(&truth[i])).collect();
            let c: Vec<usize> = (0..calls.len()).filter(|&i| calls[i].kind == kind && in_band(&calls[i])).collect();
            let found = t.iter().filter(|&&i| truth_found[i]).count();
            let good = c.iter().filter(|&&i| call_true[i]).count();
            // Breakpoint error of interval matches: mean of |start| and |end| differences.
            let mut errors: Vec<i64> = c
                .iter()
                .filter_map(|&i| {
                    let call = &calls[i];
                    truth
                        .iter()
                        .filter(|t| t.kind == call.kind && matches(call, t, a.overlap))
                        .map(|t| ((call.start - t.start).abs() + (call.end - t.end).abs()) / 2)
                        .min()
                })
                .collect();
            errors.sort_unstable();
            let ratio =
                |n: usize, d: usize| if d == 0 { "-".to_owned() } else { format!("{:.3}", n as f64 / d as f64) };
            let median = errors.get(errors.len() / 2).map_or("-".to_owned(), |e| e.to_string());
            let _ = writeln!(
                out,
                "{}\t{b}\t{}\t{found}\t{}\t{}\t{good}\t{}\t{median}",
                kind.label(),
                t.len(),
                ratio(found, t.len()),
                c.len(),
                ratio(good, c.len())
            );
        }
    }
    print!("{out}");
    if a.list {
        for (t, found) in truth.iter().zip(&truth_found) {
            if !found {
                println!("missed\t{}\t{}\t{}\t{}\t{}", t.contig, t.start, t.end, t.kind.label(), t.size);
            }
        }
        for (c, ok) in calls.iter().zip(&call_true) {
            if !ok {
                println!("false\t{}\t{}\t{}\t{}\t{}", c.contig, c.start, c.end, c.kind.label(), c.size);
            }
        }
    }
    Ok(())
}

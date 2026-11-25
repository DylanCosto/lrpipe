//! Scores tandem repeat calls against the GIAB HG002 tandem repeat benchmark (v1.0.1) by
//! allele length: each Tier 1 region's maternal and paternal length change against the
//! change of the called alleles inside it.
//!
//! The regions' haplotype length sums are reference length minus haplotype length (positive
//! for contractions, the opposite of what the benchmark's README says; checked against calls
//! where lrpipe-tr and TRGT agree), so they are negated here.
//!
//! Calls are a VCF with `END` in INFO and allele lengths `AL` in FORMAT, POS the base before
//! the repeat (lrpipe-tr, TRGT). Catalog loci (adotto) are the benchmark regions widened by a few
//! bases on each side, so a call belongs to a region when it lies inside it widened by
//! [`PAD`] bases. A region is scored from the calls lying inside it: one
//! call gives its two length changes; with several, at most one may change length (else the
//! region is "ambiguous", as the pairing of alleles across calls is unknown). Regions with no
//! call inside are "not called".
//!
//! `exact_share` and `within_1bp_share` are of the scored regions, which differ between call
//! sets; `call_rate` is the scored share of all regions, and `exact_of_all` counts the regions
//! not called or ambiguous as wrong, so call sets compare on one denominator. With `--common`,
//! only the regions both call sets score are counted.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::Parser;

/// Score tandem repeat calls against the GIAB HG002 tandem repeat benchmark.
#[derive(Parser)]
#[command(version, about)]
struct Args {
    /// Calls: VCF (plain or gzip) with INFO/END and FORMAT/AL
    calls: PathBuf,
    /// Benchmark regions: HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz (or the full bed;
    /// only Tier1 rows are scored)
    regions: PathBuf,
    /// Write each scored region (contig, start, end, truth and called length changes, whether
    /// exact) to this TSV
    #[arg(long)]
    list: Option<PathBuf>,
    /// Count only the regions these other calls score too, so two call sets are compared on
    /// the same regions
    #[arg(long)]
    common: Option<PathBuf>,
}

/// A called locus: `[start, end)` and its alleles' length changes, smallest first.
struct Call {
    start: i64,
    end: i64,
    deltas: [i64; 2],
}

/// Bases a call may extend past a region on each side and still belong to it.
const PAD: i64 = 10;

/// Size classes of the larger truth length change.
const CLASSES: [(&str, i64, i64); 4] = [("0", 0, 0), ("1-4", 1, 4), ("5-49", 5, 49), (">=50", 50, i64::MAX)];

fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// Calls per contig, sorted by start.
fn read_calls(path: &Path) -> Result<HashMap<String, Vec<Call>>> {
    let mut out: HashMap<String, Vec<Call>> = HashMap::new();
    for line in open(path)?.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            bail!("{}: a record has fewer than 10 columns", path.display());
        }
        let start: i64 = f[1].parse().context("POS")?;
        let end: i64 = f[7]
            .split(';')
            .find_map(|kv| kv.strip_prefix("END="))
            .with_context(|| format!("{}: a record lacks INFO/END", path.display()))?
            .parse()
            .context("END")?;
        let al = f[8].split(':').position(|k| k == "AL").and_then(|i| f[9].split(':').nth(i));
        let lens: Vec<i64> = match al {
            Some(v) if v != "." => v.split(',').filter_map(|x| x.parse().ok()).collect(),
            _ => continue,
        };
        let deltas: Vec<i64> = lens.iter().map(|l| l - (end - start)).collect();
        let (a, b) = match deltas[..] {
            [a] => (a, a),
            [a, b] => (a.min(b), a.max(b)),
            _ => continue,
        };
        out.entry(f[0].to_owned()).or_default().push(Call { start, end, deltas: [a, b] });
    }
    for calls in out.values_mut() {
        calls.sort_by_key(|c| c.start);
    }
    Ok(out)
}

#[derive(Default)]
struct Tally {
    regions: usize,
    not_called: usize,
    ambiguous: usize,
    exact: usize,
    within_1: usize,
}

impl Tally {
    fn row(&self, name: &str, out: &mut String) {
        let scored = self.regions - self.not_called - self.ambiguous;
        let of = |n: usize, d: usize| if d == 0 { "-".to_owned() } else { format!("{:.4}", n as f64 / d as f64) };
        let _ = writeln!(
            out,
            "{name}\t{}\t{}\t{}\t{scored}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.regions,
            self.not_called,
            self.ambiguous,
            of(scored, self.regions),
            self.exact,
            of(self.exact, scored),
            of(self.exact, self.regions),
            self.within_1,
            of(self.within_1, scored)
        );
    }
}

/// A region's called length changes: `None` when no call lies inside, `Some(None)` when the
/// calls inside are ambiguous.
fn predict(calls: &HashMap<String, Vec<Call>>, contig: &str, start: i64, end: i64) -> Option<Option<[i64; 2]>> {
    let inside: Vec<&Call> = calls.get(contig).map_or_else(Vec::new, |v| {
        let from = v.partition_point(|c| c.start < start - PAD);
        v[from..].iter().take_while(|c| c.start < end).filter(|c| c.end <= end + PAD).collect()
    });
    let changed: Vec<&&Call> = inside.iter().filter(|c| c.deltas != [0, 0]).collect();
    match (inside.len(), changed.len()) {
        (0, _) => None,
        (_, 0) => Some(Some([0, 0])),
        (_, 1) => Some(Some(changed[0].deltas)),
        _ => Some(None),
    }
}

fn main() -> Result<()> {
    let a = Args::parse();
    let calls = read_calls(&a.calls)?;
    let common = a.common.as_deref().map(read_calls).transpose()?;
    let mut list = a
        .list
        .as_deref()
        .map(|p| {
            std::fs::File::create(p).map(std::io::BufWriter::new).with_context(|| format!("creating {}", p.display()))
        })
        .transpose()?;
    let mut all = Tally::default();
    let mut by_class: Vec<Tally> = CLASSES.iter().map(|_| Tally::default()).collect();
    for line in open(&a.regions)?.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[3] != "Tier1" {
            continue;
        }
        let (start, end): (i64, i64) = (f[1].parse()?, f[2].parse()?);
        let (Ok(mat), Ok(pat)) = (f[7].parse::<i64>(), f[8].parse::<i64>()) else { continue };
        let (mat, pat) = (-mat, -pat);
        let truth = [mat.min(pat), mat.max(pat)];
        let size = mat.abs().max(pat.abs());
        let class = CLASSES.iter().position(|&(_, lo, hi)| size >= lo && size <= hi).expect("classes cover every size");
        let pred = predict(&calls, f[0], start, end);
        if let Some(other) = &common
            && !(matches!(pred, Some(Some(_))) && matches!(predict(other, f[0], start, end), Some(Some(_))))
        {
            continue;
        }
        if let (Some(w), Some(Some(p))) = (list.as_mut(), pred) {
            use std::io::Write as _;
            let [t0, t1] = truth;
            let [p0, p1] = p;
            writeln!(w, "{}\t{start}\t{end}\t{t0}\t{t1}\t{p0}\t{p1}\t{}", f[0], u8::from(p == truth))?;
        }
        for t in [&mut all, &mut by_class[class]] {
            t.regions += 1;
            match pred {
                None => t.not_called += 1,
                Some(None) => t.ambiguous += 1,
                Some(Some(p)) => {
                    t.exact += usize::from(p == truth);
                    t.within_1 += usize::from(p.iter().zip(&truth).all(|(x, y)| (x - y).abs() <= 1));
                }
            }
        }
    }
    let mut out = String::from(
        "size\tregions\tnot_called\tambiguous\tscored\tcall_rate\texact\texact_share\texact_of_all\twithin_1bp\t\
         within_1bp_share\n",
    );
    for (t, (name, ..)) in by_class.iter().zip(CLASSES) {
        t.row(name, &mut out);
    }
    all.row("all", &mut out);
    print!("{out}");
    Ok(())
}

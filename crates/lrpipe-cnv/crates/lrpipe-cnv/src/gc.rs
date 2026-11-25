//! GC content of windows and depth correction for it (an improvement over Spectre).
//!
//! Sequencing depth varies with GC content; uncorrected, GC-rich or GC-poor stretches look
//! like deletions. Each window's depth is scaled by the genome median over the median depth
//! of windows with the same GC percentage.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use pynum::sum;

use crate::Error;
use crate::coverage::ContigDepth;

/// GC fraction of each window, per contig: NaN where a window has no A, C, G or T.
pub type GcFractions = HashMap<String, Vec<f64>>;

/// GC fractions of the windows in `depths`, each measured over `span` bases centred on the
/// window (`bin_size` long from its start), from the reference FASTA. Only the contigs with
/// windows are read, through the FASTA's `.fai` index (or the whole file in order without one).
///
/// # Errors
///
/// If the reference cannot be read.
pub fn fractions(reference: &Path, depths: &[ContigDepth], bin_size: i64, span: i64) -> Result<GcFractions, Error> {
    let fai = reference.with_added_extension("fai");
    if !fai.exists() {
        return fractions_streamed(reference, depths, bin_size, span);
    }
    // Bgzipped or not by its contents, not its name.
    let mut reader = lrpipe_reads::indexed_fasta(reference).map_err(|e| Error::io(reference)(e.into_io_error()))?;
    let mut out = GcFractions::new();
    for d in depths {
        if d.positions.is_empty() {
            continue;
        }
        let in_index = reader.index().as_ref().iter().any(|r| AsRef::<[u8]>::as_ref(r.name()) == d.name.as_bytes());
        if !in_index {
            continue;
        }
        let region = pasta::core::Region::new(d.name.as_str(), ..);
        let record = reader.query(&region).map_err(Error::io(reference))?;
        out.insert(d.name.clone(), window_fractions(record.sequence().as_ref(), &d.positions, bin_size, span));
    }
    Ok(out)
}

/// Per base: whether it counts as A/T, as C/G, or neither (either case).
const BASE_CLASS: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1;
    t[b'T' as usize] = 1;
    t[b'a' as usize] = 1;
    t[b't' as usize] = 1;
    t[b'C' as usize] = 2;
    t[b'G' as usize] = 2;
    t[b'c' as usize] = 2;
    t[b'g' as usize] = 2;
    t
};

/// A/T and C/G bases in `seq`.
fn counts(seq: &[u8]) -> (u32, u32) {
    let (mut at, mut cg) = (0u32, 0u32);
    for &b in seq {
        let c = BASE_CLASS[usize::from(b)];
        at += u32::from(c == 1);
        cg += u32::from(c == 2);
    }
    (at, cg)
}

/// Each window's GC fraction over `[p - flank, p + bin_size + flank)` (clipped to the contig):
/// NaN where it has no A, C, G or T. Counts are taken once per stretch between the windows'
/// edges, which come in order, rather than per base of every window.
fn window_fractions(seq: &[u8], positions: &[i64], bin_size: i64, span: i64) -> Vec<f64> {
    let flank = (span - bin_size).max(0) / 2;
    let clamp = |x: i64| usize::try_from(x.max(0)).unwrap_or(0).min(seq.len());
    // Cumulative counts at every window edge, from one sweep over the edges in order.
    let mut edges: Vec<usize> =
        positions.iter().flat_map(|&p| [clamp(p - flank), clamp(p + bin_size + flank)]).collect();
    edges.sort_unstable();
    edges.dedup();
    let mut cumulative = Vec::with_capacity(edges.len());
    let (mut at, mut cg, mut last) = (0u32, 0u32, 0usize);
    for &e in &edges {
        let (a, c) = counts(&seq[last..e]);
        at += a;
        cg += c;
        last = e;
        cumulative.push((at, cg));
    }
    let at_edge = |x: usize| cumulative[edges.binary_search(&x).expect("every edge was counted")];
    positions
        .iter()
        .map(|&p| {
            let ((a0, c0), (a1, c1)) = (at_edge(clamp(p - flank)), at_edge(clamp(p + bin_size + flank)));
            let (a, c) = (a1 - a0, c1 - c0);
            if a + c == 0 { f64::NAN } else { f64::from(c) / f64::from(a + c) }
        })
        .collect()
}

/// [`fractions`] for a FASTA without an index: the whole file, in order.
fn fractions_streamed(
    reference: &Path,
    depths: &[ContigDepth],
    bin_size: i64,
    span: i64,
) -> Result<GcFractions, Error> {
    let file = std::fs::File::open(reference).map_err(Error::io(reference))?;
    let reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut out = GcFractions::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    let finish = |name: &str, seq: &[u8], out: &mut GcFractions| {
        let Some(d) = depths.iter().find(|d| d.name == name) else { return };
        out.insert(name.to_owned(), window_fractions(seq, &d.positions, bin_size, span));
    };
    for line in reader.split(b'\n') {
        let line = line.map_err(Error::io(reference))?;
        let line = line.strip_suffix(b"\r").unwrap_or(&line);
        if let Some(header) = line.strip_prefix(b">") {
            finish(&name, &seq, &mut out);
            name = String::from_utf8_lossy(header).split_whitespace().next().unwrap_or("").to_owned();
            seq.clear();
        } else {
            seq.extend_from_slice(line);
        }
    }
    finish(&name, &seq, &mut out);
    Ok(out)
}

/// Fewer windows than this in a GC percentage and its depths are not corrected.
const MIN_WINDOWS: usize = 100;

/// Scales each window's depth by `median / median of its GC percentage`, with the medians
/// over non-NaN windows of `reference_contigs` (the autosomes). Windows without a GC value,
/// or in a sparse GC percentage, keep their depth.
pub fn correct(depths: &mut [ContigDepth], gc: &GcFractions, reference_contigs: &[String]) {
    let pct = |g: f64| (g * 100.0).round() as usize;
    let mut by_pct: Vec<Vec<f64>> = vec![Vec::new(); 101];
    let mut all = Vec::new();
    for d in depths.iter().filter(|d| reference_contigs.contains(&d.name)) {
        let Some(g) = gc.get(&d.name) else { continue };
        for (&v, &f) in d.depth.iter().zip(g) {
            if !v.is_nan() && !f.is_nan() {
                by_pct[pct(f)].push(v);
                all.push(v);
            }
        }
    }
    if all.is_empty() {
        return;
    }
    let median = sum::median(&all);
    let factors: Vec<Option<f64>> = by_pct
        .iter()
        .map(|v| {
            let m = if v.len() >= MIN_WINDOWS { sum::median(v) } else { 0.0 };
            (m > 0.0).then(|| median / m)
        })
        .collect();
    for d in depths.iter_mut() {
        let Some(g) = gc.get(&d.name) else { continue };
        for (v, &f) in d.depth.iter_mut().zip(g) {
            if let Some(factor) = (!f.is_nan()).then(|| factors[pct(f)]).flatten() {
                *v *= factor;
            }
        }
    }
}

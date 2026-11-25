//! Multi-sample calling from `.snf` files (Sniffles2's combine mode: `CombineTask`,
//! `resolve_block_groups` and `SVGroup`).
//!
//! Each contig is one task, as Sniffles2 runs with one thread: its blocks are processed in order,
//! candidates of all samples are grouped by position, length and sequence, and each group becomes
//! one call with a genotype per sample.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::call::{Alt, Info, Phase};
use crate::config::Config;
use crate::pymath;
use crate::reads::SvType;
use crate::snf::{self, Block, Candidate, SampleGenotype, SnfBnd, SnfFile};
use crate::vcf::{self, Record, Reference};

/// The first band of the edit distance search (edlib's).
const BAND_START: usize = 64;
/// Candidates are binned by position in bins of this size (`combine_min_size`).
const BIN: i64 = 100;
/// Groups closer than this to the current bin wait for candidates of later bins.
const OVERLAP_ABS: f64 = 2500.0;
/// Bins are resolved once they hold this many candidates (at least; scaled by sample count).
const MIN_BIN_CANDIDATES: usize = 25;

/// One input sample.
#[derive(Debug)]
pub struct Sample {
    pub id: String,
    pub file: SnfFile,
}

/// Opens the inputs: several `.snf` files, or one `.tsv` listing `.snf` files and optional
/// sample IDs. A sample's ID is the listed one, else the SNF's, else the file name stem.
///
/// The files stay open unless `close_handles` (`--combine-close-handles`) or there are too many
/// to stay open under the soft limit on open files; then each block read opens its file again,
/// as Sniffles2 does with `--combine-close-handles`.
///
/// # Errors
///
/// Fails on unreadable lists or SNF files, malformed list lines, duplicate sample IDs, and
/// incompatible contig dictionaries ([`contigs`]).
pub fn open_samples(inputs: &[String], close_handles: bool) -> Result<Vec<Sample>> {
    let listed: Vec<(String, Option<String>)> = if let [list] = inputs
        && extension(list) == "tsv"
    {
        let text = std::fs::read_to_string(list).with_context(|| format!("reading {list}"))?;
        let mut v = Vec::new();
        for (k, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            match line.split('\t').collect::<Vec<_>>()[..] {
                [file] => v.push((file.to_string(), None)),
                [file, id] => v.push((file.to_string(), Some(id.to_string()))),
                _ => bail!(
                    "Invalid sample list .tsv : {list} : Line {} - expected either one or two columns (first \
                     column: .snf filename, second column: optional sample id to overrule the one specified in \
                     the .snf file)",
                    k + 1
                ),
            }
        }
        v
    } else {
        inputs.iter().map(|f| (f.clone(), None)).collect()
    };
    let mut keep_open = !close_handles;
    if keep_open && let Some(limit) = open_file_limit().filter(|&l| listed.len() as u64 + RESERVED_HANDLES > l) {
        eprintln!(
            "{} .snf inputs, with {RESERVED_HANDLES} files kept for the rest of the run, cannot all stay open \
             under the limit of {limit} open files (ulimit -n): each is opened again for every block read, as \
             with --combine-close-handles",
            listed.len()
        );
        keep_open = false;
    }
    let samples: Vec<Sample> = listed
        .into_iter()
        .map(|(path, id)| {
            let file = SnfFile::open(&path, keep_open)?;
            let id = id.or_else(|| file.sample_id.clone()).unwrap_or_else(|| {
                let name = Path::new(&path).file_name().map_or(Cow::Borrowed(""), |n| n.to_string_lossy());
                match name.rfind('.') {
                    Some(dot) if dot > 0 => name[..dot].to_string(),
                    _ => name.into_owned(),
                }
            });
            Ok(Sample { id, file })
        })
        .collect::<Result<_>>()?;
    contigs(&samples)?;
    let mut ids = HashMap::new();
    for sample in &samples {
        if let Some(previous) = ids.insert(&sample.id, &sample.file.path) {
            bail!("duplicate SNF sample ID {:?} in {} and {}", sample.id, previous, sample.file.path);
        }
    }
    Ok(samples)
}

/// The contigs to call across `samples`: every contig any SNF has, in the order they first
/// appear. An SNF lists only the contigs its run read (`--contig`, `--regions`, contigs too short
/// to call), so a panel's SNF and a whole genome's on the same reference combine; a contig a
/// sample lacks has no candidates from it. Sniffles2 takes the last file's list, which leaves
/// out the contigs only the others have.
///
/// # Errors
///
/// Fails when two SNFs give one contig different lengths: they were called on different
/// references.
pub fn contigs(samples: &[Sample]) -> Result<Vec<(String, i64)>> {
    union_of_contigs(samples.iter().map(|s| (s.file.path.as_str(), s.file.contig_lengths.as_slice())))
}

/// [`contigs`] of each file's path and contig list.
fn union_of_contigs<'a>(files: impl Iterator<Item = (&'a str, &'a [(String, i64)])>) -> Result<Vec<(String, i64)>> {
    let mut contigs: Vec<(String, i64)> = Vec::new();
    let mut seen: HashMap<&str, (i64, &str)> = HashMap::new();
    for (path, lengths) in files {
        for (name, length) in lengths {
            match seen.get(name.as_str()) {
                Some(&(other, first)) if other != *length => bail!(
                    "contig {name} is {other} bp in {first} and {length} bp in {path}: the SNF files were called on \
                     different references"
                ),
                Some(_) => {}
                None => {
                    seen.insert(name, (*length, path));
                    contigs.push((name.clone(), *length));
                }
            }
        }
    }
    Ok(contigs)
}

/// Open files left for everything but the kept-open `.snf` inputs: the standard streams, the
/// output, the reference, the thread pool's and the allocator's.
const RESERVED_HANDLES: u64 = 64;

/// The soft limit on open files (`RLIMIT_NOFILE`), or `None` if unlimited or unknown.
fn open_file_limit() -> Option<u64> {
    let mut limit = libc::rlimit { rlim_cur: 0, rlim_max: 0 };
    // SAFETY: `getrlimit` only writes the limits to `limit`, which outlives the call.
    if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &raw mut limit) } != 0 || limit.rlim_cur == libc::RLIM_INFINITY {
        return None;
    }
    Some(limit.rlim_cur)
}

/// The lowercase text after the last `.` (Sniffles2's input type test).
pub fn extension(path: &str) -> String {
    path.rsplit('.').next().unwrap_or_default().to_ascii_lowercase()
}

/// A per-sample genotype of a combined call; `id` lists the sample's merged candidate IDs.
#[derive(Debug, Clone, PartialEq)]
struct SampleCall {
    gt: SampleGenotype,
    id: String,
}

/// A combined call.
#[derive(Debug)]
pub struct CombinedCall {
    contig: String,
    svtype: SvType,
    pos: i64,
    end: i64,
    svlen: i64,
    id: String,
    alt: Alt,
    qual: Option<i64>,
    filter: String,
    info: Vec<(String, Info)>,
    genotypes: Vec<SampleCall>,
    precise: bool,
    support: i64,
    rnames: Vec<String>,
    fwd: i64,
    rev: i64,
    coverage: [Option<i64>; 5],
}

/// Candidates of several samples that describe the same SV (`SVGroup`).
#[derive(Debug)]
struct Group {
    candidates: Vec<Candidate>,
    pos_mean: f64,
    len_mean: f64,
    included: BTreeSet<usize>,
    /// Coverage of samples without a candidate, at the group's position.
    coverage_nonincluded: HashMap<usize, i64>,
    /// The first candidate's mate and orientation, and the mean mate position, for BNDs.
    mate: Option<(SnfBnd, f64)>,
}

impl Group {
    fn new(c: Candidate) -> Self {
        Group {
            pos_mean: c.pos as f64,
            len_mean: c.svlen.abs() as f64,
            included: BTreeSet::from([c.sample]),
            coverage_nonincluded: HashMap::new(),
            mate: c.bnd.as_ref().map(|b| (b.clone(), b.mate_pos as f64)),
            candidates: vec![c],
        }
    }

    /// `add_candidate`: updates the running means exactly as Sniffles2 does.
    fn add(&mut self, c: Candidate) {
        let n = self.candidates.len() as f64;
        self.pos_mean = (self.pos_mean * n + c.pos as f64) / (n + 1.0);
        self.len_mean = (self.len_mean * n + c.svlen.abs() as f64) / (n + 1.0);
        if let (Some((_, mean)), Some(b)) = (&mut self.mate, &c.bnd) {
            *mean = (*mean * n + b.mate_pos as f64) / (n + 1.0);
        }
        self.included.insert(c.sample);
        self.candidates.push(c);
    }

    /// `align_call`: whether the candidate's sequence is similar enough to the group's first.
    fn aligns(&self, c: &Candidate, limit: f64) -> bool {
        if limit == 0.0 {
            return true;
        }
        let similar = |distance: usize| (self.len_mean - distance as f64) / self.len_mean > limit;
        let (a, b) = (self.candidates[0].alt.as_bytes(), c.alt.as_bytes());
        // Any distance above this fails; search with a growing band (as edlib) up to it.
        let max_k = (self.len_mean * (1.0 - limit)).max(0.0).ceil() as usize + 1;
        if a.len().abs_diff(b.len()) > max_k {
            return false;
        }
        let query = EdQuery::new(a);
        let mut k = BAND_START.min(max_k);
        loop {
            if let Some(d) = edit_distance_within(&query, b, k) {
                return similar(d);
            }
            if k >= max_k {
                return false;
            }
            k = (k * 2).min(max_k);
        }
    }
}

/// `resolve_block_groups`: adds candidates, by decreasing support, to their closest group.
fn resolve_groups(svtype: SvType, mut candidates: Vec<Candidate>, mut groups: Vec<Group>, cfg: &Config) -> Vec<Group> {
    // A stable sort, as Python's sorted(reverse=True).
    candidates.sort_by_key(|c| std::cmp::Reverse(c.support));
    let a = &cfg.args;
    for c in candidates {
        let mut best: Option<usize> = None;
        let mut best_dist = f64::INFINITY;
        for (k, g) in groups.iter().enumerate() {
            let separate = a.combine_separate_intra && g.included.contains(&c.sample);
            if svtype == SvType::Bnd {
                let (Some(b), Some((g_bnd, g_mate))) = (&c.bnd, &g.mate) else { continue };
                let dist = (g.pos_mean - c.pos as f64).abs() + (g_mate - b.mate_pos as f64).abs();
                // Sniffles2 compares the mate contig alone, so opposite orientations, different
                // adjacencies, became one allele of the cohort (LRP-016).
                let same_adjacency = g_bnd.mate_contig == b.mate_contig
                    && g_bnd.is_first == b.is_first
                    && g_bnd.is_reverse == b.is_reverse;
                if dist < best_dist && dist <= (a.cluster_merge_bnd * 2) as f64 && same_adjacency && !separate {
                    best = Some(k);
                    best_dist = dist;
                }
            } else {
                let dist = (g.pos_mean - c.pos as f64).abs() + (g.len_mean.abs() - c.svlen.abs() as f64).abs();
                let minlen = g.len_mean.abs().min(c.svlen.abs() as f64);
                if minlen > 0.0
                    && dist < best_dist
                    && dist <= a.combine_match as f64 * minlen.sqrt()
                    && dist <= a.combine_match_max as f64
                    && !separate
                    && g.aligns(&c, a.combine_pctseq)
                {
                    best = Some(k);
                    best_dist = dist;
                }
            }
        }
        match best {
            Some(k) => groups[k].add(c),
            None => groups.push(Group::new(c)),
        }
    }
    groups
}

/// Python's `round(sum / n)` over the values present, or `None`.
fn mean_or_none_round(values: impl Iterator<Item = Option<i64>>) -> Option<i64> {
    let v: Vec<i64> = values.flatten().collect();
    (!v.is_empty()).then(|| pymath::py_round(v.iter().sum::<i64>() as f64 / v.len() as f64))
}

/// `int(statistics.median(values))`.
fn median_int(values: impl Iterator<Item = i64>) -> i64 {
    let mut v: Vec<i64> = values.collect();
    v.sort_unstable();
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { ((v[n / 2 - 1] + v[n / 2]) as f64 / 2.0) as i64 }
}

/// `SVGroup.call`: the group's combined call, or `None` if it is filtered.
fn call_group(g: &Group, cfg: &Config, n_samples: usize, task_id: usize, sv_id: &mut usize) -> Option<CombinedCall> {
    let a = &cfg.args;
    let first = &g.candidates[0];
    let total = g.included.len();
    let passed = g.candidates.iter().filter(|c| c.qc).count();
    let samples = n_samples as f64;
    let qc = (passed > 0 && passed as f64 / samples >= a.combine_high_confidence)
        || (total as f64 / samples >= a.combine_low_confidence && total >= a.combine_low_confidence_abs);
    let single_no_qc = a.no_qc && n_samples == 1;
    if !qc && !single_no_qc {
        return None;
    }
    if !a.combine_output_filtered && !g.candidates.iter().any(|c| c.qc && c.filter == "PASS") && !single_no_qc {
        return None;
    }

    let mut rnames = Vec::new();
    let mut genotypes: BTreeMap<usize, SampleCall> = BTreeMap::new();
    for c in &g.candidates {
        if let Some(r) = &c.rnames {
            rnames.extend(r.iter().cloned());
        }
        let gt = c.genotype.clone().unwrap_or(SampleGenotype {
            a: None,
            b: None,
            gq: 0,
            dr: 0,
            dv: c.support,
            phase: Phase::default(),
        });
        let id = format!("Sniffles2.{}", c.id);
        match genotypes.get_mut(&c.sample) {
            // Several candidates of one sample: keep the highest genotype, list all IDs.
            Some(cur) => {
                let id = format!("{},{id}", cur.id);
                let higher = cur.gt.a.is_none() || (gt.a.is_some() && (gt.a, gt.b) >= (cur.gt.a, cur.gt.b));
                if higher {
                    cur.gt = gt;
                }
                cur.id = id;
            }
            None => {
                genotypes.insert(c.sample, SampleCall { gt, id });
            }
        }
    }
    for s in 0..n_samples {
        if genotypes.contains_key(&s) {
            continue;
        }
        let coverage = g.coverage_nonincluded.get(&s).copied().unwrap_or(0);
        let allele = (coverage >= a.combine_null_min_coverage).then_some(0);
        let gt = SampleGenotype { a: allele, b: allele, gq: 0, dr: coverage, dv: 0, phase: Phase::default() };
        genotypes.insert(s, SampleCall { gt, id: "NULL".into() });
    }
    if a.combine_pair_relabel {
        let threshold = a.combine_pair_relabel_threshold;
        let best = genotypes
            .values()
            .filter(|s| s.gt.gq > threshold && s.gt.a.is_some())
            .map(|s| (s.gt.a, s.gt.b))
            .fold((Some(0), Some(0)), std::cmp::max);
        if best != (Some(0), Some(0)) {
            for s in genotypes.values_mut().filter(|s| s.gt.gq < threshold && s.gt.a.is_some()) {
                (s.gt.a, s.gt.b) = best;
            }
        }
    }

    let svlen_median = median_int(g.candidates.iter().map(|c| c.svlen));
    let mut alt = &first.alt;
    if first.svtype == SvType::Ins {
        let mut mindist = (first.alt.len() as i64 - svlen_median).abs();
        for c in &g.candidates {
            let dist = (c.alt.len() as i64 - svlen_median).abs();
            if dist < mindist {
                mindist = dist;
                alt = &c.alt;
            }
        }
    }
    let alt =
        if *alt == format!("<{}>", first.svtype.as_str()) { Alt::Symbolic } else { Alt::Seq(alt.clone().into_bytes()) };
    let n = g.candidates.len();
    let mut info = if n_samples == 1 { first.info.clone() } else { Vec::new() };
    if n_samples != 1 {
        let pos: Vec<i64> = g.candidates.iter().map(|c| c.pos).collect();
        let len: Vec<i64> = g.candidates.iter().map(|c| c.svlen).collect();
        info.push(("STDEV_POS".into(), Info::from(pymath::util_stdev(&pos))));
        info.push(("STDEV_LEN".into(), Info::from(pymath::util_stdev(&len))));
    }
    // Breakends have no length (SVLEN 0): as in single-sample calling, the screen is for the
    // other types. Sniffles2 applies it to them too and so drops every BND when combining
    // (LRP-007).
    if first.svtype != SvType::Bnd && first.svlen.abs() < cfg.minsvlen_screen {
        return None;
    }
    let call = CombinedCall {
        contig: first.contig.clone(),
        svtype: first.svtype,
        pos: first.pos,
        end: first.end,
        svlen: first.svlen,
        id: format!("{}.{:X}M{task_id:X}", first.svtype.as_str(), *sv_id),
        alt,
        qual: mean_or_none_round(g.candidates.iter().map(|c| c.qual)),
        filter: if n_samples != 1 { "PASS".into() } else { first.filter.clone() },
        info,
        genotypes: genotypes.into_values().collect(),
        precise: g.candidates.iter().filter(|c| c.precise).count() as f64 / n as f64 > 0.5,
        support: pymath::py_round(g.candidates.iter().map(|c| c.support).sum::<i64>() as f64 / n as f64),
        rnames,
        fwd: g.candidates.iter().map(|c| c.fwd).sum(),
        rev: g.candidates.iter().map(|c| c.rev).sum(),
        coverage: std::array::from_fn(|k| mean_or_none_round(g.candidates.iter().map(|c| c.coverage[k]))),
    };
    *sv_id += 1;
    Some(call)
}

/// The block starts of a contig or its regions (`CombineTask.generate_blocks`).
fn block_starts(len: i64, regions: &[(i64, i64)]) -> Vec<i64> {
    let step = snf::BLOCK_SIZE as usize;
    if regions.is_empty() {
        return (0..len - 1 + snf::BLOCK_SIZE).step_by(step).collect();
    }
    let mut set = BTreeSet::new();
    for &(start, end) in regions {
        set.extend((start / snf::BLOCK_SIZE * snf::BLOCK_SIZE..end + snf::BLOCK_SIZE).step_by(step));
    }
    set.into_iter().collect()
}

/// Combines one contig of all samples (`CombineTask.execute`); calls are in Sniffles2's order.
///
/// # Errors
///
/// Fails on unreadable SNF blocks.
pub fn combine_contig(
    samples: &[Sample],
    contig: &str,
    len: i64,
    regions: &[(i64, i64)],
    task_id: usize,
    cfg: &Config,
) -> Result<Vec<CombinedCall>> {
    let n_samples = samples.len();
    let bin_max_candidates = MIN_BIN_CANDIDATES.max(n_samples / 2);
    let mut sv_id = 0;
    let mut keep_by_type: [Vec<Group>; 5] = Default::default();
    let mut calls = Vec::new();
    for block_start in block_starts(len, regions) {
        let mut blocks: Vec<Option<Vec<Block>>> = samples
            .par_iter()
            .enumerate()
            .map(|(k, s)| s.file.read_blocks(contig, block_start, k))
            .collect::<Result<_>>()?;
        for svtype in SvType::ALL {
            let mut bins: BTreeMap<i64, Vec<Candidate>> = BTreeMap::new();
            for parts in blocks.iter_mut().flatten() {
                for block in parts {
                    for c in std::mem::take(&mut block.candidates[svtype as usize]) {
                        if c.support < cfg.args.combine_support_threshold {
                            continue;
                        }
                        let bin = (c.pos as f64 / BIN as f64) as i64 * BIN;
                        bins.entry(bin).or_default().push(c);
                    }
                }
            }
            let Some(&last_bin) = bins.keys().next_back() else { continue };
            let mut keep = std::mem::take(&mut keep_by_type[svtype as usize]);
            let mut size = 0;
            let mut pending: Vec<Candidate> = Vec::new();
            for (bin, cands) in bins {
                pending.extend(cands);
                size += BIN;
                if pending.len() < bin_max_candidates && bin != last_bin {
                    continue;
                }
                let groups = resolve_groups(svtype, std::mem::take(&mut pending), std::mem::take(&mut keep), cfg);
                let mut to_call = Vec::new();
                for mut g in groups {
                    let window = (g.pos_mean / snf::COVERAGE_WINDOW as f64) as i64 * snf::COVERAGE_WINDOW;
                    for (s, parts) in blocks.iter().enumerate() {
                        if g.included.contains(&s) {
                            continue;
                        }
                        let coverage =
                            parts.as_ref().and_then(|p| p.first()).and_then(|b| b.coverage.get(&window)).copied();
                        let slot = g.coverage_nonincluded.entry(s).or_insert(0);
                        *slot = (*slot).max(coverage.unwrap_or(0));
                    }
                    if (g.pos_mean - bin as f64).abs() < (size as f64 * 0.5).max(OVERLAP_ABS) {
                        keep.push(g);
                    } else {
                        to_call.push(g);
                    }
                }
                calls.extend(to_call.iter().filter_map(|g| call_group(g, cfg, n_samples, task_id, &mut sv_id)));
                size = 0;
            }
            keep_by_type[svtype as usize] = keep;
        }
    }
    for keep in &keep_by_type {
        calls.extend(keep.iter().filter_map(|g| call_group(g, cfg, n_samples, task_id, &mut sv_id)));
    }
    if cfg.sort {
        calls.sort_by_key(|c| c.pos);
    }
    Ok(calls)
}

/// Appends a combined call's record (multi-sample `VCF.write_call`).
///
/// # Errors
///
/// Fails if the reference cannot be read.
pub fn write_call(
    out: &mut String,
    call: &CombinedCall,
    n_samples: usize,
    cfg: &Config,
    reference: Reference<'_>,
) -> std::io::Result<()> {
    let mut info = call.info.clone();
    let mut filter = call.filter.as_str();
    if n_samples > 1 {
        let mut ac = 0;
        let mut supp_vec = String::with_capacity(n_samples);
        for s in &call.genotypes {
            if let (Some(a), Some(b)) = (s.gt.a, s.gt.b)
                && s.gt.dv > 0
            {
                ac += a + b;
                supp_vec.push('1');
            } else {
                supp_vec.push('0');
            }
        }
        if !supp_vec.contains('1') {
            return Ok(());
        }
        if ac == 0 {
            filter = "GT";
        }
        info.push(("AC".into(), Info::Int(ac)));
        info.push(("SUPP_VEC".into(), Info::Str(supp_vec)));
    }
    let record = Record {
        contig: &call.contig,
        svtype: call.svtype,
        pos: call.pos,
        end: call.end,
        svlen: call.svlen,
        id: &call.id,
        alt: &call.alt,
        qual: call.qual,
        filter,
        precise: call.precise,
        support: call.support,
        rnames: call.rnames.iter().map(|n| Cow::Borrowed(n.as_str())).collect(),
        coverage: call.coverage,
        fwd: call.fwd,
        rev: call.rev,
        nm: Info::Int(-1),
        svlens: None,
        info: info.iter().map(|(k, v)| (k.as_str(), v)).collect(),
    };
    if !vcf::write_record(out, record, "GT:GQ:DR:DV:PS:ID", cfg, reference)? {
        return Ok(());
    }
    for (k, s) in call.genotypes.iter().enumerate() {
        if k > 0 {
            out.push('\t');
        }
        vcf::format_sample_genotype(out, &s.gt);
        out.push(':');
        out.push_str(&s.id);
    }
    out.push('\n');
    Ok(())
}

/// Global edit distance (Levenshtein) of two byte strings, with Myers' bit-parallel algorithm in
/// 64-row blocks over the whole matrix (the reference for the banded search).
#[cfg(test)]
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let (query, target) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let m = query.len();
    if m == 0 {
        return target.len();
    }
    let n_blocks = m.div_ceil(64);
    // Match masks per byte value and block.
    let mut peq = vec![0u64; 256 * n_blocks];
    for (i, &c) in query.iter().enumerate() {
        peq[c as usize * n_blocks + i / 64] |= 1 << (i % 64);
    }
    let mut pv = vec![u64::MAX; n_blocks];
    let mut mv = vec![0u64; n_blocks];
    let last_bit = 1u64 << ((m - 1) % 64);
    let mut score = m as i64;
    for &c in target {
        let eq_row = &peq[c as usize * n_blocks..(c as usize + 1) * n_blocks];
        // D[0][j] = j: each column enters the top block with a +1 horizontal delta.
        let mut hin: i64 = 1;
        for b in 0..n_blocks {
            let (p, mvb) = (pv[b], mv[b]);
            let mut eq = eq_row[b];
            let xv = eq | mvb;
            if hin < 0 {
                eq |= 1;
            }
            let xh = ((eq & p).wrapping_add(p) ^ p) | eq;
            let mut ph = mvb | !(xh | p);
            let mut mh = p & xh;
            let high = if b + 1 == n_blocks { last_bit } else { 1 << 63 };
            let hout = i64::from(ph & high != 0) - i64::from(mh & high != 0);
            ph <<= 1;
            mh <<= 1;
            if hin < 0 {
                mh |= 1;
            } else if hin > 0 {
                ph |= 1;
            }
            pv[b] = mh | !(xv | ph);
            mv[b] = ph & xv;
            hin = hout;
        }
        score += hin;
    }
    score as usize
}

/// One 64-row block of the bit-parallel edit distance: vertical deltas and the bottom cell.
#[derive(Debug, Clone, Copy)]
struct EdBlock {
    p: u64,
    m: u64,
    score: i64,
}

/// Advances a block by one column (edlib's `calculateBlock`); returns the horizontal delta out of
/// its bottom row.
fn advance_block(b: &mut EdBlock, mut eq: u64, hin: i64) -> i64 {
    let (pv, mv) = (b.p, b.m);
    let xv = eq | mv;
    if hin < 0 {
        eq |= 1;
    }
    let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
    let mut ph = mv | !(xh | pv);
    let mut mh = pv & xh;
    let hout = (ph >> 63) as i64 - (mh >> 63) as i64;
    ph <<= 1;
    mh <<= 1;
    if hin < 0 {
        mh |= 1;
    } else if hin > 0 {
        ph |= 1;
    }
    b.p = mh | !(xv | ph);
    b.m = ph & xv;
    hout
}

/// A query prepared for edit distance searches: match masks per byte value and 64-row block,
/// padded to whole blocks with rows that match any byte (as edlib).
#[derive(Debug)]
pub struct EdQuery {
    peq: Vec<u64>,
    len: i64,
    n_blocks: usize,
}

impl EdQuery {
    pub fn new(query: &[u8]) -> Self {
        let n_blocks = query.len().div_ceil(64).max(1);
        let mut peq = vec![0u64; 256 * n_blocks];
        let pad_from = query.len() % 64;
        if pad_from != 0 || query.is_empty() {
            let pad = if query.is_empty() { u64::MAX } else { u64::MAX << pad_from };
            for c in 0..256 {
                peq[c * n_blocks + n_blocks - 1] |= pad;
            }
        }
        for (r, &c) in query.iter().enumerate() {
            peq[c as usize * n_blocks + r / 64] |= 1 << (r % 64);
        }
        EdQuery { peq, len: query.len() as i64, n_blocks }
    }
}

/// The global edit distance if it is at most `k`, else `None`: edlib's banded Myers algorithm
/// (`myersCalcEditDistanceNW`).
pub fn edit_distance_within(query: &EdQuery, target: &[u8], k: usize) -> Option<usize> {
    const W: i64 = 64;
    let (m, n) = (query.len, target.len() as i64);
    let mut k = k as i64;
    if k < (n - m).abs() {
        return None;
    }
    if m == 0 {
        return Some(n as usize);
    }
    let nb = query.n_blocks;
    let n_blocks = nb as i64;
    let padding = n_blocks * W - m;
    let peq = &query.peq;
    let mut first = 0i64;
    let mut last = (n_blocks).min((k.min((k + m - n) / 2) + 1 + W - 1) / W) - 1;
    let mut blocks = vec![EdBlock { p: u64::MAX, m: 0, score: 0 }; nb];
    for (b, block) in blocks.iter_mut().enumerate().take((last + 1).max(0) as usize) {
        block.score = (b as i64 + 1) * W;
    }
    for (c, &t) in target.iter().enumerate() {
        let c = c as i64;
        let eq_col = &peq[t as usize * nb..(t as usize + 1) * nb];
        let mut hout = 1;
        for b in first..=last {
            let block = &mut blocks[b as usize];
            hout = advance_block(block, eq_col[b as usize], hout);
            block.score += hout;
        }
        let bottom = blocks[last as usize].score;
        k = k.min(
            bottom + (n - c - 1).max(m - ((1 + last) * W - 1) - 1) + if last == n_blocks - 1 { padding } else { 0 },
        );
        if last + 1 < n_blocks && (last + 1) * W - 1 <= k - bottom + 2 * W - 2 - n + c + m {
            last += 1;
            let above = blocks[(last - 1) as usize].score;
            let block = &mut blocks[last as usize];
            *block = EdBlock { p: u64::MAX, m: 0, score: 0 };
            let new_hout = advance_block(block, eq_col[last as usize], hout);
            block.score = above - hout + W + new_hout;
        }
        while last >= first
            && (blocks[last as usize].score >= k + W
                || (last + 1) * W - 1 > k - blocks[last as usize].score + 2 * W - 2 - n + c + m + 1)
        {
            last -= 1;
        }
        while first <= last
            && (blocks[first as usize].score >= k + W
                || (first + 1) * W - 1 < blocks[first as usize].score - k - n + m + c)
        {
            first += 1;
        }
        if last < first {
            return None;
        }
    }
    if last != n_blocks - 1 {
        return None;
    }
    // The cell `padding` rows above the bottom of the last block is the query's last row.
    let block = blocks[last as usize];
    let mut score = block.score;
    let mut mask = 1u64 << 63;
    for _ in 0..padding {
        if block.p & mask != 0 {
            score -= 1;
        }
        if block.m & mask != 0 {
            score += 1;
        }
        mask >>= 1;
    }
    (score <= k).then_some(score as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snf_contig_lists_combine_when_their_lengths_agree() {
        let contig = |name: &str, length: i64| (name.to_owned(), length);
        let genome = [contig("chr1", 248), contig("chr2", 242), contig("chr20", 64)];
        let panel = [contig("chr20", 64)];
        let other = [contig("chr20", 63)];
        // A panel's SNF lists only its contigs: the union is the genome's, in its order.
        let both = union_of_contigs([("panel.snf", &panel[..]), ("genome.snf", &genome[..])].into_iter()).unwrap();
        assert_eq!(both, [contig("chr20", 64), contig("chr1", 248), contig("chr2", 242)]);
        let same = union_of_contigs([("a.snf", &genome[..]), ("b.snf", &genome[..])].into_iter()).unwrap();
        assert_eq!(same, genome);
        let error = union_of_contigs([("a.snf", &panel[..]), ("b.snf", &other[..])].into_iter()).unwrap_err();
        assert!(error.to_string().contains("contig chr20 is 64 bp in a.snf and 63 bp in b.snf"), "{error}");
    }

    fn breakend(sample: usize, alt: &str, mate_pos: i64, is_first: bool, is_reverse: bool) -> Candidate {
        Candidate {
            svtype: SvType::Bnd,
            contig: "chr1".into(),
            pos: 2500,
            end: 2500,
            svlen: 0,
            id: format!("bnd{sample}"),
            alt: alt.into(),
            qual: Some(60),
            filter: "PASS".into(),
            info: Vec::new(),
            genotype: None,
            precise: true,
            support: 10,
            rnames: None,
            qc: true,
            fwd: 5,
            rev: 5,
            coverage: [Some(20); 5],
            bnd: Some(SnfBnd { mate_contig: "chrY".into(), mate_pos, is_first, is_reverse }),
            sample,
        }
    }

    fn group_samples(candidates: Vec<Candidate>) -> Vec<Vec<usize>> {
        let args = <crate::config::Args as clap::Parser>::parse_from(["lrpipe-sv", "-i", "x.snf", "--vcf", "x.vcf"]);
        let cfg = Config::from_args(args, String::new()).expect("valid config");
        resolve_groups(SvType::Bnd, candidates, Vec::new(), &cfg)
            .iter()
            .map(|g| g.candidates.iter().map(|c| c.sample).collect())
            .collect()
    }

    /// Two samples' breakends at the same place join only if they join the same way (LRP-016).
    #[test]
    fn breakends_combine_by_adjacency() {
        let same = group_samples(vec![
            breakend(0, "N[chrY:10001[", 10_001, true, false),
            breakend(1, "N[chrY:10001[", 10_051, true, false),
        ]);
        assert_eq!(same, vec![vec![0, 1]], "the same adjacency is one allele");
        for (alt, is_first, is_reverse) in
            [("N]chrY:10001]", true, true), ("[chrY:10001[N", false, false), ("]chrY:10001]N", false, true)]
        {
            let apart = group_samples(vec![
                breakend(0, "N[chrY:10001[", 10_001, true, false),
                breakend(1, alt, 10_001, is_first, is_reverse),
            ]);
            assert_eq!(apart, vec![vec![0], vec![1]], "{alt} is another adjacency");
        }
        let far = group_samples(vec![
            breakend(0, "N[chrY:10001[", 10_001, true, false),
            breakend(1, "N[chrY:900001[", 900_001, true, false),
        ]);
        assert_eq!(far, vec![vec![0], vec![1]], "a distant mate is another adjacency");
    }

    fn naive(a: &[u8], b: &[u8]) -> usize {
        let mut prev: Vec<usize> = (0..=b.len()).collect();
        for (i, &x) in a.iter().enumerate() {
            let mut cur = vec![i + 1; b.len() + 1];
            for (j, &y) in b.iter().enumerate() {
                cur[j + 1] = (prev[j] + usize::from(x != y)).min(prev[j + 1] + 1).min(cur[j] + 1);
            }
            prev = cur;
        }
        prev[b.len()]
    }

    #[test]
    fn edit_distance_matches_dynamic_programming() {
        // A small deterministic generator over a 4-letter alphabet, lengths across block edges.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..300 {
            let la = (next() % 200) as usize;
            let lb = (next() % 200) as usize;
            let a: Vec<u8> = (0..la).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let mut b: Vec<u8> = a.iter().copied().filter(|_| next() % 10 != 0).collect();
            b.extend((0..lb % 20).map(|_| b"ACGT"[(next() % 4) as usize]));
            let d = naive(&a, &b);
            assert_eq!(edit_distance(&a, &b), d, "{a:?} {b:?}");
            for k in [0, 1, 5, d.saturating_sub(1), d, d + 1, 64, 128, 1000] {
                assert_eq!(edit_distance_within(&EdQuery::new(&a), &b, k), (d <= k).then_some(d), "k={k} {a:?} {b:?}");
                assert_eq!(edit_distance_within(&EdQuery::new(&b), &a, k), (d <= k).then_some(d), "k={k} swapped");
            }
        }
        assert_eq!(edit_distance(b"", b"ACG"), 3);
        assert_eq!(edit_distance_within(&EdQuery::new(b"<DEL>"), b"<DEL>", 64), Some(0));
        assert_eq!(edit_distance(b"<DEL>", b"<DEL>"), 0);
    }

    #[test]
    fn median_truncates_like_python() {
        assert_eq!(median_int([3, 1, 2].into_iter()), 2);
        assert_eq!(median_int([-51, -50].into_iter()), -50);
        assert_eq!(median_int([1, 2].into_iter()), 1);
    }

    #[test]
    fn block_starts_cover_the_contig() {
        assert_eq!(block_starts(200_001, &[]), vec![0, 100_000, 200_000]);
        assert_eq!(block_starts(250_001, &[]), vec![0, 100_000, 200_000, 300_000]);
        assert_eq!(block_starts(1_000_000, &[(150_000, 160_000)]), vec![100_000, 200_000]);
    }
}

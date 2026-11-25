//! `preprocess/CheckEnvs.py`: which contigs are called, in how many chunks, and the BED each
//! chunk is restricted to.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use snv_merge::Gender;
use snv_pileup::AlignmentFile;

/// Bases added on each side of `--bed_fn` regions (`param.no_of_positions`).
const BED_EXPANSION: i64 = 33;
/// `MIN_CHUNK_LENGTH` and `MAX_CHUNK_LENGTH`, for the `--chunk_num` warnings.
const MIN_CHUNK_LENGTH: f64 = 200_000.0;
const MAX_CHUNK_LENGTH: f64 = 20_000_000.0;

/// `major_contigs_order`.
fn major_contigs() -> Vec<String> {
    let names: Vec<String> = (1..=22).map(|i| i.to_string()).chain(["X".to_owned(), "Y".to_owned()]).collect();
    names.iter().map(|n| format!("chr{n}")).chain(names.iter().cloned()).collect()
}

/// Whether lrpipe-snv calls contig `name` without `--ctg-name`, `--bed-fn` or
/// `--include-all-ctgs`: chr1-22, X and Y, or 1-22, X and Y (`major_contigs_order`).
#[must_use]
pub fn is_default_contig(name: &str) -> bool {
    major_contigs().iter().any(|m| m == name)
}

/// `(name, length)` of each `.fai` row.
fn fai_contigs(reference: &Path) -> Result<Vec<(String, i64)>> {
    let mut fai = reference.as_os_str().to_owned();
    fai.push(".fai");
    let text = fs::read_to_string(&fai).with_context(|| format!("reading {}", Path::new(&fai).display()))?;
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let mut cols = line.trim().split('\t');
            let name = cols.next().context("empty .fai row")?.to_owned();
            let len = cols.next().context("short .fai row")?.parse().context("bad .fai length")?;
            Ok((name, len))
        })
        .collect()
}

/// A plain or gzipped text file (`pigz -fdc` reads both).
pub(crate) fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let n = fs::File::open(path)
        .and_then(|mut f| f.read(&mut magic))
        .with_context(|| format!("reading {}", path.display()))?;
    let file = fs::File::open(path).with_context(|| format!("reading {}", path.display()))?;
    Ok(if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// `(contig, start, end)` of every BED row, validated as `bed_tree_from` does. Comment lines and
/// UCSC's `track` and `browser` lines are skipped.
pub(crate) fn bed_rows(path: &Path) -> Result<Vec<(String, i64, i64)>> {
    let mut rows = Vec::new();
    for (row_id, line) in open_text(path)?.lines().enumerate() {
        let line = line.with_context(|| format!("reading {}", path.display()))?;
        let first = line.split_whitespace().next().unwrap_or_default();
        if line.starts_with('#') || first.is_empty() || first == "track" || first == "browser" {
            continue;
        }
        let columns: Vec<&str> = line.split_whitespace().collect();
        let bad = || format!("Invalid bed input in {}-th row {line}", row_id + 1);
        let [contig, start, end, ..] = columns[..] else { bail!("{}: fewer than 3 columns", bad()) };
        let (start, end): (i64, i64) = (
            start.parse().with_context(|| format!("{}: the start is not a number", bad()))?,
            end.parse().with_context(|| format!("{}: the end is not a number", bad()))?,
        );
        if start < 0 {
            bail!("{}: the start is negative", bad());
        }
        if end < start {
            bail!("{}: the end is before the start", bad());
        }
        rows.push((contig.to_owned(), start, end));
    }
    Ok(rows)
}

/// The regions `tensor_generator_for_unchunked` calls in the original's `None` job
/// (`--chunk_num -1` with `--bed_fn`): per contig in first-seen order, `bed_tree_from`'s intervals
/// after `merge_overlaps`, in the order Python iterates the tree's `set` of intervals.
fn unchunked_regions(rows: &[(String, i64, i64)]) -> Vec<(String, i64, i64)> {
    let mut per_contig: Vec<(&str, Vec<(i64, i64)>)> = Vec::new();
    for (contig, start, end) in rows {
        // A zero-length row becomes one base; equal intervals are added once.
        let interval = (*start, if start == end { end + 1 } else { *end });
        match per_contig.iter_mut().find(|(c, _)| c == contig) {
            Some((_, intervals)) => intervals.push(interval),
            None => per_contig.push((contig, vec![interval])),
        }
    }
    let mut regions = Vec::new();
    for (contig, mut intervals) in per_contig {
        intervals.sort_unstable();
        intervals.dedup();
        // `merge_overlaps(strict=True)`: only overlapping intervals merge, not touching ones.
        let mut merged: Vec<(i64, i64)> = Vec::with_capacity(intervals.len());
        for (start, end) in intervals {
            match merged.last_mut() {
                Some(last) if start < last.1 => last.1 = last.1.max(end),
                _ => merged.push((start, end)),
            }
        }
        regions.extend(python_set_order(&merged).into_iter().map(|(start, end)| (contig.to_owned(), start, end)));
    }
    regions
}

/// The order CPython 3.11 iterates `set(intervals)` built from `intervals` in turn, where an
/// `Interval` hashes as the tuple `(begin, end)` (`setobject.c`: linear probing, then
/// perturbation; resized to four times its size when two thirds full).
fn python_set_order(intervals: &[(i64, i64)]) -> Vec<(i64, i64)> {
    const MIN_SIZE: usize = 8;
    const LINEAR_PROBES: usize = 9;
    const PERTURB_SHIFT: u32 = 5;
    // Inserting into a table without equal keys or dummies, as `set_insert_clean` does;
    // `set_add_entry` probes the same slots when the keys are distinct.
    fn insert(table: &mut [Option<(u64, (i64, i64))>], hash: u64, key: (i64, i64)) {
        let mask = table.len() - 1;
        let mut perturb = hash;
        let mut i = hash as usize & mask;
        loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            if let Some(free) = (i..=i + probes).find(|&j| table[j].is_none()) {
                table[free] = Some((hash, key));
                return;
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }
    let mut table: Vec<Option<(u64, (i64, i64))>> = vec![None; MIN_SIZE];
    let mut used = 0usize;
    for &key in intervals {
        insert(&mut table, python_pair_hash(key), key);
        used += 1;
        let mask = table.len() - 1;
        if used * 5 >= mask * 3 {
            let min_used = if used > 50_000 { used * 2 } else { used * 4 };
            let mut size = MIN_SIZE;
            while size <= min_used {
                size <<= 1;
            }
            let old = std::mem::replace(&mut table, vec![None; size]);
            for (hash, key) in old.into_iter().flatten() {
                insert(&mut table, hash, key);
            }
        }
    }
    table.into_iter().flatten().map(|(_, key)| key).collect()
}

/// `hash((begin, end))` for non-negative ints below 2^61, which hash to themselves
/// (`tuplehash`'s xxHash-based mixing on a 64-bit build).
fn python_pair_hash((begin, end): (i64, i64)) -> u64 {
    const PRIME_1: u64 = 11_400_714_785_074_694_791;
    const PRIME_2: u64 = 14_029_467_366_897_019_727;
    const PRIME_5: u64 = 2_870_177_450_012_600_261;
    let mut acc = PRIME_5;
    for lane in [begin, end] {
        acc = acc.wrapping_add((lane as u64).wrapping_mul(PRIME_2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(PRIME_1);
    }
    acc = acc.wrapping_add(2 ^ (PRIME_5 ^ 3_527_539));
    if acc == u64::MAX { 1_546_275_796 } else { acc }
}

/// The sites of a `--vcf_fn` VCF.
#[derive(Debug, Clone, Default)]
struct KnownVcf {
    /// `split_extend_vcf`: per contig, in first-seen order, the rows written to `split_beds`.
    split_rows: Vec<(String, Vec<String>)>,
    /// `vcf_candidates_from`: each contig's sorted, distinct 1-based positions.
    sites: HashMap<String, Vec<i64>>,
}

impl KnownVcf {
    fn read(path: &Path) -> Result<Self> {
        let mut known = Self::default();
        for (row_id, line) in open_text(path)?.lines().enumerate() {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            if line.starts_with('#') {
                continue;
            }
            let mut columns = line.split_whitespace();
            let (Some(contig), Some(pos)) = (columns.next(), columns.next()) else {
                bail!("Invalid VCF input in {}-th row {line}", row_id + 1);
            };
            let pos: i64 = pos.parse().with_context(|| format!("Invalid VCF input in {}-th row {line}", row_id + 1))?;
            known.sites.entry(contig.to_owned()).or_default().push(pos);
            let (start, end) = (pos - 1, pos);
            if start < 0 {
                bail!("Invalid VCF input in {}-th row {contig} {pos}", row_id + 1);
            }
            if start - BED_EXPANSION < 0 {
                continue;
            }
            let row = format!("{contig} {} {}", start - BED_EXPANSION, end + BED_EXPANSION);
            match known.split_rows.iter_mut().find(|(c, _)| c == contig) {
                Some((_, rows)) => rows.push(row),
                None => known.split_rows.push((contig.to_owned(), vec![row])),
            }
        }
        for sites in known.sites.values_mut() {
            sites.sort_unstable();
            sites.dedup();
        }
        Ok(known)
    }
}

/// The contigs to call and where their chunk lists and BEDs went.
#[derive(Debug, Clone)]
pub(crate) struct Contigs {
    /// `(contig, chunk count)` in calling order. With `--chunk_num -1` every count is 0: the
    /// original then runs one job for a contig named `None` (`CHUNK_LIST` is `None 0 0`), which
    /// calls nothing.
    pub(crate) chunks: Vec<(String, i64)>,
    /// `tmp/split_beds`, when `--bed_fn` or `--vcf_fn` was given.
    pub(crate) split_beds: Option<PathBuf>,
    /// With `--vcf_fn`, each contig's sorted known sites.
    pub(crate) known_sites: Option<HashMap<String, Vec<i64>>>,
    /// With `--chunk_num -1` and `--bed_fn`, the regions the original's `None` job calls into
    /// `pileup_None_0.vcf`, which its sort step then ignores; otherwise empty.
    pub(crate) unchunked: Vec<(String, i64, i64)>,
}

impl Contigs {
    /// Contig names in calling order.
    pub(crate) fn names(&self) -> Vec<&str> {
        self.chunks.iter().map(|(c, _)| c.as_str()).collect()
    }

    /// The `--extend_bed` of a contig's pileup jobs, if it exists.
    pub(crate) fn extend_bed(&self, contig: &str) -> Option<PathBuf> {
        self.split_beds.as_ref().map(|dir| dir.join(contig)).filter(|p| p.is_file())
    }
}

/// What [`check_envs`] reads.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ContigRequest<'a> {
    pub(crate) bam: &'a Path,
    pub(crate) reference: &'a Path,
    /// `--bed_fn`'s rows, as [`bed_rows`] read them.
    pub(crate) bed: Option<&'a [(String, i64, i64)]>,
    pub(crate) vcf: Option<&'a Path>,
    pub(crate) ctg_name: Option<&'a str>,
    pub(crate) include_all_ctgs: bool,
    pub(crate) min_contig_size: i64,
    pub(crate) chunk_size: i64,
    pub(crate) chunk_num: i64,
    pub(crate) gender: Gender,
}

/// `CheckEnvs`: chooses the contigs, writes `tmp/CONTIGS`, `tmp/CHUNK_LIST` and, with a BED or
/// a VCF of known sites, `tmp/split_beds/<contig>` (the regions or sites widened by 33 bases on
/// each side).
///
/// # Errors
///
/// Fails if an input cannot be read or the BED is malformed.
pub(crate) fn check_envs(request: &ContigRequest<'_>, tmp: &Path) -> Result<Contigs> {
    let major = major_contigs();
    if request.bed.is_some() && request.vcf.is_some() {
        bail!("Please provide either --vcf_fn or --bed_fn only");
    }
    let bed = request.bed;
    let known = request.vcf.map(KnownVcf::read).transpose()?;
    let split_dir = tmp.join("split_beds");
    // The original creates the folder for either input, even when it writes nothing into it.
    if bed.is_some() || known.is_some() {
        fs::create_dir_all(&split_dir).with_context(|| format!("creating {}", split_dir.display()))?;
    }
    if let Some(known) = &known {
        for (contig, rows) in &known.split_rows {
            let path = split_dir.join(contig);
            fs::write(&path, rows.join("\n")).with_context(|| format!("writing {}", path.display()))?;
        }
    }
    // The contigs in the BED, or with a site in the VCF.
    let bed_contigs: Option<HashSet<&str>> = match (&bed, &known) {
        (Some(rows), _) => Some(rows.iter().map(|r| r.0.as_str()).collect()),
        (None, Some(known)) => Some(known.split_rows.iter().map(|(c, _)| c.as_str()).collect()),
        (None, None) => None,
    };
    let requested: Option<HashSet<&str>> =
        request.ctg_name.filter(|s| !s.is_empty() && *s != "EMPTY").map(|s| s.split(',').collect());
    if requested.is_some() && bed.is_some() {
        log::warn!("both --ctg_name and --bed_fn provided, will only proceed contigs in intersection");
    }
    if requested.is_some() && known.is_some() {
        log::warn!("both --ctg_name and --vcf_fn provided, will only proceed contigs in intersection");
    }
    // The contigs named by --ctg_name and/or the BED.
    let mut contig_set: HashSet<String> = match (&requested, &bed_contigs) {
        (Some(names), Some(in_bed)) => names.intersection(in_bed).map(|s| (*s).to_owned()).collect(),
        (Some(names), None) => names.iter().map(|s| (*s).to_owned()).collect(),
        (None, Some(in_bed)) => in_bed.iter().map(|s| (*s).to_owned()).collect(),
        (None, None) => HashSet::new(),
    };

    let mut chosen: Vec<(String, i64, i64)> = Vec::new();
    for (name, len) in fai_contigs(request.reference)? {
        if !request.include_all_ctgs && bed_contigs.is_none() && requested.is_none() && !major.contains(&name) {
            continue;
        }
        if bed_contigs.as_ref().is_some_and(|c| !c.contains(name.as_str())) {
            continue;
        }
        if requested.is_some() && !contig_set.contains(&name) {
            continue;
        }
        if request.min_contig_size > 0 && len < request.min_contig_size {
            log::warn!(
                "{name} contig length {len} is smaller than minimum contig size {}, will skip it!",
                request.min_contig_size
            );
            contig_set.remove(&name);
            continue;
        }
        contig_set.insert(name.clone());
        let chunks =
            if len % request.chunk_size != 0 { len / request.chunk_size + 1 } else { len / request.chunk_size };
        chosen.push((name, len, chunks.max(1)));
    }
    let lengths: Vec<f64> = chosen.iter().map(|c| c.1 as f64).collect();

    // sorted(contig_set, key=(major_contigs_order + list(contig_set)).index). Python's set order
    // is random for the other contigs; they keep their reference order here.
    let rank = |name: &str| major.iter().position(|m| m == name).unwrap_or(major.len());
    let mut sorted: Vec<String> = contig_set.iter().cloned().collect();
    let fai_order: Vec<&str> = chosen.iter().map(|c| c.0.as_str()).collect();
    sorted.sort_by_key(|name| (rank(name), fai_order.iter().position(|c| c == name).unwrap_or(usize::MAX)));
    for name in &sorted {
        if !fai_order.contains(&name.as_str()) {
            log::warn!("Contig {name} given but not found in reference fai file");
        }
    }

    // check_contig_in_bam: contigs without mapped reads (by the index) are skipped.
    let bam = AlignmentFile::open(request.bam, request.reference)?;
    let mut called: Vec<String> = Vec::new();
    for name in &sorted {
        if bam.mapped_reads(name).is_none_or(|n| n > 0) && bam.tid(name).is_some() {
            called.push(name.clone());
        } else {
            log::warn!("Contig name {name} provided but no mapped reads in BAM, skip!");
        }
    }
    if called.is_empty() && !sorted.is_empty() {
        log::warn!("No mapped reads support in BAM for provided contigs set {}", sorted.join(" "));
    }
    if request.gender == Gender::Female {
        called.retain(|c| c != "chrY" && c != "Y");
    }
    let chunking = request.chunk_num > -1;
    let chunk_count = |name: &str| {
        if chunking { chosen.iter().find(|c| c.0 == name).map_or(1, |c| c.2) } else { 0 }
    };
    let chunks: Vec<(String, i64)> = called.iter().map(|name| (name.clone(), chunk_count(name))).collect();

    if !chunks.is_empty() {
        log::info!("Call variant in contigs: {}", called.join(" "));
        if chunking {
            log::info!(
                "Chunk number for each contig: {}",
                chunks.iter().map(|(_, n)| n.to_string()).collect::<Vec<_>>().join(" ")
            );
        } else {
            log::info!("No genome chunking due to chunk_num == -1");
        }
        let (min_len, max_len) = lengths.iter().fold((f64::INFINITY, 0f64), |(lo, hi), &l| (lo.min(l), hi.max(l)));
        if request.chunk_num > 0 {
            let (min_chunk, max_chunk) = (min_len / request.chunk_num as f64, max_len / request.chunk_num as f64);
            if max_chunk > MAX_CHUNK_LENGTH {
                log::warn!(
                    "Current maximum chunk size {min_chunk} is larger than default maximum chunk size {MAX_CHUNK_LENGTH}, \
                     You may set a larger chunk_num by setting --chunk_num=$ for better parallelism."
                );
            } else if min_chunk < MIN_CHUNK_LENGTH {
                log::warn!(
                    "Current minimum chunk size {min_chunk} is smaller than default minimum chunk size \
                     {MIN_CHUNK_LENGTH}, You may set a smaller chunk_num by setting --chunk_num=$."
                );
            }
        } else if request.chunk_num == 0 && max_len < request.chunk_size as f64 / 5.0 {
            log::warn!(
                "Current maximum contig length {max_len} is much smaller than default chunk size {}, You may set a \
                 smaller chunk size by setting --chunk_size=$ for better parallelism.",
                request.chunk_size
            );
        }
    }

    // split_extend_bed, for the contigs chosen before the BAM and --gender checks.
    let split_beds = match &bed {
        Some(rows) if !chunks.is_empty() && chunking => {
            let dir = split_dir;
            let mut per_contig: Vec<(&str, Vec<String>)> = Vec::new();
            for (contig, start, end) in rows.iter().filter(|r| contig_set.contains(&r.0)) {
                let row = format!("{contig} {} {}", (start - BED_EXPANSION).max(0), (end + BED_EXPANSION).max(0));
                match per_contig.iter_mut().find(|(c, _)| c == contig) {
                    Some((_, lines)) => lines.push(row),
                    None => per_contig.push((contig, vec![row])),
                }
            }
            for (contig, lines) in per_contig {
                let path = dir.join(contig);
                fs::write(&path, lines.join("\n")).with_context(|| format!("writing {}", path.display()))?;
            }
            Some(dir)
        }
        _ if known.is_some() => Some(split_dir),
        _ => None,
    };
    let unchunked = match &bed {
        Some(rows) if !chunking && !chunks.is_empty() => unchunked_regions(rows),
        _ => Vec::new(),
    };

    fs::write(tmp.join("CONTIGS"), called.join("\n"))?;
    let chunk_list: String = if chunking {
        chunks.iter().flat_map(|(c, n)| (1..=*n).map(move |id| format!("{c} {id} {n}\n"))).collect()
    } else {
        "None 0 0\n".to_owned()
    };
    fs::write(tmp.join("CHUNK_LIST"), chunk_list)?;
    Ok(Contigs { chunks, split_beds, known_sites: known.map(|k| k.sites), unchunked })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(contig: &str, intervals: &[(i64, i64)]) -> Vec<(String, i64, i64)> {
        intervals.iter().map(|&(start, end)| (contig.to_owned(), start, end)).collect()
    }

    #[test]
    fn unchunked_regions_follow_python_set_order() {
        // Orders printed by the original's env (Python 3.11.17) for the same rows, through
        // `shared.intervaltree`'s `addi`, `merge_overlaps` and `all_intervals`.
        let fixture = rows(
            "chr20",
            &[(30_100_000, 30_600_000), (31_000_000, 31_050_000), (33_000_000, 34_000_000), (35_000_000, 35_000_000)],
        );
        let starts: Vec<i64> = unchunked_regions(&fixture).iter().map(|r| r.1).collect();
        assert_eq!(starts, [35_000_000, 31_000_000, 33_000_000, 30_100_000]);
        assert_eq!(unchunked_regions(&fixture)[0].2, 35_000_001);

        // 30 rows resize the table twice; two rows merge into others and two only touch.
        let mut many: Vec<(i64, i64)> = (0..30).map(|i| (i * 1000, i * 1000 + 10)).collect();
        many.extend([(5000, 5500), (7000, 7000), (100_000, 100_010), (100_010, 100_020), (250_000, 250_000)]);
        let regions: Vec<(i64, i64)> = unchunked_regions(&rows("chr1", &many)).iter().map(|r| (r.1, r.2)).collect();
        let expected = [
            (19000, 19010),
            (4000, 4010),
            (26000, 26010),
            (11000, 11010),
            (18000, 18010),
            (3000, 3010),
            (5000, 5500),
            (25000, 25010),
            (10000, 10010),
            (2000, 2010),
            (17000, 17010),
            (16000, 16010),
            (1000, 1010),
            (100_000, 100_010),
            (24000, 24010),
            (9000, 9010),
            (7000, 7010),
            (8000, 8010),
            (23000, 23010),
            (15000, 15010),
            (0, 10),
            (250_000, 250_001),
            (22000, 22010),
            (100_010, 100_020),
            (29000, 29010),
            (14000, 14010),
            (21000, 21010),
            (6000, 6010),
            (20000, 20010),
            (12000, 12010),
            (27000, 27010),
            (28000, 28010),
            (13000, 13010),
        ];
        assert_eq!(regions, expected);
    }

    #[test]
    fn unchunked_regions_keep_bed_contig_order() {
        let mut bed = rows("chr2", &[(10, 20)]);
        bed.extend(rows("chr1", &[(5, 9)]));
        bed.extend(rows("chr2", &[(15, 30)]));
        let regions = unchunked_regions(&bed);
        assert_eq!(regions, [("chr2".to_owned(), 10, 30), ("chr1".to_owned(), 5, 9)]);
    }
}

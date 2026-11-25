//! Port of `preprocess/CreateTensorFullAlignment.py` as the original's Python pipeline runs it
//! (`--phasing_info_in_bam`, reads haplotagged into `phase_bam/<contig>.bam`).
//!
//! For Illumina the reads are piped from `RealignReads` through `samtools mpileup -
//! --reverse-del --output-QNAME --output-MQ --output-extra HP`, and every column whose alleles
//! pass the allele-frequency thresholds is a candidate. For ONT and HiFi (`--disable_c_impl`)
//! `samtools mpileup` reads the haplotagged BAM itself, limited to the candidate BED (`-l`),
//! and the candidates are the BED rows' centres.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use snv_mpileup::{AlignmentFile, Column, EXCLUDED_FLAGS, ReadFilter, Regions, text};
use snv_pileup::Fasta;

use crate::Error;
use crate::pyrandom::sample_seed0;
use crate::reads::{RealignInput, RealignedReads};

/// `param.flankingBaseNum`.
const FLANKING: i64 = 16;
/// `no_of_positions`.
const POSITIONS: usize = 2 * FLANKING as usize + 1;
/// `param_f.channel_size`.
const CHANNELS: usize = 8;
/// The sequencing platform (`--platform`): where reads come from, which positions are
/// candidates, and the rows of each matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Platform {
    /// Oxford Nanopore.
    Ont,
    /// PacBio HiFi.
    Hifi,
    /// Illumina, with realignment.
    Ilmn,
}

impl Platform {
    /// `param_f.matrix_depth_dict`: 89 rows for ONT, 55 otherwise.
    pub fn matrix_depth(self) -> usize {
        match self {
            Self::Ont => 89,
            Self::Hifi | Self::Ilmn => 55,
        }
    }
}
/// `extend_bp_distance = no_of_positions + param.extend_bp`.
const EXTEND_DISTANCE: i64 = POSITIONS as i64 + 10;
/// `CreateTensorFullAlignment`'s `--snp_min_af` and `--indel_min_af` defaults, which CallVarBam
/// does not override.
const MIN_SNP_AF: f64 = 0.08;
const MIN_INDEL_AF: f64 = 0.15;
/// `HAP_TYPE`: unphased 60, haplotype 1 30, haplotype 2 90.
fn hap_type(hap: u8) -> i32 {
    match hap {
        1 => 30,
        2 => 90,
        _ => 60,
    }
}

/// `ACGT_NUM`.
fn acgt_num(symbol: u8) -> Option<i32> {
    Some(match symbol {
        b'A' | b'N' => 100,
        b'C' => 25,
        b'G' => 75,
        b'T' => 50,
        b'+' => -50,
        b'-' => -100,
        b'*' | b'#' => 0,
        _ => return None,
    })
}

/// One full-alignment job (`CallVarBam --full_aln_regions BED` of the Python pipeline).
#[derive(Debug, Clone, Copy)]
pub struct FullAlignmentJob<'a> {
    /// `--platform`.
    pub platform: Platform,
    /// The haplotagged BAM (`tmp/phase_output/phase_bam/<contig>.bam`).
    pub bam: &'a Path,
    /// Indexed reference FASTA.
    pub reference: &'a Path,
    /// `--ctgName`.
    pub contig: &'a str,
    /// `--full_aln_regions`: the candidate BED.
    pub regions: &'a Path,
    /// `--minMQ` of `samtools mpileup`.
    pub min_mapq: u8,
    /// `--minCoverage`.
    pub min_coverage: i64,
    /// `--vcf_fn`: the contig's known sites, sorted; without BED centres only they become
    /// candidates.
    pub known_sites: Option<&'a [i64]>,
}

/// One candidate's full-alignment tensor, before padding to [`Platform::matrix_depth`] rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FullAlignmentCandidate {
    /// `"{contig}:{pos}:{centre reference base}"`.
    pub position: String,
    /// `"{depth}-{alleles}"`.
    pub alt_info: String,
    /// `rows x 33 x 8` values.
    pub tensor: Vec<i8>,
}

impl FullAlignmentCandidate {
    /// Reads in the tensor.
    pub fn rows(&self) -> usize {
        self.tensor.len() / (POSITIONS * CHANNELS)
    }
}

/// The model inputs `tensor_generator_from` builds: each tensor centred in `depth` rows of
/// zeros (`(depth - rows) // 2` before it).
pub fn padded_inputs(candidates: &[FullAlignmentCandidate], depth: usize) -> Vec<i8> {
    let row = POSITIONS * CHANNELS;
    let mut out = Vec::with_capacity(candidates.len() * depth * row);
    for c in candidates {
        let padding = depth.saturating_sub(c.rows());
        out.extend(std::iter::repeat_n(0, padding / 2 * row));
        out.extend_from_slice(&c.tensor);
        out.extend(std::iter::repeat_n(0, (padding - padding / 2) * row));
    }
    out
}

/// A read's entry in one column (`base_list` item).
#[derive(Debug, Clone)]
struct ReadBase {
    base: u8,
    indel: Vec<u8>,
}

/// A read's channels at one column and the bases it inserts after it (`read_info`).
type ReadInfo = ([i32; CHANNELS], Vec<u8>);

/// A pileup column kept for tensors (`Position`).
struct Site {
    ref_base: u8,
    /// `read_name_list`, as ids of [`Tensors::ids`].
    names: Vec<u32>,
    bases: Vec<ReadBase>,
    base_quality: Vec<i64>,
    mapping_quality: Vec<i64>,
    /// `read_info` once `update_infos` ran.
    info: Option<HashMap<u32, ReadInfo>>,
}

/// `decode_pileup_bases`: one entry per read base, an indel attached to the base before it.
fn decode_bases(bases: &[u8]) -> Result<Vec<ReadBase>, Error> {
    let mut out: Vec<ReadBase> = Vec::new();
    let mut i = 0;
    while i < bases.len() {
        let base = bases[i];
        if base == b'+' || base == b'-' {
            i += 1;
            let mut advance = 0usize;
            loop {
                let num = *bases.get(i).ok_or_else(|| Error::Exit("IndexError: string index out of range".into()))?;
                if !num.is_ascii_digit() {
                    break;
                }
                advance = advance * 10 + usize::from(num - b'0');
                i += 1;
            }
            let last = out.last_mut().ok_or_else(|| Error::Exit("IndexError: list index out of range".into()))?;
            last.indel = [&[base], &bases[i..(i + advance).min(bases.len())]].concat();
            i = (i + advance).wrapping_sub(1);
        } else if b"ACGTNacgtn#*".contains(&base) {
            out.push(ReadBase { base, indel: Vec::new() });
        } else if base == b'^' {
            i += 1;
        }
        i = i.wrapping_add(1);
    }
    Ok(out)
}

/// `dict[key] += count` on an insertion-ordered list.
fn add<K: PartialEq + Clone>(map: &mut Vec<(K, i64)>, key: &K, count: i64) {
    match map.iter_mut().find(|(k, _)| k == key) {
        Some((_, v)) => *v += count,
        None => map.push((key.clone(), count)),
    }
}

/// The allele-frequency test of `decode_pileup_bases`: `(depth, pass_af)`.
fn passes_af(bases: &[ReadBase], reference_base: u8) -> (i64, bool) {
    // `Counter(base + indel)`, keyed by the parts.
    let mut keys: Vec<((u8, &[u8]), i64)> = Vec::new();
    for b in bases {
        add(&mut keys, &(b.base, b.indel.as_slice()), 1);
    }
    let mut pileup: Vec<(u8, i64)> = Vec::new();
    let mut depth = 0;
    for &((base, indel), count) in &keys {
        let first = base.to_ascii_uppercase();
        if b"ACGT".contains(&first) {
            add(&mut pileup, &first, count);
            depth += count;
        }
        match indel.first() {
            Some(b'+') => add(&mut pileup, &b'I', count),
            Some(b'-') => add(&mut pileup, &b'D', count),
            _ => {}
        }
    }
    let denominator = if depth > 0 { depth } else { 1 } as f64;
    let mut sorted = pileup.clone();
    sorted.sort_by_key(|item| std::cmp::Reverse(item.1));
    let mut pass_af = sorted.first().is_some_and(|(item, _)| *item != reference_base);
    let (mut pass_snp, mut pass_indel) = (false, false);
    for &(item, count) in &sorted {
        if item == reference_base {
            continue;
        }
        if item == b'I' || item == b'D' {
            pass_indel = pass_indel || count as f64 / denominator >= MIN_INDEL_AF;
            continue;
        }
        pass_snp = pass_snp || count as f64 / denominator >= MIN_SNP_AF;
    }
    pass_af = pass_af || pass_snp || pass_indel;
    (depth, pass_af)
}

/// `int(100 * min(x, cap) / cap)`.
fn normalize(x: i64, cap: f64) -> i32 {
    (100.0 * (x as f64).min(cap) / cap) as i32
}

/// `get_tensor_info`: a read's 8 channels at one column and its inserted bases.
fn tensor_info(b: &ReadBase, bq: i32, ref_base: u8, mq: i32) -> ReadInfo {
    let mut channels = [0i32; CHANNELS];
    if b.base == b'*' || b.base == b'#' {
        return (channels, Vec::new());
    }
    let strand = if b"ACGT".contains(&b.base) { 100 } else { 50 };
    let upper = b.base.to_ascii_uppercase();
    let mut alt = 0;
    if !b.indel.is_empty() {
        alt = acgt_num(b.indel[0]).expect("indels start with + or -");
    } else if upper != ref_base && b"ACGT".contains(&upper) {
        alt = acgt_num(upper).expect("ACGT have values");
    }
    let reference = acgt_num(if b"ACGT".contains(&ref_base) { ref_base } else { b'A' }).expect("ACGT have values");
    let inserted = if b.indel.first() == Some(&b'+') { b.indel[1..].to_ascii_uppercase() } else { Vec::new() };
    channels[..5].copy_from_slice(&[reference, alt, strand, mq, bq]);
    (channels, inserted)
}

impl Site {
    /// `update_infos`.
    fn read_info(&mut self) -> &HashMap<u32, ReadInfo> {
        self.info.get_or_insert_with(|| {
            let mut info = HashMap::with_capacity(self.names.len());
            for (((&name, b), &bq), &mq) in
                self.names.iter().zip(&self.bases).zip(&self.base_quality).zip(&self.mapping_quality)
            {
                info.insert(name, tensor_info(b, normalize(bq, 40.0), self.ref_base, normalize(mq, 60.0)));
            }
            info
        })
    }
}

/// The candidate BED as `CreateTensorFullAlignment` reads it: the smallest start + 1 and the
/// largest end + 1 of the contig's rows, and the centre of each row of three columns (rows with
/// more describe heterozygous SNPs).
fn bed_regions(path: &Path, contig: &str) -> Result<((i64, i64), Vec<i64>), Error> {
    let text = std::fs::read_to_string(path).map_err(|source| Error::Io { path: path.to_owned(), source })?;
    let mut span: Option<(i64, i64)> = None;
    let mut centers = Vec::new();
    for row in text.lines() {
        let cols: Vec<&str> = row.trim_end().split('\t').collect();
        if cols.first() != Some(&contig) || cols.len() < 3 {
            continue;
        }
        let parse = |s: &str| s.parse::<i64>().map_err(|_| Error::Exit(format!("bad BED row {row:?}")));
        let (start, end) = (parse(cols[1])? + 1, parse(cols[2])? + 1);
        span = Some(span.map_or((start, end), |(s, e)| (s.min(start), e.max(end))));
        if cols.len() == 3 {
            centers.push(if start == 1 { end - FLANKING - 2 } else { start + (end - start).div_euclid(2) - 1 });
        }
    }
    let span = span
        .ok_or_else(|| Error::Exit(format!("[ERROR] ctg_name({contig}) not exists in bed file({})", path.display())))?;
    Ok((span, centers))
}

/// The tensors of one Illumina full-alignment job, in the order the original writes them.
///
/// # Errors
///
/// Fails if an input cannot be read, or where the original exits.
pub fn full_alignment_candidates(job: &FullAlignmentJob<'_>) -> Result<Vec<FullAlignmentCandidate>, Error> {
    let contig = job.contig;
    let ilmn = job.platform == Platform::Ilmn;
    let ((ctg_start, ctg_end), centers) = bed_regions(job.regions, contig)?;
    // `need_phasing_pos_set`: Illumina skips the rows, and centres outside the span are dropped.
    let mut centers: Vec<i64> =
        if ilmn { Vec::new() } else { centers.into_iter().filter(|&c| c >= ctg_start && c <= ctg_end).collect() };
    centers.sort_unstable();
    centers.dedup();
    let fasta = Fasta::open(job.reference)?;
    let mut file = AlignmentFile::open(job.bam, job.reference)?;
    let regions = Regions::read(job.regions)?;
    let mut view_file = if ilmn { Some(AlignmentFile::open(job.bam, job.reference)?) } else { None };

    let reference_start = (ctg_start - 1000).max(1);
    let reference_end = ctg_end + 1000;
    let mut reference = fasta.fetch(
        contig,
        i32::try_from(reference_start - 1).map_err(|_| Error::Exit("position out of range".into()))?,
        i32::try_from(reference_end - 1).map_err(|_| Error::Exit("position out of range".into()))?,
    )?;
    reference.make_ascii_uppercase();

    let mut pileup = match view_file.as_mut() {
        Some(view_file) => {
            let reads = RealignedReads::new(
                view_file,
                &file,
                &fasta,
                RealignInput { contig, span: (ctg_start, ctg_end) },
                regions,
            )?;
            let filter = ReadFilter { excluded_flags: EXCLUDED_FLAGS, min_mapq: job.min_mapq, ..ReadFilter::default() };
            let lines = reads.map(|line| line.map_err(|e| snv_mpileup::Error::Exit(e.to_string())));
            file.pileup_sam_lines(lines, filter)?
        }
        None => {
            // `samtools mpileup BAM -r CTG:START-END -l BED`, 33 bases past the BED's span.
            let region = format!("{contig}:{}-{}", (ctg_start - POSITIONS as i64).max(1), ctg_end + POSITIONS as i64);
            let filter = ReadFilter {
                excluded_flags: EXCLUDED_FLAGS,
                min_mapq: job.min_mapq,
                bed: Some(regions),
                ..ReadFilter::default()
            };
            file.pileup(&region, filter)?
        }
    };
    let known: Option<HashSet<i64>> = job.known_sites.map(|s| s.iter().copied().collect());
    let mut tensors = Tensors {
        contig,
        reference: &reference,
        reference_start,
        matrix_depth: job.platform.matrix_depth(),
        sites: BTreeMap::new(),
        ids: HashMap::new(),
        hap: Vec::new(),
        hp_unusable: Vec::new(),
        seen: Vec::new(),
        stamp: 0,
        out: Vec::new(),
    };
    // With BED centres they are the candidates (`has_pileup_candidates`); otherwise the columns
    // that pass the thresholds, or the known sites, are.
    let has_centers = !centers.is_empty();
    let mut pending: Vec<i64> = centers;
    let mut next = 0usize;
    while let Some(column) = pileup.next_column()? {
        let pos = column.pos + 1;
        // The original's issue 105: realigned positions outside the reference are skipped.
        if ilmn && (pos < reference_start + FLANKING || pos > reference_end - FLANKING + 1) {
            continue;
        }
        let Some(site) = tensors.read_column(&column, pos)? else { continue };
        if !has_centers {
            let (depth, pass_af) = passes_af(&site.bases, site.ref_base);
            let is_candidate = match &known {
                Some(known) => known.contains(&pos),
                None => b"ACGT".contains(&site.ref_base) && pass_af && depth >= job.min_coverage,
            };
            if is_candidate {
                pending.push(pos);
            }
        }
        tensors.sites.insert(pos, site);
        if next < pending.len() && pos - pending[next] > EXTEND_DISTANCE {
            tensors.candidate(pending[next])?;
            tensors.prune(pending[next]);
            next += 1;
        }
    }
    while next < pending.len() {
        tensors.candidate(pending[next])?;
        tensors.prune(pending[next]);
        next += 1;
    }
    Ok(tensors.out)
}

/// The state of `CreateTensorFullAlignment`'s main loop.
struct Tensors<'a> {
    contig: &'a str,
    reference: &'a [u8],
    reference_start: i64,
    /// `param.matrix_depth_dict[platform]`.
    matrix_depth: usize,
    /// `pileup_dict`.
    sites: BTreeMap<i64, Site>,
    /// Read names by id, so columns hold `u32`s (Python keys its dicts by name, so reads with
    /// one name, such as mates, share an id).
    ids: HashMap<Box<str>, u32>,
    /// `hap_dict` by read id: the `HP` of the first column a read has 1 or 2 in, or 0.
    hap: Vec<u8>,
    /// Unpaired reads whose `HP` was read and is neither 1 nor 2. Their tags do not change
    /// between columns, so it is not read again (long reads carry long `MM` strings that the
    /// lookup scans); mates share a name but not their tags.
    hp_unusable: Vec<bool>,
    /// The `stamp` of the last window each read was counted in.
    seen: Vec<u32>,
    stamp: u32,
    out: Vec<FullAlignmentCandidate>,
}

impl Tensors<'_> {
    /// The id of read `name`.
    fn id(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = u32::try_from(self.hap.len()).expect("fewer than 2^32 reads");
        self.ids.insert(name.into(), id);
        self.hap.push(0);
        self.hp_unusable.push(false);
        self.seen.push(0);
        id
    }

    /// One row of `samtools_pileup_generator_from`, or `None` where it skips the row.
    fn read_column(&mut self, column: &Column<'_>, pos: i64) -> Result<Option<Site>, Error> {
        let bases = text::read_bases(column);
        let entries = column.entries().len();
        let mut names = Vec::with_capacity(entries);
        let (mut base_quality, mut mapping_quality) = (Vec::with_capacity(entries), Vec::with_capacity(entries));
        for entry in column.entries() {
            let record = entry.record();
            let q = if entry.qpos() < record.len() { record.qual(entry.qpos()) } else { 0 };
            base_quality.push(i64::from(q.saturating_add(33).min(126)) - 33);
            mapping_quality.push(i64::from(record.mapq().saturating_add(33).min(126)) - 33);
            let id = self.id(record.qname());
            names.push(id);
            // `if hap in '12' and read_name not in hap_dict: hap_dict[read_name] = int(hap)`.
            if self.hap[id as usize] != 0 || self.hp_unusable[id as usize] {
                continue;
            }
            let hap = match record.aux_integer(c"HP") {
                Some(1) => 1,
                Some(2) => 2,
                Some(12) => return Err(Error::Exit("unsupported HP value \"12\"".to_owned())),
                Some(_) => 0,
                None => match record.aux_text(c"HP").as_deref() {
                    Some("1") => 1,
                    Some("2") => 2,
                    Some(hp @ ("" | "12")) => return Err(Error::Exit(format!("unsupported HP value {hp:?}"))),
                    _ => 0,
                },
            };
            if hap == 0 {
                // 0x1: paired.
                self.hp_unusable[id as usize] = record.flag() & 0x1 == 0;
                continue;
            }
            self.hap[id as usize] = hap;
        }
        let ref_base = usize::try_from(pos - self.reference_start)
            .ok()
            .and_then(|i| self.reference.get(i))
            .copied()
            .ok_or_else(|| Error::Exit(format!("IndexError: no reference base at {}:{pos}", self.contig)))?;
        let read_bases = decode_bases(bases.as_bytes())?;
        if names.len() != read_bases.len() {
            return Ok(None);
        }
        Ok(Some(Site { ref_base, names, bases: read_bases, base_quality, mapping_quality, info: None }))
    }

    /// Forgets the columns more than 43 bp before `center`.
    fn prune(&mut self, center: i64) {
        self.sites = self.sites.split_off(&(center - EXTEND_DISTANCE));
    }

    /// `sorted_by_hap_read_name`: the window's reads in first-seen order, at most `matrix_depth`
    /// (picked with Python's `random.sample` seeded with 0), sorted by haplotype.
    fn window_reads(&mut self, center: i64) -> Vec<(u8, usize, u32)> {
        self.stamp += 1;
        let mut names: Vec<u32> = Vec::new();
        for p in center - FLANKING..=center + FLANKING {
            if let Some(site) = self.sites.get(&p) {
                for &name in &site.names {
                    if self.seen[name as usize] != self.stamp {
                        self.seen[name as usize] = self.stamp;
                        names.push(name);
                    }
                }
            }
        }
        if names.len() > self.matrix_depth {
            let mut indices = sample_seed0(names.len(), self.matrix_depth);
            indices.sort_unstable();
            names = indices.into_iter().map(|i| names[i]).collect();
        }
        // `(hap, order, name)`: the orders are distinct, so names never break ties.
        let mut sorted: Vec<(u8, usize, u32)> =
            names.into_iter().enumerate().map(|(order, name)| (self.hap[name as usize], order, name)).collect();
        sorted.sort_unstable();
        sorted
    }

    /// One candidate: `generate_tensor`.
    fn candidate(&mut self, center: i64) -> Result<(), Error> {
        if !self.sites.contains_key(&center) {
            return Ok(());
        }
        let reads = self.window_reads(center);
        if reads.is_empty() {
            return Ok(());
        }
        let flanked = flanked_sequence(self.reference, center, self.reference_start);
        let reference_base = flanked[FLANKING as usize];

        let site = &self.sites[&center];
        let mut alt: Vec<(Vec<u8>, i64)> = Vec::new();
        let (mut depth, mut ref_count) = (0i64, 0i64);
        for b in &site.bases {
            depth += 1;
            if b.base == b'#' || b.base == b'*' {
                continue;
            }
            let upper = b.base.to_ascii_uppercase();
            if let Some(&sign) = b.indel.first() {
                if sign == b'+' {
                    add(&mut alt, &[b"+".as_slice(), &[upper], &b.indel[1..].to_ascii_uppercase()].concat(), 1);
                } else {
                    add(&mut alt, &b.indel.to_ascii_uppercase(), 1);
                }
            } else if upper != reference_base {
                add(&mut alt, &vec![upper], 1);
            } else {
                ref_count += 1;
            }
        }
        // Without `--bed_fn` every candidate passes the confident-region test.
        let mut tensor = vec![[[0i32; CHANNELS]; POSITIONS]; reads.len()];
        let mut insertions: Vec<(usize, usize, Vec<u8>)> = Vec::new();
        let start = center - FLANKING;
        for p in start..=center + FLANKING {
            let Some(site) = self.sites.get_mut(&p) else { continue };
            let info = site.read_info();
            let offset = (p - start) as usize;
            for (row, (_, _, name)) in reads.iter().enumerate() {
                if let Some((channels, inserted)) = info.get(name) {
                    tensor[row][offset] = *channels;
                    if !inserted.is_empty() && p < center + FLANKING {
                        insertions.push((row, offset, inserted.clone()));
                    }
                }
            }
        }
        for (row, offset, inserted) in insertions {
            for (i, &base) in inserted.iter().take(POSITIONS - offset).enumerate() {
                tensor[row][offset + i][6] = acgt_num(base)
                    .ok_or_else(|| Error::Exit(format!("KeyError: inserted base {:?}", char::from(base))))?;
            }
        }
        let site = &self.sites[&center];
        let by_name: HashMap<u32, &ReadBase> = site.names.iter().copied().zip(&site.bases).collect();
        let lookup = |key: &[u8]| alt.iter().find(|(k, _)| k == key).map(|&(_, v)| v);
        for (row, (hap, _, name)) in reads.iter().enumerate() {
            let mut af = 0.0f64;
            if let Some(b) = by_name.get(name) {
                let upper = b.base.to_ascii_uppercase();
                let key: Vec<u8> = match b.indel.first() {
                    Some(b'+') => [b"+".as_slice(), &[upper], &b.indel[1..].to_ascii_uppercase()].concat(),
                    Some(_) => b.indel.to_ascii_uppercase(),
                    None => vec![upper],
                };
                if let Some(count) = lookup(&key) {
                    af = count as f64 / (depth as f64).max(1.0);
                }
            }
            let af = if af != 0.0 { (100.0 * af.min(1.0) / 1.0) as i32 } else { 0 };
            for position in &mut tensor[row] {
                if position[2] != 0 {
                    position[5] = af;
                    position[7] = hap_type(*hap);
                }
            }
        }

        let mut alleles: Vec<String> = Vec::new();
        for (key, count) in &alt {
            let label = match key[0] {
                b'+' => format!("I{}", String::from_utf8_lossy(&key[1..]).to_uppercase()),
                b'-' => {
                    let from = center - self.reference_start + 1;
                    let deleted = py_slice(self.reference, from, from + (key.len() - 1) as i64);
                    format!("D{}", String::from_utf8_lossy(deleted))
                }
                _ => format!("X{}", String::from_utf8_lossy(key)),
            };
            alleles.push(format!("{label} {count}"));
        }
        if ref_count > 0 {
            alleles.push(format!("R{} {ref_count}", char::from(reference_base)));
        }
        let values = tensor.iter().flatten().flatten();
        let tensor = values
            .map(|&v| i8::try_from(v).map_err(|_| Error::Exit(format!("tensor value {v} out of int8"))))
            .collect::<Result<Vec<i8>, Error>>()?;
        self.out.push(FullAlignmentCandidate {
            position: format!("{}:{center}:{}", self.contig, char::from(reference_base)),
            alt_info: format!("{depth}-{}", alleles.join(" ")),
            tensor,
        });
        Ok(())
    }
}

/// `sequence[start:stop]` with Python's slice bounds.
fn py_slice(sequence: &[u8], start: i64, stop: i64) -> &[u8] {
    let len = sequence.len() as i64;
    let bound = |i: i64| if i < 0 { (len + i).max(0) } else { i.min(len) } as usize;
    sequence.get(bound(start)..bound(stop)).unwrap_or_default()
}

/// `get_flanked_sequence`.
fn flanked_sequence(reference: &[u8], center: i64, reference_start: i64) -> Vec<u8> {
    let left = center - FLANKING - reference_start;
    let right = center + FLANKING + 1 - reference_start;
    let len = reference.len() as i64;
    let mut out = Vec::new();
    let mut from = left;
    if left < 0 {
        out.extend(std::iter::repeat_n(b'A', (-left) as usize));
        from = 0;
    }
    out.extend_from_slice(py_slice(reference, from, right));
    if right > len {
        out.extend(std::iter::repeat_n(b'A', (right - len) as usize));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indels_attach_to_the_previous_base() {
        let decoded = decode_bases(b"^5A+2CG$c-1n*#").expect("valid column");
        let got: Vec<(u8, Vec<u8>)> = decoded.into_iter().map(|b| (b.base, b.indel)).collect();
        assert_eq!(got, [(b'A', b"+CG".to_vec()), (b'c', b"-n".to_vec()), (b'*', vec![]), (b'#', vec![])]);
    }

    #[test]
    fn padding_centres_the_reads() {
        let c =
            FullAlignmentCandidate { position: String::new(), alt_info: String::new(), tensor: vec![1; 2 * 33 * 8] };
        let padded = padded_inputs(&[c], 5);
        assert_eq!(padded.len(), 5 * 33 * 8);
        assert_eq!(padded[33 * 8], 1);
        assert_eq!(padded[0], 0);
        assert_eq!(padded[3 * 33 * 8], 0);
    }
}

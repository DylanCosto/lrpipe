//! Port of `preprocess/CreateTensorPileup.py` in calling mode.
//!
//! Each `samtools mpileup` column becomes 18 counts (`generate_tensor`); a sliding window of 33
//! columns turns every candidate position into a tensor once 16 more columns have been read.
//! With `--gvcf`, `samtools mpileup -a` also writes the positions without reads, and every
//! position of the chunk goes to the non-variant blocks ([`GvcfChunk`]).

use std::collections::{HashMap, HashSet, VecDeque};
use std::path::Path;

use snv_pileup::{Bed, Candidate, CandidateAlleles, ChunkRegion, Fasta, FlankRequirement, Interval, Site as GvcfSite};

use crate::hts::{AlignmentFile, Positions, ReadFilter};
use crate::text::read_bases;
use crate::{EXCLUDED_FLAGS, Error, Regions};

/// `param.flankingBaseNum`.
const FLANKING: usize = 16;
/// `no_of_positions`.
const WINDOW: usize = 2 * FLANKING + 1;
/// `param.channel_size`.
const CHANNELS: usize = 18;
/// `param.expandReferenceRegion`.
const REFERENCE_PADDING: i64 = 1000;

/// Channels of the symbols `BASE2INDEX` maps (`param.channel`).
const INS: usize = 4;
const INS_MAX: usize = 5;
const DEL: usize = 6;
const DEL_MAX: usize = 7;
const REV_INS: usize = 13;
const REV_INS_MAX: usize = 14;
const REV_DEL: usize = 15;
const REV_DEL_MAX: usize = 16;

/// `BASE2INDEX[symbol]` for the single-character symbols.
fn channel(symbol: u8) -> Option<usize> {
    Some(match symbol {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'*' => 8,
        b'a' => 9,
        b'c' => 10,
        b'g' => 11,
        b't' => 12,
        b'#' => 17,
        _ => return None,
    })
}

/// The arguments of one `CreateTensorPileup` job.
#[derive(Debug, Clone)]
pub struct PileupJob<'a> {
    /// `--ctgName`.
    pub contig: &'a str,
    /// The chunk (`--chunk_id`/`--chunk_num`, `--ctgStart`/`--ctgEnd` or the known sites).
    pub chunk: ChunkRegion,
    /// `--extend_bed`: `samtools mpileup -l`, and the span chunked by [`ChunkRegion::BedChunk`].
    pub extend_bed: Option<&'a Path>,
    /// `--bed_fn`: candidates must overlap it.
    pub confident_bed: Option<&'a Path>,
    /// `--snp_min_af`.
    pub min_snp_af: f64,
    /// `--indel_min_af`.
    pub min_indel_af: f64,
    /// `--minCoverage`.
    pub min_coverage: i64,
    /// `--minMQ`.
    pub min_mapq: u8,
    /// `--call_snp_only`.
    pub alleles: CandidateAlleles,
    /// `--enable_variant_calling_at_sequence_head_and_tail` ([`FlankRequirement::None`]).
    pub flanks: FlankRequirement,
    /// `--fast_mode` with `--platform ont`.
    pub fast_mode: FastMode,
}

/// ONT's fast mode (`--fast_mode`, which `CreateTensorPileup` ignores on other platforms).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum FastMode {
    /// SNPs pass at `--snp_min_af`.
    #[default]
    Off,
    /// SNPs need an AF of at least `max(--snp_min_af, 0.15)` and 4 reads.
    On,
}

/// `chunk_size = span // count + 1 if span % count else span // count`.
fn chunk_size(span: i64, count: i64) -> i64 {
    if span % count != 0 { span / count + 1 } else { span / count }
}

/// One parsed read token: a base symbol, or an indel with its sign.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Token<'a> {
    /// One of `ACGTNacgtn#*`.
    Base(u8),
    /// `+` or `-` and the bases after the length.
    Indel(u8, &'a [u8]),
}

impl Token<'_> {
    /// The string Python counts (`base` or `sign + bases`).
    fn key(&self) -> Vec<u8> {
        match *self {
            Token::Base(b) => vec![b],
            Token::Indel(sign, bases) => [&[sign], bases].concat(),
        }
    }
}

/// Splits a bases column into one token per read base or indel, as the parsing loops of
/// `generate_tensor` and `decode_pileup_bases` do: `^` and the mapping quality after it are
/// skipped, and so are `$`, `.`, `,`, `>`, `<` and IUPAC letters.
///
/// # Errors
///
/// Fails where Python raises `IndexError` (a column ending inside an indel length).
pub(crate) fn tokens(bases: &[u8]) -> Result<Vec<Token<'_>>, Error> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < bases.len() {
        let base = bases[i];
        if b"ACGTNacgtn#*".contains(&base) {
            out.push(Token::Base(base));
        } else if base == b'+' || base == b'-' {
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
            out.push(Token::Indel(base, &bases[i..(i + advance).min(bases.len())]));
            // `base_idx += advance - 1`, then the loop's `base_idx += 1`.
            i = (i + advance).wrapping_sub(1);
        } else if base == b'^' {
            i += 1;
        }
        i = i.wrapping_add(1);
    }
    Ok(out)
}

/// One column's counts: what `generate_tensor` returns and calling uses.
struct Site {
    tensor: [i32; CHANNELS],
    /// `alt_dict`, in insertion order.
    alt: Vec<(String, i64)>,
    /// `pileup_dict`: reads by `A`, `C`, `G`, `T`, `I` and `D`.
    pileup: Vec<(u8, i64)>,
    depth: i64,
    pass_af: bool,
    max_del_length: usize,
}

impl Site {
    /// `n_ref` and `n_total` of the gVCF site info: the reads of `reference_base` (as read from
    /// the reference, not the `A` an `N` becomes) and the reads of every key.
    fn gvcf_counts(&self, reference_base: u8) -> (u64, u64) {
        let n_ref = self.pileup.iter().find(|(key, _)| *key == reference_base).map_or(0, |&(_, n)| n);
        let n_total: i64 = self.pileup.iter().map(|&(_, n)| n).sum();
        (n_ref.unsigned_abs(), n_total.unsigned_abs())
    }
}

/// What one chunk passes to `variantInfoCalculator` with `--gvcf`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GvcfChunk {
    /// `--ctgName`.
    pub contig: String,
    /// `ctg_start`.
    pub start: i64,
    /// `ctg_end`.
    pub end: i64,
    /// The chunk's positions, or [`GvcfSites::Empty`].
    pub sites: GvcfSites,
}

/// The sites of a [`GvcfChunk`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GvcfSites {
    /// The positions `start <= pos < end` that `samtools mpileup -a` wrote, in order.
    Sites(Vec<GvcfSite>),
    /// `samtools mpileup` wrote no column (`empty_pileup_flag`), so the original writes one
    /// block for the chunk (`write_empty_pileup`).
    Empty,
}

/// `dict[key] += count` on an insertion-ordered list.
fn add<K: PartialEq + Clone>(map: &mut Vec<(K, i64)>, key: &K, count: i64) {
    match map.iter_mut().find(|(k, _)| k == key) {
        Some((_, v)) => *v += count,
        None => map.push((key.clone(), count)),
    }
}

/// `Counter(items).items()`: counts in first-occurrence order.
pub(crate) fn counter<K: PartialEq + Clone>(items: impl IntoIterator<Item = K>) -> Vec<(K, i64)> {
    let mut counts = Vec::new();
    for item in items {
        add(&mut counts, &item, 1);
    }
    counts
}

/// `sorted(items, key=count, reverse=True)`: by count, descending, ties in order.
pub(crate) fn sorted_by_count<K: Clone>(items: &[(K, i64)]) -> Vec<(K, i64)> {
    let mut sorted = items.to_vec();
    sorted.sort_by_key(|item| std::cmp::Reverse(item.1));
    sorted
}

/// The thresholds `generate_tensor` compares against.
struct Thresholds {
    min_snp_af: f64,
    min_indel_af: f64,
    alleles: CandidateAlleles,
    fast_mode: FastMode,
}

/// `param.min_af_dict['ont']`, the SNP AF floor of fast mode.
const FAST_MODE_MIN_SNP_AF: f64 = 0.15;
/// Reads a SNP needs in fast mode.
const FAST_MODE_MIN_SNP_READS: i64 = 4;

/// `generate_tensor` for one column.
fn generate_tensor(
    pos: i64,
    bases: &[u8],
    reference: &[u8],
    reference_start: i64,
    reference_base: u8,
    thresholds: &Thresholds,
) -> Result<Site, Error> {
    // `evc_base_from` on an upper-case base.
    let reference_base = if b"ACGT".contains(&reference_base) { reference_base } else { b'A' };
    let ref_char = char::from(reference_base);
    let mut tensor = [0i32; CHANNELS];
    let mut alt: Vec<(String, i64)> = Vec::new();
    let mut pileup: Vec<(u8, i64)> = Vec::new();
    let (mut alt_count, mut ins_count, mut del_count) = (0i64, 0i64, 0i64);
    let (mut depth, mut max_ins_0, mut max_del_0, mut max_ins_1, mut max_del_1) = (0i64, 0, 0, 0, 0);
    let mut max_del_length = 0usize;

    let tokens = tokens(bases)?;
    for (token, count) in counter(tokens.iter().map(Token::key)) {
        match token[0] {
            b'+' => {
                let inserted = String::from_utf8_lossy(&token[1..]).to_uppercase();
                add(&mut alt, &format!("I{ref_char}{inserted}"), count);
                add(&mut pileup, &b'I', count);
                ins_count += count;
                if token.get(1).is_some_and(|b| b"ACGTN*".contains(b)) {
                    tensor[INS] += count as i32;
                    max_ins_0 = max_ins_0.max(count);
                } else {
                    tensor[REV_INS] += count as i32;
                    max_ins_1 = max_ins_1.max(count);
                }
            }
            b'-' => {
                let start = pos - reference_start + 1;
                let deleted = py_slice(reference, start, start + (token.len() - 1) as i64);
                add(&mut alt, &format!("D{}", String::from_utf8_lossy(deleted)), count);
                add(&mut pileup, &b'D', count);
                max_del_length = max_del_length.max(deleted.len());
                del_count += count;
                if token.get(1).is_some_and(|b| b"N*ACGT".contains(b)) {
                    tensor[DEL] += count as i32;
                    max_del_0 = max_del_0.max(count);
                } else {
                    tensor[REV_DEL] += count as i32;
                    max_del_1 = max_del_1.max(count);
                }
            }
            base => {
                let upper = base.to_ascii_uppercase();
                if b"ACGT".contains(&upper) {
                    add(&mut pileup, &upper, count);
                    depth += count;
                    if upper != reference_base {
                        add(&mut alt, &format!("X{}", char::from(upper)), count);
                        alt_count += count;
                    }
                    tensor[channel(base).expect("ACGTacgt have channels")] += count as i32;
                } else if base == b'#' || base == b'*' {
                    del_count += count;
                    tensor[channel(base).expect("#* have channels")] += count as i32;
                    depth += count;
                }
            }
        }
    }
    let ref_count = (depth - del_count - ins_count - alt_count).max(0);
    if ref_count > 0 {
        add(&mut alt, &format!("R{ref_char}"), ref_count);
    }
    tensor[INS_MAX] = max_ins_0 as i32;
    tensor[REV_INS_MAX] = max_ins_1 as i32;
    tensor[DEL_MAX] = max_del_0 as i32;
    tensor[REV_DEL_MAX] = max_del_1 as i32;
    let denominator = if depth > 0 { depth } else { 1 } as f64;
    let pileup_list = sorted_by_count(&pileup);

    let (mut pass_snp_af, mut pass_indel_af) = (false, false);
    let fast = thresholds.fast_mode == FastMode::On;
    let min_snp_af = if fast { thresholds.min_snp_af.max(FAST_MODE_MIN_SNP_AF) } else { thresholds.min_snp_af };
    let mut pass_af = pileup_list.first().is_some_and(|(item, _)| *item != reference_base);
    for &(item, count) in &pileup_list {
        if item == reference_base {
            continue;
        }
        if item == b'I' || item == b'D' {
            pass_indel_af = pass_indel_af || count as f64 / denominator >= thresholds.min_indel_af;
            continue;
        }
        pass_snp_af =
            pass_snp_af || (count as f64 / denominator >= min_snp_af && (!fast || count >= FAST_MODE_MIN_SNP_READS));
    }
    // The reference channel becomes minus the sum of all four bases, the reference included.
    let ref_channel = channel(reference_base).expect("ACGT have channels");
    let lower_channel = channel(reference_base.to_ascii_lowercase()).expect("acgt have channels");
    tensor[ref_channel] = -(tensor[0] + tensor[1] + tensor[2] + tensor[3]);
    tensor[lower_channel] = -(tensor[9] + tensor[10] + tensor[11] + tensor[12]);
    pass_af = match thresholds.alleles {
        CandidateAlleles::SnpsOnly => pass_snp_af,
        _ => pass_af || pass_snp_af || pass_indel_af,
    };
    Ok(Site { tensor, alt, pileup, depth, pass_af, max_del_length })
}

/// `sequence[start:stop]` with Python's slice bounds (negative values count from the end).
pub(crate) fn py_slice(sequence: &[u8], start: i64, stop: i64) -> &[u8] {
    let len = sequence.len() as i64;
    let bound = |i: i64| if i < 0 { (len + i).max(0) } else { i.min(len) } as usize;
    let (start, stop) = (bound(start), bound(stop));
    sequence.get(start..stop).unwrap_or_default()
}

/// `get_flanked_sequence`: the 33 reference bases around `center`, padded with `A`.
pub(crate) fn flanked_sequence(reference: &[u8], center: i64, reference_start: i64) -> Vec<u8> {
    let left = center - FLANKING as i64 - reference_start;
    let right = center + FLANKING as i64 + 1 - reference_start;
    let len = reference.len() as i64;
    if left >= 0 && right <= len {
        return reference[left as usize..right as usize].to_vec();
    }
    let mut out = Vec::new();
    let mut start = left;
    if left < 0 {
        out.extend(std::iter::repeat_n(b'A', (-left) as usize));
        start = 0;
    }
    out.extend_from_slice(py_slice(reference, start, right));
    if right > len {
        out.extend(std::iter::repeat_n(b'A', (right - len) as usize));
    }
    out
}

/// The reference bases `samtools faidx REF CONTIG:START-END` prints, upper-cased
/// (`reference_sequence_from`).
pub(crate) fn reference_sequence(fasta: &Fasta, contig: &str, start: i64, end: i64) -> Result<Vec<u8>, Error> {
    let to_i32 = |v: i64| i32::try_from(v).map_err(|_| Error::Exit(format!("position {v} out of range")));
    let mut sequence = fasta.fetch(contig, to_i32(start - 1)?, to_i32(end - 1)?)?;
    sequence.make_ascii_uppercase();
    if sequence.is_empty() {
        return Err(Error::Exit(format!("[ERROR] Failed to load reference sequence of {contig}")));
    }
    Ok(sequence)
}

/// The candidates of one chunk, like `CreateTensorPileup(args)` writes them, in order.
///
/// Positions are `"{contig}:{pos}:{center reference base}"`: the original writes the 33
/// flanking bases there, but variant output only reads the centre.
///
/// # Errors
///
/// Fails if a file cannot be read, or where the original exits.
pub fn pileup_candidates(bam: &mut AlignmentFile, fasta: &Fasta, job: &PileupJob<'_>) -> Result<Vec<Candidate>, Error> {
    Ok(candidates_and_sites(bam, fasta, job, Positions::Covered)?.0)
}

/// [`pileup_candidates`] with `--gvcf`: the candidates, and the chunk's sites for the
/// non-variant blocks, or `None` for a known-sites chunk without sites (where the original
/// returns before it opens the gVCF writer). The `samtools mpileup -a` columns without reads change
/// the candidates' flanks at the chunk's edges, as they do in the original.
///
/// # Errors
///
/// As [`pileup_candidates`].
pub fn pileup_candidates_with_sites(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    job: &PileupJob<'_>,
) -> Result<(Vec<Candidate>, Option<GvcfChunk>), Error> {
    candidates_and_sites(bam, fasta, job, Positions::All)
}

/// `CreateTensorPileup`, without `--gvcf` ([`Positions::Covered`]) or with it.
fn candidates_and_sites(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    job: &PileupJob<'_>,
    positions: Positions,
) -> Result<(Vec<Candidate>, Option<GvcfChunk>), Error> {
    let gvcf = positions == Positions::All;
    let contig = job.contig;
    let extend_bed = job.extend_bed.map(|p| Bed::read_contig(p, contig)).transpose()?;
    let (ctg_start, ctg_end) = match &job.chunk {
        ChunkRegion::Range { start, end } => (*start, *end),
        ChunkRegion::FaiChunk { id, count } => {
            let size = chunk_size(fasta.contig_len(contig)?, *count);
            (size * id, size * id + size)
        }
        ChunkRegion::BedChunk { id, count } => {
            let span = extend_bed.as_ref().and_then(Bed::span);
            let (bed_start, bed_end) =
                span.ok_or_else(|| Error::Exit(format!("[ERROR] ctg_name {contig} not exists in bed file")))?;
            let size = chunk_size(bed_end - bed_start, *count);
            let start = bed_start + 1 + size * id;
            (start, start + size)
        }
        ChunkRegion::KnownSites { sites } => match (sites.first(), sites.last()) {
            (Some(&first), Some(&last)) => (first, last),
            _ => return Ok((Vec::new(), None)),
        },
    };
    let known: Option<HashSet<i64>> = match &job.chunk {
        ChunkRegion::KnownSites { sites } => Some(sites.iter().copied().collect()),
        _ => None,
    };
    let extend_start = (ctg_start - WINDOW as i64).max(1);
    let extend_end = ctg_end + WINDOW as i64;
    let reference_start = (ctg_start - REFERENCE_PADDING).max(1);
    let reference_end = ctg_end + REFERENCE_PADDING;
    let reference = reference_sequence(fasta, contig, reference_start, reference_end)?;
    let confident: Option<Vec<Interval>> = job
        .confident_bed
        .map(|p| Bed::read_contig(p, contig).map(|bed| bed.overlapping(extend_start, extend_end)))
        .transpose()?;

    let filter = ReadFilter {
        excluded_flags: EXCLUDED_FLAGS,
        min_mapq: job.min_mapq,
        bed: job.extend_bed.map(Regions::read).transpose()?,
        positions,
    };
    let thresholds = Thresholds {
        min_snp_af: job.min_snp_af,
        min_indel_af: job.min_indel_af,
        alleles: job.alleles,
        fast_mode: job.fast_mode,
    };
    let head_and_tail = job.flanks == FlankRequirement::None;
    let mut window = Window::new(head_and_tail);
    let mut pre_pos = -1i64;
    let mut pending: VecDeque<i64> = VecDeque::new();
    let mut alleles: HashMap<i64, Alleles> = HashMap::new();
    let mut out = Vec::new();
    let mut emit = |center: i64, window: &Window, alleles: &mut HashMap<i64, Alleles>| {
        let Some(Alleles { alt, depth }) = alleles.remove(&center) else { return };
        let Some(tensor) = window.tensor() else { return };
        let flanked = flanked_sequence(&reference, center, reference_start);
        let alt: Vec<String> = alt.iter().map(|(k, v)| format!("{k} {v}")).collect();
        out.push(Candidate {
            tensor,
            position: format!("{contig}:{center}:{}", char::from(flanked[FLANKING])),
            alt_info: format!("{depth}-{}", alt.join(" ")),
        });
    };

    let mut sites: Vec<GvcfSite> = Vec::new();
    let mut any_column = false;
    let region = format!("{contig}:{extend_start}-{extend_end}");
    let mut pileup = bam.pileup(&region, filter)?;
    while let Some(column) = pileup.next_column()? {
        any_column = true;
        let pos = column.pos + 1;
        let reference_base = usize::try_from(pos - reference_start)
            .ok()
            .and_then(|i| reference.get(i))
            .copied()
            .ok_or_else(|| Error::Exit(format!("IndexError: no reference base at {contig}:{pos}")))?;
        let within = pos >= ctg_start && pos < ctg_end;
        if gvcf && within && column.entries().len() == 0 {
            // A position of the chunk without reads is a site, and nothing else: it neither
            // clears nor fills the window. Outside the chunk, its `*` is read as a deletion.
            sites.push(GvcfSite { pos, reference: reference_base, n_ref: 0, n_total: 0 });
            continue;
        }
        let bases = read_bases(&column);
        // A new region clears the window.
        if pre_pos + 1 != pos {
            window = Window::new(head_and_tail);
            pending.clear();
        }
        pre_pos = pos;
        let site = generate_tensor(pos, bases.as_bytes(), &reference, reference_start, reference_base, &thresholds)?;
        let pass_af = site.pass_af || (site.depth > 0 && (job.min_snp_af == 0.0 || job.min_indel_af == 0.0));
        if gvcf && within {
            let (n_ref, n_total) = site.gvcf_counts(reference_base);
            sites.push(GvcfSite { pos, reference: reference_base, n_ref, n_total });
        }
        let pass_confident = confident.as_ref().is_none_or(|tree| {
            let (start, end) = (pos - 1, pos + site.max_del_length as i64 + 1);
            tree.iter().any(|iv| iv.overlaps(start, end))
        });
        let is_candidate = match &known {
            Some(known) => known.contains(&pos),
            None => pass_confident && b"ACGT".contains(&reference_base) && pass_af && site.depth >= job.min_coverage,
        };
        if is_candidate {
            pending.push_back(pos);
            alleles.insert(pos, Alleles { alt: site.alt, depth: site.depth });
        }
        window.push(Some(site.tensor));
        if pending.front().is_some_and(|&first| pos - first == FLANKING as i64) {
            let center = pending.pop_front().expect("checked non-empty");
            emit(center, &window, &mut alleles);
        }
    }
    if head_and_tail {
        // The last positions are padded with zero columns.
        for pos in pre_pos + 1..=pre_pos + FLANKING as i64 {
            window.push(Some([0; CHANNELS]));
            let center = pos - FLANKING as i64;
            if pending.contains(&center) {
                emit(center, &window, &mut alleles);
            }
        }
    }
    let chunk = gvcf.then(|| GvcfChunk {
        contig: contig.to_owned(),
        start: ctg_start,
        end: ctg_end,
        sites: if any_column { GvcfSites::Sites(sites) } else { GvcfSites::Empty },
    });
    Ok((out, chunk))
}

/// A candidate's alleles (`all_alt_dict[pos]`) and depth (`depth_dict[pos]`).
struct Alleles {
    alt: Vec<(String, i64)>,
    depth: i64,
}

/// The ring of the last 33 columns (`tensor` and `pos_offset`); `None` is an empty slot.
struct Window {
    rows: [Option<[i32; CHANNELS]>; WINDOW],
    offset: usize,
}

impl Window {
    /// Empty slots, or zero columns with `--enable_variant_calling_at_sequence_head_and_tail`.
    fn new(head_and_tail: bool) -> Self {
        let rows = if head_and_tail { [Some([0; CHANNELS]); WINDOW] } else { [None; WINDOW] };
        Self { rows, offset: 0 }
    }

    fn push(&mut self, row: Option<[i32; CHANNELS]>) {
        self.rows[self.offset] = row;
        self.offset = (self.offset + 1) % WINDOW;
    }

    /// The 33 columns oldest first, unless a slot is empty.
    fn tensor(&self) -> Option<Vec<i32>> {
        let ordered = self.rows[self.offset..].iter().chain(&self.rows[..self.offset]);
        let mut tensor = Vec::with_capacity(WINDOW * CHANNELS);
        for row in ordered {
            tensor.extend_from_slice(row.as_ref()?);
        }
        Some(tensor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_follow_the_python_parser() {
        let column = b"^!A$.,c+2AG-3nnnT>*#";
        let got = tokens(column).expect("valid column");
        let want = [
            Token::Base(b'A'),
            Token::Base(b'c'),
            Token::Indel(b'+', b"AG"),
            Token::Indel(b'-', b"nnn"),
            Token::Base(b'T'),
            Token::Base(b'*'),
            Token::Base(b'#'),
        ];
        assert_eq!(got, want);
        assert!(tokens(b"A+").is_err());
    }

    #[test]
    fn counts_keep_first_occurrence_order_and_stable_ties() {
        let counts = counter(*b"CACGA");
        assert_eq!(counts, [(b'C', 2), (b'A', 2), (b'G', 1)]);
        assert_eq!(sorted_by_count(&[(b'G', 1), (b'C', 2), (b'A', 2)]), [(b'C', 2), (b'A', 2), (b'G', 1)]);
    }

    #[test]
    fn reference_channel_is_minus_the_base_sum() {
        let thresholds = Thresholds {
            min_snp_af: 0.08,
            min_indel_af: 0.08,
            alleles: CandidateAlleles::SnpsAndIndels,
            fast_mode: FastMode::Off,
        };
        let site = generate_tensor(1, b"AAaCc+1T-1n", b"AAC", 1, b'A', &thresholds).expect("valid column");
        assert_eq!(site.tensor[0], -3);
        assert_eq!(site.tensor[9], -2);
        assert_eq!(site.tensor[REV_DEL], 1);
        assert_eq!(site.tensor[INS], 1);
        let alt: Vec<String> = site.alt.iter().map(|(k, v)| format!("{k} {v}")).collect();
        assert_eq!(alt, ["XC 2", "IAT 1", "DA 1", "RA 1"]);
        assert!(site.pass_af);
    }

    #[test]
    fn gvcf_counts_use_the_reference_base_as_read() {
        let thresholds = Thresholds {
            min_snp_af: 0.08,
            min_indel_af: 0.08,
            alleles: CandidateAlleles::SnpsAndIndels,
            fast_mode: FastMode::Off,
        };
        // Deletion bases (`*`) are not counted; indels are.
        let site = generate_tensor(1, b"AAaC*+1T-1n", b"AAC", 1, b'A', &thresholds).expect("valid column");
        assert_eq!(site.gvcf_counts(b'A'), (3, 6));
        // An `N` becomes `A` for the tensor, but no read supports the `N` itself.
        let site = generate_tensor(1, b"AAaC", b"NAC", 1, b'N', &thresholds).expect("valid column");
        assert_eq!(site.gvcf_counts(b'N'), (0, 4));
    }
}

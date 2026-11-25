//! Port of `CreateTensorPileup` (`preprocess/CreateTensorPileupFromCffi.py`) for calling.
//!
//! Turns the pileup counts of a chunk into one 33 x 18 tensor per candidate, plus the position
//! and alt-info strings the original passes on to variant output.

use std::path::Path;

use crate::Error;
use crate::bed::{Bed, Interval};
use crate::counts::{CandidateParams, FEATURES, FlankRequirement, PileupCounts, pileup_counts};
use crate::gvcf::ChunkSites;
use crate::hts::{AlignmentFile, Fasta};

/// Flanking positions on each side of a candidate (`param.flankingBaseNum`).
pub const FLANKING_BASES: i64 = 16;
/// Positions in a candidate tensor (`no_of_positions`).
pub const POSITIONS: usize = 2 * FLANKING_BASES as usize + 1;
/// Values in a candidate tensor.
pub const TENSOR_LEN: usize = POSITIONS * FEATURES;
/// Longest stretch of a chunk counted at once by [`candidate_tensors`].
const PIECE_LENGTH: i64 = 1_000_000;

/// Depth above which the original scales a tensor down: `max_depth * 1.5` with `max_depth = 144`.
const MAX_DEPTH: f64 = 144.0;

/// Which part of a contig to process.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ChunkRegion {
    /// `--ctgStart`/`--ctgEnd`: 1-based, as the original receives them.
    Range {
        /// First position (1-based).
        start: i64,
        /// End position as passed to the original.
        end: i64,
    },
    /// Chunk `id` (0-based) of `count` equal chunks of the whole contig, by its `.fai` length
    /// (`--chunk_id`/`--chunk_num` without a usable `--extend_bed`).
    FaiChunk {
        /// 0-based chunk number (the original's `--chunk_id` minus 1).
        id: i64,
        /// Number of chunks (`--chunk_num`).
        count: i64,
    },
    /// Chunk `id` (0-based) of `count` equal chunks of the region spanned by the
    /// `--extend_bed` intervals, which also filter candidates.
    BedChunk {
        /// 0-based chunk number (the original's `--chunk_id` minus 1).
        id: i64,
        /// Number of chunks (`--chunk_num`).
        count: i64,
    },
    /// `--vcf_fn`: one chunk of the contig's known sites (see [`known_sites_chunk`]). Only
    /// candidates at these 1-based positions are kept, and the chunk spans the first to the last.
    KnownSites {
        /// The chunk's sites, sorted.
        sites: Vec<i64>,
    },
}

/// Chunk `id` (0-based) of `count` of a contig's sorted, distinct known sites, as
/// `CreateTensorPileup` splits `vcf_candidates_from(vcf_fn, contig)`.
pub fn known_sites_chunk(sites: &[i64], id: i64, count: i64) -> &[i64] {
    let total = sites.len() as i64;
    let size = if total % count == 0 { total / count } else { total / count + 1 };
    let start = (id * size).clamp(0, total) as usize;
    let end = (id * size + size).clamp(0, total) as usize;
    &sites[start..end]
}

/// `chunk_size = span // count + 1 if span % count else span // count`.
fn chunk_size(span: i64, count: i64) -> i64 {
    if span % count != 0 { span / count + 1 } else { span / count }
}

/// One candidate tensor with its labels.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Candidate {
    /// `TENSOR_LEN` values, position-major, as `np.int32`.
    pub tensor: Vec<i32>,
    /// `"{contig}:{pos}:{ref_base}"`.
    pub position: String,
    /// `"{depth}-{alleles}"`.
    pub alt_info: String,
}

/// Builds the candidate tensors of one chunk, like `CreateTensorPileup(args)` in calling mode.
///
/// `extend_bed` is the original's `--extend_bed` (required for [`ChunkRegion::BedChunk`]).
///
/// # Errors
///
/// Fails if a file cannot be read or the chunk does not fit the inputs.
pub fn candidate_tensors(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    contig: &str,
    chunk: &ChunkRegion,
    extend_bed: Option<&Path>,
    params: &CandidateParams,
) -> Result<Vec<Candidate>, Error> {
    if matches!(chunk, ChunkRegion::KnownSites { .. }) {
        return Ok(chunk_pileup(bam, fasta, contig, chunk, extend_bed, params)?
            .map(|p| p.candidates)
            .unwrap_or_default());
    }
    let bed = extend_bed.map(|path| Bed::read_contig(path, contig)).transpose()?;
    let Some((start, end)) = chunk_bounds(fasta, contig, chunk, bed.as_ref())? else { return Ok(Vec::new()) };
    let start = start.max(1);
    // The chunk is counted in pieces, each with its own flanks, for a fraction of the memory. A
    // candidate's tensor and its flank check only see positions within 16 of it, so each piece
    // keeps the candidates in its own range; the first and the last also keep those in the
    // chunk's flanks, which the chunk reports as well.
    let mut candidates = Vec::new();
    let mut piece = start;
    loop {
        let piece_end = if end - piece > PIECE_LENGTH { piece + PIECE_LENGTH } else { end };
        let (first, last) = (piece == start, piece_end >= end);
        let found = range_pileup(bam, fasta, contig, (piece, piece_end), bed.as_ref(), None, params)?.candidates;
        candidates.extend(found.into_iter().filter(|c| {
            let pos = known_site(c).expect("candidate positions are written as integers");
            (first || pos >= piece) && (last || pos < piece_end)
        }));
        if last {
            return Ok(candidates);
        }
        piece = piece_end;
    }
}

/// [`candidate_tensors`] with `--gvcf`: also returns the chunk's positions for the non-variant
/// gVCF, unless a [`ChunkRegion::KnownSites`] chunk has no sites (the original then writes none).
///
/// # Errors
///
/// Fails if a file cannot be read or the chunk does not fit the inputs.
pub fn candidate_tensors_with_sites(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    contig: &str,
    chunk: &ChunkRegion,
    extend_bed: Option<&Path>,
    params: &CandidateParams,
) -> Result<(Vec<Candidate>, Option<ChunkSites>), Error> {
    let Some(mut pileup) = chunk_pileup(bam, fasta, contig, chunk, extend_bed, params)? else {
        return Ok((Vec::new(), None));
    };
    let sites = ChunkSites::new(&mut pileup.counts, fasta, contig, pileup.start, pileup.end, pileup.extend_start)?;
    Ok((pileup.candidates, Some(sites)))
}

/// One chunk's candidates and the pileup they came from.
struct ChunkPileup {
    candidates: Vec<Candidate>,
    counts: PileupCounts,
    /// `ctg_start` (at least 1) and `ctg_end`.
    start: i64,
    end: i64,
    extend_start: i64,
}

fn chunk_pileup(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    contig: &str,
    chunk: &ChunkRegion,
    extend_bed: Option<&Path>,
    params: &CandidateParams,
) -> Result<Option<ChunkPileup>, Error> {
    let bed = extend_bed.map(|path| Bed::read_contig(path, contig)).transpose()?;
    let Some(bounds) = chunk_bounds(fasta, contig, chunk, bed.as_ref())? else { return Ok(None) };
    let sites = match chunk {
        ChunkRegion::KnownSites { sites } => Some(sites.as_slice()),
        _ => None,
    };
    Ok(Some(range_pileup(bam, fasta, contig, bounds, bed.as_ref(), sites, params)?))
}

/// `ctg_start` and `ctg_end` of a chunk, or `None` for a known-sites chunk without sites.
fn chunk_bounds(
    fasta: &Fasta,
    contig: &str,
    chunk: &ChunkRegion,
    bed: Option<&Bed>,
) -> Result<Option<(i64, i64)>, Error> {
    Ok(Some(match chunk {
        ChunkRegion::Range { start, end } => (*start, *end),
        ChunkRegion::FaiChunk { id, count } => {
            let size = chunk_size(fasta.contig_len(contig)?, *count);
            (size * id, size * id + size)
        }
        ChunkRegion::BedChunk { id, count } => {
            let bed = bed.ok_or(Error::MissingBed)?;
            let (bed_start, bed_end) = bed.span().ok_or_else(|| Error::EmptyBed { contig: contig.to_owned() })?;
            let size = chunk_size(bed_end - bed_start, *count);
            let start = bed_start + 1 + size * id;
            (start, start + size)
        }
        ChunkRegion::KnownSites { sites } => match (sites.first(), sites.last()) {
            (Some(&first), Some(&last)) => (first, last),
            // `return [], [], []`: no tensors and no gVCF for the chunk.
            _ => return Ok(None),
        },
    }))
}

/// The candidates of positions `ctg_start..ctg_end`, from the pileup of the range and its
/// flanks; with `sites`, only those at these 1-based positions.
fn range_pileup(
    bam: &mut AlignmentFile,
    fasta: &Fasta,
    contig: &str,
    (ctg_start, ctg_end): (i64, i64),
    bed: Option<&Bed>,
    sites: Option<&[i64]>,
    params: &CandidateParams,
) -> Result<ChunkPileup, Error> {
    let ctg_start = ctg_start.max(1);
    let extend_start = (ctg_start - POSITIONS as i64).max(1);
    let extend_end = ctg_end + POSITIONS as i64;

    // The original's `pileup_counts_*` formats the medaka Region as "{name}:{max(0, start -
    // 1)}-{end}".
    let region = format!("{contig}:{}-{extend_end}", (extend_start - 1).max(0));
    let counts = pileup_counts(bam, fasta, &region, params)?;
    let confident = bed.map(|bed| bed.overlapping(extend_start, extend_end));
    let mut candidates = slice_candidates(contig, &counts, confident.as_deref(), params.flanks);
    if let Some(sites) = sites {
        candidates.retain(|c| known_site(c).is_some_and(|pos| sites.binary_search(&pos).is_ok()));
    }
    Ok(ChunkPileup { candidates, counts, start: ctg_start, end: ctg_end, extend_start })
}

/// The 1-based position in a candidate's `"{contig}:{pos}:{ref}"`.
fn known_site(candidate: &Candidate) -> Option<i64> {
    candidate.position.rsplit(':').nth(1)?.parse().ok()
}

/// A run of rows at consecutive positions (`__enforce_pileup_chunk_contiguity`).
struct Run<'a> {
    rows: &'a [[i32; FEATURES]],
    first: i64,
    last: i64,
}

fn contiguous_runs(counts: &PileupCounts) -> Vec<Run<'_>> {
    let mut runs = Vec::new();
    let mut begin = 0;
    for i in 1..=counts.positions.len() {
        if i == counts.positions.len() || counts.positions[i] - counts.positions[i - 1] > 1 {
            if i > begin {
                runs.push(Run {
                    rows: &counts.rows[begin..i],
                    first: counts.positions[begin],
                    last: counts.positions[i - 1],
                });
            }
            begin = i;
        }
    }
    runs
}

fn slice_candidates(
    contig: &str,
    counts: &PileupCounts,
    confident: Option<&[Interval]>,
    flanks: FlankRequirement,
) -> Vec<Candidate> {
    let runs = contiguous_runs(counts);
    let mut out = Vec::new();
    for alt_info in &counts.candidates {
        // Python: `pos, depth, center_ref_base, alt = alt_info.rstrip().split('-')[:4]`.
        let fields: Vec<&str> = alt_info.trim_end().split('-').collect();
        let [pos, depth, ref_base, alleles, ..] = fields[..] else { continue };
        let pos_1: i64 = pos.parse().expect("the C code writes the position as an integer");
        if let Some(intervals) = confident
            && !intervals.iter().any(|iv| iv.overlaps(pos_1 - 1, pos_1 + 1))
        {
            continue;
        }
        let mut push = |tensor: Vec<i32>| {
            out.push(Candidate {
                tensor,
                position: format!("{contig}:{pos}:{ref_base}"),
                alt_info: format!("{depth}-{alleles}"),
            });
        };
        let (start, end) = (pos_1 - FLANKING_BASES, pos_1 + FLANKING_BASES + 1);
        for run in &runs {
            // `start - 1 >= first and end <= last`.
            if start > run.first && end <= run.last {
                let offset = (start - run.first - 1) as usize;
                let Some(window) = run.rows.get(offset..offset + POSITIONS) else { continue };
                if window.iter().any(|row| row.iter().all(|&v| v == 0)) {
                    continue;
                }
                push(window.iter().flatten().copied().collect());
            }
            if flanks == FlankRequirement::None {
                if let Some(tensor) = head_tensor(run, start, pos_1) {
                    push(tensor);
                }
                if let Some(tensor) = tail_tensor(run, start, end, pos_1) {
                    push(tensor);
                }
            }
        }
    }
    out
}

/// `rows[start:stop]` with Python's slice bounds: negative values count from the end.
fn py_slice(rows: &[[i32; FEATURES]], start: i64, stop: Option<i64>) -> &[[i32; FEATURES]] {
    let len = rows.len() as i64;
    let bound = |i: i64| if i < 0 { (len + i).max(0) } else { i.min(len) } as usize;
    let (start, stop) = (bound(start), stop.map_or(rows.len(), bound));
    rows.get(start..stop).unwrap_or_default()
}

/// The zero rows, then `rows`, flattened, if that makes a whole tensor (`tensor.shape ==
/// (no_of_positions, channel_size)`).
fn padded(before: usize, rows: &[[i32; FEATURES]], after: usize) -> Option<Vec<i32>> {
    if before + rows.len() + after != POSITIONS {
        return None;
    }
    let mut tensor = vec![0; before * FEATURES];
    tensor.extend(rows.iter().flatten());
    tensor.resize(TENSOR_LEN, 0);
    Some(tensor)
}

/// With `--enable_variant_calling_at_sequence_head_and_tail`, a candidate whose window starts
/// before the run: zero rows, then the run's first rows.
fn head_tensor(run: &Run<'_>, start: i64, pos: i64) -> Option<Vec<i32>> {
    if pos < run.first || pos > run.last || start > run.first {
        return None;
    }
    let offset = start - run.first - 1;
    padded((-offset) as usize, py_slice(run.rows, 0, Some(offset + POSITIONS as i64)), 0)
}

/// The same for a window that ends after the run: the run's last rows, then zero rows. The
/// original slices from `start - first - 1`, which counts from the end when negative.
fn tail_tensor(run: &Run<'_>, start: i64, end: i64, pos: i64) -> Option<Vec<i32>> {
    if pos < run.first || pos > run.last || end <= run.last {
        return None;
    }
    let offset = start - run.first - 1;
    let missing = end - run.last - 2;
    if missing > 0 {
        padded(0, py_slice(run.rows, offset, None), missing as usize)
    } else {
        padded(0, py_slice(run.rows, offset, Some(offset + POSITIONS as i64)), 0)
    }
}

/// Scales a tensor in place for very deep candidates, as `CallVariantsFromCffi` does before
/// pileup inference: above 1.5 x 144 reads every value is divided by `depth / 144` and truncated
/// back to `int32`.
pub fn scale_for_depth(tensor: &mut [i32], depth: u64) {
    if depth > 0 && depth as f64 > MAX_DEPTH * 1.5 {
        let scale = depth as f64 / MAX_DEPTH;
        for v in tensor {
            *v = (f64::from(*v) / scale) as i32;
        }
    }
}

/// The depth at the start of a candidate's alt info (`int(alt_info.split('-', 1)[0])`).
pub fn alt_info_depth(alt_info: &str) -> Option<u64> {
    alt_info.split('-').next()?.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(positions: &[i64]) -> PileupCounts {
        PileupCounts {
            rows: positions.iter().map(|&p| [p as i32; FEATURES]).collect(),
            positions: positions.to_vec(),
            ..PileupCounts::default()
        }
    }

    #[test]
    fn runs_split_at_gaps() {
        let c = counts(&[10, 11, 12, 20, 21, 30]);
        let runs: Vec<_> = contiguous_runs(&c).iter().map(|r| (r.first, r.last, r.rows.len())).collect();
        assert_eq!(runs, vec![(10, 12, 3), (20, 21, 2), (30, 30, 1)]);
    }

    #[test]
    fn window_is_centred_on_candidate() {
        let mut c = counts(&(100..200).collect::<Vec<_>>());
        // 1-based position 150 is 0-based row 149.
        c.candidates.push("150-30-A-XC 5 RA 25 ".to_owned());
        let got = slice_candidates("chr1", &c, None, FlankRequirement::Contiguous);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].tensor[..FEATURES], [133; FEATURES]);
        assert_eq!(got[0].tensor[16 * FEATURES], 149);
        assert_eq!((got[0].position.as_str(), got[0].alt_info.as_str()), ("chr1:150:A", "30-XC 5 RA 25"));
    }

    #[test]
    fn sequence_ends_are_padded_with_zeros() {
        let mut c = counts(&(100..140).collect::<Vec<_>>());
        // 1-based 105 is 0-based row 104, 5 rows into the run: 12 zero rows, then 21 rows.
        c.candidates.push("105-30-A-XC 5 ".to_owned());
        assert!(slice_candidates("chr1", &c, None, FlankRequirement::Contiguous).is_empty());
        let got = slice_candidates("chr1", &c, None, FlankRequirement::None);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].tensor[..12 * FEATURES], [0; 12 * FEATURES]);
        assert_eq!(got[0].tensor[12 * FEATURES], 100);
        assert_eq!(got[0].tensor[TENSOR_LEN - 1], 120);
        // 1-based 137 is row 136, near the end (139): 20 rows, then 13 zero rows.
        c.candidates[0] = "137-30-A-XC 5 ".to_owned();
        let got = slice_candidates("chr1", &c, None, FlankRequirement::None);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].tensor[0], 120);
        assert_eq!(got[0].tensor[19 * FEATURES], 139);
        assert_eq!(got[0].tensor[20 * FEATURES..], [0; 13 * FEATURES]);
    }

    #[test]
    fn known_sites_split_into_ceil_sized_chunks() {
        let sites: Vec<i64> = (1..=10).collect();
        assert_eq!(known_sites_chunk(&sites, 0, 3), &[1, 2, 3, 4]);
        assert_eq!(known_sites_chunk(&sites, 2, 3), &[9, 10]);
        assert_eq!(known_sites_chunk(&sites, 1, 5), &[3, 4]);
        assert!(known_sites_chunk(&sites[..2], 4, 5).is_empty());
    }

    #[test]
    fn deep_tensors_are_scaled_and_truncated() {
        let mut t = vec![300, -300, 7];
        scale_for_depth(&mut t, 300);
        assert_eq!(t, vec![144, -144, 3]);
        let mut shallow = vec![5];
        scale_for_depth(&mut shallow, 216);
        assert_eq!(shallow, vec![5]);
    }
}

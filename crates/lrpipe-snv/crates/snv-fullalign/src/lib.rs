//! Full-alignment feature tensors, with the original's in-memory read haplotagging.
//!
//! A port of the original's C full alignment (`calculate_*_full_alignment`, with or without the
//! dwell-time channel) and `CreateTensorFullAlignment`
//! (`preprocess/CreateTensorFullAlignmentFromCffi.py`), version 2.0.3, producing the same `int8`
//! tensors, position strings and alt-info strings. [`FullAlignmentJob::run`] does what one
//! `CallVariantsFromCffi --full_aln_regions` job does before inference; [`full_alignment_tensors`]
//! takes already parsed inputs, so a caller running many jobs can read the phased VCF once.
//!
//! Reproduced on purpose (see the module docs for details):
//!
//! - reads are chosen for deep candidates with glibc `rand()` from its default seed, drawn in
//!   candidate order through the whole job, as each job of the original is its own process;
//! - insertions and deletions appear in the alt info in khash bucket order;
//! - the allele-frequency channel divides `100 * count` by the depth in `float`, because
//!   `normalize_af` does not parenthesise its argument;
//! - bases the original never set (positions inside an `N` CIGAR operation, `'='`) index its
//!   lookup tables out of bounds; the bytes its C library has there are used;
//! - candidates below position 16 get empty matrices (their window start wraps around) but
//!   full alt info;
//! - haplotagging reproduces `cigar_prefix_length`, which returns zero flanks near read ends.

mod features;
mod haplotag;
mod inputs;
mod khash;
mod rand;
mod record;

use std::path::{Path, PathBuf};

pub use inputs::{CandidateRegions, PhasedGenotype, PhasedSnp, PhasedSnps};

/// Flanking positions on each side of a candidate (`flanking_base_num`).
pub const FLANKING_BASES: u64 = 16;
/// Positions in a candidate matrix (`no_of_positions`).
pub const POSITIONS: usize = 2 * FLANKING_BASES as usize + 1;
/// Feature channels per position without dwell time (`channel_size`): reference base,
/// alternative base, strand, mapping quality, base quality, candidate proportion, insertion
/// bases and haplotype.
pub const CHANNELS: usize = 8;
/// Feature channels per position with dwell time (`channel_size_dwell`): [`CHANNELS`], then the
/// signal length of the base.
pub const CHANNELS_WITH_DWELL_TIME: usize = 9;

/// Errors from reading inputs or building features.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A file could not be read.
    #[error("cannot read {path:?}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// A row of the candidate BED is malformed.
    #[error("invalid BED row {line} in {path:?}")]
    Bed {
        /// The BED file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
    },
    /// The candidate BED has no row on the contig (the original then fails to read the BAM).
    #[error("no full-alignment region on {contig:?} in {path:?}")]
    EmptyRegions {
        /// The BED file.
        path: PathBuf,
        /// The contig.
        contig: String,
    },
    /// A record of the phased VCF is malformed.
    #[error("invalid phased VCF record on line {line} of {path:?}")]
    Vcf {
        /// The VCF file.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
    },
    /// The phased VCF positions decrease.
    #[error("phased VCF {path:?} is not sorted at line {line}")]
    UnsortedVcf {
        /// The VCF file (empty for SNPs given directly).
        path: PathBuf,
        /// 1-based line number (or SNP number).
        line: usize,
    },
    /// A position does not fit the C types the original passes it through.
    #[error("position {position} is out of range")]
    Position {
        /// The position.
        position: i64,
    },
    /// The reference could not be opened or the contig fetched.
    #[error("cannot read the reference")]
    Reference {
        /// The underlying error.
        #[source]
        source: snv_pileup::Error,
    },
    /// The alignments could not be opened, queried or decoded.
    #[error("cannot read alignments from {path:?}")]
    Alignments {
        /// The alignment file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
    /// The contig is not in the alignment header.
    #[error("contig {contig:?} is not in the alignment header")]
    Contig {
        /// The contig name.
        contig: String,
    },
}

/// Sequencing platform, which sets the matrix depth (`param.matrix_depth_dict`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum Platform {
    /// Oxford Nanopore: 89 rows.
    #[default]
    Ont,
    /// PacBio HiFi: 55 rows.
    Hifi,
    /// Illumina: 55 rows.
    Ilmn,
}

impl Platform {
    /// Rows of a candidate matrix.
    pub fn matrix_depth(self) -> usize {
        match self {
            Self::Ont => 89,
            Self::Hifi | Self::Ilmn => 55,
        }
    }
}

/// Longest insertion or deletion written into the alt info.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum IndelLimit {
    /// `maximum_variant_length_that_need_infer`: 50 bases.
    #[default]
    Standard,
    /// `--enable_long_indel`: 100,000 bases.
    Long,
}

impl IndelLimit {
    /// The limit in bases.
    pub fn max_length(self) -> u64 {
        match self {
            Self::Standard => 50,
            Self::Long => 100_000,
        }
    }
}

/// Whether the matrices carry each base's signal length from the read's Dorado move table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum DwellTime {
    /// [`CHANNELS`] channels.
    #[default]
    Off,
    /// `--enable_dwell_time`: [`CHANNELS_WITH_DWELL_TIME`] channels, the last one the number of
    /// signal samples of the base (plus those of an insertion after it) from the `mv` tag, as
    /// `int8` (so it wraps above 127). Reads without a usable `mv` tag get 0.
    On,
}

impl DwellTime {
    /// Feature channels per position.
    pub fn channels(self) -> usize {
        match self {
            Self::Off => CHANNELS,
            Self::On => CHANNELS_WITH_DWELL_TIME,
        }
    }
}

/// Feature options of a full-alignment job.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct FeatureOptions {
    /// `--minMQ`: reads below this mapping quality are dropped.
    pub min_mapq: u64,
    /// `--platform`.
    pub platform: Platform,
    /// `--enable_long_indel`.
    pub indel_limit: IndelLimit,
    /// `--enable_dwell_time`.
    pub dwell_time: DwellTime,
    /// Order reads tied at the same start before selecting the full-alignment rows.
    pub canonical_read_order: bool,
}

impl Default for FeatureOptions {
    /// The values the original's pipeline passes: `--minMQ 5`, ONT, no long indels, no dwell time.
    fn default() -> Self {
        Self::new()
    }
}

impl FeatureOptions {
    /// The values the original's pipeline passes: `--minMQ 5`, ONT, no long indels, no dwell time.
    pub fn new() -> Self {
        Self {
            min_mapq: 5,
            platform: Platform::Ont,
            indel_limit: IndelLimit::Standard,
            dwell_time: DwellTime::Off,
            canonical_read_order: false,
        }
    }
}

/// Whether and how reads are haplotagged.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum Haplotagging<'a> {
    /// `--no_phasing_for_fa True`: every read is unphased.
    Off,
    /// Haplotag reads with mapping quality of at least 20 from these SNPs.
    Phased(&'a PhasedSnps),
}

/// The tensors of one job, as `CreateTensorFullAlignment` returns them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FullAlignment {
    matrix_depth: usize,
    channels: usize,
    tensors: Vec<i8>,
    positions: Vec<String>,
    alt_info: Vec<String>,
}

impl FullAlignment {
    /// Number of candidates.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// True if there are no candidates.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Rows of each matrix.
    pub fn matrix_depth(&self) -> usize {
        self.matrix_depth
    }

    /// Feature channels per position ([`DwellTime::channels`]).
    pub fn channels(&self) -> usize {
        self.channels
    }

    /// Shape of [`Self::tensors`]: candidates, rows, positions, channels.
    pub fn shape(&self) -> [usize; 4] {
        [self.len(), self.matrix_depth, POSITIONS, self.channels]
    }

    /// All matrices, C order (`np.int8` of [`Self::shape`]).
    pub fn tensors(&self) -> &[i8] {
        &self.tensors
    }

    /// The matrix of candidate `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is not below [`Self::len`].
    pub fn tensor(&self, i: usize) -> &[i8] {
        assert!(i < self.len(), "candidate {i} of {}", self.len());
        let n = self.matrix_depth * POSITIONS * self.channels;
        &self.tensors[i * n..(i + 1) * n]
    }

    /// `"{contig}:{pos}:{ref_base}"` per candidate (1-based position).
    pub fn positions(&self) -> &[String] {
        &self.positions
    }

    /// `"{depth}-{alleles}"` per candidate.
    pub fn alt_info(&self) -> &[String] {
        &self.alt_info
    }

    /// The tensors, positions and alt info.
    pub fn into_parts(self) -> (Vec<i8>, Vec<String>, Vec<String>) {
        (self.tensors, self.positions, self.alt_info)
    }
}

/// What the first records of an alignment file say about move tables (`mv` tags), as
/// the original's pipeline script (`run_*.py`) samples them in `check_bam_for_valid_mv_tag`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoveTables {
    /// A record has an `mv` tag with a value.
    Found,
    /// No usable tag in the first `checked` records, but some had an empty `mv` array (or
    /// string).
    Empty {
        /// Records read.
        checked: usize,
    },
    /// No `mv` tag in the first `checked` records.
    Missing {
        /// Records read.
        checked: usize,
    },
}

/// Reads up to `max_records` records of `bam` from the start and looks for an `mv` tag with a
/// value.
///
/// # Errors
///
/// Fails if the file cannot be opened or a record cannot be decoded.
pub fn sample_move_tables(
    bam: impl AsRef<Path>,
    reference: impl AsRef<Path>,
    max_records: usize,
) -> Result<MoveTables, Error> {
    use snv_hts::Aux;
    let path = bam.as_ref();
    let alignments = |e: snv_hts::Error| Error::Alignments { path: path.to_owned(), source: Box::new(e) };
    let mut reader = snv_hts::Alignments::open(path, reference).map_err(alignments)?;
    let (mut checked, mut empty) = (0, false);
    for record in reader.first_reads(max_records).map_err(alignments)? {
        checked += 1;
        for (tag, value) in snv_hts::aux_fields(&record) {
            if tag != *b"mv" {
                continue;
            }
            // `_mv_field_has_value`: a B array needs an element, any other type a value.
            let has_value = match value {
                Aux::Ints(a) => !a.is_empty(),
                Aux::Floats(a) => !a.is_empty(),
                Aux::Text(s) | Aux::Hex(s) => !s.trim_ascii().is_empty(),
                _ => true,
            };
            if has_value {
                return Ok(MoveTables::Found);
            }
            empty = true;
        }
    }
    Ok(if empty { MoveTables::Empty { checked } } else { MoveTables::Missing { checked } })
}

/// `bytes.decode('utf8', 'ignore')`.
fn decode_ignore(bytes: &[u8]) -> String {
    bytes.utf8_chunks().map(|chunk| chunk.valid()).collect()
}

/// Builds the full-alignment tensors of `regions` from `bam` (the C full-alignment
/// function plus the Python that splits its alt-info strings).
///
/// # Errors
///
/// Fails if a file cannot be read, the contig is missing, or a position does not fit a C `int`.
pub fn full_alignment_tensors(
    bam: impl AsRef<Path>,
    reference: impl AsRef<Path>,
    regions: &CandidateRegions,
    haplotagging: Haplotagging<'_>,
    options: &FeatureOptions,
) -> Result<FullAlignment, Error> {
    let job = feature_job(bam.as_ref(), reference.as_ref(), regions, haplotagging, options)?;
    let mut tensors = Vec::new();
    let (positions, alt_info) =
        full_alignment_batched(&job, regions.contig(), usize::MAX, &mut |batch| tensors.extend_from_slice(batch))?;
    let (matrix_depth, channels) = (job.matrix_depth, options.dwell_time.channels());
    Ok(FullAlignment { matrix_depth, channels, tensors, positions, alt_info })
}

/// [`full_alignment_tensors`], handing the matrices to `sink` `batch` candidates at a time (the
/// last batch may be shorter) as soon as each batch is built, so that a job never holds all of
/// them: the candidates' positions and alt info.
///
/// # Errors
///
/// Fails if a file cannot be read, the contig is missing, or a position does not fit a C `int`.
pub fn full_alignment_tensors_batched(
    bam: impl AsRef<Path>,
    reference: impl AsRef<Path>,
    regions: &CandidateRegions,
    haplotagging: Haplotagging<'_>,
    options: &FeatureOptions,
    batch: usize,
    sink: &mut dyn FnMut(&[i8]),
) -> Result<(Vec<String>, Vec<String>), Error> {
    let job = feature_job(bam.as_ref(), reference.as_ref(), regions, haplotagging, options)?;
    full_alignment_batched(&job, regions.contig(), batch, sink)
}

/// The feature builder's job for `regions`.
fn feature_job<'a>(
    bam: &'a Path,
    reference: &'a Path,
    regions: &'a CandidateRegions,
    haplotagging: Haplotagging<'a>,
    options: &FeatureOptions,
) -> Result<features::Job<'a>, Error> {
    let to_c_int = |v: i64| i32::try_from(v).map(i64::from).map_err(|_| Error::Position { position: v });
    // `hts_parse_reg("{ctg}:{start}-{end}")`: 0-based start, exclusive end.
    let beg = to_c_int(regions.start() - 1)?.max(0);
    let end = to_c_int(regions.end())?;
    let matrix_depth = options.platform.matrix_depth();
    Ok(features::Job {
        bam,
        reference,
        contig: regions.contig(),
        beg,
        end,
        candidates: regions.centers(),
        snps: match haplotagging {
            Haplotagging::Off => None,
            Haplotagging::Phased(snps) => Some(snps.as_slice()),
        },
        min_mapq: options.min_mapq,
        matrix_depth,
        max_indel_length: options.indel_limit.max_length(),
        dwell_time: options.dwell_time,
        canonical_read_order: options.canonical_read_order,
    })
}

/// Builds a job's matrices, handing them to `sink` `batch` candidates at a time: the
/// candidates' positions and alt info.
fn full_alignment_batched(
    job: &features::Job<'_>,
    contig: &str,
    batch: usize,
    sink: &mut dyn FnMut(&[i8]),
) -> Result<(Vec<String>, Vec<String>), Error> {
    let data = features::calculate_full_alignment(job, batch, sink)?;
    let mut positions = Vec::with_capacity(data.alt_info.len());
    let mut alt_info = Vec::with_capacity(data.alt_info.len());
    for raw in &data.alt_info {
        // `pos, depth, center_ref_base, alt = alt_info_string.rstrip().split('-')[:4]`.
        let text = decode_ignore(raw);
        let mut parts = inputs::py_rstrip(&text).split('-');
        let mut part = || parts.next().unwrap_or_default();
        let (pos, depth, center_ref_base, alt) = (part(), part(), part(), part());
        positions.push(format!("{contig}:{pos}:{center_ref_base}"));
        alt_info.push(format!("{depth}-{alt}"));
    }
    Ok((positions, alt_info))
}

/// One full-alignment job with the arguments `CallVariantsFromCffi` passes on to
/// `CreateTensorFullAlignment`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FullAlignmentJob {
    bam: PathBuf,
    reference: PathBuf,
    contig: String,
    full_aln_regions: PathBuf,
    phased_vcf: Option<PathBuf>,
    options: FeatureOptions,
}

impl FullAlignmentJob {
    /// A job over the candidates in `full_aln_regions` on `contig`, without haplotagging and
    /// with [`FeatureOptions::new`].
    pub fn new(
        bam: impl Into<PathBuf>,
        reference: impl Into<PathBuf>,
        contig: impl Into<String>,
        full_aln_regions: impl Into<PathBuf>,
    ) -> Self {
        Self {
            bam: bam.into(),
            reference: reference.into(),
            contig: contig.into(),
            full_aln_regions: full_aln_regions.into(),
            phased_vcf: None,
            options: FeatureOptions::new(),
        }
    }

    /// Haplotags reads from this phased VCF (`--phased_vcf_fn` with `--no_phasing_for_fa
    /// False`). As in the original, a path that does not exist leaves every read unphased.
    pub fn phased_vcf(mut self, path: impl Into<PathBuf>) -> Self {
        self.phased_vcf = Some(path.into());
        self
    }

    /// Sets the feature options.
    pub fn options(mut self, options: FeatureOptions) -> Self {
        self.options = options;
        self
    }

    /// Reads the inputs and builds the tensors, as `CreateTensorFullAlignment(args)` does with
    /// `tensor_can_fn` `PIPE`.
    ///
    /// # Errors
    ///
    /// Fails as [`CandidateRegions::from_bed`], [`PhasedSnps::from_vcf`] and
    /// [`full_alignment_tensors`] do.
    pub fn run(&self) -> Result<FullAlignment, Error> {
        let regions = CandidateRegions::from_bed(&self.full_aln_regions, &self.contig)?;
        let snps = match &self.phased_vcf {
            Some(path) if path.exists() => PhasedSnps::from_vcf(path, &self.contig)?,
            _ => PhasedSnps::default(),
        };
        let haplotagging = if self.phased_vcf.is_some() { Haplotagging::Phased(&snps) } else { Haplotagging::Off };
        full_alignment_tensors(&self.bam, &self.reference, &regions, haplotagging, &self.options)
    }

    /// [`Self::run`], handing the matrices to `sink` `batch` candidates at a time (see
    /// [`full_alignment_tensors_batched`]): the candidates' positions and alt info.
    ///
    /// # Errors
    ///
    /// As [`Self::run`].
    pub fn run_batched(&self, batch: usize, sink: &mut dyn FnMut(&[i8])) -> Result<(Vec<String>, Vec<String>), Error> {
        let regions = CandidateRegions::from_bed(&self.full_aln_regions, &self.contig)?;
        let snps = match &self.phased_vcf {
            Some(path) if path.exists() => PhasedSnps::from_vcf(path, &self.contig)?,
            _ => PhasedSnps::default(),
        };
        let haplotagging = if self.phased_vcf.is_some() { Haplotagging::Phased(&snps) } else { Haplotagging::Off };
        full_alignment_tensors_batched(&self.bam, &self.reference, &regions, haplotagging, &self.options, batch, sink)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decoding_drops_invalid_utf8() {
        assert_eq!(decode_ignore(b"12-3-A-XC 2 "), "12-3-A-XC 2 ");
        assert_eq!(decode_ignore(b"1-2-\xff-"), "1-2--");
    }

    #[test]
    fn platform_depths() {
        assert_eq!(Platform::Ont.matrix_depth(), 89);
        assert_eq!(Platform::Hifi.matrix_depth(), 55);
        assert_eq!(IndelLimit::Long.max_length(), 100_000);
        assert_eq!(DwellTime::On.channels(), 9);
    }
}

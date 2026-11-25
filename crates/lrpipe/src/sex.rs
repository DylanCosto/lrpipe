//! The sample's sex chromosomes from the reads' index alone, before anything reads the
//! alignments: the read pass needs them for chrX and chrY ploidy as it genotypes repeats.
//!
//! A BAM index (BAI or CSI) counts each reference's mapped reads; a CRAM index lists each
//! slice's reference and size in bytes. Either, over the reference's length from the file's
//! header, gives a density per contig. chrX and chrY are measured against the median
//! autosome: chrX is about 1 for XX and 0.5 for XY. chrY is low for XX but not 0, reads of
//! other chromosomes mapping to its repeats (0.083 by BAI read counts on HG004, a woman), and
//! for XY well above that (0.49 by CRAM slice bytes and 1.4 by BAI read counts on HG002, her
//! son, whose secondary and supplementary alignments pile into chrY's repeats);
//! [`Y_THRESHOLD`] lies between, about as far from each. A sample fitting neither (45,X,
//! 47,XXY) is left unclear. A sample in several files is measured over all of them
//! ([`estimate_all`]).
//!
//! The index also gives the unplaced unmapped reads, which nothing else reads, for the QC
//! summary's fraction of the reads mapped ([`unplaced_unmapped`]).

use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use flate2::read::{GzDecoder, MultiGzDecoder};

use crate::plan::Sex;

/// chrX against the autosomes below which a sample has one X.
pub const ONE_X_BELOW: f64 = 0.75;
/// chrY against the autosomes from which a sample has a Y: 2.4 times HG004's (XX, BAI) and
/// 2.4 times below HG002's (XY, CRAM).
pub const Y_THRESHOLD: f64 = 0.2;
/// A sex chromosome's depth against the autosomes' above which it fits neither XX nor XY: at
/// most one copy per autosome pair's two, so about 1 at most (an XXY or XXX sample's chrX is 1.5).
const MOST_PLAUSIBLE: f64 = 1.6;
/// An autosome's depth against the median autosome's above which the reads are not a whole
/// genome's: a trisomy is 1.5.
const MOST_EVEN: f64 = 3.0;

/// Why the index could not be read.
#[derive(Debug, thiserror::Error)]
pub enum SexError {
    /// A file could not be read.
    #[error("reading {}", path.display())]
    Read {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// The reads have no index lrpipe reads.
    #[error("{} has no .bai, .csi or .crai index", reads.display())]
    NoIndex {
        /// The reads.
        reads: PathBuf,
    },
    /// A file is not what its name says.
    #[error("{}: {reason}", path.display())]
    Format {
        /// The file.
        path: PathBuf,
        /// What is wrong.
        reason: String,
    },
    /// Most autosomes have no reads: the reads are of some contigs (a chr20 BAM), not a whole
    /// genome, whose autosomes the sex chromosomes are measured against.
    #[error(
        "{} has reads on {} of the 22 autosomes{}: the sex is estimated from a whole genome's reads, against \
         its autosomes",
        reads.display(),
        with.len(),
        listed(with)
    )]
    Partial {
        /// The reads.
        reads: PathBuf,
        /// The autosomes with reads.
        with: Vec<String>,
    },
    /// One autosome has far more reads than the others: the reads are a panel's or adaptive
    /// sampling's targets, with a few reads mismapped across the rest, not a whole genome.
    #[error(
        "{} has {contig} at {ratio:.1} times the median autosome's depth: the sex is estimated from a whole \
         genome's reads, whose autosomes are about evenly covered (a panel or adaptive sampling? give --sex)",
        reads.display()
    )]
    Uneven {
        /// The reads.
        reads: PathBuf,
        /// The deepest autosome.
        contig: String,
        /// Its depth against the median's.
        ratio: f64,
    },
    /// The reference lacks most autosomes, or chrX: reads aligned to it are not a whole genome's,
    /// which is known before they are aligned.
    #[error(
        "the reference {} has {} of the 22 autosomes{}{}: the sex is estimated from a whole genome's reads, \
         against its autosomes",
        reference.display(),
        autosomes.len(),
        listed(autosomes),
        if *x { "" } else { " and no chrX" }
    )]
    Reference {
        /// The reference.
        reference: PathBuf,
        /// Its autosomes.
        autosomes: Vec<String>,
        /// Whether it has chrX.
        x: bool,
    },
    /// The index has no chrX.
    #[error("the index of {} has no {what}", reads.display())]
    Missing {
        /// The reads.
        reads: PathBuf,
        /// What it lacks.
        what: &'static str,
    },
}

/// ` (chr1, chr2)`, or nothing for no names.
fn listed(names: &[String]) -> String {
    if names.is_empty() { String::new() } else { format!(" ({})", names.join(", ")) }
}

/// Autosomes of the 22 with reads that a sample needs for its sex to be estimated: below half,
/// the median autosome has none.
const AUTOSOMES_NEEDED: usize = 11;

/// The prefix of the human contigs' names among `names`: `chr` (`chr1`, `chrX`), or none (`1`,
/// `X`, as GRCh37 names them).
fn prefix<'a>(mut names: impl Iterator<Item = &'a str>) -> &'static str {
    let human = |n: &str| n == "X" || n == "Y" || n.parse::<u8>().is_ok_and(|i| (1..=22).contains(&i));
    if names.any(|n| n.strip_prefix("chr").is_some_and(human)) { "chr" } else { "" }
}

/// Checks, from the reference's `.fai`, that reads aligned to it can give the sex: it has chrX and
/// at least [`AUTOSOMES_NEEDED`] of the 22 autosomes. For unaligned reads, whose index is not there
/// until they are aligned: a reference of some contigs (a chr20 window) is refused before then.
///
/// # Errors
///
/// [`SexError::Reference`] if it has too few autosomes or no chrX, or if the `.fai` cannot be read.
pub fn check_reference(reference: &Path) -> Result<(), SexError> {
    let fai = crate::plan::Layout::index(reference, "fai");
    let text = std::fs::read_to_string(&fai).map_err(|source| SexError::Read { path: fai, source })?;
    let names: Vec<&str> = text.lines().filter_map(|l| l.split('\t').next()).collect();
    let prefix = prefix(names.iter().copied());
    let autosomes: Vec<String> =
        (1..=22).map(|i| format!("{prefix}{i}")).filter(|n| names.contains(&n.as_str())).collect();
    let x = names.contains(&format!("{prefix}X").as_str());
    if autosomes.len() < AUTOSOMES_NEEDED || !x {
        return Err(SexError::Reference { reference: reference.to_path_buf(), autosomes, x });
    }
    Ok(())
}

/// The sex chromosomes' density against the autosomes', and what it says.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Estimate {
    /// chrX's.
    pub x: f64,
    /// chrY's (0 when the reference has no chrY).
    pub y: f64,
    /// The sex, when the two agree: one X and a Y, or not one X and no Y.
    pub sex: Option<Sex>,
}

impl Estimate {
    fn new(x: f64, y: f64) -> Self {
        let sex = if x > MOST_PLAUSIBLE || y > MOST_PLAUSIBLE {
            None
        } else if x < ONE_X_BELOW && y >= Y_THRESHOLD {
            Some(Sex::Male)
        } else if x >= ONE_X_BELOW && y < Y_THRESHOLD {
            Some(Sex::Female)
        } else {
            None
        };
        Self { x, y, sex }
    }
}

/// Estimates the sex of the sample in `reads` (BAM or CRAM) from its index.
///
/// # Errors
///
/// If the reads have no index, or it or the header cannot be read, or they have no chrX or
/// reads on too few autosomes.
pub fn estimate(reads: &Path) -> Result<Estimate, SexError> {
    estimate_all(std::slice::from_ref(&reads.to_path_buf()))
}

/// Estimates the sex of the sample whose reads are in the files `reads` (BAM or CRAM, with the
/// same reference sequences) from their indexes: each contig's reads (BAM) or bytes (CRAM)
/// added up over the files. BAM and CRAM indexes count in other units, so with both the BAMs'
/// densities and the CRAMs' are each measured against their own autosomes, and the two averaged.
///
/// # Errors
///
/// If the reads have no index, or it or the header cannot be read, or they have no chrX or
/// reads on too few autosomes.
pub fn estimate_all(reads: &[PathBuf]) -> Result<Estimate, SexError> {
    // Per unit (BAM's reads, CRAM's bytes): whether CRAM, the first file, and its amounts.
    let mut units: Vec<(bool, &Path, Amounts)> = Vec::new();
    for path in reads {
        let cram = crate::plan::is_cram(path);
        let amounts = index_amounts(path)?;
        match units.iter_mut().find(|u| u.0 == cram) {
            Some((_, _, sums)) => {
                for ((name, _), amount) in amounts.contigs.iter().zip(&amounts.amounts) {
                    if let Some(k) = sums.contigs.iter().position(|(n, _)| n == name) {
                        sums.amounts[k] += amount;
                    }
                }
            }
            None => units.push((cram, path, amounts)),
        }
    }
    let mut ratios = Vec::with_capacity(units.len());
    for (_, path, amounts) in &units {
        ratios.push(ratios_against_autosomes(path, &amounts.contigs, &amounts.amounts)?);
    }
    let n = ratios.len().max(1) as f64;
    let (x, y) = ratios.iter().fold((0.0, 0.0), |(x, y), r| (x + r.0 / n, y + r.1 / n));
    Ok(Estimate::new(x, y))
}

/// A file's contigs and the amount of reads on each by its index.
struct Amounts {
    /// Each contig's name and length, from the header.
    contigs: Vec<(String, u64)>,
    /// Each one's reads mapped (BAI, CSI) or bytes (CRAI).
    amounts: Vec<f64>,
}

/// The contigs of `reads`' header, and the amount of reads on each by its index.
fn index_amounts(reads: &Path) -> Result<Amounts, SexError> {
    let name = reads.as_os_str().to_string_lossy();
    let beside = |ext: &str| {
        [PathBuf::from(format!("{name}{ext}")), reads.with_extension(&ext[1..])].into_iter().find(|p| p.exists())
    };
    if crate::plan::is_cram(reads) {
        let index = beside(".crai").ok_or_else(|| SexError::NoIndex { reads: reads.into() })?;
        let contigs = cram_contigs(reads)?;
        let amounts = crai_bytes(&index, contigs.len())?;
        Ok(Amounts { contigs, amounts })
    } else {
        let contigs = bam_contigs(reads)?;
        let (amounts, _) = bam_index_counts(reads, contigs.len(), &beside)?;
        Ok(Amounts { contigs, amounts })
    }
}

/// The mapped reads per reference of a BAM's index (BAI or CSI), and its unplaced unmapped
/// reads when it counts them; `beside` finds the index by its extension.
fn bam_index_counts(
    reads: &Path,
    n_contigs: usize,
    beside: &dyn Fn(&str) -> Option<PathBuf>,
) -> Result<(Vec<f64>, Option<u64>), SexError> {
    if let Some(bai) = beside(".bai") {
        bai_mapped(&bai, n_contigs)
    } else if let Some(csi) = beside(".csi") {
        csi_mapped(&csi, n_contigs)
    } else {
        Err(SexError::NoIndex { reads: reads.into() })
    }
}

/// The unplaced unmapped reads of `reads` (past every contig, with no position), from its
/// index: a BAI's or CSI's count of them (none when it has none), or for CRAM the records of
/// the containers its index lists for unmapped reads alone, from their headers, and of a
/// container shared with mapped reads (htslib's first one of them), its unplaced records,
/// decoded (their flags and reference alone, which take no reference to decode).
///
/// # Errors
///
/// If the reads have no index, or it, the header or a container cannot be read.
pub fn unplaced_unmapped(reads: &Path) -> Result<Option<u64>, SexError> {
    let name = reads.as_os_str().to_string_lossy();
    let beside = |ext: &str| {
        [PathBuf::from(format!("{name}{ext}")), reads.with_extension(&ext[1..])].into_iter().find(|p| p.exists())
    };
    if !crate::plan::is_cram(reads) {
        let contigs = bam_contigs(reads)?;
        return Ok(bam_index_counts(reads, contigs.len(), &beside)?.1);
    }
    let index = beside(".crai").ok_or_else(|| SexError::NoIndex { reads: reads.into() })?;
    let mut text = String::new();
    MultiGzDecoder::new(BufReader::new(open(&index)?))
        .read_to_string(&mut text)
        .map_err(|source| SexError::Read { path: index.clone(), source })?;
    // Each container once, with whether it holds mapped reads too: a container has a line for
    // each reference of each of its slices.
    let lines: Vec<(bool, u64)> = text
        .lines()
        .filter_map(|line| {
            let mut f = line.split('\t');
            let mapped = f.next()? != "-1";
            Some((mapped, f.nth(2)?.parse().ok()?))
        })
        .collect();
    let mut containers: std::collections::BTreeMap<u64, bool> =
        lines.iter().filter(|(mapped, _)| !mapped).map(|&(_, offset)| (offset, false)).collect();
    for (_, offset) in lines.iter().filter(|(mapped, _)| *mapped) {
        if let Some(shared) = containers.get_mut(offset) {
            *shared = true;
        }
    }
    let read = |source| SexError::Read { path: reads.into(), source };
    let mut file = BufReader::new(open(reads)?);
    let mut unmapped = 0;
    for (&offset, &shared) in &containers {
        if shared {
            unmapped += shared_container_unplaced(reads, offset).map_err(read)?;
            continue;
        }
        file.seek(SeekFrom::Start(offset)).map_err(read)?;
        // Container header: length, reference, start, span, records.
        let mut r = Bytes::new(&mut file, reads);
        r.i32()?;
        r.itf8()?;
        r.itf8()?;
        r.itf8()?;
        let records = r.itf8()?;
        unmapped += u64::try_from(records).map_err(|_| r.format(format!("{records} records in a container")))?;
    }
    Ok(Some(unmapped))
}

/// The unplaced records of the CRAM container at `offset` of `reads`, by decoding its records'
/// flags and references.
fn shared_container_unplaced(reads: &Path, offset: u64) -> io::Result<u64> {
    use lrpipe_reads::{Fields, MdNm, Tags};
    use pasta::sam::alignment::Record as _;
    let mut reader = pasta::cram::io::Reader::new(File::open(reads)?);
    let header = reader.read_header()?;
    reader.seek(SeekFrom::Start(offset))?;
    let mut container = pasta::cram::io::reader::Container::default();
    if reader.read_container(&mut container)? == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, format!("no container at byte {offset}")));
    }
    let compression = container.compression_header()?;
    let (fields, tags) = (Fields::FLAGS | Fields::REFERENCE_SEQUENCE_ID, Tags::Only(Vec::new()));
    let mut unplaced = 0;
    for slice in container.slices() {
        let slice = slice?;
        let (core, external) = slice.decode_blocks_for(&compression, fields, &tags)?;
        let records = slice.records_with_fields(
            pasta::fasta::Repository::default(),
            &header,
            &compression,
            &core,
            &external,
            MdNm::AsStored,
            fields,
            &tags,
        )?;
        unplaced += records.iter().filter(|r| r.reference_sequence_id(&header).is_none()).count() as u64;
    }
    Ok(unplaced)
}

/// chrX's and chrY's density (`amounts` over each contig's length) against the median
/// autosome's, of the 22: one the header lacks has no reads, so the reads of a chr20 reference
/// are said to be of 1 autosome, as those of a whole genome's header with reads on chr20 alone.
fn ratios_against_autosomes(reads: &Path, contigs: &[(String, u64)], amounts: &[f64]) -> Result<(f64, f64), SexError> {
    let density = |name: &str| {
        contigs.iter().zip(amounts).find(|((n, _), _)| n == name).map(|((_, len), &a)| a / (*len).max(1) as f64)
    };
    let prefix = prefix(contigs.iter().map(|(n, _)| n.as_str()));
    let named: Vec<(String, f64)> =
        (1..=22).map(|i| format!("{prefix}{i}")).map(|n| (n.clone(), density(&n).unwrap_or(0.0))).collect();
    let mut autosomes: Vec<f64> = named.iter().map(|&(_, d)| d).collect();
    autosomes.sort_by(f64::total_cmp);
    let median = autosomes[autosomes.len() / 2];
    if median <= 0.0 {
        let with = named.into_iter().filter(|&(_, d)| d > 0.0).map(|(n, _)| n).collect();
        return Err(SexError::Partial { reads: reads.into(), with });
    }
    if let Some((contig, deepest)) = named.iter().max_by(|a, b| a.1.total_cmp(&b.1))
        && deepest / median > MOST_EVEN
    {
        return Err(SexError::Uneven { reads: reads.into(), contig: contig.clone(), ratio: deepest / median });
    }
    let x = density(&format!("{prefix}X")).ok_or_else(|| SexError::Missing { reads: reads.into(), what: "chrX" })?;
    let y = density(&format!("{prefix}Y")).unwrap_or(0.0);
    Ok((x / median, y / median))
}

/// Little-endian reads, and CRAM's variable-length integers.
struct Bytes<R> {
    inner: R,
    path: PathBuf,
}

impl<R: Read> Bytes<R> {
    fn new(inner: R, path: &Path) -> Self {
        Self { inner, path: path.into() }
    }

    fn err(&self, source: io::Error) -> SexError {
        SexError::Read { path: self.path.clone(), source }
    }

    fn format(&self, reason: impl Into<String>) -> SexError {
        SexError::Format { path: self.path.clone(), reason: reason.into() }
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N], SexError> {
        let mut b = [0; N];
        self.inner.read_exact(&mut b).map_err(|e| self.err(e))?;
        Ok(b)
    }

    fn vec(&mut self, n: usize) -> Result<Vec<u8>, SexError> {
        let mut b = vec![0; n];
        self.inner.read_exact(&mut b).map_err(|e| self.err(e))?;
        Ok(b)
    }

    fn skip(&mut self, n: u64) -> Result<(), SexError> {
        let copied = io::copy(&mut (&mut self.inner).take(n), &mut io::sink()).map_err(|e| self.err(e))?;
        if copied == n { Ok(()) } else { Err(self.format("ends early")) }
    }

    fn u8(&mut self) -> Result<u8, SexError> {
        Ok(self.array::<1>()?[0])
    }

    fn i32(&mut self) -> Result<i32, SexError> {
        Ok(i32::from_le_bytes(self.array()?))
    }

    fn u64(&mut self) -> Result<u64, SexError> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    fn count(&mut self) -> Result<usize, SexError> {
        let n = self.i32()?;
        usize::try_from(n).map_err(|_| self.format(format!("negative count {n}")))
    }

    /// CRAM's ITF8.
    fn itf8(&mut self) -> Result<i32, SexError> {
        let first = self.u8()?;
        let extra = first.leading_ones().min(4);
        let b0 = u32::from(first);
        let mut v = if extra == 4 { b0 & 0x0f } else { b0 & (0xff >> (extra + 1)) };
        for i in 0..extra {
            let b = u32::from(self.u8()?);
            v = if i == 3 { (v << 4) | (b & 0x0f) } else { (v << 8) | b };
        }
        Ok(v as i32)
    }

    /// CRAM's LTF8.
    fn ltf8(&mut self) -> Result<i64, SexError> {
        let b0 = self.u8()?;
        let extra = b0.leading_ones();
        let mut v = if extra >= 7 { 0 } else { u64::from(b0 & (0xff >> (extra + 1))) };
        for _ in 0..extra {
            v = (v << 8) | u64::from(self.u8()?);
        }
        Ok(v as i64)
    }
}

fn open(path: &Path) -> Result<File, SexError> {
    File::open(path).map_err(|source| SexError::Read { path: path.into(), source })
}

/// Contig names and lengths from a BAM header.
fn bam_contigs(path: &Path) -> Result<Vec<(String, u64)>, SexError> {
    let mut r = Bytes::new(MultiGzDecoder::new(BufReader::new(open(path)?)), path);
    if &r.array::<4>()? != b"BAM\x01" {
        return Err(r.format("not BAM"));
    }
    let text = r.count()?;
    r.skip(text as u64)?;
    let n = r.count()?;
    (0..n)
        .map(|_| {
            let len = r.count()?;
            let name = r.vec(len)?;
            let name = String::from_utf8_lossy(name.strip_suffix(&[0]).unwrap_or(&name)).into_owned();
            Ok((name, r.count()? as u64))
        })
        .collect()
}

/// Mapped reads per reference, from a BAI's pseudo-bins, and its unplaced unmapped reads.
fn bai_mapped(path: &Path, n_contigs: usize) -> Result<(Vec<f64>, Option<u64>), SexError> {
    let mut r = Bytes::new(BufReader::new(open(path)?), path);
    if &r.array::<4>()? != b"BAI\x01" {
        return Err(r.format("not BAI"));
    }
    mapped_counts(&mut r, n_contigs, 37_450, Index::Bai)
}

/// Mapped reads per reference, from a CSI's pseudo-bins, and its unplaced unmapped reads.
fn csi_mapped(path: &Path, n_contigs: usize) -> Result<(Vec<f64>, Option<u64>), SexError> {
    let mut r = Bytes::new(MultiGzDecoder::new(BufReader::new(open(path)?)), path);
    if &r.array::<4>()? != b"CSI\x01" {
        return Err(r.format("not CSI"));
    }
    let _min_shift = r.i32()?;
    let depth = r.i32()?;
    let aux = r.count()?;
    r.skip(aux as u64)?;
    let pseudo = ((1u32 << ((depth + 1) * 3)) - 1) / 7 + 1;
    mapped_counts(&mut r, n_contigs, pseudo, Index::Csi)
}

/// BAM index formats, which lay out their bins differently.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Index {
    /// BAI: bins, then a linear index.
    Bai,
    /// CSI: each bin with an offset, and no linear index.
    Csi,
}

/// Reads the per-reference bins of an index, taking the mapped count from each reference's
/// pseudo-bin, and the count of unplaced unmapped reads after them, which is optional.
fn mapped_counts<R: Read>(
    r: &mut Bytes<R>,
    n_contigs: usize,
    pseudo: u32,
    index: Index,
) -> Result<(Vec<f64>, Option<u64>), SexError> {
    let n = r.count()?;
    if n != n_contigs {
        return Err(r.format(format!("indexes {n} references; the header has {n_contigs}")));
    }
    let mut mapped = vec![0.0; n];
    for m in &mut mapped {
        for _ in 0..r.count()? {
            let bin = u32::from_le_bytes(r.array()?);
            if index == Index::Csi {
                r.u64()?;
            }
            let chunks = r.count()?;
            if bin == pseudo && chunks == 2 {
                r.skip(16)?;
                *m = r.u64()? as f64;
                r.skip(8)?;
            } else {
                r.skip(16 * chunks as u64)?;
            }
        }
        if index == Index::Bai {
            let intervals = r.count()?;
            r.skip(8 * intervals as u64)?;
        }
    }
    let mut unplaced = [0; 8];
    let unplaced = match r.inner.read_exact(&mut unplaced) {
        Ok(()) => Some(u64::from_le_bytes(unplaced)),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => None,
        Err(e) => return Err(r.err(e)),
    };
    Ok((mapped, unplaced))
}

/// Contig names and lengths from a CRAM's SAM header.
fn cram_contigs(path: &Path) -> Result<Vec<(String, u64)>, SexError> {
    let mut r = Bytes::new(BufReader::new(open(path)?), path);
    let definition = r.array::<26>()?;
    if &definition[..4] != b"CRAM" {
        return Err(r.format("not CRAM"));
    }
    let major = definition[4];
    // Container header: length, reference, start, span, records, record counter, bases,
    // blocks, landmarks, then (from 3.0) a CRC.
    r.i32()?;
    for _ in 0..4 {
        r.itf8()?;
    }
    r.ltf8()?;
    r.ltf8()?;
    r.itf8()?;
    for _ in 0..r.itf8()? {
        r.itf8()?;
    }
    if major >= 3 {
        r.skip(4)?;
    }
    // The first block holds the header text.
    let method = r.u8()?;
    let _content_type = r.u8()?;
    r.itf8()?;
    let compressed = usize::try_from(r.itf8()?).map_err(|_| r.format("negative block size"))?;
    let _raw = r.itf8()?;
    let data = r.vec(compressed)?;
    let raw = match method {
        0 => data,
        1 => {
            let mut out = Vec::new();
            GzDecoder::new(&data[..]).read_to_end(&mut out).map_err(|e| r.err(e))?;
            out
        }
        m => return Err(r.format(format!("header block compressed with method {m}, which lrpipe does not read"))),
    };
    let len = raw.get(..4).map_or(0, |b| i32::from_le_bytes(b.try_into().expect("four bytes")));
    let text = raw.get(4..4 + usize::try_from(len).unwrap_or(0)).ok_or_else(|| r.format("short header block"))?;
    Ok(String::from_utf8_lossy(text)
        .lines()
        .filter(|l| l.starts_with("@SQ"))
        .filter_map(|l| {
            let field = |key: &str| l.split('\t').find_map(|f| f.strip_prefix(key));
            Some((field("SN:")?.to_owned(), field("LN:")?.parse().ok()?))
        })
        .collect())
}

/// Bytes of single-reference slices per reference, from a CRAM index.
fn crai_bytes(path: &Path, n_contigs: usize) -> Result<Vec<f64>, SexError> {
    let mut text = String::new();
    MultiGzDecoder::new(BufReader::new(open(path)?))
        .read_to_string(&mut text)
        .map_err(|source| SexError::Read { path: path.into(), source })?;
    let mut bytes = vec![0.0; n_contigs];
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let (Some(reference), Some(size)) = (f.first(), f.get(5)) else { continue };
        if let (Ok(i), Ok(size)) = (reference.parse::<usize>(), size.parse::<f64>())
            && let Some(b) = bytes.get_mut(i)
        {
            *b += size;
        }
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_panel_with_mismapped_reads_elsewhere_is_no_whole_genome() {
        // Edge-case testing: reads of one region, realigned genome-wide, put a few reads on every
        // autosome, so the median was tiny and chrY 10.5 times it called the sample male.
        let mut contigs: Vec<(String, u64)> = (1..=22).map(|i| (format!("chr{i}"), 1_000_000)).collect();
        contigs.extend([("chrX".to_owned(), 1_000_000), ("chrY".to_owned(), 1_000_000)]);
        let mut amounts = vec![10.0; 22];
        amounts[19] = 50_000.0;
        amounts.extend([4.65, 105.4]);
        let err = ratios_against_autosomes(Path::new("panel.bam"), &contigs, &amounts).expect_err("uneven");
        assert!(matches!(&err, SexError::Uneven { contig, .. } if contig == "chr20"), "{err}");
        // An even genome is measured; chrX or chrY far above an autosome pair fits no sex.
        amounts[19] = 10.0;
        let (x, y) = ratios_against_autosomes(Path::new("wgs.bam"), &contigs, &amounts).expect("even");
        assert_eq!(Estimate::new(x, y).sex, None, "chrY 10.5 times the autosomes is no Y count");
        assert_eq!(Estimate::new(0.5, 0.4).sex, Some(Sex::Male));
        assert_eq!(Estimate::new(1.0, 0.05).sex, Some(Sex::Female));
        assert_eq!(Estimate::new(1.5, 0.05).sex, Some(Sex::Female), "XXX is still no Y");
    }

    #[test]
    fn classifies_by_both_chromosomes() {
        assert_eq!(Estimate::new(0.486, 0.488).sex, Some(Sex::Male), "HG002 ONT by CRAM slice bytes");
        assert_eq!(Estimate::new(0.48, 1.4).sex, Some(Sex::Male), "HG002 HiFi by BAI read counts");
        assert_eq!(Estimate::new(0.986, 0.083).sex, Some(Sex::Female), "HG004 HiFi by BAI read counts");
        assert_eq!(Estimate::new(1.0, 0.001).sex, Some(Sex::Female));
        assert_eq!(Estimate::new(0.5, 0.0).sex, None, "one X without a Y (45,X) is neither");
        assert_eq!(Estimate::new(1.0, 0.5).sex, None, "two X with a Y (47,XXY) is neither");
    }

    #[test]
    fn reads_cram_integers() {
        let read = |bytes: &[u8]| Bytes::new(bytes, Path::new("t")).itf8().expect("readable");
        assert_eq!(read(&[0x05]), 5);
        assert_eq!(read(&[0x81, 0x00]), 256);
        assert_eq!(read(&[0xc1, 0x00, 0x00]), 65_536);
        assert_eq!(read(&[0xe1, 0x00, 0x00, 0x00]), 16_777_216);
        assert_eq!(read(&[0xff, 0xff, 0xff, 0xff, 0x0f]), -1);
        let ltf8 = |bytes: &[u8]| Bytes::new(bytes, Path::new("t")).ltf8().expect("readable");
        assert_eq!(ltf8(&[0x05]), 5);
        assert_eq!(ltf8(&[0x81, 0x00]), 256);
    }
}

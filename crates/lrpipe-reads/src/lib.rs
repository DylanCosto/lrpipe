//! Indexed BAM and CRAM reading shared by lrpipe's tools, on [pasta].
//!
//! [`Reference`] is a FASTA cache that every reader opened with it shares, so threads hold one
//! copy of each contig. [`IndexedReads`] opens a BAM (with `.bai` or `.csi`) or a CRAM (with
//! `.crai`, decoded with a [`Reference`]) and streams the records over a region to a callback,
//! as BAM records whichever the format; tools convert those into what they keep.
//! [`CramSource`] is the CRAM half on its own, for readers that keep their own files open.
//! [`FastaFile`] reads a plain or bgzipped FASTA at the `.fai`'s offsets, for tools that read
//! the reference a window at a time.
//! [`IndexedReads::sample_tag`] looks for a tag (the `HP` of haplotagged reads) in records from
//! several contigs, without reading the rest. [`haplotags`] reads the haplotypes a haplotag list
//! assigns reads, for tagging them as `haplotag` would without writing the alignments again.
//! [`regions`] reads a BED file of regions that a tool is restricted to. [`progress`] reports how
//! far a tool has got to `lrpipe run`'s display.

pub mod haplotags;
pub mod par;
pub mod progress;
pub mod regions;

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

use pasta::bam;
use pasta::core::region::Interval;
use pasta::core::{Position, Region};
use pasta::cram;
use pasta::fasta;
use pasta::sam;

pub use pasta::cram::io::reader::{Fields, MdNm, Tags};

/// Returns an uppercase VCF reference base. IUPAC ambiguity codes other than `N` use the first
/// possible nucleotide in alphabetical order, as defined by the VCF specification.
#[must_use]
pub fn vcf_reference_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' | b'D' | b'H' | b'M' | b'R' | b'V' | b'W' => b'A',
        b'B' | b'C' | b'S' | b'Y' => b'C',
        b'G' | b'K' => b'G',
        b'T' | b'U' => b'T',
        _ => b'N',
    }
}

/// What decoding a CRAM record gives a tool. BAM records are read as stored whatever this says.
///
/// The default decodes every field and adds the MD and NM tags the file left out, as htslib
/// does. A tool that reads less names only what it reads, and the data only the rest needs is
/// not decompressed: on ONT CRAMs, the `MM`/`ML` methylation tags and the base qualities are
/// most of the file. Fields left out are empty or unset in the records.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decoding {
    /// Whether the `MD` and `NM` tags the file left out are added (only those [`Self::tags`]
    /// selects, and only with the sequence and data decoded).
    pub md_nm: MdNm,
    /// The record fields decoded.
    pub fields: Fields,
    /// The tags decoded with [`Fields::DATA`].
    pub tags: Tags,
}

impl Default for Decoding {
    fn default() -> Self {
        Decoding { md_nm: MdNm::Regenerate, fields: Fields::all(), tags: Tags::All }
    }
}

/// Why a file could not be opened.
#[derive(Debug, thiserror::Error)]
pub enum OpenError {
    /// A file (or its index) could not be opened or read.
    #[error("reading {}", path.display())]
    Io {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// CRAM was given without a reference to decode it.
    #[error("{} is CRAM, which needs a reference", .0.display())]
    CramWithoutReference(PathBuf),
}

impl OpenError {
    /// The I/O error, for callers that report one with the path themselves: the cause of
    /// [`Self::Io`], or this error as one.
    #[must_use]
    pub fn into_io_error(self) -> io::Error {
        match self {
            Self::Io { source, .. } => source,
            other => io::Error::other(other),
        }
    }
}

/// The index of the BAM or CRAM at `reads` (`.bai`, `.csi` or `.crai` beside it, or `.bai` in
/// place of `.bam`) when it is older than the reads, as samtools warns: an index left from
/// before the file was written again points into the old file, and regions read through it can
/// come out empty without an error.
#[must_use]
pub fn stale_index(reads: &Path) -> Option<PathBuf> {
    let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
    let data = modified(reads)?;
    ["bai", "csi", "crai"]
        .into_iter()
        .map(|ext| reads.with_added_extension(ext))
        .chain([reads.with_extension("bai")])
        .find(|index| modified(index).is_some_and(|m| m < data))
}

/// Contigs [`FastaFile::check_index`] checks, spread over the index.
pub const FAI_CHECKS: usize = 64;

/// Reference sequences, loaded a contig at a time on first use and shared by every reader
/// opened with them (clones share the cache).
#[derive(Clone)]
pub struct Reference {
    repository: fasta::Repository,
    lengths: Arc<std::collections::HashMap<Vec<u8>, u64>>,
}

impl std::fmt::Debug for Reference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reference").finish_non_exhaustive()
    }
}

impl Reference {
    /// Opens a FASTA with a `.fai` index beside it: plain, or compressed with bgzip and with a
    /// `.gzi` index beside it too (as `samtools faidx` reads it).
    ///
    /// # Errors
    ///
    /// If the FASTA or an index cannot be read, or the FASTA is compressed but not with bgzip.
    pub fn open(path: &Path) -> Result<Self, OpenError> {
        let fai_path = path.with_added_extension("fai");
        let fai = fasta::fai::fs::read(&fai_path).map_err(|source| OpenError::Io { path: fai_path, source })?;
        Self::with_index(path, fai)
    }

    /// Opens a FASTA, plain or bgzipped (with a `.gzi` index), whose `.fai` index has already
    /// been read.
    ///
    /// # Errors
    ///
    /// If the FASTA or its `.gzi` cannot be read, or the FASTA is compressed but not with bgzip.
    pub fn with_index(path: &Path, fai: fasta::fai::Index) -> Result<Self, OpenError> {
        let lengths = fai.as_ref().iter().map(|r| (AsRef::<[u8]>::as_ref(r.name()).to_vec(), r.length())).collect();
        let file = FastaFile::open(path)?;
        file.check_index_once(&fai, path)?;
        let reader = file.indexed(fai);
        Ok(Reference {
            repository: fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader)),
            lengths: Arc::new(lengths),
        })
    }

    /// The cache itself, for readers built directly on pasta.
    #[must_use]
    pub fn repository(&self) -> &fasta::Repository {
        &self.repository
    }

    /// `contig`'s length by the `.fai` index; `None` when the reference has no such contig.
    #[must_use]
    pub fn length(&self, contig: &str) -> Option<u64> {
        self.lengths.get(contig.as_bytes()).copied()
    }

    /// `contig`'s sequence, read on first use; `None` when the reference has no such contig.
    ///
    /// # Errors
    ///
    /// If the sequence cannot be read.
    pub fn sequence(&self, contig: &str) -> io::Result<Option<Arc<fasta::record::Sequence>>> {
        self.repository.get(contig.as_bytes()).transpose()
    }

    /// Drops every contig read so far (they are kept until then, and a genome's worth would
    /// stay in memory); a contig needed again is read again.
    pub fn clear(&self) {
        self.repository.clear();
    }

    /// How many contigs are held, read and not yet dropped by [`Self::clear`].
    #[must_use]
    pub fn cached(&self) -> usize {
        self.repository.len()
    }
}

/// A FASTA file opened for reading by its `.fai` offsets, plain or compressed with bgzip.
///
/// The two are told apart by contents, not by name, as htslib does: the `.fai` of a bgzipped
/// FASTA holds offsets into the uncompressed text, which the `.gzi` maps to BGZF blocks.
#[derive(Debug)]
pub enum FastaFile {
    /// Uncompressed: the `.fai` offsets are the file's.
    Plain(File),
    /// Compressed with bgzip.
    Bgzf {
        /// The compressed file.
        file: File,
        /// Its `.gzi` index: where each BGZF block starts, compressed and uncompressed.
        gzi: pasta::bgzf::gzi::Index,
    },
}

/// The FASTA at `path`, plain or bgzipped (told apart by contents, as [`FastaFile`]), as pasta's
/// indexed reader with the `.fai` beside it: for readers built on pasta, which by itself takes a
/// FASTA for bgzipped by its name.
///
/// # Errors
///
/// If the FASTA or an index cannot be read, or the FASTA is compressed but not with bgzip.
pub fn indexed_fasta(path: &Path) -> Result<fasta::io::IndexedReader<fasta::io::BufReader<File>>, OpenError> {
    let fai_path = path.with_added_extension("fai");
    let fai = fasta::fai::fs::read(&fai_path).map_err(|source| OpenError::Io { path: fai_path, source })?;
    let file = FastaFile::open(path)?;
    file.check_index_once(&fai, path)?;
    Ok(file.indexed(fai))
}

/// Checks that the `.fai` beside the FASTA at `path` is this FASTA's (see
/// [`FastaFile::check_index`]).
///
/// # Errors
///
/// If the FASTA or its indexes cannot be read, or the `.fai` does not fit the FASTA.
pub fn check_fai(path: &Path) -> Result<(), OpenError> {
    let fai_path = path.with_added_extension("fai");
    let fai = fasta::fai::fs::read(&fai_path).map_err(|source| OpenError::Io { path: fai_path, source })?;
    FastaFile::open(path)?.check_index_once(&fai, path)
}

/// The FASTAs whose `.fai` has passed [`FastaFile::check_index`] in this process, with both
/// files' modification times: tools open the reference per contig or per thread.
static CHECKED: Mutex<BTreeSet<Checked>> = Mutex::new(BTreeSet::new());

/// A FASTA's path, and its and its `.fai`'s modification times.
type Checked = (PathBuf, Option<SystemTime>, Option<SystemTime>);

impl FastaFile {
    /// Checks that `fai` fits this FASTA (at `path`), on up to [`FAI_CHECKS`] contigs spread over
    /// it: a sequence starts after its header line's newline, and its first line ends where the
    /// `.fai` says. A `.fai` of another FASTA (other line lengths, another assembly) reads wrong
    /// bases everywhere, so every tool called thousands of false variants from it.
    ///
    /// # Errors
    ///
    /// If the FASTA cannot be read where the `.fai` points, or the bytes there do not fit.
    pub fn check_index(&self, fai: &fasta::fai::Index, path: &Path) -> Result<(), OpenError> {
        let records = fai.as_ref();
        let step = records.len().div_ceil(FAI_CHECKS).max(1);
        let mismatch = |name: &[u8], what: &str| OpenError::Io {
            path: path.with_added_extension("fai"),
            source: io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "its .fai is not this FASTA's index: {} ({what}); index it again with `samtools faidx {}`",
                    String::from_utf8_lossy(name),
                    path.display()
                ),
            ),
        };
        for record in records.iter().step_by(step) {
            let name: &[u8] = record.name().as_ref();
            let offset = record.position();
            let (bases, width) = (record.line_base_count().get(), record.line_width().get());
            // Bytes past the end of the file are not checked: a FASTA cut short fails when it is
            // read, where an index of another FASTA reads wrong bases without a word.
            if offset > 0 {
                let mut before = [0u8];
                if self.read_exact_at(&mut before, offset - 1).is_ok() && before[0] != b'\n' {
                    return Err(mismatch(name, "its sequence does not start after a header line"));
                }
            }
            if record.length() > bases && width > bases {
                let mut end = vec![0u8; (width - bases) as usize];
                if self.read_exact_at(&mut end, offset + bases).is_ok()
                    && (end.last() != Some(&b'\n') || end.iter().any(|&b| b != b'\n' && b != b'\r'))
                {
                    return Err(mismatch(name, "its lines are not as long as the index says"));
                }
            }
        }
        Ok(())
    }

    /// [`check_index`](Self::check_index), once per FASTA and `.fai` in a process.
    fn check_index_once(&self, fai: &fasta::fai::Index, path: &Path) -> Result<(), OpenError> {
        let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
        let key = (path.to_path_buf(), modified(path), modified(&path.with_added_extension("fai")));
        if CHECKED.lock().unwrap_or_else(std::sync::PoisonError::into_inner).contains(&key) {
            return Ok(());
        }
        self.check_index(fai, path)?;
        CHECKED.lock().unwrap_or_else(std::sync::PoisonError::into_inner).insert(key);
        Ok(())
    }

    /// pasta's indexed reader of this FASTA, with its `.fai` index `fai`.
    #[must_use]
    pub fn indexed(self, fai: fasta::fai::Index) -> fasta::io::IndexedReader<fasta::io::BufReader<File>> {
        let reader = match self {
            FastaFile::Plain(file) => fasta::io::BufReader::Uncompressed(io::BufReader::new(file)),
            FastaFile::Bgzf { file, gzi } => fasta::io::BufReader::Bgzf(pasta::bgzf::io::IndexedReader::new(file, gzi)),
        };
        fasta::io::indexed_reader::IndexedReader::new(reader, fai)
    }

    /// Opens the FASTA at `path`, with the `.gzi` beside it when it is bgzipped.
    ///
    /// # Errors
    ///
    /// If the FASTA or its `.gzi` cannot be read, or the FASTA is compressed but not with bgzip.
    pub fn open(path: &Path) -> Result<Self, OpenError> {
        let io_error = |path: &Path| {
            let path = path.to_owned();
            move |source| OpenError::Io { path, source }
        };
        let mut file = File::open(path).map_err(io_error(path))?;
        if !is_gzip(&mut file).map_err(io_error(path))? {
            return Ok(FastaFile::Plain(file));
        }
        if !is_bgzf(&mut file).map_err(io_error(path))? {
            return Err(OpenError::Io {
                path: path.into(),
                source: io::Error::new(
                    io::ErrorKind::InvalidData,
                    "compressed with gzip, not bgzip: recompress it with `bgzip` and index it with \
                     `samtools faidx`, or decompress it",
                ),
            });
        }
        let gzi_path = path.with_added_extension("gzi");
        let gzi = pasta::bgzf::gzi::fs::read(&gzi_path).map_err(io_error(&gzi_path))?;
        Ok(FastaFile::Bgzf { file, gzi })
    }

    /// Fills `buf` with the FASTA's text from `offset` on, an offset as the `.fai` counts it (into
    /// the uncompressed text). Reads at a position without moving a shared one, so threads may
    /// read the same `FastaFile` at once.
    ///
    /// # Errors
    ///
    /// If the file cannot be read, or ends before `buf` is full.
    pub fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> io::Result<()> {
        use std::os::unix::fs::FileExt;
        match self {
            FastaFile::Plain(file) => file.read_exact_at(buf, offset),
            FastaFile::Bgzf { file, gzi } => {
                let mut reader = pasta::bgzf::io::Reader::new(At { file, position: 0 });
                reader.seek_by_uncompressed_position(gzi, offset)?;
                reader.read_exact(buf)
            }
        }
    }
}

/// A file read from a position of its own, by positioned reads.
struct At<'a> {
    file: &'a File,
    position: u64,
}

impl Read for At<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        use std::os::unix::fs::FileExt;
        let n = self.file.read_at(buf, self.position)?;
        self.position += n as u64;
        Ok(n)
    }
}

impl Seek for At<'_> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let to = match pos {
            SeekFrom::Start(p) => Some(p),
            SeekFrom::Current(d) => self.position.checked_add_signed(d),
            SeekFrom::End(d) => self.file.metadata()?.len().checked_add_signed(d),
        };
        self.position = to.ok_or_else(|| io::Error::from(io::ErrorKind::InvalidInput))?;
        Ok(self.position)
    }
}

/// Whether `file` starts with gzip's magic bytes; leaves it at its start.
fn is_gzip(file: &mut File) -> io::Result<bool> {
    let mut magic = [0; 2];
    let n = read_up_to(file, &mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    Ok(n == 2 && magic == [0x1f, 0x8b])
}

/// Whether the gzip member `file` starts with is a BGZF block (its `BC` extra field); leaves
/// it at its start.
fn is_bgzf(file: &mut File) -> io::Result<bool> {
    let mut header = [0; 16];
    let n = read_up_to(file, &mut header)?;
    file.seek(SeekFrom::Start(0))?;
    // FLG.FEXTRA, XLEN 6, then the subfield `BC` of length 2.
    Ok(n == 16 && header[3] & 4 != 0 && header[10..12] == [6, 0] && header[12..16] == [b'B', b'C', 2, 0])
}

/// Reads until `buf` is full or the file ends; the bytes read.
fn read_up_to(file: &mut File, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match file.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// Whether `path` names a CRAM file (by extension).
#[must_use]
pub fn is_cram(path: &Path) -> bool {
    path.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"))
}

enum Inner {
    Bam(Box<bam::io::IndexedReader<pasta::bgzf::io::Reader<File>>>),
    Cram {
        reader: Box<cram::io::Reader<File>>,
        index: cram::crai::Index,
        reference: fasta::Repository,
        decoding: Decoding,
    },
}

/// What [`IndexedReads::sample_tag`] found.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TagSample {
    /// Mapped records looked at.
    pub records: usize,
    /// Whether one of them carries the tag with an integer value: the search stops at the
    /// first that does.
    pub found: bool,
    /// Records looked at whose tag has a value of another type (`HP:Z:1`), which tools that read
    /// an integer take for no tag.
    pub other_type: usize,
}

/// Contigs at least this long are sampled first by [`IndexedReads::sample_tag`]: a short one
/// (the mitochondrial genome) can be all a file starts with, and haplotagging leaves its reads
/// untagged.
pub const LONG_CONTIG: usize = 1_000_000;

/// Ends a query once a search has its answer.
#[derive(Debug)]
struct Enough;

impl std::fmt::Display for Enough {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("enough records")
    }
}

impl std::error::Error for Enough {}

/// An open, indexed alignment file. One per thread: queries move its position.
pub struct IndexedReads {
    path: PathBuf,
    header: sam::Header,
    inner: Inner,
}

impl std::fmt::Debug for IndexedReads {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndexedReads").field("path", &self.path).finish_non_exhaustive()
    }
}

impl IndexedReads {
    /// Opens `path`: BAM with `.bai` or `.csi`, or CRAM with `.crai`, decoded with `reference`.
    ///
    /// # Errors
    ///
    /// If the file or its index cannot be opened, or CRAM is given without a reference.
    pub fn open(path: &Path, reference: Option<&Reference>) -> Result<Self, OpenError> {
        let io = |source| OpenError::Io { path: path.into(), source };
        let (header, inner) = if is_cram(path) {
            let reference = reference.ok_or_else(|| OpenError::CramWithoutReference(path.into()))?;
            let index = cram::fs::read_associated_index(path).map_err(io)?;
            let mut reader = File::open(path).map(cram::io::Reader::new).map_err(io)?;
            let header = reader.read_header().map_err(io)?;
            let inner = Inner::Cram {
                reader: Box::new(reader),
                index,
                reference: reference.repository.clone(),
                decoding: Decoding::default(),
            };
            (header, inner)
        } else {
            let mut r = bam::io::indexed_reader::Builder::default().build_from_path(path).map_err(io)?;
            (r.read_header().map_err(io)?, Inner::Bam(Box::new(r)))
        };
        Ok(IndexedReads { path: path.into(), header, inner })
    }

    /// What decoding CRAM records gives ([`Decoding::default`] until set).
    pub fn set_decoding(&mut self, decoding: Decoding) {
        if let Inner::Cram { decoding: d, .. } = &mut self.inner {
            *d = decoding;
        }
    }

    /// The file's path.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The file's header.
    #[must_use]
    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// The compressed bytes the index says a query of `[start, end)` (0-based) of `contig`
    /// reads: BAM chunks' spans, or CRAM slices' lengths. A measure of how much data the
    /// range holds, which grows with depth; 0 for a contig the file does not name.
    ///
    /// # Errors
    ///
    /// If the BAM index cannot be queried.
    pub fn indexed_bytes(&self, contig: &str, start: i64, end: i64) -> io::Result<u64> {
        let Some(id) = self.header.reference_sequences().get_index_of(contig.as_bytes()) else { return Ok(0) };
        let Some((s, e)) = positions(start, end) else { return Ok(0) };
        Ok(match &self.inner {
            Inner::Bam(r) => {
                let chunks = r.index().query(id, (s..=e).into())?;
                chunks.iter().map(|c| c.end().compressed().saturating_sub(c.start().compressed())).sum()
            }
            Inner::Cram { index, .. } => index
                .iter()
                .filter(|x| x.reference_sequence_id() == Some(id))
                .filter(|x| {
                    // A slice's span counts from its start; unknown starts reach everywhere.
                    x.alignment_start().is_none_or(|a| {
                        let a = usize::from(a);
                        a <= usize::from(e) && a + x.alignment_span() > usize::from(s)
                    })
                })
                .map(|x| x.slice_length())
                .sum(),
        })
    }

    /// Calls `f` with each record overlapping `[start, end)` (0-based) of `contig`, in file
    /// order, until `f` fails. A contig the file does not name has none. CRAM records come
    /// encoded as BAM records (see [`CramSource::query`]).
    ///
    /// # Errors
    ///
    /// If records cannot be read, or `f` fails.
    pub fn query(
        &mut self,
        contig: &str,
        start: i64,
        end: i64,
        mut f: impl FnMut(&bam::Record) -> io::Result<()>,
    ) -> io::Result<()> {
        let Some(id) = self.header.reference_sequences().get_index_of(contig.as_bytes()) else { return Ok(()) };
        let Some((s, e)) = positions(start, end) else { return Ok(()) };
        match &mut self.inner {
            Inner::Bam(r) => {
                let region = Region::new(contig, s..=e);
                for rec in r.query(&self.header, &region)?.records() {
                    f(&rec?)?;
                }
            }
            Inner::Cram { reader, index, reference, decoding } => {
                let source = CramSource { header: &self.header, index, reference, decoding };
                source.query(reader, id, (s..=e).into(), f)?;
            }
        }
        Ok(())
    }
}

impl IndexedReads {
    /// Calls `f` with every record placed on `contig` (unmapped reads placed there with their
    /// mates too), in file order, until `f` fails: a whole contig, for reading a file through.
    /// A contig the file does not name has none.
    ///
    /// # Errors
    ///
    /// If records cannot be read, or `f` fails.
    pub fn query_contig(&mut self, contig: &str, mut f: impl FnMut(&bam::Record) -> io::Result<()>) -> io::Result<()> {
        let Some(id) = self.header.reference_sequences().get_index_of(contig.as_bytes()) else { return Ok(()) };
        match &mut self.inner {
            Inner::Bam(r) => {
                let region = Region::new(contig, ..);
                for rec in r.query(&self.header, &region)?.records() {
                    f(&rec?)?;
                }
            }
            Inner::Cram { reader, index, reference, decoding } => {
                let source = CramSource { header: &self.header, index, reference, decoding };
                source.query(reader, id, (..).into(), f)?;
            }
        }
        Ok(())
    }
}

impl IndexedReads {
    /// Looks for `tag`, with an integer value, in `per_contig` mapped records (secondary and
    /// supplementary alignments too) from the middle of each of the first `contigs` contigs with
    /// reads, stopping at the first record that carries it. Contigs of at least
    /// [`LONG_CONTIG`] come first, in the header's order, then the others; a contig's records
    /// are read from its middle on, and from its start if that half holds too few. Only a
    /// few BGZF blocks or CRAM containers of each contig are read, and of a CRAM only the
    /// flags, positions, CIGAR and `tag` are decoded.
    ///
    /// # Errors
    ///
    /// If records cannot be read.
    pub fn sample_tag(&mut self, tag: [u8; 2], per_contig: usize, contigs: usize) -> io::Result<TagSample> {
        let mut order: Vec<(String, usize)> = self
            .header
            .reference_sequences()
            .iter()
            .map(|(name, map)| (String::from_utf8_lossy(name).into_owned(), map.length().get()))
            .collect();
        // Stable: the header's order within each kind.
        order.sort_by_key(|&(_, len)| len < LONG_CONTIG);
        let with_reads: Option<std::collections::HashSet<usize>> = match &self.inner {
            Inner::Cram { index, .. } => {
                Some(index.iter().filter_map(cram::crai::Record::reference_sequence_id).collect())
            }
            Inner::Bam(_) => None,
        };
        let previous = match &mut self.inner {
            Inner::Cram { decoding, .. } => Some(std::mem::replace(
                decoding,
                Decoding {
                    md_nm: MdNm::AsStored,
                    fields: Fields::FLAGS
                        | Fields::REFERENCE_SEQUENCE_ID
                        | Fields::ALIGNMENT_START
                        | Fields::CIGAR
                        | Fields::DATA,
                    tags: Tags::Only(vec![sam::alignment::record::data::field::Tag::from(tag)]),
                },
            )),
            Inner::Bam(_) => None,
        };
        let mut sample = TagSample::default();
        let result = (|| {
            let mut sampled = 0;
            for (name, len) in &order {
                if sample.found || sampled == contigs {
                    break;
                }
                let id = self.header.reference_sequences().get_index_of(name.as_bytes());
                if with_reads.as_ref().is_some_and(|ids| id.is_none_or(|id| !ids.contains(&id))) {
                    continue;
                }
                let (len, mid) = (*len as i64, (*len / 2) as i64);
                let mut taken = 0;
                for (start, end) in [(mid, len), (0, mid)] {
                    if taken == per_contig || sample.found {
                        break;
                    }
                    // A record across the middle is seen by both queries, which matters only to
                    // the count.
                    let read = self.query(name, start, end, |record| {
                        if record.flags().is_unmapped() {
                            return Ok(());
                        }
                        taken += 1;
                        match record.data().get(&tag).transpose()? {
                            Some(value) if value.as_int().is_some() => sample.found = true,
                            Some(_) => sample.other_type += 1,
                            None => {}
                        }
                        if sample.found || taken == per_contig {
                            return Err(io::Error::other(Enough));
                        }
                        Ok(())
                    });
                    match read {
                        Err(e) if e.get_ref().is_some_and(|e| e.is::<Enough>()) => {}
                        other => other?,
                    }
                }
                sample.records += taken;
                if taken > 0 {
                    sampled += 1;
                }
            }
            Ok(())
        })();
        if let (Some(previous), Inner::Cram { decoding, .. }) = (previous, &mut self.inner) {
            *decoding = previous;
        }
        result.map(|()| sample)
    }
}

/// What decoding a CRAM's records takes besides the file: its header and index, and the
/// reference its bases are stored against.
#[derive(Clone, Copy)]
pub struct CramSource<'a> {
    /// The CRAM's header.
    pub header: &'a sam::Header,
    /// The CRAM's `.crai` index.
    pub index: &'a cram::crai::Index,
    /// The reference.
    pub reference: &'a fasta::Repository,
    /// What decoding gives the records.
    pub decoding: &'a Decoding,
}

impl std::fmt::Debug for CramSource<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CramSource").field("decoding", self.decoding).finish_non_exhaustive()
    }
}

impl CramSource<'_> {
    /// Calls `f` with each record of reference sequence `id` overlapping `interval` (1-based,
    /// inclusive), in file order, until `f` fails: the records pasta's CRAM query gives, but
    /// encoded straight to BAM records. The query builds a `RecordBuf` for each record, which a
    /// tool would then read through the slower generic accessors or write out again as BAM.
    ///
    /// The records hold what the query's would, as BAM stores it: a mapping quality the CRAM
    /// does not store (an unmapped read's) is 255, bases are BAM's 4-bit codes, and a CIGAR of
    /// more than 65,535 operations is moved to the `CG` tag (which [`bam::Record::cigar`] reads).
    ///
    /// # Errors
    ///
    /// If a container or slice cannot be read or decoded, a record cannot be encoded as BAM, or
    /// `f` fails.
    pub fn query<R: Read + Seek>(
        &self,
        reader: &mut cram::io::Reader<R>,
        id: usize,
        interval: Interval,
        mut f: impl FnMut(&bam::Record) -> io::Result<()>,
    ) -> io::Result<()> {
        let mut container = cram::io::reader::Container::default();
        let mut buf = Vec::new();
        for entry in self.index.iter().filter(|x| slice_overlaps(x, id, interval)) {
            reader.seek(SeekFrom::Start(entry.offset()))?;
            if reader.read_container(&mut container)? == 0 {
                break;
            }
            let compression_header = container.compression_header()?;
            let slice = container
                .header()
                .landmarks()
                .iter()
                .position(|&n| n as u64 == entry.landmark())
                .and_then(|i| container.slices().nth(i))
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid slice landmark"))??;
            let Decoding { md_nm, fields, tags } = self.decoding;
            let (core, external) = slice.decode_blocks_for(&compression_header, *fields, tags)?;
            let records = slice.records_with_fields(
                self.reference.clone(),
                self.header,
                &compression_header,
                &core,
                &external,
                *md_nm,
                *fields,
                tags,
            )?;
            for record in &records {
                record.encode_bam_record(self.header, &mut buf)?;
                let record = bam::Record::try_from(std::mem::take(&mut buf))?;
                if record_overlaps(&record, id, interval)? {
                    f(&record)?;
                }
                buf = record.into();
            }
        }
        Ok(())
    }
}

/// Whether an index entry's slice can hold records of `id` overlapping `interval` (pasta's
/// query's test).
fn slice_overlaps(entry: &cram::crai::Record, id: usize, interval: Interval) -> bool {
    if entry.reference_sequence_id() != Some(id) {
        return false;
    }
    if is_unbounded(interval) {
        return true;
    }
    let Some(start) = entry.alignment_start() else { return false };
    let span = entry.alignment_span();
    if span == 0 {
        return false;
    }
    let end = start.checked_add(span - 1).expect("slice end overflows");
    interval.intersects((start..=end).into())
}

/// Whether a record lies on `id` and overlaps `interval` (pasta's query's test: an alignment
/// that consumes no reference covers its start).
fn record_overlaps(record: &bam::Record, id: usize, interval: Interval) -> io::Result<bool> {
    if record.reference_sequence_id().transpose()? != Some(id) {
        return Ok(false);
    }
    if is_unbounded(interval) {
        return Ok(true);
    }
    let Some(start) = record.alignment_start().transpose()? else { return Ok(false) };
    let end = start.checked_add(reference_span(record.cigar().as_bytes()).max(1) - 1).expect("alignment end overflows");
    Ok(interval.intersects((start..=end).into()))
}

/// The reference bases packed BAM CIGAR operations cover (M, D, N, = and X).
fn reference_span(cigar: &[u8]) -> usize {
    cigar
        .as_chunks::<4>()
        .0
        .iter()
        .map(|op| u32::from_le_bytes(*op))
        .filter(|op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|op| (op >> 4) as usize)
        .sum()
}

fn is_unbounded(interval: Interval) -> bool {
    interval.start().is_none() && interval.end().is_none()
}

/// `[start, end)`, 0-based, as the 1-based inclusive positions of a region, clamped to start
/// at the contig's first base and to cover at least it (as the tools have always queried).
fn positions(start: i64, end: i64) -> Option<(Position, Position)> {
    Some((Position::new(start.max(0) as usize + 1)?, Position::new(end.max(1) as usize)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_fai_of_another_fasta_is_refused() {
        // Edge-case testing: a .fai with other line lengths read wrong bases everywhere.
        let dir = std::env::temp_dir().join(format!("lrpipe-reads-fai-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let fasta = dir.join("ref.fa");
        let fai = dir.join("ref.fa.fai");
        let check = |text: &str, index: &str| {
            std::fs::write(&fasta, text).expect("FASTA");
            std::fs::write(&fai, index).expect(".fai");
            check_fai(&fasta).map_err(|e| e.into_io_error().to_string())
        };
        let two = ">a\nACGT\nAC\n>b\nGGGG\nG\n";
        assert!(check(two, "a\t6\t3\t4\t5\nb\t5\t14\t4\t5\n").is_ok(), "its own index");
        let crlf = ">a\r\nACGT\r\nAC\r\n";
        assert!(check(crlf, "a\t6\t4\t4\t6\n").is_ok(), "CRLF lines");
        let err = check(two, "a\t6\t3\t3\t4\nb\t5\t14\t4\t5\n").expect_err("other line lengths");
        assert!(err.contains("its .fai is not this FASTA's index: a (its lines"), "{err}");
        let err = check(two, "a\t6\t3\t4\t5\nb\t5\t12\t4\t5\n").expect_err("another offset");
        assert!(err.contains("b (its sequence does not start after a header line)"), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn cram_is_known_by_extension() {
        assert!(is_cram(Path::new("x.cram")));
        assert!(is_cram(Path::new("x.CRAM")));
        assert!(!is_cram(Path::new("x.bam")));
    }

    #[test]
    fn regions_are_one_based_and_inclusive() {
        let (s, e) = positions(0, 10).expect("non-empty");
        assert_eq!((usize::from(s), usize::from(e)), (1, 10));
        let (s, e) = positions(-5, 0).expect("clamped");
        assert_eq!((usize::from(s), usize::from(e)), (1, 1));
    }

    #[test]
    fn missing_files_say_which() {
        let err = IndexedReads::open(Path::new("/nonexistent.bam"), None).expect_err("no such file");
        assert!(err.to_string().contains("/nonexistent.bam"));
        let err = IndexedReads::open(Path::new("/nonexistent.cram"), None).expect_err("no reference");
        assert!(matches!(err, OpenError::CramWithoutReference(_)));
    }

    /// A FASTA of two contigs (60 bases a line), its `.fai`, and its sequences.
    fn fasta() -> (String, String, Vec<(&'static str, Vec<u8>)>) {
        let contigs: Vec<(&str, Vec<u8>)> = [("chr1", 5000u32), ("chr2", 3001)]
            .into_iter()
            .map(|(name, len)| (name, (0..len).map(|i| b"ACGTTGCA"[((i * 7 + i / 13) % 8) as usize]).collect()))
            .collect();
        let (mut text, mut fai) = (String::new(), String::new());
        for (name, seq) in &contigs {
            text.push_str(&format!(">{name} test\n"));
            fai.push_str(&format!("{name}\t{}\t{}\t60\t61\n", seq.len(), text.len()));
            for line in seq.chunks(60) {
                text.push_str(std::str::from_utf8(line).expect("ASCII"));
                text.push('\n');
            }
        }
        (text, fai, contigs)
    }

    /// `text` in BGZF blocks of 1000 bytes, and their `.gzi`.
    fn bgzip(text: &[u8]) -> (Vec<u8>, pasta::bgzf::gzi::Index) {
        use std::io::Write;
        let mut w = pasta::bgzf::io::Writer::new(Vec::new());
        let mut blocks = Vec::new();
        for (i, chunk) in text.chunks(1000).enumerate() {
            if i > 0 {
                blocks.push((w.position(), (i * 1000) as u64));
            }
            w.write_all(chunk).expect("written");
            w.flush().expect("flushed");
        }
        (w.finish().expect("finished"), pasta::bgzf::gzi::Index::from(blocks))
    }

    /// LRUN-06: a bgzipped FASTA with its `.gzi` gives the plain FASTA's sequences, whatever it
    /// is named; one compressed with gzip alone, or without its `.gzi`, is an error.
    #[test]
    fn bgzipped_reference_reads_as_the_plain_one() {
        let dir = std::env::temp_dir().join(format!("lrpipe-reads-bgzf-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let (text, fai, contigs) = fasta();
        let write = |name: &str, data: &[u8]| {
            let path = dir.join(name);
            std::fs::write(&path, data).expect("written");
            std::fs::write(path.with_added_extension("fai"), &fai).expect("written");
            path
        };
        let plain = write("plain.fa", text.as_bytes());
        let (bgzf, gzi) = bgzip(text.as_bytes());
        let mut opened = vec![Reference::open(&plain).expect("plain")];
        // Named as references are distributed, and with a name that does not say it.
        for name in ["ref.fa.gz", "ref.fa"] {
            let path = write(name, &bgzf);
            pasta::bgzf::gzi::fs::write(path.with_added_extension("gzi"), &gzi).expect("written");
            opened.push(Reference::open(&path).expect("bgzipped"));
        }
        for reference in &opened {
            for (name, seq) in &contigs {
                let got = reference.sequence(name).expect("read").expect("present");
                assert_eq!(AsRef::<[u8]>::as_ref(got.as_ref()), &seq[..], "{name}");
                assert_eq!(reference.length(name), Some(seq.len() as u64));
            }
            assert!(reference.sequence("chr3").expect("read").is_none());
        }

        let no_gzi = write("no_gzi.fa.gz", &bgzf);
        let err = Reference::open(&no_gzi).expect_err("no .gzi");
        assert!(err.to_string().contains("no_gzi.fa.gz.gzi"), "{err}");

        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut gz, text.as_bytes()).expect("compressed");
        let gzip = write("gzip.fa.gz", &gz.finish().expect("compressed"));
        let err = Reference::open(&gzip).expect_err("gzip, not bgzip");
        assert!(format!("{:?}", std::error::Error::source(&err)).contains("not bgzip"), "{err:?}");
        std::fs::remove_dir_all(&dir).expect("removed");
    }

    /// A bgzipped FASTA reads at the `.fai`'s offsets as the plain one does: within a block,
    /// across blocks, to the end, and from several threads at once; past the end is an error.
    #[test]
    fn bgzipped_fasta_reads_at_offsets_as_the_plain_one() {
        let dir = std::env::temp_dir().join(format!("lrpipe-reads-at-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let (text, _, _) = fasta();
        let plain = dir.join("plain.fa");
        std::fs::write(&plain, &text).expect("written");
        let (bgzf, gzi) = bgzip(text.as_bytes());
        let packed = dir.join("packed.fa");
        std::fs::write(&packed, bgzf).expect("written");
        pasta::bgzf::gzi::fs::write(packed.with_added_extension("gzi"), &gzi).expect("written");
        let files = [FastaFile::open(&plain).expect("plain"), FastaFile::open(&packed).expect("bgzipped")];
        assert!(matches!(files[0], FastaFile::Plain(_)) && matches!(files[1], FastaFile::Bgzf { .. }));
        let len = text.len();
        std::thread::scope(|scope| {
            for file in &files {
                for (offset, n) in [(0, 10), (990, 20), (1000, 1000), (1500, 2600), (len - 7, 7), (0, len)] {
                    let text = &text;
                    scope.spawn(move || {
                        let mut buf = vec![0; n];
                        file.read_exact_at(&mut buf, offset as u64).expect("read");
                        assert_eq!(buf, text.as_bytes()[offset..offset + n], "{offset}+{n}");
                    });
                }
            }
        });
        for file in &files {
            let mut buf = [0; 2];
            assert!(file.read_exact_at(&mut buf, (len - 1) as u64).is_err(), "past the end");
        }
        std::fs::remove_dir_all(&dir).expect("removed");
    }
}

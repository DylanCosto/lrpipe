//! CRAM files read with pasta::cram; records are converted to BAM layout
//! (as htslib decodes CRAM into `bam1_t`).

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::thread::JoinHandle;

use anyhow::{Context, Result, bail};
use pasta::core::{Position, Region};
use pasta::cram::io::reader::{Fields, MdNm, Tags};
use pasta::sam::alignment::record::data::field::Tag;
use pasta::{cram, fasta, sam};

use super::bam::find_index;
use super::header::Header;
use super::record::Record;

/// Adds MD and NM tags as htslib does when decoding CRAM (`decode_md`, on by
/// default when reading) for a mapped read with a sequence that has none.
/// htslib derives them from the CRAM features; htslib's encoder uses base
/// features only at mismatches, so comparing the read with the reference
/// gives the same result: a reference N always counts as a mismatch,
/// deletions add `^` and the reference bases.
fn add_md_nm(rec: &mut Record, reference: &[u8]) {
    use super::record::{CDEL, CDIFF, CEQUAL, CINS, CMATCH, CREF_SKIP, CSOFT_CLIP};
    let need_md = rec.aux(b"MD").is_err();
    let need_nm = rec.aux(b"NM").is_err();
    if rec.is_unmapped() || rec.tid < 0 || rec.l_seq <= 0 || !(need_md || need_nm) {
        return;
    }
    let seq = rec.seq();
    let ref_base = |p: i64| -> u8 {
        usize::try_from(p).ok().and_then(|p| reference.get(p)).map_or(b'N', |b| b.to_ascii_uppercase())
    };
    let mut md = String::new();
    let mut md_dist: u32 = 0;
    let mut nm: u32 = 0;
    let (mut ref_pos, mut seq_pos) = (rec.pos(), 0usize);
    for (op, len) in rec.cigar().iter() {
        match op {
            CMATCH | CEQUAL | CDIFF => {
                for _ in 0..len {
                    let r = ref_base(ref_pos);
                    if r == b'N' || seq.get(seq_pos) != Some(&r) {
                        md.push_str(&md_dist.to_string());
                        md.push(r as char);
                        md_dist = 0;
                        nm += 1;
                    } else {
                        md_dist += 1;
                    }
                    ref_pos += 1;
                    seq_pos += 1;
                }
            }
            CINS => {
                nm += len;
                seq_pos += len as usize;
            }
            CSOFT_CLIP => seq_pos += len as usize,
            CDEL => {
                md.push_str(&md_dist.to_string());
                md.push('^');
                for _ in 0..len {
                    md.push(ref_base(ref_pos) as char);
                    ref_pos += 1;
                }
                md_dist = 0;
                nm += len;
            }
            CREF_SKIP => ref_pos += len as i64,
            _ => {}
        }
    }
    md.push_str(&md_dist.to_string());
    let mut fields = Vec::new();
    if need_md {
        fields.extend_from_slice(b"MDZ");
        fields.extend_from_slice(md.as_bytes());
        fields.push(0);
    }
    if need_nm {
        fields.extend_from_slice(b"NM");
        if nm <= u8::MAX as u32 {
            fields.push(b'C');
            fields.push(nm as u8);
        } else if nm <= u16::MAX as u32 {
            fields.push(b'S');
            fields.extend_from_slice(&(nm as u16).to_le_bytes());
        } else {
            fields.push(b'I');
            fields.extend_from_slice(&nm.to_le_bytes());
        }
    }
    rec.insert_aux_before_rg(&fields);
}

enum Iter {
    Rest,
    /// Records of `tid` overlapping [beg, end)
    Region {
        tid: i32,
        beg: i64,
        end: i64,
    },
    /// The same, decoded on worker threads (see [`Threaded`])
    Threaded(Box<Threaded>),
    Unmapped,
    Done,
}

/// Decoded batches a threaded region read keeps ready ahead of [`CramReader::read`]. With
/// [`Decoding::Alleles`], a slice of HG002 ONT reads (about 4.5 Mbases) decodes to a few MB.
const THREADED_BATCHES_AHEAD: usize = 2;

/// A region read decoded on worker threads: pasta's `fold_queries` runs on a thread of its own,
/// with a reader of its own, and sends each slice's records, in query order, through a bounded
/// channel. Dropping it stops the decoding.
struct Threaded {
    batches: mpsc::Receiver<Result<Vec<Record>>>,
    current: std::vec::IntoIter<Record>,
    worker: Option<JoinHandle<()>>,
}

impl Threaded {
    /// The next record, `None` at the end of the region.
    fn next(&mut self) -> Result<Option<Record>> {
        loop {
            if let Some(record) = self.current.next() {
                return Ok(Some(record));
            }
            match self.batches.recv() {
                Ok(batch) => self.current = batch?.into_iter(),
                // The thread is done: it sent every batch, or panicked.
                Err(mpsc::RecvError) => {
                    if self.worker.take().is_some_and(|w| w.join().is_err()) {
                        bail!("the thread decoding CRAM records panicked");
                    }
                    return Ok(None);
                }
            }
        }
    }
}

impl Drop for Threaded {
    fn drop(&mut self) {
        // Closing the channel stops the decoding at its next batch.
        let (_, closed) = mpsc::sync_channel(0);
        drop(std::mem::replace(&mut self.batches, closed));
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

pub struct CramReader {
    path: PathBuf,
    /// Threads decoding a region with [`Decoding::Alleles`]; one decodes on the caller's thread.
    threads: usize,
    reader: cram::io::Reader<File>,
    sam_header: sam::Header,
    pub header: Header,
    index: Option<cram::crai::Index>,
    first_container: u64,
    records: std::vec::IntoIter<Record>,
    iter: Iter,
    /// The record being converted, as BAM stores it (without `block_size`)
    encoded: Vec<u8>,
    repository: fasta::Repository,
    has_reference: bool,
    /// Whether `repository` is shared with a caller who clears it (see
    /// [`CramReader::share_reference`]); otherwise each fetch clears it.
    shared_reference: bool,
    /// The reference sequence of the last record (tid, bases)
    current_reference: Option<(i32, std::sync::Arc<fasta::record::Sequence>)>,
    decoding: Decoding,
}

/// What decoding a CRAM record gives.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Decoding {
    /// Every field, with the MD and NM tags htslib adds, as `samtools view`
    /// gives the records.
    #[default]
    Full,
    /// What reading alleles takes: every field but only the read group, `PS`,
    /// `HP` and `BX` of the tags, and no MD or NM. On ONT CRAMs this leaves
    /// out most of the data, the `MM`/`ML` methylation calls.
    Alleles,
}

/// `record` as BAM stores it, as htslib's `cram_to_bam` gives it (MD and NM aside), encoded
/// through `encoded`.
fn to_bam_record(header: &sam::Header, record: &cram::Record<'_>, encoded: &mut Vec<u8>) -> Result<Record> {
    record.encode_bam_record(header, encoded)?;
    let mut rec = Record::from_bam(encoded)?;
    rec.set_bin_from_alignment();
    if rec.is_unmapped() {
        rec.mapq = 0;
    }
    Ok(rec)
}

/// The tags [`Decoding::Alleles`] decodes.
const ALLELE_TAGS: [Tag; 4] = [Tag::READ_GROUP, Tag::new(b'P', b'S'), Tag::new(b'H', b'P'), Tag::new(b'B', b'X')];

/// A FASTA file as reference sequence repository.
pub(crate) fn repository(reference: &Path) -> Result<fasta::Repository> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference)
        .with_context(|| format!("opening reference {}", reference.display()))?;
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader)))
}

impl CramReader {
    /// Opens `path`, decoding with `reference`; a region read with [`Decoding::Alleles`] decodes
    /// on `threads` threads when that is more than one.
    pub fn open(path: &Path, reference: Option<&Path>, threads: usize) -> Result<CramReader> {
        let repository = match reference {
            Some(r) => repository(r)?,
            None => fasta::Repository::default(),
        };
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_path(path)
            .with_context(|| format!("opening {}", path.display()))?;
        let text = {
            let mut header_reader = reader.header_reader();
            header_reader.read_magic_number()?;
            header_reader.read_format_version()?;
            header_reader.read_file_id()?;
            let mut container = header_reader.container_reader()?;
            let mut text = String::new();
            {
                let mut raw = container.raw_sam_header_reader()?;
                raw.read_to_string(&mut text)?;
                raw.discard_to_end()?;
            }
            container.discard_to_end()?;
            text
        };
        let text = match text.find('\0') {
            Some(end) => text[..end].to_string(),
            None => text,
        };
        let sam_header: sam::Header =
            text.parse().with_context(|| format!("parsing the header of {}", path.display()))?;
        let refs = sam_header
            .reference_sequences()
            .iter()
            .map(|(name, rs)| (name.to_string(), usize::from(rs.length()) as u32))
            .collect();
        let first_container = reader.position()?;
        Ok(CramReader {
            path: path.to_path_buf(),
            threads,
            reader,
            sam_header,
            header: Header { text, refs },
            index: None,
            first_container,
            records: Vec::new().into_iter(),
            iter: Iter::Rest,
            encoded: Vec::new(),
            repository,
            has_reference: reference.is_some(),
            shared_reference: false,
            current_reference: None,
            decoding: Decoding::Full,
        })
    }

    /// Sets what decoding gives the records from now on.
    pub fn set_decoding(&mut self, decoding: Decoding) {
        self.decoding = decoding;
    }

    pub fn load_index(&mut self, path: &Path) -> Result<bool> {
        let Some(index_path) = find_index(path, &[".crai"]) else { return Ok(false) };
        self.index = Some(cram::crai::fs::read(index_path)?);
        Ok(true)
    }

    pub fn has_index(&self) -> bool {
        self.index.is_some()
    }

    /// Whether `tid` has records, by the index (`None` without one): a slice holds records of
    /// each reference sequence the index lists it under, so a region read of the whole contig
    /// gives one exactly when this is true, without decoding (and checking the reference of) a
    /// slice.
    pub fn has_records(&self, tid: usize) -> Option<bool> {
        let index = self.index.as_ref()?;
        Some(index.iter().any(|r| r.reference_sequence_id() == Some(tid)))
    }

    /// Decodes with `repository`, the reference this reader was opened with
    /// read through another cache, instead of its own: a caller that also
    /// reads whole chromosomes from it then holds them once. The caller clears
    /// it; fetches no longer do.
    pub fn share_reference(&mut self, repository: fasta::Repository) {
        self.repository = repository;
        self.shared_reference = true;
        self.current_reference = None;
    }

    fn seek(&mut self, offset: u64) -> Result<()> {
        // keep only the reference sequences needed from now on
        if !self.shared_reference {
            self.repository.clear();
        }
        self.current_reference = None;
        self.reader.get_mut().seek(SeekFrom::Start(offset))?;
        self.records = Vec::new().into_iter();
        Ok(())
    }

    pub fn fetch_all(&mut self) -> Result<()> {
        self.seek(self.first_container)?;
        self.iter = Iter::Rest;
        Ok(())
    }

    pub fn fetch(&mut self, tid: usize, beg: i64, end: i64) -> Result<()> {
        if self.index.is_none() {
            bail!("no index available");
        }
        // Stops a threaded read still going.
        self.iter = Iter::Done;
        // An empty interval keeps the serial read's htslib semantics (the reads that start before
        // `beg` and cover it), which a region query cannot express.
        if self.decoding == Decoding::Alleles && self.threads > 1 && end > beg {
            self.iter = Iter::Threaded(Box::new(self.threaded(tid, beg, end)?));
            return Ok(());
        }
        let Some(index) = &self.index else { bail!("no index available") };
        // the first container with records of `tid` that may reach `beg`
        let offset = index
            .iter()
            .filter(|r| {
                r.reference_sequence_id() == Some(tid)
                    && r.alignment_start().is_some_and(|s| usize::from(s) as i64 - 1 + r.alignment_span() as i64 > beg)
            })
            .map(|r| r.offset())
            .min();
        match offset {
            Some(o) => {
                self.seek(o)?;
                self.iter = Iter::Region { tid: tid as i32, beg, end };
            }
            None => self.iter = Iter::Done,
        }
        Ok(())
    }

    /// Starts decoding the records of `tid` overlapping [beg, end) on [`Self::threads`] threads.
    fn threaded(&self, tid: usize, beg: i64, end: i64) -> Result<Threaded> {
        let index = self.index.clone().ok_or_else(|| anyhow::anyhow!("no index available"))?;
        let (name, length) = self.header.refs[tid].clone();
        let start = Position::try_from(usize::try_from(beg.max(0))? + 1)?;
        let region = if end >= i64::from(length) {
            Region::new(name, start..)
        } else {
            Region::new(name, start..=Position::try_from(usize::try_from(end)?)?)
        };
        let workers = NonZeroUsize::new(self.threads).expect("more than one thread");
        let (tx, batches) = mpsc::sync_channel(THREADED_BATCHES_AHEAD);
        let (path, header, repository) = (self.path.clone(), self.sam_header.clone(), self.repository.clone());
        let worker = std::thread::Builder::new().name("cram-region".to_owned()).spawn(move || {
            let run = || -> Result<()> {
                let mut reader = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repository)
                    .set_required_tags(Tags::Only(ALLELE_TAGS.to_vec()))
                    .set_md_nm(MdNm::AsStored)
                    .build_from_path(&path)?;
                reader.fold_queries(
                    &header,
                    &index,
                    &[region],
                    workers,
                    |batch: &mut Vec<Record>, _, record| {
                        let rec = to_bam_record(&header, record, &mut Vec::new()).map_err(io::Error::other)?;
                        batch.push(rec);
                        Ok(())
                    },
                    // The reader is gone when the channel is closed: stop.
                    |batch| tx.send(Ok(batch)).map_err(|_| io::Error::other("the reader stopped")),
                )?;
                Ok(())
            };
            if let Err(e) = run() {
                let _ = tx.send(Err(e));
            }
        })?;
        Ok(Threaded { batches, current: Vec::new().into_iter(), worker: Some(worker) })
    }

    pub fn fetch_unmapped(&mut self) -> Result<()> {
        let Some(index) = &self.index else { bail!("no index available") };
        match index.iter().find(|r| r.reference_sequence_id().is_none()).map(|r| r.offset()) {
            Some(o) => {
                self.seek(o)?;
                self.iter = Iter::Unmapped;
            }
            None => self.iter = Iter::Done,
        }
        Ok(())
    }

    /// Decodes the next container into `records`; false at the end.
    fn next_container(&mut self) -> Result<bool> {
        let mut container = cram::io::reader::Container::default();
        if self.reader.read_container(&mut container)? == 0 {
            return Ok(false);
        }
        let header = &self.sam_header;
        let compression_header = container.compression_header()?;
        let mut out = Vec::new();
        let (md_nm, tags) = match self.decoding {
            Decoding::Full => (MdNm::Regenerate, Tags::All),
            Decoding::Alleles => (MdNm::AsStored, Tags::Only(ALLELE_TAGS.to_vec())),
        };
        let fields = Fields::all();
        for slice in container.slices() {
            let slice = slice?;
            let (core, external) = slice.decode_blocks_for(&compression_header, fields, &tags)?;
            let records = slice.records_with_fields(
                self.repository.clone(),
                header,
                &compression_header,
                &core,
                &external,
                md_nm,
                fields,
                &tags,
            )?;
            for record in records {
                let mut rec = to_bam_record(header, &record, &mut self.encoded)?;
                if self.decoding == Decoding::Full && self.has_reference && rec.tid >= 0 && !rec.is_unmapped() {
                    if self.current_reference.as_ref().is_none_or(|(t, _)| *t != rec.tid) {
                        let name = &self.header.refs[rec.tid as usize].0;
                        let seq = self
                            .repository
                            .get(name.as_bytes())
                            .ok_or_else(|| anyhow::anyhow!("reference sequence {name} not found"))??;
                        self.current_reference = Some((rec.tid, seq));
                    }
                    let (_, seq) = self.current_reference.as_ref().expect("set above");
                    add_md_nm(&mut rec, seq.as_ref().as_ref());
                }
                out.push(rec);
            }
        }
        self.records = out.into_iter();
        Ok(true)
    }

    fn read_record(&mut self) -> Result<Option<Record>> {
        loop {
            if let Some(r) = self.records.next() {
                return Ok(Some(r));
            }
            if !self.next_container()? {
                return Ok(None);
            }
        }
    }

    pub fn read(&mut self) -> Result<Option<Record>> {
        loop {
            match self.iter {
                Iter::Done => return Ok(None),
                Iter::Rest => return self.read_record(),
                Iter::Threaded(ref mut threaded) => {
                    let record = threaded.next();
                    if !matches!(record, Ok(Some(_))) {
                        self.iter = Iter::Done;
                    }
                    return record;
                }
                Iter::Unmapped => match self.read_record()? {
                    None => return Ok(None),
                    Some(r) if r.tid < 0 => return Ok(Some(r)),
                    Some(_) => {}
                },
                Iter::Region { tid, beg, end } => {
                    let Some(rec) = self.read_record()? else {
                        self.iter = Iter::Done;
                        return Ok(None);
                    };
                    if rec.tid != tid {
                        if rec.tid > tid || rec.tid < 0 {
                            self.iter = Iter::Done;
                            return Ok(None);
                        }
                        continue;
                    }
                    if rec.pos() >= end {
                        self.iter = Iter::Done;
                        return Ok(None);
                    }
                    if rec.end_pos() > beg {
                        return Ok(Some(rec));
                    }
                }
            }
        }
    }
}

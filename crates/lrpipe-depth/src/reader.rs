//! Reading the records of one contig from an indexed BAM or CRAM file.

use std::fs::File;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use pasta::bam;
use pasta::bgzf;
use pasta::core::Region;
use pasta::cram as ncram;
use pasta::csi;
use pasta::fasta;
use pasta::sam::{self, alignment::record::data::field::Tag};

use crate::alignment::{Alignment, Cigar};

mod cram;

enum BamIndex {
    Bai(bam::bai::Index),
    Csi(csi::Index),
}

enum Format {
    Bam(BamIndex),
    /// `offsets[tid]`: the file offsets of the containers with records on contig `tid`.
    Cram {
        offsets: Vec<Vec<u64>>,
        reference: PathBuf,
        fai: fasta::fai::Index,
    },
}

fn repository(reference: &Path, fai: &fasta::fai::Index) -> Result<fasta::Repository> {
    let reference = lrpipe_reads::Reference::with_index(reference, fai.clone())?;
    Ok(reference.repository().clone())
}

/// An indexed alignment file. Cheap to share between threads; each thread opens its own
/// [`ContigReader`].
pub struct Alignments {
    path: PathBuf,
    header: sam::Header,
    format: Format,
    contigs: Vec<(String, u64)>,
}

/// htslib looks for `x.bam.bai`, `x.bam.csi`, then `x.bai` and `x.csi`.
fn find_index(path: &Path, exts: &[&str]) -> Option<PathBuf> {
    exts.iter()
        .map(|e| path.with_added_extension(e))
        .chain(exts.iter().map(|e| path.with_extension(e)))
        .find(|p| p.exists())
}

impl Alignments {
    /// Opens `path`; CRAM needs the reference FASTA (with a `.fai` index).
    pub fn open(path: &Path, reference: Option<&Path>) -> Result<Self> {
        let is_cram = path.extension().is_some_and(|e| e.eq_ignore_ascii_case("cram"));
        let (header, format) = if is_cram {
            let Some(reference) = reference else {
                bail!("specify a reference file (or set REF_PATH env var) for decoding CRAM");
            };
            let fai = fasta::fai::fs::read(reference.with_added_extension("fai"))
                .with_context(|| format!("reading the .fai index of {}", reference.display()))?;
            let index_path = find_index(path, &["crai"]).context("alignment file must be indexed")?;
            let index = ncram::crai::fs::read(&index_path)?;
            let mut reader = ncram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository(reference, &fai)?)
                .build_from_path(path)
                .with_context(|| format!("opening {}", path.display()))?;
            let header = reader.read_header()?;
            // The index has a line per slice, in file order; containers hold one or more.
            let mut offsets = vec![Vec::new(); header.reference_sequences().len()];
            for r in index.iter() {
                if let Some(o) = r.reference_sequence_id().and_then(|id| offsets.get_mut(id))
                    && o.last() != Some(&r.offset())
                {
                    o.push(r.offset());
                }
            }
            (header, Format::Cram { offsets, reference: reference.to_owned(), fai })
        } else {
            let index = match find_index(path, &["bai", "csi"]) {
                Some(p) if p.extension().is_some_and(|e| e == "csi") => BamIndex::Csi(csi::fs::read(&p)?),
                Some(p) => BamIndex::Bai(bam::bai::fs::read(&p)?),
                None => bail!("alignment file must be indexed"),
            };
            let mut reader =
                File::open(path).map(bam::io::Reader::new).with_context(|| format!("opening {}", path.display()))?;
            (reader.read_header()?, Format::Bam(index))
        };
        let contigs = header
            .reference_sequences()
            .iter()
            .map(|(name, rs)| (name.to_string(), usize::from(rs.length()) as u64))
            .collect();
        Ok(Alignments { path: path.to_owned(), header, format, contigs })
    }

    /// Reference sequence names and lengths, in header order.
    pub fn contigs(&self) -> &[(String, u64)] {
        &self.contigs
    }

    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Whether the index says contig `tid` has records (CRAM: any container on it).
    pub fn has_records(&self, tid: usize) -> bool {
        // Contigs without records have no bins; metadata, when present, also counts them.
        fn bam_index<I: csi::binning_index::index::reference_sequence::Index>(
            rs: Option<&csi::binning_index::index::ReferenceSequence<I>>,
        ) -> bool {
            use csi::binning_index::ReferenceSequence as _;
            rs.is_some_and(|rs| match rs.metadata() {
                Some(m) => m.mapped_record_count() + m.unmapped_record_count() > 0,
                None => !rs.bins().is_empty(),
            })
        }
        match &self.format {
            Format::Bam(BamIndex::Bai(i)) => bam_index(i.reference_sequences().get(tid)),
            Format::Bam(BamIndex::Csi(i)) => bam_index(i.reference_sequences().get(tid)),
            Format::Cram { offsets, .. } => offsets.get(tid).is_some_and(|o| !o.is_empty()),
        }
    }

    /// A reader of its own. With more than one `decompression_threads`, BAM blocks are
    /// decompressed, or CRAM containers decoded, in parallel.
    pub fn reader(&self, decompression_threads: usize) -> Result<ContigReader<'_>> {
        let inner = match &self.format {
            Format::Bam(_) => {
                let file = File::open(&self.path)?;
                match NonZero::new(decompression_threads).filter(|n| n.get() > 1) {
                    Some(n) => {
                        Inner::BamMt(bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(n, file)))
                    }
                    None => Inner::Bam(bam::io::Reader::new(file)),
                }
            }
            Format::Cram { reference, fai, .. } => {
                Inner::Cram { repository: repository(reference, fai)?, threads: decompression_threads }
            }
        };
        Ok(ContigReader { src: self, inner })
    }
}

enum Inner {
    Bam(bam::io::Reader<bgzf::io::Reader<File>>),
    BamMt(bam::io::Reader<bgzf::io::MultithreadedReader<File>>),
    /// The reference cache is shared by the threads decoding one contig.
    Cram {
        repository: fasta::Repository,
        threads: usize,
    },
}

/// Reads whole contigs from one open file handle.
pub struct ContigReader<'a> {
    src: &'a Alignments,
    inner: Inner,
}

impl ContigReader<'_> {
    /// Calls `f` with every record placed on contig `tid`, mapped or not, in file order: what
    /// htslib returns when mosdepth queries the whole contig. `with_read_group` decodes RG.
    pub fn for_each(&mut self, tid: usize, with_read_group: bool, f: impl FnMut(&Alignment<'_>)) -> Result<()> {
        let src = self.src;
        let ctx = || read_context(src, tid);
        match (&mut self.inner, &src.format) {
            (Inner::Bam(reader), Format::Bam(index)) => read_bam(reader, src, index, tid, with_read_group, f)?,
            (Inner::BamMt(reader), Format::Bam(index)) => read_bam(reader, src, index, tid, with_read_group, f)?,
            (Inner::Cram { repository, threads }, Format::Cram { offsets, .. }) => {
                cram::for_each(&src.path, &src.header, repository, &offsets[tid], tid, with_read_group, *threads, f)
                    .with_context(ctx)?
            }
            _ => unreachable!("reader and file format differ"),
        }
        Ok(())
    }

    /// Clears this reader's cached reference sequences (CRAM), so only the contig being read
    /// stays in memory.
    pub fn release_reference(&self) {
        if let Inner::Cram { repository, .. } = &self.inner {
            repository.clear();
        }
    }
}

fn read_context(src: &Alignments, tid: usize) -> String {
    format!("reading {} from {} (truncated or corrupt file?)", src.contigs[tid].0, src.path.display())
}

fn read_bam<R>(
    reader: &mut bam::io::Reader<R>,
    src: &Alignments,
    index: &BamIndex,
    tid: usize,
    with_read_group: bool,
    mut f: impl FnMut(&Alignment<'_>),
) -> Result<()>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    let region = Region::new(src.contigs[tid].0.as_str(), ..);
    let ctx = || read_context(src, tid);
    let mut query = match index {
        BamIndex::Bai(i) => reader.query(&src.header, i, &region),
        BamIndex::Csi(i) => reader.query(&src.header, i, &region),
    }
    .with_context(ctx)?;
    let mut record = bam::Record::default();
    while query.read_record(&mut record).with_context(ctx)? != 0 {
        let read_group = match with_read_group {
            true => match record.data().get(&Tag::READ_GROUP).transpose().with_context(ctx)? {
                Some(sam::alignment::record::data::field::Value::String(s)) => Some(s.as_ref()),
                _ => None,
            },
            false => None,
        };
        let a = Alignment {
            tid,
            pos: position(record.alignment_start()).with_context(ctx)?,
            flags: record.flags().bits(),
            mapq: record.mapping_quality().map_or(255, u8::from),
            mate_tid: record.mate_reference_sequence_id().transpose().with_context(ctx)?,
            mate_pos: position(record.mate_alignment_start()).with_context(ctx)?,
            template_len: record.template_length(),
            name: record.name().map_or(&b""[..], |n| n.as_ref()),
            cigar: Cigar::Packed(record.cigar().as_bytes()),
            read_group,
        };
        f(&a);
    }
    Ok(())
}

fn position(p: Option<std::io::Result<pasta::core::Position>>) -> std::io::Result<i64> {
    Ok(p.transpose()?.map_or(-1, |p| usize::from(p) as i64 - 1))
}

/// The BAM operation code of a CIGAR operation kind.
fn cigar_code(kind: sam::alignment::record::cigar::op::Kind) -> u8 {
    use crate::alignment::op;
    use sam::alignment::record::cigar::op::Kind;
    match kind {
        Kind::Match => op::MATCH,
        Kind::Insertion => op::INS,
        Kind::Deletion => op::DEL,
        Kind::Skip => op::REF_SKIP,
        Kind::SoftClip => op::SOFT_CLIP,
        Kind::HardClip => op::HARD_CLIP,
        Kind::Pad => op::PAD,
        Kind::SequenceMatch => op::EQUAL,
        Kind::SequenceMismatch => op::DIFF,
    }
}

//! `lrpipe merge-reads`: several sorted, indexed BAMs and CRAMs of one sample merged into one
//! sorted BAM with its index, without samtools: the `merge-reads` stage of `lrpipe run` given
//! several `--reads`.
//!
//! The inputs are streamed, not loaded: a thread for each reads its records contig by contig
//! (CRAM decoded into BAM records), a few MiB ahead, and they are merged by position, k-way, a
//! contig at a time. Records at the same position come forward strand first, as samtools sorts
//! them, then in the order of the inputs, each input's in its own order; so the merge is the
//! same each time, and two files split from one sorted BAM merge back into its order but where
//! two reads of one strand start at one position. That order matters: lrpipe-snv, as Clair3,
//! reads a position's reads in file order, and its qualities move a little with it. The merged records are written
//! in pieces that worker threads compress into BGZF blocks side by side, and the index is built
//! from where each record lands, so the file is not read again to index it.
//!
//! Unplaced unmapped reads (those without a position, at a BAM's end) are left out: no stage
//! reads them. The header is [`crate::input::combine_headers`]'s.

use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, VecDeque};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::{Arc, Mutex};

use lrpipe_reads::{IndexedReads, Reference};
use pasta::bam;
use pasta::bgzf::{self, VirtualPosition};
use pasta::core::Position;
use pasta::csi::binning_index::index::reference_sequence::bin::Chunk;
use pasta::sam::alignment::Record as _;

use crate::input::{self, InputError};
use crate::plan::{BamIndex, Layout};

/// Record bytes an input's reader sends at a time.
const BATCH_BYTES: usize = 1 << 20;

/// Batches an input's reader may be ahead of the merge.
const BATCHES_AHEAD: usize = 4;

/// Uncompressed bytes compressed as one piece by a worker (about 16 BGZF blocks).
const PIECE_BYTES: usize = 1 << 20;

/// Why the reads could not be merged.
#[derive(Debug, thiserror::Error)]
pub enum MergeError {
    /// The inputs cannot be merged.
    #[error(transparent)]
    Input(#[from] InputError),
    /// An input could not be opened or read.
    #[error("reading {}", path.display())]
    Read {
        /// The input.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// An input's records are not sorted by position.
    #[error("{} is not sorted by position: a record at {contig}:{at} follows one at {contig}:{after}", path.display())]
    Unsorted {
        /// The input.
        path: PathBuf,
        /// The contig.
        contig: String,
        /// The out-of-order record's position (1-based).
        at: usize,
        /// The position of the record before it.
        after: usize,
    },
    /// The output is one of the inputs (or an input's index), which writing it would destroy.
    #[error(
        "the merged reads {} would overwrite the input {}; write them to another file (lrpipe run: another \
         --output)",
        output.display(),
        input.display()
    )]
    OutputIsInput {
        /// The merged BAM.
        output: PathBuf,
        /// The input it would overwrite.
        input: PathBuf,
    },
    /// The merged BAM or its index could not be written.
    #[error("writing {}", path.display())]
    Write {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
}

/// What `lrpipe merge-reads` merges, and where to.
#[derive(Clone, Debug)]
pub struct Options {
    /// The sorted, indexed BAMs and CRAMs, in order: records at the same position come in it.
    pub inputs: Vec<PathBuf>,
    /// The merged BAM; its index (`.bai`, or `.csi` for a reference sequence too long for a BAI)
    /// is written beside it.
    pub output: PathBuf,
    /// The reference the CRAM inputs are decoded with.
    pub reference: Option<PathBuf>,
    /// Threads compressing the output (each input is read by a thread of its own besides).
    pub threads: usize,
    /// The command line, for the `@PG` line the merge adds.
    pub command: String,
}

/// What a reader sends the merge.
enum Message {
    /// The next records of the contig.
    Records(VecDeque<bam::Record>),
    /// The contig has no more.
    End,
}

/// Merges `options.inputs` into `options.output`, with its index; returns the records written.
///
/// # Errors
///
/// If the inputs' headers cannot be merged ([`input::combine_headers`]), an input cannot be
/// read or is not sorted by position, a CRAM input is given without a reference, or the output
/// cannot be written.
pub fn merge(options: &Options) -> Result<u64, MergeError> {
    refuse_an_input_as_output(&options.inputs, &options.output)?;
    let mut headers = Vec::with_capacity(options.inputs.len());
    for path in &options.inputs {
        let header = input::read_header(path).map_err(|source| InputError::Header { path: path.clone(), source })?;
        headers.push((path.clone(), header));
    }
    let mut header = input::combine_headers(&headers)?;
    add_own_program(&mut header, &options.command);
    // A reference for each CRAM input, holding the contig its reader is on: one shared by
    // them all could not drop a contig another reader still decodes against.
    let mut references = Vec::with_capacity(options.inputs.len());
    for input in &options.inputs {
        references.push(match &options.reference {
            Some(path) if crate::plan::is_cram(input) => Some(
                Reference::open(path)
                    .map_err(|e| MergeError::Read { path: path.clone(), source: io::Error::other(e) })?,
            ),
            _ => None,
        });
    }
    let contigs: Vec<String> = header.reference_sequences().keys().map(ToString::to_string).collect();
    let longest = header.reference_sequences().values().map(|r| r.length().get()).max().unwrap_or(1);
    let write_error = |path: &Path| {
        let path = path.to_path_buf();
        move |source| MergeError::Write { path, source }
    };
    if let Some(dir) = options.output.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir).map_err(write_error(dir))?;
    }
    // Written under hidden names beside the output, renamed into place once whole: a merge that
    // fails leaves no BAM that looks complete, nor an index of another.
    let staged = Staged::new(&options.output, bam_index_for(longest));
    let result = write_merged(options, &header, &contigs, longest, &references, &staged);
    match result {
        Ok(records) => staged.commit().map(|()| records),
        Err(e) => {
            staged.discard();
            Err(e)
        }
    }
}

/// Writes the merge of `options.inputs` to `staged`'s files; the records written.
fn write_merged(
    options: &Options,
    header: &pasta::sam::Header,
    contigs: &[String],
    longest: usize,
    references: &[Option<Reference>],
    staged: &Staged,
) -> Result<u64, MergeError> {
    let write_error = |path: &Path| {
        let path = path.to_path_buf();
        move |source| MergeError::Write { path, source }
    };
    let file = File::create(&staged.bam).map_err(write_error(&staged.bam))?;

    let (written, index) = std::thread::scope(|scope| -> Result<_, MergeError> {
        let mut receivers = Vec::with_capacity(options.inputs.len());
        let mut readers = Vec::with_capacity(options.inputs.len());
        for (path, reference) in options.inputs.iter().zip(references) {
            let (send, receive) = sync_channel(BATCHES_AHEAD);
            receivers.push(receive);
            let reference = reference.as_ref();
            readers.push(scope.spawn(move || read_input(path, reference, contigs, &send)));
        }
        let (pieces, results) = compressors(scope, options.threads.max(1));
        let writer = scope.spawn(move || write_pieces(BufWriter::new(file), header, longest, &results));
        let merged = merge_records(&options.inputs, contigs, &receivers, &pieces);
        // Ends the readers waiting to send (when the merge failed) and the compressors.
        drop(receivers);
        drop(pieces);
        let written = writer.join().expect("the writer does not panic").map_err(write_error(&staged.bam));
        // A reader's failure, or the writer's, stops the merge: theirs is the error to report.
        for (path, reader) in options.inputs.iter().zip(readers) {
            reader.join().expect("a reader does not panic").map_err(|e| match e {
                ReadFailure::Unsorted { contig, at, after } => {
                    MergeError::Unsorted { path: path.clone(), contig, at, after }
                }
                ReadFailure::Io(source) => MergeError::Read { path: path.clone(), source },
            })?;
        }
        let index = written?;
        Ok((merged?, index))
    })?;
    match index {
        bam::Index::Bai(index) => bam::bai::fs::write(&staged.index, &index),
        bam::Index::Csi(index) => pasta::csi::fs::write(&staged.index, &index),
    }
    .map_err(write_error(&staged.index))?;
    Ok(written)
}

/// Adds the merge's own `@PG` line, as htslib's `sam_hdr_add_pg` does: one after the end of
/// each chain of the inputs' (`lrpipe-merge-reads`, then `lrpipe-merge-reads.1`, ...), or one
/// after the last line when no line follows another, so the merged BAM says how it was made.
fn add_own_program(header: &mut pasta::sam::Header, command: &str) {
    use pasta::sam::header::record::value::{Map, map::Program, map::program::tag};
    let program = Map::<Program>::builder()
        .insert(tag::NAME, "lrpipe")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(tag::COMMAND_LINE, command)
        .build();
    if let Ok(program) = program
        && header.add_program_record("lrpipe-merge-reads", program.clone()).is_err()
    {
        // Only a `PP` of the line's own that is not in the header fails; it has none.
        header.programs_mut().as_mut().insert("lrpipe-merge-reads".into(), program);
    }
}

/// The index a merged BAM gets for reference sequences of `longest` bases at most: a BAI where
/// one can index them, else a CSI. pasta's indexer decides, so the plan names the index the
/// merge writes.
#[must_use]
pub fn bam_index_for(longest: usize) -> BamIndex {
    let index = Position::new(longest.max(1))
        .and_then(|p| bam::index::Indexer::builder().set_max_position_hint(p).build().ok())
        .map(|indexer| indexer.build(0));
    match index {
        Some(bam::Index::Bai(_)) => BamIndex::Bai,
        _ => BamIndex::Csi,
    }
}

/// The merged BAM and its index while they are written: hidden, beside the output.
struct Staged {
    bam: PathBuf,
    index: PathBuf,
    output: PathBuf,
    output_index: PathBuf,
}

impl Staged {
    fn new(output: &Path, index: BamIndex) -> Self {
        let name = output.file_name().map_or_else(Default::default, |n| n.to_string_lossy().into_owned());
        let bam = output.with_file_name(format!(".{name}.merging"));
        Staged {
            index: Layout::index(&bam, index.extension()),
            bam,
            output: output.to_path_buf(),
            output_index: Layout::index(output, index.extension()),
        }
    }

    /// Renames the index, then the BAM, into place.
    fn commit(&self) -> Result<(), MergeError> {
        for (from, to) in [(&self.index, &self.output_index), (&self.bam, &self.output)] {
            std::fs::rename(from, to).map_err(|source| MergeError::Write { path: to.clone(), source })?;
        }
        Ok(())
    }

    /// Removes what was written.
    fn discard(&self) {
        let _ = std::fs::remove_file(&self.bam);
        let _ = std::fs::remove_file(&self.index);
    }
}

/// Refuses an `output` that is one of the `inputs` or an input's index (or would be written over
/// one, as its own index), once paths are made absolute with links resolved: `lrpipe run` checks
/// it when planning, and the merge again.
///
/// # Errors
///
/// [`MergeError::OutputIsInput`].
pub fn refuse_an_input_as_output(inputs: &[PathBuf], output: &Path) -> Result<(), MergeError> {
    let written: Vec<PathBuf> = [output.to_path_buf(), Layout::index(output, "bai"), Layout::index(output, "csi")]
        .iter()
        .map(|p| resolved(p))
        .collect();
    for input in inputs {
        let read = std::iter::once(input.clone()).chain(crate::plan::index_of(input));
        if read.map(|p| resolved(&p)).any(|p| written.contains(&p)) {
            return Err(MergeError::OutputIsInput { output: output.to_path_buf(), input: input.clone() });
        }
    }
    Ok(())
}

/// `path` absolute, with links resolved as far as it exists.
fn resolved(path: &Path) -> PathBuf {
    let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    absolute.canonicalize().unwrap_or_else(|_| crate::provenance::canonical(&absolute))
}

/// Why an input's reader stopped.
enum ReadFailure {
    Io(io::Error),
    Unsorted { contig: String, at: usize, after: usize },
}

/// Sends `path`'s records to the merge, contig by contig (each in batches, then [`Message::End`]),
/// checking they are sorted. Stops early, without an error, when the merge stops taking them.
/// A CRAM's `reference` (its own) holds only the contig being read: each is dropped as the next
/// begins, as lrpipe-sv's pass does, or a whole genome's worth would stay in memory.
fn read_input(
    path: &Path,
    reference: Option<&Reference>,
    contigs: &[String],
    send: &SyncSender<Message>,
) -> Result<(), ReadFailure> {
    let mut file = IndexedReads::open(path, reference).map_err(|e| ReadFailure::Io(io::Error::other(e)))?;
    for contig in contigs {
        if let Some(reference) = reference {
            reference.clear();
        }
        let mut batch = VecDeque::new();
        let mut bytes = 0;
        let mut last = 0;
        let mut unsorted = None;
        let mut gone = false;
        let read = file.query_contig(contig, |record| {
            let at = position(record);
            if at < last {
                unsorted = Some((at, last));
                return Err(io::Error::other("unsorted"));
            }
            last = at;
            bytes += record.as_ref().len();
            batch.push_back(record.clone());
            if bytes >= BATCH_BYTES {
                bytes = 0;
                if send.send(Message::Records(std::mem::take(&mut batch))).is_err() {
                    gone = true;
                    return Err(io::Error::other("the merge stopped"));
                }
            }
            Ok(())
        });
        if let Some((at, after)) = unsorted {
            return Err(ReadFailure::Unsorted { contig: contig.clone(), at, after });
        }
        if gone {
            return Ok(());
        }
        read.map_err(ReadFailure::Io)?;
        if (!batch.is_empty() && send.send(Message::Records(batch)).is_err()) || send.send(Message::End).is_err() {
            return Ok(());
        }
    }
    Ok(())
}

/// A record's 1-based start (0 for none), by which inputs must be sorted.
fn position(record: &bam::Record) -> usize {
    record.alignment_start().and_then(Result::ok).map_or(0, usize::from)
}

/// What records are merged by: the start, then the strand (forward first), then the input.
fn key(record: &bam::Record, input: usize) -> (usize, bool, usize) {
    (position(record), record.flags().is_reverse_complemented(), input)
}

/// A piece of the merged BAM for a worker to compress: its records' bytes as BAM stores them
/// (each after its length), and each record's end in them with what the index needs of it.
struct Piece {
    number: usize,
    bytes: Vec<u8>,
    records: Vec<(usize, Option<Context>)>,
}

/// What a BAM index takes of a record: its reference sequence, start, end and whether it is
/// mapped.
type Context = (usize, Position, Position, bool);

/// A compressed [`Piece`]: its BGZF blocks (without the end-of-file block), and each record's
/// start and end in them.
struct Compressed {
    number: usize,
    blocks: Vec<u8>,
    chunks: Vec<(VirtualPosition, VirtualPosition, Option<Context>)>,
}

/// Merges the inputs' records, contig by contig, into [`Piece`]s for the compressors; the
/// records merged. A reader that stops before the end (it failed: the caller reports why)
/// stops the merge.
fn merge_records(
    inputs: &[PathBuf],
    contigs: &[String],
    receivers: &[Receiver<Message>],
    pieces: &SyncSender<Piece>,
) -> Result<u64, MergeError> {
    let stopped =
        |i: usize| MergeError::Read { path: inputs[i].clone(), source: io::Error::other("its reader stopped") };
    let mut piece = Piece { number: 0, bytes: Vec::with_capacity(PIECE_BYTES + BATCH_BYTES), records: Vec::new() };
    let mut written = 0;
    let mut heads: Vec<VecDeque<bam::Record>> = vec![VecDeque::new(); receivers.len()];
    for _ in contigs {
        // The input whose next record comes first, by (start, reverse strand, input).
        let mut heap = BinaryHeap::new();
        let next = |i: usize, head: &mut VecDeque<bam::Record>| -> Result<bool, MergeError> {
            loop {
                match receivers[i].recv() {
                    Ok(Message::Records(records)) if records.is_empty() => {}
                    Ok(Message::Records(records)) => {
                        *head = records;
                        return Ok(true);
                    }
                    Ok(Message::End) => return Ok(false),
                    Err(_) => return Err(stopped(i)),
                }
            }
        };
        for (i, head) in heads.iter_mut().enumerate() {
            if next(i, head)? {
                heap.push(Reverse(key(&head[0], i)));
            }
        }
        while let Some(Reverse((_, _, i))) = heap.pop() {
            let record = heads[i].pop_front().expect("an input in the heap has a record");
            let bytes = record.as_ref();
            let length = u32::try_from(bytes.len()).map_err(|e| MergeError::Read {
                path: inputs[i].clone(),
                source: io::Error::new(io::ErrorKind::InvalidData, e),
            })?;
            piece.bytes.extend(length.to_le_bytes());
            piece.bytes.extend_from_slice(bytes);
            piece.records.push((piece.bytes.len(), context(&record)));
            written += 1;
            if piece.bytes.len() >= PIECE_BYTES {
                let number = piece.number + 1;
                let full = std::mem::replace(
                    &mut piece,
                    Piece { number, bytes: Vec::with_capacity(PIECE_BYTES + BATCH_BYTES), records: Vec::new() },
                );
                if pieces.send(full).is_err() {
                    return Err(MergeError::Write { path: PathBuf::new(), source: io::Error::other("stopped") });
                }
            }
            if heads[i].is_empty() && !next(i, &mut heads[i])? {
                continue;
            }
            heap.push(Reverse(key(&heads[i][0], i)));
        }
    }
    if !piece.records.is_empty() {
        let _ = pieces.send(piece);
    }
    Ok(written)
}

/// The record's reference sequence, start and end, and whether it is mapped, for the index;
/// none for a record without them.
fn context(record: &bam::Record) -> Option<Context> {
    let id = record.reference_sequence_id()?.ok()?;
    let start = record.alignment_start()?.ok()?;
    let end = record.alignment_end()?.ok()?;
    Some((id, start, end, !record.flags().is_unmapped()))
}

/// `threads` workers compressing [`Piece`]s: the sender to give them pieces, and the receiver of
/// what they compress, in any order.
fn compressors<'scope>(
    scope: &'scope std::thread::Scope<'scope, '_>,
    threads: usize,
) -> (SyncSender<Piece>, Receiver<io::Result<Compressed>>) {
    let (pieces, work) = sync_channel::<Piece>(threads);
    let (done, results) = sync_channel(threads * 2);
    let work = Arc::new(Mutex::new(work));
    for _ in 0..threads {
        let (work, done) = (Arc::clone(&work), done.clone());
        scope.spawn(move || {
            loop {
                let piece = work.lock().unwrap_or_else(std::sync::PoisonError::into_inner).recv();
                let Ok(piece) = piece else { return };
                if done.send(compress(piece)).is_err() {
                    return;
                }
            }
        });
    }
    (pieces, results)
}

/// Compresses a piece into BGZF blocks, noting where each record starts and ends in them.
fn compress(piece: Piece) -> io::Result<Compressed> {
    let mut w = bgzf::io::Writer::new(Vec::with_capacity(piece.bytes.len() / 2));
    let mut chunks = Vec::with_capacity(piece.records.len());
    let mut start = 0;
    for (end, context) in piece.records {
        let before = w.virtual_position();
        w.write_all(&piece.bytes[start..end])?;
        chunks.push((before, w.virtual_position(), context));
        start = end;
    }
    w.flush()?;
    Ok(Compressed { number: piece.number, blocks: w.into_inner(), chunks })
}

/// Writes the BAM: its header, the compressed pieces in order as they come, and the end-of-file
/// block; and indexes the records from where each landed. `longest` is the longest reference
/// sequence, which decides between a BAI and a CSI.
fn write_pieces(
    mut out: BufWriter<File>,
    header: &pasta::sam::Header,
    longest: usize,
    results: &Receiver<io::Result<Compressed>>,
) -> io::Result<bam::Index> {
    // As [`bam_index_for`] decides.
    let mut indexer = bam::index::Indexer::builder()
        .set_max_position_hint(Position::new(longest.max(1)).expect("not 0"))
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let mut head = bam::io::Writer::from(bgzf::io::Writer::new(Vec::new()));
    head.write_header(header)?;
    let mut head = head.into_inner();
    head.flush()?;
    let head = head.into_inner();
    out.write_all(&head)?;
    let mut offset = head.len() as u64;
    let mut waiting = BTreeMap::new();
    let mut next = 0;
    for compressed in results {
        let compressed = compressed?;
        waiting.insert(compressed.number, compressed);
        while let Some(piece) = waiting.remove(&next) {
            let at = |v: VirtualPosition| {
                VirtualPosition::try_from((v.compressed() + offset, v.uncompressed())).map_err(io::Error::other)
            };
            for (start, end, context) in piece.chunks {
                indexer.add_record(context, Chunk::new(at(start)?, at(end)?)).map_err(io::Error::other)?;
            }
            out.write_all(&piece.blocks)?;
            offset += piece.blocks.len() as u64;
            next += 1;
        }
    }
    if !waiting.is_empty() {
        return Err(io::Error::other("a piece of the output was lost"));
    }
    out.write_all(&bgzf::io::Writer::new(Vec::new()).finish()?)?;
    out.into_inner().map_err(io::IntoInnerError::into_error)?.sync_all()?;
    Ok(indexer.build(header.reference_sequences().len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    use pasta::sam;

    /// A CRAM of one record on each of `contigs` (each 1,000 `ACGT`s), against `reference.fa`
    /// written beside it with its `.fai`; and its index. Each record spans over half its contig,
    /// so decoding it reads the whole contig, which the reference keeps (a shorter slice's span
    /// is read alone, and not kept).
    fn small_cram(dir: &Path, contigs: &[&str]) -> (PathBuf, PathBuf) {
        use pasta::sam::alignment::io::Write as _;
        use pasta::sam::alignment::record::cigar::{Op, op::Kind};
        use pasta::sam::alignment::record::{Flags, MappingQuality};
        use pasta::sam::alignment::record_buf::{Cigar, RecordBuf, Sequence};
        use pasta::sam::header::record::value::{Map, map::ReferenceSequence};
        let seq = b"ACGT".repeat(250);
        let (mut fasta, mut fai) = (String::new(), String::new());
        for name in contigs {
            fasta.push_str(&format!(">{name}\n"));
            fai.push_str(&format!("{name}\t1000\t{}\t1000\t1001\n", fasta.len()));
            fasta.push_str(std::str::from_utf8(&seq).expect("ASCII"));
            fasta.push('\n');
        }
        let reference = dir.join("reference.fa");
        std::fs::write(&reference, fasta).expect("written");
        std::fs::write(dir.join("reference.fa.fai"), fai).expect("written");
        let mut header = sam::Header::builder();
        for name in contigs {
            header =
                header.add_reference_sequence(*name, Map::<ReferenceSequence>::new(1000.try_into().expect("1000")));
        }
        let header = header.build();
        let cram = dir.join("reads.cram");
        let opened = Reference::open(&reference).expect("the reference");
        // A slice for each record, so each contig's is decoded alone.
        let one_per_slice =
            pasta::cram::io::writer::compression::Options { records_per_slice: 1, ..Default::default() };
        let mut writer = pasta::cram::io::writer::Builder::default()
            .set_compression_options(one_per_slice)
            .set_reference_sequence_repository(opened.repository().clone())
            .build_from_path(&cram)
            .expect("created");
        writer.write_alignment_header(&header).expect("header");
        for id in 0..contigs.len() {
            let record = RecordBuf::builder()
                .set_name(format!("r{id}"))
                .set_flags(Flags::empty())
                .set_reference_sequence_id(id)
                .set_alignment_start(Position::new(11).expect("non-zero"))
                .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
                .set_cigar([Op::new(Kind::Match, 600)].into_iter().collect::<Cigar>())
                .set_sequence(Sequence::from(seq[10..610].to_vec()))
                .build();
            writer.write_alignment_record(&header, &record).expect("record");
        }
        writer.finish(&header).expect("finished");
        eprintln!("finished");
        drop(writer);
        let index = pasta::cram::fs::index(&cram).expect("indexed");
        pasta::cram::crai::fs::write(crate::plan::Layout::index(&cram, "crai"), &index).expect("index written");
        (cram, reference)
    }

    /// The merge's own `@PG` line follows each input's chain, as samtools 1.24's does: after the
    /// last line when none follows another, else one after each chain's end.
    #[test]
    fn the_merge_adds_its_program_after_each_chain() {
        let added = |lines: &str| -> Vec<(String, String)> {
            let text = format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n{lines}");
            let mut header: pasta::sam::Header = text.parse().expect("a header");
            add_own_program(&mut header, "lrpipe merge-reads --output m.bam x.bam y.bam");
            let tag = pasta::sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;
            let programs = header.programs();
            let ours = programs.as_ref().iter().filter(|(id, _)| id.starts_with(b"lrpipe-merge-reads"));
            ours.map(|(id, p)| {
                (id.to_string(), p.other_fields().get(&tag).map(ToString::to_string).unwrap_or_default())
            })
            .collect()
        };
        let pair = |id: &str, pp: &str| (id.to_string(), pp.to_string());
        assert_eq!(added("@PG\tID:a\tPN:a\n@PG\tID:b\tPN:b\n"), [pair("lrpipe-merge-reads", "b")]);
        let chains = "@PG\tID:a\tPN:a\n@PG\tID:a2\tPN:a2\tPP:a\n@PG\tID:b\tPN:b\n@PG\tID:b2\tPN:b2\tPP:b\n";
        assert_eq!(added(chains), [pair("lrpipe-merge-reads", "a2"), pair("lrpipe-merge-reads.1", "b2")]);
    }

    /// The BAI or CSI choice is pasta's indexer's: a BAI indexes positions up to 2^29.
    #[test]
    fn a_bai_indexes_up_to_2_29_bases() {
        assert_eq!(bam_index_for(248_956_422), BamIndex::Bai);
        assert_eq!(bam_index_for(1 << 29), BamIndex::Bai);
        assert_eq!(bam_index_for((1 << 29) + 1), BamIndex::Csi);
    }

    /// A CRAM input's reference holds only the contig being read: the one before is dropped
    /// as the next begins, so a merge does not end holding the whole genome.
    #[test]
    fn a_cram_readers_reference_holds_one_contig() {
        let dir = std::env::temp_dir().join(format!("lrpipe-merge-cram-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let contigs = ["chr1", "chr2", "chr3"];
        let (cram, reference) = small_cram(&dir, &contigs);
        eprintln!("written");
        let reference = Reference::open(&reference).expect("the reference");
        let names: Vec<String> = contigs.iter().map(ToString::to_string).collect();
        // A rendezvous: the reader waits in each send until the test takes it.
        let (send, receive) = sync_channel(0);
        let (mut records, mut held) = (Vec::new(), Vec::new());
        std::thread::scope(|scope| {
            let (cram, reference, names) = (&cram, &reference, &names);
            // The sender goes with the reader, so the channel closes when it is done.
            let reader = scope.spawn(move || read_input(cram, Some(reference), names, &send));
            while let Ok(message) = receive.recv() {
                if let Message::Records(batch) = message {
                    records.push(batch.len());
                    // The reader has decoded the contig and waits to send its end.
                    held.push(reference.cached());
                }
            }
            assert!(reader.join().expect("no panic").is_ok());
        });
        assert_eq!(records, [1, 1, 1], "a record on each contig");
        assert_eq!(held, [1, 1, 1], "on each contig, the earlier ones are dropped");
        std::fs::remove_dir_all(&dir).expect("removed");
    }
}

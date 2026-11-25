//! CRAM records of one contig, decoded a container per thread and handed out in file order.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc;

use anyhow::{Result, anyhow};
use pasta::cram;
use pasta::cram::io::reader::{Fields, MdNm, Tags};
use pasta::fasta;
use pasta::sam::{
    self,
    alignment::{Record as _, record::data::field::Tag},
};

use crate::alignment::{Alignment, Cigar};

/// One record's fixed fields; variable-length ones live in the batch's buffers.
struct Fixed {
    pos: i64,
    flags: u16,
    mapq: u8,
    mate_tid: Option<usize>,
    mate_pos: i64,
    template_len: i32,
    name_end: usize,
    cigar_end: usize,
    /// End of the RG value in `rgs`, if the record has one.
    rg_end: Option<usize>,
}

/// The records of one container that are on the contig being read.
#[derive(Default)]
struct Batch {
    recs: Vec<Fixed>,
    names: Vec<u8>,
    cigars: Vec<(u8, u32)>,
    rgs: Vec<u8>,
}

impl Batch {
    fn for_each(&self, tid: usize, f: &mut impl FnMut(&Alignment<'_>)) {
        let (mut name_start, mut cigar_start, mut rg_start) = (0, 0, 0);
        for r in &self.recs {
            let read_group = r.rg_end.map(|end| {
                let rg = &self.rgs[rg_start..end];
                rg_start = end;
                rg
            });
            f(&Alignment {
                tid,
                pos: r.pos,
                flags: r.flags,
                mapq: r.mapq,
                mate_tid: r.mate_tid,
                mate_pos: r.mate_pos,
                template_len: r.template_len,
                name: &self.names[name_start..r.name_end],
                cigar: Cigar::Ops(&self.cigars[cigar_start..r.cigar_end]),
                read_group,
            });
            name_start = r.name_end;
            cigar_start = r.cigar_end;
        }
    }
}

fn position(p: Option<io::Result<pasta::core::Position>>) -> io::Result<i64> {
    Ok(p.transpose()?.map_or(-1, |p| usize::from(p) as i64 - 1))
}

/// Decodes the container at `offset`, keeping records on contig `tid` (the records an htslib
/// query of the whole contig returns).
fn decode_container(
    reader: &mut cram::io::Reader<File>,
    offset: u64,
    header: &sam::Header,
    repository: &fasta::Repository,
    tid: usize,
    with_read_group: bool,
) -> io::Result<Batch> {
    reader.seek(SeekFrom::Start(offset))?;
    let mut container = cram::io::reader::Container::default();
    reader.read_container(&mut container)?;
    let compression_header = container.compression_header()?;
    // Only what depth reads: no bases (so no reference), qualities or tags but the read group.
    let mut fields = Fields::NAME
        | Fields::FLAGS
        | Fields::REFERENCE_SEQUENCE_ID
        | Fields::ALIGNMENT_START
        | Fields::MAPPING_QUALITY
        | Fields::CIGAR
        | Fields::MATE_REFERENCE_SEQUENCE_ID
        | Fields::MATE_ALIGNMENT_START
        | Fields::TEMPLATE_LENGTH;
    if with_read_group {
        fields |= Fields::READ_GROUP | Fields::DATA;
    }
    let tags = Tags::Only(vec![Tag::READ_GROUP]);
    let mut batch = Batch::default();
    for slice in container.slices() {
        let slice = slice?;
        let (core, external) = slice.decode_blocks_for(&compression_header, fields, &tags)?;
        let records = slice.records_with_fields(
            repository.clone(),
            header,
            &compression_header,
            &core,
            &external,
            MdNm::AsStored,
            fields,
            &tags,
        )?;
        for record in records {
            if record.reference_sequence_id(header).transpose()? != Some(tid) {
                continue;
            }
            if let Some(name) = record.name() {
                batch.names.extend_from_slice(name);
            }
            for op in record.cigar().iter() {
                let op = op?;
                batch.cigars.push((super::cigar_code(op.kind()), op.len() as u32));
            }
            let mut rg_end = None;
            if with_read_group
                && let Some(value) = record.data().get(&Tag::READ_GROUP).transpose()?
                && let sam::alignment::record::data::field::Value::String(s) = value
            {
                batch.rgs.extend_from_slice(s);
                rg_end = Some(batch.rgs.len());
            }
            let flags = record.flags()?.bits();
            batch.recs.push(Fixed {
                pos: position(record.alignment_start())?,
                flags,
                // CRAM stores no mapping quality for unmapped reads; htslib reports 0.
                mapq: match record.mapping_quality().transpose()? {
                    Some(q) => u8::from(q),
                    None if flags & crate::alignment::flag::UNMAPPED != 0 => 0,
                    None => 255,
                },
                mate_tid: record.mate_reference_sequence_id(header).transpose()?,
                mate_pos: position(record.mate_alignment_start())?,
                template_len: record.template_length()?,
                name_end: batch.names.len(),
                cigar_end: batch.cigars.len(),
                rg_end,
            });
        }
    }
    Ok(batch)
}

/// Calls `f` with the records on contig `tid`, in file order, decoding up to `threads`
/// containers at once. `offsets` are the contig's container offsets from the CRAI index.
#[allow(clippy::too_many_arguments)]
pub(super) fn for_each(
    path: &Path,
    header: &sam::Header,
    repository: &fasta::Repository,
    offsets: &[u64],
    tid: usize,
    with_read_group: bool,
    threads: usize,
    mut f: impl FnMut(&Alignment<'_>),
) -> Result<()> {
    let threads = threads.clamp(1, offsets.len().max(1));
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let (tx, rx) = mpsc::channel::<(usize, io::Result<Batch>)>();
    std::thread::scope(|s| -> Result<()> {
        for _ in 0..threads {
            let tx = tx.clone();
            let (next, failed) = (&next, &failed);
            // pasta's CRAM decoding needs more than the default 2 MiB of stack.
            std::thread::Builder::new().stack_size(64 << 20).spawn_scoped(s, move || {
                let mut reader = match File::open(path) {
                    Ok(file) => cram::io::Reader::new(file),
                    Err(e) => {
                        let _ = tx.send((usize::MAX, Err(e)));
                        return;
                    }
                };
                while !failed.load(Ordering::Relaxed) {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(&offset) = offsets.get(i) else { break };
                    let batch = decode_container(&mut reader, offset, header, repository, tid, with_read_group);
                    if batch.is_err() {
                        failed.store(true, Ordering::Relaxed);
                    }
                    if tx.send((i, batch)).is_err() {
                        break;
                    }
                }
            })?;
        }
        drop(tx);
        let mut pending = BTreeMap::new();
        let mut done = 0;
        for (i, batch) in rx {
            let batch = match batch {
                Ok(b) => b,
                Err(e) => {
                    failed.store(true, Ordering::Relaxed);
                    let at = offsets.get(i).map_or(String::new(), |o| format!(" at byte {o}"));
                    return Err(anyhow::Error::new(e).context(format!("decoding a CRAM container{at}")));
                }
            };
            pending.insert(i, batch);
            while let Some(b) = pending.remove(&done) {
                b.for_each(tid, &mut f);
                done += 1;
            }
        }
        if done != offsets.len() {
            return Err(anyhow!("internal error: {done} of {} containers decoded", offsets.len()));
        }
        Ok(())
    })
}

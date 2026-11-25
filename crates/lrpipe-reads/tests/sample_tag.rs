//! `IndexedReads::sample_tag` on small BAMs and CRAMs: haplotagged or not, on one contig of
//! several or only away from a contig's start, a short contig listed first, unmapped records,
//! tags of another type, and the limits.

use std::fs;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use lrpipe_reads::{IndexedReads, Reference, TagSample};
use pasta::bam;
use pasta::core::Position;
use pasta::cram;
use pasta::fasta;
use pasta::sam::{
    self,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::{Cigar, Sequence, data::field::Value},
    },
    header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    },
};

const HP: [u8; 2] = *b"HP";
/// A contig shorter than `LONG_CONTIG`, listed first, and two long ones.
const CONTIGS: [(&str, usize); 3] = [("chrM", 1_000), ("chr1", 1_200_000), ("chr2", 1_200_000)];
/// The long contigs' middle (0-based).
const MID: usize = 600_000;
const LEN: usize = 100;

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-reads-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// The reference's bases at `start` (0-based) of contig `id`.
fn bases(id: usize, start: usize) -> Vec<u8> {
    (start..start + LEN).map(|i| b"ACGT"[(i * 7 + id) % 4]).collect()
}

/// The HP tag of a read.
#[derive(Clone, Copy)]
enum Hp {
    None,
    Int(u8),
    Text(&'static str),
}

/// A read on contig `id` at `start` (0-based), tagged `hp`.
fn read(name: &str, id: usize, start: usize, hp: Hp, unmapped: bool) -> RecordBuf {
    let mut flags = Flags::empty();
    if unmapped {
        flags |= Flags::UNMAPPED;
    }
    let mut builder = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(id)
        .set_alignment_start(Position::new(start + 1).expect("non-zero"))
        .set_sequence(Sequence::from(bases(id, start)));
    if !unmapped {
        builder = builder
            .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
            .set_cigar([Op::new(Kind::Match, LEN)].into_iter().collect::<Cigar>());
    }
    let value = match hp {
        Hp::None => None,
        Hp::Int(h) => Some(Value::from(h)),
        Hp::Text(h) => Some(Value::String(h.into())),
    };
    if let Some(value) = value {
        builder = builder.set_data([(Tag::new(b'H', b'P'), value)].into_iter().collect());
    }
    builder.build()
}

/// `reads.bam` and `reads.cram` (both indexed) of `records`, and the reference `ref.fa`.
fn write(dir: &Path, records: &[RecordBuf]) -> Result<(), Box<dyn std::error::Error>> {
    let (mut text, mut fai) = (String::new(), String::new());
    for (id, (name, len)) in CONTIGS.iter().enumerate() {
        text.push_str(&format!(">{name}\n"));
        fai.push_str(&format!("{name}\t{len}\t{}\t{len}\t{}\n", text.len(), len + 1));
        let seq: Vec<u8> = (0..*len).map(|i| b"ACGT"[(i * 7 + id) % 4]).collect();
        text.push_str(std::str::from_utf8(&seq)?);
        text.push('\n');
    }
    fs::write(dir.join("ref.fa"), text)?;
    fs::write(dir.join("ref.fa.fai"), fai)?;
    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?;
    let mut header = sam::Header::builder().set_header(hd);
    for (name, len) in CONTIGS {
        header =
            header.add_reference_sequence(name, Map::<ReferenceSequence>::new(NonZeroUsize::new(len).expect("long")));
    }
    let header = header.build();

    let path = dir.join("reads.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&path)?);
    writer.write_header(&header)?;
    for record in records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish()?;
    drop(writer);
    match bam::fs::index(&path)? {
        bam::Index::Bai(index) => bam::bai::fs::write(dir.join("reads.bam.bai"), &index)?,
        bam::Index::Csi(index) => pasta::csi::fs::write(dir.join("reads.bam.csi"), &index)?,
    }

    let fasta = fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("ref.fa"))?;
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta));
    let path = dir.join("reads.cram");
    let mut writer =
        cram::io::writer::Builder::default().set_reference_sequence_repository(repository).build_from_path(&path)?;
    writer.write_header(&header)?;
    for record in records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish(&header)?;
    drop(writer);
    cram::crai::fs::write(dir.join("reads.cram.crai"), &cram::fs::index(&path)?)?;
    Ok(())
}

/// `sample_tag(HP, per_contig, contigs)` on the BAM and on the CRAM, which must agree.
fn sample(dir: &Path, per_contig: usize, contigs: usize) -> TagSample {
    let reference = Reference::open(&dir.join("ref.fa")).expect("reference");
    let mut bam = IndexedReads::open(&dir.join("reads.bam"), None).expect("BAM");
    let found = bam.sample_tag(HP, per_contig, contigs).expect("sampled");
    let mut cram = IndexedReads::open(&dir.join("reads.cram"), Some(&reference)).expect("CRAM");
    assert_eq!(cram.sample_tag(HP, per_contig, contigs).expect("sampled"), found, "the CRAM's is the BAM's");
    found
}

/// `n` reads named `prefix` on contig `id` from `start`, 100 bases apart, tagged `hp`.
fn reads(prefix: &str, id: usize, start: usize, n: usize, hp: Hp) -> Vec<RecordBuf> {
    (0..n).map(|i| read(&format!("{prefix}{i}"), id, start + i * 100, hp, false)).collect()
}

#[test]
fn a_contig_tagged_after_an_untagged_one_is_found() {
    let dir = temp_dir("second-contig");
    // chrM's reads come first in the file but are untagged, as haplotagging leaves them; chr1
    // is untagged (no phased sites, say) and chr2 tagged.
    let mut records = reads("m", 0, 0, 5, Hp::None);
    records.extend(reads("a", 1, 0, 5, Hp::None));
    records.extend(reads("b", 1, MID, 5, Hp::None));
    records.push(read("unmapped", 1, MID + 600, Hp::None, true));
    records.extend(reads("c", 2, MID, 3, Hp::None));
    records.extend(reads("t", 2, MID + 300, 2, Hp::Int(1)));
    write(&dir, &records).expect("inputs");
    // chr1's 10 (its middle on, then its start; the unmapped read not counted), then chr2's to
    // its first tagged read.
    assert_eq!(sample(&dir, 100, 8), TagSample { records: 14, found: true, other_type: 0 });
    // Only chr1.
    assert_eq!(sample(&dir, 100, 1), TagSample { records: 10, found: false, other_type: 0 });
    // Three of each long contig, then chrM's.
    assert_eq!(sample(&dir, 3, 8), TagSample { records: 9, found: false, other_type: 0 });
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_contig_is_sampled_from_its_middle() {
    let dir = temp_dir("middle");
    // Untagged at chr1's start, tagged from its middle on.
    let mut records = reads("a", 1, 0, 20, Hp::None);
    records.extend(reads("t", 1, MID, 2, Hp::Int(2)));
    write(&dir, &records).expect("inputs");
    assert_eq!(sample(&dir, 5, 8), TagSample { records: 1, found: true, other_type: 0 });

    // Reads only in a contig's first half are still sampled.
    let mut records = reads("a", 1, 0, 3, Hp::None);
    records.push(read("t", 1, 400, Hp::Int(1), false));
    write(&dir, &records).expect("inputs");
    assert_eq!(sample(&dir, 5, 8), TagSample { records: 4, found: true, other_type: 0 });
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn untagged_files_and_text_tags_are_not_haplotagged() {
    let dir = temp_dir("untagged");
    let mut records = reads("m", 0, 0, 4, Hp::None);
    records.extend(reads("a", 1, MID, 6, Hp::None));
    write(&dir, &records).expect("inputs");
    // chr1's 6, then chrM's 4.
    assert_eq!(sample(&dir, 100, 8), TagSample { records: 10, found: false, other_type: 0 });

    // HP:Z, which the tools that read haplotypes take for no tag.
    let records = reads("z", 1, MID, 3, Hp::Text("1"));
    write(&dir, &records).expect("inputs");
    assert_eq!(sample(&dir, 100, 8), TagSample { records: 3, found: false, other_type: 3 });

    // Only short contigs with reads: those are sampled.
    let mut records = reads("m", 0, 0, 2, Hp::None);
    records.push(read("t", 0, 300, Hp::Int(2), false));
    write(&dir, &records).expect("inputs");
    assert_eq!(sample(&dir, 100, 8), TagSample { records: 3, found: true, other_type: 0 });
    let _ = fs::remove_dir_all(dir);
}

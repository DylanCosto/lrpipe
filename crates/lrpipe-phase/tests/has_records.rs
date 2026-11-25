//! Which contigs of an indexed CRAM have records is read from its index
//! ([`AlignmentFile::has_records`]), and agrees with a region read of each whole contig: a
//! contig with mapped reads, one without any, and one with only an unmapped read placed on it.

use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;

use lrpipe_phase::alignment::{AlignmentFile, Fetch};
use pasta::core::Position;
use pasta::sam::{
    self,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
        },
        record_buf::{Cigar, QualityScores, Sequence},
    },
    header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    },
};
use pasta::{cram, fasta};

const LEN: usize = 300;
const CONTIGS: [&str; 3] = ["chr1", "chr2", "chr3"];

fn reference(contig: usize) -> Vec<u8> {
    (0..LEN).map(|i| b"ACGT"[(i * 7 + i / 5 + contig) % 4]).collect()
}

fn read(name: &str, contig: usize, start: usize, flags: Flags) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(contig)
        .set_alignment_start(Position::new(start + 1).expect("non-zero"))
        .set_sequence(Sequence::from(reference(contig)[start..start + 100].to_vec()))
        .set_quality_scores(QualityScores::from(vec![30; 100]));
    if !flags.is_unmapped() {
        builder = builder
            .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
            .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect::<Cigar>());
    }
    builder.build()
}

fn write_cram(dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let mut fa = Vec::new();
    let mut fai = String::new();
    for (i, name) in CONTIGS.iter().enumerate() {
        fa.extend_from_slice(format!(">{name}\n").as_bytes());
        let offset = fa.len();
        fa.extend_from_slice(&reference(i));
        fa.push(b'\n');
        fai.push_str(&format!("{name}\t{LEN}\t{offset}\t{LEN}\t{}\n", LEN + 1));
    }
    fs::write(dir.join("ref.fa"), fa)?;
    fs::write(dir.join("ref.fa.fai"), fai)?;
    let mut header = sam::Header::builder()
        .set_header(Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?);
    for name in CONTIGS {
        header = header
            .add_reference_sequence(name, Map::<ReferenceSequence>::new(NonZeroUsize::new(LEN).expect("non-empty")));
    }
    let header = header.build();
    let fasta = fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("ref.fa"))?;
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta));
    let path = dir.join("reads.cram");
    let mut writer =
        cram::io::writer::Builder::default().set_reference_sequence_repository(repository).build_from_path(&path)?;
    writer.write_header(&header)?;
    for i in 0..5 {
        writer.write_alignment_record(&header, &read(&format!("r{i}"), 0, i * 20, Flags::empty()))?;
    }
    writer.write_alignment_record(&header, &read("placed", 2, 50, Flags::UNMAPPED))?;
    writer.try_finish(&header)?;
    drop(writer);
    cram::crai::fs::write(dir.join("reads.cram.crai"), &cram::fs::index(&path)?)?;
    Ok(())
}

#[test]
fn a_crams_index_says_which_contigs_a_region_read_finds_records_on() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-phase-has-records-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_cram(&dir)?;
    let path = dir.join("reads.cram");
    let open = || -> Result<AlignmentFile, Box<dyn std::error::Error>> {
        let mut file = AlignmentFile::open(path.to_str().ok_or("path")?, Some(&dir.join("ref.fa")), 1)?;
        assert!(file.load_index()?);
        Ok(file)
    };
    let mut file = open()?;
    let by_index: Vec<bool> = (0..CONTIGS.len()).map(|tid| file.has_records(tid)).collect::<Result<_, _>>()?;
    let mut file = open()?;
    let mut by_reading = Vec::new();
    for tid in 0..CONTIGS.len() {
        file.fetch(Fetch::Region(tid, 0, i64::MAX >> 1))?;
        by_reading.push(file.read()?.is_some());
    }
    assert_eq!(by_index, [true, false, true]);
    assert_eq!(by_index, by_reading);
    fs::remove_dir_all(&dir)?;
    Ok(())
}

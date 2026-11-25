//! Reading a CRAM for alleles ([`Decoding::Alleles`]) decodes the read group and the `PS`, `HP`
//! and `BX` tags that allele detection reads, and leaves the other tags out; everything else a
//! record holds is what a full decode gives.

use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;

use lrpipe_phase::alignment::{AlignmentFile, Decoding, Fetch, Record};
use pasta::core::Position;
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
        record_buf::{
            Cigar, QualityScores, Sequence,
            data::field::{Value, value::Array},
        },
    },
    header::record::value::{
        Map,
        map::{self, ReadGroup, ReferenceSequence},
    },
};
use pasta::{cram, fasta};

const LEN: usize = 300;

fn reference() -> Vec<u8> {
    (0..LEN).map(|i| b"ACGT"[(i * 7 + i / 5) % 4]).collect()
}

fn read(i: usize) -> RecordBuf {
    let start = i * 20;
    let mut seq = reference()[start..start + 100].to_vec();
    seq[50] = if seq[50] == b'A' { b'C' } else { b'A' };
    RecordBuf::builder()
        .set_name(format!("r{i}"))
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(start + 1).expect("non-zero"))
        .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
        .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from((0..100).map(|q| q as u8 % 40).collect::<Vec<_>>()))
        .set_data(
            [
                (Tag::READ_GROUP, Value::String("rg1".into())),
                (Tag::new(b'H', b'P'), Value::from(1 + (i % 2) as u8)),
                (Tag::new(b'P', b'S'), Value::from(1000 + i as i32)),
                (Tag::new(b'B', b'X'), Value::String(format!("BARCODE{i}").into())),
                (Tag::BASE_MODIFICATIONS, Value::String("C+m?,0;".into())),
                (Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(vec![200]))),
                (Tag::new(b'X', b'X'), Value::String("unread".into())),
            ]
            .into_iter()
            .collect(),
        )
        .build()
}

fn write_cram(dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(dir.join("ref.fa"), [&b">chr1\n"[..], &reference(), b"\n"].concat())?;
    fs::write(dir.join("ref.fa.fai"), format!("chr1\t{LEN}\t6\t{LEN}\t{}\n", LEN + 1))?;
    let header = sam::Header::builder()
        .set_header(Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(NonZeroUsize::new(LEN).expect("non-empty")))
        .add_read_group("rg1", Map::<ReadGroup>::builder().insert(map::read_group::tag::SAMPLE, "S").build()?)
        .build();
    let fasta = fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("ref.fa"))?;
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta));
    let path = dir.join("reads.cram");
    let mut writer =
        cram::io::writer::Builder::default().set_reference_sequence_repository(repository).build_from_path(&path)?;
    writer.write_header(&header)?;
    for i in 0..10 {
        writer.write_alignment_record(&header, &read(i))?;
    }
    writer.try_finish(&header)?;
    drop(writer);
    cram::crai::fs::write(dir.join("reads.crai"), &cram::fs::index(&path)?)?;
    Ok(())
}

fn records(dir: &Path, decoding: Decoding) -> Result<Vec<Record>, Box<dyn std::error::Error>> {
    let mut file = AlignmentFile::open(dir.join("reads.cram").to_str().ok_or("path")?, Some(&dir.join("ref.fa")), 1)?;
    file.set_decoding(decoding);
    file.fetch(Fetch::All)?;
    let mut out = Vec::new();
    while let Some(record) = file.read()? {
        out.push(record);
    }
    Ok(out)
}

#[test]
fn alleles_decoding_keeps_what_allele_detection_reads() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-phase-cram-decoding-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_cram(&dir)?;
    let (full, alleles) = (records(&dir, Decoding::Full)?, records(&dir, Decoding::Alleles)?);
    assert_eq!(full.len(), 10);
    assert_eq!(alleles.len(), full.len());
    for (f, a) in full.iter().zip(&alleles) {
        let name = String::from_utf8_lossy(f.qname()).into_owned();
        assert_eq!((a.qname(), a.pos(), a.mapq()), (f.qname(), f.pos(), f.mapq()), "{name}");
        assert_eq!(a.cigar().iter().collect::<Vec<_>>(), f.cigar().iter().collect::<Vec<_>>(), "{name}");
        assert_eq!((a.seq(), a.qual()), (f.seq(), f.qual()), "{name}");
        for tag in [b"RG", b"HP", b"PS", b"BX"] {
            let value = f.aux(tag).map_err(|e| format!("{name}: full decode: {e}"))?;
            assert_eq!(a.aux(tag).ok(), Some(value), "{name}: {}", String::from_utf8_lossy(tag));
        }
        for tag in [b"MM", b"ML", b"XX", b"MD", b"NM"] {
            assert!(f.aux(tag).is_ok(), "{name}: full decode lacks {}", String::from_utf8_lossy(tag));
            assert!(a.aux(tag).is_err(), "{name}: alleles decode has {}", String::from_utf8_lossy(tag));
        }
    }
    fs::remove_dir_all(&dir)?;
    Ok(())
}

/// A CRAM of many small slices: 600 reads of 100 bases across a 20 kb contig, some reverse and
/// some unmapped but placed, 16 records per slice.
fn write_sliced_cram(dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    use pasta::cram::file_definition::Version;
    use pasta::cram::io::writer::compression::{Options, Profile};
    const CONTIG: usize = 20_000;
    let reference: Vec<u8> = (0..CONTIG).map(|i| b"ACGT"[(i * 7 + i / 5 + i / 97) % 4]).collect();
    fs::write(dir.join("long.fa"), [&b">chr1\n"[..], &reference, b"\n"].concat())?;
    fs::write(dir.join("long.fa.fai"), format!("chr1\t{CONTIG}\t6\t{CONTIG}\t{}\n", CONTIG + 1))?;
    let header = sam::Header::builder()
        .set_header(Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(NonZeroUsize::new(CONTIG).expect("non-empty")))
        .add_read_group("rg1", Map::<ReadGroup>::builder().insert(map::read_group::tag::SAMPLE, "S").build()?)
        .build();
    let fasta = fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("long.fa"))?;
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta));
    let mut options = Options::new(Version::new(3, 1), Profile::Normal);
    options.records_per_slice = 16;
    let path = dir.join("long.cram");
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .set_compression_options(options)
        .build_from_path(&path)?;
    writer.write_header(&header)?;
    for i in 0..600 {
        let start = i * 33;
        let mut seq = reference[start..start + 100].to_vec();
        seq[i % 100] = if seq[i % 100] == b'A' { b'C' } else { b'A' };
        let flags = match i % 7 {
            3 => Flags::REVERSE_COMPLEMENTED,
            5 => Flags::UNMAPPED,
            _ => Flags::empty(),
        };
        let mut builder = RecordBuf::builder()
            .set_name(format!("r{i}"))
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(start + 1).expect("non-zero"))
            .set_sequence(Sequence::from(seq))
            .set_quality_scores(QualityScores::from(vec![(i % 40) as u8; 100]))
            .set_data(
                [
                    (Tag::READ_GROUP, Value::String("rg1".into())),
                    (Tag::new(b'H', b'P'), Value::from(1 + (i % 2) as u8)),
                    (Tag::BASE_MODIFICATIONS, Value::String("C+m?;".into())),
                ]
                .into_iter()
                .collect(),
            );
        if !flags.is_unmapped() {
            builder = builder
                .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
                .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect::<Cigar>());
        }
        writer.write_alignment_record(&header, &builder.build())?;
    }
    writer.try_finish(&header)?;
    drop(writer);
    cram::crai::fs::write(dir.join("long.cram.crai"), &cram::fs::index(&path)?)?;
    Ok(())
}

/// The records of `[beg, end)` of chr1, read for alleles on `threads` threads.
fn region(dir: &Path, threads: usize, beg: i64, end: i64) -> Result<Vec<Record>, Box<dyn std::error::Error>> {
    let mut file =
        AlignmentFile::open(dir.join("long.cram").to_str().ok_or("path")?, Some(&dir.join("long.fa")), threads)?;
    file.set_decoding(Decoding::Alleles);
    assert!(file.load_index()?, "the index is found");
    file.fetch(Fetch::Region(0, beg, end))?;
    let mut out = Vec::new();
    while let Some(record) = file.read()? {
        out.push(record);
    }
    Ok(out)
}

#[test]
fn threaded_region_reads_give_the_serial_records() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-phase-cram-threaded-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_sliced_cram(&dir)?;
    // The whole contig (and past its end, as pysam's open intervals), a stretch inside it,
    // single bases, and nothing.
    let regions = [(0, (1 << 31) - 1), (0, 20_000), (5_000, 9_000), (6_532, 6_533), (19_000, 30_000), (300, 300)];
    for (beg, end) in regions {
        let serial = region(&dir, 1, beg, end)?;
        for threads in [2, 4] {
            assert_eq!(region(&dir, threads, beg, end)?, serial, "[{beg}, {end}) on {threads} threads");
        }
        if end > beg + 1_000 {
            assert!(serial.len() > 20, "[{beg}, {end}): {} records", serial.len());
        }
    }
    // A read stopped early, then another fetch on the same reader.
    let mut file = AlignmentFile::open(dir.join("long.cram").to_str().ok_or("path")?, Some(&dir.join("long.fa")), 4)?;
    file.set_decoding(Decoding::Alleles);
    file.load_index()?;
    file.fetch(Fetch::Region(0, 0, 20_000))?;
    file.read()?.ok_or("a record")?;
    file.fetch(Fetch::Region(0, 5_000, 9_000))?;
    let mut again = Vec::new();
    while let Some(record) = file.read()? {
        again.push(record);
    }
    assert_eq!(again, region(&dir, 1, 5_000, 9_000)?);
    fs::remove_dir_all(&dir)?;
    Ok(())
}

//! Counting within another pass's pieces writes the table a run writes: reads at staggered
//! starts, some with 5hmC calls alone, CpGs across window ends, pieces merged out of order.

use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;

use pasta::bam;
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
            Cigar, Sequence,
            data::field::{Value, value::Array},
        },
    },
    header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    },
};

use lrpipe_methyl::input::ReadFilter;
use lrpipe_methyl::mods::Code;
use lrpipe_methyl::pileup::{Calling, Mods, Motif, Partition, Rows, Strands, Thresholds};
use lrpipe_methyl::read::Convention;
use lrpipe_methyl::run::{Coordinates, Job, Report, run};
use lrpipe_methyl::stream::Pass;

const LEN: usize = 200;
const WINDOW: i64 = 50;

/// `AT` repeats with CpGs here and there, two of them across a window's end (49-50, 99-100).
fn reference() -> Vec<u8> {
    let mut seq: Vec<u8> = (0..LEN).map(|i| b"AT"[i % 2]).collect();
    for s in [3, 17, 30, 49, 64, 99, 120, 133, 160, 181, 197] {
        seq[s] = b'C';
        seq[s + 1] = b'G';
    }
    seq
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    }
}

/// Read `i`: 60 bases from `15 i`, reverse when `i` is odd, with 5hmC calls alone when `i` is a
/// multiple of 3 and 5mC and 5hmC calls otherwise, at every C as sequenced.
fn read(i: usize) -> RecordBuf {
    let start = i * 15;
    let seq = reference()[start..(start + 60).min(LEN)].to_vec();
    let reverse = i % 2 == 1;
    let sequenced: Vec<u8> = if reverse { seq.iter().rev().map(|&b| complement(b)).collect() } else { seq.clone() };
    let cs = sequenced.iter().filter(|&&b| b == b'C').count();
    let group = |code: char| format!("C+{code}?{};", ",0".repeat(cs));
    let value = |k: usize| u8::try_from((i * 37 + k * 53) % 256).expect("a byte");
    let (mm, ml): (String, Vec<u8>) = if i.is_multiple_of(3) {
        (group('h'), (0..cs).map(value).collect())
    } else {
        (group('m') + &group('h'), (0..2 * cs).map(value).collect())
    };
    RecordBuf::builder()
        .set_name(format!("r{i}"))
        .set_flags(if reverse { Flags::REVERSE_COMPLEMENTED } else { Flags::empty() })
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(start + 1).expect("non-zero"))
        .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
        .set_cigar([Op::new(Kind::Match, seq.len())].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(seq))
        .set_data(
            [
                (Tag::BASE_MODIFICATIONS, Value::String(mm.into())),
                (Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml))),
            ]
            .into_iter()
            .collect(),
        )
        .build()
}

const READS: usize = 13;

fn write_inputs(dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let reference = reference();
    fs::write(dir.join("ref.fa"), [&b">chr1\n"[..], &reference, b"\n"].concat())?;
    fs::write(dir.join("ref.fa.fai"), format!("chr1\t{LEN}\t6\t{LEN}\t{}\n", LEN + 1))?;
    let length = NonZeroUsize::new(LEN).expect("non-empty");
    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?;
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(length))
        .build();
    let path = dir.join("reads.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&path)?);
    writer.write_header(&header)?;
    for i in 0..READS {
        writer.write_alignment_record(&header, &read(i))?;
    }
    writer.try_finish()?;
    drop(writer);
    match bam::fs::index(&path)? {
        bam::Index::Bai(index) => bam::bai::fs::write(dir.join("reads.bam.bai"), &index)?,
        bam::Index::Csi(index) => pasta::csi::fs::write(dir.join("reads.bam.csi"), &index)?,
    }
    Ok(())
}

/// The table counted in pieces from these starts, merged last piece first.
fn in_pieces(dir: &Path, job: &Job, starts: &[i64]) -> Result<String, Box<dyn std::error::Error>> {
    let reference = reference();
    let fetch = |s: i64, e: i64| -> std::io::Result<Vec<u8>> {
        let (s, e) = (s.clamp(0, LEN as i64) as usize, e.clamp(0, LEN as i64) as usize);
        Ok(reference[s..e.max(s)].to_vec())
    };
    let pass = Pass::new(job)?;
    let contig = pass.contig("chr1", LEN as i64, fetch)?;
    let mut pieces: Vec<_> = starts.iter().map(|&s| contig.piece(s)).collect();
    let mut reader = bam::io::Reader::new(fs::File::open(dir.join("reads.bam"))?);
    reader.read_header()?;
    for record in reader.records() {
        let record = record?;
        let start = record.alignment_start().transpose()?.map_or(0, |p| usize::from(p) as i64 - 1);
        let k = starts.iter().rposition(|&s| s <= start).expect("a piece");
        pieces[k].add(&record)?;
    }
    for (k, piece) in pieces.into_iter().enumerate().rev() {
        contig.merge(k, piece);
    }
    let counted = contig.finish();
    let mut text = String::new();
    for w in 0..counted.windows() {
        counted.write_window(w, fetch, &mut text)?;
    }
    Ok(text)
}

#[test]
fn pieces_write_the_table_a_run_writes() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-stream-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let mut thresholds = Thresholds::NONE;
    thresholds.set(b'C', 0.6);
    for strands in [Strands::Combined, Strands::Separate] {
        let rows = Rows {
            modified_bases: Some(vec![(b'C', Some(Code::Letter(b'm'))), (b'C', Some(Code::Letter(b'h')))]),
            motif: Motif::CpG,
            strands,
            mods: Mods::Separate,
        };
        let job = Job {
            reads: dir.join("reads.bam"),
            reference: Some(dir.join("ref.fa")),
            regions: Vec::new(),
            coordinates: Coordinates::ZeroBased,
            window: WINDOW,
            threads: 2,
            filter: ReadFilter::PRIMARY,
            report: Report::BedMethyl { calling: Calling { thresholds, convention: Convention::Reference }, rows },
            outputs: vec![(Partition::All, dir.join("run.bed"))],
            haplotags: None,
        };
        run(&job)?;
        let whole = fs::read_to_string(dir.join("run.bed"))?;
        assert!(whole.lines().count() > 20, "{whole}");
        for starts in [&[0][..], &[0, 37, 90, 140], &[0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]] {
            assert_eq!(in_pieces(&dir, &job, starts)?, whole, "pieces from {starts:?}, {strands:?}");
        }
    }
    fs::remove_dir_all(&dir)?;
    Ok(())
}

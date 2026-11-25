//! A whole run on a small reference and BAM with known repeat alleles.

use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;

use lrpipe_tr::catalog;
use lrpipe_tr::extract::RefWindow;
use lrpipe_tr::reads::Reads;
use lrpipe_tr::run::{Job, RunOptions, records, run};
use lrpipe_tr::stream::{Collector, Observations};
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
        record_buf::{Cigar, QualityScores, Sequence, data::field::Value},
    },
    header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    },
};

const REPEAT_START: usize = 1000;
const REPEAT_COPIES: usize = 10;

/// A reference without tandem repeats except `CAG` x 10 at 1000..1030.
fn reference() -> Vec<u8> {
    let mut x: u64 = 12345;
    let mut seq: Vec<u8> = (0..2000)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 33) as usize % 4]
        })
        .collect();
    seq[REPEAT_START..REPEAT_START + 3 * REPEAT_COPIES].copy_from_slice(&b"CAG".repeat(REPEAT_COPIES));
    seq
}

/// A read over 400..1700 carrying `copies` CAGs, the extra ones aligned as an insertion at
/// the repeat start.
fn read(name: &str, reference: &[u8], copies: usize, haplotype: u8) -> RecordBuf {
    let end = REPEAT_START + 3 * REPEAT_COPIES;
    let seq = [&reference[400..REPEAT_START], &b"CAG".repeat(copies)[..], &reference[end..1700]].concat();
    let extra = 3 * (copies - REPEAT_COPIES);
    let cigar: Cigar = if extra == 0 {
        [Op::new(Kind::Match, 1300)].into_iter().collect()
    } else {
        [
            Op::new(Kind::Match, REPEAT_START - 400),
            Op::new(Kind::Insertion, extra),
            Op::new(Kind::Match, 1700 - REPEAT_START),
        ]
        .into_iter()
        .collect()
    };
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(401).expect("non-zero"))
        .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
        .set_cigar(cigar)
        // Qualities and a tag nothing here reads, which CRAM decoding leaves out.
        .set_quality_scores(QualityScores::from(vec![30; seq.len()]))
        .set_sequence(Sequence::from(seq))
        .set_data(
            [(Tag::new(b'H', b'P'), Value::from(haplotype)), (Tag::new(b'X', b'X'), Value::String("unread".into()))]
                .into_iter()
                .collect(),
        )
        .build()
}

fn write_inputs(dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let reference = reference();
    fs::write(dir.join("ref.fa"), [&b">chr1\n"[..], &reference, b"\n"].concat())?;
    fs::write(dir.join("ref.fa.fai"), format!("chr1\t{0}\t6\t{0}\t{1}\n", reference.len(), reference.len() + 1))?;
    fs::write(dir.join("catalog.bed"), "chr1\t1000\t1030\tID=TEST;MOTIFS=CAG;STRUC=(CAG)n\n")?;

    let length = NonZeroUsize::new(reference.len()).expect("non-empty reference");
    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?;
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(length))
        .build();
    let records: Vec<RecordBuf> = (0..6)
        .flat_map(|i| [read(&format!("ref{i}"), &reference, 10, 1), read(&format!("alt{i}"), &reference, 14, 2)])
        .collect();
    let path = dir.join("reads.bam");
    let mut writer = bam::io::Writer::new(fs::File::create(&path)?);
    writer.write_header(&header)?;
    for record in &records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish()?;
    drop(writer);
    match bam::fs::index(&path)? {
        bam::Index::Bai(index) => bam::bai::fs::write(dir.join("reads.bam.bai"), &index)?,
        bam::Index::Csi(index) => pasta::csi::fs::write(dir.join("reads.bam.csi"), &index)?,
    }
    // The same records as CRAM.
    let fasta = fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("ref.fa"))?;
    let repository = fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(fasta));
    let path = dir.join("reads.cram");
    let mut writer =
        cram::io::writer::Builder::default().set_reference_sequence_repository(repository).build_from_path(&path)?;
    writer.write_header(&header)?;
    for record in &records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish(&header)?;
    drop(writer);
    cram::crai::fs::write(dir.join("reads.cram.crai"), &cram::fs::index(&path)?)?;
    Ok(())
}

#[test]
fn calls_a_heterozygous_expansion() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-tr-e2e-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let job = Job {
        reads: dir.join("reads.bam"),
        genome: dir.join("ref.fa"),
        output: dir.join("out.vcf"),
        sample: "S".into(),
        haplotags: None,
    };
    let loci = catalog::read(&dir.join("catalog.bed"))?;
    let summary = run(&job, &loci, &RunOptions { threads: 2, ..RunOptions::default() })?;
    assert_eq!((summary.loci, summary.called), (1, 1));

    let vcf = fs::read_to_string(&job.output)?;
    let record = vcf.lines().find(|l| !l.starts_with('#')).ok_or("no record")?;
    let f: Vec<&str> = record.split('\t').collect();
    assert_eq!(f[1], "1000", "POS is the base before the repeat");
    assert_eq!(f[3].len(), 31, "REF is the anchor base and the repeat");
    assert_eq!(f[4].len(), 43, "ALT carries 14 CAGs");
    // Haplotagged reads: one allele per haplotype.
    assert_eq!(f[9], "0/1:30,42:30-30,42-42:6,6:10,14", "record {record}");
    fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn streamed_reads_give_the_same_records() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-tr-stream-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let (bam, genome) = (dir.join("reads.bam"), dir.join("ref.fa"));
    let loci = catalog::read(&dir.join("catalog.bed"))?;
    let opts = RunOptions { threads: 1, ..RunOptions::default() };
    let (want, _) = records(&bam, &genome, &loci, &opts)?;

    // Another program's pass: all reads of the contig, offered in two pieces.
    let mut reads = Vec::new();
    Reads::open(&bam, None)?.fetch("chr1", 0, 2000, opts.min_mapq, &mut reads)?;
    let reference = RefWindow::contig(&genome, "chr1")?;
    let (first, second) = reads.split_at(reads.len() / 2);
    let pieces: Vec<_> = [first, second]
        .into_iter()
        .map(|piece| {
            let mut c = Collector::new(&loci, &reference, &opts);
            for r in piece {
                if c.spans_any(r.start, r.end, lrpipe_tr::extract::soft_clips(&r.cigar)) {
                    c.push(r);
                }
            }
            c.into_observations()
        })
        .collect();
    let obs = Observations::merge(pieces, loci.len(), opts.max_reads);
    let (got, called) = obs.records(&loci, 0..loci.len(), &reference, &opts);
    assert_eq!(called, 1);
    assert_eq!(got, want);
    fs::remove_dir_all(&dir)?;
    Ok(())
}

/// Reading a CRAM decodes only what genotyping reads; the records must be those the same reads
/// give as BAM.
#[test]
fn cram_gives_the_records_bam_gives() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-tr-cram-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let genome = dir.join("ref.fa");
    let loci = catalog::read(&dir.join("catalog.bed"))?;
    let opts = RunOptions { threads: 1, ..RunOptions::default() };
    let (bam, _) = records(&dir.join("reads.bam"), &genome, &loci, &opts)?;
    let (cram, called) = records(&dir.join("reads.cram"), &genome, &loci, &opts)?;
    assert_eq!(called, 1);
    assert_eq!(cram, bam);
    // The reads themselves, haplotypes included (these alleles genotype alike without them).
    let fetch = |path: &Path| -> Result<_, Box<dyn std::error::Error>> {
        let mut reads = Vec::new();
        let reference = lrpipe_tr::reads::SharedReference::open(&genome)?;
        Reads::open(path, Some(&reference))?.fetch("chr1", 0, 2000, opts.min_mapq, &mut reads)?;
        Ok(reads)
    };
    let (from_bam, from_cram) = (fetch(&dir.join("reads.bam"))?, fetch(&dir.join("reads.cram"))?);
    assert!(from_bam.iter().all(|r| r.haplotype.is_some()), "{from_bam:?}");
    assert_eq!(from_cram, from_bam);
    fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn a_haplotag_list_replaces_the_hp_tags() -> Result<(), Box<dyn std::error::Error>> {
    // F1b: `lrpipe run` genotypes the repeats with its own haplotags, a list that wins over the
    // reads' HP tags (ref reads HP 1, alt reads HP 2): ref0 on haplotype 2, alt0 on 1, alt1 on
    // none, the rest unlisted.
    let dir = std::env::temp_dir().join(format!("lrpipe-tr-haplotags-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let list = dir.join("list.tsv");
    fs::write(&list, "ref0\tH2\t7\tchr1\nalt0\tH1\t7\tchr1\nalt1\tnone\tnone\tchr1\n")?;
    let haplotags = std::sync::Arc::new(lrpipe_reads::haplotags::Haplotags::read(&list)?);
    let mut reads = Reads::open(&dir.join("reads.bam"), None)?;
    reads.set_haplotags(std::sync::Arc::clone(&haplotags));
    let mut out = Vec::new();
    reads.fetch("chr1", 1000, 1030, 1, &mut out)?;
    let haplotypes: Vec<Option<u8>> = out.iter().map(|r| r.haplotype).collect();
    assert_eq!(&haplotypes[..4], [Some(2), Some(1), None, None]);
    assert!(haplotypes[4..].iter().all(Option::is_none), "{haplotypes:?}");
    assert_eq!(haplotags.matches().counts(), (12, 3));

    // The list of the reads' own tags gives the calls the tags give.
    let own: String = (0..6).map(|i| format!("ref{i}\tH1\t7\tchr1\nalt{i}\tH2\t7\tchr1\n")).collect();
    fs::write(&list, own)?;
    let opts = RunOptions { threads: 1, ..RunOptions::default() };
    let loci = catalog::read(&dir.join("catalog.bed"))?;
    let genome = dir.join("ref.fa");
    let listed = std::sync::Arc::new(lrpipe_reads::haplotags::Haplotags::read(&list)?);
    let (with_list, _) = lrpipe_tr::run::records_with(&dir.join("reads.bam"), &genome, &loci, &opts, Some(&listed))?;
    assert_eq!(with_list, records(&dir.join("reads.bam"), &genome, &loci, &opts)?.0);
    fs::remove_dir_all(&dir)?;
    Ok(())
}

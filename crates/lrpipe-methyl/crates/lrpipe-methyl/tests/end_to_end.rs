//! Whole runs on a small reference and BAM: forward reads methylated at every CpG, reverse
//! reads unmethylated, each strand on its own haplotype.

use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;

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
        record_buf::{
            Cigar, QualityScores, Sequence,
            data::field::{Value, value::Array},
        },
    },
    header::record::value::{
        Map,
        map::{self, ReferenceSequence},
    },
};

use lrpipe_methyl::cpg::Sites;
use lrpipe_methyl::input::ReadFilter;
use lrpipe_methyl::mods::Code;
use lrpipe_methyl::pileup::{Calling, Mods, Motif, Partition, Rows, Strands, Thresholds};
use lrpipe_methyl::read::Convention;
use lrpipe_methyl::run::{Coordinates, Estimate, Job, Report, run};

/// CpG sites of the reference (the C), one near the start where no read counts for
/// pb-CpG-tools' end rule.
const SITES: [usize; 4] = [10, 60, 120, 170];
const LEN: usize = 200;

fn reference() -> Vec<u8> {
    let mut seq: Vec<u8> = (0..LEN).map(|i| b"AT"[i % 2]).collect();
    for s in SITES {
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

/// An `MM` group listing the read-as-sequenced C's at `wanted` (indices into it), `?` mode.
fn mm(sequenced: &[u8], wanted: &[usize]) -> String {
    let mut out = String::from("C+m?");
    let mut skipped = 0;
    for (i, &b) in sequenced.iter().enumerate() {
        if b != b'C' {
            continue;
        }
        if wanted.contains(&i) {
            out.push_str(&format!(",{skipped}"));
            skipped = 0;
        } else {
            skipped += 1;
        }
    }
    out.push(';');
    out
}

/// A read over the whole reference: forward with 5mC value `value` at each CpG's C, or
/// reverse with it at each CpG's G (its own C).
fn read(name: &str, reverse: bool, value: u8, haplotype: u8) -> RecordBuf {
    let seq = reference();
    let (sequenced, wanted): (Vec<u8>, Vec<usize>) = if reverse {
        let rc: Vec<u8> = seq.iter().rev().map(|&b| complement(b)).collect();
        (rc, SITES.iter().map(|&s| LEN - 1 - (s + 1)).collect())
    } else {
        (seq.clone(), SITES.to_vec())
    };
    let mut wanted = wanted;
    wanted.sort_unstable();
    let ml = vec![value; wanted.len()];
    RecordBuf::builder()
        .set_name(name)
        .set_flags(if reverse { Flags::REVERSE_COMPLEMENTED } else { Flags::empty() })
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::new(1).expect("non-zero"))
        .set_mapping_quality(MappingQuality::new(60).expect("not 255"))
        .set_cigar([Op::new(Kind::Match, LEN)].into_iter().collect::<Cigar>())
        .set_sequence(Sequence::from(seq))
        // Qualities and a tag nothing here reads, which CRAM decoding leaves out.
        .set_quality_scores(QualityScores::from(vec![30; LEN]))
        .set_data(
            [
                (Tag::BASE_MODIFICATIONS, Value::String(mm(&sequenced, &wanted).into())),
                (Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml))),
                (Tag::new(b'H', b'P'), Value::from(haplotype)),
                (Tag::new(b'X', b'X'), Value::String("unread".into())),
            ]
            .into_iter()
            .collect(),
        )
        .build()
}

fn write_inputs(dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let reference = reference();
    fs::write(dir.join("ref.fa"), [&b">chr1\n"[..], &reference, b"\n"].concat())?;
    fs::write(dir.join("ref.fa.fai"), format!("chr1\t{LEN}\t6\t{LEN}\t{}\n", LEN + 1))?;
    let length = NonZeroUsize::new(LEN).expect("non-empty");
    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?;
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(length))
        // A contig without reads, which the reference lacks (as a chr20 reference serves reads
        // of chr20 with a whole genome's header).
        .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(length))
        .build();
    let records: Vec<RecordBuf> =
        (0..4).flat_map(|i| [read(&format!("f{i}"), false, 230, 1), read(&format!("r{i}"), true, 20, 2)]).collect();
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

fn job(dir: &Path, report: Report, outputs: Vec<(Partition, &str)>) -> Job {
    Job {
        reads: dir.join("reads.bam"),
        reference: Some(dir.join("ref.fa")),
        regions: Vec::new(),
        coordinates: Coordinates::OneBased,
        window: 50,
        threads: 2,
        filter: ReadFilter::PRIMARY,
        report,
        outputs: outputs.into_iter().map(|(p, name)| (p, dir.join(name))).collect(),
        haplotags: None,
    }
}

fn rows(path: &Path) -> Result<Vec<Vec<String>>, Box<dyn std::error::Error>> {
    Ok(fs::read_to_string(path)?
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect())
}

#[test]
fn bedmethyl_combines_a_cpgs_strands() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-bedmethyl-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let rows_wanted = Rows {
        modified_bases: Some(vec![(b'C', Some(Code::Letter(b'm')))]),
        motif: Motif::CpG,
        strands: Strands::Combined,
        mods: Mods::Separate,
    };
    let calling = Calling { thresholds: Thresholds::NONE, convention: Convention::Reference };
    run(&job(&dir, Report::BedMethyl { calling, rows: rows_wanted }, vec![(Partition::All, "out.bed")]))?;
    let got = rows(&dir.join("out.bed"))?;
    let positions: Vec<&str> = got.iter().map(|r| r[1].as_str()).collect();
    assert_eq!(positions, ["10", "60", "120", "170"], "one row per CpG");
    for r in &got {
        // Valid coverage 8 (four reads per strand), four modified, four canonical, strand '.'.
        assert_eq!((r[3].as_str(), r[5].as_str(), r[9].as_str()), ("m", ".", "8"), "row {r:?}");
        assert_eq!((r[10].as_str(), r[11].as_str(), r[12].as_str()), ("50.00", "4", "4"), "row {r:?}");
    }
    fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn cpg_counts_per_haplotype_and_leaves_out_read_ends() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-cpg-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let report = Report::CpgCount { sites: Sites::Reference, estimate: Estimate::Count, min_coverage: 4, min_mapq: 1 };
    let outputs =
        vec![(Partition::All, "c.bed"), (Partition::Haplotype(1), "h1.bed"), (Partition::Haplotype(2), "h2.bed")];
    run(&job(&dir, report, outputs))?;
    let combined = rows(&dir.join("c.bed"))?;
    // The CpG at 10 is within 20 aligned bases of every read's start.
    let positions: Vec<&str> = combined.iter().map(|r| r[1].as_str()).collect();
    assert_eq!(positions, ["60", "120", "170"]);
    for r in &combined {
        assert_eq!(&r[3..], ["50.0", "Total", "8", "4", "4", "0.900", "0.080"], "row {r:?}");
    }
    let hap1 = rows(&dir.join("h1.bed"))?;
    assert!(hap1.iter().all(|r| r[3..6] == ["100.0", "hap1", "4"]), "{hap1:?}");
    let hap2 = rows(&dir.join("h2.bed"))?;
    assert!(hap2.iter().all(|r| r[3..6] == ["0.0", "hap2", "4"]), "{hap2:?}");
    fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn a_haplotag_list_replaces_the_hp_tags() -> Result<(), Box<dyn std::error::Error>> {
    // F1b: `lrpipe run` haplotags the reads from its own phasing and gives the CpG tables the
    // list, which wins over the reads' HP tags (forward reads HP 1, reverse reads HP 2): f0 and
    // f1 on haplotype 2, r0 on haplotype 1, the rest on none (f2 is listed on another contig).
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-haplotags-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let list = dir.join("haplotags.tsv");
    fs::write(
        &list,
        "#readname\thaplotype\tphaseset\tchromosome\nf0\tH2\t5\tchr1\nf1\tH2\t5\tchr1\nr0\tH1\t5\tchr1\n\
         r1\tnone\tnone\tchr1\nf2\tH1\t5\tchr2\n",
    )?;
    let report = Report::CpgCount { sites: Sites::Reference, estimate: Estimate::Count, min_coverage: 1, min_mapq: 1 };
    let outputs =
        vec![(Partition::All, "c.bed"), (Partition::Haplotype(1), "h1.bed"), (Partition::Haplotype(2), "h2.bed")];
    let haplotags = lrpipe_methyl::input::Haplotags::read(&list)?;
    run(&Job { haplotags: Some(std::sync::Arc::new(haplotags)), ..job(&dir, report, outputs) })?;
    let combined = rows(&dir.join("c.bed"))?;
    assert!(combined.iter().all(|r| r[3..6] == ["50.0", "Total", "8"]), "{combined:?}");
    let hap1 = rows(&dir.join("h1.bed"))?;
    assert_eq!(hap1.len(), 3);
    assert!(hap1.iter().all(|r| r[3..6] == ["0.0", "hap1", "1"]), "{hap1:?}");
    let hap2 = rows(&dir.join("h2.bed"))?;
    assert!(hap2.iter().all(|r| r[3..6] == ["100.0", "hap2", "2"]), "{hap2:?}");

    // On the command line.
    let args = || -> Vec<std::ffi::OsString> {
        let mut args: Vec<std::ffi::OsString> =
            ["lrpipe-methyl", "cpg", "--min-coverage", "1", "--bam"].map(Into::into).to_vec();
        args.extend([dir.join("reads.bam").into(), "--ref".into(), dir.join("ref.fa").into()]);
        args.extend(["--output-prefix".into(), dir.join("out").into(), "--haplotags".into(), list.clone().into()]);
        args
    };
    lrpipe_methyl::cli::main(args())?;
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut pasta::bgzf::io::Reader::new(fs::File::open(dir.join("out.hap2.bed.gz"))?),
        &mut text,
    )?;
    let hap2: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(hap2.len() == 3 && hap2.iter().all(|l| l.contains("\t100.0\thap2\t2\t")), "{text}");
    fs::remove_dir_all(&dir)?;
    Ok(())
}

/// `records`, as `name` in `dir`, indexed, with [`write_inputs`]'s header.
fn write_bam(dir: &Path, name: &str, records: &[RecordBuf]) -> Result<(), Box<dyn std::error::Error>> {
    let length = NonZeroUsize::new(LEN).expect("non-empty");
    let hd = Map::<map::Header>::builder().insert(map::header::tag::SORT_ORDER, "coordinate").build()?;
    let header = sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(length))
        .build();
    let path = dir.join(name);
    let mut writer = bam::io::Writer::new(fs::File::create(&path)?);
    writer.write_header(&header)?;
    for record in records {
        writer.write_alignment_record(&header, record)?;
    }
    writer.try_finish()?;
    drop(writer);
    let bam::Index::Bai(index) = bam::fs::index(&path)? else { panic!("a BAI") };
    bam::bai::fs::write(path.with_added_extension("bai"), &index)?;
    Ok(())
}

/// `record` with `flags` and the haplotype `hp` (none: no HP tag).
fn with(mut record: RecordBuf, flags: Flags, hp: Option<u8>) -> RecordBuf {
    *record.flags_mut() = flags;
    let data = record.data_mut();
    data.remove(&Tag::new(b'H', b'P'));
    if let Some(hp) = hp {
        data.insert(Tag::new(b'H', b'P'), Value::from(hp));
    }
    record
}

#[test]
fn a_list_tags_primary_alignments_alone_as_haplotag_does() -> Result<(), Box<dyn std::error::Error>> {
    // F1b: a read's supplementary and secondary alignments get no haplotype from the list, as
    // `haplotag` (`--tag-supplementary skip`) leaves them untagged: the tables of reads tagged
    // elsewhere, given the list, are those of the reads `haplotag` writes. pb-CpG-tools counts the
    // supplementary alignment among all reads, not the secondary one.
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-primary-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let f0 = read("f0", false, 230, 1);
    let r0 = read("r0", true, 20, 2);
    let supplementary = Flags::SUPPLEMENTARY;
    let secondary = Flags::SECONDARY;
    // As tagged by another pipeline, and as `haplotag` tags them from the list.
    let elsewhere =
        [f0.clone(), with(f0.clone(), supplementary, Some(1)), with(f0.clone(), secondary, Some(1)), r0.clone()];
    let haplotagged = [
        with(f0.clone(), Flags::empty(), Some(2)),
        with(f0.clone(), supplementary, None),
        with(f0.clone(), secondary, None),
        with(r0, Flags::REVERSE_COMPLEMENTED, None),
    ];
    write_bam(&dir, "elsewhere.bam", &elsewhere)?;
    write_bam(&dir, "haplotagged.bam", &haplotagged)?;
    let list = dir.join("list.tsv");
    fs::write(&list, "f0\tH2\t5\tchr1\nr0\tnone\tnone\tchr1\n")?;
    let report = Report::CpgCount { sites: Sites::Reference, estimate: Estimate::Count, min_coverage: 1, min_mapq: 1 };
    let tables =
        |reads: &str, prefix: &str, list: Option<&Path>| -> Result<Vec<Vec<Vec<String>>>, Box<dyn std::error::Error>> {
            let outputs = [(Partition::All, "c"), (Partition::Haplotype(1), "h1"), (Partition::Haplotype(2), "h2")]
                .map(|(p, t)| (p, format!("{prefix}.{t}.bed")));
            let mut job = job(&dir, report.clone(), outputs.iter().map(|(p, t)| (*p, t.as_str())).collect());
            job.reads = dir.join(reads);
            // As `cpg` reads them: supplementary alignments too.
            job.filter = ReadFilter { skip: Flags::UNMAPPED.union(Flags::SECONDARY), ..ReadFilter::PRIMARY };
            job.haplotags = list.map(lrpipe_methyl::input::Haplotags::read).transpose()?.map(std::sync::Arc::new);
            run(&job)?;
            outputs.iter().map(|(_, t)| rows(&dir.join(t))).collect()
        };
    let listed = tables("elsewhere.bam", "listed", Some(&list))?;
    assert_eq!(listed, tables("haplotagged.bam", "tagged", None)?);
    assert!(listed[0].iter().all(|r| r[5] == "3"), "f0, its supplementary alignment and r0: {:?}", listed[0]);
    assert!(listed[1].is_empty(), "{:?}", listed[1]);
    assert!(listed[2].len() == 3 && listed[2].iter().all(|r| r[5] == "1"), "f0's primary alignment: {:?}", listed[2]);
    fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn cpg_combined_only_writes_no_haplotype_tables() -> Result<(), Box<dyn std::error::Error>> {
    // F1a: `lrpipe run` on reads without HP tags asks for the combined table alone, rather than
    // empty hap1 and hap2 tables.
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-combined-only-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let prefix = dir.join("out");
    let run_cpg = |extra: &[&str]| {
        let mut args: Vec<std::ffi::OsString> = ["lrpipe-methyl", "cpg", "--bam"].map(Into::into).to_vec();
        args.extend([dir.join("reads.bam").into(), "--ref".into(), dir.join("ref.fa").into()]);
        args.extend(["--output-prefix".into(), prefix.clone().into_os_string()]);
        args.extend(extra.iter().map(Into::into));
        lrpipe_methyl::cli::main(args)
    };
    let table = |name: &str| dir.join(format!("out.{name}.bed.gz"));
    run_cpg(&["--combined-only"])?;
    assert!(table("combined").exists());
    assert!(!table("hap1").exists() && !table("hap2").exists());
    let combined = fs::read(table("combined"))?;
    run_cpg(&[])?;
    assert!(table("hap1").exists() && table("hap2").exists());
    assert_eq!(fs::read(table("combined"))?, combined, "the combined table is the same either way");
    fs::remove_dir_all(&dir)?;
    Ok(())
}

/// Reading a CRAM decodes only what methylation reads ([`lrpipe_methyl::input::TAGS`]); the
/// tables must be those the same reads give as BAM.
#[test]
fn cram_gives_the_tables_bam_gives() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-cram-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let bedmethyl = || Report::BedMethyl {
        calling: Calling { thresholds: Thresholds::NONE, convention: Convention::Reference },
        rows: Rows {
            modified_bases: Some(vec![(b'C', Some(Code::Letter(b'm')))]),
            motif: Motif::CpG,
            strands: Strands::Separate,
            mods: Mods::Separate,
        },
    };
    let cpg = || Report::CpgCount { sites: Sites::Reference, estimate: Estimate::Count, min_coverage: 1, min_mapq: 1 };
    for (name, report) in [("bedmethyl", bedmethyl()), ("cpg", cpg())] {
        let mut tables = Vec::new();
        for reads in ["reads.bam", "reads.cram"] {
            let outputs = vec![(Partition::All, "all.bed"), (Partition::Haplotype(1), "h1.bed")];
            let outputs = match report {
                Report::BedMethyl { .. } => outputs[..1].to_vec(),
                _ => outputs,
            };
            let mut job = job(&dir, report.clone(), outputs);
            job.reads = dir.join(reads);
            run(&job)?;
            let mut text = fs::read_to_string(dir.join("all.bed"))?;
            if dir.join("h1.bed").exists() {
                text += &fs::read_to_string(dir.join("h1.bed"))?;
            }
            tables.push(text);
        }
        assert!(tables[0].lines().count() > 2, "{name}: {}", tables[0]);
        assert_eq!(tables[0], tables[1], "{name}: the CRAM's table differs from the BAM's");
    }
    fs::remove_dir_all(&dir)?;
    Ok(())
}

/// LPH-004: a reference that lacks a contig with reads, or has it of another length, is an error for the
/// tables that need the reference, not tables without that contig's rows. (Every run here has
/// a contig without reads that the reference lacks, which is not an error.)
#[test]
fn a_reference_that_does_not_match_the_reads_is_an_error() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-mismatch-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let reference = reference();
    // Named as Ensembl names it, against the reads' UCSC name.
    fs::write(dir.join("ensembl.fa"), [&b">1\n"[..], &reference, b"\n"].concat())?;
    fs::write(dir.join("ensembl.fa.fai"), format!("1\t{LEN}\t3\t{LEN}\t{}\n", LEN + 1))?;
    let short = LEN - 50;
    fs::write(dir.join("short.fa"), [&b">chr1\n"[..], &reference[..short], b"\n"].concat())?;
    fs::write(dir.join("short.fa.fai"), format!("chr1\t{short}\t6\t{short}\t{}\n", short + 1))?;
    // Another assembly's chr1, longer than the reads' (as CHM13's chr2 is than GRCh38's).
    let long = [&reference[..], &reference[..50]].concat();
    fs::write(dir.join("long.fa"), [&b">chr1\n"[..], &long, b"\n"].concat())?;
    fs::write(dir.join("long.fa.fai"), format!("chr1\t{}\t6\t{}\t{}\n", LEN + 50, LEN + 50, LEN + 51))?;
    let cpg = || Report::CpgCount { sites: Sites::Reference, estimate: Estimate::Count, min_coverage: 1, min_mapq: 1 };
    let bedmethyl = |motif| Report::BedMethyl {
        calling: Calling { thresholds: Thresholds::NONE, convention: Convention::Read },
        rows: Rows { modified_bases: None, motif, strands: Strands::Separate, mods: Mods::Separate },
    };
    let cases = [
        ("ensembl.fa", "no contig chr1"),
        ("short.fa", "contig chr1 has 150 bases, but 200"),
        ("long.fa", "contig chr1 has 250 bases, but 200"),
    ];
    for (fa, wanted) in cases {
        for report in [cpg(), bedmethyl(Motif::CpG)] {
            let mut job = job(&dir, report, vec![(Partition::All, "out.bed")]);
            job.reference = Some(dir.join(fa));
            let e = run(&job).expect_err("mismatched reference").to_string();
            assert!(e.contains(wanted) && e.contains(fa), "{fa}: {e}");
        }
        // A table that does not read the reference does not need it to match.
        let mut job = job(&dir, bedmethyl(Motif::Any), vec![(Partition::All, "out.bed")]);
        job.reference = Some(dir.join(fa));
        run(&job)?;
        assert!(!rows(&dir.join("out.bed"))?.is_empty());
    }
    fs::remove_dir_all(&dir)?;
    Ok(())
}

/// `pileup`'s output `-` or `stdout` is the standard output, as modkit takes them; it wrote a
/// file of that name.
#[test]
fn pileup_writes_to_stdout_for_a_dash_or_stdout() -> Result<(), Box<dyn std::error::Error>> {
    let dir = std::env::temp_dir().join(format!("lrpipe-methyl-stdout-{}", std::process::id()));
    fs::create_dir_all(&dir)?;
    write_inputs(&dir)?;
    let pileup = |out: &str| {
        std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe-methyl"))
            .current_dir(&dir)
            .args(["pileup", "reads.bam", out, "--ref", "ref.fa", "--cpg", "--no-filtering"])
            .output()
    };
    let file = pileup("out.bed")?;
    assert!(file.status.success(), "{}", String::from_utf8_lossy(&file.stderr));
    let wanted = fs::read(dir.join("out.bed"))?;
    assert!(!wanted.is_empty());
    for out in ["-", "stdout"] {
        let piped = pileup(out)?;
        assert!(piped.status.success(), "{}", String::from_utf8_lossy(&piped.stderr));
        assert_eq!(piped.stdout, wanted, "{out}");
        assert!(!dir.join(out).exists(), "{out} is not a file");
    }
    fs::remove_dir_all(&dir)?;
    Ok(())
}

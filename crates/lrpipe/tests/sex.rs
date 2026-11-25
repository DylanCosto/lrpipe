//! Sex from hand-built BAM and CRAM headers and indexes.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;

use lrpipe::plan::Sex;
use lrpipe::sex::{check_reference, estimate, estimate_all, unplaced_unmapped};

const LEN: u32 = 1_000_000;

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-sex-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn contigs() -> Vec<String> {
    (1..=22).map(|i| format!("chr{i}")).chain(["chrX".into(), "chrY".into(), "chrM".into()]).collect()
}

/// Reads (or bytes) on each contig: autosomes 1000, chrX and chrY as given, chrM 5.
fn amounts(x: u64, y: u64) -> Vec<u64> {
    std::iter::repeat_n(1000, 22).chain([x, y, 5]).collect()
}

fn gzip(path: &Path, bytes: &[u8]) {
    let mut e = GzEncoder::new(fs::File::create(path).expect("created"), Compression::fast());
    e.write_all(bytes).expect("written");
    e.finish().expect("finished");
}

fn bam(dir: &Path, x: u64, y: u64) -> PathBuf {
    bam_of(dir, &amounts(x, y))
}

/// A BAM and its BAI with `amounts` of mapped reads on [`contigs`].
fn bam_of(dir: &Path, amounts: &[u64]) -> PathBuf {
    bam_on(dir, &contigs(), amounts)
}

/// A BAM and its BAI with `amounts` of mapped reads on `names`.
fn bam_on(dir: &Path, names: &[String], amounts: &[u64]) -> PathBuf {
    let path = dir.join("s.bam");
    let mut h = b"BAM\x01".to_vec();
    h.extend(0i32.to_le_bytes());
    h.extend((names.len() as i32).to_le_bytes());
    for name in names {
        h.extend((name.len() as i32 + 1).to_le_bytes());
        h.extend(name.as_bytes());
        h.push(0);
        h.extend(LEN.to_le_bytes());
    }
    gzip(&path, &h);

    let mut bai = b"BAI\x01".to_vec();
    bai.extend((names.len() as i32).to_le_bytes());
    for &mapped in amounts {
        bai.extend(2i32.to_le_bytes());
        // An ordinary bin with one chunk, then the pseudo-bin.
        bai.extend(4681u32.to_le_bytes());
        bai.extend(1i32.to_le_bytes());
        bai.extend([0u8; 16]);
        bai.extend(37_450u32.to_le_bytes());
        bai.extend(2i32.to_le_bytes());
        bai.extend([0u8; 16]);
        bai.extend(mapped.to_le_bytes());
        bai.extend(0u64.to_le_bytes());
        bai.extend(1i32.to_le_bytes());
        bai.extend(0u64.to_le_bytes());
    }
    fs::write(dir.join("s.bam.bai"), bai).expect("written");
    path
}

fn itf8(v: u32) -> Vec<u8> {
    match v {
        0..0x80 => vec![v as u8],
        0x80..0x4000 => vec![0x80 | (v >> 8) as u8, v as u8],
        _ => vec![0xc0 | (v >> 16) as u8, (v >> 8) as u8, v as u8],
    }
}

/// A CRAM container header of `records` records of `reference` (-1 for unmapped, -2 for several).
fn container(reference: i32, records: u32) -> Vec<u8> {
    let mut c = 0i32.to_le_bytes().to_vec();
    c.extend(match reference {
        -1 => vec![0xff, 0xff, 0xff, 0xff, 0x0f],
        -2 => vec![0xff, 0xff, 0xff, 0xff, 0x0e],
        r => itf8(r as u32),
    });
    c.extend(itf8(0)); // start
    c.extend(itf8(0)); // span
    c.extend(itf8(records));
    c
}

/// A CRAM of [`contigs`] and its index, with [`amounts`]' bytes of slices on them, and containers
/// of 7 and 300 unmapped reads (the second, as htslib writes them, of several references).
fn cram(dir: &Path, x: u64, y: u64) -> PathBuf {
    let path = dir.join("s.cram");
    let text: String = contigs().iter().map(|n| format!("@SQ\tSN:{n}\tLN:{LEN}\n")).collect();
    let mut data = (text.len() as i32).to_le_bytes().to_vec();
    data.extend(text.as_bytes());
    let mut block = vec![0u8, 0];
    block.extend(itf8(0));
    block.extend(itf8(data.len() as u32));
    block.extend(itf8(data.len() as u32));
    block.extend(&data);
    block.extend([0u8; 4]);

    let mut f = b"CRAM\x03\x00".to_vec();
    f.extend([0u8; 20]);
    f.extend((block.len() as i32).to_le_bytes());
    for _ in 0..4 {
        f.extend(itf8(0)); // reference, start, span, records
    }
    f.extend([0u8, 0]); // record counter, bases (LTF8)
    f.extend(itf8(1)); // blocks
    f.extend(itf8(0)); // landmarks
    f.extend([0u8; 4]); // CRC
    f.extend(block);
    let unmapped = [f.len(), f.len() + container(-1, 7).len()];
    f.extend(container(-1, 7));
    f.extend(container(-2, 300));
    fs::write(&path, f).expect("written");

    let crai: String = amounts(x, y)
        .iter()
        .enumerate()
        .flat_map(|(i, &bytes)| {
            [format!("{i}\t1\t500\t0\t0\t{}\n", bytes / 2), format!("{i}\t500\t500\t0\t0\t{}\n", bytes / 2)]
        })
        .chain(unmapped.iter().map(|o| format!("-1\t0\t0\t{o}\t0\t99999\n")))
        .collect();
    gzip(&dir.join("s.cram.crai"), crai.as_bytes());
    path
}

#[test]
fn bam_index() {
    let dir = temp_dir("bam");
    let e = estimate(&bam(&dir, 500, 300)).expect("estimated");
    assert_eq!((e.x, e.y, e.sex), (0.5, 0.3, Some(Sex::Male)));
    let e = estimate(&bam(&dir, 1000, 10)).expect("estimated");
    assert_eq!(e.sex, Some(Sex::Female));
    // F9: a BAI need not count the unplaced unmapped reads, and this one does not.
    assert_eq!(unplaced_unmapped(&dir.join("s.bam")).expect("read"), None);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn cram_index() {
    let dir = temp_dir("cram");
    let e = estimate(&cram(&dir, 500, 300)).expect("estimated");
    assert_eq!((e.x, e.y, e.sex), (0.5, 0.3, Some(Sex::Male)));
    let e = estimate(&cram(&dir, 980, 0)).expect("estimated");
    assert_eq!(e.sex, Some(Sex::Female));
    // F9: the unmapped reads are counted from the headers of their containers, which hold no
    // mapped reads (tests/qc.rs decodes one that does).
    assert_eq!(unplaced_unmapped(&dir.join("s.cram")).expect("read"), Some(307));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn reads_without_an_index_are_an_error() {
    let dir = temp_dir("none");
    let reads = bam(&dir, 500, 300);
    fs::remove_file(dir.join("s.bam.bai")).expect("removed");
    let err = estimate(&reads).expect_err("no index");
    assert!(err.to_string().contains("no .bai"), "{err}");
    let _ = fs::remove_dir_all(dir);
}

/// LUX-10: reads of one autosome (a chr20 BAM) are said to be so, not to have no reads on the
/// autosomes.
#[test]
fn reads_of_some_autosomes_are_said_to_be_too_few() {
    let dir = temp_dir("partial");
    let mut amounts = vec![0; 25];
    amounts[19] = 1000;
    let err = estimate(&bam_of(&dir, &amounts)).expect_err("too few autosomes").to_string();
    assert!(err.contains("has reads on 1 of the 22 autosomes (chr20)"), "{err}");
    assert!(err.contains("whole genome"), "{err}");
    // Aligned to a reference of chr20 alone (bonefish's BAM of unaligned reads): the same, not
    // "the index has no autosomes" (its chrX-less header was taken for one naming them `20`).
    let err = estimate(&bam_on(&dir, &["chr20".into()], &[1000])).expect_err("one autosome").to_string();
    assert!(err.contains("has reads on 1 of the 22 autosomes (chr20)"), "{err}");

    // Known from the reference before the reads are aligned.
    let reference = dir.join("ref.fa");
    let fai = |names: &[String]| {
        let text: String = names.iter().map(|n| format!("{n}\t{LEN}\t6\t60\t61\n")).collect();
        fs::write(dir.join("ref.fa.fai"), text).expect("written");
    };
    fai(&["chr20".into()]);
    let err = check_reference(&reference).expect_err("one autosome").to_string();
    assert!(err.contains("ref.fa has 1 of the 22 autosomes (chr20) and no chrX: the sex is estimated"), "{err}");
    fai(&contigs());
    check_reference(&reference).expect("a whole genome");
    let grch37: Vec<String> = (1..=22).map(|i| i.to_string()).chain(["X".into(), "MT".into()]).collect();
    fai(&grch37);
    check_reference(&reference).expect("GRCh37's names");
    let _ = fs::remove_dir_all(dir);
}

/// F5: a sample in several files is measured over all of them: their BAM indexes' reads added
/// up, which neither file alone tells; BAM's reads and CRAM's bytes each measured on their own
/// and averaged.
#[test]
fn several_files_add_up() {
    let (first, second) = (temp_dir("first"), temp_dir("second"));
    let (x_only, y_only) = (bam(&first, 0, 0), bam(&second, 1000, 600));
    assert_eq!(estimate(&x_only).expect("estimated").sex, None, "no X at all");
    assert_eq!(estimate(&y_only).expect("estimated").sex, None, "two Xs and a Y");
    let e = estimate_all(&[x_only.clone(), y_only]).expect("estimated");
    assert_eq!((e.x, e.y, e.sex), (0.5, 0.3, Some(Sex::Male)));

    let mixed = (bam(&first, 400, 200), cram(&second, 600, 400));
    let e = estimate_all(&[mixed.0, mixed.1]).expect("estimated");
    assert!((e.x - 0.5).abs() < 1e-12 && (e.y - 0.3).abs() < 1e-12, "{e:?}");
    let _ = fs::remove_dir_all(first);
    let _ = fs::remove_dir_all(second);
}

//! A BAI without its optional per-contig counts (the pseudo-bin 37450), for every contig or
//! some, is not one of a BAM without reads (LRP-005).

use std::fs;
use std::path::{Path, PathBuf};

use snv_hts::{Alignments, ReadSource};

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../lrpipe-sv/tests/data")
}

/// `bai` with the pseudo-bins of the contigs `strip` holds for left out.
fn without_counts(bai: &[u8], strip: impl Fn(usize) -> bool) -> Vec<u8> {
    let int = |o: usize| u32::from_le_bytes(bai[o..o + 4].try_into().expect("four bytes")) as usize;
    let mut out = bai[..8].to_vec();
    let mut o = 8;
    for tid in 0..int(4) {
        let (bins, mut kept) = (int(o), Vec::new());
        o += 4;
        let mut n = 0u32;
        for _ in 0..bins {
            let len = 8 + 16 * int(o + 4);
            if int(o) != 37450 || !strip(tid) {
                kept.extend_from_slice(&bai[o..o + len]);
                n += 1;
            }
            o += len;
        }
        out.extend(n.to_le_bytes());
        out.extend(kept);
        let len = 4 + 8 * int(o);
        out.extend_from_slice(&bai[o..o + len]);
        o += len;
    }
    out.extend_from_slice(&bai[o..]);
    out
}

#[test]
fn an_index_without_counts_does_not_mean_no_reads() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("index-counts");
    fs::create_dir_all(&dir).expect("scratch directory");
    let (bam, reference) = (dir.join("sample1.bam"), data().join("reference.window.fa"));
    fs::copy(data().join("sample1.bam"), &bam).expect("BAM copied");
    let bai = fs::read(data().join("sample1.bam.bai")).expect("BAI read");
    fs::write(dir.join("sample1.bam.bai"), &bai).expect("BAI written");

    let with = Alignments::open(&bam, &reference).expect("opens");
    let names = with.target_names().to_vec();
    let (reads, empty) = (names[16].as_str(), names[0].as_str());
    assert_eq!(with.mapped_reads(reads), Some(562));
    assert_eq!(with.mapped_reads(empty), Some(0), "counted elsewhere, so none here");

    fs::write(dir.join("sample1.bam.bai"), without_counts(&bai, |_| true)).expect("BAI written");
    let without = Alignments::open(&bam, &reference).expect("opens");
    assert_eq!(without.mapped_reads(reads), None, "unknown, not zero");
    assert_eq!(without.mapped_reads(empty), Some(0), "no bins, so no reads");

    // Counts for every contig but the one with reads: still unknown there.
    fs::write(dir.join("sample1.bam.bai"), without_counts(&bai, |tid| tid == 16)).expect("BAI written");
    let mut partial = Alignments::open(&bam, &reference).expect("opens");
    assert_eq!(partial.mapped_reads(reads), None, "unknown, though other contigs are counted");
    assert_eq!(partial.mapped_reads(empty), Some(0));
    let (mut records, mut n) = (partial.query(reads).expect("queries"), 0);
    while records.next_read().expect("reads").is_some() {
        n += 1;
    }
    assert_eq!(n, 562, "the reads are there");
}

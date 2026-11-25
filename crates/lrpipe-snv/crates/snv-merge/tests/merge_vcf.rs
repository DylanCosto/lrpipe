//! `merge_vcf` on small synthetic inputs. The expected outputs are what the
//! original's `MergeVcf.py` writes for the same files.

use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use flate2::Compression;
use flate2::write::GzEncoder;
use snv_merge::{Error, Gender, MergeMode, MergeOptions, RefCalls, merge_vcf};

const FA_HEADER: &str =
    "##fileformat=VCFv4.2\n##INFO=<ID=F>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
const PILEUP_HEADER: &str =
    "##fileformat=VCFv4.2\n##INFO=<ID=P>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("merge_vcf").join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory can be created");
    dir
}

fn record(contig: &str, pos: u32, alt: &str, qual: &str, gt: &str, info: &str) -> String {
    let filter = if alt == "." { "RefCall" } else { "PASS" };
    format!("{contig}\t{pos}\t.\tA\t{alt}\t{qual}\t{filter}\t{info}\tGT:GQ\t{gt}:9\n")
}

fn gzip(path: &Path, text: &str) {
    let mut encoder = GzEncoder::new(fs::File::create(path).expect("file can be created"), Compression::default());
    encoder.write_all(text.as_bytes()).expect("gzip written");
    encoder.finish().expect("gzip finished");
}

/// Pileup calls (gzipped) and full-alignment calls (plain) for chrX.
fn inputs(dir: &Path) -> MergeOptions {
    let pileup = [
        record("chrX", 10, "C", "20", "0/1", "P"),
        record("chrX", 20, "C", "20", "1/1", "P"),
        record("chrX", 30, "C", "1.5", "0/1", "P"),
        record("chrX", 40, ".", "9", "0/0", "P"),
        record("chrX", 50, "C", "20", "1/1", "P"),
        record("chr1", 15, "C", "20", "1/1", "P"),
    ]
    .concat();
    gzip(&dir.join("pileup.vcf.gz"), &format!("{PILEUP_HEADER}{pileup}"));
    // Unsorted, with a RefCall at 20 that hides the pileup call there.
    let full_alignment = [
        record("chrX", 45, "G", "2", "1/1", "F"),
        record("chrX", 20, ".", "7", "0/0", "F"),
        record("chrX", 5, "G", "30", "0/1", "F"),
    ]
    .concat();
    fs::write(dir.join("fa.vcf"), format!("{FA_HEADER}{full_alignment}")).expect("written");
    let mut options = MergeOptions::new(dir.join("pileup.vcf.gz"), dir.join("fa.vcf"), dir.join("merge_chrX.vcf"));
    options.contig = Some("chrX".to_owned());
    options
}

fn output(options: &MergeOptions) -> String {
    fs::read_to_string(&options.output).expect("output is readable")
}

#[test]
fn full_alignment_calls_replace_pileup_calls() {
    let dir = scratch("candidates");
    let options = inputs(&dir);
    let summary = merge_vcf(&options).expect("merge succeeds");
    let want = [
        record("chrX", 10, "C", "20", "0/1", "P"),
        record("chrX", 30, "C", "1.5", "0/1", "P").replace("PASS", "LowQual"),
        record("chrX", 45, "G", "2", "1/1", "F").replace("PASS", "LowQual"),
        record("chrX", 5, "G", "30", "0/1", "F"),
        record("chrX", 50, "C", "20", "1/1", "P"),
    ]
    .concat();
    assert_eq!(output(&options), format!("{FA_HEADER}{want}"));
    assert_eq!((summary.records, summary.pileup, summary.full_alignment), (5, 3, 2));
}

#[test]
fn male_haploid_with_ref_calls() {
    let dir = scratch("male");
    let mut options = inputs(&dir);
    options.gender = Gender::Male;
    options.ref_calls = RefCalls::Show;
    merge_vcf(&options).expect("merge succeeds");
    // Heterozygous calls are emptied; 1/1 and 0/0 become 1 and 0. The full-alignment rows
    // stay in file order, interleaved with the pileup rows as heapq.merge compares them.
    let want = [
        record("chrX", 40, ".", "9", "0", "P"),
        record("chrX", 45, "G", "2", "1", "F").replace("PASS", "LowQual"),
        record("chrX", 20, ".", "7", "0", "F"),
        record("chrX", 50, "C", "20", "1", "P"),
    ]
    .concat();
    assert_eq!(output(&options), format!("{FA_HEADER}{want}"));
}

#[test]
fn candidate_regions_choose_the_source() {
    let dir = scratch("regions");
    let mut options = inputs(&dir);
    let beds = dir.join("candidate_bed");
    fs::create_dir_all(&beds).expect("created");
    // Padded by 33 at both ends: [33, 46) covers 40 and 45.
    fs::write(beds.join("chrX.0_1"), "chrX\t0\t79\n").expect("written");
    fs::write(beds.join("chr1.0_1"), "chr1\t0\t100\n").expect("written");
    options.mode = MergeMode::Regions { candidate_beds: beds.clone() };
    options.ref_calls = RefCalls::Show;
    let summary = merge_vcf(&options).expect("merge succeeds");
    // Pileup calls outside the region, full-alignment calls inside it, by position.
    let want = [
        record("chrX", 10, "C", "20", "0/1", "P"),
        record("chrX", 20, "C", "20", "1/1", "P"),
        record("chrX", 30, "C", "1.5", "0/1", "P").replace("PASS", "LowQual"),
        record("chrX", 45, "G", "2", "1/1", "F").replace("PASS", "LowQual"),
        record("chrX", 50, "C", "20", "1/1", "P"),
    ]
    .concat();
    assert_eq!(output(&options), format!("{PILEUP_HEADER}{want}"));
    let bed = beds.join("full_aln_regions_chrX");
    assert_eq!(summary.candidate_bed.as_deref(), Some(bed.as_path()));
    assert_eq!(fs::read_to_string(bed).expect("written"), "chrX\t0\t79\n");
}

#[test]
fn a_header_only_input_without_its_last_newline_keeps_the_records_apart() {
    // SortVcf's header-only VCF (`get_header`) has no newline after its `#CHROM` line: lrpipe-snv
    // merges a contig without full-alignment calls with one (merging each contig as its calls
    // finish, on the GPU with final phasing), and its first pileup record joined the `#CHROM`
    // line, which lrpipe-cnv then refused (HG002 ONT chr20).
    let dir = scratch("unterminated");
    let options = inputs(&dir);
    fs::write(dir.join("fa.vcf"), FA_HEADER.trim_end_matches('\n')).expect("written");
    merge_vcf(&options).expect("merges");
    let text = output(&options);
    let chrom = text.lines().find(|l| l.starts_with("#CHROM")).expect("a #CHROM line");
    assert_eq!(chrom, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE");
    assert!(text.lines().filter(|l| !l.starts_with('#')).all(|l| l.split('\t').count() == 10), "{text}");
    assert_eq!(text.lines().filter(|l| l.starts_with("chrX\t")).count(), 4, "{text}");

    // The pileup VCF's header in the regions mode too.
    let mut options = inputs(&dir);
    gzip(&dir.join("pileup.vcf.gz"), PILEUP_HEADER.trim_end_matches('\n'));
    fs::write(dir.join("fa.vcf"), format!("{FA_HEADER}{}", record("chrX", 45, "G", "30", "0/1", "F")))
        .expect("written");
    let beds = dir.join("beds");
    fs::create_dir_all(&beds).expect("directory");
    // Padded by 33 at both ends: [33, 67) covers 45.
    fs::write(beds.join("chrX.0_1"), "chrX\t0\t100\n").expect("written");
    options.mode = MergeMode::Regions { candidate_beds: beds };
    merge_vcf(&options).expect("merges");
    let text = output(&options);
    assert!(text.contains("SAMPLE\nchrX\t45\t"), "{text}");
}

#[test]
fn missing_inputs_read_as_empty() {
    let dir = scratch("missing");
    let mut options = inputs(&dir);
    options.full_alignment_vcf = dir.join("absent.vcf.gz");
    let summary = merge_vcf(&options).expect("merge succeeds");
    // No header: the original takes it from the full-alignment VCF.
    let want = [
        record("chrX", 10, "C", "20", "0/1", "P"),
        record("chrX", 20, "C", "20", "1/1", "P"),
        record("chrX", 30, "C", "1.5", "0/1", "P").replace("PASS", "LowQual"),
        record("chrX", 50, "C", "20", "1/1", "P"),
    ]
    .concat();
    assert_eq!(output(&options), want);
    assert_eq!(summary.full_alignment, 0);
}

#[test]
fn invalid_modes_fail() {
    let dir = scratch("errors");
    let mut options = inputs(&dir);
    options.mode = MergeMode::Regions { candidate_beds: dir.join("absent") };
    assert!(matches!(merge_vcf(&options), Err(Error::InputDir { .. })));
    options.mode = MergeMode::Regions { candidate_beds: dir.clone() };
    options.contig = None;
    assert!(matches!(merge_vcf(&options), Err(Error::ContigUnset)));
    let mut options = inputs(&dir);
    fs::write(dir.join("fa.vcf"), format!("{FA_HEADER}chrX\t1x\t.\tA\tC\t5\tPASS\n")).expect("written");
    match merge_vcf(&options) {
        Err(Error::Line { line: 4, .. }) => {}
        other => panic!("expected a line error, got {other:?}"),
    }
    options.contig = Some("chrY".to_owned());
    assert!(merge_vcf(&options).is_ok(), "rows of other contigs are not parsed");
}

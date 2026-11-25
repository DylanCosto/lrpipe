//! `sort_vcf` and `compress_and_index` on small synthetic inputs.

use std::fs;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use snv_merge::{
    ChunkLayout, EmptyReason, Error, Haploid, IndexStatus, Postprocess, SortOptions, compress_and_index, sort_vcf,
};

const HEADER: &str = "##fileformat=VCFv4.2\n##source=Clair3\n##clair3_version=2.0.3\n##cmdline=run\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("merge").join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("in")).expect("scratch directory can be created");
    dir
}

fn gunzip(path: &Path) -> String {
    let mut out = String::new();
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    MultiGzDecoder::new(file).read_to_string(&mut out).expect("valid gzip text");
    out
}

fn record(contig: &str, pos: u32, alt: &str, gt: &str) -> String {
    let filter = if alt == "." { "RefCall" } else { "PASS" };
    format!("{contig}\t{pos}\t.\tA\t{alt}\t5.5\t{filter}\tP\tGT:GQ:DP:AD:AF\t{gt}:5:3:1,2:0.6667\n")
}

fn options(dir: &Path, contigs: &str) -> SortOptions {
    fs::write(dir.join("CONTIGS"), contigs).expect("contigs written");
    let mut options = SortOptions::new(dir.join("in"), dir.join("pileup.vcf"));
    options.prefix = Some("pileup".to_owned());
    options.contigs = Some(dir.join("CONTIGS"));
    options
}

#[test]
fn merges_chunks_in_contig_order() {
    let dir = scratch("order");
    let input = dir.join("in");
    // A later row at the same position replaces an earlier one.
    let chunk1 = [
        record("chr1", 30, "C", "0/1"),
        record("chr1", 10, "C", "1/1"),
        record("chr1", 20, ".", "0/0"),
        record("chr1", 30, "G", "1/1"),
    ];
    fs::write(input.join("pileup_chr1_1.vcf"), format!("{HEADER}{}", chunk1.concat())).expect("written");
    // Another chunk's extra header line is kept.
    let chunk2 = record("chr1", 40, "G", "1/1");
    fs::write(input.join("pileup_chr1_2.vcf"), format!("##extra=1\n{HEADER}{chunk2}")).expect("written");
    fs::write(input.join("pileup_chr11_1.vcf"), format!("{HEADER}{}", record("chr11", 5, "C", "0/1")))
        .expect("written");
    fs::write(input.join("pileup_chrUn_1.vcf"), format!("{HEADER}{}", record("chrUn", 7, "C", "0/1")))
        .expect("written");
    fs::write(input.join("full_alignment_chr1_1.vcf"), format!("{HEADER}{}", record("chr1", 1, "C", "0/1")))
        .expect("written");
    let summary = sort_vcf(&options(&dir, "chrUn\nchr11\nchr1\n")).expect("sort_vcf succeeds");

    assert_eq!(summary.index, IndexStatus::Built);
    assert_eq!(summary.records, 6);
    assert!(!dir.join("pileup.vcf").exists());
    assert!(dir.join("pileup.vcf.gz.tbi").exists());
    let text = gunzip(&summary.vcf);
    let (header, body) = text.split_at(text.find("chr1\t").expect("a chr1 record"));
    assert!(header.contains("##extra=1\n") && header.matches("#CHROM").count() == 1, "{header}");
    let want = [
        record("chr1", 10, "C", "1/1"),
        record("chr1", 20, ".", "0/0"),
        record("chr1", 30, "G", "1/1"),
        record("chr1", 40, "G", "1/1"),
        record("chr11", 5, "C", "0/1"),
        record("chrUn", 7, "C", "0/1"),
    ];
    assert_eq!(body, want.concat());
}

#[test]
fn excludes_appledouble_sidecars_when_sorting_by_suffix() {
    let dir = scratch("appledouble");
    let input = dir.join("in");
    let row = record("chr1", 10, "C", "0/1");
    fs::write(input.join("pileup_chr1_1.vcf"), format!("{HEADER}{row}")).expect("write VCF");
    fs::write(input.join("._pileup_chr1_1.vcf"), [0x00, 0x05, 0x16, 0x07, 0xff]).expect("write sidecar");
    let hidden_row = record("chr1", 20, "G", "0/1");
    fs::write(input.join("._pileup_chr1_2.vcf"), format!("{HEADER}{hidden_row}")).expect("write hidden VCF");
    let mut opts = options(&dir, "chr1\n");
    opts.prefix = None;
    opts.layout = ChunkLayout::Unordered;
    let summary = sort_vcf(&opts).expect("sort VCFs without reading AppleDouble metadata");
    assert_eq!(summary.index, IndexStatus::Built);
    assert_eq!(summary.records, 2);
    assert_eq!(gunzip(&summary.vcf), format!("{HEADER}{row}{hidden_row}"));
}

#[test]
fn pileup_only_processes_rows_in_both_layouts() {
    let dir = scratch("pileup_only");
    let rows = [record("chrX", 10, "C", "0/1"), record("chrX", 20, ".", "0/0"), record("chrX", 30, "C", "1/1")];
    fs::write(dir.join("in/pileup_chrX_1.vcf"), format!("{HEADER}{}", rows.concat())).expect("written");
    let mut opts = options(&dir, "chrX\n");
    opts.postprocess = Postprocess::PileupOnly;
    opts.haploid = Haploid::Precise;
    opts.qual = 10;

    // Both layouts write the processed rows: the RefCall and the het call (empty when haploid)
    // go, the hom call becomes haploid and LowQual. Per contig, the original wrote the rows
    // unchanged, dropping only the RefCall (LRP-012).
    let want = rows[2].replace("1/1:", "1:").replace("PASS", "LowQual");
    for layout in [ChunkLayout::ByContig, ChunkLayout::Unordered] {
        opts.layout = layout;
        let text = gunzip(&sort_vcf(&opts).expect("sort_vcf succeeds").vcf);
        assert!(text.ends_with(&format!("SAMPLE\n{want}")), "{layout:?}: {text}");
    }
}

#[test]
fn empty_inputs_give_the_reference_header() {
    let dir = scratch("empty");
    fs::write(dir.join("ref.fa"), ">chr1\nACGT\n").expect("written");
    fs::write(dir.join("ref.fa.fai"), "chr1\t4\t6\t4\t5\n").expect("written");
    fs::write(dir.join("in/other.vcf"), "x\n").expect("written");
    let mut opts = options(&dir, "chr1\n");
    opts.reference = Some(dir.join("ref.fa"));
    let summary = sort_vcf(&opts).expect("sort_vcf succeeds");
    assert_eq!(summary.empty, Some(EmptyReason::NoPrefixMatch));
    let text = gunzip(&summary.vcf);
    assert!(
        text.ends_with("##contig=<ID=chr1,length=4>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    );
    assert_eq!(summary.copied_to, Some(dir.join("merge_output.vcf.gz")));
    assert_eq!(gunzip(&dir.join("merge_output.vcf.gz")), text);

    fs::write(dir.join("in/pileup_chr1_1.vcf"), HEADER).expect("written");
    assert_eq!(sort_vcf(&opts).expect("sort_vcf succeeds").empty, Some(EmptyReason::NoRecords));
    fs::write(dir.join("in/pileup_chr1_1.vcf"), "").expect("written");
    assert_eq!(sort_vcf(&opts).expect("sort_vcf succeeds").empty, Some(EmptyReason::NoLines));
}

#[test]
fn original_failures_are_errors() {
    let dir = scratch("errors");
    fs::write(dir.join("in/pileup_chr1_1.vcf"), HEADER.replace("##cmdline=run\n", "")).expect("written");
    let opts = options(&dir, "chr1\n");
    assert!(matches!(sort_vcf(&opts), Err(Error::CmdFileUnset)));

    let mut no_contigs = opts.clone();
    no_contigs.contigs = Some(dir.join("missing"));
    assert!(matches!(sort_vcf(&no_contigs), Err(Error::ContigsFile { .. })));
    let mut no_dir = opts.clone();
    no_dir.input_dir = dir.join("missing");
    assert!(matches!(sort_vcf(&no_dir), Err(Error::InputDir { .. })));

    fs::write(dir.join("in/pileup_chr1_1.vcf"), format!("{HEADER}chr1\tx\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n"))
        .expect("written");
    assert!(matches!(sort_vcf(&opts), Err(Error::Line { line: 6, .. })));
}

/// Whether the original's `tabix` (htslib 1.24) indexes a VCF with `middle` between two records,
/// as observed with `tabix -f -p vcf`. Without a REF column it does (htslib 1.19 did not).
const TABIX_CASES: [(&str, bool); 16] = [
    ("chr2\n", false),
    ("chr20\t200\n", true),
    ("chr20\t200\t.\n", true),
    ("chr20\t200\t.\tA\n", true),
    ("chr20\tabc\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", false),
    ("chr20\t200x\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", true),
    ("chr20\t0\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", false),
    ("\n", false),
    ("chr20\t200\t.\tA\tC\t5\tPASS\tEND=abc\tGT\t0/1\n", true),
    ("chr20\t50\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", false),
    ("chr20\t\t150\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", true),
    ("chr20\t\tx\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", false),
    ("chr20\t0x96\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", true),
    ("chr20\t200\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\r\n", true),
    ("chr21\t200\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", false),
    ("#comment\n", true),
];

#[test]
fn indexing_fails_where_the_original_tabix_fails() {
    let dir = scratch("tabix");
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
    let (first, last) = ("chr20\t100\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n", "chr20\t300\t.\tA\tC\t5\tPASS\tP\tGT\t0/1\n");
    for (i, (middle, indexed)) in TABIX_CASES.iter().enumerate() {
        let vcf = dir.join(format!("case{i}.vcf"));
        fs::write(&vcf, format!("{header}{first}{middle}{last}")).expect("written");
        let (gz, status) = compress_and_index(&vcf).expect("compresses");
        assert_eq!(status == IndexStatus::Built, *indexed, "{middle:?}");
        assert_eq!(dir.join(format!("case{i}.vcf.gz.tbi")).exists(), *indexed, "{middle:?}");
        assert_eq!(gunzip(&gz), format!("{header}{first}{middle}{last}"));
    }
}

//! The `lrpipe-phase` command line, also `lrpipe phase` and `lrpipe haplotag`.

use std::ffi::OsString;
use std::io::Write;

use anyhow::{Result, bail};
use clap::{Args, Parser, Subcommand};

use crate::compare::{CompareOptions, run_compare};
use crate::genotype::{GenotypeOptions, load_kmerald_costs, run_genotype};
use crate::hapcut2vcf::run_hapcut2vcf;
use crate::haplotag::{HaplotagOptions, SupplementaryStrategy, run_haplotag};
use crate::haplotagphase::{HaplotagphaseOptions, run_haplotagphase};
use crate::phase::{Algorithm, PhaseOptions, run_phase};
use crate::split::{SplitOptions, run_split};
use crate::stats::{StatsOptions, run_stats};
use crate::unphase::run_unphase;
use crate::variants::Realignment;
use crate::vcf::PhaseTag;

#[derive(Parser)]
#[command(
    name = "lrpipe-phase",
    version,
    propagate_version = true,
    about = "Read-based phasing (Rust port of WhatsHap 2.8)"
)]
struct Cli {
    /// Print debug messages
    #[arg(long, global = true)]
    debug: bool,
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Phase variants in a VCF with the WhatsHap algorithm
    #[command(
        after_help = "Other commands: haplotag, genotype, stats, compare, split, unphase, hapcut2vcf and haplotagphase; \
                      `help` lists them."
    )]
    Phase(PhaseArgs),
    /// Tag reads by haplotype
    Haplotag(HaplotagArgs),
    /// Print phasing statistics of a single VCF file
    Stats(StatsArgs),
    /// Compare two or more phased variant files
    Compare(CompareArgs),
    /// Split reads by haplotype
    Split(SplitArgs),
    /// Remove phasing information from a VCF file
    Unphase(UnphaseArgs),
    /// Convert hapCUT output format to VCF
    Hapcut2vcf(Hapcut2vcfArgs),
    /// Phase variants in VCF based on information from haplotagged reads
    Haplotagphase(HaplotagphaseArgs),
    /// Genotype variants
    Genotype(GenotypeArgs),
}

#[derive(Args)]
#[command(args_override_self = true)]
struct GenotypeArgs {
    /// VCF file with variants to be genotyped (can be gzip-compressed)
    #[arg(value_name = "VCF")]
    variant_file: String,
    /// BAM or VCF file(s) with phase information, either through sequencing reads (BAM) or through phased
    /// blocks (VCF)
    #[arg(value_name = "PHASEINPUT")]
    phase_input_files: Vec<String>,
    /// Output VCF file. Add .gz to the file name to get compressed output. If omitted, use standard output.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    /// Reference file. Provide this to detect alleles through re-alignment.
    #[arg(short = 'r', long, value_name = "FASTA")]
    reference: Option<String>,
    /// Reduce coverage to at most MAXCOV.
    #[arg(short = 'H', long = "max-coverage", value_name = "MAXCOV", default_value_t = 15)]
    max_coverage: usize,
    /// Minimum mapping quality
    #[arg(long = "mapping-quality", visible_alias = "mapq", value_name = "QUAL", default_value_t = 20)]
    mapping_quality: i32,
    #[arg(long = "indels", hide = true)]
    indels_used: bool,
    /// Genotype only SNVs
    #[arg(long = "only-snvs")]
    only_snvs: bool,
    /// Ignore read groups in BAM header and assume all reads come from the same sample.
    #[arg(long = "ignore-read-groups")]
    ignore_read_groups: bool,
    /// Name of a sample to genotype. If not given, all samples in the input VCF are genotyped. Can be used
    /// multiple times.
    #[arg(long = "sample", value_name = "SAMPLE", action = clap::ArgAction::Append)]
    samples: Vec<String>,
    /// Name of chromosome to genotyped. If not given, all chromosomes in the input VCF are genotyped. Can be
    /// used multiple times.
    #[arg(long = "chromosome", value_name = "CHROMOSOME", action = clap::ArgAction::Append)]
    chromosomes: Vec<String>,
    /// Name of chromosome not to be genotyped.
    #[arg(long = "exclude-chromosome", action = clap::ArgAction::Append)]
    excluded_chromosomes: Vec<String>,
    /// Phred scaled error probability threshold used for genotyping. Must be at least 0. If error probability of
    /// genotype is higher, genotype ./. is output.
    #[arg(
        long = "gt-qual-threshold",
        value_name = "GTQUALTHRESHOLD",
        default_value_t = 0.0,
        allow_negative_numbers = true
    )]
    gt_qual_threshold: f64,
    /// Skip initial prior genotyping and use uniform priors.
    #[arg(long = "no-priors")]
    nopriors: bool,
    /// output prior genotype likelihoods to the given file.
    #[arg(short = 'p', long)]
    prioroutput: Option<String>,
    /// When --reference is used, extend alignment by this many bases to left and right when realigning.
    #[arg(long, value_name = "OVERHANG", default_value_t = 10)]
    overhang: i64,
    /// This constant is used to regularize the priors.
    #[arg(long, value_name = "CONSTANT", default_value_t = 0.0, allow_negative_numbers = true)]
    constant: f64,
    /// When detecting alleles through re-alignment, use affine gap costs (EXPERIMENTAL).
    #[arg(long = "affine-gap")]
    affine_gap: bool,
    /// gap starting penalty in case affine gap costs are used (default: 10).
    #[arg(long = "gap-start", value_name = "GAPSTART")]
    gap_start: Option<f64>,
    /// gap extend penalty in case affine gap costs are used (default: 7).
    #[arg(long = "gap-extend", value_name = "GAPEXTEND")]
    gap_extend: Option<f64>,
    /// mismatch cost in case affine gap costs are used (default: 15)
    #[arg(long, value_name = "MISMATCH", default_value_t = 15.0)]
    mismatch: f64,
    /// Use pedigree information in PED file to improve genotyping (EXPERIMENTAL).
    #[arg(long, value_name = "PED/FAM")]
    ped: Option<String>,
    /// Recombination rate in cM/Mb (used with --ped).
    #[arg(long, value_name = "RECOMBRATE", default_value_t = 1.26)]
    recombrate: f64,
    /// File with genetic map (used with --ped) to be used instead of constant recombination rate.
    #[arg(long, value_name = "FILE")]
    genmap: Option<String>,
    /// Only work on samples mentioned in the provided PED file.
    #[arg(long = "use-ped-samples")]
    use_ped_samples: bool,
    /// Use kmerald for detecting alleles through re-alignment.
    #[arg(long = "use-kmerald")]
    use_kmerald: bool,
    /// Error model based costs used by kmerald during re-alignment.
    #[arg(long = "kmeralign-costs", value_name = "COSTS")]
    kmeralign_costs_path: Option<String>,
    /// kmer size used by kmerald during re-alignment.
    #[arg(long = "kmer-size", value_name = "KMER")]
    kmer_size: Option<u32>,
    /// Gap penalty used by kmerald during re-alignment.
    #[arg(long = "kmerald-gappenalty", value_name = "GAP")]
    kmerald_gappenalty: Option<f64>,
    /// Consider this many bases on the left and side of a variant position for kmerald based re-alignment.
    #[arg(long = "kmerald-window", value_name = "WINDOW")]
    kmerald_window: Option<f64>,
    /// Number of chromosomes to genotype at a time (lrpipe-phase only; the output is identical for any value)
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,
}

#[derive(Args)]
struct UnphaseArgs {
    /// VCF file. Use "-" to read from standard input
    #[arg(value_name = "VCF")]
    vcf: String,
}

#[derive(Args)]
#[command(args_override_self = true)]
struct Hapcut2vcfArgs {
    /// Output VCF file. If omitted, use standard output.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    /// VCF file
    #[arg(value_name = "VCF")]
    vcf: String,
    /// hapCUT result file
    #[arg(value_name = "HAPCUT-RESULT")]
    hapcut: String,
}

#[derive(Args)]
#[command(args_override_self = true)]
struct HaplotagphaseArgs {
    /// Output file. If omitted, use standard output.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    /// Reference file. Must be accompanied by .fai index (create with samtools faidx)
    #[arg(short = 'r', long, value_name = "FASTA")]
    reference: Option<String>,
    /// Threshold percentage for qualities. If the percentage of votes for the variant is less than this value,
    /// the algorithm does not assign any information to the variant.
    #[arg(
        short = 'g',
        long = "gap-threshold",
        value_name = "PERCENT",
        default_value_t = 70,
        allow_negative_numbers = true
    )]
    gap_threshold: i64,
    /// Ignore variants within homopolymers longer than the cut value.
    #[arg(short = 'c', long = "cut-poly", value_name = "LENGTH", default_value_t = 10, allow_negative_numbers = true)]
    cut_poly: i64,
    /// Add phasing information only to indels.
    #[arg(short = 'i', long = "only-indels")]
    only_indels: bool,
    /// Name of a sample to phase. If not given, all samples in the input VCF are phased. Can be used multiple times.
    #[arg(long = "sample", value_name = "SAMPLE", action = clap::ArgAction::Append)]
    samples: Vec<String>,
    /// Ignore read groups in BAM/CRAM header and assume all reads come from the same sample.
    #[arg(long = "ignore-read-groups")]
    ignore_read_groups: bool,
    /// Name of chromosome to phase. If not given, all chromosomes in the input VCF are phased. Can be used
    /// multiple times.
    #[arg(long = "chromosome", value_name = "CHROMOSOME", action = clap::ArgAction::Append)]
    chromosomes: Vec<String>,
    /// Ignore multiallelic variants.
    #[arg(long = "no-mav")]
    no_mav: bool,
    /// Name of chromosome not to phase.
    #[arg(long = "exclude-chromosome", action = clap::ArgAction::Append)]
    excluded_chromosomes: Vec<String>,
    /// VCF file with variants to phase (must be gzip-compressed and indexed)
    #[arg(value_name = "VCF")]
    variant_file: String,
    /// BAM/CRAM file with alignments tagged by haplotype and phase set
    #[arg(value_name = "ALIGNMENTS")]
    alignment_file: String,
}

#[derive(Args)]
#[command(args_override_self = true)]
struct SplitArgs {
    /// Output haplotype 1 reads to FILE (.gz supported)
    #[arg(long = "output-h1", value_name = "FILE")]
    output_h1: Option<String>,
    /// Output haplotype 2 reads to FILE (.gz supported)
    #[arg(long = "output-h2", value_name = "FILE")]
    output_h2: Option<String>,
    /// Output haplotype reads to FILE. Use this option as many times as there are haplotypes in the input. The
    /// first -o is used for H1, second for H2 etc.
    #[arg(short = 'o', long = "output", value_name = "FILE", action = clap::ArgAction::Append)]
    outputs: Option<Vec<String>>,
    /// Output file to write untagged reads to (.gz supported)
    #[arg(long = "output-untagged")]
    output_untagged: Option<String>,
    /// Add reads without tag to all (H1, H2, ...) outputs.
    #[arg(long = "add-untagged")]
    add_untagged: bool,
    /// Only consider reads to be tagged if they belong to the largest phased block (in terms of read count) on
    /// their respective chromosome
    #[arg(long = "only-largest-block")]
    only_largest_block: bool,
    /// Only check the haplotype of reads listed in the haplotag list file. Reads (read names) not contained in
    /// this file will be discarded.
    #[arg(long = "discard-unknown-reads")]
    discard_unknown_reads: bool,
    /// Output file to write read lengths histogram to in tab-separated format.
    #[arg(long = "read-lengths-histogram")]
    read_lengths_histogram: Option<String>,
    /// Input FASTQ/BAM file with reads (FASTQ can be gzipped)
    #[arg(value_name = "READS")]
    reads_file: String,
    /// Tab-separated list with (at least) two columns <readname> and <haplotype> (can be gzipped). The
    /// haplotypes have to be named H1, H2, etc. (or none). The output of the "haplotag" command can be used (4
    /// columns), and this is required for using the "--only-largest-block" option.
    #[arg(value_name = "LIST")]
    list_file: String,
}

#[derive(Args)]
#[command(args_override_self = true)]
struct CompareArgs {
    /// Name of the sample to process. If not given, use first sample found in VCF.
    #[arg(long, value_name = "SAMPLE")]
    sample: Option<String>,
    /// Comma-separated list of data set names to be used in the report (in same order as VCFs).
    #[arg(long, value_name = "NAMES")]
    names: Option<String>,
    /// For single sample VCFs, ignore sample name and assume all samples are the same.
    #[arg(long = "ignore-sample-name")]
    ignore_sample_name: bool,
    /// Filename to write comparison results from pair-wise comparison to (tab-separated).
    #[arg(long = "tsv-pairwise", value_name = "TSVPAIRWISE")]
    tsv_pairwise: Option<String>,
    /// Filename to write comparison results from multiway comparison to (tab-separated). Only for diploid VCFs.
    #[arg(long = "tsv-multiway", value_name = "TSVMULTIWAY")]
    tsv_multiway: Option<String>,
    /// Only process SNVs and ignore all other variants.
    #[arg(long = "only-snvs")]
    only_snvs: bool,
    /// Write BED file with switch error positions to given filename. Only for diploid VCFs.
    #[arg(long = "switch-error-bed")]
    switch_error_bed: Option<String>,
    /// Not supported by lrpipe-phase (WhatsHap writes a PDF with matplotlib)
    #[arg(long = "plot-blocksizes", hide = true)]
    plot_blocksizes: Option<String>,
    /// Not supported by lrpipe-phase (WhatsHap writes a PDF with matplotlib)
    #[arg(long = "plot-sum-of-blocksizes", hide = true)]
    plot_sum_of_blocksizes: Option<String>,
    /// Write position-wise agreement of longest joint blocks in each chromosome to tab-separated file. Only for
    /// diploid VCFs.
    #[arg(long = "longest-block-tsv")]
    longest_block_tsv: Option<String>,
    /// The ploidy of the sample(s).
    #[arg(short = 'p', long, value_name = "PLOIDY", default_value_t = 2, allow_negative_numbers = true)]
    ploidy: i64,
    /// At least two phased variant files (VCF or BCF) to be compared.
    #[arg(value_name = "VCF/BCF", required = true)]
    vcf: Vec<String>,
}

#[derive(Args)]
#[command(args_override_self = true)]
struct StatsArgs {
    /// Write phased blocks as GTF with each block represented as a 'gene'. If blocks are interleaved or
    /// nested, they are split into multiple 'exons'.
    #[arg(long, value_name = "FILE")]
    gtf: Option<String>,
    /// Write list of all blocks to FILE (one block per line). Nested/interleaved blocks are not split.
    #[arg(long = "block-list", value_name = "FILE")]
    block_list: Option<String>,
    /// Name of the sample to process. If not given, use first sample found in VCF.
    #[arg(long, value_name = "SAMPLE")]
    sample: Option<String>,
    /// Override chromosome lengths in VCF with those from FILE (one line per chromosome, tab separated
    /// '<chr> <length>'). Lengths are used to compute NG50 values.
    #[arg(long = "chr-lengths", value_name = "FILE")]
    chr_lengths: Option<String>,
    /// Write statistics in tab-separated value format to FILE
    #[arg(long, value_name = "FILE")]
    tsv: Option<String>,
    /// Only process SNVs and ignore all other variants.
    #[arg(long = "only-snvs")]
    only_snvs: bool,
    /// Name of chromosome(s) to process. If not given, all chromosomes in the input VCF are considered.
    /// Can be used multiple times and accepts a comma-separated list.
    #[arg(long = "chromosome", value_name = "CHROMOSOME", action = clap::ArgAction::Append)]
    chromosomes: Vec<String>,
    /// Phased VCF file
    #[arg(value_name = "VCF")]
    vcf: String,
}

#[derive(Args)]
struct HaplotagArgs {
    /// VCF file with phased variants (must be gzip-compressed and indexed)
    #[arg(value_name = "VCF")]
    variant_file: String,
    /// BAM/CRAM file with alignments to be tagged by haplotype
    #[arg(value_name = "ALIGNMENTS")]
    alignment_file: String,
    /// Output file. If omitted, use standard output.
    #[arg(short = 'o', long)]
    output: Option<String>,
    /// Reference file. Must be accompanied by .fai index (create with samtools faidx)
    #[arg(short = 'r', long, value_name = "FASTA")]
    reference: Option<String>,
    /// Detect alleles without requiring a reference, at the expense of phasing quality
    #[arg(long)]
    no_reference: bool,
    /// Region(s) of interest (chrom:start-end, chrom, or chrom:start)
    #[arg(long = "regions", value_name = "REGION")]
    regions: Option<Vec<String>>,
    /// Region(s) of interest from a BED file in place of --regions, sorted and merged: a read
    /// reaching several is written once (lrpipe-phase's own option)
    #[arg(long = "regions-bed", value_name = "BED", conflicts_with = "regions")]
    regions_bed: Option<String>,
    /// Ignore linkage information stored in BX tags of the reads.
    #[arg(long)]
    ignore_linked_read: bool,
    /// Assume reads with identical BX tags belong to different read clouds if their distance is larger than this
    #[arg(
        short = 'd',
        long = "linked-read-distance-cutoff",
        value_name = "LINKEDREADDISTANCE",
        default_value_t = 50000
    )]
    linked_read_distance_cutoff: i64,
    /// Ignore read groups in BAM/CRAM header and assume all reads come from the same sample.
    #[arg(long)]
    ignore_read_groups: bool,
    /// Name of a sample to phase. Can be used multiple times.
    #[arg(long = "sample", value_name = "SAMPLE")]
    given_samples: Option<Vec<String>>,
    /// Write assignments of read names to haplotypes (tab separated) to given file (.gz: gzipped)
    #[arg(long = "output-haplotag-list", value_name = "HAPLOTAG_LIST")]
    haplotag_list: Option<String>,
    /// How to tag supplementary alignments: skip, copy-primary (if value omitted), independent-or-skip, independent-or-copy-primary
    #[arg(long = "tag-supplementary", num_args = 0..=1, default_missing_value = "copy-primary", default_value = "skip",
          value_parser = ["skip", "copy-primary", "independent-or-skip", "independent-or-copy-primary"])]
    tag_supplementary: String,
    /// Maximum distance between supplementary and primary alignment for copying tags
    #[arg(long = "supplementary-distance", default_value_t = 100_000)]
    supplementary_distance_threshold: i64,
    /// Allow strand mismatch between supplementary and primary alignment when copying tags
    #[arg(long = "no-supplementary-strand-match")]
    no_supplementary_strand_match: bool,
    /// Ploidy
    #[arg(long, default_value_t = 2)]
    ploidy: usize,
    /// Skip reads that map to a contig that does not exist in the VCF
    #[arg(long)]
    skip_missing_contigs: bool,
    /// Write reads that map to a contig that does not exist in the VCF untagged, as on a contig
    /// without variants, instead of failing (lrpipe-phase's own option; --skip-missing-contigs
    /// leaves them out of the output)
    #[arg(long, conflicts_with = "skip_missing_contigs")]
    untag_missing_contigs: bool,
    /// Write only the haplotag list (--output-haplotag-list), no alignments (lrpipe-phase's own
    /// option)
    #[arg(long, requires = "haplotag_list", conflicts_with_all = ["output", "write_index"])]
    list_only: bool,
    /// Index the output BAM (--output FILE.bam) as `samtools index` does: FILE.bam.bai, or
    /// FILE.bam.csi for contigs longer than a BAI allows (lrpipe-phase's own option)
    #[arg(long, requires = "output")]
    write_index: bool,
    /// Number of threads for compressing the output file (default: the value of --threads)
    #[arg(long = "output-threads", visible_alias = "out-threads")]
    output_threads: Option<usize>,
    /// Compression level of the output BAM/CRAM (0-9; default: htslib's default). Only the
    /// compressed bytes change, not the records.
    #[arg(long = "compression-level", value_name = "LEVEL", value_parser = clap::value_parser!(u32).range(0..=9))]
    compression_level: Option<u32>,
    /// Number of threads for allele detection and input decompression; also the default for
    /// --output-threads (output records are identical for any value)
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,
}

#[derive(Args)]
struct PhaseArgs {
    /// VCF or BCF file with variants to be phased (can be gzip-compressed)
    #[arg(value_name = "VCF")]
    variant_file: String,
    /// BAM or CRAM file(s) with phase information
    #[arg(value_name = "PHASEINPUT")]
    phase_input_files: Vec<String>,
    /// Output VCF file. Add .gz to the file name to get compressed output. If omitted, use standard output.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    /// Reference file. Must be accompanied by .fai index (create with samtools faidx)
    #[arg(short = 'r', long, value_name = "FASTA")]
    reference: Option<String>,
    /// Detect alleles without requiring a reference, at the expense of phasing quality
    #[arg(long)]
    no_reference: bool,
    /// Store phasing information with PS tag (standardized) or HP tag
    #[arg(long, default_value = "PS", value_parser = ["PS", "HP"])]
    tag: String,
    /// Write reads that have been used for phasing to FILE.
    #[arg(long = "output-read-list", value_name = "FILE")]
    read_list_filename: Option<String>,
    /// Phasing algorithm to use
    #[arg(long, default_value = "whatshap", value_parser = ["whatshap", "hapchat", "heuristic"])]
    algorithm: String,
    /// Merge reads which are likely to come from the same haplotype
    #[arg(long = "merge-reads")]
    read_merging: bool,
    /// For the heuristic: Maximum number of memorized intermediate solutions. Larger values increase runtime and memory consumption, but can improve phasing quality. (default: 256)
    #[arg(long = "row-limit", short = 'L', value_name = "ROWLIMIT", allow_negative_numbers = true)]
    row_limit: Option<i32>,
    /// Coverage reduction parameter in the internal core phasing algorithm
    #[arg(long = "internal-downsampling", value_name = "COVERAGE", default_value_t = 15)]
    max_coverage: usize,
    /// Minimum mapping quality
    #[arg(long = "mapping-quality", visible_alias = "mapq", value_name = "QUAL", default_value_t = 20)]
    mapping_quality: i32,
    /// Phase only SNVs
    #[arg(long)]
    only_snvs: bool,
    /// Ignore read groups in BAM/CRAM header and assume all reads come from the same sample.
    #[arg(long)]
    ignore_read_groups: bool,
    /// Name of a sample to phase. Can be used multiple times.
    #[arg(long = "sample", value_name = "SAMPLE")]
    samples: Vec<String>,
    /// Name of chromosome to phase. Can be used multiple times.
    #[arg(long = "chromosome", value_name = "CHROMOSOME")]
    chromosomes: Vec<String>,
    /// Name of chromosome not to phase.
    #[arg(long = "exclude-chromosome", value_name = "CHROMOSOME")]
    excluded_chromosomes: Vec<String>,
    /// Allow switching variants from hetero- to homozygous in an optimal solution
    #[arg(long)]
    distrust_genotypes: bool,
    /// Also work on homozygous variants, which might be turned to heterozygous
    #[arg(long)]
    include_homozygous: bool,
    /// Default genotype quality used as cost of changing a genotype when no genotype likelihoods are available
    #[arg(long, default_value_t = 30)]
    default_gq: i32,
    /// Constant (float) to be used to regularize genotype likelihoods read from input VCF
    #[arg(long)]
    gl_regularizer: Option<f64>,
    /// Write list of changed genotypes to FILE.
    #[arg(long = "changed-genotype-list", value_name = "FILE")]
    gtchange_list_filename: Option<String>,
    /// Use pedigree information in PED file to improve phasing (switches to PedMEC algorithm)
    #[arg(long, value_name = "PED/FAM")]
    ped: Option<String>,
    /// Write putative recombination events to FILE.
    #[arg(long = "recombination-list", value_name = "FILE")]
    recombination_list_filename: Option<String>,
    /// Recombination rate in cM/Mb (used with --ped)
    #[arg(long, value_name = "RECOMBRATE", default_value_t = 1.26)]
    recombrate: f64,
    /// File with genetic map (used with --ped) to be used instead of constant recombination rate
    #[arg(long, value_name = "FILE")]
    genmap: Option<String>,
    /// Do not merge blocks that are not connected by reads (i.e. solely based on genotype status)
    #[arg(long)]
    no_genetic_haplotyping: bool,
    /// Only work on samples mentioned in the provided PED file.
    #[arg(long)]
    use_ped_samples: bool,
    /// Use also supplementary alignments
    #[arg(long)]
    use_supplementary: bool,
    /// Skip supplementary alignments further than DIST bp away from the primary alignment
    #[arg(long = "supplementary-distance", value_name = "DIST", default_value_t = 100_000)]
    supplementary_distance_threshold: i64,
    /// Number of threads; up to half as many chromosomes are phased at a time, fewer if they would exceed
    /// --max-memory (output is identical for any value)
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,
    /// Memory in MB the phasing algorithm may use to keep backtrace tables instead of recomputing them
    /// (output is identical for any value)
    #[arg(long = "dp-memory", value_name = "MB", default_value_t = 1024)]
    dp_memory_mb: u64,
    /// Approximate memory limit in MB: with --threads, fewer chromosomes are phased at a time if their
    /// estimated memory does not fit next to --dp-memory and 96 MB per thread (default: half of the system
    /// memory or cgroup limit; output is identical for any value)
    #[arg(long = "max-memory", value_name = "MB")]
    max_memory_mb: Option<u64>,
}

/// Runs the command line `args` (the program name first) as WhatsHap would.
///
/// # Errors
///
/// When the command fails, as WhatsHap would.
pub fn main(args: Vec<OsString>) -> Result<()> {
    let cli = Cli::parse_from(&args);
    // As WhatsHap records it: the arguments after the program name.
    let args: Vec<String> = args.iter().skip(1).map(|a| a.to_string_lossy().into_owned()).collect();
    let level = if cli.debug { log::LevelFilter::Debug } else { log::LevelFilter::Info };
    env_logger::Builder::new()
        .filter_level(level)
        .format(|buf, record| writeln!(buf, "{}", record.args()))
        .target(env_logger::Target::Stderr)
        .try_init()
        .ok();
    let whatshap_command_line = || format!("(lrpipe-phase {}) {}", env!("CARGO_PKG_VERSION"), args.join(" "));
    match cli.command {
        Command::Genotype(a) => {
            if a.ignore_read_groups && a.ped.is_some() {
                bail!("Option --ignore-read-groups cannot be used together with --ped");
            }
            if a.genmap.is_some() && a.ped.is_none() {
                bail!("Option --genmap can only be used together with --ped");
            }
            if a.genmap.is_some() && a.chromosomes.len() != 1 {
                bail!("Option --genmap can only be used when working on exactly one chromosome (use --chromosome)");
            }
            if a.phase_input_files.is_empty() {
                bail!("Not providing any PHASEINPUT files not allowed for genotyping.");
            }
            if a.gt_qual_threshold < 0.0 {
                bail!("Genotype quality threshold (gt-qual-threshold) must be at least 0.");
            }
            if a.prioroutput.is_some() && a.nopriors {
                bail!("Genotype priors are only computed if --no-priors is NOT set.");
            }
            if a.constant != 0.0 && a.nopriors {
                bail!("--constant can only be used if --no-priors is NOT set..");
            }
            if a.affine_gap && a.reference.is_none() {
                bail!("Option --affine-gap can only be used together with --reference.");
            }
            if a.use_ped_samples && a.ped.is_none() {
                bail!("Option --use-ped-samples can only be used when PED file is provided (--ped).");
            }
            if a.use_ped_samples && !a.samples.is_empty() {
                bail!("Option --use-ped-samples cannot be used together with --samples");
            }
            if a.indels_used {
                log::warn!("Ignoring --indels as indel genotyping is default in WhatsHap 2.0+");
            }
            if a.use_kmerald && a.kmeralign_costs_path.is_none() {
                bail!(
                    "Option --use-kmerald can only be used when the costs to be used for kmer alignment --kmeralign-costs are provided."
                );
            }
            let realignment = if a.use_kmerald {
                // WhatsHap needs all kmerald parameters on the command line
                let (Some(kmer_size), Some(gap_penalty), Some(window)) =
                    (a.kmer_size, a.kmerald_gappenalty, a.kmerald_window)
                else {
                    bail!(
                        "TypeError: --kmer-size, --kmerald-gappenalty and --kmerald-window must be given with --use-kmerald"
                    );
                };
                Realignment::Kmerald {
                    costs: load_kmerald_costs(a.kmeralign_costs_path.as_deref().expect("checked above"))?,
                    kmer_size,
                    gap_penalty: gap_penalty as f32,
                    window: window as i64,
                }
            } else if a.affine_gap {
                // WhatsHap's affine gap alignment takes integers and fails on
                // the float values of --gap-start/--gap-extend
                if a.gap_start.is_some() || a.gap_extend.is_some() {
                    bail!("TypeError: 'float' object cannot be interpreted as an integer (--gap-start/--gap-extend)");
                }
                Realignment::AffineGap { gap_start: 10, gap_extend: 7, mismatch: a.mismatch as f32 }
            } else {
                Realignment::EditDistance
            };
            run_genotype(&GenotypeOptions {
                variant_file: a.variant_file,
                phase_input_files: a.phase_input_files,
                reference: a.reference,
                output: a.output,
                samples: a.samples,
                chromosomes: a.chromosomes,
                excluded_chromosomes: a.excluded_chromosomes,
                ignore_read_groups: a.ignore_read_groups,
                only_snvs: a.only_snvs,
                mapping_quality: a.mapping_quality,
                max_coverage: a.max_coverage,
                nopriors: a.nopriors,
                ped: a.ped,
                recombrate: a.recombrate,
                genmap: a.genmap,
                gt_qual_threshold: a.gt_qual_threshold,
                prioroutput: a.prioroutput,
                constant: a.constant,
                overhang: a.overhang,
                realignment,
                use_ped_samples: a.use_ped_samples,
                command_line: whatshap_command_line(),
                threads: a.threads,
            })
        }
        Command::Unphase(a) => run_unphase(&a.vcf),
        Command::Hapcut2vcf(a) => run_hapcut2vcf(&a.hapcut, &a.vcf, &a.output, &whatshap_command_line()),
        Command::Haplotagphase(a) => run_haplotagphase(&HaplotagphaseOptions {
            variant_file: a.variant_file,
            alignment_file: a.alignment_file,
            output: a.output,
            reference: a.reference,
            gap_threshold: a.gap_threshold,
            cut_poly: a.cut_poly,
            only_indels: a.only_indels,
            samples: a.samples,
            ignore_read_groups: a.ignore_read_groups,
            chromosomes: a.chromosomes,
            excluded_chromosomes: a.excluded_chromosomes,
            mav: !a.no_mav,
            command_line: whatshap_command_line(),
        }),
        Command::Split(a) => {
            if a.output_h1.is_none() && a.output_h2.is_none() && a.outputs.is_none() && a.output_untagged.is_some() {
                bail!("Nothing to be done since neither --output-h1/h2, --outputs/-o nor --output-untagged are given.");
            }
            if (a.output_h1.is_some() || a.output_h2.is_some()) && a.outputs.is_some() {
                bail!("--output-h1/-h2 cannot be used together with --outputs/-o");
            }
            run_split(&SplitOptions {
                reads_file: a.reads_file,
                list_file: a.list_file,
                output_h1: a.output_h1,
                output_h2: a.output_h2,
                outputs: a.outputs,
                output_untagged: a.output_untagged,
                add_untagged: a.add_untagged,
                only_largest_block: a.only_largest_block,
                discard_unknown_reads: a.discard_unknown_reads,
                read_lengths_histogram: a.read_lengths_histogram,
            })
        }
        Command::Compare(a) => {
            if a.vcf.len() < 2 {
                bail!("At least two VCFs need to be given.");
            }
            if a.ploidy < 2 {
                bail!("Ploidy must be > 1.");
            }
            let stdout = std::io::stdout();
            let mut out = std::io::BufWriter::new(stdout.lock());
            let result = run_compare(
                &CompareOptions {
                    vcf: a.vcf,
                    ploidy: a.ploidy as usize,
                    names: a.names,
                    sample: a.sample,
                    ignore_sample_name: a.ignore_sample_name,
                    tsv_pairwise: a.tsv_pairwise,
                    tsv_multiway: a.tsv_multiway,
                    only_snvs: a.only_snvs,
                    switch_error_bed: a.switch_error_bed,
                    plot_blocksizes: a.plot_blocksizes,
                    plot_sum_of_blocksizes: a.plot_sum_of_blocksizes,
                    longest_block_tsv: a.longest_block_tsv,
                },
                &mut out,
            );
            out.flush()?;
            result
        }
        Command::Stats(a) => {
            let stdout = std::io::stdout();
            let mut out = std::io::BufWriter::new(stdout.lock());
            let result = run_stats(
                &StatsOptions {
                    vcf: a.vcf,
                    sample: a.sample,
                    gtf: a.gtf,
                    tsv: a.tsv,
                    block_list: a.block_list,
                    only_snvs: a.only_snvs,
                    chromosomes: a.chromosomes,
                    chr_lengths: a.chr_lengths,
                },
                &mut out,
            );
            // what was printed before an error is still written, as in Python
            out.flush()?;
            result
        }
        Command::Haplotag(a) => {
            if a.reference.is_some() && a.no_reference {
                bail!("Options --reference and --no-reference cannot be used together");
            }
            run_haplotag(HaplotagOptions {
                variant_file: a.variant_file,
                alignment_file: a.alignment_file,
                output: a.output,
                reference: a.reference,
                no_reference: a.no_reference,
                regions: a.regions,
                regions_bed: a.regions_bed,
                ignore_linked_read: a.ignore_linked_read,
                linked_read_distance_cutoff: a.linked_read_distance_cutoff,
                ignore_read_groups: a.ignore_read_groups,
                given_samples: a.given_samples,
                haplotag_list: a.haplotag_list,
                supplementary_strategy: SupplementaryStrategy::parse(&a.tag_supplementary)?,
                supplementary_distance_threshold: a.supplementary_distance_threshold,
                supplementary_strand_match: !a.no_supplementary_strand_match,
                ploidy: a.ploidy,
                skip_missing_contigs: a.skip_missing_contigs,
                untag_missing_contigs: a.untag_missing_contigs,
                list_only: a.list_only,
                write_index: a.write_index,
                output_threads: a.output_threads,
                compression_level: a.compression_level,
                threads: a.threads,
                command_line: std::iter::once("lrpipe-phase")
                    .chain(args.iter().map(String::as_str))
                    .collect::<Vec<_>>()
                    .join(" "),
            })
        }
        Command::Phase(a) => {
            if a.reference.is_some() && a.no_reference {
                bail!("Options --reference and --no-reference cannot be used together");
            }
            if a.phase_input_files.is_empty() && a.ped.is_none() {
                bail!("Not providing any PHASEINPUT files only allowed in --ped mode.");
            }
            if a.max_coverage > 23 {
                bail!("Coverage downsampling parameter must not exceed 23.");
            }
            if a.ignore_read_groups && a.ped.is_some() {
                bail!("Option --ignore-read-groups cannot be used together with --ped");
            }
            if a.genmap.is_some() && a.ped.is_none() {
                bail!("Option --genmap can only be used together with --ped");
            }
            if a.genmap.is_some() && a.chromosomes.len() != 1 {
                bail!("Option --genmap can only be used when working on exactly one chromosome (use --chromosome)");
            }
            if a.include_homozygous && !a.distrust_genotypes {
                bail!("Option --include-homozygous can only be used with --distrust-genotypes.");
            }
            if a.use_ped_samples && a.ped.is_none() {
                bail!("Option --use-ped-samples can only be used when PED file is provided (--ped).");
            }
            if a.use_ped_samples && !a.samples.is_empty() {
                bail!("Option --use-ped-samples cannot be used together with --samples");
            }
            let row_limit = match a.row_limit {
                None => 256,
                Some(l) => {
                    if a.algorithm != "heuristic" {
                        log::warn!("Ignoring --row-limit as heuristic is not used as algorithm.");
                    } else if l > 65535 {
                        bail!("Row limit parameter must not exceed 65535.");
                    }
                    l
                }
            };
            let algorithm = match a.algorithm.as_str() {
                // a negative limit wraps around in WhatsHap's unsigned conversion
                "heuristic" => Algorithm::Heuristic { row_limit: row_limit as u32 },
                "hapchat" => Algorithm::HapChat,
                _ => Algorithm::WhatsHap,
            };
            if algorithm == Algorithm::HapChat && a.ped.is_some() {
                bail!("The hapchat algorithm cannot do pedigree phasing");
            }
            if algorithm == Algorithm::HapChat && a.read_list_filename.is_some() {
                log::warn!(
                    "On which haplotype a read occurs in the inferred solution is not yet implemented in hapchat, and so the corresponding column in the read list file contains no information about this"
                );
            }
            let tag = if a.tag == "HP" { PhaseTag::HP } else { PhaseTag::PS };
            let command_line = format!("(lrpipe-phase {}) {}", env!("CARGO_PKG_VERSION"), args.join(" "));
            run_phase(PhaseOptions {
                variant_file: a.variant_file,
                phase_input_files: a.phase_input_files,
                reference: a.reference,
                no_reference: a.no_reference,
                output: a.output,
                samples: a.samples,
                chromosomes: a.chromosomes,
                excluded_chromosomes: a.excluded_chromosomes,
                ignore_read_groups: a.ignore_read_groups,
                only_snvs: a.only_snvs,
                mapping_quality: a.mapping_quality,
                max_coverage: a.max_coverage,
                tag,
                read_list_filename: a.read_list_filename,
                use_supplementary: a.use_supplementary,
                supplementary_distance_threshold: a.supplementary_distance_threshold,
                command_line: Some(command_line),
                threads: a.threads,
                dp_memory_mb: a.dp_memory_mb,
                max_memory_mb: a.max_memory_mb,
                read_merging: a.read_merging,
                algorithm,
                ped: a.ped,
                use_ped_samples: a.use_ped_samples,
                recombrate: a.recombrate,
                genmap: a.genmap,
                genetic_haplotyping: !a.no_genetic_haplotyping,
                recombination_list_filename: a.recombination_list_filename,
                distrust_genotypes: a.distrust_genotypes,
                include_homozygous: a.include_homozygous,
                default_gq: a.default_gq,
                gl_regularizer: a.gl_regularizer,
                gtchange_list_filename: a.gtchange_list_filename,
            })
        }
    }
}

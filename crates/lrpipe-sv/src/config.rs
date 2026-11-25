//! Command line and derived settings, with Sniffles2 2.8.1's defaults (`config.py`).

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, bail};
use clap::Parser;

fn tobool(v: &str) -> Result<bool, String> {
    match v.trim().to_ascii_lowercase().as_str() {
        "true" | "1" => Ok(true),
        "false" | "0" => Ok(false),
        _ => Err("Boolean value (True | False) required for argument".into()),
    }
}

#[derive(Debug, Parser)]
#[command(version, about = "Fast structural variant calling from long reads (Sniffles2 2.8.1 port)")]
pub struct Args {
    /// Sorted and indexed BAM or CRAM input (CRAM needs --reference); or, for multi-sample
    /// calling, one or more .snf files or a .tsv list of .snf files (with optional sample IDs)
    #[arg(short, long, num_args = 1.., required = true)]
    pub input: Vec<String>,
    /// VCF output filename
    #[arg(short, long)]
    pub vcf: Option<String>,
    /// SNF output filename: SV candidates for later multi-sample calling
    #[arg(long, value_name = "OUT.snf")]
    pub snf: Option<String>,
    /// Reference sequence (FASTA with .fai): DEL sequences and the base before each SV
    #[arg(long)]
    pub reference: Option<String>,
    /// Accepted for compatibility; phase information is always determined (as in Sniffles2)
    #[arg(long)]
    pub phase: bool,
    /// Number of parallel threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,
    /// Only process these contigs (may be given several times)
    #[arg(short, long)]
    pub contig: Vec<String>,
    /// Call SVs on contigs shorter than 1 Mbp too
    #[arg(long)]
    pub all_contigs: bool,
    /// Only process the regions in this BED file
    #[arg(long, value_name = "REGIONS.bed")]
    pub regions: Option<String>,
    /// Only process this region (may be given several times)
    #[arg(long, value_name = "contig:start-end")]
    pub region: Vec<String>,

    /// Minimum number of supporting reads ("auto" for coverage-based)
    #[arg(long, default_value = "3")]
    pub minsupport: String,
    /// Coverage based minimum support multiplier (auto minsupport only)
    #[arg(long)]
    pub minsupport_auto_mult: Option<f64>,
    /// Minimum SV length; a ~ prefix allows slightly smaller SVs if strongly supported
    #[arg(long, default_value = "~50")]
    pub minsvlen: String,
    #[arg(long, default_value_t = 0.9)]
    pub minsvlen_screen_ratio: f64,
    /// Alignments with mapping quality below this are ignored
    #[arg(long, default_value_t = 20)]
    pub mapq: u8,
    /// Output all SV candidates, including those failing QC
    #[arg(long, alias = "qc-output-all")]
    pub no_qc: bool,
    /// Only output calls with FILTER PASS
    #[arg(long)]
    pub pass_only: bool,
    #[arg(long, default_value = "True", value_parser = tobool, action = clap::ArgAction::Set)]
    pub qc_stdev: bool,
    #[arg(long, default_value_t = 500.0)]
    pub qc_stdev_abs_max: f64,
    #[arg(long, default_value = "False", value_parser = tobool, action = clap::ArgAction::Set)]
    pub qc_strand: bool,
    #[arg(long, default_value_t = 1)]
    pub qc_coverage: i64,
    #[arg(long, default_value_t = 2500)]
    pub long_ins_length: i64,
    #[arg(long, default_value_t = 50000)]
    pub long_del_length: i64,
    #[arg(long, default_value_t = 10000)]
    pub long_inv_length: i64,
    #[arg(long, default_value_t = 0.66)]
    pub long_del_coverage: f64,
    #[arg(long, default_value_t = 50000)]
    pub long_dup_length: i64,
    #[arg(long, default_value_t = 1.33)]
    pub long_dup_coverage: f64,
    #[arg(long, default_value = "True", value_parser = tobool, action = clap::ArgAction::Set)]
    pub qc_bnd_filter_strand: bool,
    #[arg(long, default_value_t = 1000)]
    pub bnd_min_split_length: i64,
    #[arg(long, default_value_t = 0.1)]
    pub max_splits_kb: f64,
    #[arg(long, default_value_t = 3)]
    pub max_splits_base: i64,
    #[arg(long, default_value_t = 1000)]
    pub min_alignment_length: i64,
    #[arg(long, default_value_t = 0.1)]
    pub phase_conflict_threshold: f64,
    #[arg(long, default_value = "True", value_parser = tobool, action = clap::ArgAction::Set)]
    pub detect_large_ins: bool,
    #[arg(long, default_value_t = 0.5)]
    pub max_unknown_pct: f64,

    #[arg(long, default_value_t = 100)]
    pub cluster_binsize: i64,
    #[arg(long, default_value_t = 2.5)]
    pub cluster_r: f64,
    #[arg(long, default_value_t = 1.5)]
    pub cluster_repeat_h: f64,
    #[arg(long, default_value_t = 1000.0)]
    pub cluster_repeat_h_max: f64,
    /// Maximum distance between two indels of a read that are joined into one lead
    /// [default: 150, or 50 with --accuracy]
    #[arg(long, allow_hyphen_values = true)]
    pub cluster_merge_pos: Option<i64>,
    #[arg(long)]
    pub cluster_merge_len: Option<f64>,
    #[arg(long, default_value_t = 1000)]
    pub cluster_merge_bnd: i64,

    #[arg(long, default_value_t = 2)]
    pub genotype_ploidy: i64,
    #[arg(long, default_value_t = 0.05)]
    pub genotype_error: f64,
    /// Sample ID in the VCF output
    #[arg(long)]
    pub sample_id: Option<String>,
    /// Determine genotypes of the SVs in this VCF instead of calling SVs
    #[arg(long, value_name = "IN.vcf")]
    pub genotype_vcf: Option<String>,
    /// Maximum distance (scaled by the square root of the SV length) for matching SVs
    #[arg(long, default_value_t = 250)]
    pub combine_match: i64,
    /// Maximum distance for matching SVs
    #[arg(long, default_value_t = 1000)]
    pub combine_match_max: i64,
    /// Minimum fraction of samples in which an SV must have passed QC to be reported (multi-sample)
    #[arg(long, default_value_t = 0.0)]
    pub combine_high_confidence: f64,
    /// Minimum fraction of samples in which an SV must be present (failed QC) to be reported
    #[arg(long, default_value_t = 0.2)]
    pub combine_low_confidence: f64,
    /// Minimum number of samples in which an SV must be present (failed QC) to be reported
    #[arg(long, default_value_t = 2)]
    pub combine_low_confidence_abs: usize,
    /// Minimum coverage for an absent sample's genotype to be 0/0 rather than ./.
    #[arg(long, default_value_t = 5)]
    pub combine_null_min_coverage: i64,
    /// Do not combine SVs of the same sample
    #[arg(long)]
    pub combine_separate_intra: bool,
    /// Include low-confidence / mosaic SVs in multi-sample calling
    #[arg(long)]
    pub combine_output_filtered: bool,
    /// Override low-quality genotypes with the best one (e.g. for tumor-normal pairs)
    #[arg(long)]
    pub combine_pair_relabel: bool,
    /// Genotype quality below which --combine-pair-relabel relabels a genotype
    #[arg(long, default_value_t = 20)]
    pub combine_pair_relabel_threshold: i64,
    /// Open each .snf file again for every block read instead of keeping all open, for more
    /// inputs than the open-file limit (ulimit -n) allows; done anyway when they would not fit
    #[arg(long)]
    pub combine_close_handles: bool,
    /// Minimum sequence similarity (1 - edit distance / length) for merging insertions; 0 disables
    #[arg(long, default_value_t = 0.7)]
    pub combine_pctseq: f64,
    /// Accepted for compatibility (output is always sorted)
    #[arg(long, default_value_t = 20)]
    pub combine_max_inmemory_results: usize,
    /// Minimum support for SVs to be considered for multi-sample calling
    #[arg(long, default_value_t = 3)]
    pub combine_support_threshold: i64,

    /// Output names of supporting reads
    #[arg(long)]
    pub output_rnames: bool,
    /// Do not compute consensus sequences for insertions
    #[arg(long)]
    pub no_consensus: bool,
    /// Do not sort the output VCF by position
    #[arg(long)]
    pub no_sort: bool,
    #[arg(long, default_value_t = 50000)]
    pub max_del_seq_len: i64,
    /// Symbolic ALTs (<INS>, <DEL>, ...) instead of sequences
    #[arg(long)]
    pub symbolic: bool,
    /// Overwrite existing output files
    #[arg(long)]
    pub allow_overwrite: bool,

    /// Mosaic SV calling mode
    #[arg(long)]
    pub mosaic: bool,
    #[arg(long, default_value_t = 0.218)]
    pub mosaic_af_max: f64,
    #[arg(long, default_value_t = 0.05)]
    pub mosaic_af_min: f64,
    #[arg(long, default_value_t = 500)]
    pub mosaic_qc_invdup_min_length: i64,
    #[arg(long, default_value_t = 1.66)]
    pub mosaic_qc_nm_mult: f64,
    #[arg(long, default_value_t = -1.0)]
    pub mosaic_qc_coverage_max_change_frac: f64,
    #[arg(long, default_value = "True", value_parser = tobool, action = clap::ArgAction::Set)]
    pub mosaic_qc_strand: bool,
    #[arg(long)]
    pub mosaic_include_germline: bool,
    #[arg(long, default_value_t = 50000)]
    pub max_svlen_mosaic: i64,
    /// Not a command-line option in Sniffles2 (class attribute)
    #[arg(skip = 3usize)]
    pub mosaic_min_reads: usize,
    #[arg(skip = 10usize)]
    pub mosaic_use_strand_thresholds: usize,

    /// Tandem repeat annotations (BED) for improved calling in repeats
    #[arg(long)]
    pub tandem_repeats: Option<String>,
    #[arg(long, default_value_t = 50000)]
    pub dev_seq_cache_maxlen: i64,
    #[arg(long, default_value_t = 10)]
    pub consensus_max_reads_bin: usize,
    #[arg(long)]
    pub repeat: bool,
    #[arg(long)]
    pub qc_nm: bool,
    #[arg(long, default_value_t = 1.66)]
    pub qc_nm_mult: f64,
    #[arg(long, default_value_t = -1.0)]
    pub qc_coverage_max_change_frac: f64,
    #[arg(long, default_value_t = 5)]
    pub coverage_updown_bins: i64,
    #[arg(long, default_value_t = 20)]
    pub cluster_resplit_binsize: i64,
    #[arg(long, alias = "excl-flags", short = 'F')]
    pub exclude_flags: Option<u16>,
    #[arg(long)]
    pub dev_no_qc: bool,
    #[arg(long)]
    pub dev_no_resplit: bool,
    #[arg(long)]
    pub dev_no_resplit_repeat: bool,
    #[arg(long)]
    pub dev_keep_lowqual_splits: bool,
    #[arg(long, default_value_t = -1)]
    pub dev_min_leads_cluster: i64,
    #[arg(long, default_value_t = 1.0 / 6.0)]
    pub dev_min_dup_vaf: f64,
    #[arg(long, default_value_t = 200000)]
    pub dev_longer_del: i64,
    #[arg(long, default_value_t = 200000)]
    pub dev_longer_dup: i64,
    #[arg(long, default_value_t = 5)]
    pub dev_minreads_extra: i64,
    #[arg(long, default_value_t = 10000)]
    pub dev_maxsvlen_extra: i64,
    #[arg(long, default_value_t = 0.80)]
    pub dev_inline_sa_support_max: f64,
    #[arg(long, default_value_t = 500)]
    pub dev_min_close_edge_dist: i64,
    #[arg(long, default_value_t = 0.75)]
    pub dev_min_read_close_edge_prop: f64,
    #[arg(long)]
    pub dev_emit_sv_lengths: bool,

    /// lrpipe-sv improvement (changes output): best-accuracy settings, --split-alleles,
    /// --read-support-genotypes and --cluster-merge-pos 50 (an explicit --cluster-merge-pos
    /// still applies)
    #[arg(long, help_heading = "lrpipe-sv improvements")]
    pub accuracy: bool,
    /// lrpipe-sv improvement (changes output): genotype insertions of --long-ins-length or more
    /// from the reads that support them. Sniffles2 scales their support up with their length,
    /// for the reads clipped inside them, but already counts those reads (SUPPORT_LONG), so its
    /// heterozygous long insertions come out 1/1. Only GT, GQ, DR and DV change: VAF and every
    /// filter keep the scaled count, so the same calls pass
    #[arg(long, help_heading = "lrpipe-sv improvements")]
    pub read_support_genotypes: bool,
    /// lrpipe-sv improvement (changes output): report two alleles at a locus when the
    /// supporting reads show two clearly separated SV lengths
    #[arg(long, help_heading = "lrpipe-sv improvements")]
    pub split_alleles: bool,
    /// --split-alleles: minimum share of a cluster's reads in each allele
    #[arg(long, default_value_t = 0.2, help_heading = "lrpipe-sv improvements")]
    pub split_min_frac: f64,
    /// --split-alleles: minimum reads in each allele
    #[arg(long, default_value_t = 3, help_heading = "lrpipe-sv improvements")]
    pub split_min_support: u32,
    /// --split-alleles: minimum empty gap between the alleles' lengths, in within-allele standard deviations
    #[arg(long, default_value_t = 1.0, help_heading = "lrpipe-sv improvements")]
    pub split_min_gap_sd: f64,
    /// --split-alleles: length noise floor, relative to the longer allele
    #[arg(long, default_value_t = 0.01, help_heading = "lrpipe-sv improvements")]
    pub split_noise: f64,

    /// Also write read depth, counted in the same pass over the reads, with output identical to
    /// `mosdepth -n` (and lrpipe-depth): PREFIX.mosdepth.summary.txt and .global.dist.txt, and with
    /// --depth-by PREFIX.regions.bed.gz and .region.dist.txt. Covers all contigs (only --contig
    /// restricts it)
    #[arg(long, value_name = "PREFIX", help_heading = "Depth output")]
    pub depth: Option<String>,
    /// --depth: windows of this many bases, or the regions of a BED file (mosdepth -b)
    #[arg(long, value_name = "bed|window", help_heading = "Depth output")]
    pub depth_by: Option<String>,
    /// --depth: ignore reads with a lower mapping quality (mosdepth -Q)
    #[arg(long, default_value_t = 0, help_heading = "Depth output")]
    pub depth_mapq: u8,
    /// --depth: ignore reads with any of these flag bits (mosdepth -F)
    #[arg(long, default_value_t = 1796, help_heading = "Depth output")]
    pub depth_flag: u16,
    /// --depth: count each read's whole span, not only its aligned bases (mosdepth -x)
    #[arg(long, help_heading = "Depth output")]
    pub depth_fast: bool,
    /// --depth: report each region's median depth instead of its mean (mosdepth -m)
    #[arg(long, help_heading = "Depth output")]
    pub depth_median: bool,
    /// --depth: per region, count the bases with at least these depths (mosdepth -T)
    #[arg(long, value_delimiter = ',', help_heading = "Depth output")]
    pub depth_thresholds: Vec<i32>,
    /// --depth: count only the reads of --regions (sorted, apart), in the pass, for the
    /// --depth-by BED's regions, which must lie in them: those regions' rows are mosdepth's,
    /// but the summary has only the regions' rows (CONTIG_region, total_region) and no global
    /// distribution is written, as the contigs' depth is not known (lrpipe's own option)
    #[arg(long, requires_all = ["depth", "regions"], help_heading = "Depth output")]
    pub depth_in_regions: bool,

    /// Also genotype tandem repeats (lrpipe-tr) from the same pass over the reads: the loci of this
    /// catalog (BED with TRGT-style `ID=..;MOTIFS=..`, or a motif list in column 4). Needs
    /// --reference and --repeats-vcf
    #[arg(long, value_name = "BED", help_heading = "Tandem repeat output")]
    pub repeats: Option<String>,
    /// --repeats: output VCF (.vcf.gz is bgzipped), identical to lrpipe-tr's on the same reads
    #[arg(long, value_name = "VCF", help_heading = "Tandem repeat output")]
    pub repeats_vcf: Option<String>,
    /// --repeats: sequencing platform, which sets lrpipe-tr's genotyping noise model
    #[arg(long, value_enum, default_value = "hifi", help_heading = "Tandem repeat output")]
    pub repeats_platform: RepeatPlatform,
    /// --repeats: sex chromosomes, XX, or XY for one copy of chrX (outside the PARs) and chrY
    #[arg(long, value_enum, default_value = "xx", help_heading = "Tandem repeat output")]
    pub karyotype: KaryotypeArg,
    /// --repeats: only the catalog's loci that overlap these regions (BED), each genotyped from
    /// all its reads, as `lrpipe tr --regions` (lrpipe's own option)
    #[arg(long, value_name = "BED", help_heading = "Tandem repeat output")]
    pub repeats_regions: Option<String>,

    /// Also call copy number variants (lrpipe-cnv) into this directory, from the --depth windows
    /// of the same pass, which must be mosdepth's `-x -Q 20 -b 1000` (--depth-fast --depth-mapq
    /// 20 --depth-by 1000). Needs --reference
    #[arg(long, value_name = "DIR", help_heading = "Copy number output")]
    pub cnv: Option<String>,
    /// --cnv: SNV VCF of the sample (lrpipe-snv, Clair3), for copy-neutral depth and LoH calls
    #[arg(long, value_name = "VCF", help_heading = "Copy number output")]
    pub cnv_snv: Option<String>,
    /// --cnv: blacklist BED (Spectre's data/grch38_blacklist_spectre.bed)
    #[arg(long, value_name = "BED", help_heading = "Copy number output")]
    pub cnv_blacklist: Option<String>,
    /// --cnv: the reference's .mdr (N runs); made in the --cnv directory if not given
    #[arg(long, value_name = "MDR", help_heading = "Copy number output")]
    pub cnv_metadata: Option<String>,
    /// --cnv: lrpipe-cnv's improved method (sex chromosomes, GC correction, amplifications,
    /// copy numbers from the final calls, breakpoints from this run's SVs, valid VCF 4.2), or
    /// Spectre 0.2.1's exactly
    #[arg(long, value_enum, default_value = "improved", help_heading = "Copy number output")]
    pub cnv_mode: CnvMode,

    /// Also write CpG methylation (lrpipe-methyl), counted in the same pass over the reads, as
    /// bedMethyl (.gz for BGZF): identical to `lrpipe methyl pileup INPUT BEDMETHYL --ref
    /// REFERENCE --modified-bases 5mC 5hmC --cpg --combine-strands` (modkit's CpG recipe),
    /// with its pass thresholds. Covers all contigs (only --contig restricts it). Needs
    /// --reference
    #[arg(long, value_name = "BEDMETHYL", help_heading = "Methylation output")]
    pub bedmethyl: Option<String>,
    /// --bedmethyl: only the positions in these regions (BED), with the pass thresholds
    /// estimated from their reads: identical to `lrpipe methyl pileup` with `--include-bed`, as
    /// modkit's (lrpipe's own option)
    #[arg(long, value_name = "BED", requires = "bedmethyl", help_heading = "Methylation output")]
    pub bedmethyl_regions: Option<String>,

    /// Also write the reads' statistics, counted in the same pass, as JSON: the primary
    /// alignments of mapped reads (their count, bases, length N50, mean length and mean base
    /// quality) and the secondary, supplementary and placed unmapped records. Covers all contigs
    /// (only --contig and --regions restrict it). From CRAM, without the mean quality: the pass
    /// does not decode CRAM's qualities, which would make it two thirds slower
    #[arg(long, value_name = "JSON", help_heading = "Read statistics output")]
    pub read_stats: Option<String>,

    /// The reads' haplotypes from a haplotag list (WhatsHap's `haplotag --output-haplotag-list`)
    /// in place of their HP and PS tags, for the SVs' phases and the repeats, as `haplotag`
    /// would have tagged them: a listed read's primary alignment gets its haplotype and phase
    /// set, every other alignment none
    #[arg(long, value_name = "LIST", help_heading = "Phasing")]
    pub haplotags: Option<String>,
    /// Phase the SVs later, from a haplotag list (`lrpipe sv phase`): call them as if the reads
    /// had no HP and PS tags, and write what phasing them needs (their reads' names) to FILE.
    /// The repeats still read the tags
    #[arg(long, value_name = "FILE", conflicts_with_all = ["haplotags", "genotype_vcf"], help_heading = "Phasing")]
    pub phasing_evidence: Option<String>,
}

/// Where the SV calls take the reads' haplotypes from.
#[derive(Debug, Default)]
pub enum Haplotypes {
    /// Their alignments' HP and PS tags.
    #[default]
    Tags,
    /// A haplotag list (`--haplotags`).
    List(lrpipe_reads::haplotags::Haplotags),
    /// Nowhere: the SVs are phased later (`--phasing-evidence`).
    Later,
}

/// `--repeats-platform`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum RepeatPlatform {
    Hifi,
    Ont,
}

/// `--karyotype`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum KaryotypeArg {
    Xx,
    Xy,
}

/// `--cnv-mode`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum CnvMode {
    Improved,
    Spectre,
}

/// `--repeats`: the catalog and how to genotype it.
#[derive(Debug)]
pub struct RepeatConfig {
    pub output: String,
    /// The reference, for flank alignment.
    pub genome: std::path::PathBuf,
    /// Loci by contig, each list sorted by start; contigs in the reference's order.
    pub loci: Vec<(String, Vec<lrpipe_tr::catalog::Locus>)>,
    pub opts: lrpipe_tr::run::RunOptions,
}

impl RepeatConfig {
    fn from_args(args: &Args) -> Result<Option<Self>> {
        let Some(catalog) = args.repeats.as_deref() else {
            if args.repeats_vcf.is_some() || args.repeats_regions.is_some() {
                bail!("--repeats-vcf and --repeats-regions need --repeats");
            }
            return Ok(None);
        };
        let Some(output) = args.repeats_vcf.clone() else { bail!("--repeats needs --repeats-vcf") };
        let Some(reference) = args.reference.as_deref() else { bail!("--repeats needs --reference") };
        let genome = std::path::PathBuf::from(reference);
        let fai_path = genome.with_added_extension("fai");
        let fai = pasta::fasta::fai::fs::read(&fai_path).with_context(|| format!("reading {}", fai_path.display()))?;
        let contigs: Vec<(String, u64)> = fai.as_ref().iter().map(|r| (r.name().to_string(), r.length())).collect();
        let mut all =
            lrpipe_tr::catalog::read(std::path::Path::new(catalog)).with_context(|| format!("reading {catalog}"))?;
        lrpipe_tr::run::sort_loci(&mut all, &contigs);
        if let Some(path) = args.repeats_regions.as_deref() {
            let regions = lrpipe_reads::regions::Regions::read(Path::new(path))?;
            lrpipe_tr::catalog::retain_in(&mut all, &regions);
        }
        let mut loci: Vec<(String, Vec<lrpipe_tr::catalog::Locus>)> = Vec::new();
        for l in all {
            match loci.last_mut() {
                Some((c, v)) if *c == l.contig => v.push(l),
                _ => loci.push((l.contig.clone(), vec![l])),
            }
        }
        let platform = match args.repeats_platform {
            RepeatPlatform::Hifi => lrpipe_tr::genotype::Platform::Hifi,
            RepeatPlatform::Ont => lrpipe_tr::genotype::Platform::Ont,
        };
        let karyotype = match args.karyotype {
            KaryotypeArg::Xx => lrpipe_tr::run::Karyotype::Xx,
            KaryotypeArg::Xy => lrpipe_tr::run::Karyotype::Xy,
        };
        let opts = lrpipe_tr::run::RunOptions {
            threads: args.threads.max(1),
            karyotype,
            genotype: lrpipe_tr::genotype::GenotypeOptions::for_platform(platform),
            x_pars: lrpipe_tr::run::x_pars(&genome, karyotype)?,
            ..lrpipe_tr::run::RunOptions::default()
        };
        Ok(Some(RepeatConfig { output, genome, loci, opts }))
    }

    /// The loci of `contig`, sorted by start.
    pub fn contig_loci(&self, contig: &str) -> &[lrpipe_tr::catalog::Locus] {
        self.loci.iter().find(|(c, _)| c == contig).map_or(&[], |(_, v)| v.as_slice())
    }
}

/// `--cnv`: where to write and how to call.
#[derive(Debug)]
pub struct CnvConfig {
    pub dir: String,
    pub mode: CnvMode,
    /// The contigs to call with `--contig`: those of them lrpipe-cnv calls (longer than its
    /// minimum contig length); `None` for every such contig of the reference.
    pub only_chr: Option<Vec<String>>,
}

impl CnvConfig {
    fn from_args(args: &Args, depth: Option<&DepthConfig>) -> Result<Option<Self>> {
        let Some(dir) = args.cnv.clone() else {
            if args.cnv_snv.is_some() || args.cnv_blacklist.is_some() || args.cnv_metadata.is_some() {
                bail!("--cnv-snv, --cnv-blacklist and --cnv-metadata need --cnv");
            }
            return Ok(None);
        };
        if args.reference.is_none() {
            bail!("--cnv needs --reference");
        }
        let spectre_windows = depth.is_some_and(|d| {
            d.opts.mode == lrpipe_depth::Mode::Fast
                && d.opts.filter.min_mapq == 20
                && d.opts.filter.exclude_flags == 1796
                && matches!(d.opts.report.by, Some(lrpipe_depth::By::Windows(1000)))
                && !d.opts.report.use_median
        });
        if !spectre_windows {
            bail!(
                "--cnv calls from the --depth windows, which must be Spectre's: add --depth PREFIX --depth-by 1000 \
                 --depth-fast --depth-mapq 20 (and no --depth-median or --depth-flag)"
            );
        }
        // The depth windows cover only the --contig contigs, so only they are called; checked
        // here, as lrpipe-cnv would fail on the others' missing coverage after the whole pass.
        let only_chr = if args.contig.is_empty() {
            None
        } else {
            let reference = args.reference.as_deref().unwrap_or_default();
            let min_len = lrpipe_cnv::pipeline::DEFAULT_MIN_CHR_LEN;
            let genome = lrpipe_cnv::genome::Genome::from_fasta_index(reference.as_ref(), min_len, 2, &[])
                .with_context(|| format!("--cnv: reading the index of {reference}"))?;
            // In the reference's order, each once, as lrpipe-cnv takes the contigs without
            // --contig: its thresholds are drawn from them in order.
            let called: Vec<String> =
                genome.contigs.iter().filter(|c| args.contig.contains(&c.name)).map(|c| c.name.clone()).collect();
            if called.is_empty() {
                bail!(
                    "--cnv calls reference contigs longer than {min_len} bp, and no --contig is one ({})",
                    args.contig.join(", ")
                );
            }
            Some(called)
        };
        Ok(Some(CnvConfig { dir, mode: args.cnv_mode, only_chr }))
    }
}

/// `--depth`: how to count read depth and what to write.
#[derive(Debug)]
pub struct DepthConfig {
    pub prefix: String,
    /// Also used for contigs the SV pass does not read whole (see `main::write_depth`).
    pub opts: lrpipe_depth::run::RunOptions,
}

impl DepthConfig {
    fn from_args(args: &Args) -> Result<Option<Self>> {
        let Some(prefix) = args.depth.clone() else { return Ok(None) };
        if !args.depth_thresholds.is_empty() && args.depth_by.is_none() {
            bail!("--depth-thresholds can only be used with --depth-by");
        }
        let mut thresholds = args.depth_thresholds.clone();
        thresholds.sort_unstable();
        // As mosdepth (and lrpipe-depth), MOSDEPTH_PRECISION sets the decimal places.
        let precision = match std::env::var("MOSDEPTH_PRECISION") {
            Ok(p) => p.parse().context("invalid MOSDEPTH_PRECISION")?,
            Err(_) => 2,
        };
        let opts = lrpipe_depth::run::RunOptions {
            mode: if args.depth_fast { lrpipe_depth::Mode::Fast } else { lrpipe_depth::Mode::Cigar },
            filter: lrpipe_depth::Filter {
                min_mapq: args.depth_mapq,
                exclude_flags: args.depth_flag,
                ..lrpipe_depth::Filter::default()
            },
            report: lrpipe_depth::ReportOptions {
                by: args.depth_by.as_deref().map(lrpipe_depth::By::parse).transpose()?,
                use_median: args.depth_median,
                thresholds,
                precision,
                regions_only: args.depth_in_regions,
            },
            chrom: None,
            threads: args.threads,
        };
        Ok(Some(DepthConfig { prefix, opts }))
    }
}

/// Regions from a BED file, or else from `--region` strings, grouped by contig in input order.
/// Malformed entries are skipped with a warning, as in Sniffles2.
fn read_regions(bed: Option<&str>, strings: &[String]) -> Result<HashMap<String, Vec<(i64, i64)>>> {
    let mut by_contig: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    let parse = |s: &str| s.trim().parse::<i64>().ok();
    let mut bed_lines = 0usize;
    if let Some(path) = bed {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            bed_lines += 1;
            let f: Vec<&str> = line.split('\t').collect();
            match (f.first(), f.get(1).and_then(|s| parse(s)), f.get(2).and_then(|s| parse(s))) {
                (Some(c), Some(s), Some(e)) => by_contig.entry((*c).to_string()).or_default().push((s, e)),
                _ => eprintln!("Invalid region line: {line}"),
            }
        }
    } else {
        for r in strings {
            let parsed = r.split_once(':').and_then(|(c, se)| {
                let (s, e) = se.split_once('-')?;
                Some((c, parse(s)?, parse(e)?))
            });
            match parsed {
                Some((c, s, e)) if !c.contains(':') => by_contig.entry(c.to_string()).or_default().push((s, e)),
                _ => eprintln!("Invalid region string: {r}"),
            }
        }
    }
    // Sniffles2 then calls the whole genome: regions asked for but none valid is a mistake.
    if by_contig.is_empty() {
        if let Some(path) = bed {
            if bed_lines > 0 {
                bail!(
                    "--regions {path} has no valid region: each line needs a contig, start and end, separated by tabs"
                );
            }
        } else if !strings.is_empty() {
            bail!(
                "no valid --region among {strings:?}: give each as CONTIG:START-END, such as chr20:1000-2000, without commas"
            );
        }
    }
    Ok(by_contig)
}

/// Checks `--depth-in-regions`: `regions` (by contig) sorted and apart, which counting each read
/// once needs, and every region of the `--depth-by` BED inside one of them, so that its depth is
/// of all the reads that reach it.
fn depth_regions_within(regions: &HashMap<String, Vec<(i64, i64)>>, by: Option<&str>) -> Result<()> {
    for (contig, spans) in regions {
        if spans.windows(2).any(|w| w[1].0 < w[0].1) {
            bail!("--depth-in-regions needs --regions sorted and apart (bedtools merge), as {contig}'s are not");
        }
    }
    let Some(lrpipe_depth::By::Bed(bed)) = by.map(lrpipe_depth::By::parse).transpose()? else {
        bail!("--depth-in-regions needs --depth-by with a BED of regions inside --regions");
    };
    let names: Vec<String> = regions.keys().cloned().collect();
    if let Some((contig, _)) = bed.missing(&names).next() {
        bail!("--depth-in-regions: the --depth-by BED has {contig}, which --regions does not");
    }
    for (contig, spans) in regions {
        for r in bed.get(contig) {
            let (start, stop) = (i64::from(r.start), i64::from(r.stop));
            if !spans.iter().any(|&(s, e)| s <= start && stop <= e) {
                bail!(
                    "--depth-in-regions: {contig}:{start}-{stop} of the --depth-by BED is not inside one of --regions"
                );
            }
        }
    }
    Ok(())
}

/// Whether a VCF path asks for bgzipped output (`.gz` or `.bgz`).
pub fn is_bgzip_name(path: &str) -> bool {
    path.ends_with(".gz") || path.ends_with(".bgz")
}

/// How `--minsupport` is applied.
#[derive(Debug, Clone, Copy)]
pub enum MinSupport {
    Auto,
    Const(usize),
}

/// All settings used by the caller.
#[derive(Debug)]
pub struct Config {
    pub args: Args,
    pub command: String,
    pub minsvlen: i64,
    pub minsvlen_hard_cap: bool,
    pub minsvlen_screen: i64,
    pub minsupport: MinSupport,
    pub minsupport_const: usize,
    pub minsupport_auto_mult: f64,
    pub minsupport_auto_base: f64,
    pub minsupport_auto_regional_coverage_weight: f64,
    pub coverage_binsize: i64,
    pub consensus_min_reads: usize,
    pub consensus_kmer_len: usize,
    pub consensus_kmer_skip_base: usize,
    pub consensus_kmer_skip_seqlen_mult: f64,
    pub long_ins_rescale_base: f64,
    pub long_ins_rescale_mult: f64,
    pub dev_longer_dup: i64,
    pub dev_longer_del: i64,
    pub genotype_min_z_score: i64,
    pub precise: f64,
    pub tandem_repeat_region_pad: i64,
    pub qc_nm_measure: bool,
    pub cluster_merge_len: f64,
    pub dev_min_leads_cluster: usize,
    /// Effective `--cluster-merge-pos`.
    pub cluster_merge_pos: i64,
    /// `--regions` / `--region` by contig, in the order given ([start, end), 0-based).
    pub regions_by_contig: HashMap<String, Vec<(i64, i64)>>,
    pub sort: bool,
    /// `--depth` settings.
    pub depth: Option<DepthConfig>,
    /// `--repeats` settings.
    pub repeats: Option<RepeatConfig>,
    /// `--cnv` settings.
    pub cnv: Option<CnvConfig>,
    pub sample_id: String,
    /// Where the SV calls take the reads' haplotypes from.
    pub haplotypes: Haplotypes,
}

const DEFAULT_CLUSTER_MERGE_LEN: f64 = 0.22;
/// Sniffles2's `--cluster-merge-pos` default.
const DEFAULT_CLUSTER_MERGE_POS: i64 = 150;
/// `--cluster-merge-pos` under `--accuracy`: joining a read's indels across up to 150 bp merges
/// separate events in tandem repeats (0-50 bp score the same on HG002; 100 bp is worse).
const ACCURACY_CLUSTER_MERGE_POS: i64 = 50;
const DEFAULT_CLUSTER_MERGE_LEN_MOSAIC: f64 = 0.27;

impl Config {
    /// Whether the run calls SVs: with an SV output (`--vcf`, `--snf`) or `--cnv`, which uses
    /// them. Without, the pass only counts depth, repeats and methylation, and reads no SV leads.
    #[must_use]
    pub fn calls_svs(&self) -> bool {
        self.args.vcf.is_some() || self.args.snf.is_some() || self.cnv.is_some()
    }

    pub fn from_args(mut args: Args, command: String) -> Result<Self> {
        let (minsvlen_hard_cap, minsvlen) = match args.minsvlen.strip_prefix('~') {
            Some(v) => (false, v.parse()?),
            None => (true, args.minsvlen.parse()?),
        };
        let minsvlen_screen = (args.minsvlen_screen_ratio * minsvlen as f64) as i64;
        let (minsupport, minsupport_const) = if args.minsupport == "auto" {
            (MinSupport::Auto, 3)
        } else {
            let v: usize = args.minsupport.parse()?;
            (MinSupport::Const(v), v)
        };
        if args.accuracy {
            args.split_alleles = true;
            args.read_support_genotypes = true;
        }
        if !args.contig.is_empty() && args.regions.is_some() {
            bail!("Please provide either --contig or --regions, not both.");
        }
        let regions_by_contig = read_regions(args.regions.as_deref(), &args.region)?;
        if args.depth_in_regions {
            depth_regions_within(&regions_by_contig, args.depth_by.as_deref())?;
        }
        let cluster_merge_pos = args.cluster_merge_pos.unwrap_or(if args.accuracy {
            ACCURACY_CLUSTER_MERGE_POS
        } else {
            DEFAULT_CLUSTER_MERGE_POS
        });
        if args.dev_no_qc {
            args.no_qc = true;
        }
        if args.genotype_ploidy != 2 {
            bail!("Currently only --genotype-ploidy 2 is supported");
        }
        if args.mosaic_include_germline {
            args.mosaic = true;
        }
        let qc_nm_measure = args.qc_nm || args.mosaic;
        let mut cluster_merge_len = args.cluster_merge_len.unwrap_or(DEFAULT_CLUSTER_MERGE_LEN);
        if args.mosaic && cluster_merge_len == DEFAULT_CLUSTER_MERGE_LEN {
            cluster_merge_len = DEFAULT_CLUSTER_MERGE_LEN_MOSAIC;
        }
        let dev_min_leads_cluster = match args.dev_min_leads_cluster {
            -1 if args.no_qc => 1,
            -1 => 2,
            v => v.max(0) as usize,
        };
        let depth = DepthConfig::from_args(&args)?;
        if args.bedmethyl.is_some() && args.reference.is_none() {
            bail!("--bedmethyl needs --reference");
        }
        let cnv = CnvConfig::from_args(&args, depth.as_ref())?;
        let haplotypes = match (&args.haplotags, &args.phasing_evidence) {
            (Some(list), _) => Haplotypes::List(lrpipe_reads::haplotags::Haplotags::read(Path::new(list))?),
            (None, Some(_)) => Haplotypes::Later,
            (None, None) => Haplotypes::Tags,
        };
        Ok(Config {
            command,
            minsvlen,
            minsvlen_hard_cap,
            minsvlen_screen,
            minsupport,
            minsupport_const,
            minsupport_auto_mult: args.minsupport_auto_mult.unwrap_or(0.1),
            minsupport_auto_base: 1.5,
            minsupport_auto_regional_coverage_weight: 0.75,
            coverage_binsize: args.cluster_binsize,
            consensus_min_reads: 4,
            consensus_kmer_len: 6,
            consensus_kmer_skip_base: 3,
            consensus_kmer_skip_seqlen_mult: 1.0 / 500.0,
            long_ins_rescale_base: 1.66,
            long_ins_rescale_mult: 0.33,
            dev_longer_dup: (args.long_dup_length * 4).min(args.dev_longer_dup),
            dev_longer_del: (args.long_del_length * 4).min(args.dev_longer_del),
            genotype_min_z_score: 5,
            precise: 25.0,
            tandem_repeat_region_pad: 500,
            qc_nm_measure,
            cluster_merge_len,
            dev_min_leads_cluster,
            cluster_merge_pos,
            regions_by_contig,
            // Sniffles2 always sorts bgzipped output.
            sort: !args.no_sort || args.vcf.as_deref().is_some_and(is_bgzip_name),
            sample_id: args.sample_id.clone().unwrap_or_else(|| "SAMPLE".into()),
            depth,
            repeats: RepeatConfig::from_args(&args)?,
            cnv,
            haplotypes,
            args,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg(extra: &[&str]) -> Config {
        let mut argv = vec!["lrpipe-sv", "-i", "x.bam"];
        argv.extend_from_slice(extra);
        Config::from_args(Args::parse_from(argv), String::new()).expect("valid options")
    }

    #[test]
    fn defaults_match_sniffles2() {
        let c = cfg(&[]);
        assert!(!c.args.split_alleles);
        assert!(!c.args.read_support_genotypes);
        assert_eq!(c.cluster_merge_pos, 150);
    }

    #[test]
    fn accuracy_preset_turns_on_every_setting() {
        let c = cfg(&["--accuracy"]);
        assert!(c.args.split_alleles);
        assert!(c.args.read_support_genotypes);
        assert_eq!(c.cluster_merge_pos, 50);
    }

    #[test]
    fn cnv_calls_the_long_contig_options_once_in_reference_order() {
        // LSV-002: lrpipe-cnv draws its thresholds from the contigs in order, and a contig
        // given twice was called twice.
        let dir = std::env::temp_dir().join(format!("lrpipe-sv-cnv-contigs-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let reference = dir.join("ref.fa");
        std::fs::write(
            reference.with_added_extension("fai"),
            "c1\t2000000\t4\t60\t61\nshort\t5000\t1\t60\t61\nc2\t3000000\t1\t60\t61\n",
        )
        .expect("index");
        let mut argv = vec!["lrpipe-sv", "-i", "x.bam", "--reference", reference.to_str().expect("path")];
        argv.extend(["--depth", "d", "--depth-by", "1000", "--depth-fast", "--depth-mapq", "20", "--cnv", "out"]);
        argv.extend(["--contig", "c2", "--contig", "short", "--contig", "c1", "--contig", "c2"]);
        let c = cfg(&argv[3..]);
        let cnv = c.cnv.expect("--cnv");
        assert_eq!(cnv.only_chr, Some(vec!["c1".to_owned(), "c2".to_owned()]));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn explicit_cluster_merge_pos_wins() {
        assert_eq!(cfg(&["--accuracy", "--cluster-merge-pos", "20"]).cluster_merge_pos, 20);
        assert_eq!(cfg(&["--cluster-merge-pos", "-1"]).cluster_merge_pos, -1);
    }
}

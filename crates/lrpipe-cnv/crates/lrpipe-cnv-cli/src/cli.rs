use std::ffi::OsString;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};

use clap::{Parser, Subcommand};
use lrpipe_cnv::extensions::{Caller, CopyNumbers, Extensions, GcCorrection, Mappability, SexChromosomes};
use lrpipe_cnv::hmm::HmmOptions;
use lrpipe_cnv::normalize::OutlierRemoval;
use lrpipe_cnv::output::vcf;
use lrpipe_cnv::pipeline::{Options, call_cnvs, call_cnvs_with_depths, write_outputs};
use lrpipe_cnv::sex::Sex;
use lrpipe_cnv::thresholds::Clamping;

/// Long-read CNV calling, ported from Spectre 0.2.1 with identical output.
#[derive(Parser)]
#[command(name = "lrpipe-cnv", version, propagate_version = true, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Call CNVs, and with SNVs loss of heterozygosity (Spectre's CNVCaller)
    #[command(name = "CNVCaller", alias = "cnv-caller")]
    CnvCaller(CnvArgs),
}

#[derive(clap::Args)]
struct CnvArgs {
    /// mosdepth regions.bed.gz (from `mosdepth -x -b 1000 -Q 20`)
    #[arg(short = 'c', long, required_unless_present = "bam", conflicts_with = "bam")]
    coverage: Option<PathBuf>,
    /// BAM or CRAM to count depth from directly (as `mosdepth -x -b 1000 -Q 20`); CRAM is
    /// decoded with --cram-reference or -r
    #[arg(long)]
    bam: Option<PathBuf>,
    /// Reference the CRAM was written against, if not -r
    #[arg(long)]
    cram_reference: Option<PathBuf>,
    /// Threads for counting depth from --bam
    #[arg(short = 't', long, default_value_t = 8)]
    threads: usize,
    /// Sample ID (names the outputs and seeds the window sample)
    #[arg(short = 's', long)]
    sample_id: String,
    /// Output directory
    #[arg(short = 'd', long)]
    output_dir: PathBuf,
    /// Reference FASTA, with a .fai index
    #[arg(short = 'r', long)]
    reference: PathBuf,
    /// .mdr file of reference N runs; made from the reference if missing (default:
    /// `metadata.mdr` in the output directory)
    #[arg(short = 'm', long)]
    metadata: Option<PathBuf>,
    /// Single-sample SNV VCF (Clair3 or lrpipe-snv): copy-neutral depth from heterozygous
    /// SNVs, and LoH calls
    #[arg(short = 'v', long)]
    snv: Option<PathBuf>,
    /// Not supported yet
    #[arg(long)]
    snfj: Option<PathBuf>,
    /// Blacklist BED
    #[arg(short = 'l', long)]
    blacklist: Option<PathBuf>,
    /// Comma-separated contigs to call [default: all kept contigs]
    #[arg(short = 'o', long)]
    only_chr: Option<String>,
    /// Ploidy
    #[arg(short = 'p', long, default_value_t = 2)]
    ploidy: i64,
    /// Per-contig ploidy, `chr:N,...` (a coverage multiplier of ploidy²/N, as in Spectre)
    #[arg(long)]
    ploidy_chr: Option<String>,
    /// Shortest N run in the .mdr
    #[arg(short = 'n', long, default_value_t = 5)]
    n_size: i64,
    /// Shortest CNV reported (Spectre only accepts 10001..=100000)
    #[arg(long, default_value_t = 100_000)]
    min_cnv_len: i64,
    /// Contigs this long or shorter are ignored
    #[arg(long, default_value_t = 1_000_000)]
    min_chr_len: u64,
    /// Use this depth instead of the genome-wide median
    #[arg(long)]
    sample_coverage_overwrite: Option<f64>,
    /// Do not mask windows far above the genome-wide mean
    #[arg(long)]
    disable_max_coverage: bool,
    /// Do not clamp thresholds to ploidy ± 1
    #[arg(short = 'a', long)]
    cancer: bool,
    /// Share of a call's windows that may lie inside the strict thresholds
    #[arg(long, default_value_t = 0.5)]
    qc_noise_allow: f64,
    /// Not supported (always fails in Spectre 0.2.1)
    #[arg(long)]
    loh_only: bool,
    /// LoH calls need more than this many SNVs per kb ...
    #[arg(long, default_value_t = 10)]
    loh_min_snv_perkb: i64,
    /// ... or more than this many SNVs in total
    #[arg(long, default_value_t = 1000)]
    loh_min_snv_total: i64,
    /// LoH calls must span more than this many bases
    #[arg(long, default_value_t = 1_000_000.0)]
    loh_min_region_size: f64,
    /// Seed for the dip test's random draws, which Spectre leaves unseeded; it only matters
    /// when heterozygous SNV depths are not clearly unimodal
    #[arg(long, default_value_t = 0)]
    dip_seed: u32,
    /// Outlier cutoff, in standard deviations above the genome-wide mean
    #[arg(long, default_value_t = 5)]
    dev_max_std_outlier_rm: i64,
    /// Recorded in the .spc only
    #[arg(long, default_value_t = 0.10)]
    mosdepth_cov_genome_chr_diff: f64,
    /// Largest difference of median depth between candidates merged
    #[arg(long, default_value_t = 0.80)]
    cov_diff_threshold: f64,
    /// Recorded in the .spc only
    #[arg(long, default_value_t = 0.25)]
    dist_proportion: f64,

    /// Turn on every change below that needs no extra input (sex detection, GC correction,
    /// keeping amplifications, copy numbers from the final calls)
    #[arg(long, help_heading = "Changes to Spectre's method (off by default)")]
    improved: bool,
    /// Copies of chrX and chrY: from depth (detect), or the sample's sex
    #[arg(long, value_enum, help_heading = "Changes to Spectre's method (off by default)")]
    sex: Option<SexArg>,
    /// Correct depth for GC content
    #[arg(long, help_heading = "Changes to Spectre's method (off by default)")]
    gc_correction: bool,
    /// Bases over which a window's GC content is measured, centred on it (about a read's
    /// length)
    #[arg(long, default_value_t = lrpipe_cnv::extensions::DEFAULT_GC_SPAN, help_heading = "Changes to Spectre's method (off by default)")]
    gc_span: i64,
    /// Keep runs of windows above the outlier cutoff (amplifications); mask only short ones
    #[arg(long, help_heading = "Changes to Spectre's method (off by default)")]
    keep_amplifications: bool,
    /// Report CN, GT and CD from each call's final windows, not its first raw run
    #[arg(long, help_heading = "Changes to Spectre's method (off by default)")]
    refresh_copy_numbers: bool,
    /// Depth over all reads, any MAPQ (`mosdepth -x -b 1000 -Q 0`): windows whose coverage
    /// falls below --min-mapq-ratio of it are masked as ambiguously mapped
    #[arg(long, help_heading = "Changes to Spectre's method (off by default)")]
    all_mapq_coverage: Option<PathBuf>,
    /// See --all-mapq-coverage
    #[arg(long, default_value_t = 0.5, help_heading = "Changes to Spectre's method (off by default)")]
    min_mapq_ratio: f64,
    /// How windows become calls: Spectre's threshold runs, or a copy-number HMM (which also
    /// uses heterozygous SNV allele fractions with -v)
    #[arg(
        long,
        value_enum,
        default_value = "thresholds",
        help_heading = "Changes to Spectre's method (off by default)"
    )]
    caller: CallerArg,
    /// SV calls (lrpipe-sv, Sniffles2 VCF): CNV calls a deletion or duplication overlaps
    /// reciprocally by half take its base-pair breakpoints
    #[arg(long, help_heading = "Changes to Spectre's method (off by default)")]
    sv_breakpoints: Option<PathBuf>,
    /// VCF as Spectre writes it, or valid VCF 4.2 for strict readers (fields declared with the
    /// number of values they hold, END the event's last base, negative SVLEN for deletions)
    #[arg(long, value_enum, default_value = "spectre")]
    vcf_style: VcfStyleArg,
}

#[derive(Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
enum VcfStyleArg {
    Spectre,
    Standard,
}

#[derive(Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
enum CallerArg {
    Thresholds,
    Hmm,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum SexArg {
    Detect,
    Male,
    Female,
}

/// Spectre's short options of more than one letter, which clap cannot parse as short options.
const SPECTRE_SHORT: &[(&str, &str)] = &[
    ("-sj", "--snfj"),
    ("-pc", "--ploidy-chr"),
    ("-mcl", "--min-cnv-len"),
    ("-mchl", "--min-chr-len"),
    ("-sco", "--sample-coverage-overwrite"),
    ("-dmc", "--disable-max-coverage"),
    ("-qcna", "--qc-noise-allow"),
    ("-lohonly", "--loh-only"),
    ("-lohkb", "--loh-min-snv-perkb"),
    ("-lohsnv", "--loh-min-snv-total"),
    ("-lohsize", "--loh-min-region-size"),
    ("-01", "--dev-max-std-outlier-rm"),
    ("-02", "--mosdepth-cov-genome-chr-diff"),
    ("-05", "--cov-diff-threshold"),
    ("-06", "--dist-proportion"),
];

/// Runs the command line `args` (the program name first), as the `lrpipe-cnv` binary does.
///
/// # Errors
///
/// If the arguments are invalid or the run fails.
pub fn main(args: Vec<OsString>) -> Result<()> {
    let args = args.into_iter().map(|a| a.to_string_lossy().into_owned()).map(|a| {
        match SPECTRE_SHORT.iter().find(|(s, _)| *s == a) {
            Some((_, long)) => (*long).to_owned(),
            None => a,
        }
    });
    let Command::CnvCaller(a) = Cli::parse_from(args).command;
    if a.snfj.is_some() {
        bail!("--snfj is not ported yet");
    }
    if a.loh_only {
        bail!("--loh-only always fails in Spectre 0.2.1 and is not supported");
    }
    let coverage = a.coverage.clone().or_else(|| a.bam.clone()).context("-c or --bam is required")?;
    let mut opts = Options::new(coverage, a.sample_id, a.output_dir, a.reference);
    opts.metadata = a.metadata;
    opts.blacklist = a.blacklist;
    opts.snv = a.snv;
    opts.dip_seed = a.dip_seed;
    opts.loh.min_snv_perkb = a.loh_min_snv_perkb;
    opts.loh.min_snv_total = a.loh_min_snv_total;
    opts.loh.min_region_size = a.loh_min_region_size;
    opts.only_chr = a.only_chr.map(|s| s.split(',').map(str::to_owned).collect());
    opts.ploidy = a.ploidy;
    opts.ploidy_chr = match a.ploidy_chr {
        Some(list) => list
            .split(',')
            .map(|item| {
                let (chr, p) = item.split_once(':').with_context(|| format!("--ploidy-chr: {item} is not chr:N"))?;
                let p: i64 = p.parse().with_context(|| format!("--ploidy-chr: {item} is not chr:N"))?;
                Ok((chr.to_owned(), p))
            })
            .collect::<Result<_>>()?,
        None => Vec::new(),
    };
    opts.n_size = a.n_size;
    opts.min_chr_len = a.min_chr_len;
    opts.clamping = if a.cancer { Clamping::Unclamped } else { Clamping::PloidyBounds };
    opts.mosdepth_cov_genome_chr_diff = a.mosdepth_cov_genome_chr_diff;
    opts.dist_proportion = a.dist_proportion;
    opts.normalize.sample_coverage_overwrite = a.sample_coverage_overwrite;
    opts.normalize.max_std_outlier_rm = a.dev_max_std_outlier_rm as f64;
    opts.normalize.outlier_removal =
        if a.disable_max_coverage { OutlierRemoval::Disabled } else { OutlierRemoval::Enabled };
    opts.refine.min_cnv_len = a.min_cnv_len;
    opts.refine.qc_noise_allow = a.qc_noise_allow;
    opts.refine.cov_diff_threshold = a.cov_diff_threshold;

    let mut ext = if a.improved { Extensions::all() } else { Extensions::default() };
    if let Some(sex) = a.sex {
        ext.sex_chromosomes = match sex {
            SexArg::Detect => SexChromosomes::Detect,
            SexArg::Male => SexChromosomes::Given(Sex::Male),
            SexArg::Female => SexChromosomes::Given(Sex::Female),
        };
    }
    if a.gc_correction || a.improved {
        ext.gc = GcCorrection::Median { span: a.gc_span };
    }
    if a.refresh_copy_numbers {
        ext.copy_numbers = CopyNumbers::FromCall;
    }
    if a.keep_amplifications || (a.improved && !a.disable_max_coverage) {
        opts.normalize.outlier_removal = OutlierRemoval::IsolatedWindows;
    }
    ext.mappability =
        a.all_mapq_coverage.map(|path| Mappability { all_mapq_coverage: path, min_ratio: a.min_mapq_ratio });
    if a.caller == CallerArg::Hmm {
        ext.caller = Caller::Hmm(HmmOptions::default());
    }
    ext.sv_breakpoints = a.sv_breakpoints;
    opts.extensions = ext;
    opts.vcf_style = match a.vcf_style {
        VcfStyleArg::Spectre => vcf::Style::Spectre,
        VcfStyleArg::Standard => vcf::Style::Standard,
    };

    let start = Instant::now();
    opts.check_reference()?;
    std::fs::create_dir_all(&opts.output_dir).with_context(|| format!("creating {}", opts.output_dir.display()))?;
    let run = match &a.bam {
        Some(bam) => {
            let t = Instant::now();
            let reference = a.cram_reference.as_deref().unwrap_or(&opts.reference);
            let depths = crate::bam::window_depths(bam, Some(reference), 1000, 20, a.threads)
                .with_context(|| format!("counting depth in {}", bam.display()))?;
            eprintln!("lrpipe-cnv: counted depth in {} ({:.1} s)", bam.display(), t.elapsed().as_secs_f64());
            call_cnvs_with_depths(&opts, depths)?
        }
        None => call_cnvs(&opts)?,
    };
    let paths = write_outputs(&run, &opts)?;
    let p = &run.prepared;
    if let Some(d) = p.sex_depth {
        let sex = p.sex.map_or("unclear, sex chromosomes treated as autosomes", |s| match s {
            Sex::Male => "male",
            Sex::Female => "female",
        });
        eprintln!("lrpipe-cnv: chrX at {:.2} and chrY at {:.2} of autosome depth: {sex}", d.x, d.y);
    }
    if opts.extensions.mappability.is_some() {
        eprintln!("lrpipe-cnv: {} windows masked as ambiguously mapped", p.masked_windows);
    }
    if let Some(n) = run.snapped {
        eprintln!("lrpipe-cnv: {n} calls took SV breakpoints");
    }
    let n: usize = run.calls.iter().map(|(_, c)| c.len()).sum();
    let loh: usize = run.loh.as_ref().map_or(0, |l| l.passing.iter().map(|(_, c)| c.len()).sum());
    eprintln!(
        "lrpipe-cnv: {n} CNVs and {loh} LoH regions in {} ({:.1} s); also {} and {}",
        paths.vcf.display(),
        start.elapsed().as_secs_f64(),
        paths.bed.display(),
        paths.spc.display()
    );
    Ok(())
}

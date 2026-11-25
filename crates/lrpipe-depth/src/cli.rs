use std::ffi::OsString;
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Parser;

use crate::counter::{Filter, Mode};
use crate::regions::By;
use crate::report::ReportOptions;
use crate::run::{RunOptions, run};

/// Binned read depth from an indexed BAM or CRAM, with the same options and output files as
/// mosdepth 0.3.14 (per-base and quantized output are not supported).
#[derive(Parser)]
#[command(version, about)]
struct Args {
    /// Output prefix: writes {prefix}.mosdepth.summary.txt, {prefix}.mosdepth.global.dist.txt,
    /// and with --by {prefix}.regions.bed.gz and {prefix}.mosdepth.region.dist.txt
    prefix: String,
    /// The indexed alignment file
    #[arg(value_name = "BAM-or-CRAM")]
    alignments: PathBuf,

    /// Contigs counted at once (each needs 4 bytes per base of its contig)
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    /// Only this chromosome
    #[arg(short, long)]
    chrom: Option<String>,
    /// A BED file or an integer window size
    #[arg(short, long, value_name = "bed|window")]
    by: Option<String>,
    /// Don't output per-base depth (required: per-base output is not implemented)
    #[arg(short, long)]
    no_per_base: bool,
    /// Reference FASTA for CRAM (default: $REF_PATH)
    #[arg(short, long)]
    fasta: Option<PathBuf>,

    /// Exclude reads with any of these flag bits
    #[arg(short = 'F', long = "flag", default_value_t = 1796)]
    flag: u16,
    /// Only include reads with any of these flag bits (0: no requirement)
    #[arg(short, long = "include-flag", default_value_t = 0)]
    include_flag: u16,
    /// Count each read's whole span; don't look inside CIGARs or correct mate overlaps
    #[arg(short = 'x', long)]
    fast_mode: bool,
    /// Count the whole fragment of proper pairs (read 1 only)
    #[arg(short = 'a', long)]
    fragment_mode: bool,
    /// Not supported
    #[arg(short, long, hide = true)]
    quantize: Option<String>,
    /// Ignore reads with a lower mapping quality
    #[arg(short = 'Q', long, default_value_t = 0)]
    mapq: u8,
    /// Ignore reads with a smaller |insert size|
    #[arg(short = 'l', long = "min-frag-len", default_value_t = -1, allow_negative_numbers = true)]
    min_frag_len: i64,
    /// Ignore reads with a larger |insert size| (-1: no limit)
    #[arg(short = 'u', long = "max-frag-len", default_value_t = -1, allow_negative_numbers = true)]
    max_frag_len: i64,
    /// For each region, count bases with at least each of these depths (comma-separated)
    #[arg(short = 'T', long, value_delimiter = ',')]
    thresholds: Vec<i32>,
    /// Report each region's median instead of its mean
    #[arg(short = 'm', long)]
    use_median: bool,
    /// Only count these read groups (comma-separated IDs)
    #[arg(short = 'R', long, value_delimiter = ',')]
    read_groups: Vec<String>,
}

/// Runs the command line `args` (the program name first), as the `lrpipe-depth` binary does.
///
/// # Errors
///
/// If the arguments are invalid or the run fails.
pub fn main(args: Vec<OsString>) -> Result<()> {
    let a = Args::parse_from(args);
    if a.quantize.is_some() {
        bail!("--quantize is not supported");
    }
    if !a.no_per_base {
        bail!("per-base output is not supported; pass -n/--no-per-base");
    }
    if a.fast_mode && a.fragment_mode {
        bail!("only one of --fast-mode and --fragment-mode can be specified");
    }
    let max_frag_len = if a.max_frag_len < 0 { i64::MAX } else { a.max_frag_len };
    if max_frag_len < a.min_frag_len {
        bail!("--max-frag-len was lower than --min-frag-len");
    }
    if !a.thresholds.is_empty() && a.by.is_none() {
        bail!("--thresholds can only be used when --by is specified");
    }
    let precision = match std::env::var("MOSDEPTH_PRECISION") {
        Ok(p) => p.parse().context("invalid MOSDEPTH_PRECISION")?,
        Err(_) => 2,
    };
    let mut thresholds = a.thresholds;
    thresholds.sort();
    let reference = a.fasta.or_else(|| std::env::var_os("REF_PATH").map(PathBuf::from).filter(|p| p.is_file()));

    let opts = RunOptions {
        mode: if a.fast_mode {
            Mode::Fast
        } else if a.fragment_mode {
            Mode::Fragment
        } else {
            Mode::Cigar
        },
        filter: Filter {
            min_mapq: a.mapq,
            min_frag_len: a.min_frag_len,
            max_frag_len,
            exclude_flags: a.flag,
            include_flags: a.include_flag,
            read_groups: a.read_groups.into_iter().map(String::into_bytes).collect(),
        },
        report: ReportOptions {
            by: a.by.as_deref().map(By::parse).transpose()?,
            use_median: a.use_median,
            thresholds,
            precision,
            regions_only: false,
        },
        chrom: a.chrom.as_deref().map(contig_of),
        threads: a.threads,
    };
    run(&a.alignments, reference.as_deref(), &a.prefix, &opts)
}

/// mosdepth's `region_line_to_region` (mosdepth.nim:211-227): only the contig of a
/// `chrom:start-end` region is used.
fn contig_of(region: &str) -> String {
    if !region.contains(':') {
        return region.to_owned();
    }
    region.rsplitn(3, [':', '-']).nth(2).unwrap_or(region).to_owned()
}

//! lrpipe-methyl's command line.

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};

use crate::cpg::Sites;
use crate::input::{Fetched, Haplotags, ReadFilter};
use crate::mods::{Bases, Code};
use crate::pileup::{Calling, Mods, Motif, Partition, Rows, Strands, Thresholds};
use crate::read::Convention;
use crate::run::{self, Coordinates, Estimate, Job, Report};
use pasta::sam::alignment::record::Flags;

#[derive(Debug, Parser)]
#[command(version, propagate_version = true, about = "Base modification pileups from long reads (ONT and PacBio)")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Counts modification calls per reference position into bedMethyl (modkit's columns)
    Pileup(PileupArgs),
    /// CpG methylation from 5mC calls, combined and per haplotype (pb-CpG-tools' columns)
    Cpg(CpgArgs),
}

#[derive(Debug, clap::Parser)]
struct PileupArgs {
    /// Indexed BAM or CRAM with MM/ML tags
    reads: PathBuf,
    /// Output bedMethyl (.gz for BGZF; `-` or `stdout` for the standard output, as modkit)
    output: PathBuf,
    /// Reference FASTA with .fai (needed for CRAM)
    #[arg(long = "ref")]
    reference: Option<PathBuf>,
    /// Regions (contig or contig:start-end, the start 0-based as modkit reads it); every contig
    /// by default
    #[arg(long)]
    region: Vec<String>,
    /// Only the positions in these regions (BED), from which the pass thresholds are estimated
    /// too, as modkit's --include-bed
    #[arg(long, value_name = "BED", conflicts_with = "region")]
    include_bed: Option<PathBuf>,
    /// Worker threads
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Bases per window processed by one thread at a time [default: 1000000 for CRAM, whose
    /// slices span 100 kb or more and are decoded again by each window they reach; 250000
    /// otherwise]
    #[arg(long)]
    window: Option<i64>,
    /// Pass threshold, a probability, for every base (0.8) or one base (C:0.8); may repeat
    #[arg(long, conflicts_with = "filter_percentile")]
    filter_threshold: Vec<String>,
    /// Count every call, failing none
    #[arg(long, conflicts_with_all = ["filter_threshold", "filter_percentile"])]
    no_filtering: bool,
    /// Without --filter-threshold: the percentile of call confidences that fail
    #[arg(long, default_value_t = 0.1)]
    filter_percentile: f64,
    /// Without --filter-threshold: reads sampled to estimate thresholds
    #[arg(long, default_value_t = 10_042)]
    num_reads: usize,
    /// OUTPUT is a directory with combined.bedmethyl and, from the HP tag, hp1.bedmethyl
    /// and hp2.bedmethyl
    #[arg(long)]
    phased: bool,
    /// The reads' haplotypes from a haplotag list (WhatsHap's `haplotag --output-haplotag-list`)
    /// in place of the HP tags, as `haplotag` would tag them: a listed read's primary alignment
    /// gets its haplotype, every other alignment none (lrpipe's own option)
    #[arg(long, value_name = "LIST")]
    haplotags: Option<PathBuf>,
    /// Only these modifications, where the reference has their base: 5mC, 5hmC, 6mA, 4mC,
    /// BASE:CODE (C:m), or a base alone (C) with --combine-mods
    #[arg(long, num_args = 1..)]
    modified_bases: Vec<String>,
    /// Only CpG sites (needs --ref)
    #[arg(long)]
    cpg: bool,
    /// Sum each CpG's two strands onto its C (needs --cpg)
    #[arg(long)]
    combine_strands: bool,
    /// Count every modification of a base as modified, in one row per base
    #[arg(long)]
    combine_mods: bool,
}

#[derive(Debug, clap::Args)]
struct CpgArgs {
    /// Indexed BAM or CRAM with MM/ML tags (5mC)
    #[arg(long)]
    bam: PathBuf,
    /// Reference FASTA with .fai
    #[arg(long = "ref")]
    reference: PathBuf,
    /// Prefix of the outputs: PREFIX.combined.bed.gz, PREFIX.hap1.bed.gz, PREFIX.hap2.bed.gz
    #[arg(long, default_value = "lrpipe_cpg")]
    output_prefix: String,
    /// Write only PREFIX.combined.bed.gz, for reads without HP tags, whose haplotype tables
    /// would be empty (lrpipe's own option; pb-CpG-tools writes all three)
    #[arg(long, conflicts_with = "haplotags")]
    combined_only: bool,
    /// The reads' haplotypes from a haplotag list (WhatsHap's `haplotag --output-haplotag-list`)
    /// in place of the HP tags, as `haplotag` would tag them: a listed read's primary alignment
    /// gets its haplotype, every other alignment none (lrpipe's own option)
    #[arg(long, value_name = "LIST")]
    haplotags: Option<PathBuf>,
    /// How a site's methylation is estimated
    #[arg(long, value_enum, default_value_t = PileupMode::Count)]
    pileup_mode: PileupMode,
    /// Which CpG sites are scored
    #[arg(long, value_enum, default_value_t = SitesArg::Reference)]
    modsites_mode: SitesArg,
    /// Sites with fewer calls are left out
    #[arg(long, default_value_t = 4)]
    min_coverage: u32,
    /// Reads with a lower mapping quality are left out
    #[arg(long, default_value_t = 1)]
    min_mapq: u8,
    /// Regions (contig or contig:start-end, 1-based); every contig by default
    #[arg(long)]
    region: Vec<String>,
    /// Only the CpG sites in these regions (BED) (lrpipe's own option)
    #[arg(long, value_name = "BED", conflicts_with = "region")]
    include_bed: Option<PathBuf>,
    /// Worker threads
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// Bases per window processed by one thread at a time [default: 1000000 for CRAM, whose
    /// slices span 100 kb or more and are decoded again by each window they reach; 250000
    /// otherwise]
    #[arg(long)]
    window: Option<i64>,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
enum PileupMode {
    /// Share of calls with probability over one half
    Count,
    /// Expected level under a fixed model of methylated and unmethylated reads' values,
    /// fitted on PacBio Revio reads against EM-seq
    Model,
}

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
enum SitesArg {
    /// CpG in the reference
    Reference,
    /// CpG in most reads
    Denovo,
}

/// Says what was left out, and warns when no alignment had modification calls.
fn report(fetched: &Fetched) {
    if fetched.unplaceable > 0 {
        eprintln!(
            "lrpipe-methyl: {} alignments left out: MM/ML tags that do not fit the sequence",
            fetched.unplaceable
        );
    }
    if fetched.kept == 0 && fetched.untagged > 0 {
        eprintln!("lrpipe-methyl: warning: no alignment has MM/ML tags, so every output is empty");
    }
}

/// The window size asked for, or the default for the reads' format.
fn window(asked: Option<i64>, reads: &Path) -> i64 {
    asked.unwrap_or(if crate::input::is_cram(reads) { 1_000_000 } else { 250_000 })
}

fn cpg(args: &CpgArgs) -> Result<()> {
    let sites = match args.modsites_mode {
        SitesArg::Reference => Sites::Reference,
        SitesArg::Denovo => Sites::Denovo,
    };
    let path = |suffix: &str| PathBuf::from(format!("{}.{suffix}.bed.gz", args.output_prefix));
    let job = Job {
        reads: args.bam.clone(),
        reference: Some(args.reference.clone()),
        regions: regions(
            &args.region,
            args.include_bed.as_deref(),
            Coordinates::OneBased,
            &args.bam,
            Some(&args.reference),
        )?,
        coordinates: Coordinates::OneBased,
        window: window(args.window, &args.bam),
        threads: args.threads,
        // pb-CpG-tools counts supplementary alignments too.
        filter: ReadFilter {
            skip: Flags::UNMAPPED.union(Flags::SECONDARY),
            min_mapq: args.min_mapq,
            bases: Bases::only(b"C"),
        },
        report: Report::CpgCount {
            sites,
            estimate: match args.pileup_mode {
                PileupMode::Count => Estimate::Count,
                PileupMode::Model => Estimate::Model,
            },
            min_coverage: args.min_coverage,
            min_mapq: args.min_mapq,
        },
        outputs: if args.combined_only {
            vec![(Partition::All, path("combined"))]
        } else {
            vec![
                (Partition::All, path("combined")),
                (Partition::Haplotype(1), path("hap1")),
                (Partition::Haplotype(2), path("hap2")),
            ]
        },
        haplotags: read_haplotags(args.haplotags.as_deref())?,
    };
    report(&run::run(&job).context("cpg")?);
    report_haplotags(&job, args.haplotags.as_deref());
    Ok(())
}

/// The regions a job counts, as `coordinates` spell them: `--region`'s, or those of the BED of
/// `--include-bed`, sorted and merged, with their contigs in the order of the header of `reads`
/// (decoded with `reference`), so that the output is in its order.
///
/// # Errors
///
/// If the BED or the reads' header cannot be read.
fn regions(
    region: &[String],
    include_bed: Option<&Path>,
    coordinates: Coordinates,
    reads: &Path,
    reference: Option<&Path>,
) -> Result<Vec<String>> {
    let Some(path) = include_bed else { return Ok(region.to_vec()) };
    let first = match coordinates {
        Coordinates::OneBased => 1,
        Coordinates::ZeroBased => 0,
    };
    let reference = reference.map(crate::input::SharedReference::open).transpose()?;
    let header = crate::input::Reads::open(reads, reference.as_ref())?.header().clone();
    let names: Vec<String> =
        header.reference_sequences().keys().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
    let bed = lrpipe_reads::regions::Regions::read(path)?.in_order(names.iter().map(String::as_str));
    let spans: Vec<String> = bed
        .iter()
        .flat_map(|(contig, intervals)| intervals.iter().map(move |(s, e)| format!("{contig}:{}-{e}", s + first)))
        .collect();
    if spans.is_empty() {
        bail!("--include-bed {} has no regions", path.display());
    }
    Ok(spans)
}

/// `--haplotags`, read whole.
fn read_haplotags(path: Option<&Path>) -> Result<Option<Arc<Haplotags>>> {
    let Some(path) = path else { return Ok(None) };
    let list = Haplotags::read(path)?;
    eprintln!(
        "lrpipe-methyl: {} reads assigned a haplotype by {} ({:.0} MiB)",
        list.len(),
        path.display(),
        list.heap_bytes() as f64 / f64::from(1 << 20)
    );
    if list.ambiguous() > 0 {
        eprintln!(
            "lrpipe-methyl: warning: {} reads listed more than once on a contig with different haplotypes (a list \
             written with --tag-supplementary?) are left untagged",
            list.ambiguous()
        );
    }
    Ok(Some(Arc::new(list)))
}

/// A `--modified-bases` name as a canonical base and code (`None` for every code).
fn modified_base(name: &str) -> Result<(u8, Option<Code>)> {
    Ok(match name {
        "5mC" => (b'C', Some(Code::Letter(b'm'))),
        "5hmC" => (b'C', Some(Code::Letter(b'h'))),
        "5fC" => (b'C', Some(Code::Letter(b'f'))),
        "5caC" => (b'C', Some(Code::Letter(b'c'))),
        "4mC" => (b'C', Some(Code::Chebi(21839))),
        "6mA" | "m6A" => (b'A', Some(Code::Letter(b'a'))),
        "inosine" => (b'A', Some(Code::Chebi(17596))),
        _ => match name.split_once(':') {
            Some((b, c)) if b.len() == 1 => {
                let code = match c.parse::<u32>() {
                    Ok(n) => Code::Chebi(n),
                    Err(_) if c.len() == 1 => Code::Letter(c.as_bytes()[0]),
                    Err(_) => bail!("unknown modified base {name:?}"),
                };
                (b.as_bytes()[0].to_ascii_uppercase(), Some(code))
            }
            None if name.len() == 1 => (name.as_bytes()[0].to_ascii_uppercase(), None),
            _ => bail!("unknown modified base {name:?}"),
        },
    })
}

/// The pass thresholds of `--filter-threshold`, in the forms modkit takes: a probability for
/// every base (`0.8`), or one for a primary base (`C:0.8`, with `A`, `C`, `G` or `T`). A base's
/// own threshold wins over the one for every base, whichever comes first.
fn filter_thresholds(specs: &[String]) -> Result<Thresholds> {
    let probability = |spec: &str, p: &str| -> Result<f32> {
        let p: f32 = p.parse().ok().filter(|p| (0.0..=1.0).contains(p)).with_context(|| {
            format!("--filter-threshold {spec:?}: the threshold is not a probability between 0 and 1")
        })?;
        Ok(p)
    };
    let mut every = None;
    let mut per_base = Vec::new();
    for spec in specs {
        match spec.split_once(':') {
            None if every.is_some() => bail!("--filter-threshold {spec:?}: repeated threshold for every base"),
            None => every = Some(probability(spec, spec)?),
            Some((base @ ("A" | "C" | "G" | "T"), _)) if per_base.iter().any(|&(b, _)| b == base.as_bytes()[0]) => {
                bail!("--filter-threshold {spec:?}: repeated threshold for base {base}")
            }
            Some((base @ ("A" | "C" | "G" | "T"), p)) => per_base.push((base.as_bytes()[0], probability(spec, p)?)),
            Some((base, _)) => bail!(
                "--filter-threshold {spec:?}: {base:?} is not a primary base (A, C, G or T); give a base's \
                 threshold as C:0.8, or one for every base as 0.8"
            ),
        }
    }
    let mut t = Thresholds::NONE;
    if let Some(p) = every {
        t.set(b'N', p);
    }
    for (base, p) in per_base {
        t.set(base, p);
    }
    Ok(t)
}

/// Runs the command line `args` (the program name first), as the `lrpipe-methyl` binary does.
///
/// # Errors
///
/// If the arguments are invalid or the run fails.
pub fn main(args: Vec<OsString>) -> Result<()> {
    let cli = Cli::parse_from(args);
    match cli.command {
        Commands::Pileup(args) => pileup(&args),
        Commands::Cpg(args) => cpg(&args),
    }
}

fn pileup(args: &PileupArgs) -> Result<()> {
    let job = pileup_job(args)?;
    report(&run::run(&job).context("pileup")?);
    report_haplotags(&job, args.haplotags.as_deref());
    Ok(())
}

/// How many of the reads `--haplotags` named (an alignment is counted by each window it reaches),
/// with a warning when few.
fn report_haplotags(job: &Job, path: Option<&Path>) {
    if let Some(report) = job.haplotags.as_ref().zip(path).and_then(|(list, path)| list.matches().report(path)) {
        eprintln!("lrpipe-methyl: {report}");
    }
}

/// The job `pileup` runs for `args` (the words after `pileup`), with its pass thresholds
/// estimated as `pileup` estimates them: for counting the same table within another pass over
/// the reads ([`crate::stream`]).
///
/// # Errors
///
/// If the arguments are invalid or the thresholds cannot be estimated.
pub fn pileup_job_from(args: &[OsString]) -> Result<Job> {
    let args = PileupArgs::try_parse_from(std::iter::once(OsString::from("pileup")).chain(args.iter().cloned()))?;
    pileup_job(&args)
}

fn pileup_job(args: &PileupArgs) -> Result<Job> {
    // Checked even when unused, so a mistyped option fails before any output is made.
    let thresholds = filter_thresholds(&args.filter_threshold)?;
    if !(0.0..1.0).contains(&args.filter_percentile) {
        bail!(
            "--filter-percentile {} is not in [0, 1): 0.1 fails the 10% least confident calls",
            args.filter_percentile
        );
    }
    let outputs = if args.phased {
        if crate::output::is_stdout(&args.output) {
            bail!("--phased writes a directory of tables, so OUT cannot be the standard output");
        }
        std::fs::create_dir_all(&args.output).with_context(|| format!("creating {}", args.output.display()))?;
        vec![
            (Partition::All, args.output.join("combined.bedmethyl")),
            (Partition::Haplotype(1), args.output.join("hp1.bedmethyl")),
            (Partition::Haplotype(2), args.output.join("hp2.bedmethyl")),
        ]
    } else {
        vec![(Partition::All, args.output.clone())]
    };
    if args.combine_strands && !args.cpg {
        bail!("--combine-strands needs --cpg");
    }
    if (args.cpg || !args.modified_bases.is_empty()) && args.reference.is_none() {
        bail!("--cpg and --modified-bases need --ref");
    }
    let rows = Rows {
        modified_bases: if args.modified_bases.is_empty() {
            None
        } else {
            let bases: Vec<(u8, Option<Code>)> =
                args.modified_bases.iter().map(|m| modified_base(m)).collect::<Result<_>>()?;
            if let Some((base, _)) = bases.iter().find(|(_, code)| code.is_none())
                && !args.combine_mods
            {
                bail!(
                    "--modified-bases {}: a base alone counts all its modifications as one, which needs --combine-mods; \
                     or name them (5mC 5hmC)",
                    char::from(*base)
                );
            }
            Some(bases)
        },
        motif: if args.cpg { Motif::CpG } else { Motif::Any },
        strands: if args.combine_strands { Strands::Combined } else { Strands::Separate },
        mods: if args.combine_mods { Mods::Combined } else { Mods::Separate },
    };
    // Only the bases asked for are decoded.
    let bases = rows
        .modified_bases
        .as_ref()
        .map_or(Bases::ALL, |m| Bases::only(&m.iter().map(|&(b, _)| b).collect::<Vec<u8>>()));
    // modkit tallies by its reference convention when given --modified-bases.
    let convention = if rows.modified_bases.is_some() { Convention::Reference } else { Convention::Read };
    let mut calling = Calling { thresholds: Thresholds::NONE, convention };
    let mut job = Job {
        reads: args.reads.clone(),
        reference: args.reference.clone(),
        regions: regions(
            &args.region,
            args.include_bed.as_deref(),
            Coordinates::ZeroBased,
            &args.reads,
            args.reference.as_deref(),
        )?,
        // As modkit reads them.
        coordinates: Coordinates::ZeroBased,
        window: window(args.window, &args.reads),
        threads: args.threads,
        filter: ReadFilter { bases, ..ReadFilter::PRIMARY },
        report: Report::BedMethyl { calling, rows: rows.clone() },
        outputs,
        haplotags: read_haplotags(args.haplotags.as_deref())?,
    };
    if !args.no_filtering {
        if args.filter_threshold.is_empty() {
            const WINDOWS: usize = 256;
            let per_window = args.num_reads.div_ceil(WINDOWS);
            let (t, totals) = run::estimate_thresholds(&job, args.filter_percentile, WINDOWS, per_window)
                .context("estimating pass thresholds")?;
            for (base, total) in b"ACGT".iter().copied().zip(totals) {
                if total > 0 {
                    eprintln!(
                        "lrpipe-methyl: threshold {:.7} for {}, {} calls sampled",
                        t.get(base),
                        char::from(base),
                        total
                    );
                }
            }
            calling.thresholds = t;
        } else {
            calling.thresholds = thresholds;
        }
    }
    job.report = Report::BedMethyl { calling, rows };
    Ok(job)
}

/// Says what the alignments' calls left out, and warns when none had any.
pub fn report_fetched(fetched: &Fetched) {
    report(fetched);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn specs(s: &[&str]) -> Vec<String> {
        s.iter().map(|&s| s.to_owned()).collect()
    }

    #[test]
    fn include_bed_regions_are_spelled_in_each_commands_coordinates() {
        let bed = std::env::temp_dir().join(format!("lrpipe-methyl-include-{}.bed", std::process::id()));
        std::fs::write(&bed, "chr2\t5\t9\nchr1\t10\t20\nchr1\t15\t30\n").expect("written");
        // The ONT fixture's header names chr1 before chr2.
        let reads = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../lrpipe-snv/ci/fixture/ont.bam");
        let regions = |coordinates| regions(&[], Some(&bed), coordinates, &reads, None);
        // Sorted in the header's order and merged; modkit's 0-based start, or pb-CpG-tools' 1-based.
        let zero = regions(Coordinates::ZeroBased).expect("read");
        assert_eq!(zero, specs(&["chr1:10-30", "chr2:5-9"]));
        let one = regions(Coordinates::OneBased).expect("read");
        assert_eq!(one, specs(&["chr1:11-30", "chr2:6-9"]));
        let header: pasta::sam::Header = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n".parse().expect("a header");
        for (spans, coordinates) in [(zero, Coordinates::ZeroBased), (one, Coordinates::OneBased)] {
            let windows = crate::run::windows(&header, &spans, coordinates, 1000).expect("windows");
            let got: Vec<(String, i64, i64)> = windows.into_iter().map(|w| (w.contig, w.start, w.end)).collect();
            assert_eq!(got, [("chr1".to_owned(), 10, 30), ("chr2".to_owned(), 5, 9)]);
        }
        std::fs::write(&bed, "# none\n").expect("written");
        assert!(regions(Coordinates::ZeroBased).is_err(), "no regions is not every contig");
        let _ = std::fs::remove_file(bed);
    }

    fn job_error(args: &[&str]) -> String {
        let args: Vec<OsString> = args.iter().map(OsString::from).collect();
        format!("{:#}", pileup_job_from(&args).expect_err("rejected"))
    }

    #[test]
    fn filter_thresholds_take_modkits_forms() {
        let t = filter_thresholds(&specs(&["C:0.8"])).expect("parses");
        assert_eq!((t.get(b'C'), t.get(b'A')), (0.8, 0.0));
        // A base's own threshold wins over the one for every base, in either order.
        for order in [["A:0.7", "0.9"], ["0.9", "A:0.7"]] {
            let t = filter_thresholds(&specs(&order)).expect("parses");
            assert_eq!((t.get(b'A'), t.get(b'C'), t.get(b'T')), (0.7, 0.9, 0.9));
        }
    }

    /// LPH-005: a modification name, a lowercase base (a modification code in modkit), or a
    /// value that is not a probability is an error, not a threshold set for no base.
    #[test]
    fn malformed_filter_thresholds_are_errors() {
        for bad in ["5mC:0.9", "m:0.9", "c:0.8", "N:0.8", "CG:0.8", ":0.8", "C:", "C:1.5", "-0.1", "high", "C:0.5:1"] {
            assert!(filter_thresholds(&specs(&[bad])).is_err(), "{bad} accepted");
        }
        let e = job_error(&["in.bam", "out.bed", "--filter-threshold", "5mC:0.99"]);
        assert!(e.contains("\"5mC\" is not a primary base"), "{e}");
    }

    /// LPH-005: a percentile that no call reaches is an error, not a threshold of 0.
    #[test]
    fn a_repeated_filter_threshold_is_an_error() {
        // As modkit: one threshold per base, and one for every base.
        for repeated in [&["C:0.8", "C:0.6"][..], &["0.8", "0.7"], &["C:0.8", "A:0.7", "C:0.8"]] {
            let err = filter_thresholds(&specs(repeated)).expect_err("repeated").to_string();
            assert!(err.contains("repeated threshold"), "{repeated:?}: {err}");
        }
        assert!(filter_thresholds(&specs(&["0.8", "C:0.6", "A:0.7"])).is_ok(), "every base and two bases");
    }

    #[test]
    fn filter_percentile_outside_0_to_1_is_an_error() {
        for bad in ["1.5", "1", "-0.1", "NaN"] {
            let e = job_error(&["in.bam", "out.bed", &format!("--filter-percentile={bad}")]);
            assert!(e.contains("is not in [0, 1)"), "{bad}: {e}");
        }
    }
}

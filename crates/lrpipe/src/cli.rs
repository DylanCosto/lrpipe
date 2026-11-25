//! lrpipe's own commands: `run`, `setup`, `sex`, `merge`, `qc` and `compare`. (The binary runs the tools' subcommands
//! before this sees them.)

use std::collections::BTreeMap;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::SystemTime;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand, ValueEnum};
use signal_hook::consts::{SIGHUP, SIGINT, SIGQUIT, SIGTERM};

use crate::checks;
use crate::compare::Run;
use crate::config::{MethylationReport, Site};
use crate::cores::{self, Cores};
use crate::input::{self, Format, Input as Reads};
use crate::manifest::{
    self, HAPLOTYPES_INPUT, HAPLOTYPES_NONE, HAPLOTYPES_RUN, Manifest, ModelRecord, RunStatus, SampleRecord,
    SexEvidence, StageRecord,
};
use crate::merge::{self, Input, Source};
use crate::plan::{
    self, Analyses, Analysis, Gvcf, Haplotag, HaplotaggedReads, Haplotypes, Intermediates, Layout, Merged, MergedVcf,
    Options, Pars, PlanError, Platform, Sample, Sex, Stage, StageName, SvMode, Threads,
};
use crate::provenance::{FileRecord, ToolRecord, utc};
use crate::run::{self, Event, Rerun, RunError, Settings};
use crate::sex;
use crate::snv_model::{self, Choice, Given, Mismatch, ModelError, Request};
use crate::tui;
use lrpipe_reads::par::Build;
use lrpipe_reads::{IndexedReads, Reference};

/// `eprintln!` that never panics: a run's progress goes to stderr, which may stop taking it (a
/// pipe whose reader exited, a terminal closed after `disown`) while hours of stages are still
/// to finish. With `panic = "abort"` a failed write would kill lrpipe and every stage with it
/// and leave no manifest (LRUN-05), so the line is dropped instead.
macro_rules! say {
    ($($arg:tt)*) => {{
        use std::io::Write as _;
        let _ = writeln!(std::io::stderr(), $($arg)*);
    }};
}

#[derive(Debug, Parser)]
#[command(
    version,
    propagate_version = true,
    about = "Long-read calling: the pipeline on a sample, and each of its tools",
    after_help = "Tools (each takes its own tool's options; `lrpipe <TOOL> --help` lists them):
  snv     SNVs and indels with Clair3 2.0.3's run_clair3.sh options and output (`run`, implied),
          and `pileup` and `phase` for single steps [lrpipe-snv]
  sv      structural variants with Sniffles 2.8.1's options and output, and depth, tandem repeats
          and copy number from the same pass [lrpipe-sv]
  cnv     copy number variants with Spectre 0.2.1's CNVCaller options and output [lrpipe-cnv]
  tr      tandem repeat genotypes [lrpipe-tr]
  depth   binned read depth with mosdepth 0.3.14's options and output [lrpipe-depth]
  methyl  base modification pileups: `pileup` as modkit's, `cpg` as pb-CpG-tools' [lrpipe-methyl]
  phase   read-based phasing with WhatsHap 2.8's options and output (`phase`, implied), and its
          `haplotag`, `genotype`, `stats`, `compare`, `split`, ... [lrpipe-phase]
  haplotag  `lrpipe phase haplotag`: tag reads by haplotype as WhatsHap's haplotag"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Run the pipeline on one sample
    Run(RunArgs),
    /// Fetch a genome's reference, Clair3's models and the catalogues from lrpipe's manifest,
    /// check them, build what is missing and write DIR/site.toml
    Setup(SetupArgs),
    /// Estimate a sample's sex chromosomes from its BAM or CRAM index
    Sex {
        /// Sorted, indexed BAM or CRAM
        reads: PathBuf,
    },
    /// One sorted, indexed VCF from a run's variant VCFs, each record tagged with its source
    /// (INFO/SOURCE); `lrpipe run` writes it as SAMPLE.vcf.gz
    Merge(MergeArgs),
    /// One sorted, indexed BAM from several sorted, indexed BAMs and CRAMs of a sample, with the
    /// same @SQ lines; `lrpipe run` given several aligned --reads writes it as reads/SAMPLE.bam
    MergeReads(MergeReadsArgs),
    /// A run's QC summary as JSON and a table: the reads, coverage, sex, SNVs and indels, SVs,
    /// phasing and CpGs, from the run's files; `lrpipe run` writes it as qc/SAMPLE.qc.json
    Qc(crate::qc::QcArgs),
    /// Compare two runs' outputs record by record, and their tools, times and peak memory;
    /// exits with 1 if any output's records differ
    Compare {
        /// The earlier run's output directory
        old: PathBuf,
        /// The later run's
        new: PathBuf,
    },
}

#[derive(Debug, clap::Args)]
struct RunArgs {
    /// Site configuration (TOML): reference, models, catalogues and tool binaries
    #[arg(long)]
    config: PathBuf,
    /// The reads: sorted, indexed BAMs or CRAMs (several are merged into OUTPUT/reads/SAMPLE.bam
    /// first); or unaligned reads (an unaligned BAM, as dorado and PacBio write them, FASTQ or
    /// FASTQ.gz), which bonefish aligns into OUTPUT/align/SAMPLE.bam before the other stages run.
    /// Give --reads once for each file
    #[arg(long, required = true)]
    reads: Vec<PathBuf>,
    /// Sample name, in every output's sample column and file names
    #[arg(long)]
    sample: String,
    /// Alternate @RG SM label for the same donor in one aligned BAM or CRAM. Repeat for each
    /// label other than --sample, after verifying that all the reads are from that donor
    #[arg(long = "sample-alias", value_name = "NAME")]
    sample_aliases: Vec<String>,
    /// Sequencing platform
    #[arg(long, value_enum)]
    platform: PlatformArg,
    /// Sex chromosomes: one copy of chrX (outside the PARs) and chrY for male; auto reads them
    /// from the index (chrX and chrY reads against the autosomes' per base)
    #[arg(long, value_enum, default_value_t = SexArg::Auto)]
    sex: SexArg,
    /// Call haploid SNVs and small indels across a bacterial reference, without sex inference
    #[arg(long, conflicts_with = "sex")]
    bacterial: bool,
    /// Output directory
    #[arg(long)]
    output: PathBuf,
    /// Cores the run uses in all: lrpipe-snv and the read pass, which run at the same time, each
    /// get this many threads and share them [default: all cores]
    #[arg(long, value_parser = checks::threads)]
    threads: Option<usize>,
    /// Which stages run: those whose command or inputs changed since they last succeeded, or all
    #[arg(long, value_enum, default_value_t = RerunArg::Changed)]
    rerun: RerunArg,
    /// GiB the stages running together may use (32, 32G or 32GiB; at least 1); a stage waits
    /// while its expected peak would not fit beside the running ones' [default: memory available
    /// at the start, within the cgroup's limit (a SLURM job's, a container's)]
    #[arg(long, value_parser = checks::memory_gib)]
    memory: Option<f64>,
    /// Analyses to run, comma-separated; what they need runs too (cnv needs snv's calls and the
    /// read pass's depth and SVs) [default: all, tr when the configuration has a repeat catalog]
    #[arg(long, value_enum, value_delimiter = ',')]
    analyses: Option<Vec<AnalysisArg>>,
    /// Also write SAMPLE.vcf.gz: the chosen variant analyses' VCFs (SNVs and indels, phased with
    /// `phase`, SVs, CNVs and tandem repeats) as one; each tool's own files are written either way
    #[arg(long, value_enum, default_value_t = OnOff::On)]
    merged_vcf: OnOff,
    /// The Clair3 model for the reads, in place of the site configuration's or the one chosen
    /// from the reads' header: a model's name in the configuration's snv.models directory (as
    /// r1041_e82_400bps_sup_v500), or a model directory (a path with a /)
    #[arg(long)]
    model: Option<PathBuf>,
    /// What a model for another platform or chemistry than the reads' header says does: stops
    /// the run, or runs with a warning
    #[arg(long, value_enum, default_value_t = MismatchArg::Error)]
    model_mismatch: MismatchArg,
    /// Also write lrpipe-snv's gVCF, snv/merge_output.gvcf.gz, for joint genotyping (GLnexus)
    #[arg(long)]
    gvcf: bool,
    /// Also write the reads haplotagged from the run's phasing, phase/SAMPLE.haplotagged.bam with
    /// its index, for IGV and other tools (needs phase among the analyses); the SVs' phases and
    /// the CpG tables per haplotype take the haplotag list either way
    #[arg(long)]
    haplotagged_reads: bool,
    /// Haplotag the reads from the run's phased SNVs, so the SVs' phases, the tandem repeats'
    /// alleles and the CpG tables per haplotype follow the SNVs' phasing. Off, they follow the
    /// reads' own HP tags, if any, and run beside lrpipe-snv instead of after it (a few minutes
    /// less on a whole genome); the SNVs are phased either way
    #[arg(long, value_enum, default_value_t = OnOff::On)]
    haplotag: OnOff,
    /// How SVs are called: accuracy (lrpipe-sv's --accuracy: a repeat expansion's two alleles
    /// called apart and long insertions genotyped from their reads, which finds more of GIAB's
    /// SVs and genotypes them better) or sniffles2
    /// (Sniffles2 2.8.1's calls, record for record)
    #[arg(long, value_enum, default_value_t = SvModeArg::Accuracy)]
    sv_mode: SvModeArg,
    /// Restrict every analysis to the regions of this BED (a targeted panel, adaptive sampling's
    /// targets, one chromosome for a test run): SNVs, SVs, depth, tandem repeats and methylation
    /// in them, each by its tool's own region option; not cnv, which normalises over the whole
    /// genome. The reads are still all aligned (unaligned input) or merged (several files)
    #[arg(long, value_name = "BED")]
    regions: Option<PathBuf>,
    /// Keep lrpipe-snv's intermediate files
    #[arg(long)]
    keep_intermediates: bool,
    /// Print the stages' commands without running them
    #[arg(long)]
    dry_run: bool,
    /// How the run shows its stages: `tui` a full-screen display of every stage with its
    /// progress, `plain` a line for each event, `auto` the display on a terminal and lines
    /// elsewhere (a pipe, a file, a batch job)
    #[arg(long, value_enum, default_value_t = DisplayArg::Auto)]
    display: DisplayArg,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum DisplayArg {
    /// The display on a terminal, lines elsewhere
    Auto,
    /// The full-screen display (lines when stdin and stderr are not a terminal)
    Tui,
    /// A line for each event
    Plain,
}

#[derive(Debug, clap::Args)]
struct SetupArgs {
    /// The genome, whose manifest is built into lrpipe [required unless --manifest, --pin or
    /// --verify of a directory set up before]
    #[arg(long, value_enum, ignore_case = true)]
    genome: Option<GenomeArg>,
    /// Platforms whose Clair3 models are fetched, comma-separated [required unless --pin, or
    /// --verify of a directory set up before]
    #[arg(long, value_enum, value_delimiter = ',')]
    platform: Option<Vec<PlatformArg>>,
    /// The setup directory: the files, and site.toml with their absolute paths (with --pin,
    /// where the files are fetched to)
    #[arg(long)]
    dir: PathBuf,
    /// Print the URLs, checksums and curl/wget commands that fetch the files into DIR, for a
    /// machine without the internet, and fetch nothing
    #[arg(long, conflicts_with_all = ["verify", "pin"])]
    offline: bool,
    /// Check DIR's files against the manifest, reading every one again (by size where it has no
    /// checksum); exits with 1 if any is missing or wrong
    #[arg(long, conflicts_with = "pin")]
    verify: bool,
    /// Accept the terms of the files that need them (printed when they are not accepted)
    #[arg(long)]
    accept_licences: bool,
    /// A manifest of your own (TOML, with its models) in place of the genome's
    #[arg(long)]
    manifest: Option<PathBuf>,
    /// Maintainers: fetch every file of --manifest that has no SHA-256 into DIR and write
    /// their SHA-256s into the manifest
    #[arg(long, requires = "manifest")]
    pin: bool,
    /// Build the reference's index with this bonefish (DIR/REFERENCE.bonefish.raw, 8.5 GB for
    /// GRCh38), for aligning unaligned reads, and name both in site.toml ([tools] bonefish,
    /// [align] index); a later setup without it keeps them. Without it, the bonefish beside
    /// lrpipe (the release bundle's) builds it. With --offline, the setup printed to run at the
    /// end builds it
    #[arg(long, conflicts_with_all = ["verify", "pin"])]
    bonefish: Option<PathBuf>,
    /// Leave out bonefish's index of the reference, which the bundle's bonefish otherwise builds:
    /// unaligned reads are then aligned against the reference, building its index each run
    #[arg(long, conflicts_with = "bonefish")]
    no_bonefish_index: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum GenomeArg {
    /// GRCh38's no-alt analysis set
    #[value(name = "GRCh38")]
    Grch38,
    /// T2T-CHM13 v2.0's analysis set
    #[value(name = "CHM13")]
    Chm13,
    /// GRCh37 as hs37d5
    #[value(name = "GRCh37")]
    Grch37,
}

#[derive(Debug, clap::Args)]
struct MergeArgs {
    /// Reference FASTA (with its .fai): the contigs and their order
    #[arg(long)]
    reference: PathBuf,
    /// Sample name, which every input's sample column must be
    #[arg(long)]
    sample: String,
    /// Output VCF (.vcf.gz; a .tbi index is written beside it)
    #[arg(long)]
    output: PathBuf,
    /// SNVs and indels (lrpipe snv)
    #[arg(long)]
    snv: Option<PathBuf>,
    /// Structural variants (lrpipe sv)
    #[arg(long)]
    sv: Option<PathBuf>,
    /// Copy number variants (lrpipe cnv)
    #[arg(long)]
    cnv: Option<PathBuf>,
    /// Tandem repeat genotypes (lrpipe tr)
    #[arg(long)]
    tr: Option<PathBuf>,
}

#[derive(Debug, clap::Args)]
struct MergeReadsArgs {
    /// The merged BAM; its index (.bai, or .csi for a reference sequence over 2^29 - 1 bases)
    /// is written beside it
    #[arg(long)]
    output: PathBuf,
    /// Reference FASTA (with its .fai), for decoding CRAM inputs
    #[arg(long)]
    reference: Option<PathBuf>,
    /// Threads compressing the output; each input is read by a thread of its own besides
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// Sorted, indexed BAMs and CRAMs; records at the same position are written in this order
    #[arg(required = true, num_args = 2..)]
    inputs: Vec<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum AnalysisArg {
    /// Every analysis, as without --analyses (tr with a repeat catalogue configured)
    All,
    /// SNVs and indels
    Snv,
    /// The SNVs and indels phased (runs snv's calling)
    Phase,
    /// Structural variants
    Sv,
    /// Copy number variants and loss of heterozygosity (runs snv's calling, depth and sv's too)
    Cnv,
    /// Tandem repeat genotypes
    Tr,
    /// Read depth in 1 kb windows
    Depth,
    /// CpG methylation
    Methyl,
}

impl AnalysisArg {
    /// The analysis, or none for [`AnalysisArg::All`].
    fn analysis(self) -> Option<Analysis> {
        Some(match self {
            Self::All => return None,
            Self::Snv => Analysis::Snv,
            Self::Phase => Analysis::Phase,
            Self::Sv => Analysis::Sv,
            Self::Cnv => Analysis::Cnv,
            Self::Tr => Analysis::Tr,
            Self::Depth => Analysis::Depth,
            Self::Methyl => Analysis::Methyl,
        })
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum OnOff {
    On,
    Off,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SvModeArg {
    Accuracy,
    Sniffles2,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum PlatformArg {
    Hifi,
    Ont,
}

impl PlatformArg {
    fn platform(self) -> Platform {
        match self {
            Self::Hifi => Platform::Hifi,
            Self::Ont => Platform::Ont,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SexArg {
    Auto,
    Female,
    Male,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum MismatchArg {
    Error,
    Warn,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum RerunArg {
    Changed,
    All,
}

/// The names of lrpipe's own commands (`run`, `setup`, ...), and their aliases.
#[must_use]
pub fn commands() -> Vec<String> {
    use clap::CommandFactory;
    Cli::command()
        .get_subcommands()
        .flat_map(|c| std::iter::once(c.get_name().to_owned()).chain(c.get_all_aliases().map(str::to_owned)))
        .collect()
}

/// Runs the command line `args` (the program name first), as the `lrpipe` binary does. Messages
/// are prefixed with the program's name.
#[must_use]
pub fn main(args: Vec<OsString>) -> ExitCode {
    let name = args
        .first()
        .and_then(|a| Path::new(a).file_stem())
        .map_or_else(|| "lrpipe".to_owned(), |s| s.to_string_lossy().into_owned());
    let command: Vec<String> = args.iter().map(|a| a.to_string_lossy().into_owned()).collect();
    let cli = Cli::parse_from(args);
    let result = match cli.command {
        Commands::Run(args) => run(&args, &name, command).map(|()| ExitCode::SUCCESS),
        Commands::Setup(args) => setup(&args),
        Commands::Merge(args) => merge_vcfs(&args).map(|()| ExitCode::SUCCESS),
        Commands::MergeReads(args) => merge_reads(&args).map(|()| ExitCode::SUCCESS),
        Commands::Qc(args) => crate::qc::main(&args).map(|()| ExitCode::SUCCESS),
        Commands::Compare { old, new } => compare(&old, &new),
        Commands::Sex { reads } => sex::estimate(&reads).map_err(Into::into).map(|e| {
            println!("chrX {:.3}\tchrY {:.3}\t{}", e.x, e.y, e.sex.map_or("unclear", Sex::as_str));
            ExitCode::SUCCESS
        }),
    };
    result.unwrap_or_else(|e| {
        if matches!(e.downcast_ref::<RunError>(), Some(RunError::Interrupted)) {
            say!("{name}: interrupted; running it again resumes");
            return ExitCode::from(130);
        }
        say!("{name}: error: {e:#}");
        if let Some(RunError::Stage { log, .. }) = e.downcast_ref::<RunError>() {
            let tail = last_lines(log, LOG_TAIL);
            if !tail.is_empty() {
                say!("{name}: the end of {}:", log.display());
                for line in tail {
                    say!("  {line}");
                }
            }
        }
        ExitCode::from(2)
    })
}

/// Lines of a failed stage's log shown with the error: its tools end with the reason, which the
/// terminal did not show before (LUX-03).
const LOG_TAIL: usize = 10;

/// The last `n` lines of the file at `path` that are not blank, from its last 64 KiB; none if
/// it cannot be read.
fn last_lines(path: &Path, n: usize) -> Vec<String> {
    use std::io::{Read, Seek, SeekFrom};
    let Ok(mut file) = std::fs::File::open(path) else { return Vec::new() };
    let len = file.metadata().map_or(0, |m| m.len());
    let mut text = Vec::new();
    if file.seek(SeekFrom::Start(len.saturating_sub(64 << 10))).is_err() || file.read_to_end(&mut text).is_err() {
        return Vec::new();
    }
    let text = String::from_utf8_lossy(&text);
    let lines: Vec<&str> = text.lines().map(str::trim_end).filter(|l| !l.is_empty()).collect();
    lines[lines.len().saturating_sub(n)..].iter().map(|l| (*l).to_owned()).collect()
}

fn merge_vcfs(args: &MergeArgs) -> Result<()> {
    let inputs: Vec<Input> =
        [(Source::Snv, &args.snv), (Source::Sv, &args.sv), (Source::Cnv, &args.cnv), (Source::Tr, &args.tr)]
            .into_iter()
            .filter_map(|(source, path)| Some(Input { source, path: path.clone()? }))
            .collect();
    if inputs.is_empty() {
        bail!("give at least one of --snv, --sv, --cnv and --tr");
    }
    merge::merge(&inputs, &args.reference, &args.sample, &args.output)
}

fn merge_reads(args: &MergeReadsArgs) -> Result<()> {
    let options = crate::merge_reads::Options {
        inputs: args.inputs.clone(),
        output: args.output.clone(),
        reference: args.reference.clone(),
        threads: args.threads,
        command: std::env::args_os().map(|a| a.to_string_lossy().into_owned()).collect::<Vec<_>>().join(" "),
    };
    let records = crate::merge_reads::merge(&options)?;
    say!("merge-reads: {records} records from {} files into {}", args.inputs.len(), args.output.display());
    Ok(())
}

fn setup(args: &SetupArgs) -> Result<ExitCode> {
    use crate::setup::{self, Options, manifest::Genome};
    let fetcher = setup::fetch::Fetcher::new();
    if args.pin {
        let manifest = args.manifest.as_deref().expect("required by clap");
        let pinned = setup::pin(manifest, &args.dir, args.accept_licences, &fetcher)?;
        println!("pinned {} files in {}", pinned.len(), manifest.display());
        return Ok(ExitCode::SUCCESS);
    }
    let genome = args.genome.map(|g| match g {
        GenomeArg::Grch38 => Genome::Grch38,
        GenomeArg::Chm13 => Genome::Chm13,
        GenomeArg::Grch37 => Genome::Grch37,
    });
    let platforms: Option<Vec<Platform>> = args.platform.as_ref().map(|ps| {
        ps.iter().map(|p| p.platform()).fold(Vec::new(), |mut v, p| {
            if !v.contains(&p) {
                v.push(p);
            }
            v
        })
    });
    let given = match (&args.manifest, genome) {
        (Some(path), genome) => {
            let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
            let named =
                setup::manifest::Manifest::parse(&text).with_context(|| format!("in {}", path.display()))?.genome;
            if let (Some(genome), Some(named)) = (genome, &named)
                && !genome.as_str().eq_ignore_ascii_case(named)
            {
                bail!("{} is a manifest of {named}, not {}", path.display(), genome.as_str());
            }
            Some(text)
        }
        (None, Some(genome)) => Some(genome.manifest_text()),
        (None, None) => None,
    };
    let (text, platforms) = match (given, platforms) {
        (Some(text), Some(platforms)) => (text, platforms),
        (given, platforms) if args.verify => {
            let (text, recorded) = setup::recorded_setup(&args.dir)?;
            (given.unwrap_or(text), platforms.unwrap_or(recorded))
        }
        (None, _) => bail!("give --genome (or --manifest)"),
        (Some(_), None) => bail!("give --platform"),
    };
    if args.offline {
        // The manifest and bonefish as the setup to run on the other machine names them.
        let manifest = args.manifest.as_deref().map(absolute).transpose()?;
        let rerun = setup::Rerun { manifest: manifest.as_deref(), bonefish: args.bonefish.as_deref() };
        print_out(&setup::offline(&text, &absolute(&args.dir)?, &platforms, &rerun)?)?;
        return Ok(ExitCode::SUCCESS);
    }
    if args.verify {
        let verification = setup::verify(&text, &args.dir, &platforms)?;
        print_out(&verification.to_string())?;
        return Ok(if verification.passed() { ExitCode::SUCCESS } else { ExitCode::from(1) });
    }
    let bonefish = match &args.bonefish {
        Some(bonefish) => Some(absolute_tool(bonefish)?),
        None if args.no_bonefish_index => None,
        None => bundled_bonefish(),
    };
    // Before anything is fetched.
    if let Some(bonefish) = &bonefish {
        checks::bonefish(bonefish)?;
        if args.bonefish.is_none() {
            say!(
                "lrpipe setup: bonefish's index of the reference is made with {} (--no-bonefish-index leaves it out)",
                bonefish.display()
            );
        }
    }
    let options = Options { dir: args.dir.clone(), platforms, accept_licences: args.accept_licences, bonefish };
    let site = setup::setup(&text, &options, &fetcher)?;
    println!("{}", site.display());
    Ok(ExitCode::SUCCESS)
}

fn compare(old: &Path, new: &Path) -> Result<ExitCode> {
    let comparison = crate::compare::compare(&Run::open(old)?, &Run::open(new)?)?;
    print_out(&comparison.to_string())?;
    Ok(if comparison.outputs_match() { ExitCode::SUCCESS } else { ExitCode::from(1) })
}

/// Writes `text` to stdout. A reader that stopped reading (`lrpipe compare ... | head`) ends the
/// output quietly, as other Unix tools' does, where `print!` panicked.
///
/// # Errors
///
/// If stdout cannot be written for another reason.
fn print_out(text: &str) -> Result<()> {
    use std::io::Write;
    let mut out = std::io::stdout().lock();
    match out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        written => written.context("writing to stdout"),
    }
}

fn run(args: &RunArgs, name: &str, command: Vec<String>) -> Result<()> {
    checks::sample_name(&args.sample)?;
    if args.bacterial && args.analyses.iter().flatten().any(|a| !matches!(a, AnalysisArg::All | AnalysisArg::Snv)) {
        return Err(PlanError::BacterialAnalyses.into());
    }
    let requested_sex = if args.bacterial { Some(Sex::NotApplicable) } else { sex_given(args.sex) };
    let ploidy = args.bacterial.then_some(1);
    let mut site = Site::load(&args.config)?;
    if site.tools.lrpipe.is_none() {
        site.tools.lrpipe = std::env::current_exe().ok();
    }
    // Before anything reads the reference's index or the reads' (classify checks each aligned
    // file and its index).
    checks::reference(&site.reference)?;
    let given: Vec<PathBuf> = args.reads.iter().map(|r| absolute(r)).collect::<Result<_>>()?;
    let aliases = (!args.sample_aliases.is_empty()).then_some((args.sample.as_str(), args.sample_aliases.as_slice()));
    let classified = input::classify_with_sample_aliases(&given, &site.reference, aliases)?;
    let input = &classified.input;
    let layout = Layout::new(absolute(&args.output)?, &args.sample);
    if layout.root().exists() && !layout.root().is_dir() {
        bail!("--output {} is a file, not a directory", layout.root().display());
    }
    // What taking the directory would refuse, which a dry run, taking nothing, checks here.
    if args.dry_run {
        run::check(&layout)?;
    }
    let platform = args.platform.platform();
    // Checked, and sorted in the reference's order, before anything is planned from them.
    let mut region_warnings = Vec::new();
    let regions = match &args.regions {
        None => None,
        Some(bed) => {
            let contigs = fai_contigs(&site.reference)?;
            let checked = checks::regions(&absolute(bed)?, &contigs)?;
            let restriction =
                plan::restrict_contigs(&checked, &contigs, |contig, len| args.bacterial || plan::called(contig, len));
            if !restriction.left_out.is_empty() {
                say!(
                    "{name}: warning: --regions on {} are left out: a run of the whole genome calls nothing there (SNVs \
                     are called on chr1-22, X and Y, SVs on contigs of 1 Mbp or more)",
                    restriction.left_out.join(", ")
                );
            }
            if restriction.regions.is_empty() {
                bail!(
                    "--regions {}: no region is on a contig a run calls (chr1-22, X and Y, or of 1 Mbp or more)",
                    bed.display()
                );
            }
            match restriction.targets {
                None => {
                    say!(
                        "{name}: --regions cover every contig a run of the whole genome calls: the run is of the whole genome"
                    );
                    None
                }
                Some(targets) => {
                    region_warnings = checks::regions_in_reads(&restriction.regions, bed, &classified)?;
                    Some((targets, restriction.regions, contigs))
                }
            }
        }
    };
    // Several aligned files are merged, and unaligned reads aligned, into the run's own BAM, which
    // the other stages read.
    let (reads, unaligned, unaligned_fastq, merged) = match input {
        Reads::Aligned(reads) => (reads.clone(), Vec::new(), false, None),
        Reads::SeveralAligned(inputs, index) => {
            crate::merge_reads::refuse_an_input_as_output(inputs, &layout.merged_reads())?;
            (layout.merged_reads(), Vec::new(), false, Some(Merged { inputs: inputs.clone(), index: *index }))
        }
        Reads::Unaligned(files) => (
            layout.aligned_bam(),
            files.iter().map(|(p, _)| p.clone()).collect(),
            files.iter().all(|(_, f)| *f == Format::Fastq),
            None,
        ),
    };
    // The aligned reads the sex and HP tags are read from before planning: the files given (all
    // of several, which are not merged yet); for unaligned reads, the align stage's BAM, once it
    // is written.
    let aligned_reads: Vec<PathBuf> = match input {
        Reads::Aligned(reads) => vec![reads.clone()],
        Reads::SeveralAligned(inputs, _) => inputs.clone(),
        Reads::Unaligned(_) => vec![layout.aligned_bam()],
    };
    if !unaligned.is_empty() {
        let bonefish = site.tools.bonefish.clone().unwrap_or_else(default_bonefish);
        checks::bonefish(&bonefish).context(
            "bonefish aligns the unaligned --reads; set [tools] bonefish in the site configuration, or put bonefish \
             beside lrpipe or on PATH",
        )?;
        site.tools.bonefish = Some(bonefish);
    }
    let (threads, cores) = match args.threads {
        None => (Threads::each(cores::available()), Cores::All),
        Some(n) if n > cores::available() => {
            say!(
                "{name}: warning: --threads {n} is more than the {} lrpipe may use: the stages share those",
                count(cores::available(), "core")
            );
            (Threads::each(n), Cores::All)
        }
        Some(n) => match cores::confine(n).context("confining the run to --threads cores")? {
            Cores::Unconfinable => (Threads::split(n), Cores::Unconfinable),
            confined => (Threads::each(n), confined),
        },
    };
    let mut warnings = classified.warnings.clone();
    warnings.extend(region_warnings);
    // Of each aligned file: the reference has its contigs with reads. Unaligned reads' aligned BAM
    // is checked once it is written: it must have the reference's contigs exactly.
    for aligned in &classified.aligned {
        warnings.extend(checks::contigs(aligned, &site.reference)?);
    }
    // Over every file's read groups, as the merged or aligned BAM has them.
    warnings.extend(checks::sample_warning(&classified.header, &args.sample));
    let default = || {
        Analysis::ALL
            .into_iter()
            .filter(|&a| if args.bacterial { a == Analysis::Snv } else { a != Analysis::Tr || site.repeats.is_some() })
    };
    let all = args.analyses.as_ref().is_some_and(|a| a.contains(&AnalysisArg::All));
    let mut analyses: Vec<Analysis> = match &args.analyses {
        Some(chosen) if all => default().chain(chosen.iter().filter_map(|a| a.analysis())).collect(),
        Some(chosen) => chosen.iter().filter_map(|a| a.analysis()).collect(),
        None => default().collect(),
    };
    let explicit = args.analyses.is_some() && !all;
    if let Reads::Unaligned(files) = input {
        unaligned_methylation(&site, files, explicit, &mut analyses, &mut warnings)?;
    } else {
        aligned_methylation(&site, &aligned_reads, explicit, &mut analyses, &mut warnings)?;
    }
    if let Some((_, chosen, contigs)) = &regions {
        // lrpipe-cnv's copy numbers are the depth against the whole genome's (its median window),
        // which the regions' depth alone would not give: in a panel, each target's depth against
        // the targets'. Asked for, that is an error rather than calls that mislead.
        if analyses.contains(&Analysis::Cnv) {
            if args.analyses.as_ref().is_some_and(|a| a.contains(&AnalysisArg::Cnv)) {
                bail!("{}: call the CNVs in a run of the whole genome", PlanError::CnvInRegions);
            }
            analyses.retain(|&a| a != Analysis::Cnv);
            say!(
                "{name}: cnv is not called in --regions: lrpipe-cnv normalises the depth over the whole genome, not \
                 the regions' (call it in a run of the whole genome)"
            );
        }
        let reference: u64 = contigs.iter().map(|&(_, l)| l).sum();
        say!(
            "{name}: --regions: {} on {}, {} bases ({:.2}% of the reference)",
            count(chosen.len(), "region"),
            count(chosen.iter().count(), "contig"),
            chosen.bases(),
            chosen.bases() as f64 * 100.0 / reference.max(1) as f64
        );
    }
    let analyses = Analyses::of(&analyses);
    if args.haplotagged_reads && !analyses.contains(Analysis::Phase) {
        bail!("--haplotagged-reads haplotags the reads from the run's phasing: give phase among --analyses");
    }
    if args.haplotagged_reads && matches!(args.haplotag, OnOff::Off) {
        bail!("--haplotagged-reads haplotags the reads from the run's phasing, which --haplotag off leaves out");
    }
    let memory_gib = match args.memory {
        Some(gib) => gib,
        None => available_memory_gib().context("reading the memory available (or give --memory)")?,
    };
    // Chosen before the reads are aligned, from the unaligned reads' header, whose `@RG` lines
    // bonefish keeps: a run whose reads it cannot call stops before hours of aligning.
    let model = if analyses.calls_snvs() {
        checks::device(&site.snv.device)?;
        Some(snv_model(args, &site, platform, &classified)?)
    } else {
        if let Some(model) = &args.model {
            say!(
                "{name}: warning: --model {} is not used: no analysis chosen calls SNVs (snv, phase, cnv)",
                model.display()
            );
        }
        None
    };
    for warning in model.iter().flat_map(|c| &c.warnings) {
        say!("{name}: warning: {warning}");
    }
    for warning in &warnings {
        say!("{name}: warning: {warning}");
    }
    let rerun = match args.rerun {
        RerunArg::Changed => Rerun::Changed,
        RerunArg::All => Rerun::All,
    };
    let par_bed = if args.bacterial { None } else { par_bed(&site.reference)? };
    let options = |haplotypes| Options {
        threads,
        intermediates: if args.keep_intermediates { Intermediates::Keep } else { Intermediates::Remove },
        analyses,
        merged_vcf: match args.merged_vcf {
            OnOff::On => MergedVcf::On,
            OnOff::Off => MergedVcf::Off,
        },
        gvcf: if args.gvcf { Gvcf::On } else { Gvcf::Off },
        haplotypes,
        haplotag: match args.haplotag {
            OnOff::On => Haplotag::On,
            OnOff::Off => Haplotag::Off,
        },
        haplotagged_reads: if args.haplotagged_reads { HaplotaggedReads::On } else { HaplotaggedReads::Off },
        sv_mode: match args.sv_mode {
            SvModeArg::Accuracy => SvMode::Accuracy,
            SvModeArg::Sniffles2 => SvMode::Sniffles2,
        },
        pars: if par_bed.is_some() { Pars::Known } else { Pars::Unknown },
        memory_gib,
    };
    // The sex is the index's unless given, which for unaligned reads is not there yet: until then
    // the sample is planned as XX.
    let mut sample = Sample {
        name: args.sample.clone(),
        reads,
        unaligned,
        unaligned_fastq,
        merged,
        platform,
        sex: requested_sex.unwrap_or(Sex::Female),
        model: model.as_ref().map(|c| c.path.clone()),
        regions: regions.as_ref().map(|(targets, ..)| targets.clone()),
    };
    let unaligned_reads: Vec<FileRecord> =
        sample.unaligned.iter().map(|p| record("reads", p)).collect::<Result<_>>()?;
    let started = SystemTime::now();
    let mut stop = None;
    let mut lock = None;
    // The manifest of this run's align stage, and the stage records of the run before it, when
    // this run aligned the reads before planning the rest from them.
    let mut aligned: Option<(Manifest, Vec<StageRecord>)> = None;

    // The align stage, run first and alone: the sex (from the BAM's index) and the reads' HP tags,
    // which the other stages' commands depend on, are read from what it writes.
    if !sample.unaligned.is_empty() {
        // What would stop the run after aligning stops it now: an index of other contigs than
        // the reference's, the plan (with the sex given or XX, and no HP tags) and its files.
        if let Some(index) = &site.align.index {
            if !index.exists() {
                bail!(
                    "the site configuration's [align] index {} does not exist: build it with `lrpipe setup --bonefish` \
                     (or `bonefish -d INDEX --idx-raw REFERENCE`), or leave [align] index out for bonefish to index the \
                     reference on every run",
                    index.display()
                );
            }
            let bonefish = site.tools.bonefish.as_deref().expect("found above");
            same_contigs(index, &index_contigs(bonefish, index)?, &site.reference, INDEX_HINT)?;
        }
        let planned = plan::plan(&site, &sample, options(Haplotypes::Untagged), &layout)?;
        check_resources(
            &site,
            &planned,
            sample.model.as_deref().map(|m| (sample.platform, m)),
            args.regions.as_deref(),
        )?;
        // The sex is read from the aligned reads' index, which a reference of some contigs (a
        // chr20 window) cannot give: known now, not after hours of aligning. After the files are
        // checked, so a file missing is said to be.
        if requested_sex.is_none() {
            sex::check_reference(&site.reference)
                .context("estimating the sex from the aligned reads' index (or give --sex)")?;
        }
        let align = planned.into_iter().next().expect("the align stage first");
        if rerun == Rerun::All || !run::is_done(&align, &layout) {
            if let Some(warning) = memory_warning(std::slice::from_ref(&align), memory_gib, args.memory.is_some()) {
                say!("{name}: warning: {warning}");
                warnings.push(warning);
            }
            if args.dry_run {
                print_model(model.as_ref());
                print_stages(std::slice::from_ref(&align), &layout, rerun);
                if args.bacterial {
                    println!("# then haploid small-variant calling from {}", layout.aligned_bam().display());
                } else {
                    println!(
                        "# then the other stages, planned from {} once the align stage has written it (its index \
                         gives the sex, and without phasing its reads' HP tags the haplotypes)",
                        layout.aligned_bam().display()
                    );
                }
                return Ok(());
            }
            std::fs::create_dir_all(layout.root()).with_context(|| format!("creating {}", layout.root().display()))?;
            let held = run::lock(&layout)?;
            let mut earlier = manifest::earlier_stages(&layout.manifest());
            let stages = [align];
            let mut manifest = Manifest {
                lrpipe: env!("CARGO_PKG_VERSION").into(),
                command: command.clone(),
                output: layout.root().to_path_buf(),
                memory_budget_gib: round(memory_gib),
                started: utc(started),
                finished: None,
                status: RunStatus::Running,
                wall_seconds: None,
                error: None,
                sample: SampleRecord {
                    name: args.sample.clone(),
                    platform: platform.as_str().into(),
                    sex: requested_sex.map_or("auto: from the aligned reads' index", Sex::as_str).into(),
                    sex_from_index: None,
                    ploidy,
                    reads: Vec::new(),
                    unaligned_reads: unaligned_reads.clone(),
                    haplotypes: if args.bacterial { "haploid sample" } else { "not known until the reads are aligned" }
                        .into(),
                },
                snv_model: model.as_ref().map(model_record),
                resources: resources(&site, &stages, None, None)?,
                tools: tools(&stages)?,
                stages: Vec::new(),
                outputs: outputs(&stages, &layout),
                warnings: warnings.clone(),
            };
            say!(
                "{name}: {}: aligning {} with bonefish ({:.1} GiB) into {}",
                args.sample,
                match sample.unaligned.len() {
                    1 => sample.unaligned[0].display().to_string(),
                    n => format!("{n} files of reads"),
                },
                memory_gib,
                layout.aligned_bam().display()
            );
            let flag = stop_on_signals()?;
            let settings = Settings { rerun, memory_gib };
            let header = tui::Header {
                sample: args.sample.clone(),
                platform: platform.as_str().into(),
                output: layout.root().to_path_buf(),
                threads: "every core, for bonefish".into(),
                memory_gib,
                device: site.snv.device.clone(),
                model: model.as_ref().map(model_name),
                analyses: tui::analyses_line(Analysis::ALL.into_iter().filter(|&a| analyses.contains(a))),
                then: Some(format!(
                    "then, once the reads are aligned: {}, merged into one VCF, and the QC summary",
                    tui::analyses_line(Analysis::ALL.into_iter().filter(|&a| analyses.contains(a)))
                )),
            };
            shown(args.display, header, &stages, &earlier, &[], &layout, &flag, |report| {
                run::execute_after(&stages, &layout, &held, settings, &flag, &mut manifest, &earlier, report)
            })?;
            // Signalled as the align stage finished: the rest does not start.
            if flag.load(std::sync::atomic::Ordering::Relaxed) {
                finish(&mut manifest, &layout, started, RunStatus::Interrupted, None)?;
                return Err(RunError::Interrupted.into());
            }
            // Not complete: until the other stages run, the manifest is of a run under way.
            manifest.status = RunStatus::Running;
            manifest.finished = None;
            manifest.wall_seconds = None;
            manifest.write(&layout.manifest()).with_context(|| format!("writing {}", layout.manifest().display()))?;
            // The other stages' records are the run's before, whose measurements they keep; the
            // align stage's is this run's.
            earlier.retain(|r| r.name != StageName::Align);
            (stop, lock, aligned) = (Some(flag), Some(held), Some((manifest, earlier)));
        }
    }

    // Everything from here to the stages' start, which a run that aligned its reads above records
    // in its manifest if it fails.
    let calling = (|| -> Result<Option<Calling>> {
        let estimate = match requested_sex {
            Some(Sex::NotApplicable) => None,
            // Checked against the index, which a wrong --sex contradicts (as the QC summary says).
            Some(given) => {
                let index = sex::estimate_all(&aligned_reads).ok();
                if let Some(warning) = index.and_then(|e| crate::qc::sex_warning(&e, given)) {
                    say!("{name}: warning: {warning}");
                    warnings.push(warning);
                }
                None
            }
            None => {
                let e =
                    sex::estimate_all(&aligned_reads).context("estimating the sex from the index (or give --sex)")?;
                let Some(sex) = e.sex else {
                    bail!(
                        "chrX at {:.3} and chrY at {:.3} of the autosomes' reads per base fit neither XX nor XY; give \
                         --sex",
                        e.x,
                        e.y
                    );
                };
                say!("{name}: sex {} from the index (chrX {:.3}, chrY {:.3} of autosomes)", sex.as_str(), e.x, e.y);
                sample.sex = sex;
                Some(e)
            }
        };
        if !sample.unaligned.is_empty() {
            same_contigs(&sample.reads, &bam_contigs(&sample.reads)?, &site.reference, BAM_HINT)?;
            // The aligned BAM's index, as each aligned --reads was checked.
            let header = input::read_header(&sample.reads)
                .with_context(|| format!("reading the header of {}", sample.reads.display()))?;
            checks::aligned(&sample.reads, header)?;
        }
        let tags = if args.bacterial { Vec::new() } else { haplotag_samples(&aligned_reads, &site.reference)? };
        let haplotypes = if tags.iter().any(|t| t.found) { Haplotypes::Tagged } else { Haplotypes::Untagged };
        let options = options(haplotypes);
        if sample.sex == Sex::Male && par_bed.is_none() {
            say!(
                "{name}: warning: {} is not GRCh38, GRCh37 or CHM13 (by chrX's length), whose pseudo-autosomal regions \
                 lrpipe knows: chrX and chrY SNVs and indels are called diploid, and all of chrX is one copy for \
                 tandem repeats and copy number",
                site.reference.display()
            );
        }
        let stages = plan::plan(&site, &sample, options, &layout)?;
        let haplotagged = stages.iter().any(|s| s.name == StageName::Haplotag);
        // The align stage, if it runs, was warned of before it ran.
        let calling: Vec<Stage> = stages.iter().filter(|s| s.name != StageName::Align).cloned().collect();
        let memory = memory_warning(&calling, memory_gib, args.memory.is_some());
        for warning in haplotype_warning(&aligned_reads, &tags, options, &stages, &layout).into_iter().chain(memory) {
            say!("{name}: warning: {warning}");
            warnings.push(warning);
        }
        if args.dry_run {
            // What the run checks before its stages start: the files they read and the programs
            // they run.
            check_resources(
                &site,
                &stages,
                sample.model.as_deref().map(|m| (sample.platform, m)),
                args.regions.as_deref(),
            )?;
            tools(&stages)?;
            print_model(model.as_ref());
            print_stages(&stages, &layout, rerun);
            return Ok(None);
        }

        // The align stage ran above as --rerun said; it is not run again with the others.
        let rerun = if aligned.is_some() && rerun == Rerun::All {
            for stage in stages.iter().filter(|s| s.name != StageName::Align) {
                let _ = std::fs::remove_file(layout.stamp(stage.name));
            }
            Rerun::Changed
        } else {
            rerun
        };
        std::fs::create_dir_all(layout.root()).with_context(|| format!("creating {}", layout.root().display()))?;
        let lock = match lock.take() {
            Some(lock) => lock,
            None => run::lock(&layout)?,
        };
        if let Some(bed) = &par_bed {
            let path = layout.par_bed();
            std::fs::create_dir_all(path.parent().expect("in .lrpipe"))
                .with_context(|| format!("creating {}", path.display()))?;
            std::fs::write(&path, bed).with_context(|| format!("writing {}", path.display()))?;
        }
        // Rewritten only when they change, so that a large BED's time (read for a fingerprint in
        // place of its content) does not rerun its stages.
        if let Some((_, chosen, _)) = &regions {
            write_if_changed(&layout.regions_bed(), &chosen.bed())?;
            if stages.iter().any(|s| s.inputs.contains(&layout.region_windows())) {
                write_if_changed(&layout.region_windows(), &plan::windows_in(chosen))?;
            }
        }
        let reads = aligned_reads.iter().map(|r| record("reads", r)).collect::<Result<Vec<_>>>()?;
        remove_haplotype_tables(&layout, &stages)?;
        // After the lock and the reads, before the resources it is one of: a run turned away, or
        // about to fail on its reads, does not first spend minutes on it.
        if stages.iter().any(|s| s.name == StageName::Cnv)
            && let Some(mdr) = &site.cnv.metadata
        {
            make_cnv_metadata(mdr, &site.reference, name)?;
        }
        let manifest = Manifest {
            lrpipe: env!("CARGO_PKG_VERSION").into(),
            command: command.clone(),
            output: layout.root().to_path_buf(),
            memory_budget_gib: round(memory_gib),
            started: utc(started),
            finished: None,
            status: RunStatus::Running,
            wall_seconds: None,
            error: None,
            sample: SampleRecord {
                name: sample.name.clone(),
                platform: sample.platform.as_str().into(),
                sex: sample.sex.as_str().into(),
                sex_from_index: estimate.map(|e| SexEvidence { chrx: round3(e.x), chry: round3(e.y) }),
                ploidy,
                reads,
                unaligned_reads: unaligned_reads.clone(),
                haplotypes: match haplotypes {
                    _ if args.bacterial => "haploid sample",
                    _ if haplotagged => HAPLOTYPES_RUN,
                    Haplotypes::Tagged => HAPLOTYPES_INPUT,
                    Haplotypes::Untagged => HAPLOTYPES_NONE,
                }
                .into(),
            },
            snv_model: model.as_ref().map(model_record),
            resources: resources(
                &site,
                &stages,
                sample.model.as_deref().map(|m| (sample.platform, m)),
                args.regions.as_deref(),
            )?,
            tools: tools(&stages)?,
            // The align stage as this run ran it, if it did.
            stages: aligned.as_ref().map(|(m, _)| m.stages.clone()).unwrap_or_default(),
            outputs: outputs(&stages, &layout),
            warnings: warnings.clone(),
        };
        Ok(Some(Calling { options, stages, lock, manifest, settings: Settings { rerun, memory_gib } }))
    })();
    let Calling { options, stages, lock, mut manifest, settings } = match calling {
        Ok(Some(calling)) => calling,
        Ok(None) => return Ok(()),
        Err(e) => {
            if let Some((manifest, _)) = &mut aligned {
                finish(manifest, &layout, started, RunStatus::Failed, Some(format!("{e:#}")))?;
            }
            return Err(e);
        }
    };
    // A signal while the aligned reads were being read: no stage starts.
    if let Some(stop) = &stop
        && stop.load(std::sync::atomic::Ordering::Relaxed)
    {
        if let Some((manifest, _)) = &mut aligned {
            finish(manifest, &layout, started, RunStatus::Interrupted, None)?;
        }
        return Err(RunError::Interrupted.into());
    }

    say!(
        "{name}: {} ({} reads, {} for lrpipe-snv and {} for the read pass{}, {:.1} GiB) into {}",
        sample.name,
        sample.platform.as_str(),
        count(options.threads.snv, "thread"),
        options.threads.reads,
        match cores {
            Cores::Confined(n) => format!(" sharing {}", count(n, "core")),
            Cores::All | Cores::Unconfinable => String::new(),
        },
        memory_gib,
        layout.root().display()
    );
    if let Some(choice) = &model {
        say!("{name}: snv model {} ({}: {})", choice.path.display(), choice.by, choice.reason);
    }
    let stop = match stop {
        Some(stop) => stop,
        None => stop_on_signals()?,
    };
    let earlier = match &aligned {
        Some((_, earlier)) => earlier.clone(),
        None => manifest::earlier_stages(&layout.manifest()),
    };
    let header = tui::Header {
        sample: sample.name.clone(),
        platform: sample.platform.as_str().into(),
        output: layout.root().to_path_buf(),
        threads: format!("{} for lrpipe-snv, {} for the read pass", options.threads.snv, options.threads.reads),
        memory_gib,
        device: site.snv.device.clone(),
        model: model.as_ref().map(model_name),
        analyses: tui::analyses_line(Analysis::ALL.into_iter().filter(|&a| options.analyses.contains(a))),
        then: None,
    };
    let ran = manifest.stages.clone();
    let result = shown(args.display, header, &stages, &earlier, &ran, &layout, &stop, |report| {
        run::execute_after(&stages, &layout, &lock, settings, &stop, &mut manifest, &earlier, report)
    });
    if aligned.is_some() {
        // The run's time includes the alignment, which ran before these stages.
        manifest.wall_seconds = Some(round(started.elapsed().unwrap_or_default().as_secs_f64()));
        manifest.write(&layout.manifest()).with_context(|| format!("writing {}", layout.manifest().display()))?;
    }
    result?;
    if let Some(warning) = checks::leftovers(&layout, &stages) {
        say!("{name}: warning: {warning}");
    }
    // The QC summary's main values, and its warnings with the run's (and in its manifest); or,
    // as its failure does not fail the run, why there is none.
    if let Some(qc) = manifest.stages.iter().find(|s| s.name == StageName::Qc) {
        let report = match &qc.error {
            Some(error) => Err(anyhow::anyhow!("its stage failed ({error}); see {}", qc.log.display())),
            None => crate::qc::report(&layout.qc_json()),
        };
        let new = match report {
            Ok(report) => {
                for line in &report.lines {
                    say!("{name}: QC: {line}");
                }
                say!("{name}: QC: {} and {}", layout.qc_json().display(), layout.qc_tsv().display());
                report.warnings
            }
            Err(e) => vec![format!("no QC summary: {e:#}")],
        };
        let new: Vec<String> = new.into_iter().filter(|w| !warnings.contains(w)).collect();
        if !new.is_empty() {
            manifest.warnings.extend(new.iter().cloned());
            manifest.write(&layout.manifest()).with_context(|| format!("writing {}", layout.manifest().display()))?;
            warnings.extend(new);
        }
    }
    say!(
        "{name}: done in {:.0} s; manifest {}",
        manifest.wall_seconds.unwrap_or_default(),
        layout.manifest().display()
    );
    for warning in &warnings {
        say!("{name}: warning: {warning}");
    }
    Ok(())
}

/// What a run plans from its aligned reads, before its stages start.
struct Calling {
    options: Options,
    stages: Vec<Stage>,
    lock: run::OutputLock,
    manifest: Manifest,
    settings: Settings,
}

/// A warning when the run's memory (`--memory`, if `given`, else what is available) is less than
/// some of `stages` are expected to use at their peak: each of those runs with no other stage
/// beside it, and may run out of memory.
fn memory_warning(stages: &[Stage], memory_gib: f64, given: bool) -> Option<String> {
    let mut over: Vec<&Stage> = stages.iter().filter(|s| s.memory_gib > memory_gib).collect();
    if over.is_empty() {
        return None;
    }
    over.sort_by(|a, b| b.memory_gib.total_cmp(&a.memory_gib));
    let peaks: Vec<String> = over.iter().map(|s| format!("{} {:.1} GiB", s.name, s.memory_gib)).collect();
    Some(format!(
        "{} is less than {} expected to use at {} peak ({}): {} with no other stage beside it, and may run out of \
         memory",
        if given {
            format!("--memory {memory_gib:.1} GiB")
        } else {
            format!("the {memory_gib:.1} GiB of memory available")
        },
        if over.len() == 1 { "a stage is" } else { "stages are" },
        if over.len() == 1 { "its" } else { "their" },
        peaks.join(", "),
        if over.len() == 1 { "it runs" } else { "each runs" },
    ))
}

/// Finishes the manifest of a run that stopped after its align stage and before the others
/// started, as `status`, with why (`error`).
fn finish(
    manifest: &mut Manifest,
    layout: &Layout,
    started: SystemTime,
    status: RunStatus,
    error: Option<String>,
) -> Result<()> {
    manifest.status = status;
    manifest.error = error;
    manifest.finished = Some(utc(SystemTime::now()));
    manifest.wall_seconds = Some(round(started.elapsed().unwrap_or_default().as_secs_f64()));
    manifest.write(&layout.manifest()).with_context(|| format!("writing {}", layout.manifest().display()))
}

/// What to do about an index of other contigs than the reference's.
const INDEX_HINT: &str = "build bonefish's index again from the site's reference (lrpipe setup --bonefish), or set the \
                          reference it was built from";

/// What to do about aligned reads of other contigs than the reference's.
const BAM_HINT: &str = "bonefish mapped to other contigs than the reference's .fai lists: index the reference again \
                        (samtools faidx), or set the reference bonefish's index was built from";

/// The contigs, each name and length, of the reference's `.fai`.
fn fai_contigs(reference: &Path) -> Result<Vec<(String, u64)>> {
    let fai = Layout::index(reference, "fai");
    let text = std::fs::read_to_string(&fai).with_context(|| format!("reading {}", fai.display()))?;
    text.lines()
        .map(|l| {
            let mut f = l.split('\t');
            let name = f.next().unwrap_or_default().to_owned();
            let length = f.next().and_then(|n| n.parse().ok());
            length.map(|n| (name, n)).with_context(|| format!("{}: a line without a length: {l}", fai.display()))
        })
        .collect()
}

/// The contigs bonefish's `index` maps to: the `@SQ` lines bonefish prints mapping no reads to it
/// (`bonefish -a -t 1 INDEX -`, its reads from an empty input), which for a raw index takes under
/// a second.
fn index_contigs(bonefish: &Path, index: &Path) -> Result<Vec<(String, u64)>> {
    let out = std::process::Command::new(bonefish)
        .args(["-a", "-t", "1"])
        .arg(index)
        .arg("-")
        .stdin(std::process::Stdio::null())
        .output()
        .with_context(|| format!("running {}", bonefish.display()))?;
    if !out.status.success() {
        bail!(
            "reading the contigs of bonefish's index {}: {} {}",
            index.display(),
            out.status,
            String::from_utf8_lossy(&out.stderr).trim()
        );
    }
    Ok(sam_contigs(&String::from_utf8_lossy(&out.stdout)))
}

/// The contigs of a BAM's header.
fn bam_contigs(bam: &Path) -> Result<Vec<(String, u64)>> {
    let header = snv_model::header_text(bam)?;
    Ok(sam_contigs(&header))
}

/// The contigs of the `@SQ` lines of a SAM header.
fn sam_contigs(header: &str) -> Vec<(String, u64)> {
    header
        .lines()
        .filter(|l| l.starts_with("@SQ\t"))
        .filter_map(|l| {
            let field = |tag: &str| l.split('\t').skip(1).find_map(|f| f.strip_prefix(tag));
            Some((field("SN:")?.to_owned(), field("LN:")?.parse().ok()?))
        })
        .collect()
}

/// Refuses `what`, whose contigs are `theirs`, unless they are the reference's, in its order,
/// saying how they differ and what to do (`hint`).
fn same_contigs(what: &Path, theirs: &[(String, u64)], reference: &Path, hint: &str) -> Result<()> {
    let ours = fai_contigs(reference)?;
    if theirs == ours {
        return Ok(());
    }
    let contig = |c: Option<&(String, u64)>| c.map_or_else(|| "none".to_owned(), |(n, l)| format!("{n} of {l} bases"));
    let at = theirs.iter().zip(&ours).position(|(a, b)| a != b).unwrap_or(theirs.len().min(ours.len()));
    bail!(
        "{} has other contigs than {} ({} contigs against {}; the {} is {} against {}): {hint}",
        what.display(),
        Layout::index(reference, "fai").display(),
        theirs.len(),
        ours.len(),
        ordinal(at + 1),
        contig(theirs.get(at)),
        contig(ours.get(at))
    )
}

/// `n` as an English ordinal: 1st, 2nd, 3rd, 4th, 11th, 21st, ...
fn ordinal(n: usize) -> String {
    let suffix = match (n % 10, n % 100) {
        (_, 11..=13) => "th",
        (1, _) => "st",
        (2, _) => "nd",
        (3, _) => "rd",
        _ => "th",
    };
    format!("{n}{suffix}")
}

/// The sex `--sex` gives, unless it is to be read from the index.
fn sex_given(sex: SexArg) -> Option<Sex> {
    match sex {
        SexArg::Female => Some(Sex::Female),
        SexArg::Male => Some(Sex::Male),
        SexArg::Auto => None,
    }
}

/// SIGINT, SIGTERM, SIGQUIT and SIGHUP (a terminal closed, an ssh session dropped: the stages, in
/// process groups of their own, do not get it) stop the stages and finish the manifest, in
/// either phase of a run (aligning, then the rest). Another one changes nothing: the stages are
/// killed at once, and the run ends within seconds of it however they end (`timeout -s INT` and
/// a terminal's Ctrl-C can each deliver two). SIGHUP is left alone when it is ignored, as under
/// `nohup`, whose run must outlive the terminal (LRUN-05).
fn stop_on_signals() -> Result<Arc<AtomicBool>> {
    let stop = Arc::new(AtomicBool::new(false));
    let hangup = (!crate::process::ignored(SIGHUP)).then_some(SIGHUP);
    for signal in [SIGINT, SIGTERM, SIGQUIT].into_iter().chain(hangup) {
        signal_hook::flag::register(signal, Arc::clone(&stop)).context("handling signals")?;
    }
    Ok(stop)
}

/// Prints lrpipe-snv's model and why, as `--dry-run` does.
fn print_model(model: Option<&Choice>) {
    if let Some(choice) = model {
        println!("# snv model: {} ({}: {})", choice.path.display(), choice.by, choice.reason);
    }
}

/// Prints each stage's name, the stages it follows, whether it is done (a run would skip it: it
/// succeeded with the same command and inputs, as did every stage it follows) and its command, as
/// a shell would take it, as `--dry-run` does.
fn print_stages(stages: &[Stage], layout: &Layout, rerun: Rerun) {
    let mut done: Vec<StageName> = Vec::new();
    for stage in stages {
        let after: Vec<&str> = stage.after.iter().map(|s| s.as_str()).collect();
        let skipped =
            rerun == Rerun::Changed && stage.after.iter().all(|a| done.contains(a)) && run::is_done(stage, layout);
        if skipped {
            done.push(stage.name);
        }
        println!(
            "# {} (after: {}){}",
            stage.name,
            if after.is_empty() { "-".into() } else { after.join(", ") },
            if skipped { ": done, would be skipped" } else { "" }
        );
        let words: Vec<String> = std::iter::once(stage.program.as_os_str())
            .chain(stage.args.iter().map(|a| a.as_os_str()))
            .map(|w| crate::setup::shell(&w.to_string_lossy()))
            .collect();
        println!("{}", words.join(" "));
    }
}

fn model_record(choice: &Choice) -> ModelRecord {
    ModelRecord {
        path: choice.path.clone(),
        by: choice.by.clone(),
        reason: choice.reason.clone(),
        warnings: choice.warnings.clone(),
    }
}

/// The bonefish beside the lrpipe binary, where lrpipe's release bundle puts it, if there is one.
fn bundled_bonefish() -> Option<PathBuf> {
    std::env::current_exe().ok().and_then(|exe| Some(exe.parent()?.join("bonefish"))).filter(|path| path.is_file())
}

/// The bonefish the align stage runs when the configuration names none: the one beside the
/// lrpipe binary (where lrpipe's release bundle puts it), else `bonefish` on `PATH`.
fn default_bonefish() -> PathBuf {
    bundled_bonefish().unwrap_or_else(|| PathBuf::from("bonefish"))
}

/// Methylation of unaligned input. It comes from the reads' `MM`/`ML` tags, which bonefish
/// carries over from an unaligned BAM's reads and a FASTQ's comments (dorado's, or `samtools
/// fastq -T MM,ML`'s); reads basecalled without modified bases have none, and every methylation
/// output of them would be empty. When no input has them, `methyl` asked for in `--analyses`
/// (`explicit`) is an error, and the default analyses leave it out with a warning, as they leave
/// out tandem repeats without a catalogue; when some inputs have them, the run warns that the
/// others give no methylation.
fn unaligned_methylation(
    site: &Site,
    files: &[(PathBuf, Format)],
    explicit: bool,
    analyses: &mut Vec<Analysis>,
    warnings: &mut Vec<String>,
) -> Result<()> {
    if !analyses.contains(&Analysis::Methyl) || site.methylation.report == MethylationReport::None {
        return Ok(());
    }
    let mut without = Vec::new();
    for (path, format) in files {
        if !input::has_modifications(path, *format)? {
            without.push(path.display().to_string());
        }
    }
    if without.is_empty() {
        return Ok(());
    }
    let (has, its) = if without.len() == 1 { ("has", "its") } else { ("have", "their") };
    let what = format!(
        "{} {has} no base modification (MM/ML) tags on {its} first {} reads",
        without.join(", "),
        input::MODIFICATION_SAMPLE
    );
    let instead = "reads basecalled with modified bases (dorado's unaligned BAM, or FASTQ with its tags: `samtools \
                   fastq -T MM,ML`)";
    if without.len() < files.len() {
        warnings.push(format!("{what}: their reads give no methylation; the other --reads' do"));
    } else if explicit {
        bail!("{what}, so methyl would find no methylation: give {instead}, or leave methyl out of --analyses");
    } else {
        analyses.retain(|a| *a != Analysis::Methyl);
        warnings.push(format!("{what}: no methylation is called (for it, give {instead})"));
    }
    Ok(())
}

/// Methylation of aligned input, as [`unaligned_methylation`] has it of unaligned input: from
/// the reads' `MM` tags, looked for in the reads [`haplotag_samples`] looks at for `HP` (a
/// whole genome's few BGZF blocks or CRAM containers). Reads aligned without them (basecalled
/// without modified bases, or aligned from FASTQ that dropped them) would give empty methylation
/// tables and an exit status of 0.
fn aligned_methylation(
    site: &Site,
    files: &[PathBuf],
    explicit: bool,
    analyses: &mut Vec<Analysis>,
    warnings: &mut Vec<String>,
) -> Result<()> {
    if !analyses.contains(&Analysis::Methyl) || site.methylation.report == MethylationReport::None {
        return Ok(());
    }
    let samples = tag_samples(files, &site.reference, *b"MM", "base modification (MM)")?;
    // MM is a string: a record with it counts among those whose tag is not an integer.
    let without: Vec<String> = files
        .iter()
        .zip(&samples)
        .filter(|(_, s)| !s.found && s.other_type == 0)
        .map(|(path, _)| path.display().to_string())
        .collect();
    if without.is_empty() {
        return Ok(());
    }
    let looked: usize = samples.iter().map(|s| s.records).sum();
    let (has, its) = if without.len() == 1 { ("has", "its") } else { ("have", "their") };
    let what = format!(
        "{} {has} no base modification (MM/ML) tags on {its} {looked} mapped reads sampled",
        without.join(", ")
    );
    let instead = "reads basecalled with modified bases and aligned with their MM/ML tags kept (dorado's BAM aligned by \
                   dorado, minimap2 -y or pbmm2, or lrpipe run from the unaligned BAM)";
    if without.len() < files.len() {
        warnings.push(format!("{what}: their reads give no methylation; the other --reads' do"));
    } else if explicit {
        bail!("{what}, so methyl would find no methylation: give {instead}, or leave methyl out of --analyses");
    } else {
        analyses.retain(|a| *a != Analysis::Methyl);
        warnings.push(format!("{what}: no methylation is called (for it, give {instead})"));
    }
    Ok(())
}

/// Mapped reads [`haplotag_samples`] looks at for an `HP` tag in each contig it samples.
const HAPLOTAG_PER_CONTIG: usize = 200;

/// Contigs with reads [`haplotag_samples`] samples: the first long ones (1 Mb or more), so a file
/// haplotagged only on some contigs, or with few phased sites where a contig starts, is found to
/// be. On a whole genome it reads a few BGZF blocks or CRAM containers of each, however large
/// the file.
const HAPLOTAG_CONTIGS: usize = 8;

/// The reads' `HP` tags: an integer one (as lrpipe-sv and lrpipe-methyl read it) on any of
/// [`HAPLOTAG_PER_CONTIG`] mapped reads from the middle of each of the first [`HAPLOTAG_CONTIGS`]
/// contigs with reads; of each file, when the sample is in several.
fn haplotag_samples(reads: &[PathBuf], reference: &Path) -> Result<Vec<lrpipe_reads::TagSample>> {
    tag_samples(reads, reference, *b"HP", "haplotype (HP)")
}

/// `tag` in the reads [`haplotag_samples`] looks at, of each of `reads`; `what` names it in an
/// error.
fn tag_samples(reads: &[PathBuf], reference: &Path, tag: [u8; 2], what: &str) -> Result<Vec<lrpipe_reads::TagSample>> {
    let mut opened = None;
    let mut samples = Vec::with_capacity(reads.len());
    for path in reads {
        let context = || format!("checking {} for {what} tags", path.display());
        let reference = if lrpipe_reads::is_cram(path) {
            if opened.is_none() {
                opened = Some(Reference::open(reference).with_context(context)?);
            }
            opened.as_ref()
        } else {
            None
        };
        let sample = IndexedReads::open(path, reference)
            .map_err(anyhow::Error::from)
            .and_then(|mut file| Ok(file.sample_tag(tag, HAPLOTAG_PER_CONTIG, HAPLOTAG_CONTIGS)?))
            .with_context(context)?;
        samples.push(sample);
    }
    Ok(samples)
}

/// What the run warns of about haplotypes, if any planned stage reads them and the run does not
/// haplotag the reads from its own phasing: reads without `HP` tags (`tags`, of each of `reads`)
/// give unphased SVs, no CpG tables per haplotype and tandem repeats split into alleles by length
/// alone; and of several files only some haplotagged give haplotypes of their reads alone.
fn haplotype_warning(
    reads: &[PathBuf],
    tags: &[lrpipe_reads::TagSample],
    options: Options,
    stages: &[Stage],
    layout: &Layout,
) -> Vec<String> {
    // With phasing, the run haplotags the reads itself, whatever tags they carry.
    if stages.iter().any(|s| s.name == StageName::Haplotag) {
        return Vec::new();
    }
    let svs = stages.iter().any(|s| s.outputs.contains(&layout.sv_vcf()));
    let cpg = stages.iter().any(|s| s.name == StageName::Methylation);
    let repeats = stages.iter().any(|s| s.outputs.contains(&layout.repeats_vcf()));
    let readers: Vec<&str> =
        [(svs, "the SVs' phases"), (repeats, "the tandem repeats' alleles"), (cpg, "the CpG tables per haplotype")]
            .into_iter()
            .filter_map(|(planned, reader)| planned.then_some(reader))
            .collect();
    let Some((last, rest)) = readers.split_last() else { return Vec::new() };
    let readers = if rest.is_empty() { (*last).to_owned() } else { format!("{} and {last}", rest.join(", ")) };
    let names = |which: &dyn Fn(&lrpipe_reads::TagSample) -> bool| -> String {
        let named: Vec<String> =
            reads.iter().zip(tags).filter(|(_, t)| which(t)).map(|(r, _)| r.display().to_string()).collect();
        named.join(", ")
    };
    let (records, other_type) = tags.iter().fold((0, 0), |(r, o), t| (r + t.records, o + t.other_type));
    let are = if reads.len() == 1 { "is" } else { "are" };
    // The SNVs phased, but the reads not haplotagged from them (--haplotag off).
    let phased = options.analyses.contains(Analysis::Phase);
    let mut warnings = Vec::new();
    match options.haplotypes {
        Haplotypes::Untagged => warnings.push(format!(
            "{} {are} not haplotagged ({}), and `lrpipe run` haplotags reads only when it phases the SNVs{}: {}; for \
             them, {}",
            names(&|_| true),
            if other_type > 0 {
                format!(
                    "no integer HP tag in any of {records} mapped reads sampled, but HP tags of another type, read as \
                     none, on {other_type}"
                )
            } else {
                format!("no HP tag in any of {records} mapped reads sampled")
            },
            if phased { " and --haplotag is not off" } else { "" },
            {
                let effects: Vec<&str> = [
                    (svs, "the SVs are unphased"),
                    (cpg, "no CpG tables per haplotype are written"),
                    (repeats, "the tandem repeats' alleles are told apart by length alone"),
                ]
                .into_iter()
                .filter_map(|(planned, effect)| planned.then_some(effect))
                .collect();
                match effects.split_last() {
                    Some((last, [])) => (*last).to_owned(),
                    Some((last, rest)) => format!("{}, and {last}", rest.join(", ")),
                    None => unreachable!("one of them is planned"),
                }
            },
            if phased { "leave out --haplotag off" } else { "add phase to --analyses" }
        )),
        // Without phasing (with it the run haplotags the reads itself, above), the input's own
        // tags are the haplotypes, of the files that have them.
        Haplotypes::Tagged => {
            if phased {
                warnings.push(format!(
                    "{readers} follow the reads' own HP tags (--haplotag off), not the run's phasing: their \
                     haplotypes and phase sets need not be the phased SNVs'"
                ));
            }
            if tags.iter().any(|t| !t.found) {
                warnings.push(format!(
                    "{} {} haplotagged and {} not (no HP tag in their sampled reads): {readers} come from the \
                     haplotagged reads alone",
                    names(&|t| t.found),
                    if tags.iter().filter(|t| t.found).count() == 1 { "is" } else { "are" },
                    names(&|t| !t.found)
                ));
            }
        }
    }
    warnings
}

/// Removes CpG tables per haplotype that an earlier run of haplotagged reads left in the output
/// directory, when this run's methylation stage writes only the combined table (its reads are
/// not haplotagged): they would be taken for this run's.
fn remove_haplotype_tables(layout: &Layout, stages: &[Stage]) -> Result<()> {
    let Some(cpg) = stages.iter().find(|s| s.name == StageName::Methylation) else { return Ok(()) };
    if cpg.outputs.contains(&layout.cpg_table("hap1")) {
        return Ok(());
    }
    for table in ["hap1", "hap2"].map(|t| layout.cpg_table(t)) {
        match std::fs::remove_file(&table) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e).with_context(|| format!("removing {}", table.display())),
        }
    }
    Ok(())
}

/// lrpipe-snv's model: `--model` (a name in `[snv] models`, or a directory), else the site's for
/// the platform, either checked against the reads' header; else the one for the header from
/// `[snv] models`. Unaligned reads' header is that of their unaligned BAMs, whose read groups
/// bonefish keeps; FASTQ has none.
fn snv_model(args: &RunArgs, site: &Site, platform: Platform, classified: &input::Classified) -> Result<Choice> {
    let configured = match platform {
        Platform::Hifi => &site.snv.hifi_model,
        Platform::Ont => &site.snv.ont_model,
    };
    let given = match &args.model {
        Some(model) if !model.as_os_str().as_encoded_bytes().contains(&b'/') => {
            let Some(dir) = &site.snv.models else {
                bail!(
                    "--model {} names a model in the site configuration's [snv] models, which it does not set; give \
                     the model's directory (a path with a /)",
                    model.display()
                );
            };
            Some(Given { path: dir.join(model), by: "--model".into() })
        }
        Some(model) => Some(Given { path: absolute(model)?, by: "--model".into() }),
        None => {
            configured.as_ref().map(|path| Given { path: path.clone(), by: format!("snv.{}_model", platform.as_str()) })
        }
    };
    // Of several files, their read groups together, as the merged or aligned BAM has them.
    let mut header = snv_model::text_of(&classified.header).context("writing the reads' header as text")?;
    let fastq: Vec<&PathBuf> = match &classified.input {
        Reads::Unaligned(files) if files.iter().all(|(_, f)| *f == Format::Fastq) => {
            files.iter().map(|(p, _)| p).collect()
        }
        _ => Vec::new(),
    };
    // dorado's FASTQ names each read's read group, its run and basecall model, in an `RG:Z` tag.
    let mut from_tags = false;
    for path in &fastq {
        for (id, model) in input::fastq_read_groups(path)? {
            header.push_str(&format!("@RG\tID:{id}\tPL:ONT\tDS:basecall_model={model}\n"));
            from_tags = true;
        }
    }
    let headerless = !fastq.is_empty() && !from_tags;
    let mismatch = match args.model_mismatch {
        MismatchArg::Error => Mismatch::Error,
        MismatchArg::Warn => Mismatch::Warn,
    };
    let request = Request { platform, given, models: site.snv.models.as_deref(), mismatch, fastq: headerless };
    let mut choice = snv_model::choose(&header, &request).map_err(|e| match e {
        ModelError::Unknown { platform, .. } if headerless => anyhow::anyhow!(
            "FASTQ has no header saying how its reads were basecalled, which their Clair3 model depends on, and its \
             reads have no read group tags naming the basecall model (dorado's RG:Z): give the model with --model or \
             snv.{platform}_model"
        ),
        e => e.into(),
    })?;
    if from_tags {
        choice.reason.push_str(" (the FASTQ's header as its reads' RG:Z tags give it)");
    }
    checks::model(&choice.path, &choice.by, site.snv.models.as_deref())?;
    Ok(choice)
}

fn report(event: Event) {
    say!("{}", event_line(&event));
}

/// The line a run prints for `event`, with the time.
fn event_line(event: &Event) -> String {
    let now = utc(SystemTime::now());
    match event {
        Event::Skipped(s) => format!("[{now}] {s}: already done, skipped"),
        Event::Waiting(s, needs, in_use) => {
            format!("[{now}] {s}: waiting for memory (expects {needs:.1} GiB; running stages expect {in_use:.1})")
        }
        Event::Started(s, _) => format!("[{now}] {s}: started"),
        Event::Succeeded(s, secs) => format!("[{now}] {s}: done in {secs:.0} s"),
        Event::Failed(s, why) => format!("[{now}] {s}: FAILED: {why}"),
        Event::Cancelled(s) => format!("[{now}] {s}: stopped"),
        Event::Interrupted => format!("[{now}] interrupted: stopping the running stages"),
        Event::Abandoned(s, pid) => format!(
            "[{now}] {s}: killed, but still running as process {pid}; do not run again into this directory until \
             it has ended"
        ),
    }
}

/// Runs `execute` with the run's events shown as `display` asks: on the full-screen display of
/// `stages` (only the run's stages listed), or as lines.
#[expect(clippy::too_many_arguments, reason = "what the display shows, and the run it shows")]
fn shown<T>(
    display: DisplayArg,
    header: tui::Header,
    stages: &[Stage],
    earlier: &[StageRecord],
    ran: &[StageRecord],
    layout: &Layout,
    stop: &Arc<AtomicBool>,
    execute: impl FnOnce(&mut dyn FnMut(Event)) -> T,
) -> T {
    let wanted = match display {
        DisplayArg::Plain => false,
        DisplayArg::Auto | DisplayArg::Tui => tui::available(),
    };
    if display == DisplayArg::Tui && !wanted {
        say!("lrpipe: warning: --display tui needs stdin and stderr on a terminal: printing lines instead");
    }
    let screen = if wanted {
        match tui::Display::start(header, stages, earlier, ran, layout, Arc::clone(stop)) {
            Ok(screen) => Some(screen),
            Err(e) => {
                say!("lrpipe: warning: the display could not start ({e}): printing lines instead");
                None
            }
        }
    } else {
        None
    };
    match screen {
        Some(screen) => {
            let result = execute(&mut |event: Event| {
                let line = event_line(&event);
                screen.report(event, line);
            });
            screen.finish();
            result
        }
        None => execute(&mut report),
    }
}

/// lrpipe-snv's model as the display names it: its directory's name.
fn model_name(choice: &Choice) -> String {
    choice.path.file_name().map_or_else(|| choice.path.display().to_string(), |n| n.to_string_lossy().into_owned())
}

/// `n` of `noun`, plural but for one: "1 thread", "4 threads".
fn count(n: usize, noun: &str) -> String {
    format!("{n} {noun}{}", if n == 1 { "" } else { "s" })
}

fn round3(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

fn round(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}

/// Memory the kernel can give without swapping, in GiB: `MemAvailable`, or what is left under
/// the limits of lrpipe's cgroup if that is less. A SLURM job's `--mem`, `docker run -m` and a
/// Kubernetes pod's limit are cgroup limits, and `MemAvailable` is the whole machine's (LRUN-03).
#[cfg(not(target_os = "macos"))]
fn available_memory_gib() -> Result<f64> {
    let meminfo = std::fs::read_to_string("/proc/meminfo")?;
    let kib = meminfo
        .lines()
        .find_map(|l| l.strip_prefix("MemAvailable:"))
        .and_then(|v| v.trim().trim_end_matches("kB").trim().parse::<f64>().ok())
        .context("no MemAvailable line")?;
    let cgroup = std::fs::read_to_string("/proc/self/cgroup")
        .ok()
        .and_then(|cgroups| cgroup_headroom(&cgroups, Path::new("/sys/fs/cgroup")));
    let bytes = cgroup.map_or(kib * 1024.0, |left| (kib * 1024.0).min(left as f64));
    Ok(bytes / 1_073_741_824.0)
}

/// Bytes left under the tightest memory limit of the cgroups that `cgroups` (as
/// `/proc/self/cgroup` lists them) names and of their ancestors, with the cgroup file systems
/// mounted at `root` (`/sys/fs/cgroup`); none if none has a limit. A cgroup v2 limit is
/// `memory.max` or `memory.high`, less `memory.current`; a v1 limit is `memory.limit_in_bytes`,
/// less `memory.usage_in_bytes`. Either usage counts page cache, of which the inactive part
/// (`memory.stat`'s `inactive_file`, v1 `total_inactive_file`) is given back before the limit is
/// reached, so that part is not counted as used. Walking up from the process's cgroup also finds
/// the limit of a container that sees only its own cgroup at `root`.
#[cfg_attr(all(target_os = "macos", not(test)), expect(dead_code, reason = "cgroups are Linux's"))]
fn cgroup_headroom(cgroups: &str, root: &Path) -> Option<u64> {
    let mut least: Option<u64> = None;
    for line in cgroups.lines() {
        // hierarchy-id:controllers:path
        let mut fields = line.splitn(3, ':');
        let (Some(id), Some(controllers), Some(path)) = (fields.next(), fields.next(), fields.next()) else {
            continue;
        };
        let (mount, limits, usage, inactive): (PathBuf, &[&str], &str, &str) = if id == "0" && controllers.is_empty() {
            (root.to_path_buf(), &["memory.max", "memory.high"], "memory.current", "inactive_file")
        } else if controllers.split(',').any(|c| c == "memory") {
            (root.join("memory"), &["memory.limit_in_bytes"], "memory.usage_in_bytes", "total_inactive_file")
        } else {
            continue;
        };
        let mut dir = mount.join(path.trim_start_matches('/'));
        while dir.starts_with(&mount) {
            let read = |file: &str| std::fs::read_to_string(dir.join(file)).ok();
            let number = |file: &str| read(file).and_then(|v| v.trim().parse::<u64>().ok());
            // "max" is no limit.
            if let Some(limit) = limits.iter().filter_map(|f| number(f)).min() {
                let reclaimable = read("memory.stat").and_then(|stat| {
                    stat.lines().find_map(|l| l.strip_prefix(inactive)?.strip_prefix(' ')?.trim().parse::<u64>().ok())
                });
                let used = number(usage).unwrap_or(0).saturating_sub(reclaimable.unwrap_or(0));
                let left = limit.saturating_sub(used);
                least = Some(least.map_or(left, |l| l.min(left)));
            }
            if !dir.pop() {
                break;
            }
        }
    }
    least
}

/// Memory the kernel can give without swapping, in GiB: the share of physical memory that
/// macOS's memory pressure counts as available (`kern.memorystatus_level`, a percentage, which
/// counts compressible pages as well as free ones).
#[cfg(target_os = "macos")]
fn available_memory_gib() -> Result<f64> {
    let out = std::process::Command::new("/usr/sbin/sysctl")
        .args(["-n", "hw.memsize", "kern.memorystatus_level"])
        .output()?;
    if !out.status.success() {
        bail!("sysctl failed: {}", String::from_utf8_lossy(&out.stderr).trim());
    }
    let values: Vec<f64> = String::from_utf8(out.stdout)?
        .lines()
        .map(|l| l.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .context("sysctl printed a value that is not a number")?;
    let [bytes, percent] = values[..] else { bail!("sysctl printed {} values, not 2", values.len()) };
    Ok(bytes * percent / 100.0 / 1_073_741_824.0)
}

/// A tool's path made absolute, for site.toml: a bare name is looked up on `PATH`.
fn absolute_tool(tool: &Path) -> Result<PathBuf> {
    if tool.components().count() > 1 {
        return absolute(tool);
    }
    crate::provenance::resolve(tool).with_context(|| format!("{} is not on PATH", tool.display()))
}

fn absolute(path: &Path) -> Result<PathBuf> {
    std::path::absolute(path).with_context(|| format!("resolving {}", path.display()))
}

fn record(role: &str, path: &Path) -> Result<FileRecord> {
    FileRecord::cached(path).with_context(|| format!("{role} {}", path.display()))
}

/// The reference, every other data file the site configuration names that a planned stage
/// reads, lrpipe-snv's model (on its platform) and `--regions`' BED, after checking each exists.
/// A file only another analysis or platform would read (the repeat catalogue of a run without
/// `tr`, the HiFi model of an ONT run) is neither needed nor checked (LRUN-12).
fn resources(
    site: &Site,
    stages: &[Stage],
    model: Option<(Platform, &Path)>,
    regions: Option<&Path>,
) -> Result<BTreeMap<String, FileRecord>> {
    let read = |path: &Path| stages.iter().any(|s| s.inputs.iter().any(|i| i == path));
    // Each file's role (its name in the manifest) and what the site configuration calls it.
    let mut files: Vec<(&str, &str, &Path)> = vec![("reference", "reference", &site.reference)];
    let fai = PathBuf::from(format!("{}.fai", site.reference.display()));
    files.push(("reference_index", "reference's index", &fai));
    // A bgzipped reference is read through its `.gzi` too; checked here rather than by a stage.
    let gzi = PathBuf::from(format!("{}.gzi", site.reference.display()));
    if is_gzip(&site.reference).with_context(|| format!("reference {}", site.reference.display()))? {
        files.push(("reference_gzi", "reference's .gzi index", &gzi));
    }
    let optional = [
        ("tandem_repeats", "[sv] tandem_repeats", site.sv.tandem_repeats.as_deref()),
        ("repeat_catalog", "[repeats] catalog", site.repeats.as_ref().map(|r| r.catalog.as_path())),
        ("cnv_blacklist", "[cnv] blacklist", site.cnv.blacklist.as_deref()),
        ("cnv_metadata", "[cnv] metadata", site.cnv.metadata.as_deref()),
        ("align_index", "[align] index", site.align.index.as_deref()),
    ];
    files.extend(optional.into_iter().filter_map(|(role, key, p)| Some((role, key, p?))).filter(|(_, _, p)| read(p)));
    // `--regions`' BED as given (the stages read it sorted and merged, from the output directory).
    files.extend(regions.map(|p| ("regions", "--regions", p)));
    let mut out = BTreeMap::new();
    for (role, key, path) in files {
        // What names the file: the site configuration's key, or the option.
        let named = if key.starts_with("--") { key.to_owned() } else { format!("the site configuration's {key}") };
        if !path.exists() {
            bail!("{named} {} does not exist", path.display());
        }
        out.insert(role.to_owned(), record(&named, path)?);
    }
    if let Some((platform, dir)) = model.filter(|(_, dir)| read(&snv_model::network(dir, "pileup"))) {
        let role = format!("{}_model", platform.as_str());
        out.insert(role.clone(), record(&role, &snv_model::network(dir, "pileup"))?);
        let full_alignment = snv_model::network(dir, "full_alignment");
        out.insert(format!("{role}_full_alignment"), record(&role, &full_alignment)?);
    }
    Ok(out)
}

/// Checks the files [`resources`] records, before the run makes the configured cnv metadata (an
/// `.mdr` of the reference's N runs, which need not be there yet).
fn check_resources(
    site: &Site,
    stages: &[Stage],
    model: Option<(Platform, &Path)>,
    regions: Option<&Path>,
) -> Result<()> {
    let mut checked = site.clone();
    checked.cnv.metadata = site.cnv.metadata.clone().filter(|m| m.exists());
    resources(&checked, stages, model, regions).map(drop)
}

/// Writes `text` to `path` (in `.lrpipe`) unless it holds it already.
fn write_if_changed(path: &Path, text: &str) -> Result<()> {
    if std::fs::read(path).is_ok_and(|held| held == text.as_bytes()) {
        return Ok(());
    }
    std::fs::create_dir_all(path.parent().expect("in .lrpipe"))
        .with_context(|| format!("creating {}", path.display()))?;
    std::fs::write(path, text).with_context(|| format!("writing {}", path.display()))
}

/// Whether the file at `path` is gzip-compressed (by its first two bytes, as htslib tells).
fn is_gzip(path: &Path) -> std::io::Result<bool> {
    use std::io::Read;
    let mut magic = Vec::with_capacity(2);
    std::fs::File::open(path)?.take(2).read_to_end(&mut magic)?;
    Ok(magic == [0x1f, 0x8b])
}

/// Makes the configured `.mdr` of the reference's N runs if it is not there yet, as lrpipe-cnv
/// would (LRUN-12). Made before the run, it is the CNV stage's input from the first run on, so
/// a resumed run does not redo the CNVs for finding it there the second time.
fn make_cnv_metadata(mdr: &Path, reference: &Path, name: &str) -> Result<()> {
    if mdr.exists() {
        return Ok(());
    }
    say!("{name}: making {} from {} (its N runs, for lrpipe-cnv)", mdr.display(), reference.display());
    if let Some(dir) = mdr.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    }
    let defaults = lrpipe_cnv::pipeline::Options::new(PathBuf::new(), String::new(), PathBuf::new(), PathBuf::new());
    lrpipe_cnv::metadata::write_mdr(reference, defaults.n_size, mdr)
        .with_context(|| format!("making the cnv metadata {} from {}", mdr.display(), reference.display()))
}

/// The programs the stages run: bonefish for the align stage, and a tool of the same lrpipe
/// binary for each of the others.
fn tools(stages: &[Stage]) -> Result<BTreeMap<String, ToolRecord>> {
    let mut tools = BTreeMap::new();
    for stage in stages {
        let name = if stage.name == StageName::Align { "bonefish" } else { "lrpipe" };
        if !tools.contains_key(name) {
            let tool =
                ToolRecord::of(&stage.program).with_context(|| format!("{name} ({})", stage.program.display()))?;
            tools.insert(name.to_owned(), tool);
        }
    }
    Ok(tools)
}

fn outputs(stages: &[Stage], layout: &Layout) -> BTreeMap<String, PathBuf> {
    let mut out = BTreeMap::new();
    let mut add = |kind: &str, path: PathBuf| {
        out.insert(kind.to_owned(), path);
    };
    let written = |path: &PathBuf| stages.iter().any(|s| s.outputs.contains(path));
    let files = [
        ("aligned_bam", layout.aligned_bam()),
        ("merged_reads", layout.merged_reads()),
        ("vcf", layout.merged_vcf()),
        ("snv_vcf", layout.snv_vcf()),
        ("phased_snv_vcf", layout.phased_snv_vcf()),
        ("gvcf", layout.snv_gvcf()),
        ("sv_vcf", layout.sv_vcf()),
        ("snf", layout.snf()),
        ("haplotags", layout.haplotags()),
        ("haplotagged_reads", layout.haplotagged_reads()),
        ("depth_windows", layout.depth_windows()),
        ("depth_summary", layout.depth_summary()),
        ("cnv_vcf", layout.cnv_vcf()),
        ("cnv_bed", layout.cnv_bed()),
        ("repeats_vcf", layout.repeats_vcf()),
        ("bedmethyl", layout.bedmethyl()),
        ("read_stats", layout.read_stats()),
        ("qc_json", layout.qc_json()),
        ("qc_tsv", layout.qc_tsv()),
    ];
    for (kind, path) in files {
        if written(&path) {
            add(kind, path);
        }
    }
    for table in ["combined", "hap1", "hap2"] {
        if written(&layout.cpg_table(table)) {
            add(&format!("cpg_{table}"), layout.cpg_table(table));
        }
    }
    out
}

/// The PARs of `reference` as a BED, if it is a build whose PARs are known (by chrX's length in
/// its `.fai`), named as its sex chromosomes are (`chrX` or `X`).
fn par_bed(reference: &Path) -> Result<Option<String>> {
    let fai = PathBuf::from(format!("{}.fai", reference.display()));
    let text = std::fs::read_to_string(&fai).with_context(|| format!("reading {}", fai.display()))?;
    let mut x = text.lines().filter_map(|l| {
        let mut f = l.split('\t');
        Some((f.next()?, f.next()?.parse::<u64>().ok()?))
    });
    Ok(x.find(|(name, _)| matches!(*name, "chrX" | "X"))
        .and_then(|(name, len)| Some(Build::from_x_length(len)?.bed(name.strip_suffix('X').unwrap_or("")))))
}

#[cfg(test)]
mod par_tests {
    use super::*;

    #[test]
    fn par_beds_follow_the_reference() {
        let dir = std::env::temp_dir().join(format!("lrpipe-par-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let reference = dir.join("ref.fa");
        let fai =
            |x: &str| std::fs::write(dir.join("ref.fa.fai"), format!("chr1\t248956422\t6\t60\t61\n{x}\t0\t60\t61\n"));
        fai("chrX\t156040895").expect("fai");
        let bed = par_bed(&reference).expect("reads").expect("GRCh38");
        assert_eq!(bed.lines().next(), Some("chrX\t10000\t2781479"));
        assert_eq!(bed.lines().nth(2), Some("chrY\t10000\t2781479"));
        fai("X\t155270560").expect("fai");
        assert_eq!(par_bed(&reference).expect("reads").expect("GRCh37").lines().next(), Some("X\t60000\t2699520"));
        fai("chrX\t1000").expect("fai");
        assert_eq!(par_bed(&reference).expect("reads"), None);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// LRUN-06: a bgzipped reference needs its `.gzi`, checked before any stage runs; a plain
    /// one does not.
    #[test]
    fn a_bgzipped_reference_needs_its_gzi() {
        let dir = std::env::temp_dir().join(format!("lrpipe-gzi-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let site =
            |name: &str| Site::parse(&format!("reference = \"{name}\"\n[snv]\n"), &dir).expect("valid configuration");
        for name in ["plain.fa", "ref.fa.gz"] {
            let data: &[u8] = if name == "plain.fa" { b">chr1\nACGT\n" } else { &[0x1f, 0x8b, 8, 4] };
            std::fs::write(dir.join(name), data).expect("written");
            std::fs::write(dir.join(format!("{name}.fai")), "chr1\t4\t6\t4\t5\n").expect("written");
        }
        assert!(!resources(&site("plain.fa"), &[], None, None).expect("plain").contains_key("reference_gzi"));
        let e = format!("{:#}", resources(&site("ref.fa.gz"), &[], None, None).expect_err("no .gzi"));
        assert!(e.contains("reference's .gzi index") && e.contains("ref.fa.gz.gzi"), "{e}");
        std::fs::write(dir.join("ref.fa.gz.gzi"), [0u8; 8]).expect("written");
        assert!(resources(&site("ref.fa.gz"), &[], None, None).expect("with .gzi").contains_key("reference_gzi"));
        std::fs::remove_dir_all(&dir).expect("removed");
    }

    #[test]
    fn reads_the_memory_available() {
        let gib = available_memory_gib().expect("memory available");
        assert!(gib > 0.0 && gib < 1e6, "{gib} GiB");
    }

    /// A scratch `/sys/fs/cgroup` with `files` (path, contents) written under it.
    fn cgroup_root(name: &str, files: &[(&str, String)]) -> PathBuf {
        let root = std::env::temp_dir().join(format!("lrpipe-cgroup-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        for (path, text) in files {
            let path = root.join(path);
            std::fs::create_dir_all(path.parent().expect("a directory")).expect("scratch directory");
            std::fs::write(path, text).expect("written");
        }
        root
    }

    const GIB: u64 = 1 << 30;

    #[test]
    fn cgroup_v2_limits_are_read_up_the_hierarchy() {
        // LRUN-03: a SLURM job's cgroup on the unified hierarchy: the job's limit is on a cgroup
        // above the step's, which has none ("max").
        let job = "system.slice/slurmstepd.scope/job_7";
        let root = cgroup_root(
            "v2",
            &[
                (&format!("{job}/memory.max"), format!("{}\n", 16 * GIB)),
                (&format!("{job}/memory.current"), format!("{}\n", 2 * GIB)),
                (&format!("{job}/memory.stat"), format!("anon 1\ninactive_file {GIB}\n")),
                (&format!("{job}/step_0/memory.max"), "max\n".into()),
                (&format!("{job}/step_0/memory.high"), "max\n".into()),
                (&format!("{job}/step_0/memory.current"), "4096\n".into()),
                ("system.slice/memory.max", "max\n".into()),
            ],
        );
        let cgroups = format!("0::/{job}/step_0\n");
        // 16 GiB less the 2 GiB used, of which 1 GiB is inactive page cache.
        assert_eq!(cgroup_headroom(&cgroups, &root), Some(15 * GIB));

        // memory.high is a limit too; the tighter one counts.
        std::fs::write(root.join(format!("{job}/step_0/memory.high")), format!("{}", 4 * GIB)).expect("written");
        assert_eq!(cgroup_headroom(&cgroups, &root), Some(4 * GIB - 4096));

        // No limit anywhere.
        let free = cgroup_root(
            "v2-free",
            &[("user.slice/memory.max", "max\n".into()), ("user.slice/memory.current", "1\n".into())],
        );
        assert_eq!(cgroup_headroom("0::/user.slice\n", &free), None);
        let _ = std::fs::remove_dir_all(root);
        let _ = std::fs::remove_dir_all(free);
    }

    #[test]
    fn cgroup_v1_limits_are_read_from_the_memory_controller() {
        // LRUN-03: `docker run -m 2g` on cgroup v1 without a cgroup namespace: the container
        // sees its own cgroup at the controller's mount, not at the path /proc/self/cgroup names.
        let root = cgroup_root(
            "v1",
            &[
                ("memory/memory.limit_in_bytes", format!("{}\n", 2 * GIB)),
                ("memory/memory.usage_in_bytes", format!("{}\n", GIB / 2)),
                ("memory/memory.stat", format!("inactive_file 99\ntotal_inactive_file {}\n", GIB / 4)),
                ("cpu/cpu.shares", "1024\n".into()),
            ],
        );
        let cgroups = "12:cpu,cpuacct:/docker/abc\n4:memory:/docker/abc\n1:name=systemd:/docker/abc\n0::/\n";
        assert_eq!(cgroup_headroom(cgroups, &root), Some(2 * GIB - GIB / 4));

        // Nested: a job with a tighter limit than its parent's, under a root without one (v1's
        // "no limit" is a huge number).
        let nested = cgroup_root(
            "v1-nested",
            &[
                ("memory/memory.limit_in_bytes", "9223372036854771712\n".into()),
                ("memory/memory.usage_in_bytes", format!("{}\n", 100 * GIB)),
                ("memory/batch/memory.limit_in_bytes", format!("{}\n", 64 * GIB)),
                ("memory/batch/memory.usage_in_bytes", format!("{}\n", 60 * GIB)),
                ("memory/batch/job/memory.limit_in_bytes", format!("{}\n", 32 * GIB)),
                ("memory/batch/job/memory.usage_in_bytes", format!("{}\n", GIB)),
            ],
        );
        // The job has 31 GiB left under its own limit, but its parent has only 4.
        assert_eq!(cgroup_headroom("4:memory:/batch/job\n", &nested), Some(4 * GIB));
        // Usage over a limit leaves nothing.
        std::fs::write(nested.join("memory/batch/job/memory.usage_in_bytes"), format!("{}", 40 * GIB))
            .expect("written");
        assert_eq!(cgroup_headroom("4:memory:/batch/job\n", &nested), Some(0));
        assert_eq!(cgroup_headroom("", &nested), None, "no cgroups listed");
        let _ = std::fs::remove_dir_all(root);
        let _ = std::fs::remove_dir_all(nested);
    }
}

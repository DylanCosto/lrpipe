//! lrpipe-snv's command line: `run` (the whole pipeline, as `run_clair3.sh`), and `pileup` and
//! `phase` for single steps.

use std::ffi::OsString;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::{chunk, pipeline, settings};
use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use ndarray::{Array2, Array3};
use ndarray_npy::write_npy;
use snv_infer::{ModelConfig, PILEUP_OUTPUTS};
use snv_phase::PhaseJob;
use snv_pileup::{AlignmentFile, ChunkRegion, Fasta, POSITIONS, TENSOR_LEN, candidate_tensors};

#[derive(Debug, Parser)]
#[command(version, propagate_version = true, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Build pileup candidate tensors for one chunk and optionally score them, writing
    /// `PREFIX.tensors.npy`, `PREFIX.inputs.npy`, `PREFIX.info.tsv` and `PREFIX.probs.npy`.
    Pileup(PileupArgs),
    /// Phase one contig's selected het SNPs, as the original's `whatshap phase` step.
    Phase(PhaseArgs),
    /// Run the whole pipeline from the BAM to `merge_output.vcf.gz` (pileup calling, het-SNP
    /// selection, phasing, full-alignment calling and merging), in the original's output layout.
    ///
    /// Also takes the original run script's own spellings (`--bam_fn`, `-b`, `--gvcf=True`,
    /// unique prefixes), and runs when lrpipe-snv is invoked through a link named `run_clair3.sh`.
    #[command(
        infer_long_args = true,
        after_help = "Other commands: pileup and phase, single steps of the pipeline; `help` lists them."
    )]
    Run(Box<settings::RunArgs>),
}

/// Arguments of the original's intermediate `whatshap phase` call.
#[derive(Debug, clap::Args)]
struct PhaseArgs {
    /// Het SNPs selected for phasing (`<contig>.vcf`).
    #[arg(long)]
    variants: PathBuf,
    /// Sorted, indexed BAM/CRAM.
    #[arg(long)]
    bam: PathBuf,
    /// Indexed reference FASTA.
    #[arg(long = "ref")]
    reference: PathBuf,
    /// Contig to phase.
    #[arg(long)]
    ctg_name: String,
    /// Output VCF (`.vcf.gz` for bgzip).
    #[arg(long)]
    output: PathBuf,
    /// Threads (the output does not depend on it).
    #[arg(long, default_value_t = 1)]
    threads: usize,
}

/// Arguments mirroring `CallVariantsFromCffi --pileup` for one chunk.
#[derive(Debug, clap::Args)]
struct PileupArgs {
    /// Sorted, indexed BAM/CRAM.
    #[arg(long)]
    bam: PathBuf,
    /// Indexed reference FASTA.
    #[arg(long = "ref")]
    reference: PathBuf,
    /// Contig to process.
    #[arg(long)]
    ctg_name: String,
    /// 1-based chunk number (`--chunk_id`).
    #[arg(long, requires = "chunk_num")]
    chunk_id: Option<i64>,
    /// Number of chunks (`--chunk_num`).
    #[arg(long)]
    chunk_num: Option<i64>,
    /// First position (`--ctgStart`), when not chunking.
    #[arg(long, conflicts_with = "chunk_id", requires = "ctg_end")]
    ctg_start: Option<i64>,
    /// End position (`--ctgEnd`), when not chunking.
    #[arg(long)]
    ctg_end: Option<i64>,
    /// `--extend_bed`; like the original, a path that does not exist counts as not given.
    #[arg(long)]
    extend_bed: Option<PathBuf>,
    /// Minimum SNP allele frequency.
    #[arg(long, default_value_t = 0.08)]
    snp_min_af: f64,
    /// Minimum indel allele frequency.
    #[arg(long, default_value_t = 0.15)]
    indel_min_af: f64,
    /// Minimum depth.
    #[arg(long, default_value_t = 2)]
    min_coverage: usize,
    /// Minimum mapping quality.
    #[arg(long, default_value_t = 5)]
    min_mq: u8,
    /// Pileup checkpoint (`pileup.pt`); without it only tensors are written.
    #[arg(long)]
    model: Option<PathBuf>,
    /// Candidates per inference batch (`predictBatchSize * threads` in the original).
    #[arg(long, default_value_t = 200)]
    batch_size: usize,
    /// Accepted for the original's job command lines (`--threads`); the networks score on one
    /// thread.
    #[arg(long, default_value_t = 1)]
    threads: i32,
    /// Write the chunk's pileup VCF here, as `CallVariantsFromCffi --call_fn` does (needs --model).
    #[arg(long, requires = "model")]
    call_fn: Option<PathBuf>,
    /// File whose contents become the `##cmdline` header (`--cmd_fn`).
    #[arg(long)]
    cmd_fn: Option<PathBuf>,
    /// Sample name in the VCF.
    #[arg(long, default_value = "SAMPLE")]
    sample_name: String,
    /// Calls below this QUAL are LowQual (`--qual`).
    #[arg(long, default_value_t = 2)]
    qual: i64,
    /// Output prefix.
    #[arg(long)]
    out: PathBuf,
}

/// Runs lrpipe-snv with `args`, the program name first, as its binary and `lrpipe snv` do.
///
/// # Errors
///
/// Fails if the command fails.
pub fn main(args: Vec<OsString>) -> Result<()> {
    // `lrpipe snv` may have set up logging already.
    let _ = env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).try_init();
    match Cli::parse_from(settings::run_script_spellings(args)).command {
        Command::Pileup(args) => run_pileup(&args),
        Command::Run(args) => pipeline::run(&args),
        Command::Phase(args) => Ok(snv_phase::phase_contig(&PhaseJob {
            variants: &args.variants,
            alignments: &args.bam,
            reference: &args.reference,
            contig: &args.ctg_name,
            output: &args.output,
            threads: args.threads,
        })?),
    }
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    s.into()
}

fn run_pileup(args: &PileupArgs) -> Result<()> {
    let extend_bed = args.extend_bed.as_deref().filter(|p| p.is_file());
    let chunk = match (args.chunk_id, args.chunk_num, args.ctg_start, args.ctg_end) {
        (Some(id), Some(count), _, _) if extend_bed.is_some() => ChunkRegion::BedChunk { id: id - 1, count },
        (Some(id), Some(count), _, _) => ChunkRegion::FaiChunk { id: id - 1, count },
        (None, _, Some(start), Some(end)) => ChunkRegion::Range { start, end },
        _ => bail!("give --chunk-id and --chunk-num, or --ctg-start and --ctg-end"),
    };
    let params = chunk::candidate_params(args.snp_min_af, args.indel_min_af, args.min_coverage, args.min_mq);

    let mut bam = AlignmentFile::open(&args.bam, &args.reference)?;
    let fasta = Fasta::open(&args.reference)?;
    let started = Instant::now();
    let candidates = candidate_tensors(&mut bam, &fasta, &args.ctg_name, &chunk, extend_bed, &params)
        .with_context(|| format!("building tensors for {}", args.ctg_name))?;
    log::info!("{} candidates in {:.2?}", candidates.len(), started.elapsed());

    let n = candidates.len();
    let raw: Vec<i32> = candidates.iter().flat_map(|c| c.tensor.iter().copied()).collect();
    let inputs = chunk::model_inputs(&candidates)?;
    let shape = (n, POSITIONS, TENSOR_LEN / POSITIONS);
    write_npy(with_suffix(&args.out, ".tensors.npy"), &Array3::from_shape_vec(shape, raw)?)?;
    write_npy(with_suffix(&args.out, ".inputs.npy"), &Array3::from_shape_vec(shape, inputs.clone())?)?;
    let mut info = BufWriter::new(File::create(with_suffix(&args.out, ".info.tsv"))?);
    for c in &candidates {
        writeln!(info, "{}\t{}", c.position, c.alt_info)?;
    }
    info.flush()?;

    if let Some(model_path) = &args.model {
        let mut model = snv_infer::load_pileup_model(model_path, ModelConfig::default())?;
        let started = Instant::now();
        let probs = chunk::score(&mut model, &inputs, args.batch_size)?;
        log::info!("scored {n} candidates in {:.2?}", started.elapsed());
        if let Some(call_fn) = &args.call_fn {
            let output = chunk::VcfOutput {
                reference: &args.reference,
                cmd_file: args.cmd_fn.as_deref(),
                sample_name: &args.sample_name,
                qual: args.qual,
                iupac: snv_call::IupacBases::Vcf,
                max_indel_length: 50,
                likelihoods: snv_call::Likelihoods::Omit,
            };
            output.write(call_fn, &candidates, &probs)?;
        }
        write_npy(with_suffix(&args.out, ".probs.npy"), &Array2::from_shape_vec((n, PILEUP_OUTPUTS), probs)?)?;
    }
    Ok(())
}

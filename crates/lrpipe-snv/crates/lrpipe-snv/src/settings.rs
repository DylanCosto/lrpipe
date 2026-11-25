//! `lrpipe-snv run`'s arguments, named like those of the original's pipeline script (`run_*.py`),
//! and the values it resolves from them before it starts the pipeline.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use snv_call::{Haploid, IupacBases, RefCalls};
use snv_fullalign::{DwellTime, IndelLimit, MoveTables, Platform as FaPlatform};
use snv_infer::{InferenceDevice, ModelConfig, Precision};
use snv_merge::Gender;
use snv_pileup::{CandidateAlleles, CandidateParams, FlankRequirement};

use crate::gvcf::GvcfArgs;

/// `str2bool` in `shared/utils.py`: the original's boolean options take `--flag` or `--flag=VALUE`.
pub(crate) fn str2bool(value: &str) -> Result<bool, String> {
    match value.to_lowercase().as_str() {
        "yes" | "true" | "t" | "y" | "1" => Ok(true),
        "no" | "false" | "f" | "n" | "0" => Ok(false),
        _ => Err("Boolean value expected.".to_owned()),
    }
}

/// Sequencing platform (`--platform`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum Platform {
    /// Oxford Nanopore.
    Ont,
    /// PacBio HiFi.
    Hifi,
    /// Illumina.
    Ilmn,
}

impl Platform {
    /// The value the original's scripts compare against.
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Ont => "ont",
            Self::Hifi => "hifi",
            Self::Ilmn => "ilmn",
        }
    }

    /// The full-alignment matrix depth that goes with it.
    pub(crate) fn full_alignment(self) -> FaPlatform {
        match self {
            Self::Ont => FaPlatform::Ont,
            Self::Hifi => FaPlatform::Hifi,
            Self::Ilmn => FaPlatform::Ilmn,
        }
    }
}

/// How the threads that wait for the GPU wait (`--gpu-wait`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GpuWait {
    /// On a core, as CUDA does by default, so the GPU's results are taken up soonest.
    Spin,
    /// Asleep, leaving the core to phasing and building tensors for slightly later wake-ups (the
    /// default): on the whole HG002 HiFi genome 4 min 55 s instead of 5 min 10 s, on HG003 ONT
    /// 7 min 22 s instead of 7 min 25 s, with 3-6% less CPU time.
    Block,
}

/// Sample sex (`--gender`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GenderArg {
    /// No sex-aware calling.
    Unknown,
    /// Haploid-precise calling on chrX/chrY outside the PAR regions.
    Male,
    /// chrY is not called.
    Female,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum ReadOrder {
    /// Order reads tied at a coordinate by strand, name and alignment fields.
    Canonical,
    /// Use file order for ties, as Clair3 does.
    Input,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum AmbiguousBases {
    /// The VCF specification's: an ambiguous REF base as its first possible nucleotide.
    Vcf,
    /// Clair3's: ambiguous REF and ALT bases as `N`.
    N,
}

/// Arguments of `lrpipe-snv run`, named like those of the original's pipeline script (with
/// dashes; its own spellings work too, see [`run_script_spellings`]). Boolean options take
/// `--flag` or `--flag=true/false`, as the original's do.
#[derive(Debug, clap::Args)]
pub(crate) struct RunArgs {
    /// Sorted, indexed BAM/CRAM (`--bam_fn`).
    #[arg(long, short = 'b', visible_alias = "bam-fn")]
    pub(crate) bam: PathBuf,
    /// Indexed reference FASTA (`--ref_fn`).
    #[arg(long = "ref", short = 'f', visible_alias = "ref-fn")]
    pub(crate) reference: PathBuf,
    /// Model directory holding `pileup.pt` and `full_alignment.pt` (`--model_path`).
    #[arg(long, short = 'm')]
    pub(crate) model_path: PathBuf,
    /// Output directory (`--output`).
    #[arg(long, short = 'o')]
    pub(crate) output: PathBuf,
    /// Worker threads (`--threads`).
    #[arg(long, short = 't', default_value_t = 8)]
    pub(crate) threads: usize,
    /// Sequencing platform (`--platform`).
    #[arg(long, short = 'p', value_enum, default_value = "ont")]
    pub(crate) platform: Platform,
    /// Order ties in the native ONT and HiFi pileup and full-alignment stages.
    #[arg(long, value_enum, default_value = "canonical")]
    pub(crate) read_order: ReadOrder,
    /// How ambiguous (IUPAC) reference bases are written, unless `--keep-iupac-bases` keeps them.
    #[arg(long, value_enum, default_value = "vcf")]
    pub(crate) ambiguous_bases: AmbiguousBases,
    /// Only call variants in these regions (`--bed_fn`).
    #[arg(long)]
    pub(crate) bed_fn: Option<PathBuf>,
    /// Only genotype the sites of this VCF (`--vcf_fn`); sites without a call are added back as
    /// `./.`.
    #[arg(long)]
    pub(crate) vcf_fn: Option<PathBuf>,
    /// Comma-separated contigs (`--ctg_name`); default: chr1-22, X, Y (or 1-22, X, Y) present in
    /// the reference.
    #[arg(long)]
    pub(crate) ctg_name: Option<String>,
    /// Sample name (`--sample_name`).
    #[arg(long, default_value = "SAMPLE")]
    pub(crate) sample_name: String,
    /// Calls with QUAL at or below this are LowQual (`--qual`).
    #[arg(long, default_value_t = 2)]
    pub(crate) qual: i64,
    /// Chunks per contig (`--chunk_num`): only checked, the chunks come from `--chunk-size`;
    /// `-1` turns chunking off, with which the original (2.0.3) calls nothing.
    #[arg(long, default_value_t = 0, allow_negative_numbers = true)]
    pub(crate) chunk_num: i64,
    /// Contig chunk size (`--chunk_size`).
    #[arg(long, default_value_t = 5_000_000)]
    pub(crate) chunk_size: i64,
    /// Minimum SNP allele frequency of a candidate (`--snp_min_af`; 0 is the platform default).
    #[arg(long, default_value_t = 0.0)]
    pub(crate) snp_min_af: f64,
    /// Minimum indel allele frequency of a candidate (`--indel_min_af`; 0 is the platform
    /// default).
    #[arg(long, default_value_t = 0.0)]
    pub(crate) indel_min_af: f64,
    /// Proportion of variant calls re-called in full alignment (`--var_pct_full`; 0 is the
    /// platform default).
    #[arg(long, default_value_t = 0.0)]
    pub(crate) var_pct_full: f64,
    /// Proportion of reference calls re-called in full alignment (`--ref_pct_full`; 0 is the
    /// platform default).
    #[arg(long, default_value_t = 0.0)]
    pub(crate) ref_pct_full: f64,
    /// Proportion of het SNPs kept for phasing (`--var_pct_phasing`; 0 is 0.7, or 0.8 for the
    /// guppy5 models).
    #[arg(long, default_value_t = 0.0)]
    pub(crate) var_pct_phasing: f64,
    /// Pileup checkpoint name in the model directory (`--pileup_model_prefix`).
    #[arg(long, default_value = "pileup")]
    pub(crate) pileup_model_prefix: String,
    /// Full-alignment checkpoint name in the model directory (`--fa_model_prefix`).
    #[arg(long, default_value = "full_alignment")]
    pub(crate) fa_model_prefix: String,
    /// Minimum mapping quality (`--min_mq`, at least 5).
    #[arg(long, default_value_t = 5)]
    pub(crate) min_mq: i64,
    /// Minimum depth of a candidate (`--min_coverage`, at least 2).
    #[arg(long, default_value_t = 2)]
    pub(crate) min_coverage: i64,
    /// Skip contigs shorter than this (`--min_contig_size`).
    #[arg(long, default_value_t = 0)]
    pub(crate) min_contig_size: i64,
    /// ONT only: skip candidates with AF below 0.15 (`--fast_mode`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) fast_mode: bool,
    /// Stop after pileup calling (`--pileup_only`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) pileup_only: bool,
    /// Keep reference calls in the final output (`--print_ref_calls`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) print_ref_calls: bool,
    /// Call every contig of the reference (`--include_all_ctgs`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) include_all_ctgs: bool,
    /// Phase the final calls with lrpipe-phase (`--use_whatshap_for_final_output_phasing`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) use_whatshap_for_final_output_phasing: bool,
    /// The same (`--enable_phasing`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) enable_phasing: bool,
    /// Tag the reads by haplotype into `phased_output.bam`; implies final-output phasing
    /// (`--use_whatshap_for_final_output_haplotagging`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) use_whatshap_for_final_output_haplotagging: bool,
    /// Delete `tmp/` at the end (`--remove_intermediate_dir`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) remove_intermediate_dir: bool,
    /// Only SNP candidates (`--call_snp_only`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) call_snp_only: bool,
    /// Call candidates whose flanks are not covered (`--enable_variant_calling_at_sequence_head_and_tail`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) enable_variant_calling_at_sequence_head_and_tail: bool,
    /// Indels up to 100,000 bases instead of 50 (`--enable_long_indel`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) enable_long_indel: bool,
    /// Keep IUPAC reference bases instead of writing N (`--keep_iupac_bases`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) keep_iupac_bases: bool,
    /// Skip phasing; full-alignment reads are not haplotagged (`--no_phasing_for_fa`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) no_phasing_for_fa: bool,
    /// Haploid calling that drops heterozygous calls (`--haploid_precise`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) haploid_precise: bool,
    /// Haploid calling that keeps heterozygous calls as 1/1 (`--haploid_sensitive`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) haploid_sensitive: bool,
    /// Sample sex (`--gender`).
    #[arg(long, value_enum, default_value = "unknown")]
    pub(crate) gender: GenderArg,
    /// Pseudo-autosomal regions kept diploid with `--gender male` (`--par_regions_bed`).
    #[arg(long)]
    pub(crate) par_regions_bed: Option<PathBuf>,
    /// Where the networks run: `cpu` (output identical to the original's CPU run), `cuda:N` (the
    /// original's `--use_gpu --device cuda:N` calling, with the CPU run's pileup inputs and TF32
    /// as libtorch's cuDNN rounds) or `metal` (an Apple GPU, with the CPU run's inputs: `float32`
    /// calls within the last bits of its).
    #[arg(long, default_value = "cpu")]
    pub(crate) device: InferenceDevice,
    /// LongPhase binary for the `--use-longphase-*` options (`--longphase`; default: `longphase`
    /// on `PATH`).
    #[arg(long)]
    pub(crate) longphase: Option<PathBuf>,
    /// Phase het SNPs with LongPhase instead of lrpipe-phase (`--use_longphase_for_intermediate_phasing`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) use_longphase_for_intermediate_phasing: bool,
    /// The same (`--longphase_for_phasing`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) longphase_for_phasing: bool,
    /// Phase the final calls with LongPhase (`--use_longphase_for_final_output_phasing`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) use_longphase_for_final_output_phasing: bool,
    /// Run the original's Python pipeline (`samtools mpileup` text tensors, reads haplotagged into
    /// `phase_bam/`) instead of the C implementation (`--disable_c_impl`); always on for
    /// Illumina, and turned off again by `--enable-dwell-time`.
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) disable_c_impl: bool,
    /// Add each base's signal length from the Dorado move table (`mv` tag) to the
    /// full-alignment tensors, for the `_with_mv` models (`--enable_dwell_time`; ONT only).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) enable_dwell_time: bool,
    /// Run the networks on the GPU (`--use_gpu`); `--device cpu` then means `cuda:0`.
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) use_gpu: bool,
    /// lrpipe-snv only: give up output identical to the original's for speed on the CPU; the calls
    /// agree with an exact run on more than 99%. With AVX-VNNI (x86-64) or `sdot` (aarch64) the
    /// full-alignment convolutions run in 8 bits, and so does the pileup network's arithmetic if a check on the run's own
    /// candidates shows it close enough; close calls, and pileup calls near a quality cut-off, are
    /// scored again in float32. Half of the ONT variant calls go to full alignment
    /// (`--var_pct_full 0.5`). A GPU runs float32 either way. On a CPU without these instructions
    /// it has no effect, with a warning.
    #[arg(long)]
    pub(crate) approximate: bool,
    /// lrpipe-snv only: how threads wait for the GPU. `block` (the default) sleeps, leaving the core
    /// to the rest of the run; `spin` keeps a core busy, as CUDA does by default, for the soonest
    /// wake-ups.
    #[arg(long, value_enum, default_value = "block")]
    pub(crate) gpu_wait: GpuWait,
    /// lrpipe-snv only: NVRTC (`libnvrtc.so.12`), which compiles the GPU kernels, when it is not on
    /// the library path, `$CUDA_HOME`, `$CUDA_PATH` or `/usr/local/cuda`.
    #[arg(long)]
    pub(crate) nvrtc: Option<PathBuf>,
    #[command(flatten)]
    pub(crate) gvcf: GvcfArgs,
    /// Accepted and ignored, as in the original, whose pipeline does not read it
    /// (`--use_whatshap_for_intermediate_phasing`).
    #[arg(long, hide = true, num_args = 0..=1, default_value = "true", default_missing_value = "true", value_parser = str2bool)]
    pub(crate) use_whatshap_for_intermediate_phasing: bool,
    /// Accepted for the original's command lines and ignored: lrpipe-snv runs `samtools` in
    /// process.
    #[arg(long, hide = true)]
    pub(crate) samtools: Option<String>,
    /// Accepted and ignored: no step runs Python.
    #[arg(long, hide = true)]
    pub(crate) python: Option<String>,
    /// Accepted and ignored: no step runs PyPy.
    #[arg(long, hide = true)]
    pub(crate) pypy: Option<String>,
    /// Accepted and ignored: lrpipe-snv schedules the jobs GNU parallel would run.
    #[arg(long, hide = true)]
    pub(crate) parallel: Option<String>,
    /// Accepted and ignored: WhatsHap runs in process (lrpipe-phase).
    #[arg(long, hide = true)]
    pub(crate) whatshap: Option<String>,
}

/// Options for which the original's pipeline script takes the value `EMPTY` to mean "not given".
const EMPTY_MEANS_NOT_GIVEN: [&str; 6] = ["bed-fn", "vcf-fn", "ctg-name", "longphase", "device", "par-regions-bed"];

/// The original's command lines for `lrpipe-snv run`: underscores in option names become dashes
/// (`--bam_fn` is `--bam-fn`) and an option given `EMPTY` (its "not given") is dropped; values
/// are left as they are. Invoked through a link to lrpipe-snv named like the original's run
/// script (`run_clair3.sh` or `run_clair3.py`), the arguments are those of `lrpipe-snv run` and
/// read order and ambiguous bases default to Clair3's. Other command lines are returned unchanged.
pub(crate) fn run_script_spellings(args: Vec<OsString>) -> Vec<OsString> {
    let mut args = args.into_iter();
    let Some(program) = args.next() else { return Vec::new() };
    let invoked_as_run_script =
        Path::new(&program).file_name().and_then(|n| n.to_str()).is_some_and(|n| n.starts_with("run_clair3"));
    let mut out = vec![program];
    if invoked_as_run_script {
        out.push("run".into());
    } else {
        match args.next() {
            Some(sub) if sub == "run" => out.push(sub),
            Some(other) => {
                out.push(other);
                out.extend(args);
                return out;
            }
            None => return out,
        }
    }
    let mut args = args.peekable();
    let (mut has_read_order, mut has_ambiguous_bases) = (false, false);
    while let Some(arg) = args.next() {
        let Some(option) = arg.to_str().and_then(|a| a.strip_prefix("--")).filter(|o| !o.is_empty()) else {
            out.push(arg);
            continue;
        };
        let (name, value) = match option.split_once('=') {
            Some((name, value)) => (name.replace('_', "-"), Some(value.to_owned())),
            None => (option.replace('_', "-"), None),
        };
        has_read_order |= "read-order".starts_with(&name);
        has_ambiguous_bases |= "ambiguous-bases".starts_with(&name);
        if EMPTY_MEANS_NOT_GIVEN.contains(&name.as_str()) {
            match &value {
                Some(v) if v == "EMPTY" => continue,
                None if args.peek().is_some_and(|next| next == "EMPTY") => {
                    args.next();
                    continue;
                }
                _ => {}
            }
        }
        out.push(match value {
            Some(value) => format!("--{name}={value}").into(),
            None => format!("--{name}").into(),
        });
    }
    // Clair3's own: its order of tied reads and its N for ambiguous bases.
    if invoked_as_run_script && !has_ambiguous_bases {
        out.insert(2, "--ambiguous-bases=n".into());
    }
    if invoked_as_run_script && !has_read_order {
        out.insert(2, "--read-order=input".into());
    }
    out
}

/// The models the original's pipeline script gives a phasing proportion of 0.8.
const GUPPY5_MODELS: [&str; 3] = ["r941_prom_sup_g5014", "r941_prom_hac_g5014", "ont_guppy5"];

/// A phasing tool.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Phaser {
    /// lrpipe-phase, in process.
    WhatsHap,
    /// This LongPhase binary.
    LongPhase(PathBuf),
}

/// Values the original's pipeline scripts (`run_*.py` and `*_c_impl_pipeline.py`) derive from the
/// arguments.
#[derive(Debug, Clone)]
pub(crate) struct Settings {
    pub(crate) var_pct_full: f64,
    pub(crate) ref_pct_full: f64,
    pub(crate) var_pct_phasing: f64,
    pub(crate) candidates: CandidateParams,
    /// `--snp_min_af` and `--indel_min_af` as the original's pipeline script passes
    /// them on (Illumina's Python tensors compare in `float64`).
    pub(crate) min_af: (f64, f64),
    pub(crate) min_coverage: i64,
    pub(crate) min_mq: i64,
    pub(crate) min_contig_size: i64,
    /// Which phaser the final calls go through, if any.
    pub(crate) final_phasing: Option<Phaser>,
    /// The phaser for the het SNPs that haplotag full-alignment reads.
    pub(crate) intermediate_phasing: Phaser,
    pub(crate) haplotag: bool,
    pub(crate) haploid: Haploid,
    pub(crate) gender: Gender,
    pub(crate) par_regions: Option<PathBuf>,
    pub(crate) ref_calls: RefCalls,
    pub(crate) iupac: IupacBases,
    pub(crate) indel_limit: IndelLimit,
    pub(crate) dwell_time: DwellTime,
    pub(crate) pileup_model: PathBuf,
    pub(crate) fa_model: PathBuf,
    pub(crate) device: InferenceDevice,
    pub(crate) gpu_wait: GpuWait,
    #[cfg_attr(not(feature = "gpu"), expect(dead_code, reason = "only GPU calling loads NVRTC"))]
    pub(crate) nvrtc: Option<PathBuf>,
    pub(crate) precision: Precision,
}

impl Settings {
    /// How to run a network on `device` at `precision`.
    pub(crate) fn model_config(&self, device: InferenceDevice, precision: Precision) -> ModelConfig {
        ModelConfig { device, precision }
    }
}

/// The original pipeline script's dwell-time checks: move-table (`_with_mv`) models are ONT-only
/// and need `--enable_dwell_time`, which needs an `mv` tag with values in the first 50 records.
fn dwell_time(args: &RunArgs, model_name: &str) -> Result<DwellTime> {
    let move_table_model = [model_name, &args.pileup_model_prefix, &args.fa_model_prefix]
        .iter()
        .any(|name| name.to_lowercase().contains("with_mv") || name.to_lowercase().contains("with_move"));
    if move_table_model {
        if args.platform != Platform::Ont {
            bail!(
                "Model '{model_name}' is a move-table (signal-aware) model and is ONT-only, but --platform is '{}'. \
                 Use '--platform ont' with ONT data, or choose a non move-table model for this platform.",
                args.platform.as_str()
            );
        }
        if !args.enable_dwell_time {
            bail!(
                "Model '{model_name}' is a move-table (signal-aware) model and requires the dwell-time feature. Please \
                 add '--enable_dwell_time'. The input BAM must also contain Dorado 'mv' tags (basecalled with \
                 '--emit-moves'); if it has none, use a non move-table ONT model instead."
            );
        }
    }
    if !args.enable_dwell_time {
        return Ok(DwellTime::Off);
    }
    if args.platform != Platform::Ont {
        bail!("Dwell time is not supported for non-ONT platform, disabling option! exit!");
    }
    match snv_fullalign::sample_move_tables(&args.bam, &args.reference, 50)
        .map_err(|e| anyhow::anyhow!("Failed to check BAM file for 'mv' tag: {e}"))?
    {
        MoveTables::Found => Ok(DwellTime::On),
        MoveTables::Empty { checked } => bail!(
            "Dwell time is enabled but within the first {checked} alignments, an 'mv' tag was found without a valid \
             value. Please ensure 'mv' is populated (e.g. mv:B:* with values) or disable --enable_dwell_time."
        ),
        MoveTables::Missing { checked } => bail!(
            "Dwell time is enabled but no valid 'mv' tag was found in the first {checked} alignments. The 'mv' tag \
             (move table) is required for dwell time analysis. Please provide a BAM file containing 'mv' tags \
             earlier in the file or disable --enable_dwell_time."
        ),
    }
}

/// `model_path / prefix`, or with `.pt` when that file exists (`Path.with_suffix(".pt")`), or
/// when neither does, so the error on reading it names the file a model directory holds
/// (`cannot read /models/x/pileup.pt`, not `.../pileup`).
fn checkpoint(model_path: &Path, prefix: &str) -> PathBuf {
    let path = model_path.join(prefix);
    let with_pt = path.with_extension("pt");
    if !prefix.ends_with(".pt") && (with_pt.exists() || !path.exists()) { with_pt } else { path }
}

impl Settings {
    /// Resolves the arguments as the original's pipeline script does, warning where it does.
    ///
    /// # Errors
    ///
    /// Fails for options that are not ported and for values `CheckEnvs` rejects.
    pub(crate) fn resolve(args: &RunArgs) -> Result<Self> {
        Self::resolve_with(args, snv_infer::eight_bit_supported())
    }

    /// [`Self::resolve`], with whether this CPU has the 8-bit kernels `--approximate` uses given.
    fn resolve_with(args: &RunArgs, eight_bit_supported: bool) -> Result<Self> {
        let tools = [
            ("samtools", &args.samtools),
            ("python", &args.python),
            ("pypy", &args.pypy),
            ("parallel", &args.parallel),
            ("whatshap", &args.whatshap),
        ];
        for (name, _) in tools.iter().filter(|(_, path)| path.is_some()) {
            log::info!("--{name} is ignored: lrpipe-snv runs that step in process");
        }
        if let Some(vcf) = args.vcf_fn.as_ref().filter(|v| !v.is_file()) {
            bail!("VCF file {} provides but not found", vcf.display());
        }
        // The original stops after the pileup without a gVCF, and so did lrpipe-snv, without a
        // word, where the pileup had calls (lrpipe#54).
        if args.pileup_only && args.gvcf.on() {
            bail!(
                "--pileup-only stops before the gVCF is written: drop --gvcf, or drop --pileup-only to write \
                 merge_output.gvcf.gz"
            );
        }
        let ont = args.platform == Platform::Ont;
        let or_default = |value: f64, ont_default: f64, other: f64| {
            if value != 0.0 {
                value
            } else if ont {
                ont_default
            } else {
                other
            }
        };
        let device = match args.device {
            InferenceDevice::Cpu if args.use_gpu => InferenceDevice::Cuda(0),
            device => device,
        };
        if matches!(device, InferenceDevice::Cuda(_)) && !cfg!(feature = "gpu") {
            bail!("this build has no GPU calling (--device cuda:N needs the `gpu` feature)");
        }
        if device == InferenceDevice::Metal && !cfg!(all(feature = "metal", target_os = "macos")) {
            bail!("this build has no Metal calling (--device metal needs the `metal` feature, on macOS)");
        }
        // `--approximate` on a CPU without the 8-bit kernels would gain nothing in the networks,
        // so it changes nothing there: no fewer full-alignment calls, no 8-bit check.
        let approximate = if args.approximate && device == InferenceDevice::Cpu && !eight_bit_supported {
            log::warn!(
                "--approximate has no effect: this CPU has none of the 8-bit instructions it uses (AVX-VNNI or \
                 AVX-512 VNNI on x86-64, dotprod on aarch64), so the networks run in float32 and full alignment \
                 gets its usual share of the calls"
            );
            false
        } else {
            args.approximate
        };
        // On the CPU full alignment is most of the time; `--approximate` sends half of the ONT
        // variant calls there instead of 70%, which changes fewer than 0.5% of the final calls.
        let ont_var_pct_full = if approximate && device == InferenceDevice::Cpu { 0.5 } else { 0.7 };
        let var_pct_full = or_default(args.var_pct_full, ont_var_pct_full, 0.3);
        let ref_pct_full = or_default(args.ref_pct_full, 0.1, 0.3);
        let mut snp_min_af = or_default(args.snp_min_af, 0.08, 0.08);
        let mut indel_min_af = or_default(args.indel_min_af, 0.15, 0.08);
        if args.vcf_fn.is_some() {
            // Genotyping: every covered known site is a candidate.
            (snp_min_af, indel_min_af) = (0.0, 0.0);
        }
        let model_name = args.model_path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
        let dwell_time = dwell_time(args, model_name)?;
        let var_pct_phasing = if GUPPY5_MODELS.contains(&model_name) {
            0.8
        } else if args.var_pct_phasing != 0.0 {
            args.var_pct_phasing
        } else {
            0.7
        };
        let mut min_mq = args.min_mq;
        if min_mq < 5 {
            log::warn!("Invalid minimum mapping quality input --min_mq>=5");
            min_mq = 5;
        }
        let mut min_coverage = args.min_coverage;
        if min_coverage < 2 {
            log::warn!("Invalid minimum coverage input --min_coverage>=2");
            min_coverage = 2;
        }
        let mut min_contig_size = args.min_contig_size;
        if min_contig_size < 0 {
            log::warn!("Invalid minimum contig size --min_contig_size>=0");
            min_contig_size = 0;
        }
        for (name, value) in [
            ("var_pct_full", var_pct_full),
            ("ref_pct_full", ref_pct_full),
            ("snp_min_af", snp_min_af),
            ("indel_min_af", indel_min_af),
        ] {
            if !(0.0..=1.0).contains(&value) {
                bail!("parameter --{name}={value} out of range [0, 1]");
            }
        }
        if args.qual < 0 {
            bail!("parameter --qual={} out of range [0, inf)", args.qual);
        }
        if ref_pct_full > 0.3 {
            log::warn!("For efficiency, we use a maximum 30% reference candidates for full-alignment calling");
        }
        // `CreateTensorPileup`: fast mode only applies to ONT.
        let fast_mode = ont && args.fast_mode;
        let candidates = CandidateParams {
            min_depth: usize::try_from(if fast_mode { min_coverage.max(4) } else { min_coverage })?,
            // Python passes the floats through cffi as C `float`.
            min_snp_af: (if fast_mode { snp_min_af.max(0.15) } else { snp_min_af }) as f32,
            min_indel_af: indel_min_af as f32,
            min_mapq: u8::try_from(min_mq.min(255))?,
            max_indel_length: if args.enable_long_indel { 100_000 } else { 50 },
            alleles: if args.call_snp_only { CandidateAlleles::SnpsOnly } else { CandidateAlleles::SnpsAndIndels },
            flanks: if args.enable_variant_calling_at_sequence_head_and_tail {
                FlankRequirement::None
            } else {
                FlankRequirement::Contiguous
            },
        };
        let longphase_bin = || -> Result<PathBuf> {
            let bin = args.longphase.clone().unwrap_or_else(|| PathBuf::from("longphase"));
            let found = if bin.components().count() > 1 {
                bin.is_file().then_some(bin.clone())
            } else {
                std::env::var_os("PATH")
                    .and_then(|paths| std::env::split_paths(&paths).map(|dir| dir.join(&bin)).find(|p| p.is_file()))
            };
            found.with_context(|| format!("Cannot find LongPhase binary {}", bin.display()))
        };
        let haplotag = args.use_whatshap_for_final_output_haplotagging;
        let whatshap_final = args.use_whatshap_for_final_output_phasing
            || args.enable_phasing
            || (haplotag && !args.use_longphase_for_final_output_phasing);
        let final_phasing = if whatshap_final {
            Some(Phaser::WhatsHap)
        } else if args.use_longphase_for_final_output_phasing {
            Some(Phaser::LongPhase(longphase_bin()?))
        } else {
            None
        };
        let use_longphase = args.use_longphase_for_intermediate_phasing || args.longphase_for_phasing;
        if use_longphase && args.platform == Platform::Ilmn {
            log::warn!("Illumina platform do not support longphase phasing, will enable whatshap phasing!");
        }
        let intermediate_phasing = if use_longphase && args.platform != Platform::Ilmn {
            Phaser::LongPhase(longphase_bin()?)
        } else {
            Phaser::WhatsHap
        };
        let gender = match args.gender {
            GenderArg::Unknown => Gender::Unknown,
            GenderArg::Male => Gender::Male,
            GenderArg::Female => Gender::Female,
        };
        let mut par_regions = args.par_regions_bed.clone();
        if gender == Gender::Male {
            match &par_regions {
                Some(bed) if !bed.is_file() || bed.metadata().map(|m| m.len()).unwrap_or(0) == 0 => {
                    log::warn!(
                        "--par_regions_bed {} is missing or empty; proceeding WITHOUT PAR (chrX/chrY PAR regions will be \
                         treated as haploid)",
                        bed.display()
                    );
                    par_regions = None;
                }
                Some(_) => {}
                None => log::info!(
                    "--gender male without --par_regions_bed: chrX/chrY PAR regions will be treated as haploid. \
                     Consider providing --par_regions_bed for pseudo-autosomal regions."
                ),
            }
        }
        Ok(Self {
            var_pct_full,
            ref_pct_full,
            var_pct_phasing,
            candidates,
            min_af: (snp_min_af, indel_min_af),
            min_coverage,
            min_mq,
            min_contig_size,
            final_phasing,
            intermediate_phasing,
            haplotag,
            haploid: if args.haploid_precise {
                Haploid::Precise
            } else if args.haploid_sensitive {
                Haploid::Sensitive
            } else {
                Haploid::Off
            },
            gender,
            par_regions,
            ref_calls: if args.print_ref_calls { RefCalls::Show } else { RefCalls::Hide },
            iupac: match (args.keep_iupac_bases, args.ambiguous_bases) {
                (true, _) => IupacBases::Keep,
                (false, AmbiguousBases::Vcf) => IupacBases::Vcf,
                (false, AmbiguousBases::N) => IupacBases::ToN,
            },
            indel_limit: if args.enable_long_indel { IndelLimit::Long } else { IndelLimit::Standard },
            dwell_time,
            pileup_model: checkpoint(&args.model_path, &args.pileup_model_prefix),
            fa_model: checkpoint(&args.model_path, &args.fa_model_prefix),
            device,
            gpu_wait: args.gpu_wait,
            nvrtc: args.nvrtc.clone(),
            precision: if approximate { Precision::Mixed } else { Precision::Float32 },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn booleans_parse_like_str2bool() {
        for yes in ["yes", "True", "t", "Y", "1"] {
            assert_eq!(str2bool(yes), Ok(true));
        }
        for no in ["no", "FALSE", "f", "n", "0"] {
            assert_eq!(str2bool(no), Ok(false));
        }
        assert!(str2bool("maybe").is_err());
    }

    fn spelled(args: &[&str]) -> Vec<String> {
        let args = args.iter().map(OsString::from).collect();
        run_script_spellings(args).into_iter().map(|a| a.into_string().expect("UTF-8")).collect()
    }

    #[test]
    fn run_script_spellings_become_lrpipe_snv_options() {
        let line = [
            "lrpipe-snv",
            "run",
            "--bam_fn=/data/a_b.bam",
            "-f",
            "/ref/x_y.fa",
            "--chunk_size",
            "1000000",
            "--bed_fn=EMPTY",
            "--ctg_name",
            "EMPTY",
            "--gvcf=True",
            "--print_ref_calls",
            "False",
        ];
        let want = [
            "lrpipe-snv",
            "run",
            "--bam-fn=/data/a_b.bam",
            "-f",
            "/ref/x_y.fa",
            "--chunk-size",
            "1000000",
            "--gvcf=True",
            "--print-ref-calls",
            "False",
        ];
        assert_eq!(spelled(&line), want);
    }

    #[test]
    fn a_link_named_like_the_run_script_runs_lrpipe_snv_run() {
        assert_eq!(
            spelled(&["/opt/bin/run_clair3.sh", "--model_path=m"]),
            ["/opt/bin/run_clair3.sh", "run", "--read-order=input", "--ambiguous-bases=n", "--model-path=m"]
        );
        for option in ["--read_order=canonical", "--read-o=canonical"] {
            assert_eq!(
                spelled(&["/opt/bin/run_clair3.py", option]),
                ["/opt/bin/run_clair3.py", "run", "--ambiguous-bases=n", &option.replace('_', "-")]
            );
        }
        // Other subcommands keep their arguments.
        assert_eq!(spelled(&["lrpipe-snv", "phase", "--ctg_name", "x"]), ["lrpipe-snv", "phase", "--ctg_name", "x"]);
    }

    /// LSNV-004: `--approximate` keeps ONT's `--var_pct_full` at 0.7, and calls in float32, on a
    /// CPU without the 8-bit kernels; with them it sends half of the variant calls.
    #[test]
    fn approximate_without_eight_bit_kernels_changes_nothing() {
        use clap::Parser;
        #[derive(Debug, Parser)]
        struct Cli {
            #[command(flatten)]
            run: RunArgs,
        }
        let parse = |extra: &[&str]| {
            let mut line = vec!["lrpipe-snv", "-b", "in.bam", "-f", "ref.fa", "-m", "model", "-t", "4", "-o", "out"];
            line.extend_from_slice(extra);
            Cli::try_parse_from(line).expect("parses").run
        };
        let ont = parse(&["-p", "ont", "--approximate"]);
        let without = Settings::resolve_with(&ont, false).expect("resolves");
        assert_eq!((without.var_pct_full, without.precision), (0.7, Precision::Float32));
        let with = Settings::resolve_with(&ont, true).expect("resolves");
        assert_eq!((with.var_pct_full, with.precision), (0.5, Precision::Mixed));
        // Without `--approximate`, or with an explicit proportion, the flag does not matter.
        for supported in [false, true] {
            let exact = Settings::resolve_with(&parse(&["-p", "ont"]), supported).expect("resolves");
            assert_eq!((exact.var_pct_full, exact.precision), (0.7, Precision::Float32));
            let set = parse(&["-p", "ont", "--approximate", "--var-pct-full", "0.6"]);
            assert_eq!(Settings::resolve_with(&set, supported).expect("resolves").var_pct_full, 0.6);
        }
    }

    #[test]
    fn a_run_script_command_line_parses() {
        use clap::Parser;
        #[derive(Debug, Parser)]
        struct Cli {
            #[command(flatten)]
            run: RunArgs,
        }
        let line = spelled(&[
            "run_clair3.sh",
            "-b",
            "in.bam",
            "-f",
            "ref.fa",
            "-m",
            "model",
            "-t",
            "4",
            "-p",
            "hifi",
            "-o",
            "out",
            "--ctg_name=chr20",
            "--gvcf",
            "--enable_phasing=True",
            "--samtools=samtools",
            "--longphase=EMPTY",
            "--use_whatshap_for_intermediate_phasing",
            "True",
            "--chunk_num=-1",
        ]);
        let args = Cli::try_parse_from(line.iter().filter(|a| *a != "run")).expect("parses").run;
        assert_eq!(args.bam, Path::new("in.bam"));
        assert_eq!((args.threads, args.platform, args.chunk_num), (4, Platform::Hifi, -1));
        assert_eq!(args.ctg_name.as_deref(), Some("chr20"));
        assert!(args.gvcf.steps(Path::new("tmp"), None).is_some() && args.enable_phasing && args.longphase.is_none());
    }

    #[test]
    fn pileup_only_with_a_gvcf_is_refused() {
        // It exited 0 without the gVCF wherever the pileup had calls (lrpipe#54).
        use clap::Parser;
        #[derive(Debug, Parser)]
        struct Cli {
            #[command(flatten)]
            run: RunArgs,
        }
        let parse = |extra: &[&str]| {
            let mut line = vec!["lrpipe-snv", "-b", "in.bam", "-f", "ref.fa", "-m", "model", "-t", "4", "-o", "out"];
            line.extend_from_slice(&["-p", "hifi"]);
            line.extend_from_slice(extra);
            Cli::try_parse_from(line).expect("parses").run
        };
        let refused = Settings::resolve_with(&parse(&["--pileup-only", "--gvcf"]), false).expect_err("refused");
        assert!(refused.to_string().contains("--pileup-only stops before the gVCF"), "{refused}");
        for alone in [&["--pileup-only"][..], &["--gvcf"][..]] {
            Settings::resolve_with(&parse(alone), false).unwrap_or_else(|e| panic!("{alone:?}: {e}"));
        }
    }

    #[test]
    fn a_missing_checkpoint_is_named_with_its_pt() {
        // `-m /nope` failed with "cannot read /nope/pileup": the file a model holds is pileup.pt.
        assert_eq!(checkpoint(Path::new("/nope"), "pileup"), Path::new("/nope/pileup.pt"));
        let dir = std::env::temp_dir().join(format!("lrpipe-snv-checkpoint-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        // A file named as the prefix, without .pt, is the original's, and is taken.
        std::fs::write(dir.join("pileup"), "").expect("written");
        assert_eq!(checkpoint(&dir, "pileup"), dir.join("pileup"));
        std::fs::write(dir.join("pileup.pt"), "").expect("written");
        assert_eq!(checkpoint(&dir, "pileup"), dir.join("pileup.pt"));
        assert_eq!(checkpoint(&dir, "custom.pt"), dir.join("custom.pt"));
        let _ = std::fs::remove_dir_all(dir);
    }
}

//! The `lrpipe-tr` command: genotypes a repeat catalog from a BAM or CRAM into a VCF.

use std::ffi::OsString;
use std::path::PathBuf;
use std::time::Instant;

use crate::catalog;
use crate::genotype::{GenotypeOptions, Platform};
use crate::run::{Job, Karyotype, RunOptions, run, sort_loci, x_pars};
use anyhow::{Context, Result, bail};
use clap::Parser;
use pasta::fasta;

/// Genotype tandem repeats from long reads (HiFi or ONT).
#[derive(Parser)]
#[command(version, about)]
struct Args {
    /// Aligned reads: indexed BAM or CRAM
    #[arg(long)]
    reads: PathBuf,
    /// Reference FASTA with a .fai index (also decodes CRAM)
    #[arg(long, visible_alias = "reference")]
    genome: PathBuf,
    /// Repeat catalog: BED with TRGT-style `ID=..;MOTIFS=..` or a motif list in column 4
    #[arg(long)]
    repeats: PathBuf,
    /// Only the catalog's loci that overlap these regions (BED); each is genotyped from all its
    /// reads, as without (lrpipe-tr's own option)
    #[arg(long, value_name = "BED")]
    regions: Option<PathBuf>,
    /// Output VCF (.vcf.gz is bgzipped, with a .tbi index)
    #[arg(long)]
    output: PathBuf,
    /// Sample name in the VCF
    #[arg(long, default_value = "SAMPLE")]
    sample: String,
    /// Sex chromosomes: XX, or XY for one copy of chrX (outside the PARs) and chrY
    #[arg(long, value_enum, ignore_case = true, default_value = "xx")]
    karyotype: KaryotypeArg,
    /// Worker threads
    #[arg(short = 't', long, default_value_t = 8)]
    threads: usize,
    /// Lowest mapping quality of reads used
    #[arg(long, default_value_t = 1)]
    min_mapq: u8,
    /// Sequencing platform: sets the genotyping parameters below
    #[arg(long, value_enum, ignore_case = true, default_value = "hifi")]
    platform: PlatformArg,
    /// Per-read repeat length variance floor, in squared bases [default: by platform]
    #[arg(long, help_heading = "Genotyping")]
    noise_floor: Option<f64>,
    /// Per-read repeat length standard deviation per base of allele length [default: by platform]
    #[arg(long, help_heading = "Genotyping")]
    noise_per_base: Option<f64>,
    /// Penalty, in units of ln(reads), a two-allele length split must beat [default: by platform]
    #[arg(long, help_heading = "Genotyping")]
    split_penalty: Option<f64>,
    /// Net share of tagged reads that must prefer their own haplotype's allele [default: by platform]
    #[arg(long, help_heading = "Genotyping")]
    min_haplotag_support: Option<f64>,
    /// Share of the reads each allele of a heterozygous call needs [default: by platform]
    #[arg(long, help_heading = "Genotyping")]
    min_allele_share: Option<f64>,
    /// Share of the reads each of two exact-length peaks needs [default: by platform]
    #[arg(long, help_heading = "Genotyping")]
    min_peak_share: Option<f64>,
    /// The reads' haplotypes from a haplotag list (WhatsHap's `haplotag --output-haplotag-list`)
    /// in place of their HP tags, as `haplotag` would tag them: a listed read's primary
    /// alignment gets its haplotype, every other alignment none
    #[arg(long, value_name = "LIST")]
    haplotags: Option<PathBuf>,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum PlatformArg {
    Hifi,
    Ont,
}

#[derive(Clone, Copy, clap::ValueEnum)]
enum KaryotypeArg {
    Xx,
    Xy,
}

/// Refuses a reference whose contigs differ in length from the reads' header: the reads were
/// aligned to another assembly, so their positions do not match its loci.
fn check_reads_contigs(reads: &std::path::Path, genome: &std::path::Path, contigs: &[(String, u64)]) -> Result<()> {
    let reference = if lrpipe_reads::is_cram(reads) {
        Some(lrpipe_reads::Reference::open(genome).map_err(|e| e.into_io_error())?)
    } else {
        None
    };
    let file = lrpipe_reads::IndexedReads::open(reads, reference.as_ref())
        .map_err(|e| e.into_io_error())
        .with_context(|| format!("opening {}", reads.display()))?;
    for (name, sequence) in file.header().reference_sequences() {
        let len = sequence.length().get() as u64;
        if let Some(&(_, other)) = contigs.iter().find(|(c, _)| c.as_bytes() == name.as_slice())
            && other != len
        {
            bail!(
                "contig {name} is {len} bp in {} but {other} bp in the reference {}: give the reference the reads were \
                 aligned to",
                reads.display(),
                genome.display()
            );
        }
    }
    Ok(())
}

/// Runs the command line `args` (the program name first), as the `lrpipe-tr` binary does.
///
/// # Errors
///
/// If the arguments are invalid or the run fails.
pub fn main(args: Vec<OsString>) -> Result<()> {
    let a = Args::parse_from(args);
    let start = Instant::now();
    let mut loci = catalog::read(&a.repeats).with_context(|| format!("reading {}", a.repeats.display()))?;
    let fai_path = a.genome.with_added_extension("fai");
    let fai = fasta::fai::fs::read(&fai_path).with_context(|| format!("reading {}", fai_path.display()))?;
    let contigs: Vec<(String, u64)> = fai.as_ref().iter().map(|r| (r.name().to_string(), r.length())).collect();
    let total = loci.len();
    sort_loci(&mut loci, &contigs);
    if loci.len() < total {
        eprintln!("lrpipe-tr: {} loci on contigs not in the reference skipped", total - loci.len());
    }
    let length = |contig: &str| contigs.iter().find(|(c, _)| c == contig).map(|&(_, l)| l);
    if let Some(l) = loci.iter().find(|l| length(&l.contig).is_some_and(|len| l.end > len as i64)) {
        bail!(
            "the catalog's locus {} ({}:{}-{}) ends past {}'s end in the reference ({} bp): is the catalog for this \
             reference?",
            l.id,
            l.contig,
            l.start,
            l.end,
            l.contig,
            length(&l.contig).unwrap_or_default()
        );
    }
    check_reads_contigs(&a.reads, &a.genome, &contigs)?;
    if let Some(path) = &a.regions {
        let regions = lrpipe_reads::regions::Regions::read(path)?;
        let before = loci.len();
        catalog::retain_in(&mut loci, &regions);
        eprintln!("lrpipe-tr: {} of {before} loci overlap the regions of {}", loci.len(), path.display());
    }
    let platform = match a.platform {
        PlatformArg::Hifi => Platform::Hifi,
        PlatformArg::Ont => Platform::Ont,
    };
    let preset = GenotypeOptions::for_platform(platform);
    let genotype = GenotypeOptions {
        noise_floor: a.noise_floor.unwrap_or(preset.noise_floor),
        noise_per_base: a.noise_per_base.unwrap_or(preset.noise_per_base),
        split_penalty: a.split_penalty.unwrap_or(preset.split_penalty),
        min_haplotag_support: a.min_haplotag_support.unwrap_or(preset.min_haplotag_support),
        min_allele_share: a.min_allele_share.unwrap_or(preset.min_allele_share),
        min_peak_share: a.min_peak_share.unwrap_or(preset.min_peak_share),
        ..preset
    };
    let karyotype = match a.karyotype {
        KaryotypeArg::Xx => Karyotype::Xx,
        KaryotypeArg::Xy => Karyotype::Xy,
    };
    let opts = RunOptions {
        threads: a.threads,
        min_mapq: a.min_mapq,
        karyotype,
        genotype,
        x_pars: x_pars(&a.genome, karyotype)?,
        ..RunOptions::default()
    };
    let haplotags = match &a.haplotags {
        Some(path) => {
            let list = lrpipe_reads::haplotags::Haplotags::read(path)?;
            if list.ambiguous() > 0 {
                eprintln!(
                    "lrpipe-tr: warning: {} reads listed more than once on a contig with different haplotypes (a \
                     list written with --tag-supplementary?) are left untagged",
                    list.ambiguous()
                );
            }
            Some(std::sync::Arc::new(list))
        }
        None => None,
    };
    let job = Job { reads: a.reads, genome: a.genome, output: a.output, sample: a.sample, haplotags };
    let s = run(&job, &loci, &opts).context("genotyping")?;
    // An alignment is counted by each chunk of loci it reaches.
    if let Some(report) = job.haplotags.as_ref().zip(a.haplotags.as_ref()).and_then(|(l, p)| l.matches().report(p)) {
        eprintln!("lrpipe-tr: {report}");
    }
    eprintln!(
        "lrpipe-tr: {} of {} loci called in {:.1} s; {}",
        s.called,
        s.loci,
        start.elapsed().as_secs_f64(),
        job.output.display()
    );
    Ok(())
}

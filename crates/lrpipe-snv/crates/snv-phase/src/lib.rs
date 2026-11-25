//! The original's intermediate phasing step, run in-process with lrpipe-phase.
//!
//! The original phases the selected het SNPs of each contig with
//! `whatshap phase --output OUT --reference REF --chromosome CTG --distrust-genotypes
//! --ignore-read-groups VCF BAM`. [`phase_contig`] runs the same phasing through lrpipe-phase,
//! whose output matches WhatsHap 2.8 byte for byte apart from the `##commandline` header.
//!
//! With `--use_whatshap_for_final_output_phasing` the original also phases each contig's final
//! calls ([`phase_output`]), and with `--use_whatshap_for_final_output_haplotagging` it tags the
//! reads by haplotype ([`haplotag_output`]).

use std::path::{Path, PathBuf};

use lrpipe_phase::haplotag::{HaplotagOptions, SupplementaryStrategy, run_haplotag};
use lrpipe_phase::phase::{Algorithm, PhaseOptions, run_phase};
use lrpipe_phase::vcf::PhaseTag;

/// Errors from phasing.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A path is not valid UTF-8 (lrpipe-phase takes paths as strings).
    #[error("unsupported path {path:?}")]
    Path {
        /// The path.
        path: PathBuf,
    },
    /// Haplotagging or indexing its output failed.
    #[error("haplotagging {path:?} failed")]
    Haplotag {
        /// The output BAM.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
    /// lrpipe-phase failed.
    #[error("phasing {contig} failed")]
    Phase {
        /// The contig being phased.
        contig: String,
        /// The underlying error.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

/// One `whatshap phase` job of the original's pipeline.
#[derive(Debug, Clone)]
pub struct PhaseJob<'a> {
    /// Het SNPs to phase (`<phase_vcf_dir>/<contig>.vcf` from SelectHetSnp).
    pub variants: &'a Path,
    /// Sorted, indexed BAM/CRAM.
    pub alignments: &'a Path,
    /// Indexed reference FASTA.
    pub reference: &'a Path,
    /// Contig to phase (`--chromosome`).
    pub contig: &'a str,
    /// Output VCF (`phased_<contig>.vcf.gz`).
    pub output: &'a Path,
    /// lrpipe-phase threads; the output does not depend on it.
    pub threads: usize,
}

fn path_str(path: &Path) -> Result<String, Error> {
    path.to_str().map(str::to_owned).ok_or_else(|| Error::Path { path: path.to_owned() })
}

/// The arguments of `whatshap phase` that the original sets, with WhatsHap 2.8's defaults for
/// the rest.
struct PhaseCall {
    variants: String,
    alignments: String,
    reference: String,
    output: String,
    chromosomes: Vec<String>,
    distrust_genotypes: bool,
    threads: usize,
    command_line: String,
}

impl PhaseCall {
    fn run(self) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        let options = PhaseOptions {
            variant_file: self.variants,
            phase_input_files: vec![self.alignments],
            reference: Some(self.reference),
            no_reference: false,
            output: self.output,
            samples: Vec::new(),
            chromosomes: self.chromosomes,
            excluded_chromosomes: Vec::new(),
            ignore_read_groups: true,
            only_snvs: false,
            mapping_quality: 20,
            max_coverage: 15,
            tag: PhaseTag::PS,
            read_list_filename: None,
            use_supplementary: false,
            supplementary_distance_threshold: 100_000,
            command_line: Some(self.command_line),
            threads: self.threads.max(1),
            dp_memory_mb: 1024,
            max_memory_mb: None,
            read_merging: false,
            algorithm: Algorithm::WhatsHap,
            ped: None,
            use_ped_samples: false,
            recombrate: 1.26,
            genmap: None,
            genetic_haplotyping: true,
            recombination_list_filename: None,
            distrust_genotypes: self.distrust_genotypes,
            include_homozygous: false,
            default_gq: 30,
            gl_regularizer: None,
            gtchange_list_filename: None,
        };
        Ok(run_phase(options)?)
    }
}

/// Phases one contig like the original's `whatshap phase` call.
///
/// # Errors
///
/// Fails if a path is not UTF-8 or lrpipe-phase reports an error.
pub fn phase_contig(job: &PhaseJob<'_>) -> Result<(), Error> {
    let (variants, alignments) = (path_str(job.variants)?, path_str(job.alignments)?);
    let (reference, output) = (path_str(job.reference)?, path_str(job.output)?);
    let command_line = format!(
        "(lrpipe-snv {} via lrpipe-phase) phase --output {output} --reference {reference} --chromosome {} \
         --distrust-genotypes --ignore-read-groups {variants} {alignments}",
        env!("CARGO_PKG_VERSION"),
        job.contig,
    );
    PhaseCall {
        variants,
        alignments,
        reference,
        output,
        chromosomes: vec![job.contig.to_owned()],
        distrust_genotypes: true,
        threads: job.threads,
        command_line,
    }
    .run()
    .map_err(|source| Error::Phase { contig: job.contig.to_owned(), source })
}

/// Phases one contig's final calls like `--use_whatshap_for_final_output_phasing`:
/// `whatshap phase --output OUT --reference REF --ignore-read-groups merge_<contig>.vcf BAM`,
/// which trusts the genotypes and phases every contig in the VCF.
///
/// # Errors
///
/// Fails if a path is not UTF-8 or lrpipe-phase reports an error.
pub fn phase_output(job: &PhaseJob<'_>) -> Result<(), Error> {
    let (variants, alignments) = (path_str(job.variants)?, path_str(job.alignments)?);
    let (reference, output) = (path_str(job.reference)?, path_str(job.output)?);
    let command_line = format!(
        "(lrpipe-snv {} via lrpipe-phase) phase --output {output} --reference {reference} --ignore-read-groups \
         {variants} {alignments}",
        env!("CARGO_PKG_VERSION"),
    );
    PhaseCall {
        variants,
        alignments,
        reference,
        output,
        chromosomes: Vec::new(),
        distrust_genotypes: false,
        threads: job.threads,
        command_line,
    }
    .run()
    .map_err(|source| Error::Phase { contig: job.contig.to_owned(), source })
}

/// Tags reads by haplotype like `--use_whatshap_for_final_output_haplotagging`: `whatshap
/// haplotag --output OUT --reference REF --ignore-read-groups phased_merge_output.vcf.gz BAM`
/// with WhatsHap 2.8's defaults, then `samtools index` (a BAI next to `output`).
///
/// # Errors
///
/// Fails if a path is not UTF-8, or lrpipe-phase or the indexing reports an error.
pub fn haplotag_output(
    phased_vcf: &Path,
    alignments: &Path,
    reference: &Path,
    output: &Path,
    threads: usize,
) -> Result<(), Error> {
    haplotag(phased_vcf, alignments, reference, None, output, threads)
}

/// Tags one contig's reads like the original's Illumina step 4/7: `whatshap haplotag --output
/// phase_bam/<contig>.bam --reference REF --ignore-read-groups --regions CONTIG
/// phased_<contig>.vcf.gz BAM`, then `samtools index`.
///
/// # Errors
///
/// Fails if a path is not UTF-8, or lrpipe-phase or the indexing reports an error.
pub fn haplotag_contig(
    phased_vcf: &Path,
    alignments: &Path,
    reference: &Path,
    contig: &str,
    output: &Path,
    threads: usize,
) -> Result<(), Error> {
    haplotag(phased_vcf, alignments, reference, Some(contig), output, threads)
}

fn haplotag(
    phased_vcf: &Path,
    alignments: &Path,
    reference: &Path,
    region: Option<&str>,
    output: &Path,
    threads: usize,
) -> Result<(), Error> {
    let (variant_file, alignment_file) = (path_str(phased_vcf)?, path_str(alignments)?);
    let (reference, out) = (path_str(reference)?, path_str(output)?);
    let regions = region.map_or_else(String::new, |r| format!(" --regions {r}"));
    let command_line = format!(
        "(lrpipe-snv {} via lrpipe-phase) haplotag --output {out} --reference {reference} --ignore-read-groups\
         {regions} {variant_file} {alignment_file}",
        env!("CARGO_PKG_VERSION"),
    );
    let failed = |source: Box<dyn std::error::Error + Send + Sync>| Error::Haplotag { path: output.to_owned(), source };
    run_haplotag(HaplotagOptions {
        variant_file,
        alignment_file,
        output: Some(out),
        reference: Some(reference),
        no_reference: false,
        regions: region.map(|r| vec![r.to_owned()]),
        regions_bed: None,
        ignore_linked_read: false,
        linked_read_distance_cutoff: 50_000,
        ignore_read_groups: true,
        given_samples: None,
        haplotag_list: None,
        supplementary_strategy: SupplementaryStrategy::Skip,
        supplementary_distance_threshold: 100_000,
        supplementary_strand_match: true,
        ploidy: 2,
        skip_missing_contigs: false,
        untag_missing_contigs: false,
        list_only: false,
        write_index: false,
        output_threads: None,
        compression_level: None,
        threads: threads.max(1),
        command_line,
    })
    .map_err(|e| failed(e.into()))?;
    snv_hts::index_bam(output).map_err(|e| failed(e.into()))
}

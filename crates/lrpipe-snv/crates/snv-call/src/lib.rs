//! Variant output: genotypes, alleles and VCF lines from network probabilities.
//!
//! A port of `CallVariants.py` (`output_with`, `output_from`, `compute_PL`, `get_header`) from
//! version 2.0.3 of the original, producing the same bytes: [`vcf_record`] for models without
//! indel-length outputs (pileup) and [`vcf_record_with_indel_lengths`] for models with them (full
//! alignment), both without long-indel calling (`enable_long_indel=False`).
//!
//! For gVCF output (`--gvcf`), records carry PL ([`Likelihoods`]), [`NonVariantWriter`] writes
//! a pileup chunk's non-variant blocks (`variantInfoCalculator` in `preprocess/utils.py`), and
//! [`Lz4Writer`] and [`open_lz4`] handle the lz4 frames the original keeps its intermediate gVCFs
//! in.

mod alt_info;
mod decode;
mod gvcf;
mod lz4;
mod record;
mod variant_length;

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

pub use decode::{GENOTYPE_CLASSES, GT21_CLASSES, VARIANT_LENGTH_CLASSES};
pub use gvcf::{NonVariantOptions, NonVariantWriter, Resolution, chunk_gvcf_name};
pub use lz4::{Lz4Writer, open_lz4};
pub use record::{
    CallSource, Haploid, IupacBases, Likelihoods, OutputConfig, RefCalls, vcf_record, vcf_record_with_indel_lengths,
};

/// The version of the original this crate ports, written to the header's version line.
pub const PORTED_VERSION: &str = "2.0.3";

/// Errors from building VCF output.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A candidate's position, alt info or probabilities cannot be decoded.
    #[error("cannot decode candidate {position:?}")]
    Candidate {
        /// The candidate's position string.
        position: String,
    },
    /// A file could not be read.
    #[error("cannot read {path:?}")]
    Io {
        /// The file.
        path: PathBuf,
        /// The underlying error.
        #[source]
        source: std::io::Error,
    },
    /// Output could not be written.
    #[error("cannot write the output")]
    Write(#[source] std::io::Error),
    /// `--gq_bin_size` is 0, which the original divides by.
    #[error("the GQ bin size must not be 0")]
    GqBinSize,
}

/// `get_header`'s template up to the INFO lines; [`FORMAT_START`] and [`FORMAT_END`] follow, with
/// the gVCF template's extra lines between the parts.
const HEADER_START: &str = r#"##fileformat=VCFv4.2
##source=Clair3
##clair3_version={version}
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##FILTER=<ID=RefCall,Description="Reference call">
##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">"#;

/// The gVCF template's extra lines after the INFO lines.
const GVCF_ALT_END: &str = r#"##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">"#;

/// The FORMAT lines before the gVCF template's `MIN_DP`.
const FORMAT_START: &str = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#;

/// The gVCF template's `MIN_DP` line.
const GVCF_MIN_DP: &str =
    r#"##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">"#;

/// The FORMAT lines after the gVCF template's `MIN_DP`.
const FORMAT_END: &str = r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">"#;

/// Which template `get_header` fills.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Template {
    Vcf,
    Gvcf,
}

/// `get_header`: the template, `##reference` and `##cmdline` inserted at line 3, then the
/// `##contig` lines and the `#CHROM` line from the reference's `.fai` (without a trailing
/// newline), and the contig lengths.
fn header(
    template: Template,
    reference: Option<&Path>,
    cmd_file: Option<&Path>,
    sample: &str,
) -> Result<(String, HashMap<String, i64>), Error> {
    let mut parts = vec![HEADER_START];
    if template == Template::Gvcf {
        parts.push(GVCF_ALT_END);
    }
    parts.push(FORMAT_START);
    if template == Template::Gvcf {
        parts.push(GVCF_MIN_DP);
    }
    parts.push(FORMAT_END);
    let text = parts.join("\n").replace("{version}", PORTED_VERSION);
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    if let Some(r) = reference.filter(|r| r.exists()) {
        lines.insert(3, format!("##reference={}", r.display()));
    }
    if let Some(cmd) = cmd_file.and_then(|c| fs::read_to_string(c).ok()) {
        let cmd = cmd.trim_end();
        if !cmd.is_empty() {
            lines.insert(3, format!("##cmdline={cmd}"));
        }
    }
    let mut header = lines.join("\n") + "\n";
    let mut lengths = HashMap::new();
    if let Some(r) = reference {
        let mut fai = r.as_os_str().to_owned();
        fai.push(".fai");
        let fai = PathBuf::from(fai);
        let text = fs::read_to_string(&fai).map_err(|source| Error::Io { path: fai.clone(), source })?;
        for row in text.lines() {
            let mut cols = row.trim().split('\t');
            let (Some(name), Some(len)) = (cols.next(), cols.next()) else { continue };
            writeln!(header, "##contig=<ID={name},length={len}>").expect("writing to a String cannot fail");
            if let Ok(len) = len.parse() {
                lengths.insert(name.to_owned(), len);
            }
        }
        write!(header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")
            .expect("writing to a String cannot fail");
    }
    Ok((header, lengths))
}

/// The VCF header `get_header` writes (non-gVCF), ending with the `#CHROM` line and a newline
/// as `CallVariantsFromCffi` adds it.
///
/// `reference` adds `##reference` (if the file exists) and the `##contig` lines from its `.fai`;
/// `cmd_file` is the original's `--cmd_fn`, whose contents become `##cmdline` when non-empty.
///
/// # Errors
///
/// Fails if the reference's `.fai` cannot be read.
pub fn vcf_header(reference: Option<&Path>, cmd_file: Option<&Path>, sample: &str) -> Result<String, Error> {
    let (mut header, _) = header(Template::Vcf, reference, cmd_file, sample)?;
    // Without a reference there is no #CHROM line, and the caller's newline leaves a blank one.
    header.push('\n');
    Ok(header)
}

/// The header `variantInfoCalculator` writes at the top of a chunk's non-variant gVCF
/// (`get_header(gvcf=True)` and a newline), and the contig lengths from the reference's `.fai`.
///
/// # Errors
///
/// Fails if the reference's `.fai` cannot be read.
pub fn gvcf_header(reference: &Path, sample: &str) -> Result<(String, HashMap<String, i64>), Error> {
    let (mut header, lengths) = header(Template::Gvcf, Some(reference), None, sample)?;
    header.push('\n');
    Ok((header, lengths))
}

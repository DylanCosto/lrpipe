//! One pileup job: candidates, model inputs, probabilities and the chunk VCF.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use rayon::prelude::*;
use snv_call::{
    CallSource, IupacBases, Likelihoods, OutputConfig, vcf_header, vcf_record, vcf_record_with_indel_lengths,
};
use snv_infer::{FULL_ALIGNMENT_OUTPUTS, PILEUP_OUTPUTS, PileupModel};
use snv_pileup::{
    Candidate, CandidateAlleles, CandidateParams, FlankRequirement, TENSOR_LEN, alt_info_depth, scale_for_depth,
};

/// The original's candidate thresholds for one pileup job with the pipeline's defaults.
pub(crate) fn candidate_params(snp_min_af: f64, indel_min_af: f64, min_coverage: usize, min_mq: u8) -> CandidateParams {
    CandidateParams {
        min_depth: min_coverage,
        // Python passes the floats through cffi as C `float`.
        min_snp_af: snp_min_af as f32,
        min_indel_af: indel_min_af as f32,
        min_mapq: min_mq,
        // `param.maximum_variant_length_that_need_infer` (long indels off).
        max_indel_length: 50,
        alleles: CandidateAlleles::SnpsAndIndels,
        flanks: FlankRequirement::Contiguous,
    }
}

/// The tensors of all candidates, with `CallVariantsFromCffi`'s high-depth scaling applied.
pub(crate) fn model_inputs(candidates: &[Candidate]) -> Result<Vec<i32>> {
    let mut inputs: Vec<i32> = candidates.iter().flat_map(|c| c.tensor.iter().copied()).collect();
    for (c, tensor) in candidates.iter().zip(inputs.chunks_mut(TENSOR_LEN)) {
        let depth = alt_info_depth(&c.alt_info).context("alt info without a depth")?;
        scale_for_depth(tensor, depth);
    }
    Ok(inputs)
}

/// Scores `inputs` in batches of `batch_size` candidates, as the original batches them.
pub(crate) fn score(model: &mut impl PileupModel, inputs: &[i32], batch_size: usize) -> Result<Vec<f32>> {
    let mut probs = Vec::with_capacity(inputs.len() / TENSOR_LEN * PILEUP_OUTPUTS);
    for batch in inputs.chunks(batch_size * TENSOR_LEN) {
        probs.extend(model.predict(batch)?);
    }
    Ok(probs)
}

/// Header and records of a chunk VCF, as `call_variants_from_cffi` writes them.
#[derive(Debug, Clone)]
pub(crate) struct VcfOutput<'a> {
    pub(crate) reference: &'a Path,
    pub(crate) cmd_file: Option<&'a Path>,
    pub(crate) sample_name: &'a str,
    /// `--qual`; the original's pipeline leaves it at 2 for these jobs.
    pub(crate) qual: i64,
    /// `--keep_iupac_bases`.
    pub(crate) iupac: IupacBases,
    /// `maximum_variant_length_that_need_infer` (`--enable_long_indel`).
    pub(crate) max_indel_length: usize,
    /// PL in every record (`--gvcf`).
    pub(crate) likelihoods: Likelihoods,
}

/// One pileup job's candidates and probabilities.
pub(crate) type PileupCalls<'a> = (&'a [Candidate], &'a [f32]);

/// One full-alignment job's positions, alt infos and probabilities.
pub(crate) type FullAlignmentCalls<'a> = (&'a [String], &'a [String], &'a [f32]);

impl VcfOutput<'_> {
    fn config(&self, source: CallSource) -> OutputConfig {
        OutputConfig {
            min_pass_quality: Some(self.qual as f64),
            iupac: self.iupac,
            max_indel_length: self.max_indel_length,
            likelihoods: self.likelihoods,
            source,
            ..OutputConfig::default()
        }
    }

    fn create(&self, path: &Path) -> Result<BufWriter<File>> {
        let mut out = BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
        out.write_all(vcf_header(Some(self.reference), self.cmd_file, self.sample_name)?.as_bytes())?;
        Ok(out)
    }

    /// One job's VCF; like `CallVariantsFromCffi`, a VCF without records is removed.
    pub(crate) fn write(&self, path: &Path, candidates: &[Candidate], probs: &[f32]) -> Result<()> {
        self.write_each(path, &[(candidates, probs)])
    }

    /// One VCF of several jobs' records in turn, removed if it has none.
    pub(crate) fn write_each(&self, path: &Path, jobs: &[PileupCalls<'_>]) -> Result<()> {
        let records = self.write_jobs(path, jobs)?;
        remove_if_empty(path, records)
    }

    /// The record `vcf_record` makes of each candidate, `None` where it makes none.
    pub(crate) fn pileup_records(&self, candidates: &[Candidate], probs: &[f32]) -> Result<Vec<Option<String>>> {
        let config = self.config(CallSource::Pileup);
        let records = candidates.iter().zip(probs.chunks(PILEUP_OUTPUTS));
        Ok(records.map(|(c, p)| vcf_record(&c.position, &c.alt_info, p, &config)).collect::<Result<_, _>>()?)
    }

    /// One header, then the records of each job in turn, as the original's GPU calling writes all
    /// of a GPU's jobs into one VCF. Returns the number of records.
    pub(crate) fn write_jobs(&self, path: &Path, jobs: &[PileupCalls<'_>]) -> Result<u64> {
        let config = self.config(CallSource::Pileup);
        let mut out = self.create(path)?;
        let records = write_records(&mut out, jobs, |(candidates, probs), text| {
            let mut records = 0;
            for (c, p) in candidates.iter().zip(probs.chunks(PILEUP_OUTPUTS)) {
                if let Some(line) = vcf_record(&c.position, &c.alt_info, p, &config)? {
                    text.push_str(&line);
                    records += 1;
                }
            }
            Ok(records)
        })?;
        out.flush()?;
        Ok(records)
    }

    /// The same for a full-alignment job: INFO `F` and indel-length heads.
    pub(crate) fn write_full_alignment(
        &self,
        path: &Path,
        positions: &[String],
        alt_info: &[String],
        probs: &[f32],
    ) -> Result<()> {
        let records = self.write_full_alignment_jobs(path, &[(positions, alt_info, probs)])?;
        remove_if_empty(path, records)
    }

    /// One header, then the records of each full-alignment job in turn. Returns the number of
    /// records.
    pub(crate) fn write_full_alignment_jobs(&self, path: &Path, jobs: &[FullAlignmentCalls<'_>]) -> Result<u64> {
        let config = self.config(CallSource::FullAlignment);
        let mut out = self.create(path)?;
        let records = write_records(&mut out, jobs, |(positions, alt_info, probs), text| {
            let mut records = 0;
            for ((position, alt), p) in positions.iter().zip(*alt_info).zip(probs.chunks(FULL_ALIGNMENT_OUTPUTS)) {
                if let Some(line) = vcf_record_with_indel_lengths(position, alt, p, &config)? {
                    text.push_str(&line);
                    records += 1;
                }
            }
            Ok(records)
        })?;
        out.flush()?;
        Ok(records)
    }
}

/// Jobs formatted at once by [`write_records`].
const JOBS_PER_GROUP: usize = 64;

/// Writes the records `format` appends for each job, in job order. A genome's calls take a
/// minute to format on one thread, so groups of jobs are formatted in parallel.
fn write_records<J: Sync>(
    out: &mut impl Write,
    jobs: &[J],
    format: impl Fn(&J, &mut String) -> Result<u64> + Sync,
) -> Result<u64> {
    let mut records = 0;
    for group in jobs.chunks(JOBS_PER_GROUP) {
        let texts = group
            .par_iter()
            .map(|job| {
                let mut text = String::new();
                let n = format(job, &mut text)?;
                Ok((text, n))
            })
            .collect::<Result<Vec<_>>>()?;
        for (text, n) in texts {
            out.write_all(text.as_bytes())?;
            records += n;
        }
    }
    Ok(records)
}

/// `call_variants_from_cffi` removes a job's VCF when it holds no record.
fn remove_if_empty(path: &Path, records: u64) -> Result<()> {
    if records == 0 {
        log::info!("[INFO] No vcf output for file {}, remove empty file", path.display());
        std::fs::remove_file(path).with_context(|| format!("removing {}", path.display()))?;
    }
    Ok(())
}

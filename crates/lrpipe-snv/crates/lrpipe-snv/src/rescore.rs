//! `--approximate` pileup calling on the CPU: every candidate is scored with the network's 8-bit
//! products (which already scores close calls again in `float32`), then those whose QUAL lies
//! near a cut-off it is compared with are scored again in `float32` too.
//!
//! On ONT's model the 8-bit network moves a pileup QUAL by 0.08 at the median and at most 1.4
//! (HG002 ONT chr20, 326k records), but the cut-offs that pick het SNPs for phasing and candidates for full
//! alignment are quantiles with thousands of calls near them. Moving a few het SNPs across the
//! phasing cut-off re-phases whole blocks and so changes the haplotags of many full-alignment
//! tensors. Rescoring the candidates near a cut-off puts each on the side its `float32` QUAL does.
//!
//! How far 8 bits move a model's outputs depends on the model: on HiFi's, the recurrent products
//! moved probabilities by up to 0.96. So a run first checks the 8-bit network on candidates of
//! its own ([`check_eight_bit`]) and calls in `float32` when it moves them more than the
//! rescoring allows for.

use anyhow::Result;
use rayon::prelude::*;
use snv_infer::{PILEUP_OUTPUTS, PileupModel};
use snv_pileup::Candidate;

use crate::chunk::VcfOutput;
use crate::settings::Settings;

/// Rescored when the 8-bit QUAL is within this of a cut-off whose sides call differently: the
/// `--qual` thresholds and the phasing cut-off (moving het SNPs across it re-phases blocks).
/// [`check_eight_bit`] allows QUAL changes up to [`CHECK_QUAL`] on a sample; over all of HG002
/// ONT chr20 the largest was 15% more than on the check's sample (1.4 and 1.23).
const QUAL_WINDOW: f64 = 2.0;

/// Rescored when the 8-bit QUAL is within this of a full-alignment cut-off, which only decides
/// whether full alignment or pileup calls the candidate. On HG002 ONT chr20, a window of 0.75
/// (which 99.9% of QUAL changes are within) agreed with exact calls as well as rescoring every
/// candidate within 1.5 did (99.21% and 99.19% of PASS calls, F1 the same to three decimals), in
/// 56 s instead of 59.
const FULL_ALIGNMENT_WINDOW: f64 = 0.75;

/// Rescored when the 8-bit QUAL is below this: the call itself is close. On HG002 ONT chr20 every
/// call the 8-bit network changed (155 of 326k, before rescoring) had a QUAL below 5.4, so with
/// QUAL changes of up to [`QUAL_WINDOW`] the 8-bit QUAL is below 7.4.
const LOW_QUAL: f64 = 7.5;

/// The most the 8-bit network may move a probability on the checked candidates: half the class
/// gap the network itself rescores at (0.16), so no call it leaves alone can change.
const CHECK_PROBABILITY: f32 = 0.08;

/// The most the 8-bit network may move a QUAL on the checked candidates, inside [`QUAL_WINDOW`]
/// with room for the larger changes elsewhere in the run.
const CHECK_QUAL: f64 = 1.5;

/// Whether the 8-bit network (`mixed`) keeps a run's calls within what the rescoring catches,
/// judged on `candidates` (with their model inputs `inputs`) against `exact`: `Err` says how far it
/// moved them.
pub(crate) fn check_eight_bit(
    candidates: &[Candidate],
    inputs: &[i32],
    output: &VcfOutput<'_>,
    mixed: &mut impl PileupModel,
    exact: &mut impl PileupModel,
    batch: usize,
) -> Result<std::result::Result<(), String>> {
    let score = |model: &mut dyn FnMut(&[i32]) -> Result<Vec<f32>>| -> Result<Vec<f32>> {
        let mut probs = Vec::with_capacity(candidates.len() * PILEUP_OUTPUTS);
        for b in inputs.chunks(batch * snv_pileup::TENSOR_LEN) {
            probs.extend(model(b)?);
        }
        Ok(probs)
    };
    let fast = score(&mut |b| Ok(mixed.predict(b)?))?;
    let want = score(&mut |b| Ok(exact.predict(b)?))?;
    let probability = fast.iter().zip(&want).fold(0f32, |m, (a, b)| m.max((a - b).abs()));
    let quals = |probs: &[f32]| -> Result<Vec<Option<f64>>> {
        let lines = output.pileup_records(candidates, probs)?;
        Ok(lines.iter().map(|l| l.as_deref().and_then(parse_record).map(|(q, _)| q)).collect())
    };
    let qual = quals(&fast)?
        .iter()
        .zip(quals(&want)?)
        .filter_map(|(a, b)| Some((a.as_ref()? - b?).abs()))
        .fold(0f64, f64::max);
    log::info!(
        "the 8-bit pileup network moved {} candidates' probabilities by up to {probability:.3} and QUALs by up to {qual:.2}",
        candidates.len()
    );
    if probability > CHECK_PROBABILITY || qual > CHECK_QUAL {
        return Ok(Err(format!(
            "on {} of its candidates it moved a probability by up to {probability:.3} and a QUAL by up to {qual:.2} \
             (allowed: {CHECK_PROBABILITY} and {CHECK_QUAL})",
            candidates.len()
        )));
    }
    Ok(Ok(()))
}

/// A candidate's model inputs, kept small until the cut-offs are known: each run of zeros (most
/// of a tensor) as a zero byte and its length, every other value as a zigzag varint, which
/// never starts with a zero byte.
#[derive(Debug, Clone)]
pub(crate) struct PackedInputs(Vec<u8>);

impl PackedInputs {
    pub(crate) fn new(values: &[i32]) -> Self {
        let mut bytes = Vec::new();
        let mut zeros = 0u8;
        for &v in values {
            if v == 0 {
                if zeros == u8::MAX {
                    bytes.extend([0, zeros]);
                    zeros = 0;
                }
                zeros += 1;
                continue;
            }
            if zeros > 0 {
                bytes.extend([0, zeros]);
                zeros = 0;
            }
            let mut z = ((v << 1) ^ (v >> 31)) as u32;
            while z >= 0x80 {
                bytes.push((z as u8 & 0x7f) | 0x80);
                z >>= 7;
            }
            bytes.push(z as u8);
        }
        if zeros > 0 {
            bytes.extend([0, zeros]);
        }
        Self(bytes)
    }

    fn unpack_into(&self, out: &mut Vec<i32>) {
        let mut bytes = self.0.iter();
        while let Some(&b) = bytes.next() {
            if b == 0 {
                let run = *bytes.next().expect("a zero byte is followed by its run length");
                out.extend(std::iter::repeat_n(0, usize::from(run)));
                continue;
            }
            let (mut z, mut shift, mut byte) = (0u32, 0, b);
            loop {
                z |= u32::from(byte & 0x7f) << shift;
                if byte < 0x80 {
                    break;
                }
                shift += 7;
                byte = *bytes.next().expect("a varint ends");
            }
            out.push((z >> 1) as i32 ^ -((z & 1) as i32));
        }
    }
}

/// One pileup chunk as `--approximate` CPU calling keeps it until its VCF is written.
#[derive(Debug)]
pub(crate) struct ScoredChunk {
    pub(crate) contig: String,
    pub(crate) id: i64,
    /// Without their tensors, which `packed` holds as model inputs.
    pub(crate) candidates: Vec<Candidate>,
    pub(crate) probs: Vec<f32>,
    pub(crate) packed: Vec<PackedInputs>,
}

/// What a pileup record is compared with: reference calls with the reference cut-off, het SNPs
/// with the phasing and variant cut-offs, other variants with the variant cut-off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Reference,
    HetSnp,
    Variant,
}

/// QUAL and kind of a record written by `vcf_record`.
fn parse_record(line: &str) -> Option<(f64, Kind)> {
    let columns: Vec<&str> = line.splitn(11, '\t').collect();
    let qual = columns.get(5)?.parse().ok()?;
    let genotype = columns.get(9)?.split(':').next()?;
    let snp = columns.get(3)?.chars().count() == 1 && columns.get(4)?.chars().count() == 1;
    let kind = match genotype {
        "0/0" => Kind::Reference,
        "0/1" if snp => Kind::HetSnp,
        _ => Kind::Variant,
    };
    Some((qual, kind))
}

/// The highest of the lowest `pct` of `quals`, as `SelectQual` picks a cut-off.
fn lowest_cutoff(mut quals: Vec<f64>, pct: f64) -> Option<f64> {
    quals.sort_by(f64::total_cmp);
    let len = (pct * quals.len() as f64) as usize;
    len.checked_sub(1).map(|last| quals[last])
}

/// Rescores in `float32` (with networks from `load`, one per thread) every candidate of `chunks`
/// whose QUAL is near a cut-off, `batch` candidates at a time, and returns how many it rescored.
pub(crate) fn rescore_near_cutoffs<M: PileupModel>(
    pool: &rayon::ThreadPool,
    chunks: &mut [ScoredChunk],
    output: &VcfOutput<'_>,
    settings: &Settings,
    final_qual: f64,
    load: impl Fn() -> Result<M> + Sync,
    batch: usize,
) -> Result<usize> {
    let records = pool.install(|| {
        chunks
            .par_iter()
            .map(|c| {
                let lines = output.pileup_records(&c.candidates, &c.probs)?;
                Ok(lines.iter().map(|line| line.as_deref().and_then(parse_record)).collect())
            })
            .collect::<Result<Vec<Vec<Option<(f64, Kind)>>>>>()
    })?;
    let (mut het_snps, mut variants, mut references) = (Vec::new(), Vec::new(), Vec::new());
    for &(qual, kind) in records.iter().flatten().flatten() {
        match kind {
            Kind::Reference => references.push(qual),
            Kind::HetSnp => {
                het_snps.push(qual);
                variants.push(qual);
            }
            Kind::Variant => variants.push(qual),
        }
    }
    let ref_pct = if settings.ref_pct_full == 0.0 { settings.var_pct_full } else { settings.ref_pct_full };
    let phase = lowest_cutoff(het_snps, 1.0 - settings.var_pct_phasing);
    let variant = lowest_cutoff(variants, settings.var_pct_full);
    let reference = lowest_cutoff(references, ref_pct);
    let thresholds = [output.qual as f64, final_qual];
    let near = |qual: f64, cutoff: Option<f64>| cutoff.is_some_and(|c| (qual - c).abs() <= QUAL_WINDOW);
    let near_fa = |qual: f64, cutoff: Option<f64>| cutoff.is_some_and(|c| (qual - c).abs() <= FULL_ALIGNMENT_WINDOW);
    let mut chosen: Vec<(usize, usize)> = Vec::new();
    for (j, parsed) in records.iter().enumerate() {
        for (i, (qual, kind)) in parsed.iter().enumerate().filter_map(|(i, r)| r.map(|r| (i, r))) {
            let close = qual < LOW_QUAL
                || thresholds.iter().any(|&t| near(qual, Some(t)))
                || match kind {
                    Kind::HetSnp => near(qual, phase) || near_fa(qual, variant),
                    Kind::Variant => near_fa(qual, variant),
                    Kind::Reference => near_fa(qual, reference),
                };
            if close {
                chosen.push((j, i));
            }
        }
    }
    drop(records);

    // Batches handed out as threads free up, each worker with its own network.
    let chunks_ref = &*chunks;
    let rescored: Vec<Vec<f32>> = pool.install(|| {
        chosen
            .par_chunks(batch)
            .with_min_len(1)
            .map_init(
                || load().map(|model| (model, Vec::new())),
                |state, group| {
                    let (model, inputs) = state.as_mut().map_err(|e| anyhow::anyhow!("{e:#}"))?;
                    inputs.clear();
                    for &(j, i) in group {
                        chunks_ref[j].packed[i].unpack_into(inputs);
                    }
                    Ok(model.predict(inputs)?)
                },
            )
            .collect::<Result<_>>()
    })?;
    for (&(j, i), p) in chosen.iter().zip(rescored.iter().flat_map(|p| p.as_chunks::<PILEUP_OUTPUTS>().0)) {
        chunks[j].probs[i * PILEUP_OUTPUTS..(i + 1) * PILEUP_OUTPUTS].copy_from_slice(p);
    }
    Ok(chosen.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_inputs_round_trip() {
        let mut values: Vec<i32> = vec![0; 600];
        values[0] = 5;
        values[300] = -128;
        values[301] = 127;
        values[302] = 40_000;
        values[303] = -1;
        values.extend([1, 0, 0, i32::MIN, i32::MAX, 64, -64, 63, -65]);
        let mut out = Vec::new();
        PackedInputs::new(&values).unpack_into(&mut out);
        assert_eq!(out, values);
    }

    #[test]
    fn records_are_classified_as_select_qual_reads_them() {
        let record = |r: &str, a: &str, gt: &str| format!("chr1\t10\t.\t{r}\t{a}\t12.5\tPASS\tP\tGT:GQ\t{gt}:9\n");
        assert_eq!(parse_record(&record("A", "C", "0/1")), Some((12.5, Kind::HetSnp)));
        assert_eq!(parse_record(&record("A", "AC", "0/1")), Some((12.5, Kind::Variant)));
        assert_eq!(parse_record(&record("A", "C", "1/1")), Some((12.5, Kind::Variant)));
        assert_eq!(parse_record(&record("A", ".", "0/0")), Some((12.5, Kind::Reference)));
    }

    #[test]
    fn cutoffs_take_the_highest_of_the_lowest_share() {
        assert_eq!(lowest_cutoff(vec![5.0, 1.0, 3.0, 2.0], 0.5), Some(2.0));
        assert_eq!(lowest_cutoff(vec![5.0], 0.5), None);
    }
}

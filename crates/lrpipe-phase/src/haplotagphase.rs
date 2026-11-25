//! Port of `cli/haplotagphase.py`: phase variants in a VCF using the HP and
//! PS tags of haplotagged reads.

use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result, bail};

use crate::core::read::{Read, ReadSet};
use crate::phase::{ReferenceCache, add_phased_blocks_and_sort, detect_file_format, read_alignments};
use crate::variants::ReadSetReader;
use crate::vcf::{PhaseTag, PhasedVcfWriter, VariantCallPhase, VcfReader, VcfVariant};

pub struct HaplotagphaseOptions {
    pub variant_file: String,
    pub alignment_file: String,
    /// "-": standard output
    pub output: String,
    pub reference: Option<String>,
    pub gap_threshold: i64,
    pub cut_poly: i64,
    pub only_indels: bool,
    pub samples: Vec<String>,
    pub ignore_read_groups: bool,
    pub chromosomes: Vec<String>,
    pub excluded_chromosomes: Vec<String>,
    pub mav: bool,
    pub command_line: String,
}

/// An insertion-ordered dict
struct OrderedMap<K, V> {
    index: HashMap<K, usize>,
    entries: Vec<(K, V)>,
}

impl<K: std::hash::Hash + Eq + Clone, V> OrderedMap<K, V> {
    fn new() -> Self {
        OrderedMap { index: HashMap::new(), entries: Vec::new() }
    }

    fn contains_key(&self, k: &K) -> bool {
        self.index.contains_key(k)
    }

    fn get_mut(&mut self, k: &K) -> Option<&mut V> {
        self.index.get(k).map(|&i| &mut self.entries[i].1)
    }

    fn entry_or_insert_with(&mut self, k: K, f: impl FnOnce() -> V) -> &mut V {
        let n = self.entries.len();
        let i = *self.index.entry(k.clone()).or_insert(n);
        if i == n {
            self.entries.push((k, f()));
        }
        &mut self.entries[i].1
    }
}

/// Votes per variant position: (phase set index, haplotype) -> summed quality
type Votes = OrderedMap<i32, OrderedMap<(i32, i64), i64>>;

/// `compute_votes`
fn compute_votes(
    is_homozygous: &HashMap<i32, bool>,
    reads: &ReadSet,
    allele_to_id: &HashMap<i32, HashMap<i32, i64>>,
) -> Result<Votes> {
    let mut votes: Votes = OrderedMap::new();
    let mut number_of_skipped = 0;
    for read in reads.iter() {
        let (ps, ht) = (read.ps_tag - 1, i64::from(read.hp_tag) - 1);
        if ht < 0 || ps < 0 {
            continue;
        }
        if ht > 1 {
            number_of_skipped += 1;
            continue;
        }
        for variant in &read.variants {
            let pos = variant.position;
            if *is_homozygous.get(&pos).with_context(|| format!("KeyError: {pos}"))? {
                continue;
            }
            let vote = votes.entry_or_insert_with(pos, OrderedMap::new);
            if !vote.contains_key(&(ps, 0)) {
                vote.entry_or_insert_with((ps, 0), || 0);
                vote.entry_or_insert_with((ps, 1), || 0);
            }
            let id = allele_to_id
                .get(&pos)
                .and_then(|m| m.get(&variant.allele))
                .with_context(|| format!("KeyError: {}", variant.allele))?;
            let key = (ps, ht ^ id);
            *vote.get_mut(&key).with_context(|| format!("KeyError: {key:?}"))? += i64::from(variant.quality);
        }
    }
    if number_of_skipped > 0 {
        log::warn!(
            "{number_of_skipped} reads were skipped due incorrect HP. The haplotagphase command supports only a diploid input"
        );
    }
    Ok(votes)
}

/// `best_candidate`: (allele, phase set, fraction of the total score, score)
/// of the first candidate with the highest score
fn best_candidate(var: &OrderedMap<(i32, i64), i64>) -> Result<(i64, i32, f64, i64)> {
    let mut best: Option<&((i32, i64), i64)> = None;
    for e in &var.entries {
        if best.is_none_or(|b| e.1 > b.1) {
            best = Some(e);
        }
    }
    let &((phase_set, allele), score) = best.expect("a vote has candidates");
    let total: i64 = var.entries.iter().map(|e| e.1).sum();
    if total == 0 {
        bail!("ZeroDivisionError: division by zero");
    }
    Ok((allele, phase_set, score as f64 / total as f64, score))
}

/// `length_of_homopolymer`: run of `ref[start]` from `start` in steps of
/// `step`, at most `threshold`
fn length_of_homopolymer(reference: &[u8], start: i64, step: i64, threshold: i64) -> i64 {
    let mut res = 0;
    let mut i = start;
    while res < threshold
        && 0 <= i
        && (i as usize) < reference.len()
        && reference[i as usize] == reference[start as usize]
    {
        res += 1;
        i += step;
    }
    res
}

/// `consensus`: the two super reads and the phase set of each voted position
#[allow(clippy::too_many_arguments)]
fn consensus(
    only_indels: bool,
    gap_threshold: i64,
    cut_homopolymers: i64,
    refseq: &[u8],
    change: &HashMap<i32, &VcfVariant>,
    phased: &HashMap<i32, &Option<VariantCallPhase>>,
    votes: &Votes,
    id_to_allele: &HashMap<i32, HashMap<i64, u32>>,
) -> Result<(ReadSet, HashMap<i32, i32>)> {
    let mut super_reads = [Read::new("superread_0", 0, 0, 0), Read::new("superread_1", 0, 0, 0)];
    let mut components = HashMap::new();
    for (pos, vote) in &votes.entries {
        let pos = *pos;
        let (best_allele, phase_set, fraction, score) = best_candidate(vote)?;
        components.insert(pos, phase_set);
        if phased[&pos].is_none() {
            if 100.0 * fraction < gap_threshold as f64 {
                continue;
            }
            if only_indels && change[&pos].is_snv() {
                continue;
            }
            if cut_homopolymers > 0 {
                let p = i64::from(pos);
                let max_length = length_of_homopolymer(refseq, p + 1, 1, cut_homopolymers).max(length_of_homopolymer(
                    refseq,
                    p,
                    -1,
                    cut_homopolymers,
                ));
                if max_length > cut_homopolymers {
                    continue;
                }
            }
        }
        let allele = |a: i64| id_to_allele[&pos].get(&a).copied().with_context(|| format!("KeyError: {a}"));
        let quality = i32::try_from(score).unwrap_or(i32::MAX);
        super_reads[0].add_variant(pos, allele(best_allele)? as i32, quality);
        super_reads[1].add_variant(pos, allele(1 - best_allele)? as i32, quality);
    }
    let mut readset = ReadSet::new();
    for mut read in super_reads {
        read.variants.sort_by_key(|v| v.position);
        readset.add(read)?;
    }
    Ok((readset, components))
}

pub fn run_haplotagphase(opts: &HaplotagphaseOptions) -> Result<()> {
    let Some(reference) = &opts.reference else { bail!("Option --reference should be specified") };
    crate::fasta::check_reference(reference)?;
    // PhasedInputReader
    match detect_file_format(&opts.alignment_file)? {
        Some("BAM" | "CRAM") => {}
        Some(_) => bail!("ValueError: call PhasedInputReader.read_vcfs() first"),
        None => bail!("Unable to determine type of input file '{}'", opts.alignment_file),
    }
    // (WhatsHap opens the reference twice, for the reads and for homopolymers)
    let mut fasta = ReferenceCache::open(reference)?;
    let mut readset_reader =
        ReadSetReader::new(std::slice::from_ref(&opts.alignment_file), Some(reference), 20, false, 100_000, 1)?;
    let mut writer =
        PhasedVcfWriter::new(&opts.variant_file, Some(&opts.command_line), &opts.output, PhaseTag::PS, false)?
            .with_mav(opts.mav);
    let mut vcf_reader = VcfReader::open_with(&opts.variant_file, false, true, None)?.with_mav(opts.mav);
    if opts.ignore_read_groups && opts.samples.is_empty() && vcf_reader.samples.len() > 1 {
        bail!("When using --ignore-read-groups on a VCF with multiple samples, --sample must also be used.");
    }
    let samples: Vec<String> = if opts.samples.is_empty() { vcf_reader.samples.clone() } else { opts.samples.clone() };
    let vcf_samples: HashSet<&String> = vcf_reader.samples.iter().collect();
    for sample in &samples {
        if !vcf_samples.contains(sample) {
            bail!("Sample '{sample}' requested on command-line not found in VCF");
        }
    }
    let mut numeric_ids: HashMap<String, i32> = HashMap::new();
    while let Some(table) = vcf_reader.next_table()? {
        let chromosome = table.chromosome.clone();
        let reference_seq: &[u8] = fasta.get(&chromosome)?;
        log::info!("Processing chromosome {chromosome}...");
        let included = (opts.chromosomes.is_empty() || opts.chromosomes.contains(&chromosome))
            && !opts.excluded_chromosomes.contains(&chromosome);
        if !included {
            log::info!("Leaving chromosome {chromosome} unchanged");
            writer.write_unchanged(&chromosome)?;
            continue;
        }
        let mut sample_to_super_reads: Vec<(String, ReadSet)> = Vec::new();
        let mut sample_to_components: HashMap<String, HashMap<i32, i32>> = HashMap::new();
        for sample in &vcf_reader.samples {
            log::info!("Processing sample {sample}");
            let genotypes = table.genotypes_of(sample);
            // PhasedInputReader.read
            let bam_sample = if opts.ignore_read_groups { None } else { Some(sample.as_str()) };
            let n = numeric_ids.len() as i32;
            let numeric_id = match bam_sample {
                None => 0,
                Some(s) => *numeric_ids.entry(s.to_string()).or_insert(n),
            };
            let mut reads = read_alignments(
                &mut readset_reader,
                &chromosome,
                &table.variants,
                bam_sample,
                numeric_id,
                Some(reference_seq),
                Some(genotypes),
            )?;
            // (no VCF inputs: the reads are only sorted)
            add_phased_blocks_and_sort(&mut reads, &[], 0, &chromosome, &table.variants, sample, numeric_id)?;
            let phases = table.phases_of(sample);
            if !samples.contains(sample) {
                log::info!("Skipping sample {sample}");
                continue;
            }
            let mut homozygous = HashMap::new();
            let mut change = HashMap::new();
            let mut phased = HashMap::new();
            // alleles of the genotype -> 0/1 and back
            let mut allele_to_id: HashMap<i32, HashMap<i32, i64>> = HashMap::new();
            let mut id_to_allele: HashMap<i32, HashMap<i64, u32>> = HashMap::new();
            let (mut homozygous_number, mut phased_number) = (0, 0);
            for ((variant, phase), genotype) in table.variants.iter().zip(phases).zip(genotypes) {
                let pos = variant.position;
                for (i, v) in genotype.as_vector().into_iter().enumerate() {
                    allele_to_id.entry(pos).or_default().insert(v as i32, i as i64);
                    id_to_allele.entry(pos).or_default().insert(i as i64, v);
                }
                homozygous.insert(pos, genotype.is_homozygous());
                phased.insert(pos, phase);
                phased_number += usize::from(phase.is_some());
                homozygous_number += usize::from(genotype.is_homozygous());
                change.insert(pos, variant);
            }
            log::info!("Number of homozygous variants is {homozygous_number}");
            log::info!("Number of already phased variants is {phased_number}");
            let votes = compute_votes(&homozygous, &reads, &allele_to_id)?;
            let (super_reads, components) = consensus(
                opts.only_indels,
                opts.gap_threshold,
                opts.cut_poly,
                reference_seq,
                &change,
                &phased,
                &votes,
                &id_to_allele,
            )?;
            sample_to_super_reads.push((sample.clone(), super_reads));
            sample_to_components.insert(sample.clone(), components);
        }
        writer.write(&chromosome, &sample_to_super_reads, &sample_to_components)?;
    }
    writer.close()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homopolymer_length() {
        assert_eq!(length_of_homopolymer(b"AAABBBCCC", 0, 1, 10), 3);
        assert_eq!(length_of_homopolymer(b"AAABBBCCC", 2, -1, 10), 3);
        assert_eq!(length_of_homopolymer(b"AAABBBCCC", 3, 1, 2), 2);
        assert_eq!(length_of_homopolymer(b"AABBBCCCC", 5, 1, 5), 4);
        assert_eq!(length_of_homopolymer(b"", 0, 1, 10), 0);
    }

    #[test]
    fn best_candidate_first_of_the_highest() {
        let mut m = OrderedMap::new();
        for (k, v) in [((1, 2), 100), ((2, 2), 100), ((3, 3), 100)] {
            m.entry_or_insert_with(k, || v);
        }
        let (allele, ps, q, score) = best_candidate(&m).unwrap();
        assert_eq!((allele, ps, score), (2, 1, 100));
        assert!((q - 1.0 / 3.0).abs() < 1e-15);
    }
}

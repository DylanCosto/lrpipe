//! Port of `cli/genotype.py`: genotype variants with the forward-backward
//! algorithm (priors from all reads, then a DP over phase-informative reads).

use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::compat::pystr::{py_int, py_strip};
use crate::core::genotype::{Genotype, PhredGenotypeLikelihoods};
use crate::core::genotypedp::{GenotypeDPTable, compute_genotypes};
use crate::core::pedigree::Pedigree;
use crate::core::read::ReadSet;
use crate::pedigree::{Trio, load_genetic_map, ped_samples, read_ped};
use crate::phase::{
    RecombinationCostComputer, ReferenceCache, add_phased_blocks_and_sort, detect_file_format, read_alignments,
    setup_families,
};
use crate::readselect::readselection;
use crate::variants::{ReadSetReader, Realignment};
use crate::vcf::{GenotypeVcfWriter, SampleGenotypes, VariantTable, VcfReader, VcfVariant};

pub struct GenotypeOptions {
    pub variant_file: String,
    pub phase_input_files: Vec<String>,
    pub reference: Option<String>,
    /// "-": standard output
    pub output: String,
    pub samples: Vec<String>,
    pub chromosomes: Vec<String>,
    pub excluded_chromosomes: Vec<String>,
    pub ignore_read_groups: bool,
    pub only_snvs: bool,
    pub mapping_quality: i32,
    pub max_coverage: usize,
    pub nopriors: bool,
    pub ped: Option<String>,
    pub recombrate: f64,
    pub genmap: Option<String>,
    pub gt_qual_threshold: f64,
    pub prioroutput: Option<String>,
    pub constant: f64,
    pub overhang: i64,
    pub realignment: Realignment,
    pub use_ped_samples: bool,
    pub command_line: String,
    /// chromosomes genotyped in parallel (lrpipe-phase only; the output does
    /// not change)
    pub threads: usize,
}

/// `PhasedInputReader` for genotyping (one per thread)
struct InputReader {
    readset_reader: Option<ReadSetReader>,
    fasta: Option<ReferenceCache>,
}

/// Everything shared by the chromosomes
struct Shared<'a> {
    opts: &'a GenotypeOptions,
    samples: Vec<String>,
    families: BTreeMap<String, Vec<String>>,
    family_trios: HashMap<String, Vec<Trio>>,
    /// `NumericSampleIds`: WhatsHap assigns ids in order of first use, which
    /// is the same on every chromosome; computed in advance
    ids: HashMap<String, u32>,
    bam_paths: Vec<String>,
    phase_vcfs: Vec<HashMap<String, VariantTable>>,
    recombination_cost_computer: RecombinationCostComputer,
    gt_prob: f64,
}

/// Genotyping results of one chromosome
struct ChromosomeCalls {
    calls: HashMap<String, SampleGenotypes>,
    /// the priors, if `--prioroutput` is given
    priors: Option<HashMap<String, SampleGenotypes>>,
}

fn open_reference(opts: &GenotypeOptions) -> Result<Option<ReferenceCache>> {
    opts.reference.as_deref().map(ReferenceCache::open).transpose()
}

fn open_readset_reader(opts: &GenotypeOptions, bam_paths: &[String]) -> Result<Option<ReadSetReader>> {
    if bam_paths.is_empty() {
        return Ok(None);
    }
    let mut r = ReadSetReader::new(bam_paths, opts.reference.as_deref(), opts.mapping_quality, false, 100_000, 1)?;
    r.set_overhang(opts.overhang);
    r.set_realignment(opts.realignment.clone());
    Ok(Some(r))
}

impl Shared<'_> {
    fn open_input(&self) -> Result<InputReader> {
        Ok(InputReader {
            fasta: open_reference(self.opts)?,
            readset_reader: open_readset_reader(self.opts, &self.bam_paths)?,
        })
    }

    /// `PhasedInputReader.read`: reads of `sample`, with the phased blocks of
    /// the VCF inputs if `read_vcf`, sorted; and the source ids of the VCFs
    fn read(
        &self,
        input: &mut InputReader,
        chromosome: &str,
        variants: &[VcfVariant],
        sample: &str,
        read_vcf: bool,
    ) -> Result<(ReadSet, HashSet<i32>)> {
        let reference = match input.fasta.as_mut() {
            Some(f) => Some(f.get(chromosome)?),
            None => None,
        };
        let bam_sample = if self.opts.ignore_read_groups { None } else { Some(sample) };
        // (WhatsHap's hard-coded zero without read groups)
        let numeric_id = bam_sample.map_or(0, |s| self.ids[s]);
        let mut readset = match &mut input.readset_reader {
            None => ReadSet::new(),
            Some(reader) => {
                read_alignments(reader, chromosome, variants, bam_sample, numeric_id as i32, reference, None)?
            }
        };
        let phase_vcfs: &[HashMap<String, VariantTable>] = if read_vcf { &self.phase_vcfs } else { &[] };
        let sample_id = if read_vcf { self.ids[sample] as i32 } else { 0 };
        let vcf_source_ids = add_phased_blocks_and_sort(
            &mut readset,
            phase_vcfs,
            self.bam_paths.len(),
            chromosome,
            variants,
            sample,
            sample_id,
        )?;
        Ok((readset, vcf_source_ids))
    }

    /// Priors and genotypes of one chromosome, by sample
    fn genotype_table(&self, input: &mut InputReader, table: &VariantTable) -> Result<ChromosomeCalls> {
        let opts = self.opts;
        let chromosome = &table.chromosome;
        let var_to_pos: HashMap<i32, usize> = table.variants.iter().enumerate().map(|(i, v)| (v.position, i)).collect();
        log::info!("======== Working on chromosome '{chromosome}'");
        let n = table.variants.len();
        let positions: Vec<i32> = table.variants.iter().map(|v| v.position).collect();
        let mut calls: HashMap<String, SampleGenotypes> = table
            .samples
            .iter()
            .map(|s| (s.clone(), SampleGenotypes { genotypes: vec![Genotype::none(); n], likelihoods: vec![None; n] }))
            .collect();
        if !opts.nopriors {
            // prior genotype likelihoods from all reads
            for sample in &self.samples {
                log::info!("---- Initial genotyping of {sample}");
                let (mut readset, _) = self.read(input, chromosome, &table.variants, sample, false)?;
                readset.sort();
                let (mut genotypes, genotype_likelihoods) = compute_genotypes(&readset, &positions)?;
                let mut regularized = Vec::with_capacity(n);
                for (gl, genotype) in genotype_likelihoods.iter().zip(genotypes.iter_mut()) {
                    let c = opts.constant;
                    let norm_sum = gl[0] + gl[1] + gl[2] + 3.0 * c;
                    let r = [(gl[0] + c) / norm_sum, (gl[1] + c) / norm_sum, (gl[2] + c) / norm_sum];
                    *genotype = determine_genotype(&r, self.gt_prob);
                    regularized.push(Some(r));
                }
                let entry = calls.get_mut(sample).expect("samples are in the VCF");
                entry.likelihoods = regularized;
                entry.genotypes = genotypes;
            }
        } else {
            for sample in &self.samples {
                calls.get_mut(sample).expect("samples are in the VCF").likelihoods = vec![Some([1.0 / 3.0; 3]); n];
            }
        }
        let priors = opts.prioroutput.as_ref().map(|_| calls.clone());
        let no_trios = Vec::new();
        for (representative, family) in &self.families {
            if family.len() == 1 {
                log::info!("---- Processing individual {representative}");
            } else {
                log::info!("---- Processing family with individuals: {}", family.join(","));
            }
            let max_coverage_per_sample = (opts.max_coverage / family.len()).max(1);
            log::info!("Using maximum coverage per sample of {max_coverage_per_sample}X");
            let trios = self.family_trios.get(representative).unwrap_or(&no_trios);
            let mut all_reads = ReadSet::new();
            for sample in family {
                let (readset, vcf_source_ids) = self.read(input, chromosome, &table.variants, sample, true)?;
                let keep: Vec<usize> = (0..readset.len()).filter(|&i| readset.get(i).len() >= 2).collect();
                let readset = readset.into_subset(keep);
                log::info!("Kept {} reads that cover at least two variants each", readset.len());
                let selected = readselection(&readset, max_coverage_per_sample, Some(&vcf_source_ids), true)?;
                for read in readset.into_subset(selected).iter() {
                    all_reads.add(read.clone())?;
                }
            }
            all_reads.sort();
            let accessible_positions = all_reads.get_positions();
            log::info!(
                "Variants covered by at least one phase-informative read in at least one individual after read selection: {}",
                accessible_positions.len()
            );
            let mut pedigree = Pedigree::new();
            for sample in family {
                let likelihoods = &calls[sample].likelihoods;
                let genotype_l = accessible_positions
                    .iter()
                    .map(|p| {
                        likelihoods[var_to_pos[p]]
                            .map(|gl| PhredGenotypeLikelihoods::new(gl.to_vec(), 2, 2))
                            .transpose()
                    })
                    .collect::<Result<Vec<_>>>()?;
                pedigree.add_individual(
                    self.ids[sample],
                    vec![Genotype::none(); accessible_positions.len()],
                    genotype_l,
                );
            }
            for trio in trios {
                let (Some(father), Some(mother)) = (&trio.father, &trio.mother) else { continue };
                pedigree.add_relationship(self.ids[father], self.ids[mother], self.ids[&trio.child])?;
            }
            let recombination_costs = self.recombination_cost_computer.compute(&accessible_positions);
            log::info!(
                "Genotype {} sample{} by solving the genotyping problem ...",
                family.len(),
                if family.len() > 1 { "s" } else { "" }
            );
            let dp = GenotypeDPTable::new(&mut all_reads, &recombination_costs, &pedigree, &accessible_positions)?;
            for s in family {
                let id = self.ids[s];
                let entry = calls.get_mut(s).context("sample not in the VCF")?;
                for (k, p) in accessible_positions.iter().enumerate() {
                    let likelihoods = dp.genotype_likelihoods(id, k);
                    entry.genotypes[var_to_pos[p]] = determine_genotype(&likelihoods, self.gt_prob);
                    entry.likelihoods[var_to_pos[p]] = Some(likelihoods);
                }
            }
        }
        Ok(ChromosomeCalls { calls, priors })
    }
}

/// `int_to_diploid_biallelic_gt`
fn diploid_genotype(numeric: i64) -> Genotype {
    match numeric {
        0 => Genotype::new(&[0, 0]).expect("valid genotype"),
        1 => Genotype::new(&[0, 1]).expect("valid genotype"),
        2 => Genotype::new(&[1, 1]).expect("valid genotype"),
        _ => Genotype::none(),
    }
}

/// `determine_genotype`: the likeliest genotype if it is unique and above
/// the threshold
fn determine_genotype(likelihoods: &[f64; 3], threshold_prob: f64) -> Genotype {
    let mut to_sort = [(likelihoods[0], 0), (likelihoods[1], 1), (likelihoods[2], 2)];
    // Python's sort by likelihood (also where a likelihood is NaN)
    crate::compat::pysort::sort_small(&mut to_sort, |a, b| a.0 < b.0);
    if to_sort[2].0 > to_sort[1].0 && to_sort[2].0 > threshold_prob {
        diploid_genotype(to_sort[2].1)
    } else {
        diploid_genotype(-1)
    }
}

/// The result of a chromosome, to be written in VCF order
enum ChromosomeResult {
    Unchanged(String),
    Genotyped { chromosome: String, positions: Vec<i32>, calls: ChromosomeCalls },
}

pub fn run_genotype(opts: &GenotypeOptions) -> Result<()> {
    log::info!("This is lrpipe-phase {} (genotyping; port of WhatsHap 2.8)", env!("CARGO_PKG_VERSION"));
    if let Some(reference) = opts.reference.as_deref() {
        crate::fasta::check_reference(reference)?;
    }
    // PhasedInputReader
    let mut bam_paths = Vec::new();
    let mut vcf_paths = Vec::new();
    for path in &opts.phase_input_files {
        match detect_file_format(path)? {
            Some("BAM" | "CRAM") => bam_paths.push(path.clone()),
            Some("VCF") => vcf_paths.push(path.clone()),
            _ => bail!("Unable to determine type of input file '{path}'"),
        }
    }
    // opened in the order of WhatsHap's PhasedInputReader (used by the first
    // or only thread)
    let fasta = open_reference(opts)?;
    let mut phase_vcf_readers = Vec::new();
    for path in &vcf_paths {
        phase_vcf_readers.push(VcfReader::open_with(path, opts.only_snvs, true, None)?);
    }
    let first_input = InputReader { fasta, readset_reader: open_readset_reader(opts, &bam_paths)? };
    let mut vcf_writer = GenotypeVcfWriter::new(&opts.variant_file, Some(&opts.command_line), &opts.output, false)?;
    let mut prior_vcf_writer = match &opts.prioroutput {
        Some(path) => Some(GenotypeVcfWriter::new(&opts.variant_file, Some(&opts.command_line), path, true)?),
        None => None,
    };
    let mut vcf_reader = VcfReader::open(&opts.variant_file, opts.only_snvs)?.with_ignore_genotypes(true);
    if opts.ignore_read_groups && opts.samples.is_empty() && vcf_reader.samples.len() > 1 {
        bail!("When using --ignore-read-groups on a VCF with multiple samples, --sample must also be used.");
    }
    let ped_trios = match &opts.ped {
        Some(path) => Some(read_ped(path)?),
        None => None,
    };
    let mut samples: Vec<String> =
        if opts.samples.is_empty() { vcf_reader.samples.clone() } else { opts.samples.clone() };
    if let (Some(trios), true) = (&ped_trios, opts.use_ped_samples) {
        samples = ped_samples(trios);
    }
    let vcf_samples: HashSet<&String> = vcf_reader.samples.iter().collect();
    for s in &samples {
        if !vcf_samples.contains(s) {
            bail!("Sample '{s}' requested on command-line not found in VCF");
        }
    }
    let recombination_cost_computer = match (&opts.ped, &opts.genmap) {
        (Some(_), Some(genmap)) => {
            log::info!("Using region-specific recombination rates from genetic map {genmap}.");
            RecombinationCostComputer::GeneticMap(load_genetic_map(genmap)?)
        }
        _ => {
            if opts.ped.is_some() {
                log::info!("Using uniform recombination rate of {} cM/Mb.", opts.recombrate);
            }
            RecombinationCostComputer::Uniform(opts.recombrate)
        }
    };
    // WhatsHap uses a frozenset (hash-seed-dependent order); here: first
    // occurrence order
    let mut unique = HashSet::new();
    samples.retain(|s| unique.insert(s.clone()));
    let (families, family_trios) = setup_families(&samples, ped_trios.as_deref(), opts.max_coverage);
    // numeric sample ids in the order WhatsHap assigns them: children of
    // trios, the samples while computing priors (with read groups), then
    // the members of each family
    let mut ids: HashMap<String, u32> = HashMap::new();
    let mut assign = |s: &str| {
        let n = ids.len() as u32;
        ids.entry(s.to_string()).or_insert(n);
    };
    for trios in family_trios.values() {
        for trio in trios {
            assign(&trio.child);
        }
    }
    if !opts.nopriors && !opts.ignore_read_groups {
        for s in &samples {
            assign(s);
        }
    }
    for family in families.values() {
        for s in family {
            assign(s);
        }
    }
    let mut phase_vcfs = Vec::new();
    for mut reader in phase_vcf_readers {
        let mut m = HashMap::new();
        log::info!("Reading phased blocks from '{}'", reader.path());
        while let Some(table) = reader.next_table()? {
            m.insert(table.chromosome.clone(), table);
        }
        phase_vcfs.push(m);
    }
    let ctx = Shared {
        opts,
        samples,
        families,
        family_trios,
        ids,
        bam_paths,
        phase_vcfs,
        recombination_cost_computer,
        gt_prob: 1.0 - 10f64.powf(-opts.gt_qual_threshold / 10.0),
    };
    let mut write = |result: ChromosomeResult| -> Result<()> {
        match result {
            ChromosomeResult::Unchanged(chromosome) => {
                vcf_writer.write_unchanged(&chromosome)?;
                if let Some(w) = &mut prior_vcf_writer {
                    w.write_unchanged(&chromosome)?;
                }
            }
            ChromosomeResult::Genotyped { chromosome, positions, calls } => {
                if let (Some(w), Some(priors)) = (&mut prior_vcf_writer, &calls.priors) {
                    w.write_genotypes(&chromosome, &positions, priors)?;
                }
                log::info!("======== Writing VCF");
                vcf_writer.write_genotypes(&chromosome, &positions, &calls.calls)?;
                log::info!("Done writing VCF");
            }
        }
        Ok(())
    };
    let included = |chromosome: &str| {
        (opts.chromosomes.is_empty() || opts.chromosomes.iter().any(|c| c == chromosome))
            && !opts.excluded_chromosomes.iter().any(|c| c == chromosome)
    };
    let leave_unchanged = |chromosome: &str| {
        log::info!(
            "Leaving chromosome '{chromosome}' unchanged (present in VCF but not requested by option --chromosome)"
        );
    };
    let threads = opts.threads.max(1);
    if threads == 1 {
        let mut input = first_input;
        while let Some(table) = vcf_reader.next_table()? {
            if !included(&table.chromosome) {
                leave_unchanged(&table.chromosome);
                write(ChromosomeResult::Unchanged(table.chromosome))?;
                continue;
            }
            let calls = ctx.genotype_table(&mut input, &table)?;
            let positions = table.variants.iter().map(|v| v.position).collect();
            write(ChromosomeResult::Genotyped { chromosome: table.chromosome, positions, calls })?;
        }
    } else {
        // chromosomes in parallel, written in VCF order
        let (sender, receiver) = std::sync::mpsc::channel::<(usize, Result<ChromosomeResult>)>();
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build()?;
        let ctx = &ctx;
        let mut outcome: Result<()> = Ok(());
        pool.in_place_scope(|scope| {
            let mut pending: BTreeMap<usize, Result<ChromosomeResult>> = BTreeMap::new();
            let (mut next_index, mut next_to_write, mut in_flight) = (0usize, 0usize, 0usize);
            let mut eof = false;
            loop {
                while outcome.is_ok() && !eof && in_flight < threads {
                    match vcf_reader.next_table() {
                        Ok(Some(table)) => {
                            let index = next_index;
                            next_index += 1;
                            if !included(&table.chromosome) {
                                leave_unchanged(&table.chromosome);
                                pending.insert(index, Ok(ChromosomeResult::Unchanged(table.chromosome)));
                                continue;
                            }
                            in_flight += 1;
                            let sender = sender.clone();
                            scope.spawn(move |_| {
                                let result = ctx.open_input().and_then(|mut input| {
                                    let calls = ctx.genotype_table(&mut input, &table)?;
                                    let positions = table.variants.iter().map(|v| v.position).collect();
                                    Ok(ChromosomeResult::Genotyped { chromosome: table.chromosome, positions, calls })
                                });
                                let _ = sender.send((index, result));
                            });
                        }
                        Ok(None) => eof = true,
                        Err(e) => outcome = Err(e),
                    }
                }
                while outcome.is_ok()
                    && let Some(result) = pending.remove(&next_to_write)
                {
                    next_to_write += 1;
                    if let Err(e) = result.and_then(&mut write) {
                        outcome = Err(e);
                    }
                }
                if in_flight == 0 {
                    break;
                }
                let (index, result) = receiver.recv().expect("workers send results");
                in_flight -= 1;
                pending.insert(index, result);
            }
        });
        outcome?;
    }
    vcf_writer.close()?;
    if let Some(w) = prior_vcf_writer {
        w.close()?;
    }
    Ok(())
}

/// The k-mer alignment costs of `--kmeralign-costs` (tab-separated: k-mer,
/// k-mer, cost)
pub fn load_kmerald_costs(path: &str) -> Result<Arc<HashMap<(i64, i64), f32>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("cannot read {path}"))?;
    let mut costs = HashMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            bail!("IndexError: list index out of range in {path}");
        }
        let int = |s: &str| {
            py_int(s)
                .and_then(|v| i64::try_from(v).ok())
                .with_context(|| format!("invalid literal for int() with base 10: '{s}'"))
        };
        let (k1, k2) = (int(f[0])?, int(f[1])?);
        let cost =
            py_strip(f[2]).parse::<f64>().with_context(|| format!("could not convert string to float: '{}'", f[2]))?;
        costs.insert((k1, k2), cost as f32);
    }
    Ok(Arc::new(costs))
}

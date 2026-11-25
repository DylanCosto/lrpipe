//! Port of `cli/phase.py` (`run_whatshap`): read-based and pedigree-based
//! phasing.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail};

use crate::core::genotype::{Genotype, PhredGenotypeLikelihoods};
use crate::core::pedigree::Pedigree;
use crate::core::read::ReadSet;
use crate::core::{dp, hapchat, heuristic};
use crate::fasta::Fasta;
use crate::graph::ComponentFinder;
use crate::memory::{self, ChromosomeUsage, FinishOnDrop, MemoryScheduler};
use crate::merge::ReadMerger;
use crate::pedigree::{
    RecombinationMapEntry, Trio, find_recombination, load_genetic_map, mendelian_conflict, ped_samples, read_ped,
    recombination_cost_map, uniform_recombination_map,
};
use crate::readselect::readselection;
use crate::variants::{ReadSetError, ReadSetReader};
use crate::vcf::{GenotypeChange, PhaseTag, PhasedVcfWriter, VariantTable, VcfReader};

/// Threads decompressing each BAM input file (each reader decompresses up
/// to this many BGZF blocks ahead). With htslib's default of two queued
/// blocks per thread, each CRAM reader held up to 48 decoded containers with
/// 24 threads (about 2 GB with HG002 ONT reads) without making whole-genome
/// runs faster.
const DECOMPRESSION_THREADS: usize = 8;

/// `--algorithm`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Algorithm {
    /// The exact (Ped)MEC dynamic program
    WhatsHap,
    /// `PedMecHeuristic`, keeping up to `row_limit` partial solutions
    Heuristic { row_limit: u32 },
    /// HapCHAT
    HapChat,
}

pub struct PhaseOptions {
    pub variant_file: String,
    pub phase_input_files: Vec<String>,
    pub reference: Option<String>,
    pub no_reference: bool,
    pub output: String,
    pub samples: Vec<String>,
    pub chromosomes: Vec<String>,
    pub excluded_chromosomes: Vec<String>,
    pub ignore_read_groups: bool,
    pub only_snvs: bool,
    pub mapping_quality: i32,
    pub max_coverage: usize,
    pub tag: PhaseTag,
    pub read_list_filename: Option<String>,
    pub use_supplementary: bool,
    pub supplementary_distance_threshold: i64,
    pub command_line: Option<String>,
    pub threads: usize,
    /// Memory (MB) for keeping DP backtrace tables instead of recomputing them
    pub dp_memory_mb: u64,
    /// Approximate memory limit for phasing chromosomes in parallel (None:
    /// half of the system memory)
    pub max_memory_mb: Option<u64>,
    /// `--merge-reads`
    pub read_merging: bool,
    pub algorithm: Algorithm,
    /// PED file for pedigree phasing
    pub ped: Option<String>,
    pub use_ped_samples: bool,
    /// Uniform recombination rate in cM/Mb
    pub recombrate: f64,
    pub genmap: Option<String>,
    pub genetic_haplotyping: bool,
    pub recombination_list_filename: Option<String>,
    pub distrust_genotypes: bool,
    pub include_homozygous: bool,
    pub default_gq: i32,
    pub gl_regularizer: Option<f64>,
    pub gtchange_list_filename: Option<String>,
}

/// `NumericSampleIds`: ids are assigned in order of first use.
#[derive(Default)]
struct NumericSampleIds {
    mapping: HashMap<String, i32>,
    order: Vec<String>,
}

impl NumericSampleIds {
    fn get(&mut self, sample: &str) -> i32 {
        if let Some(&id) = self.mapping.get(sample) {
            return id;
        }
        let id = self.mapping.len() as i32;
        self.mapping.insert(sample.to_string(), id);
        self.order.push(sample.to_string());
        id
    }

    fn name(&self, id: i32) -> &str {
        &self.order[id as usize]
    }

    /// Id of a sample that has already been registered.
    fn id(&self, sample: &str) -> i32 {
        self.mapping[sample]
    }
}

/// Detects the file type as `utils.detect_file_format` does.
pub(crate) fn detect_file_format(path: &str) -> Result<Option<&'static str>> {
    use std::io::Read;
    let mut f = std::fs::File::open(path).with_context(|| format!("[Errno 2] No such file or directory: '{path}'"))?;
    let mut first = [0u8; 16];
    let n = f.read(&mut first)?;
    let first = &first[..n];
    if first.starts_with(b"CRAM") {
        return Ok(Some("CRAM"));
    }
    if first.starts_with(b"##fileformat=VCF") {
        return Ok(Some("VCF"));
    }
    if first.starts_with(b"\x1f\x8b") {
        let mut d = flate2::read::MultiGzDecoder::new(std::fs::File::open(path)?);
        let mut buf = [0u8; 16];
        let mut got = 0;
        while got < 16 {
            let k = d.read(&mut buf[got..])?;
            if k == 0 {
                break;
            }
            got += k;
        }
        let b = &buf[..got];
        if b.starts_with(b"BAM\x01") {
            return Ok(Some("BAM"));
        }
        if b.starts_with(b"##fileformat=VCF") {
            return Ok(Some("VCF"));
        }
    }
    Ok(None)
}

pub struct ReferenceCache {
    /// Read through pasta's cache, which CRAM readers can decode with too
    /// (see [`ReadSetReader::share_reference`]).
    repository: pasta::fasta::Repository,
    path: String,
    names: HashSet<String>,
    cached: Option<(String, Bases)>,
}

/// A chromosome's bases, uppercase.
enum Bases {
    /// The repository's own copy, which is already uppercase.
    Shared(Arc<pasta::fasta::record::Sequence>),
    /// An uppercase copy of a soft-masked sequence.
    Uppercased(Vec<u8>),
}

impl ReferenceCache {
    pub fn open(path: &str) -> Result<Self> {
        if !std::path::Path::new(&format!("{path}.fai")).exists() {
            bail!(
                "An index file (.fai) for the reference FASTA '{path}' could not be found. Please create one with 'samtools faidx'."
            );
        }
        let fasta = Fasta::open(path).map_err(|e| anyhow!("Error while opening FASTA reference file: {e}"))?;
        let names: HashSet<String> = fasta.names().into_iter().collect();
        let repository = crate::alignment::cram::repository(std::path::Path::new(path))
            .map_err(|e| anyhow!("Error while opening FASTA reference file: {e}"))?;
        Ok(ReferenceCache { repository, path: path.to_string(), names, cached: None })
    }

    /// The cache the sequences are read through.
    pub fn repository(&self) -> &pasta::fasta::Repository {
        &self.repository
    }

    /// Whether the bases [`Self::get`] gave are the repository's own copy
    /// (else a CRAM reader sharing the repository holds a second one).
    pub fn shares_bases(&self) -> bool {
        matches!(self.cached, Some((_, Bases::Shared(_))))
    }

    pub fn get(&mut self, chromosome: &str) -> Result<&[u8]> {
        if self.cached.as_ref().map(|(c, _)| c.as_str()) != Some(chromosome) {
            if !self.names.contains(chromosome) {
                bail!("Chromosome '{chromosome}' present in VCF file, but not in the reference FASTA '{}'", self.path);
            }
            // The last chromosome's sequence goes, here and in readers sharing the cache.
            self.cached = None;
            self.repository.clear();
            let seq = self
                .repository
                .get(chromosome.as_bytes())
                .ok_or_else(|| anyhow!("Chromosome '{chromosome}' not found in '{}'", self.path))??;
            let bases: &[u8] = seq.as_ref().as_ref();
            let bases = if bases.iter().any(u8::is_ascii_lowercase) {
                Bases::Uppercased(bases.to_ascii_uppercase())
            } else {
                Bases::Shared(seq)
            };
            self.cached = Some((chromosome.to_string(), bases));
        }
        Ok(match &self.cached.as_ref().expect("cached above").1 {
            Bases::Shared(seq) => seq.as_ref().as_ref(),
            Bases::Uppercased(bases) => bases,
        })
    }
}

/// `find_mendelian_conflicts`
fn find_mendelian_conflicts(trios: &[Trio], table: &VariantTable) -> HashSet<usize> {
    let mut conflicts = HashSet::new();
    for trio in trios {
        let (Some(mother), Some(father)) = (&trio.mother, &trio.father) else { continue };
        let gm = table.genotypes_of(mother);
        let gf = table.genotypes_of(father);
        let gc = table.genotypes_of(&trio.child);
        for (index, ((m, f), c)) in gm.iter().zip(gf).zip(gc).enumerate() {
            if !m.is_none() && !f.is_none() && !c.is_none() && mendelian_conflict(m, f, c) {
                conflicts.insert(index);
            }
        }
    }
    conflicts
}

/// `find_phaseable_variants`: (positions homozygous in some individual, table
/// of variants to phase)
fn find_phaseable_variants(
    family: &[String],
    include_homozygous: bool,
    trios: &[Trio],
    table: &VariantTable,
) -> Result<(Vec<i32>, VariantTable)> {
    let mut missing = HashSet::new();
    let mut heterozygous = HashSet::new();
    let mut homozygous = HashSet::new();
    for sample in family {
        for (index, gt) in table.genotypes_of(sample).iter().enumerate() {
            if gt.is_none() {
                missing.insert(index);
            } else if !gt.is_homozygous() {
                heterozygous.insert(index);
            } else {
                if !gt.is_diploid_and_biallelic() {
                    bail!("AssertionError: homozygous genotype {gt} is not diploid and biallelic");
                }
                homozygous.insert(index);
            }
        }
    }
    let mendelian_conflicts = find_mendelian_conflicts(trios, table);
    let to_retain: HashSet<usize> = if include_homozygous { (0..table.len()).collect() } else { heterozygous };
    let to_retain: HashSet<usize> =
        to_retain.into_iter().filter(|i| !missing.contains(i) && !mendelian_conflicts.contains(i)).collect();
    let to_discard: HashSet<usize> = (0..table.len()).filter(|i| !to_retain.contains(i)).collect();
    let mut homozygous_positions: Vec<i32> =
        to_retain.intersection(&homozygous).map(|&i| table.variants[i].position).collect();
    homozygous_positions.sort_unstable();
    let mut phasable = table.clone();
    phasable.remove_rows_by_index(&to_discard);
    if family.len() == 1 {
        log::info!(
            "Found {} usable{} variants ({} skipped due to missing genotypes)",
            phasable.len(),
            if include_homozygous { "" } else { " heterozygous" },
            missing.len()
        );
    } else {
        log::info!(
            "Found {} usable variants ({} skipped due to Mendelian conflicts)",
            phasable.len(),
            mendelian_conflicts.len()
        );
    }
    Ok((homozygous_positions, phasable))
}

/// `find_components`: maps each position to its component, identified by the
/// leftmost position.
fn find_components(
    phased_positions: &[i32],
    reads: &ReadSet,
    master_block: Option<&[i32]>,
    heterozygous_positions: Option<&HashMap<i32, HashSet<i32>>>,
) -> HashMap<i32, i32> {
    let mut cf = ComponentFinder::new(phased_positions.iter().copied());
    let set: HashSet<i32> = phased_positions.iter().copied().collect();
    for read in reads.iter() {
        let positions: Vec<i32> = match heterozygous_positions {
            None => read.variants.iter().map(|v| v.position).filter(|p| set.contains(p)).collect(),
            Some(het) => {
                let het = &het[&read.sample_id];
                read.variants.iter().map(|v| v.position).filter(|p| set.contains(p) && het.contains(p)).collect()
            }
        };
        for &p in positions.iter().skip(1) {
            cf.merge(positions[0], p);
        }
    }
    if let Some(master) = master_block {
        for &p in master.iter().skip(1) {
            cf.merge(master[0], p);
        }
    }
    set.iter().map(|&p| (p, cf.find(p))).collect()
}

/// `compute_overall_components`
#[allow(clippy::too_many_arguments)]
fn compute_overall_components(
    accessible_positions: &[i32],
    all_reads: &ReadSet,
    distrust_genotypes: bool,
    family: &[String],
    genetic_haplotyping: bool,
    homozygous_positions: &[i32],
    ids: &NumericSampleIds,
    superreads_list: &[ReadSet],
) -> HashMap<i32, i32> {
    let accessible_set: HashSet<i32> = accessible_positions.iter().copied().collect();
    let mut master_block: Option<Vec<i32>> = None;
    let mut heterozygous_by_sample: Option<HashMap<i32, HashSet<i32>>> = None;
    if distrust_genotypes {
        let mut hom_in_any_sample = HashSet::new();
        let mut het_by_sample = HashMap::new();
        for (sample, sr) in family.iter().zip(superreads_list) {
            let mut hets = HashSet::new();
            for (v1, v2) in sr.get(0).variants.iter().zip(&sr.get(1).variants) {
                assert_eq!(v1.position, v2.position);
                if !accessible_set.contains(&v1.position) {
                    continue;
                }
                match (v1.allele, v2.allele) {
                    (0, 1) | (1, 0) => {
                        hets.insert(v1.position);
                    }
                    (0, 0) | (1, 1) => {
                        hom_in_any_sample.insert(v1.position);
                    }
                    _ => {}
                }
            }
            het_by_sample.insert(ids.id(sample), hets);
        }
        heterozygous_by_sample = Some(het_by_sample);
        if family.len() > 1 && genetic_haplotyping {
            let mut m: Vec<i32> = hom_in_any_sample.into_iter().collect();
            m.sort_unstable();
            master_block = Some(m);
        }
    } else if family.len() > 1 && genetic_haplotyping {
        let mut m: Vec<i32> = homozygous_positions.iter().copied().filter(|p| accessible_set.contains(p)).collect();
        m.sort_unstable();
        m.dedup();
        master_block = Some(m);
    }
    find_components(accessible_positions, all_reads, master_block.as_deref(), heterozygous_by_sample.as_ref())
}

/// `create_pedigree`
fn create_pedigree(
    opts: &PhaseOptions,
    family: &[String],
    ids: &NumericSampleIds,
    phasable: &VariantTable,
    trios: &[Trio],
) -> Result<Pedigree> {
    let mut pedigree = Pedigree::new();
    for sample in family {
        let genotypes: Vec<Genotype> = phasable.genotypes_of(sample).to_vec();
        let gls: Vec<Option<PhredGenotypeLikelihoods>> = if opts.distrust_genotypes {
            let mut v = Vec::with_capacity(genotypes.len());
            for (gt, gl) in genotypes.iter().zip(phasable.genotype_likelihoods_of(sample)) {
                if !gt.is_diploid_and_biallelic() {
                    bail!("AssertionError: genotype {gt} is not diploid and biallelic");
                }
                v.push(Some(match gl {
                    None => {
                        // all genotypes get default_gq, except the called one, which gets 0
                        let mut x = vec![opts.default_gq as f64; 3];
                        x[gt.get_index() as usize] = 0.0;
                        PhredGenotypeLikelihoods::new(x, 2, 2)?
                    }
                    Some(gl) => gl.as_phred(2, opts.gl_regularizer)?,
                }));
            }
            v
        } else {
            vec![None; genotypes.len()]
        };
        pedigree.add_individual(ids.id(sample) as u32, genotypes, gls);
    }
    for trio in trios {
        let (Some(father), Some(mother)) = (&trio.father, &trio.mother) else { continue };
        pedigree.add_relationship(ids.id(father) as u32, ids.id(mother) as u32, ids.id(&trio.child) as u32)?;
    }
    Ok(pedigree)
}

/// `setup_families`: (family representative -> members, representative -> trios)
pub(crate) fn setup_families(
    samples: &[String],
    trios_in_ped: Option<&[Trio]>,
    max_coverage: usize,
) -> (BTreeMap<String, Vec<String>>, HashMap<String, Vec<Trio>>) {
    let mut family_finder = ComponentFinder::new(samples.iter().cloned());
    let mut all_trios = Vec::new();
    if let Some(ped_trios) = trios_in_ped {
        let sample_set: HashSet<&String> = samples.iter().collect();
        let mut warned_unknown = false;
        let mut warned_missing = false;
        for trio in ped_trios {
            let (Some(father), Some(mother)) = (&trio.father, &trio.mother) else {
                if !warned_unknown {
                    log::warn!(
                        "Relationship {}/{}/{} ignored because at least one of the individuals is unknown. Hiding further warnings of this type, use --debug to show",
                        trio.child,
                        trio.father.as_deref().unwrap_or("None"),
                        trio.mother.as_deref().unwrap_or("None")
                    );
                    warned_unknown = true;
                }
                continue;
            };
            if !sample_set.contains(mother) || !sample_set.contains(father) || !sample_set.contains(&trio.child) {
                if !warned_missing {
                    log::warn!(
                        "Relationship {}/{}/{} ignored because at least one of the individuals was not among the samples to be phased (either not in the input VCF or restricted by --sample). Hiding further warnings of this type, use --debug to show",
                        trio.child,
                        father,
                        mother
                    );
                    warned_missing = true;
                }
                continue;
            }
            all_trios.push(trio.clone());
        }
        for trio in &all_trios {
            family_finder.merge(trio.father.clone().unwrap(), trio.child.clone());
            family_finder.merge(trio.mother.clone().unwrap(), trio.child.clone());
        }
    }
    let mut families: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for sample in samples {
        families.entry(family_finder.find(sample.clone())).or_default().push(sample.clone());
    }
    let mut family_trios: HashMap<String, Vec<Trio>> = HashMap::new();
    for trio in &all_trios {
        family_trios.entry(family_finder.find(trio.child.clone())).or_default().push(trio.clone());
    }
    log::info!(
        "Working on {} sample{} from {} famil{}",
        samples.len(),
        if samples.len() == 1 { "" } else { "s" },
        families.len(),
        if families.len() == 1 { "y" } else { "ies" }
    );
    let largest_trio_count = family_trios.values().map(|t| t.len()).max().unwrap_or(0);
    if max_coverage + 2 * largest_trio_count > 23 {
        log::warn!(
            "The maximum coverage is too high! WhatsHap may take a long time to finish and require a huge amount of memory."
        );
    }
    (families, family_trios)
}

/// Phred-scaled recombination costs between consecutive positions.
pub(crate) enum RecombinationCostComputer {
    Uniform(f64),
    GeneticMap(Vec<RecombinationMapEntry>),
}

impl RecombinationCostComputer {
    pub(crate) fn compute(&self, positions: &[i32]) -> Vec<u32> {
        match self {
            RecombinationCostComputer::Uniform(rate) => uniform_recombination_map(*rate, positions),
            RecombinationCostComputer::GeneticMap(map) => recombination_cost_map(map, positions),
        }
    }
}

fn best_case_blocks(reads: &ReadSet) -> (usize, usize) {
    let positions: HashSet<i32> = reads.iter().flat_map(|r| r.variants.iter().map(|v| v.position)).collect();
    let mut cf = ComponentFinder::new(positions.iter().copied());
    for read in reads.iter() {
        for v in read.variants.iter().skip(1) {
            cf.merge(read.variants[0].position, v.position);
        }
    }
    let mut sizes: HashMap<i32, usize> = HashMap::new();
    for &p in &positions {
        *sizes.entry(cf.find(p)).or_default() += 1;
    }
    (sizes.len(), sizes.values().filter(|&&s| s > 1).count())
}

fn log_component_stats(components: &HashMap<i32, i32>, n_accessible: usize) {
    let mut blocks: HashMap<i32, Vec<i32>> = HashMap::new();
    for (&pos, &b) in components {
        blocks.entry(b).or_default().push(pos);
    }
    let mut largest: Vec<i32> = blocks.values().max_by_key(|b| b.len()).cloned().unwrap_or_default();
    largest.sort_unstable();
    if !largest.is_empty() {
        log::info!(
            "Largest block contains {} variants ({:.1}% of accessible variants) between position {} and {}",
            largest.len(),
            100.0 * largest.len() as f64 / n_accessible as f64,
            largest[0] + 1,
            largest[largest.len() - 1] + 1
        );
    } else {
        log::info!("No. of phased blocks: {}", blocks.len());
    }
}

/// Accumulates wall-clock time per stage (`timer.StageTimer`).
#[derive(Default)]
struct StageTimer {
    elapsed: HashMap<&'static str, std::time::Duration>,
}

impl StageTimer {
    fn time<T>(&mut self, stage: &'static str, f: impl FnOnce() -> T) -> T {
        let t = std::time::Instant::now();
        let r = f();
        *self.elapsed.entry(stage).or_default() += t.elapsed();
        r
    }

    fn get(&self, stage: &str) -> f64 {
        self.elapsed.get(stage).map(|d| d.as_secs_f64()).unwrap_or(0.0)
    }
}

/// Immutable state shared by all workers.
struct Shared {
    opts: PhaseOptions,
    bam_paths: Vec<String>,
    families: BTreeMap<String, Vec<String>>,
    family_trios: HashMap<String, Vec<Trio>>,
    ids: NumericSampleIds,
    /// Phased VCFs given as phasing input: chromosome -> variant table
    phase_vcfs: Vec<HashMap<String, VariantTable>>,
    recombination_cost_computer: RecombinationCostComputer,
    /// Threads for allele detection and for solving independent parts of a
    /// chromosome
    pool: Arc<rayon::ThreadPool>,
    /// Threads decompressing each input file
    decompression_threads: usize,
}

/// Per-thread resources (file handles).
struct Worker {
    readset_reader: ReadSetReader,
    fasta: Option<ReferenceCache>,
    timers: StageTimer,
    /// Memory use of the chromosome being phased (if chromosomes are phased
    /// in parallel)
    usage: Option<Arc<ChromosomeUsage>>,
}

impl Worker {
    fn new(shared: &Shared) -> Result<Self> {
        let opts = &shared.opts;
        let fasta = match (&opts.reference, !opts.no_reference) {
            (Some(path), true) => Some(ReferenceCache::open(path)?),
            _ => None,
        };
        let mut readset_reader = ReadSetReader::new(
            &shared.bam_paths,
            opts.reference.as_deref(),
            opts.mapping_quality,
            opts.use_supplementary,
            opts.supplementary_distance_threshold,
            shared.decompression_threads,
        )?;
        if opts.threads > 1 {
            readset_reader.set_detection_pool(Arc::clone(&shared.pool));
        }
        if let Some(f) = &fasta {
            readset_reader.share_reference(f.repository());
        }
        Ok(Worker { readset_reader, fasta, timers: StageTimer::default(), usage: None })
    }

    fn set_usage(&mut self, usage: Option<Arc<ChromosomeUsage>>) {
        self.readset_reader.set_progress(usage.as_ref().map(|u| Arc::clone(&u.progress)));
        self.usage = usage;
    }

    /// `PhasedInputReader.read`: reads of `sample` (from BAM/CRAM and phased
    /// VCF blocks), sorted, and the source ids of the VCF inputs.
    fn read(
        &mut self,
        shared: &Shared,
        chromosome: &str,
        variants: &[crate::vcf::VcfVariant],
        sample: &str,
    ) -> Result<(ReadSet, HashSet<i32>)> {
        let opts = &shared.opts;
        let (reference, copies) = match self.fasta.as_mut() {
            Some(f) => {
                f.get(chromosome)?;
                // CRAM files decode with the cache's copy, which is these bases unless they
                // had to be uppercased.
                let copies = 1 + usize::from(self.readset_reader.cram_files() > 0 && !f.shares_bases());
                (Some(f.get(chromosome)?), copies)
            }
            None => (None, 0),
        };
        if let (Some(usage), Some(reference)) = (&self.usage, reference) {
            usage.set_reference_bytes(reference.len() * copies);
        }
        let bam_sample = if opts.ignore_read_groups { None } else { Some(sample) };
        let numeric_id = match bam_sample {
            None => 0,
            Some(s) => shared.ids.id(s),
        };
        let readset_reader = &mut self.readset_reader;
        let mut readset = self.timers.time("read_bam", || {
            read_alignments(readset_reader, chromosome, variants, bam_sample, numeric_id, reference, None)
        })?;
        let sample_id = shared.ids.id(sample);
        let vcf_source_ids = add_phased_blocks_and_sort(
            &mut readset,
            &shared.phase_vcfs,
            shared.bam_paths.len(),
            chromosome,
            variants,
            sample,
            sample_id,
        )?;
        Ok((readset, vcf_source_ids))
    }
}

/// The alignment part of `PhasedInputReader.read`: the reads of `bam_sample`
/// (all reads if None) with alleles detected; a sample that is not in the
/// files gives no reads, a chromosome that is not in them an error.
/// `restricted_genotypes`: only these alleles are considered per variant.
pub(crate) fn read_alignments(
    readset_reader: &mut ReadSetReader,
    chromosome: &str,
    variants: &[crate::vcf::VcfVariant],
    bam_sample: Option<&str>,
    numeric_id: i32,
    reference: Option<&[u8]>,
    restricted_genotypes: Option<&[Genotype]>,
) -> Result<ReadSet> {
    let result = match restricted_genotypes {
        Some(g) => readset_reader.read_restricted(chromosome, variants, bam_sample, numeric_id, reference, g),
        None => readset_reader.read(chromosome, variants, bam_sample, numeric_id, reference),
    };
    match result {
        Ok(rs) => Ok(rs),
        Err(ReadSetError::SampleNotFound) => {
            log::warn!("Sample '{}' not found in any BAM/CRAM file.", bam_sample.unwrap_or("None"));
            Ok(ReadSet::new())
        }
        Err(ReadSetError::ReferenceNotFound(_)) => {
            let alternative = match chromosome.strip_prefix("chr") {
                Some(c) => c.to_string(),
                None => format!("chr{chromosome}"),
            };
            let mut message = format!("The chromosome '{chromosome}' was not found in the BAM/CRAM file.");
            if readset_reader.has_reference(&alternative) {
                message += &format!(" Found '{alternative}' instead");
            }
            bail!(message);
        }
        Err(ReadSetError::Other(e)) => Err(e),
    }
}

/// The rest of `PhasedInputReader.read`: adds the phased blocks of the VCF
/// inputs for `sample` and sorts the reads; returns the VCF source ids.
pub(crate) fn add_phased_blocks_and_sort(
    readset: &mut ReadSet,
    phase_vcfs: &[HashMap<String, VariantTable>],
    n_bam_paths: usize,
    chromosome: &str,
    variants: &[crate::vcf::VcfVariant],
    sample: &str,
    sample_id: i32,
) -> Result<HashSet<i32>> {
    let mut vcf_source_ids = HashSet::new();
    for (i, vcf) in phase_vcfs.iter().enumerate() {
        if let Some(table) = vcf.get(chromosome) {
            let source_id = (n_bam_paths + i) as i32;
            vcf_source_ids.insert(source_id);
            for read in table.phased_blocks_as_reads(sample, variants, source_id, sample_id)? {
                readset.add(read)?;
            }
        }
    }
    for read in readset.iter_mut() {
        read.sort()?;
    }
    readset.sort();
    log::info!("Found {} reads covering {} variants", readset.len(), readset.get_positions().len());
    Ok(vcf_source_ids)
}

/// Phasing result for one chromosome (one VCF block).
struct ChromResult {
    chromosome: String,
    superreads: Vec<(String, ReadSet)>,
    components: HashMap<String, HashMap<i32, i32>>,
    read_list_lines: Vec<String>,
    /// Contents of the recombination list file, once per family (each
    /// overwrites the previous one, as in WhatsHap)
    recombination_lists: Vec<String>,
}

fn phase_chromosome(shared: &Shared, worker: &mut Worker, variant_table: VariantTable) -> Result<ChromResult> {
    let opts = &shared.opts;
    let ids = &shared.ids;
    let chromosome = variant_table.chromosome.clone();
    let mut superreads: Vec<(String, ReadSet)> = Vec::new();
    let mut components: HashMap<String, HashMap<i32, i32>> = HashMap::new();
    let mut read_list_lines = Vec::new();
    let mut recombination_lists = Vec::new();
    let no_trios: Vec<Trio> = Vec::new();
    for (representative, family) in &shared.families {
        log::info!("");
        if family.len() == 1 {
            log::info!("# Working on contig {chromosome} in individual {representative}");
        } else {
            log::info!("# Working on contig {chromosome} in family individuals {}", family.join(","));
        }
        let max_coverage_per_sample = 1.max(opts.max_coverage / family.len());
        let trios = shared.family_trios.get(representative).unwrap_or(&no_trios);
        assert!(family.len() == 1 || !trios.is_empty());
        let (homozygous_positions, mut phasable) =
            find_phaseable_variants(family, opts.include_homozygous, trios, &variant_table)?;
        let mut readsets: Vec<ReadSet> = Vec::new();
        for sample in family {
            let (readset, vcf_source_ids) = worker.read(shared, &chromosome, &phasable.variants, sample)?;
            let keep: Vec<usize> = (0..readset.len()).filter(|&i| readset.get(i).len() >= 2).collect();
            let readset = readset.into_subset(keep);
            log::info!("Kept {} reads that cover at least two variants each", readset.len());
            let merged;
            let to_select = if opts.read_merging {
                merged = worker.timers.time("select", || ReadMerger::default().merge(&readset))?;
                &merged
            } else {
                &readset
            };
            let selected = if matches!(opts.algorithm, Algorithm::Heuristic { .. }) {
                // the heuristic phases all reads
                to_select.clone()
            } else {
                let selected_indices = worker.timers.time("select", || {
                    readselection(to_select, max_coverage_per_sample, Some(&vcf_source_ids), true)
                })?;
                to_select.subset(selected_indices)
            };
            log::info!(
                "Selected {} most phase-informative reads covering {} variants",
                selected.len(),
                selected.get_positions().len()
            );
            if family.len() == 1 && !opts.distrust_genotypes {
                let (_, nonsingleton_all) = best_case_blocks(&readset);
                let (n_cov, nonsingleton_cov) = best_case_blocks(&selected);
                log::info!(
                    "Best-case phasing would result in {} non-singleton phased block{} ({} singletons). ",
                    nonsingleton_cov,
                    if nonsingleton_cov == 1 { "" } else { "s" },
                    n_cov - nonsingleton_cov
                );
                log::debug!("... would be {nonsingleton_all} non-singleton phased blocks without read selection");
            }
            readsets.push(selected);
        }
        // merge_readsets
        let mut all_reads = ReadSet::new();
        for rs in readsets {
            for r in rs.into_reads() {
                all_reads.add(r)?;
            }
        }
        all_reads.sort();
        let mut accessible_positions = all_reads.get_positions();
        if family.len() > 1 && opts.genetic_haplotyping {
            // also retain positions homozygous in at least one individual
            accessible_positions.extend(homozygous_positions.iter().copied());
            accessible_positions.sort_unstable();
            accessible_positions.dedup();
            log::info!(
                "Variants either covered by phase-informative read or homozygous in at least one individual: {}",
                accessible_positions.len()
            );
        }
        let accessible_set: HashSet<i32> = accessible_positions.iter().copied().collect();
        phasable.subset_rows_by_position(&accessible_set);
        assert_eq!(phasable.len(), accessible_positions.len());

        let pedigree = create_pedigree(opts, family, ids, &phasable, trios)?;
        let recombination_costs = shared.recombination_cost_computer.compute(&accessible_positions);
        let problem_name = if family.len() == 1 { "MEC" } else { "PedMEC" };
        log::info!(
            "Phasing {} sample{} by solving the {problem_name} problem ...",
            family.len(),
            if family.len() == 1 { "" } else { "s" }
        );
        all_reads.reassign_read_ids();
        let solution = worker.timers.time("phase", || match opts.algorithm {
            Algorithm::HapChat => hapchat::solve(&all_reads),
            Algorithm::Heuristic { row_limit } => heuristic::solve(
                &all_reads,
                &recombination_costs,
                &pedigree,
                opts.distrust_genotypes,
                &accessible_positions,
                row_limit,
            ),
            Algorithm::WhatsHap => shared.pool.install(|| {
                dp::solve(
                    &all_reads,
                    &recombination_costs,
                    &pedigree,
                    opts.distrust_genotypes,
                    Some(&accessible_positions),
                )
            }),
        })?;
        let dp::Solution { optimal_score, superreads: superreads_list, transmission_vector, partitioning } = solution;
        log::debug!("{problem_name} cost: {optimal_score}");
        let overall_components = worker.timers.time("components", || {
            compute_overall_components(
                &accessible_positions,
                &all_reads,
                opts.distrust_genotypes,
                family,
                opts.genetic_haplotyping,
                &homozygous_positions,
                ids,
                &superreads_list,
            )
        });
        log_component_stats(&overall_components, accessible_positions.len());

        if opts.recombination_list_filename.is_some() {
            if opts.algorithm == Algorithm::HapChat {
                // `assert transmission_vector is not None` in WhatsHap
                bail!("HapChat computes no transmission vector for --recombination-list");
            }
            let (text, n) = recombination_list(
                &chromosome,
                &accessible_positions,
                &overall_components,
                &recombination_costs,
                &transmission_vector,
                trios,
            )?;
            log::info!("Total no. of detected recombination events: {n}");
            recombination_lists.push(text);
        }

        if opts.read_list_filename.is_some() {
            for (i, read) in all_reads.iter().enumerate() {
                // HapChat does not compute a partitioning (WhatsHap writes "*")
                let haplotype = partitioning.get(i).map_or_else(|| "*".to_string(), |h| h.to_string());
                let phaseset = overall_components[&read.variants[0].position] + 1;
                read_list_lines.push(format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read.name,
                    read.source_id,
                    ids.name(read.sample_id),
                    phaseset,
                    haplotype,
                    read.len(),
                    read.variants[0].position + 1,
                    read.variants[read.len() - 1].position + 1
                ));
            }
        }
        for (sample, sr) in family.iter().zip(superreads_list) {
            // (an assertion in WhatsHap; the heuristic orders samples by id)
            if sr.get(0).sample_id != ids.id(sample) || sr.get(1).sample_id != ids.id(sample) {
                bail!("The super reads of sample {sample} belong to another sample");
            }
            superreads.push((sample.clone(), sr));
            components.insert(sample.clone(), overall_components.clone());
        }
    }
    Ok(ChromResult { chromosome, superreads, components, read_list_lines, recombination_lists })
}

/// `write_recombination_list`: file contents and number of events. (Python's
/// `print` without `sep` separates fields by spaces.)
fn recombination_list(
    chromosome: &str,
    accessible_positions: &[i32],
    overall_components: &HashMap<i32, i32>,
    recombination_costs: &[u32],
    transmission_vector: &[u32],
    trios: &[Trio],
) -> Result<(String, usize)> {
    let mut per_child: HashMap<&str, Vec<u32>> = HashMap::new();
    for &value in transmission_vector {
        let mut v = value;
        for trio in trios {
            per_child.entry(trio.child.as_str()).or_default().push(v % 4);
            v /= 4;
        }
    }
    let mut text = String::from(
        "#child_id chromosome position1 position2 transmitted_hap_father1 transmitted_hap_father2 transmitted_hap_mother1 transmitted_hap_mother2 recombination_cost\n",
    );
    let mut n = 0;
    for trio in trios {
        let tv = per_child.get(trio.child.as_str()).map(|v| v.as_slice()).unwrap_or(&[]);
        let events = find_recombination(tv, overall_components, accessible_positions, recombination_costs)?;
        for e in &events {
            text += &format!(
                "{} {} {} {} {} {} {} {} {}\n",
                trio.child,
                chromosome,
                e.position1 + 1,
                e.position2 + 1,
                e.transmitted_hap_father1,
                e.transmitted_hap_father2,
                e.transmitted_hap_mother1,
                e.transmitted_hap_mother2,
                e.recombination_cost
            );
        }
        n += events.len();
    }
    Ok((text, n))
}

/// Writes per-chromosome results in VCF order.
struct OrderedOutput {
    vcf_writer: PhasedVcfWriter,
    read_list: Option<std::io::BufWriter<std::fs::File>>,
    recombination_list_filename: Option<String>,
    gtchange_list_filename: Option<String>,
    distrust_genotypes: bool,
    timers: StageTimer,
}

/// `write_changed_genotypes`
fn write_changed_genotypes(path: &str, changes: &[GenotypeChange]) -> Result<()> {
    let mut text = String::from("#sample\tchromosome\tposition\tREF\tALT\told_gt\tnew_gt\n");
    for c in changes {
        text += &format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            c.sample, c.chromosome, c.position, c.reference, c.alternative, c.old_gt, c.new_gt
        );
    }
    std::fs::write(path, text)?;
    Ok(())
}

impl OrderedOutput {
    /// `result` is None for chromosomes that are left unchanged.
    fn write(&mut self, chromosome: &str, result: Option<ChromResult>) -> Result<()> {
        let vcf_writer = &mut self.vcf_writer;
        match result {
            None => {
                self.timers.time("write_vcf", || vcf_writer.write(chromosome, &[], &HashMap::new()))?;
            }
            Some(r) => {
                if let Some(path) = &self.recombination_list_filename {
                    for text in &r.recombination_lists {
                        std::fs::write(path, text)?;
                    }
                }
                log::debug!("Writing phasing result to output VCF");
                let changes =
                    self.timers.time("write_vcf", || vcf_writer.write(&r.chromosome, &r.superreads, &r.components))?;
                if !changes.is_empty() {
                    // `assert distrust_genotypes` in WhatsHap (HapChat can change
                    // genotypes after its first block of reads)
                    if !self.distrust_genotypes {
                        bail!("Genotypes changed while writing the VCF without --distrust-genotypes");
                    }
                    log::info!("Changed {} genotypes while writing VCF", changes.len());
                }
                if let Some(path) = &self.gtchange_list_filename {
                    log::info!("Writing list of changed genotypes to '{path}'");
                    write_changed_genotypes(path, &changes)?;
                }
                if let Some(f) = self.read_list.as_mut() {
                    for line in &r.read_list_lines {
                        writeln!(f, "{line}")?;
                    }
                }
                log::debug!("Chromosome '{}' finished", r.chromosome);
            }
        }
        Ok(())
    }
}

pub fn run_phase(opts: PhaseOptions) -> Result<()> {
    let start_time = std::time::Instant::now();
    log::info!("This is lrpipe-phase {} (port of WhatsHap 2.8)", env!("CARGO_PKG_VERSION"));
    if let Some(reference) = opts.reference.as_deref().filter(|_| !opts.no_reference) {
        crate::fasta::check_reference(reference)?;
    }

    let mut bam_paths = Vec::new();
    let mut vcf_paths = Vec::new();
    for path in &opts.phase_input_files {
        match detect_file_format(path)? {
            Some("BAM") | Some("CRAM") => bam_paths.push(path.clone()),
            Some("VCF") => vcf_paths.push(path.clone()),
            _ => bail!("Unable to determine type of input file '{path}'"),
        }
    }
    let mut phase_vcf_readers = Vec::new();
    for path in &vcf_paths {
        phase_vcf_readers.push(VcfReader::open_with(path, opts.only_snvs, true, None)?);
    }
    if !bam_paths.is_empty() && opts.reference.is_none() && !opts.no_reference {
        bail!(
            "A reference FASTA needs to be provided with -r/--reference; or use --no-reference at the expense of phasing quality."
        );
    }
    let vcf_writer =
        PhasedVcfWriter::new(&opts.variant_file, opts.command_line.as_deref(), &opts.output, opts.tag, opts.only_snvs)?;
    let mut vcf_reader =
        VcfReader::open(&opts.variant_file, opts.only_snvs)?.with_genotype_likelihoods(opts.distrust_genotypes);
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
    let (families, family_trios) = setup_families(&samples, ped_trios.as_deref(), opts.max_coverage);
    // WhatsHap assigns numeric sample ids on first use: children of trios
    // first, then in the order in which families are processed.
    let mut ids = NumericSampleIds::default();
    for trios in family_trios.values() {
        for t in trios {
            ids.get(&t.child);
        }
    }
    for family in families.values() {
        for s in family {
            ids.get(s);
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
    let read_list = match &opts.read_list_filename {
        Some(path) => {
            let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
            writeln!(
                f,
                "#readname\tsource_id\tsample\tphaseset\thaplotype\tcovered_variants\tfirst_variant_pos\tlast_variant_pos"
            )?;
            Some(f)
        }
        None => None,
    };
    let threads = opts.threads.max(1);
    let recombination_list_filename = opts.recombination_list_filename.clone();
    let gtchange_list_filename = opts.gtchange_list_filename.clone();
    let distrust_genotypes = opts.distrust_genotypes;
    dp::set_table_memory_budget(opts.dp_memory_mb << 20);
    let pool = Arc::new(rayon::ThreadPoolBuilder::new().num_threads(threads).build()?);
    let decompression_threads = if threads > 1 { threads.min(DECOMPRESSION_THREADS) } else { 1 };
    let shared = Shared {
        opts,
        bam_paths,
        families,
        family_trios,
        ids,
        phase_vcfs,
        recombination_cost_computer,
        pool,
        decompression_threads,
    };
    let mut out = OrderedOutput {
        vcf_writer,
        read_list,
        recombination_list_filename,
        gtchange_list_filename,
        distrust_genotypes,
        timers: StageTimer::default(),
    };
    let mut parse_timer = StageTimer::default();
    let included = |chromosome: &str| {
        (shared.opts.chromosomes.is_empty() || shared.opts.chromosomes.iter().any(|c| c == chromosome))
            && !shared.opts.excluded_chromosomes.iter().any(|c| c == chromosome)
    };

    let mut worker_timers: Vec<StageTimer> = Vec::new();
    if threads == 1 {
        let mut worker = Worker::new(&shared)?;
        while let Some(table) = parse_timer.time("parse_vcf", || vcf_reader.next_table())? {
            let chromosome = table.chromosome.clone();
            if !included(&chromosome) {
                log::info!(
                    "Leaving chromosome '{chromosome}' unchanged (present in VCF but not requested by --chromosome)"
                );
                out.write(&chromosome, None)?;
                continue;
            }
            let result = phase_chromosome(&shared, &mut worker, table)?;
            out.write(&chromosome, Some(result))?;
        }
        worker_timers.push(worker.timers);
    } else {
        // Chromosomes are phased in parallel; results are written in VCF order.
        use std::sync::Mutex;
        use std::sync::mpsc::{channel, sync_channel};
        // Chromosomes phased at the same time. Allele detection, the DP and
        // decompression use all threads, so half as many chromosomes as
        // threads keeps the cores busy (on HG003 it was as fast as one per
        // thread) while holding fewer chromosomes in memory.
        let workers = threads.div_ceil(2);
        // Fewer are phased at a time if their estimated memory does not fit
        // next to the DP tables (which have their own budget) and the memory
        // that grows with the number of threads (stacks, allocator caches,
        // alignments queued for allele detection, decompression buffers:
        // about 1.4 GB with 16 threads on HG003 ONT).
        let max_memory = match shared.opts.max_memory_mb {
            Some(mb) => mb << 20,
            None => memory::system_memory().unwrap_or(16 << 30) / 2,
        };
        let budget = max_memory
            .saturating_sub(shared.opts.dp_memory_mb << 20)
            .saturating_sub(threads as u64 * memory::THREAD_MEMORY);
        let mut scheduler = MemoryScheduler::new(budget);
        log::info!(
            "Phasing up to {workers} chromosomes at a time while their estimated memory fits into {:.1} GB",
            max_memory as f64 / f64::from(1 << 30)
        );
        // Opened here so that an error (such as a missing reference index)
        // is reported before any thread waits for chromosomes.
        let pending_workers = (0..workers).map(|_| Worker::new(&shared)).collect::<Result<Vec<_>>>()?;
        let (task_tx, task_rx) = sync_channel::<(usize, VariantTable, Arc<ChromosomeUsage>)>(workers);
        let task_rx = Arc::new(Mutex::new(task_rx));
        let (res_tx, res_rx) = channel::<(usize, Result<ChromResult>)>();
        let shared_ref = &shared;
        std::thread::scope(|scope| -> Result<()> {
            let mut handles = Vec::new();
            for mut worker in pending_workers {
                let task_rx = Arc::clone(&task_rx);
                let res_tx = res_tx.clone();
                handles.push(scope.spawn(move || -> StageTimer {
                    loop {
                        let job = task_rx.lock().unwrap().recv();
                        let Ok((i, table, usage)) = job else { break };
                        let finished = FinishOnDrop(Arc::clone(&usage));
                        worker.set_usage(Some(usage));
                        let r = phase_chromosome(shared_ref, &mut worker, table);
                        worker.set_usage(None);
                        drop(finished);
                        let failed = r.is_err();
                        if res_tx.send((i, r)).is_err() || failed {
                            break;
                        }
                    }
                    worker.timers
                }));
            }
            // Sending fails once all workers have stopped.
            drop(task_rx);
            drop(res_tx);
            // order of blocks: (chromosome, phased?)
            let mut blocks: Vec<(String, bool)> = Vec::new();
            let mut pending: BTreeMap<usize, ChromResult> = BTreeMap::new();
            let mut next_to_write = 0usize;
            let mut flush = |blocks: &Vec<(String, bool)>,
                             pending: &mut BTreeMap<usize, ChromResult>,
                             next: &mut usize|
             -> Result<()> {
                while *next < blocks.len() {
                    let (chrom, phased) = &blocks[*next];
                    if !*phased {
                        out.write(chrom, None)?;
                    } else if let Some(r) = pending.remove(next) {
                        out.write(chrom, Some(r))?;
                    } else {
                        break;
                    }
                    *next += 1;
                }
                Ok(())
            };
            let mut error: Option<anyhow::Error> = None;
            fn handle(
                (i, r): (usize, Result<ChromResult>),
                pending: &mut BTreeMap<usize, ChromResult>,
                error: &mut Option<anyhow::Error>,
            ) {
                match r {
                    Ok(res) => {
                        pending.insert(i, res);
                    }
                    Err(e) => {
                        if error.is_none() {
                            *error = Some(e);
                        }
                    }
                }
            }
            loop {
                let table = match parse_timer.time("parse_vcf", || vcf_reader.next_table()) {
                    Ok(Some(t)) => t,
                    Ok(None) => break,
                    Err(e) => {
                        error = Some(e);
                        break;
                    }
                };
                let i = blocks.len();
                let chromosome = table.chromosome.clone();
                if !included(&chromosome) {
                    log::info!(
                        "Leaving chromosome '{chromosome}' unchanged (present in VCF but not requested by --chromosome)"
                    );
                    blocks.push((chromosome, false));
                } else {
                    blocks.push((chromosome, true));
                    // Write results while waiting for memory.
                    while !scheduler.fits(table.variants.len()) && error.is_none() {
                        while let Ok(msg) = res_rx.try_recv() {
                            handle(msg, &mut pending, &mut error);
                        }
                        flush(&blocks, &mut pending, &mut next_to_write)?;
                        std::thread::sleep(std::time::Duration::from_millis(10));
                    }
                    if error.is_some() {
                        break;
                    }
                    let usage = Arc::new(ChromosomeUsage::new(table.variants.len()));
                    if task_tx.send((i, table, Arc::clone(&usage))).is_err() {
                        break;
                    }
                    scheduler.start(usage);
                }
                while let Ok(msg) = res_rx.try_recv() {
                    handle(msg, &mut pending, &mut error);
                }
                if error.is_some() {
                    break;
                }
                flush(&blocks, &mut pending, &mut next_to_write)?;
            }
            drop(task_tx);
            if error.is_none() {
                for msg in res_rx.iter() {
                    handle(msg, &mut pending, &mut error);
                    if error.is_some() {
                        break;
                    }
                    flush(&blocks, &mut pending, &mut next_to_write)?;
                }
                if error.is_none() {
                    flush(&blocks, &mut pending, &mut next_to_write)?;
                }
            }
            for h in handles {
                worker_timers.push(h.join().expect("worker thread panicked"));
            }
            match error {
                Some(e) => Err(e),
                None => Ok(()),
            }
        })?;
    }
    let OrderedOutput { vcf_writer, read_list, timers: out_timers, .. } = out;
    vcf_writer.close()?;
    if let Some(mut f) = read_list {
        f.flush()?;
    }
    let total = start_time.elapsed().as_secs_f64();
    let get = |stage: &str| -> f64 {
        worker_timers.iter().map(|t| t.get(stage)).sum::<f64>() + out_timers.get(stage) + parse_timer.get(stage)
    };
    log::info!("\n# Resource usage");
    if threads > 1 {
        log::info!("(stage times are summed over {threads} threads)");
    }
    log::info!("Time spent reading BAM/CRAM:                 {:6.1} s", get("read_bam"));
    log::info!("Time spent parsing VCF:                      {:6.1} s", get("parse_vcf"));
    log::info!("Time spent selecting reads:                  {:6.1} s", get("select"));
    log::info!("Time spent phasing:                          {:6.1} s", get("phase"));
    log::info!("Time spent writing VCF:                      {:6.1} s", get("write_vcf"));
    log::info!("Time spent finding components:               {:6.1} s", get("components"));
    log::info!("Total elapsed time:                          {:6.1} s", total);
    Ok(())
}

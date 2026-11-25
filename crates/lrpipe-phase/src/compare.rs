//! Port of `cli/compare.py`: compare two or more phased variant files.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result, bail};

use crate::compat::pyfloat::PyNum;
use crate::core::genotype::Genotype;
use crate::core::switchflip::SwitchFlipCalculator;
use crate::vcf::{BlockKey, PloidyError, VariantCallPhase, VariantTable, VcfReader, VcfVariant, block_key};

const COUNT_WIDTH: usize = 9;

pub struct CompareOptions {
    pub vcf: Vec<String>,
    pub ploidy: usize,
    pub names: Option<String>,
    pub sample: Option<String>,
    pub ignore_sample_name: bool,
    pub tsv_pairwise: Option<String>,
    pub tsv_multiway: Option<String>,
    pub only_snvs: bool,
    pub switch_error_bed: Option<String>,
    pub plot_blocksizes: Option<String>,
    pub plot_sum_of_blocksizes: Option<String>,
    pub longest_block_tsv: Option<String>,
}

/// `SwitchFlips`: counts are ints for diploid comparisons and floats
/// (divided by the ploidy) for polyploid ones.
#[derive(Clone, Copy, Debug)]
struct SwitchFlips {
    switches: PyNum,
    flips: PyNum,
}

impl Default for SwitchFlips {
    fn default() -> Self {
        SwitchFlips { switches: PyNum::Int(0), flips: PyNum::Int(0) }
    }
}

impl std::fmt::Display for SwitchFlips {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}", self.switches, self.flips)
    }
}

#[derive(Clone, Copy, Debug)]
struct PhasingErrors {
    switches: PyNum,
    hamming: PyNum,
    switch_flips: SwitchFlips,
    diff_genotypes: i64,
}

impl Default for PhasingErrors {
    fn default() -> Self {
        PhasingErrors {
            switches: PyNum::Int(0),
            hamming: PyNum::Int(0),
            switch_flips: SwitchFlips::default(),
            diff_genotypes: 0,
        }
    }
}

impl std::ops::AddAssign for PhasingErrors {
    fn add_assign(&mut self, other: PhasingErrors) {
        self.switches += other.switches;
        self.hamming += other.hamming;
        self.switch_flips.switches += other.switch_flips.switches;
        self.switch_flips.flips += other.switch_flips.flips;
        self.diff_genotypes += other.diff_genotypes;
    }
}

/// A haplotype of a block as Python's string of allele numbers
type Hap = Vec<u8>;

fn complement(s: &[u8]) -> Result<Vec<u8>> {
    s.iter()
        .map(|&c| match c {
            b'0' => Ok(b'1'),
            b'1' => Ok(b'0'),
            _ => bail!("KeyError: '{}'", c as char),
        })
        .collect()
}

/// Number of differing positions of two equally long sequences
fn hamming<T: PartialEq>(s0: &[T], s1: &[T]) -> Result<i64> {
    if s0.len() != s1.len() {
        bail!("AssertionError: hamming distance of sequences of different lengths");
    }
    Ok(s0.iter().zip(s1).filter(|(a, b)| a != b).count() as i64)
}

fn switch_encoding(phasing: &[u8]) -> Vec<u8> {
    phasing.windows(2).map(|w| if w[0] == w[1] { b'0' } else { b'1' }).collect()
}

fn compute_switch_flips(phasing0: &[u8], phasing1: &[u8]) -> Result<SwitchFlips> {
    if phasing0.len() != phasing1.len() {
        bail!("AssertionError: phasings of different lengths");
    }
    let s0 = switch_encoding(phasing0);
    let s1 = switch_encoding(phasing1);
    let (mut switches, mut flips) = (0, 0);
    let mut switches_in_a_row = 0;
    for (i, (p0, p1)) in s0.iter().zip(&s1).enumerate() {
        if p0 != p1 {
            switches_in_a_row += 1;
        }
        if i + 1 == s0.len() || p0 == p1 {
            flips += switches_in_a_row / 2;
            switches += switches_in_a_row % 2;
            switches_in_a_row = 0;
        }
    }
    Ok(SwitchFlips { switches: PyNum::Int(switches), flips: PyNum::Int(flips) })
}

/// `hap[i]` of a Python string
fn char_at(hap: &[u8], i: usize) -> Result<u8> {
    hap.get(i).copied().context("IndexError: string index out of range")
}

/// The digit `int(hap[i])`
fn digit_at(hap: &[u8], i: usize) -> Result<u32> {
    Ok(u32::from(char_at(hap, i)? - b'0'))
}

/// Positions at which both phasings agree on the genotype
fn compute_matching_genotype_pos(phasing0: &[Hap], phasing1: &[Hap]) -> Result<Vec<usize>> {
    if phasing0.len() != phasing1.len()
        || phasing0.len() < 2
        || phasing0[0].len() != phasing1[0].len()
        || phasing0.iter().any(|h| h.len() != phasing0[0].len())
    {
        bail!("AssertionError: incompatible phasings");
    }
    let mut matching = Vec::new();
    for i in 0..phasing0[0].len() {
        let g0: Vec<u32> = phasing0.iter().map(|h| digit_at(h, i)).collect::<Result<_>>()?;
        let g1: Vec<u32> = phasing1.iter().map(|h| digit_at(h, i)).collect::<Result<_>>()?;
        if Genotype::new(&g0)? == Genotype::new(&g1)? {
            matching.push(i);
        }
    }
    Ok(matching)
}

/// `compute_switch_flips_poly`
fn compute_switch_flips_poly(
    phasing0: &[Hap],
    phasing1: &[Hap],
    switch_cost: f64,
    flip_cost: f64,
) -> Result<SwitchFlips> {
    if phasing0.len() != phasing1.len() {
        log::error!(
            "Incompatible phasings. Number of haplotypes is not equal ({} != {}).",
            phasing0.len(),
            phasing1.len()
        );
        bail!("AssertionError: incompatible phasings");
    }
    let num_pos = phasing0[0].len();
    if num_pos == 0 {
        return Ok(SwitchFlips::default());
    }
    let ploidy = phasing0.len();
    for (i, hap) in phasing1.iter().enumerate() {
        if hap.len() != num_pos {
            log::error!(
                "Inconsistent input for phasing. Haplotypes have different lengths ( len(phasing1[0]={num_pos} != len(phasing1[{i}]={}.",
                hap.len()
            );
            bail!("AssertionError: haplotypes of different lengths");
        }
        if phasing0[i].len() != num_pos {
            log::error!(
                "Inconsistent input for phasing. Haplotypes have different lengths ( len(phasing1[0]={num_pos} != len(phasing0[{i}]={}.",
                phasing0[i].len()
            );
        }
    }
    if ploidy > 6 {
        log::warn!("Computing vector error with more than 6 haplotypes. This may take very long ...");
    }
    let position_wise = |phasing: &[Hap]| -> Result<Vec<Vec<u32>>> {
        (0..num_pos).map(|i| phasing.iter().map(|h| digit_at(h, i)).collect()).collect()
    };
    let calc = SwitchFlipCalculator::new(ploidy as u32, switch_cost, flip_cost);
    let r = calc.compare(&position_wise(phasing0)?, &position_wise(phasing1)?);
    Ok(SwitchFlips { switches: PyNum::Float(r.switches / ploidy as f64), flips: PyNum::Float(r.flips / ploidy as f64) })
}

/// `compute_switch_errors_poly`: switches needed, ignoring positions with
/// different genotypes
fn compute_switch_errors_poly(phasing0: &[Hap], phasing1: &[Hap], matching_pos: &[usize]) -> Result<PyNum> {
    let num_vars = phasing0[0].len();
    let matched = |phasing: &[Hap]| -> Result<Vec<Hap>> {
        phasing.iter().map(|h| matching_pos.iter().map(|&i| char_at(h, i)).collect()).collect()
    };
    let vector_error = compute_switch_flips_poly(
        &matched(phasing0)?,
        &matched(phasing1)?,
        1.0,
        (2 * num_vars * phasing0.len() + 1) as f64,
    )?;
    if vector_error.flips.as_f64() != 0.0 {
        bail!("AssertionError: flips in switch error computation");
    }
    Ok(vector_error.switches)
}

/// Itertools' `permutations(range(n))`, in lexicographic order
fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut cur: Vec<usize> = (0..n).collect();
    let mut result = vec![cur.clone()];
    while let Some(i) = (1..n).rev().find(|&i| cur[i - 1] < cur[i]) {
        let j = (i..n).rev().find(|&j| cur[j] > cur[i - 1]).expect("a larger element exists");
        cur.swap(i - 1, j);
        cur[i..].reverse();
        result.push(cur.clone());
    }
    result
}

/// `compare_block`: phasings are lists of haplotype strings
fn compare_block(phasing0: &[Hap], phasing1: &[Hap]) -> Result<PhasingErrors> {
    let ploidy = phasing0.len();
    let mut minimum_hamming_distance = f64::INFINITY;
    for perm in permutations(ploidy) {
        let mut total_hamming = 0;
        for i in 0..ploidy {
            total_hamming += hamming(&phasing1[i], &phasing0[perm[i]])?;
        }
        let total_hamming = total_hamming as f64 / ploidy as f64;
        minimum_hamming_distance = minimum_hamming_distance.min(total_hamming);
    }
    let matching_pos = compute_matching_genotype_pos(phasing0, phasing1)?;
    let (switches, switch_flips, hamming_distance) = if ploidy == 2 {
        (
            PyNum::Int(hamming(&switch_encoding(&phasing0[0]), &switch_encoding(&phasing1[0]))?),
            compute_switch_flips(&phasing0[0], &phasing1[0])?,
            PyNum::Int(minimum_hamming_distance as i64),
        )
    } else {
        (
            compute_switch_errors_poly(phasing0, phasing1, &matching_pos)?,
            compute_switch_flips_poly(phasing0, phasing1, 1.0, 1.0)?,
            PyNum::Float(minimum_hamming_distance),
        )
    };
    Ok(PhasingErrors {
        switches,
        hamming: hamming_distance,
        switch_flips,
        diff_genotypes: (phasing0[0].len() - matching_pos.len()) as i64,
    })
}

fn fraction2percentstr(nominator: PyNum, denominator: i64) -> String {
    if denominator == 0 { "--".to_string() } else { format!("{:.2}%", nominator.as_f64() * 100.0 / denominator as f64) }
}

fn safefraction(nominator: PyNum, denominator: i64) -> PyNum {
    PyNum::Float(if denominator == 0 { f64::NAN } else { nominator.as_f64() / denominator as f64 })
}

/// Python's `s.rjust(width)`
fn rjust(s: &str, width: usize) -> String {
    let n = s.chars().count();
    if n >= width { s.to_string() } else { format!("{}{s}", " ".repeat(width - n)) }
}

/// A value printed by `print_stat`
enum Stat {
    None,
    Dash,
    Value(String),
}

impl From<i64> for Stat {
    fn from(v: i64) -> Stat {
        Stat::Value(v.to_string())
    }
}

impl From<usize> for Stat {
    fn from(v: usize) -> Stat {
        Stat::Value(v.to_string())
    }
}

impl From<PyNum> for Stat {
    fn from(v: PyNum) -> Stat {
        Stat::Value(v.py_str())
    }
}

impl From<String> for Stat {
    fn from(v: String) -> Stat {
        Stat::Value(v)
    }
}

impl From<SwitchFlips> for Stat {
    fn from(v: SwitchFlips) -> Stat {
        Stat::Value(v.to_string())
    }
}

/// `print_stat`: "text: value [value2]" with the text right-aligned
fn print_stat(out: &mut impl Write, text: &str, value: Stat, value2: Option<String>) -> Result<()> {
    let text = rjust(text, 37);
    match value {
        Stat::None => writeln!(out, "{text}")?,
        Stat::Dash => writeln!(out, "{text}: {}", "-".repeat(COUNT_WIDTH))?,
        Stat::Value(v) => match value2 {
            None => writeln!(out, "{text}: {}", rjust(&v, COUNT_WIDTH))?,
            Some(v2) => writeln!(out, "{text}: {} {}", rjust(&v, COUNT_WIDTH), rjust(&v2, COUNT_WIDTH))?,
        },
    }
    Ok(())
}

fn stat(out: &mut impl Write, text: &str, value: impl Into<Stat>) -> Result<()> {
    print_stat(out, text, value.into(), None)
}

fn print_errors(out: &mut impl Write, errors: &PhasingErrors, phased_pairs: i64) -> Result<()> {
    stat(out, "phased pairs of variants assessed", phased_pairs)?;
    stat(out, "switch errors", errors.switches)?;
    stat(out, "switch error rate", fraction2percentstr(errors.switches, phased_pairs))?;
    stat(out, "switch/flip decomposition", errors.switch_flips)?;
    stat(
        out,
        "switch/flip rate",
        fraction2percentstr(errors.switch_flips.switches + errors.switch_flips.flips, phased_pairs),
    )
}

/// `PairwiseComparisonResults`, as printed to the TSV
struct PairwiseComparisonResults {
    values: Vec<String>,
}

const PAIRWISE_FIELDS: [&str; 20] = [
    "intersection_blocks",
    "covered_variants",
    "all_assessed_pairs",
    "all_switches",
    "all_switch_rate",
    "all_switchflips",
    "all_switchflip_rate",
    "blockwise_hamming",
    "blockwise_hamming_rate",
    "blockwise_diff_genotypes",
    "blockwise_diff_genotypes_rate",
    "largestblock_assessed_pairs",
    "largestblock_switches",
    "largestblock_switch_rate",
    "largestblock_switchflips",
    "largestblock_switchflip_rate",
    "largestblock_hamming",
    "largestblock_hamming_rate",
    "largestblock_diff_genotypes",
    "largestblock_diff_genotypes_rate",
];

type BedRecord = (String, i64, i64, String);

/// (dataset_list0, dataset_list1) --> count, in insertion order
type MultiwayResults = Vec<((String, String), i64)>;

/// The result of comparing two variant tables
struct PairResult {
    results: PairwiseComparisonResults,
    bed_records: Vec<BedRecord>,
    longest_block_positions: Vec<i32>,
    longest_block_agreement: Vec<u8>,
}

/// An insertion-ordered dict from joint block id to variant indices
#[derive(Default)]
struct BlockIntersection {
    index: HashMap<Vec<BlockKey>, usize>,
    blocks: Vec<Vec<usize>>,
}

/// `compare`: prints the comparison of two or more variant tables; returns
/// the pairwise results for two tables, the multiway results otherwise.
fn compare(
    out: &mut impl Write,
    tables: &[&VariantTable],
    sample_names: &[&str],
    dataset_names: &[&str],
    ploidy: usize,
) -> Result<(Option<PairResult>, MultiwayResults)> {
    // common heterozygous variants
    let mut common: Option<HashSet<&VcfVariant>> = None;
    for (table, sample) in tables.iter().zip(sample_names) {
        let het: Vec<&VcfVariant> = table
            .variants
            .iter()
            .zip(table.genotypes_of(sample))
            .filter(|(_, gt)| !gt.is_homozygous())
            .map(|(v, _)| v)
            .collect();
        match &mut common {
            None => common = Some(het.into_iter().collect()),
            Some(c) => {
                let het: HashSet<&VcfVariant> = het.into_iter().collect();
                c.retain(|v| het.contains(v));
            }
        }
    }
    let common = common.expect("at least two tables");
    stat(out, "common heterozygous variants", common.len())?;
    print_stat(out, "(restricting to these below)", Stat::None, None)?;
    let mut sorted_variants: Vec<&VcfVariant> = common.iter().copied().collect();
    sorted_variants.sort_by_key(|v| v.position);
    let phases: Vec<Vec<&Option<VariantCallPhase>>> = tables
        .iter()
        .zip(sample_names)
        .map(|(table, sample)| {
            table
                .variants
                .iter()
                .zip(table.phases_of(sample))
                .filter(|(v, _)| common.contains(v))
                .map(|(_, p)| p)
                .collect()
        })
        .collect();

    let mut nan_counter = 0;
    let mut block_sizes: Vec<HashMap<BlockKey, usize>> = vec![HashMap::new(); tables.len()];
    let mut intersection = BlockIntersection::default();
    for variant_index in 0..common.len() {
        let mut any_none = false;
        let mut joint = Vec::with_capacity(phases.len());
        for (i, table_phases) in phases.iter().enumerate() {
            match table_phases[variant_index] {
                Some(phase) if phase.phase.iter().all(|a| a.is_some()) => {
                    let key = block_key(&phase.block_id, &mut nan_counter);
                    *block_sizes[i].entry(key.clone()).or_default() += 1;
                    joint.push(key);
                }
                _ => any_none = true,
            }
        }
        if !any_none {
            let n = intersection.blocks.len();
            let i = *intersection.index.entry(joint).or_insert(n);
            if i == n {
                intersection.blocks.push(Vec::new());
            }
            intersection.blocks[i].push(variant_index);
        }
    }
    for (name, sizes) in dataset_names.iter().zip(&block_sizes) {
        stat(out, &format!("non-singleton blocks in {name}"), sizes.values().filter(|&&n| n > 1).count())?;
        stat(out, "--> covered variants", sizes.values().filter(|&&n| n > 1).sum::<usize>())?;
    }
    let intersection_block_count = intersection.blocks.iter().filter(|b| b.len() > 1).count() as i64;
    let intersection_block_variants: i64 =
        intersection.blocks.iter().filter(|b| b.len() > 1).map(|b| b.len() as i64).sum();
    stat(out, "non-singleton intersection blocks", intersection_block_count)?;
    stat(out, "--> covered variants", intersection_block_variants)?;
    if tables.len() == 2 {
        let annotation = format!("{}<-->{}", dataset_names[0], dataset_names[1]);
        let pair = compare_pair(
            out,
            &intersection.blocks,
            intersection_block_count,
            intersection_block_variants,
            &phases,
            ploidy,
            &sorted_variants,
            &tables[0].chromosome,
            &annotation,
        )?;
        Ok((Some(pair), Vec::new()))
    } else {
        if ploidy != 2 {
            bail!("AssertionError: multiway comparison requires ploidy 2");
        }
        let multiway = compare_multiway(out, &intersection.blocks, dataset_names, &phases)?;
        Ok((None, multiway))
    }
}

/// The allele string of haplotype `j` of a block in one data set
fn block_haplotype(phases: &[&Option<VariantCallPhase>], block: &[usize], j: usize) -> Result<Hap> {
    let mut hap = Vec::new();
    for &i in block {
        let phase = phases[i].as_ref().expect("phases in intersection blocks are present");
        let allele = phase.phase.get(j).context("IndexError: tuple index out of range")?;
        hap.extend_from_slice(allele.expect("alleles in intersection blocks are present").to_string().as_bytes());
    }
    Ok(hap)
}

#[allow(clippy::too_many_arguments)]
fn compare_pair(
    out: &mut impl Write,
    blocks: &[Vec<usize>],
    intersection_block_count: i64,
    intersection_block_variants: i64,
    phases: &[Vec<&Option<VariantCallPhase>>],
    ploidy: usize,
    sorted_variants: &[&VcfVariant],
    chromosome: &str,
    annotation: &str,
) -> Result<PairResult> {
    let mut longest_block = 0;
    let mut longest_block_errors = PhasingErrors::default();
    let mut longest_block_positions = Vec::new();
    let mut longest_block_agreement = Vec::new();
    let mut phased_pairs = 0;
    let mut bed_records = Vec::new();
    let mut total_errors = PhasingErrors::default();
    let mut total_compared_variants = 0;
    for block in blocks.iter().filter(|b| b.len() >= 2) {
        let mut phasing0 = Vec::new();
        let mut phasing1 = Vec::new();
        for j in 0..ploidy {
            phasing0.push(block_haplotype(&phases[0], block, j)?);
            phasing1.push(block_haplotype(&phases[1], block, j)?);
        }
        let block_positions: Vec<i32> = block.iter().map(|&i| sorted_variants[i].position).collect();
        let errors = compare_block(&phasing0, &phasing1)?;
        if ploidy == 2 {
            // BedCreator.records: one record per switch error
            let (s0, s1) = (switch_encoding(&phasing0[0]), switch_encoding(&phasing1[0]));
            if phasing0[0].len() != phasing1[0].len() || phasing0[0].len() != block_positions.len() {
                bail!("AssertionError: phasings and positions of different lengths");
            }
            for (i, (sw0, sw1)) in s0.iter().zip(&s1).enumerate() {
                if sw0 != sw1 {
                    bed_records.push((
                        chromosome.to_string(),
                        i64::from(block_positions[i]) + 1,
                        i64::from(block_positions[i + 1]) + 1,
                        annotation.to_string(),
                    ));
                }
            }
        }
        total_errors += errors;
        phased_pairs += block.len() as i64 - 1;
        total_compared_variants += block.len() as i64;
        if block.len() > longest_block {
            longest_block = block.len();
            longest_block_errors = errors;
            longest_block_positions = block_positions;
            if ploidy == 2 {
                // compares the lists of haplotypes, not the haplotypes
                let agree = hamming(&phasing0, &phasing1)? < hamming(&phasing0[0], &complement(&phasing1[0])?)?;
                longest_block_agreement = phasing0[0]
                    .iter()
                    .zip(&phasing1[0])
                    .map(|(p0, p1)| u8::from(if agree { p0 == p1 } else { p0 != p1 }))
                    .collect();
            }
        }
    }
    let longest_block = longest_block as i64;
    let longest_block_assessed_pairs = (longest_block - 1).max(0);
    print_stat(out, "ALL INTERSECTION BLOCKS", Stat::Dash, None)?;
    print_errors(out, &total_errors, phased_pairs)?;
    stat(out, "Block-wise Hamming distance", total_errors.hamming)?;
    stat(out, "Block-wise Hamming distance [%]", fraction2percentstr(total_errors.hamming, total_compared_variants))?;
    stat(out, "Different genotypes", total_errors.diff_genotypes)?;
    stat(
        out,
        "Different genotypes [%]",
        fraction2percentstr(PyNum::Int(total_errors.diff_genotypes), total_compared_variants),
    )?;
    print_stat(out, "LARGEST INTERSECTION BLOCK", Stat::Dash, None)?;
    print_errors(out, &longest_block_errors, longest_block_assessed_pairs)?;
    stat(out, "Hamming distance", longest_block_errors.hamming)?;
    stat(out, "Hamming distance [%]", fraction2percentstr(longest_block_errors.hamming, longest_block))?;
    stat(out, "Different genotypes", longest_block_errors.diff_genotypes)?;
    stat(
        out,
        "Different genotypes [%]",
        fraction2percentstr(PyNum::Int(longest_block_errors.diff_genotypes), longest_block),
    )?;
    let sf_sum = |e: &PhasingErrors| e.switch_flips.switches + e.switch_flips.flips;
    let s = |v: PyNum| v.py_str();
    let values = vec![
        intersection_block_count.to_string(),
        intersection_block_variants.to_string(),
        phased_pairs.to_string(),
        s(total_errors.switches),
        s(safefraction(total_errors.switches, phased_pairs)),
        total_errors.switch_flips.to_string(),
        s(safefraction(sf_sum(&total_errors), phased_pairs)),
        s(total_errors.hamming),
        s(safefraction(total_errors.hamming, total_compared_variants)),
        total_errors.diff_genotypes.to_string(),
        s(safefraction(PyNum::Int(total_errors.diff_genotypes), total_compared_variants)),
        longest_block_assessed_pairs.to_string(),
        s(longest_block_errors.switches),
        s(safefraction(longest_block_errors.switches, longest_block_assessed_pairs)),
        longest_block_errors.switch_flips.to_string(),
        s(safefraction(sf_sum(&longest_block_errors), longest_block_assessed_pairs)),
        s(longest_block_errors.hamming),
        s(safefraction(longest_block_errors.hamming, longest_block)),
        longest_block_errors.diff_genotypes.to_string(),
        s(safefraction(PyNum::Int(longest_block_errors.diff_genotypes), longest_block)),
    ];
    Ok(PairResult {
        results: PairwiseComparisonResults { values },
        bed_records,
        longest_block_positions,
        longest_block_agreement,
    })
}

/// `compare_multiway`: how the data sets split at each pair of variants
fn compare_multiway(
    out: &mut impl Write,
    blocks: &[Vec<usize>],
    dataset_names: &[&str],
    phases: &[Vec<&Option<VariantCallPhase>>],
) -> Result<MultiwayResults> {
    let mut histogram: HashMap<Vec<u8>, i64> = HashMap::new();
    let mut total_compared = 0;
    for block in blocks.iter().filter(|b| b.len() >= 2) {
        total_compared += block.len() as i64 - 1;
        let switch_encodings: Vec<Vec<u8>> =
            phases.iter().map(|p| block_haplotype(p, block, 0).map(|h| switch_encoding(&h))).collect::<Result<_>>()?;
        for i in 0..block.len() - 1 {
            let s: Vec<u8> = switch_encodings.iter().map(|e| char_at(e, i)).collect::<Result<_>>()?;
            let c = complement(&s)?;
            *histogram.entry(s.min(c)).or_default() += 1;
        }
    }
    stat(out, "Compared pairs of variants", total_compared)?;
    let mut bipartitions: Vec<&Vec<u8>> = histogram.keys().collect();
    bipartitions.sort();
    let mut results: MultiwayResults = Vec::new();
    for (i, s) in bipartitions.into_iter().enumerate() {
        let count = histogram[s];
        if i == 0 {
            if s.iter().any(|&c| c != b'0') {
                bail!("AssertionError: no pairs of variants on which all phasings agree");
            }
            writeln!(out, "ALL AGREE")?;
        } else if i == 1 {
            writeln!(out, "DISAGREEMENT")?;
        }
        let (mut left, mut right) = (Vec::new(), Vec::new());
        for (name, &leftright) in dataset_names.iter().zip(s) {
            if leftright == b'0' {
                left.push(*name);
            } else {
                right.push(*name);
            }
        }
        let (left, right) = (left.join(","), right.join(","));
        print_stat(
            out,
            &format!("{{{left}}} vs. {{{right}}}"),
            Stat::Value(count.to_string()),
            Some(fraction2percentstr(PyNum::Int(count), total_compared)),
        )?;
        let key = (left, right);
        match results.iter_mut().find(|(k, _)| *k == key) {
            Some(entry) => entry.1 = count,
            None => results.push((key, count)),
        }
    }
    Ok(results)
}

/// `get_sample_names`
fn sample_names(readers: &[VcfReader], requested: Option<&str>, ignore_name: bool) -> Result<Vec<String>> {
    let mut first_samples = Vec::new();
    let mut intersection: Option<Vec<String>> = None;
    for reader in readers {
        intersection = Some(match intersection {
            None => reader.samples.clone(),
            Some(i) => i.into_iter().filter(|s| reader.samples.contains(s)).collect(),
        });
        if ignore_name && reader.samples.len() > 1 {
            bail!("File '{}' contains multiple samples, option --ignore-sample-name not available.", reader.path());
        }
        first_samples.push(reader.samples.first().context("IndexError: list index out of range")?.clone());
    }
    let mut intersection = intersection.expect("at least two readers");
    intersection.sort();
    intersection.dedup();
    if let Some(requested) = requested.filter(|s| !s.is_empty()) {
        if !intersection.iter().any(|s| s == requested) {
            bail!("Sample '{requested}' requested on command-line not found in all VCFs");
        }
        Ok(vec![requested.to_string(); readers.len()])
    } else if ignore_name {
        Ok(first_samples)
    } else {
        match intersection.len() {
            0 => bail!("None of the samples is present in all VCFs"),
            1 => Ok(vec![intersection[0].clone(); readers.len()]),
            _ => bail!(
                "More than one sample is present in all VCFs, please use --sample to specify which sample to work on."
            ),
        }
    }
}

fn create(path: &str) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path).with_context(|| format!("cannot create {path}"))?))
}

/// `run_compare`, printing the report to `out`
pub fn run_compare(opts: &CompareOptions, out: &mut impl Write) -> Result<()> {
    if opts.vcf.len() < 2 {
        bail!("At least two VCFs need to be given.");
    }
    if opts.ploidy < 2 {
        bail!("Ploidy must be > 1.");
    }
    if opts.ploidy > 2 && opts.tsv_multiway.is_some() {
        bail!("Option --tsv-multiway can only be used if ploidy=2.");
    }
    if opts.ploidy > 2 && opts.switch_error_bed.is_some() {
        bail!("Option --switch-error-bed can only be used if ploidy=2.");
    }
    if opts.ploidy > 2 && opts.longest_block_tsv.is_some() {
        bail!("Option --longest-block-tsv can only be used if ploidy=2.");
    }
    if opts.plot_blocksizes.is_some() || opts.plot_sum_of_blocksizes.is_some() {
        bail!("Options --plot-blocksizes and --plot-sum-of-blocksizes are not supported by lrpipe-phase");
    }
    let ploidy = opts.ploidy;
    let mut readers = Vec::new();
    for f in &opts.vcf {
        readers.push(VcfReader::open_with(f, opts.only_snvs, true, Some(ploidy))?.with_mav(ploidy > 2));
    }
    let dataset_names: Vec<String> = match opts.names.as_deref().filter(|n| !n.is_empty()) {
        Some(names) => {
            let names: Vec<String> = names.split(',').map(String::from).collect();
            if names.len() != opts.vcf.len() {
                bail!("Number of names given with --names does not equal number of VCFs.");
            }
            names
        }
        None => (0..opts.vcf.len()).map(|i| format!("file{i}")).collect(),
    };
    let sample_names = sample_names(&readers, opts.sample.as_deref(), opts.ignore_sample_name)?;

    let mut tsv_pairwise = opts.tsv_pairwise.as_deref().map(create).transpose()?;
    let mut tsv_multiway = match &opts.tsv_multiway {
        Some(path) => {
            let mut f = create(path)?;
            writeln!(f, "#sample\tchromosome\tdataset_list0\tdataset_list1\tcount")?;
            Some(f)
        }
        None => None,
    };
    let mut longest_block_tsv = match &opts.longest_block_tsv {
        Some(path) => {
            let mut f = create(path)?;
            writeln!(f, "#dataset_name0\tdataset_name1\t#sample\tchromosome\tposition\tphase_agreeing")?;
            Some(f)
        }
        None => None,
    };
    if let Some(f) = &mut tsv_pairwise {
        writeln!(
            f,
            "#sample\tchromosome\tdataset_name0\tdataset_name1\tfile_name0\tfile_name1\t{}\thet_variants0\tonly_snvs",
            PAIRWISE_FIELDS.join("\t")
        )?;
    }
    let mut switch_error_bed = opts.switch_error_bed.as_deref().map(create).transpose()?;

    let unique_samples: Vec<&String> = {
        let mut seen = Vec::new();
        for s in &sample_names {
            if !seen.contains(&s) {
                seen.push(s);
            }
        }
        seen
    };
    if unique_samples.len() > 1 && opts.ignore_sample_name {
        writeln!(out, "Comparing phasings for samples: {}  (--ignore-sample-names selected)", sample_names.join(", "))?;
    } else {
        writeln!(out, "Comparing phasings for sample {}", sample_names[0])?;
    }

    // get_variant_tables
    let mut vcfs: Vec<HashMap<String, VariantTable>> = Vec::new();
    for (reader, filename) in readers.iter_mut().zip(&opts.vcf) {
        log::info!("Reading phasing from '{filename}'");
        let mut m = HashMap::new();
        loop {
            match reader.next_table() {
                Ok(Some(table)) => {
                    m.insert(table.chromosome.clone(), table);
                }
                Ok(None) => break,
                Err(e) => match e.downcast_ref::<PloidyError>() {
                    Some(p) => bail!("Provided ploidy is invalid: {p}. Aborting."),
                    None => return Err(e),
                },
            }
        }
        vcfs.push(m);
    }
    let mut chromosomes: Vec<&String> =
        vcfs[0].keys().filter(|c| vcfs[1..].iter().all(|m| m.contains_key(*c))).collect();
    chromosomes.sort();
    if chromosomes.is_empty() {
        bail!("No chromosome is contained in all VCFs. Aborting.");
    }
    log::info!(
        "Chromosomes present in all VCFs: {}",
        chromosomes.iter().map(|c| c.as_str()).collect::<Vec<_>>().join(", ")
    );

    writeln!(out, "FILENAMES")?;
    let longest_name = dataset_names.iter().map(|n| n.chars().count()).max().expect("at least two names");
    for (name, filename) in dataset_names.iter().zip(&opts.vcf) {
        writeln!(out, "{} = {filename}", rjust(name, longest_name + 2))?;
    }
    let width = longest_name.max(15) + 5;
    let names: Vec<&str> = dataset_names.iter().map(String::as_str).collect();
    let samples: Vec<&str> = sample_names.iter().map(String::as_str).collect();

    for chromosome in chromosomes {
        writeln!(out, "---------------- Chromosome {chromosome} ----------------")?;
        let mut all_bed_records = Vec::new();
        let tables: Vec<&VariantTable> = vcfs.iter().map(|m| &m[chromosome]).collect();
        let mut all_union: HashSet<&VcfVariant> = HashSet::new();
        let mut all_intersection: Option<HashSet<&VcfVariant>> = None;
        let mut het_union: HashSet<&VcfVariant> = HashSet::new();
        let mut het_intersection: Option<HashSet<&VcfVariant>> = None;
        let mut het_variants0 = None;
        writeln!(out, "VARIANT COUNTS (heterozygous / all): ")?;
        for ((table, name), sample) in tables.iter().zip(&names).zip(&samples) {
            all_union.extend(table.variants.iter());
            let het: Vec<&VcfVariant> = table
                .variants
                .iter()
                .zip(table.genotypes_of(sample))
                .filter(|(_, gt)| !gt.is_homozygous())
                .map(|(v, _)| v)
                .collect();
            het_variants0.get_or_insert(het.len());
            het_union.extend(het.iter().copied());
            let all: HashSet<&VcfVariant> = table.variants.iter().collect();
            let het_set: HashSet<&VcfVariant> = het.iter().copied().collect();
            match (&mut all_intersection, &mut het_intersection) {
                (Some(a), Some(h)) => {
                    a.retain(|v| all.contains(v));
                    h.retain(|v| het_set.contains(v));
                }
                _ => {
                    all_intersection = Some(all);
                    het_intersection = Some(het_set);
                }
            }
            writeln!(
                out,
                "{} {} / {}",
                rjust(&format!("{name}:"), width),
                rjust(&het.len().to_string(), COUNT_WIDTH),
                rjust(&table.variants.len().to_string(), COUNT_WIDTH)
            )?;
        }
        let count = |s: &HashSet<&VcfVariant>| rjust(&s.len().to_string(), COUNT_WIDTH);
        writeln!(out, "{} {} / {}", rjust("UNION:", width), count(&het_union), count(&all_union))?;
        writeln!(
            out,
            "{} {} / {}",
            rjust("INTERSECTION:", width),
            count(het_intersection.as_ref().expect("at least two tables")),
            count(all_intersection.as_ref().expect("at least two tables"))
        )?;

        for i in 0..tables.len() {
            for j in i + 1..tables.len() {
                writeln!(out, "PAIRWISE COMPARISON: {} <--> {}:", names[i], names[j])?;
                let (pair, _) =
                    compare(out, &[tables[i], tables[j]], &[samples[i], samples[j]], &[names[i], names[j]], ploidy)?;
                let pair = pair.expect("two tables give pairwise results");
                all_bed_records.extend(pair.bed_records);
                let sample_name = if opts.ignore_sample_name {
                    format!("{}_{}", samples[i], samples[j])
                } else {
                    samples[i].to_string()
                };
                if let Some(f) = &mut tsv_pairwise {
                    writeln!(
                        f,
                        "{sample_name}\t{chromosome}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        names[i],
                        names[j],
                        opts.vcf[i],
                        opts.vcf[j],
                        pair.results.values.join("\t"),
                        het_variants0.expect("at least one table"),
                        u8::from(opts.only_snvs)
                    )?;
                }
                if let Some(f) = &mut longest_block_tsv {
                    for (position, agreeing) in pair.longest_block_positions.iter().zip(&pair.longest_block_agreement) {
                        writeln!(f, "{}\t{}\t{sample_name}\t{chromosome}\t{position}\t{agreeing}", names[i], names[j])?;
                    }
                }
            }
        }

        // all switch errors found in the current chromosome
        if let Some(f) = &mut switch_error_bed {
            all_bed_records.sort();
            for (chrom, start, end, annotation) in &all_bed_records {
                writeln!(f, "{chrom}\t{start}\t{end}\t{annotation}")?;
            }
        }

        if tables.len() > 2 {
            if ploidy != 2 {
                bail!("AssertionError: multiway comparison requires ploidy 2");
            }
            writeln!(out, "MULTIWAY COMPARISON OF ALL PHASINGS:")?;
            let (_, multiway) = compare(out, &tables, &samples, &names, ploidy)?;
            if let Some(f) = &mut tsv_multiway {
                // WhatsHap joins a set of the sample names, whose order varies
                // between runs; here they are in order of first appearance
                let sample_name = if opts.ignore_sample_name {
                    unique_samples.iter().map(|s| s.as_str()).collect::<Vec<_>>().join("_")
                } else {
                    samples[0].to_string()
                };
                for ((list0, list1), count) in multiway {
                    writeln!(f, "{sample_name}\t{chromosome}\t{{{list0}}}\t{{{list1}}}\t{count}")?;
                }
            }
        }
    }
    for mut f in [tsv_pairwise, tsv_multiway, longest_block_tsv, switch_error_bed].into_iter().flatten() {
        f.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn switch_flip_decomposition() {
        let sf = compute_switch_flips(b"00011", b"00100").unwrap();
        assert_eq!((sf.switches, sf.flips), (PyNum::Int(1), PyNum::Int(0)));
        let sf = compute_switch_flips(b"00011", b"00111").unwrap();
        assert_eq!((sf.switches, sf.flips), (PyNum::Int(0), PyNum::Int(1)));
        let sf = compute_switch_flips(b"000", b"001").unwrap();
        assert_eq!((sf.switches, sf.flips), (PyNum::Int(1), PyNum::Int(0)));
        assert_eq!(switch_encoding(b"0001011"), b"001110");
        assert_eq!(complement(b"01100").unwrap(), b"10011");
    }

    #[test]
    fn itertools_permutations() {
        assert_eq!(
            permutations(3),
            vec![vec![0, 1, 2], vec![0, 2, 1], vec![1, 0, 2], vec![1, 2, 0], vec![2, 0, 1], vec![2, 1, 0]]
        );
    }
}

//! Port of `cli/stats.py`: phasing statistics of a single VCF file.

use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use crate::compat::pyfloat::PyNum as Num;
use crate::compat::pystr::py_int;
use crate::hts::{IndexedVcf, VcfFetchError};
use crate::vcf::{BlockId, BlockKey, VariantTable, VcfReader, block_key};

pub struct StatsOptions {
    pub vcf: String,
    pub sample: Option<String>,
    pub gtf: Option<String>,
    pub tsv: Option<String>,
    pub block_list: Option<String>,
    pub only_snvs: bool,
    pub chromosomes: Vec<String>,
    pub chr_lengths: Option<String>,
}

/// A phase block: positions of its variants, in increasing order (variants in
/// a variant table have distinct, sorted positions).
#[derive(Clone, Debug, Default)]
struct PhasedBlock {
    positions: Vec<i32>,
    snvs: usize,
}

impl PhasedBlock {
    fn len(&self) -> usize {
        self.positions.len()
    }

    fn leftmost(&self) -> i32 {
        self.positions[0]
    }

    fn rightmost(&self) -> i32 {
        self.positions[self.positions.len() - 1]
    }

    fn span(&self) -> i64 {
        i64::from(self.rightmost()) - i64::from(self.leftmost())
    }

    /// `PhasedBlock.split`: the variants left of `split_left` and those right
    /// of `split_right` (the ones in between are dropped).
    fn split(&self, split_left: i32, split_right: i32) -> (PhasedBlock, PhasedBlock) {
        let left = self.positions.iter().copied().filter(|&p| p < split_left).collect();
        let right = self.positions.iter().copied().filter(|&p| p > split_right).collect();
        (PhasedBlock { positions: left, snvs: 0 }, PhasedBlock { positions: right, snvs: 0 })
    }
}

/// A block after splitting interleaved blocks: only its span and chromosome
/// are used.
struct SplitBlock {
    chromosome: Arc<str>,
    span: i64,
}

/// A block as needed for the block-size statistics.
struct BlockSummary {
    len: usize,
    snvs: usize,
}

const NAN: Num = Num::Float(f64::NAN);

const FIELD_NAMES: [&str; 21] = [
    "variants",
    "phased",
    "unphased",
    "singletons",
    "blocks",
    "variant_per_block_median",
    "variant_per_block_avg",
    "variant_per_block_min",
    "variant_per_block_max",
    "variant_per_block_sum",
    "bp_per_block_median",
    "bp_per_block_avg",
    "bp_per_block_min",
    "bp_per_block_max",
    "bp_per_block_sum",
    "heterozygous_variants",
    "heterozygous_snvs",
    "phased_snvs",
    "phased_fraction",
    "phased_snvs_fraction",
    "block_n50",
];

struct DetailedStats {
    variants: Num,
    phased: Num,
    unphased: Num,
    singletons: Num,
    blocks: Num,
    variant_per_block_median: Num,
    variant_per_block_avg: Num,
    variant_per_block_min: Num,
    variant_per_block_max: Num,
    variant_per_block_sum: Num,
    bp_per_block_median: Num,
    bp_per_block_avg: Num,
    bp_per_block_min: Num,
    bp_per_block_max: Num,
    bp_per_block_sum: Num,
    heterozygous_variants: Num,
    heterozygous_snvs: Num,
    phased_snvs: Num,
    phased_fraction: Num,
    phased_snvs_fraction: Num,
    block_n50: Num,
}

impl DetailedStats {
    /// `dataclasses.astuple`, in the order of `FIELD_NAMES`
    fn values(&self) -> [Num; 21] {
        [
            self.variants,
            self.phased,
            self.unphased,
            self.singletons,
            self.blocks,
            self.variant_per_block_median,
            self.variant_per_block_avg,
            self.variant_per_block_min,
            self.variant_per_block_max,
            self.variant_per_block_sum,
            self.bp_per_block_median,
            self.bp_per_block_avg,
            self.bp_per_block_min,
            self.bp_per_block_max,
            self.bp_per_block_sum,
            self.heterozygous_variants,
            self.heterozygous_snvs,
            self.phased_snvs,
            self.phased_fraction,
            self.phased_snvs_fraction,
            self.block_n50,
        ]
    }

    /// `DetailedStats.print`
    fn print(&self, out: &mut impl Write) -> Result<()> {
        let max_integer_width = self
            .values()
            .iter()
            .filter(|v| !v.is_nan())
            .map(|v| match v {
                Num::Int(i) => i.to_string().len(),
                Num::Float(f) => (f.trunc() as i64).to_string().len(),
            })
            .max()
            .expect("there are values");
        let w = max_integer_width.max(8);
        let fi = |v: Num| fmt_fixed(v.as_f64(), 0, w);
        let ff = |v: Num| fmt_fixed(v.as_f64(), 2, w + 3);
        let fp = |v: Num| {
            let x = v.as_f64() * 100.0;
            let s = if x.is_nan() { "nan%".to_string() } else { format!("{x:.1}%") };
            format!("{s:>width$}", width = w + 3)
        };
        let lines = [
            format!("{:>21}: {}", "Variants in VCF", fi(self.variants)),
            format!(
                "{:>21}: {}    ({}    SNVs)",
                "Heterozygous",
                fi(self.heterozygous_variants),
                fi(self.heterozygous_snvs)
            ),
            format!("{:>21}: {}    ({}    SNVs)", "Phased", fi(self.phased), fi(self.phased_snvs)),
            format!(
                "{:>21}: {} ({} SNVs)",
                "Heterozygous phased",
                fp(self.phased_fraction),
                fp(self.phased_snvs_fraction)
            ),
            format!("{:>21}: {}    (not considered below)", "Unphased", fi(self.unphased)),
            format!("{:>21}: {}    (not considered below)", "Singletons", fi(self.singletons)),
            format!("{:>21}: {}", "Blocks", fi(self.blocks)),
            String::new(),
            "Block sizes (no. of variants)".to_string(),
            format!("{:>21}: {}    variants", "Sum of sizes", fi(self.variant_per_block_sum)),
            format!("{:>21}: {} variants", "Median block size", ff(self.variant_per_block_median)),
            format!("{:>21}: {} variants", "Average block size", ff(self.variant_per_block_avg)),
            format!("{:>21}: {}    variants", "Largest block", fi(self.variant_per_block_max)),
            format!("{:>21}: {}    variants", "Smallest block", fi(self.variant_per_block_min)),
            String::new(),
            "Block lengths (basepairs)".to_string(),
            format!("{:>21}: {}    bp", "Sum of lengths", fi(self.bp_per_block_sum)),
            format!("{:>21}: {} bp", "Median block length", ff(self.bp_per_block_median)),
            format!("{:>21}: {} bp", "Average block length", ff(self.bp_per_block_avg)),
            format!("{:>21}: {}    bp", "Longest block", fi(self.bp_per_block_max)),
            format!("{:>21}: {}    bp", "Shortest block", fi(self.bp_per_block_min)),
            format!("{:>21}: {}    bp", "Block NG50", fi(self.block_n50)),
        ];
        for line in lines {
            writeln!(out, "{line}")?;
        }
        Ok(())
    }
}

/// Python's `format(x, f"{width}.{precision}f")`
fn fmt_fixed(x: f64, precision: usize, width: usize) -> String {
    let s = if x.is_nan() { "nan".to_string() } else { format!("{x:.precision$}") };
    format!("{s:>width$}")
}

/// `statistics.median` of a sorted, non-empty list of ints
fn median(sorted: &[i64]) -> Num {
    let n = sorted.len();
    if n % 2 == 1 { Num::Int(sorted[n / 2]) } else { Num::Float((sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0) }
}

/// `n50`: the length at which the lengths, longest first, add up to half of
/// the target length (0 if they never do)
fn n50(lengths: &mut [i64], target_length: i128) -> i64 {
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let mut total: i128 = 0;
    for &length in lengths.iter() {
        total += i128::from(length);
        // total >= 0.5 * target_length
        if 2 * total >= target_length {
            return length;
        }
    }
    0
}

/// `compute_ng50`
fn compute_ng50(blocks: &[SplitBlock], chr_lengths: &HashMap<String, i128>) -> Num {
    let mut chromosomes: Vec<&str> = blocks.iter().map(|b| &*b.chromosome).collect();
    chromosomes.sort_unstable();
    chromosomes.dedup();
    let mut target_length = 0;
    for chromosome in chromosomes {
        match chr_lengths.get(chromosome) {
            Some(length) => target_length += length,
            None => {
                log::warn!("Not able to compute NG50 because length of contig '{chromosome}' not available");
                return NAN;
            }
        }
    }
    let mut block_lengths: Vec<i64> = blocks.iter().map(|b| b.span).collect();
    Num::Int(n50(&mut block_lengths, target_length))
}

#[derive(Default)]
struct PhasingStats {
    blocks: Vec<BlockSummary>,
    split_blocks: Vec<SplitBlock>,
    unphased: i64,
    variants: i64,
    heterozygous_variants: i64,
    heterozygous_snvs: i64,
}

impl PhasingStats {
    fn extend(&mut self, other: PhasingStats) {
        self.blocks.extend(other.blocks);
        self.split_blocks.extend(other.split_blocks);
        self.unphased += other.unphased;
        self.variants += other.variants;
        self.heterozygous_variants += other.heterozygous_variants;
        self.heterozygous_snvs += other.heterozygous_snvs;
    }

    /// `add_blocks` for the blocks of one chromosome
    fn add_blocks(&mut self, chromosome: &Arc<str>, blocks: &[PhasedBlock]) {
        self.blocks.extend(blocks.iter().map(|b| BlockSummary { len: b.len(), snvs: b.snvs }));
        self.split_blocks.extend(
            nonoverlapping_blocks(blocks)
                .into_iter()
                .map(|b| SplitBlock { chromosome: Arc::clone(chromosome), span: b.span() }),
        );
    }

    /// `get_detailed_stats`
    fn detailed_stats(&self, chr_lengths: &HashMap<String, i128>) -> DetailedStats {
        let mut block_sizes: Vec<i64> = self.blocks.iter().filter(|b| b.len > 1).map(|b| b.len as i64).collect();
        block_sizes.sort_unstable();
        let n_singletons = self.blocks.iter().filter(|b| b.len == 1).count() as i64;
        // Block length stats calculated from split interleaved blocks to avoid inflating values
        let mut block_lengths: Vec<i64> = self.split_blocks.iter().map(|b| b.span).collect();
        block_lengths.sort_unstable();
        let phased_snvs: i64 = self.blocks.iter().filter(|b| b.len > 1).map(|b| b.snvs as i64).sum();
        if block_sizes.is_empty() {
            return DetailedStats {
                variants: Num::Int(self.variants),
                phased: Num::Int(0),
                unphased: Num::Int(self.unphased),
                singletons: Num::Int(n_singletons),
                blocks: Num::Int(0),
                variant_per_block_median: NAN,
                variant_per_block_avg: NAN,
                variant_per_block_min: Num::Int(0),
                variant_per_block_max: Num::Int(0),
                variant_per_block_sum: Num::Int(0),
                bp_per_block_median: NAN,
                bp_per_block_avg: NAN,
                bp_per_block_min: Num::Int(0),
                bp_per_block_max: Num::Int(0),
                bp_per_block_sum: Num::Int(0),
                heterozygous_variants: Num::Int(self.heterozygous_variants),
                heterozygous_snvs: Num::Int(self.heterozygous_snvs),
                phased_snvs: Num::Int(0),
                phased_fraction: Num::Float(0.0),
                phased_snvs_fraction: Num::Float(0.0),
                block_n50: NAN,
            };
        }
        // a block of two or more variants leaves at least one split block
        let sizes_sum: i64 = block_sizes.iter().sum();
        let lengths_sum: i64 = block_lengths.iter().sum();
        DetailedStats {
            variants: Num::Int(self.variants),
            phased: Num::Int(sizes_sum),
            unphased: Num::Int(self.unphased),
            singletons: Num::Int(n_singletons),
            blocks: Num::Int(block_sizes.len() as i64),
            variant_per_block_median: median(&block_sizes),
            variant_per_block_avg: Num::Float(sizes_sum as f64 / block_sizes.len() as f64),
            variant_per_block_min: Num::Int(block_sizes[0]),
            variant_per_block_max: Num::Int(block_sizes[block_sizes.len() - 1]),
            variant_per_block_sum: Num::Int(sizes_sum),
            bp_per_block_median: median(&block_lengths),
            bp_per_block_avg: Num::Float(lengths_sum as f64 / block_lengths.len() as f64),
            bp_per_block_min: Num::Int(block_lengths[0]),
            bp_per_block_max: Num::Int(block_lengths[block_lengths.len() - 1]),
            bp_per_block_sum: Num::Int(lengths_sum),
            heterozygous_variants: Num::Int(self.heterozygous_variants),
            heterozygous_snvs: Num::Int(self.heterozygous_snvs),
            phased_snvs: Num::Int(phased_snvs),
            phased_fraction: if self.heterozygous_variants != 0 {
                Num::Float(sizes_sum as f64 / self.heterozygous_variants as f64)
            } else {
                NAN
            },
            phased_snvs_fraction: if self.heterozygous_snvs != 0 {
                Num::Float(phased_snvs as f64 / self.heterozygous_snvs as f64)
            } else {
                NAN
            },
            block_n50: compute_ng50(&self.split_blocks, chr_lengths),
        }
    }
}

/// A block in the queue of `get_nonoverlapping_blocks`, ordered so that the
/// block with the leftmost start is popped first.
struct Queued(PhasedBlock);

impl PartialEq for Queued {
    fn eq(&self, other: &Self) -> bool {
        self.0.leftmost() == other.0.leftmost()
    }
}

impl Eq for Queued {}

impl PartialOrd for Queued {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Queued {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.0.leftmost().cmp(&self.0.leftmost())
    }
}

/// `get_nonoverlapping_blocks` for the blocks of one chromosome: split
/// blocks where another block starts before they end. The blocks' variants
/// are disjoint, so no two blocks start at the same position.
fn nonoverlapping_blocks(blocks: &[PhasedBlock]) -> Vec<PhasedBlock> {
    let mut queue: BinaryHeap<Queued> = blocks.iter().filter(|b| b.len() > 1).cloned().map(Queued).collect();
    let mut split_blocks = Vec::new();
    while let Some(Queued(mut block)) = queue.pop() {
        if let Some(Queued(next_block)) = queue.peek() {
            // Check if next block overlaps current. If so split the current block.
            if block.rightmost() > next_block.leftmost() {
                let (left, right) = block.split(next_block.leftmost(), next_block.rightmost());
                if right.len() > 1 {
                    queue.push(Queued(right));
                }
                // Skip the left-side block if is it too short after splitting
                if left.len() < 2 {
                    continue;
                }
                block = left;
            }
        }
        split_blocks.push(block);
    }
    split_blocks
}

/// Python's `a != b` for two block ids
fn block_ids_differ(a: &BlockId, b: &BlockId) -> bool {
    match (a, b) {
        (BlockId::Int(x), BlockId::Int(y)) => x != y,
        (BlockId::Int(x), BlockId::Float(y)) | (BlockId::Float(y), BlockId::Int(x)) => *x as f64 != *y,
        (BlockId::Float(x), BlockId::Float(y)) => x != y,
        (BlockId::Str(x), BlockId::Str(y)) => x != y,
        _ => true,
    }
}

/// `str(block_id)`
fn block_id_str(id: &Option<BlockId>) -> String {
    match id {
        None => "None".to_string(),
        Some(id) => id.to_string(),
    }
}

fn py_type_name(id: &Option<BlockId>) -> &'static str {
    match id {
        None => "NoneType",
        Some(BlockId::Int(_)) => "int",
        Some(BlockId::Float(_)) => "float",
        Some(BlockId::Str(_)) => "str",
    }
}

/// The phase blocks of a chromosome by block id, in the order in which the
/// ids first occur (a Python dict). Keys are the first-seen id objects.
struct Blocks {
    keys: Vec<Option<BlockId>>,
    blocks: Vec<PhasedBlock>,
}

impl Blocks {
    /// `sorted(blocks.keys())`, failing like Python where ids cannot be
    /// compared
    fn sorted_indices(&self) -> Result<Vec<usize>> {
        let mut indices: Vec<usize> = (0..self.keys.len()).collect();
        if indices.len() < 2 {
            return Ok(indices);
        }
        let all_str = self.keys.iter().all(|k| matches!(k, Some(BlockId::Str(_))));
        let all_num = self.keys.iter().all(|k| matches!(k, Some(BlockId::Int(_) | BlockId::Float(_))));
        if !all_str && !all_num {
            // the first pair of incomparable types
            let first = &self.keys[0];
            let other =
                self.keys.iter().find(|k| py_type_name(k) != py_type_name(first) || k.is_none()).unwrap_or(first);
            bail!("'<' not supported between instances of '{}' and '{}'", py_type_name(other), py_type_name(first));
        }
        let value = |k: &Option<BlockId>| match k {
            Some(BlockId::Int(v)) => *v as f64,
            Some(BlockId::Float(v)) => *v,
            _ => 0.0,
        };
        if all_str {
            indices.sort_by(|&a, &b| block_id_str(&self.keys[a]).cmp(&block_id_str(&self.keys[b])));
        } else {
            indices.sort_by(|&a, &b| {
                value(&self.keys[a]).partial_cmp(&value(&self.keys[b])).unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        Ok(indices)
    }
}

struct GtfWriter {
    file: BufWriter<File>,
}

impl GtfWriter {
    /// Write a feature to the GTF. start is 0-based.
    fn write(&mut self, chromosome: &str, start: i64, stop: i64, name: &str) -> Result<()> {
        assert!(start < stop);
        writeln!(
            self.file,
            "{chromosome}\tPhasing\texon\t{}\t{stop}\t.\t+\t.\tgene_id \"{name}\"; transcript_id \"{name}.1\";",
            start + 1
        )?;
        Ok(())
    }
}

/// A run of consecutive phased variants of the same block, for the GTF
struct GtfBlock {
    start: i64,
    end: i64,
    id: Option<BlockId>,
}

/// `get_phase_blocks`: count variants into `stats` and return the phase
/// blocks of `sample`.
fn phase_blocks(
    chromosome: &str,
    gtfwriter: &mut Option<GtfWriter>,
    sample: &str,
    stats: &mut PhasingStats,
    table: &VariantTable,
) -> Result<Blocks> {
    let genotypes = table.genotypes_of(sample);
    let phases = table.phases_of(sample);
    let mut index: HashMap<BlockKey, usize> = HashMap::new();
    let mut nan_counter = 0;
    let mut blocks = Blocks { keys: Vec::new(), blocks: Vec::new() };
    let mut prev_block = GtfBlock { start: 0, end: 0, id: None };
    for ((variant, genotype), phase) in table.variants.iter().zip(genotypes).zip(phases) {
        stats.variants += 1;
        if genotype.is_homozygous() {
            continue;
        }
        stats.heterozygous_variants += 1;
        let is_snv = variant.is_snv();
        if is_snv {
            stats.heterozygous_snvs += 1;
        }
        let Some(phase) = phase else {
            stats.unphased += 1;
            continue;
        };
        let i = *index.entry(block_key(&phase.block_id, &mut nan_counter)).or_insert_with(|| {
            blocks.keys.push(phase.block_id.clone());
            blocks.blocks.push(PhasedBlock::default());
            blocks.blocks.len() - 1
        });
        let block = &mut blocks.blocks[i];
        block.positions.push(variant.position);
        block.snvs += usize::from(is_snv);
        if let Some(gtf) = gtfwriter {
            let pos = i64::from(variant.position);
            let new_block = GtfBlock { start: pos, end: pos + 1, id: phase.block_id.clone() };
            match &prev_block.id {
                None => prev_block = new_block,
                Some(prev_id) => {
                    if phase.block_id.as_ref().is_none_or(|id| block_ids_differ(prev_id, id)) {
                        gtf.write(chromosome, prev_block.start, prev_block.end, &prev_id.to_string())?;
                        prev_block = new_block;
                    }
                    prev_block.end = pos + 1;
                }
            }
        }
    }
    if let Some(gtf) = gtfwriter
        && let Some(id) = &prev_block.id
    {
        gtf.write(chromosome, prev_block.start, prev_block.end, &id.to_string())?;
    }
    Ok(blocks)
}

/// `unpack_chromosomes`: split comma-separated entries
fn unpack_chromosomes(chromosomes: &[String]) -> Vec<String> {
    chromosomes.iter().flat_map(|entry| entry.split(',')).filter(|c| !c.is_empty()).map(String::from).collect()
}

/// `parse_chr_lengths`: two tab-separated columns, name and length
fn parse_chr_lengths(filename: &str) -> Result<HashMap<String, i128>> {
    let content = std::fs::read(filename).with_context(|| format!("cannot open {filename}"))?;
    let content = String::from_utf8(content).with_context(|| format!("{filename} is not UTF-8"))?;
    let mut chr_lengths = HashMap::new();
    for line in content.split_inclusive('\n') {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            bail!("AssertionError: line in {filename} does not have two tab-separated fields: {line:?}");
        }
        let length =
            py_int(fields[1]).with_context(|| format!("invalid literal for int() with base 10: {:?}", fields[1]))?;
        chr_lengths.insert(fields[0].to_string(), length);
    }
    Ok(chr_lengths)
}

/// `get_chr_lengths`
fn chr_lengths(reader: &VcfReader, chr_lengths_file: Option<&str>) -> Result<HashMap<String, i128>> {
    if let Some(file) = chr_lengths_file {
        let chr_lengths = parse_chr_lengths(file)?;
        log::info!("Read length of {} chromosomes from {}", chr_lengths.len(), file);
        Ok(chr_lengths)
    } else {
        let chr_lengths: HashMap<String, i128> = reader
            .header()
            .contig_lengths()
            .into_iter()
            .filter(|&(_, length)| length != 0)
            .map(|(name, length)| (name, i128::from(length)))
            .collect();
        if chr_lengths.is_empty() {
            log::warn!("VCF header does not contain contig lengths, cannot compute NG50. Consider using --chr-lengths");
        }
        Ok(chr_lengths)
    }
}

/// Variant tables of the requested chromosomes by index lookup, or of all
/// chromosomes in file order (`parse_variant_tables`).
struct VariantTables {
    reader: VcfReader,
    fetch: Option<(IndexedVcf, std::vec::IntoIter<String>)>,
}

impl VariantTables {
    fn next(&mut self) -> Result<Option<VariantTable>> {
        let Some((indexed, chromosomes)) = &mut self.fetch else { return self.reader.next_table() };
        let Some(chromosome) = chromosomes.next() else { return Ok(None) };
        match self.reader.fetch_regions(indexed, &chromosome, &[(0, None)]) {
            Ok(table) => Ok(Some(table)),
            Err(VcfFetchError::InvalidContig) => bail!("invalid contig `{chromosome}`"),
            Err(VcfFetchError::IndexMissing) => bail!("{} is missing an index (.tbi or .csi)", self.reader.path()),
            Err(VcfFetchError::Other(e)) => Err(e),
        }
    }
}

fn create(path: &str) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path).with_context(|| format!("cannot create {path}"))?))
}

/// The phasing of a VCF's first sample over all its chromosomes, as `stats` reports it for
/// them all (its `ALL` line, or its one chromosome's).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Totals {
    /// Heterozygous SNVs (`heterozygous_snvs`).
    pub heterozygous_snvs: i64,
    /// Those in a block of two or more variants (`phased_snvs`).
    pub phased_snvs: i64,
    /// Blocks of two or more variants (`blocks`).
    pub blocks: i64,
    /// The blocks' NG50 (`block_n50`): the span of the shortest of the longest blocks, split
    /// where blocks interleave, that together span at least half the chromosomes with blocks
    /// (by the header's contig lengths); none without blocks or lengths.
    pub block_n50: Option<i64>,
}

/// [`Totals`] of the VCF at `path`, multi-allelic sites left out as `stats` leaves them.
///
/// # Errors
///
/// If the VCF cannot be read, or has no samples.
pub fn totals(path: &str) -> Result<Totals> {
    let reader = VcfReader::open_with(path, false, true, None)?;
    let Some(sample) = reader.samples.first().cloned() else { bail!("{path} has no samples") };
    let chr_lengths = chr_lengths(&reader, None)?;
    let mut tables = VariantTables { reader, fetch: None };
    let mut total = PhasingStats::default();
    while let Some(table) = tables.next()? {
        let chromosome: Arc<str> = Arc::from(table.chromosome.as_str());
        let mut stats = PhasingStats::default();
        let blocks = phase_blocks(&chromosome, &mut None, &sample, &mut stats, &table)?;
        stats.add_blocks(&chromosome, &blocks.blocks);
        total.extend(stats);
    }
    let detailed = total.detailed_stats(&chr_lengths);
    let int = |n: Num| match n {
        Num::Int(n) => Some(n),
        Num::Float(_) => None,
    };
    Ok(Totals {
        heterozygous_snvs: total.heterozygous_snvs,
        phased_snvs: int(detailed.phased_snvs).unwrap_or(0),
        blocks: int(detailed.blocks).unwrap_or(0),
        block_n50: int(detailed.block_n50),
    })
}

/// `run_stats`, printing the statistics to `out`
pub fn run_stats(opts: &StatsOptions, out: &mut impl Write) -> Result<()> {
    let chromosomes = unpack_chromosomes(&opts.chromosomes);
    let mut gtfwriter = match &opts.gtf {
        Some(path) => Some(GtfWriter { file: create(path)? }),
        None => None,
    };
    let reader = VcfReader::open_with(&opts.vcf, opts.only_snvs, true, None)?;
    if reader.samples.is_empty() {
        log::error!("Input VCF does not contain any sample");
        return Ok(());
    }
    log::info!("Found {} sample(s) in input VCF", reader.samples.len());
    let sample = match &opts.sample {
        Some(sample) if reader.samples.contains(sample) => sample.clone(),
        Some(sample) => {
            log::error!("Requested sample ({sample}) not found");
            return Ok(());
        }
        None => {
            let sample = reader.samples[0].clone();
            log::info!("Reporting results for sample {sample}");
            sample
        }
    };
    let chr_lengths = chr_lengths(&reader, opts.chr_lengths.as_deref())?;

    let mut tsv_file = match &opts.tsv {
        Some(path) => {
            let mut f = create(path)?;
            writeln!(f, "#sample\tchromosome\tfile_name\t{}", FIELD_NAMES.join("\t"))?;
            Some(f)
        }
        None => None,
    };
    let mut block_list_file = match &opts.block_list {
        Some(path) => {
            let mut f = create(path)?;
            writeln!(f, "#sample\tchromosome\tphase_set\tfrom\tto\tvariants")?;
            Some(f)
        }
        None => None,
    };
    let write_tsv = |f: &mut Option<BufWriter<File>>, chromosome: &str, stats: &DetailedStats| -> Result<()> {
        if let Some(f) = f {
            let values: Vec<String> = stats.values().iter().map(|v| v.py_str()).collect();
            writeln!(f, "{sample}\t{chromosome}\t{}\t{}", opts.vcf, values.join("\t"))?;
        }
        Ok(())
    };

    writeln!(out, "Phasing statistics for sample {sample} from file {}", opts.vcf)?;
    let fetch = if !chromosomes.is_empty() {
        let indexed = IndexedVcf::open(&opts.vcf)?;
        indexed.has_index().then(|| (indexed, chromosomes.clone().into_iter()))
    } else {
        None
    };
    let mut tables = VariantTables { reader, fetch };
    let mut total_stats = PhasingStats::default();
    let given: HashSet<&str> = chromosomes.iter().map(String::as_str).collect();
    let mut seen_chromosomes: HashSet<String> = HashSet::new();
    while let Some(table) = tables.next()? {
        let chromosome: Arc<str> = Arc::from(table.chromosome.as_str());
        seen_chromosomes.insert(table.chromosome.clone());
        if !given.is_empty() && !given.contains(&*chromosome) {
            continue;
        }
        let mut stats = PhasingStats::default();
        writeln!(out, "---------------- Chromosome {chromosome} ----------------")?;
        let blocks = phase_blocks(&chromosome, &mut gtfwriter, &sample, &mut stats, &table)?;
        if let Some(f) = &mut block_list_file {
            for i in blocks.sorted_indices()? {
                let block = &blocks.blocks[i];
                writeln!(
                    f,
                    "{sample}\t{chromosome}\t{}\t{}\t{}\t{}",
                    block_id_str(&blocks.keys[i]),
                    i64::from(block.leftmost()) + 1,
                    i64::from(block.rightmost()) + 1,
                    block.len()
                )?;
            }
        }
        stats.add_blocks(&chromosome, &blocks.blocks);
        let detailed_stats = stats.detailed_stats(&chr_lengths);
        detailed_stats.print(out)?;
        write_tsv(&mut tsv_file, &chromosome, &detailed_stats)?;
        total_stats.extend(stats);
        if !given.is_empty() && given.iter().all(|c| seen_chromosomes.contains(*c)) {
            break;
        }
    }
    if seen_chromosomes.len() > 1 {
        writeln!(out, "---------------- ALL chromosomes (aggregated) ----------------")?;
        let detailed_stats = total_stats.detailed_stats(&chr_lengths);
        detailed_stats.print(out)?;
        write_tsv(&mut tsv_file, "ALL", &detailed_stats)?;
    }
    for mut f in [gtfwriter.map(|g| g.file), tsv_file, block_list_file].into_iter().flatten() {
        f.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn block(positions: &[i32]) -> PhasedBlock {
        PhasedBlock { positions: positions.to_vec(), snvs: 0 }
    }

    fn spans(blocks: &[PhasedBlock]) -> Vec<(i32, i32)> {
        blocks.iter().map(|b| (b.leftmost(), b.rightmost())).collect()
    }

    #[test]
    fn interleaved_blocks_are_split() {
        // A: 1 3 20 30, B: 5 10 -> A splits into [1 3] and [20 30]
        let split = nonoverlapping_blocks(&[block(&[1, 3, 20, 30]), block(&[5, 10])]);
        assert_eq!(spans(&split), vec![(1, 3), (5, 10), (20, 30)]);
        // singleton remainders are dropped
        let split = nonoverlapping_blocks(&[block(&[1, 10]), block(&[5, 6]), block(&[7])]);
        assert_eq!(spans(&split), vec![(5, 6)]);
    }

    #[test]
    fn totals_count_the_phased_het_snvs_and_the_blocks_ng50() {
        // Blocks spanning 100 (with an indel), 400 and, on chr2, 100 bases; an unphased het SNV,
        // a homozygous one and a singleton block. WhatsHap 2.8's `stats` gives the same.
        let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=950>\n##contig=<ID=chr2,length=200>\n\
                   ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                   ##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
                   chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:10\n\
                   chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT:PS\t1|0:10\n\
                   chr1\t30\t.\tA\tAT\t.\tPASS\t.\tGT:PS\t0|1:10\n\
                   chr1\t110\t.\tG\tA\t.\tPASS\t.\tGT:PS\t0|1:10\n\
                   chr1\t200\t.\tT\tC\t.\tPASS\t.\tGT:PS\t0|1:200\n\
                   chr1\t600\t.\tA\tC\t.\tPASS\t.\tGT:PS\t1|0:200\n\
                   chr1\t700\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
                   chr1\t800\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n\
                   chr1\t900\t.\tC\tG\t.\tPASS\t.\tGT:PS\t0|1:900\n\
                   chr2\t50\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:50\n\
                   chr2\t150\t.\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:50\n";
        let path = std::env::temp_dir().join(format!("lrpipe-phase-totals-{}.vcf", std::process::id()));
        std::fs::write(&path, vcf).expect("written");
        let totals = totals(path.to_str().expect("UTF-8")).expect("read");
        let _ = std::fs::remove_file(&path);
        // Spans 400, 100 and 100 against 1150 bases: the third reaches half.
        assert_eq!(totals, Totals { heterozygous_snvs: 9, phased_snvs: 7, blocks: 3, block_n50: Some(100) });
    }

    #[test]
    fn n50_and_median() {
        assert_eq!(n50(&mut [10, 5, 30], 100), 0);
        assert_eq!(n50(&mut [10, 5, 30], 60), 30);
        assert_eq!(n50(&mut [10, 5, 30], 80), 10);
        assert_eq!(n50(&mut [10, 5, 30], 81), 5);
        assert_eq!(median(&[1, 2, 7]), Num::Int(2));
        assert_eq!(median(&[1, 2]), Num::Float(1.5));
    }
}

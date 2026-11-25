//! gVCF output: the PL of variant records (`compute_PL`, `CallVariants.py`) and the
//! non-variant blocks of a pileup chunk (`variantInfoCalculator`, `preprocess/utils.py`).

use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;

use crate::{Error, gvcf_header};

/// `IUPAC_base_to_ACGT_base_dict`.
fn iupac_to_acgt(base: char) -> Option<char> {
    let i = "ACGTURYSWKMBDHVN".find(base)?;
    "ACGTTACCAGACAAAA".chars().nth(i)
}

/// `GT21_LABELS`, in the order of the network's GT21 outputs.
const GT21_LABELS: [&str; 21] = [
    "AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT", "DelDel", "ADel", "CDel", "GDel", "TDel", "InsIns",
    "AIns", "CIns", "GIns", "TIns", "InsDel",
];

/// `partial_label_from`.
fn partial_label(reference: &str, alternate: &str) -> Option<String> {
    let (r, a) = (reference.chars().count(), alternate.chars().count());
    Some(if r > a {
        "Del".to_owned()
    } else if r < a {
        "Ins".to_owned()
    } else {
        alternate.chars().next()?.to_string()
    })
}

/// `mix_two_partial_labels`.
fn mix_labels(label1: &str, label2: &str) -> String {
    let (n1, n2) = (label1.chars().count(), label2.chars().count());
    if n1 == 1 && n2 == 1 {
        return if label1 <= label2 { format!("{label1}{label2}") } else { format!("{label2}{label1}") };
    }
    let (single, multi) = if n1 > 1 && n2 == 1 { (label2, label1) } else { (label1, label2) };
    if multi.chars().count() > 1 && single.chars().count() == 1 {
        return format!("{single}{multi}");
    }
    if n1 > 0 && n2 > 0 && label1 == label2 {
        return format!("{label1}{label2}");
    }
    "InsDel".to_owned()
}

/// The genotype head index of a genotype (`genotype_enum_for_task(genotype_enum_from(..))`):
/// 0/0, 1/1, then every heterozygous genotype.
fn zygosity(g1: usize, g2: usize) -> usize {
    match (g1, g2) {
        (0, 0) => 0,
        _ if g1 == g2 => 1,
        _ => 2,
    }
}

/// `compute_PL`: the PL of a record, from the raw GT21 and genotype probabilities.
///
/// Each genotype's likelihood is its GT21 probability times its zygosity probability, in
/// `float32` as numpy computes it; they are normalised to their (`float32`) sum, offset by
/// `1e-8`, Phred-scaled in `float64` and shifted so the smallest is 0 (rounded up). A genotype
/// without a GT21 class (a RefCall's `.`, an `N` allele) gives PL 990 for every genotype, or a
/// single 990 for a RefCall. `None` is where the original raises: a reference base outside the
/// IUPAC table or more than two alternate alleles.
pub(crate) fn genotype_likelihoods(
    reference: &str,
    alternate: &str,
    gt21: &[f32],
    genotype: &[f32],
) -> Option<Vec<i64>> {
    let alternates: Vec<&str> = alternate.split(',').collect();
    let genotypes: &[(usize, usize)] = match alternates.len() {
        1 => &[(0, 0), (0, 1), (1, 1)],
        2 => &[(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)],
        _ => return None,
    };
    let mut chars = reference.chars();
    let reference = match (chars.next(), chars.next()) {
        (Some(base), None) => iupac_to_acgt(base)?.to_string(),
        _ => reference.to_owned(),
    };
    let alleles: Vec<&str> = std::iter::once(reference.as_str()).chain(alternates).collect();
    let mut likelihoods: Vec<f32> = Vec::with_capacity(genotypes.len());
    for &(g1, g2) in genotypes {
        let label = mix_labels(&partial_label(&reference, alleles[g1])?, &partial_label(&reference, alleles[g2])?);
        let Some(index) = GT21_LABELS.iter().position(|&l| l == label) else {
            return Some(if alternate == "." { vec![990] } else { vec![990; genotypes.len()] });
        };
        likelihoods.push(gt21[index] * genotype[zygosity(g1, g2)]);
    }
    let sum = likelihoods.iter().fold(0f32, |acc, &p| acc + p);
    let ln_10 = 10f64.ln();
    let phred: Vec<f64> = likelihoods.iter().map(|&p| -10.0 * (f64::from(p / sum + 1e-8f32).ln() / ln_10)).collect();
    let min = phred.iter().copied().fold(f64::INFINITY, f64::min);
    Some(phred.iter().map(|&p| (p - min).ceil() as i64).collect())
}

/// `LOG_10` in `preprocess/utils.py`: an approximation of ln 10 that the original divides by.
#[expect(clippy::approx_constant, reason = "the original's own approximation, which its output depends on")]
const LOG_10: f64 = 2.3025;
/// `LOG_2`: an approximation of log10(2).
#[expect(clippy::approx_constant, reason = "the original's own approximation, which its output depends on")]
const LOG_2: f64 = 0.3010;
/// `self.max_gq`.
const MAX_GQ: i64 = 50;

/// Python's `round(x, 6)`: correctly rounded, like Rust's fixed-precision formatting.
fn round6(x: f64) -> f64 {
    format!("{x:.6}").parse().expect("a formatted float parses")
}

/// `log10sumexp` (C, through cffi): `m + log(sum(10^(x - m))) / LOG_10`.
///
/// `getMyMaxItem` reads one element past the array (`i <= n_list`); in the original that value is
/// never a valid maximum and behaves as 0 against the negative log-probabilities, so `m` is
/// taken as `max(x..., 0)`. Because `LOG_10` is not ln 10, the result depends on `m`.
fn log10_sum_exp(values: &[f64; 3]) -> f64 {
    let m = values.iter().copied().fold(0.0, |m: f64, x| if x > m { x } else { m });
    let sum: f64 = values.iter().map(|&x| 10f64.powf(x - m)).sum();
    m + sum.ln() / LOG_10
}

/// `normalize_log10_prob`: each value minus the rounded log-sum, capped at 0.
fn normalize(values: [f64; 3]) -> [f64; 3] {
    let lse = round6(log10_sum_exp(&values));
    values.map(|x| {
        let x = x - lse;
        if 0.0 < x { 0.0 } else { x }
    })
}

/// `log10p_to_phred` (C, through cffi, then `round(.., 6)`): -10 log10(1 - 10^p), with the
/// original's `LOG_10`, or 50 when 10^p is exactly 1.
fn log10p_to_phred(log10p: f64) -> f64 {
    let p = 10f64.powf(log10p);
    round6(if p == 1.0 { 50.0 } else { -10.0 * ((1.0 - p).ln() / LOG_10) })
}

/// Python's `a // b`.
fn floor_div(a: i64, b: i64) -> i64 {
    let q = a / b;
    if (a % b != 0) && ((a < 0) != (b < 0)) { q - 1 } else { q }
}

/// The genotype quality and likelihoods of a non-variant site (`_cal_reference_likelihood`).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Likelihood {
    /// Whether 0/0 is the most likely genotype (`validPL`).
    valid: bool,
    gq: i64,
    binned_gq: i64,
    log10_probs: [f64; 3],
}

/// The constants of `variantInfoCalculator`.
#[derive(Debug, Clone, PartialEq)]
struct Model {
    /// `math.log(p_err) / LOG_10` (the original computes it with mpmath, correctly rounded).
    log_p: f64,
    /// `math.log1p(-p_err) / LOG_10`.
    log_1mp: f64,
    gq_bin_size: i64,
    /// `normalize_log10_prob([-1.0, -1.0, -1.0])`, used at zero coverage.
    no_coverage: [f64; 3],
}

impl Model {
    fn new(base_error: f64, gq_bin_size: i64) -> Self {
        Self {
            log_p: base_error.ln() / LOG_10,
            log_1mp: (-base_error).ln_1p() / LOG_10,
            gq_bin_size,
            no_coverage: normalize([-1.0; 3]),
        }
    }

    fn likelihood(&self, n_ref: u64, n_total: u64) -> Likelihood {
        let log10_probs = if n_total == 0 {
            self.no_coverage
        } else {
            // Counts are numpy int64; mpmath and numpy multiply them exactly before rounding.
            let (n_ref, n_alt, n_total) = (n_ref as f64, n_total as f64 - n_ref as f64, n_total as f64);
            let p_ref = n_ref * self.log_1mp + n_alt * self.log_p;
            let p_het = -n_total * LOG_2;
            let p_hom_alt = n_ref * self.log_p + n_alt * self.log_1mp;
            normalize([p_ref, p_het, p_hom_alt])
        };
        let gq = (log10p_to_phred(log10_probs[0]) as i64).min(MAX_GQ);
        let binned_gq = if gq >= 1 { floor_div(gq - 1, self.gq_bin_size) * self.gq_bin_size + 1 } else { 0 };
        let max = log10_probs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Likelihood { valid: log10_probs[0] == max, gq, binned_gq, log10_probs }
    }
}

/// How non-variant sites are written (`--bp_resolution`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Resolution {
    /// Runs of similar sites become one block (the default).
    #[default]
    Blocks,
    /// Every site is its own record, except runs of `N` reference bases.
    PerBase,
}

/// Parameters of non-variant gVCF output.
#[derive(Debug, Clone, PartialEq)]
pub struct NonVariantOptions {
    /// Estimated base error rate (`--base_err`, default 0.001).
    pub base_error: f64,
    /// Width of the GQ bins that split blocks (`--gq_bin_size`, default 5).
    pub gq_bin_size: i64,
    /// Blocks or one record per site (`--bp_resolution`).
    pub resolution: Resolution,
}

impl Default for NonVariantOptions {
    fn default() -> Self {
        Self { base_error: 0.001, gq_bin_size: 5, resolution: Resolution::Blocks }
    }
}

/// One non-variant record (`non_variant_info`).
#[derive(Debug, Clone, PartialEq)]
struct Record {
    pos: i64,
    end: i64,
    reference: char,
    genotype: &'static str,
    /// Written in the GQ column (`binned_gq`).
    gq: i64,
    /// The unbinned GQ (`gq`).
    raw_gq: i64,
    min_dp: u64,
    pl: [i64; 3],
}

/// The block being extended (`current_block` and the `cur_*` fields).
#[derive(Debug)]
struct Block {
    records: Vec<Record>,
    gq_bin: i64,
    genotype: &'static str,
    min_dp: u64,
    max_dp: u64,
    min_raw_gq: i64,
}

impl Block {
    fn new(record: Record) -> Self {
        Self {
            gq_bin: record.gq,
            genotype: record.genotype,
            min_dp: record.min_dp,
            max_dp: record.min_dp,
            min_raw_gq: record.raw_gq,
            records: vec![record],
        }
    }

    fn push(&mut self, record: Record) {
        self.min_raw_gq = self.min_raw_gq.min(record.raw_gq);
        self.records.push(record);
    }
}

/// `ceil(dp + dp * 0.3)`: how far a block's depth may spread above its minimum.
fn depth_limit(min_dp: u64) -> f64 {
    let dp = min_dp as f64;
    (dp + dp * 0.3).ceil()
}

/// Writes the non-variant records of one pileup chunk, as `variantInfoCalculator` builds them
/// from the sites `CreateTensorPileup` feeds it (`make_gvcf_online`).
///
/// Sites are pushed in position order. Consecutive sites join a block while they share the GQ
/// bin and the genotype (0/0 when it is the most likely genotype, `./.` otherwise), neither
/// side of a reference base change is `N`, and the block's depths stay within 30% of its
/// minimum. A block is written as one record from its first to its last site, with its
/// smallest GQ; blocks of `./.` sites are written site by site.
pub struct NonVariantWriter<W: Write> {
    out: W,
    contig: String,
    model: Model,
    resolution: Resolution,
    contig_lengths: HashMap<String, i64>,
    interval_end: Option<i64>,
    block: Option<Block>,
    /// `cur_ref`: the reference base as given (not turned into `N`), which the original does not
    /// update when a depth change starts a new block.
    block_reference: u8,
}

impl<W: Write> std::fmt::Debug for NonVariantWriter<W> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NonVariantWriter")
            .field("contig", &self.contig)
            .field("model", &self.model)
            .field("resolution", &self.resolution)
            .field("block", &self.block)
            .finish_non_exhaustive()
    }
}

/// The file name stem the original gives a pileup chunk's non-variant gVCF:
/// `{sample}.{contig}.{start}.{end}`, with `.tmp.gvcf` appended.
pub fn chunk_gvcf_name(sample: &str, contig: &str, start: i64, end: i64) -> String {
    format!("{sample}.{contig}.{start}.{end}")
}

impl<W: Write> NonVariantWriter<W> {
    /// Writes the gVCF header to `out` and returns a writer for `contig`'s sites.
    ///
    /// `chunk_name` is the [`chunk_gvcf_name`] of the chunk. Its last two fields and the known
    /// contig suffix are removed to recover the sample, even when the sample or contig has dots.
    /// `reference` provides `##reference`, the `##contig` lines and the contig lengths.
    ///
    /// # Errors
    ///
    /// Fails if the reference's `.fai` cannot be read, `gq_bin_size` is 0 (the original divides by
    /// it at the first site), or writing fails.
    pub fn new(
        mut out: W,
        reference: &Path,
        chunk_name: &str,
        contig: &str,
        options: &NonVariantOptions,
    ) -> Result<Self, Error> {
        if options.gq_bin_size == 0 {
            return Err(Error::GqBinSize);
        }
        let sample = chunk_name
            .rsplitn(3, '.')
            .nth(2)
            .and_then(|stem| stem.strip_suffix(contig))
            .and_then(|stem| stem.strip_suffix('.'))
            .unwrap_or_else(|| chunk_name.rsplitn(4, '.').last().unwrap_or(chunk_name));
        let (header, contig_lengths) = gvcf_header(reference, sample)?;
        out.write_all(header.as_bytes()).map_err(Error::Write)?;
        Ok(Self {
            out,
            contig: contig.to_owned(),
            model: Model::new(options.base_error, options.gq_bin_size),
            resolution: options.resolution,
            contig_lengths,
            interval_end: None,
            block: None,
            block_reference: b'N',
        })
    }

    /// Flushes the previous block before writing sites in a new interval. The last position
    /// (1-based, inclusive) also limits extension to the end of the contig.
    ///
    /// # Errors
    /// Fails if writing the previous block fails.
    pub fn start_interval(&mut self, last_pos: i64) -> io::Result<()> {
        self.flush_block()?;
        self.interval_end = Some(last_pos);
        Ok(())
    }

    /// `reference_likelihood`: the record of one site.
    fn site_record(&self, pos: i64, reference: u8, n_ref: u64, n_total: u64) -> Record {
        let likelihood = self.model.likelihood(n_ref, n_total);
        let genotype = if likelihood.valid { "0/0" } else { "./." };
        let phred = likelihood.log10_probs.map(|x| -10.0 * x);
        let min = phred.iter().copied().fold(f64::INFINITY, f64::min);
        let (reference, gq, binned_gq, pl) = match reference {
            b'A' | b'C' | b'G' | b'T' => {
                (char::from(reference), likelihood.gq, likelihood.binned_gq, phred.map(|x| (x - min) as i64))
            }
            _ => ('N', 1, 1, [0; 3]),
        };
        Record { pos, end: pos, reference, genotype, gq: binned_gq, raw_gq: gq, min_dp: n_total, pl }
    }

    /// Adds the site at 1-based `pos` with reference base `reference` (uppercase), `n_ref`
    /// reads supporting it and `n_total` reads in all.
    ///
    /// # Errors
    ///
    /// Fails if writing a finished block fails.
    pub fn push(&mut self, pos: i64, reference: u8, n_ref: u64, n_total: u64) -> io::Result<()> {
        let record = self.site_record(pos, reference, n_ref, n_total);
        let Some(block) = &mut self.block else {
            self.block_reference = reference;
            self.block = Some(Block::new(record));
            return Ok(());
        };
        // Only a literal `N` on either side splits: an IUPAC code such as `R` is written as `N`
        // but joins the neighbouring block.
        let n_boundary = reference != self.block_reference && (reference == b'N' || self.block_reference == b'N');
        if record.gq != block.gq_bin || record.genotype != block.genotype || n_boundary {
            self.block_reference = reference;
            return self.start_block(record);
        }
        let dp = record.min_dp;
        if dp < block.min_dp {
            if block.max_dp as f64 > depth_limit(dp) {
                return self.start_block(record);
            }
            block.min_dp = dp;
            block.push(record);
        } else if dp > block.max_dp {
            if dp as f64 > depth_limit(block.min_dp) {
                return self.start_block(record);
            }
            block.max_dp = dp;
            block.push(record);
        } else {
            block.push(record);
        }
        Ok(())
    }

    /// Writes the current block and starts a new one at `record`.
    fn start_block(&mut self, record: Record) -> io::Result<()> {
        if let Some(block) = self.block.replace(Block::new(record)) {
            self.write_block(&block)?;
        }
        Ok(())
    }

    /// `write_to_gvcf_batch`.
    fn write_block(&mut self, block: &Block) -> io::Result<()> {
        let first = &block.records[0];
        let per_site = self.resolution == Resolution::PerBase || first.genotype == "./.";
        if per_site && first.reference != 'N' {
            for record in &block.records {
                self.write_record(record)?;
            }
            return Ok(());
        }
        let end = block.records.last().map_or(first.pos, |r| r.pos);
        let record = if first.reference == 'N' {
            Record { end, genotype: "./.", gq: 1, raw_gq: 1, min_dp: block.min_dp, pl: [0; 3], ..first.clone() }
        } else {
            Record { end, gq: block.min_raw_gq, min_dp: block.min_dp, ..first.clone() }
        };
        self.write_record(&record)
    }

    /// `write_to_gvcf`. An END one short of the contig length is written as the length
    /// if it stays within the current interval.
    fn write_record(&mut self, r: &Record) -> io::Result<()> {
        let end = match self.contig_lengths.get(&self.contig) {
            Some(&len) if r.end == len - 1 && self.interval_end.is_none_or(|end| len <= end) => len,
            _ => r.end,
        };
        writeln!(
            self.out,
            "{}\t{}\t.\t{}\t<NON_REF>\t0\t.\tEND={end}\tGT:GQ:MIN_DP:PL\t{}:{}:{}:{},{},{}",
            self.contig, r.pos, r.reference, r.genotype, r.gq, r.min_dp, r.pl[0], r.pl[1], r.pl[2]
        )
    }

    /// `write_empty_pileup`: one `N` record for a chunk without any read, whose sites are
    /// `start..end` (`start` at least 1): to `end - 1`, or the contig's end if that comes first.
    /// The original's runs to `end`, so it overlaps the next chunk's first site, and on the last
    /// chunk past the contig's end.
    ///
    /// # Errors
    ///
    /// Fails if writing fails.
    pub fn write_empty(&mut self, start: i64, end: i64) -> io::Result<()> {
        let start = start.max(1);
        let end = self.contig_lengths.get(&self.contig).map_or(end - 1, |&len| (end - 1).min(len));
        if end < start {
            return Ok(());
        }
        let record =
            Record { pos: start, end, reference: 'N', genotype: "./.", gq: 1, raw_gq: 1, min_dp: 0, pl: [0; 3] };
        self.write_record(&record)
    }

    /// Writes the last block and returns the inner writer.
    ///
    /// # Errors
    ///
    /// Fails if writing fails.
    pub fn finish(mut self) -> io::Result<W> {
        self.flush_block()?;
        Ok(self.out)
    }

    fn flush_block(&mut self) -> io::Result<()> {
        if let Some(block) = self.block.take() {
            self.write_block(&block)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_mix_like_the_original() {
        assert_eq!(mix_labels("G", "A"), "AG");
        assert_eq!(mix_labels("Del", "A"), "ADel");
        assert_eq!(mix_labels("A", "Ins"), "AIns");
        assert_eq!(mix_labels("Ins", "Ins"), "InsIns");
        assert_eq!(mix_labels("Ins", "Del"), "InsDel");
        assert_eq!(mix_labels(".", "A"), ".A");
    }

    #[test]
    fn genotype_likelihoods_match_compute_pl() {
        // The original's `compute_PL` on random float32 probabilities.
        let rows = include_str!("../tests/data/compute_pl.tsv");
        for row in rows.lines().filter(|r| !r.starts_with('#')) {
            let f: Vec<&str> = row.split('\t').collect();
            let floats = |s: &str| -> Vec<f32> { s.split(',').map(|x| x.parse::<f64>().unwrap() as f32).collect() };
            let want: Vec<i64> = f[4].split(',').map(|x| x.parse().unwrap()).collect();
            assert_eq!(genotype_likelihoods(f[0], f[1], &floats(f[2]), &floats(f[3])), Some(want), "{row}");
        }
    }

    #[test]
    fn a_chunk_without_reads_is_one_block_of_its_own_sites() {
        // The sites of chunks of chrT (length 10) are 1..4, 4..8 and 8..12: the blocks end
        // before the next chunk's first site, and at the contig's end.
        let blocks = |chunks: &[(i64, i64)]| {
            let mut writer = NonVariantWriter {
                out: Vec::new(),
                contig: "chrT".to_owned(),
                model: Model::new(0.001, 5),
                resolution: Resolution::Blocks,
                contig_lengths: HashMap::from([("chrT".to_owned(), 10)]),
                interval_end: None,
                block: None,
                block_reference: b'N',
            };
            for &(start, end) in chunks {
                writer.write_empty(start, end).unwrap();
            }
            let text = String::from_utf8(writer.finish().unwrap()).unwrap();
            text.lines().map(|l| l.split('\t').take(8).collect::<Vec<_>>().join(" ")).collect::<Vec<_>>()
        };
        assert_eq!(
            blocks(&[(0, 4), (4, 8), (8, 12)]),
            ["chrT 1 . N <NON_REF> 0 . END=3", "chrT 4 . N <NON_REF> 0 . END=7", "chrT 8 . N <NON_REF> 0 . END=10"]
        );
        assert!(blocks(&[(11, 15)]).is_empty(), "a chunk past the contig's end has no sites");
    }

    #[test]
    fn python_floor_division() {
        assert_eq!([floor_div(7, 5), floor_div(-1, 5), floor_div(7, -5), floor_div(0, 3)], [1, -1, -2, 0]);
    }

    #[test]
    fn round6_is_half_even_on_the_exact_value() {
        // 1/128 = 0.0078125 exactly, a tie at six digits.
        assert_eq!(round6(0.0078125), 0.007812);
        assert_eq!(round6(-0.5228977), -0.522898);
    }

    /// Compares with rows of `p_err n_ref n_total gt gq binned_gq pl0 pl1 pl2` from the original's
    /// `variantInfoCalculator.reference_likelihood`.
    fn check_likelihoods(rows: &str) {
        let mut models: HashMap<String, Model> = HashMap::new();
        for row in rows.lines().filter(|r| !r.starts_with('#')) {
            let f: Vec<&str> = row.split('\t').collect();
            let model = models.entry(f[0].to_owned()).or_insert_with(|| Model::new(f[0].parse().unwrap(), 5));
            let int = |i: usize| f[i].parse::<i64>().unwrap();
            let (n_ref, n_total) = (int(1) as u64, int(2) as u64);
            let writer = NonVariantWriter {
                out: Vec::new(),
                contig: "chr1".to_owned(),
                model: model.clone(),
                resolution: Resolution::Blocks,
                contig_lengths: HashMap::new(),
                interval_end: None,
                block: None,
                block_reference: b'N',
            };
            let record = writer.site_record(1, b'A', n_ref, n_total);
            let want = (f[3], int(4), int(5), [int(6), int(7), int(8)]);
            assert_eq!((record.genotype, record.raw_gq, record.gq, record.pl), want, "{row}");
        }
    }

    #[test]
    fn site_likelihoods_match_the_original() {
        check_likelihoods(include_str!("../tests/data/reference_likelihoods.tsv"));
        // A larger grid (`SNV_GVCF_GRID`), from the same script.
        if let Ok(path) = std::env::var("SNV_GVCF_GRID") {
            check_likelihoods(&std::fs::read_to_string(path).unwrap());
        }
    }

    /// `tests/data/non_variant/<case>.sites.tsv` lists `pos ref n_ref n_total` after a
    /// `#name bp_resolution p_err gq_bin_size` line; `<case>.expected.gvcf` is what the original's
    /// `variantInfoCalculator` writes for them (`make_gvcf_online`, then the last block).
    #[test]
    fn blocks_match_variant_info_calculator() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/non_variant");
        for case in ["blocks", "per_base", "other_parameters"] {
            let sites = std::fs::read_to_string(data.join(format!("{case}.sites.tsv"))).unwrap();
            let mut lines = sites.lines();
            let params: Vec<&str> = lines.next().unwrap().trim_start_matches('#').split('\t').collect();
            let options = NonVariantOptions {
                base_error: params[2].parse().unwrap(),
                gq_bin_size: params[3].parse().unwrap(),
                resolution: if params[1] == "1" { Resolution::PerBase } else { Resolution::Blocks },
            };
            let mut writer =
                NonVariantWriter::new(Vec::new(), &data.join("ref.fa"), params[0], "chrT", &options).unwrap();
            for line in lines {
                let f: Vec<&str> = line.split('\t').collect();
                writer
                    .push(f[0].parse().unwrap(), f[1].as_bytes()[0], f[2].parse().unwrap(), f[3].parse().unwrap())
                    .unwrap();
            }
            // `##reference` names the path the fixture was made with.
            let without_reference = |text: String| -> String {
                text.split_inclusive('\n').filter(|l| !l.starts_with("##reference=")).collect()
            };
            let got = without_reference(String::from_utf8(writer.finish().unwrap()).unwrap());
            let want = without_reference(std::fs::read_to_string(data.join(format!("{case}.expected.gvcf"))).unwrap());
            assert_eq!(got, want, "case {case}");
        }
    }

    #[test]
    fn a_dotted_contig_name_does_not_change_the_gvcf_sample() {
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/non_variant/ref.fa");
        let sample = "isolate.v1";
        let contig = "NC_000913.3";
        let chunk = chunk_gvcf_name(sample, contig, 1, 10);
        let options = NonVariantOptions { base_error: 0.001, gq_bin_size: 5, resolution: Resolution::Blocks };
        let writer = NonVariantWriter::new(Vec::new(), &reference, &chunk, contig, &options).expect("gVCF header");
        let header = String::from_utf8(writer.finish().expect("finished")).expect("UTF-8");
        let columns = header.lines().find(|line| line.starts_with("#CHROM\t")).expect("column names");
        assert_eq!(columns.split('\t').next_back(), Some(sample));
    }

    #[test]
    fn zero_coverage_constants() {
        // `variantInfoCalculator(...).constant_log10_probs` in the original: [-0.477102] * 3.
        let model = Model::new(0.001, 5);
        assert_eq!(model.no_coverage, [-0.477102; 3]);
        assert_eq!(
            model.likelihood(0, 0),
            Likelihood { valid: true, gq: 1, binned_gq: 1, log10_probs: [-0.477102; 3] }
        );
    }
}

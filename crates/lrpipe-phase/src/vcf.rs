//! Port of `vcf.py` (the parts used by `whatshap phase`).

use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result, bail};

use crate::compat::pystr::py_int;
use crate::core::genotype::{Genotype, PhredGenotypeLikelihoods};
use crate::core::read::{Read, ReadSet};
use crate::hts::{
    self, BCF_HL_FMT, BCF_HL_INFO, Header, IndexedVcf, PyFormatValue, Record, VcfFetchError, VcfFile, VcfWriter,
};

/// A variant in a VCF file (not to be confused with a variant on a read).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct VcfVariant {
    pub position: i32,
    pub reference_allele: String,
    pub alternative_alleles: Vec<String>,
}

impl VcfVariant {
    pub fn get_allele(&self, a: usize) -> &str {
        if a == 0 { &self.reference_allele } else { &self.alternative_alleles[a - 1] }
    }

    pub fn is_snv(&self) -> bool {
        self.alternative_alleles.iter().any(|alt| *alt != self.reference_allele)
            && self.reference_allele.len() == 1
            && self.alternative_alleles.iter().all(|alt| alt.len() == 1)
    }

    /// Remove common suffix and prefix (`BiallelicVcfVariant.normalized` and
    /// the multiallelic variant of it).
    pub fn normalized(&self) -> VcfVariant {
        let mut pos = self.position;
        let mut r: &str = &self.reference_allele;
        let mut alts: Vec<&str> = self.alternative_alleles.iter().map(|s| s.as_str()).collect();
        if alts.len() == 1 {
            let mut a: &str = alts[0];
            while !r.is_empty() && !a.is_empty() && r.as_bytes()[r.len() - 1] == a.as_bytes()[a.len() - 1] {
                r = &r[..r.len() - 1];
                a = &a[..a.len() - 1];
            }
            while !r.is_empty() && !a.is_empty() && r.as_bytes()[0] == a.as_bytes()[0] {
                r = &r[1..];
                a = &a[1..];
                pos += 1;
            }
            alts = vec![a];
        } else {
            while !r.is_empty()
                && alts.iter().all(|a| !a.is_empty())
                && alts.iter().all(|a| r.as_bytes()[r.len() - 1] == a.as_bytes()[a.len() - 1])
            {
                r = &r[..r.len() - 1];
                alts = alts.iter().map(|a| &a[..a.len() - 1]).collect();
            }
            while !r.is_empty()
                && alts.iter().all(|a| !a.is_empty())
                && alts.iter().all(|a| r.as_bytes()[0] == a.as_bytes()[0])
            {
                r = &r[1..];
                alts = alts.iter().map(|a| &a[1..]).collect();
                pos += 1;
            }
        }
        VcfVariant {
            position: pos,
            reference_allele: r.to_string(),
            alternative_alleles: alts.into_iter().map(|s| s.to_string()).collect(),
        }
    }
}

/// A ploidy that is too high or inconsistent (`vcf.PloidyError`).
#[derive(Debug)]
pub struct PloidyError(pub String);

impl std::fmt::Display for PloidyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for PloidyError {}

/// Identity of a block id as a key of a Python dict: ints and floats with
/// equal values are the same key, and every NaN is a key of its own (pysam
/// creates a new float object for each value).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum BlockKey {
    None,
    Num(i64),
    Fraction(u64),
    Str(String),
    Nan(usize),
}

/// The dict key of `block_id`; `nan_counter` numbers the NaNs.
pub fn block_key(id: &Option<BlockId>, nan_counter: &mut usize) -> BlockKey {
    match id {
        None => BlockKey::None,
        Some(BlockId::Int(v)) => BlockKey::Num(*v),
        Some(BlockId::Float(v)) if v.is_nan() => {
            *nan_counter += 1;
            BlockKey::Nan(*nan_counter)
        }
        Some(BlockId::Float(v)) if v.fract() == 0.0 && v.abs() < 9e18 => BlockKey::Num(*v as i64),
        Some(BlockId::Float(v)) => BlockKey::Fraction(v.to_bits()),
        Some(BlockId::Str(s)) => BlockKey::Str(s.clone()),
    }
}

/// For a single chromosome, variants (rows) and genotypes per sample (columns).
/// A phase set id as pysam returns it: its Python type follows the header
/// type of PS (Integer, Float, or String if PS is undefined in the header).
#[derive(Clone, Debug, PartialEq)]
pub enum BlockId {
    Int(i64),
    Float(f64),
    Str(String),
}

impl BlockId {
    /// Python's `int(block_id)`
    pub fn to_int(&self) -> Result<i64> {
        match self {
            BlockId::Int(v) => Ok(*v),
            BlockId::Float(v) => Ok(v.trunc() as i64),
            BlockId::Str(s) => py_int(s)
                .and_then(|v| i64::try_from(v).ok())
                .ok_or_else(|| anyhow::anyhow!("invalid literal for int() with base 10: '{s}'")),
        }
    }
}

impl std::fmt::Display for BlockId {
    /// Python's `str(block_id)`
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BlockId::Int(v) => write!(f, "{v}"),
            BlockId::Float(v) => write!(f, "{}", crate::compat::pyfloat::repr(*v)),
            BlockId::Str(s) => write!(f, "{s}"),
        }
    }
}

/// Phase of a variant call: block id (PS, None if missing) and the alleles
/// in haplotype order (`VariantCallPhase`).
#[derive(Clone, Debug, PartialEq)]
pub struct VariantCallPhase {
    pub block_id: Option<BlockId>,
    pub phase: Vec<Option<u32>>,
    /// PQ value, if present
    pub quality: Option<f64>,
}

#[derive(Clone)]
pub struct VariantTable {
    pub chromosome: String,
    pub samples: Vec<String>,
    pub genotypes: Vec<Vec<Genotype>>,
    pub phases: Vec<Vec<Option<VariantCallPhase>>>,
    pub genotype_likelihoods: Vec<Vec<Option<GenotypeLikelihoods>>>,
    pub variants: Vec<VcfVariant>,
    sample_to_index: HashMap<String, usize>,
}

/// Genotype likelihoods as log10 probabilities (`vcf.GenotypeLikelihoods`);
/// entries can be missing (None).
#[derive(Clone, Debug, PartialEq)]
pub struct GenotypeLikelihoods {
    pub log_prob_genotypes: Vec<Option<f64>>,
}

/// Python's `sum()` of floats (CPython 3.12 uses Neumaier summation).
fn py_float_sum(xs: &[f64]) -> f64 {
    let Some((&first, rest)) = xs.split_first() else { return 0.0 };
    let mut f = 0.0 + first;
    let mut c = 0.0;
    for &x in rest {
        let t = f + x;
        if f.abs() >= x.abs() {
            c += (f - t) + x;
        } else {
            c += (x - t) + f;
        }
        f = t;
    }
    if c != 0.0 && c.is_finite() {
        f += c;
    }
    f
}

impl GenotypeLikelihoods {
    /// `GenotypeLikelihoods.as_phred`
    pub fn as_phred(&self, ploidy: u32, regularizer: Option<f64>) -> Result<PhredGenotypeLikelihoods> {
        let mut probs = Vec::with_capacity(self.log_prob_genotypes.len());
        for p in &self.log_prob_genotypes {
            probs.push(p.ok_or_else(|| anyhow::anyhow!("TypeError: missing value in genotype likelihoods"))?);
        }
        let phred: Vec<f64> = match regularizer {
            None => {
                let m = probs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                probs.iter().map(|&p| ((p - m) * -10.0).round_ties_even()).collect()
            }
            Some(r) => {
                let p: Vec<f64> = probs.iter().map(|&x| 10f64.powf(x)).collect();
                let s = py_float_sum(&p);
                let p: Vec<f64> = p.iter().map(|&x| x / s + r).collect();
                let m = p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                p.iter().map(|&x| (-10.0 * (x / m).log10()).round_ties_even()).collect()
            }
        };
        PhredGenotypeLikelihoods::new(phred, ploidy, 2)
    }
}

impl VariantTable {
    pub fn new(chromosome: &str, samples: &[String]) -> Self {
        VariantTable {
            chromosome: chromosome.to_string(),
            samples: samples.to_vec(),
            genotypes: vec![Vec::new(); samples.len()],
            phases: vec![Vec::new(); samples.len()],
            genotype_likelihoods: vec![Vec::new(); samples.len()],
            variants: Vec::new(),
            sample_to_index: samples.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.variants.len()
    }

    pub fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    pub fn add_variant(
        &mut self,
        variant: VcfVariant,
        genotypes: Vec<Genotype>,
        phases: Vec<Option<VariantCallPhase>>,
        genotype_likelihoods: Vec<Option<GenotypeLikelihoods>>,
    ) {
        self.variants.push(variant);
        for (i, gl) in genotype_likelihoods.into_iter().enumerate() {
            self.genotype_likelihoods[i].push(gl);
        }
        for (i, g) in genotypes.into_iter().enumerate() {
            self.genotypes[i].push(g);
        }
        for (i, p) in phases.into_iter().enumerate() {
            self.phases[i].push(p);
        }
    }

    pub fn genotypes_of(&self, sample: &str) -> &[Genotype] {
        &self.genotypes[self.sample_to_index[sample]]
    }

    pub fn phases_of(&self, sample: &str) -> &[Option<VariantCallPhase>] {
        &self.phases[self.sample_to_index[sample]]
    }

    pub fn genotype_likelihoods_of(&self, sample: &str) -> &[Option<GenotypeLikelihoods>] {
        &self.genotype_likelihoods[self.sample_to_index[sample]]
    }

    /// `VariantTable.phased_blocks_as_reads`: one read per haplotype of each
    /// phased block of `sample`, restricted to `input_variants`.
    pub fn phased_blocks_as_reads(
        &self,
        sample: &str,
        input_variants: &[VcfVariant],
        source_id: i32,
        numeric_sample_id: i32,
    ) -> Result<Vec<Read>> {
        const DEFAULT_QUALITY: i32 = 20;
        const MAPQ: i32 = 100;
        const TARGET_PLOIDY: usize = 2;
        let Some(&sample_index) = self.sample_to_index.get(sample) else { return Ok(Vec::new()) };
        let input_variant_set: HashSet<&VcfVariant> = input_variants.iter().collect();
        // block_id -> reads, in insertion order
        let mut blocks: Vec<(Option<BlockId>, Vec<Read>)> = Vec::new();
        for ((variant, genotype), phase) in
            self.variants.iter().zip(&self.genotypes[sample_index]).zip(&self.phases[sample_index])
        {
            if genotype.as_vector().len() != TARGET_PLOIDY
                || !input_variant_set.contains(variant)
                || genotype.is_homozygous()
            {
                continue;
            }
            let Some(phase) = phase else { continue };
            if phase.phase[0].is_none() {
                continue;
            }
            // pysam returns PQ as float for Float fields; Cython truncates it to int
            let quality = phase.quality.map(|q| q as i32).unwrap_or(DEFAULT_QUALITY);
            let alleles: Vec<i32> = phase
                .phase
                .iter()
                .map(|a| {
                    a.map(|a| a as i32).ok_or_else(|| anyhow::anyhow!("TypeError: missing allele in phased genotype"))
                })
                .collect::<Result<_>>()?;
            match blocks.iter_mut().find(|(b, _)| *b == phase.block_id) {
                Some((_, reads)) => {
                    for (i, &allele) in alleles.iter().enumerate() {
                        reads[i].add_variant(variant.position, allele, quality);
                    }
                }
                None => {
                    let block_name = phase.block_id.as_ref().map(|b| b.to_string()).unwrap_or_else(|| "None".into());
                    let mut reads = Vec::new();
                    for (i, &allele) in alleles.iter().enumerate() {
                        let mut r = Read::new(
                            &format!("{sample}_phase_{i}_block_{block_name}"),
                            MAPQ,
                            source_id,
                            numeric_sample_id,
                        );
                        r.add_variant(variant.position, allele, quality);
                        reads.push(r);
                    }
                    blocks.push((phase.block_id.clone(), reads));
                }
            }
        }
        let mut out = Vec::new();
        for (_, reads) in blocks {
            for mut r in reads {
                if r.len() > 1 {
                    r.sort()?;
                    out.push(r);
                }
            }
        }
        Ok(out)
    }

    pub fn remove_rows_by_index(&mut self, indices: &HashSet<usize>) {
        let keep: Vec<bool> = (0..self.variants.len()).map(|i| !indices.contains(&i)).collect();
        let mut k = keep.iter();
        self.variants.retain(|_| *k.next().unwrap());
        for g in self.genotypes.iter_mut() {
            let mut k = keep.iter();
            g.retain(|_| *k.next().unwrap());
        }
        for p in self.phases.iter_mut() {
            let mut k = keep.iter();
            p.retain(|_| *k.next().unwrap());
        }
        for g in self.genotype_likelihoods.iter_mut() {
            let mut k = keep.iter();
            g.retain(|_| *k.next().unwrap());
        }
    }

    pub fn subset_rows_by_position(&mut self, positions: &HashSet<i32>) {
        let discard: HashSet<usize> = self
            .variants
            .iter()
            .enumerate()
            .filter(|(_, v)| !positions.contains(&v.position))
            .map(|(i, _)| i)
            .collect();
        self.remove_rows_by_index(&discard);
    }
}

/// `vcf.genotype_code`
fn genotype_code(gt: &Option<Vec<Option<u32>>>) -> Result<Genotype> {
    match gt {
        None => Ok(Genotype::none()),
        Some(alleles) => {
            if alleles.iter().any(|a| a.is_none()) {
                Ok(Genotype::none())
            } else {
                let a: Vec<u32> = alleles.iter().map(|a| a.unwrap()).collect();
                Genotype::new(&a)
            }
        }
    }
}

/// Parses records into a `VariantTable` (`VcfReader._process_single_chromosome`).
struct TableBuilder<'a> {
    only_snvs: bool,
    phases: bool,
    genotype_likelihoods: bool,
    mav: bool,
    ignore_genotypes: bool,
    ploidy: &'a mut Option<usize>,
    warned_duplicate: &'a mut bool,
    table: VariantTable,
    prev_position: Option<i64>,
    phase_detected: Option<&'static str>,
    n_snvs: usize,
    n_other: usize,
    n_multi: usize,
    /// Haploid calls left out of phasing.
    n_haploid: usize,
}

impl TableBuilder<'_> {
    fn process(&mut self, header: &Header, record: &Record) -> Result<()> {
        let chromosome = self.table.chromosome.clone();
        let Some(alts) = record.alts()? else { return Ok(()) };
        if alts.len() > 1 {
            self.n_multi += 1;
            if !self.mav || alts.len() >= crate::core::genotype::MAX_ALLELES as usize {
                return Ok(());
            }
        }
        let pos = record.pos();
        let reference = record.reference()?.unwrap_or_default();
        if reference.len() == 1 && alts.iter().all(|a| a.len() == 1) {
            self.n_snvs += 1;
        } else {
            self.n_other += 1;
            if self.only_snvs {
                return Ok(());
            }
        }
        if let Some(prev) = self.prev_position {
            if prev > pos {
                bail!("VCF not ordered: {}:{} appears before {}:{}", chromosome, prev + 1, chromosome, pos + 1);
            }
            if prev == pos {
                if !*self.warned_duplicate {
                    log::warn!(
                        "Skipping duplicated position {} on chromosome '{}' Hiding further warnings of this type, use --debug to show",
                        pos + 1,
                        chromosome
                    );
                    *self.warned_duplicate = true;
                } else {
                    log::debug!("Skipping duplicated position {} on chromosome '{}'", pos + 1, chromosome);
                }
                return Ok(());
            }
        }
        self.prev_position = Some(pos);
        let n_samples = self.table.samples.len();
        let mut phases = vec![None; n_samples];
        if self.phases {
            for (s, phase_slot) in phases.iter_mut().enumerate() {
                let mut phase = None;
                for (p, name) in
                    [(extract_hp_phase(header, record, s)?, "HP"), (extract_gt_ps_phase(header, record, s)?, "GT_PS")]
                {
                    let Some(p) = p else { continue };
                    match self.phase_detected {
                        None => self.phase_detected = Some(name),
                        Some(d) if d != name => {
                            bail!("Mixed phasing information in input VCF (e.g. mixing PS and HP fields)")
                        }
                        _ => {}
                    }
                    let phase_ploidy = p.phase.len();
                    if phase_ploidy > 15 {
                        return Err(PloidyError("Ploidies higher than 15 are not supported.".into()).into());
                    } else if p.block_id.is_none() {
                    } else if self.ploidy.is_none() {
                        *self.ploidy = Some(phase_ploidy);
                    } else if Some(phase_ploidy) != *self.ploidy {
                        return Err(PloidyError(format!(
                            "Phasing information contains inconsistent ploidy ({} and {})",
                            self.ploidy.unwrap(),
                            phase_ploidy
                        ))
                        .into());
                    }
                    phase = Some(p);
                }
                *phase_slot = phase;
            }
        }
        let mut gls = vec![None; n_samples];
        if self.genotype_likelihoods {
            for (s, slot) in gls.iter_mut().enumerate() {
                // Prefer GLs (floats) over PLs (ints) if both should be present
                if let Some(gl) = record.format_numbers(header, s, "GL")? {
                    *slot = Some(GenotypeLikelihoods { log_prob_genotypes: gl });
                } else if let Some(pl) = record.format_numbers(header, s, "PL")? {
                    *slot = Some(GenotypeLikelihoods {
                        log_prob_genotypes: pl.iter().map(|p| p.map(|p| p / -10.0)).collect(),
                    });
                }
            }
        }
        if self.ignore_genotypes {
            let variant = VcfVariant { position: pos as i32, reference_allele: reference, alternative_alleles: alts };
            self.table.add_variant(variant, vec![Genotype::none(); n_samples], vec![None; n_samples], gls);
            return Ok(());
        }
        let mut gts = Vec::with_capacity(n_samples);
        for s in 0..n_samples {
            gts.push(record.gt(header, s)?);
        }
        // A haploid call (chrX or chrY outside the PARs of an XY sample) has nothing to phase:
        // it is left out and written as it is, beside the chromosome's diploid calls. WhatsHap
        // stops there with "Inconsistent ploidy (2 and 1)" (LRP-003). It is left out per
        // sample: another sample's diploid call at the same site is phased (LRP-017).
        let haploid = |g: &Option<Vec<Option<u32>>>| g.as_ref().is_some_and(|g| g.len() == 1 && g[0].is_some());
        if gts.iter().any(haploid) {
            self.n_haploid += 1;
            let diploid =
                |g: &Option<Vec<Option<u32>>>| g.as_ref().is_some_and(|g| g.len() > 1 && g.iter().all(Option::is_some));
            if !gts.iter().any(diploid) {
                return Ok(());
            }
            for (g, gl) in gts.iter_mut().zip(gls.iter_mut()) {
                if haploid(g) {
                    *g = None;
                    *gl = None;
                }
            }
        }
        for g in gts.iter().flatten() {
            if g.iter().any(|a| a.is_none()) {
                continue;
            }
            let ploidy = g.len();
            if ploidy > 15 {
                return Err(PloidyError("Ploidies higher than 15 are not supported.".into()).into());
            }
            match *self.ploidy {
                None => *self.ploidy = Some(ploidy),
                Some(p) if p != ploidy => {
                    return Err(PloidyError(format!("Inconsistent ploidy ({p} and {ploidy})")).into());
                }
                _ => {}
            }
        }
        let genotypes = gts.iter().map(genotype_code).collect::<Result<Vec<_>>>()?;
        let variant = VcfVariant { position: pos as i32, reference_allele: reference, alternative_alleles: alts };
        self.table.add_variant(variant, genotypes, phases, gls);
        Ok(())
    }

    fn finish(self) -> VariantTable {
        log::debug!(
            "Parsed {} SNVs and {} non-SNVs. Also found {} multi-ALTs.",
            self.n_snvs,
            self.n_other,
            self.n_multi
        );
        if self.n_haploid > 0 {
            log::info!("Left {} haploid calls unphased", self.n_haploid);
        }
        self.table
    }
}

/// `VcfReader._extract_HP_phase`
fn extract_hp_phase(header: &Header, record: &Record, sample: usize) -> Result<Option<VariantCallPhase>> {
    if !header.is_defined(BCF_HL_FMT, "HP") {
        return Ok(None);
    }
    let Some(hp) = record.format_strings(header, sample, "HP")? else { return Ok(None) };
    if hp.is_empty() || hp == [Some(".".to_string())] || hp.iter().any(|v| v.is_none()) {
        if hp.is_empty() || hp == [Some(".".to_string())] {
            return Ok(None);
        }
        bail!("invalid HP value");
    }
    let mut fields: Vec<Vec<i64>> = Vec::new();
    for v in hp.iter().flatten() {
        let parts = v.split('-').map(|x| py_int(x).and_then(|x| i64::try_from(x).ok())).collect::<Option<Vec<_>>>();
        fields.push(parts.ok_or_else(|| anyhow::anyhow!("invalid literal for int() in HP tag '{v}'"))?);
    }
    if fields.iter().any(|f| f[0] != fields[0][0]) {
        bail!(
            "AssertionError: HP tag values of different blocks: {}",
            hp.iter().flatten().cloned().collect::<Vec<_>>().join(",")
        );
    }
    let block_id = fields[0][0];
    let order = fields
        .iter()
        .map(|f| f.get(1).map(|o| o - 1))
        .collect::<Option<Vec<i64>>>()
        .ok_or_else(|| anyhow::anyhow!("IndexError in HP tag: list index out of range"))?;
    let gt = record.gt(header, sample)?.unwrap_or_default();
    let mut phase = Vec::new();
    for i in 0..order.len() {
        let k = order.iter().position(|&o| o == i as i64).ok_or_else(|| anyhow::anyhow!("ValueError in HP tag"))?;
        phase.push(*gt.get(k).ok_or_else(|| anyhow::anyhow!("IndexError in HP tag"))?);
    }
    Ok(Some(VariantCallPhase {
        block_id: Some(BlockId::Int(block_id)),
        phase,
        quality: phase_quality(header, record, sample)?,
    }))
}

/// `VcfReader._extract_GT_PS_phase`
fn extract_gt_ps_phase(header: &Header, record: &Record, sample: usize) -> Result<Option<VariantCallPhase>> {
    if !record.sample_phased(header, sample)? {
        return Ok(None);
    }
    let gt = record.gt(header, sample)?.unwrap_or_default();
    let is_het = !gt.iter().all(|a| *a == gt[0]);
    if !is_het {
        return Ok(None);
    }
    let block_id = if !record.has_format(header, "PS")? {
        Some(BlockId::Int(0))
    } else if header.id2type(BCF_HL_FMT, header.id("PS")) == hts::BCF_HT_STR {
        match record.format_strings(header, sample, "PS")? {
            Some(v) if v.len() == 1 && v[0].is_some() => Some(BlockId::Str(v[0].clone().unwrap())),
            Some(v) if v.len() == 1 => None,
            _ => bail!("unsupported PS value"),
        }
    } else {
        match record.format_scalar(header, sample, "PS")? {
            None => Some(BlockId::Int(0)),
            Some(None) => None,
            Some(Some(hts::FormatScalar::Int(v))) => Some(BlockId::Int(v)),
            Some(Some(hts::FormatScalar::Float(v))) => Some(BlockId::Float(v as f64)),
        }
    };
    Ok(Some(VariantCallPhase { block_id, phase: gt, quality: phase_quality(header, record, sample)? }))
}

/// `call.get("PQ", None)`
fn phase_quality(header: &Header, record: &Record, sample: usize) -> Result<Option<f64>> {
    Ok(match record.format_scalar(header, sample, "PQ")? {
        Some(Some(hts::FormatScalar::Int(v))) => Some(v as f64),
        Some(Some(hts::FormatScalar::Float(v))) => Some(v as f64),
        _ => None,
    })
}

/// Reads a VCF chromosome by chromosome (`VcfReader.__iter__`).
pub struct VcfReader {
    file: VcfFile,
    path: String,
    only_snvs: bool,
    phases: bool,
    genotype_likelihoods: bool,
    mav: bool,
    ignore_genotypes: bool,
    pub samples: Vec<String>,
    ploidy: Option<usize>,
    pending: Option<Record>,
    eof: bool,
    warned_duplicate: bool,
}

impl VcfReader {
    pub fn open(path: &str, only_snvs: bool) -> Result<Self> {
        Self::open_with(path, only_snvs, false, None)
    }

    pub fn open_with(path: &str, only_snvs: bool, phases: bool, ploidy: Option<usize>) -> Result<Self> {
        let file = VcfFile::open(path)?;
        let samples = file.header.samples();
        log::debug!("Found {} sample(s) in the VCF file.", samples.len());
        Ok(VcfReader {
            file,
            path: path.to_string(),
            only_snvs,
            phases,
            genotype_likelihoods: false,
            mav: false,
            ignore_genotypes: false,
            samples,
            ploidy,
            pending: None,
            eof: false,
            warned_duplicate: false,
        })
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    pub fn header(&self) -> &Header {
        &self.file.header
    }

    /// Also read genotype likelihoods (GL, or PL if there is no GL).
    /// Do not read genotypes (all are none) or phases.
    pub fn with_ignore_genotypes(mut self, ignore_genotypes: bool) -> Self {
        self.ignore_genotypes = ignore_genotypes;
        self
    }

    /// Also read multi-ALT sites (with fewer than 16 alleles).
    pub fn with_mav(mut self, mav: bool) -> Self {
        self.mav = mav;
        self
    }

    pub fn with_genotype_likelihoods(mut self, genotype_likelihoods: bool) -> Self {
        self.genotype_likelihoods = genotype_likelihoods;
        self
    }

    /// Next chromosome's variant table, or None at the end.
    pub fn next_table(&mut self) -> Result<Option<VariantTable>> {
        if self.pending.is_none() && !self.eof {
            let mut rec = Record::new();
            if self.file.read(&mut rec)? {
                self.pending = Some(rec);
            } else {
                self.eof = true;
            }
        }
        let Some(first) = self.pending.take() else { return Ok(None) };
        let chromosome = first.chrom(&self.file.header);
        let mut b = TableBuilder {
            only_snvs: self.only_snvs,
            phases: self.phases,
            genotype_likelihoods: self.genotype_likelihoods,
            mav: self.mav,
            ignore_genotypes: self.ignore_genotypes,
            ploidy: &mut self.ploidy,
            warned_duplicate: &mut self.warned_duplicate,
            table: VariantTable::new(&chromosome, &self.samples),
            prev_position: None,
            phase_detected: None,
            n_snvs: 0,
            n_other: 0,
            n_multi: 0,
            n_haploid: 0,
        };
        b.process(&self.file.header, &first)?;
        let mut next = Record::new();
        loop {
            if !self.file.read(&mut next)? {
                self.eof = true;
                break;
            }
            if next.chrom(&self.file.header) != chromosome {
                self.pending = Some(next);
                break;
            }
            b.process(&self.file.header, &next)?;
        }
        Ok(Some(b.finish()))
    }

    /// `VcfReader.fetch_regions`: records of `chromosome` overlapping the
    /// regions, concatenated in the given order.
    pub fn fetch_regions(
        &mut self,
        indexed: &mut IndexedVcf,
        chromosome: &str,
        regions: &[(i64, Option<i64>)],
    ) -> Result<VariantTable, VcfFetchError> {
        let samples = self.samples.clone();
        fetch_variant_table(
            indexed,
            &samples,
            self.only_snvs,
            self.phases,
            &mut self.ploidy,
            &mut self.warned_duplicate,
            chromosome,
            regions,
        )
    }
}

/// `VcfReader.fetch_regions` without a `VcfReader` (so that it can be used
/// from several threads, each with its own `IndexedVcf`).
#[allow(clippy::too_many_arguments)]
pub fn fetch_variant_table(
    indexed: &mut IndexedVcf,
    samples: &[String],
    only_snvs: bool,
    phases: bool,
    ploidy: &mut Option<usize>,
    warned_duplicate: &mut bool,
    chromosome: &str,
    regions: &[(i64, Option<i64>)],
) -> Result<VariantTable, VcfFetchError> {
    let mut records = Vec::new();
    for &(start, end) in regions {
        records.extend(indexed.fetch(chromosome, start, end)?);
    }
    let header = indexed.header();
    let mut b = TableBuilder {
        only_snvs,
        phases,
        genotype_likelihoods: false,
        mav: false,
        ignore_genotypes: false,
        ploidy,
        warned_duplicate,
        table: VariantTable::new(chromosome, samples),
        prev_position: None,
        phase_detected: None,
        n_snvs: 0,
        n_other: 0,
        n_multi: 0,
        n_haploid: 0,
    };
    for r in &records {
        b.process(header, r).map_err(VcfFetchError::Other)?;
    }
    Ok(b.finish())
}

#[derive(Clone)]
struct PredefinedHeader {
    format_or_info: &'static str,
    id: &'static str,
    number: &'static str,
    typ: &'static str,
    description: &'static str,
}

impl PredefinedHeader {
    fn line(&self) -> String {
        format!(
            "##{}=<ID={},Number={},Type={},Description=\"{}\">",
            self.format_or_info, self.id, self.number, self.typ, self.description
        )
    }
}

fn predefined_format(id: &str) -> Option<PredefinedHeader> {
    let h = |id, number, typ, description| PredefinedHeader { format_or_info: "FORMAT", id, number, typ, description };
    Some(match id {
        "GL" => h(
            "GL",
            "G",
            "Float",
            "Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy",
        ),
        "GQ" => h("GQ", "1", "Integer", "Phred-scaled genotype quality"),
        "GT" => h("GT", "1", "String", "Genotype"),
        "HP" => h("HP", ".", "String", "Phasing haplotype identifier"),
        "PQ" => h("PQ", "1", "Float", "Phasing quality"),
        "PS" => h("PS", "1", "Integer", "Phase set identifier"),
        "HS" => h("HS", ".", "Integer", "Haploid phase set identifier"),
        "AD" => h("AD", ".", "Integer", "Observed allele depths"),
        _ => return None,
    })
}

fn predefined_info(id: &str) -> Option<PredefinedHeader> {
    let h = |id, number, typ, description| PredefinedHeader { format_or_info: "INFO", id, number, typ, description };
    Some(match id {
        "AC" => h("AC", "A", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed"),
        "AN" => h("AN", "A", "Integer", "Total number of alleles in called genotypes"),
        "END" => h("END", "1", "Integer", "Stop position of the interval"),
        "SVLEN" => h("SVLEN", ".", "Integer", "Difference in length between REF and ALT alleles"),
        "SVTYPE" => h("SVTYPE", "1", "String", "Type of structural variant"),
        _ => return None,
    })
}

/// `vcf.missing_headers`: (missing contigs, formats to (re)define, missing infos)
fn missing_headers(path: &str) -> Result<(Vec<String>, Vec<String>, Vec<String>)> {
    let mut file = VcfFile::open(path)?;
    let header = file.header.dup();
    let mut incorrect_formats = Vec::new();
    for (fmt, v) in header.defined(BCF_HL_FMT) {
        let Some(h) = predefined_format(&fmt) else { continue };
        if v.number != h.number || (v.typ != h.typ && !(v.typ == "Float" && h.typ == "Integer")) {
            if fmt == "PS" && v.typ != h.typ {
                bail!(
                    "The input VCF/BCF contains phase set ('PS') tags that are of the non-standard type '{}' instead of 'Integer'. WhatsHap cannot overwrite these as it could produce inconsistent files. To proceed, you can use 'whatshap unphase' to remove phasing information from the input file",
                    v.typ
                );
            }
            incorrect_formats.push(fmt);
        }
    }
    let mut contigs: Vec<String> = Vec::new();
    let mut contig_set: HashSet<String> = HashSet::new();
    let mut formats: Vec<String> = Vec::new();
    let mut format_set: HashSet<String> = HashSet::new();
    let mut seen_infos: Vec<String> = Vec::new();
    let mut seen_info_set: HashSet<String> = HashSet::new();
    let mut rec = Record::new();
    while file.read(&mut rec)? {
        for k in rec.info_keys(&file.header)? {
            if seen_info_set.insert(k.clone()) {
                seen_infos.push(k);
            }
        }
        if rec.alts()?.unwrap_or_default().iter().any(|a| a.starts_with('<')) && seen_info_set.insert("END".into()) {
            seen_infos.push("END".into());
        }
        let c = rec.chrom(&file.header);
        if contig_set.insert(c.clone()) {
            contigs.push(c);
        }
        for f in rec.format_keys(&file.header)? {
            if format_set.insert(f.clone()) {
                formats.push(f);
            }
        }
    }
    let header_contigs: HashSet<String> = header.contigs().into_iter().collect();
    let header_formats: HashSet<String> = header.defined(BCF_HL_FMT).into_iter().map(|(k, _)| k).collect();
    let header_infos: HashSet<String> = header.defined(BCF_HL_INFO).into_iter().map(|(k, _)| k).collect();
    let missing_contigs = contigs.into_iter().filter(|c| !header_contigs.contains(c)).collect();
    let missing_formats: Vec<String> = formats.into_iter().filter(|f| !header_formats.contains(f)).collect();
    // `list(set(seen_infos) - set(header.info))`: order follows Python's
    // string-set iteration order, which is randomized per process; we use
    // first-seen order.
    let missing_infos = seen_infos.into_iter().filter(|i| !header_infos.contains(i)).collect();
    incorrect_formats.extend(missing_formats);
    Ok((missing_contigs, incorrect_formats, missing_infos))
}

/// `vcf.augment_header`
fn augment_header(header: &mut Header, contigs: &[String], formats: &[String], infos: &[String]) -> Result<()> {
    for contig in contigs {
        header.add_meta_items("contig", &[("ID", contig, false)])?;
    }
    for fmt in formats {
        if header.is_defined(BCF_HL_FMT, fmt) {
            header.remove(BCF_HL_FMT, fmt);
        }
        let Some(h) = predefined_format(fmt) else { bail!("FORMAT '{fmt}' not defined in VCF header") };
        header.add_line(&h.line())?;
    }
    for info in infos {
        let Some(h) = predefined_info(info) else { bail!("INFO '{info}' not defined in VCF header") };
        header.add_line(&h.line())?;
    }
    Ok(())
}

/// A genotype changed while writing (with --distrust-genotypes).
#[derive(Clone, Debug)]
pub struct GenotypeChange {
    pub sample: String,
    pub chromosome: String,
    /// 0-based position
    pub position: i64,
    pub reference: String,
    pub alternative: String,
    pub old_gt: Genotype,
    pub new_gt: Genotype,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhaseTag {
    PS,
    HP,
}

/// Opens the input VCF with the header of the phased output VCF: missing
/// definitions added, the command line, and the phase tag definition
/// (`VcfAugmenter.__init__` + `PhasedVcfWriter.setup_header`).
/// `VcfAugmenter.__init__` up to `setup_header`: the input VCF with missing
/// header lines added and the command line recorded
fn open_augmented(in_path: &str, command_line: Option<&str>) -> Result<VcfFile> {
    log::debug!("Reading the input VCF to find possibly missing headers");
    let (contigs, formats, infos) = missing_headers(in_path)?;
    let mut reader = VcfFile::open(in_path)?;
    augment_header(&mut reader.header, &contigs, &formats, &infos)?;
    if let Some(cl) = command_line {
        let quoted = format!("\"{}\"", cl.replace('"', ""));
        reader.header.add_meta_value("commandline", &quoted)?;
    }
    Ok(reader)
}

/// `VcfAugmenter`: reads the input VCF chromosome by chromosome and writes
/// each record after it was modified
struct VcfAugmenter {
    reader: VcfFile,
    writer: Option<VcfWriter>,
    unprocessed: Option<Record>,
    reader_eof: bool,
}

impl VcfAugmenter {
    fn new(reader: VcfFile, out_path: &str, plain: bool) -> Result<Self> {
        let writer = if plain {
            VcfWriter::create_plain(out_path, &reader.header)?
        } else {
            VcfWriter::create(out_path, &reader.header)?
        };
        Ok(VcfAugmenter { reader, writer: Some(writer), unprocessed: None, reader_eof: false })
    }

    /// `_iterrecords`: the next record of the chromosome, which must come
    /// next; `n` counts the records of this pass (0 at its start). As in
    /// WhatsHap, the saved first record of the next chromosome is kept after
    /// it was returned, so a later pass for another chromosome fails once the
    /// input is exhausted.
    fn next_record(&mut self, chromosome: &str, n: &mut usize) -> Result<Option<Record>> {
        if *n == 0
            && let Some(r) = &self.unprocessed
        {
            if r.chrom(&self.reader.header) != chromosome {
                bail!("AssertionError: VCF records of chromosome {chromosome} expected");
            }
            *n += 1;
            return Ok(Some(r.clone()));
        }
        if self.reader_eof {
            return Ok(None);
        }
        let mut rec = Record::new();
        if !self.reader.read(&mut rec)? {
            self.reader_eof = true;
            return Ok(None);
        }
        *n += 1;
        if rec.chrom(&self.reader.header) != chromosome {
            if *n == 1 {
                bail!("AssertionError: chromosome {chromosome} is not next in the VCF");
            }
            self.unprocessed = Some(rec);
            return Ok(None);
        }
        Ok(Some(rec))
    }

    fn write_record(&mut self, rec: &mut Record) -> Result<()> {
        self.writer.as_mut().unwrap().write(rec, &self.reader.header)
    }

    fn write_unchanged(&mut self, chromosome: &str) -> Result<()> {
        let mut n = 0;
        while let Some(mut r) = self.next_record(chromosome, &mut n)? {
            self.write_record(&mut r)?;
        }
        Ok(())
    }

    fn close(mut self) -> Result<()> {
        self.writer.take().unwrap().close()
    }
}

pub fn open_for_phased_output(in_path: &str, command_line: Option<&str>, tag: PhaseTag) -> Result<VcfFile> {
    let mut reader = open_augmented(in_path, command_line)?;
    // setup_header: remove FreeBayes' "phasing=none" and add the tag definition
    reader.header.remove_first_record_with_key("phasing");
    let tag_name = match tag {
        PhaseTag::PS => "PS",
        PhaseTag::HP => "HP",
    };
    reader.header.add_line(&predefined_format(tag_name).unwrap().line())?;
    Ok(reader)
}

/// Reads the input VCF and writes it back out with phasing information
/// (`VcfAugmenter` + `PhasedVcfWriter`).
pub struct PhasedVcfWriter {
    aug: VcfAugmenter,
    tag: PhaseTag,
    only_snvs: bool,
    /// also phase multi-ALT sites, with any allele
    mav: bool,
    phase_tag_found_warned: bool,
    pub samples: Vec<String>,
}

impl PhasedVcfWriter {
    pub fn new(
        in_path: &str,
        command_line: Option<&str>,
        out_path: &str,
        tag: PhaseTag,
        only_snvs: bool,
    ) -> Result<Self> {
        Self::new_as(in_path, command_line, out_path, tag, only_snvs, false)
    }

    /// Like `new`, with the output opened as a text file (uncompressed VCF
    /// whatever the file name), as `hapcut2vcf` does.
    pub fn new_plain_output(in_path: &str, command_line: Option<&str>, out_path: &str) -> Result<Self> {
        Self::new_as(in_path, command_line, out_path, PhaseTag::PS, false, true)
    }

    fn new_as(
        in_path: &str,
        command_line: Option<&str>,
        out_path: &str,
        tag: PhaseTag,
        only_snvs: bool,
        plain: bool,
    ) -> Result<Self> {
        let reader = open_for_phased_output(in_path, command_line, tag)?;
        let samples = reader.header.samples();
        Ok(PhasedVcfWriter {
            aug: VcfAugmenter::new(reader, out_path, plain)?,
            tag,
            only_snvs,
            mav: false,
            phase_tag_found_warned: false,
            samples,
        })
    }

    /// Also phase multi-ALT sites and alleles other than 0 and 1 (`mav`).
    pub fn with_mav(mut self, mav: bool) -> Self {
        self.mav = mav;
        self
    }

    fn next_record(&mut self, chromosome: &str, n: &mut usize) -> Result<Option<Record>> {
        self.aug.next_record(chromosome, n)
    }

    fn write_record(&mut self, rec: &mut Record) -> Result<()> {
        self.aug.write_record(rec)
    }

    /// Add phasing information to all variants of one chromosome.
    /// `sample_superreads` maps a sample to its two super reads,
    /// `sample_components` maps a sample to position -> component.
    pub fn write(
        &mut self,
        chromosome: &str,
        sample_superreads: &[(String, ReadSet)],
        sample_components: &HashMap<String, HashMap<i32, i32>>,
    ) -> Result<Vec<GenotypeChange>> {
        let mut sample_phases: HashMap<&str, HashMap<i32, Vec<u32>>> = HashMap::new();
        for (sample, superreads) in sample_superreads {
            let mut phases = HashMap::new();
            let (r0, r1) = (superreads.get(0), superreads.get(1));
            for (v0, v1) in r0.variants.iter().zip(r1.variants.iter()) {
                let ok = [v0.allele, v1.allele].iter().all(|&a| a == 0 || a == 1 || self.mav);
                if ok {
                    phases.insert(v0.position, vec![v0.allele as u32, v1.allele as u32]);
                }
            }
            sample_phases.insert(sample.as_str(), phases);
        }
        let mut changes = Vec::new();
        let phased_samples: Vec<&str> = sample_superreads.iter().map(|(s, _)| s.as_str()).collect();
        let sample_index: HashMap<String, usize> =
            self.samples.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect();
        let mut prev_pos: Option<i64> = None;
        let mut n = 0;
        while let Some(mut record) = self.next_record(chromosome, &mut n)? {
            self.modify_record(
                &mut record,
                chromosome,
                &phased_samples,
                &sample_index,
                &sample_phases,
                sample_components,
                &mut prev_pos,
                &mut changes,
            )?;
            self.write_record(&mut record)?;
        }
        Ok(changes)
    }

    #[allow(clippy::too_many_arguments)]
    fn modify_record(
        &mut self,
        record: &mut Record,
        chromosome: &str,
        phased_samples: &[&str],
        sample_index: &HashMap<String, usize>,
        sample_phases: &HashMap<&str, HashMap<i32, Vec<u32>>>,
        sample_components: &HashMap<String, HashMap<i32, i32>>,
        prev_pos: &mut Option<i64>,
        changes: &mut Vec<GenotypeChange>,
    ) -> Result<()> {
        let header = &self.aug.reader.header;
        // _remove_existing_phasing
        if matches!(self.tag, PhaseTag::PS) {
            for s in phased_samples {
                let si = sample_index[*s];
                if !record.has_format(header, "GT")? {
                    continue;
                }
                record.set_phased(header, si, false)?;
                if let Some(gt) = record.gt(header, si)?
                    && gt.iter().all(|a| a.is_some())
                {
                    let mut a: Vec<u32> = gt.iter().map(|a| a.unwrap()).collect();
                    a.sort_unstable();
                    record.set_gt(header, si, &a)?;
                }
            }
        }
        let pos = record.pos();
        let Some(alts) = record.alts()? else { return Ok(()) };
        if alts.len() > 1 && !self.mav {
            // multi-ALT sites are not phased unless requested
            return Ok(());
        }
        if *prev_pos == Some(pos) {
            return Ok(());
        }
        let reference = record.reference()?.unwrap_or_default();
        let is_snv = reference.len() == 1 && alts[0].len() == 1;
        if self.only_snvs && !is_snv {
            return Ok(());
        }
        let pos32 = pos as i32;
        // Is the variant phased in any sample?
        let mut phased_any = false;
        for sample in &self.samples {
            if let Some(phases) = sample_phases.get(sample.as_str())
                && sample_components[sample].contains_key(&pos32)
                && phases.contains_key(&pos32)
            {
                phased_any = true;
                break;
            }
        }
        if !phased_any {
            return Ok(());
        }
        let tag = match self.tag {
            PhaseTag::PS => "PS",
            PhaseTag::HP => "HP",
        };
        for s in phased_samples {
            let si = sample_index[*s];
            let components = &sample_components[*s];
            let phases = &sample_phases[s];
            if !self.phase_tag_found_warned
                && record.has_format(header, tag)?
                && self.tag_value_present(record, si, tag)?
            {
                log::warn!("Ignoring existing phasing information found in input VCF ({tag} tag exists).");
                self.phase_tag_found_warned = true;
            }
            let gt = record.gt(header, si)?;
            let gt_type = genotype_code(&gt)?;
            let mut is_het = !gt_type.is_homozygous();
            // is the genotype to be changed?
            if let Some(phase) = phases.get(&pos32) {
                let new_gt = Genotype::new(phase)?;
                if new_gt != gt_type {
                    record.set_gt(header, si, &new_gt.as_vector())?;
                    changes.push(GenotypeChange {
                        sample: s.to_string(),
                        chromosome: chromosome.to_string(),
                        position: pos,
                        reference: reference.clone(),
                        alternative: alts[0].clone(),
                        old_gt: gt_type,
                        new_gt,
                    });
                    is_het = !new_gt.is_homozygous();
                }
            }
            match (components.get(&pos32), phases.get(&pos32), is_het, self.tag) {
                (Some(&component), Some(phase), true, PhaseTag::PS) => {
                    record.format_set_int(header, si, "PS", &[Some(component as i64 + 1)])?;
                    record.set_gt(header, si, phase)?;
                    record.set_phased(header, si, true)?;
                }
                (Some(&component), Some(phase), true, PhaseTag::HP) => {
                    let value: Vec<String> = phase.iter().map(|a| format!("{}-{}", component + 1, a + 1)).collect();
                    record.format_set_str(header, si, "HP", Some(&value.join(",")))?;
                }
                (_, _, _, PhaseTag::PS) => record.format_set_int(header, si, "PS", &[None])?,
                (_, _, _, PhaseTag::HP) => record.format_set_str(header, si, "HP", None)?,
            }
        }
        *prev_pos = Some(pos);
        Ok(())
    }

    /// pysam `call[tag] is not None` (only used for a warning).
    fn tag_value_present(&self, record: &Record, sample: usize, tag: &str) -> Result<bool> {
        record.format_value_present(&self.aug.reader.header, sample, tag)
    }

    /// Write all records of a chromosome unchanged.
    pub fn write_unchanged(&mut self, chromosome: &str) -> Result<()> {
        self.aug.write_unchanged(chromosome)
    }

    pub fn close(self) -> Result<()> {
        self.aug.close()
    }
}

/// Genotypes and likelihoods (0/0, 0/1, 1/1) of one sample, per variant of
/// a variant table; `None`: not genotyped
#[derive(Clone)]
pub struct SampleGenotypes {
    pub genotypes: Vec<Genotype>,
    pub likelihoods: Vec<Option<[f64; 3]>>,
}

/// Reads the input VCF and writes it back out with genotypes, genotype
/// likelihoods and qualities (`GenotypeVcfWriter`).
pub struct GenotypeVcfWriter {
    aug: VcfAugmenter,
}

impl GenotypeVcfWriter {
    /// `plain`: write uncompressed VCF whatever the file name (output opened
    /// as a text file)
    pub fn new(in_path: &str, command_line: Option<&str>, out_path: &str, plain: bool) -> Result<Self> {
        let mut reader = open_augmented(in_path, command_line)?;
        // setup_header
        for line in [
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype computed by WhatsHap genotyping algorithm\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled genotype quality computed by WhatsHap genotyping algorithm\">",
            "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled likelihoods for genotypes: 0/0, 0/1, 1/1, computed by WhatsHap genotyping algorithm\">",
        ] {
            reader.header.add_line(line)?;
        }
        Ok(GenotypeVcfWriter { aug: VcfAugmenter::new(reader, out_path, plain)? })
    }

    pub fn write_unchanged(&mut self, chromosome: &str) -> Result<()> {
        self.aug.write_unchanged(chromosome)
    }

    pub fn close(self) -> Result<()> {
        self.aug.close()
    }

    /// `write_genotypes` (diploid): `positions` are the variants of the
    /// table, `calls` the genotyping results by sample name
    pub fn write_genotypes(
        &mut self,
        chromosome: &str,
        positions: &[i32],
        calls: &HashMap<String, SampleGenotypes>,
    ) -> Result<()> {
        const PLOIDY: usize = 2;
        let genotyped: HashMap<i32, usize> = positions.iter().enumerate().map(|(i, &p)| (p, i)).collect();
        let samples = self.aug.reader.header.samples();
        let mut n = 0;
        while let Some(mut record) = self.aug.next_record(chromosome, &mut n)? {
            let header = &self.aug.reader.header;
            if let Some(alts) = record.alts()? {
                let pos = record.pos() as i32;
                let n_alleles = 1 + alts.len() as i32;
                let n_genotypes =
                    crate::core::genotype::binomial_coefficient(PLOIDY as i32 + n_alleles - 1, n_alleles - 1);
                for (s, sample) in samples.iter().enumerate() {
                    let mut geno = Genotype::none();
                    let mut geno_l: Vec<f64> = vec![1.0 / n_genotypes as f64; n_genotypes.max(0) as usize];
                    if let Some(&idx) = genotyped.get(&pos)
                        && let Some(sg) = calls.get(sample)
                        && let Some(likelihoods) = sg.likelihoods[idx]
                    {
                        geno_l = likelihoods.to_vec();
                        geno = sg.genotypes[idx];
                    }
                    let geno_index = geno.get_index() as usize;
                    let mut geno_q = 0.0;
                    for i in 0..n_genotypes.max(0) as usize {
                        let l = *geno_l.get(i).context("IndexError: list index out of range")?;
                        if i != geno_index {
                            geno_q += l;
                        }
                    }
                    record.format_set_value(header, s, "GT", &PyFormatValue::Alleles(vec![Some(0); PLOIDY]))?;
                    let gl: Vec<Option<f64>> =
                        geno_l.iter().map(|&j| Some(if j > 0.0 { j.log10().max(-1000.0) } else { -1000.0 })).collect();
                    record.format_set_value(header, s, "GL", &PyFormatValue::Floats(gl))?;
                    let alleles: Vec<Option<u32>> = geno.as_vector().into_iter().map(Some).collect();
                    record.format_set_value(header, s, "GT", &PyFormatValue::Alleles(alleles))?;
                    let gq = if geno.is_none() {
                        PyFormatValue::Nulls(1)
                    } else if geno_q > 0.0 {
                        PyFormatValue::Ints(vec![Some(((-10.0 * geno_q.log10()).round_ties_even() as i64).min(10000))])
                    } else {
                        PyFormatValue::Ints(vec![Some(10000)])
                    };
                    record.format_set_value(header, s, "GQ", &gq)?;
                    record.set_qual_missing();
                    // delete all other genotype information
                    for tag in record.format_keys_present(header) {
                        if !matches!(tag.as_str(), "GT" | "GL" | "GQ") {
                            let null = record.format_null_value(header, s, &tag)?;
                            record.format_set_value(header, s, &tag, &null)?;
                        }
                    }
                }
            }
            self.aug.write_record(&mut record)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod ploidy_tests {
    use super::*;

    #[test]
    fn haploid_calls_beside_diploid_ones_are_left_unphased() {
        // LRP-003: chrX of an XY sample, diploid in a PAR and haploid outside it.
        let dir = std::env::temp_dir().join(format!("lrpipe-phase-ploidy-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let path = dir.join("x.vcf");
        let header = "##fileformat=VCFv4.2\n##contig=<ID=chrX,length=1000>\n\
                      ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
        let records = [
            "chrX\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/1",
            "chrX\t200\t.\tA\tG\t30\tPASS\t.\tGT\t1",
            "chrX\t300\t.\tC\tT\t30\tPASS\t.\tGT\t0/1",
        ];
        std::fs::write(&path, format!("{header}{}\n", records.join("\n"))).expect("VCF written");
        let mut reader = VcfReader::open(path.to_str().expect("UTF-8"), false).expect("opens");
        let table = reader.next_table().expect("no ploidy error").expect("chrX");
        let positions: Vec<i32> = table.variants.iter().map(|v| v.position).collect();
        assert_eq!(positions, [99, 299], "the haploid call is left out");
    }

    #[test]
    fn a_haploid_sample_leaves_another_samples_calls_in() {
        // LRP-017: F diploid, M haploid at the same chrX sites; M's calls alone stay out.
        let dir = std::env::temp_dir().join(format!("lrpipe-phase-samples-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let path = dir.join("x.vcf");
        let header = "##fileformat=VCFv4.2\n##contig=<ID=chrX,length=1000>\n\
                      ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\n";
        let records = [
            "chrX\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/1\t1",
            "chrX\t200\t.\tA\tG\t30\tPASS\t.\tGT\t./.\t1",
            "chrX\t300\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\t0",
            "chrX\t400\t.\tG\tA\t30\tPASS\t.\tGT\t1/1\t.",
        ];
        std::fs::write(&path, format!("{header}{}\n", records.join("\n"))).expect("VCF written");
        let mut reader = VcfReader::open(path.to_str().expect("UTF-8"), false).expect("opens");
        let table = reader.next_table().expect("no ploidy error").expect("chrX");
        let positions: Vec<i32> = table.variants.iter().map(|v| v.position).collect();
        assert_eq!(positions, [99, 299, 399], "only the site no sample is diploid at is left out");
        let (f, m) = (table.genotypes_of("F"), table.genotypes_of("M"));
        assert_eq!(f[0], Genotype::new(&[0, 1]).expect("0/1"), "F's call is kept");
        assert_eq!(f[1], Genotype::new(&[0, 1]).expect("0/1"));
        assert_eq!((m[0], m[1]), (Genotype::none(), Genotype::none()), "M's haploid calls are not");
    }
}

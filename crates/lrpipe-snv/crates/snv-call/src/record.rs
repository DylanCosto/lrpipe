//! Port of `output_with` (`CallVariants.py`): one VCF line per candidate.

use std::fmt::Write as _;

use crate::Error;
use crate::alt_info::{AltInfo, count_of, first_count_with_len, long_insertion_reads};
use crate::decode::{GENOTYPE_CLASSES, GT21_CLASSES, Heads, VARIANT_LENGTH_CLASSES, decode};
use crate::gvcf::genotype_likelihoods;

/// Where the probabilities come from; written as the INFO flag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallSource {
    /// Pileup model (`P`).
    Pileup,
    /// Full-alignment model (`F`).
    FullAlignment,
}

/// Haploid calling modes (`--haploid_precise`, `--haploid_sensitive`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Haploid {
    /// Diploid calls (the default).
    Off,
    /// Drop heterozygous calls.
    Precise,
    /// Drop multi-allelic calls.
    Sensitive,
}

/// Whether RefCall records are written (`--showRef`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefCalls {
    /// Write them (the original's default for its intermediate VCFs).
    Show,
    /// Drop them.
    Hide,
}

/// Whether non-ACGTN bases are kept in REF/ALT (`--keep_iupac_bases`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IupacBases {
    /// Use the VCF reference representation and replace unresolved ALT bases with `N` (default).
    Vcf,
    /// Replace ambiguous REF and ALT bases with `N`, as in Clair3 2.0.3.
    ToN,
    /// Keep them.
    Keep,
}

/// Whether records carry genotype likelihoods (`--gvcf`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Likelihoods {
    /// `GT:GQ:DP:AD:AF` (the default).
    #[default]
    Omit,
    /// `GT:GQ:DP:AD:AF:PL`, as for gVCF output.
    Write,
}

/// The `OutputConfig` fields that affect a VCF record.
#[derive(Debug, Clone, PartialEq)]
pub struct OutputConfig {
    /// Calls below this QUAL get `LowQual` (`--qual`, default 2); `None` passes everything.
    pub min_pass_quality: Option<f64>,
    /// RefCall records.
    pub ref_calls: RefCalls,
    /// Haploid mode.
    pub haploid: Haploid,
    /// IUPAC handling.
    pub iupac: IupacBases,
    /// Longest indel allele considered (`maximum_variant_length_that_need_infer`).
    pub max_indel_length: usize,
    /// INFO flag.
    pub source: CallSource,
    /// PL field (`gvcf`).
    pub likelihoods: Likelihoods,
}

impl Default for OutputConfig {
    /// `CallVariantsFromCffi --pileup` defaults.
    fn default() -> Self {
        Self {
            min_pass_quality: Some(2.0),
            ref_calls: RefCalls::Show,
            haploid: Haploid::Off,
            iupac: IupacBases::Vcf,
            max_indel_length: 50,
            source: CallSource::Pileup,
            likelihoods: Likelihoods::Omit,
        }
    }
}

/// `Phred_Trans = -10 * log(e, 10)`, with `math.log(x, base)` computed as `log(x) / log(base)`.
fn phred_trans() -> f64 {
    -10.0 * (std::f64::consts::E.ln() / 10f64.ln())
}

/// `quality_score_from`: `(1 - p + 1e-10) / (p + 1e-10)` in `float32` (numpy scalar arithmetic),
/// then a Phred transform in `float64`, clamped at 0 and rounded to 2 decimals.
pub(crate) fn quality_score(p: f32) -> f64 {
    let ratio = ((1.0f32 - p) + 1e-10f32) / (p + 1e-10f32);
    let q = phred_trans() * f64::from(ratio).ln() + 10.0;
    let q = if 0.0 > q { 0.0 } else { q };
    // Python's round(q, 2) is correctly rounded, like Rust's fixed-precision formatting.
    format!("{q:.2}").parse().expect("a formatted float parses")
}

/// `convert_iupac_to_n`.
fn iupac_to_n(s: &str) -> String {
    if s == "." {
        return s.to_owned();
    }
    s.chars().map(|c| if "ACGTN,.".contains(c.to_ascii_uppercase()) { c } else { 'N' }).collect()
}

/// The ALT alleles as written beside `written`, the REF written for `reference`: IUPAC codes
/// become N, except an indel's or MNP's first base, which is the REF's own first base (the anchor
/// both alleles share), so it is written as `written`'s is. Otherwise an indel anchored on an
/// ambiguous base would read `REF=ATT ALT=N`.
fn anchored_alternates(alternates: &str, reference: &str, written: &str) -> String {
    if alternates == "." {
        return alternates.to_owned();
    }
    let (Some(anchor), Some(written_anchor)) = (reference.chars().next(), written.chars().next()) else {
        return iupac_to_n(alternates);
    };
    alternates
        .split(',')
        .map(|allele| {
            let mut chars = allele.chars();
            match chars.next() {
                Some(first) if (allele.len() > 1 || reference.len() > 1) && first.eq_ignore_ascii_case(&anchor) => {
                    format!("{written_anchor}{}", iupac_to_n(chars.as_str()))
                }
                _ => iupac_to_n(allele),
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// Builds the VCF line the original writes for one candidate, or `None` if it writes none, for a
/// model without indel-length heads (`add_indel_length=False`, the pileup model).
///
/// `position` is `"{contig}:{pos}:{ref_base}"`, `alt_info` is `"{depth}-{alleles}"` and
/// `probabilities` holds 21 GT21 then 3 genotype probabilities (any further values are ignored).
///
/// # Errors
///
/// Fails where the original would raise: malformed inputs, a non-IUPAC reference base, or an
/// outcome it has no genotype for.
pub fn vcf_record(
    position: &str,
    alt_info: &str,
    probabilities: &[f32],
    config: &OutputConfig,
) -> Result<Option<String>, Error> {
    if probabilities.len() < GT21_CLASSES + GENOTYPE_CLASSES {
        return Err(Error::Candidate { position: position.to_owned() });
    }
    let (gt21, rest) = probabilities.split_at(GT21_CLASSES);
    let heads = Heads { gt21, genotype: &rest[..GENOTYPE_CLASSES], variant_lengths: None };
    record(position, alt_info, &heads, config)
}

/// Builds the VCF line the original writes for one candidate, or `None` if it writes none, for a
/// model with indel-length heads (`add_indel_length=True`, the full-alignment model).
///
/// As [`vcf_record`], but `probabilities` holds 21 GT21, 3 genotype and two times
/// [`VARIANT_LENGTH_CLASSES`] indel-length probabilities (90 values, as the model writes them;
/// any further values are ignored). The indel lengths weight every outcome and choose which
/// insertion or deletion of the alt info is called.
///
/// # Errors
///
/// Fails where the original would raise: malformed inputs, a non-IUPAC reference base, or an
/// outcome it has no genotype for.
pub fn vcf_record_with_indel_lengths(
    position: &str,
    alt_info: &str,
    probabilities: &[f32],
    config: &OutputConfig,
) -> Result<Option<String>, Error> {
    if probabilities.len() < GT21_CLASSES + GENOTYPE_CLASSES + 2 * VARIANT_LENGTH_CLASSES {
        return Err(Error::Candidate { position: position.to_owned() });
    }
    let (gt21, rest) = probabilities.split_at(GT21_CLASSES);
    let (genotype, rest) = rest.split_at(GENOTYPE_CLASSES);
    let (length_1, rest) = rest.split_at(VARIANT_LENGTH_CLASSES);
    let heads = Heads { gt21, genotype, variant_lengths: Some([length_1, &rest[..VARIANT_LENGTH_CLASSES]]) };
    record(position, alt_info, &heads, config)
}

/// `output_with` (non-gVCF, `enable_long_indel=False`).
fn record(position: &str, alt_info: &str, heads: &Heads<'_>, config: &OutputConfig) -> Result<Option<String>, Error> {
    let bad = || Error::Candidate { position: position.to_owned() };
    let fields: Vec<&str> = position.trim_end().split(':').collect();
    let (chrom, pos, reference_sequence) = match fields[..] {
        [c, p, r] => (c.to_owned(), p, r),
        [.., p, r] if fields.len() > 3 => (fields[..fields.len() - 2].join(":"), p, r),
        _ => return Err(bad()),
    };
    let pos: i64 = pos.parse().map_err(|_| bad())?;
    let alt = AltInfo::parse(alt_info).ok_or_else(bad)?;
    let call = decode(reference_sequence, heads, &alt, config.max_indel_length).ok_or_else(bad)?;
    let o = call.outcome;

    if (config.ref_calls == RefCalls::Hide && o.reference) || (!o.reference && call.reference == call.alternate) {
        return Ok(None);
    }
    let reference = call.reference;
    let mut alternate = call.alternate;
    let is_multi = alternate.contains(',');
    let any_hetero =
        o.hetero_snp || o.hetero_acgt_ins || o.hetero_ins_ins || o.hetero_acgt_del || o.hetero_del_del || o.ins_del;
    match config.haploid {
        Haploid::Precise if any_hetero => return Ok(None),
        Haploid::Sensitive if is_multi => return Ok(None),
        _ => {}
    }

    let mut genotype_string = if o.reference {
        Some("0/0")
    } else if o.homo_snp || o.homo_ins || o.homo_del {
        Some("1/1")
    } else if o.hetero_snp || o.hetero_acgt_ins || o.hetero_ins_ins || o.hetero_acgt_del || o.hetero_del_del {
        Some("0/1")
    } else {
        None
    };
    if is_multi {
        genotype_string = Some("1/2");
    }
    // Python would raise UnboundLocalError here.
    let genotype_string = genotype_string.ok_or_else(bad)?;

    let alleles = alt.decode();
    let ref_count = alleles.ref_count;
    // `supported_reads_for_long_ins`, with `--enable_long_indel` only.
    let long_ins = |ins: &str| if config.max_indel_length > 50 { long_insertion_reads(&alleles.ins, ins) } else { 0 };
    let mut supported = 0i64;
    let mut ad: Vec<i64> = Vec::new();
    if o.reference {
        supported = ref_count;
        alternate = ".".to_owned();
    } else if o.homo_snp || o.hetero_snp {
        for base in alternate.chars().filter(|&c| c != ',') {
            let n = count_of(&alleles.snp, &base.to_string());
            supported += n;
            ad.push(n);
        }
    } else if o.homo_ins || o.hetero_ins_ins {
        for ins in alternate.split(',') {
            let n = count_of(&alleles.ins, ins) + long_ins(ins);
            supported += n;
            ad.push(n);
        }
    } else if o.hetero_acgt_ins {
        let parts: Vec<&str> = alternate.split(',').collect();
        let (snp_base, ins) = if is_multi { (parts[0].get(..1), parts[1]) } else { (None, alternate.as_str()) };
        let snp_reads = snp_base.map_or(0, |b| count_of(&alleles.snp, b));
        let ins_reads = count_of(&alleles.ins, ins) + long_ins(ins);
        supported = ins_reads + snp_reads;
        if snp_base.is_some_and(|b| !b.is_empty()) {
            ad.push(snp_reads);
        }
        ad.push(ins_reads);
    } else if o.homo_del || o.hetero_del_del {
        if !alleles.del.is_empty() {
            if o.homo_del {
                let n = reference.get(1..).filter(|d| !d.is_empty()).map_or(0, |d| count_of(&alleles.del, d));
                supported = n;
                ad.push(n);
            } else if alleles.del.len() > 1 {
                for allele in alternate.split(',') {
                    let n = first_count_with_len(&alleles.del, reference.len() as i64 - allele.len() as i64);
                    ad.push(n);
                    supported += n;
                }
            }
        }
    } else if o.hetero_acgt_del {
        let parts: Vec<&str> = alternate.split(',').collect();
        let snp_base = if is_multi { parts.get(1).and_then(|p| p.get(..1)) } else { None };
        let snp_reads = snp_base.map_or(0, |b| count_of(&alleles.snp, b));
        let del_reads = reference.get(1..).filter(|d| !d.is_empty()).map_or(0, |d| count_of(&alleles.del, d));
        supported = del_reads + snp_reads;
        if snp_base.is_some_and(|b| !b.is_empty()) {
            ad.push(snp_reads);
        }
        ad.push(del_reads);
    } else if o.ins_del {
        for allele in alternate.split(',') {
            let diff = reference.len() as i64 - allele.len() as i64;
            let n = if diff < 0 {
                let ins = if reference.len() > 1 { &allele[..allele.len() - (reference.len() - 1)] } else { allele };
                count_of(&alleles.ins, ins) + long_ins(ins)
            } else {
                first_count_with_len(&alleles.del, diff)
            };
            ad.push(n);
            supported += n;
        }
    }

    let depth = alt.depth;
    let af = if depth != 0 { supported as f64 / depth as f64 } else { 0.0 };
    let af = if af > 1.0 { 1.0 } else { af };
    let quality = quality_score(call.probability);
    let genotype_string = match config.haploid {
        Haploid::Off => genotype_string,
        Haploid::Precise | Haploid::Sensitive if genotype_string.contains('1') => "1",
        Haploid::Precise | Haploid::Sensitive => "0",
    };
    let filter = if o.reference {
        "RefCall"
    } else if config.min_pass_quality.is_none_or(|min| quality >= min) {
        "PASS"
    } else {
        "LowQual"
    };
    let (reference, alternate) = match config.iupac {
        IupacBases::Vcf => {
            let mapped = reference.bytes().map(lrpipe_reads::vcf_reference_base).map(char::from).collect::<String>();
            let same = |allele: &str| allele.eq_ignore_ascii_case(&mapped);
            let written = match (o.reference, iupac_to_n(&alternate).split(',').filter(|allele| same(allele)).count()) {
                (true, _) | (false, 0) => mapped,
                // The one ALT is the base that stands for the ambiguous REF: no variant.
                (false, _) if !alternate.contains(',') => return Ok(None),
                // One of several ALTs is: the others are still variants, so the record is kept
                // with the ambiguous REF as N, as Clair3 writes it, rather than dropped.
                (false, _) => iupac_to_n(&reference),
            };
            let alternate = anchored_alternates(&alternate, &reference, &written);
            (written, alternate)
        }
        IupacBases::ToN => (iupac_to_n(&reference), iupac_to_n(&alternate)),
        IupacBases::Keep => (reference, alternate),
    };

    let mut allele_depth = ref_count.to_string();
    for n in &ad {
        write!(allele_depth, ",{n}").expect("writing to a String cannot fail");
    }
    let af_string = if ad.len() <= 1 {
        format!("{af:.4}")
    } else if depth == 0 {
        // `1.0 * item / read_depth` raises ZeroDivisionError.
        return Err(bad());
    } else {
        ad.iter().map(|&n| format!("{:.4}", f64::min(1.0, n as f64 / depth as f64))).collect::<Vec<_>>().join(",")
    };
    let info = match config.source {
        CallSource::Pileup => "P",
        CallSource::FullAlignment => "F",
    };
    let mut line = format!(
        "{chrom}\t{pos}\t.\t{reference}\t{alternate}\t{quality:.2}\t{filter}\t{info}\tGT:GQ:DP:AD:AF{pl_key}\t{genotype_string}:{gq}:{depth}:{allele_depth}:{af_string}",
        gq = quality as i64,
        pl_key = if config.likelihoods == Likelihoods::Write { ":PL" } else { "" },
    );
    if config.likelihoods == Likelihoods::Write {
        // Computed from the bases as written, after the IUPAC replacement.
        let pl = genotype_likelihoods(&reference, &alternate, heads.gt21, heads.genotype).ok_or_else(bad)?;
        let pl: Vec<String> = pl.iter().map(i64::to_string).collect();
        write!(line, ":{}", pl.join(",")).expect("writing to a String cannot fail");
    }
    line.push('\n');
    Ok(Some(line))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quality_rounds_and_clamps() {
        assert_eq!(quality_score(0.0), 0.0);
        // p = 0.9: 10 * log10(0.9 / 0.1) + 10 = 19.54
        assert_eq!(quality_score(0.9), 19.54);
    }

    #[test]
    fn an_indel_keeps_its_anchor_base_beside_an_ambiguous_ref() {
        // A deletion and an insertion anchored on R, written with R as A.
        assert_eq!(anchored_alternates("R", "RTT", "ATT"), "A");
        assert_eq!(anchored_alternates("RACG", "R", "A"), "AACG");
        // With the REF as N (one of several ALTs is the ambiguous base's), the anchors are N too.
        assert_eq!(anchored_alternates("A,RTT", "R", "N"), "A,NTT");
        // A SNV's ALT is a base of its own, not an anchor; other IUPAC codes still become N.
        assert_eq!(anchored_alternates("G", "R", "A"), "G");
        assert_eq!(anchored_alternates("RMG", "R", "A"), "ANG");
        assert_eq!(anchored_alternates(".", "R", "A"), ".");
        // An unambiguous REF is unchanged.
        assert_eq!(anchored_alternates("A", "ATT", "ATT"), "A");
    }

    #[test]
    fn iupac_bases_become_n() {
        assert_eq!(iupac_to_n("AMG,T"), "ANG,T");
        assert_eq!(iupac_to_n("."), ".");
        assert_eq!(iupac_to_n("acg"), "acg");
    }

    fn probs(gt21_idx: usize, gt21_p: f32, genotype_idx: usize, genotype_p: f32) -> Vec<f32> {
        let mut p = vec![(1.0 - gt21_p) / 20.0; 21];
        p[gt21_idx] = gt21_p;
        let mut g = vec![(1.0 - genotype_p) / 2.0; 3];
        g[genotype_idx] = genotype_p;
        p.extend(g);
        p
    }

    // Expected lines are the original's `output_with` on the same float32 inputs.
    #[test]
    fn het_snp_record() {
        // AG (index 2) with 0/1 (index 2).
        let line = vcf_record("chr1:100:A", "30-XG 14 RA 16", &probs(2, 0.9, 2, 0.9), &OutputConfig::default());
        assert_eq!(
            line.unwrap().unwrap(),
            "chr1\t100\t.\tA\tG\t16.30\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:16:30:16,14:0.4667\n"
        );
    }

    #[test]
    fn ref_call_record() {
        let line = vcf_record("chr1:100:A", "30-XG 2 RA 28", &probs(0, 0.9, 0, 0.9), &OutputConfig::default());
        assert_eq!(
            line.unwrap().unwrap(),
            "chr1\t100\t.\tA\t.\t16.30\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:16:30:28:0.9333\n"
        );
    }
}

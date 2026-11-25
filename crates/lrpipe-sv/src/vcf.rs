//! VCF output in Sniffles2's format (`vcf.py`).

use std::borrow::Cow;
use std::fmt::Write as _;

use crate::call::{Alt, Genotype, Info, Phase, SvCall};
use crate::config::Config;
use crate::fasta::Contig;
use crate::pymath;
use crate::reads::SvType;
use crate::snf::SampleGenotype;

const HEADER_BODY: &str = r#"##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend; Translocation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">
##FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=GT,Description="Genotype filter">
##FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">
##FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">
##FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">
##FILTER=<ID=COV_MIN,Description="Minimum coverage filter">
##FILTER=<ID=COV_MIN_GT,Description="Minimum coverage filter (missing genotype)">
##FILTER=<ID=COV_CHANGE_DEL,Description="Coverage change filter for DEL">
##FILTER=<ID=COV_CHANGE_DUP,Description="Coverage change filter for DUP">
##FILTER=<ID=COV_CHANGE_INS,Description="Coverage change filter for INS">
##FILTER=<ID=COV_CHANGE_FRAC_US,Description="Coverage fractional change filter: upstream-start">
##FILTER=<ID=COV_CHANGE_FRAC_SC,Description="Coverage fractional change filter: start-center">
##FILTER=<ID=COV_CHANGE_FRAC_CE,Description="Coverage fractional change filter: center-end">
##FILTER=<ID=COV_CHANGE_FRAC_ED,Description="Coverage fractional change filter: end-downstream">
##FILTER=<ID=COV_VAR,Description="Coverage variance exceeded">
##FILTER=<ID=MOSAIC_VAF,Description="Mosaic variant allele fraction filter">
##FILTER=<ID=NOT_MOSAIC_VAF,Description="Variant allele fraction filter for non-mosaic">
##FILTER=<ID=ALN_NM,Description="Length adjusted mismatch filter">
##FILTER=<ID=STRAND_BND,Description="Strand support filter for BNDs">
##FILTER=<ID=STRAND,Description="Strand support filter for germline SVs">
##FILTER=<ID=STRAND_MOSAIC,Description="Strand support filter for mosaic SVs">
##FILTER=<ID=SVLEN_MIN,Description="SV length filter">
##FILTER=<ID=SVLEN_MIN_MOSAIC,Description="SV length filter for mosaic SVs (min)">
##FILTER=<ID=SVLEN_MAX_MOSAIC,Description="SV length filter for mosaic SVs (max)">
##FILTER=<ID=SINGLE_BREAK,Description="A single break point was detected but not classified as an SV.">
##FILTER=<ID=INLINE_SA,Description="INLINE/CIGAR-based SV is mostly supported by SA reads">
##FILTER=<ID=MOSAIC_SV_CLOSE_EDGE,Description="For mosaic SVs, the location is close to the end of the read (either end)">
##FILTER=<ID=GT_FAILED,Description="Sniffles was unable to genotype this call.">
##FILTER=<ID=UNRESOLVED,Description="Call is within an unresolved region.">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
##INFO=<ID=MOSAIC,Number=0,Type=Flag,Description="Structural variation classified as putative mosaic">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
"#;

const HEADER_BODY2: &str = r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
##INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only), or reads for a DUP SV reported as INS by the mapper.">
##INFO=<ID=SUPPORT_SA,Number=1,Type=Integer,Description="Number of reads supporting a DEL SV through supplementary alignments (split events)">
##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">
##INFO=<ID=SUPPORT_UNSCALED,Number=1,Type=Float,Description="Base support for long insertions, before adjustment">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">
##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">
##INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">
##INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">
##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">
##INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">
##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction">
##INFO=<ID=COVERAGE_VAR,Number=1,Type=Float,Description="Variance of coverage across large events">
##INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">
##INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">
##INFO=<ID=LASM,Number=0,Type=Flag,Description="Local assembly used to detect the structural variant">
"#;

/// QUAL is clamped to this.
const MAX_QUAL: i64 = 60;

/// The VCF header (Sniffles2's, with this program in `source` and `command`).
pub fn header(cfg: &Config, contigs: &[(String, i64)], samples: &[String], date: &str) -> String {
    let mut h = String::new();
    h.push_str("##fileformat=VCFv4.2\n");
    let _ = writeln!(h, "##source=lrpipe-sv_{}", env!("CARGO_PKG_VERSION"));
    let _ = writeln!(h, "##command=\"{}\"", cfg.command);
    let _ = writeln!(h, "##fileDate=\"{date}\"");
    for (name, len) in contigs {
        let _ = writeln!(h, "##contig=<ID={name},length={len}>");
    }
    h.push_str(HEADER_BODY);
    if cfg.args.dev_emit_sv_lengths {
        h.push_str(
            "##INFO=<ID=SVLENGTHS,Number=.,Type=Integer,Description=\"Lengths of structural variation (all)\">\n",
        );
    }
    h.push_str(HEADER_BODY2);
    let _ = writeln!(h, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", samples.join("\t"));
    h
}

/// The reference sequence available for a contig.
#[derive(Debug, Clone, Copy)]
pub enum Reference<'a> {
    /// No `--reference` given.
    Absent,
    /// The reference has no sequence of this name.
    Missing,
    Present(Contig<'a>),
}

fn format_info(out: &mut String, k: &str, v: &Info) {
    match v {
        Info::Int(i) => {
            let _ = write!(out, "{k}={i}");
        }
        Info::Float(f) => {
            let _ = write!(out, "{k}={}", pymath::fmt3(*f));
        }
        Info::Str(s) => {
            let _ = write!(out, "{k}={s}");
        }
        Info::Bool(true) => out.push_str(k),
        Info::Bool(false) => {
            let _ = write!(out, "{k}=False");
        }
        Info::None => {
            let _ = write!(out, "{k}=.");
        }
    }
}

pub fn format_genotype(out: &mut String, gt: Option<&Genotype>) {
    match gt {
        Some(g) => format_gt(out, Some(g.a), Some(g.b), g.gq, g.dr, g.dv, &g.phase),
        None => out.push_str("./.:0:0:0:."),
    }
}

/// A multi-sample genotype without its ID column.
pub fn format_sample_genotype(out: &mut String, g: &SampleGenotype) {
    format_gt(out, g.a, g.b, g.gq, g.dr, g.dv, &g.phase);
}

/// `format_genotype` with phasing: `.` alleles are `None`.
fn format_gt(out: &mut String, a: Option<i64>, b: Option<i64>, gq: i64, dr: i64, dv: i64, phase: &Phase) {
    let ps = match phase.ps.as_deref() {
        None | Some("NULL") => ".",
        Some(p) => p,
    };
    let (mut a, mut b) = (a, b);
    let sep = if phase.hp.is_some() && matches!((a, b), (Some(0), Some(1)) | (Some(1), Some(1))) {
        if phase.hp.as_deref() == Some("1") {
            (a, b) = (b, a);
        }
        '|'
    } else {
        '/'
    };
    let allele = |x: Option<i64>| x.map_or_else(|| ".".to_string(), |v| v.to_string());
    let _ = write!(out, "{}{sep}{}:{gq}:{dr}:{dv}:{ps}", allele(a), allele(b));
}

fn iupac_to_n(s: &mut [u8]) {
    for c in s {
        if b"RYSWKMBDHV".contains(c) {
            *c = b'N';
        }
    }
}

/// [`iupac_to_n`] on an ALT's bases only: a breakend's mate (`t[p[`, `t]p]`, `]p]t`, `[p[t`) is
/// left as it is, so `chrY` and `chrM` stay themselves. Sniffles2 rewrites the whole ALT, which
/// turns a mate on chrY into `chrN` while `CHR2` still says chrY (LRP-004).
fn iupac_to_n_in_bases(alt: &mut [u8]) {
    let brackets = |c: &u8| *c == b'[' || *c == b']';
    match (alt.iter().position(brackets), alt.iter().rposition(brackets)) {
        (Some(first), Some(last)) => {
            let (before, rest) = alt.split_at_mut(first);
            iupac_to_n(before);
            iupac_to_n(&mut rest[last + 1 - first..]);
        }
        _ => iupac_to_n(alt),
    }
}

/// The fields of one output record, from a single-sample call or a combined call.
#[derive(Debug)]
pub struct Record<'a> {
    pub contig: &'a str,
    pub svtype: SvType,
    pub pos: i64,
    pub end: i64,
    pub svlen: i64,
    pub id: &'a str,
    pub alt: &'a Alt,
    pub qual: Option<i64>,
    pub filter: &'a str,
    pub precise: bool,
    pub support: i64,
    /// Supporting read names (only written with `--output-rnames`).
    pub rnames: Vec<Cow<'a, str>>,
    /// Coverage upstream, at start, center, end and downstream (`None` is written as Python does).
    pub coverage: [Option<i64>; 5],
    pub fwd: i64,
    pub rev: i64,
    pub nm: Info,
    pub svlens: Option<&'a [i64]>,
    /// Additional INFO fields; written sorted by key.
    pub info: Vec<(&'a str, &'a Info)>,
}

impl<'a> Record<'a> {
    /// A single-sample call's record.
    pub fn from_call(call: &'a SvCall, contig: &'a str, cfg: &Config) -> Self {
        let rnames = if cfg.args.output_rnames {
            call.rnames.iter().map(|n| String::from_utf8_lossy(n)).collect()
        } else {
            Vec::new()
        };
        Record {
            contig,
            svtype: call.svtype,
            pos: call.pos,
            end: call.end,
            svlen: call.svlen,
            id: &call.id,
            alt: &call.alt,
            qual: Some(call.qual),
            filter: call.filter,
            precise: call.precise,
            support: call.support as i64,
            rnames,
            coverage: [call.cov_up, call.cov_start, call.cov_center, call.cov_end, call.cov_down].map(Some),
            fwd: call.fwd as i64,
            rev: call.rev as i64,
            nm: Info::Float(call.nm),
            svlens: call.svlens.as_deref(),
            info: call.info.iter().map(|(k, v)| (*k, v)).collect(),
        }
    }
}

/// `VCF.write_call` for a single-sample call: appends the record line, or nothing if Sniffles2
/// would skip the call.
///
/// # Errors
///
/// Fails if the reference cannot be read.
pub fn write_call(
    out: &mut String,
    call: &SvCall,
    contig: &str,
    cfg: &Config,
    reference: Reference<'_>,
) -> std::io::Result<()> {
    if write_record(out, Record::from_call(call, contig, cfg), "GT:GQ:DR:DV:PS", cfg, reference)? {
        format_genotype(out, call.genotype.as_ref());
        out.push('\n');
    }
    Ok(())
}

/// `VCF.write_call` up to the sample columns: appends the fixed columns and FORMAT, each followed
/// by a tab, and returns true, or appends nothing and returns false if Sniffles2 skips the record.
///
/// # Errors
///
/// Fails if the reference cannot be read.
pub fn write_record(
    out: &mut String,
    mut rec: Record<'_>,
    format: &str,
    cfg: &Config,
    reference: Reference<'_>,
) -> std::io::Result<bool> {
    let pos = if rec.pos > 0 { rec.pos } else { 1 };
    let mut end = rec.end;
    if rec.precise && rec.svtype == SvType::Del {
        end = pos + rec.svlen.abs();
    }
    if rec.svtype == SvType::Ins {
        if let Alt::Seq(s) = rec.alt
            && rec.svlen != s.len() as i64
            && !cfg.args.symbolic
        {
            rec.svlen = s.len() as i64;
        }
        if rec.svlen < cfg.minsvlen {
            return Ok(false);
        }
    }
    let info = info_field(&mut rec, cfg, end);
    let Some((ref_seq, alt)) = ref_and_alt(&rec, cfg, reference)? else { return Ok(false) };
    let qual = rec.qual.map_or_else(|| ".".to_string(), |q| q.clamp(0, MAX_QUAL).to_string());
    let _ = write!(
        out,
        "{}\t{pos}\tSniffles2.{}\t{}\t{}\t{qual}\t{}\t{info}\t{format}\t",
        rec.contig,
        rec.id,
        String::from_utf8_lossy(&ref_seq),
        String::from_utf8_lossy(&alt),
        rec.filter
    );
    Ok(true)
}

/// Python's `str()` of an optional integer.
fn py_opt(v: Option<i64>) -> String {
    v.map_or_else(|| "None".to_string(), |x| x.to_string())
}

/// The INFO column: fixed fields in Sniffles2's order, then the record's extra fields by name.
fn info_field(rec: &mut Record<'_>, cfg: &Config, end: i64) -> String {
    let mut info = String::new();
    info.push_str(if rec.precise { "PRECISE" } else { "IMPRECISE" });
    rec.info.sort_by_key(|(k, _)| *k);
    let af = match rec.info.iter().find(|(k, _)| *k == "VAF").map(|(_, v)| *v) {
        Some(Info::Float(f)) => *f,
        Some(Info::Int(i)) => *i as f64,
        _ => 0.0,
    };
    if af <= cfg.args.mosaic_af_max && cfg.args.mosaic {
        info.push_str(";MOSAIC");
    }
    let _ = write!(info, ";SVTYPE={}", rec.svtype.as_str());
    if rec.svtype != SvType::Bnd {
        let _ = write!(info, ";SVLEN={};END={end}", rec.svlen);
    }
    let _ = write!(info, ";SUPPORT={}", rec.support);
    if cfg.args.output_rnames {
        let _ = write!(info, ";RNAMES={}", rec.rnames.join(","));
    }
    let [up, start, center, cend, down] = rec.coverage.map(py_opt);
    let _ = write!(info, ";COVERAGE={up},{start},{center},{cend},{down}");
    let _ = write!(info, ";STRAND={}{}", if rec.fwd > 0 { "+" } else { "" }, if rec.rev > 0 { "-" } else { "" });
    if cfg.qc_nm_measure {
        info.push(';');
        format_info(&mut info, "NM", &rec.nm);
    }
    if cfg.args.dev_emit_sv_lengths
        && rec.svtype != SvType::Bnd
        && let Some(l) = rec.svlens.filter(|l| !l.is_empty())
    {
        let s: Vec<String> = l.iter().map(i64::to_string).collect();
        let _ = write!(info, ";SVLENGTHS={}", s.join(","));
    }
    for (k, v) in &rec.info {
        if **v == Info::None {
            continue;
        }
        info.push(';');
        format_info(&mut info, k, v);
    }
    info
}

/// REF and ALT as Sniffles2 writes them, or `None` if it skips the call (a deletion whose
/// reference sequence is mostly `N`).
fn ref_and_alt(
    call: &Record<'_>,
    cfg: &Config,
    reference: Reference<'_>,
) -> std::io::Result<Option<(Vec<u8>, Vec<u8>)>> {
    let symbolic = || format!("<{}>", call.svtype.as_str()).into_bytes();
    let mut ref_seq: Vec<u8> = b"N".to_vec();
    let mut alt: Vec<u8> = match call.alt {
        Alt::Seq(s) => s.clone(),
        Alt::Symbolic => symbolic(),
    };
    if !cfg.args.symbolic
        && call.svtype == SvType::Del
        && let Reference::Present(fa) = reference
        && call.svlen.abs() <= cfg.args.max_del_seq_len
    {
        // pysam rejects a negative start (ValueError) and an unknown contig (KeyError).
        let fetched = if call.pos > 0 { Some(fa.fetch(call.pos - 1, call.pos - call.svlen)?) } else { None };
        match fetched {
            Some(r) if !r.is_empty() => {
                let n = r.iter().filter(|&&c| c == b'N').count();
                if n > 0 && n as f64 / r.len() as f64 > cfg.args.max_unknown_pct {
                    return Ok(None);
                }
                alt = vec![r[0]];
                ref_seq = r;
            }
            _ => {
                ref_seq = b"N".to_vec();
                alt = symbolic();
            }
        }
    }
    if cfg.args.symbolic {
        ref_seq = b"N".to_vec();
        if call.svtype != SvType::Bnd {
            alt = symbolic();
        }
    } else if !matches!(reference, Reference::Absent) && ref_seq == b"N" {
        // Sniffles2 only cleans up IUPAC codes in this branch (and in BND mate names, which
        // lrpipe-sv leaves alone).
        if let Reference::Present(fa) = reference {
            let start = (call.pos - 1).max(0);
            ref_seq = fa.fetch(start, start + 1)?;
            let is_symbolic = alt.first() == Some(&b'<');
            if call.svtype == SvType::Ins && !is_symbolic {
                alt = [&ref_seq[..], &alt[..]].concat();
            } else if call.svtype == SvType::Bnd && !is_symbolic {
                alt = if alt.first() == Some(&b'N') {
                    [&ref_seq[..], &alt[1..]].concat()
                } else {
                    [&alt[..alt.len() - 1], &ref_seq[..]].concat()
                };
            }
        }
        if ref_seq.first() != Some(&b'<') {
            iupac_to_n(&mut ref_seq);
        }
        if alt.first() != Some(&b'<') {
            iupac_to_n_in_bases(&mut alt);
        }
    }
    Ok(Some((ref_seq, alt)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gt(a: i64, b: i64, hp: Option<&str>, ps: Option<&str>) -> Genotype {
        Genotype {
            a,
            b,
            gq: 40,
            dr: 10,
            dv: 12,
            phase: Phase { hp: hp.map(str::to_string), ps: ps.map(str::to_string) },
        }
    }

    fn format(g: Option<&Genotype>) -> String {
        let mut s = String::new();
        format_genotype(&mut s, g);
        s
    }

    #[test]
    fn formats_genotypes_like_sniffles2() {
        assert_eq!(format(None), "./.:0:0:0:.");
        assert_eq!(format(Some(&gt(0, 1, None, None))), "0/1:40:10:12:.");
        // HP 1 puts the alternative allele first.
        assert_eq!(format(Some(&gt(0, 1, Some("1"), Some("777")))), "1|0:40:10:12:777");
        assert_eq!(format(Some(&gt(0, 1, Some("2"), Some("NULL")))), "0|1:40:10:12:.");
        // Hom-ref is never phased.
        assert_eq!(format(Some(&gt(0, 0, Some("1"), Some("5")))), "0/0:40:10:12:5");
    }

    #[test]
    fn replaces_iupac_codes() {
        let mut s = b"ACGTRYKMN".to_vec();
        iupac_to_n(&mut s);
        assert_eq!(s, b"ACGTNNNNN");
    }

    #[test]
    fn breakend_mates_keep_their_names() {
        // LRP-004: all four orientations, mates whose names hold IUPAC letters.
        for (alt, want) in [
            ("R[chrY:10001[", "N[chrY:10001["),
            ("A]chrM:5]", "A]chrM:5]"),
            ("]MT:7]K", "]MT:7]N"),
            ("[HSCHRUN_RANDOM_CTG1:9[S", "[HSCHRUN_RANDOM_CTG1:9[N"),
            ("ACGTRY", "ACGTNN"),
        ] {
            let mut s = alt.as_bytes().to_vec();
            iupac_to_n_in_bases(&mut s);
            assert_eq!(String::from_utf8(s).expect("ascii"), want, "{alt}");
        }
    }
}

//! Genotyping known SVs (`--genotype-vcf`), as Sniffles2's `GenotypeTask` and VCF rewriting.

use std::collections::HashMap;
use std::io::Read;

use anyhow::{Context, Result, bail};

use crate::call::{Coverage, Genotype, Phase, SvCall};
use crate::config::Config;
use crate::reads::SvType;
use crate::vcf;

/// Input SVs are binned by position in bins of this size.
const BIN: i64 = 5000;
/// SVs this close to a bin edge are also placed in the neighbouring bin.
const BIN_EDGE: i64 = BIN / 10;

/// One SV from the input VCF.
#[derive(Debug, Clone)]
pub struct InputSv {
    /// The input line, trimmed.
    pub raw_line: String,
    pub contig: String,
    /// 0-based position (VCF POS - 1).
    pub pos: i64,
    pub svtype: String,
    pub svlen: i64,
    /// The breakend's mate and orientation, for BNDs.
    pub bnd: Option<InputBnd>,
}

/// A breakend of the input: its mate and both sides' orientation, as a BND ALT gives them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputBnd {
    pub mate_contig: String,
    /// The mate's position as the ALT writes it.
    pub mate_pos: i64,
    /// The local sequence comes first (`t[p[`, `t]p]`), not after the mate (`]p]t`, `[p[t`).
    pub is_first: bool,
    /// The mate's side is `]p]`, not `[p[`.
    pub is_reverse: bool,
}

impl InputBnd {
    /// Parses a BND ALT in any of the VCF's four forms. The local side comes from where the
    /// brackets are, not from an `N`: Sniffles2 reads `alt[0] == 'N'`, which a reference base in
    /// its place (as Sniffles2 writes with `--reference`) turns into the other orientation.
    fn parse(alt: &str) -> Result<Self> {
        let parts: Vec<&str> = alt.split(['[', ']']).collect();
        if parts.len() <= 2 {
            bail!("BND ALT not formatted according to VCF 4.2 specifications");
        }
        let (mate_contig, mate_pos) =
            parts[1].split_once(':').context("BND ALT not formatted according to VCF 4.2 specifications")?;
        Ok(Self {
            mate_contig: mate_contig.to_string(),
            mate_pos: mate_pos.parse()?,
            is_first: !alt.starts_with(['[', ']']),
            is_reverse: alt.contains(']'),
        })
    }

    /// Whether a called breakend is this one: the same mate contig, the mate within `tolerance`,
    /// and the same orientation on both sides. Sniffles2 compares the mate contig alone, so a
    /// target took the support of any breakend near it to that contig, however far its mate or
    /// whichever way it joins (LRP-015).
    fn matches(&self, call: &crate::reads::BndInfo, tolerance: i64) -> bool {
        *self.mate_contig == *call.mate_contig
            && (self.mate_pos - call.mate_ref_start).abs() <= tolerance
            && self.is_first == call.is_first
            && self.is_reverse == call.is_reverse
    }
}

/// The parsed `--genotype-vcf` input.
#[derive(Debug)]
pub struct Input {
    /// Header lines (each trimmed, newline-terminated).
    pub header: String,
    pub svs_by_contig: HashMap<String, Vec<InputSv>>,
}

/// Reads the input VCF (plain or gzipped), as Sniffles2's `VCF.read_svs_iter`.
///
/// # Errors
///
/// Fails on unreadable files and malformed records (as Sniffles2's fatal parse errors).
pub fn read_input(path: &str) -> Result<Input> {
    let mut bytes = Vec::new();
    std::fs::File::open(path).with_context(|| format!("opening {path}"))?.read_to_end(&mut bytes)?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut text = Vec::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_end(&mut text)
            .with_context(|| format!("reading {path}"))?;
        bytes = text;
    }
    let text = String::from_utf8(bytes).with_context(|| format!("{path} is not UTF-8 text"))?;
    let mut header = String::new();
    let mut svs_by_contig: HashMap<String, Vec<InputSv>> = HashMap::new();
    for (k, line) in text.lines().enumerate() {
        let line_index = k + 1;
        let strip = line.trim();
        if strip.starts_with('#') {
            header.push_str(strip);
            header.push('\n');
            continue;
        }
        if strip.is_empty() {
            bail!("Error parsing input VCF: Line {line_index}: empty line");
        }
        // Sniffles2 splits the line with its newline, which an 8-column record's INFO keeps: its
        // last field's value ends in "\n", so a sites-only record ending in SVTYPE=BND was an
        // unsupported "BND\n" and went ungenotyped. The line is split without it.
        let sv = parse_record(line.trim_end_matches(['\r', '\n']), strip)
            .with_context(|| format!("Error parsing input VCF: Line {line_index}"))?;
        svs_by_contig.entry(sv.contig.clone()).or_default().push(sv);
    }
    Ok(Input { header, svs_by_contig })
}

fn parse_record(line: &str, stripped: &str) -> Result<InputSv> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        bail!("expected at least 8 columns");
    }
    let (contig, pos, reference, alt, info) = (f[0], f[1].trim().parse::<i64>()? - 1, f[3], f[4], f[7]);
    let mut fields: HashMap<&str, Option<&str>> = HashMap::new();
    for item in info.split(';') {
        let mut kv = item.split('=');
        let key = kv.next().unwrap_or_default();
        let value = kv.next();
        if kv.next().is_some() {
            bail!("malformed INFO field '{item}'");
        }
        fields.insert(key, value);
    }
    let (mut svtype, mut svlen) = if alt.len() > reference.len() {
        ("INS".to_string(), alt.len() as i64)
    } else {
        ("DEL".to_string(), -(reference.len() as i64))
    };
    if let Some(&Some(t)) = fields.get("SVTYPE") {
        svtype = if t == "TRA" { "BND".to_string() } else { t.to_string() };
    }
    if let Some(&Some(l)) = fields.get("SVLEN") {
        svlen = l.trim().parse()?;
    }
    let bnd = if svtype == "BND" { Some(InputBnd::parse(alt)?) } else { None };
    Ok(InputSv { raw_line: stripped.to_string(), contig: contig.to_string(), pos, svtype, svlen, bnd })
}

/// The output header: the input's, with Sniffles2's `##genotype*` lines, any missing FORMAT
/// definition of the fields each record has (GT, GQ, DR, DV and PS), and a `#CHROM` line naming
/// the one sample each record has: the genotyped one (`--sample-id`). Sniffles2 keeps the
/// input's `#CHROM` line and leaves out PS, so a sites-only input gave 8 header columns over 10
/// record columns, a multi-sample one more than the records have, and PS was undeclared
/// (LRP-008).
pub fn header(input: &Input, cfg: &Config, date: &str) -> String {
    let mut lines: Vec<String> = input.header.split('\n').map(str::to_string).collect();
    let chrom = format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", cfg.sample_id);
    match lines.iter().position(|l| l.starts_with("#CHROM")) {
        Some(at) => lines[at] = chrom,
        None => {
            let at = lines.iter().rposition(|l| !l.is_empty()).map_or(0, |i| i + 1);
            lines.insert(at, chrom);
        }
    }
    let at = 1.min(lines.len());
    lines.insert(at, format!("##genotypeFileDate=\"{date}\""));
    lines.insert(at, format!("##genotypeCommand=\"{}\"", cfg.command));
    lines.insert(at, format!("##genotypeSource=lrpipe-sv_{}", env!("CARGO_PKG_VERSION")));
    let definitions = [
        ("GT", r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#),
        ("GQ", r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">"#),
        ("DR", r#"##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">"#),
        ("DV", r#"##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">"#),
        ("PS", r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">"#),
    ];
    for (id, def) in definitions {
        let tag = format!("##FORMAT=<ID={id},");
        if !lines.iter().any(|l| l.contains(&tag)) {
            let at = lines.len().saturating_sub(2);
            lines.insert(at, def.to_string());
        }
    }
    lines.join("\n")
}

/// The genotyped records of one contig's input SVs, in input order (`GenotypeTask.execute` and
/// `VCF.rewrite_genotype`).
pub fn write_records(svs: &[InputSv], candidates: &[SvCall], coverage: &Coverage<'_>, cfg: &Config) -> String {
    // Input SVs by (type, bin), each also in the neighbouring bin when near an edge.
    let mut bins: HashMap<(SvType, i64), Vec<usize>> = HashMap::new();
    for (k, sv) in svs.iter().enumerate() {
        let Some(t) = svtype_of(&sv.svtype) else {
            eprintln!("Unsupported SVTYPE: {}", sv.svtype);
            continue;
        };
        let b = sv.pos / BIN * BIN;
        let mut targets = vec![b];
        if sv.pos % BIN < BIN_EDGE {
            targets.push(b - BIN);
        }
        if sv.pos % BIN > BIN - BIN_EDGE {
            targets.push(b + BIN);
        }
        for bin in targets {
            bins.entry((t, bin)).or_default().push(k);
        }
    }
    let mut best: Vec<Option<(usize, f64)>> = vec![None; svs.len()];
    for (ci, cand) in candidates.iter().enumerate() {
        let Some(list) = bins.get(&(cand.svtype, cand.pos / BIN * BIN)) else { continue };
        for &k in list {
            let sv = &svs[k];
            let best_dist = best[k].map_or(f64::INFINITY, |(_, d)| d);
            if cand.svtype == SvType::Bnd {
                let dist = (sv.pos - cand.pos).abs() as f64;
                let tolerance = cfg.args.cluster_merge_bnd;
                let same = matches!((&sv.bnd, &cand.bnd), (Some(t), Some(c)) if t.matches(c, tolerance));
                if dist < best_dist && dist <= tolerance as f64 && same {
                    best[k] = Some((ci, dist));
                }
            } else {
                let dist = ((sv.pos - cand.pos).abs() + (sv.svlen.abs() - cand.svlen.abs()).abs()) as f64;
                let minlen = sv.svlen.abs().min(cand.svlen.abs()) as f64;
                if minlen > 0.0
                    && dist < best_dist
                    && dist <= cfg.args.combine_match as f64 * minlen.sqrt()
                    && dist <= cfg.args.combine_match_max as f64
                {
                    best[k] = Some((ci, dist));
                }
            }
        }
    }
    let mut out = String::new();
    for (k, sv) in svs.iter().enumerate() {
        let matched = best[k].and_then(|(ci, _)| candidates[ci].genotype.as_ref());
        let fallback = input_genotype(sv, coverage, cfg);
        let gt = matched.or(fallback.as_ref());
        let cols: Vec<&str> = sv.raw_line.split('\t').take(8).collect();
        out.push_str(&cols.join("\t"));
        // Each genotype ends with its PS, so FORMAT names it; Sniffles2 leaves PS out of FORMAT,
        // which makes its genotyped VCFs unreadable to bcftools (LRP-008).
        out.push_str("\tGT:GQ:DR:DV:PS\t");
        vcf::format_genotype(&mut out, gt);
        out.push('\n');
    }
    out
}

fn svtype_of(s: &str) -> Option<SvType> {
    SvType::ALL.into_iter().find(|t| t.as_str() == s)
}

/// 0/0 with the mean coverage at the SV as DR, or `None` (./.) without coverage.
fn input_genotype(sv: &InputSv, cv: &Coverage<'_>, cfg: &Config) -> Option<Genotype> {
    let bs = cfg.coverage_binsize;
    let mut start = sv.pos;
    let (s, c, e) = match sv.svtype.as_str() {
        "INS" | "BND" => {
            if sv.svtype == "BND" && sv.bnd.as_ref().is_some_and(|b| b.is_first) {
                start -= 1;
            }
            let end = start + 1;
            (cv.get(start - bs), cv.get(start), cv.get(end + bs))
        }
        _ => {
            let end = sv.pos + sv.svlen.abs();
            (cv.get(start), cv.get((start + end) / 2), cv.get(end - bs))
        }
    };
    let covs = [s.unwrap_or(0), c.unwrap_or(0), e.unwrap_or(0)];
    let coverage = crate::pymath::py_round(covs.iter().sum::<i64>() as f64 / 3.0);
    (coverage > 0).then(|| Genotype { a: 0, b: 0, gq: 0, dr: coverage, dv: 0, phase: Phase::default() })
}

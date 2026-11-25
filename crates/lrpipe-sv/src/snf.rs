//! Sniffles2 `.snf` files (`snf.py`): SV candidates for multi-sample calling.
//!
//! An SNF file is one JSON header line followed by gzip-compressed pickles, one per 100 kbp block
//! of a contig. Each block holds the candidates (`sniffles.sv.SVCall` objects) by SV type and the
//! sample's mean coverage in 500 bp windows. Offsets in the header's index are relative to the end
//! of the header line.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use serde_json::{Value as Json, json};

use crate::call::{Alt, Coverage, Info, Phase, SvCall};
use crate::config::Config;
use crate::pickle::{self, Value};
use crate::reads::SvType;

/// Candidates are stored in blocks of this many reference bases.
pub const BLOCK_SIZE: i64 = 100_000;
/// Coverage is stored as the mean over windows of this many bases (`coverage_binsize_combine`).
pub const COVERAGE_WINDOW: i64 = 500;
/// `snf_format_version` of Sniffles2 2.8.1; files of other versions are rejected.
const FORMAT_VERSION: &str = "S2_rc4";
/// The Sniffles2 build this writer's files are compatible with.
const COMPATIBLE_BUILD: &str = "2.8.1";
/// Files written before this Sniffles2 build need their QC redone, which is not supported.
const MIN_BUILD: &str = "2.5.3";

/// The candidate's block (`int(pos / snf_block_size) * snf_block_size`).
pub fn block_of(pos: i64) -> i64 {
    (pos as f64 / BLOCK_SIZE as f64) as i64 * BLOCK_SIZE
}

fn s(v: &str) -> Value {
    Value::Str(v.into())
}

fn object(class: &str, state: Vec<(&str, Value)>) -> Value {
    Value::Object {
        class: Box::new(Value::Global("sniffles.sv".into(), class.into())),
        args: Vec::new(),
        state: Box::new(Value::Dict(state.into_iter().map(|(k, v)| (s(k), v)).collect())),
    }
}

fn info_value(v: &Info) -> Value {
    match v {
        Info::Int(i) => Value::Int(*i),
        Info::Float(f) => Value::Float(*f),
        Info::Str(t) => s(t),
        Info::Bool(b) => Value::Bool(*b),
        Info::None => Value::None,
    }
}

/// A phase block ID is the integer PS tag where it is one.
fn phase_value(v: Option<&str>) -> Value {
    match v {
        None => Value::None,
        Some(p) => p.parse().map_or_else(|_| s(p), Value::Int),
    }
}

/// A genotype's phase as Sniffles2 pickles it: (HP, PS).
pub fn phase_tuple(phase: &Phase) -> Value {
    Value::Tuple(vec![phase.hp.as_deref().map_or(Value::None, s), phase_value(phase.ps.as_deref())])
}

/// The `SVCall` object Sniffles2 pickles for a candidate (its dataclass fields, in order).
fn call_value(c: &SvCall, contig: &str, cfg: &Config) -> Value {
    let alt = match &c.alt {
        Alt::Seq(seq) => String::from_utf8_lossy(seq).into_owned(),
        Alt::Symbolic => format!("<{}>", c.svtype.as_str()),
    };
    let genotypes = match &c.genotype {
        Some(g) => vec![(
            Value::Int(0),
            Value::Tuple(vec![
                Value::Int(g.a),
                Value::Int(g.b),
                Value::Int(g.gq),
                Value::Int(g.dr),
                Value::Int(g.dv),
                phase_tuple(&g.phase),
            ]),
        )],
        None => Vec::new(),
    };
    let rnames = if cfg.args.output_rnames {
        Value::List(c.rnames.iter().map(|n| s(&String::from_utf8_lossy(n))).collect())
    } else {
        Value::None
    };
    let bnd_info = c.bnd.as_ref().map_or(Value::None, |b| {
        object(
            "SVCallBNDInfo",
            vec![
                ("mate_contig", s(&b.mate_contig)),
                ("mate_ref_start", Value::Int(b.mate_ref_start)),
                ("is_first", Value::Bool(b.is_first)),
                ("is_reverse", Value::Bool(b.is_reverse)),
            ],
        )
    });
    let sampler = object(
        "ForwardDifferenceWelford",
        vec![("n", Value::Int(0)), ("m1", Value::Int(0)), ("m2", Value::Int(0)), ("last", Value::None)],
    );
    object(
        "SVCall",
        vec![
            ("contig", s(contig)),
            ("pos", Value::Int(c.pos)),
            ("id", s(&c.id)),
            ("ref", s("N")),
            ("alt", Value::Str(alt.into())),
            ("qual", Value::Int(c.qual)),
            ("filter", s(c.filter)),
            ("info", Value::Dict(c.info.iter().map(|(k, v)| (s(k), info_value(v))).collect())),
            ("svtype", s(c.svtype.as_str())),
            ("svlen", Value::Int(c.svlen)),
            ("end", Value::Int(c.end)),
            ("genotypes", Value::Dict(genotypes)),
            ("precise", Value::Bool(c.precise)),
            ("support", Value::Int(c.support as i64)),
            ("rnames", rnames),
            ("qc", Value::Bool(c.qc)),
            ("nm", if cfg.qc_nm_measure { Value::Float(c.nm) } else { Value::Int(-1) }),
            ("postprocess", Value::None),
            (
                "svlens",
                c.svlens.as_ref().map_or(Value::None, |l| Value::List(l.iter().map(|&x| Value::Int(x)).collect())),
            ),
            ("fwd", Value::Int(c.fwd as i64)),
            ("rev", Value::Int(c.rev as i64)),
            ("forward_difference_sampler", sampler),
            ("coverage_upstream", Value::Int(c.cov_up)),
            ("coverage_downstream", Value::Int(c.cov_down)),
            ("coverage_start", Value::Int(c.cov_start)),
            ("coverage_center", Value::Int(c.cov_center)),
            ("coverage_end", Value::Int(c.cov_end)),
            ("sample_internal_id", Value::None),
            ("bnd_info", bnd_info),
            ("support_inline", Value::None),
            ("support_splits", Value::None),
            ("raw_vcf_line", Value::None),
            ("raw_vcf_line_index", Value::None),
            ("is_single_break", Value::Bool(false)),
        ],
    )
}

/// One contig's SNF blocks: (block start, gzip-compressed pickle), by block start.
#[derive(Debug, Default)]
pub struct ContigBlocks {
    pub blocks: Vec<(i64, Vec<u8>)>,
    pub candidate_count: usize,
    /// Mean coverage of the contig (`coverage_average_total`).
    pub mean_coverage: f64,
}

/// Packs a contig's candidates (all of them, QC failures included, in calling order) and the
/// positions of its single break candidates into blocks with the coverage of each block (`SNFile.store`, `annotate_block_coverages`, `write_and_index`).
///
/// # Errors
///
/// Fails if compression fails.
pub fn encode_contig(
    calls: &[SvCall],
    single_breaks: &[i64],
    contig: &str,
    coverage: &Coverage<'_>,
    cfg: &Config,
) -> Result<ContigBlocks> {
    let mut by_block: std::collections::BTreeMap<i64, [Vec<Value>; 5]> = std::collections::BTreeMap::new();
    for c in calls {
        by_block.entry(block_of(c.pos)).or_default()[c.svtype as usize].push(call_value(c, contig, cfg));
    }
    // Single break candidates are not stored, but create their block.
    for &pos in single_breaks {
        by_block.entry(block_of(pos)).or_default();
    }
    let windows = coverage.window_sums(COVERAGE_WINDOW);
    let mut blocks = Vec::with_capacity(by_block.len());
    for (start, lists) in by_block {
        let mut items: Vec<(Value, Value)> =
            SvType::ALL.into_iter().zip(lists).map(|(t, l)| (s(t.as_str()), Value::List(l))).collect();
        let first = (start / COVERAGE_WINDOW).max(0);
        let cov: Vec<(Value, Value)> = (first..first + BLOCK_SIZE / COVERAGE_WINDOW)
            .filter_map(|w| {
                let sum = *windows.get(usize::try_from(w).ok()?)?;
                let mean = crate::pymath::py_round(sum as f64 / COVERAGE_WINDOW as f64);
                Some((Value::Int(start + (w - first) * COVERAGE_WINDOW), Value::Int(mean)))
            })
            .collect();
        items.push((s("_COVERAGE"), Value::Dict(cov)));
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(&pickle::encode(&Value::Dict(items)))?;
        blocks.push((start, gz.finish()?));
    }
    Ok(ContigBlocks { blocks, candidate_count: calls.len() + single_breaks.len(), mean_coverage: coverage.mean() })
}

/// Writes an SNF file from the contigs' blocks (`SNFileBase.write_results`).
///
/// # Errors
///
/// Fails on I/O errors.
pub fn write(path: &str, cfg: &Config, contigs: &[(String, i64)], parts: &[(String, ContigBlocks)]) -> Result<()> {
    let mut index = serde_json::Map::new();
    let mut offset = 0usize;
    for (contig, part) in parts {
        let mut blocks = serde_json::Map::new();
        for (start, data) in &part.blocks {
            blocks.insert(start.to_string(), json!([[offset, data.len()]]));
            offset += data.len();
        }
        index.insert(contig.clone(), Json::Object(blocks));
    }
    let coverages: serde_json::Map<String, Json> = contigs
        .iter()
        .map(|(name, _)| {
            let cov = parts.iter().find(|(c, _)| c == name).map_or(0.0, |(_, p)| p.mean_coverage);
            (name.clone(), json!(cov))
        })
        .collect();
    let header = json!({
        "config": {
            "version": "lrpipe-sv",
            "build": COMPATIBLE_BUILD,
            "lrpipe_sv_version": env!("CARGO_PKG_VERSION"),
            "command": cfg.command,
            "snf_format_version": FORMAT_VERSION,
            "snf_block_size": BLOCK_SIZE,
            "sample_id": cfg.args.sample_id,
            "contig_lengths": contigs,
            "contig_coverages": coverages,
            "minsvlen": cfg.minsvlen,
            "output_rnames": cfg.args.output_rnames,
        },
        "index": index,
        "snf_candidate_count": parts.iter().map(|(_, p)| p.candidate_count).sum::<usize>(),
    });
    let mut out = std::io::BufWriter::new(File::create(path).with_context(|| format!("creating {path}"))?);
    serde_json::to_writer(&mut out, &header)?;
    out.write_all(b"\n")?;
    for (_, part) in parts {
        for (_, data) in &part.blocks {
            out.write_all(data)?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Writes `input` again as `output`, each candidate (a pickled `SVCall` object) changed by
/// `patch`: the same header (but for the blocks' places) and blocks, in the same order.
///
/// # Errors
///
/// If `input` cannot be read or is not an SNF file this writer wrote, `patch` fails, or `output`
/// cannot be written.
pub fn rewrite(input: &Path, output: &str, mut patch: impl FnMut(&mut Value) -> Result<()>) -> Result<()> {
    let context = || format!("reading {}", input.display());
    let data = std::fs::read(input).with_context(context)?;
    let header_end = data.iter().position(|&b| b == b'\n').with_context(context)? + 1;
    let header = std::str::from_utf8(&data[..header_end]).with_context(context)?;
    // The index sits between the config and the candidate count; nothing in the config holds
    // these keys' unescaped quotes.
    let (index_at, count_at) = match (header.rfind(",\"index\":"), header.rfind(",\"snf_candidate_count\":")) {
        (Some(i), Some(c)) if i < c => (i, c),
        _ => bail!("{} is not an SNF file lrpipe-sv wrote", input.display()),
    };
    let index: serde_json::Map<String, Json> =
        serde_json::from_str(&header[index_at + ",\"index\":".len()..count_at]).with_context(context)?;
    let blocks = &data[header_end..];
    let mut new_index = serde_json::Map::new();
    let mut out = Vec::with_capacity(blocks.len());
    for (contig, starts) in &index {
        let mut new_starts = serde_json::Map::new();
        for (start, parts) in starts.as_object().with_context(context)? {
            let mut new_parts = Vec::new();
            for part in parts.as_array().with_context(context)? {
                let (offset, len) = part
                    .as_array()
                    .and_then(|p| {
                        Some((usize::try_from(p.first()?.as_u64()?).ok()?, usize::try_from(p.get(1)?.as_u64()?).ok()?))
                    })
                    .with_context(context)?;
                let gz = blocks.get(offset..offset + len).with_context(context)?;
                let mut pickled = Vec::new();
                MultiGzDecoder::new(gz).read_to_end(&mut pickled).with_context(context)?;
                let mut block = pickle::decode(&pickled).with_context(context)?;
                if let Value::Dict(items) = &mut block {
                    for (_, list) in items.iter_mut().filter(|(k, _)| k.as_str().is_some_and(|k| k != "_COVERAGE")) {
                        if let Value::List(candidates) = list {
                            candidates.iter_mut().try_for_each(&mut patch)?;
                        }
                    }
                }
                let mut gz = GzEncoder::new(Vec::new(), Compression::default());
                gz.write_all(&pickle::encode(&block))?;
                let gz = gz.finish()?;
                new_parts.push(json!([out.len(), gz.len()]));
                out.extend_from_slice(&gz);
            }
            new_starts.insert(start.clone(), Json::Array(new_parts));
        }
        new_index.insert(contig.clone(), Json::Object(new_starts));
    }
    let mut file = std::io::BufWriter::new(File::create(output).with_context(|| format!("creating {output}"))?);
    file.write_all(&header.as_bytes()[..index_at])?;
    write!(file, ",\"index\":{}", Json::Object(new_index))?;
    file.write_all(&header.as_bytes()[count_at..])?;
    file.write_all(&out)?;
    file.flush()?;
    Ok(())
}

/// One sample genotype: `.` alleles are `None`.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleGenotype {
    pub a: Option<i64>,
    pub b: Option<i64>,
    pub gq: i64,
    pub dr: i64,
    pub dv: i64,
    pub phase: Phase,
}

/// A BND candidate's mate and its orientation on both sides (`SVCallBNDInfo`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnfBnd {
    pub mate_contig: String,
    pub mate_pos: i64,
    /// The local sequence comes before the mate (`N[p[`, `N]p]`).
    pub is_first: bool,
    /// The mate's side is `]p]`, not `[p[`.
    pub is_reverse: bool,
}

/// A candidate read from an SNF block (the `SVCall` fields multi-sample calling uses).
#[derive(Debug, Clone)]
pub struct Candidate {
    pub svtype: SvType,
    pub contig: String,
    pub pos: i64,
    pub end: i64,
    pub svlen: i64,
    pub id: String,
    pub alt: String,
    pub qual: Option<i64>,
    pub filter: String,
    pub info: Vec<(String, Info)>,
    pub genotype: Option<SampleGenotype>,
    pub precise: bool,
    pub support: i64,
    pub rnames: Option<Vec<String>>,
    pub qc: bool,
    pub fwd: i64,
    pub rev: i64,
    /// Coverage upstream, at start, center, end and downstream.
    pub coverage: [Option<i64>; 5],
    /// The mate and orientation, for BNDs.
    pub bnd: Option<SnfBnd>,
    /// The input file's index.
    pub sample: usize,
}

/// A decoded block: candidates by SV type and coverage by window start.
#[derive(Debug, Default)]
pub struct Block {
    pub candidates: [Vec<Candidate>; 5],
    pub coverage: HashMap<i64, i64>,
}

fn field<'a>(obj: &'a Value, name: &str) -> Result<&'a Value> {
    obj.attr(name).with_context(|| format!("SVCall without '{name}'"))
}

fn int_field(obj: &Value, name: &str) -> Result<i64> {
    field(obj, name)?.as_int().with_context(|| format!("SVCall '{name}' is not an integer"))
}

fn bool_field(obj: &Value, name: &str) -> Result<bool> {
    field(obj, name)?.as_bool().with_context(|| format!("SVCall '{name}' is not a bool"))
}

fn str_field(obj: &Value, name: &str) -> Result<String> {
    Ok(field(obj, name)?.as_str().with_context(|| format!("SVCall '{name}' is not a string"))?.to_string())
}

/// A genotype allele, `.` or an integer.
fn allele(v: &Value) -> Result<Option<i64>> {
    match v {
        Value::Str(t) if &**t == "." => Ok(None),
        _ => v.as_int().map(Some).context("invalid genotype allele"),
    }
}

/// A phase tuple entry as Sniffles2 formats it.
fn phase_text(v: &Value) -> Option<String> {
    match v {
        Value::None => None,
        Value::Str(t) => Some(t.to_string()),
        Value::Int(i) => Some(i.to_string()),
        Value::Bool(b) => Some(if *b { "True" } else { "False" }.to_string()),
        _ => None,
    }
}

fn genotype(v: &Value) -> Result<Option<SampleGenotype>> {
    let Value::Dict(items) = v else { bail!("SVCall genotypes is not a dict") };
    let Some((_, gt)) = items.iter().find(|(k, _)| k.as_int() == Some(0)) else { return Ok(None) };
    let Value::Tuple(f) = gt else { bail!("genotype is not a tuple") };
    if f.len() < 5 {
        bail!("genotype has {} fields", f.len());
    }
    let phase = match f.get(5) {
        Some(Value::Tuple(p)) if p.len() == 2 => Phase { hp: phase_text(&p[0]), ps: phase_text(&p[1]) },
        _ => Phase::default(),
    };
    let int = |k: usize| f[k].as_int().context("invalid genotype field");
    Ok(Some(SampleGenotype { a: allele(&f[0])?, b: allele(&f[1])?, gq: int(2)?, dr: int(3)?, dv: int(4)?, phase }))
}

fn info(v: &Value) -> Vec<(String, Info)> {
    let Value::Dict(items) = v else { return Vec::new() };
    items
        .iter()
        .filter_map(|(k, v)| {
            let value = match v {
                Value::Int(i) => Info::Int(*i),
                Value::Float(f) => Info::Float(*f),
                Value::Str(t) => Info::Str(t.to_string()),
                Value::Bool(b) => Info::Bool(*b),
                Value::None => Info::None,
                _ => return None,
            };
            Some((k.as_str()?.to_string(), value))
        })
        .collect()
}

/// The `SVCall` fields read, by position in [`candidate`]'s lookup.
const FIELDS: [&str; 23] = [
    "svtype",
    "contig",
    "pos",
    "end",
    "svlen",
    "id",
    "alt",
    "qual",
    "filter",
    "info",
    "genotypes",
    "precise",
    "support",
    "rnames",
    "qc",
    "fwd",
    "rev",
    "coverage_upstream",
    "coverage_start",
    "coverage_center",
    "coverage_end",
    "coverage_downstream",
    "bnd_info",
];

fn candidate(v: &Value, sample: usize) -> Result<Candidate> {
    let Value::Object { state, .. } = v else { bail!("SNF candidate is not an object") };
    let Value::Dict(items) = &**state else { bail!("SNF candidate without state") };
    // One pass over the object's attributes.
    let mut slots: [Option<&Value>; FIELDS.len()] = [None; FIELDS.len()];
    for (k, val) in items {
        if let Some(i) = k.as_str().and_then(|k| FIELDS.iter().position(|f| *f == k)) {
            slots[i] = Some(val);
        }
    }
    let get = |i: usize| slots[i].with_context(|| format!("SVCall without '{}'", FIELDS[i]));
    let int = |i: usize| get(i)?.as_int().with_context(|| format!("SVCall '{}' is not an integer", FIELDS[i]));
    let text = |i: usize| -> Result<String> {
        Ok(get(i)?.as_str().with_context(|| format!("SVCall '{}' is not a string", FIELDS[i]))?.to_string())
    };
    let flag = |i: usize| get(i)?.as_bool().with_context(|| format!("SVCall '{}' is not a bool", FIELDS[i]));
    let svtype_name = text(0)?;
    let svtype = SvType::ALL
        .into_iter()
        .find(|t| t.as_str() == svtype_name)
        .with_context(|| format!("unknown SV type {svtype_name}"))?;
    let bnd = match slots[22] {
        Some(b @ Value::Object { .. }) => Some(SnfBnd {
            mate_contig: str_field(b, "mate_contig")?,
            mate_pos: int_field(b, "mate_ref_start")?,
            is_first: bool_field(b, "is_first")?,
            is_reverse: bool_field(b, "is_reverse")?,
        }),
        _ => None,
    };
    let rnames = match slots[13] {
        Some(Value::List(l)) => Some(l.iter().filter_map(|x| x.as_str().map(str::to_string)).collect()),
        _ => None,
    };
    let opt = |i: usize| slots[i].and_then(Value::as_int);
    Ok(Candidate {
        svtype,
        contig: text(1)?,
        pos: int(2)?,
        end: int(3)?,
        svlen: int(4)?,
        id: text(5)?,
        alt: text(6)?,
        qual: match get(7)? {
            Value::None => None,
            q => Some(q.as_int().or_else(|| q.as_float().map(|f| f as i64)).context("invalid SVCall qual")?),
        },
        filter: text(8)?,
        info: info(get(9)?),
        genotype: genotype(get(10)?)?,
        precise: flag(11)?,
        support: int(12)?,
        rnames,
        qc: flag(14)?,
        fwd: opt(15).unwrap_or(0),
        rev: opt(16).unwrap_or(0),
        coverage: [opt(17), opt(18), opt(19), opt(20), opt(21)],
        bnd,
        sample,
    })
}

/// Decodes a block's pickle.
fn decode_block(data: &[u8], sample: usize) -> Result<Block> {
    let v = pickle::decode(data)?;
    let mut block = Block::default();
    for t in SvType::ALL {
        if let Some(Value::List(items)) = v.get(t.as_str()) {
            block.candidates[t as usize] = items.iter().map(|c| candidate(c, sample)).collect::<Result<_>>()?;
        }
    }
    if let Some(Value::Dict(items)) = v.get("_COVERAGE") {
        block.coverage = items.iter().filter_map(|(k, v)| Some((k.as_int()?, v.as_int()?))).collect();
    }
    Ok(block)
}

/// An SNF file opened for reading.
#[derive(Debug)]
pub struct SnfFile {
    pub path: String,
    /// The file, kept open; `None` when each block read opens it again
    /// (`--combine-close-handles`), so that any number of inputs stays under the open-file limit.
    file: Option<File>,
    header_len: u64,
    pub sample_id: Option<String>,
    pub contig_lengths: Vec<(String, i64)>,
    pub candidate_count: u64,
    /// Contig -> block start -> (offset, length) of each part.
    index: HashMap<String, HashMap<i64, Vec<(u64, u64)>>>,
}

impl SnfFile {
    /// Opens an SNF file and reads its header; the file stays open only with `keep_open`.
    ///
    /// # Errors
    ///
    /// Fails on unreadable or invalid files, other block sizes and format versions, and files
    /// written by Sniffles2 before 2.5.3.
    pub fn open(path: &str, keep_open: bool) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {path}"))?;
        let mut line = String::new();
        BufReader::new(&file).read_line(&mut line).with_context(|| format!("reading the SNF header of {path}"))?;
        let header_len = line.len() as u64;
        // Python's json module writes NaN and Infinity, which JSON does not allow.
        let text =
            line.trim().replace(": NaN", ": null").replace(": Infinity", ": null").replace(": -Infinity", ": null");
        let h: Json = serde_json::from_str(&text)
            .with_context(|| format!("{path}: the file may not be a valid .snf file or could have been corrupted"))?;
        let config = &h["config"];
        if config["snf_block_size"].as_i64() != Some(BLOCK_SIZE) {
            bail!("SNF block size differs for {path}");
        }
        if config["snf_format_version"].as_str() != Some(FORMAT_VERSION) {
            bail!("SNF format version for {path} is not supported");
        }
        let build = config["build"].as_str().unwrap_or_default();
        if build.split('-').next().unwrap_or_default() < MIN_BUILD {
            bail!("{path} was written by Sniffles2 {build}; files from builds before {MIN_BUILD} are not supported");
        }
        let contig_lengths = config["contig_lengths"]
            .as_array()
            .context("SNF header without contig_lengths")?
            .iter()
            .map(|c| Some((c.get(0)?.as_str()?.to_string(), c.get(1)?.as_i64()?)))
            .collect::<Option<Vec<_>>>()
            .context("invalid contig_lengths in the SNF header")?;
        let mut index: HashMap<String, HashMap<i64, Vec<(u64, u64)>>> = HashMap::new();
        for (contig, blocks) in h["index"].as_object().context("SNF header without index")? {
            let entry = index.entry(contig.clone()).or_default();
            for (start, parts) in blocks.as_object().context("invalid SNF index")? {
                let parts = parts
                    .as_array()
                    .context("invalid SNF index")?
                    .iter()
                    .map(|p| Some((p.get(0)?.as_u64()?, p.get(1)?.as_u64()?)))
                    .collect::<Option<Vec<_>>>()
                    .context("invalid SNF index")?;
                entry.insert(start.parse().context("invalid SNF block start")?, parts);
            }
        }
        Ok(SnfFile {
            path: path.to_string(),
            file: keep_open.then_some(file),
            header_len,
            sample_id: config["sample_id"].as_str().map(str::to_string),
            contig_lengths,
            candidate_count: h["snf_candidate_count"].as_u64().unwrap_or(0),
            index,
        })
    }

    /// The parts of one block (`read_blocks`), or `None` if the file has none there.
    ///
    /// # Errors
    ///
    /// Fails on unreadable or corrupt blocks.
    pub fn read_blocks(&self, contig: &str, start: i64, sample: usize) -> Result<Option<Vec<Block>>> {
        let Some(parts) = self.index.get(contig).and_then(|b| b.get(&start)) else { return Ok(None) };
        let reopened;
        let file = match &self.file {
            Some(f) => f,
            None => {
                reopened = File::open(&self.path).with_context(|| format!("opening {}", self.path))?;
                &reopened
            }
        };
        let mut blocks = Vec::with_capacity(parts.len());
        for &(offset, len) in parts {
            let block = self.read_block(file, offset, len, sample).with_context(|| {
                format!("reading block '{contig}.{start}' from '{}' (truncated or corrupt file?)", self.path)
            })?;
            blocks.push(block);
        }
        Ok(Some(blocks))
    }

    fn read_block(&self, file: &File, offset: u64, len: u64, sample: usize) -> Result<Block> {
        let mut gz = vec![0u8; usize::try_from(len)?];
        file.read_exact_at(&mut gz, self.header_len + offset)?;
        let mut data = Vec::new();
        MultiGzDecoder::new(&gz[..]).read_to_end(&mut data)?;
        decode_block(&data, sample)
    }
}

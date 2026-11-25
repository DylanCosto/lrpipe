//! SVs phased after calling, from a haplotag list: the read pass calls them as if the reads had
//! no HP and PS tags and writes what phasing them needs (`--phasing-evidence`), and `lrpipe sv
//! phase` gives them the phases Sniffles2 2.8.1 gives them from reads tagged as the list says.
//!
//! The reads' haplotypes change three things in a single-sample call, all from its reads
//! alone: its PHASE INFO and the genotype's phase (`postprocessing.phase_sv` and the phase of a
//! 1/1 genotype in `annotate_sv`), and whether a MOSAIC_VAF call is re-passed as heterozygous
//! (`Task.rescue_phasing`), which also counts the haplotypes of the other leads and reads at its
//! first bin. The evidence holds, for each call in calling order, the reads behind each of
//! these (a primary alignment by its read's name, any other alignment as one a haplotag list
//! never tags), the genotype, and for a call the rescue may re-pass its VCF record. The phases
//! are then worked out by the same functions the read pass uses ([`call::phase_of`],
//! [`call::phase_homozygous`], [`call::rescued_by_phase`]), and the read pass's VCF and SNF
//! rewritten with them: the PHASE values, genotypes and PS, and the filter, QC and genotype of
//! the calls re-passed. So the output is the one `lrpipe sv` writes from the alignments
//! `haplotag` tags with the list, which is Sniffles2's.
//!
//! An evidence file is gzipped text. Its first line is [`MAGIC`], its second the options that
//! matter (`threshold=`, `sort=`, `no_qc=`), then each contig's lines: `C` and the contig, `R`
//! and a read name (the contig's reads are numbered from 0 in that order), and one `V` line per
//! call with its ID, position, QC (0 or 1), whether the N-region check filtered it, whether it
//! is a BND, its genotype (`a,b` or `.`), its alignments (read numbers, `-` for an alignment no
//! list tags), what the rescue counts (`-`, or the read numbers of the leads at its bin and of
//! the reads over it, as `leads;reads`) and its VCF record (empty if the rescue cannot pass it,
//! or when the record is skipped).

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{BufRead, BufReader, Read as _, Write as _};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, bail};
use clap::Parser;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use lrpipe_reads::haplotags::Haplotags;

use crate::call::{self, Genotype, Phase, ReadPhase, SvCall};
use crate::cluster::Leadtab;
use crate::config::Config;
use crate::pickle::Value;
use crate::reads::SvType;
use crate::vcf;

/// The first line of an evidence file.
pub const MAGIC: &str = "#lrpipe-sv phasing evidence 1";

/// An HP or PS tag's value as Sniffles2 has it when the tag is absent.
const NO_HAPLOTYPE: &str = "0";
const NO_PHASE_SET: &str = "NULL";

/// A gzip member of `text`.
fn gzip(text: &str) -> Result<Vec<u8>> {
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(text.as_bytes())?;
    Ok(gz.finish()?)
}

/// The evidence file's first lines, as a gzip member.
///
/// # Errors
///
/// If compression fails.
pub fn header(cfg: &Config) -> Result<Vec<u8>> {
    gzip(&format!(
        "{MAGIC}\nthreshold={:?}\tsort={}\tno_qc={}\n",
        cfg.args.phase_conflict_threshold,
        u8::from(cfg.sort),
        u8::from(cfg.args.no_qc)
    ))
}

/// The reads one contig's evidence names, numbered in the order first named.
#[derive(Default)]
struct Names {
    numbers: HashMap<Arc<[u8]>, usize>,
    text: String,
}

impl Names {
    /// `name`'s number, naming it on first use.
    fn number(&mut self, name: &Arc<[u8]>) -> usize {
        let next = self.numbers.len();
        *self.numbers.entry(Arc::clone(name)).or_insert_with(|| {
            let _ = writeln!(self.text, "R\t{}", String::from_utf8_lossy(name));
            next
        })
    }
}

/// The phasing evidence of one contig's calls (all of them, QC failures included, in calling
/// order, after QC and the N-region check), as a gzip member.
///
/// `hap_bins` and `hap_reads` are the contig's alignments' bins with their read names
/// (`ChunkLeads::hap_reads`); `reference` gives the VCF records.
///
/// # Errors
///
/// If the reference cannot be read or compression fails.
pub fn contig_evidence(
    contig: &str,
    calls: &[SvCall],
    tabs: &[Leadtab; 5],
    hap_bins: &[(i64, i64, u8)],
    hap_reads: &[Option<Arc<[u8]>>],
    cfg: &Config,
    reference: vcf::Reference<'_>,
) -> Result<Vec<u8>> {
    let mut names = Names::default();
    // The reads over each bin the rescue counts at.
    let mut bins: Vec<i64> = calls.iter().filter(|c| c.rescuable).map(|c| c.cluster.start).collect();
    bins.sort_unstable();
    bins.dedup();
    let mut over: Vec<Vec<usize>> = vec![Vec::new(); bins.len()];
    for (&(first, end, _), name) in hap_bins.iter().zip(hap_reads) {
        // `HapRef` counts the reads whose bins [first, end) hold the bin, once they span one.
        let Some(name) = name.as_ref().filter(|_| end > first) else { continue };
        let from = bins.partition_point(|&b| b < first);
        let to = bins.partition_point(|&b| b < end);
        if from < to {
            let n = names.number(name);
            for reads in &mut over[from..to] {
                reads.push(n);
            }
        }
    }
    let mut lines = String::new();
    for c in calls {
        let tab = &tabs[c.svtype as usize];
        // One per alignment, as `phase_sv` takes them.
        let mut seen = std::collections::HashSet::new();
        let mut reads = Vec::new();
        for &i in &c.cluster.leads {
            let lead = &tab.leads[i];
            if seen.insert(lead.read_id) {
                reads.push(lead_token(lead, &mut names));
            }
        }
        let rescue = if c.rescuable {
            let at = tab.bin(c.cluster.start).unwrap_or_default();
            let leads: Vec<String> =
                at.iter().map(|&i| lead_token(&tab.leads[i], &mut names)).filter(|t| t != "-").collect();
            let k = bins.binary_search(&c.cluster.start).expect("a rescuable call's bin");
            let reads: Vec<String> = over[k].iter().map(usize::to_string).collect();
            format!("{};{}", leads.join(","), reads.join(","))
        } else {
            "-".to_owned()
        };
        let mut record = String::new();
        if c.rescuable {
            vcf::write_call(&mut record, c, contig, cfg, reference)?;
        }
        let genotype = c.genotype.as_ref().map_or_else(|| ".".to_owned(), |g| format!("{},{}", g.a, g.b));
        let _ = writeln!(
            lines,
            "V\t{}\t{}\t{}\t{}\t{}\t{genotype}\t{}\t{rescue}\t{}",
            c.id,
            c.pos,
            u8::from(c.qc),
            u8::from(c.filter == "UNRESOLVED"),
            u8::from(c.svtype == SvType::Bnd),
            reads.join(","),
            record.trim_end_matches('\n'),
        );
    }
    gzip(&format!("C\t{contig}\n{}{lines}", names.text))
}

/// A lead's alignment in the evidence: its read's number when a haplotag list tags it (the
/// primary alignment's leads but BND ones, which Sniffles2 never phases), else `-`.
fn lead_token(lead: &crate::reads::Lead, names: &mut Names) -> String {
    if lead.is_sa || lead.bnd.is_some() { "-".to_owned() } else { names.number(&lead.qname).to_string() }
}

/// One call of the evidence.
#[derive(Debug)]
struct Evidence {
    id: String,
    pos: i64,
    qc: bool,
    unresolved: bool,
    bnd: bool,
    genotype: Option<(i64, i64)>,
    reads: Vec<Option<usize>>,
    /// The read numbers of the leads at the call's bin and of the reads over it.
    rescue: Option<(Vec<usize>, Vec<usize>)>,
    record: String,
}

/// One contig of the evidence.
#[derive(Debug, Default)]
struct ContigEvidence {
    name: String,
    reads: Vec<Vec<u8>>,
    calls: Vec<Evidence>,
}

/// The options the evidence was written with.
#[derive(Debug)]
struct Written {
    threshold: f64,
    sort: bool,
    no_qc: bool,
}

/// Reads an evidence file.
fn read_evidence(path: &Path) -> Result<(Written, Vec<ContigEvidence>)> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(MultiGzDecoder::new(BufReader::new(file)));
    let mut lines = reader.lines().enumerate();
    let mut next = || -> Result<Option<(usize, String)>> {
        let line = lines.next().map(|(i, l)| l.map(|l| (i + 1, l))).transpose();
        line.with_context(|| format!("reading {}", path.display()))
    };
    let bad = |line: usize, what: &str| anyhow::anyhow!("{}, line {line}: {what}", path.display());
    if next()?.map(|(_, l)| l).as_deref() != Some(MAGIC) {
        bail!("{} is not lrpipe-sv's phasing evidence (--phasing-evidence)", path.display());
    }
    let (n, options) = next()?.ok_or_else(|| bad(2, "no options"))?;
    let option =
        |key: &str| options.split('\t').find_map(|o| o.strip_prefix(key)?.strip_prefix('=')).ok_or_else(|| bad(n, key));
    let written = Written {
        threshold: option("threshold")?.parse().map_err(|_| bad(n, "threshold"))?,
        sort: option("sort")? == "1",
        no_qc: option("no_qc")? == "1",
    };
    let mut contigs: Vec<ContigEvidence> = Vec::new();
    while let Some((n, line)) = next()? {
        let fields: Vec<&str> = line.splitn(10, '\t').collect();
        match fields[..] {
            ["C", name] => contigs.push(ContigEvidence { name: name.to_owned(), ..ContigEvidence::default() }),
            ["R", name] => contigs.last_mut().ok_or_else(|| bad(n, "a read before a contig"))?.reads.push(name.into()),
            ["V", id, pos, qc, unresolved, bnd, genotype, reads, rescue, record] => {
                let contig = contigs.last_mut().ok_or_else(|| bad(n, "a call before a contig"))?;
                let number = |t: &str| -> Result<usize> {
                    let k: usize = t.parse().map_err(|_| bad(n, "a read number"))?;
                    if k >= contig.reads.len() {
                        return Err(bad(n, "an unnamed read"));
                    }
                    Ok(k)
                };
                let numbers = |list: &str| -> Result<Vec<usize>> {
                    list.split(',').filter(|t| !t.is_empty()).map(number).collect()
                };
                let genotype = match genotype.split_once(',') {
                    None => None,
                    Some((a, b)) => Some((
                        a.parse().map_err(|_| bad(n, "a genotype"))?,
                        b.parse().map_err(|_| bad(n, "a genotype"))?,
                    )),
                };
                let reads = reads
                    .split(',')
                    .filter(|t| !t.is_empty())
                    .map(|t| if t == "-" { Ok(None) } else { number(t).map(Some) })
                    .collect::<Result<_>>()?;
                let rescue = match rescue.split_once(';') {
                    None => None,
                    Some((leads, over)) => Some((numbers(leads)?, numbers(over)?)),
                };
                contig.calls.push(Evidence {
                    id: id.to_owned(),
                    pos: pos.parse().map_err(|_| bad(n, "a position"))?,
                    qc: qc == "1",
                    unresolved: unresolved == "1",
                    bnd: bnd == "1",
                    genotype,
                    reads,
                    rescue,
                    record: record.to_owned(),
                });
            }
            _ => return Err(bad(n, "not a C, R or V line")),
        }
    }
    Ok((written, contigs))
}

/// A call's phasing.
#[derive(Debug, Clone, PartialEq)]
struct Phased {
    /// Its PHASE INFO value.
    info: String,
    /// Its genotype's phase.
    phase: Phase,
    /// Whether the phasing rescue re-passed it: its filter PASS (unless the N-region check
    /// filters it), its QC passed and its genotype's second allele 1.
    rescued: bool,
}

/// `call`'s phasing with the haplotypes `list` gives `contig`'s reads.
fn phase(call: &Evidence, contig: &ContigEvidence, list: &Haplotags, threshold: f64) -> Result<Phased> {
    let listed = list.contig(&contig.name);
    let haplotype = |read: Option<usize>| -> Result<(String, String)> {
        let assigned = read.and_then(|k| listed?.get(&contig.reads[k]));
        Ok(match assigned {
            None => (NO_HAPLOTYPE.to_owned(), NO_PHASE_SET.to_owned()),
            Some(a) if a.haplotype > 2 => {
                bail!("haplotype H{} is not supported (Sniffles2 expects HP 0, 1 or 2)", a.haplotype)
            }
            Some(a) => (a.haplotype.to_string(), a.phase_set.to_string()),
        })
    };
    let reads: Vec<ReadPhase> = call
        .reads
        .iter()
        .map(|&r| {
            // A BND's leads have no phase set.
            haplotype(r).map(|(hp, ps)| (hp, (!call.bnd).then_some(ps)))
        })
        .collect::<Result<_>>()?;
    let (mut info, mut phase) = call::phase_of(&reads, threshold);
    if let Some((a, b)) = call.genotype {
        let mut g = Genotype { a, b, gq: 0, dr: 0, dv: 0, phase: phase.clone() };
        if let Some(homozygous) = call::phase_homozygous(&mut g, &info) {
            info = homozygous;
            phase = g.phase;
        }
    }
    let rescued = match &call.rescue {
        None => false,
        Some((leads, over)) => {
            let mut counts = [0i64; 6];
            for (k, reads) in [(0, leads), (3, over)] {
                for &r in reads {
                    let (hp, _) = haplotype(Some(r))?;
                    let slot = &mut counts[k + usize::from(hp.as_bytes()[0] - b'0')];
                    // Sniffles2 counts in 16-bit arrays.
                    *slot = (*slot + 1).min(i64::from(u16::MAX));
                }
            }
            call::rescued_by_phase(&info, counts)
        }
    };
    Ok(Phased { info, phase, rescued })
}

#[derive(Debug, Parser)]
#[command(
    name = "lrpipe-sv phase",
    about = "Phase the SVs of a run with --phasing-evidence from a haplotag list: write its VCF (and \
             SNF) as from reads haplotagged with the list"
)]
struct PhaseArgs {
    /// The run's --phasing-evidence
    #[arg(long, value_name = "FILE")]
    evidence: PathBuf,
    /// Haplotag list (WhatsHap's `haplotag --output-haplotag-list`): a listed read's primary
    /// alignment has its haplotype and phase set, every other alignment none
    #[arg(long, value_name = "LIST")]
    haplotags: PathBuf,
    /// The run's VCF (--vcf)
    #[arg(long, value_name = "IN.vcf")]
    input_vcf: PathBuf,
    /// The phased VCF (.gz for bgzip, with a tabix index)
    #[arg(long, value_name = "OUT.vcf")]
    vcf: String,
    /// The run's SNF (--snf)
    #[arg(long, value_name = "IN.snf", requires = "snf")]
    input_snf: Option<PathBuf>,
    /// The phased SNF
    #[arg(long, value_name = "OUT.snf", requires = "input_snf")]
    snf: Option<String>,
}

/// `lrpipe sv phase`: `args` are the program name, `phase` and its options.
///
/// # Errors
///
/// If an input cannot be read or does not fit the others, or an output cannot be written.
pub fn main(args: Vec<std::ffi::OsString>) -> Result<()> {
    let args = PhaseArgs::parse_from(args.into_iter().skip(1));
    let started = std::time::Instant::now();
    let list = Haplotags::read(&args.haplotags)?;
    let (written, contigs) = read_evidence(&args.evidence)?;
    // The calls' reads the list names: a list of other reads names none.
    for contig in &contigs {
        let listed = list.contig(&contig.name);
        for name in &contig.reads {
            list.matches().count(listed.is_some_and(|l| l.lists(name)));
        }
    }
    if let Some(report) = list.matches().report(&args.haplotags) {
        eprintln!("lrpipe-sv phase: {report} (the calls' reads)");
    }
    if list.ambiguous() > 0 {
        eprintln!(
            "lrpipe-sv phase: warning: {} reads listed more than once on a contig with different haplotypes (a list \
             written with --tag-supplementary?) are left untagged",
            list.ambiguous()
        );
    }
    let mut phased: HashMap<&str, (Phased, &Evidence)> = HashMap::new();
    let mut n_phased = 0;
    let mut n_rescued = 0;
    for contig in &contigs {
        for call in &contig.calls {
            let p = phase(call, contig, &list, written.threshold)
                .with_context(|| format!("phasing {} on {}", call.id, contig.name))?;
            n_phased += usize::from(p.phase.hp.is_some());
            n_rescued += usize::from(p.rescued);
            if phased.insert(&call.id, (p, call)).is_some() {
                bail!("{} names call {} twice", args.evidence.display(), call.id);
            }
        }
    }
    write_vcf(&args, &written, &contigs, &phased)?;
    if let (Some(input), Some(output)) = (&args.input_snf, &args.snf) {
        crate::snf::rewrite(input, output, |c| patch_candidate(c, &phased))?;
    }
    eprintln!(
        "Phased {n_phased} of {} SVs with {} reads' haplotypes from {}, {n_rescued} re-passed ({:.2}s)",
        phased.len(),
        list.len(),
        args.haplotags.display(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// The input VCF's header and records by ID (without `Sniffles2.`).
fn read_vcf(path: &Path) -> Result<(String, HashMap<String, String>)> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut text = String::new();
    let mut reader = BufReader::new(file);
    if reader.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        MultiGzDecoder::new(reader).read_to_string(&mut text)?;
    } else {
        reader.read_to_string(&mut text)?;
    }
    let mut header = String::new();
    let mut records = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            header.push_str(line);
            header.push('\n');
        } else if let Some(id) = line.split('\t').nth(2) {
            let id = id.strip_prefix("Sniffles2.").unwrap_or(id).to_owned();
            records.insert(id, line.to_owned());
        }
    }
    Ok((header, records))
}

/// The phased VCF: the calls that pass QC once phased (all with `--no-qc`), in the input's
/// order (by position with `sort`), each record with its PHASE, genotype and, once re-passed,
/// filter.
fn write_vcf(
    args: &PhaseArgs,
    written: &Written,
    contigs: &[ContigEvidence],
    phased: &HashMap<&str, (Phased, &Evidence)>,
) -> Result<()> {
    let (header, mut records) = read_vcf(&args.input_vcf)?;
    let mut out = Vec::with_capacity(contigs.len());
    for contig in contigs {
        let mut calls: Vec<&Evidence> = contig
            .calls
            .iter()
            .filter(|c| written.no_qc || c.qc || (phased[c.id.as_str()].0.rescued && !c.unresolved))
            .collect();
        if written.sort {
            calls.sort_by_key(|c| c.pos);
        }
        let mut text = String::new();
        for call in calls {
            let p = &phased[call.id.as_str()].0;
            let record = match records.remove(&call.id) {
                Some(r) => r,
                None if p.rescued && !call.record.is_empty() => call.record.clone(),
                // A call whose record Sniffles2 skips.
                None => continue,
            };
            text.push_str(&patch_record(&record, call, p).with_context(|| format!("phasing {}", call.id))?);
            text.push('\n');
        }
        out.push(text);
    }
    if let Some(id) = records.keys().next() {
        bail!("{} has call {id}, which {} does not", args.input_vcf.display(), args.evidence.display());
    }
    crate::cli::write_vcf(&args.vcf, &header, &out)
}

/// `record` with `p`'s PHASE, genotype phase and, once re-passed, filter and second allele.
fn patch_record(record: &str, call: &Evidence, p: &Phased) -> Result<String> {
    let mut columns: Vec<String> = record.split('\t').map(str::to_owned).collect();
    if columns.len() != 10 {
        bail!("a VCF record of {} columns, not 10", columns.len());
    }
    if p.rescued && !call.unresolved {
        columns[6] = "PASS".to_owned();
    }
    columns[7] = columns[7]
        .split(';')
        .map(|kv| if kv.starts_with("PHASE=") { format!("PHASE={}", p.info) } else { kv.to_owned() })
        .collect::<Vec<_>>()
        .join(";");
    if let Some((a, b)) = call.genotype {
        // GT:GQ:DR:DV:PS, unphased as the input's reads have no haplotypes.
        let f: Vec<&str> = columns[9].split(':').collect();
        let int = |s: &str| s.parse::<i64>().context("a genotype field is not an integer");
        let [_, gq, dr, dv, _] = f[..] else { bail!("a genotype of {} fields, not 5", f.len()) };
        let g = Genotype {
            a,
            b: if p.rescued { 1 } else { b },
            gq: int(gq)?,
            dr: int(dr)?,
            dv: int(dv)?,
            phase: p.phase.clone(),
        };
        let mut sample = String::new();
        vcf::format_genotype(&mut sample, Some(&g));
        columns[9] = sample;
    }
    Ok(columns.join("\t"))
}

/// An SNF candidate (a pickled `SVCall`) with its phasing.
fn patch_candidate(candidate: &mut Value, phased: &HashMap<&str, (Phased, &Evidence)>) -> Result<()> {
    let id = candidate.attr("id").and_then(Value::as_str).context("an SNF candidate without an ID")?;
    let (p, call) = phased.get(id).with_context(|| format!("SNF candidate {id} is not in the evidence"))?;
    let state = candidate.state_mut().context("an SNF candidate without state")?;
    if p.rescued {
        let qc = !call.unresolved;
        if qc {
            state.set("filter", Value::Str("PASS".into()));
        }
        state.set("qc", Value::Bool(qc));
    }
    if let Some(info) = state.get_mut("info") {
        info.set("PHASE", Value::Str(p.info.as_str().into()));
    }
    if let Some(Value::Dict(genotypes)) = state.get_mut("genotypes")
        && let Some((_, Value::Tuple(g))) = genotypes.first_mut()
        && g.len() == 6
    {
        if p.rescued {
            g[1] = Value::Int(1);
        }
        g[5] = crate::snf::phase_tuple(&p.phase);
    }
    Ok(())
}

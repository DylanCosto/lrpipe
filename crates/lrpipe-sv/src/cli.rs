//! The command line: argument handling and the run over contigs.

use std::collections::HashMap;
use std::ffi::OsString;
use std::io::{BufRead as _, Write as _};
use std::path::Path;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use anyhow::{Context, Result, bail};
use clap::Parser;
use lrpipe_reads::progress::Progress;
use lrpipe_reads::regions::Regions;
use pasta::variant::index::IndexFormat;
use pasta::variant::index::tabix::{Tabix, TabixBuilder, TabixConfig};
use rayon::prelude::*;

use crate::call::{ContigCtx, Coverage, SvCall};
use crate::cluster::{HapRef, Leadtab};
use crate::config::{Args, Config, DepthConfig, MinSupport};
use crate::fasta::{Fasta, NRuns};
use crate::methyl::{self, Methylation};
use crate::read_stats::ReadStats;
use crate::reads::{Alignments, ChunkLeads, SvType};
use crate::{alleles, call, cluster, combine, config, genotype, phasing, pymath, reads, snf, vcf};

/// Reads are fetched in pieces of this many reference bases, processed in parallel.
const CHUNK_LEN: i64 = 2_000_000;

/// Reading the contigs' share of the pass, for its progress; writing the calls takes the rest.
const READ_SHARE: f64 = 0.95;

/// Stack size of worker threads.
const WORKER_STACK: usize = 16 << 20;

/// Contigs shorter than this are skipped unless requested (`--all-contigs`, `--contig`).
pub const MIN_CONTIG_LEN: i64 = 1_000_000;

/// One contig to call (a Sniffles2 task; Sniffles2 makes one per contig).
#[derive(Debug)]
struct Task {
    id: usize,
    tid: usize,
    name: String,
    len: i64,
    /// `--regions` on this contig; empty means the whole contig.
    regions: Vec<(i64, i64)>,
    /// `--genotype-vcf` SVs on this contig, in input order.
    genotype_svs: Vec<genotype::InputSv>,
}

impl Task {
    /// The Sniffles2 regions to read: the given ones, or [0, len - 1) for the whole contig.
    fn read_regions(&self) -> Vec<(i64, i64)> {
        if self.regions.is_empty() { vec![(0, self.len - 1)] } else { self.regions.clone() }
    }
}

/// Runs the command line `args` (the program name first), as the `lrpipe-sv` binary does.
///
/// # Errors
///
/// If the arguments are invalid or the run fails.
pub fn main(args: Vec<OsString>) -> Result<()> {
    // `lrpipe-sv phase`: Sniffles2's command line takes no words but its options' values.
    if args.get(1).is_some_and(|a| a == "phase") {
        return crate::phasing::main(args);
    }
    let command = args.iter().map(|a| a.to_string_lossy()).collect::<Vec<_>>().join(" ");
    let cfg = Config::from_args(Args::parse_from(args), command)?;
    // CRAM decoding (pasta) needs more than the default 2 MB of stack per thread.
    rayon::ThreadPoolBuilder::new().num_threads(cfg.args.threads.max(1)).stack_size(WORKER_STACK).build_global()?;
    run(&cfg)
}

fn run(cfg: &Config) -> Result<()> {
    let start_time = Instant::now();
    let exts: Vec<String> = cfg.args.input.iter().map(|f| combine::extension(f)).collect();
    if exts.iter().any(|e| *e != exts[0]) {
        bail!("Please specify either a single .bam/.cram file, one or more .snf files, or a single .tsv file as input");
    }
    match exts[0].as_str() {
        "snf" | "tsv" if cfg.args.haplotags.is_some() || cfg.args.phasing_evidence.is_some() => bail!(
            "--haplotags and --phasing-evidence phase the SVs of a BAM or CRAM, not multi-sample calling from .snf \
             files (phase each sample's .snf with lrpipe sv phase first)"
        ),
        "snf" | "tsv" => return run_combine(cfg, start_time),
        "bam" | "cram" if cfg.args.input.len() == 1 => {}
        "bam" | "cram" => bail!("Please specify max 1 .bam/.cram file as input (got {})", cfg.args.input.len()),
        _ => bail!(
            "Failed to determine run mode from input. Please specify either: a single .bam/.cram file - OR - one \
             or more .snf files - OR - a single .tsv file containing a list of .snf files and optional sample ids"
        ),
    }
    if cfg.args.genotype_vcf.is_some() && cfg.args.snf.is_some() {
        bail!("--snf cannot be used with run mode genotype_vcf");
    }
    // Without SV output the pass may still count depth, repeats, methylation or the reads'
    // statistics (lrpipe only).
    if !cfg.calls_svs()
        && cfg.depth.is_none()
        && cfg.repeats.is_none()
        && cfg.args.bedmethyl.is_none()
        && cfg.args.read_stats.is_none()
    {
        bail!("Please specify at least one of: --vcf or --snf for output (both may be used at the same time)");
    }
    let depth_summary = cfg.depth.as_ref().map(|d| format!("{}.mosdepth.summary.txt", d.prefix));
    let outputs =
        [cfg.args.vcf.as_deref(), cfg.args.snf.as_deref(), depth_summary.as_deref(), cfg.args.read_stats.as_deref()];
    for path in outputs.into_iter().flatten() {
        if Path::new(path).exists() && !cfg.args.allow_overwrite {
            bail!("Output file '{path}' already exists! Use --allow-overwrite to ignore this check and overwrite.");
        }
    }
    let input = &cfg.args.input[0];
    let mut bam = Alignments::open(input, cfg.args.reference.as_deref())?;
    if cfg.args.bedmethyl.is_some() {
        bam.decode_tags(&lrpipe_methyl::input::TAGS);
    }
    let bam = bam;
    let reference = cfg.args.reference.as_deref().map(Fasta::open).transpose()?;
    let tandem_repeats =
        cfg.args.tandem_repeats.as_deref().map(|p| load_tandem_repeats(p, cfg.tandem_repeat_region_pad)).transpose()?;

    let genotype_input = cfg.args.genotype_vcf.as_deref().map(genotype::read_input).transpose()?;
    // `util.should_process_contig`
    let contigs: Vec<(usize, String, i64)> = bam
        .contigs
        .iter()
        .enumerate()
        .filter(|(_, (name, len))| {
            let requested = cfg.args.contig.iter().any(|c| c == name);
            let in_regions = cfg.regions_by_contig.contains_key(name);
            (cfg.args.contig.is_empty() || requested)
                && (cfg.regions_by_contig.is_empty() || in_regions)
                && (cfg.args.all_contigs || *len >= MIN_CONTIG_LEN || requested || in_regions)
        })
        .map(|(tid, (name, len))| (tid, name.clone(), *len))
        .collect();
    let tasks: Vec<Task> = contigs
        .iter()
        .filter(|(_, _, len)| *len > 1)
        .enumerate()
        .map(|(id, (tid, name, len))| Task {
            id,
            tid: *tid,
            name: name.clone(),
            len: *len,
            regions: cfg.regions_by_contig.get(name).cloned().unwrap_or_default(),
            // Sniffles2's task covers [0, len - 1).
            genotype_svs: genotype_input
                .as_ref()
                .and_then(|g| g.svs_by_contig.get(name))
                .map(|v| v.iter().filter(|sv| sv.pos >= 0 && sv.pos < *len - 1).cloned().collect())
                .unwrap_or_default(),
        })
        // Contigs without input SVs have no output when genotyping.
        .filter(|t| genotype_input.is_none() || !t.genotype_svs.is_empty())
        .collect();
    let contigs: Vec<(String, i64)> = contigs.into_iter().map(|(_, name, len)| (name, len)).collect();
    if contigs.is_empty() && !cfg.args.contig.is_empty() {
        bail!("--contig {} names no contig of {input} (chr20 or 20?)", cfg.args.contig.join(","));
    }
    if contigs.is_empty() && !cfg.regions_by_contig.is_empty() {
        bail!("no --regions contig is a contig of {input} (chr20 or 20?)");
    }
    if let Some(fasta) = &reference {
        let missing = fasta.check_contigs(&contigs, &format!("the header of {input}"))?;
        if let Some(warning) = crate::fasta::missing_contigs_warning(&missing) {
            eprintln!("Warning: {warning}");
        }
    }
    let action = if cfg.calls_svs() { "calling SVs on" } else { "reading, without SV calls," };
    eprintln!("lrpipe-sv {}: {action} {} contigs of {}", env!("CARGO_PKG_VERSION"), tasks.len(), input);
    let mut methylation = match (cfg.args.bedmethyl.as_deref(), cfg.args.reference.as_deref()) {
        (Some(path), Some(fasta)) => {
            // Only --contig and --bedmethyl-regions restrict it, not --regions, as for --depth.
            let regions = cfg.args.bedmethyl_regions.as_deref().map(|p| Regions::read(Path::new(p))).transpose()?;
            let wanted: Vec<(usize, String, Vec<String>)> = bam
                .contigs
                .iter()
                .enumerate()
                .filter(|(_, (name, _))| cfg.args.contig.is_empty() || cfg.args.contig.contains(name))
                .filter_map(|(t, (name, _))| match &regions {
                    None => Some((t, name.clone(), vec![name.clone()])),
                    // As lrpipe-methyl's regions: a 0-based start, the end excluded.
                    Some(r) if r.has(name) => {
                        let spans = r.contig(name).iter().map(|(s, e)| format!("{name}:{s}-{e}")).collect();
                        Some((t, name.clone(), spans))
                    }
                    Some(_) => None,
                })
                .collect();
            // In regions, the thresholds are sampled from their reads, and the pass counts none
            // of them: they are read on their own.
            let (sampled, counted): (Vec<String>, Vec<usize>) = match &regions {
                None => {
                    (cfg.args.contig.clone(), tasks.iter().filter(|t| t.regions.is_empty()).map(|t| t.tid).collect())
                }
                Some(_) => (wanted.iter().flat_map(|(.., spans)| spans.iter().cloned()).collect(), Vec::new()),
            };
            Some(Methylation::start(input, fasta, path, cfg.args.threads, &sampled, &wanted, &counted)?)
        }
        _ => None,
    };
    let (methyl_pass, mut methyl_writer) = methylation.as_mut().map(Methylation::parts).unzip();
    let methyl_pass = methyl_pass.filter(|_| cfg.args.bedmethyl_regions.is_none());

    // Two-stage pipeline: contig i is called while contig i + 1 is read, so at most two
    // contigs' leads are in memory. Each stage is parallel within its contig.
    let mut results: Vec<ContigOutput> = Vec::with_capacity(tasks.len());
    let mut depth_reports: HashMap<usize, lrpipe_depth::ContigReport> = HashMap::new();
    let mut repeat_texts: HashMap<String, String> = HashMap::new();
    let mut pending: Option<(&Task, ContigLeads)> = None;
    let mut read_stats = cfg.args.read_stats.as_ref().map(|_| ReadStats::default());
    // How far the pass has got, for `lrpipe run`'s display: the bases of the contigs read, of
    // all of them, which is most of its time.
    let progress = Progress::from_env();
    let bases: i64 = tasks.iter().map(|t| t.len).sum::<i64>().max(1);
    let mut read_bases = 0;
    progress.mark(0.0, &format!("reading 0 of {} contigs", tasks.len()));
    for (contigs_read, next) in tasks.iter().map(Some).chain(std::iter::once(None)).enumerate() {
        let prev = pending.take();
        let (called, read) = rayon::join(
            || {
                prev.map(|(t, mut leads)| {
                    // The contig's depth, repeats and methylation are finished alongside its SV
                    // calls.
                    let depth = leads.depth.take();
                    let repeats = leads.repeats.take();
                    let methylation = leads.methylation.take();
                    let writer = methyl_writer.as_deref_mut();
                    let tr = tandem_repeats.as_ref().map(|m| m.get(&t.name).map_or(&[][..], Vec::as_slice));
                    let ((called, written), (report, repeat_text)) = rayon::join(
                        || {
                            rayon::join(
                                || {
                                    if cfg.calls_svs() {
                                        call_contig(cfg, reference.as_ref(), tr, t, leads)
                                    } else {
                                        Ok(ContigOutput {
                                            text: String::new(),
                                            n_records: 0,
                                            snf: None,
                                            evidence: None,
                                        })
                                    }
                                },
                                || {
                                    writer
                                        .zip(reference.as_ref())
                                        .map(|(w, fasta)| w.write(t.tid, &t.name, methylation, fasta))
                                },
                            )
                        },
                        || {
                            rayon::join(
                                || depth.zip(cfg.depth.as_ref()).map(|(chunks, d)| depth_report(&d.opts, t, chunks)),
                                || {
                                    repeats
                                        .zip(cfg.repeats.as_ref())
                                        .map(|((reference, obs), r)| repeat_records(r, t, &reference, &obs))
                                },
                            )
                        },
                    );
                    (t, called, written, report, repeat_text)
                })
            },
            || next.map(|t| read_contig(cfg, &bam, reference.as_ref(), methyl_pass, t)),
        );
        if let Some((t, c, written, report, repeat_text)) = called {
            results.push(c?);
            written.transpose()?;
            if let Some(report) = report {
                depth_reports.insert(t.tid, report);
            }
            if let Some(text) = repeat_text {
                repeat_texts.insert(t.name.clone(), text);
            }
        }
        if let (Some(t), Some(leads)) = (next, read) {
            let mut leads = leads?;
            if let (Some(total), Some(stats)) = (read_stats.as_mut(), leads.read_stats.take()) {
                total.merge(stats);
            }
            read_bases += t.len;
            progress.report(
                READ_SHARE * read_bases as f64 / bases as f64,
                &format!("read {} ({} of {} contigs)", t.name, contigs_read + 1, tasks.len()),
            );
            pending = Some((t, leads));
        }
    }
    progress.mark(READ_SHARE, "writing the calls");
    if let (config::Haplotypes::List(list), Some(path)) = (&cfg.haplotypes, cfg.args.haplotags.as_deref()) {
        if let Some(report) = list.matches().report(Path::new(path)) {
            eprintln!("lrpipe-sv: {report}");
        }
        if list.ambiguous() > 0 {
            eprintln!(
                "lrpipe-sv: warning: {} reads listed more than once on a contig with different haplotypes (a list \
                 written with --tag-supplementary?) are left untagged",
                list.ambiguous()
            );
        }
    }
    if let Some(m) = methylation {
        m.finish()?;
        eprintln!(
            "Wrote methylation to {} ({:.2}s)",
            cfg.args.bedmethyl.as_deref().unwrap_or_default(),
            start_time.elapsed().as_secs_f64()
        );
    }

    let mut total = 0;
    let mut records = Vec::with_capacity(results.len());
    let mut snf_parts = Vec::new();
    let mut evidence = Vec::new();
    for (task, out) in tasks.iter().zip(results) {
        total += out.n_records;
        records.push(out.text);
        if let Some(blocks) = out.snf {
            snf_parts.push((task.name.clone(), blocks));
        }
        evidence.extend(out.evidence);
    }
    if let Some(path) = cfg.args.phasing_evidence.as_deref() {
        let mut file =
            std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {path}"))?);
        file.write_all(&phasing::header(cfg)?)?;
        for contig in &evidence {
            file.write_all(contig)?;
        }
        file.flush().with_context(|| format!("writing {path}"))?;
        eprintln!("Wrote the SVs' phasing evidence to {path} (for lrpipe sv phase)");
    }
    if let Some(path) = cfg.args.snf.as_deref() {
        snf::write(path, cfg, &contigs, &snf_parts)?;
        let n: usize = snf_parts.iter().map(|(_, p)| p.candidate_count).sum();
        eprintln!("Wrote {n} SV candidates to {path} (for multi-sample calling)");
    }
    if let Some(vcf_path) = cfg.args.vcf.as_deref() {
        let header = match &genotype_input {
            Some(input) => genotype::header(input, cfg, &now_string()),
            None => vcf::header(cfg, &contigs, std::slice::from_ref(&cfg.sample_id), &now_string()),
        };
        write_vcf(vcf_path, &header, &records)?;
        eprintln!("Wrote {total} called SVs to {vcf_path} ({:.2}s)", start_time.elapsed().as_secs_f64());
    }
    // Depth, repeats and the reads' statistics each read the contigs the pass left before
    // writing, side by side.
    let (stats, (depth, repeats)) = rayon::join(
        || -> Result<()> {
            let (Some(path), Some(mut stats)) = (cfg.args.read_stats.as_deref(), read_stats) else { return Ok(()) };
            // The contigs --contig and --regions leave that the pass did not read (the short ones).
            let read: std::collections::HashSet<usize> = tasks.iter().map(|t| t.tid).collect();
            let rest: Vec<usize> = (0..bam.contigs.len())
                .filter(|t| !read.contains(t) && cfg.regions_by_contig.is_empty())
                .filter(|&t| cfg.args.contig.is_empty() || cfg.args.contig.contains(&bam.contigs[t].0))
                .collect();
            bam.count_contigs(&rest, &mut stats)?;
            stats.restricted = !cfg.args.contig.is_empty() || !cfg.regions_by_contig.is_empty();
            stats.write(Path::new(path))?;
            eprintln!("Wrote the reads' statistics to {path} ({:.2}s)", start_time.elapsed().as_secs_f64());
            Ok(())
        },
        || {
            rayon::join(
                || -> Result<()> {
                    if let Some(depth) = &cfg.depth {
                        write_depth(cfg, depth, input, &bam.contigs, depth_reports)?;
                        eprintln!(
                            "Wrote read depth to {}.* ({:.2}s)",
                            depth.prefix,
                            start_time.elapsed().as_secs_f64()
                        );
                    }
                    if let Some(c) = &cfg.cnv {
                        call_cnvs(cfg, c)?;
                        eprintln!("Wrote copy number calls to {} ({:.2}s)", c.dir, start_time.elapsed().as_secs_f64());
                    }
                    Ok(())
                },
                || -> Result<()> {
                    if let Some(r) = &cfg.repeats {
                        let (n, called) = write_repeats(cfg, r, input, repeat_texts)?;
                        eprintln!(
                            "Wrote {called} of {n} tandem repeat genotypes to {} ({:.2}s)",
                            r.output,
                            start_time.elapsed().as_secs_f64()
                        );
                    }
                    Ok(())
                },
            )
        },
    );
    stats?;
    depth?;
    repeats?;
    progress.report(1.0, "done");
    Ok(())
}

/// `--cnv`: copy number (and with --cnv-snv, LoH) calls by lrpipe-cnv from the `--depth`
/// windows just written. In the improved mode the calls take base-pair breakpoints from this
/// run's SV calls.
fn call_cnvs(cfg: &Config, c: &config::CnvConfig) -> Result<()> {
    let depth = cfg.depth.as_ref().context("--cnv needs --depth (checked with the options)")?;
    let reference = cfg.args.reference.as_deref().context("--cnv needs --reference (checked with the options)")?;
    let regions = format!("{}.regions.bed.gz", depth.prefix);
    let mut o = lrpipe_cnv::pipeline::Options::new(
        regions.into(),
        cfg.sample_id.clone(),
        c.dir.clone().into(),
        reference.into(),
    );
    o.metadata = cfg.args.cnv_metadata.as_deref().map(Into::into);
    o.blacklist = cfg.args.cnv_blacklist.as_deref().map(Into::into);
    o.snv = cfg.args.cnv_snv.as_deref().map(Into::into);
    // With --contig the depth windows cover only those contigs.
    o.only_chr.clone_from(&c.only_chr);
    if c.mode == config::CnvMode::Improved {
        o.extensions = lrpipe_cnv::extensions::Extensions::all();
        o.extensions.sv_breakpoints = cfg.args.vcf.as_deref().map(Into::into);
        o.normalize.outlier_removal = lrpipe_cnv::normalize::OutlierRemoval::IsolatedWindows;
        o.vcf_style = lrpipe_cnv::output::vcf::Style::Standard;
    }
    std::fs::create_dir_all(&c.dir).with_context(|| format!("creating {}", c.dir))?;
    let run = lrpipe_cnv::pipeline::call_cnvs(&o).context("calling copy number variants")?;
    lrpipe_cnv::pipeline::write_outputs(&run, &o).context("writing copy number calls")?;
    Ok(())
}

/// Writes `--repeats` output: the pass's records for the contigs it read whole, and for the
/// catalog's other contigs records from lrpipe-tr's own reading, all in reference order; as
/// lrpipe-tr writes the same loci from the same reads. Returns the loci written and how many
/// were called.
fn write_repeats(
    cfg: &Config,
    r: &config::RepeatConfig,
    input: &str,
    mut counted: HashMap<String, String>,
) -> Result<(usize, usize)> {
    let fai_path = r.genome.with_added_extension("fai");
    let fai = pasta::fasta::fai::fs::read(&fai_path).with_context(|| format!("reading {}", fai_path.display()))?;
    let contigs: Vec<(String, u64)> = fai.as_ref().iter().map(|c| (c.name().to_string(), c.length())).collect();
    // Only --contig restricts the output, as for --depth.
    let wanted = |c: &str| cfg.args.contig.is_empty() || cfg.args.contig.iter().any(|x| x == c);
    let rest: Vec<lrpipe_tr::catalog::Locus> = r
        .loci
        .iter()
        .filter(|(c, _)| wanted(c) && !counted.contains_key(c))
        .flat_map(|(_, v)| v.iter().cloned())
        .collect();
    if !rest.is_empty() {
        let (text, _) = lrpipe_tr::run::records(Path::new(input), &r.genome, &rest, &r.opts)?;
        for line in text.lines() {
            let contig = line.split('\t').next().unwrap_or_default();
            let t = counted.entry(contig.to_owned()).or_default();
            t.push_str(line);
            t.push('\n');
        }
    }
    let file = std::fs::File::create(&r.output).with_context(|| format!("creating {}", r.output))?;
    let mut w: Box<dyn std::io::Write> = if r.output.ends_with(".gz") {
        Box::new(pasta::bgzf::io::Writer::new(file))
    } else {
        Box::new(std::io::BufWriter::new(file))
    };
    let header = lrpipe_tr::vcf::header(&contigs, &cfg.sample_id);
    w.write_all(header.as_bytes())?;
    let (mut n, mut called) = (0, 0);
    for (contig, _) in r.loci.iter().filter(|(c, _)| wanted(c)) {
        let text = counted.remove(contig).unwrap_or_default();
        for line in text.lines() {
            n += 1;
            // A no-call's FORMAT values are all missing.
            called += usize::from(!line.ends_with(":.:.:.:."));
        }
        w.write_all(text.as_bytes())?;
    }
    w.flush()?;
    // Dropping the writer ends the BGZF stream, which the index must cover.
    drop(w);
    if r.output.ends_with(".gz") {
        // A fresh index: one left by an earlier run (lrpipe-tr's own stage) would miss records.
        index_bgzipped_vcf(&r.output, &header)?;
    }
    Ok((n, called))
}

/// Writes `--depth` output for all contigs (or the `--contig` ones), as lrpipe-depth would for the
/// whole file. `counted` has the reports of the contigs the SV pass read whole; the others
/// (short contigs, `--regions`, contigs without SVs to genotype) are counted here.
fn write_depth(
    cfg: &Config,
    depth: &DepthConfig,
    input: &str,
    contigs: &[(String, i64)],
    mut counted: HashMap<usize, lrpipe_depth::ContigReport>,
) -> Result<()> {
    let opts = &depth.opts;
    let wanted: Vec<usize> = (0..contigs.len())
        .filter(|&t| cfg.args.contig.is_empty() || cfg.args.contig.contains(&contigs[t].0))
        .filter(|&t| !opts.report.skips(&contigs[t].0))
        .collect();
    let rest: Vec<usize> = wanted.iter().copied().filter(|t| !counted.contains_key(t)).collect();
    let mut rest_reports = HashMap::new();
    if !rest.is_empty() {
        let aln =
            lrpipe_depth::reader::Alignments::open(Path::new(input), cfg.args.reference.as_deref().map(Path::new))?;
        lrpipe_depth::run::count_contigs(&aln, &rest, opts, |tid, report| {
            rest_reports.insert(tid, report);
            Ok(())
        })?;
    }
    let max_len = contigs.iter().map(|c| c.1).max().unwrap_or(0) as u64;
    let mut writer = lrpipe_depth::Writer::new(&depth.prefix, &opts.report, max_len)?;
    for t in wanted {
        let report = counted.remove(&t).or_else(|| rest_reports.remove(&t));
        writer.write_contig(&report.with_context(|| format!("no depth for {}", contigs[t].0))?)?;
    }
    writer.finish()?;
    if let (Some(lrpipe_depth::By::Bed(bed)), true) = (&opts.report.by, cfg.args.contig.is_empty()) {
        let names: Vec<String> = contigs.iter().map(|c| c.0.clone()).collect();
        for (chrom, n) in bed.missing(&names) {
            eprintln!("warning: --depth-by: contig {chrom} from the BED file ({n} regions) is not in the input");
        }
    }
    Ok(())
}

/// Multi-sample calling from .snf files: one task per contig, in parallel.
fn run_combine(cfg: &Config, start_time: Instant) -> Result<()> {
    if cfg.args.snf.is_some() {
        bail!("--snf cannot be used with run mode combine");
    }
    if cfg.depth.is_some() {
        bail!("--depth needs a BAM or CRAM input");
    }
    let Some(vcf_path) = cfg.args.vcf.as_deref() else { bail!("Please specify --vcf for output") };
    if Path::new(vcf_path).exists() && !cfg.args.allow_overwrite {
        bail!("Output file '{vcf_path}' already exists! Use --allow-overwrite to ignore this check and overwrite.");
    }
    let samples = combine::open_samples(&cfg.args.input, cfg.args.combine_close_handles)?;
    let mut contigs = combine::contigs(&samples)?;
    if !cfg.args.contig.is_empty() {
        contigs.retain(|(name, _)| cfg.args.contig.contains(name));
    } else if !cfg.regions_by_contig.is_empty() {
        contigs.retain(|(name, _)| cfg.regions_by_contig.contains_key(name));
    }
    eprintln!(
        "lrpipe-sv {}: calling SVs across {} samples ({} candidates total) on {} contigs",
        env!("CARGO_PKG_VERSION"),
        samples.len(),
        samples.iter().map(|s| s.file.candidate_count).sum::<u64>(),
        contigs.len()
    );
    let reference = cfg.args.reference.as_deref().map(Fasta::open).transpose()?;
    if let Some(fasta) = &reference {
        let missing = fasta.check_contigs(&contigs, "the SNF files")?;
        if let Some(warning) = crate::fasta::missing_contigs_warning(&missing) {
            eprintln!("Warning: {warning}");
        }
    }
    let results: Vec<(String, usize)> = contigs
        .par_iter()
        .enumerate()
        .map(|(task_id, (name, len))| {
            let regions = cfg.regions_by_contig.get(name).map_or(&[][..], Vec::as_slice);
            let calls = combine::combine_contig(&samples, name, *len, regions, task_id, cfg)?;
            let fasta_contig = reference.as_ref().and_then(|fa| fa.contig(name));
            let reference = match (&reference, fasta_contig) {
                (None, _) => vcf::Reference::Absent,
                (Some(_), None) => vcf::Reference::Missing,
                (Some(_), Some(c)) => vcf::Reference::Present(c),
            };
            let mut text = String::new();
            let mut n = 0;
            for c in &calls {
                let before = text.len();
                combine::write_call(&mut text, c, samples.len(), cfg, reference)?;
                n += usize::from(text.len() > before);
            }
            Ok((text, n))
        })
        .collect::<Result<_>>()?;
    let names: Vec<String> = samples.iter().map(|s| s.id.clone()).collect();
    let header = vcf::header(cfg, &contigs, &names, &now_string());
    let total: usize = results.iter().map(|(_, n)| n).sum();
    let records: Vec<String> = results.into_iter().map(|(t, _)| t).collect();
    write_vcf(vcf_path, &header, &records)?;
    eprintln!("Wrote {total} called SVs to {vcf_path} ({:.2}s)", start_time.elapsed().as_secs_f64());
    Ok(())
}

/// Writes the VCF, bgzipped with a tabix index if the name ends in `.gz` or `.bgz`.
/// Writes a VCF: plain, or bgzipped with a tabix index (a CSI one for long contigs) when `path`
/// ends in `.gz` or `.bgz`.
pub(crate) fn write_vcf(path: &str, header: &str, records: &[String]) -> Result<()> {
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    if !config::is_bgzip_name(path) {
        let mut out = std::io::BufWriter::new(file);
        out.write_all(header.as_bytes())?;
        for r in records {
            out.write_all(r.as_bytes())?;
        }
        out.flush()?;
        return Ok(());
    }
    let mut out = pasta::bgzf::io::Writer::new(std::io::BufWriter::new(file));
    out.write_all(header.as_bytes())?;
    for r in records {
        out.write_all(r.as_bytes())?;
    }
    out.finish()?.flush()?;
    index_bgzipped_vcf(path, header)
}

/// Indexes the bgzipped VCF at `path`, whose header is `header`: a tabix index, or a CSI one for
/// long contigs, removing any index of the other kind an earlier run left beside it.
fn index_bgzipped_vcf(path: &str, header: &str) -> Result<()> {
    // The index is the VCF's name plus `.tbi` or `.csi`, where bcftools and htslib look for it.
    // Sniffles2 (pysam.tabix_index) names a `.vcf.bgz`'s index `.vcf.gz.tbi`, which bcftools does
    // not find and which overwrote the index of a `.vcf.gz` of the same stem (lrpipe#52).
    // Contigs too long for tabix (> 2^29 bp) need a CSI index.
    let longest = longest_contig(header).with_context(|| format!("indexing {path}"))?;
    let (min_shift, index_path, other) = if longest > TBI_MAX_LEN {
        (CSI_MIN_SHIFT, format!("{path}.csi"), format!("{path}.tbi"))
    } else {
        (0, format!("{path}.tbi"), format!("{path}.csi"))
    };
    let (format, tabix) = index_vcf(Path::new(path), min_shift).with_context(|| format!("indexing {path}"))?;
    // `save` finishes the BGZF stream and flushes the file, so write errors surface here.
    tabix.index.save(Path::new(&index_path), format).with_context(|| format!("writing {index_path}"))?;
    match std::fs::remove_file(&other) {
        Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
            Err(e).with_context(|| format!("removing the stale index {other}"))
        }
        _ => Ok(()),
    }
}

/// The longest contig a tabix (TBI) index holds: 2^29 bp, from its 14-bit minimum shift and five
/// levels.
const TBI_MAX_LEN: u64 = 1 << 29;

/// The CSI minimum shift of HTSlib's `tabix --csi` and `bcftools index`; HTSlib adds levels to fit
/// the longest `##contig` length.
const CSI_MIN_SHIFT: u32 = 14;

/// The longest `##contig` length in a VCF header, or 0 when none has a length.
fn longest_contig(header: &str) -> Result<u64> {
    let mut longest = 0;
    for line in header.lines() {
        let Some(fields) = line.strip_prefix("##contig=<").and_then(|l| l.strip_suffix('>')) else {
            continue;
        };
        // From the right: an ID may contain commas, and lrpipe-sv writes the length last.
        let Some(len) = fields.rsplit(',').find_map(|f| f.strip_prefix("length=")) else {
            continue;
        };
        match len.parse::<u64>() {
            Ok(len) if len > 0 => longest = longest.max(len),
            _ => bail!("invalid contig length in the VCF header: {line}"),
        }
    }
    Ok(longest)
}

/// Indexes a bgzipped VCF as HTSlib 1.24's `tabix -p vcf` does: TBI when `min_shift` is 0,
/// otherwise CSI with that minimum shift (`tabix --min-shift`).
///
/// # Errors
///
/// If the VCF cannot be read, or a data line cannot be parsed or is out of order (as in `tabix`).
fn index_vcf(path: &Path, min_shift: u32) -> std::io::Result<(IndexFormat, Tabix)> {
    let mut reader = pasta::bgzf::io::Reader::new(std::fs::File::open(path)?);
    let mut builder = TabixBuilder::new(TabixConfig::VCF, min_shift);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        // HTSlib's `bgzf_getline` drops the newline and a carriage return before it.
        if line.last() == Some(&b'\n') {
            line.pop();
            if line.last() == Some(&b'\r') {
                line.pop();
            }
        }
        builder.push_line(&line, u64::from(reader.virtual_position()))?;
    }
    let format = builder.format();
    Ok((format, builder.finish(u64::from(reader.virtual_position()))))
}

/// A contig's leads, coverage spans and mismatch rates, gathered from all its chunks.
#[derive(Debug)]
struct ContigLeads {
    per_type: [Vec<(i64, reads::Lead)>; 5],
    /// Single break leads (`--snf` only).
    single: [Vec<(i64, reads::Lead)>; 2],
    spans: Vec<(i64, i64)>,
    hap_bins: Vec<(i64, i64, u8)>,
    /// `--phasing-evidence`: the read names beside `hap_bins` ([`ChunkLeads::hap_reads`]).
    hap_reads: Vec<Option<std::sync::Arc<[u8]>>>,
    nms: Vec<f64>,
    /// `--depth`, when the whole contig was read: the depth counters of its pieces, in order.
    depth: Option<Vec<lrpipe_depth::ContigCounter>>,
    /// `--repeats`, when the whole contig was read: its reference and the loci's observations.
    repeats: Option<(lrpipe_tr::extract::RefWindow, lrpipe_tr::stream::Observations)>,
    /// `--bedmethyl`, when the whole contig was read and the reference has it: its counts.
    methylation: Option<lrpipe_methyl::stream::Counted>,
    /// `--read-stats`: the records read.
    read_stats: Option<ReadStats>,
}

/// Reads one contig's alignments in parallel chunks and extracts its leads, and with `--depth`,
/// `--repeats` and `--bedmethyl` counts the contig's depth, repeats and methylation from the
/// same reads.
fn read_contig(
    cfg: &Config,
    bam: &Alignments,
    reference: Option<&Fasta>,
    methylation: Option<&lrpipe_methyl::stream::Pass>,
    task: &Task,
) -> Result<ContigLeads> {
    bam.start_contig();
    // Depth needs every read of the contig, but for `--depth-in-regions`: only whole-contig tasks
    // count it here, or then the regions' tasks.
    let depth = cfg
        .depth
        .as_ref()
        .map(|d| &d.opts)
        .filter(|o| (task.regions.is_empty() || o.report.regions_only) && !o.report.skips(&task.name));
    // (region index, region, chunk start, chunk end): each region is read in pieces, in order.
    let regions = task.read_regions();
    let pieces: Vec<(usize, (i64, i64), i64, i64)> = regions
        .iter()
        .enumerate()
        .flat_map(|(ri, &(rs, re))| {
            let n = ((re - rs + CHUNK_LEN - 1) / CHUNK_LEN).max(1);
            (0..n).map(move |k| (ri, (rs, re), rs + k * CHUNK_LEN, (rs + (k + 1) * CHUNK_LEN).min(re)))
        })
        .collect();
    let len = task.len as usize;
    // Repeats too need every read of the contig.
    let repeat_loci =
        cfg.repeats.as_ref().map(|r| r.contig_loci(&task.name)).filter(|l| !l.is_empty() && task.regions.is_empty());
    let repeat_ref = match (repeat_loci, cfg.repeats.as_ref()) {
        (Some(_), Some(r)) => Some(lrpipe_tr::extract::RefWindow::contig(&r.genome, &task.name)?),
        _ => None,
    };
    let repeat_piece = match (repeat_loci, repeat_ref.as_ref(), cfg.repeats.as_ref()) {
        (Some(loci), Some(reference), Some(r)) => Some(reads::RepeatPiece { loci, reference, opts: &r.opts }),
        _ => None,
    };
    // Methylation too needs every read of the contig.
    let methyl_contig = match (methylation, reference) {
        (Some(pass), Some(fasta)) if task.regions.is_empty() => methyl::contig(pass, fasta, &task.name, task.len)?,
        _ => None,
    };
    let chunks: Vec<ChunkLeads> = pieces
        .par_iter()
        .enumerate()
        .map(|(k, &(ri, region, s, e))| {
            // The task region ends a base short of the contig (Sniffles2's [0, len - 1)), so the
            // last piece counts depth and methylation to the contig end.
            let to_end = if k + 1 == pieces.len() { task.len } else { e };
            // A read is counted by the piece it starts in. In regions, the first piece of a region
            // counts those starting before it that reach it but not the region before (those it
            // counted), so that each read reaching the regions is counted once.
            let (depth_start, depth_end) = match task.regions.is_empty() {
                true => (s, to_end),
                false if s == region.0 => (if ri == 0 { 0 } else { regions[ri - 1].1 }, e),
                false => (s, e),
            };
            let depth = depth
                .map(|o| (depth_start, depth_end, lrpipe_depth::ContigCounter::chunk(len, o.mode, o.filter.clone())));
            let mut piece = methyl_contig.as_ref().map(|c| c.piece(s));
            let methylation = piece.as_mut().map(|p| (to_end, p));
            // The reads' statistics of the same reads as the depth: in regions, each read reaching
            // them once.
            let stats = cfg.args.read_stats.as_ref().map(|_| (depth_start, depth_end));
            // Read ids only need to be unique within the contig.
            let id = (k as u64) << 28;
            let leads = bam.chunk_leads(cfg, task.tid, s, e, region, id, depth, repeat_piece, methylation, stats)?;
            if let (Some(c), Some(p)) = (methyl_contig.as_ref(), piece) {
                c.merge(k, p);
            }
            Ok(leads)
        })
        .collect::<Result<_>>()?;
    let mut out = ContigLeads {
        per_type: Default::default(),
        single: Default::default(),
        spans: Vec::new(),
        hap_bins: Vec::new(),
        hap_reads: Vec::new(),
        nms: Vec::new(),
        depth: None,
        repeats: None,
        methylation: methyl_contig.map(lrpipe_methyl::stream::Contig::finish),
        read_stats: None,
    };
    let mut counters = depth.map(|_| Vec::with_capacity(pieces.len()));
    let mut repeat_pieces = Vec::with_capacity(pieces.len());
    let last_region = regions.len() - 1;
    for (c, &(ri, ..)) in chunks.into_iter().zip(&pieces) {
        for (dst, src) in out.per_type.iter_mut().zip(c.leads) {
            dst.extend(src);
        }
        for (dst, src) in out.single.iter_mut().zip(c.single) {
            dst.extend(src);
        }
        out.spans.extend(c.spans);
        out.hap_bins.extend(c.hap_bins);
        out.hap_reads.extend(c.hap_reads);
        if let (Some(counters), Some(piece)) = (counters.as_mut(), c.depth) {
            counters.push(piece);
        }
        repeat_pieces.push(c.repeats);
        if let Some(stats) = c.read_stats {
            out.read_stats.get_or_insert_default().merge(stats);
        }
        // Sniffles2 recomputes the average NM per region, so the last region's wins.
        if ri == last_region {
            out.nms.extend(c.nms);
        }
    }
    out.depth = counters;
    if let (Some(loci), Some(reference), Some(r)) = (repeat_loci, repeat_ref, cfg.repeats.as_ref()) {
        let obs = lrpipe_tr::stream::Observations::merge(repeat_pieces, loci.len(), r.opts.max_reads);
        out.repeats = Some((reference, obs));
    }
    Ok(out)
}

/// Genotypes a contig's repeat loci from the pass's observations, in parallel: its VCF records.
/// Each locus is a job of its own, as a few long loci take seconds each and lie together in
/// the catalog.
fn repeat_records(
    r: &config::RepeatConfig,
    task: &Task,
    reference: &lrpipe_tr::extract::RefWindow,
    obs: &lrpipe_tr::stream::Observations,
) -> String {
    let loci = r.contig_loci(&task.name);
    let records: Vec<String> = (0..loci.len())
        .into_par_iter()
        .with_max_len(1)
        .map(|i| obs.records(loci, i..i + 1, reference, &r.opts).0)
        .collect();
    records.concat()
}

/// Combines a contig's per-piece depth counters, in read order, into its depth report.
fn depth_report(
    opts: &lrpipe_depth::run::RunOptions,
    task: &Task,
    pieces: Vec<lrpipe_depth::ContigCounter>,
) -> lrpipe_depth::ContigReport {
    let mut counter = lrpipe_depth::ContigCounter::chunk(task.len as usize, opts.mode, opts.filter.clone());
    for piece in pieces {
        counter.append(piece);
    }
    let coverage = counter.finish();
    lrpipe_depth::ContigReport::new(
        &task.name,
        task.len as u32,
        &coverage,
        &opts.report,
        &mut lrpipe_depth::Scratch::new(&opts.report),
    )
}

/// One contig's output: VCF records and, with `--snf`, its SNF blocks, and with
/// `--phasing-evidence` its evidence.
#[derive(Debug)]
struct ContigOutput {
    text: String,
    n_records: usize,
    snf: Option<snf::ContigBlocks>,
    /// A gzip member.
    evidence: Option<Vec<u8>>,
}

/// Calls one contig from its leads and returns its VCF records.
fn call_contig(
    cfg: &Config,
    reference: Option<&Fasta>,
    tandem_repeats: Option<&[(i64, i64)]>,
    task: &Task,
    leads: ContigLeads,
) -> Result<ContigOutput> {
    let ContigLeads {
        per_type,
        single,
        spans,
        hap_bins,
        hap_reads,
        nms,
        depth: _,
        repeats: _,
        methylation: _,
        read_stats: _,
    } = leads;
    let mut tabs: [Leadtab; 5] = per_type.map(|v| Leadtab::new(v, cfg.args.consensus_max_reads_bin));
    let hap_ref = HapRef::new(&hap_bins);
    // The phasing evidence names the reads at some bins.
    let hap_bins = if cfg.args.phasing_evidence.is_some() { hap_bins } else { Vec::new() };
    // Single breaks are clustered like inversions (Sniffles2's generic path) and become
    // candidates at the center of their leads' starts; only their positions are kept.
    let single_positions: Vec<i64> = single
        .into_iter()
        .flat_map(|v| {
            let mut tab = Leadtab::new(v, cfg.args.consensus_max_reads_bin);
            let clusters = cluster::resolve(SvType::Inv, &mut tab, &hap_ref, tandem_repeats, cfg);
            clusters
                .into_iter()
                .filter(|c| !c.leads.is_empty())
                .map(|c| {
                    let starts: Vec<i64> = c.leads.iter().map(|&i| tab.leads[i].ref_start).collect();
                    pymath::center(&starts)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    // average_regional_nm: a sequential sum over reads, as in Sniffles2.
    let nm_sum = nms.iter().fold(0.0, |s, &v| s + v);
    let nm_threshold = nm_sum / nms.len().max(1) as f64;

    // SV types are clustered independently; IDs are numbered afterwards in Sniffles2's order.
    let per_type_calls: Vec<Vec<SvCall>> = tabs
        .par_iter_mut()
        .zip(SvType::ALL)
        .map(|(tab, svtype)| {
            let mut clusters = cluster::resolve(svtype, tab, &hap_ref, tandem_repeats, cfg);
            if cfg.args.split_alleles && matches!(svtype, SvType::Ins | SvType::Del) {
                clusters = clusters.into_iter().flat_map(|c| alleles::split(c, tab, cfg)).collect();
            }
            clusters
                .into_iter()
                .filter_map(|mut c| {
                    c.compute_sa_counts(tab);
                    call::call_from(c, tab, cfg)
                })
                .collect()
        })
        .collect();
    let mut candidates: Vec<SvCall> = per_type_calls.into_iter().flatten().collect();
    for (sv_id, c) in candidates.iter_mut().enumerate() {
        c.id = format!("{}.{sv_id:X}S{:X}", c.svtype.as_str(), task.id);
    }

    let fasta_contig = reference.and_then(|fa| fa.contig(&task.name));
    // Only needed for coverage and DEL checks of candidate calls.
    let n_runs = fasta_contig
        .filter(|_| !candidates.is_empty())
        .map(|fa| fa.n_runs())
        .transpose()
        .with_context(|| format!("reading {}", task.name))?;
    // Sniffles2 masks N bases only inside --regions; for a whole contig its numpy masking fails
    // (and is skipped) if the reference and BAM lengths differ.
    let region_mask = n_runs.as_ref().filter(|_| !task.regions.is_empty()).map(|r| r.clipped_to(&task.regions));
    let n_mask = match &region_mask {
        Some(m) => Some(m),
        None => n_runs.as_ref().filter(|_| fasta_contig.is_some_and(|fa| fa.len() == task.len)),
    };
    let coverage = Coverage::new(&spans, task.len, n_mask);
    drop(spans);
    let coverage_total =
        if matches!(cfg.minsupport, MinSupport::Auto) || cfg.args.mosaic { coverage.mean() } else { 0.0 };
    call::annotate_coverage(&mut candidates, &coverage, cfg);

    let ctx = ContigCtx { cfg, tabs: &tabs };
    candidates.par_iter_mut().for_each(|c| call::finalize(c, &ctx, coverage_total, nm_threshold));
    if !task.genotype_svs.is_empty() {
        // Genotyping uses all candidates, without the QC filter, sorting or the N-region check.
        let text = genotype::write_records(&task.genotype_svs, &candidates, &coverage, cfg);
        return Ok(ContigOutput { text, n_records: task.genotype_svs.len(), snf: None, evidence: None });
    }

    if let Some(runs) = &n_runs
        && cfg.args.max_unknown_pct < 1.0
    {
        for c in candidates.iter_mut().filter(|c| c.svtype == SvType::Del) {
            if n_percentage(runs, c.pos, c.end) > cfg.args.max_unknown_pct {
                c.filter = "UNRESOLVED";
                c.qc = false;
            }
        }
    }
    let reference_of = || match (reference, fasta_contig) {
        (None, _) => vcf::Reference::Absent,
        (Some(_), None) => vcf::Reference::Missing,
        (Some(_), Some(c)) => vcf::Reference::Present(c),
    };
    let evidence = cfg
        .args
        .phasing_evidence
        .is_some()
        .then(|| phasing::contig_evidence(&task.name, &candidates, &tabs, &hap_bins, &hap_reads, cfg, reference_of()))
        .transpose()?;
    drop((hap_bins, hap_reads));
    // The SNF keeps all candidates, QC failures included, in calling order.
    let snf = cfg
        .args
        .snf
        .is_some()
        .then(|| snf::encode_contig(&candidates, &single_positions, &task.name, &coverage, cfg))
        .transpose()?;
    if !cfg.args.no_qc {
        candidates.retain(|c| c.qc);
    }
    if cfg.sort {
        candidates.sort_by_key(|c| c.pos);
    }
    let mut text = String::new();
    let mut n = 0;
    for c in &mut candidates {
        let before = text.len();
        vcf::write_call(&mut text, c, &task.name, cfg, reference_of())?;
        n += usize::from(text.len() > before);
    }
    Ok(ContigOutput { text, n_records: n, snf, evidence })
}

/// `UnresolvedReferenceRegions.get_N_percentage`: percent of [pos, end) in runs of 100+ `N`s.
fn n_percentage(runs: &NRuns, pos: i64, end: i64) -> f64 {
    if end <= pos {
        return 0.0;
    }
    runs.long_overlap(pos, end, 100) as f64 / (end - pos) as f64 * 100.0
}

/// `util.load_tandem_repeats`: padded intervals per contig, sorted.
fn load_tandem_repeats(path: &str, pad: i64) -> Result<HashMap<String, Vec<(i64, i64)>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut map: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    let mut unsorted = false;
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let (start, end): (i64, i64) = (parts[1].trim().parse()?, parts[2].trim().parse()?);
        let v = map.entry(parts[0].to_string()).or_default();
        if v.last().is_some_and(|&(last_start, _)| start < last_start) {
            unsorted = true;
        }
        v.push(((start - pad).max(0), end + pad));
    }
    if unsorted {
        for v in map.values_mut() {
            v.sort_unstable();
        }
    }
    Ok(map)
}

/// The current UTC time as `YYYY/MM/DD HH:MM:SS`.
fn now_string() -> String {
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).map_or(0, |d| d.as_secs()) as i64;
    let (days, rem) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    // Civil date from days since 1970-01-01 (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}/{m:02}/{d:02} {:02}:{:02}:{:02}", rem / 3600, rem % 3600 / 60, rem % 60)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A scratch directory for one test, removed on drop.
    struct ScratchDir(std::path::PathBuf);

    impl ScratchDir {
        fn new(test: &str) -> Self {
            let dir = std::env::temp_dir().join(format!("lrpipe-sv-cli-{}-{test}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("creating the scratch directory");
            Self(dir)
        }

        fn path(&self, name: &str) -> String {
            self.0.join(name).to_str().expect("a UTF-8 temporary directory").to_owned()
        }
    }

    impl Drop for ScratchDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn header(contigs: &[(&str, u64)]) -> String {
        let mut h = String::from("##fileformat=VCFv4.2\n");
        for (name, len) in contigs {
            h.push_str(&format!("##contig=<ID={name},length={len}>\n"));
        }
        h.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        h
    }

    fn record(contig: &str, pos: u64) -> String {
        format!("{contig}\t{pos}\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-100;END={}\n", pos + 100)
    }

    #[test]
    fn longest_contig_reads_the_contig_lengths() {
        assert_eq!(longest_contig("##fileformat=VCFv4.2\n").expect("no contigs"), 0);
        assert_eq!(longest_contig(&header(&[("a", 5), ("b", 13), ("c", 8)])).expect("valid lengths"), 13);
        let h = "##contig=<ID=chr1,length=248956422,assembly=GRCh38>\n##contig=<ID=chrUn>\n";
        assert_eq!(longest_contig(h).expect("valid lengths"), 248_956_422);
        assert!(longest_contig("##contig=<ID=a,length=0>\n").is_err());
        assert!(longest_contig("##contig=<ID=a,length=x>\n").is_err());
    }

    #[test]
    fn short_contigs_get_a_tabix_index_named_after_the_vcf() {
        let dir = ScratchDir::new("tbi");
        let vcf = dir.path("calls.vcf.bgz");
        let records = [record("chr1", 100), record("chr1", 5000), record("chr2", 70)];
        write_vcf(&vcf, &header(&[("chr1", TBI_MAX_LEN), ("chr2", 1000)]), &records).expect("writing the VCF");
        let tabix = Tabix::load(Path::new(&dir.path("calls.vcf.bgz.tbi"))).expect("loading the index");
        assert_eq!(tabix.index.format(), IndexFormat::Tbi);
        assert_eq!(tabix.names(), [b"chr1".to_vec(), b"chr2".to_vec()]);
        assert!(!Path::new(&dir.path("calls.vcf.bgz.csi")).exists());
        assert!(!Path::new(&dir.path("calls.vcf.gz.tbi")).exists(), "Sniffles2's name for it");
    }

    #[test]
    fn a_bgz_vcf_leaves_the_index_of_a_gz_vcf_of_the_same_stem() {
        // A `.vcf.bgz` named its index `.vcf.gz.tbi`, replacing that of a `.vcf.gz` beside it,
        // whose indexed lookups then found nothing (lrpipe#52); with CSI indexes alike.
        for (kind, len) in [("tbi", TBI_MAX_LEN), ("csi", 600_000_000)] {
            let dir = ScratchDir::new(&format!("stem-{kind}"));
            let h = header(&[("chr1", len)]);
            write_vcf(&dir.path("calls.vcf.gz"), &h, &[record("chr1", 100), record("chr1", 5000)]).expect("the .gz");
            let gz_index = dir.path(&format!("calls.vcf.gz.{kind}"));
            let before = std::fs::read(&gz_index).expect("the .gz's index");
            write_vcf(&dir.path("calls.vcf.bgz"), &h, &[record("chr1", 9000)]).expect("the .bgz");
            assert_eq!(std::fs::read(&gz_index).expect("the .gz's index"), before, "{kind}");
            assert!(Path::new(&dir.path(&format!("calls.vcf.bgz.{kind}"))).exists(), "{kind}");
        }
    }

    #[test]
    fn contigs_longer_than_tabix_allows_get_a_csi_index() {
        let dir = ScratchDir::new("csi");
        let vcf = dir.path("calls.vcf.gz");
        let records = [record("chr1", 100), record("chr1", 550_000_000)];
        write_vcf(&vcf, &header(&[("chr1", 600_000_000)]), &records).expect("writing the VCF");
        let tabix = Tabix::load(Path::new(&dir.path("calls.vcf.gz.csi"))).expect("loading the index");
        assert_eq!(tabix.index.format(), IndexFormat::Csi);
        // HTSlib's `tabix --csi`: a 14-bit minimum shift, with six levels to cover 600 Mbp.
        assert_eq!((tabix.index.min_shift(), tabix.index.n_lvls()), (14, 6));
        assert!(tabix.index.query(0, 549_999_999, 550_000_000).is_some());
        assert!(!Path::new(&dir.path("calls.vcf.gz.tbi")).exists());
    }
}

//! The run's QC summary: `lrpipe qc`, the qc stage, reads what the other stages wrote and writes
//! `qc/SAMPLE.qc.json` and `qc/SAMPLE.qc.tsv`, and `lrpipe run` prints its main values at the end.
//!
//! Nothing is read again for it but the reads' index. The reads (the primary alignments of
//! mapped reads: count, yield, length N50, mean length and base quality, and the secondary and
//! supplementary alignments) are counted by the read pass as it reads them
//! ([`Layout::read_stats`]); the reads' index gives the unplaced unmapped reads, for the
//! fraction mapped, and the sex, checked against the sex the run called the sample with. The
//! coverage is the depth summary's (MAPQ 20 or more), the SNVs and indels are lrpipe-snv's PASS
//! calls, the SVs the read pass's PASS calls, the phasing WhatsHap's `stats` of the phased SNVs
//! (lrpipe-phase's), and the CpGs the rows of the bedMethyl or CpG tables. What the run did not
//! write is said to be not measured, rather than read from the reads again.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use serde::Serialize;
use serde::ser::SerializeMap;

use crate::config::Site;
use crate::plan::{Layout, Priority, Sample, Sex, Stage, StageName};

/// `lrpipe qc`'s arguments: the run's files to summarise, each optional but the reads.
#[derive(Debug, clap::Args)]
pub struct QcArgs {
    /// Sample name
    #[arg(long)]
    pub sample: String,
    /// The sex checked against the reads' index, or not-applicable for bacterial calling
    #[arg(long, value_parser = ["female", "male", "not-applicable"])]
    pub sex: String,
    /// The reads: a sorted, indexed BAM or CRAM, whose index gives the sex and the unplaced
    /// unmapped reads
    #[arg(long)]
    pub reads: PathBuf,
    /// The files whose indexes give the sex, as `lrpipe run` estimates it as it starts: each of
    /// several aligned --reads, which the run merges into --reads, so that the summary's
    /// estimate is the same [default: --reads]
    #[arg(long)]
    pub sex_reads: Vec<PathBuf>,
    /// The regions of a run in regions (`lrpipe run --regions`, as it merged them), which the
    /// files are of
    #[arg(long)]
    pub regions: Option<PathBuf>,
    /// The read pass's statistics of the reads (`lrpipe sv --read-stats`)
    #[arg(long)]
    pub read_stats: Option<PathBuf>,
    /// The depth summary (mosdepth's, `lrpipe sv --depth`)
    #[arg(long)]
    pub depth_summary: Option<PathBuf>,
    /// SNVs and indels (lrpipe-snv's merge_output.vcf.gz)
    #[arg(long)]
    pub snv_vcf: Option<PathBuf>,
    /// The phased SNVs and indels (lrpipe-snv's phased_merge_output.vcf.gz)
    #[arg(long)]
    pub phased_vcf: Option<PathBuf>,
    /// Structural variants (`lrpipe sv`)
    #[arg(long)]
    pub sv_vcf: Option<PathBuf>,
    /// CpG methylation as bedMethyl (`lrpipe sv --bedmethyl`)
    #[arg(long)]
    pub bedmethyl: Option<PathBuf>,
    /// pb-CpG-tools' table of all reads (`lrpipe methyl cpg`)
    #[arg(long)]
    pub cpg_combined: Option<PathBuf>,
    /// The table of haplotype 1's reads
    #[arg(long)]
    pub cpg_hap1: Option<PathBuf>,
    /// The table of haplotype 2's reads
    #[arg(long)]
    pub cpg_hap2: Option<PathBuf>,
    /// The summary as JSON
    #[arg(long)]
    pub json: PathBuf,
    /// The summary as a table: one value per line, its name and the value (NA when not measured)
    #[arg(long)]
    pub tsv: PathBuf,
}

/// The qc stage: `lrpipe qc` on the files of `stages` it summarises, last, after every one of
/// `stages`. Its inputs are those files and the reads with their index, so it reruns when any
/// changes. Its failure does not fail the run ([`StageName::fails_the_run`]).
#[must_use]
pub fn stage(site: &Site, sample: &Sample, layout: &Layout, stages: &[Stage]) -> Stage {
    let writer = |path: &PathBuf| stages.iter().find(|s| s.outputs.contains(path)).map(|s| s.name);
    let reads = &sample.reads;
    // The reads' index: the one the stage writing the reads writes, else the one beside them.
    let index = ["bai", "csi", "crai"]
        .into_iter()
        .map(|e| Layout::index(reads, e))
        .find(|i| writer(i).is_some())
        .or_else(|| crate::plan::index_of(reads));
    let files = [
        ("--read-stats", layout.read_stats()),
        ("--depth-summary", layout.depth_summary()),
        ("--snv-vcf", layout.snv_vcf()),
        ("--phased-vcf", layout.phased_snv_vcf()),
        ("--sv-vcf", layout.sv_vcf()),
        ("--bedmethyl", layout.bedmethyl()),
        ("--cpg-combined", layout.cpg_table("combined")),
        ("--cpg-hap1", layout.cpg_table("hap1")),
        ("--cpg-hap2", layout.cpg_table("hap2")),
    ];
    let mut args: Vec<std::ffi::OsString> = vec![
        "qc".into(),
        "--sample".into(),
        sample.name.clone().into(),
        "--sex".into(),
        sample.sex.as_str().into(),
        "--reads".into(),
        reads.into(),
    ];
    let mut inputs: Vec<PathBuf> = std::iter::once(reads.clone()).chain(index).collect();
    // A run in regions: what is summarised is of them, and other regions rerun it.
    if sample.regions.is_some() {
        args.extend(["--regions".into(), layout.regions_bed().into()]);
        inputs.push(layout.regions_bed());
    }
    // The sex as the run estimated it: over the several files merged into the reads.
    if sample.sex != Sex::NotApplicable {
        for input in sample.merged.iter().flat_map(|m| &m.inputs) {
            args.extend(["--sex-reads".into(), input.into()]);
            inputs.push(input.clone());
            inputs.extend(crate::plan::index_of(input));
        }
    }
    for (flag, path) in files {
        if writer(&path).is_some() {
            args.extend([flag.into(), path.clone().into()]);
            inputs.push(path);
        }
    }
    // Last: after every other stage, so that it never runs beside one.
    let after: Vec<StageName> = stages.iter().map(|s| s.name).collect();
    args.extend(["--json".into(), layout.qc_json().into(), "--tsv".into(), layout.qc_tsv().into()]);
    Stage {
        name: StageName::Qc,
        program: site.tools.lrpipe.clone().unwrap_or_else(|| PathBuf::from("lrpipe")),
        args,
        inputs,
        outputs: vec![layout.qc_json(), layout.qc_tsv()],
        after,
        memory_gib: crate::plan::expected_memory(StageName::Qc, sample.platform),
        priority: Priority::Critical,
        scratch: None,
    }
}

/// The summary.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Summary {
    /// The sample.
    pub sample: String,
    /// The regions of a run in regions (`lrpipe run --regions`), which every value is of but the
    /// sex; none for a whole genome.
    pub regions: Option<Within>,
    /// The reads, from the read pass; none without one.
    pub reads: Option<Reads>,
    /// Mean depth, from the depth summary; none without one.
    pub coverage: Option<Coverage>,
    /// The sex the run called the sample with, against the index's.
    pub sex: SexCheck,
    /// PASS SNVs; none without SNV calls.
    pub snvs: Option<Snvs>,
    /// PASS indels; none without SNV calls.
    pub indels: Option<Indels>,
    /// PASS SVs; none without SV calls.
    pub svs: Option<Svs>,
    /// The SNVs' phasing; none without phasing.
    pub phasing: Option<Phasing>,
    /// CpGs reported; none without methylation.
    pub methylation: Option<Methylation>,
    /// What was not measured, and why.
    pub not_measured: Vec<String>,
    /// What looks wrong.
    pub warnings: Vec<String>,
}

/// The regions a summary is of.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Within {
    /// How many there are (merged).
    pub regions: usize,
    /// Their bases.
    pub bases: u64,
}

/// The reads: the primary alignments of mapped reads (duplicates and QC failures included, as
/// they are reads sequenced: long-read pipelines rarely mark them), and the other records.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Reads {
    /// Mapped reads (primary alignments).
    pub mapped_reads: u64,
    /// Their bases: the mapped yield (the unmapped reads' lengths are not read).
    pub mapped_bases: u64,
    /// The mapped reads' length N50.
    pub n50: u64,
    /// Their mean length.
    pub mean_length: Option<f64>,
    /// The arithmetic mean of their bases' Phred qualities, over the reads that have qualities:
    /// samtools stats' "average quality". NanoPlot and cramino average the bases' error
    /// probabilities instead, and report lower values.
    pub mean_quality: Option<f64>,
    /// Unmapped reads: those placed at a position, which the read pass reads, and the unplaced
    /// ones by the index; none when the index does not count them.
    pub unmapped: Option<u64>,
    /// Mapped reads, as a percentage of all reads.
    pub mapped_percent: Option<f64>,
    /// Secondary alignments.
    pub secondary: u64,
    /// Secondary alignments, as a percentage of all alignments (primary, secondary and
    /// supplementary).
    pub secondary_percent: Option<f64>,
    /// Supplementary alignments.
    pub supplementary: u64,
    /// Supplementary alignments, as a percentage of all alignments.
    pub supplementary_percent: Option<f64>,
}

/// Mean depth (MAPQ 20 or more, each read's span).
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Coverage {
    /// Over all the autosomes (1 to 22) there are.
    pub autosomes: Option<f64>,
    /// chrX's.
    pub chrx: Option<f64>,
    /// chrY's.
    pub chry: Option<f64>,
    /// Each autosome's, in the reference's order.
    pub per_autosome: Means,
}

/// Contigs' mean depths, in order; a JSON object.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Means(pub Vec<(String, f64)>);

impl Serialize for Means {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut map = serializer.serialize_map(Some(self.0.len()))?;
        for (contig, mean) in &self.0 {
            map.serialize_entry(contig, mean)?;
        }
        map.end()
    }
}

/// The sex the run called the sample with, against the one its reads' index gives
/// ([`crate::sex`]).
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SexCheck {
    /// The run's inferred or specified sex, or `not-applicable` for a bacterial run.
    pub called: String,
    /// The index's: female, male or unclear (fitting neither); none when it cannot tell.
    pub index: Option<String>,
    /// chrX's reads per base against the autosomes'.
    pub chrx_ratio: Option<f64>,
    /// chrY's.
    pub chry_ratio: Option<f64>,
    /// Whether the two agree; none when the index cannot tell.
    pub agrees: Option<bool>,
    /// Why the index cannot tell.
    pub note: Option<String>,
}

/// PASS SNVs.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Snvs {
    /// Records.
    pub pass: u64,
    /// Heterozygous (0/1, 1/2) calls.
    pub het: u64,
    /// Homozygous ALT calls.
    pub hom_alt: u64,
    /// Haploid ALT calls.
    pub hemizygous: u64,
    /// `het` over `hom_alt`.
    pub het_hom_ratio: Option<f64>,
    /// ALT alleles that are transitions (A-G, C-T).
    pub transitions: u64,
    /// ALT alleles that are transversions.
    pub transversions: u64,
    /// `transitions` over `transversions`.
    pub ti_tv: Option<f64>,
}

/// PASS indels: the records that are not SNVs, MNPs (REF and ALT of more than one base, the same
/// length) included, which bcftools counts apart (Clair3 calls none).
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Indels {
    /// Records.
    pub pass: u64,
    /// Heterozygous calls.
    pub het: u64,
    /// Homozygous ALT calls.
    pub hom_alt: u64,
    /// Haploid ALT calls.
    pub hemizygous: u64,
    /// `het` over `hom_alt`.
    pub het_hom_ratio: Option<f64>,
}

/// PASS SVs.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Svs {
    /// Records.
    pub pass: u64,
    /// Records of each SVTYPE, by name.
    pub by_type: std::collections::BTreeMap<String, u64>,
}

/// The SNVs' phasing, as WhatsHap's `stats` gives it over all chromosomes: of every record, not
/// only the PASS ones, and multi-allelic sites left out, so its heterozygous SNVs are not
/// [`Snvs::het`].
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Phasing {
    /// Heterozygous SNVs, of any FILTER, multi-allelic sites left out.
    pub het_snvs_any_filter: i64,
    /// Those in a phase block of two or more variants.
    pub phased_snvs: i64,
    /// `phased_snvs` as a percentage of `het_snvs_any_filter`.
    pub phased_percent: Option<f64>,
    /// Phase blocks of two or more variants.
    pub blocks: i64,
    /// The blocks' NG50 (`stats`' `block_n50`): the span of the shortest of the longest blocks,
    /// split where blocks interleave, that together span at least half the chromosomes with
    /// blocks.
    pub block_n50: Option<i64>,
}

/// CpGs reported.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Methylation {
    /// Of all the reads (bedMethyl, or the combined CpG table).
    pub combined: Option<u64>,
    /// Of haplotype 1's reads (the CpG table).
    pub hap1: Option<u64>,
    /// Of haplotype 2's.
    pub hap2: Option<u64>,
}

/// Runs `lrpipe qc`: summarises the files of `args` and writes the summary, saying its warnings.
///
/// # Errors
///
/// If a file given cannot be read or is not what it should be, or the summary cannot be written.
pub fn main(args: &QcArgs) -> Result<()> {
    check_files(args)?;
    let summary = summarise(args)?;
    for path in [&args.json, &args.tsv] {
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
        }
    }
    let json = serde_json::to_string_pretty(&summary).expect("a summary serialises") + "\n";
    std::fs::write(&args.json, json).with_context(|| format!("writing {}", args.json.display()))?;
    std::fs::write(&args.tsv, table(&summary)).with_context(|| format!("writing {}", args.tsv.display()))?;
    for warning in &summary.warnings {
        eprintln!("lrpipe qc: warning: {warning}");
    }
    Ok(())
}

/// Refuses an output that is one of the inputs or the other output, before anything is read or
/// written, and VCFs of another sample than `--sample`: their records would be summarised as
/// its.
fn check_files(args: &QcArgs) -> Result<()> {
    // Following links: an output that is a link to an input would be written through it.
    let spelled = |path: &Path| {
        path.canonicalize().unwrap_or_else(|_| {
            crate::provenance::canonical(&std::path::absolute(path).unwrap_or_else(|_| path.to_owned()))
        })
    };
    let inputs: Vec<(&str, &PathBuf)> = [
        ("--reads", Some(&args.reads)),
        ("--regions", args.regions.as_ref()),
        ("--read-stats", args.read_stats.as_ref()),
        ("--depth-summary", args.depth_summary.as_ref()),
        ("--snv-vcf", args.snv_vcf.as_ref()),
        ("--phased-vcf", args.phased_vcf.as_ref()),
        ("--sv-vcf", args.sv_vcf.as_ref()),
        ("--bedmethyl", args.bedmethyl.as_ref()),
        ("--cpg-combined", args.cpg_combined.as_ref()),
        ("--cpg-hap1", args.cpg_hap1.as_ref()),
        ("--cpg-hap2", args.cpg_hap2.as_ref()),
    ]
    .into_iter()
    .filter_map(|(option, path)| path.map(|p| (option, p)))
    .chain(args.sex_reads.iter().map(|p| ("--sex-reads", p)))
    .collect();
    if spelled(&args.json) == spelled(&args.tsv) {
        bail!("--json and --tsv are the same file, {}: give each its own", args.json.display());
    }
    for (output, path) in [("--json", &args.json), ("--tsv", &args.tsv)] {
        if let Some((input, _)) = inputs.iter().find(|(_, input)| spelled(input) == spelled(path)) {
            bail!("{output} {} would overwrite the {input} input: write the summary to another file", path.display());
        }
    }
    for (option, vcf) in [("--snv-vcf", &args.snv_vcf), ("--phased-vcf", &args.phased_vcf), ("--sv-vcf", &args.sv_vcf)]
    {
        let Some(vcf) = vcf else { continue };
        let mut samples = None;
        for line in lines(vcf)?.lines() {
            let line = line.with_context(|| format!("reading {}", vcf.display()))?;
            if let Some(columns) = line.strip_prefix("#CHROM") {
                samples = Some(columns.split('\t').skip(9).map(str::to_owned).collect::<Vec<_>>());
                break;
            }
            if !line.starts_with('#') {
                break;
            }
        }
        match samples.as_deref() {
            Some([one]) if one != &args.sample => {
                bail!("{option} {} is sample {one}'s, not --sample {}", vcf.display(), args.sample)
            }
            Some([]) | Some([_]) | None => {}
            Some(many) => {
                bail!("{option} {} has several samples ({}): give one sample's VCF", vcf.display(), many.join(", "))
            }
        }
    }
    Ok(())
}

/// The summary of `args`' files.
///
/// # Errors
///
/// If a file given cannot be read or is not what it should be.
pub fn summarise(args: &QcArgs) -> Result<Summary> {
    let mut not_measured = Vec::new();
    let mut warnings = Vec::new();
    let regions = match &args.regions {
        Some(path) => {
            let r =
                lrpipe_reads::regions::Regions::read(path).with_context(|| format!("reading {}", path.display()))?;
            Some(Within { regions: r.iter().map(|(_, spans)| spans.len()).sum(), bases: r.bases() })
        }
        None => None,
    };
    let reads = match &args.read_stats {
        Some(path) => {
            let (r, error, restricted) = reads(path, &args.reads)?;
            if r.mean_quality.is_none() {
                not_measured.push(if crate::plan::is_cram(&args.reads) {
                    "reads.mean_quality: the read pass does not decode a CRAM's qualities (it would take two thirds \
                     longer)"
                        .into()
                } else {
                    "reads.mean_quality: the reads have no base qualities".into()
                });
            }
            match error {
                _ if restricted => not_measured.push(
                    "reads.unmapped: the reads are counted within the regions, the unplaced unmapped reads only in \
                     the whole file"
                        .into(),
                ),
                Some(error) => {
                    not_measured.push(format!("reads.unmapped: {error}"));
                    warnings.push(error);
                }
                None if r.unmapped.is_none() => not_measured
                    .push("reads.unmapped: the reads' index does not count the unplaced unmapped reads".into()),
                None => {}
            }
            Some(r)
        }
        None => {
            not_measured.push(
                "reads: the read pass, which counts them, did not run (it runs for sv, cnv, depth and tr, and for \
                 methyl on ONT)"
                    .into(),
            );
            None
        }
    };
    let coverage = match &args.depth_summary {
        Some(path) => Some(coverage(path)?),
        None => {
            not_measured.push("coverage: depth was not measured (depth or cnv among --analyses)".into());
            None
        }
    };
    let sex_reads = if args.sex_reads.is_empty() { std::slice::from_ref(&args.reads) } else { &args.sex_reads[..] };
    let (sex, estimate) = sex_check(sex_reads, &args.sex);
    warnings.extend(estimate.and_then(|e| {
        let called = if args.sex == "male" { Sex::Male } else { Sex::Female };
        sex_warning(&e, called)
    }));
    // The variants and the CpGs, each file on a thread of its own: a whole genome's bedMethyl
    // alone takes half a minute.
    let (small, svs, phasing, [bedmethyl, combined, hap1, hap2]) = std::thread::scope(|scope| {
        let small = args.snv_vcf.as_deref().map(|p| scope.spawn(move || small_variants(p)));
        let svs = args.sv_vcf.as_deref().map(|p| scope.spawn(move || svs(p)));
        let phasing = args.phased_vcf.as_deref().map(|p| scope.spawn(move || phasing(p)));
        let tables = [&args.bedmethyl, &args.cpg_combined, &args.cpg_hap1, &args.cpg_hap2]
            .map(|p| p.as_deref().map(|p| scope.spawn(move || cpgs(p))));
        Ok::<_, anyhow::Error>((joined(small)?, joined(svs)?, joined(phasing)?, {
            let [a, b, c, d] = tables;
            [joined(a)?, joined(b)?, joined(c)?, joined(d)?]
        }))
    })?;
    let (snvs, indels) = small.unzip();
    if snvs.is_none() {
        not_measured.push("snvs and indels: no SNVs were called (snv, phase or cnv among --analyses)".into());
    }
    if svs.is_none() {
        not_measured.push("svs: no SVs were called (sv or cnv among --analyses)".into());
    }
    if phasing.is_none() {
        not_measured.push("phasing: the SNVs were not phased (phase among --analyses)".into());
    }
    let methylation = Methylation { combined: bedmethyl.or(combined), hap1, hap2 };
    let methylation = if methylation == Methylation::default() {
        not_measured.push("methylation: none was reported (methyl among --analyses)".into());
        None
    } else {
        Some(methylation)
    };
    Ok(Summary {
        sample: args.sample.clone(),
        regions,
        reads,
        coverage,
        sex,
        snvs,
        indels,
        svs,
        phasing,
        methylation,
        not_measured,
        warnings,
    })
}

/// What a thread of [`summarise`] gave, if it ran.
fn joined<T>(thread: Option<std::thread::ScopedJoinHandle<'_, Result<T>>>) -> Result<Option<T>> {
    thread.map(|t| t.join().expect("the QC's threads do not panic")).transpose()
}

/// A text file's lines, gzipped (and BGZF) or not.
fn lines(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let gzipped = file.fill_buf().with_context(|| format!("reading {}", path.display()))?.starts_with(&[0x1f, 0x8b]);
    Ok(if gzipped { Box::new(BufReader::new(MultiGzDecoder::new(file))) } else { Box::new(file) })
}

/// Calls `each` on every line of `path` that is not a header (`#`), without its newline.
fn each_line(path: &Path, mut each: impl FnMut(&str) -> Result<()>) -> Result<()> {
    let mut reader = lines(path)?;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line).with_context(|| format!("reading {}", path.display()))? == 0 {
            return Ok(());
        }
        let l = line.trim_end_matches(['\n', '\r']);
        if !l.is_empty() && !l.starts_with('#') {
            each(l).with_context(|| format!("in {}", path.display()))?;
        }
    }
}

/// The read pass's statistics, with the unplaced unmapped reads of the index of `reads`, why the
/// index could not be read for them, and whether the reads are of some contigs or regions alone.
fn reads(path: &Path, reads: &Path) -> Result<(Reads, Option<String>, bool)> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let v: serde_json::Value = serde_json::from_str(&text).with_context(|| format!("reading {}", path.display()))?;
    let int = |key: &str| v[key].as_u64().with_context(|| format!("{}: no {key}", path.display()));
    let (n, secondary, supplementary) = (int("reads")?, int("secondary")?, int("supplementary")?);
    // Unknown when the index does not count the unplaced ones, or cannot be read; and not
    // measured for reads counted within regions, as the unplaced ones are not.
    let restricted = v["restricted"].as_bool() == Some(true);
    let counted = if restricted { Ok(None) } else { crate::sex::unplaced_unmapped(reads) };
    let (unplaced, error) = match counted {
        Ok(unplaced) => (unplaced, None),
        Err(e) => (None, Some(format!("the reads' unplaced unmapped reads cannot be counted: {e}"))),
    };
    let unmapped = unplaced.map(|u| u + int("unmapped_placed").unwrap_or(0));
    let alignments = n + secondary + supplementary;
    let reads = Reads {
        mapped_reads: n,
        mapped_bases: int("bases")?,
        n50: int("n50")?,
        mean_length: v["mean_length"].as_f64(),
        mean_quality: v["mean_quality"].as_f64(),
        unmapped,
        mapped_percent: unmapped.and_then(|u| percent(n, n + u)),
        secondary,
        secondary_percent: percent(secondary, alignments),
        supplementary,
        supplementary_percent: percent(supplementary, alignments),
    };
    Ok((reads, error, restricted))
}

/// `part` as a percentage of `whole`, to three decimals; none of nothing.
fn percent(part: u64, whole: u64) -> Option<f64> {
    (whole > 0).then(|| round(100.0 * part as f64 / whole as f64))
}

/// `a` over `b`, to three decimals; none over 0.
fn ratio(a: u64, b: u64) -> Option<f64> {
    (b > 0).then(|| round(a as f64 / b as f64))
}

fn round(x: f64) -> f64 {
    (x * 1000.0).round() / 1000.0
}

/// The autosome `contig` is, 1 to 22, with or without `chr`.
fn autosome(contig: &str) -> Option<u8> {
    let n: u8 = contig.strip_prefix("chr").unwrap_or(contig).parse().ok()?;
    (1..=22).contains(&n).then_some(n)
}

/// The mean depths of a mosdepth summary: its contig rows (not `_region` nor `total`), of
/// `chrom length bases mean min max`.
///
/// A summary of the regions alone (`lrpipe sv --depth-in-regions`, a run in regions) has only the
/// `CONTIG_region` rows, which are then each contig's depth within the regions.
fn coverage(path: &Path) -> Result<Coverage> {
    // (contig, length, bases, mean) of the contig rows, and of the region rows.
    let (mut contigs, mut regions) = (Vec::new(), Vec::new());
    each_line(path, |l| {
        let f: Vec<&str> = l.split('\t').collect();
        let [contig, len, sum, mean, ..] = f[..] else { bail!("not a mosdepth summary row: {l}") };
        if contig == "chrom" || contig.starts_with("total") {
            return Ok(());
        }
        let number = |s: &str| s.parse::<f64>().with_context(|| format!("not a number: {s}"));
        let row = (number(len)? as u64, number(sum)? as u64, number(mean)?);
        match contig.strip_suffix("_region") {
            Some(contig) => regions.push((contig.to_owned(), row)),
            None => contigs.push((contig.to_owned(), row)),
        }
        Ok(())
    })?;
    let rows = if contigs.is_empty() { regions } else { contigs };
    let (mut length, mut bases) = (0u64, 0u64);
    let mut coverage = Coverage { autosomes: None, chrx: None, chry: None, per_autosome: Means::default() };
    for (contig, (len, sum, mean)) in rows {
        if autosome(&contig).is_some() {
            length += len;
            bases += sum;
            coverage.per_autosome.0.push((contig, mean));
        } else if matches!(contig.as_str(), "chrX" | "X") {
            coverage.chrx = Some(mean);
        } else if matches!(contig.as_str(), "chrY" | "Y") {
            coverage.chry = Some(mean);
        }
    }
    coverage.autosomes = ratio(bases, length);
    Ok(coverage)
}

/// The sex the run called the sample with (`called`) against its reads' index's (over all of
/// `reads`), and the index's estimate.
fn sex_check(reads: &[PathBuf], called: &str) -> (SexCheck, Option<crate::sex::Estimate>) {
    if called == Sex::NotApplicable.as_str() {
        return (
            SexCheck {
                called: called.to_owned(),
                index: None,
                chrx_ratio: None,
                chry_ratio: None,
                agrees: None,
                note: None,
            },
            None,
        );
    }
    match crate::sex::estimate_all(reads) {
        Ok(e) => (
            SexCheck {
                called: called.to_owned(),
                index: Some(e.sex.map_or("unclear", Sex::as_str).to_owned()),
                chrx_ratio: Some(round(e.x)),
                chry_ratio: Some(round(e.y)),
                agrees: e.sex.map(|s| s.as_str() == called),
                note: e.sex.is_none().then(|| "chrX and chrY fit neither XX nor XY".to_owned()),
            },
            Some(e),
        ),
        Err(e) => (
            SexCheck {
                called: called.to_owned(),
                index: None,
                chrx_ratio: None,
                chry_ratio: None,
                agrees: None,
                // Without the reads' path, which is the run's own.
                note: Some(match e {
                    crate::sex::SexError::Partial { with, .. } => format!(
                        "the reads are on {} of the 22 autosomes{}: the index tells the sex from a whole genome's reads",
                        with.len(),
                        if with.is_empty() { String::new() } else { format!(" ({})", with.join(", ")) }
                    ),
                    crate::sex::SexError::Uneven { contig, ratio, .. } => format!(
                        "the reads cover {contig} at {ratio:.1} times the median autosome's depth: the index tells the \
                         sex from a whole genome's reads"
                    ),
                    other => other.to_string(),
                }),
            },
            None,
        ),
    }
}

/// The warning that the reads' index (`estimate`) gives another sex than the one the run calls
/// the sample with; none when they agree or the index fits neither.
#[must_use]
pub fn sex_warning(estimate: &crate::sex::Estimate, called: Sex) -> Option<String> {
    let index = estimate.sex.filter(|&s| s != called)?.as_str();
    Some(format!(
        "the reads' index says the sample is {index} (chrX at {:.3} and chrY at {:.3} of the autosomes' reads per \
         base), but it is called {}: if the index is right, run again with --sex {index} (or auto)",
        estimate.x,
        estimate.y,
        called.as_str()
    ))
}

/// What a genotype says of a call.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Zygosity {
    /// Two different alleles.
    Het,
    /// Two of one ALT allele.
    HomAlt,
    /// One ALT allele, haploid.
    Hemizygous,
    /// No ALT allele, or a missing one.
    Other,
}

/// The zygosity of the genotype `gt` (`0/1`, `1|1`, `1`, `./.`, ...).
fn zygosity(gt: &str) -> Zygosity {
    let alleles: Vec<Option<u32>> = gt.split(['/', '|']).map(|a| a.parse().ok()).collect();
    match alleles[..] {
        [Some(a)] if a > 0 => Zygosity::Hemizygous,
        [Some(a), Some(b)] if a != b => Zygosity::Het,
        [Some(a), Some(b)] if a > 0 && a == b => Zygosity::HomAlt,
        _ => Zygosity::Other,
    }
}

/// Whether the SNV `reference` > `alt` is a transition (A-G or C-T).
fn transition(reference: u8, alt: u8) -> bool {
    let pair = [reference.to_ascii_uppercase(), alt.to_ascii_uppercase()];
    matches!(pair, [b'A', b'G'] | [b'G', b'A'] | [b'C', b'T'] | [b'T', b'C'])
}

/// The first sample's genotype of a VCF record's fields.
fn genotype<'a>(fields: &[&'a str]) -> &'a str {
    let at = fields.get(8).and_then(|f| f.split(':').position(|k| k == "GT"));
    at.and_then(|i| fields.get(9)?.split(':').nth(i)).unwrap_or(".")
}

/// The PASS SNVs and indels of a VCF: a record is an SNV when its REF and every ALT are one base.
fn small_variants(path: &Path) -> Result<(Snvs, Indels)> {
    let (mut snvs, mut indels) = (Snvs::default(), Indels::default());
    each_line(path, |l| {
        let f: Vec<&str> = l.splitn(11, '\t').collect();
        if f.len() < 8 {
            bail!("not a VCF record: {l}");
        }
        if f[6] != "PASS" {
            return Ok(());
        }
        let alts: Vec<&str> = f[4].split(',').collect();
        let base = |a: &str| a.len() == 1 && a != "*" && a != ".";
        let zygosity = zygosity(genotype(&f));
        if base(f[3]) && alts.iter().all(|a| base(a)) {
            snvs.pass += 1;
            let r = f[3].as_bytes()[0];
            for alt in alts {
                if transition(r, alt.as_bytes()[0]) {
                    snvs.transitions += 1;
                } else {
                    snvs.transversions += 1;
                }
            }
            match zygosity {
                Zygosity::Het => snvs.het += 1,
                Zygosity::HomAlt => snvs.hom_alt += 1,
                Zygosity::Hemizygous => snvs.hemizygous += 1,
                Zygosity::Other => {}
            }
        } else {
            indels.pass += 1;
            match zygosity {
                Zygosity::Het => indels.het += 1,
                Zygosity::HomAlt => indels.hom_alt += 1,
                Zygosity::Hemizygous => indels.hemizygous += 1,
                Zygosity::Other => {}
            }
        }
        Ok(())
    })?;
    snvs.het_hom_ratio = ratio(snvs.het, snvs.hom_alt);
    snvs.ti_tv = ratio(snvs.transitions, snvs.transversions);
    indels.het_hom_ratio = ratio(indels.het, indels.hom_alt);
    Ok((snvs, indels))
}

/// The PASS SVs of a VCF, by their INFO's SVTYPE.
fn svs(path: &Path) -> Result<Svs> {
    let mut svs = Svs::default();
    each_line(path, |l| {
        let f: Vec<&str> = l.splitn(9, '\t').collect();
        if f.len() < 8 {
            bail!("not a VCF record: {l}");
        }
        if f[6] == "PASS" {
            svs.pass += 1;
            let svtype = f[7].split(';').find_map(|i| i.strip_prefix("SVTYPE=")).unwrap_or(".");
            *svs.by_type.entry(svtype.to_owned()).or_default() += 1;
        }
        Ok(())
    })?;
    Ok(svs)
}

/// The phasing of a phased VCF, by lrpipe-phase's `stats`.
fn phasing(path: &Path) -> Result<Phasing> {
    let name = path.to_str().with_context(|| format!("{} is not UTF-8", path.display()))?;
    let t = lrpipe_phase::stats::totals(name).with_context(|| format!("the phasing of {}", path.display()))?;
    let (het, phased) = (u64::try_from(t.heterozygous_snvs).unwrap_or(0), u64::try_from(t.phased_snvs).unwrap_or(0));
    Ok(Phasing {
        het_snvs_any_filter: t.heterozygous_snvs,
        phased_snvs: t.phased_snvs,
        phased_percent: percent(phased, het),
        blocks: t.blocks,
        block_n50: t.block_n50,
    })
}

/// The CpGs of a bedMethyl or CpG table: its rows' distinct positions (a bedMethyl row is of a
/// modification, 5mC or 5hmC, at a CpG; its rows are sorted).
fn cpgs(path: &Path) -> Result<u64> {
    let mut n = 0;
    let mut last = String::new();
    each_line(path, |l| {
        let mut f = l.splitn(3, '\t');
        let (contig, start) = (f.next().unwrap_or_default(), f.next().unwrap_or_default());
        let at = l.get(..contig.len() + 1 + start.len()).unwrap_or(l);
        if at != last {
            n += 1;
            last.clear();
            last.push_str(at);
        }
        Ok(())
    })?;
    Ok(n)
}

/// The summary as a table, `name<TAB>value` (NA when not measured): its values by section and
/// name, and a line for each note and warning.
#[must_use]
pub fn table(summary: &Summary) -> String {
    let mut rows: Vec<(String, String)> = vec![("sample".into(), summary.sample.clone())];
    let value = serde_json::to_value(summary).expect("a summary serialises");
    for key in ["regions", "bases"] {
        rows.push((format!("regions.{key}"), cell(&value["regions"][key])));
    }
    let mut add = |section: &str, keys: &[&str]| {
        for key in keys {
            let v = &value[section][key];
            rows.push((format!("{section}.{key}"), cell(v)));
        }
    };
    add(
        "reads",
        &[
            "mapped_reads",
            "mapped_bases",
            "n50",
            "mean_length",
            "mean_quality",
            "unmapped",
            "mapped_percent",
            "secondary",
            "secondary_percent",
            "supplementary",
            "supplementary_percent",
        ],
    );
    add("coverage", &["autosomes", "chrx", "chry"]);
    add("sex", &["called", "index", "chrx_ratio", "chry_ratio", "agrees"]);
    add("snvs", &["pass", "het", "hom_alt", "hemizygous", "het_hom_ratio", "transitions", "transversions", "ti_tv"]);
    add("indels", &["pass", "het", "hom_alt", "hemizygous", "het_hom_ratio"]);
    add("svs", &["pass"]);
    add("phasing", &["het_snvs_any_filter", "phased_snvs", "phased_percent", "blocks", "block_n50"]);
    add("methylation", &["combined", "hap1", "hap2"]);
    if let Some(coverage) = &summary.coverage {
        for (contig, mean) in &coverage.per_autosome.0 {
            rows.push((format!("coverage.{contig}"), mean.to_string()));
        }
    }
    if let Some(svs) = &summary.svs {
        for (svtype, n) in &svs.by_type {
            rows.push((format!("svs.{svtype}"), n.to_string()));
        }
    }
    rows.extend(summary.sex.note.iter().map(|n| ("sex.note".into(), n.clone())));
    rows.extend(summary.not_measured.iter().map(|n| ("not_measured".into(), n.clone())));
    rows.extend(summary.warnings.iter().map(|w| ("warning".into(), w.clone())));
    rows.into_iter().map(|(k, v)| format!("{k}\t{}\n", v.replace(['\t', '\n'], " "))).collect()
}

/// A JSON value as a table cell.
fn cell(v: &serde_json::Value) -> String {
    match v {
        serde_json::Value::Null => "NA".into(),
        serde_json::Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// What `lrpipe run` says of the summary at the end: its main values, a line or a few, and its
/// warnings.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Report {
    /// The values.
    pub lines: Vec<String>,
    /// The warnings.
    pub warnings: Vec<String>,
}

/// The [`Report`] of the summary at `path` (the qc stage's JSON).
///
/// # Errors
///
/// If it cannot be read, or is not JSON.
pub fn report(path: &Path) -> Result<Report> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let v: serde_json::Value = serde_json::from_str(&text).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = Vec::new();
    let int = |v: &serde_json::Value| v.as_u64().map(thousands);
    let sex = &v["sex"];
    let bacterial = sex["called"].as_str() == Some(Sex::NotApplicable.as_str());
    if let Some(n) = v["regions"]["regions"].as_u64() {
        let bases = bases(v["regions"]["bases"].as_u64().unwrap_or(0));
        let regions = if n == 1 { "region".to_owned() } else { format!("{} regions", thousands(n)) };
        let scope = if bacterial { "" } else { ", every value but the sex (the reads' index's)" };
        lines.push(format!("within the {regions} ({bases}){scope}"));
    }
    let float = |v: &serde_json::Value, digits: usize| v.as_f64().map(|x| format!("{x:.digits$}"));
    let r = &v["reads"];
    if let Some(reads) = int(&r["mapped_reads"]) {
        let yield_ = bases(r["mapped_bases"].as_u64().unwrap_or(0));
        let mut parts = vec![format!("{reads} mapped reads ({yield_})")];
        parts.push(format!("read N50 {}", bases(r["n50"].as_u64().unwrap_or(0))));
        parts.extend(r["mean_length"].as_f64().map(|m| format!("mean length {}", bases(m as u64))));
        parts.extend(float(&r["mean_quality"], 1).map(|q| format!("mean base quality {q} (Phred, per base)")));
        parts.extend(float(&r["mapped_percent"], 2).map(|p| format!("{p}% mapped")));
        let others: Vec<String> = [("supplementary", "supplementary_percent"), ("secondary", "secondary_percent")]
            .iter()
            .filter_map(|(name, key)| float(&r[*key], 2).map(|p| format!("{p}% {name}")))
            .collect();
        if !others.is_empty() {
            parts.push(format!("of the alignments {}", others.join(", ")));
        }
        lines.push(parts.join(", "));
    }
    let c = &v["coverage"];
    let mut parts = Vec::new();
    // A tenth of a read deep, or a hundredth below 10 (a targeted run's off-target depth).
    let depth =
        |v: &serde_json::Value| v.as_f64().map(|x| if x < 10.0 { format!("{x:.2}x") } else { format!("{x:.1}x") });
    if let Some(autosomes) = depth(&c["autosomes"]) {
        let mut coverage = format!("coverage {autosomes} over the autosomes");
        for (name, key) in [("chrX", "chrx"), ("chrY", "chry")] {
            coverage.extend(depth(&c[key]).map(|m| format!(", {name} {m}")));
        }
        parts.push(coverage);
    }
    if let Some(called) = sex["called"].as_str().filter(|_| !bacterial) {
        parts.push(match (sex["index"].as_str(), sex["agrees"].as_bool()) {
            (Some(_), Some(true)) => format!("sex {called}, as the index says"),
            (Some(index), Some(false)) => format!("sex {called}, but the index says {index}"),
            (Some(_), None) => format!("sex {called}; the index fits neither XX nor XY"),
            (None, _) => format!("sex {called}; the index cannot tell"),
        });
    }
    if !parts.is_empty() {
        lines.push(parts.join("; "));
    }
    let mut parts = Vec::new();
    let (snvs, indels) = (&v["snvs"], &v["indels"]);
    if let Some(n) = int(&snvs["pass"]) {
        let mut snv = format!("{n} SNVs");
        let detail: Vec<String> = [
            float(&snvs["ti_tv"], 2).map(|t| format!("Ti/Tv {t}")),
            float(&snvs["het_hom_ratio"], 2).map(|h| format!("het/hom {h}")),
        ]
        .into_iter()
        .flatten()
        .collect();
        if !detail.is_empty() {
            snv += &format!(" ({})", detail.join(", "));
        }
        parts.push(snv);
    }
    if let Some(n) = int(&indels["pass"]) {
        parts.push(match float(&indels["het_hom_ratio"], 2) {
            Some(h) => format!("{n} indels (het/hom {h})"),
            None => format!("{n} indels"),
        });
    }
    if let Some(n) = int(&v["svs"]["pass"]) {
        let types: Vec<String> = v["svs"]["by_type"]
            .as_object()
            .map(|t| t.iter().filter_map(|(k, n)| Some(format!("{k} {}", int(n)?))).collect())
            .unwrap_or_default();
        parts.push(if types.is_empty() { format!("{n} SVs") } else { format!("{n} SVs ({})", types.join(", ")) });
    }
    if !parts.is_empty() {
        lines.push(format!("PASS: {}", parts.join(", ")));
    }
    let mut parts = Vec::new();
    let p = &v["phasing"];
    if let Some(het) = int(&p["het_snvs_any_filter"]) {
        let mut phasing = match float(&p["phased_percent"], 1) {
            Some(percent) => format!("{percent}% of {het} het SNVs (any FILTER) phased"),
            None => format!("{het} het SNVs (any FILTER), none phased"),
        };
        phasing.extend(p["block_n50"].as_u64().map(|n| format!(", block N50 {}", bases(n))));
        parts.push(phasing);
    }
    let m = &v["methylation"];
    if let Some(cpgs) = int(&m["combined"]) {
        let haplotypes: Vec<String> =
            ["hap1", "hap2"].iter().filter_map(|h| Some(format!("{h} {}", int(&m[*h])?))).collect();
        parts.push(if haplotypes.is_empty() {
            format!("{cpgs} CpGs")
        } else {
            format!("{cpgs} CpGs ({})", haplotypes.join(", "))
        });
    }
    if !parts.is_empty() {
        lines.push(parts.join("; "));
    }
    let warnings = v["warnings"].as_array().into_iter().flatten().filter_map(|w| w.as_str().map(str::to_owned));
    Ok(Report { lines, warnings: warnings.collect() })
}

/// `n` with commas between thousands.
fn thousands(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::new();
    for (i, d) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(d);
    }
    out
}

/// A number of bases in b, kb, Mb or Gb.
fn bases(n: u64) -> String {
    match n {
        0..1_000 => format!("{n} b"),
        1_000..1_000_000 => format!("{:.1} kb", n as f64 / 1e3),
        1_000_000..1_000_000_000 => format!("{:.2} Mb", n as f64 / 1e6),
        _ => format!("{:.2} Gb", n as f64 / 1e9),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A scratch file, removed when dropped.
    struct Temp(PathBuf);

    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    /// A scratch file named `name` holding `bytes`.
    fn temp(name: &str, bytes: impl AsRef<[u8]>) -> Temp {
        let path = std::env::temp_dir().join(format!("lrpipe-qc-{}-{name}", std::process::id()));
        std::fs::write(&path, bytes).expect("written");
        Temp(path)
    }

    #[test]
    fn ti_tv_and_het_hom_count_pass_records_by_their_alleles() {
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
        let records = [
            // Transitions: A>G (het), C>T (hom), T>C (haploid).
            "chr1\t1\t.\tA\tG\t30\tPASS\t.\tGT:GQ\t0/1:20",
            "chr1\t2\t.\tc\tt\t30\tPASS\t.\tGT\t1|1",
            "chrX\t3\t.\tT\tC\t30\tPASS\t.\tGT\t1",
            // A transversion, and a site of two ALTs: one transition, one transversion.
            "chr1\t4\t.\tA\tC\t30\tPASS\t.\tGT\t0|1",
            "chr1\t5\t.\tG\tA,T\t30\tPASS\t.\tGT\t1/2",
            // Not PASS: left out.
            "chr1\t6\t.\tA\tG\t3\tLowQual\t.\tGT\t0/1",
            "chr1\t7\t.\tA\t.\t3\tRefCall\t.\tGT\t0/0",
            // Indels: an insertion (hom), a deletion (het), and one beside an SNV allele.
            "chr1\t8\t.\tA\tAT\t30\tPASS\t.\tGT\t1/1",
            "chr1\t9\t.\tAT\tA\t30\tPASS\t.\tGT\t0/1",
            "chr1\t10\t.\tA\tG,AT\t30\tPASS\t.\tGT\t1/2",
            // A missing genotype counts as a record only.
            "chr1\t11\t.\tA\tT\t30\tPASS\t.\tGT\t./.",
        ];
        let file = temp("small.vcf", format!("{header}{}\n", records.join("\n")));
        let (snvs, indels) = small_variants(&file.0).expect("read");
        assert_eq!((snvs.pass, snvs.het, snvs.hom_alt, snvs.hemizygous), (6, 3, 1, 1));
        assert_eq!((snvs.transitions, snvs.transversions), (4, 3));
        assert_eq!(snvs.ti_tv, Some(1.333));
        assert_eq!(snvs.het_hom_ratio, Some(3.0));
        assert_eq!((indels.pass, indels.het, indels.hom_alt), (3, 2, 1));
        assert_eq!(indels.het_hom_ratio, Some(2.0));
        // No transversions, no ratio.
        let file = temp("ti.vcf", format!("{header}{}\n", records[0]));
        assert_eq!(small_variants(&file.0).expect("read").0.ti_tv, None);
    }

    #[test]
    fn genotypes_are_read_from_their_format_field() {
        assert_eq!(zygosity("0/1"), Zygosity::Het);
        assert_eq!(zygosity("1|2"), Zygosity::Het);
        assert_eq!(zygosity("2/2"), Zygosity::HomAlt);
        assert_eq!(zygosity("1"), Zygosity::Hemizygous);
        assert_eq!(zygosity("0"), Zygosity::Other);
        assert_eq!(zygosity("0/0"), Zygosity::Other);
        assert_eq!(zygosity("./1"), Zygosity::Other);
        assert_eq!(genotype(&["c", "1", ".", "A", "G", ".", "PASS", ".", "GQ:GT", "40:1/1"]), "1/1");
        assert_eq!(genotype(&["c", "1", ".", "A", "G", ".", "PASS", "."]), ".");
    }

    #[test]
    fn coverage_is_of_the_contig_rows_with_the_autosomes_by_their_bases() {
        let summary = "chrom\tlength\tbases\tmean\tmin\tmax\n\
                       chr1\t100\t3000\t30.00\t0\t50\nchr1_region\t100\t3000\t30.00\t0\t50\n\
                       chr2\t300\t3000\t10.00\t0\t50\nchrX\t100\t1500\t15.00\t0\t40\n\
                       chrY\t50\t500\t10.00\t0\t20\nchrM\t10\t10000\t1000.00\t0\t2000\n\
                       total\t560\t21000\t37.5\t0\t2000\n";
        let c = coverage(&temp("depth.txt", summary).0).expect("read");
        assert_eq!(c.autosomes, Some(15.0), "6000 bases over 400");
        assert_eq!((c.chrx, c.chry), (Some(15.0), Some(10.0)));
        assert_eq!(c.per_autosome.0, [("chr1".to_owned(), 30.0), ("chr2".to_owned(), 10.0)]);
        assert_eq!(serde_json::to_string(&c.per_autosome).expect("json"), r#"{"chr1":30.0,"chr2":10.0}"#);
    }

    #[test]
    fn a_summary_of_regions_alone_gives_their_coverage() {
        // `lrpipe sv --depth-in-regions`: each contig's depth within the regions.
        let summary = "chrom\tlength\tbases\tmean\tmin\tmax\n\
                       chr1_region\t100\t3000\t30.00\t0\t50\nchr2_region\t300\t3000\t10.00\t0\t50\n\
                       chrX_region\t100\t1500\t15.00\t0\t40\ntotal_region\t500\t7500\t15.00\t0\t50\n";
        let c = coverage(&temp("regions.txt", summary).0).expect("read");
        assert_eq!(c.autosomes, Some(15.0), "6000 bases over 400");
        assert_eq!((c.chrx, c.chry), (Some(15.0), None));
        assert_eq!(c.per_autosome.0, [("chr1".to_owned(), 30.0), ("chr2".to_owned(), 10.0)]);
    }

    #[test]
    fn svs_and_cpgs_are_counted() {
        let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   chr1\t1\ta\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;SVLEN=-50\n\
                   chr1\t2\tb\tN\t<INS>\t60\tPASS\tSVTYPE=INS\n\
                   chr1\t3\tc\tN\t<DEL>\t60\tGT\tSVTYPE=DEL\n\
                   chr1\t4\td\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL\n";
        let svs = svs(&temp("sv.vcf", vcf).0).expect("read");
        assert_eq!(svs.pass, 3);
        assert_eq!(svs.by_type.into_iter().collect::<Vec<_>>(), [("DEL".to_owned(), 2), ("INS".to_owned(), 1)]);
        // bedMethyl: a row for each modification at a CpG.
        let bed = "chr1\t10\t11\tm\t5\t+\nchr1\t10\t11\th\t5\t+\nchr1\t100\t101\tm\t3\t+\nchr2\t10\t11\tm\t3\t+\n";
        assert_eq!(cpgs(&temp("m.bed", bed).0).expect("read"), 3);
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut e, format!("#header\n{bed}").as_bytes()).expect("written");
        let gz = temp("m.bed.gz", e.finish().expect("finished"));
        assert_eq!(cpgs(&gz.0).expect("read"), 3, "gzipped, with a header");
    }

    #[test]
    fn several_reads_give_the_sex_from_each_as_the_run_estimates_it() {
        // One estimate, the run's at its start: over the files merged, not the merged BAM.
        let site = Site::parse("reference = \"/ref.fa\"\n[snv]\n", Path::new("/")).expect("a site");
        let inputs: Vec<PathBuf> = ["/reads/a.bam", "/reads/b.cram"].map(PathBuf::from).to_vec();
        let sample = Sample {
            name: "S".into(),
            reads: "/out/reads/S.bam".into(),
            unaligned: Vec::new(),
            unaligned_fastq: false,
            merged: Some(crate::plan::Merged { inputs: inputs.clone(), index: crate::plan::BamIndex::Bai }),
            regions: None,
            platform: crate::plan::Platform::Ont,
            sex: Sex::Male,
            model: None,
        };
        let qc = stage(&site, &sample, &Layout::new("/out", "S"), &[]);
        let words: Vec<String> = qc.args.iter().map(|a| a.to_string_lossy().into_owned()).collect();
        let words = words.join(" ");
        assert!(
            words.contains("--reads /out/reads/S.bam --sex-reads /reads/a.bam --sex-reads /reads/b.cram"),
            "{words}"
        );
        assert!(inputs.iter().all(|i| qc.inputs.contains(i)), "{:?}", qc.inputs);

        // In regions: of them, which are its input, so that other regions rerun it.
        let regions = Sample { regions: Some(crate::plan::Targets::Bed), merged: None, ..sample };
        let qc = stage(&site, &regions, &Layout::new("/out", "S"), &[]);
        let words: Vec<String> = qc.args.iter().map(|a| a.to_string_lossy().into_owned()).collect();
        assert!(words.join(" ").contains("--regions /out/.lrpipe/regions.bed"), "{words:?}");
        assert!(qc.inputs.contains(&PathBuf::from("/out/.lrpipe/regions.bed")), "{:?}", qc.inputs);
    }

    #[test]
    fn numbers_read_as_people_write_them() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_000), "1,000");
        assert_eq!(thousands(4_512_331), "4,512,331");
        assert_eq!(bases(950), "950 b");
        assert_eq!(bases(23_456), "23.5 kb");
        assert_eq!(bases(1_200_000), "1.20 Mb");
        assert_eq!(bases(123_400_000_000), "123.40 Gb");
    }
}

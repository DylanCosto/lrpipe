//! Port of `cli/haplotag.py`: tag reads by haplotype.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Write;
use std::sync::mpsc::{channel, sync_channel};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use crate::alignment::{self, AlignmentFile, Aux, Fetch, Record};
use anyhow::{Context, Result, anyhow, bail};

use crate::compat::pystr::py_int;
use crate::core::read::{Read, ReadSet};
use crate::hts::{IndexedVcf, VcfFetchError};
use crate::phase::ReferenceCache;
use crate::variants::{PRIMARY_DEFAULT_SUB_ALIGNMENT_ID, ReadSetError, ReadSetReader, is_primary, sub_alignment_id};
use crate::vcf::{VariantCallPhase, VariantTable, VcfReader, VcfVariant};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SupplementaryStrategy {
    Skip,
    CopyPrimary,
    IndependentOrSkip,
    IndependentOrCopyPrimary,
}

impl SupplementaryStrategy {
    pub fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "skip" => Self::Skip,
            "copy-primary" => Self::CopyPrimary,
            "independent-or-skip" => Self::IndependentOrSkip,
            "independent-or-copy-primary" => Self::IndependentOrCopyPrimary,
            _ => bail!("invalid supplementary strategy '{s}'"),
        })
    }

    fn consider_supplementary(self) -> bool {
        self != Self::Skip
    }

    fn attempt_to_haplotag_each_supplementary(self) -> bool {
        matches!(self, Self::IndependentOrSkip | Self::IndependentOrCopyPrimary)
    }

    fn attempt_to_copy_primary(self) -> bool {
        matches!(self, Self::CopyPrimary | Self::IndependentOrCopyPrimary)
    }
}

pub struct HaplotagOptions {
    pub variant_file: String,
    pub alignment_file: String,
    pub output: Option<String>,
    pub reference: Option<String>,
    pub no_reference: bool,
    pub regions: Option<Vec<String>>,
    /// The regions of a BED file in place of `regions`, sorted and merged; a read reaching
    /// several is written once (lrpipe-phase's own option)
    pub regions_bed: Option<String>,
    pub ignore_linked_read: bool,
    pub linked_read_distance_cutoff: i64,
    pub ignore_read_groups: bool,
    pub given_samples: Option<Vec<String>>,
    pub haplotag_list: Option<String>,
    pub supplementary_strategy: SupplementaryStrategy,
    pub supplementary_distance_threshold: i64,
    pub supplementary_strand_match: bool,
    pub ploidy: usize,
    pub skip_missing_contigs: bool,
    /// Reads on contigs the VCF does not have are written (and listed) untagged, as on a contig
    /// without variants, instead of failing the run (lrpipe-phase's own option)
    pub untag_missing_contigs: bool,
    /// Only the haplotag list is written, no alignments (lrpipe-phase's own option)
    pub list_only: bool,
    /// The output BAM is indexed (`.bai`, or `.csi` for contigs too long for it), as `samtools
    /// index` would (lrpipe-phase's own option)
    pub write_index: bool,
    /// Threads for compressing the output (None: use `threads`)
    pub output_threads: Option<usize>,
    /// Output compression level (None: htslib default)
    pub compression_level: Option<u32>,
    /// Number of chromosomes to process in parallel
    pub threads: usize,
    /// Value of the CL field of the @PG header line
    pub command_line: String,
}

/// `ReadAlignmentRepresentation`
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct ReadRepr {
    read_name: String,
    chromosome: String,
    is_supplementary: bool,
    sub_alignment_id: String,
}

#[derive(Clone, Copy, Debug)]
struct PrimaryInfo {
    reference_start: i64,
    reference_end: i64,
    is_reverse: bool,
}

/// `read_representation`
fn read_representation(read: &Read, as_primary: bool) -> ReadRepr {
    let is_supplementary = if as_primary { false } else { read.is_supplementary };
    let sub_alignment_id =
        if as_primary { PRIMARY_DEFAULT_SUB_ALIGNMENT_ID.to_string() } else { read.sub_alignment_id.clone() };
    let mut read_name = read.name.clone();
    if read_name.ends_with(&sub_alignment_id) {
        read_name.truncate(read_name.len() - read.sub_alignment_id.len());
    }
    ReadRepr { read_name, chromosome: read.chromosome.clone(), is_supplementary, sub_alignment_id }
}

/// `alignment_representation`
fn alignment_representation(rec: &Record, chromosome: &str, as_primary: bool) -> ReadRepr {
    let primary = as_primary || is_primary(rec);
    ReadRepr {
        read_name: String::from_utf8_lossy(rec.qname()).into_owned(),
        chromosome: chromosome.to_string(),
        is_supplementary: !primary,
        sub_alignment_id: if primary { PRIMARY_DEFAULT_SUB_ALIGNMENT_ID.to_string() } else { sub_alignment_id(rec) },
    }
}

fn min_alignment_distance(first_start: i64, first_end: i64, second_start: i64, second_end: i64) -> i64 {
    let (a, b) = if first_start < first_end { (first_start, first_end) } else { (first_end, first_start) };
    let (c, d) = if second_start < second_end { (second_start, second_end) } else { (second_end, second_start) };
    if a > d {
        return a - d;
    }
    if c > b {
        return c - b;
    }
    0
}

/// Variant position -> (phase set, alleles in haplotype order)
type PhaseInfo = HashMap<i32, (i64, Vec<Option<u32>>)>;

/// Regions of one chromosome: (start, end or None for "to the end")
type Regions = Vec<(i64, Option<i64>)>;

/// `get_variant_information`
fn get_variant_information(table: &VariantTable, sample: &str) -> Result<(PhaseInfo, Vec<VcfVariant>)> {
    let mut vpos_to_phase_info = HashMap::new();
    let mut variants = Vec::new();
    for ((v, gt), phase) in table.variants.iter().zip(table.genotypes_of(sample)).zip(table.phases_of(sample)) {
        let Some(VariantCallPhase { block_id: Some(block_id), phase, .. }) = phase else { continue };
        vpos_to_phase_info.insert(v.position, (block_id.to_int()?, phase.clone()));
        if !gt.is_homozygous() {
            variants.push(v.clone());
        }
    }
    Ok((vpos_to_phase_info, variants))
}

struct HaplotagInfo {
    bx_tag_to_haplotype: HashMap<String, Vec<(i32, usize, i64)>>,
    read_to_haplotype: HashMap<ReadRepr, (usize, i64, i64)>,
    n_multiple_phase_sets: usize,
    primary_info_by_repr: HashMap<ReadRepr, PrimaryInfo>,
}

impl HaplotagInfo {
    /// No read of the chromosome has a haplotype (`--untag-missing-contigs`).
    fn untagged() -> Self {
        HaplotagInfo {
            bx_tag_to_haplotype: HashMap::new(),
            read_to_haplotype: HashMap::new(),
            n_multiple_phase_sets: 0,
            primary_info_by_repr: HashMap::new(),
        }
    }
}

/// `utils.Region.parse`
fn parse_region(spec: &str) -> Result<(String, i64, Option<i64>)> {
    let mut parts = spec.splitn(2, ':');
    let chromosome = parts.next().unwrap().to_string();
    let rest = parts.next().unwrap_or("");
    if rest.is_empty() {
        return Ok((chromosome, 0, None));
    }
    let sep = if rest.contains(':') { ':' } else { '-' };
    let mut se = rest.splitn(2, sep);
    let invalid = || anyhow!("Region must be specified as chrom[:start[-end]])");
    let int = |s: &str| py_int(s).and_then(|v| i64::try_from(v).ok()).ok_or_else(invalid);
    let start = int(se.next().unwrap())? - 1;
    let end = match se.next() {
        None | Some("") => None,
        Some(e) => {
            let e = int(e)?;
            if e <= start {
                bail!("end is before start in specified region");
            }
            Some(e)
        }
    };
    Ok((chromosome, start, end))
}

/// Python's `str.isupper()` (ASCII)
fn py_isupper(s: &str) -> bool {
    let mut cased = false;
    for c in s.chars() {
        if c.is_lowercase() {
            return false;
        }
        if c.is_uppercase() {
            cased = true;
        }
    }
    cased
}

const VALID_HEADERS: [&str; 5] = ["HD", "SQ", "RG", "PG", "CO"];

fn valid_header_order(record: &str) -> &'static [&'static str] {
    match record {
        "HD" => &["VN", "SO", "SS", "GO"],
        "SQ" => &["SN", "LN", "AS", "M5", "UR", "SP", "AH", "TP", "DS", "AN"],
        "RG" => &["ID", "CN", "SM", "LB", "PU", "PI", "DT", "DS", "PL", "FO", "KS", "PG", "PM", "BC"],
        "PG" => &["PN", "ID", "VN", "PP", "DS", "CL"],
        _ => &[],
    }
}

fn known_header_field(record: &str, key: &str) -> bool {
    let fields: &[&str] = match record {
        "HD" => &["VN", "SO", "GO", "SS"],
        "SQ" => &["SN", "LN", "AS", "M5", "SP", "UR", "AH", "TP", "DS", "AN"],
        "RG" => &["ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", "PL", "PM", "PU", "SM", "BC"],
        "PG" => &["ID", "PN", "CL", "PP", "DS", "VN"],
        _ => &[],
    };
    fields.contains(&key)
}

/// An insertion-ordered dict of header fields (Python dict semantics).
type Fields = Vec<(String, String)>;

fn fields_set(fields: &mut Fields, key: &str, value: String) {
    if let Some(e) = fields.iter_mut().find(|(k, _)| k == key) {
        e.1 = value;
    } else {
        fields.push((key.to_string(), value));
    }
}

/// pysam `AlignmentHeader.to_dict()`: record type -> list of field dicts
/// (CO lines are stored as a single ("", text) entry).
fn header_to_dict(text: &str, references: &[(String, u64)]) -> Result<Vec<(String, Vec<Fields>)>> {
    let mut text = text.to_string();
    if !text.contains("@SQ") {
        text.push('\n');
        for (name, len) in references {
            text.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
        }
    }
    let mut result: Vec<(String, Vec<Fields>)> = Vec::new();
    for line in text.split('\n') {
        let line = line.trim_matches(|c| c == ' ' || c == '\0');
        if line.is_empty() {
            continue;
        }
        if !line.starts_with('@') {
            bail!("header line without '@': '{line}'");
        }
        let fields: Vec<&str> = line[1..].split('\t').collect();
        let record = fields[0];
        if !VALID_HEADERS.contains(&record) {
            bail!("header line with invalid type '{record}': '{line}'");
        }
        let entry_index = match result.iter().position(|(r, _)| r == record) {
            Some(i) => i,
            None => {
                result.push((record.to_string(), Vec::new()));
                result.len() - 1
            }
        };
        if record == "CO" {
            result[entry_index].1.push(vec![(String::new(), fields[1..].join("\t"))]);
            continue;
        }
        let mut x: Fields = Vec::new();
        for (idx, field) in fields[1..].iter().enumerate() {
            let Some((key, value)) = field.split_once(':') else { bail!("malformatted header: no ':' in field") };
            if key == "CL" {
                let rest = fields[idx + 1..].join("\t");
                let (key, value) = rest.split_once(':').unwrap();
                if !known_header_field(record, key) {
                    bail!("KeyError: '{key}'");
                }
                fields_set(&mut x, key, value.to_string());
                break;
            }
            let value = if record == "SQ" && key == "LN" {
                value
                    .trim()
                    .parse::<i64>()
                    .map_err(|_| anyhow!("invalid literal for int() with base 10: '{value}'"))?
                    .to_string()
            } else {
                value.to_string()
            };
            fields_set(&mut x, key, value);
        }
        if record == "HD" && !result[entry_index].1.is_empty() {
            bail!("multiple 'HD' lines are not permitted");
        }
        result[entry_index].1.push(x);
    }
    if !result.iter().any(|(r, _)| r == "SQ") {
        let sq = references
            .iter()
            .map(|(n, l)| vec![("LN".to_string(), l.to_string()), ("SN".to_string(), n.clone())])
            .collect();
        result.push(("SQ".to_string(), sq));
    }
    Ok(result)
}

/// pysam `build_header_line`
fn build_header_line(fields: &Fields, record: &str) -> String {
    let mut line = vec![format!("@{record}")];
    if record == "CO" {
        line.push(fields[0].1.clone());
    } else {
        for key in valid_header_order(record) {
            if let Some((_, v)) = fields.iter().find(|(k, _)| k == key) {
                line.push(format!("{key}:{v}"));
            }
        }
        for (k, v) in fields {
            if !py_isupper(k) {
                line.push(format!("{k}:{v}"));
            }
        }
    }
    line.join("\t")
}

/// pysam `AlignmentHeader.from_dict()`: header text and references.
fn header_from_dict(dict: &[(String, Vec<Fields>)]) -> (String, Vec<(String, u32)>) {
    let mut lines = Vec::new();
    for record in VALID_HEADERS {
        if let Some((_, entries)) = dict.iter().find(|(r, _)| r == record) {
            for fields in entries {
                lines.push(build_header_line(fields, record));
            }
        }
    }
    let text = lines.join("\n") + "\n";
    let mut refs = Vec::new();
    if let Some((_, sq)) = dict.iter().find(|(r, _)| r == "SQ") {
        for fields in sq {
            let get = |k: &str| fields.iter().find(|(key, _)| key == k).map(|(_, v)| v.clone());
            refs.push((get("SN").unwrap_or_default(), get("LN").and_then(|l| l.parse().ok()).unwrap_or(0)));
        }
    }
    (text, refs)
}

/// Opens the output like pysam's `AlignmentFile(..., header=AlignmentHeader.from_dict(...))`:
/// `mode` is "wb", "wb<level>" or "wc<level>", `path` "-" is standard output.
fn create_output(
    path: &str,
    mode: &str,
    text: &str,
    refs: &[(String, u32)],
    reference: Option<&str>,
    threads: usize,
) -> Result<alignment::AlignmentWriter> {
    let header = alignment::Header { text: text.to_string(), refs: refs.to_vec() };
    let level = mode[2..].parse::<u32>().ok();
    if mode.starts_with("wc") {
        let reference = reference.ok_or_else(|| anyhow!("CRAM output needs a reference"))?;
        return Ok(alignment::AlignmentWriter::Cram(Box::new(alignment::CramWriter::create(
            std::path::Path::new(path),
            &header,
            std::path::Path::new(reference),
        )?)));
    }
    let out: Box<dyn Write + Send> = if path == "-" {
        Box::new(std::io::stdout())
    } else {
        Box::new(
            std::fs::File::create(path)
                .with_context(|| format!("Error while initializing alignment output file at path: {path}"))?,
        )
    };
    Ok(alignment::AlignmentWriter::Bam(alignment::bam::BamWriter::create(out, &header, level, threads)?))
}

/// pysam `AlignedSegment.set_tag(tag, value)` for integers (or None to delete).
fn set_int_tag(rec: &mut Record, tag: &[u8; 2], value: Option<i64>) -> Result<()> {
    rec.remove_aux(tag);
    if let Some(v) = value {
        rec.push_aux_int(tag, v)?;
    }
    Ok(())
}

fn bx_tag(rec: &Record) -> Option<String> {
    match rec.aux(b"BX") {
        Ok(Aux::String(s)) => Some(s.to_string()),
        Ok(other) => Some(format!("{other:?}")),
        Err(_) => None,
    }
}

struct Haplotagger<'a> {
    opts: &'a HaplotagOptions,
    readset_reader: ReadSetReader,
    fasta: Option<ReferenceCache>,
}

impl Haplotagger<'_> {
    /// `PhasedInputReader.read` (without phased VCF inputs)
    fn read(
        &mut self,
        chromosome: &str,
        variants: &[VcfVariant],
        sample: &str,
        numeric_sample_id: i32,
        regions: &[(i64, Option<i64>)],
    ) -> Result<ReadSet> {
        let reference = match self.fasta.as_mut() {
            Some(f) => Some(f.get(chromosome)?),
            None => None,
        };
        let bam_sample = if self.opts.ignore_read_groups { None } else { Some(sample) };
        let mut readset = match self.readset_reader.read_regions(
            chromosome,
            variants,
            bam_sample,
            numeric_sample_id,
            reference,
            regions,
        ) {
            Ok(rs) => rs,
            Err(ReadSetError::SampleNotFound) => {
                log::warn!("Sample '{}' not found in any BAM/CRAM file.", bam_sample.unwrap_or("None"));
                ReadSet::new()
            }
            Err(ReadSetError::ReferenceNotFound(_)) => {
                bail!("The chromosome '{chromosome}' was not found in the BAM/CRAM file.")
            }
            Err(ReadSetError::Other(e)) => return Err(e),
        };
        for read in readset.iter_mut() {
            read.sort()?;
        }
        readset.sort();
        log::info!("Found {} reads covering {} variants", readset.len(), readset.get_positions().len());
        Ok(readset)
    }

    /// `prepare_haplotag_information`
    fn prepare(
        &mut self,
        table: &VariantTable,
        shared_samples: &[String],
        sample_ids: &mut HashMap<String, i32>,
        regions: &[(i64, Option<i64>)],
    ) -> Result<HaplotagInfo> {
        let ploidy = self.opts.ploidy;
        let ignore_linked_read = self.opts.ignore_linked_read;
        let cutoff = self.opts.linked_read_distance_cutoff;
        let mut info = HaplotagInfo {
            bx_tag_to_haplotype: HashMap::new(),
            read_to_haplotype: HashMap::new(),
            n_multiple_phase_sets: 0,
            primary_info_by_repr: HashMap::new(),
        };
        for sample in shared_samples {
            let (vpos_to_phase_info, variants) = get_variant_information(table, sample)?;
            let numeric_id = if self.opts.ignore_read_groups {
                0
            } else {
                let n = sample_ids.len() as i32;
                *sample_ids.entry(sample.clone()).or_insert(n)
            };
            let read_set = self.read(&table.chromosome, &variants, sample, numeric_id, regions)?;
            let mut bx_tag_to_readlist: HashMap<&str, Vec<usize>> = HashMap::new();
            if !ignore_linked_read {
                for (i, read) in read_set.iter().enumerate() {
                    if !read.bx_tag.is_empty() {
                        bx_tag_to_readlist.entry(read.bx_tag.as_str()).or_default().push(i);
                    }
                }
            }
            for read in read_set.iter() {
                if !read.is_supplementary {
                    info.primary_info_by_repr.insert(
                        read_representation(read, true),
                        PrimaryInfo {
                            reference_start: read.reference_start as i64,
                            reference_end: read.reference_end as i64,
                            is_reverse: read.is_reverse,
                        },
                    );
                }
            }
            let mut processed: HashSet<ReadRepr> = HashSet::new();
            for (ri, read) in read_set.iter().enumerate() {
                let repr = read_representation(read, false);
                if processed.contains(&repr) {
                    continue;
                }
                processed.insert(repr);
                let mut reads_to_consider = vec![ri];
                if !ignore_linked_read && !read.bx_tag.is_empty() {
                    for &j in &bx_tag_to_readlist[read.bx_tag.as_str()] {
                        let r = read_set.get(j);
                        if !processed.contains(&read_representation(r, false))
                            && (read.reference_start as i64 - r.reference_start as i64).abs() <= cutoff
                            && !reads_to_consider.contains(&j)
                        {
                            reads_to_consider.push(j);
                        }
                    }
                }
                // phaseset -> costs of assigning the read(s) to each haplotype,
                // in insertion order (Python dict)
                let mut haplotype_costs: Vec<(i64, Vec<i64>)> = Vec::new();
                for &j in &reads_to_consider {
                    let r = read_set.get(j);
                    processed.insert(read_representation(r, false));
                    for v in &r.variants {
                        if v.allele != 0 && v.allele != 1 {
                            bail!("AssertionError: allele {} not in [0, 1]", v.allele);
                        }
                        let (phaseset, phasing) = &vpos_to_phase_info[&v.position];
                        let k = match haplotype_costs.iter().position(|(ps, _)| ps == phaseset) {
                            Some(k) => k,
                            None => {
                                haplotype_costs.push((*phaseset, vec![0; ploidy]));
                                haplotype_costs.len() - 1
                            }
                        };
                        for (hap_index, hap_allele) in phasing.iter().enumerate() {
                            if Some(v.allele as u32) == *hap_allele {
                                haplotype_costs[k].1[hap_index] += v.quality as i64;
                            }
                        }
                    }
                }
                // stable sort by maximum score, descending
                haplotype_costs.sort_by(|a, b| b.1.iter().max().cmp(&a.1.iter().max()));
                if haplotype_costs.is_empty() {
                    continue;
                }
                if haplotype_costs.len() > 1 {
                    info.n_multiple_phase_sets += 1;
                }
                let (phaseset, scores) = &haplotype_costs[0];
                let mut scores_list: Vec<(usize, i64)> = scores.iter().copied().enumerate().collect();
                scores_list.sort_by_key(|t| std::cmp::Reverse(t.1));
                let (first_ht, first_score) = scores_list[0];
                let (_second_ht, second_score) = scores_list[1];
                let quality = first_score - second_score;
                if quality == 0 {
                    continue;
                }
                if !ignore_linked_read && !read.bx_tag.is_empty() {
                    info.bx_tag_to_haplotype.entry(read.bx_tag.clone()).or_default().push((
                        read.reference_start,
                        first_ht,
                        *phaseset,
                    ));
                }
                for &j in &reads_to_consider {
                    let r = read_set.get(j);
                    info.read_to_haplotype.insert(read_representation(r, false), (first_ht, quality, *phaseset));
                    log::debug!(
                        "Assigned read {} to haplotype {} with a quality of {} based on {} covered variants",
                        r.name,
                        first_ht,
                        quality,
                        r.len()
                    );
                }
            }
        }
        Ok(info)
    }
}

/// `attempt_add_phase_information`: returns (is_tagged, haplotype_name, phaseset)
fn attempt_add_phase_information(
    rec: &mut Record,
    chromosome: &str,
    info: &HaplotagInfo,
    opts: &HaplotagOptions,
) -> Result<(bool, String, String)> {
    let strategy = opts.supplementary_strategy;
    let mut is_tagged = false;
    let mut haplotype_name = "none".to_string();
    let mut phaseset = "none".to_string();
    let is_supplementary = rec.is_supplementary();
    if is_supplementary && strategy == SupplementaryStrategy::Skip {
        return Ok((is_tagged, haplotype_name, phaseset));
    }
    let mut representations =
        vec![alignment_representation(rec, chromosome, false), alignment_representation(rec, chromosome, true)];
    match strategy {
        SupplementaryStrategy::CopyPrimary => {
            representations.remove(0);
        }
        SupplementaryStrategy::IndependentOrSkip => {
            representations.remove(1);
        }
        _ => {}
    }
    if is_supplementary
        && strategy.attempt_to_copy_primary()
        && let Some(pi) = info.primary_info_by_repr.get(&alignment_representation(rec, chromosome, true))
    {
        let end = rec.cigar().end_pos();
        let mut remove = opts.supplementary_distance_threshold
            < min_alignment_distance(pi.reference_start, pi.reference_end, rec.pos(), end);
        if opts.supplementary_strand_match {
            remove |= pi.is_reverse != rec.is_reverse();
        }
        if remove {
            representations.pop();
        }
    }
    let mut found = false;
    for repr in &representations {
        if let Some(&(haplotype, quality, ps)) = info.read_to_haplotype.get(repr) {
            haplotype_name = format!("H{}", haplotype + 1);
            set_int_tag(rec, b"HP", Some(haplotype as i64 + 1))?;
            set_int_tag(rec, b"PC", Some(quality))?;
            set_int_tag(rec, b"PS", Some(ps))?;
            phaseset = ps.to_string();
            is_tagged = true;
            found = true;
            break;
        }
    }
    if !found && !opts.ignore_linked_read {
        let read_clouds: &[(i32, usize, i64)] = match bx_tag(rec) {
            None => &[],
            Some(tag) => info.bx_tag_to_haplotype.get(&tag).map(|v| v.as_slice()).unwrap_or(&[]),
        };
        for &(reference_start, haplotype, ps) in read_clouds {
            // Python's loop variable `phaseset` keeps the last cloud's value
            // even if no cloud is close enough.
            phaseset = ps.to_string();
            if (reference_start as i64 - rec.pos()).abs() <= opts.linked_read_distance_cutoff {
                haplotype_name = format!("H{}", haplotype + 1);
                set_int_tag(rec, b"HP", Some(haplotype as i64 + 1))?;
                set_int_tag(rec, b"PC", None)?;
                set_int_tag(rec, b"PS", Some(ps))?;
                phaseset = ps.to_string();
                is_tagged = true;
                break;
            }
        }
    }
    Ok((is_tagged, haplotype_name, phaseset))
}

/// `ignore_read`
fn ignore_read(rec: &Record, include_supplementary: bool) -> bool {
    if rec.is_unmapped() || rec.is_secondary() {
        true
    } else if include_supplementary && rec.is_supplementary() {
        false
    } else {
        rec.is_supplementary()
    }
}

fn md5_of(path: &str) -> Result<String> {
    use md5::{Digest, Md5};
    let data = std::fs::read(path)?;
    Ok(Md5::digest(&data).iter().map(|b| format!("{b:02x}")).collect())
}

/// Writes `path`'s index beside it: `.bai`, or `.csi` when a contig is too long for a BAI.
fn write_bam_index(path: &std::path::Path) -> Result<()> {
    let context = || format!("indexing {}", path.display());
    match pasta::bam::fs::index(path).with_context(context)? {
        pasta::bam::Index::Bai(index) => {
            pasta::bam::bai::fs::write(path.with_added_extension("bai"), &index).with_context(context)
        }
        pasta::bam::Index::Csi(index) => {
            pasta::csi::fs::write(path.with_added_extension("csi"), &index).with_context(context)
        }
    }
}

/// `utils.stdout_is_regular_file`
fn stdout_is_regular_file() -> bool {
    unsafe {
        let mut st: libc::stat = std::mem::zeroed();
        libc::fstat(1, &mut st) == 0 && (st.st_mode & libc::S_IFMT) == libc::S_IFREG
    }
}

#[derive(Default)]
struct Counts {
    n_alignments: u64,
    n_tagged: u64,
    n_multiple_phase_sets: usize,
}

/// Computes the read-to-haplotype assignments of the chromosomes, in
/// parallel with `threads` > 1, and hands them to the writer in order.
struct Preparation<'a> {
    opts: &'a HaplotagOptions,
    tasks: &'a [(String, Regions)],
    samples: &'a [String],
    shared_samples: &'a [String],
    pools: Option<&'a Pools>,
    threads: usize,
}

impl Preparation<'_> {
    fn run(&self, mut write: impl FnMut(&str, &Regions, HaplotagInfo) -> Result<()>) -> Result<()> {
        let (opts, tasks, samples, shared_samples) = (self.opts, self.tasks, self.samples, self.shared_samples);
        if self.threads == 1 {
            let mut worker = PrepWorker::new(opts, None)?;
            for (chrom, regions) in tasks {
                log::debug!("Processing chromosome {chrom}");
                let t = Instant::now();
                let info = worker.prepare(samples, shared_samples, chrom, regions)?;
                log::debug!("Prepared {chrom} in {:.1} s", t.elapsed().as_secs_f64());
                if let Some(info) = info {
                    write(chrom, regions, info)?;
                }
            }
            return Ok(());
        }
        // Chromosomes prepared at the same time (and ahead of the one being
        // written). Writing is the bottleneck and preparing uses all threads
        // for allele detection, so two suffice; more only cost memory and
        // compete with the writer.
        let workers = self.threads.min(2);
        let (task_tx, task_rx) = sync_channel::<usize>(workers);
        let task_rx = Arc::new(Mutex::new(task_rx));
        let (res_tx, res_rx) = channel::<(usize, Result<Option<HaplotagInfo>>)>();
        let pools = self.pools;
        std::thread::scope(|scope| -> Result<()> {
            for _ in 0..workers {
                let task_rx = Arc::clone(&task_rx);
                let res_tx = res_tx.clone();
                scope.spawn(move || {
                    let mut worker = match PrepWorker::new(opts, pools) {
                        Ok(w) => w,
                        Err(e) => {
                            let _ = res_tx.send((usize::MAX, Err(e)));
                            return;
                        }
                    };
                    loop {
                        let job = task_rx.lock().unwrap().recv();
                        let Ok(i) = job else { break };
                        let (chrom, regions) = &tasks[i];
                        let t = Instant::now();
                        let r = worker.prepare(samples, shared_samples, chrom, regions);
                        log::debug!("Prepared {chrom} in {:.1} s", t.elapsed().as_secs_f64());
                        if res_tx.send((i, r)).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(res_tx);
            let mut pending: HashMap<usize, Result<Option<HaplotagInfo>>> = HashMap::new();
            let mut next_task = 0;
            // keep up to `workers` chromosomes in flight beyond the one being written
            while next_task < tasks.len().min(workers) {
                task_tx.send(next_task).ok();
                next_task += 1;
            }
            let mut result = Ok(());
            'outer: for i in 0..tasks.len() {
                let waiting = Instant::now();
                let info = loop {
                    if let Some(r) = pending.remove(&i) {
                        break r;
                    }
                    match res_rx.recv() {
                        Ok((usize::MAX, Err(e))) => {
                            result = Err(e);
                            break 'outer;
                        }
                        Ok((j, r)) => {
                            pending.insert(j, r);
                        }
                        Err(_) => {
                            result = Err(anyhow!("worker threads terminated unexpectedly"));
                            break 'outer;
                        }
                    }
                };
                if next_task < tasks.len() {
                    task_tx.send(next_task).ok();
                    next_task += 1;
                }
                let (chrom, regions) = &tasks[i];
                log::debug!(
                    "Processing chromosome {chrom} (waited {:.1} s for its haplotype assignments)",
                    waiting.elapsed().as_secs_f64()
                );
                match info {
                    Ok(Some(info)) => {
                        if let Err(e) = write(chrom, regions, info) {
                            result = Err(e);
                            break;
                        }
                    }
                    Ok(None) => {}
                    Err(e) => {
                        result = Err(e);
                        break;
                    }
                }
            }
            drop(task_tx);
            result
        })
    }
}

/// Thread pools used with --threads > 1
struct Pools {
    /// Threads decompressing each input file
    decompression_threads: usize,
    /// Allele detection
    detection: Arc<rayon::ThreadPool>,
}

/// Per-thread resources for computing read-to-haplotype assignments.
struct PrepWorker<'a> {
    tagger: Haplotagger<'a>,
    indexed: IndexedVcf,
    ploidy: Option<usize>,
    warned_duplicate: bool,
    sample_ids: HashMap<String, i32>,
}

impl<'a> PrepWorker<'a> {
    fn new(opts: &'a HaplotagOptions, pools: Option<&Pools>) -> Result<Self> {
        let mut readset_reader = ReadSetReader::new(
            std::slice::from_ref(&opts.alignment_file),
            opts.reference.as_deref(),
            20,
            opts.supplementary_strategy.consider_supplementary(),
            100_000,
            pools.map_or(1, |p| p.decompression_threads),
        )?;
        readset_reader.set_haplotag_options(true, opts.supplementary_strategy.attempt_to_haplotag_each_supplementary());
        if let Some(pools) = pools {
            readset_reader.set_detection_pool(Arc::clone(&pools.detection));
        }
        let fasta = match (&opts.reference, opts.no_reference) {
            (Some(path), false) => Some(ReferenceCache::open(path)?),
            _ => None,
        };
        Ok(PrepWorker {
            tagger: Haplotagger { opts, readset_reader, fasta },
            indexed: IndexedVcf::open(&opts.variant_file)?,
            ploidy: Some(opts.ploidy),
            warned_duplicate: false,
            sample_ids: HashMap::new(),
        })
    }

    /// Variant loading + `prepare_haplotag_information` for one chromosome;
    /// None if the chromosome is skipped (--skip-missing-contigs).
    fn prepare(
        &mut self,
        samples: &[String],
        shared_samples: &[String],
        chrom: &str,
        regions: &[(i64, Option<i64>)],
    ) -> Result<Option<HaplotagInfo>> {
        let opts = self.tagger.opts;
        let table = match crate::vcf::fetch_variant_table(
            &mut self.indexed,
            samples,
            false,
            true,
            &mut self.ploidy,
            &mut self.warned_duplicate,
            chrom,
            regions,
        ) {
            Ok(t) => t,
            Err(VcfFetchError::InvalidContig) => {
                if opts.untag_missing_contigs {
                    log::info!("Writing reads on '{chrom}' untagged because the contig does not exist in the VCF");
                    return Ok(Some(HaplotagInfo::untagged()));
                }
                if opts.skip_missing_contigs {
                    log::info!("Skipping reads on '{chrom}' because the contig does not exist in the VCF");
                    return Ok(None);
                }
                bail!(
                    "Input BAM/CRAM contains reads on contig '{chrom}', but that contig does not exist in the VCF header. To bypass this check, use --skip-missing-contigs"
                );
            }
            Err(VcfFetchError::IndexMissing) => bail!("{} is missing an index (.tbi or .csi)", opts.variant_file),
            Err(VcfFetchError::Other(e)) => return Err(e),
        };
        log::debug!("Loaded {} variants for chromosome {chrom} in VCF", table.len());
        Ok(Some(self.tagger.prepare(&table, shared_samples, &mut self.sample_ids, regions)?))
    }
}

/// The haplotag list (`--output-haplotag-list`), gzip-compressed if the
/// file name ends in .gz.
enum ListWriter {
    Plain(std::io::BufWriter<std::fs::File>),
    Gzip(Box<std::io::BufWriter<flate2::write::GzEncoder<std::fs::File>>>),
}

impl Write for ListWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            ListWriter::Plain(w) => w.write(buf),
            ListWriter::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            ListWriter::Plain(w) => w.flush(),
            ListWriter::Gzip(w) => w.flush(),
        }
    }
}

impl ListWriter {
    /// Writes the remaining data (and the gzip trailer), reporting errors.
    fn finish(self) -> std::io::Result<()> {
        match self {
            ListWriter::Plain(mut w) => w.flush(),
            ListWriter::Gzip(w) => (*w).into_inner().map_err(|e| e.into_error())?.finish()?.flush(),
        }
    }
}

/// The alignment-writing side of haplotag.
struct TagOutput<'a> {
    bam_reader: &'a mut AlignmentFile,
    /// None with `--list-only`.
    writer: Option<&'a mut alignment::AlignmentWriter>,
    list: &'a mut Option<ListWriter>,
    counts: &'a mut Counts,
}

/// Tag (or untag) one alignment of `chrom`; returns the haplotype name and
/// phase set for the haplotag list and whether it was tagged.
fn tag_alignment(
    rec: &mut Record,
    chrom: &str,
    info: &HaplotagInfo,
    opts: &HaplotagOptions,
) -> Result<(bool, String, String)> {
    let consider_supp = opts.supplementary_strategy.consider_supplementary();
    let (mut tagged, mut haplotype_name, mut phaseset) = (false, "none".to_string(), "none".to_string());
    if !ignore_read(rec, consider_supp) {
        (tagged, haplotype_name, phaseset) = attempt_add_phase_information(rec, chrom, info, opts)?;
    }
    if !tagged {
        set_int_tag(rec, b"HP", None)?;
        set_int_tag(rec, b"PC", None)?;
        set_int_tag(rec, b"PS", None)?;
    }
    Ok((tagged, haplotype_name, phaseset))
}

/// Whether an alignment appears in the haplotag list.
fn listed(rec: &Record, opts: &HaplotagOptions) -> bool {
    !(rec.is_secondary() || (rec.is_supplementary() && !opts.supplementary_strategy.consider_supplementary()))
}

impl TagOutput<'_> {
    fn write_chromosome(
        &mut self,
        chrom: &str,
        regions: &[(i64, Option<i64>)],
        info: &HaplotagInfo,
        opts: &HaplotagOptions,
    ) -> Result<()> {
        self.counts.n_multiple_phase_sets += info.n_multiple_phase_sets;
        let tid = self.bam_reader.header().tid(chrom).ok_or_else(|| anyhow!("unknown reference {chrom}"))?;
        // `--regions-bed`'s are sorted and apart: an alignment starting before the end of the
        // region before was written with it. WhatsHap writes it again for each of its
        // `--regions` it reaches, unsorted.
        let once = opts.regions_bed.is_some();
        let mut previous_end = None;
        for &(start, end) in regions {
            let fetch = match (start, end) {
                (0, None) => Fetch::Region(tid, 0, WHOLE_CONTIG_END),
                (s, e) => Fetch::Region(tid, s, e.unwrap_or(MAX_REGION_END)),
            };
            self.bam_reader.fetch(fetch)?;
            let written_before = previous_end.filter(|_| once);
            previous_end = end;
            while let Some(mut rec) = self.bam_reader.read()? {
                if written_before.is_some_and(|e| rec.pos() < e) {
                    continue;
                }
                self.counts.n_alignments += 1;
                let (tagged, haplotype_name, phaseset) = tag_alignment(&mut rec, chrom, info, opts)?;
                self.counts.n_tagged += tagged as u64;
                if let Some(w) = self.writer.as_mut() {
                    w.write(&rec)?;
                }
                if let Some(w) = self.list.as_mut()
                    && listed(&rec, opts)
                {
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}",
                        String::from_utf8_lossy(rec.qname()),
                        haplotype_name,
                        phaseset,
                        chrom
                    )?;
                }
            }
        }
        Ok(())
    }
}

/// End coordinate used for regions without an end
const MAX_REGION_END: i64 = (1 << 31) - 1;

/// End of a whole-contig fetch (`HTS_POS_MAX`)
const WHOLE_CONTIG_END: i64 = i64::MAX >> 1;

pub fn run_haplotag(opts: HaplotagOptions) -> Result<()> {
    let start_time = Instant::now();
    if let Some(reference) = opts.reference.as_deref().filter(|_| !opts.no_reference) {
        crate::fasta::check_reference(reference)?;
    }
    if opts.list_only {
        if opts.haplotag_list.is_none() {
            bail!("--list-only writes only the haplotag list: give --output-haplotag-list");
        }
        if opts.output.is_some() || opts.write_index {
            bail!("--list-only writes no alignments: leave out --output and --write-index");
        }
    }
    if opts.write_index && !opts.output.as_deref().is_some_and(|o| o.ends_with(".bam")) {
        bail!("--write-index indexes a BAM output file: give --output FILE.bam");
    }
    if opts.output.is_none() && !opts.list_only && unsafe { libc::isatty(1) } == 1 {
        bail!(
            "Refusing to write BAM to the terminal. Either use the '-o' option or redirect standard output with '>'."
        );
    }
    let vcf_reader = VcfReader::open_with(&opts.variant_file, false, true, Some(opts.ploidy))
        .with_context(|| format!("Error while loading variant file {}", opts.variant_file))?;
    let indexed_vcf = IndexedVcf::open(&opts.variant_file)?;

    // compute_variant_file_samples_to_use
    let samples_in_vcf: BTreeSet<String> = vcf_reader.samples.iter().cloned().collect();
    if samples_in_vcf.is_empty() {
        bail!("No samples detected in VCF file; cannot perform haplotagging");
    }
    log::info!("Found {} sample(s) in input VCF", samples_in_vcf.len());
    if opts.ignore_read_groups && opts.given_samples.is_none() && samples_in_vcf.len() > 1 {
        bail!(
            "When setting \"--ignore-read-groups\" on a multi-sample VCF, samples to be used must be specified via the \"--sample\" parameter."
        );
    }
    let use_vcf_samples: BTreeSet<String> = match &opts.given_samples {
        None => samples_in_vcf.clone(),
        Some(given) => {
            let missing: BTreeSet<&String> = given.iter().filter(|s| !samples_in_vcf.contains(*s)).collect();
            if !missing.is_empty() {
                bail!(
                    "The following samples were specified via the \"--sample\" parameter, but are not part of the input VCF: {:?}",
                    missing
                );
            }
            let keep: BTreeSet<String> = given.iter().cloned().collect();
            log::info!("Keeping {} sample(s) for haplo-tagging", keep.len());
            keep
        }
    };

    let threads = opts.threads.max(1);
    // Shared by all input files (must outlive them)
    let pools = if threads > 1 {
        Some(Pools {
            decompression_threads: threads,
            detection: Arc::new(rayon::ThreadPoolBuilder::new().num_threads(threads).build()?),
        })
    } else {
        None
    };
    let decompression_threads = if threads > 1 { threads } else { 1 };
    let reference = match (&opts.reference, opts.no_reference) {
        (Some(r), false) => Some(std::path::PathBuf::from(r)),
        _ => None,
    };
    let mut bam_reader = AlignmentFile::open(&opts.alignment_file, reference.as_deref(), decompression_threads)
        .with_context(|| format!("Error while loading alignment file {}", opts.alignment_file))?;
    if !bam_reader.load_index()? {
        bail!("Error while loading alignment file {}: the file is not indexed", opts.alignment_file);
    }
    let references: Vec<(String, u64)> = bam_reader.header().refs.iter().map(|(n, l)| (n.clone(), *l as u64)).collect();
    let header_text = bam_reader.header().text.clone();
    let mut header_dict = header_to_dict(&header_text, &references)?;

    // compute_shared_samples
    let mut bam_samples: BTreeSet<String> = BTreeSet::new();
    if let Some((_, rgs)) = header_dict.iter().find(|(r, _)| r == "RG") {
        for rg in rgs {
            bam_samples.insert(rg.iter().find(|(k, _)| k == "SM").map(|(_, v)| v.clone()).unwrap_or_default());
        }
    }
    log::info!("Found {} sample(s) in BAM file", bam_samples.len());
    let shared_samples: Vec<String> = if !opts.ignore_read_groups {
        let shared: Vec<String> = bam_samples.intersection(&use_vcf_samples).cloned().collect();
        if shared.is_empty() {
            bail!(
                "No common samples between VCF and BAM file detected. You may restart the analysis setting \"--ignore-read-groups\" (if appropriate) to avoid this error."
            );
        } else if shared.len() < bam_samples.len() {
            let missing: Vec<&String> = bam_samples.iter().filter(|s| !shared.contains(s)).collect();
            log::warn!(
                "Ignoring the following sample(s) for haplo-tagging because they are not part of the VCF or were not requested via \"--sample\": {}",
                missing.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(" | ")
            );
        }
        shared
    } else {
        use_vcf_samples.iter().cloned().collect()
    };

    // normalize_user_regions
    let mut user_regions: Vec<(String, Regions)> = Vec::new();
    match &opts.regions {
        None if opts.regions_bed.is_some() => {
            let path = opts.regions_bed.as_deref().unwrap_or_default();
            // In the header's order, as the output is.
            let bed = lrpipe_reads::regions::Regions::read(std::path::Path::new(path))?
                .in_order(references.iter().map(|(n, _)| n.as_str()));
            let refnames: HashSet<&str> = references.iter().map(|(n, _)| n.as_str()).collect();
            for (chrom, intervals) in bed.iter() {
                if !refnames.contains(chrom) {
                    bail!("Requested reference '{chrom}' of {path} not found in input BAM/CRAM");
                }
                user_regions.push((chrom.to_owned(), intervals.iter().map(|&(s, e)| (s, Some(e))).collect()));
            }
        }
        None => {
            for (name, _) in &references {
                user_regions.push((name.clone(), vec![(0, None)]));
            }
        }
        Some(specs) => {
            let refnames: HashSet<&String> = references.iter().map(|(n, _)| n).collect();
            for spec in specs {
                let (chrom, start, end) = parse_region(spec)?;
                if !refnames.contains(&chrom) {
                    bail!("Requested reference '{{region.chromosome}}' not found in input BAM/CRAM");
                }
                match user_regions.iter_mut().find(|(c, _)| *c == chrom) {
                    Some((_, v)) => v.push((start, end)),
                    None => user_regions.push((chrom, vec![(start, end)])),
                }
            }
        }
    }
    let include_unmapped = opts.regions.is_none() && opts.regions_bed.is_none();

    if opts.reference.is_none() && !opts.no_reference {
        bail!(
            "A reference FASTA needs to be provided with -r/--reference; or use --no-reference at the expense of phasing quality."
        );
    }

    // open_output_alignment_file (none with --list-only)
    let mut writer = if opts.list_only {
        None
    } else {
        let pg_id_prefix = "lrpipe-phase";
        let mut pg: Fields = vec![
            ("ID".into(), pg_id_prefix.into()),
            ("PN".into(), "lrpipe-phase".into()),
            ("VN".into(), env!("CARGO_PKG_VERSION").into()),
            ("CL".into(), opts.command_line.clone()),
            ("m5".into(), md5_of(&opts.variant_file)?),
        ];
        match header_dict.iter_mut().find(|(r, _)| r == "PG") {
            Some((_, pgs)) => {
                let ids: HashSet<String> = pgs
                    .iter()
                    .filter_map(|f| f.iter().find(|(k, _)| k == "ID").map(|(_, v)| v.clone()))
                    .filter(|id| id.starts_with(pg_id_prefix))
                    .collect();
                let mut nr = 1;
                while ids.contains(&pg[0].1) {
                    pg[0].1 = format!("{pg_id_prefix}.{nr}");
                    nr += 1;
                }
                pgs.push(pg);
            }
            None => header_dict.push(("PG".into(), vec![pg])),
        }
        let (text, refs) = header_from_dict(&header_dict);
        let (out_path, mode, out_ref) = match &opts.output {
            // Like WhatsHap, write uncompressed BAM to a pipe unless a level is requested
            None => match opts.compression_level {
                Some(level) => ("-".to_string(), format!("wb{level}"), None),
                None => ("-".to_string(), if stdout_is_regular_file() { "wb" } else { "wb0" }.to_string(), None),
            },
            Some(p) if p.ends_with(".cram") => {
                let Some(r) = &opts.reference else {
                    bail!("Writing CRAM output requires FASTA reference file given via \"--reference\"")
                };
                let level = opts.compression_level.map(|l| l.to_string()).unwrap_or_default();
                (p.clone(), format!("wc{level}"), Some(r.as_str()))
            }
            Some(p) => {
                let level = opts.compression_level.map(|l| l.to_string()).unwrap_or_default();
                (p.clone(), format!("wb{level}"), None)
            }
        };
        let output_threads = opts.output_threads.unwrap_or(opts.threads);
        Some(create_output(&out_path, &mode, &text, &refs, out_ref, output_threads)?)
    };
    let mut haplotag_writer: Option<ListWriter> = match &opts.haplotag_list {
        None => None,
        Some(path) => {
            let f = std::fs::File::create(path)
                .with_context(|| format!("Error while initializing haplotag list output at path: {path}"))?;
            let mut w = if path.ends_with(".gz") {
                ListWriter::Gzip(Box::new(std::io::BufWriter::new(flate2::write::GzEncoder::new(
                    f,
                    flate2::Compression::default(),
                ))))
            } else {
                ListWriter::Plain(std::io::BufWriter::new(f))
            };
            writeln!(w, "#readname\thaplotype\tphaseset\tchromosome")?;
            Some(w)
        }
    };

    // Without alignments to write, the CRAM records being listed need no more than reading
    // alleles does (no `MM`/`ML`, most of an ONT CRAM, nor MD and NM), decoded on the threads:
    // HG002 ONT chr20 took 44 s with 8 threads, writing the list as it decoded every field.
    if opts.list_only {
        bam_reader.set_decoding(alignment::Decoding::Alleles);
    }
    let mut counts = Counts::default();
    // contigs_with_alignments
    let mut has_alignments: HashSet<String> = HashSet::new();
    for (tid, (name, _)) in references.iter().enumerate() {
        if bam_reader.has_records(tid)? {
            has_alignments.insert(name.clone());
        }
    }
    let samples = vcf_reader.samples.clone();
    drop(indexed_vcf);
    let tasks: Vec<(String, Regions)> =
        user_regions.iter().filter(|(c, _)| has_alignments.contains(c)).cloned().collect();
    let preparation = Preparation {
        opts: &opts,
        tasks: &tasks,
        samples: &samples,
        shared_samples: &shared_samples,
        pools: pools.as_ref(),
        threads,
    };
    let mut out = TagOutput {
        bam_reader: &mut bam_reader,
        writer: writer.as_mut(),
        list: &mut haplotag_writer,
        counts: &mut counts,
    };
    preparation.run(|chrom, regions, info| {
        let t = Instant::now();
        out.write_chromosome(chrom, regions, &info, &opts)?;
        log::debug!("Wrote {chrom} in {:.1} s", t.elapsed().as_secs_f64());
        Ok(())
    })?;
    if let Some(writer) = writer.as_mut().filter(|_| include_unmapped) {
        log::debug!("Copying unmapped reads to output");
        let t = Instant::now();
        bam_reader.fetch(Fetch::Unmapped)?;
        while let Some(rec) = bam_reader.read()? {
            writer.write(&rec)?;
        }
        log::debug!("Copied unmapped reads in {:.1} s", t.elapsed().as_secs_f64());
    }
    if let Some(writer) = writer {
        writer.finish()?;
    }
    if opts.write_index
        && let Some(path) = &opts.output
    {
        let t = Instant::now();
        write_bam_index(std::path::Path::new(path))?;
        log::debug!("Indexed {path} in {:.1} s", t.elapsed().as_secs_f64());
    }
    if let Some(w) = haplotag_writer {
        w.finish().context("error writing the haplotag list")?;
    }
    log::info!("\n== SUMMARY ==");
    log::info!("Total alignments processed:              {:12}", counts.n_alignments);
    log::info!("Alignments that could be tagged:         {:12}", counts.n_tagged);
    log::info!("Alignments spanning multiple phase sets: {:12}", counts.n_multiple_phase_sets);
    log::info!("Finished in {:.1} s", start_time.elapsed().as_secs_f64());
    Ok(())
}

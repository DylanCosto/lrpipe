//! Port of `cli/split.py`: split reads (FASTQ or BAM) by haplotype, using a
//! haplotype assignment list such as `haplotag --output-haplotag-list`.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use anyhow::{Context, Result, bail};

use crate::alignment::bam::{BamReader, BamWriter, sanitise_header_text};
use crate::alignment::record::Record;
use crate::compat::pystr::py_strip;
use crate::fastx::FastxFile;
use crate::phase::detect_file_format;
use crate::xopen::{self, TextLines};

pub struct SplitOptions {
    pub reads_file: String,
    pub list_file: String,
    pub output_h1: Option<String>,
    pub output_h2: Option<String>,
    pub outputs: Option<Vec<String>>,
    pub output_untagged: Option<String>,
    pub add_untagged: bool,
    pub only_largest_block: bool,
    pub discard_unknown_reads: bool,
    pub read_lengths_histogram: Option<String>,
}

/// A parsed line of the haplotag list: read name, haplotype, phase set and
/// chromosome (the last two only with four columns)
type ListLine = (String, String, Option<String>, Option<String>);

/// The lines of the haplotag list
fn list_lines(path: &str) -> Result<TextLines<Box<dyn std::io::BufRead>>> {
    Ok(TextLines::new(xopen::open(path)?, path))
}

/// `check_haplotag_list_information`: whether the list has four columns
/// (phase set and chromosome)
fn has_four_columns(path: &str) -> Result<bool> {
    let first = list_lines(path)?.next().transpose()?.unwrap_or_default();
    let first_line = py_strip(&first);
    let n = first_line.split('\t').count();
    if n < 2 {
        bail!(
            "First line of haplotag list file does not have at least 2 columns, or it is not tab-separated: {first_line}"
        );
    }
    Ok(n >= 4)
}

fn parse_line(line: &str, four_columns: bool) -> Result<ListLine> {
    let cols: Vec<&str> = py_strip(line).split('\t').collect();
    if four_columns {
        if cols.len() < 4 {
            bail!("ValueError: not enough values to unpack (expected 4, got {})", cols.len());
        }
        Ok((cols[0].to_string(), cols[1].to_string(), Some(cols[2].to_string()), Some(cols[3].to_string())))
    } else {
        if cols.len() < 2 {
            bail!("IndexError: list index out of range");
        }
        Ok((cols[0].to_string(), cols[1].to_string(), None, None))
    }
}

/// Read counts per phase set of one chromosome, in order of first occurrence
/// (a `Counter`)
#[derive(Default)]
struct BlockCounts {
    order: Vec<String>,
    counts: HashMap<String, usize>,
}

/// `process_haplotag_list_file`: read name -> haplotype, and the read names
/// in the list (only with `discard_unknown_reads`)
fn process_haplotag_list(
    path: &str,
    four_columns: bool,
    only_largest_blocks: bool,
    discard_unknown_reads: bool,
    ploidy: usize,
) -> Result<(HashMap<String, usize>, HashSet<String>)> {
    let mut lines = list_lines(path)?.peekable();
    // a header line starting with '#' is skipped
    if let Some(Ok(first)) = lines.peek()
        && first.starts_with('#')
    {
        lines.next();
    }
    let mut block_sizes: HashMap<String, BlockCounts> = HashMap::new();
    let mut blocks_to_readnames: HashMap<(String, String), HashSet<String>> = HashMap::new();
    let mut known_reads = HashSet::new();
    let mut readname_to_haplotype: HashMap<String, usize> = HashMap::new();
    let mut total_reads = 0;
    for line in lines {
        let (readname, haplo_name, phaseset, chromosome) = parse_line(&line?, four_columns)?;
        total_reads += 1;
        let haplo_num = if haplo_name == "none" {
            0
        } else {
            match haplo_name.strip_prefix('H').and_then(|n| n.parse::<usize>().ok()) {
                Some(n) if (1..=ploidy).contains(&n) && haplo_name == format!("H{n}") => n,
                _ => {
                    log::error!(
                        "Haplotype name '{{haplo_name}}' in haplotype list file not recognized; must be one of 'none', 'H1', ..., 'H{ploidy}'"
                    );
                    bail!("KeyError: '{haplo_name}'");
                }
            }
        };
        if haplo_num == 0 {
            if discard_unknown_reads {
                known_reads.insert(readname);
            }
            continue;
        }
        readname_to_haplotype.insert(readname.clone(), haplo_num);
        if only_largest_blocks {
            let (chromosome, phaseset) = (chromosome.unwrap_or_default(), phaseset.unwrap_or_default());
            let counts = block_sizes.entry(chromosome.clone()).or_default();
            if !counts.counts.contains_key(&phaseset) {
                counts.order.push(phaseset.clone());
            }
            *counts.counts.entry(phaseset.clone()).or_default() += 1;
            blocks_to_readnames.entry((chromosome, phaseset)).or_default().insert(readname);
        }
    }
    let tagged_reads = readname_to_haplotype.len() as i64;
    log::info!("Total number of reads in haplotag list: {total_reads}");
    log::info!("Total number of haplo-tagged reads: {tagged_reads}");
    log::info!("Total number of untagged reads: {}", total_reads - tagged_reads);
    if discard_unknown_reads {
        known_reads.extend(readname_to_haplotype.keys().cloned());
        if total_reads != known_reads.len() as i64 {
            bail!(
                "AssertionError: Total number of reads is not equal to number of known reads: {total_reads} != {}",
                known_reads.len()
            );
        }
    }
    if only_largest_blocks {
        // select_reads_in_largest_phased_blocks
        let mut selected: HashSet<String> = HashSet::new();
        log::info!("Determining largest blocks/phasesets per chromosome");
        for (chromosome, counts) in &block_sizes {
            // Counter.most_common(1): the first of the largest
            let mut best: Option<(&String, usize)> = None;
            for name in &counts.order {
                let n = counts.counts[name];
                if best.is_none_or(|(_, m)| n > m) {
                    best = Some((name, n));
                }
            }
            let (block_name, reads_in_block) = best.expect("a chromosome has at least one block");
            log::info!("Chromosome: {chromosome} - Phaseset: {block_name} - Tagged reads: {reads_in_block}");
            selected.extend(blocks_to_readnames[&(chromosome.clone(), block_name.clone())].iter().cloned());
        }
        log::info!("Total number of haplo-tagged reads in all largest phased blocks: {}", selected.len());
        readname_to_haplotype =
            selected.into_iter().map(|k| (k.clone(), readname_to_haplotype.get(&k).copied().unwrap_or(0))).collect();
        log::info!(
            "Number of reads removed / reads not overlapping largest phased blocks: {}",
            total_reads - readname_to_haplotype.len() as i64
        );
    }
    Ok((readname_to_haplotype, known_reads))
}

#[derive(Clone, Copy, PartialEq)]
enum InputFormat {
    Bam,
    Fastq,
}

/// `utils.detect_file_format` and the FASTQ file name heuristic of
/// `initialize_io_files`
fn input_format(reads_file: &str) -> Result<InputFormat> {
    let detected = detect_file_format(reads_file)?;
    match detected {
        Some("BAM") => Ok(InputFormat::Bam),
        Some(f) => bail!("Input file format detected as: {f}. Currently, only BAM and FASTQ is supported."),
        None => {
            // the list in WhatsHap lacks two commas
            let extensions = ["fastq", "fastq.gz", "fastq.gzipfq", "fq.gzfq.gzip"];
            if extensions.iter().any(|ext| reads_file.ends_with(ext)) {
                Ok(InputFormat::Fastq)
            } else {
                bail!("Undetected file format for input reads. Expecting BAM or FASTQ (gzipped): {reads_file}")
            }
        }
    }
}

/// An output file for one haplotype (or untagged reads); None for /dev/null
enum Output {
    Bam(Option<BamWriter>),
    Fastq(Option<xopen::Writer>),
}

/// `write_read_length_histogram`
fn write_read_length_histogram(length_counts: &[HashMap<usize, u64>], path: &str) -> Result<()> {
    let mut all_read_lengths: Vec<usize> = length_counts.iter().flat_map(|lc| lc.keys().copied()).collect();
    all_read_lengths.sort_unstable();
    let mut f = xopen::create(path, None)?;
    let columns: Vec<String> = (1..length_counts.len()).map(|i| format!("count-h{i}")).collect();
    write!(f, "#length\tcount-untagged")?;
    for c in &columns {
        write!(f, "\t{c}")?;
    }
    writeln!(f)?;
    for rlen in all_read_lengths {
        write!(f, "{rlen}")?;
        for lc in length_counts {
            write!(f, "\t{}", lc.get(&rlen).copied().unwrap_or(0))?;
        }
        writeln!(f)?;
    }
    f.finish()
}

/// `_bam_iterator`: the read length, from the sequence or else the CIGAR
fn bam_read_length(rec: &Record) -> usize {
    if rec.l_seq > 0 {
        return rec.l_seq as usize;
    }
    // calculateQueryLengthWithoutHardClipping: M, I, S, = and X
    let inferred: u32 = rec
        .cigar()
        .iter()
        .filter(|(op, _)| matches!(op, 0 | 1 | 4 | 7 | 8))
        .fold(0u32, |acc, (_, n)| acc.wrapping_add(n));
    let inferred = inferred as i32;
    if inferred > 0 { inferred as usize } else { 0 }
}

/// Per-haplotype counters of processed reads
#[derive(Default)]
struct ReadCounter {
    total_reads: u64,
    unknown_reads: u64,
    skipped_reads: u64,
    by_haplotype: HashMap<usize, u64>,
}

pub fn run_split(opts: &SplitOptions) -> Result<()> {
    let (outputs, ploidy): (Vec<Option<String>>, usize) = if opts.output_h1.as_deref().is_some_and(|s| !s.is_empty())
        || opts.output_h2.as_deref().is_some_and(|s| !s.is_empty())
    {
        if opts.outputs.as_ref().is_some_and(|o| !o.is_empty()) {
            bail!("Cannot use output_h1/output_h2 and outputs at the same time");
        }
        (vec![opts.output_untagged.clone(), opts.output_h1.clone(), opts.output_h2.clone()], 2)
    } else {
        let Some(outputs) = &opts.outputs else { bail!("TypeError: object of type 'NoneType' has no len()") };
        let mut all = vec![opts.output_untagged.clone()];
        all.extend(outputs.iter().cloned().map(Some));
        (all, outputs.len())
    };
    let t_run = Instant::now();
    let t_init = Instant::now();

    let has_chrom_info = has_four_columns(&opts.list_file)?;
    if opts.only_largest_block {
        log::debug!(
            "User selected \"--only-largest-block\", this requires chromosome and phaseset information to be present in the haplotag list file."
        );
        if !has_chrom_info {
            bail!(
                "The haplotag list file does not contain phaseset and chromosome information, which is required to select only reads from the largest phased block. Columns 3 and 4 are missing."
            );
        }
    }
    let t_list = Instant::now();
    let (readname_to_haplotype, known_reads) = process_haplotag_list(
        &opts.list_file,
        has_chrom_info,
        opts.only_largest_block,
        opts.discard_unknown_reads,
        ploidy,
    )?;
    let mut missing_reads: i64 = if opts.discard_unknown_reads {
        log::debug!(
            "User selected to discard unknown reads, i.e., ignore all reads that are not part of the haplotag list input file."
        );
        if known_reads.is_empty() {
            bail!("AssertionError: No known reads in input set - would discard everything, this is probably wrong");
        }
        known_reads.len() as i64
    } else {
        -1
    };
    let list_time = t_list.elapsed();

    // initialize_io_files
    let format = input_format(&opts.reads_file)?;
    let mut bam_reader = None;
    let mut fastx_reader = None;
    let mut writers: Vec<Output> = Vec::new();
    match format {
        InputFormat::Bam => {
            let reader = BamReader::open(Path::new(&opts.reads_file), 1)?;
            let text = sanitise_header_text(&reader.raw_text)?;
            for path in &outputs {
                writers.push(Output::Bam(match path {
                    None => None,
                    Some(p) => {
                        // pysam's AlignmentFile("-", "wb") writes to standard output
                        let out: Box<dyn Write + Send> = if p == "-" {
                            Box::new(std::io::stdout())
                        } else {
                            Box::new(std::fs::File::create(p).with_context(|| format!("cannot create {p}"))?)
                        };
                        Some(BamWriter::create_with_text(out, &text, &reader.header.refs, None, 1)?)
                    }
                }));
            }
            bam_reader = Some(reader);
        }
        InputFormat::Fastq => {
            fastx_reader = Some(FastxFile::open(&opts.reads_file)?);
            for path in &outputs {
                writers.push(Output::Fastq(match path {
                    None => None,
                    Some(p) => Some(xopen::create(p, Some(1))?),
                }));
            }
        }
    }
    let init_time = t_init.elapsed();

    let mut histogram_data: Vec<HashMap<usize, u64>> = vec![HashMap::new(); outputs.len()];
    let mut counter = ReadCounter::default();
    let mut process_haplotype: Vec<bool> = outputs.iter().map(|o| o.is_some()).collect();
    process_haplotype[0] = process_haplotype[0] || opts.add_untagged;

    let t_iter = Instant::now();
    loop {
        // the next read: name, length and the record to write
        let (read_name, read_length, bam_record, fastq_text) = if let Some(reader) = &mut bam_reader {
            let Some(rec) = reader.read()? else { break };
            let name = String::from_utf8(rec.qname().to_vec()).context("UnicodeDecodeError in a read name")?;
            let len = bam_read_length(&rec);
            (name, len, Some(rec), None)
        } else {
            let reader = fastx_reader.as_mut().expect("FASTQ input");
            let Some(rec) = reader.next_record()? else { break };
            let len = rec.sequence.chars().count();
            let text = format!("{rec}\n");
            (rec.name, len, None, Some(text))
        };
        counter.total_reads += 1;
        if opts.discard_unknown_reads && !known_reads.contains(&read_name) {
            counter.unknown_reads += 1;
            continue;
        }
        let read_haplotype = readname_to_haplotype.get(&read_name).copied().unwrap_or(0);
        if !process_haplotype[read_haplotype] {
            counter.skipped_reads += 1;
            continue;
        }
        *histogram_data[read_haplotype].entry(read_length).or_default() += 1;
        *counter.by_haplotype.entry(read_haplotype).or_default() += 1;
        let targets: Vec<usize> =
            if read_haplotype == 0 && opts.add_untagged { (0..writers.len()).collect() } else { vec![read_haplotype] };
        for i in targets {
            match &mut writers[i] {
                Output::Bam(Some(w)) => w.write(bam_record.as_ref().expect("BAM input"))?,
                Output::Fastq(Some(w)) => w.write_all(fastq_text.as_ref().expect("FASTQ input").as_bytes())?,
                _ => {}
            }
        }
        if opts.discard_unknown_reads {
            missing_reads -= 1;
            if missing_reads == 0 {
                log::info!("All known reads processed - cancel processing...");
                break;
            }
        }
    }
    let iter_time = t_iter.elapsed();

    let mut histogram_time = None;
    if let Some(path) = &opts.read_lengths_histogram {
        let t = Instant::now();
        write_read_length_histogram(&histogram_data, path)?;
        histogram_time = Some(t.elapsed());
    }
    for w in writers {
        match w {
            Output::Bam(Some(w)) => w.finish()?,
            Output::Fastq(Some(w)) => w.finish()?,
            _ => {}
        }
    }

    let by = |h: usize| counter.by_haplotype.get(&h).copied().unwrap_or(0);
    log::info!("\n== SUMMARY ==");
    log::info!("Total reads processed: {}", counter.total_reads);
    log::info!("Number of output reads \"untagged\": {}", by(0));
    for h in 1..=ploidy {
        log::info!("Number of output reads haplotype {h}: {}", by(h));
    }
    log::info!("Number of unknown (dropped) reads: {}", counter.unknown_reads);
    log::info!("Number of skipped reads (per user request): {}", counter.skipped_reads);
    log::info!("Time for processing haplotag list: {:.3} sec", list_time.as_secs_f64());
    log::info!("Time for total initial setup: {:.3} sec", init_time.as_secs_f64());
    log::info!("Time for iterating input reads: {:.3} sec", iter_time.as_secs_f64());
    if let Some(t) = histogram_time {
        log::info!("Time for creating histogram output: {:.3} sec", t.as_secs_f64());
    }
    log::info!("Total run time: {:.3} sec", t_run.elapsed().as_secs_f64());
    Ok(())
}

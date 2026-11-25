//! Checks of a run's inputs before it is planned, each saying what is wrong in the user's terms
//! and what to do. Without them the mistake reaches a stage, which fails with a message only its
//! log holds ("invalid BAI header", "cannot read .../pileup"), or does not fail at all: a
//! reference naming its contigs `20` against reads on `chr20` gave an empty VCF and exit code 0.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use md5::{Digest, Md5};
use pasta::csi::binning_index::ReferenceSequence as _;
use pasta::csi::binning_index::index::ReferenceSequence;
use pasta::{bam, cram, csi, sam};

use crate::plan::{Layout, Stage};
use lrpipe_reads::regions::Regions;
use lrpipe_reads::{IndexedReads, Reference};

/// Refuses a sample name that is not letters, digits, `.`, `-` and `_`, or that starts with `.`
/// or `-`. The name makes the outputs' file names (`sv/NAME.sv.vcf`) and fills the VCFs' sample
/// column: a `/` or `..` put files outside the output directory, an empty name or one starting
/// with `.` made hidden files, a tab broke the merged VCF, a space made the sample column two
/// words, and one starting with `-` reads as an option to the tools. It is at most
/// [`SAMPLE_NAME_MAX`] characters, which with the longest suffix fits a file name.
///
/// # Errors
///
/// If the name is not one of those.
pub fn sample_name(name: &str) -> Result<()> {
    if name.len() > SAMPLE_NAME_MAX {
        bail!("--sample is {} characters long: a sample name is at most {SAMPLE_NAME_MAX}", name.len());
    }
    let allowed = |c: char| c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_');
    let wrong = if name.is_empty() {
        Some("it is empty".to_owned())
    } else if let Some(c) = name.chars().find(|&c| !allowed(c)) {
        Some(format!("it has {c:?}"))
    } else if name.starts_with(['.', '-']) {
        Some(format!("it starts with {:?}", &name[..1]))
    } else {
        None
    };
    if let Some(wrong) = wrong {
        bail!(
            "--sample {name:?}: {wrong}; a sample name is letters, digits, '.', '-' and '_', not starting with '.' or \
             '-' (it names the outputs' files and fills the VCFs' sample column)"
        );
    }
    Ok(())
}

/// The longest sample name: file names are at most 255 bytes, and the outputs add up to 22
/// (`.cpg.combined.bed.gz.tbi`).
pub const SAMPLE_NAME_MAX: usize = 200;

/// Checks that the reference exists and that its indexes fit the file.
///
/// # Errors
///
/// If it or an index is missing, unreadable, or does not match the reference.
pub fn reference(reference: &Path) -> Result<()> {
    let shown = reference.display();
    match std::fs::metadata(reference) {
        Ok(m) if m.is_dir() => bail!("the reference {shown} (the site configuration's `reference`) is a directory"),
        Ok(_) => {}
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            bail!("the reference {shown} (the site configuration's `reference`) does not exist")
        }
        Err(e) => return Err(e).with_context(|| format!("reading the reference {shown}")),
    }
    let faidx = format!("make the indexes with `samtools faidx {shown}`");
    let fai = Layout::index(reference, "fai");
    if !fai.exists() {
        bail!("the reference {shown} has no index {}: {faidx}", fai.display());
    }
    let mut magic = [0; 2];
    let n = File::open(reference).and_then(|mut f| read_up_to(&mut f, &mut magic)).context("reading the reference")?;
    let gzi = Layout::index(reference, "gzi");
    if n == 2 && magic == [0x1f, 0x8b] && !gzi.exists() {
        bail!("the reference {shown} is compressed and has no index {}: {faidx}", gzi.display());
    }
    lrpipe_reads::check_fai(reference).with_context(|| format!("checking the reference index for {shown}"))
}

/// Reads until `buf` is full or the reader ends; the bytes read.
fn read_up_to(reader: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match reader.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// An aligned file of reads as [`aligned`] found it.
#[derive(Clone, Debug)]
pub struct Aligned {
    /// The file.
    pub path: PathBuf,
    /// Its header.
    pub header: sam::Header,
    /// Mapped reads on each of the header's contigs, by a BAM index's counts; none for CRAM, or
    /// an index without them.
    pub mapped: Option<Vec<u64>>,
    /// What the run should be warned of.
    pub warnings: Vec<String>,
}

/// Checks an aligned BAM or CRAM, whose format and name [`crate::input`] has checked, with its
/// `header`, before anything reads its index: that it is not sorted by read name, and that its
/// index is beside it under the name the tools look for (`READS.bai` or `READS.csi`,
/// `READS.crai`) and indexes a file with its header's contigs. A file whose header says
/// `SO:unsorted` but that has such an index is taken as sorted, with a warning: many tools
/// write `unsorted` for what they did not sort themselves.
///
/// # Errors
///
/// If any of that is not so, saying how to make the file or index the tools need.
pub fn aligned(reads: &Path, header: sam::Header) -> Result<Aligned> {
    let shown = reads.display();
    let order =
        header.header().and_then(|hd| hd.other_fields().get(&sam::header::record::value::map::header::tag::SORT_ORDER));
    let sort = format!("`samtools sort -o sorted.bam {shown} && samtools index sorted.bam`");
    let mut warnings = Vec::new();
    if let Some(index) = lrpipe_reads::stale_index(reads) {
        warnings.push(format!(
            "the index {} is older than --reads {shown}: if the reads were written again since, index them again \
             (`samtools index`), as regions read through an old index can come out empty",
            index.display()
        ));
    }
    match order.map(|o| &o[..]) {
        Some(b"queryname") => bail!(
            "--reads {shown} is sorted by read name (@HD SO:queryname), not by position as its index needs: sort and \
             index it, as with {sort}"
        ),
        Some(b"unsorted") => warnings.push(format!(
            "--reads {shown} says it is unsorted (@HD SO:unsorted) but has an index, so it is taken as sorted by \
             position; if it is not, sort and index it, as with {sort}"
        )),
        _ => {}
    }
    let mapped = index(reads, lrpipe_reads::is_cram(reads), header.reference_sequences().len())?;
    Ok(Aligned { path: reads.to_path_buf(), header, mapped, warnings })
}

/// Checks the reads' index: there under the name the tools look for, readable, and of as many
/// contigs as the header has (one made for another file is not, as a rule). Returns a BAM
/// index's mapped reads per contig, when it counts them.
fn index(reads: &Path, cram: bool, contigs: usize) -> Result<Option<Vec<u64>>> {
    let shown = reads.display();
    let samtools = format!("index it with `samtools index {shown}`");
    let names: &[&str] = if cram { &["crai"] } else { &["bai", "csi"] };
    let Some((ext, path)) = names.iter().map(|&ext| (ext, Layout::index(reads, ext))).find(|(_, p)| p.exists()) else {
        let elsewhere = names
            .iter()
            .map(|ext| reads.with_extension(ext))
            .find(|p| p.exists())
            .map(|p| {
                format!(
                    " ({} is there, but the tools look for the index under the reads' name and the index's extension)",
                    p.display()
                )
            })
            .unwrap_or_default();
        let wanted: Vec<String> = names.iter().map(|ext| Layout::index(reads, ext).display().to_string()).collect();
        bail!("--reads {shown} has no index {}{elsewhere}: {samtools}", wanted.join(" or "));
    };
    let not_its =
        |what: String| anyhow::anyhow!("{} is not the index of --reads {shown} ({what}): {samtools}", path.display());
    let count = |n: usize| {
        (n != contigs)
            .then(|| format!("it indexes {n} contig{}, the reads' header has {contigs}", if n == 1 { "" } else { "s" }))
    };
    let (problem, mapped) = match ext {
        "bai" => {
            bam::bai::fs::read(&path).map(|i| (count(i.reference_sequences().len()), counted(i.reference_sequences())))
        }
        "csi" => csi::fs::read(&path).map(|i| (count(i.reference_sequences().len()), counted(i.reference_sequences()))),
        _ => cram::crai::fs::read(&path).and_then(|i| {
            let last = i.iter().filter_map(|r| r.reference_sequence_id()).max();
            let problem = match last.filter(|&id| id >= contigs) {
                Some(id) => Some(format!("it indexes contig {}, the reads' header has {contigs}", id + 1)),
                None => misplaced_container(reads, &i)?.map(|offset| {
                    format!(
                        "it puts a container at byte {offset}, where none starts, as an index of another version of \
                         the file does"
                    )
                }),
            };
            Ok((problem, None))
        }),
    }
    .map_err(|e| not_its(format!("reading it: {e}")))?;
    problem.map_or(Ok(mapped), |what| Err(not_its(what)))
}

/// The first of a CRAM index's first, middle and last entries whose container offset is not
/// where a container of `reads` starts, by the container header's CRC32 (CRAM 3). Without this,
/// an index of another version of the file (ONT's published HG002 and HG004 indexes are 63 and
/// 68 bytes off theirs) reads as a corrupt file: "container header checksum mismatch".
fn misplaced_container(reads: &Path, index: &cram::crai::Index) -> io::Result<Option<u64>> {
    use std::io::{Seek, SeekFrom};
    let mut file = File::open(reads)?;
    let mut magic = [0; 5];
    file.read_exact(&mut magic)?;
    if &magic[..4] != b"CRAM" || magic[4] < 3 {
        return Ok(None);
    }
    let mut offsets: Vec<u64> = [0, index.len() / 2, index.len().saturating_sub(1)]
        .iter()
        .filter_map(|&i| index.get(i).map(cram::crai::Record::offset))
        .collect();
    offsets.dedup();
    for offset in offsets {
        let mut header = Vec::with_capacity(CONTAINER_HEADER_READ);
        file.seek(SeekFrom::Start(offset))?;
        (&mut file).take(CONTAINER_HEADER_READ as u64).read_to_end(&mut header)?;
        if container_header_crc(&header) == Some(false) {
            return Ok(Some(offset));
        }
    }
    Ok(None)
}

/// Bytes read for a container header: one with as many slices as this allows is not checked.
const CONTAINER_HEADER_READ: usize = 4096;

/// Whether `bytes` start with a CRAM 3 container header whose CRC32 holds; `None` when they end
/// before the header does.
fn container_header_crc(bytes: &[u8]) -> Option<bool> {
    fn itf8(b: &[u8], p: &mut usize) -> Option<()> {
        let n = b.get(*p)?.leading_ones().min(4) as usize;
        *p += 1 + n;
        (*p <= b.len()).then_some(())
    }
    fn ltf8(b: &[u8], p: &mut usize) -> Option<()> {
        let n = b.get(*p)?.leading_ones() as usize;
        *p += 1 + n;
        (*p <= b.len()).then_some(())
    }
    let mut p = 4;
    for _ in 0..4 {
        itf8(bytes, &mut p)?;
    }
    ltf8(bytes, &mut p)?;
    ltf8(bytes, &mut p)?;
    itf8(bytes, &mut p)?;
    let landmarks_at = p;
    itf8(bytes, &mut p)?;
    let landmarks = itf8_value(&bytes[landmarks_at..p]);
    for _ in 0..landmarks {
        itf8(bytes, &mut p)?;
    }
    let stored = u32::from_le_bytes(bytes.get(p..p + 4)?.try_into().ok()?);
    let mut crc = flate2::Crc::new();
    crc.update(&bytes[..p]);
    Some(crc.sum() == stored)
}

/// The value of one ITF8 integer, `bytes` long (1 to 5).
fn itf8_value(bytes: &[u8]) -> u32 {
    let b = |i: usize| u32::from(bytes[i]);
    match bytes.len() {
        1 => b(0),
        2 => ((b(0) & 0x3F) << 8) | b(1),
        3 => ((b(0) & 0x1F) << 16) | (b(1) << 8) | b(2),
        4 => ((b(0) & 0x0F) << 24) | (b(1) << 16) | (b(2) << 8) | b(3),
        _ => ((b(0) & 0x0F) << 28) | (b(1) << 20) | (b(2) << 12) | (b(3) << 4) | (b(4) & 0x0F),
    }
}

/// Each contig's mapped reads, if the index counts them for every one: a contig without bins
/// has none, as `samtools idxstats` counts it (an index has no metadata for it).
fn counted<I>(contigs: &[ReferenceSequence<I>]) -> Option<Vec<u64>>
where
    I: pasta::csi::binning_index::index::reference_sequence::Index,
    ReferenceSequence<I>: pasta::csi::binning_index::ReferenceSequence,
{
    contigs
        .iter()
        .map(|r| r.metadata().map(|m| m.mapped_record_count()).or_else(|| r.bins().is_empty().then_some(0)))
        .collect()
}

/// The share of a BAM's mapped reads that may be on contigs the reference lacks, with a warning:
/// HLA, decoy and EBV contigs of a whole-genome reference that the no-alt analysis set leaves
/// out hold some reads, which then get no calls.
pub const MISSING_SHARE: f64 = 0.01;

/// Checks that the reference matches the contigs with reads by length and, when present in the
/// header, MD5. A contig it lacks gets no calls (lrpipe-snv wrote an empty VCF and exited 0 for a
/// reference naming its contigs `20` against reads on `chr20`), and one of another length is
/// another assembly's, an error. Contigs without reads are not checked, as lrpipe-methyl checks
/// (LPH-004), so a reference of some contigs serves reads of those whose header names a whole
/// genome. Lacking contigs with reads is an error for CRAM, whose reads there cannot be decoded,
/// and for BAM when they hold every read or more than [`MISSING_SHARE`] of them (by the index's
/// counts, else its bytes); fewer are a warning, which is returned.
///
/// # Errors
///
/// If the reference or the reads' index cannot be read, or the reference does not match.
pub fn contigs(reads: &Aligned, reference: &Path) -> Result<Option<String>> {
    let path = &reads.path;
    let fai = Layout::index(reference, "fai");
    let text = std::fs::read_to_string(&fai).with_context(|| format!("reading {}", fai.display()))?;
    let lengths: HashMap<&str, u64> = text
        .lines()
        .filter_map(|l| {
            let mut f = l.split('\t');
            Some((f.next()?, f.next()?.parse().ok()?))
        })
        .collect();
    let cram = lrpipe_reads::is_cram(path);
    // Reads on each contig: the index's counts, else the bytes it gives the contig.
    let amounts: Vec<u64> = match &reads.mapped {
        Some(counts) => counts.clone(),
        None => {
            let decoder = if cram { Some(Reference::open(reference)?) } else { None };
            let file = IndexedReads::open(path, decoder.as_ref())?;
            reads
                .header
                .reference_sequences()
                .iter()
                .map(|(name, c)| file.indexed_bytes(&String::from_utf8_lossy(name), 0, usize::from(c.length()) as i64))
                .collect::<io::Result<_>>()?
        }
    };
    let unit = if reads.mapped.is_some() { "reads" } else { "bytes of reads" };
    let (mut missing, mut other) = (Vec::new(), Vec::new());
    let mut declared: Vec<(String, Vec<u8>)> = Vec::new();
    for ((name, contig), &amount) in reads.header.reference_sequences().iter().zip(&amounts) {
        if amount == 0 {
            continue;
        }
        let name = String::from_utf8_lossy(name);
        let length = usize::from(contig.length()) as u64;
        let theirs = lengths.get(name.as_ref()).copied();
        if theirs == Some(length) {
            let tag = &sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
            if let Some(checksum) = contig.other_fields().get(tag) {
                let checksum: &[u8] = checksum.as_ref();
                declared.push((name.into_owned(), checksum.to_vec()));
            }
            continue;
        }
        match theirs {
            None => missing.push((name.into_owned(), amount)),
            Some(theirs) => other.push((name.into_owned(), theirs, length)),
        }
    }
    check_sequence_md5s(reference, path, &declared)?;
    let shown = reference.display();
    let file = path.display();
    let more = |n: usize| if n > 1 { format!(" and {} more with reads", n - 1) } else { String::new() };
    if let Some((name, theirs, length)) = other.first() {
        bail!(
            "contig {name} has {theirs} bases in the reference {shown}, but {length} in the header of {file}{}: is \
             it another assembly's? Give the reference the reads were aligned to",
            more(other.len())
        );
    }
    let total: u64 = amounts.iter().sum();
    let lacked: u64 = missing.iter().map(|&(_, n)| n).sum();
    let share = lacked as f64 / total.max(1) as f64;
    if !missing.is_empty() && !cram && lacked < total && share <= MISSING_SHARE {
        let mut named: Vec<String> = missing.iter().take(5).map(|(c, n)| format!("{c} {n}")).collect();
        if missing.len() > 5 {
            named.push(format!("and {} more", missing.len() - 5));
        }
        return Ok(Some(format!(
            "the reference {shown} lacks contigs that {:.2}% of the reads of {file} are on ({}, by {unit}): their \
             reads get no calls",
            share * 100.0,
            named.join(", ")
        )));
    }
    if let Some((first, _)) = missing.first() {
        let other_name = first.strip_prefix("chr").map_or_else(|| format!("chr{first}"), str::to_owned);
        let hint = if lengths.contains_key(other_name.as_str()) {
            format!(
                ": it has {other_name}, so it names its contigs {}, and the reads' header {}",
                if first.starts_with("chr") { "without `chr`" } else { "with `chr`" },
                if first.starts_with("chr") { "with it" } else { "without it" }
            )
        } else {
            String::new()
        };
        let why = if cram {
            "CRAM reads there cannot be decoded without it".to_owned()
        } else {
            format!("{:.1}% of its reads are on contigs it lacks", share * 100.0)
        };
        bail!(
            "the reference {shown} has no contig {first}{}, which {file} has reads on{hint} ({why}); give the \
             reference the reads were aligned to (contigs without reads need not be in it)",
            more(missing.len())
        );
    }
    Ok(None)
}

/// Contigs hashed at a time: each holds its contig's bases while it is, so more would hold
/// more of the reference than a run should before its stages start.
const MD5_THREADS: usize = 4;

/// Checks the reference's contigs against the checksums the reads' header declares for them
/// (`declared`: name and M5), several at a time: a whole genome's is 3 Gbp to read and hash
/// before every run, a resumed one too. The error is the first contig's that does not match, in
/// the header's order.
fn check_sequence_md5s(reference: &Path, reads: &Path, declared: &[(String, Vec<u8>)]) -> Result<()> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    if declared.is_empty() {
        return Ok(());
    }
    let next = AtomicUsize::new(0);
    let threads = MD5_THREADS.min(declared.len()).min(crate::cores::available().max(1));
    let mut failures: Vec<(usize, anyhow::Error)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| -> Option<(usize, anyhow::Error)> {
                    let fasta = match Reference::open(reference) {
                        Ok(fasta) => fasta,
                        Err(e) => return Some((0, e.into())),
                    };
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let (name, expected) = declared.get(i)?;
                        if let Err(e) = check_sequence_md5(&fasta, reference, reads, name, expected) {
                            return Some((i, e));
                        }
                    }
                })
            })
            .collect();
        workers.into_iter().filter_map(|w| w.join().expect("a checksum thread panicked")).collect()
    });
    failures.sort_by_key(|&(i, _)| i);
    match failures.into_iter().next() {
        Some((_, e)) => Err(e),
        None => Ok(()),
    }
}

fn check_sequence_md5(fasta: &Reference, reference: &Path, reads: &Path, name: &str, expected: &[u8]) -> Result<()> {
    let sequence = fasta.sequence(name)?.with_context(|| format!("reading {name} from {}", reference.display()))?;
    let mut md5 = Md5::new();
    let mut buffer = [0u8; 8192];
    for bases in sequence.as_ref().as_ref().chunks(buffer.len()) {
        let bytes = &mut buffer[..bases.len()];
        bytes.copy_from_slice(bases);
        bytes.make_ascii_uppercase();
        md5.update(&*bytes);
    }
    fasta.clear();
    let actual = crate::provenance::hex(&md5.finalize());
    if !expected.eq_ignore_ascii_case(actual.as_bytes()) {
        bail!(
            "contig {name} in the reference {} has MD5 {actual}, but {} declares M5 {}; give the reference the reads were aligned to",
            reference.display(),
            reads.display(),
            String::from_utf8_lossy(expected)
        );
    }
    Ok(())
}

/// The bonefish at `program` (a path, or a name on `PATH`), checked to be bonefish by its
/// `--version`, which prints `bonefish VERSION`: another program given in its place (`/bin/true`)
/// ran, wrote nothing, and failed later on what it had not written.
///
/// # Errors
///
/// If it cannot be found or run, or its `--version` is not bonefish's.
pub fn bonefish(program: &Path) -> Result<crate::provenance::ToolRecord> {
    let tool = crate::provenance::ToolRecord::of(program)
        .with_context(|| format!("running {} --version", program.display()))?;
    if !tool.version.starts_with("bonefish ") {
        bail!(
            "{} is not bonefish: its --version printed {:?}, not \"bonefish VERSION\"",
            tool.path.display(),
            tool.version
        );
    }
    Ok(tool)
}

/// Reads `--regions`' BED and checks it against the reference's `contigs` (each name and length,
/// in its order): every line a contig, start and end with `0 <= start < end`, on a contig of the
/// reference, ending within it, and at least one. Returns the regions sorted in the reference's
/// order, and merged where they overlap or touch, as every tool is given them.
///
/// # Errors
///
/// If the file cannot be read, or any of that is not so.
pub fn regions(bed: &Path, contigs: &[(String, u64)]) -> Result<Regions> {
    let read = Regions::read(bed).context("--regions")?;
    let mut sorted = Vec::with_capacity(read.len());
    for (contig, intervals) in read.iter() {
        let Some(&(_, length)) = contigs.iter().find(|(c, _)| c == contig) else {
            let other = contig.strip_prefix("chr").map_or_else(|| format!("chr{contig}"), str::to_owned);
            let hint = if contigs.iter().any(|(c, _)| *c == other) {
                format!(": it has {other}, so name the regions' contigs as it does")
            } else {
                ": are the regions of another assembly?".to_owned()
            };
            bail!(
                "--regions {}: contig {contig} is not in the reference ({} in its .fai){hint}",
                bed.display(),
                count_contigs(contigs.len())
            );
        };
        if let Some(&(start, end)) = intervals.last().filter(|&&(_, e)| e.unsigned_abs() > length) {
            bail!(
                "--regions {}: {contig}:{start}-{end} ends past {contig}'s {length} bases in the reference: are the \
                 regions of another assembly?",
                bed.display()
            );
        }
    }
    for (contig, _) in contigs {
        sorted.extend(read.contig(contig).iter().map(|&(s, e)| (contig.clone(), s, e)));
    }
    let regions = Regions::from_intervals(sorted);
    if regions.is_empty() {
        bail!("--regions {} has no regions (lines of a contig, a start and an end)", bed.display());
    }
    Ok(regions)
}

/// Checks `--regions`' contigs (`regions`, of the BED `bed`) against the aligned reads of
/// `classified` (unaligned reads are aligned to the reference's contigs): each must be in every
/// file's header, which the tools read them by. One without reads is not an error, as Clair3
/// takes it: its regions get no calls, which the returned warning says (by a BAM index's
/// counts; a CRAM's index does not count them).
///
/// # Errors
///
/// If a file's header lacks a region's contig.
pub fn regions_in_reads(regions: &Regions, bed: &Path, classified: &crate::input::Classified) -> Result<Vec<String>> {
    let mut without_reads = Vec::new();
    for (contig, _) in regions.iter() {
        let mut reads = Some(0u64);
        for file in &classified.aligned {
            let Some(index) = file.header.reference_sequences().get_index_of(contig.as_bytes()) else {
                bail!(
                    "--regions {}: {} has no contig {contig} (its header does not name it), which the tools cannot read",
                    bed.display(),
                    file.path.display()
                );
            };
            reads = reads.zip(file.mapped.as_ref().and_then(|m| m.get(index).copied())).map(|(a, b)| a + b);
        }
        if !classified.aligned.is_empty() && reads == Some(0) {
            without_reads.push(contig.to_owned());
        }
    }
    Ok(match without_reads.as_slice() {
        [] => Vec::new(),
        contigs => vec![format!(
            "--regions on {} get no calls: the reads have none there (by their index)",
            contigs.join(", ")
        )],
    })
}

/// "1 contig", "24 contigs".
fn count_contigs(n: usize) -> String {
    format!("{n} contig{}", if n == 1 { "" } else { "s" })
}

/// Checks `[snv] device` as lrpipe-snv reads it (`cpu`, `cuda`, `cuda:N` or `metal`).
///
/// # Errors
///
/// If lrpipe-snv would refuse it: a name it does not take, or a device this build cannot use
/// (CUDA without the `gpu` feature, Metal outside a macOS build with the `metal` feature), which
/// lrpipe-snv refuses only when its stage starts, after the reads are aligned.
pub fn device(device: &str) -> Result<()> {
    let Ok(parsed) = device.parse::<snv_infer::InferenceDevice>() else {
        bail!(
            "the site configuration's [snv] device = {device:?}: lrpipe snv takes \"cpu\", \"cuda\" (the first GPU), \"cuda:N\" or \"metal\" (an Apple GPU)"
        );
    };
    snv_infer::ModelConfig { device: parsed, precision: snv_infer::Precision::Float32 }
        .check()
        .with_context(|| format!("the site configuration's [snv] device = {device:?}"))?;
    Ok(())
}

/// Checks that lrpipe-snv's model, given by `by` (`--model`, `snv.ont_model`, ...), is a
/// directory holding its two networks; `models` is the site's directory of Clair3's models, whose
/// models are listed when a name in it is not one.
///
/// # Errors
///
/// If the directory or a network is missing.
pub fn model(path: &Path, by: &str, models: Option<&Path>) -> Result<()> {
    let shown = path.display();
    if !path.is_dir() {
        let held = models.filter(|m| path.parent() == Some(*m)).map(|m| {
            let mut names: Vec<String> = std::fs::read_dir(m)
                .into_iter()
                .flatten()
                .flatten()
                .filter(|e| crate::snv_model::network(&e.path(), "pileup").is_file())
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .collect();
            names.sort();
            format!("; {} holds {}", m.display(), if names.is_empty() { "none".into() } else { names.join(", ") })
        });
        bail!("the snv model {shown} (given by {by}) does not exist{}", held.unwrap_or_default());
    }
    for network in ["pileup", "full_alignment"] {
        let file = crate::snv_model::network(path, network);
        if !file.is_file() {
            bail!(
                "the snv model {shown} (given by {by}) has no {network}.pt: a Clair3 model directory holds pileup.pt and full_alignment.pt"
            );
        }
        if !complete_checkpoint(&file).with_context(|| format!("reading {}", file.display()))? {
            bail!(
                "the snv model {shown} (given by {by}) has a {} cut short (a zip archive, as PyTorch saves one, without \
                 its end): fetch the model again (`lrpipe setup --verify` checks a setup's)",
                file.file_name().map_or_else(|| network.into(), |n| n.to_string_lossy()),
            );
        }
    }
    Ok(())
}

/// Bytes from the end of a zip archive its end record (22 bytes, and a comment of up to 64 KiB)
/// starts within.
const ZIP_END_WITHIN: u64 = 22 + 65_535;

/// Whether the network file at `path` is whole: a zip archive (as PyTorch saves a checkpoint, and
/// Clair3's models are) ending with its end record, which a copy or download cut short lacks. A file
/// that is not a zip archive is left to lrpipe-snv, which says what it is.
fn complete_checkpoint(path: &Path) -> io::Result<bool> {
    use std::io::{Seek, SeekFrom};
    let mut file = File::open(path)?;
    let mut magic = [0; 4];
    if read_up_to(&mut file, &mut magic)? < 4 || magic != *b"PK\x03\x04" {
        return Ok(true);
    }
    let len = file.metadata()?.len();
    file.seek(SeekFrom::Start(len.saturating_sub(ZIP_END_WITHIN)))?;
    let mut tail = Vec::new();
    file.read_to_end(&mut tail)?;
    // The end record: its signature, 18 bytes, then a comment of the length its last two give,
    // which ends the file.
    Ok((0..tail.len().saturating_sub(21)).rev().any(|at| {
        let end = &tail[at..];
        end.starts_with(b"PK\x05\x06") && end.len() == 22 + usize::from(u16::from_le_bytes([end[20], end[21]]))
    }))
}

/// `--memory`: GiB, as `32`, `32G`, `32GB` or `32GiB` (all GiB), or with `M`/`MiB` or `T`/`TiB`.
///
/// # Errors
///
/// If it is not a number with one of those units, or less than 1 GiB, which no run fits in.
pub fn memory_gib(text: &str) -> Result<f64, String> {
    let t = text.trim();
    let split = t.find(|c: char| !(c.is_ascii_digit() || c == '.')).unwrap_or(t.len());
    let (number, unit) = t.split_at(split);
    let scale = match unit.trim().to_ascii_lowercase().as_str() {
        "" | "g" | "gb" | "gib" => 1.0,
        "m" | "mb" | "mib" => 1.0 / 1024.0,
        "t" | "tb" | "tib" => 1024.0,
        _ => return Err(format!("{text:?} is not an amount of memory: give GiB, as 32, 32G or 32GiB")),
    };
    let gib = number
        .parse::<f64>()
        .map_err(|_| format!("{text:?} is not an amount of memory: give GiB, as 32, 32G or 32GiB"))?
        * scale;
    if gib < 1.0 || !gib.is_finite() {
        return Err(format!("{text} is less than the 1 GiB a run needs at the least (the unit is GiB: 32 is 32 GiB)"));
    }
    Ok(gib)
}

/// `--threads`: at least one.
///
/// # Errors
///
/// If it is not a whole number, or is 0.
pub fn threads(text: &str) -> Result<usize, String> {
    match text.parse::<usize>() {
        Ok(0) => Err("a run needs at least 1 thread".into()),
        Ok(n) => Ok(n),
        Err(e) => Err(format!("{text:?}: {e}")),
    }
}

/// A warning when `sample` is not the sample the reads' read groups name (their `SM`, of every
/// file of `header`, the reads' merged): the outputs name `sample`, which tools reading the
/// reads' own header (WhatsHap's `haplotag`) do not match.
#[must_use]
pub fn sample_warning(header: &sam::Header, sample: &str) -> Option<String> {
    use sam::header::record::value::map::read_group::tag::SAMPLE;
    let mut theirs: Vec<String> = Vec::new();
    for group in header.read_groups().values() {
        if let Some(sm) = group.other_fields().get(&SAMPLE).map(ToString::to_string)
            && !theirs.contains(&sm)
        {
            theirs.push(sm);
        }
    }
    if theirs.is_empty() || theirs.iter().any(|s| s == sample) {
        return None;
    }
    let named = theirs.iter().map(|s| format!("SM:{s}")).collect::<Vec<_>>().join(", ");
    Some(format!(
        "--sample {sample} is not the sample the reads' read groups name ({named}): the outputs name {sample}, which \
         tools that read the reads' header do not take for theirs"
    ))
}

/// A warning naming the outputs a run of this layout can write that are in its output directory
/// but that none of `stages` writes: an earlier run's with other analyses (a full run, then
/// `--analyses snv`), other input (unaligned reads aligned into `align/`, then aligned reads) or
/// options (`--haplotagged-reads`, then a run without it), which are not this run's and are not
/// in its manifest. They are left, not removed: they may be all that is left of the earlier run,
/// and a BAM of reads may be one that other work reads. Haplotagged reads older than this run's
/// haplotag list are said to be from an earlier phasing.
#[must_use]
pub fn leftovers(layout: &Layout, stages: &[Stage]) -> Option<String> {
    // What this run reads is its own too: an earlier run's aligned BAM given as --reads.
    let written = |path: &PathBuf| stages.iter().any(|s| s.outputs.contains(path) || s.inputs.contains(path));
    let left: Vec<PathBuf> = [
        layout.aligned_bam(),
        layout.merged_reads(),
        layout.haplotagged_reads(),
        layout.merged_vcf(),
        layout.phased_snv_vcf(),
        layout.snv_gvcf(),
        layout.sv_vcf(),
        layout.snf(),
        layout.depth_windows(),
        layout.depth_summary(),
        layout.cnv_vcf(),
        layout.cnv_bed(),
        layout.repeats_vcf(),
        layout.bedmethyl(),
        layout.cpg_table("combined"),
        layout.cpg_table("hap1"),
        layout.cpg_table("hap2"),
        layout.read_stats(),
        layout.qc_json(),
        layout.qc_tsv(),
    ]
    .into_iter()
    .filter(|p| !written(p) && p.exists())
    .collect();
    if left.is_empty() {
        return None;
    }
    let names: Vec<String> =
        left.iter().map(|p| p.strip_prefix(layout.root()).unwrap_or(p).display().to_string()).collect();
    let modified = |path: &Path| std::fs::metadata(path).and_then(|m| m.modified()).ok();
    let stale = left.contains(&layout.haplotagged_reads())
        && written(&layout.haplotags())
        && modified(&layout.haplotagged_reads()) < modified(&layout.haplotags());
    Some(format!(
        "{} holds {} of an earlier run that this one does not write ({}): {} not this run's results and not in its \
         manifest{}",
        layout.root().display(),
        if left.len() == 1 { "an output" } else { "outputs" },
        names.join(", "),
        if left.len() == 1 { "it is" } else { "they are" },
        if stale {
            format!(
                "; {} is older than this run's haplotag list, so its HP tags are an earlier phasing's (give \
                 --haplotagged-reads to write it anew)",
                names[left.iter().position(|p| *p == layout.haplotagged_reads()).expect("left")]
            )
        } else {
            String::new()
        }
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_device_this_build_cannot_use_is_refused_before_the_run() {
        assert!(device("cpu").is_ok());
        assert!(device("gpu").unwrap_err().to_string().contains("lrpipe snv takes"));
        let metal = device("metal");
        if cfg!(all(feature = "metal", target_os = "macos")) {
            assert!(metal.is_ok());
        } else {
            assert!(format!("{:#}", metal.unwrap_err()).contains("no Metal backend"), "metal refused at once");
        }
        assert_eq!(device("cuda:0").is_ok(), cfg!(feature = "gpu"));
    }

    #[test]
    fn reads_and_haplotagged_reads_of_an_earlier_run_are_leftovers() {
        // An earlier run of unaligned reads left align/S.bam, and one with --haplotagged-reads its
        // BAM, older than the list this run's phasing wrote.
        let root = std::env::temp_dir().join(format!("lrpipe-checks-leftovers-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let layout = Layout::new(&root, "S");
        for path in [layout.aligned_bam(), layout.merged_reads(), layout.haplotagged_reads(), layout.qc_json()] {
            std::fs::create_dir_all(path.parent().expect("a directory")).expect("directory");
            std::fs::write(path, "reads").expect("written");
        }
        std::fs::File::options()
            .write(true)
            .open(layout.haplotagged_reads())
            .expect("earlier BAM")
            .set_modified(std::time::SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(1_000_000_000))
            .expect("set the earlier timestamp");
        std::fs::write(layout.haplotags(), "list").expect("written");
        let haplotag = Stage {
            name: crate::plan::StageName::Haplotag,
            program: "lrpipe".into(),
            args: Vec::new(),
            inputs: Vec::new(),
            outputs: vec![layout.haplotags()],
            after: Vec::new(),
            memory_gib: 1.0,
            priority: crate::plan::Priority::Critical,
            scratch: None,
        };
        let w = leftovers(&layout, std::slice::from_ref(&haplotag)).expect("a warning");
        assert!(w.contains("(align/S.bam, reads/S.bam, phase/S.haplotagged.bam, qc/S.qc.json)"), "{w}");
        assert!(w.contains("phase/S.haplotagged.bam is older than this run's haplotag list"), "{w}");
        // Without this run's list, no word of its age.
        let w = leftovers(&layout, &[]).expect("a warning");
        assert!(!w.contains("is older than this run's haplotag list"), "{w}");
        // An earlier run's aligned BAM that this run reads (given as --reads) is not a leftover.
        let reading = Stage { inputs: vec![layout.aligned_bam()], ..haplotag.clone() };
        let w = leftovers(&layout, &[reading]).expect("a warning");
        assert!(!w.contains("align/S.bam"), "{w}");
        let _ = std::fs::remove_dir_all(root);
    }

    #[test]
    fn a_cram_index_of_another_version_of_the_file_is_found() {
        // ONT's published indexes of their HG002 and HG004 CRAMs are 63 and 68 bytes off, as the
        // index of a version with a longer header is.
        use std::io::Write as _;
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-sv/tests/data");
        let cram = data.join("sample1.cram");
        let crai = data.join("sample1.cram.crai");
        let own = cram::crai::fs::read(&crai).expect("the fixture's index");
        assert_eq!(misplaced_container(&cram, &own).expect("read"), None);
        let mut text = String::new();
        flate2::read::GzDecoder::new(File::open(&crai).expect("index")).read_to_string(&mut text).expect("gzip");
        let shifted: String = text
            .lines()
            .map(|line| {
                let mut f: Vec<String> = line.split('\t').map(str::to_owned).collect();
                f[3] = (f[3].parse::<u64>().expect("an offset") + 68).to_string();
                f.join("\t") + "\n"
            })
            .collect();
        let dir = std::env::temp_dir().join(format!("lrpipe-checks-crai-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let other = dir.join("shifted.crai");
        let mut gz =
            flate2::write::GzEncoder::new(File::create(&other).expect("created"), flate2::Compression::default());
        gz.write_all(shifted.as_bytes()).expect("written");
        gz.finish().expect("finished");
        let shifted = cram::crai::fs::read(&other).expect("the shifted index");
        assert_eq!(misplaced_container(&cram, &shifted).expect("read"), Some(own[0].offset() + 68));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn regions_are_checked_against_the_reference_and_sorted_in_its_order() {
        let dir = std::env::temp_dir().join(format!("lrpipe-checks-regions-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let bed = dir.join("regions.bed");
        let contigs: Vec<(String, u64)> =
            [("chr1", 1000), ("chr2", 500), ("chrX", 800)].into_iter().map(|(c, l)| (c.to_owned(), l)).collect();
        let check = |text: &str| {
            std::fs::write(&bed, text).expect("written");
            regions(&bed, &contigs)
        };
        // In the reference's order, merged, whatever the file's order.
        let sorted = check("chrX\t0\t800\nchr2\t100\t200\nchr1\t5\t10\nchr2\t150\t300\n").expect("valid");
        assert_eq!(sorted.bed(), "chr1\t5\t10\nchr2\t100\t300\nchrX\t0\t800\n");
        for (text, wanted) in [
            ("chr3\t0\t10\n", "contig chr3 is not in the reference (3 contigs in its .fai): are the regions of"),
            ("X\t0\t10\n", "it has chrX, so name the regions' contigs as it does"),
            ("chr2\t400\t501\n", "chr2:400-501 ends past chr2's 500 bases"),
            ("chr1\t10\t5\n", "line 1: the interval 10-5 ends before it starts"),
            ("# nothing\n", "has no regions"),
        ] {
            let e = format!("{:#}", check(text).expect_err(text));
            assert!(e.contains(wanted), "{text}: {e}");
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn sample_names_are_safe_file_names() {
        for good in ["HG003", "hg003", "NA12878_rep-2.1"] {
            sample_name(good).expect(good);
        }
        sample_name(&"a".repeat(SAMPLE_NAME_MAX)).expect("the longest");
        let e = sample_name(&"a".repeat(SAMPLE_NAME_MAX + 1)).expect_err("too long").to_string();
        assert!(e.contains("at most 200"), "{e}");
        for (bad, why) in [
            ("../evil", "'/'"),
            ("a/b", "'/'"),
            ("", "empty"),
            ("HG 003", "' '"),
            ("HG\t003", "'\\t'"),
            ("..", "starts with \".\""),
            (".hidden", "starts with \".\""),
            ("-x", "starts with \"-\""),
            ("HG003\u{7}", "'\\u{7}'"),
        ] {
            let e = sample_name(bad).expect_err(bad).to_string();
            assert!(e.contains(why) && e.contains("letters, digits"), "{bad:?}: {e}");
        }
    }

    /// `path`, a BAM, checked as [`crate::input::classify`] checks one.
    fn checked(path: &Path) -> Result<Aligned> {
        aligned(path, crate::input::read_header(path).expect("a header"))
    }

    /// A BAM of `contigs` (name, length) and no reads, with `order` as its `@HD SO`, and a BAI
    /// counting `mapped` reads on each; in a scratch directory named for `test`.
    fn bam(test: &str, order: &str, contigs: &[(&str, u32)], mapped: &[u64]) -> PathBuf {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("lrpipe-checks-{test}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let path = dir.join("reads.bam");
        let text = format!("@HD\tVN:1.6\tSO:{order}\n");
        let mut bam = pasta::bgzf::io::Writer::new(File::create(&path).expect("BAM"));
        let mut h = b"BAM\x01".to_vec();
        h.extend(i32::try_from(text.len()).expect("short").to_le_bytes());
        h.extend(text.as_bytes());
        h.extend(i32::try_from(contigs.len()).expect("few").to_le_bytes());
        for (name, len) in contigs {
            h.extend(i32::try_from(name.len() + 1).expect("short").to_le_bytes());
            h.extend(name.as_bytes());
            h.push(0);
            h.extend(len.to_le_bytes());
        }
        bam.write_all(&h).expect("written");
        bam.finish().expect("BAM written");
        let mut bai = b"BAI\x01".to_vec();
        bai.extend(i32::try_from(contigs.len()).expect("few").to_le_bytes());
        for &n in mapped {
            // The metadata pseudo-bin alone: the reads' span, then mapped and unmapped counts.
            bai.extend(1i32.to_le_bytes());
            bai.extend(37_450u32.to_le_bytes());
            bai.extend(2i32.to_le_bytes());
            bai.extend([0u8; 16]);
            bai.extend(n.to_le_bytes());
            bai.extend(0u64.to_le_bytes());
            bai.extend(0i32.to_le_bytes());
        }
        std::fs::write(dir.join("reads.bam.bai"), bai).expect("written");
        path
    }

    /// A reference `.fai` of `contigs` beside `reads` (the FASTA is not read).
    fn reference(reads: &Path, contigs: &[(&str, u32)]) -> PathBuf {
        let path = reads.with_file_name("ref.fa");
        std::fs::write(&path, ">x\n").expect("written");
        let fai: String = contigs.iter().map(|(n, l)| format!("{n}\t{l}\t3\t60\t61\n")).collect();
        std::fs::write(Layout::index(&path, "fai"), fai).expect("written");
        path
    }

    #[test]
    fn a_few_reads_on_contigs_the_reference_lacks_are_a_warning_most_an_error() {
        // HLA, decoy and EBV contigs of the reads' header, which the no-alt analysis set lacks,
        // hold a few reads: a warning. Most of the reads there, or all, is the wrong reference.
        let header = [("chr20", 1000), ("chrEBV", 100), ("chrY", 500)];
        let check = |test: &str, mapped: &[u64], ours: &[(&str, u32)]| {
            let path = bam(test, "coordinate", &header, mapped);
            let reads = checked(&path).expect("good reads");
            let found = contigs(&reads, &reference(&path, ours));
            let _ = std::fs::remove_dir_all(path.parent().expect("a directory"));
            found
        };
        let w = check("few", &[1000, 5, 0], &[("chr20", 1000)]).expect("a warning").expect("not silent");
        assert!(w.contains("0.50% of the reads of") && w.contains("chrEBV 5, by reads"), "{w}");
        assert_eq!(check("none", &[1000, 0, 0], &[("chr20", 1000)]).expect("fine"), None);
        let e = check("many", &[1000, 100, 0], &[("chr20", 1000)]).expect_err("too many").to_string();
        assert!(e.contains("has no contig chrEBV") && e.contains("9.1% of its reads"), "{e}");
        let e = check("named", &[1000, 0, 0], &[("20", 1000)]).expect_err("all missing").to_string();
        assert!(e.contains("names its contigs without `chr`"), "{e}");
        let e = check("length", &[1000, 0, 0], &[("chr20", 999)]).expect_err("another assembly").to_string();
        assert!(e.contains("contig chr20 has 999 bases in the reference"), "{e}");
    }

    #[test]
    fn unsorted_reads_with_an_index_are_a_warning_by_name_an_error() {
        let header = [("chr20", 1000)];
        let path = bam("unsorted", "unsorted", &header, &[10]);
        let found = checked(&path).expect("taken as sorted");
        assert!(found.warnings.iter().any(|w| w.contains("SO:unsorted")), "{:?}", found.warnings);
        assert_eq!(found.mapped, Some(vec![10]));
        let path = bam("by-name", "queryname", &header, &[10]);
        let e = checked(&path).expect_err("by name").to_string();
        assert!(e.contains("sorted by read name"), "{e}");
        let _ = std::fs::remove_dir_all(path.parent().expect("a directory"));
    }

    /// A zip archive of one file, `name` holding `data` stored uncompressed, as PyTorch saves a
    /// checkpoint: its local header and data, its central directory and its end record.
    fn stored_zip(name: &str, data: &[u8]) -> Vec<u8> {
        let mut crc = flate2::Crc::new();
        crc.update(data);
        let (crc, size, name_len) = (crc.sum(), data.len() as u32, name.len() as u16);
        let mut zip = Vec::new();
        zip.extend(b"PK\x03\x04");
        zip.extend([20, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        zip.extend(crc.to_le_bytes());
        zip.extend(size.to_le_bytes());
        zip.extend(size.to_le_bytes());
        zip.extend(name_len.to_le_bytes());
        zip.extend(0u16.to_le_bytes());
        zip.extend(name.as_bytes());
        zip.extend(data);
        let directory = zip.len() as u32;
        zip.extend(b"PK\x01\x02");
        zip.extend([20, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        zip.extend(crc.to_le_bytes());
        zip.extend(size.to_le_bytes());
        zip.extend(size.to_le_bytes());
        zip.extend(name_len.to_le_bytes());
        zip.extend([0; 12]);
        zip.extend(0u32.to_le_bytes());
        zip.extend(name.as_bytes());
        let directory_len = zip.len() as u32 - directory;
        zip.extend(b"PK\x05\x06");
        zip.extend([0, 0, 0, 0, 1, 0, 1, 0]);
        zip.extend(directory_len.to_le_bytes());
        zip.extend(directory.to_le_bytes());
        zip.extend(0u16.to_le_bytes());
        zip
    }

    #[test]
    fn a_model_takes_its_networks_as_lrpipe_snv_finds_them() {
        // lrpipe-snv takes `pileup` when there is no `pileup.pt`.
        let dir = std::env::temp_dir().join(format!("lrpipe-checks-model-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        std::fs::write(dir.join("pileup"), "").expect("written");
        std::fs::write(dir.join("full_alignment.pt"), "").expect("written");
        model(&dir, "--model", None).expect("bare pileup taken");
        std::fs::remove_file(dir.join("pileup")).expect("removed");
        let e = model(&dir, "--model", None).expect_err("no pileup").to_string();
        assert!(e.contains("has no pileup.pt"), "{e}");
        // A checkpoint cut short (a zip archive without its end record), as a copy stopped part
        // way leaves it, passed planning and failed in lrpipe-snv.
        let zip = stored_zip("archive/data.pkl", &[7; 300]);
        std::fs::write(dir.join("pileup.pt"), &zip).expect("written");
        model(&dir, "--model", None).expect("whole");
        for cut in [zip.len() / 2, zip.len() - 1] {
            std::fs::write(dir.join("pileup.pt"), &zip[..cut]).expect("written");
            let e = model(&dir, "--model", None).expect_err("cut short").to_string();
            assert!(e.contains("has a pileup.pt cut short"), "{e}");
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn memory_takes_gib_with_or_without_a_unit() {
        assert_eq!(memory_gib("8"), Ok(8.0));
        assert_eq!(memory_gib("8G"), Ok(8.0));
        assert_eq!(memory_gib("8GiB"), Ok(8.0));
        assert_eq!(memory_gib("8gb"), Ok(8.0));
        assert_eq!(memory_gib("1.5"), Ok(1.5));
        assert_eq!(memory_gib("2048M"), Ok(2.0));
        assert_eq!(memory_gib("1T"), Ok(1024.0));
        for small in ["0", "0.5", "512M"] {
            assert!(memory_gib(small).expect_err(small).contains("1 GiB"), "{small}");
        }
        for bad in ["eight", "8 parsecs", "", "-4"] {
            assert!(memory_gib(bad).is_err(), "{bad}");
        }
        assert!(threads("0").is_err());
        assert_eq!(threads("4"), Ok(4));
    }

    #[test]
    fn devices_are_read_as_lrpipe_snv_reads_them() {
        // Names lrpipe-snv takes; whether this build can use a GPU is the next test's.
        for good in ["cpu", "cuda", "cuda:1", "metal"] {
            good.parse::<snv_infer::InferenceDevice>().expect(good);
        }
        device("cpu").expect("cpu");
        let e = device("gpu").expect_err("gpu").to_string();
        assert!(e.contains("\"gpu\"") && e.contains("cuda:N"), "{e}");
    }

    #[test]
    fn a_sample_other_than_the_read_groups_is_warned_of() {
        let header = "@HD\tVN:1.6\n@RG\tID:a\tPL:ONT\tSM:hg003\n";
        let parse = |text: &str| text.parse::<sam::Header>().expect("a header");
        let w = sample_warning(&parse(header), "HG003").expect("a warning");
        assert!(w.contains("SM:hg003") && w.contains("--sample HG003"), "{w}");
        assert_eq!(sample_warning(&parse(header), "hg003"), None);
        let none = parse("@HD\tVN:1.6\n@RG\tID:a\n");
        assert_eq!(sample_warning(&none, "HG003"), None, "no SM, nothing to differ from");
    }
}

//! Genotyping a whole catalog: loci grouped into chunks of nearby loci, chunks handed to
//! threads (each with its own alignment and reference readers), records written in catalog
//! order.

use std::fs::File;
use std::io::{BufRead, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use lrpipe_reads::par::{self, Build};
use pasta::core::{Position, Region};
use pasta::fasta;

use crate::Error;
use crate::catalog::Locus;
use crate::extract::{ExtractOptions, RefWindow, allele};
use crate::genotype::{GenotypeOptions, Observation, Ploidy, call, one_per_read};
use crate::reads::{Read, Reads, SharedReference, is_cram};
use crate::vcf;

/// The sample's sex chromosomes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Karyotype {
    /// Two copies of chrX, no chrY.
    #[default]
    Xx,
    /// One copy of chrX outside its pseudoautosomal regions, and of chrY.
    Xy,
}

/// The files of a run and the sample's name.
#[derive(Clone, Debug)]
pub struct Job {
    /// Indexed BAM or CRAM.
    pub reads: PathBuf,
    /// Reference FASTA with a `.fai` index; also decodes CRAM.
    pub genome: PathBuf,
    /// Output VCF, bgzipped when it ends in `.gz` (with a tabix index beside it, `.tbi`, or
    /// `.csi` for contigs longer than tabix's 2^29 bases).
    pub output: PathBuf,
    /// Sample name in the VCF header.
    pub sample: String,
    /// The reads' haplotypes from a haplotag list, in place of the alignments' `HP` tags.
    pub haplotags: Option<std::sync::Arc<lrpipe_reads::haplotags::Haplotags>>,
}

/// Settings of a run.
#[derive(Clone, Debug, PartialEq)]
pub struct RunOptions {
    /// Worker threads, each with its own readers.
    pub threads: usize,
    /// Reads with a lower mapping quality are not used.
    pub min_mapq: u8,
    /// Copies of chrX and chrY.
    pub karyotype: Karyotype,
    /// At most this many spanning reads per locus are used, spread evenly over them in read
    /// order.
    pub max_reads: usize,
    /// How reads' repeat sequences are cut out.
    pub extract: ExtractOptions,
    /// How alleles are called.
    pub genotype: GenotypeOptions,
    /// chrX's pseudoautosomal regions, two copies in XY samples too ([`x_pars`] reads them off
    /// the reference); GRCh38's by default.
    pub x_pars: Vec<(i64, i64)>,
}

impl Default for RunOptions {
    fn default() -> Self {
        RunOptions {
            threads: 8,
            min_mapq: 1,
            karyotype: Karyotype::Xx,
            max_reads: 200,
            extract: ExtractOptions::default(),
            genotype: GenotypeOptions::default(),
            x_pars: Build::Grch38.x().to_vec(),
        }
    }
}

/// What a run did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    /// Loci genotyped.
    pub loci: usize,
    /// Loci with a call (enough spanning reads).
    pub called: usize,
}

/// chrX's pseudoautosomal regions in the reference `genome` (by chrX's length in its `.fai`):
/// GRCh38's, GRCh37's or CHM13's. They were GRCh38's whatever the reference (LRP-010).
///
/// Another chrX has no known PARs, so in an XY sample all of it is one copy, as lrpipe-cnv
/// counts it; a warning says so. It was an error, which failed `lrpipe run` at its repeats after
/// its dry run had passed (lrpipe#44).
///
/// # Errors
///
/// If the `.fai` cannot be read.
pub fn x_pars(genome: &Path, karyotype: Karyotype) -> Result<Vec<(i64, i64)>, Error> {
    let fai_path = genome.with_added_extension("fai");
    let fai = fasta::fai::fs::read(&fai_path).map_err(Error::io(&fai_path))?;
    let x = fai.as_ref().iter().find(|r| matches!(r.name().as_ref() as &[u8], b"chrX" | b"X")).map(|r| r.length());
    Ok(match (x, x.and_then(Build::from_x_length), karyotype) {
        (_, Some(build), _) => build.x().to_vec(),
        (Some(len), None, Karyotype::Xy) => {
            eprintln!(
                "warning: chrX of {} is {len} bp, which is not GRCh38's, GRCh37's or CHM13's: its pseudoautosomal \
                 regions are unknown, so all of chrX is one copy in this XY sample",
                genome.display()
            );
            Vec::new()
        }
        // No chrX, or an XX sample, which never uses them.
        _ => Build::Grch38.x().to_vec(),
    })
}

/// A chunk ends before a locus this far past the previous one: reads are fetched once per
/// chunk, so a larger gap would fetch reads no locus uses.
const CHUNK_MAX_GAP: i64 = 10_000;
/// Loci per chunk at most, so work spreads over threads evenly.
const CHUNK_MAX_LOCI: usize = 500;
/// A chunk spans at most this many bases, which bounds the reads held per thread.
const CHUNK_MAX_SPAN: i64 = 1_000_000;

/// Copies of `locus` expected (0 for chrY in XX samples), with chrX's PARs `x_pars`.
fn copies(locus: &Locus, karyotype: Karyotype, x_pars: &[(i64, i64)]) -> usize {
    let c = locus.contig.trim_start_matches("chr");
    match (c, karyotype) {
        ("M" | "MT", _) => 1,
        ("Y", Karyotype::Xx) => 0,
        ("Y", Karyotype::Xy) => 1,
        ("X", Karyotype::Xy) if !par::within(x_pars, locus.start, locus.end) => 1,
        _ => 2,
    }
}

/// Consecutive loci on one contig, close together.
fn chunks(loci: &[Locus]) -> Vec<Range<usize>> {
    let mut out = Vec::new();
    let mut first = 0;
    for i in 1..=loci.len() {
        let split = i == loci.len()
            || loci[i].contig != loci[first].contig
            || loci[i].start - loci[i - 1].end > CHUNK_MAX_GAP
            || i - first >= CHUNK_MAX_LOCI
            || loci[i].end - loci[first].start > CHUNK_MAX_SPAN;
        if split {
            out.push(first..i);
            first = i;
        }
    }
    out
}

/// Sorts loci by the reference's contig order, then start; drops loci on contigs it lacks.
pub fn sort_loci(loci: &mut Vec<Locus>, contigs: &[(String, u64)]) {
    loci.retain(|l| contigs.iter().any(|(c, _)| c == &l.contig));
    let rank = |c: &str| contigs.iter().position(|(n, _)| n == c).unwrap_or(usize::MAX);
    loci.sort_by_key(|a| (rank(&a.contig), a.start, a.end));
}

/// Genotypes `loci` (sorted by contig and start, as [`sort_loci`] leaves them) and writes the
/// VCF.
///
/// # Errors
///
/// If an input cannot be read or the output written.
///
/// # Panics
///
/// If a worker thread panics.
pub fn run(job: &Job, loci: &[Locus], opts: &RunOptions) -> Result<Summary, Error> {
    let fai_path = job.genome.with_added_extension("fai");
    let fai = fasta::fai::fs::read(&fai_path).map_err(Error::io(&fai_path))?;
    let contigs: Vec<(String, u64)> = fai.as_ref().iter().map(|r| (r.name().to_string(), r.length())).collect();
    let (text, called) = records_with(&job.reads, &job.genome, loci, opts, job.haplotags.as_ref())?;
    let out = &job.output;
    let file = File::create(out).map_err(Error::io(out))?;
    let header = vcf::header(&contigs, &job.sample);
    if out.extension().is_some_and(|e| e == "gz") {
        let mut w = pasta::bgzf::io::Writer::new(std::io::BufWriter::new(file));
        w.write_all(header.as_bytes()).map_err(Error::io(out))?;
        w.write_all(text.as_bytes()).map_err(Error::io(out))?;
        w.finish().and_then(|mut f| f.flush()).map_err(Error::io(out))?;
        let longest = contigs.iter().map(|&(_, len)| len).max().unwrap_or(0);
        index_vcf(out, longest)?;
    } else {
        let mut w = std::io::BufWriter::new(file);
        w.write_all(header.as_bytes()).map_err(Error::io(out))?;
        w.write_all(text.as_bytes()).map_err(Error::io(out))?;
        w.flush().map_err(Error::io(out))?;
    }
    Ok(Summary { loci: loci.len(), called })
}

/// The longest contig a tabix (TBI) index holds: 2^29 bases.
const TBI_MAX_LEN: u64 = 1 << 29;

/// Indexes the bgzipped VCF at `path`, of contigs up to `longest` bases, as `tabix -p vcf` does:
/// into `PATH.tbi`, or `PATH.csi` (with `tabix --csi`'s minimum shift, 14) for a longer contig,
/// removing an index of the other kind left beside it.
fn index_vcf(path: &Path, longest: u64) -> Result<(), Error> {
    use pasta::variant::index::tabix::{TabixBuilder, TabixConfig};
    let (min_shift, suffix, other) = if longest > TBI_MAX_LEN { (14, "csi", "tbi") } else { (0, "tbi", "csi") };
    let mut reader = pasta::bgzf::io::Reader::new(File::open(path).map_err(Error::io(path))?);
    let mut builder = TabixBuilder::new(TabixConfig::VCF, min_shift);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line).map_err(Error::io(path))? == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        builder.push_line(&line, u64::from(reader.virtual_position())).map_err(Error::io(path))?;
    }
    let format = builder.format();
    let tabix = builder.finish(u64::from(reader.virtual_position()));
    let index = path.with_added_extension(suffix);
    tabix.index.save(&index, format).map_err(Error::io(&index))?;
    // An index of the other kind from an earlier run would describe an older VCF.
    let other = path.with_added_extension(other);
    match std::fs::remove_file(&other) {
        Err(e) if e.kind() != std::io::ErrorKind::NotFound => Err(Error::io(&other)(e)),
        _ => Ok(()),
    }
}

/// The VCF records of `loci` (sorted as [`sort_loci`] leaves them) from the reads of `reads`,
/// decoded (when CRAM) with `genome`, and how many loci were called.
///
/// # Errors
///
/// If an input cannot be read.
///
/// # Panics
///
/// If a worker thread panics.
pub fn records(
    reads: &std::path::Path,
    genome: &std::path::Path,
    loci: &[Locus],
    opts: &RunOptions,
) -> Result<(String, usize), Error> {
    records_with(reads, genome, loci, opts, None)
}

/// [`records`], with the reads' haplotypes from `haplotags` in place of their `HP` tags when
/// given.
///
/// # Errors
///
/// If an input cannot be read.
///
/// # Panics
///
/// If a worker thread panics.
pub fn records_with(
    reads: &std::path::Path,
    genome: &std::path::Path,
    loci: &[Locus],
    opts: &RunOptions,
    haplotags: Option<&std::sync::Arc<lrpipe_reads::haplotags::Haplotags>>,
) -> Result<(String, usize), Error> {
    let chunks = chunks(loci);
    let results: Vec<Mutex<Option<String>>> = chunks.iter().map(|_| Mutex::new(None)).collect();
    let next = AtomicUsize::new(0);
    let called = AtomicUsize::new(0);
    let threads = opts.threads.max(1).min(chunks.len().max(1));
    let cram_reference = if is_cram(reads) { Some(SharedReference::open(genome)?) } else { None };
    std::thread::scope(|s| -> Result<(), Error> {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| -> Result<(), Error> {
                    let mut worker = Worker::open(reads, genome, cram_reference.as_ref(), haplotags)?;
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(range) = chunks.get(i) else { break };
                        let (text, n) = worker.chunk(&loci[range.clone()], opts)?;
                        called.fetch_add(n, Ordering::Relaxed);
                        *results[i].lock().expect("lock poisoned") = Some(text);
                    }
                    Ok(())
                })
            })
            .collect();
        for w in workers {
            w.join().expect("a worker thread panicked")?;
        }
        Ok(())
    })?;
    let text: String =
        results.into_iter().map(|r| r.into_inner().expect("lock poisoned").unwrap_or_default()).collect();
    Ok((text, called.into_inner()))
}

/// At most `max` of `obs`, spread evenly over them in their order: all of them when there are
/// no more. Reads that span a locus in a soft clip come after those that span it in their
/// alignment, so taking the first ones would leave them out of deep data.
pub fn evenly<T>(obs: Vec<T>, max: usize) -> Vec<T> {
    let n = obs.len();
    if n <= max {
        return obs;
    }
    let keep: Vec<usize> = (0..max).map(|i| i * n / max).collect();
    obs.into_iter().enumerate().filter(|(i, _)| keep.binary_search(i).is_ok()).map(|(_, o)| o).collect()
}

/// A read's observation of `locus`, when it spans it with both flanks.
pub fn observe(read: &Read, locus: &Locus, reference: &RefWindow, opts: &RunOptions) -> Option<Observation> {
    allele(read, locus, reference, &opts.extract).map(|seq| Observation {
        seq,
        haplotype: read.haplotype,
        phase_set: read.phase_set,
        read: read.name,
        supplementary: read.supplementary,
    })
}

/// One locus's VCF record from its reads' observations (in read order; the first
/// [`RunOptions::max_reads`] count), and whether it was called. `reference` must cover the
/// locus and the base before it.
pub fn genotype_locus(locus: &Locus, obs: &[Observation], reference: &RefWindow, opts: &RunOptions) -> (String, bool) {
    let copies = copies(locus, opts.karyotype, &opts.x_pars);
    let obs = &obs[..obs.len().min(opts.max_reads)];
    let ploidy = if copies == 1 { Ploidy::One } else { Ploidy::Two };
    let alleles = if copies == 0 { Vec::new() } else { call(obs, ploidy, &opts.genotype) };
    let anchor = reference.slice(locus.start - 1, locus.start).first().copied().unwrap_or(b'N');
    let ref_seq = reference.slice(locus.start, locus.end);
    (vcf::record(locus, anchor, ref_seq, &alleles, copies.max(1)), !alleles.is_empty())
}

/// One thread's readers and read buffer.
struct Worker {
    reads: Reads,
    fasta: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    buf: Vec<Read>,
}

impl Worker {
    fn open(
        reads: &std::path::Path,
        genome: &std::path::Path,
        cram_reference: Option<&SharedReference>,
        haplotags: Option<&std::sync::Arc<lrpipe_reads::haplotags::Haplotags>>,
    ) -> Result<Self, Error> {
        let mut reads = Reads::open(reads, cram_reference)?;
        if let Some(haplotags) = haplotags {
            reads.set_haplotags(std::sync::Arc::clone(haplotags));
        }
        // Bgzipped or not by its contents, not its name.
        let fasta = lrpipe_reads::indexed_fasta(genome).map_err(|e| Error::io(genome)(e.into_io_error()))?;
        Ok(Worker { reads, fasta, buf: Vec::new() })
    }

    /// The VCF records of a chunk of loci, and how many were called.
    fn chunk(&mut self, loci: &[Locus], opts: &RunOptions) -> Result<(String, usize), Error> {
        let contig = &loci[0].contig;
        let from = loci[0].start;
        let to = loci.iter().map(|l| l.end).max().unwrap_or(from);
        // Enough reference beyond each locus for refining a boundary, and for the farthest flank
        // fitted into a soft clip (three of its lengths out).
        let pad =
            (opts.extract.flank + (opts.extract.anchor + opts.extract.search) as i64).max(6 * opts.extract.flank) + 10;
        let win_start = (from - pad).max(0);
        let (Ok(first), Ok(last)) = (usize::try_from(win_start + 1), usize::try_from(to + pad)) else {
            return Err(Error::Input(format!("{contig}:{from}-{to}: bad interval")));
        };
        let (Some(first), Some(last)) = (Position::new(first), Position::new(last)) else {
            return Err(Error::Input(format!("{contig}:{from}-{to}: bad interval")));
        };
        let seq = match self.fasta.query(&Region::new(contig.as_str(), first..=last)) {
            Ok(rec) => rec.sequence().as_ref().to_ascii_uppercase(),
            Err(e) => return Err(Error::Input(format!("reference {contig}:{from}-{to}: {e}"))),
        };
        let reference = RefWindow { start: win_start, seq };
        self.buf.clear();
        self.reads.fetch(contig, from, to, opts.min_mapq, &mut self.buf)?;
        let reads = &self.buf;
        let longest = reads.iter().map(|r| r.end - r.start).max().unwrap_or(0);

        let mut text = String::new();
        let mut called = 0;
        for locus in loci {
            // Reads come sorted by start: those that can reach the locus start between
            // `start - longest` and its end (a read starting in the repeat has its left flank
            // in a soft clip).
            let hi = reads.partition_point(|r| r.start <= locus.end);
            let lo = reads[..hi].partition_point(|r| r.start < locus.start - longest);
            let obs: Vec<Observation> =
                reads[lo..hi].iter().filter_map(|r| observe(r, locus, &reference, opts)).collect();
            let obs = evenly(one_per_read(obs), opts.max_reads);
            let (record, was_called) = genotype_locus(locus, &obs, &reference, opts);
            called += usize::from(was_called);
            text.push_str(&record);
        }
        Ok((text, called))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locus(contig: &str, start: i64, end: i64) -> Locus {
        Locus {
            id: "t".into(),
            contig: contig.into(),
            start,
            end,
            motifs: vec![b"CAG".to_vec()],
            info: String::new(),
            shared_start: None,
            shared_end: None,
        }
    }

    fn genome(dir: &str, x_len: u64) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("lrpipe-tr-{dir}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let fa = dir.join("genome.fa");
        std::fs::write(fa.with_added_extension("fai"), format!("chr1\t1000\t6\t60\t61\nchrX\t{x_len}\t1030\t60\t61\n"))
            .expect("fai written");
        fa
    }

    #[test]
    fn at_most_max_reads_spread_evenly() {
        assert_eq!(evenly((0..10).collect::<Vec<u32>>(), 4), vec![0, 2, 5, 7]);
        assert_eq!(evenly((0..3).collect::<Vec<u32>>(), 4), vec![0, 1, 2]);
    }

    #[test]
    fn the_pars_are_the_references() {
        // LRP-010: GRCh37's PAR2 holds chrX:155,000,000, GRCh38's does not.
        let par2 = locus("chrX", 155_000_000, 155_000_030);
        let grch37 = x_pars(&genome("grch37", 155_270_560), Karyotype::Xy).expect("known build");
        assert_eq!(copies(&par2, Karyotype::Xy, &grch37), 2);
        let grch38 = x_pars(&genome("grch38", 156_040_895), Karyotype::Xy).expect("known build");
        assert_eq!(copies(&par2, Karyotype::Xy, &grch38), 1);
        assert_eq!(copies(&locus("chrX", 20_000, 20_030), Karyotype::Xy, &grch38), 2, "PAR1");
        assert_eq!(copies(&par2, Karyotype::Xx, &grch38), 2);
        // An unknown chrX (lrpipe#44): no PARs, so all of it is one copy in XY samples and two
        // in XX ones.
        let other = x_pars(&genome("other", 1_000_000), Karyotype::Xy).expect("an unknown chrX is no error");
        assert!(other.is_empty());
        assert_eq!(copies(&locus("chrX", 20_000, 20_030), Karyotype::Xy, &other), 1, "no PAR1");
        assert_eq!(copies(&locus("chrX", 20_000, 20_030), Karyotype::Xx, &other), 2);
    }
}

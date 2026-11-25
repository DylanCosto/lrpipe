//! Port of `variants.py` and `_variants.pyx`: detect which alleles of the
//! VCF variants each alignment supports.

use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};

use anyhow::{Result, anyhow, bail};

use crate::align::{edit_distance, edit_distance_affine_gap, enumerate_all_kmers, kmer_align};
use crate::alignment::{Aux, Record};
use crate::bam::{BamReader, FetchError};
use crate::compat::pystr::py_int;
use crate::core::genotype::Genotype;
use crate::core::read::{Read, ReadSet, ReadVariant};
use crate::vcf::VcfVariant;

pub enum ReadSetError {
    SampleNotFound,
    ReferenceNotFound(String),
    Other(anyhow::Error),
}

struct AlignedRead {
    read: Read,
    is_supplementary: bool,
    is_reverse: bool,
    reference_start: i64,
    reference_end: i64,
}

impl AlignedRead {
    fn distance(&self, other: &AlignedRead) -> i64 {
        0.max(self.reference_start - other.reference_end).max(other.reference_start - self.reference_end)
    }
}

#[derive(Clone, Default)]
struct AlleleProgress {
    progress: i64,
    length: i64,
    quality: i64,
    matched: i64,
    match_target: i64,
    inserted: i64,
    insert_target: i64,
    deleted: i64,
    delete_target: i64,
}

#[derive(Clone)]
struct VariantProgress {
    variant_id: usize,
    query_start: i64,
    alleles: Vec<AlleleProgress>,
}

impl VariantProgress {
    fn add_allele(&mut self, matches: i64, insertions: i64, deletions: i64) {
        self.alleles.push(AlleleProgress {
            length: matches + insertions + deletions,
            match_target: matches,
            insert_target: insertions,
            delete_target: deletions,
            ..Default::default()
        });
    }

    fn reset(&mut self, query_start: i64) {
        self.query_start = query_start;
        for a in self.alleles.iter_mut() {
            a.progress = 0;
            a.matched = 0;
            a.inserted = 0;
            a.deleted = 0;
            a.quality = 0;
        }
    }

    fn resolved(&self) -> Vec<usize> {
        self.alleles.iter().enumerate().filter(|(_, a)| a.progress == a.length).map(|(i, _)| i).collect()
    }

    fn n_pending(&self) -> usize {
        self.alleles.iter().filter(|a| 0 <= a.progress && a.progress < a.length).count()
    }
}

/// How alleles are detected by re-alignment against the reference
#[derive(Clone, Debug)]
pub enum Realignment {
    /// unit-cost edit distance (the default)
    EditDistance,
    /// `--affine-gap`: affine gap costs, `mismatch` for every mismatch
    AffineGap { gap_start: i64, gap_extend: i64, mismatch: f32 },
    /// `--use-kmerald`: alignment of k-mers with an error model
    Kmerald { costs: Arc<HashMap<(i64, i64), f32>>, kmer_size: u32, gap_penalty: f32, window: i64 },
}

pub struct ReadSetReader {
    reader: BamReader,
    mapq_threshold: i32,
    overhang: i64,
    realignment: Realignment,
    duplicates: bool,
    use_supplementary: bool,
    supplementary_distance_threshold: i64,
    allow_supplementary_only_read_groups: bool,
    /// Threads for allele detection (alignments are processed in order if None)
    detection_pool: Option<Arc<rayon::ThreadPool>>,
    progress: Option<Arc<ReadProgress>>,
}

/// How far reading a chromosome has got, updated while reading so that other
/// threads can estimate how much memory the whole chromosome will need.
#[derive(Debug, Default)]
pub struct ReadProgress {
    /// Approximate memory of the reads detected so far (summed over all reads
    /// with this progress)
    bytes: AtomicU64,
    /// `bytes` when the current read started (reads of earlier samples)
    bytes_at_start: AtomicU64,
    /// Variants that the alignments read so far have moved past, and all
    /// variants of the current read
    variants_done: AtomicU64,
    variants_total: AtomicU64,
    /// The reads of the current read have been detected and grouped (the
    /// memory peak of reading is over)
    loaded: AtomicBool,
}

impl ReadProgress {
    pub fn bytes(&self) -> u64 {
        self.bytes.load(Ordering::Relaxed)
    }

    /// `bytes` split into the reads of earlier reads and of the current one.
    pub fn bytes_split(&self) -> (u64, u64) {
        let before = self.bytes_at_start.load(Ordering::Relaxed);
        (before, self.bytes().saturating_sub(before))
    }

    /// Variants passed and total variants of the current read.
    pub fn variants(&self) -> (u64, u64) {
        (self.variants_done.load(Ordering::Relaxed), self.variants_total.load(Ordering::Relaxed))
    }

    /// Whether the reads of the current read are complete.
    pub fn loaded(&self) -> bool {
        self.loaded.load(Ordering::Relaxed)
    }

    pub(crate) fn start(&self, variants_total: usize) {
        self.bytes_at_start.store(self.bytes(), Ordering::Relaxed);
        self.variants_done.store(0, Ordering::Relaxed);
        self.variants_total.store(variants_total as u64, Ordering::Relaxed);
        self.loaded.store(false, Ordering::Relaxed);
    }

    #[cfg(test)]
    pub(crate) fn set_for_test(&self, variants_done: u64, bytes: u64) {
        self.variants_done.store(variants_done, Ordering::Relaxed);
        self.bytes.store(bytes, Ordering::Relaxed);
    }

    #[cfg(test)]
    pub(crate) fn set_loaded_for_test(&self) {
        self.loaded.store(true, Ordering::Relaxed);
    }

    /// Reads up to the variant with index `first` have been detected.
    fn passed(&self, first: usize) {
        self.variants_done.store(first as u64, Ordering::Relaxed);
    }

    fn add(&self, read: &AlignedRead) {
        let size = read.read.memory_size() + std::mem::size_of::<AlignedRead>() - std::mem::size_of::<Read>();
        self.bytes.fetch_add(size as u64, Ordering::Relaxed);
    }
}

/// Default sub-alignment id of primary alignments (`PRIMARY_DEFAULT_SUB_ALIGNMENT_ID`).
pub const PRIMARY_DEFAULT_SUB_ALIGNMENT_ID: &str = "____1";

/// `get_sub_alignment_id`. WhatsHap uses `str(hash((cigarstring, flag)))`
/// for non-primary alignments, which is randomized per Python process and
/// only compared within one run; any deterministic injective id works.
pub fn sub_alignment_id(rec: &Record) -> String {
    if is_primary(rec) {
        PRIMARY_DEFAULT_SUB_ALIGNMENT_ID.to_string()
    } else {
        format!("#{}#{}", rec.cigar(), rec.flags())
    }
}

/// Data about an alignment needed for allele detection (borrowed view).
struct Aln<'a> {
    rec: &'a Record,
    cigar: Vec<(u8, u32)>,
    seq: Vec<u8>,
}

impl ReadSetReader {
    pub fn new(
        paths: &[String],
        reference: Option<&str>,
        mapq_threshold: i32,
        use_supplementary: bool,
        supplementary_distance_threshold: i64,
        decompression_threads: usize,
    ) -> Result<Self> {
        Ok(ReadSetReader {
            reader: BamReader::open(paths, reference, decompression_threads)?,
            mapq_threshold,
            overhang: 10,
            realignment: Realignment::EditDistance,
            duplicates: false,
            use_supplementary,
            supplementary_distance_threshold,
            allow_supplementary_only_read_groups: false,
            detection_pool: None,
            progress: None,
        })
    }

    /// Report the progress of the following reads to `progress` (the result
    /// does not change).
    pub fn set_progress(&mut self, progress: Option<Arc<ReadProgress>>) {
        self.progress = progress;
    }

    /// Extend the re-alignment by this many bases to the left and right.
    pub fn set_overhang(&mut self, overhang: i64) {
        self.overhang = overhang;
    }

    pub fn set_realignment(&mut self, realignment: Realignment) {
        self.realignment = realignment;
    }

    /// Settings used by `haplotag`.
    pub fn set_haplotag_options(&mut self, duplicates: bool, allow_supplementary_only_read_groups: bool) {
        self.duplicates = duplicates;
        self.allow_supplementary_only_read_groups = allow_supplementary_only_read_groups;
    }

    /// Detect alleles in parallel on `pool` (the result does not change).
    pub fn set_detection_pool(&mut self, pool: Arc<rayon::ThreadPool>) {
        self.detection_pool = Some(pool);
    }

    pub fn n_paths(&self) -> usize {
        self.reader.n_paths()
    }

    pub fn has_reference(&self, chromosome: &str) -> bool {
        self.reader.has_reference(chromosome)
    }

    pub fn cram_files(&self) -> usize {
        self.reader.cram_files()
    }

    /// Decodes every CRAM file with `repository` (see
    /// [`crate::alignment::cram::CramReader::share_reference`]).
    pub fn share_reference(&mut self, repository: &pasta::fasta::Repository) {
        self.reader.share_reference(repository);
    }

    /// Detect alleles and return a ReadSet with the given variants.
    /// `reference` is the full chromosome sequence (uppercase) or None to
    /// detect alleles from the CIGAR only.
    pub fn read(
        &mut self,
        chromosome: &str,
        variants: &[VcfVariant],
        sample: Option<&str>,
        numeric_sample_id: i32,
        reference: Option<&[u8]>,
    ) -> Result<ReadSet, ReadSetError> {
        self.read_regions(chromosome, variants, sample, numeric_sample_id, reference, &[(0, None)])
    }

    /// `read` with `restricted_genotypes`: with a reference, only the alleles
    /// of `genotypes[i]` are considered for variant i.
    pub fn read_restricted(
        &mut self,
        chromosome: &str,
        variants: &[VcfVariant],
        sample: Option<&str>,
        numeric_sample_id: i32,
        reference: Option<&[u8]>,
        genotypes: &[Genotype],
    ) -> Result<ReadSet, ReadSetError> {
        assert_eq!(genotypes.len(), variants.len());
        self.read_impl(chromosome, variants, sample, numeric_sample_id, reference, &[(0, None)], Some(genotypes))
    }

    pub fn read_regions(
        &mut self,
        chromosome: &str,
        variants: &[VcfVariant],
        sample: Option<&str>,
        numeric_sample_id: i32,
        reference: Option<&[u8]>,
        regions: &[(i64, Option<i64>)],
    ) -> Result<ReadSet, ReadSetError> {
        self.read_impl(chromosome, variants, sample, numeric_sample_id, reference, regions, None)
    }

    #[allow(clippy::too_many_arguments)]
    fn read_impl(
        &mut self,
        chromosome: &str,
        variants: &[VcfVariant],
        sample: Option<&str>,
        numeric_sample_id: i32,
        reference: Option<&[u8]>,
        regions: &[(i64, Option<i64>)],
        restricted_genotypes: Option<&[Genotype]>,
    ) -> Result<ReadSet, ReadSetError> {
        let normalized: Vec<VcfVariant>;
        let variant_positions: Vec<i32> = variants.iter().map(|v| v.position).collect();
        let mut valid_positions: Vec<i32> = Vec::new();
        let mut var_progress: Vec<VariantProgress> = Vec::new();
        if reference.is_none() {
            normalized = variants.iter().map(|v| v.normalized()).collect();
            let valid_ids = detect_non_overlapping_variants(&normalized);
            valid_positions = valid_ids.iter().map(|&j| normalized[j].position).collect();
            var_progress = valid_ids.iter().map(|&j| build_var_progress(&normalized, j)).collect();
        } else {
            normalized = Vec::new();
        }
        let progress = self.progress.as_deref();
        if let Some(p) = progress {
            p.start(if reference.is_none() { valid_positions.len() } else { variant_positions.len() });
        }
        let detector = AlleleDetector {
            chromosome,
            variants,
            normalized: &normalized,
            var_progress: &var_progress,
            reference,
            restricted_genotypes,
            realignment: &self.realignment,
            overhang: self.overhang,
            numeric_sample_id,
            progress,
        };
        let mut i: usize = 0;
        let (mapq_threshold, use_supp, duplicates) = (self.mapq_threshold, self.use_supplementary, self.duplicates);
        // _usable_alignments; also returns the index of the first variant that
        // can overlap the alignment (WhatsHap advances this index over all
        // usable alignments in order, so it depends on the earlier ones)
        let mut usable = |rec: &Record| -> Option<usize> {
            if (!use_supp && rec.is_supplementary())
                || (rec.mapq() as i32) < mapq_threshold
                || rec.is_secondary()
                || rec.is_unmapped()
                || (!duplicates && rec.is_duplicate())
            {
                return None;
            }
            let positions = if reference.is_none() { &valid_positions[..] } else { &variant_positions[..] };
            while i < positions.len() && (positions[i] as i64) < rec.pos() {
                i += 1;
            }
            Some(i)
        };
        let mut aligned_reads: Vec<AlignedRead> = Vec::new();
        let res = match &self.detection_pool {
            None => self.reader.for_each_in_regions(chromosome, regions, sample, |source_id, rec| {
                if let Some(first) = usable(rec) {
                    if let Some(read) = detector.detect(source_id, rec, first)? {
                        aligned_reads.push(read);
                    }
                    if let Some(p) = progress {
                        p.passed(first);
                    }
                }
                Ok(())
            }),
            Some(pool) => {
                let reader = &mut self.reader;
                let mut result = Ok(());
                pool.in_place_scope(|scope| {
                    let mut batches = OrderedBatches::new(scope, &detector, pool.current_num_threads());
                    result = reader.for_each_in_regions(chromosome, regions, sample, |source_id, rec| {
                        if let Some(first) = usable(rec) {
                            batches.push(source_id, rec, first, &mut aligned_reads);
                        }
                        Ok(())
                    });
                    batches.finish(&mut aligned_reads);
                    if let Some(e) = batches.error {
                        // An alignment before the point where reading failed
                        // (if it did) failed first.
                        result = Err(FetchError::Other(e));
                    }
                });
                result
            }
        };
        let n_supplementary = aligned_reads.iter().filter(|r| r.is_supplementary).count();
        match res {
            Ok(()) => {}
            Err(FetchError::SampleNotFound) => return Err(ReadSetError::SampleNotFound),
            Err(FetchError::ReferenceNotFound(r)) => return Err(ReadSetError::ReferenceNotFound(r)),
            Err(FetchError::Other(e)) => return Err(ReadSetError::Other(e)),
        }
        log::info!("Number of supplementary alignments: {n_supplementary}");
        let readset = group_reads(
            aligned_reads,
            self.supplementary_distance_threshold,
            self.allow_supplementary_only_read_groups,
        )
        .map_err(ReadSetError::Other)?;
        if let Some(p) = &self.progress {
            p.loaded.store(true, Ordering::Relaxed);
        }
        Ok(readset)
    }
}

/// Everything needed to turn one usable alignment into a `Read`.
struct AlleleDetector<'a> {
    chromosome: &'a str,
    variants: &'a [VcfVariant],
    normalized: &'a [VcfVariant],
    var_progress: &'a [VariantProgress],
    reference: Option<&'a [u8]>,
    restricted_genotypes: Option<&'a [Genotype]>,
    realignment: &'a Realignment,
    overhang: i64,
    numeric_sample_id: i32,
    progress: Option<&'a ReadProgress>,
}

impl AlleleDetector<'_> {
    /// The read for `rec` if it covers any variant. `first` is the index of
    /// the first variant (of the valid variants without reference) to consider.
    fn detect(&self, source_id: usize, rec: &Record, first: usize) -> Result<Option<AlignedRead>> {
        let ps = match rec.aux(b"PS") {
            Ok(Aux::String(s)) => py_int(s)
                .and_then(|v| i32::try_from(v).ok())
                .ok_or_else(|| anyhow!("Invalid PS tag value ({s}) in read {}. PS must be an integer.", qname(rec)))?,
            Ok(a) => aux_int(&a).unwrap_or(-1),
            Err(_) => -1,
        };
        let hp = match rec.aux(b"HP") {
            Ok(a) => aux_int(&a).unwrap_or(-1),
            Err(_) => -1,
        };
        let bx = match rec.aux(b"BX") {
            Ok(Aux::String(s)) => s.to_string(),
            _ => String::new(),
        };
        let cigar: Vec<(u8, u32)> = rec.cigar().iter().collect();
        let reference_start = rec.pos();
        let reference_end = rec.cigar().end_pos();
        let mut read = Read::new(&qname(rec), rec.mapq() as i32, source_id as i32, self.numeric_sample_id);
        read.reference_start = reference_start as i32;
        read.bx_tag = bx;
        read.hp_tag = hp;
        read.ps_tag = ps;
        read.chromosome = self.chromosome.to_string();
        read.sub_alignment_id = sub_alignment_id(rec);
        read.is_supplementary = rec.is_supplementary();
        read.is_reverse = rec.is_reverse();
        read.reference_end = reference_end as i32;
        let aln = Aln { rec, cigar, seq: rec.seq() };
        let detected: Vec<(usize, usize, i32)> = match self.reference {
            None => detect_alleles_no_reference(self.normalized, self.var_progress, first, &aln)?,
            Some(refseq) => detect_alleles_by_alignment(
                self.variants,
                self.restricted_genotypes,
                first,
                &aln,
                refseq,
                self.overhang,
                self.realignment,
            )?,
        };
        for (j, allele, quality) in detected {
            read.add_variant(self.variants[j].position, allele as i32, quality);
        }
        if read.is_empty() {
            return Ok(None);
        }
        let read = AlignedRead {
            read,
            is_supplementary: rec.is_supplementary(),
            is_reverse: rec.is_reverse(),
            reference_start,
            reference_end,
        };
        if let Some(p) = self.progress {
            p.add(&read);
        }
        Ok(Some(read))
    }
}

/// Alignments per batch handed to the thread pool
const BATCH_SIZE: usize = 64;

type BatchResult = Result<Vec<AlignedRead>>;

/// Runs allele detection for batches of alignments on a rayon scope and
/// collects the reads in alignment order.
struct OrderedBatches<'scope, 'env> {
    scope: &'scope rayon::Scope<'env>,
    detector: &'env AlleleDetector<'env>,
    batch: Vec<(usize, Record, usize)>,
    sender: Sender<(usize, BatchResult)>,
    receiver: Receiver<(usize, BatchResult)>,
    /// Number of batches spawned and of batches whose reads were collected
    spawned: usize,
    collected: usize,
    max_in_flight: usize,
    finished: BTreeMap<usize, BatchResult>,
    /// For each batch not yet collected, the first variant index of its last
    /// alignment (for the progress)
    batch_firsts: VecDeque<usize>,
    /// The first error, in alignment order
    error: Option<anyhow::Error>,
}

impl<'scope, 'env> OrderedBatches<'scope, 'env> {
    fn new(scope: &'scope rayon::Scope<'env>, detector: &'env AlleleDetector<'env>, threads: usize) -> Self {
        let (sender, receiver) = channel();
        OrderedBatches {
            scope,
            detector,
            batch: Vec::with_capacity(BATCH_SIZE),
            sender,
            receiver,
            spawned: 0,
            collected: 0,
            max_in_flight: 4 * threads,
            finished: BTreeMap::new(),
            batch_firsts: VecDeque::new(),
            error: None,
        }
    }

    fn push(&mut self, source_id: usize, rec: &Record, first: usize, out: &mut Vec<AlignedRead>) {
        if self.error.is_some() {
            return;
        }
        self.batch.push((source_id, rec.clone(), first));
        if self.batch.len() == BATCH_SIZE {
            while self.spawned - self.collected >= self.max_in_flight {
                self.receive(out);
            }
            self.spawn();
        }
    }

    fn spawn(&mut self) {
        let batch = std::mem::replace(&mut self.batch, Vec::with_capacity(BATCH_SIZE));
        self.batch_firsts.push_back(batch.last().map_or(0, |&(_, _, first)| first));
        let (id, sender, detector) = (self.spawned, self.sender.clone(), self.detector);
        self.scope.spawn(move |_| {
            let mut reads = Vec::new();
            let mut result = Ok(());
            for (source_id, rec, first) in &batch {
                match detector.detect(*source_id, rec, *first) {
                    Ok(Some(read)) => reads.push(read),
                    Ok(None) => {}
                    Err(e) => {
                        result = Err(e);
                        break;
                    }
                }
            }
            // The receiver outlives the scope.
            let _ = sender.send((id, result.map(|()| reads)));
        });
        self.spawned += 1;
    }

    /// Wait for one batch, then collect all batches that are next in order.
    fn receive(&mut self, out: &mut Vec<AlignedRead>) {
        let (id, result) = self.receiver.recv().expect("batch sender alive");
        self.finished.insert(id, result);
        while let Some(result) = self.finished.remove(&self.collected) {
            self.collected += 1;
            let first = self.batch_firsts.pop_front().expect("one entry per spawned batch");
            if let Some(p) = self.detector.progress {
                p.passed(first);
            }
            match result {
                Ok(reads) if self.error.is_none() => out.extend(reads),
                Ok(_) => {}
                Err(e) => {
                    self.error.get_or_insert(e);
                }
            }
        }
    }

    fn finish(&mut self, out: &mut Vec<AlignedRead>) {
        if !self.batch.is_empty() && self.error.is_none() {
            self.spawn();
        }
        while self.collected < self.spawned {
            self.receive(out);
        }
    }
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn aux_int(a: &Aux) -> Option<i32> {
    Some(match *a {
        Aux::I8(v) => v as i32,
        Aux::U8(v) => v as i32,
        Aux::I16(v) => v as i32,
        Aux::U16(v) => v as i32,
        Aux::I32(v) => v,
        Aux::U32(v) => v as i32,
        _ => return None,
    })
}

pub fn is_primary(rec: &Record) -> bool {
    !(rec.is_supplementary() || rec.is_secondary() || rec.is_unmapped())
}

/// `ReadSetReader._group_reads` + `create_read_from_group` +
/// `_make_readset_from_grouped_reads`.
fn group_reads(
    reads: Vec<AlignedRead>,
    distance_threshold: i64,
    allow_supplementary_only_groups: bool,
) -> Result<ReadSet> {
    // dict keyed by (source_id, name, sub_alignment_id or None, sample_id), in
    // insertion order. Groups are lists of indices into `reads` (`first` and
    // `next`), so that no read or name is copied while most reads are in
    // memory twice.
    const END: u32 = u32::MAX;
    let mut first: Vec<u32> = Vec::new();
    let mut last: Vec<u32> = Vec::new();
    let mut next: Vec<u32> = vec![END; reads.len()];
    {
        type Key<'a> = (i32, &'a str, Option<&'a str>, i32);
        let mut groups: HashMap<Key, u32> = HashMap::with_capacity(reads.len());
        for (i, r) in reads.iter().enumerate() {
            let sub = allow_supplementary_only_groups.then_some(r.read.sub_alignment_id.as_str());
            let i = u32::try_from(i).expect("fewer than 2^32 alignments per chromosome");
            match groups.entry((r.read.source_id, r.read.name.as_str(), sub, r.read.sample_id)) {
                Entry::Occupied(e) => {
                    let g = *e.get() as usize;
                    next[last[g] as usize] = i;
                    last[g] = i;
                }
                Entry::Vacant(e) => {
                    e.insert(first.len() as u32);
                    first.push(i);
                    last.push(i);
                }
            }
        }
    }
    drop(last);
    let mut reads: Vec<Option<AlignedRead>> = reads.into_iter().map(Some).collect();
    let mut readset = ReadSet::with_capacity(first.len());
    let mut n_skipped = 0;
    let mut n_non_singleton = 0;
    let mut group = Vec::new();
    for head in first {
        let mut i = head;
        while i != END {
            group.push(reads[i as usize].take().expect("each read is in one group"));
            i = next[i as usize];
        }
        let read = if group.len() == 1 {
            create_read_from_single(group.pop().unwrap(), allow_supplementary_only_groups)?
        } else {
            n_non_singleton += 1;
            let read = create_read_from_group(&group, distance_threshold, allow_supplementary_only_groups)?;
            group.clear();
            read
        };
        match read {
            Some(read) => readset.add(read)?,
            None => n_skipped += 1,
        }
    }
    log::info!("Number of non-singleton groups: {n_non_singleton}");
    log::info!("Skipped {n_skipped} groups");
    Ok(readset)
}

/// `create_read_from_group` for a group of one alignment, reusing its read
/// (the result is the same as building a new read from it).
fn create_read_from_single(aligned: AlignedRead, allow_supplementary_only_groups: bool) -> Result<Option<Read>> {
    if aligned.is_supplementary && !allow_supplementary_only_groups {
        return Ok(None);
    }
    let mut read = aligned.read;
    let n_variants = read.len();
    // Variants at the same position (from VCF records at one position) keep
    // the first one unless their alleles differ, then neither is kept.
    let mut variants = std::mem::take(&mut read.variants);
    let mut seen: HashMap<i32, usize> = HashMap::new();
    let mut skip: HashSet<i32> = HashSet::new();
    let mut kept: Vec<ReadVariant> = Vec::with_capacity(variants.len());
    for v in variants.drain(..) {
        if let Some(&k) = seen.get(&v.position) {
            if kept[k].allele != v.allele {
                skip.insert(v.position);
            }
        } else {
            seen.insert(v.position, kept.len());
            kept.push(v);
        }
    }
    if !skip.is_empty() {
        kept.retain(|v| !skip.contains(&v.position));
    }
    read.variants = kept;
    read.mapqs.truncate(1);
    read.sort()?;
    if read.len() != n_variants {
        log::debug!("Converted read {} with {} variants to read with {} variants.", read.name, n_variants, read.len());
    }
    if allow_supplementary_only_groups {
        let sub = read.sub_alignment_id.clone();
        read.name.push_str(&sub);
    }
    Ok(Some(read))
}

fn create_read_from_group(
    group: &[AlignedRead],
    distance_threshold: i64,
    allow_supplementary_only_groups: bool,
) -> Result<Option<Read>> {
    if group.len() > 1 {
        log::debug!("Group of read '{}' has {} items.", group[0].read.name, group.len());
    }
    let mut primary: Option<&AlignedRead> = None;
    let mut n_primary = 0;
    for r in group {
        if !r.is_supplementary {
            n_primary += 1;
            primary = Some(r);
        }
    }
    if primary.is_none() && !allow_supplementary_only_groups {
        return Ok(None);
    }
    if n_primary > 2 {
        log::warn!("Read name '{}' has more than two primary alignments.", group[0].read.name);
        return Ok(None);
    }
    // With supplementary-only groups allowed, a group without primary has a
    // single read, which then serves as "primary" (the last read in Python).
    let primary = primary.unwrap_or(&group[group.len() - 1]);
    let mut reference_start = primary.reference_start;
    let mut variants: Vec<(i32, i32, u32)> = Vec::new();
    let mut variant_index: HashMap<i32, usize> = HashMap::new();
    let mut skip: HashSet<i32> = HashSet::new();
    for r in group {
        if r.is_supplementary {
            if r.is_reverse != primary.is_reverse {
                continue;
            }
            if primary.distance(r) > distance_threshold {
                continue;
            }
        }
        reference_start = reference_start.min(r.reference_start);
        for v in &r.read.variants {
            if let Some(&k) = variant_index.get(&v.position) {
                if variants[k].1 != v.allele {
                    skip.insert(v.position);
                }
            } else {
                variant_index.insert(v.position, variants.len());
                variants.push((v.position, v.allele, v.quality));
            }
        }
    }
    let p = &primary.read;
    let name =
        if allow_supplementary_only_groups { format!("{}{}", p.name, p.sub_alignment_id) } else { p.name.clone() };
    let mut union_read = Read::new(&name, p.mapqs[0], p.source_id, p.sample_id);
    union_read.reference_start = reference_start as i32;
    union_read.bx_tag = p.bx_tag.clone();
    union_read.hp_tag = p.hp_tag;
    union_read.ps_tag = p.ps_tag;
    union_read.chromosome = p.chromosome.clone();
    union_read.sub_alignment_id = p.sub_alignment_id.clone();
    union_read.is_supplementary = p.is_supplementary;
    union_read.is_reverse = primary.is_reverse;
    union_read.reference_end = primary.reference_end as i32;
    for (pos, allele, quality) in variants {
        if !skip.contains(&pos) {
            union_read.add_variant(pos, allele, quality as i32);
        }
    }
    union_read.sort()?;
    if union_read.len() != p.len() {
        log::debug!("Converted read {} with {} variants to read with {} variants.", p.name, p.len(), union_read.len());
    }
    Ok(Some(union_read))
}

/// `ReadSetReader.detect_non_overlapping_variants`
fn detect_non_overlapping_variants(variants: &[VcfVariant]) -> Vec<usize> {
    let mut j = 0;
    let mut conflicting = HashSet::new();
    let mut seen_pos = HashSet::new();
    while j < variants.len() {
        let v = &variants[j];
        if seen_pos.contains(&v.position) {
            conflicting.insert(j);
            j += 1;
            continue;
        }
        seen_pos.insert(v.position);
        let r = v.reference_allele.len() as i64;
        let max_del = v.alternative_alleles.iter().map(|a| r - a.len() as i64).max().unwrap();
        if max_del > 0 {
            let deletion_end = v.position as i64 + r;
            if j + 1 < variants.len() && (variants[j + 1].position as i64) < deletion_end {
                conflicting.insert(j);
                while j + 1 < variants.len() && (variants[j + 1].position as i64) < deletion_end {
                    j += 1;
                    conflicting.insert(j);
                }
            }
        }
        j += 1;
    }
    (0..variants.len()).filter(|j| !conflicting.contains(j)).collect()
}

fn build_var_progress(variants: &[VcfVariant], j: usize) -> VariantProgress {
    let mut v = VariantProgress { variant_id: j, query_start: 0, alleles: Vec::new() };
    let ref_len = variants[j].reference_allele.len() as i64;
    v.add_allele(ref_len, 0, 0);
    for alt in &variants[j].alternative_alleles {
        let alt_len = alt.len() as i64;
        v.add_allele(ref_len.min(alt_len), 0.max(alt_len - ref_len), 0.max(ref_len - alt_len));
    }
    v
}

/// `_iterate_cigar`: yields (variant index, cigar index, consumed, query_pos)
fn iterate_cigar(
    variants: &[VcfVariant],
    j: usize,
    reference_start: i64,
    cigar: &[(u8, u32)],
) -> Result<Vec<(usize, usize, i64, i64)>> {
    let mut out = Vec::new();
    let mut ref_pos = reference_start;
    let mut query_pos: i64 = 0;
    let n = variants.len();
    let mut j = j;
    while j < n && (variants[j].position as i64) < ref_pos {
        j += 1;
    }
    let mut v_position: i64 = 0;
    for (i, &(op, length)) in cigar.iter().enumerate() {
        let length = length as i64;
        if j < n {
            v_position = variants[j].position as i64;
        }
        match op {
            0 | 7 | 8 => {
                while j < n && v_position < ref_pos + length {
                    assert!(v_position >= ref_pos);
                    out.push((j, i, v_position - ref_pos, query_pos + v_position - ref_pos));
                    j += 1;
                    if j < n {
                        v_position = variants[j].position as i64;
                    }
                }
                query_pos += length;
                ref_pos += length;
            }
            1 => {
                if j < n && v_position == ref_pos {
                    out.push((j, i, 0, query_pos));
                    j += 1;
                    if j < n {
                        v_position = variants[j].position as i64;
                    }
                }
                query_pos += length;
            }
            2 => {
                while j < n && v_position < ref_pos + length {
                    assert!(v_position >= ref_pos);
                    out.push((j, i, v_position - ref_pos, query_pos));
                    j += 1;
                    if j < n {
                        v_position = variants[j].position as i64;
                    }
                }
                ref_pos += length;
            }
            3 => {
                while j < n && v_position < ref_pos + length {
                    assert!(v_position >= ref_pos);
                    j += 1;
                    if j < n {
                        v_position = variants[j].position as i64;
                    }
                }
                ref_pos += length;
            }
            4 => query_pos += length,
            5 | 6 => {}
            _ => bail!("Unsupported CIGAR operation: {op}"),
        }
    }
    Ok(out)
}

/// `split_cigar_left`: cigar ops from the split point towards the start.
fn split_cigar_left(cigar: &[(u8, u32)], i: usize, consumed: i64) -> Vec<(u8, i64)> {
    let (middle_op, middle_length) = cigar[i];
    assert!(consumed <= middle_length as i64);
    let mut out = Vec::new();
    if consumed > 0 {
        out.push((middle_op, consumed));
    }
    for j in (0..i).rev() {
        out.push((cigar[j].0, cigar[j].1 as i64));
    }
    out
}

/// `split_cigar_right`: cigar ops from the split point towards the end.
fn split_cigar_right(cigar: &[(u8, u32)], i: usize, consumed: i64) -> Vec<(u8, i64)> {
    let (middle_op, middle_length) = cigar[i];
    let mut out = Vec::new();
    if consumed < middle_length as i64 {
        out.push((middle_op, middle_length as i64 - consumed));
    }
    for &(op, l) in &cigar[i + 1..] {
        out.push((op, l as i64));
    }
    out
}

/// `cigar_prefix_length`
fn cigar_prefix_length(cigar: &[(u8, i64)], reference_bases: i64) -> Result<(i64, i64)> {
    let mut ref_pos = 0;
    let mut query_pos = 0;
    for &(op, length) in cigar {
        match op {
            0 | 7 | 8 => {
                ref_pos += length;
                query_pos += length;
                if ref_pos >= reference_bases {
                    return Ok((reference_bases, query_pos + reference_bases - ref_pos));
                }
            }
            2 => {
                ref_pos += length;
                if ref_pos >= reference_bases {
                    return Ok((reference_bases, query_pos));
                }
            }
            1 => query_pos += length,
            4 | 5 => {}
            3 => return Ok((reference_bases, query_pos)),
            _ => bail!("unknown CIGAR operator"),
        }
    }
    assert!(ref_pos < reference_bases);
    Ok((ref_pos, query_pos))
}

/// Python-style clamped slice `s[a:b]` (non-negative indices).
fn slice(s: &[u8], a: i64, b: i64) -> &[u8] {
    let len = s.len() as i64;
    let a = a.clamp(0, len) as usize;
    let b = b.clamp(0, len) as usize;
    if a >= b { &[] } else { &s[a..b] }
}

/// `ReadSetReader.realign` (default settings: unit-cost edit distance).
/// With `restricted`, only the alleles of that genotype are candidates.
#[allow(clippy::too_many_arguments)]
fn realign(
    variant: &VcfVariant,
    restricted: Option<&Genotype>,
    aln: &Aln,
    i: usize,
    consumed: i64,
    query_pos: i64,
    reference: &[u8],
    overhang: i64,
    realignment: &Realignment,
) -> Result<Option<(usize, i32)>> {
    if variant.alternative_alleles.iter().any(|a| a.starts_with('<')) {
        return Ok(None);
    }
    let left = split_cigar_left(&aln.cigar, i, consumed);
    let right = split_cigar_right(&aln.cigar, i, consumed);
    if let Realignment::Kmerald { costs, kmer_size, gap_penalty, window } = realignment {
        let (left_ref_bases, left_query_bases) = cigar_prefix_length(&left, *window)?;
        let (right_ref_bases, right_query_bases) =
            cigar_prefix_length(&right, variant.reference_allele.len() as i64 + window)?;
        let pos = variant.position as i64;
        assert!(pos - left_ref_bases >= 0);
        assert!(pos + right_ref_bases <= reference.len() as i64);
        if aln.seq.is_empty() {
            bail!("'NoneType' object is not subscriptable (read {} has no sequence)", qname(aln.rec));
        }
        let [alternative] = &variant.alternative_alleles[..] else {
            bail!("'MultiallelicVcfVariant' object has no attribute 'alternative_allele'");
        };
        let query = enumerate_all_kmers(
            slice(&aln.seq, query_pos - left_query_bases, query_pos + right_query_bases),
            *kmer_size,
        );
        let ref_kmers = enumerate_all_kmers(slice(reference, pos - left_ref_bases, pos + right_ref_bases), *kmer_size);
        let mut alt_seq = slice(reference, pos - left_ref_bases, pos).to_vec();
        alt_seq.extend_from_slice(alternative.as_bytes());
        alt_seq.extend_from_slice(slice(reference, pos + variant.reference_allele.len() as i64, pos + right_ref_bases));
        let alt_kmers = enumerate_all_kmers(&alt_seq, *kmer_size);
        let distance_ref = kmer_align(&ref_kmers, &query, costs, *gap_penalty);
        let distance_alt = kmer_align(&alt_kmers, &query, costs, *gap_penalty);
        return Ok(if distance_ref < distance_alt {
            Some((0, 30))
        } else if distance_ref > distance_alt {
            Some((1, 30))
        } else {
            None
        });
    }
    let (left_ref_bases, left_query_bases) = cigar_prefix_length(&left, overhang)?;
    let (right_ref_bases, right_query_bases) =
        cigar_prefix_length(&right, variant.reference_allele.len() as i64 + overhang)?;
    let pos = variant.position as i64;
    assert!(pos - left_ref_bases >= 0);
    assert!(pos + right_ref_bases <= reference.len() as i64);
    if aln.seq.is_empty() {
        bail!("'NoneType' object is not subscriptable (read {} has no sequence)", qname(aln.rec));
    }
    let query = slice(&aln.seq, query_pos - left_query_bases, query_pos + right_query_bases);
    let left_pad = slice(reference, pos - left_ref_bases, pos);
    let right_pad = slice(reference, pos + variant.reference_allele.len() as i64, pos + right_ref_bases);
    let mut padded: Vec<Vec<u8>> = vec![slice(reference, pos - left_ref_bases, pos + right_ref_bases).to_vec()];
    for alt in &variant.alternative_alleles {
        let mut a = left_pad.to_vec();
        a.extend_from_slice(alt.as_bytes());
        a.extend_from_slice(right_pad);
        padded.push(a);
    }
    let allowed = restricted.map(|g| g.as_vector());
    let distance = |a: &[u8]| match realignment {
        Realignment::AffineGap { gap_start, gap_extend, mismatch } => {
            edit_distance_affine_gap(query, a, *mismatch, *gap_start, *gap_extend)
        }
        _ => edit_distance(query, a),
    };
    let mut distances: Vec<(usize, i64)> = padded
        .iter()
        .enumerate()
        .filter(|(i, _)| allowed.as_ref().is_none_or(|a| a.contains(&(*i as u32))))
        .map(|(i, a)| (i, distance(a)))
        .collect();
    if distances.is_empty() {
        // a missing genotype leaves no candidate
        bail!("IndexError: list index out of range");
    }
    distances.sort_by_key(|d| d.1); // stable, like Python's sort
    let quality = match realignment {
        // (at most zero: WhatsHap subtracts the larger distance)
        Realignment::AffineGap { .. } => {
            (if distances.len() > 1 { distances[0].1 - distances[1].1 } else { distances[0].1 }) as i32
        }
        _ => 30,
    };
    if distances.len() == 1 || distances[0].1 < distances[1].1 { Ok(Some((distances[0].0, quality))) } else { Ok(None) }
}

fn detect_alleles_by_alignment(
    variants: &[VcfVariant],
    restricted_genotypes: Option<&[Genotype]>,
    j: usize,
    aln: &Aln,
    reference: &[u8],
    overhang: i64,
    realignment: &Realignment,
) -> Result<Vec<(usize, usize, i32)>> {
    let mut out = Vec::new();
    if aln.cigar.is_empty() {
        return Ok(out);
    }
    for (index, i, consumed, query_pos) in iterate_cigar(variants, j, aln.rec.pos(), &aln.cigar)? {
        let restricted = restricted_genotypes.map(|g| &g[index]);
        if let Some((allele, quality)) =
            realign(&variants[index], restricted, aln, i, consumed, query_pos, reference, overhang, realignment)?
            && allele <= variants[index].alternative_alleles.len()
        {
            out.push((index, allele, quality));
        }
    }
    Ok(out)
}

/// `_detect_alleles` (allele detection without reference, from the CIGAR)
/// Entries of `var_progress` serve as templates: every entry is reset before
/// use, so working on fresh copies is equivalent to WhatsHap's reuse of one
/// table and lets reads be processed independently.
fn detect_alleles_no_reference(
    variants: &[VcfVariant],
    var_progress: &[VariantProgress],
    first: usize,
    aln: &Aln,
) -> Result<Vec<(usize, usize, i32)>> {
    let rec = aln.rec;
    let qual = rec.qual();
    let has_qualities = !qual.is_empty() && qual[0] != 255;
    let mut out = Vec::new();
    let mut ref_pos = rec.pos();
    let mut query_pos: i64 = 0;
    let n = var_progress.len();
    let mut j = first;
    while j < n {
        let var_pos = variants[var_progress[j].variant_id].position as i64;
        if var_pos >= ref_pos {
            break;
        }
        j += 1;
    }
    let mut vqueue: VecDeque<VariantProgress> = VecDeque::new();
    let yield_resolved = |vp: &VariantProgress, out: &mut Vec<(usize, usize, i32)>| {
        let resolved = vp.resolved();
        let lengths: Vec<i64> = resolved.iter().map(|&r| vp.alleles[r].length).collect();
        let max = *lengths.iter().max().unwrap();
        let i = resolved[lengths.iter().position(|&l| l == max).unwrap()];
        let a = &vp.alleles[i];
        let q = if a.length > 0 { a.quality.div_euclid(a.length) } else { 30 };
        out.push((vp.variant_id, i, q as i32));
    };
    for &(op, length) in &aln.cigar {
        let length = length as i64;
        while j < n {
            let var_pos = variants[var_progress[j].variant_id].position as i64;
            if var_pos >= ref_pos {
                break;
            }
            j += 1;
        }
        match op {
            3 => {
                ref_pos += length;
                continue;
            }
            4 => {
                query_pos += length;
                continue;
            }
            5 | 6 => continue,
            _ => {}
        }
        let ref_end = ref_pos + length;
        while j < n {
            let var_id = var_progress[j].variant_id;
            let var_pos = variants[var_id].position as i64;
            if var_pos >= ref_end {
                break;
            }
            let ref_len = variants[var_id].reference_allele.len();
            if op == 1 && ref_len > 0 {
                break;
            }
            if op == 2 && ref_len == 0 {
                j += 1;
                continue;
            }
            let query_start = if op != 2 { query_pos + var_pos - ref_pos } else { query_pos };
            let mut entry = var_progress[j].clone();
            entry.reset(query_start);
            vqueue.push_back(entry);
            j += 1;
        }
        let (mut ref_end, mut query_end) = (ref_pos, query_pos);
        match op {
            0 | 7 | 8 => {
                ref_end += length;
                query_end += length;
            }
            1 => query_end += length,
            2 => ref_end += length,
            _ => bail!("Unsupported CIGAR operation: {op}"),
        }
        for entry in vqueue.iter_mut() {
            let variant = &variants[entry.variant_id];
            match op {
                0 | 7 | 8 => detect_match(variant, entry, aln, qual, has_qualities, query_pos, length)?,
                1 => detect_insertion(variant, entry, aln, length)?,
                2 => detect_deletion(entry, length),
                _ => unreachable!(),
            }
        }
        ref_pos = ref_end;
        query_pos = query_end;
        while let Some(vp) = vqueue.pop_front() {
            let num_resolved = vp.resolved().len();
            let num_pending = vp.n_pending();
            if num_resolved >= 1 && num_pending == 0 {
                yield_resolved(&vp, &mut out);
            } else if num_pending > 0 {
                vqueue.push_front(vp);
                break;
            }
        }
    }
    for vp in &vqueue {
        if !vp.resolved().is_empty() && vp.n_pending() == 0 {
            yield_resolved(vp, &mut out);
        }
    }
    Ok(out)
}

fn seq_at(aln: &Aln, i: i64) -> Result<u8> {
    if i < 0 || i as usize >= aln.seq.len() {
        bail!("string index out of range (read {})", qname(aln.rec));
    }
    Ok(aln.seq[i as usize])
}

fn detect_match(
    variant: &VcfVariant,
    entry: &mut VariantProgress,
    aln: &Aln,
    qual: &[u8],
    has_qualities: bool,
    query_pos: i64,
    length: i64,
) -> Result<()> {
    let query_start = entry.query_start;
    let op_start = 0.max(entry.query_start - query_pos);
    for (i, a) in entry.alleles.iter_mut().enumerate() {
        if a.progress < 0 {
            continue;
        }
        let mut ops_consumed = op_start;
        let allele_seq = variant.get_allele(i).as_bytes();
        // Faithful to WhatsHap: the query index is computed once before the
        // loop and not advanced within it.
        let query_index = query_start + a.matched + a.inserted;
        while a.matched < a.match_target && ops_consumed < length {
            let qbase = seq_at(aln, query_index)?;
            let vbase = allele_seq[(a.matched + a.inserted) as usize];
            if qbase == vbase {
                ops_consumed += 1;
                if has_qualities {
                    a.quality += qual[query_index as usize] as i64;
                } else {
                    a.quality += 30;
                }
                a.matched += 1;
                a.progress += 1;
            } else {
                break;
            }
        }
        if ops_consumed < length && a.progress < a.length {
            a.progress = -1;
        }
    }
    Ok(())
}

fn detect_insertion(variant: &VcfVariant, entry: &mut VariantProgress, aln: &Aln, length: i64) -> Result<()> {
    let query_start = entry.query_start;
    for (i, a) in entry.alleles.iter_mut().enumerate() {
        if a.progress < 0 {
            continue;
        }
        let mut ops_consumed = 0;
        let allele_seq = variant.get_allele(i).as_bytes();
        while a.inserted < a.insert_target && ops_consumed < length {
            ops_consumed += 1;
            let qbase = seq_at(aln, query_start + a.matched + a.inserted)?;
            let vbase = allele_seq[(a.matched + a.inserted) as usize];
            if qbase == vbase {
                a.inserted += 1;
                a.progress += 1;
                a.quality += 30;
            } else {
                break;
            }
        }
        if ops_consumed < length && 0 < a.progress && a.progress < a.length {
            a.progress = -1;
        }
    }
    Ok(())
}

fn detect_deletion(entry: &mut VariantProgress, length: i64) {
    for a in entry.alleles.iter_mut() {
        if a.progress < 0 {
            continue;
        }
        let mut ops_consumed = 0;
        while a.deleted < a.delete_target && ops_consumed < length {
            ops_consumed += 1;
            a.deleted += 1;
            a.progress += 1;
            a.quality += 30;
        }
        if ops_consumed < length && a.progress < a.length {
            a.progress = -1;
        }
    }
}

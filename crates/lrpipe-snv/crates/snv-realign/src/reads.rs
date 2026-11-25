//! Port of `preprocess/RealignReads.py`: realigns the reads of a full-alignment job around its
//! candidate sites and writes them as SAM lines, in the order the original pipes them
//! into `samtools mpileup -`.
//!
//! Reads are taken in 5 kb chunks. Positions where at least two good reads disagree with the
//! reference become candidates; candidates closer than 160 bp form windows, and each window's
//! reads are realigned to the haplotypes of a de Bruijn graph. A read is written once the
//! chunks have moved more than 1,020 bp past its (realigned) start.

use std::collections::{BTreeSet, HashMap, VecDeque};

use snv_mpileup::{AlignmentFile, Record, Regions, View};

use crate::Error;
use crate::dbg::candidate_haplotypes;
use crate::realigner::{WindowRead, realign_reads};
use rayon::prelude::*;

/// `min_dbg_mapping_quality` and `min_dbg_base_quality`.
const MIN_DBG_QUALITY: i64 = 20;
/// `region_expansion_in_bp` (also `expand_align_ref_region`).
const EXPANSION: i64 = 20;
/// `min_windows_distance`.
const MIN_WINDOW_DISTANCE: i64 = EXPANSION * 4;
/// `max_window_size` (also `max_region_reads_num`).
const MAX_WINDOW: i64 = 1000;
/// `param.realign_chunk_size`.
const CHUNK: i64 = 5000;
/// `graph_min_mapping_quality`.
const GRAPH_MIN_MAPQ: i64 = 14;
/// Base qualities below this are hidden from the graph.
const GRAPH_MIN_BASE_QUALITY: i64 = 15;
/// `--minCoverage` (CallVarBam does not pass it to RealignReads).
const MIN_COVERAGE: f64 = 2.0;

/// `CHAR_STR`: the 87 digits of the five-character read names.
const NAME_DIGITS: &[u8] = b"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$%&()*+./:;<=>?[]^`{|}~";
/// `EXP`.
const NAME_LEN: u32 = 5;

/// `simplfy_read_name(rs_idx)` after incrementing it.
fn read_name(index: u64) -> String {
    let base = NAME_DIGITS.len() as u64;
    let mut out = String::with_capacity(NAME_LEN as usize);
    let mut rest = index;
    for exp in (1..NAME_LEN).rev() {
        let div = base.pow(exp);
        out.push(char::from(NAME_DIGITS[(rest / div) as usize]));
        rest %= div;
    }
    out.push(char::from(NAME_DIGITS[(rest % base) as usize]));
    out
}

/// `seq[i]` with Python's negative indices; `None` where Python raises `IndexError`.
fn py_index(seq: &[u8], i: i64) -> Option<u8> {
    let len = seq.len() as i64;
    let i = if i < 0 { i + len } else { i };
    (0..len).contains(&i).then(|| seq[i as usize])
}

/// `seq[start:end]` with Python's slice bounds.
fn py_slice(seq: &[u8], start: i64, end: i64) -> &[u8] {
    let len = seq.len() as i64;
    let bound = |i: i64| if i < 0 { (len + i).max(0) } else { i.min(len) } as usize;
    let (start, end) = (bound(start), bound(end));
    seq.get(start..end).unwrap_or_default()
}

/// `(length, op)` of each `(\d+)([MIDNSHP=X])` in a CIGAR string.
fn cigar_ops(cigar: &str) -> Vec<(i64, u8)> {
    let mut ops = Vec::new();
    let mut advance = 0i64;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            advance = advance * 10 + i64::from(c - b'0');
        } else {
            ops.push((advance, c));
            advance = 0;
        }
    }
    ops
}

/// A read as `RealignReads.py` keeps it (`Read`).
#[derive(Debug, Clone)]
struct Read {
    start: i64,
    end: i64,
    seq: Vec<u8>,
    cigar: String,
    mapq: i64,
    quals: Vec<i64>,
    raw_quals: String,
    name: String,
    flag: String,
    rnext: String,
    pnext: String,
    phasing: Option<char>,
    best_cigar: String,
    best_pos: i64,
    /// `best_align_score`; Python treats 0 like "not set".
    best_score: Option<i64>,
}

impl Read {
    /// `count_align_score`: inserted and deleted bases.
    fn align_score(cigar: &str) -> i64 {
        cigar_ops(cigar).iter().filter(|(_, op)| *op == b'I' || *op == b'D').map(|(l, _)| l).sum()
    }

    /// `set_realignment_info`.
    fn set_realignment(&mut self, cigar: &str, start: i64) {
        let cigar = cigar.replace('X', "M");
        if cigar == self.cigar && start == self.start {
            return;
        }
        let best = self.best_score.filter(|&s| s != 0);
        if best.is_some() && cigar == self.best_cigar && start == self.best_pos {
            return;
        }
        let score = Self::align_score(&cigar);
        if best.is_none_or(|b| score >= b) {
            self.best_cigar = cigar;
            self.best_pos = start;
            self.best_score = Some(score);
        }
    }

    /// The SAM line `RealignReads.py` writes: the TLEN column repeats PNEXT (`self.TLEN =
    /// PNEXT`), and the last column is `HP:i:N` or empty.
    fn sam_line(&self, contig: &str) -> String {
        let phasing = self.phasing.map_or_else(String::new, |hp| format!("HP:i:{hp}"));
        format!(
            "{}\t{}\t{contig}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{phasing}",
            self.name,
            self.flag,
            self.best_pos + 1,
            self.mapq,
            self.best_cigar,
            self.rnext,
            self.pnext,
            self.pnext,
            String::from_utf8_lossy(&self.seq),
            self.raw_quals,
        )
    }
}

/// `is_too_many_soft_clipped_bases_for_a_read_from`: less than 55% aligned.
fn mostly_soft_clipped(cigar: &str) -> bool {
    let (mut soft, mut total) = (0i64, 0i64);
    for (len, op) in cigar_ops(cigar) {
        if op == b'S' {
            soft += len;
        }
        total += len;
    }
    1.0 - (soft as f64) / ((total + 1) as f64) < 0.55
}

/// `get_halpotype_tag`: the first digit of the first tag containing `HP:i:`.
fn haplotype_tag(tags: &[&str]) -> Option<char> {
    let tag = tags.iter().find(|t| t.contains("HP:i:"))?;
    tag.chars().nth(5).filter(char::is_ascii_digit)
}

/// The realigned reads of one job, produced chunk by chunk as SAM lines.
pub struct RealignedReads<'a> {
    file: &'a AlignmentFile,
    view: View<'a>,
    contig: String,
    reference: Vec<u8>,
    /// 0-based position of `reference[0]`.
    reference_offset: i64,
    chunk: Option<(i64, i64)>,
    /// `aligned_reads`, in insertion order.
    reads: Vec<Read>,
    pileup: Vec<(i64, i64)>,
    pileup_index: HashMap<i64, usize>,
    next_id: u64,
    out: VecDeque<String>,
    finished: bool,
}

impl std::fmt::Debug for RealignedReads<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RealignedReads").field("contig", &self.contig).finish_non_exhaustive()
    }
}

/// Where a job's reads and reference come from.
#[derive(Debug, Clone, Copy)]
pub struct RealignInput<'a> {
    /// The contig.
    pub contig: &'a str,
    /// `(ctg_start, ctg_end)` from the candidate BED: the smallest start + 1 and the largest
    /// end + 1.
    pub span: (i64, i64),
}

impl<'a> RealignedReads<'a> {
    /// `reads_realignment` for a `--full_aln_regions` job: `samtools view -h BAM
    /// CTG:START-END -q 5 -L BED` over the BED's span widened by 1,000 bp, and the reference
    /// widened by 5,000 bp. `view_file` and `file` are the same BAM, opened twice.
    ///
    /// # Errors
    ///
    /// Fails if the inputs cannot be read.
    pub fn new(
        view_file: &'a mut AlignmentFile,
        file: &'a AlignmentFile,
        fasta: &snv_pileup::Fasta,
        input: RealignInput<'_>,
        regions: Regions,
    ) -> Result<Self, Error> {
        let (ctg_start, ctg_end) = input.span;
        let extend_start = (ctg_start - MAX_WINDOW).max(1);
        let extend_end = ctg_end + MAX_WINDOW;
        let reference_start = (ctg_start - 5 * 1000).max(1);
        let reference_end = ctg_end + 5 * 1000;
        let to_i32 = |v: i64| i32::try_from(v).map_err(|_| Error::Exit(format!("position {v} out of range")));
        let mut reference = fasta.fetch(input.contig, to_i32(reference_start - 1)?, to_i32(reference_end - 1)?)?;
        reference.make_ascii_uppercase();
        if reference.is_empty() {
            return Err(Error::Exit("[ERROR] Failed to load reference sequence".into()));
        }
        let region = format!("{}:{extend_start}-{extend_end}", input.contig);
        let view = view_file.view(&region, 5, Some(regions))?;
        Ok(Self {
            file,
            view,
            contig: input.contig.to_owned(),
            reference,
            reference_offset: reference_start - 1,
            chunk: None,
            reads: Vec::new(),
            pileup: Vec::new(),
            pileup_index: HashMap::new(),
            next_id: 0,
            out: VecDeque::new(),
            finished: false,
        })
    }

    fn reference_at(&self, pos: i64) -> Result<u8, Error> {
        py_index(&self.reference, pos - self.reference_offset)
            .ok_or_else(|| Error::Exit(format!("IndexError: no reference base at {}:{pos}", self.contig)))
    }

    fn count(&mut self, pos: i64) {
        match self.pileup_index.get(&pos) {
            Some(&i) => self.pileup[i].1 += 1,
            None => {
                self.pileup_index.insert(pos, self.pileup.len());
                self.pileup.push((pos, 1));
            }
        }
    }

    /// The body of `samtools_view_generator_from` for one record.
    fn add_record(&mut self, record: &Record) -> Result<(), Error> {
        let line = self.file.format(record)?;
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns.len() < 11 {
            return Err(Error::Exit(format!("short SAM line {line:?}")));
        }
        if columns[2] != self.contig {
            return Ok(());
        }
        let parse = |s: &str| s.parse::<i64>().map_err(|_| Error::Exit(format!("bad SAM field {s:?}")));
        let flag = parse(columns[1])?;
        let pos = parse(columns[3])? - 1;
        let mapq = parse(columns[4])?;
        let cigar = columns[5];
        let seq = columns[9].to_ascii_uppercase().into_bytes();
        let raw_quals = columns[10];
        let quals: Vec<i64> = raw_quals.bytes().map(|q| i64::from(q) - 33).collect();
        let phasing = haplotype_tag(&columns[11..]);
        self.next_id += 1;
        let name = read_name(self.next_id - 1);
        let (_, chunk_end) = *self.chunk.get_or_insert((pos, pos + CHUNK));
        if pos >= chunk_end + EXPANSION {
            self.process_chunk()?;
            if let Some((start, end)) = self.chunk.as_mut() {
                *start += CHUNK;
                *end += CHUNK;
            }
        }
        let (chunk_start, chunk_end) = self.chunk.expect("set above");

        let deleted: i64 = cigar_ops(cigar).iter().filter(|(_, op)| *op == b'D').map(|(l, _)| l).sum();
        let read = Read {
            start: pos,
            end: pos + seq.len() as i64 + deleted,
            seq,
            cigar: cigar.to_owned(),
            mapq,
            quals,
            raw_quals: raw_quals.to_owned(),
            name,
            flag: flag.to_string(),
            rnext: columns[6].to_owned(),
            pnext: columns[7].to_owned(),
            phasing,
            best_cigar: cigar.to_owned(),
            best_pos: pos,
            best_score: None,
        };
        if cigar == "*" || mostly_soft_clipped(cigar) {
            return Ok(());
        }
        if mapq >= MIN_DBG_QUALITY {
            self.count_differences(&read, chunk_start, chunk_end)?;
        }
        self.reads.push(read);
        Ok(())
    }

    /// The pileup of positions where the read disagrees with the reference.
    fn count_differences(&mut self, read: &Read, chunk_start: i64, chunk_end: i64) -> Result<(), Error> {
        let quality = |i: i64| {
            usize::try_from(i)
                .ok()
                .and_then(|i| read.quals.get(i))
                .copied()
                .ok_or_else(|| Error::Exit("IndexError: list index out of range".into()))
        };
        let (mut ref_pos, mut query_pos) = (read.start, 0i64);
        let out_of_region = |p: i64| p < chunk_start - EXPANSION || p > chunk_end + EXPANSION;
        for (advance, op) in cigar_ops(&read.cigar) {
            match op {
                b'=' => {
                    ref_pos += advance;
                    query_pos += advance;
                }
                b'M' | b'X' => {
                    for _ in 0..advance {
                        if quality(query_pos)? >= MIN_DBG_QUALITY {
                            let reference_base = self.reference_at(ref_pos)?;
                            let query_base = usize::try_from(query_pos)
                                .ok()
                                .and_then(|i| read.seq.get(i))
                                .copied()
                                .ok_or_else(|| Error::Exit("IndexError: string index out of range".into()))?;
                            if b"ACGT".contains(&reference_base) && query_base != reference_base {
                                self.count(ref_pos);
                            }
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                b'I' | b'S' => {
                    let pre_base = self.reference_at(ref_pos - 1)?;
                    let start = usize::try_from(query_pos).unwrap_or(0).min(read.quals.len());
                    let end = (start + advance as usize).min(read.quals.len());
                    let good = read.quals[start..end].iter().all(|&q| q >= MIN_DBG_QUALITY);
                    if !out_of_region(ref_pos) && b"ACGT".contains(&pre_base) && good {
                        for p in ref_pos - advance..ref_pos + advance {
                            self.count(p);
                        }
                    }
                    query_pos += advance;
                }
                b'D' => {
                    let pre_base = self.reference_at(ref_pos - 1)?;
                    if !out_of_region(ref_pos) && b"ACGT".contains(&pre_base) {
                        for p in ref_pos..ref_pos + advance {
                            self.count(p);
                        }
                    }
                    ref_pos += advance;
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// One iteration of `reads_realignment`'s chunk loop: realign the windows of the current
    /// chunk's candidates, then write the reads the chunks have moved past.
    fn process_chunk(&mut self) -> Result<(), Error> {
        let Some((chunk_start, chunk_end)) = self.chunk else { return Ok(()) };
        let mut candidates: Vec<i64> = self
            .pileup
            .iter()
            .filter(|&&(pos, count)| {
                count as f64 >= MIN_COVERAGE && pos >= chunk_start - EXPANSION - 1 && pos < chunk_end + EXPANSION
            })
            .map(|&(pos, _)| pos)
            .collect();
        candidates.sort_unstable();
        if self.reads.is_empty() || candidates.is_empty() {
            return Ok(());
        }
        let splits = (chunk_end - chunk_start) / MAX_WINDOW;
        let mut windows: Vec<Window> = Vec::new();
        for split in 0..splits {
            let split_start = chunk_start + split * MAX_WINDOW - EXPANSION - 1;
            let split_end = split_start + MAX_WINDOW + EXPANSION * 2 + 1;
            let in_split: Vec<i64> =
                candidates.iter().copied().filter(|&p| split_start <= p && p < split_end).collect();
            windows.extend(self.split_windows(&in_split));
        }
        // A window reads only what realignment leaves alone (the reads' bases, positions and
        // CIGARs, the reference), so the windows run in parallel; each read's best alignment is
        // then updated in window order, as one window after the other would.
        let input =
            WindowInput { reference: &self.reference, reference_offset: self.reference_offset, reads: &self.reads };
        let realigned: Vec<Result<Vec<Realigned>, Error>> = windows.par_iter().map(|w| input.realign(w)).collect();
        for result in realigned {
            for (r, cigar, position) in result? {
                self.reads[r].set_realignment(&cigar, position);
            }
        }
        self.flush(chunk_start - EXPANSION - MAX_WINDOW);
        Ok(())
    }

    /// The windows of one 1 kb split, each with its reads.
    fn split_windows(&self, positions: &[i64]) -> Vec<Window> {
        let mut windows: Vec<(i64, i64)> = Vec::new();
        let mut span: Option<(i64, i64)> = None;
        for &pos in positions {
            span = Some(match span {
                None => (pos, pos),
                Some((start, end)) if pos > end + 2 * MIN_WINDOW_DISTANCE => {
                    windows.push((start - MIN_WINDOW_DISTANCE, end + MIN_WINDOW_DISTANCE));
                    (pos, pos)
                }
                Some((start, _)) => (start, pos),
            });
        }
        if let Some((start, end)) = span {
            windows.push((start - MIN_WINDOW_DISTANCE, end + MIN_WINDOW_DISTANCE));
        }
        if windows.is_empty() {
            return Vec::new();
        }
        windows.sort_by_key(|w| w.0);
        let max_window_end = windows.iter().map(|w| w.1).max().expect("non-empty");
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); windows.len()];
        for (r, read) in self.reads.iter().enumerate() {
            if read.start > max_window_end {
                continue;
            }
            // `find_max_overlap_index`: the first window with the largest overlap.
            let overlaps = windows.iter().map(|&(s, e)| (read.end.min(e) - read.start.max(s)).max(0));
            let (best, overlap) =
                overlaps.enumerate().fold((0, -1), |(bi, bo), (i, o)| if o > bo { (i, o) } else { (bi, bo) });
            if overlap > 0 {
                members[best].push(r);
            }
        }
        windows
            .into_iter()
            .zip(members)
            .filter(|&((start, end), _)| end - start <= MAX_WINDOW)
            .map(|((start, end), members)| Window { start, end, members })
            .collect()
    }
}

/// A read's new alignment from one window: its index, CIGAR (matches as `X`) and position.
type Realigned = (usize, String, i64);

/// A realignment window: its span and the indices of its reads.
#[derive(Debug)]
struct Window {
    start: i64,
    end: i64,
    members: Vec<usize>,
}

/// What realigning a window reads.
#[derive(Debug, Clone, Copy)]
struct WindowInput<'a> {
    reference: &'a [u8],
    reference_offset: i64,
    reads: &'a [Read],
}

impl WindowInput<'_> {
    /// The new `(read, cigar, position)` of each read of the window whose realignment differs
    /// from its alignment, in read order, for `set_realignment_info`.
    fn realign(&self, window: &Window) -> Result<Vec<Realigned>, Error> {
        let (start, end, members) = (window.start, window.end, window.members.as_slice());
        let offset = self.reference_offset;
        let window_ref = py_slice(self.reference, start - offset, end - offset).to_vec();
        let mut graph_reads: Vec<&[u8]> = Vec::new();
        let mut low_quality: Vec<BTreeSet<usize>> = Vec::new();
        for &r in members {
            let read = &self.reads[r];
            if read.mapq < GRAPH_MIN_MAPQ || read.start > end || read.end < start {
                continue;
            }
            graph_reads.push(&read.seq);
            low_quality.push(
                read.quals.iter().enumerate().filter(|&(_, &q)| q < GRAPH_MIN_BASE_QUALITY).map(|(i, _)| i).collect(),
            );
        }
        // The empty read list still splits into one empty read.
        if graph_reads.is_empty() {
            graph_reads.push(&[]);
            low_quality.push(BTreeSet::new());
        }
        let mut consensus = candidate_haplotypes(&window_ref, &graph_reads, &low_quality);
        // Python reads at most 200 of the C array's entries.
        consensus.truncate(200);
        if consensus.is_empty() || (consensus.len() == 1 && consensus[0] == window_ref) || members.is_empty() {
            return Ok(Vec::new());
        }
        let min_start = members.iter().map(|&r| self.reads[r].start).min().expect("non-empty");
        let max_end = members.iter().map(|&r| self.reads[r].end).max().expect("non-empty");
        let tmp_ref_start = (min_start.min(start) - EXPANSION).max(0);
        let tmp_ref_end = max_end.max(end) + EXPANSION;
        let get = |a: i64, b: i64| {
            let (a, b) = if b < a { (b, a) } else { (a, b) };
            py_slice(self.reference, a - offset, b - offset).to_vec()
        };
        let prefix = get(tmp_ref_start, start);
        let center = get(start, end);
        if tmp_ref_end < end {
            return Ok(Vec::new());
        }
        let suffix = get(end, tmp_ref_end);
        let ref_seq = [prefix.as_slice(), &center, &suffix].concat();
        let chosen: Vec<usize> = members.iter().copied().take(MAX_WINDOW as usize).collect();
        let window_reads: Vec<WindowRead<'_>> = chosen
            .iter()
            .map(|&r| {
                let read = &self.reads[r];
                WindowRead { seq: &read.seq, position: read.start as i32, cigar: &read.cigar }
            })
            .collect();
        let haplotypes: Vec<Vec<u8>> =
            consensus.iter().map(|c| [prefix.as_slice(), c.as_slice(), suffix.as_slice()].concat()).collect();
        let realigned = realign_reads(
            &window_reads,
            &ref_seq,
            &haplotypes,
            tmp_ref_start as i32,
            prefix.len() as i32,
            suffix.len() as i32,
        );
        let mut changed = Vec::new();
        for (&r, (position, cigar)) in chosen.iter().zip(realigned) {
            let read = &self.reads[r];
            let position = i64::from(position);
            if cigar.is_empty() || (read.cigar == cigar && read.start == position) {
                continue;
            }
            if cigar.replace('X', "M") == read.cigar && position == read.start {
                continue;
            }
            changed.push((r, cigar, position));
        }
        Ok(changed)
    }
}

impl RealignedReads<'_> {
    /// Writes, sorted by realigned start, the reads that start before `limit`, and forgets the
    /// pileup before it.
    fn flush(&mut self, limit: i64) {
        let mut order: Vec<usize> = (0..self.reads.len()).collect();
        order.sort_by_key(|&r| self.reads[r].best_pos);
        let mut written = vec![false; self.reads.len()];
        for &r in &order {
            if self.reads[r].best_pos < limit {
                self.out.push_back(self.reads[r].sam_line(&self.contig));
                written[r] = true;
            }
        }
        if written.iter().any(|&w| w) {
            let reads = std::mem::take(&mut self.reads);
            self.reads = reads.into_iter().zip(written).filter(|(_, w)| !w).map(|(r, _)| r).collect();
        }
        if !order.is_empty() {
            self.pileup.retain(|&(pos, _)| pos >= limit);
            self.pileup_index = self.pileup.iter().enumerate().map(|(i, &(pos, _))| (pos, i)).collect();
        }
    }

    /// Runs until there are lines to hand out or the reads are exhausted.
    fn fill(&mut self) -> Result<(), Error> {
        while self.out.is_empty() && !self.finished {
            match self.view.next() {
                Some(record) => self.add_record(&record?)?,
                None => {
                    self.finished = true;
                    if self.chunk.is_some() {
                        self.process_chunk()?;
                        // After the last chunk every remaining read is written.
                        self.flush(i64::MAX);
                    }
                }
            }
        }
        Ok(())
    }
}

impl Iterator for RealignedReads<'_> {
    type Item = Result<String, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Err(e) = self.fill() {
            self.finished = true;
            self.out.clear();
            return Some(Err(e));
        }
        self.out.pop_front().map(Ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_names_count_in_base_87() {
        assert_eq!(NAME_DIGITS.len(), 87);
        assert_eq!(read_name(0), "00000");
        assert_eq!(read_name(1), "00001");
        assert_eq!(read_name(87), "00010");
        assert_eq!(read_name(87 * 87 + 5), "00105");
    }

    #[test]
    fn soft_clip_threshold() {
        assert!(!mostly_soft_clipped("100M"));
        assert!(mostly_soft_clipped("60S40M"));
        assert!(!mostly_soft_clipped("40S60M"));
    }

    #[test]
    fn haplotype_tags() {
        assert_eq!(haplotype_tag(&["NM:i:0", "HP:i:2", "PS:i:5"]), Some('2'));
        assert_eq!(haplotype_tag(&["NM:i:0"]), None);
    }
}

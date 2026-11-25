//! mosdepth's outputs from per-contig coverage: region depths, thresholds, summary and
//! distributions (mosdepth.nim:589-840).
//!
//! [`ContigReport::new`] does the per-contig work and can run on any thread; [`Writer`]
//! takes reports in reference order and writes the files.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use pasta::bgzf;
use pasta::core::Position;
use pasta::csi::{self, binning_index::index::reference_sequence::bin::Chunk};

use crate::counter::Coverage;
use crate::regions::By;
use crate::stats::{CountStat, DepthStat, Distribution, SUMMARY_HEADER, format_float, nim_to_int};

/// What to report, besides the always-written summary and global distribution.
#[derive(Debug)]
pub struct ReportOptions {
    /// `--by`: windows or BED regions for `regions.bed.gz`.
    pub by: Option<By>,
    /// `-m`: report each region's median instead of its mean.
    pub use_median: bool,
    /// `-T`: for each region, count bases with at least these depths. Sorted ascending.
    pub thresholds: Vec<i32>,
    /// Decimal places (mosdepth's `MOSDEPTH_PRECISION`, default 2).
    pub precision: usize,
    /// The coverage counts only the reads of `by`'s regions (lrpipe's own mode, not mosdepth's):
    /// their rows are mosdepth's, but a contig's depth is not known, so the summary has only the
    /// regions' rows (`CONTIG_region`, `total_region`) and the global distribution is not written.
    pub regions_only: bool,
}

impl Default for ReportOptions {
    fn default() -> Self {
        ReportOptions { by: None, use_median: false, thresholds: Vec::new(), precision: 2, regions_only: false }
    }
}

impl ReportOptions {
    /// mosdepth skips contigs without BED regions entirely, when nothing per-base (per-base
    /// depth, thresholds, quantized) is written (mosdepth.nim:703-705): no summary line and no
    /// contribution to the totals.
    pub fn skips(&self, contig: &str) -> bool {
        self.thresholds.is_empty() && matches!(&self.by, Some(By::Bed(bed)) if !bed.contains(contig))
    }
}

/// Per-worker buffers reused across contigs.
pub struct Scratch {
    median: Option<CountStat>,
}

impl Scratch {
    pub fn new(opts: &ReportOptions) -> Self {
        // mosdepth.nim:686: medians are exact up to 65535.
        Scratch { median: opts.use_median.then(|| CountStat::new(65536)) }
    }
}

/// One line of a BED-like output, with its interval for the index.
struct Line {
    start: u32,
    stop: u32,
    end_byte: usize,
}

/// Everything mosdepth writes for one contig.
pub struct ContigReport {
    name: String,
    found: bool,
    region_text: Vec<u8>,
    region_lines: Vec<Line>,
    threshold_text: Vec<u8>,
    threshold_lines: Vec<Line>,
    chrom_stat: DepthStat,
    region_stat: DepthStat,
    global_dist: Distribution,
    region_dist: Distribution,
}

impl std::fmt::Debug for ContigReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ContigReport")
            .field("name", &self.name)
            .field("found", &self.found)
            .field("regions", &self.region_lines.len())
            .field("chrom_stat", &self.chrom_stat)
            .finish_non_exhaustive()
    }
}

impl ContigReport {
    /// This contig's `regions.bed.gz` lines, as mosdepth writes them (empty without `--by`).
    pub fn region_text(&self) -> &[u8] {
        &self.region_text
    }

    /// Computes the report for contig `name` of `len` bases from its coverage.
    pub fn new(name: &str, len: u32, cov: &Coverage, opts: &ReportOptions, scratch: &mut Scratch) -> Self {
        let mut r = ContigReport {
            name: name.to_owned(),
            found: cov.found,
            region_text: Vec::new(),
            region_lines: Vec::new(),
            threshold_text: Vec::new(),
            threshold_lines: Vec::new(),
            chrom_stat: DepthStat::default(),
            region_stat: DepthStat::default(),
            global_dist: Distribution::default(),
            region_dist: Distribution::default(),
        };
        let arr = &cov.depth[..];
        let n = arr.len() as u32;

        if let Some(by) = &opts.by {
            let windows = by.is_windows();
            let mut value = String::new();
            for reg in by.regions(name, len) {
                // mosdepth.nim:718-722: the value stays 0 on contigs without reads.
                let mut me = 0.0;
                if cov.found {
                    me = region_value(arr, reg.start, reg.stop, scratch.median.as_mut());
                    let slice = &arr[reg.start.min(n) as usize..n.min(reg.stop) as usize];
                    r.region_stat.add(&DepthStat::of(slice));
                }
                value.clear();
                format_float(&mut value, me, opts.precision);
                let t = &mut r.region_text;
                write!(t, "{name}\t{}\t{}\t", reg.start, reg.stop).unwrap();
                if let Some(rn) = &reg.name {
                    t.extend_from_slice(rn.as_bytes());
                    t.push(b'\t');
                }
                t.extend_from_slice(value.as_bytes());
                t.push(b'\n');
                r.region_lines.push(Line { start: reg.start, stop: reg.stop, end_byte: t.len() });

                if cov.found {
                    if windows {
                        r.region_dist.add_capped(nim_to_int(me));
                    } else if (reg.start as usize) < arr.len() {
                        let stop = (reg.stop as usize).min(arr.len());
                        r.region_dist.add_depths(&arr[reg.start as usize..stop]);
                    } else {
                        eprintln!(
                            "[lrpipe-depth] warning requested interval outside of chromosome range:{}..{}",
                            reg.start, reg.stop
                        );
                    }
                }

                if !opts.thresholds.is_empty() {
                    write_thresholds(
                        &mut r,
                        arr,
                        cov.found,
                        &opts.thresholds,
                        reg.start,
                        reg.stop,
                        reg.name.as_deref(),
                    );
                }
            }
        }

        if cov.found {
            let bases = cov.bases();
            r.global_dist.add_depths(bases);
            r.chrom_stat = DepthStat::of(bases);
        }
        r
    }
}

/// mosdepth's `imean` (mosdepth.nim:399-413): the mean adds each depth divided by the
/// region length, in order, so the floating-point result matches bit for bit.
fn region_value(arr: &[i32], start: u32, stop: u32, median: Option<&mut CountStat>) -> f64 {
    if start as usize > arr.len() {
        return 0.0;
    }
    let vals = &arr[start as usize..(stop as usize).min(arr.len())];
    match median {
        Some(ms) => {
            ms.clear();
            for &v in vals {
                ms.add(v);
            }
            ms.median() as f64
        }
        None => {
            let l = f64::from(stop - start);
            let mut acc = 0.0f64;
            for &v in vals {
                acc += f64::from(v) / l;
            }
            acc
        }
    }
}

/// mosdepth.nim:522-557.
fn write_thresholds(
    r: &mut ContigReport,
    arr: &[i32],
    found: bool,
    thresholds: &[i32],
    start: u32,
    stop: u32,
    name: Option<&str>,
) {
    let t = &mut r.threshold_text;
    write!(t, "{}\t{start}\t{stop}\t{}", r.name, name.unwrap_or("unknown")).unwrap();
    if !found {
        for _ in thresholds {
            t.extend_from_slice(b"\t0");
        }
    } else {
        let mut counts = vec![0u64; thresholds.len()];
        let vals = &arr[(start as usize).min(arr.len())..(stop as usize).min(arr.len())];
        for &v in vals {
            for (c, &th) in counts.iter_mut().zip(thresholds) {
                if v < th {
                    break;
                }
                *c += 1;
            }
        }
        for c in counts {
            write!(t, "\t{c}").unwrap();
        }
    }
    t.push(b'\n');
    r.threshold_lines.push(Line { start, stop, end_byte: t.len() });
}

/// A BGZF-compressed BED-like file with a CSI index, like hts-nim's `wopen_bgzi` makes.
struct IndexedBed {
    path: PathBuf,
    writer: bgzf::io::Writer<BufWriter<File>>,
    indexer: csi::binning_index::Indexer<csi::binning_index::index::reference_sequence::index::BinnedIndex>,
    names: Vec<String>,
}

impl IndexedBed {
    fn create(path: PathBuf, levels: u8) -> Result<Self> {
        let file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
        let writer = bgzf::io::Writer::new(BufWriter::new(file));
        let indexer = csi::binning_index::Indexer::new(14, levels).context("contig too long to index")?;
        Ok(IndexedBed { path, writer, indexer, names: Vec::new() })
    }

    fn write_lines(&mut self, contig: &str, text: &[u8], lines: &[Line]) -> Result<()> {
        if lines.is_empty() {
            return Ok(());
        }
        if self.names.last().map(String::as_str) != Some(contig) {
            self.names.push(contig.to_owned());
        }
        let id = self.names.len() - 1;
        let mut from = 0;
        for l in lines {
            let start_vpos = self.writer.virtual_position();
            self.writer.write_all(&text[from..l.end_byte])?;
            let end_vpos = self.writer.virtual_position();
            from = l.end_byte;
            // tabix indexes [start, stop) as 1-based [start + 1, stop], at least one base.
            let start = Position::try_from(l.start as usize + 1)?;
            let end = Position::try_from((l.stop as usize).max(l.start as usize + 1))?;
            self.indexer
                .add_record(Some((id, start, end, true)), Chunk::new(start_vpos, end_vpos))
                .map_err(|e| anyhow::anyhow!("indexing {}: {e:?}", self.path.display()))?;
        }
        Ok(())
    }

    fn finish(self) -> Result<()> {
        let IndexedBed { path, writer, indexer, names } = self;
        writer.finish()?.flush()?;
        let header = csi::binning_index::index::header::Builder::bed()
            .set_reference_sequence_names(names.iter().map(|n| n.as_str().into()).collect())
            .build();
        let index = indexer.set_header(header).build(names.len());
        let mut csi_path = path.into_os_string();
        csi_path.push(".csi");
        csi::fs::write(&csi_path, &index).with_context(|| format!("writing {}", Path::new(&csi_path).display()))?;
        Ok(())
    }
}

/// The CSI depth needed for the longest contig (mosdepth.nim:565-577).
pub fn csi_levels(max_contig_len: u64) -> u8 {
    let mut levels = 0;
    let mut s = 1u64 << 14;
    while max_contig_len > s {
        levels += 1;
        s <<= 3;
    }
    levels
}

/// Writes mosdepth's output files under `prefix`, from reports given in reference order.
pub struct Writer {
    precision: usize,
    has_regions: bool,
    regions: Option<IndexedBed>,
    thresholds: Option<IndexedBed>,
    summary: BufWriter<File>,
    /// None with [`ReportOptions::regions_only`].
    global_dist_out: Option<BufWriter<File>>,
    region_dist_out: Option<BufWriter<File>>,
    global_stat: DepthStat,
    global_region_stat: DepthStat,
    global_dist: Distribution,
    region_dist: Distribution,
    summary_header: bool,
    regions_only: bool,
}

fn create(path: String) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(&path).with_context(|| format!("creating {path}"))?))
}

impl Writer {
    /// Creates the output files. `max_contig_len` sizes the CSI indexes.
    pub fn new(prefix: &str, opts: &ReportOptions, max_contig_len: u64) -> Result<Self> {
        let levels = csi_levels(max_contig_len);
        let has_regions = opts.by.is_some();
        let regions = match has_regions {
            true => Some(IndexedBed::create(format!("{prefix}.regions.bed.gz").into(), levels)?),
            false => None,
        };
        let thresholds = match opts.thresholds.is_empty() {
            false => {
                let mut t = IndexedBed::create(format!("{prefix}.thresholds.bed.gz").into(), levels)?;
                let mut header = String::from("#chrom\tstart\tend\tregion");
                for th in &opts.thresholds {
                    header.push_str(&format!("\t{th}X"));
                }
                header.push('\n');
                t.writer.write_all(header.as_bytes())?;
                Some(t)
            }
            true => None,
        };
        Ok(Writer {
            precision: opts.precision,
            has_regions,
            regions,
            thresholds,
            global_dist_out: match opts.regions_only {
                false => Some(create(format!("{prefix}.mosdepth.global.dist.txt"))?),
                true => None,
            },
            summary: create(format!("{prefix}.mosdepth.summary.txt"))?,
            region_dist_out: match has_regions {
                true => Some(create(format!("{prefix}.mosdepth.region.dist.txt"))?),
                false => None,
            },
            global_stat: DepthStat::default(),
            global_region_stat: DepthStat::default(),
            global_dist: Distribution::default(),
            region_dist: Distribution::default(),
            summary_header: true,
            regions_only: opts.regions_only,
        })
    }

    fn summary_line(&mut self, name: &str, stat: &DepthStat) -> io::Result<()> {
        if self.summary_header {
            self.summary.write_all(SUMMARY_HEADER.as_bytes())?;
            self.summary_header = false;
        }
        self.summary.write_all(stat.summary_line(name, self.precision).as_bytes())
    }

    /// Writes one contig's output (mosdepth.nim:733-764).
    pub fn write_contig(&mut self, r: &ContigReport) -> Result<()> {
        if let Some(f) = &mut self.regions {
            f.write_lines(&r.name, &r.region_text, &r.region_lines)?;
        }
        if let Some(f) = &mut self.thresholds {
            f.write_lines(&r.name, &r.threshold_text, &r.threshold_lines)?;
        }
        if r.found {
            self.global_stat.add(&r.chrom_stat);
            if !self.regions_only {
                self.summary_line(&r.name, &r.chrom_stat.clone())?;
            }
            if self.has_regions {
                self.summary_line(&format!("{}_region", r.name), &r.region_stat.clone())?;
            }
            self.global_region_stat.add(&r.region_stat);
        }
        if let Some(out) = &mut self.global_dist_out {
            r.global_dist.write(&r.name, self.precision, out)?;
        }
        if let Some(out) = &mut self.region_dist_out {
            r.region_dist.write(&r.name, self.precision, out)?;
        }
        if r.found {
            if self.has_regions {
                self.region_dist.merge(&r.region_dist);
            }
            self.global_dist.merge(&r.global_dist);
        }
        Ok(())
    }

    /// Writes the totals and closes the files (mosdepth.nim:807-840).
    pub fn finish(mut self) -> Result<()> {
        let (g, gr) = (self.global_stat, self.global_region_stat);
        if !self.regions_only {
            self.summary_line("total", &g)?;
        }
        if self.has_regions {
            self.summary_line("total_region", &gr)?;
        }
        if let Some(out) = &mut self.global_dist_out {
            self.global_dist.write("total", self.precision, out)?;
            out.flush()?;
        }
        if let Some(out) = &mut self.region_dist_out {
            self.region_dist.write("total", self.precision, out)?;
            out.flush()?;
        }
        self.summary.flush()?;
        if let Some(f) = self.regions {
            f.finish()?;
        }
        if let Some(f) = self.thresholds {
            f.finish()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_like_mosdepth() {
        assert_eq!(csi_levels(16384), 0);
        assert_eq!(csi_levels(248_956_422), 5);
        assert_eq!(csi_levels(1 << 29), 5);
        assert_eq!(csi_levels((1 << 29) + 1), 6);
    }

    #[test]
    fn mean_divides_before_adding() {
        let arr = [1, 2, 3, 0];
        let m = region_value(&arr, 0, 3, None);
        assert_eq!(m, 1.0 / 3.0 + 2.0 / 3.0 + 3.0 / 3.0);
        assert_eq!(region_value(&arr, 5, 10, None), 0.0);
    }
}

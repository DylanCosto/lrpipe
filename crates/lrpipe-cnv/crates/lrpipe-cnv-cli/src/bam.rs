//! Window depths from a BAM or CRAM with lrpipe-depth, identical to
//! `mosdepth -n -x -Q 20 -b 1000`'s `regions.bed.gz`.

use std::path::Path;

use anyhow::Result;
use lrpipe_cnv::coverage::{ContigDepth, parse_regions};
use lrpipe_depth::counter::{Filter, Mode};
use lrpipe_depth::reader::Alignments;
use lrpipe_depth::regions::By;
use lrpipe_depth::report::ReportOptions;
use lrpipe_depth::run::{RunOptions, count_contigs};

/// Counts `path`'s depth in windows of `bin_size` over reads of MAPQ ≥ `min_mapq`, contigs in
/// parallel on `threads` threads. `reference` decodes CRAM.
pub fn window_depths(
    path: &Path,
    reference: Option<&Path>,
    bin_size: u32,
    min_mapq: u8,
    threads: usize,
) -> Result<Vec<ContigDepth>> {
    let aln = Alignments::open(path, reference)?;
    let opts = RunOptions {
        mode: Mode::Fast,
        filter: Filter { min_mapq, ..Filter::default() },
        report: ReportOptions { by: Some(By::Windows(bin_size)), ..ReportOptions::default() },
        chrom: None,
        threads,
    };
    let work: Vec<usize> = (0..aln.contigs().len()).collect();
    let mut text = Vec::new();
    count_contigs(&aln, &work, &opts, |_, report| {
        text.extend_from_slice(report.region_text());
        Ok(())
    })?;
    Ok(parse_regions(text.as_slice(), path)?)
}

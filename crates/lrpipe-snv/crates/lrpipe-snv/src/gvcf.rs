//! `lrpipe-snv run --gvcf`: the gVCF steps of the original's pipeline.
//!
//! Each pileup chunk writes its non-variant blocks to `tmp/gvcf_tmp_output/<sample>.<contig>.
//! <start>.<end>.tmp.gvcf`; after full-alignment calling they are sorted into
//! `tmp/gvcf_tmp_output/non_var.gvcf`, merged with each contig's calls into
//! `tmp/merge_output/merge_<contig>.gvcf`, and those into `merge_output.gvcf.gz`. Like the
//! original with `lz4` installed, the intermediate files are lz4 frames.

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use lrpipe_reads::regions::Regions;
use snv_call::{Likelihoods, Lz4Writer, NonVariantOptions, NonVariantWriter, Resolution, chunk_gvcf_name};
use snv_merge::{GvcfContigs, GvcfOutput, SortOptions, sort_vcf};
use snv_mpileup::{GvcfChunk, GvcfSites};
use snv_pileup::{ChunkSites, Site};

use crate::settings::str2bool;

/// The gVCF options of the original's run script.
#[derive(Debug, Clone, clap::Args)]
pub(crate) struct GvcfArgs {
    /// Also write `merge_output.gvcf.gz`, with reference blocks confined to any `--bed-fn` targets.
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    gvcf: bool,
    /// Estimated base error rate for the non-variant likelihoods (`--base_err`).
    #[arg(long, default_value_t = 0.001)]
    base_err: f64,
    /// GQ bin width that splits non-variant blocks (`--gq_bin_size`).
    #[arg(long, default_value_t = 5)]
    gq_bin_size: i64,
    /// Keep every reference contig in the gVCF header, not only the called ones
    /// (`--output_all_contigs_in_gvcf_header`).
    #[arg(long, num_args = 0..=1, default_value = "false", default_missing_value = "true", value_parser = str2bool)]
    output_all_contigs_in_gvcf_header: bool,
}

/// Where the gVCF steps read and write, when `--gvcf` is on.
#[derive(Debug, Clone)]
pub(crate) struct Gvcf {
    /// `tmp/gvcf_tmp_output`.
    dir: PathBuf,
    options: NonVariantOptions,
    contigs: GvcfContigs,
    targets: Option<Regions>,
}

/// `tmp/gvcf_tmp_output`, which `CheckEnvs` creates with or without `--gvcf`.
fn tmp_dir(tmp: &Path) -> PathBuf {
    tmp.join("gvcf_tmp_output")
}

/// Creates `tmp/gvcf_tmp_output`, as `CheckEnvs` does for every run.
pub(crate) fn create_tmp_dir(tmp: &Path) -> Result<()> {
    let dir = tmp_dir(tmp);
    std::fs::create_dir_all(&dir).with_context(|| format!("creating {}", dir.display()))
}

impl GvcfArgs {
    /// Whether `--gvcf` is on.
    pub(crate) fn on(&self) -> bool {
        self.gvcf
    }

    /// The PL field of the pileup and full-alignment VCFs.
    pub(crate) fn likelihoods(&self) -> Likelihoods {
        if self.gvcf { Likelihoods::Write } else { Likelihoods::Omit }
    }

    /// The gVCF steps under `tmp`, or `None` without `--gvcf`.
    pub(crate) fn steps(&self, tmp: &Path, bed: Option<&[(String, i64, i64)]>) -> Option<Gvcf> {
        if self.gq_bin_size < 0 {
            log::warn!("[WARNING] Invalid gq bin size --gq_bin_size>=0");
        }
        self.gvcf.then(|| Gvcf {
            dir: tmp_dir(tmp),
            options: NonVariantOptions {
                base_error: self.base_err,
                // The original's pipeline script raises a negative size to 0, which then fails as
                // the original does.
                gq_bin_size: self.gq_bin_size.max(0),
                resolution: Resolution::Blocks,
            },
            contigs: if self.output_all_contigs_in_gvcf_header { GvcfContigs::All } else { GvcfContigs::Processed },
            targets: bed.map(|rows| {
                // A zero-length BED row names one base, as in `Regions::read` and `--regions`.
                Regions::from_intervals(
                    rows.iter()
                        .map(|(contig, start, end)| (contig.clone(), *start, (*end).max(start.saturating_add(1)))),
                )
            }),
        })
    }
}

impl Gvcf {
    /// `tmp/gvcf_tmp_output/non_var.gvcf`.
    fn non_variant(&self) -> PathBuf {
        self.dir.join("non_var.gvcf")
    }

    /// Writes one pileup chunk's non-variant blocks, as `CreateTensorPileup --gvcf` does.
    pub(crate) fn write_chunk(&self, reference: &Path, sample: &str, sites: &ChunkSites) -> Result<()> {
        let blocks = (!sites.is_empty()).then(|| sites.sites());
        self.write_blocks(reference, sample, (sites.contig(), sites.start(), sites.end()), blocks)
    }

    /// [`write_chunk`](Self::write_chunk) for an Illumina chunk (`CreateTensorPileup.py`).
    pub(crate) fn write_mpileup_chunk(&self, reference: &Path, sample: &str, chunk: &GvcfChunk) -> Result<()> {
        let blocks = match &chunk.sites {
            GvcfSites::Sites(sites) => Some(sites.iter().copied()),
            GvcfSites::Empty => None,
        };
        self.write_blocks(reference, sample, (&chunk.contig, chunk.start, chunk.end), blocks)
    }

    /// Writes `<sample>.<contig>.<start>.<end>.tmp.gvcf`: the blocks of `sites`, or with `None`
    /// (no reads in the chunk) the one block of `write_empty_pileup`.
    fn write_blocks(
        &self,
        reference: &Path,
        sample: &str,
        (contig, start, end): (&str, i64, i64),
        sites: Option<impl Iterator<Item = Site>>,
    ) -> Result<()> {
        let name = chunk_gvcf_name(sample, contig, start, end);
        let path = self.dir.join(format!("{name}.tmp.gvcf"));
        let file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
        let out = Lz4Writer::new(BufWriter::new(file));
        let mut writer = NonVariantWriter::new(out, reference, &name, contig, &self.options)?;
        match (&self.targets, sites) {
            (Some(targets), sites) => {
                let targets = targets.contig(contig);
                let mut sites = sites.map(Iterator::peekable);
                // Convert the chunk's 1-based, end-exclusive range to BED coordinates.
                let (start, end) = (start.max(1) - 1, end - 1);
                let first = targets.partition_point(|&(_, stop)| stop <= start);
                for &(left, right) in targets[first..].iter().take_while(|&&(left, _)| left < end) {
                    let first_pos = start.max(left) + 1;
                    let last_pos = end.min(right);
                    writer.start_interval(last_pos)?;
                    if let Some(sites) = &mut sites {
                        while sites.next_if(|site| site.pos < first_pos).is_some() {}
                        while let Some(site) = sites.next_if(|site| site.pos <= last_pos) {
                            writer.push(site.pos, site.reference, site.n_ref, site.n_total)?;
                        }
                    } else {
                        writer.write_empty(first_pos, last_pos + 1)?;
                    }
                }
            }
            (None, Some(sites)) => {
                for site in sites {
                    writer.push(site.pos, site.reference, site.n_ref, site.n_total)?;
                }
            }
            (None, None) => writer.write_empty(start, end)?,
        }
        writer.finish()?.finish().with_context(|| format!("writing {}", path.display()))?;
        Ok(())
    }

    /// `SortVcf --vcf_fn_suffix .tmp.gvcf`: all chunks' blocks into `non_var.gvcf`.
    pub(crate) fn sort_chunks(&self, reference: &Path, contigs_file: &Path, cmd_file: &Path) -> Result<()> {
        let mut sort = SortOptions::new(&self.dir, self.non_variant());
        sort.suffix = ".tmp.gvcf".to_owned();
        sort.reference = Some(reference.to_owned());
        sort.contigs = Some(contigs_file.to_owned());
        sort.cmd_file = Some(cmd_file.to_owned());
        sort.gvcf_contigs = self.contigs;
        let summary = sort_vcf(&sort)?;
        log::info!("merged {} non-variant records into {}", summary.records, summary.vcf.display());
        Ok(())
    }

    /// `MergeVcf --gvcf True`: `merge_<contig>.gvcf` next to the contig's merged VCF.
    pub(crate) fn merge_output(&self, merged_vcf: &Path, reference: &Path) -> GvcfOutput {
        GvcfOutput::On {
            gvcf: merged_vcf.with_extension("gvcf"),
            non_variant: self.non_variant(),
            reference: reference.to_owned(),
        }
    }

    /// `SortVcf --vcf_fn_prefix merge --vcf_fn_suffix .gvcf`: `merge_output.gvcf.gz`.
    pub(crate) fn sort_merged(&self, merge_dir: &Path, output: &Path, sort: &SortOptions) -> Result<()> {
        let mut sort = sort.clone();
        sort.input_dir = merge_dir.to_owned();
        sort.output = output.join("merge_output.gvcf");
        sort.prefix = Some("merge".to_owned());
        sort.suffix = ".gvcf".to_owned();
        sort.gvcf_contigs = self.contigs;
        let summary = sort_vcf(&sort)?;
        log::info!("{} gVCF records in {}", summary.records, summary.vcf.display());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use snv_call::open_lz4;
    use std::io::BufRead;

    #[test]
    fn reference_blocks_stop_at_bed_boundaries_with_or_without_reads() {
        let tmp = std::env::temp_dir().join(format!("lrpipe-snv-gvcf-targets-{}", std::process::id()));
        create_tmp_dir(&tmp).expect("gVCF scratch directory");
        let reference = tmp.join("ref.fa");
        std::fs::write(&reference, ">chrT\nACGTACGTAC\n").expect("reference");
        std::fs::write(reference.with_extension("fa.fai"), "chrT\t10\t6\t10\t11\n").expect("FASTA index");
        let bed = [
            ("chrT".to_owned(), 6, 9),
            ("chrT".to_owned(), 1, 3),
            ("chrT".to_owned(), 0, 2),
            ("chrT".to_owned(), 3, 4),
            ("chrT".to_owned(), 5, 5),
        ];
        let args = GvcfArgs { gvcf: true, base_err: 0.001, gq_bin_size: 100, output_all_contigs_in_gvcf_header: false };
        let gvcf = args.steps(&tmp, Some(&bed)).expect("gVCF enabled");
        let read_blocks = |(start, end), sites| {
            let chunk = GvcfChunk { contig: "chrT".to_owned(), start, end, sites };
            gvcf.write_mpileup_chunk(&reference, "sample", &chunk).expect("gVCF chunk");
            let name = chunk_gvcf_name("sample", "chrT", start, end);
            open_lz4(&gvcf.dir.join(format!("{name}.tmp.gvcf")))
                .expect("compressed gVCF")
                .lines()
                .map(|line| line.expect("record"))
                .filter(|line| !line.starts_with('#'))
                .map(|line| {
                    let fields: Vec<_> = line.split('\t').collect();
                    let start = fields[1].parse::<i64>().expect("POS");
                    let end = fields[7].strip_prefix("END=").expect("END").parse::<i64>().expect("integer");
                    let min_depth = fields[9].split(':').nth(2).expect("MIN_DP").parse::<u64>().expect("integer");
                    (start, end, min_depth)
                })
                .collect::<Vec<_>>()
        };

        // Unsorted, overlapping and adjacent rows give positions 1–4 and 6–9.
        // The two excluded bases would otherwise join nearby sites and lower their MIN_DP.
        let sites = (1..=10)
            .zip(*b"ACGTACGTAC")
            .map(|(pos, reference)| {
                let depth = if matches!(pos, 5 | 10) { 8 } else { 9 };
                Site { pos, reference, n_ref: depth, n_total: depth }
            })
            .collect();
        assert_eq!(read_blocks((1, 11), GvcfSites::Sites(sites)), [(1, 4, 9), (6, 9, 9)]);
        assert_eq!(read_blocks((1, 11), GvcfSites::Empty), [(1, 4, 0), (6, 9, 0)]);
        assert_eq!(read_blocks((3, 8), GvcfSites::Empty), [(3, 4, 0), (6, 7, 0)]);
        assert!(read_blocks((5, 6), GvcfSites::Empty).is_empty(), "a chunk in the BED gap");
        std::fs::remove_dir_all(tmp).expect("remove scratch directory");
    }
}

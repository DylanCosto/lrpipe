# Parity test fixture

Used by `tests/parity.rs`; regenerate with `tools/make_fixture.sh` (needs samtools, pysam and
Sniffles 2.8.1).

- `sample1.bam`, `sample1.cram`: the primary and supplementary alignments overlapping
  chr17:46,200,000-46,400,000 from ONT's HG002 giab_2025.01 haplotagged CRAM (R10.4.1, sup,
  GRCh38 without alts). Base qualities are removed and only the SA, NM, HP and PS tags kept.
  The CRAM is the same reads, written by samtools against the padded reference below (it has no
  NM tags, as is usual for CRAM, so sniffles-rs computes NM as htslib does).
- `sample2.cram`: the reads not on haplotype 1, as a second sample for multi-sample calling.
- `reference.window.fa`: GRCh38 chr17:46,158,127-46,442,419, the bases the reads cover. The
  tests rebuild the full-length chr17 from it, padded with N (its MD5 is in the BAM header).
- `regions.bed`: a sub-window for `--regions`.
- `run_regions.bed`: two regions for `lrpipe run --regions` (`lrpipe/tests/run.rs`): the first
  holds an insertion and the end of a duplication that starts before it, the second a deletion;
  `expected/regions_run.vcf` is Sniffles 2.8.1's `--regions` calls with the reference and
  `--sample-id S`, as `lrpipe run` calls them.
- `split.bam`, `split.foreign.bam`, `split.haplotags.tsv`: the reads of `sample1.bam` around
  chr17:46,237,501, some split at its 724 bp deletion into a primary and a supplementary
  alignment (with SA tags) and some with a secondary alignment, tagged as `whatshap haplotag`
  tags them (primary alignments only), tagged HP 1 throughout as by another pipeline, and the
  haplotag list of the first (`tools/make_split_fixture.sh`, which also writes
  `expected/split.vcf`: Sniffles 2.8.1 `--phase` on `split.bam`).
- `expected/*.vcf`: what Sniffles 2.8.1 writes for these inputs (see `tools/make_fixture.sh` for
  each command), except `accuracy.vcf`, a snapshot of `sniffles-rs --accuracy`.

The reads are from Oxford Nanopore Technologies' open data release for the GIAB HG002 sample
(`s3://ont-open-data/giab_2025.01/`).

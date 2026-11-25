# Run modes, standalone outputs and targeted gVCFs

> **Status.** This records lrpipe 0.2.1 as released. Since then: the stale repeat index
> (#42, fixed in #50), damaged outputs accepted on resume (#46, #49), the standalone phasing
> `.fai` check (#43, #47), late XY repeat failures with unknown PARs (#44, #48), CUDA calls at
> high depth (#45, #51) and bonefish's MM/ML tags ([bonefish#46](https://github.com/DylanCosto/bonefish/issues/46)) are fixed on `main`, and so are the SV `.vcf.bgz` index
> (#52) and `--pileup-only --gvcf`
> (#54), both in 0.2.2, and gVCF blocks between
> disjoint targets (#53), after 0.2.2. Code links point at the 0.2.1 commit.

The follow-up uses the published lrpipe 0.2.1 Linux x86-64 package. The core matrix crosses two technologies, 12 read layouts, four region types, ten analysis and output presets, three sex settings, two SV modes and both settings for intermediate-file retention. The presets cover the default run, gVCF output, phased SV calling, both routes to repeat genotyping, and standalone CNV, methylation and depth selections. [Earlier option matrix](README.md).

The [coverage table](run-mode-coverage.tsv) records the outcomes of all six matrices. The [run-mode breakdown](run-modes-by-input.tsv) groups cases by platform, input layout and requested sex.

The 11,520 compact runs yielded 9,024 passes, 960 expected preflight refusals and 1,536 late unknown-PAR refusals. All 14,976 common VCF streams agreed across 4,512 pairs of retention settings. Automatic sex and explicit XX agreed across 12,672 VCF streams for 3,520 pairs on the compact fixture with an XX signal.

## Genotype ploidy and resource settings

The compact test reference contains a short custom chrX. Forced XY repeat calling reaches a runtime refusal for unknown pseudoautosomal regions. The separate repeat stage failed in 1,152 cases and the shared read pass failed in 384. The small-variant stage uses diploid sex-chromosome calling on that custom reference and the pipeline logs a warning. Small-variant planning, repeat reference handling.

All 48 GRCh38 controls passed with HG002 HiFi BAM and HG003 ONT CRAM. They crossed automatic sex, forced XX and forced XY with gVCF on and off, two target BEDs, and both routes to repeat genotyping. The target BEDs contain chrX PAR and non-PAR sequence and a chrY interval. One BED crosses the PAR1 boundary. The tests checked the emitted sex, PAR bounds, repeat karyotype and called genotype arity. Indexed queries verified the diploid-to-haploid boundary in all 16 XY gVCFs. XX repeat files include uncalled chrY catalog rows.

The scheduler tests ran 108 combinations across BAM, CRAM and tagged FASTQ. The test wrapper held the start of each real calling command for 350 ms to observe overlap. SNV and read-pass processes serialized in all 36 runs requesting 1 GiB. They overlapped in all 36 runs requesting 32 GiB and all 36 runs using the default budget. Explicit one-thread and two-thread settings and the affinity-aware default reached both stage commands. These are scheduler measurements.

## gVCF blocks between disjoint targets

For both real datasets, two 30 kb targets on chrX led to positive-depth reference genotype blocks across more than 100 Mb outside the BED. The indexed HG002 HiFi lookup at `chrX:54215050` returned a reference block with `GT=0`, `GQ=99` and `MIN_DP=28`. The queried site lies about 53.8 Mb from the nearest target edge. The count table uses off-target portions of each record, including partial overhangs. [Per-contig counts](targeted-gvcf.tsv).

All eight compact BAM and CRAM checks with two disjoint chr20 targets reproduced reference genotypes in the gap. The ordinary and phased variant calls stayed within the BED. In the matched compact ONT comparison, Clair3 2.0.3 agreed with lrpipe on coordinates, alleles, genotypes and minimum depth for all 2,393 gVCF rows, including 2,350 rows wholly outside the BED. A second comparison used the HG003 ONT CRAM and the chrX PAR1 boundary panel. Both callers wrote the same eight wholly off-target chrX block intervals.

To reproduce the disjoint behavior, use an indexed BAM and a reference with reads across the intervening span. Write the run to a scratch directory and query a position in the gap.

```bash
printf 'chr20\t511000\t512000\nchr20\t530000\t531000\n' > "$TARGETS"
"$LRPIPE" run --config "$SITE" --reads "$BAM" --platform ont --sample "$SAMPLE" \
  --regions "$TARGETS" --analyses snv --gvcf --sex female --output "$RUN_ROOT"
bcftools query -r chr20:520000-520000 -f '%CHROM\t%POS\t%END[\t%GT\t%MIN_DP]\n' \
  "$RUN_ROOT/snv/merge_output.gvcf.gz"
```

The published caller derives pileup chunk bounds from the span of the BED. Variant candidate selection checks the interval union, and the nonvariant gVCF writer receives the larger pileup span. Chunk and candidate code, gVCF writer.

## Standalone small variants

All 512 covered-panel runs and all 128 runs with either no reads or no variants completed. The covered-panel matrix crosses HiFi and ONT, BAM and CRAM, CPU and CUDA, read order, approximation, pileup-only, final phasing, gVCF and cleanup. Across 256 matched CPU and CUDA pairs on that panel, all 512 common VCF streams agreed on coordinates, alleles, genotypes and phasing. Four device replays also appeared as active SNV processes on the NVIDIA L40S.

With both `--pileup-only` and `--gvcf`, all 128 tested runs with pileup records stopped without a final gVCF. All 32 tested runs with those switches and zero pileup records emitted a gVCF. Pileup exit paths.

## Standalone SV `.vcf.bgz` indexing

Of 512 output combinations, 320 passed the output checks, 64 refused requests that selected no output, and 128 wrote the index under the wrong extension for `.vcf.bgz`. All 256 emitted SNFs could be used as inputs to a subsequent SV cohort call. Within each of 64 matching input and calling groups, direct VCF records agreed across output format and SNF settings. [Per-format counts](standalone-sv-output.tsv).

The standalone SV command writes `calls.vcf.gz.tbi` for `calls.vcf.bgz`. `bcftools view -r` could not discover that index for `.bgz`. Explicitly loading it returned all four test sites. A second run with the same stem and a `.bgz` filename replaced an existing `.gz` file's index. Its sequential reads still returned four records. Indexed queries returned four before the second write and zero after it. The same collision affected CSI indexes when a controlled BAM header declared a 600 Mb contig. Index path construction.

This reproduction uses a scratch directory and a BAM with chr20 SV calls. The second interval has no SV calls on the tested fixture.

```bash
call=("$LRPIPE" sv --input "$BAM" --reference "$REFERENCE" --sample-id "$SAMPLE")
"${call[@]}" --vcf "$SCRATCH/calls.vcf.gz"
bcftools view -H -r chr20 "$SCRATCH/calls.vcf.gz" | wc -l
"${call[@]}" --region chr20:1000-2000 --vcf "$SCRATCH/calls.vcf.bgz"
bcftools view -H -r chr20 "$SCRATCH/calls.vcf.gz" | wc -l
bcftools view -H "$SCRATCH/calls.vcf.gz" | wc -l
```

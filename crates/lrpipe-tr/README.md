# lrpipe-tr

Tandem repeat genotyping from long reads (PacBio HiFi and ONT), in pure Rust: allele
sequences, lengths and motif counts at every locus of a repeat catalog, from pathogenic
expansions (STRchive) to genome-wide catalogs (adotto, 937,204 loci).

Written from first principles. No code of TRGT (PacBio licence, PacBio data only, no
redistribution), LongTR (GPLv2), Straglr or STRkit (GPLv3) is read or used; TRGT's output is
used only as a black-box comparison.

```bash
lrpipe-tr --reads sample.bam --genome ref.fa --repeats catalog.bed --output sample.vcf.gz \
    --sample HG002 --karyotype xy -t 16      # --platform ont for Nanopore reads
```

The catalog is a BED with TRGT-style `ID=..;MOTIFS=..` in column 4, or just a motif list
(`CAG,CAA`). CRAM is decoded with `--genome` (`--reference` is the same option, as the other
tools name it). `--karyotype` and `--platform` take any case (`XY`, `ont`). A `.vcf.gz` output is
bgzipped and gets a tabix index beside it (`.tbi`; `.csi` for contigs over 2^29 bases).
`--regions BED` (lrpipe-tr's own option, which `lrpipe run --regions` passes) genotypes only the
catalog's loci that a region reaches (shares a base with), each whole and from all its reads, so
its record is the one a run of the whole catalog writes.

## Method

1. **Reads.** Loci are grouped into chunks of nearby loci; each thread queries a chunk's
   reads once (primary and supplementary alignments, MAPQ ≥ 1, no duplicates) and its
   reference window.
2. **Repeat sequence per read.** A read counts when it aligns at least 50 bases past both
   boundaries. The boundaries are mapped into the read through the CIGAR, an insertion right
   at a boundary counting as repeat. Where the alignment has an indel within 10 bases of a
   boundary on the flank side, the boundary is found again by fitting 30 flank bases into the
   read (edit distance, ties broken toward the alignment's boundary). A boundary two catalog
   loci share (ExpansionHunter's catalogs split FXN, ATXN7, NOP56 and others into adjacent
   repeats) is not refined: an insertion at it, or within 10 bases and sliding to it along its
   own period, is split once between the two loci, where it best continues each one's repeat,
   instead of counting in both. An expansion longer than the aligner aligns through comes as
   an alignment that spans one flank and stops in the repeat, the rest of the read
   soft-clipped (most reads of a 1 kb expansion or more). Such a read counts when its clipped
   bases start as the repeat (they repeat with a motif's period at three in four positions)
   and hold the other flank: 100 reference bases fitted at most 20 edits away, next to the
   repeat or 100 or 200 bases further out (an expanded copy's flank often differs right next
   to the repeat), the allele running up to that fit. At most 200 spanning reads per locus
   are used, spread evenly over them (clipped reads come after the others in read order).
3. **Alleles.** With haplotags (`HP`) on most reads, the two haplotypes, when their reads back
   them: each tagged read is aligned to both haplotypes' consensuses and prefers the closer;
   the tags stand when both haplotypes' reads lean to their own allele and, net, 30% of the
   reads do. Otherwise the two consensuses differ only by read noise and the call is
   homozygous. Untagged reads first go through an expansion test: some reads much longer
   than the rest (the widest gap between consecutive lengths at least doubling the length
   below it and adding 100 bases), however few (3 reads and 5%; an expansion's long molecules
   are sampled less than the other allele's), are the second allele. Then two tests on their
   lengths. The best two-group split of the sorted lengths is kept when it fits
   enough better than one group: `n ln((SS1 + n v) / (SS2 + n v)) > 5 ln n`, with per-read
   length variance `v = 0.5 + (0.006 x length)^2` for HiFi (`8 + (0.03 x length)^2` for
   `--platform ont`), each group holding 20% of the reads. Failing that, two exact lengths each
   holding 30% of the reads (50% for ONT) (the second at least half the first)
   are two alleles: single-base differences in homopolymers, where read noise spreads
   symmetrically around one length. Other reads count for the nearest allele.
   `--haplotags LIST` (lrpipe's own option) takes the haplotypes from a haplotag list
   (WhatsHap's `haplotag --output-haplotag-list`) in place of the `HP` tags, as `haplotag`
   would tag the reads (a listed read's primary alignment, no other), which `lrpipe run` passes
   from its own phased SNVs. On HG002 ONT chr20 (21,959 loci, GIAB TR Tier 1, as below) the
   haplotags count: with the CRAM's `HP` tags 84.4% of the called regions' alleles are exactly
   right, with the tags removed 69.2%.
4. **Allele sequence.** A sequence most reads share; otherwise a consensus: up to 30 reads
   (those nearest the median length: an expansion's reads vary by hundreds of bases from cell
   to cell) aligned to a backbone (banded edit distance, gaps placed consistently), each base,
   deletion and insertion voted, twice. Where the alignment has an indel near a boundary, the
   flank search spans it, so expansions an aligner places just inside the flank are kept.
   With `--platform ont`, an allele whose longest single-base run has 10 bases or more (and
   at least a tenth of the allele) takes the length its reads make most likely, and that run
   is lengthened or shortened to it: Nanopore reads miss bases in long runs more often than
   they add them, so the median read is often a base short. A read's length error `d` is
   modelled as `exp(d / (0.5 + 0.04 h))` short and `exp(-d / (0.5 + 0.02 h))` long, for a run of
   `h` bases, with 5% of reads erring at random (fitted on HG002 chr1).
5. **Output.** VCF with REF and ALT as the repeat plus the base before it, and per sample
   `GT:AL:ALLR:SD:MC` (lengths, supporting read length range, supporting reads, copies of each
   motif).

### Pseudoautosomal regions

In XY samples chrX has one copy outside its pseudoautosomal regions (PARs) and two inside.
The PARs are read off the reference by chrX's length: GRCh38's, GRCh37's or T2T-CHM13's. An
XY sample on a reference with another chrX has no known PARs, so all of its chrX is one copy
(with a warning), as lrpipe-cnv counts it; this was an error, which failed `lrpipe run` late
(lrpipe#44). They used to
be GRCh38's whatever the reference (review finding LRP-010). A `.gz` catalog may be gzip or
BGZF (LRP-014).

## Layout

- `crates/lrpipe-tr`: the library (catalog, reads, extraction, genotyping, VCF, the parallel
  run) and the `lrpipe-tr` command. `tests/end_to_end.rs` calls a heterozygous expansion from a
  small BAM it writes.
- `crates/lrpipe-tr-bench`: scores calls (lrpipe-tr's or TRGT's VCF) against the GIAB HG002 tandem
  repeat benchmark v1.0.1 by allele length change per Tier 1 region.

## Pathogenic expansions (Coriell cell lines, PacBio PureTarget)

PacBio's public PureTarget repeat expansion data (`downloads.pacbcloud.com/public/dataset/PureTargetRE`,
the Revio SPRQ "Coriell24" and Sequel IIe "Coriell" sets): 45 HiFi BAMs of 33 Coriell cell lines,
each holding the reads of its expanded locus, with PacBio's TRGT calls beside them. On the regions
TRGT's catalog gives (lrpipe's STRchive regions of AR and PABPN1 run 34 and 13 bases further),
each allele as its length over the motif's, against the repeat counts in Coriell's catalog
(25 lines give one; within max(3 units, 10%)):

| | Coriell counts matched | agrees with TRGT (max(10 bp, 5%)) |
|---|---|---|
| lrpipe-tr | 24 of 38 samples | 41 of 45 |
| TRGT (PacBio's calls) | 22 of 38 | |
| lrpipe-tr 0.1.7 | 11 of 34 (no call at 4) | 32 of 41 |

Up to 0.1.7 lrpipe-tr missed every expansion of more than about 1 kb (DMPK, FMR1 full
mutations, FXN, RFC1): their reads come soft-clipped (above). It now calls them all, e.g. DMPK
439 CTG (TRGT 439), FMR1 NA09237 884-923 CGG (Coriell 931-940, TRGT 888-913), FXN NA15850 627
and 1,037 GAA (Coriell 650 and 1,030; TRGT 604 and 738). Most of both tools' misses are the
catalog's: Coriell counts HTT's and ATN1's CAG alone, the regions include the adjacent units;
and some Coriell counts are old (NA06905 FMR1 70 CGG, both tools 79-80; NA16215 FXN 830 GAA,
both about 1,000). The C9orf72 and RFC1 lines have no Coriell size; there both tools agree.

On whole genomes (HG002 against GIAB TR v1.0.1 Tier 1, the pipeline's haplotags) the same
share of regions is exactly right as before (HiFi 94.24%, ONT 85.79%), at the same speed, and
two expansions of 2.5 and 10.5 kb (chr4:189.8 Mb, chr5:178.6 Mb) are now called as GIAB has them.

## Results (HG002 HiFi, 31x, GRCh38, 16 threads)

Speed, on 937,204 loci (adotto + STRchive), against TRGT 5.1.0 on the same BAM:

| | Wall | CPU | Peak memory |
|---|---|---|---|
| lrpipe-tr | 78 s | 895 s | 3.3 GB |
| TRGT 5.1.0 | 645 s | 7,820 s | 5.2 GB |

(Measured on 2026-11-03 with lrpipe 0.1.3 plus the review fixes; [BENCHMARKS.md](../../BENCHMARKS.md)
has every tool's run. Before reads were decoded with pasta, lrpipe-tr took 161 s.)

The 82 STRchive pathogenic loci alone take 0.3 s; 80 of them agree with TRGT's allele lengths
within max(1 bp, 5%).

Accuracy against GIAB HG002 TR v1.0.1 Tier 1 (`lrpipe-tr-bench`): both alleles' length change
exactly right, by the truth's largest change. A region is scored when the caller made a call
there; the catalog covers only part of Tier 1, so both callers call about a quarter of its
1,638,106 regions, and lrpipe-tr slightly more of them than TRGT:

| | Called | Not called | Ambiguous (several calls) |
|---|---|---|---|
| lrpipe-tr | 393,498 (24.02%) | 1,240,074 | 4,534 |
| TRGT 5.1.0 | 392,120 (23.94%) | 1,240,126 | 5,860 |

On the 391,984 regions both call, so that neither is scored where the other made no call
(`lrpipe-tr-bench --common`):

| Truth change | Regions | lrpipe-tr | TRGT 5.1.0 | lrpipe-tr ±1 bp | TRGT ±1 bp |
|---|---|---|---|---|---|
| none | 202,367 | 96.5% | 91.7% | 99.2% | 98.7% |
| 1-4 bp | 110,246 | 83.3% | 82.0% | 94.5% | 94.0% |
| 5-49 bp | 70,974 | 87.0% | 85.3% | 91.3% | 90.4% |
| >= 50 bp | 8,397 | 83.3% | 81.3% | 86.8% | 85.8% |
| all | 391,984 | 90.8% | 87.6% | 96.2% | 95.6% |

Each on all the regions it calls, lrpipe-tr is exactly right at 90.6% (21.8% of all Tier 1
regions) and TRGT at 87.6% (21.0%). The genotyping parameters were tuned on chr1; on the other
chromosomes alone lrpipe-tr is right at 90.8% of regions and TRGT at 87.6%.

HG002 ONT (haplotagged CRAM, its `HP` tags, same catalog and scoring, `--platform ont`): 235 s
and 5.9 GB on 16 threads, mostly CRAM decoding. The ONT settings were tuned on chr1; on the
other chromosomes alone 86.0% of regions are exactly right (85.1% before the length model).

| Truth change | ONT | ONT ±1 bp | before the length model |
|---|---|---|---|
| none | 93.2% | 98.7% | 92.3% |
| 1-4 bp | 72.7% | 92.7% | 71.5% |
| 5-49 bp | 85.7% | 91.7% | 85.4% |
| >= 50 bp | 86.2% | 90.8% | 86.3% |
| all | 85.9% | 95.6% | 85.0% |

Most remaining ONT errors are still single bases in homopolymers: a 1-base difference
between the haplotypes in a run of 10 or more, which read lengths cannot confirm (a likelihood
test on them does worse than the reads' preference for their own haplotype's allele), and
minisatellites whose units hold several runs, where the median read is already off.

Long runs are near what the reads allow. At chr1's 1,725 homozygous homopolymer regions (median
48 reads), the share of reads of the true length falls from 55% for runs under 15 bases to 22%
at 15-19, 17% at 20-24 and 12% at 25 or more, and the errors lean short (mean -0.6 base at
20-24, -1.3 at 25 or more). The length model gets 75.3% of these regions exactly right from
the reads' lengths. A table of the reads' observed errors per 3-base run length, fitted on half
the regions and scored on the other half, gets 77.6%, and fitted on the very regions it scores
78.3%, the most the read lengths support: about 0.2% of all regions, for a table fitted to one
basecaller on one sample, so it is not used. Refitting the model's widths by likelihood (74.7%)
and moving its peak short of the allele's length (67.1% and lower) do worse. Across the
HG002, HG003, HG004 trio, Tier 1 homopolymers are 42.9% exactly consistent against 88.8% for
2-6 base motifs, and on chrY, which a son inherits whole, 32% of his homopolymers equal his
father's (92% of 2-6 base motifs), nearly all of the rest 1 base apart.

```bash
lrpipe-tr-bench calls.vcf.gz HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz
```

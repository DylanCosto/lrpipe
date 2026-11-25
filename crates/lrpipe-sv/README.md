# lrpipe-sv

An unofficial Rust port of [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
(structural variant calling from long reads) by Moritz Smolka, Hermann Romanek, Fritz
Sedlazeck and contributors. The first goal is **VCF records identical to Sniffles 2.8.1**,
just faster; accuracy changes come after that, behind a flag until they beat Sniffles2 on
GIAB. It is not affiliated with or endorsed by the Sniffles project; Sniffles2's license is
in `LICENSE-SNIFFLES`, lrpipe-sv's own (MIT) in `LICENSE`, and the 2.8.1 source it is ported
from is in `upstream/`.

It is part of [lrpipe](..) (`lrpipe sv`), and named so as not to suggest an official Sniffles
project (it was sniffles-rs while that was a working name).

## Status

`lrpipe-sv` calls SVs from one sorted, indexed BAM or CRAM into a VCF, writes `.snf`
candidate files, and calls several samples together from `.snf` files. The VCF records are
identical to Sniffles 2.8.1's (only the `##source`, `##command` and `##fileDate` header lines
differ) in every comparison run so far, but for the records of the Sniffles2 bugs lrpipe fixes
([below](#sniffles2-bugs-fixed-on-purpose)):

- HG003 chr20 (ONT R10.4.1 sup, minimap2): default options, no `--reference`, `--symbolic`,
  `--no-qc` (11,038 records), `--mosaic`, `--minsupport auto`, `--tandem-repeats`,
  `--regions` and `--region`, CRAM input
- HG003 chr20 `--genotype-vcf`: its own calls, and GIAB v5.0q chr20 (208,757 records)
- HG002 chr20 (ONT, bonefish): default options and `--no-qc` (15,758 records)
- HG003 whole genome (ONT sup, 110 GB BAM): 27,046 records; HG002 whole genome (HiFi):
  26,287 records; HG002 whole genome (ONT R10.4.1 sup, haplotagged 99 GB CRAM from ONT's
  giab_2025.01 release): 27,420 records
- Multi-sample calling of 8, 32 and 128 Sniffles2 `.snf` files (HG003 chr20 subsamples)

`--snf` files hold the same candidates and coverage as Sniffles2's (checked by decoding both),
and Sniffles2 combines them exactly as it combines its own.

Output is the same for any `--threads`. `--vcf out.vcf.gz` writes bgzipped output with a tabix
index, as Sniffles2 does.

Written entirely in Rust: BAM and CRAM reading use [pasta](https://github.com/DylanCosto/pasta)
(begun from noodles) with the pure-Rust zlib-rs inflate; `.snf` files are read and written with a small pickle
codec.

Not ported: `--combine-population` and population SNFs, `--combine-consensus`, `--dev-*`
debugging options.

Known differences from Sniffles2:

- `RNAMES` (`--output-rnames`) lists read names in first-seen order. Sniffles2 lists them
  in Python `set` order, which changes from run to run.
- Unmapped records and records without a CIGAR are skipped.
- `##fileDate` is in UTC; Sniffles2 writes local time.
- Multi-sample calling with more than 20 inputs: Sniffles2 2.8.1 sets calls that arrive out of
  position order aside in `result-*-unsorted.part.vcf` files in its temp directory and leaves
  them out of the VCF (105 of 1,025 calls for 32 samples). lrpipe-sv keeps them and sorts the
  output; the records are otherwise identical.
- Multi-sample calling runs each contig as one task, which is what Sniffles2 does with one
  thread. With several threads Sniffles2 splits long contigs into pieces and its output depends
  on the thread count.
- Multi-sample calling from more `.snf` files than stay open under the open-file limit
  (`ulimit -n`) works without `--combine-close-handles`: each file is then opened again for
  every block read, as the flag does here and in Sniffles2, which fails without it.

### Sniffles2 bugs fixed on purpose

lrpipe fixes real bugs of the tools it reproduces rather than copying them; the VCF differs
from Sniffles2's only in these records (review findings LRP-004, 007, 008, 009, 015 and 016, and
LSV-004 of the pre-release bug hunt):

- **Breakend mates.** With a reference, Sniffles2 replaces IUPAC codes with `N` in the whole ALT
  of a breakend, mate included, so a mate on chrY or chrM became `chrN`, and one on an unplaced
  contig such as `chrUn_KI270746v1` became `chrUn_NI270746v1`, while `CHR2` still named it.
  lrpipe-sv replaces them in the bases only. On the whole genomes this changes the BND records
  with such a mate and nothing else: 29 of HG003 ONT's 27,046 records, 23 of HG002 ONT's 27,420
  and 33 of HG002 HiFi's 26,287 ([BENCHMARKS.md](../../BENCHMARKS.md)).
- **Breakends when combining.** Combining `.snf` files screened every SV by length, and a BND's
  is 0, so Sniffles2 dropped every breakend from multi-sample VCFs. As in single-sample calling,
  breakends are no longer screened by length.
- **Breakends genotyped by their mate and orientation.** `--genotype-vcf` gave a breakend the
  support of the called breakend nearest it with a mate on the same contig, however far that
  mate or whichever way it joined: a target `N[chrY:900001[` took the reads of `N[chrY:10001[`.
  A target now takes a called breakend's support only with the same mate contig, the mate
  within `--cluster-merge-bnd` and the same orientation on both sides; its orientation comes from
  the ALT's brackets, as Sniffles2 read it from a leading `N`, which a reference base replaces.
  On the parity fixture this changes one record: sample2's breakend at chr17:46292138 had the
  1 read of another candidate (0/0), and now has its own 8 (0/1), as calling sample2 finds it.
- **Breakends combined by orientation.** Combining `.snf` files grouped breakends by position
  and mate contig, so `N[chrY:10001[` of one sample and `N]chrY:10001]` of another, different
  adjacencies, became one allele of the cohort. Breakends now combine only with the same
  orientation on both sides.
- **`--genotype-vcf` output that describes itself.** Every genotype ends with its phase set,
  but Sniffles2's FORMAT is `GT:GQ:DR:DV`, which bcftools rejects; it is `GT:GQ:DR:DV:PS` now,
  with PS declared. The `#CHROM` line names one sample, the genotyped one (`--sample-id`), as
  every record has one: Sniffles2 copies the input's, which for a sites-only VCF has no sample
  column and for a multi-sample one more than the records have. A sites-only record whose INFO
  ends with SVTYPE, as catalogs often write it (`...;SVTYPE=INS`, `SVTYPE=BND`), is genotyped:
  Sniffles2 read its type with the line's newline, `INS\n`, skipped it as unsupported and wrote
  it 0/0 whatever its reads.
- **CSI indexes.** A BAM indexed with a `.csi` and no `.bai` is read (Sniffles2 reads it through
  pysam; lrpipe-sv's reader took only a BAI).

The parity tests compare with Sniffles2's output after exactly these changes.

## Speed

16 threads on a 16-core machine; all lrpipe-sv outputs below are identical to Sniffles2's (measured
before the breakend fix above, which changes a few dozen BND records per genome).

| input (ONT R10.4.1 sup, BAM on NVMe) | Sniffles 2.8.1 | lrpipe-sv |
|---|---:|---:|
| HG003 whole genome (110 GB BAM) | 126.3 s, 2.00 GB | 77.5 s, 0.38 GB |
| HG003 chr1-chr8 (53 GB BAM) | 84.8 s, 1.95 GB | 38.1 s, 0.30 GB |
| HG003 chr20 (2.5 GB BAM) | 24.4 s, 0.83 GB | 2.0 s, 0.06 GB |

Most of the run is spent inflating BGZF blocks (about two thirds of lrpipe-sv's CPU time,
with the pure-Rust zlib-rs). Sniffles2's own time is mostly in htslib's decoder too, so the
whole-genome gain is smaller than on single chromosomes, where Sniffles2 uses one core.

Sniffles2 makes one task per contig, so a single chromosome runs on one core. lrpipe-sv
reads each contig in 2 Mbp pieces in parallel and merges them so that the result does not
depend on how the contig was split. Contigs go through a two-stage pipeline (one is called
while the next is read), so memory stays at about two contigs' worth of SV evidence.

Multi-sample calling (HG003 chr20 `.snf` files; Sniffles2's best time of 1 and 16 threads):

| samples | Sniffles 2.8.1 | lrpipe-sv, 16 threads |
|---:|---:|---:|
| 8 | 5.4 s | 0.14 s |
| 32 | 6.2 s | 0.53 s |
| 128 | 10.9 s | 1.58 s |

CRAM input (16 threads): HG003 chr20 Sniffles2 45.6 s, lrpipe-sv 13.3 s (BAM: 2.0 s); HG002
ONT whole genome (99 GB) Sniffles2 573 s, lrpipe-sv 294 s. CRAM decoding is the slow part.
For CRAM records without a stored `NM` tag, lrpipe-sv computes it from the reference, as
htslib does when Sniffles2 reads them (it changes which phased calls are rescued).

## Comparison with sawfish (HiFi)

HG002 chr20, PacBio Revio HiFi reads aligned to GRCh38 with soft-clipped supplementary
alignments (`-Y`, which sawfish requires), 16 threads. Accuracy is Truvari 5.5 against GIAB
v5.0q (chr20 benchmark regions, PASS calls, SVs >= 50 bp). "strict" uses `--pctseq 0`,
"multi" adds `--pick multi`; with the default `--pctseq 0.7` every caller's matches were the
same as in "strict". sawfish 2.3.1 ran `discover` + `joint-call` with default options (no
CNV excluded-regions file).

| aligner | caller | wall | peak memory | F1 strict (precision, recall) | F1 multi |
|---|---|---:|---:|---|---:|
| minimap2 2.31 `map-hifi` | Sniffles 2.8.1 | 14.5 s | 0.81 GB | 0.811 (0.938, 0.715) | 0.892 |
| | lrpipe-sv | 1.6 s | 0.05 GB | 0.811 (identical VCF) | 0.892 |
| | sawfish 2.3.1 | 25.2 s | 1.90 GB | 0.820 (0.895, 0.756) | 0.864 |
| bonefish `-x hifi` | Sniffles 2.8.1 | 14.3 s | 0.81 GB | 0.810 (0.932, 0.716) | 0.891 |
| | lrpipe-sv | 1.1 s | 0.05 GB | 0.810 (identical VCF) | 0.891 |
| | sawfish 2.3.1 | 16.3 s | 1.86 GB | 0.819 (0.896, 0.754) | 0.864 |

sawfish assembles each candidate locus, which finds more of the truth SVs (recall 0.756 vs
0.715) at lower precision; with `--pick multi` (one call may match several truth records, as
for SVs split differently in the truth) Sniffles2's calls score higher. lrpipe-sv output is
identical to Sniffles2's on these HiFi BAMs as well.

### Whole genome

HG002, 6.53 M Revio HiFi reads aligned with bonefish `-x hifi -Y` (34.9 min, 16 threads), scored
genome-wide in the GIAB v5.0q benchmark regions (28,123 truth SVs), otherwise as above.

| caller | wall | peak memory | F1 strict (precision, recall) | F1 multi |
|---|---:|---:|---|---:|
| Sniffles 2.8.1 | 72.6 s | 1.71 GB | 0.806 (0.926, 0.714) | 0.870 |
| lrpipe-sv | 31.2 s | 0.38 GB | 0.806 (identical VCF, 26,287 records) | 0.870 |
| sawfish 2.3.1 | 410.1 s | 4.25 GB | 0.802 (0.893, 0.728) | 0.847 |

Genome-wide, sawfish's extra recall (0.728 vs 0.714) does not make up for its extra false
positives (2,463 vs 1,611), so Sniffles2's calls score slightly higher in both modes.

## Accuracy improvements (opt-in)

Options under "lrpipe-sv improvements" in `--help` change the output; without them the
output is identical to Sniffles2's. For the best accuracy use the preset:

```bash
lrpipe-sv --input sample.bam --vcf sample.vcf --reference ref.fa --threads 16 --accuracy
```

`--accuracy` turns on `--split-alleles` and `--read-support-genotypes` and sets
`--cluster-merge-pos` to 50 (an explicit `--cluster-merge-pos` still applies). The sections
below describe and measure each.

Whole genome, HG002 ONT R10.4.1 sup (ONT's haplotagged giab_2025.01 CRAM), Truvari 5.5
against GIAB v5.0q in its benchmark regions (`--passonly --pctseq 0`; "multi" adds
`--pick multi`); GT is genotype concordance of strict matches:

| caller | strict P | strict R | strict F1 | multi F1 | GT |
|---|---:|---:|---:|---:|---:|
| Sniffles 2.8.1 (= lrpipe-sv default) | 0.917 | 0.725 | 0.810 | 0.875 | 0.834 |
| lrpipe-sv `--accuracy` | 0.898 | 0.782 | **0.836** | **0.890** | **0.896** |
| ONT wf-human-variation calls (Sniffles 2.0.7, `--minsvlen 30`, tandem repeats) | 0.910 | 0.669 | 0.772 | 0.833 | 0.809 |

### `--split-alleles`

Genome-wide on HG002 HiFi, 81% of Sniffles2's false negatives against GIAB v5.0q are in
tandem repeats of 101 bp or more, and 71% have another truth SV of the same type within
500 bp: the two haplotypes of a repeat expansion with different lengths (for example +300
and +330 bp). Sniffles2 re-splits a cluster by length only when lengths differ by more than
max(50 bp, 22%), so it reports one call and the second allele is missed.

`--split-alleles` looks for the best two-way split of each INS/DEL cluster's lead lengths and
reports two calls (each genotyped from its own reads) when each side has at least 3 reads and
20% of the cluster, and the empty gap between the two groups' lengths is at least
`--split-min-gap-sd` (default 1) times the larger within-group standard deviation (at least 1%
of the length).

HG002, Truvari against GIAB v5.0q as above:

| data | F1 strict, default (= Sniffles2) | F1 strict, `--split-alleles` | F1 multi, default | F1 multi, `--split-alleles` |
|---|---:|---:|---:|---:|
| chr20 HiFi (minimap2 `-Y`) | 0.811 | 0.843 | 0.892 | 0.894 |
| chr20 ONT R10.4.1 (bonefish) | 0.814 | 0.847 | 0.896 | 0.898 |
| whole genome HiFi (bonefish `-Y`) | 0.806 | 0.828 | 0.870 | 0.874 |

Whole genome, `--split-alleles` finds 1,139 more truth SVs (recall 0.714 -> 0.754) for 326
more false positives (precision 0.926 -> 0.916). For comparison, sawfish scores 0.802 strict
and 0.847 multi on the same data.

### `--cluster-merge-pos 50` (an existing Sniffles2 option)

Sniffles2 joins two indels of the same read into one lead when they are within
`--cluster-merge-pos` (default 150 bp) of each other on both the read and the reference. In
tandem repeats this often joins pieces that belong to different events, so the call ends up
larger than either. Lowering it to 50 bp helps on its own and together with
`--split-alleles` (values from 0 to 50 bp score the same; 100 bp is worse):

| data | default | `--cluster-merge-pos 50` | `--split-alleles` | both |
|---|---:|---:|---:|---:|
| chr20 HiFi, F1 strict / multi | 0.811 / 0.892 | 0.825 / 0.909 | 0.843 / 0.894 | 0.856 / 0.912 |
| chr20 ONT, F1 strict / multi | 0.814 / 0.896 | 0.827 / 0.914 | 0.847 / 0.898 | 0.864 / 0.920 |
| whole genome HiFi, F1 strict / multi | 0.806 / 0.870 | 0.814 / 0.880 | 0.828 / 0.874 | 0.837 / 0.884 |

`--accuracy` (both settings) is the recommended choice. Whole genome on HiFi it
finds 1,528 more truth SVs than Sniffles2 (recall 0.714 -> 0.768) for 294 more false positives
(precision 0.926 -> 0.919), at the same speed.

Remaining false negatives are mostly not fixable in the caller: of those with a passing call of
the wrong size nearby (the largest group), 78% have no read whose insertion or deletion length
is within 30% of the truth, because the assembly-based truth represents the repeat change
differently (for example as several smaller records). `tools/wrongsize_reads.py` checks this.

Tried and not adopted: Sniffles2's `--tandem-repeats` mode (merges neighbouring clusters in
annotated repeats) lowers strict F1 to 0.789 on chr20 HiFi, and using the annotation only to
join each read's indels does not help.

Base qualities are not used, as in Sniffles2. On the ONT whole genome with `--accuracy`, the
mean base quality of the supporting reads' inserted bases (or of the bases flanking a deletion)
barely separates true from false calls (AUC 0.556 for INS, 0.522 for DEL; about 0.5 among calls
with the same support), and every quality threshold tried lowered multi F1 (best -0.0016).

### `--read-support-genotypes`

Sniffles2 genotypes an insertion of `--long-ins-length` (2,500 bp) or more from its support
scaled up by its length, 1.66 + 0.33 x length / 2,500 (about 2x at 2.7 kb), to make up for
reads that end inside it. But it already counts those reads among the supporters
(`SUPPORT_LONG`), and long reads mostly span such insertions, so the scaled support comes close
to the whole coverage and a heterozygous insertion is called 1/1. On HG002 ONT, a 2,748 bp
insertion GIAB has as 0/1 has 18 supporting reads of 33, scaled to 36: `DV=36`, `DR=0`, `1/1`.
In the HG002, HG003, HG004 trio this was the largest source of Mendelian errors among long
insertions: a heterozygous insertion called 1/1 in a parent or in the son.

`--read-support-genotypes` computes GT, GQ, DR and DV of these insertions from the reads that
support them. `VAF` and every filter keep the scaled support, so the same calls pass with the
same positions, lengths and filters; only the genotype and the `PHASE` field's homozygous flag
change. Whole genome, Truvari against GIAB v5.0q through `lrpipe run`'s phasing:

| data | `--accuracy` before | with `--read-support-genotypes` |
|---|---:|---:|
| HG002 HiFi, GT concordance (strict) | 0.8946 | 0.9000 |
| HG002 ONT, GT concordance (strict) | 0.8876 | 0.8958 |

Precision, recall and F1 are unchanged, call for call. On the ONT trio (`lrpipe sv --accuracy`
on each sample), Mendelian errors among PASS insertions of 1 kb or more fall from 15.1% to 6.6%
(all SVs: 5.29% to 4.78%).

### Per-technology tuning (not adopted)

`--accuracy` is one preset for both HiFi and ONT. To check whether separate presets would pay
off, each `--split-*` option and `--cluster-merge-pos` was varied one at a time around
`--accuracy` on whole-genome HG002 (HiFi bonefish `-Y` BAM and the ONT CRAM above), choosing on
the odd chromosomes and X and scoring on the even chromosomes and Y. The best combination per
technology, as the change from `--accuracy` on the validation half:

| tech | options added to `--accuracy` | strict P | strict R | strict F1 | multi F1 | GT |
|---|---|---:|---:|---:|---:|---:|
| HiFi | `--split-min-support 4 --split-min-gap-sd 0.75 --cluster-merge-pos 25` | -0.0029 | +0.0054 | +0.0020 | +0.0019 | +0.0013 |
| ONT | `--split-min-frac 0.25 --split-min-gap-sd 0.75` | -0.0017 | +0.0027 | +0.0008 | -0.0002 | +0.0033 |

HiFi gains on both halves but only by 0.002; for ONT the strict gain on the tuning half (+0.0020)
mostly disappears on validation. Every setting that helped traded precision for recall. Gains
this small, measured on one sample and one aligner per technology, don't justify two more
presets, so `--accuracy` stays the single recommendation. Single-option results, for reference
(validation-half strict F1 change): `--split-min-gap-sd 0.75` +0.0007 HiFi / -0.0010 ONT,
`--split-min-gap-sd 1.5` -0.0012 / -0.0011, `--split-noise 0.005` or `0.02` within ±0.0014,
`--cluster-merge-pos 100` -0.0043 / -0.0027.

## Usage

```bash
lrpipe-sv --input sample.bam --vcf sample.vcf --reference ref.fa --threads 16
# candidates for multi-sample calling (with or without --vcf)
lrpipe-sv --input sample.bam --snf sample.snf --reference ref.fa --threads 16
# multi-sample calling: .snf files, or a .tsv of "file.snf[<TAB>sample id]" lines
lrpipe-sv --input a.snf b.snf c.snf --vcf multi.vcf --reference ref.fa
# genotype known SVs
lrpipe-sv --input sample.bam --genotype-vcf known.vcf --vcf genotyped.vcf --reference ref.fa
```

Options use Sniffles2's names and defaults; see `lrpipe-sv --help`. `.snf` files from
Sniffles2 (2.5.3 or later) and from lrpipe-sv can be mixed.

`--reference` is a FASTA with its `.fai`: plain, or compressed with bgzip with its `.gzi` beside
it too (`ref.fa.gz`, as references are distributed, and as Sniffles2 reads them through pysam),
told apart by contents. Every part of the pass reads it so: the SV calls, CRAM decoding, depth,
copy numbers, repeats and bedMethyl, with the same output as from the plain FASTA. A FASTA
compressed with gzip alone is an error that says to recompress it with `bgzip`.

### Read depth in the same pass

Depth, tandem repeats (`--repeats`) and methylation (`--bedmethyl`) come from the same pass as
the SV calls, and also without them: given none of `--vcf`, `--snf` and `--cnv`, the pass reads
no SV evidence and calls nothing (on HG002 HiFi chr20 with depth and repeats, the other outputs
were the same either way). Sniffles2 needs `--vcf` or `--snf`.

`--depth PREFIX` also writes read depth, counted from the reads the SV calling already decodes,
with output identical to `mosdepth -n` (through lrpipe-depth):
`PREFIX.mosdepth.summary.txt` and `.mosdepth.global.dist.txt`, and with `--depth-by`
`PREFIX.regions.bed.gz` and `.mosdepth.region.dist.txt`.

```bash
# the options mosdepth -x -Q 20 -b 1000 takes, for CNV calling
lrpipe-sv --input sample.bam --vcf sample.vcf --threads 16 \
    --depth sample --depth-by 1000 --depth-fast --depth-mapq 20
```

| Option | mosdepth | |
|---|---|---|
| `--depth-by` | `-b` | window size or BED file |
| `--depth-mapq` | `-Q` | minimum mapping quality (default 0) |
| `--depth-flag` | `-F` | flag bits to exclude (default 1796) |
| `--depth-fast` | `-x` | count whole read spans |
| `--depth-median` | `-m` | region medians instead of means |
| `--depth-thresholds` | `-T` | bases at or above each depth, per region |

Depth covers every contig in the input, or only the `--contig` ones; `--regions` and the 1 Mbp
contig minimum do not restrict it (contigs the SV pass does not read whole are counted
separately), but for `--depth-in-regions` (lrpipe's own option, which `lrpipe run --regions`
passes): the depth is then counted from the reads of `--regions` alone (sorted and apart), in
the pass, each read once, for the `--depth-by` BED's regions, which must lie in them. Their rows
and region distribution are mosdepth's from the whole contig, but the summary has only the
regions' rows (`CONTIG_region`, `total_region`) and no global distribution is written, as the
contigs' depth is not counted (a panel of 20 regions of 2 kb on chr20 of the HG003 ONT BAM: 0.39 s
against 2.50 s reading chr20 whole). The SV calls are unchanged. On HG002 HiFi (39 GB BAM,
31×, 16 threads) the run takes 37.3 s instead of 33.6 s and 2.2 GB of memory instead of 0.4 GB; running lrpipe-depth
afterwards would take 28.5 s and 10.1 GB, mosdepth 44.3 s.

### Tandem repeats and copy number in the same pass

The same pass can also genotype tandem repeats ([lrpipe-tr](../lrpipe-tr))
and call copy number variants ([lrpipe-cnv](../lrpipe-cnv)), so one
read of the alignments gives SVs, depth, repeats and CNVs:

```bash
lrpipe-sv --input sample.bam --vcf sample.sv.vcf.gz --reference ref.fa --threads 16 \
    --depth sample --depth-by 1000 --depth-fast --depth-mapq 20 \
    --repeats catalog.bed --repeats-vcf sample.repeats.vcf.gz --karyotype xy \
    --cnv sample.cnv --cnv-blacklist blacklist.bed --cnv-snv sample.lrpipe-snv.vcf.gz
```

- `--repeats CATALOG --repeats-vcf OUT`: each read the SV calling decodes is offered to lrpipe-tr,
  which records its sequence at every catalog locus it spans; a contig's loci are genotyped in
  parallel with its SV calls. Contigs the pass does not read whole are genotyped by lrpipe-tr's
  own reader. The VCF is the one `lrpipe-tr` writes from the same reads.
  `--repeats-platform hifi|ont` picks lrpipe-tr's noise model, `--karyotype xx|xy` the copies of
  chrX and chrY. `--repeats-regions BED` (lrpipe's own option, as `lrpipe tr --regions`)
  genotypes only the loci a region reaches, each whole, from all its reads.
- `--cnv DIR`: lrpipe-cnv calls from the `--depth` windows, which must be Spectre's (the
  `--depth` options above; anything else is refused). `--cnv-mode improved` (the default) is
  lrpipe-cnv's `--improved` with breakpoints from this run's SV calls and valid VCF 4.2;
  `--cnv-mode spectre` is Spectre 0.2.1's method exactly. `--cnv-snv` adds copy-neutral depth
  from heterozygous SNVs and LoH calls; `--cnv-metadata` reuses a reference `.mdr`. With
  `--contig`, only those contigs are called (those of them longer than 1 Mbp, as lrpipe-cnv
  calls no shorter ones).

`tools/joined_vs_separate.sh` runs the four tools separately and in one pass and checks that
every output is the same. On HG002 HiFi (31x, 16 threads) all of them are (SV VCF, depth
windows and summary, repeat VCF over 937,204 loci, CNV VCF and BED):

| | Wall | Peak memory |
|---|---|---|
| lrpipe-sv (SVs) | 32.5 s | 0.4 GB |
| lrpipe-depth | 48.3 s | 9.6 GB |
| lrpipe-tr | 162.6 s | 3.4 GB |
| lrpipe-cnv | 26.2 s | 2.2 GB |
| separately, in total | 269.5 s | |
| one pass | 118.5 s | 5.2 GB |

On HG002 ONT (haplotagged CRAM, `--repeats-platform ont`) the outputs are the same too, and the
pass takes 713 s against 1,445 s for the four tools (lrpipe-sv 297 s, lrpipe-depth 244 s, lrpipe-tr
879 s, lrpipe-cnv 26 s), at 7.2 GB against lrpipe-depth's 11.9 GB. Most of the pass beyond the SV
calling was lrpipe-tr's genotyping (consensus and haplotag checks on noisy reads), not decoding;
with lrpipe-tr 8daa56a's bit-parallel haplotag check the pass without `--cnv` takes 494 s
(GIAB TR accuracy unchanged).

### Methylation in the same pass

`--bedmethyl OUT` (needs `--reference`) also writes CpG methylation, identical to
`lrpipe methyl pileup INPUT OUT --ref REFERENCE --modified-bases 5mC 5hmC --cpg
--combine-strands` (modkit's CpG recipe; with `--contig`, as that command with `--region` for
each). The pass thresholds are estimated first, from a sample of the reads, as that command
estimates them. Each piece of a contig then counts its reads' calls at the contig's CpG sites
alone, and the contig's table is written window by window, alongside its SV calls, from the
merged counts. Contigs the pass does not read whole are counted by lrpipe-methyl's own reader.
`--bedmethyl-regions BED` (lrpipe's own option) writes only the positions in the regions, with
the pass thresholds estimated from their reads, identical to that command with `--include-bed
BED` (modkit's option); the regions are counted by lrpipe-methyl's own reader. `--regions`
(Sniffles2's) restricts the SV calls alone, not the depth, repeats or methylation.

On HG002 ONT (haplotagged CRAM, 16 threads, with depth and repeats) the bedMethyl is the same as
`lrpipe methyl pileup`'s, and the pass takes 5,366 CPU-seconds against 3,774 without it and
3,338 for `lrpipe methyl pileup` on its own: methylation costs half the CPU when it shares the
decoding. Peak memory goes from 5.7 GB to 6.2 GB.

### The reads' statistics in the same pass

`--read-stats FILE` also writes the reads' statistics as JSON, counted from every record the
pass reads (each once, by its start), with the contigs it skips (those under 1 Mb) read for them
alone: from CRAM in one pass over their slices, each decoded once for its records' flags,
reference and CIGAR, without their reference sequences (on the HG002 ONT CRAM, 170 such
contigs, the pass took 71.4 and 71.7 s against 78.4 and 87.9 s querying each contig with its
bases; on hs38DH's some 3,000 contigs each also kept its reference sequence, which no longer
happens). It gives the primary alignments of mapped reads, duplicates and QC failures included
(`reads`, `bases`, length `n50`, `mean_length`, and `mean_quality`, the arithmetic mean of the
bases' Phred qualities, as samtools stats' "average quality"), and the `secondary`,
`supplementary` and placed unmapped records. With `--contig` or `--regions` they are of those
alone (`restricted`): in regions, of each read that reaches them, as `--depth-in-regions`
counts them (once when the regions are sorted and apart), and no short contig is read. `lrpipe run`'s
QC summary reads it. From CRAM there is no mean quality: the pass does not decode CRAM's
qualities, which would make it 65% slower (HG002 ONT chr20: 3.8 s against 2.3 s). On HG003 ONT
chr20 (a 2.5 GB BAM, 16 threads) it takes no measurable time (a median of 3.04 s against
3.01 s over 12 runs), 2% more CPU and 5 MB more memory.

### Phasing from a haplotag list

The SVs' phases come from the reads' `HP` and `PS` tags, as Sniffles2's do. Two lrpipe-only
options take them from a haplotag list instead (WhatsHap's `haplotag --output-haplotag-list`,
or `lrpipe haplotag --list-only`'s), so reads need not be written again haplotagged:

```bash
# from the list in the read pass
lrpipe-sv --input reads.bam --vcf calls.vcf --haplotags haplotags.tsv.gz
# or calls first, beside the phasing of the SNVs, and phases once the list is made
lrpipe-sv --input reads.bam --vcf unphased.vcf --snf unphased.snf --phasing-evidence evidence.gz
lrpipe-sv phase --evidence evidence.gz --haplotags haplotags.tsv.gz \
  --input-vcf unphased.vcf --vcf calls.vcf --input-snf unphased.snf --snf calls.snf
```

Either way a read is tagged as `haplotag` would tag it: its primary alignment gets the list's
haplotype and phase set, every other alignment none, whatever tags the input carries (with
`--haplotags` the tandem repeats read the list too). The output is then the one from the
alignments `haplotag` writes with the list: on reads split into primary and supplementary
alignments, and with secondary ones, tagged otherwise by another pipeline, Sniffles 2.8.1's
from the haplotagged reads (`tests/data/split*.bam`). A read listed more than once on a contig
with different haplotypes (a list written with `--tag-supplementary`, whose lines do not say
which is the primary alignment's) is left untagged, with a warning, and the log says how many
of the reads the list names, with a warning when it names fewer than half (a list of other
reads, or of another alignment of them). Neither option is taken with `.snf` inputs.

`--phasing-evidence FILE` calls the SVs as if the reads had no `HP` and `PS` tags and writes
what phasing them needs: for each call the reads behind its phase (`phase_sv`), its genotype,
and, for a MOSAIC_VAF call that Sniffles2's phasing rescue may pass as heterozygous, the reads
the rescue counts at its first bin and its VCF record. `lrpipe-sv phase` works out each call's
PHASE, genotype phase and PS, and the rescue, with the read pass's own functions, and rewrites
the read pass's VCF and SNF with them; nothing else in them depends on the haplotypes. On the
parity fixture, phased from its own tags as a list, the VCF is Sniffles 2.8.1's from the tags
and the SNF lrpipe-sv's, and a list that sets off the rescue gives the read pass's output from
the same list (`tests/phasing.rs`). `lrpipe run` phases its SVs this way from its own phased
SNVs. `--output-rnames` is not needed for it and changes nothing.

## Building

```bash
cargo build --release
```

Needs Rust 1.91+ and the system linker (`build-essential` on Debian/Ubuntu, the Xcode command
line tools on macOS). There are no C dependencies, so the same build works on x86-64 and ARM
(aarch64) Linux and on macOS. Windows is not supported (`.snf` reading uses Unix file APIs).

## Testing

```bash
cargo test --release
```

Besides the unit tests, `tests/parity.rs` runs lrpipe-sv on a bundled fixture (`tests/data`:
562 HG002 ONT alignments around chr17:46.3 Mb, 3.8 MB) and requires the output to match what
Sniffles 2.8.1 wrote for the same input, for default options at 1, 2 and 4 threads,
`--reference`, CRAM input, `--phase`, `--mosaic`, `--no-qc`, `--symbolic`, `--regions`,
`--genotype-vcf`, multi-sample calling from `.snf` files and bgzipped output; and that a
bgzipped reference gives what the plain one gives. Sniffles2 is only
needed to regenerate the fixture (`tools/make_fixture.sh`). `tools/check.sh` runs the tests,
also with Rust 1.91.

## Validating on another machine

To check a new platform (a VM, an ARM server, a Mac); step-by-step setup for a Linux VM and a
Mac is in [docs/VALIDATION.md](docs/VALIDATION.md):

```bash
tools/validate.sh
```

It records the platform, builds, runs the unit and parity tests, times the fixture and writes
`validation/validation-HOST-ARCH.txt`; it exits non-zero if anything fails. To also compare a
full-size run, point it at your data and either a Sniffles 2.8.1 VCF of the same input or a
Sniffles 2.8.1 install to run on that machine:

```bash
INPUT=sample.bam ORACLE=sniffles2.vcf tools/validate.sh
INPUT=sample.cram REFERENCE=ref.fa SNIFFLES=$(command -v sniffles) tools/validate.sh
```

Genotype likelihoods use the platform's `pow` and `log`, as Sniffles2 does through Python, so a
last-bit difference between C libraries could in principle change a rounded QUAL or GQ; running
Sniffles2 on the same machine (`SNIFFLES=`) tells such a difference apart from a bug.

## Checking parity

`tools/parity.sh OUTDIR INPUT.bam [options...]` runs Sniffles2 and lrpipe-sv with the same
options and diffs the VCFs. `tools/combine_parity.sh OUTDIR ORACLE.vcf [UNSORTED.part.vcf] --
INPUT.tsv [options...]` compares multi-sample calling with a Sniffles2 combine output.

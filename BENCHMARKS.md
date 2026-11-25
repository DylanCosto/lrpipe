# Benchmarks

Every lrpipe tool rerun on whole genomes on 2026-11-03, against the tool it replaces on the same
inputs and machine; `lrpipe snv`'s GPU runs again on 2026-11-25 with 0.2.2, whose GPU calling
changed ([GPU calling in 0.2.2](#gpu-calling-in-022)).
Each tool's README has more runs, its parity tests and its history.

**Machine:** Intel Core i5-12600K (16 threads), 31 GB RAM, NVIDIA RTX 5060 Ti 16 GB, inputs on
NVMe. Every run used 16 threads. Wall time, CPU time (user + system) and peak resident memory
come from `/usr/bin/time -v`. Each step ran alone, one after another.

**Inputs:**

| Sample | Reads | Used by |
|---|---|---|
| HG003 ONT R10.4.1 sup v5.0.0, minimap2, GRCh38 | 110 GB BAM | snv, sv |
| HG002 ONT R10.4.1 sup, haplotagged, GRCh38 | 99 GB CRAM | sv, methyl pileup |
| HG004 ONT R10.4.1 sup v5.0.0, haplotagged, GRCh38 (ONT GIAB 2025.01, PAY87778) ⁴ | 58 GB CRAM | trio check |
| HG002 PacBio Revio HiFi, 6.53 M reads, bonefish `-x hifi -Y`, GRCh38 | 39 GB BAM | snv, sv, depth, cnv, tr, phase |
| HG002 PacBio Revio HiFi, PacBio's alignment, with 5mC calls | 45 GB BAM | methyl cpg |

## Speed

| Tool | Job | Original | lrpipe | Speed-up |
|---|---|---|---|---|
| `lrpipe snv` | HG003 ONT, CPU, exact | Clair3 2.0.3: 1 h 23 min, 5.8 GB ¹ | 1 h 10 min 16 s, 60,882 CPU-s, 11.5 GB | 1.2x |
| `lrpipe snv` | HG003 ONT, CPU, `--approximate` | | 28 min 25 s, 25,093 CPU-s, 11.3 GB | 2.9x |
| `lrpipe snv` | HG003 ONT, GPU (0.2.2) | | 7 min 36 s, 4,052 CPU-s, 11.0 GB | |
| `lrpipe snv` | HG002 HiFi, CPU, exact | | 36 min 41 s, 32,797 CPU-s, 7.6 GB | |
| `lrpipe snv` | HG002 HiFi, CPU, `--approximate` | | 26 min 20 s, 21,493 CPU-s, 9.0 GB | |
| `lrpipe snv` | HG002 HiFi, GPU (0.2.2) | Clair3 2.0.3: 33 min 45 s, 11.4 GB ¹ | 4 min 18 s, 2,307 CPU-s, 10.6 GB | 7.9x |
| `lrpipe sv` | HG003 ONT BAM | Sniffles 2.8.1: 126 s, 2.0 GB ¹ | 67.5 s, 554 CPU-s, 0.2 GB | 1.9x |
| `lrpipe sv` | HG002 ONT CRAM | Sniffles 2.8.1: 573 s ¹ | 71.7 s, 526 CPU-s, 0.3 GB | 8.0x |
| `lrpipe sv` | HG002 HiFi BAM | Sniffles 2.8.1: 72.6 s, 1.7 GB ¹ | 24.1 s, 224 CPU-s, 0.2 GB | 3.0x |
| `lrpipe depth` | HG002 HiFi, `-n -x -Q 20 -b 1000` | mosdepth 0.3.14: 63.3 s, 273 CPU-s, 1.9 GB | 24.2 s, 242 CPU-s, 9.6 GB | 2.6x |
| `lrpipe cnv` | HG002 HiFi coverage | Spectre 0.2.1: 148 s, 156 CPU-s, 1.0 GB | 3.2 s, 3.2 CPU-s, 0.3 GB | 46x |
| `lrpipe tr` | HG002 HiFi, 937,204 loci | TRGT 5.1.0: 645 s, 7,820 CPU-s, 5.2 GB | 77.7 s, 895 CPU-s, 3.3 GB | 8.3x |
| `lrpipe methyl pileup` | HG002 ONT CRAM, 5mC and 5hmC at CpGs | modkit 0.6.4: 799 s, 8,325 CPU-s, 13.5 GB | 251 s, 2,773 CPU-s, 2.1 GB ² | 3.2x |
| `lrpipe methyl cpg` | HG002 Revio, count mode, reference sites | pb-CpG-tools 3.0.0: 1,932 s, 26,746 CPU-s, 12.5 GB | 217 s, 2,317 CPU-s, 1.3 GB ² | 8.9x |
| `lrpipe phase` | HG002 HiFi calls, every contig | WhatsHap 2.8: 1,326 s, 6,974 CPU-s, 1.6 GB, 6 contigs at a time ³ | 125 s, 1,483 CPU-s, 8.8 GB | 10.6x |

¹ Not rerun: the original's earlier run on the same input and machine (Clair3 from
[lrpipe-snv's README](crates/lrpipe-snv/README.md), Sniffles2 from [lrpipe-sv's](crates/lrpipe-sv/README.md)).
Clair3 was run once per platform: on the CPU for ONT and on the GPU for HiFi.

² Rerun after a fix found by this benchmark. `lrpipe` fixed glibc's mmap threshold at 128 KiB for
every tool, to keep `lrpipe sv`'s read pass's memory down, which made `methyl cpg` map and unmap
a buffer per read: 389 s, 900 s of it in the kernel, and `methyl pileup` 303 s. The threshold is
now fixed for `lrpipe sv` only. The other `lrpipe` rows ran with it too; the standalone
`lrpipe-snv`, without it, took 4 min 12 s on HiFi on the GPU (11.3 GB) against `lrpipe snv`'s
4 min 42 s (10.2 GB) on 2026-11-03; earlier runs, before the threshold was fixed, took 3:57-4:16.

³ WhatsHap phases one contig per process, so it ran as 24 processes, 6 at a time. Its wall time
is the slowest contig's batch.

⁴ ONT's published `.crai` for this CRAM points every container 68 bytes past where it starts (it
indexes a version with a longer header; HG002's was 63 bytes off), which reads as a container
checksum mismatch. Every offset moved back 68 bytes lands on a container, and that index was used.

### `lrpipe snv` on an Apple GPU (Metal)

`--device metal` (a build with `--features metal`) on an Apple M4 Pro (12 CPU cores, 16 GPU
cores, 24 GB), HG002 SPRQ chr20, 12 threads, the machine otherwise idle. The 10x rows are the
SPRQ BAM downsampled (`samtools view -s 42.22`), with the `hifi_sequel2` model, as for
low-coverage Sequel II genomes; there `--approximate` keeps the pileup network in `float32`, as
its 8-bit check fails (2,546 candidates moved by up to 0.087).

| Job | CPU, exact | CPU, `--approximate` | Metal |
|---|---|---|---|
| chr20 10x, whole run | 134.6 s, 1,497 CPU-s | 84.0 s, 902 CPU-s | 52.3 s, 19 CPU-s |
| pileup calling | 50.5 s | 54.4 s | 20.8 s |
| full-alignment calling | 80.1 s | 25.8 s | 27.7 s |
| chr20 46x (`hifi_revio`), whole run | 53.2 s, 561 CPU-s | 30.2 s, 297 CPU-s | 23.4 s, 48 CPU-s |

Metal runs the networks in `float32`, as the CPU's exact mode does, so its calls are the CPU's:
on both chr20 runs every record, filter and genotype matched, QUALs within 0.01, and the GIAB
v4.2.1 F1 on chr20 is the same (10x: SNPs 0.9881, indels 0.8626; 46x: 0.9988, 0.9901). Its
kernels run at about 2 TFLOPS, near the M4 Pro's measured `float32` multiply-add rate (2.6),
so full-alignment calling is no faster than the CPU's 8-bit networks; the gain is in pileup
calling, which `--approximate` cannot speed up on low-coverage genomes, and in leaving the CPU
free (19 CPU-s against 902) for the stages `lrpipe run` runs alongside.

#### Regional HiFi and ONT calling against Clair3

On the M4 Pro, three callers processed the same 2 Mb of HG002 on GRCh38 (BED
`chr20 5000000 7000000`). Times are mean wall seconds over two runs with four threads
requested. Each run started from an indexed BAM and wrote an indexed small-variant VCF,
phased VCF and gVCF. Speedup is the Clair3 mean divided by the lrpipe mean.

| Reads | lrpipe Metal | Clair3 2.0.3 CPU | Speedup vs CPU | Local Clair3 MPS | Speedup vs MPS |
|---|---:|---:|---:|---:|---:|
| PacBio Revio HiFi, 51.0× | 3.16 s | 52.60 s | 16.7× | 54.53 s | 17.3× |
| ONT R10.4.1 SUP v5.0, 51.3× | 1.93 s | 46.26 s | 23.9× | 42.40 s | 21.9× |

The HiFi callers used the `hifi_revio` model and the ONT callers used
`r1041_e82_400bps_sup_v500`. The ONT BAM included reads overlapping chr20:4,000,001–8,000,000,
extracted from the HG002 PAW70337 CRAM before timing. Each caller used the same BAM within
each platform. On macOS 27.0.1, upstream Clair3 2.0.3 uses the CPU. The local comparison
build selected PyTorch 2.2.2's Metal Performance Shaders (MPS) backend on top of Clair3
`0f461d9`. The Apple GPU driver recorded command execution for both GPU callers.

The timed lrpipe build was `1c54b67` with `--no-default-features --features metal`.
Each platform ran in order lrpipe, upstream Clair3, local Clair3 MPS, then in reverse order.
Timing included process startup, phasing and compression. These speedups describe the 2 Mb
interval. [Individual timings and indexed record counts](tools/validation/results/apple-metal-hg002-2026-11-24.tsv)

On ONT, all three callers agreed on all 3,640 PASS calls by chromosome, position, REF, ALT
and genotype. On HiFi, lrpipe shared 3,769 PASS calls with each Clair3 build, with three
unique to lrpipe and four unique to each comparison build. In a separate exact CPU replay
of lrpipe on HiFi, all 3,775 merged records matched the Metal run in every field except
three QUAL values that differed by 0.01.

## Results against the original tools

Outputs compared after the run, without the header lines that name the program, its command
line or the date (`compare_all.sh` in the validation directory).

| Tool | Against | Result |
|---|---|---|
| `lrpipe sv` | Sniffles 2.8.1 | identical VCF records on all three genomes (27,046, 27,420 and 26,287), but for the breakend fix: with it, the 29, 23 and 33 BND records whose mate is on chrY or a contig named with an IUPAC letter (`chrUn_KI270746v1`) keep the mate's name, where Sniffles2 writes `chrN` or `chrUn_NI270746v1` |
| `lrpipe depth` | mosdepth 0.3.14 | every output file identical |
| `lrpipe cnv` | Spectre 0.2.1 | identical VCF and BED (29 CNVs) |
| `lrpipe methyl pileup` | modkit 0.6.4 | identical bedMethyl (4.7 GB) |
| `lrpipe methyl cpg` | pb-CpG-tools 3.0.0 | identical records in all three tables (28.9 million combined, 25.3 and 25.1 million per haplotype). The header lacks pb-CpG-tools' `##basemod-source`, which it copies from the BAM's primrose `@PG` record |
| `lrpipe phase` | WhatsHap 2.8 | identical records on all 24 contigs; also identical to lrpipe-snv's own phasing of the same calls (`phased_merge_output.vcf.gz`) |
| `lrpipe snv` | Clair3 2.0.3, CPU, HG003 ONT | PASS calls (position, alleles, genotype) agree on 99.984% exact and 99.780% with `--approximate` |
| `lrpipe snv` | Clair3 2.0.3, GPU, HG002 HiFi | 99.340% with 0.2.2 (99.707% before; [GPU calling in 0.2.2](#gpu-calling-in-022)) |
| `lrpipe snv` | its own exact CPU run | `--approximate` agrees on 99.795% (ONT) and 99.967% (HiFi) |

Rerun with the fixed binary, `lrpipe cnv`, `lrpipe tr` and `lrpipe methyl pileup` wrote the
same files as before the fixes; `lrpipe sv` changed exactly the BND records above.

### GPU calling in 0.2.2

lrpipe 0.2.2 changed GPU calling: the pileup network now takes the CPU path's inputs, with their
high-depth scaling, instead of Clair3's GPU-mode inputs, whose counts wrap past 127 in 8 bits.
On deep samples the wrapping made calls collapse (E. coli K-12 HiFi at 295x: SNP F1 42.5% and
indel F1 1.4% on the GPU, against 96.4% and 97.6% on the CPU; #45).

Both whole genomes were run again on the GPU with the published 0.2.2 bundle on 2026-11-25, on
the same machine and inputs. Against each run's own exact CPU calls, PASS calls (position, alleles,
genotype) agree on 99.575% on HiFi (99.470% before) and 99.683% on ONT (99.602% before); the rest
comes from TF32 rounding tipping close calls, and phasing, where one such call can re-phase a
chromosome and change full alignment's haplotags. Against Clair3's own GPU run the HiFi calls
agree less, 99.340% (99.707% before), as lrpipe no longer reproduces its wrapping. GIAB F1 is
unchanged within 0.0002 on both genomes (the tables below). Speed is unchanged: the GPU does the
same work on other numbers.

The 0.2.2 GPU calls were scored with touchstone
0.1.0's `vcfeval`, a Rust port of RTG Tools 3.13's, on the same split, regions and
stratifications as `score2.sh`. On the earlier HiFi GPU calls it wrote the same stratified table
and the same TP, FP and FN records as RTG, in 26 s and 4.4 GB at most against RTG's 16 GB heap.
The runs, scores and scripts are in
[tools/validation/results/gpu-0.2.2-2026-11-25](tools/validation/results/gpu-0.2.2-2026-11-25).

## Accuracy against truth sets

### SNVs and indels: GIAB v4.2.1

`rtg vcfeval` (RTG Tools 3.13), PASS calls, SNPs and other variants scored apart, within the
benchmark regions of chr1-22 (with the regions GIAB marks inconsistent removed). F1 overall and
in GIAB's stratifications:

HG003 ONT:

| | Clair3 CPU | lrpipe exact | lrpipe `--approximate` | lrpipe GPU |
|---|---|---|---|---|
| SNP, all | 0.99704 | 0.99704 | 0.99705 | 0.99703 |
| SNP, difficult regions | 0.98967 | 0.98967 | 0.98967 | 0.98963 |
| SNP, not difficult | 0.99865 | 0.99865 | 0.99866 | 0.99865 |
| SNP, homopolymers >= 7 bp | 0.98547 | 0.98548 | 0.98549 | 0.98543 |
| SNP, low mappability and segdups | 0.98448 | 0.98450 | 0.98449 | 0.98438 |
| SNP, segdups | 0.96941 | 0.96942 | 0.96955 | 0.96923 |
| SNP, MHC | 0.99195 | 0.99195 | 0.99187 | 0.99200 |
| INDEL, all | 0.89236 | 0.89231 | 0.89222 | 0.89245 |
| INDEL, difficult regions | 0.84468 | 0.84460 | 0.84447 | 0.84479 |
| INDEL, not difficult | 0.99630 | 0.99630 | 0.99630 | 0.99630 |
| INDEL, homopolymers >= 7 bp | 0.76592 | 0.76580 | 0.76566 | 0.76610 |
| INDEL, low mappability and segdups | 0.97353 | 0.97353 | 0.97341 | 0.97337 |
| INDEL, segdups | 0.96481 | 0.96481 | 0.96486 | 0.96486 |
| INDEL, MHC | 0.90870 | 0.90870 | 0.90839 | 0.90896 |

HG002 HiFi:

| | Clair3 GPU | lrpipe exact | lrpipe `--approximate` | lrpipe GPU |
|---|---|---|---|---|
| SNP, all | 0.99857 | 0.99857 | 0.99857 | 0.99858 |
| SNP, difficult regions | 0.99311 | 0.99313 | 0.99314 | 0.99315 |
| SNP, not difficult | 0.99979 | 0.99979 | 0.99978 | 0.99979 |
| SNP, homopolymers >= 7 bp | 0.99474 | 0.99474 | 0.99475 | 0.99472 |
| SNP, low mappability and segdups | 0.98724 | 0.98733 | 0.98734 | 0.98735 |
| SNP, segdups | 0.97343 | 0.97357 | 0.97358 | 0.97364 |
| SNP, MHC | 0.99605 | 0.99605 | 0.99605 | 0.99623 |
| INDEL, all | 0.96875 | 0.96879 | 0.96880 | 0.96875 |
| INDEL, difficult regions | 0.95726 | 0.95733 | 0.95735 | 0.95728 |
| INDEL, not difficult | 0.99592 | 0.99588 | 0.99588 | 0.99586 |
| INDEL, homopolymers >= 7 bp | 0.93782 | 0.93791 | 0.93794 | 0.93783 |
| INDEL, low mappability and segdups | 0.97151 | 0.97171 | 0.97171 | 0.97140 |
| INDEL, segdups | 0.96178 | 0.96213 | 0.96213 | 0.96164 |
| INDEL, MHC | 0.91945 | 0.91945 | 0.91945 | 0.91945 |

Every lrpipe mode is within 0.00015 of Clair3's F1 overall, and within 0.0004 in every
stratification; the largest gaps are indels in segmental duplications. The GPU columns are from
0.2.2 (see [GPU calling in 0.2.2](#gpu-calling-in-022)); the other lrpipe columns are from
2026-11-03. The ONT indel F1 is lower here than
in [lrpipe-snv's VALIDATION.md](crates/lrpipe-snv/VALIDATION.md) (0.9042, hap.py), which scored with
hap.py's matching; this table compares the callers to each other under one method.

### Tandem repeats: GIAB HG002 TR v1.0.1 Tier 1

`lrpipe-tr-bench`: both alleles' length change exactly right, by the truth's largest change. The
catalog covers only part of Tier 1, so each caller calls about a quarter of its 1,638,106
regions:

| | Called | Exactly right, of called | of all Tier 1 | ±1 bp, of called |
|---|---|---|---|---|
| lrpipe-tr | 393,498 (24.02%) | 90.62% | 21.77% | 96.14% |
| TRGT 5.1.0 | 392,120 (23.94%) | 87.57% | 20.96% | 95.59% |

On the 391,984 regions both call (`--common`): lrpipe-tr 90.78%, TRGT 87.59% exactly right.
[lrpipe-tr's README](crates/lrpipe-tr/README.md#results-hg002-hifi-31x-grch38-16-threads) breaks this
down by the size of the change.

### Structural variants: GIAB HG002 v5.0q

`lrpipe run`'s phased SVs (`sv/SAMPLE.sv.vcf`: the read pass, then `sv phase` from the run's
haplotag list), HG002 whole genome, Truvari 5.5 in GIAB v5.0q's benchmark regions (`--passonly
--pctseq 0`; multi adds `--pick multi`; GT is genotype concordance of strict matches), 2026-11-06:

| | strict P | strict R | strict F1 | multi F1 | GT |
|---|---:|---:|---:|---:|---:|
| HiFi, `--sv-mode sniffles2` (= Sniffles 2.8.1) | 0.923 | 0.716 | 0.807 | 0.871 | 0.839 |
| HiFi, `--sv-mode accuracy` (the default) | 0.907 | 0.773 | **0.835** | **0.886** | **0.900** |
| ONT CRAM, `--sv-mode sniffles2` | 0.917 | 0.725 | 0.810 | 0.875 | 0.834 |
| ONT CRAM, `--sv-mode accuracy` | 0.898 | 0.782 | **0.836** | **0.890** | **0.896** |

Accuracy mode's GT is with `--read-support-genotypes` (2026-11-10; 0.1.8: 0.895 HiFi, 0.888
ONT), which changes genotypes only.

Accuracy mode's lower strict precision is mostly how strict matching counts: most of the
calls it adds lie within 500 bp of another call of the same type, the second allele of a repeat
expansion, and a second call matched to a truth SV another call already matched is a false
positive under one-to-one matching. Under `--pick multi` its precision is the same or higher
(HiFi 0.934 against 0.930, ONT 0.924 against 0.923; 56 and 154 more false positives). The SV
calls take the same time either way (HiFi 29 and 26 s, ONT 66 and 66 s). Sniffles2's own
tandem repeat annotation (`[sv] tandem_repeats`, its `human_GRCh38_no_alt_analysis_set.trf.bed`)
lowers strict F1 to 0.777 (HiFi) and 0.778 (ONT), 0.796 and 0.794 with accuracy mode: it merges
calls in a repeat that GIAB records apart.

### Tandem repeat expansions: Coriell cell lines

PacBio's public PureTarget data (45 HiFi BAMs of 33 Coriell lines with expansions at DMPK, FMR1,
FXN, HTT, ATXN1, ATXN3, ATN1, AR, PABPN1, C9orf72 and RFC1), each allele against the repeat
count in Coriell's catalog (within max(3 units, 10%)), on TRGT's regions: lrpipe-tr matches 24
of 38 samples, PacBio's TRGT calls 22, and lrpipe-tr agrees with TRGT on 41 of 45 (within
max(10 bp, 5%)). [lrpipe-tr's README](crates/lrpipe-tr/README.md#pathogenic-expansions-coriell-cell-lines-pacbio-puretarget)
lists what both miss and why.

### Phasing

Against GIAB's HG002 Q100 benchmark (v5.0q, GRCh38 small variants), whose genotypes are phased
chromosome-wide from the diploid assembly. (v4.2.1's benchmark VCF is unphased; its separate
phased file, phased by trio and Strand-seq for about 70% of heterozygous variants, scored other
call sets in [lrpipe-phase's benchmarks](crates/lrpipe-phase/benchmarks/hg002/README.md#phasing-accuracy).) Scored are `lrpipe run`'s
phased calls (`phased_merge_output.vcf.gz`, from its 2026-11-02 whole-genome runs) on the HiFi
BAM and the ONT CRAM: PASS calls on
chr1-22 inside the benchmark regions, SNVs and indels. Switch errors are from `whatshap compare`,
over neighbouring variants heterozygous in both and in one block; block lengths are from
`whatshap stats`:

| | HG002 HiFi | HG002 ONT |
|---|---|---|
| Heterozygous calls phased | 99.9% | 99.9% |
| Variant pairs assessed | 2,554,138 | 2,477,064 |
| Switch errors | 7,909 (0.31%) | 9,867 (0.40%) |
| of them single switches, and flips (one variant wrong, two switches) | 2,331 and 2,789 | 727 and 4,570 |
| Hamming errors, per block | 6.1% | 6.8% |
| Blocks | 7,567 | 2,175 |
| Block NG50 | 606 kb | 1.99 Mb |

`lrpipe phase` on the same HiFi calls scores the same, as it writes the same records. As those are
WhatsHap 2.8's, so are these scores (`score_phasing.sh` in the validation directory).

### Mendelian consistency: the HG002, HG003, HG004 trio (ONT)

Every truth set above but HG003's small variants is HG002's, the genome the defaults were tuned
on. Inheritance checks the other two without truth: HG002 has one allele from his father HG003
and one from his mother HG004 at each autosomal locus, his mother's on chrX and his father's on
chrY, so a call that breaks this is wrong in at least one of the three. `lrpipe run` 0.1.8 on
ONT's GIAB 2025.01 data for all three (R10.4.1 sup v5.0.0, wf-human-variation's CRAMs; HG003 as
the BAM above), 2026-11-10: 1,658 s, 966 s and 883 s. `trio.py` and `duo.py` in the validation
directory.

| | Calls checked | Mendelian errors |
|---|---:|---:|
| SNPs ¹ | 4,447,650 | 0.20% |
| Indels ¹ | 718,437 | 9.38% |
| SVs, `--sv-mode sniffles2` ² | 36,018 | 6.61% |
| SVs, `--sv-mode accuracy`, 0.1.8 ² | 37,083 | 5.05% |
| SVs, `--sv-mode accuracy` with `--read-support-genotypes` ² | 37,083 | **4.54%** |

¹ PASS, chr1-22, inside both HG002's and HG003's GIAB v4.2.1 regions; a site missing from a
sample's VCF counts as reference there. On the father-son pair alone, where GIAB's truth exists
for both, lrpipe has 1,854 SNPs (0.047%) where one is homozygous for an allele the other lacks,
and GIAB's truth 98 (0.002%). HG003 alone scores as in the table above (SNP F1 0.99704, indel
0.89213). The three VCFs are normalised first (`bcftools norm -m -any`, left-aligned against
GRCh38), so an indel written differently in two of them is one indel; compared as written,
without normalising, the indels' rate is 11.6% (85,415 of 738,681) and the SNPs' 0.20%. The
indel errors left include the same change written as different sets of alleles, which
normalising does not join, besides the genotype errors.

² PASS, chr1-22, as `lrpipe run` writes them (the read pass, then `sv phase` from the run's
haplotag list); an SV matches another sample's of the same type within 500 bp and a size ratio
of 0.7 (looser matching, 2 kb and 0.3, gives 4.29% for 0.1.8). Accuracy mode, the default, is
tuned on HG002 alone, so its lower rate on the trio than Sniffles2's is its first check on other
genomes; the gain is in calls under 300 bp (deletions 5.68% to 3.94%, insertions 6.97% to
4.15%). In 0.1.8 insertions of 1 kb or more were the weak spot, 14.7% (15.1% in Sniffles2
mode), with confident calls (1/1 calls in an error: median 19 reads, VAF 0.95). HG002's GIAB
truth showed most were heterozygous insertions called 1/1, in a parent or in the son, nearly all
of 2.5 kb or more: Sniffles2 scales those insertions' support up by their length, for reads
clipped inside them, though it already counts those reads. `--read-support-genotypes`, now part
of `--accuracy` ([lrpipe-sv's README](crates/lrpipe-sv/README.md#--read-support-genotypes)), genotypes
them from their reads, with the same calls passing: 6.26% of them break inheritance. Of all
1,683 errors left, 52% are in a repeat of lrpipe-tr's catalog whose lrpipe-tr genotypes are
consistent across the trio (within a motif unit or 5%): the SV calls write the same change in
different pieces or sizes in each sample (a 568 bp insertion as 198 + 376 bp in the mother, as
198, 580 and 376 bp in the son), or leave one sample's below a filter. 41% are outside the
catalog, and 6% in repeats lrpipe-tr's genotypes break inheritance at too.

Tandem repeats, by allele length: a locus is consistent when the son's alleles can come one
from each parent (on chrX from his mother, on chrY from his father):

| Loci | Number | Exact | ±1 bp | ±1 motif or 5% |
|---|---:|---:|---:|---:|
| chr1-22 | 905,698 | 65.85% | 93.71% | 99.85% |
| chr1-22, in GIAB's HG002 Tier 1 regions | 265,016 | 73.60% | 95.80% | 99.96% |
| … homopolymers | 59,983 | 42.86% | | |
| … motifs of 2-6 bp | 107,139 | 88.82% | | |
| … longer motifs | 97,894 | 75.78% | | |
| chrX | 25,256 | 83.89% | 94.33% | 99.42% |
| chrY (son and father) | 3,528 | 82.48% | 92.29% | 99.26% |

If each sample is exactly right at about 86% of loci, as HG002 ONT is against Tier 1, all three
are right at about 64%, and an error does not always break consistency; 73.6% fits that, so the
other two genomes look no worse than HG002. Nearly
every inconsistency is within a motif unit, and most are homopolymer lengths: on chrY, which a
son inherits whole, his alleles equal his father's at 92% of loci with 2-6 bp motifs and 32% of
homopolymers, almost all of those 1 bp apart. Father and son alone, Tier 1 loci: 84.3% share an
exact allele with 0.1.8 against 77.6% with the earlier runs (0.1.0 for HG002, 0.1.2 for HG003),
from the ONT length model.

## Haplotagging

`lrpipe run` haplotags the reads from its phased SNVs (the haplotag, repeats and sv-phase
stages, and the CpG tables after them). Measured 2026-11-03 on this machine while other builds
were running (load average 4 to 9, and 18 to 38 during the whole-genome steps), so the times
are upper bounds.

**The haplotag stage alone, whole genomes** (`lrpipe haplotag --list-only`, 16 threads, on the
2026-11-02 runs' `phased_merge_output.vcf.gz`):

| Input | Wall | Peak memory | Alignments listed | Assigned a haplotype | The list in memory (`--haplotags`) |
|---|---|---|---|---|---|
| HG002 HiFi, 39 GB BAM | 96 s (1,052% CPU), again 108 s | 2.8 GiB | 6.53 M | 5.48 M | 378 MiB (400 at the peak of reading it) |
| HG002 ONT, 99 GB CRAM | 253 s (1,194% CPU), again 311 s | 5.6 GiB, again 5.8 | 8.16 M | 6.20 M | 467 MiB |

**The repeats stage alone, whole genomes** (`lrpipe tr --haplotags` with those lists, 937,204
loci, 16 threads): HG002 HiFi BAM 88 s and 3.9 GiB, HG002 ONT CRAM 231 s and 6.7 GiB; every
primary alignment was in the list. In the read pass, without phasing, the repeats cost less
(HG002 ONT chr20: the pass took 6.0 s with them and 1.9 s without, `lrpipe tr` alone 5.8 s).
Their haplotypes matter: on HG002 ONT chr20 (GIAB TR Tier 1, `lrpipe-tr-bench`) 84.4% of the
called regions' alleles are exactly right with the CRAM's HP tags, 69.2% with them removed.

**The read pass's phasing evidence** (`--phasing-evidence`, HG002 chr1, 16 threads): peak memory
189 MB against 171 on HiFi, 332 against 315 on the ONT CRAM, and no measurable time.

**`lrpipe run` on chr20** (the reads of chr20 taken from the inputs above, the whole GRCh38
reference and repeat catalogue, 16 threads, GPU), before (`bd9c011`), with haplotagging but the
repeats still in the read pass, and with the repeats stage, two runs each:

| Input | Analyses | Before | Haplotagging | With the repeats stage | haplotag, repeats, sv-phase stages |
|---|---|---|---|---|---|
| HG002 HiFi BAM (no 5mC calls) | all | 39.9 s, 34.3 s | 43.0 s, 42.9 s | 48.3 s, 47.3 s | 3 s, 9 s, under 1 s |
| HG002 ONT CRAM | all but cnv ¹ | 49.2 s, 47.8 s | 56.0 s, 57.1 s | 68.3 s, 69.8 s | 7-8 s, 10 s, under 1 s |

On a whole genome the run waits for the haplotag stage after lrpipe-snv (96-108 s on HiFi,
253-311 s on the ONT CRAM, unless the read pass is still running), then the repeats (88 s and
231 s) and sv-phase (0.1 s on chr20) before the merge; the HiFi CpG tables (217 s,
[above](#speed)) start after the haplotag stage instead of with the run, beside the repeats.
So on ONT the repeats stage adds about its own time to the run, and on HiFi the CpG tables
about theirs.

**`lrpipe run` on whole genomes, 0.1.5 against this release** (`b00fa56`), on 2026-11-04: the
inputs above, the same site configuration (GPU), every analysis, one run after another on an
otherwise idle machine:

| Input | 0.1.5 | This release | haplotag | repeats | Peak memory (lrpipe-snv) |
|---|---|---|---|---|---|
| HG002 HiFi BAM (bonefish) | 6 min 57 s | 10 min 2 s | 96 s | 105 s | 15.7 GB, 13.4 GB |
| HG002 ONT CRAM | 19 min 38 s | 24 min 19 s | 269 s | 227 s | 13.9 GB, 14.0 GB |

The time added is the haplotag stage and then the repeats (and on HiFi the CpG tables, beside
them) after lrpipe-snv: work that 0.1.5 did beside lrpipe-snv from the input's own tags, or
without haplotypes. The read pass, without the repeats, needs 1.9 GB rather than 3.6 (HiFi) and
6.5 (ONT), so on ONT it now runs beside lrpipe-snv where 0.1.5's waited 12 min for memory.
This HiFi BAM has no `MM`/`ML` tags: 0.1.5 and this release alike wrote empty CpG tables from
it (as in the timings above); `lrpipe run` now leaves methylation out for such reads, with a
warning.

What changed, by `lrpipe compare` and the scorers above:
- SNVs and indels, phased and not, depth, CNVs (but for one END now clipped to its contig's
  end) and the ONT bedMethyl: identical.
- SVs: every 0.1.5 record is unchanged apart from its phase; the HiFi BAM, without `HP` tags,
  gains 196 PASS calls (Sniffles2's rescue of phased calls), and on the ONT CRAM, tagged by
  another pipeline, 17 rescued calls go and 14 come with the run's own haplotypes.
- Tandem repeats, GIAB HG002 TR v1.0.1 Tier 1 (`lrpipe-tr-bench`), exactly right of called:
  HiFi 90.62% → 94.24% (±1 bp 96.14% → 97.49%), since 0.1.5 split the untagged reads' alleles
  by length alone; ONT 85.04% → 84.90% (the CRAM's own tags against the run's).

¹ lrpipe-cnv stopped on these chr20 ONT runs, before and with haplotagging alike: lrpipe-snv's
`merge_output.vcf.gz` lacked the newline after its `#CHROM` line (fixed since).

**Results.** Phased from the list, the SVs, their SNF and the CpG tables are those of the
reads haplotagged into a BAM by `lrpipe haplotag` with the same phased SNVs, on whole chr20
(16 threads): HG002 HiFi (bonefish BAM) 876 SVs, 743 phased, 4 passed by the phasing rescue;
HG002 Revio BAM 854 SVs and 708,358 and 703,841 rows in the `hap1` and `hap2` tables; HG002 ONT
CRAM 933 SVs, 8 rescued, and `hap1`/`hap2` 734,394 and 731,483 rows: the same VCF records, SNF
and tables. The list without the BAM is the list written beside it. `lrpipe run` on the Revio
chr20 reads with their HP and PS tags removed wrote 709,229 and 705,180 `hap1`/`hap2` rows and
677 phased SVs, 674 with a phase set of the phased SNVs and 3 phased without one (their reads'
phase sets disagree, which Sniffles2 leaves out).

## Sex chromosomes (LRP-003)

HG002 HiFi chrX and chrY, called as `lrpipe run` now calls an XY sample (`--gender male` and the
GRCh38 PARs as `--par-regions-bed`, GPU, `--approximate`, `--enable-phasing`), then phased with
`lrpipe phase`. All 13,113 diploid calls lie inside chrX's PARs; the 111,068 chrX and 13,027 chrY
calls outside them are haploid. Both phasing steps completed, leaving the haploid calls
unphased; before the fix, phasing stopped at the first haploid call.

## Reproducing

The scripts are in [tools/validation](tools/validation), with the inputs, versions, model
hashes and truth sets they used, and their raw results (TP, FN and FP per stratification, and
the phasing comparisons per chromosome): `run_all.sh` runs every step above, timed, and skips the
steps already done; `compare_all.sh` compares the outputs and scores them; `rerun_fixed.sh`
reran the tools whose output the fixes could change. GIAB scoring uses `score2.sh` (rtg vcfeval)
and `giab_strat.py`, and phasing `score_phasing.sh`. [lrpipe-snv's BENCHMARKING.md](crates/lrpipe-snv/BENCHMARKING.md)
describes its checks in more detail.

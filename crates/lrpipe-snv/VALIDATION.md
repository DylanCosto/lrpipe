# Validation

lrpipe-snv was validated on GIAB samples on three platforms in two ways:

- **Parity**: Clair3 2.0.3 and lrpipe-snv run with the same options, and every output file both
  write is compared (`tools/option_parity.sh`). This was done on HG002 for ONT, PacBio HiFi and
  Illumina on chr20, and on the whole genome for HiFi and ONT on the CPU and ONT on the GPU.
- **Accuracy**: hap.py against GIAB v4.2.1. This was done on HG003 for ONT and Illumina (HG003
  is in none of Clair3's documented training sets), whole genome and chr20 (left out of every
  Clair3 model's training), and on HG002 chr20 for all three platforms.

Run on 2026-10-27/28, at commit `282d716`, when lrpipe-snv ran the networks through libtorch, the
same build as Clair3's. Its exact mode then wrote the same VCFs as Clair3, so its accuracy was
Clair3's on the same reads.

libtorch has since been removed: the networks run in lrpipe-snv's own Rust code, whose
probabilities differ from libtorch's in the last digits, so the files below are no longer
byte-identical to Clair3's. Against the former libtorch backend on HG002 ONT (16 threads,
i5-12600K, RTX 5060 Ti), the PASS calls agree on 100% of chr20 on the CPU, with the same GIAB
recall, 99.945% of chr20 on the GPU, with GIAB F1 the same to five decimals, and 99.827% of the
whole genome on the GPU (see the README's [Networks](README.md#networks)). The runs below have
not been repeated without libtorch on the VM. On the development machine, without libtorch, the
PASS calls agree with Clair3's own runs on 99.984% of the whole HG003 ONT genome on the CPU, on
100% of HG002 HiFi chr20 on the CPU and 99.469% on the GPU, and on 100% of every option case
on ONT, HiFi and Illumina regions (README, [Against Clair3](README.md#against-clair3)).

## Setup

| | |
|---|---|
| machine | Intel Xeon Gold 5418Y (24 cores), 121 GB RAM, NVIDIA L40S-24Q (driver 580.178.04) |
| Clair3 | bioconda `clair3=2.0.3`, `pytorch=2.13.0` (CPU, MKL), `longphase=2.1.1`; models from the package |
| Clair3 on the GPU | the same scripts with `--python` set to a conda-forge `pytorch=2.13.0=cuda129*` env (Python 3.11, for `libclair3`'s cffi module) |
| lrpipe-snv | `282d716`: `cargo build --release` against the Clair3 env's libtorch; for the GPU, `target-cuda/` against the CUDA env's |
| Illumina oracle | Clair3 2.0.3 with its realigner built (see [Illumina](#illumina-oracle)) |
| hap.py | hap.py 0.3.15, rtg-tools 3.13, `--engine vcfeval`, GIAB stratifications v3.3 |

### Data

| sample | platform | reads | model |
|---|---|---|---|
| HG002 | ONT R10.4.1, dorado sup v5.0.0 | ONT Open Data `giab_2025.01`, wf-human-variation `sup/HG002/PAW70337` haplotagged CRAM | `r1041_e82_400bps_sup_v500` |
| HG002 | PacBio Revio (SPRQ) | haplotagged GRCh38 CRAM, converted to BAM (below) | `hifi_revio` |
| HG002 | Illumina NovaSeq PCR-free 35x | DeepVariant case-study test data, `grch38_no_alt` BAM | `ilmn` |
| HG003 | ONT R10.4.1, dorado sup v5.0.0 | ONT Open Data `giab_2025.01`, `sup/HG003/PAY87794` haplotagged CRAM | `r1041_e82_400bps_sup_v500` |
| HG003 | Illumina NovaSeq PCR-free 35x | DeepVariant case-study `HG003.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam` | `ilmn` |

- Reference: GIAB's `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`, except for the Revio
  CRAM, which only decodes with UCSC's `hg38.analysisSet.fa`. The two differ only where UCSC
  turned GRCh38's IUPAC codes into `N`.
- Truth: GIAB v4.2.1 `HG00{2,3}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` with the
  `_noinconsistent.bed` regions.

Notes on the inputs:

- The `.crai` files next to the ONT CRAMs on S3 do not match the CRAMs: every region query fails
  with `Container header CRC32 failure`. Re-index them with `samtools index`.
- The Revio CRAM is CRAM 3.1. Clair3 2.0.3's `libclair3.so` bundles htslib 1.15.1, which fails
  on it with `cram_next_slice: Failure to decode slice`. For Clair3 the reads were converted to
  BAM with samtools 1.24 (63 GB, against the CRAM's 20 GB). lrpipe-snv reads the CRAM directly
  and writes the same VCFs as from the BAM.

## Parity

Each run below is `tools/option_parity.sh` with default options, with libtorch. A file counts as
the same if it is byte-identical apart from the `##cmdline` header lines. The files compared are
every top-level VCF and BAM, and in `tmp/` every chunk and job VCF, phased VCF, candidate BED,
split BED and chunk list.

### CPU

| HG002 | contigs | threads | files compared | differ | `merge_output` records | Clair3 | lrpipe-snv |
|---|---|---|---|---|---|---|---|
| ONT sup v5.0.0 | chr20 | 16 | 48 | 0 | 148,111 | 717 s | 69 s |
| Revio HiFi | chr20 | 16 | 40 | 0 | 137,188 | 392 s | 36 s |
| Illumina 35x | chr20 | 16 | 43 | 0 | 123,718 | 1,317 s | 502 s |
| Revio HiFi | chr1-22, X, Y | 20 | 1,236 | 0 | 5,286,484 | 2,247 s | 800 s |
| ONT sup v5.0.0 | chr1-22, X, Y | 20 | 1,633 | 0 | 5,794,669 | 4,068 s | 1,693 s |

`option_parity.sh` runs Clair3 and then lrpipe-snv, each under `nice -n 19`.

### GPU

The GPU runs went through `option_parity.sh` twice:

- `ONLY_CLAIR3=1` with `--use_gpu --device=cuda:0 --python=<CUDA env python3>`;
- then `SKIP_CLAIR3=1 LRPIPE_SNV=target-cuda/release/lrpipe-snv` with `--use_gpu --device=cuda:0`.

| HG002 ONT | files compared | differ | `merge_output` records | Clair3 | lrpipe-snv |
|---|---|---|---|---|---|
| chr20 | 30 | 6 (same records) | 148,886 | 590 s | 40 s |
| chr1-22, X, Y | 573 | 6 (same records) | 5,803,783 | 3,328 s | 672 s |

On both runs the six differing files were the stage intermediates, not the results:

- Clair3's GPU mode splits them across its GPU workers (`pileup_0..3.vcf`,
  `full_alignment_0..1.vcf`), while lrpipe-snv writes one per stage (`pileup_0.vcf`,
  `full_alignment_0.vcf`).
- Concatenated and sorted, both sides held the same records: 10,463,031 pileup and 4,515,474
  full-alignment records on the whole genome.
- `pileup.vcf.gz`, `full_alignment.vcf.gz`, `merge_output.vcf.gz` and every phasing and
  candidate file were identical.

GPU calls differ from CPU calls in Clair3 as in lrpipe-snv: 5,803,783 against 5,794,669 on the
whole genome.

GPU use on the whole genome was sampled from `nvidia-smi` every 5 s; Clair3's first 9 minutes
were not sampled.

| | sampled | samples with GPU > 0% | GPU ≥ 50% | mean utilization | peak memory |
|---|---|---|---|---|---|
| Clair3 | 46 min | 34% | 9% | 11% | 7.0 GB |
| lrpipe-snv | 12 min | 67% | 43% | 37% | 1.8 GB |

The two schedule the GPU differently:

- **Clair3** builds every chunk's tensors on the CPU into `.npy` files, then scores them in a
  burst with 4 GPU workers. Pileup took 644 s of tensors, then 254 s on the GPU.
- **lrpipe-snv** scores tensors as the CPU threads build them. Pileup calling took 293 s, and
  phasing and full alignment, run together, finished at 633 s.

### Illumina oracle

Stock bioconda Clair3 2.0.3 has no compiled realigner, so for Illumina its full-alignment stage
writes nothing, and `merge_output` is a copy of the pileup VCF. Its Dockerfile builds the
realigner, but `free_memory` is `static` in `realigner.cpp` (since commit 583cf9a), so the
Python side cannot find it.

The oracle is the bioconda scripts with `preprocess/realign/realigner` and `debruijn_graph` built
by the Dockerfile's commands, after removing `static` from `free_memory` (`realigner.cpp:862`):

```
g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 -I$BOOST/include debruijn_graph.cpp
```

It is passed to `option_parity.sh` as `CLAIR3_RUN`. On HG002 chr20 both write 282,765 pileup,
237,680 full-alignment and 123,718 merged records.

## Accuracy

The command was:

```
hap.py TRUTH.vcf.gz merge_output.vcf.gz -f TRUTH_noinconsistent.bed \
  -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --engine vcfeval --engine-vcfeval-template GRCh38_no_alt.sdf \
  --stratification GRCh38@all/GRCh38-all-stratifications.tsv -o happy
```

with `-l chr20` for chr20. The tables count PASS calls only. `tools/happy_summary.py
happy.extended.csv` prints the rows below.

### HG003, whole genome (chr1-22)

lrpipe-snv exact CPU mode, with libtorch. ONT ran on 22 threads in 1,403 s; Illumina on 20
threads in 10,331 s, 8,953 s of it in full alignment with realignment.

| F1 | ONT SNP | ONT INDEL | Illumina SNP | Illumina INDEL |
|---|---|---|---|---|
| all | **0.9973** | **0.9042** | **0.9949** | **0.9942** |
| all difficult regions | 0.9914 | 0.8620 | 0.9743 | 0.9923 |
| not in difficult regions | 0.9986 | 0.9977 | 0.9994 | 0.9986 |
| homopolymers ≥7 bp | 0.9886 | 0.7870 | 0.9983 | 0.9963 |
| low mappability and segdups | 0.9865 | 0.9813 | 0.9414 | 0.9359 |
| segdups | 0.9737 | 0.9746 | 0.9243 | 0.9337 |
| MHC | 0.9967 | 0.9678 | 0.9742 | 0.9553 |

| HG003, all | TP | FN | FP | recall | precision |
|---|---|---|---|---|---|
| ONT SNP | 3,316,200 | 11,295 | 6,390 | 0.9966 | 0.9981 |
| ONT INDEL | 446,677 | 57,824 | 37,989 | 0.8854 | 0.9239 |
| Illumina SNP | 3,306,273 | 21,223 | 12,494 | 0.9936 | 0.9962 |
| Illumina INDEL | 500,561 | 3,940 | 1,927 | 0.9922 | 0.9963 |

Notes:

- **ONT** matches ONT's own hap.py result for the same HG003 run (wf-human-variation's Clair3,
  `giab_2025.01/analysis/happy-benchmark/sup/`): SNP 0.9973, INDEL 0.9047. Its INDEL F1 comes
  from the reads: homopolymers ≥7 bp hold 96% of its INDEL false negatives (55,302 of 57,824)
  and 95% of its false positives.
- **Illumina** loses its SNPs where short reads do not map uniquely: low-mappability regions
  and segmental duplications hold 19,896 of its 21,223 SNP false negatives.
- No other caller (DeepVariant, GATK) was run on these reads.

### HG002 HiFi, whole genome (chr1-22), GPU, against Clair3

lrpipe-snv (lrpipe 0.1.3's code, `merge_output.vcf.gz`, unphased) and Clair3 2.0.3 on the same
HG002 Revio HiFi BAM on the GPU (RTX 5060 Ti; their PASS calls agree on 99.707%), scored with
rtg vcfeval 3.13 (hap.py's engine) against GIAB v4.2.1 within
`HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`, SNPs and other variants apart, PASS
calls only (no stratifications):

| | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|
| SNP, lrpipe-snv | 3,362,233 | 6,643 | 2,894 | 0.99803 | 0.99914 | **0.99858** |
| SNP, Clair3 | 3,362,238 | 6,709 | 2,889 | 0.99801 | 0.99914 | 0.99857 |
| INDEL, lrpipe-snv | 515,812 | 23,632 | 9,585 | 0.95620 | 0.98176 | **0.96881** |
| INDEL, Clair3 | 515,786 | 23,665 | 9,611 | 0.95614 | 0.98171 | 0.96875 |

The 0.3% of PASS calls on which the two differ leave accuracy the same: lrpipe-snv has 66 fewer
false-positive SNPs and 5 more false-negative ones, and 26 more true and 33 fewer false-positive
indels. Most of the differences come from Clair3's GPU path (pileup counts wrapping in 8 bits
at deep sites, which lrpipe-snv reproduced in this run and no longer does) and from phasing, where one het SNP across the
cut-off re-phases a chromosome ([README](README.md#against-clair3)).

### chr20

| F1 | SNP | INDEL | INDEL, homopolymers ≥7 bp | SNP, all difficult |
|---|---|---|---|---|
| HG003 ONT sup v5.0.0 | 0.9984 | 0.8956 | 0.7669 | 0.9939 |
| HG003 Illumina 35x | 0.9961 | 0.9954 | 0.9968 | 0.9755 |
| HG003 Illumina 35x, pileup only (stock Clair3 2.0.3) | 0.9943 | 0.9828 | | |
| HG002 Revio HiFi | 0.9994 | 0.9963 | 0.9942 | 0.9969 |
| HG002 ONT sup v5.0.0 | 0.9991 | 0.9027 | | |
| HG002 Illumina 35x | 0.9964 | 0.9947 | | |

HG002 is a training sample, on chromosomes other than chr20, for `hifi_revio` and `ilmn`. The
pileup-only row is what stock bioconda Clair3 2.0.3 reports for Illumina without a realigner;
full alignment with realignment adds 0.0126 INDEL F1.

# Whole-genome validation

The [0.2.2 GPU report](results/gpu-0.2.2-2026-11-25/README.md) reruns `lrpipe snv`'s GPU calling
on both whole genomes after 0.2.2 changed its inputs, scored with touchstone's `vcfeval`.

The [0.2.1 option and resume report](results/release-0.2.1-options-2026-11-25/README.md)
records pipeline combinations on Linux, repeat VCF indexing, checkpoint recovery, and CPU and CUDA comparisons.

The [2026-11-25 CNV and copy-neutral LOH report](results/human-cnv-loh-2026-11-25/README.md) covers independently
confirmed events in HG008-T and HCC1395-T, the matched-normal and GIAB controls,
and runtimes from indexed whole-genome alignments.

The [HG006 and HG007 ONT depth study](results/ont-indel-depth-2026-11-24/README.md)
compares lrpipe and Clair3 on chr20 at native depth and across three seeds per downsampled
depth. It includes indel accuracy, CUDA runtimes, and scores inside and outside homopolymers.

The scripts behind [BENCHMARKS.md](../../BENCHMARKS.md), and their raw results. Every path comes
from [`paths.sh`](paths.sh), which defaults to where each file was on the machine that ran
BENCHMARKS.md; set the variables it names to run them elsewhere.

| Script | What it does |
|---|---|
| `run_all.sh` | Every lrpipe tool and the original tools without a kept output, on whole genomes, each step timed (`/usr/bin/time -v`) into `logs/`. A step whose `done/` marker exists is skipped, so it can be run again after a failure. |
| `rerun_fixed.sh` | The tools whose output the review fixes could change, rerun into `out-fixed/`. |
| `compare_all.sh` | Each lrpipe output against the original tool's (ignoring only header lines naming the program, its command or the date), then the accuracy scores: `score2.sh`, `lrpipe-tr-bench` and `score_phasing.sh`. |
| `score2.sh` | SNVs and indels against GIAB v4.2.1 with rtg vcfeval, SNPs and other variants apart, PASS calls in the benchmark regions of chr1-22 without the ones GIAB marks inconsistent; `giab_strat.py` then counts TP, FN and FP in each GIAB stratification. |
| `score_phasing.sh` | Phasing against GIAB's HG002 Q100 v5.0q: `whatshap compare` and `whatshap stats` on PASS calls of chr1-22 in the benchmark regions. |
| `male_chrx.sh` | HG002 HiFi chrX and chrY called and phased as an XY sample (LRP-003). |

`results/` holds what they wrote for BENCHMARKS.md: `compare-identity.txt` and
`compare-scores.txt` (the output of `compare_all.sh`), `giab/` (TP, FN, FP, recall, precision and
F1 per stratification for each call set) and `phasing/` (`whatshap compare --tsv-pairwise` and
`whatshap stats --tsv`, per chromosome).

## What was run

- **lrpipe**: commit `485d01b`, built with `cargo build --release`; the phasing scores used
  `lrpipe run`'s outputs of 2026-11-02. [BENCHMARKS.md](../../BENCHMARKS.md) gives the machine.
- **Reads**: HG003 ONT R10.4.1 sup v5.0.0 aligned with minimap2 (110 GB BAM); HG002 ONT R10.4.1
  sup, ONT's `giab_2025.01` haplotagged CRAM (99 GB); HG002 PacBio Revio HiFi aligned with
  bonefish `-x hifi -Y` (39 GB BAM); HG002 PacBio Revio with PacBio's alignment and 5mC calls
  (`HG002.m84011_220902_175841_s1.GRCh38.bam`, 45 GB). Reference: GRCh38 without ALT contigs.
- **Clair3 models** (SHA-256):

  | Model | `pileup.pt` | `full_alignment.pt` |
  |---|---|---|
  | `r1041_e82_400bps_sup_v500` (ONT) | `d84d9f1a1da0a68204404cac1817724ada4d15363d8a62d87703b6ed56c4a4ff` | `320f7083ae1b8eeba390471fbbeb6fbefb93dab079300ce2605a7222d7d65388` |
  | `hifi` (HiFi) | `71a762cb98e17f163496174bcc35b2ff5606db760383d91c9a885d6637fdc6a0` | `9ae2f5c75184f40649dde5c5bc82d39a65c3bfa5738f516fe20353008f0b4c26` |

- **Original tools**: Clair3 2.0.3, Sniffles 2.8.1, mosdepth 0.3.14, Spectre 0.2.1, TRGT 5.1.0,
  modkit 0.6.4, pb-CpG-tools 3.0.0, WhatsHap 2.8.
- **Truth sets**: GIAB v4.2.1 GRCh38 HG002 and HG003 (`*_1_22_v4.2.1_benchmark.vcf.gz`, with
  `*_benchmark_noinconsistent.bed`); GIAB GRCh38 stratifications (`GRCh38@all`); GIAB HG002 TR
  v1.0.1 Tier 1; GIAB HG002 Q100 v5.0q (`HG002_GRCh38_v5.0q_smvar.vcf.gz` and its
  `benchmark.bed`, from `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/`).
  RTG Tools 3.13.

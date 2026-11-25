# lrpipe 0.2.2 GPU calling, whole genomes

25 November 2026 · the published lrpipe 0.2.2 Linux x86-64 bundle · i5-12600K, 16 threads,
RTX 5060 Ti 16 GB · the inputs and models of [BENCHMARKS.md](../../../../BENCHMARKS.md).

0.2.2 changed GPU calling: the pileup network takes the CPU path's inputs, with their
high-depth scaling, instead of Clair3's 8-bit GPU-mode inputs, which wrap past 127
(#45). These runs check what that did on
human whole genomes.

## Runs

| Run | Wall | CPU time | Peak memory |
|---|---|---|---|
| HG003 ONT, `lrpipe snv --device cuda:0` | 7 min 36 s | 4,052 CPU-s | 11.0 GB |
| HG002 HiFi, `lrpipe snv --device cuda:0` | 4 min 18 s | 2,307 CPU-s | 10.6 GB |

From `/usr/bin/time -v` (`*.time`); [`run.sh`](run.sh) ran them.

## Call agreement

PASS calls (position, alleles, unphased genotype), `lrpipe-snv/tools/concordance.py`:

| Comparison | Before (2026-11-03) | 0.2.2 |
|---|---|---|
| HiFi GPU against its exact CPU run | 99.470% | 99.575% |
| ONT GPU against its exact CPU run | 99.602% | 99.683% |
| HiFi GPU against Clair3 2.0.3's GPU run | 99.707% | 99.340% |
| HiFi GPU, 0.2.2 against before | | 99.343% |
| ONT GPU, 0.2.2 against before | | 99.412% |

## GIAB v4.2.1

[`hg003-ont-gpu.strat.tsv`](hg003-ont-gpu.strat.tsv) and
[`hg002-hifi-gpu.strat.tsv`](hg002-hifi-gpu.strat.tsv): TP, FN, FP, recall, precision and F1 by
stratification, as `score2.sh` writes them. Overall F1:

| | SNP before | SNP 0.2.2 | INDEL before | INDEL 0.2.2 |
|---|---|---|---|---|
| HG003 ONT | 0.99704 | 0.99703 | 0.89213 | 0.89245 |
| HG002 HiFi | 0.99858 | 0.99858 | 0.96881 | 0.96875 |

## Scoring with touchstone

[`score_ts.sh`](score_ts.sh) is `score2.sh` with [touchstone](https://github.com/DylanCosto/touchstone)
0.1.0's `vcfeval` (`-T 8`), a Rust port of RTG Tools 3.13's, in place of RTG's. The split into
SNPs and other variants is still RTG's `vcffilter`, with a 2 GB heap. As a check, touchstone
scored the 2026-11-03 HiFi GPU calls: the stratified table, and the TP, FP, FN and
TP-baseline records of both kinds, were identical to RTG's. It took 26 s and 4.4 GB at most for
the SNPs, and 5 s and 2.0 GB for the rest; RTG ran with a 16 GB heap and reached 19.5 GB.

# HG006 and HG007 ONT indel accuracy by depth

On 2026-11-24, lrpipe and Clair3 completed 20 paired CUDA comparisons on two ONT chr20 panels. The largest paired indel F1 difference between lrpipe and Clair3 was 0.097 percentage points. The median per-case Clair3-to-lrpipe wall-time ratio after GPU admission waits was 15.02×.

## Indel precision, recall, and F1

Downsampled rows give the mean of three replicate scores. The F1 range gives the lowest and highest replicate value. Each native row uses the full staged panel once.

| Donor | Target | Measured depth | lrpipe P / R / F1 | lrpipe F1 range | Clair3 P / R / F1 | Clair3 F1 range | lrpipe mean FP / FN |
|---|---|---:|---:|---:|---:|---:|---:|
| HG006 | Native | 26.44× | 90.86% / 90.74% / 90.80% | 90.80% | 90.91% / 90.67% / 90.79% | 90.79% | 285.0 / 289.0 |
| HG006 | 15× | 15.05–15.13× | 79.72% / 88.37% / 83.82% | 83.55%–84.33% | 79.69% / 88.37% / 83.81% | 83.55%–84.29% | 702.0 / 363.0 |
| HG006 | 20× | 19.91–20.03× | 86.67% / 90.10% / 88.35% | 88.01%–88.53% | 86.52% / 90.07% / 88.26% | 87.92%–88.45% | 432.7 / 309.0 |
| HG006 | 25× | 24.99–25.05× | 89.92% / 90.66% / 90.29% | 90.16%–90.40% | 89.93% / 90.68% / 90.31% | 90.16%–90.45% | 317.3 / 291.3 |
| HG007 | Native | 33.77× | 93.54% / 90.73% / 92.11% | 92.11% | 93.51% / 90.76% / 92.11% | 92.11% | 190.0 / 281.0 |
| HG007 | 15× | 15.03–15.11× | 79.52% / 87.43% / 83.29% | 82.88%–83.77% | 79.52% / 87.43% / 83.29% | 82.88%–83.77% | 682.7 / 381.0 |
| HG007 | 20× | 19.97–20.07× | 86.74% / 89.54% / 88.12% | 87.49%–88.55% | 86.74% / 89.54% / 88.11% | 87.46%–88.56% | 414.7 / 317.0 |
| HG007 | 25× | 24.91–25.05× | 90.35% / 90.53% / 90.44% | 90.25%–90.76% | 90.33% / 90.54% / 90.44% | 90.25%–90.74% | 293.0 / 287.0 |

- HG006 has 3,120 truth indels across 19.253 Mb. Relative to native depth, lrpipe indel F1 changed by -6.98 points at 15× and -0.51 points at 25×.
- HG007 has 3,030 truth indels across 19.232 Mb. Relative to native depth, lrpipe indel F1 changed by -8.82 points at 15× and -1.67 points at 25×.

## Homopolymers

The GIAB homopolymer mask includes exact runs of at least seven bases, imperfect runs of at least eleven bases, and 5 bp of padding. An indel overlaps this mask if its VCF reference span intersects it. [Machine-readable metrics](metrics.json) include per-seed precision, recall, F1, and error counts inside and outside homopolymers, in tandem repeats, and in the GIAB difficult-region union. [GIAB stratifications](https://github.com/usnistgov/giab-stratifications)

| Donor | Target | lrpipe indel precision inside / outside | Clair3 indel precision inside / outside | lrpipe mean indel FPs inside / total |
|---|---|---:|---:|---:|
| HG006 | Native | 76.39% / 98.71% | 76.49% / 98.71% | 259.0 / 285.0 |
| HG006 | 15× | 57.30% / 94.80% | 57.21% / 94.82% | 594.3 / 702.0 |
| HG006 | 20× | 68.24% / 97.60% | 67.98% / 97.55% | 383.7 / 432.7 |
| HG006 | 25× | 74.56% / 98.42% | 74.60% / 98.42% | 285.3 / 317.3 |
| HG007 | Native | 82.87% / 99.32% | 82.80% / 99.32% | 177.0 / 190.0 |
| HG007 | 15× | 57.59% / 95.09% | 57.59% / 95.09% | 587.0 / 682.7 |
| HG007 | 20× | 68.86% / 98.00% | 68.84% / 98.00% | 376.3 / 414.7 |
| HG007 | 25× | 76.24% / 98.66% | 76.18% / 98.67% | 267.3 / 293.0 |

Across the two panels, homopolymer-overlapping loci accounted for 81.9% of the increase in mean lrpipe indel false positives at 15×.

## Indel F1 by length change

Each cell gives lrpipe / Clair3. [The metrics](metrics.json) preserve precision, recall, error counts, and the mixed-size category for multiallelic records spanning more than one length bin.

| Donor | Target | 1 bp | 2–5 bp | 6–15 bp | 16–49 bp |
|---|---|---:|---:|---:|---:|
| HG006 | Native | 87.91% / 87.87% | 94.16% / 94.11% | 96.99% / 97.34% | 100.00% / 100.00% |
| HG006 | 15× | 78.54% / 78.50% | 89.86% / 89.80% | 94.38% / 94.61% | 97.55% / 97.28% |
| HG006 | 20× | 84.30% / 84.21% | 93.07% / 93.00% | 96.09% / 95.98% | 98.91% / 98.64% |
| HG006 | 25× | 87.09% / 87.12% | 93.93% / 93.93% | 97.05% / 97.05% | 99.46% / 99.46% |
| HG007 | Native | 89.83% / 89.92% | 95.47% / 95.38% | 98.18% / 98.18% | 100.00% / 100.00% |
| HG007 | 15× | 77.99% / 77.98% | 90.10% / 90.11% | 94.75% / 94.75% | 99.05% / 99.05% |
| HG007 | 20× | 84.02% / 84.02% | 93.67% / 93.66% | 96.14% / 96.14% | 99.37% / 99.37% |
| HG007 | 25× | 87.09% / 87.10% | 95.09% / 95.04% | 97.23% / 97.29% | 100.00% / 100.00% |

## GPU runtime comparison

Successful-run wall time starts at the indexed BAM and ends with the indexed VCF. Processing time includes each caller’s preparation and phasing. Downsampling, scoring, work from failed attempts, and recorded GPU allocation waits are excluded. The per-run table lists waits separately.

| Donor | Target | Median lrpipe seconds | Median Clair3 seconds | Median paired Clair3 / lrpipe |
|---|---|---:|---:|---:|
| HG006 | Native | 5.95 | 95.14 | 15.99× |
| HG006 | 15× | 6.28 | 76.66 | 12.14× |
| HG006 | 20× | 5.47 | 86.24 | 15.61× |
| HG006 | 25× | 5.84 | 94.62 | 15.99× |
| HG007 | Native | 6.63 | 99.56 | 15.02× |
| HG007 | 15× | 5.53 | 75.29 | 13.54× |
| HG007 | 20× | 5.72 | 85.90 | 15.02× |
| HG007 | 25× | 5.57 | 95.05 | 16.80× |

## Per-run CUDA results and runtimes

| Donor | Arm | Measured depth | lrpipe indel P / R / F1 | Clair3 indel P / R / F1 | lrpipe s | Clair3 s | Clair3 GPU wait s |
|---|---|---:|---:|---:|---:|---:|---:|
| HG006 | Native | 26.44× | 90.86% / 90.74% / 90.80% | 90.91% / 90.67% / 90.79% | 5.95 | 95.14 | 0.00 |
| HG006 | 15× seed 29 | 15.05× | 79.59% / 88.01% / 83.59% | 79.57% / 88.01% / 83.58% | 6.35 | 76.72 | 0.00 |
| HG006 | 15× seed 73 | 15.07× | 79.39% / 88.17% / 83.55% | 79.39% / 88.17% / 83.55% | 6.27 | 76.66 | 0.00 |
| HG006 | 15× seed 113 | 15.13× | 80.19% / 88.91% / 84.33% | 80.12% / 88.91% / 84.29% | 6.28 | 76.27 | 0.00 |
| HG006 | 20× seed 29 | 19.91× | 87.02% / 90.10% / 88.53% | 86.90% / 90.06% / 88.45% | 5.45 | 85.05 | 0.00 |
| HG006 | 20× seed 73 | 20.02× | 86.77% / 90.32% / 88.51% | 86.58% / 90.32% / 88.41% | 5.47 | 87.37 | 0.00 |
| HG006 | 20× seed 113 | 20.03× | 86.22% / 89.87% / 88.01% | 86.09% / 89.84% / 87.92% | 6.38 | 86.24 | 0.00 |
| HG006 | 25× seed 29 | 24.99× | 89.65% / 90.67% / 90.16% | 89.65% / 90.67% / 90.16% | 5.84 | 93.37 | 0.00 |
| HG006 | 25× seed 73 | 25.03× | 90.10% / 90.71% / 90.40% | 90.14% / 90.77% / 90.45% | 5.65 | 95.06 | 94.16 |
| HG006 | 25× seed 113 | 25.05× | 90.01% / 90.61% / 90.31% | 90.01% / 90.61% / 90.31% | 6.47 | 94.62 | 0.00 |
| HG007 | Native | 33.77× | 93.54% / 90.73% / 92.11% | 93.51% / 90.76% / 92.11% | 6.63 | 99.56 | 0.00 |
| HG007 | 15× seed 29 | 15.03× | 79.31% / 87.49% / 83.20% | 79.31% / 87.49% / 83.20% | 5.53 | 74.88 | 0.00 |
| HG007 | 15× seed 73 | 15.11× | 79.00% / 87.16% / 82.88% | 79.00% / 87.16% / 82.88% | 5.44 | 75.29 | 0.00 |
| HG007 | 15× seed 113 | 15.07× | 80.25% / 87.62% / 83.77% | 80.25% / 87.62% / 83.77% | 5.60 | 75.80 | 73.05 |
| HG007 | 20× seed 29 | 19.97× | 87.31% / 89.83% / 88.55% | 87.28% / 89.87% / 88.56% | 5.12 | 84.19 | 0.00 |
| HG007 | 20× seed 73 | 19.99× | 86.27% / 88.75% / 87.49% | 86.24% / 88.71% / 87.46% | 5.93 | 86.48 | 0.00 |
| HG007 | 20× seed 113 | 20.07× | 86.66% / 90.03% / 88.31% | 86.69% / 90.03% / 88.33% | 5.72 | 85.90 | 0.00 |
| HG007 | 25× seed 29 | 24.97× | 90.20% / 90.43% / 90.31% | 90.20% / 90.43% / 90.31% | 5.57 | 93.60 | 0.00 |
| HG007 | 25× seed 73 | 24.91× | 90.23% / 90.26% / 90.25% | 90.21% / 90.30% / 90.25% | 5.46 | 95.46 | 289.41 |
| HG007 | 25× seed 113 | 25.05× | 90.62% / 90.89% / 90.76% | 90.59% / 90.89% / 90.74% | 6.20 | 95.05 | 0.00 |

## Methods and scope

Both callers used an NVIDIA L40S-24Q with eight CPU threads, the same sorted
BAM per pair, and the same `r1041_e82_400bps_sup_v500` model weights. lrpipe used
`snv run --device cuda:0 --read-order canonical`. Clair3 2.0.3 used
`--use_gpu=True --device cuda:0`.
lrpipe source revision

The panels came from flow cells `PAY77227` and `PBA20413` in the R10.4.1 SUP v5.0.0 ONT
release and used `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`. Calls cover chr20:6,000,001–26,000,000.
Native depth and each subsampled depth were measured on that interval with
`samtools coverage --ff 0x904`. For each seed, samtools 1.24 selected read names
with `--subsample-seed` and a fraction equal to target depth divided by native
depth. The three replicate seeds were 29, 73, and 113.
[ONT dataset](https://epi2me.nanoporetech.com/giab-2025.01/)

RTG vcfeval 3.13 compared diploid PASS genotypes with `--decompose`, using
the GIAB v4.2.1 hg38 VCF and confidence BED for each donor. Scoring removed
5 kb from both panel edges and excluded ambiguous reference sequence. The
raw metrics retain separate call-side and truth-side true-positive counts.
Indels have called allele length changes from 1 to 49 bp. Seed variation
measures resampling within each donor. These results describe two chr20
panels.
[GIAB benchmarks](https://www.nist.gov/programs-projects/genome-bottle)
[RTG vcfeval](https://realtimegenomics.github.io/rtg-core/rtg_command_reference.html#vcfeval)

Shared GPU memory contention required one Clair3 retry. A study-local
wrapper waited at Clair3’s allocation checks by polling the actual free
GPU memory. Three completed calls recorded GPU waits totaling 456.62 s.
The sampled peak scratch footprint was 1.55 GiB.

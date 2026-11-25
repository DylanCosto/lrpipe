# Whole-genome HG002 benchmarks

lrpipe-phase was named braid when these were run; its binary was `target/release/braid`.

WhatsHap 2.8 and braid `phase` and `haplotag` on whole-genome HG002 data
from three sequencing technologies, run with `tests/compare_wgs.sh` in October
2026. braid ran with 1, 8 and 24 threads.

**Every output was identical to WhatsHap's**: the phased VCF (except the
`##commandline` line) and read list from `phase`, and the BAM records
(including tag types) and haplotag list from `haplotag`, for all three
technologies and all thread counts.

Raw results are in `pacbio/`, `illumina/` and `ont/`: `compare.txt` (identity
checks) and `timings.txt` (`/usr/bin/time` output). These runs used commit
`ae1f539`, before `phase --max-memory` existed; [Memory limit](#memory-limit)
has `phase` results with it, [noodles I/O](#noodles-io-d057a1c) the checks of
commit `d057a1c`, and [Phasing accuracy](#phasing-accuracy) the comparison
with the GIAB phased truth set.

Memory is the peak resident set size in GiB (2^30 bytes; `/usr/bin/time`
reports KiB).

## Results

Speed-up is WhatsHap's wall time divided by braid's.

### phase

| Data | Usable het. variants | WhatsHap | rs -t 1 | rs -t 8 | rs -t 24 |
|---|---:|---:|---:|---:|---:|
| PacBio Revio | 2,957,781 | 13,386 s (3 h 43 min) | 1,827 s (7.3×) | 241 s (56×) | 132 s (101×) |
| Illumina NovaSeq | 2,443,745 | 13,936 s (3 h 52 min) | 2,676 s (5.2×) | 626 s (22×) | 275 s (51×) |
| ONT R10.4.1 | 2,443,745 | 13,799 s (3 h 50 min) | 2,899 s (4.8×) | 317 s (44×) | 174 s (79×) |

Peak memory:

| Data | WhatsHap | rs -t 1 | rs -t 8 | rs -t 24 |
|---|---:|---:|---:|---:|
| PacBio Revio | 3.4 GiB | 2.0 GiB | 7.4 GiB | 19.0 GiB |
| Illumina NovaSeq | 8.3 GiB | 7.7 GiB | 23.1 GiB | 42.6 GiB |
| ONT R10.4.1 | 2.9 GiB | 2.0 GiB | 6.7 GiB | 22.6 GiB |

### haplotag

The input is WhatsHap's phased VCF from the `phase` run.

| Data | Record data hashed | WhatsHap | rs -t 1 | rs -t 8 | rs -t 24 |
|---|---:|---:|---:|---:|---:|
| PacBio Revio | 268 GB | 10,156 s (2 h 49 min) | 6,604 s (1.5×) | 718 s (14×) | 383 s (27×) |
| Illumina NovaSeq | 314 GB | 21,839 s (6 h 4 min) | 7,227 s (3.0×) | 1,135 s (19×) | 1,112 s (20×) |
| ONT R10.4.1 | 394 GB | 18,048 s (5 h 1 min) | 12,458 s (1.5×) | 1,181 s (15×) | 725 s (25×) |

Peak memory:

| Data | WhatsHap | rs -t 1 | rs -t 8 | rs -t 24 |
|---|---:|---:|---:|---:|
| PacBio Revio | 2.8 GiB | 2.5 GiB | 4.2 GiB | 5.4 GiB |
| Illumina NovaSeq | 11.4 GiB | 7.6 GiB | 13.6 GiB | 14.1 GiB |
| ONT R10.4.1 | 2.5 GiB | 2.6 GiB | 4.8 GiB | 6.2 GiB |

All six WhatsHap runs together took 25.3 hours; braid with 24 threads
took 47 minutes (33×).

### Notes

- `phase` with 1 thread is 5–7× faster than WhatsHap, with similar or lower
  memory. With `-t N`, braid phases up to N/2 chromosomes at a time, so
  memory grows with the thread count; it is highest for Illumina, which has
  the most reads per chromosome. `--max-memory` now limits this (see below).
- `haplotag` memory also grows with the thread count, and it has no limit.
- `haplotag` with 1 thread is 1.5–3× faster. Decoding the input and
  compressing the output BAM likely take most of the time, and both tools have
  to do that work. With more threads, braid also compresses with N
  threads; WhatsHap ran with its default single compression thread, so only
  the 1-thread `haplotag` numbers compare like with like.
- Illumina `haplotag` did not get faster from 8 to 24 threads (below 800% CPU
  in both runs), while PacBio and ONT did. The cause was not investigated;
  writing the larger output (many short records) or disk throughput are
  candidates.
- PacBio used DeepVariant calls, while Illumina and ONT used the GIAB v4.2.1
  benchmark VCF. Speeds across technologies are therefore not strictly
  comparable; the identity results do not depend on this.

## Memory limit

`phase --max-memory MB` (default: half of the system memory, here 60.8 GiB)
phases fewer chromosomes at a time when their estimated memory would not fit.
Results with 24 threads at commit `036632d`; raw numbers are in
`max-memory/timings.txt`.

| Data | Default limit | `--max-memory 12000` (11.7 GiB) | No limit |
|---|---:|---:|---:|
| PacBio Revio | 139 s, 17.6 GiB | 157 s, 12.1 GiB | 136 s, 17.6 GiB |
| Illumina NovaSeq | 275–279 s, 32.9–33.9 GiB | 1,885 s, 8.8 GiB | 255–259 s, 41.1 GiB |
| ONT R10.4.1 | 186 s, 15.9 GiB | 193–200 s, 13.2–13.5 GiB | 185 s, 16.1 GiB |

- Output was identical to WhatsHap in every run that was compared (all but
  the no-limit PacBio and ONT runs), also at `-t 8`.
- The default limit costs Illumina about 8% against no limit and saves 7–8
  GiB; for PacBio and ONT it is not reached.
- With long reads, `--max-memory 12000` is exceeded by up to 15%: the queues
  of CRAM containers being decoded and memory held by the allocator are not
  part of the estimate.
- Illumina chromosomes need up to 7.8 GiB each while their reads are grouped
  by name, so with 12000 MB only one is phased at a time, and one chromosome
  keeps only about two cores busy.
- The "no limit" Illumina runs used commit `aa47266`, whose input files
  queued more blocks for decompression; the other runs share `036632d`'s
  behaviour. Some runs used a build of `036632d` with a thread that logged
  memory once a second (marked in the raw file).

## noodles I/O (`d057a1c`)

Commit `d057a1c` reads and writes all files with noodles instead of htslib and
adds the `stats`, `compare`, `split`, `genotype`, `unphase`, `hapcut2vcf` and
`haplotagphase` subcommands. It was checked against the WhatsHap outputs of the
runs above (WhatsHap did not change). Raw results are in `main-d057a1c/`.

**Every output was identical to WhatsHap's:**

- `phase` with 24 and 8 threads: VCF and read list, all three data sets
- `haplotag` with 24 threads: BAM records and haplotag list, all three
- `stats` (stdout, `--tsv`, `--block-list`) on WhatsHap's phased VCF of each
  data set
- `compare` of the three phased VCFs (stdout, `--tsv-pairwise`,
  `--tsv-multiway`, `--switch-error-bed`, `--longest-block-tsv`); 41,422
  switch-error intervals and 101,893 longest-block rows
- `compare` of each phasing with the GIAB phased truth set (see
  [Phasing accuracy](#phasing-accuracy))

`compare` with the unphased GIAB v4.2.1 VCF was also identical, but assesses no
phased pairs, so it is not a meaningful check. The test suites passed as well
(`phase` at `-t 1`, `-t 4` and `-t 4 --max-memory 1`, `haplotag`, `stats`,
`compare`, `split`, `genotype`, `unphase`/`hapcut2vcf` and the VCF I/O fuzzer).
The 8-thread `haplotag` and `--max-memory 12000` runs were not repeated.

`phase` with 24 threads, each run alone, compared with `036632d` (default
memory limit, [Memory limit](#memory-limit)):

| Data | Input | `036632d` | `d057a1c` |
|---|---|---:|---:|
| PacBio Revio | CRAM | 139 s, 17.6 GiB | 472 s, 13.1 GiB |
| Illumina NovaSeq | BAM | 275–279 s, 32.9–33.9 GiB | 264 s, 24.9 GiB |
| ONT R10.4.1 | CRAM | 186 s, 15.9 GiB | 526 s, 12.3 GiB |

- Peak memory is 23–27% lower (reads are now grouped by name without copying
  them).
- CRAM input is 2.8–3.4× slower. Time spent reading the input (summed over
  threads) went from 423 s to 3,844 s for PacBio and from 999 s to 4,729 s for
  ONT, while the other stages took about as long as before. noodles decodes
  each CRAM file on one thread; htslib decoded containers with a thread pool.
  BAM input is decompressed with several threads and got slightly faster.
- `haplotag` with CRAM input is slower still. PacBio needed 36 minutes for
  chr1–8 while running alone (383 s for the whole genome at `ae1f539`), using
  about 3–4 cores. The `haplotag` and 8-thread `phase` times in
  `main-d057a1c/timings.txt` come from runs that overlapped each other and are
  not comparable with the tables above.
- `stats` and `compare` on these whole-genome VCFs were 3–4× faster than
  WhatsHap's (for example 82 s against 292 s for `compare` of PacBio with the
  GIAB VCF, 122 s against 505 s for the three-way comparison).

## Phasing accuracy

Since the output is identical, braid phases exactly as accurately as
WhatsHap. To show how accurate that is, each data set's phasing (from
braid `phase -t 8` at `d057a1c`, then named whatshap-rs, identical to WhatsHap's) was compared
with the
[GIAB v4.2.1 phased benchmark](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/)
(`HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz`,
phased with trio and Strand-seq data):

```
braid compare --names truth,DATA --tsv-pairwise DATA.tsv TRUTH.vcf.gz DATA.phased.vcf.gz
```

WhatsHap's `compare` gave byte-identical output for each phasing. Totals over
chr1–22, from `main-d057a1c/truth/*.tsv`:

| Data | Pairs assessed | Switch errors | Switch / flip | Block-wise Hamming |
|---|---:|---:|---:|---:|
| PacBio Revio | 1,745,817 | 2,269 (0.13%) | 663 / 803 | 43,551 (2.5%) |
| ONT R10.4.1 | 1,751,380 | 8,810 (0.50%) | 328 / 4,241 | 32,117 (1.8%) |
| Illumina NovaSeq | 1,077,175 | 8,246 (0.77%) | 6,308 / 969 | 14,361 (1.0%) |

- Only variants phased in both files are assessed: the truth set phases about
  70% of its heterozygous variants. The blocks shared by both cover 1.75
  million variants for PacBio and ONT and 1.40 million for Illumina.
- The block-wise Hamming distance counts variants on the wrong haplotype
  within each block, after choosing the block's better orientation; a single
  switch in a long block can move many variants, which is why long reads have
  more of them despite fewer switch errors.
- The PacBio phasing uses its own DeepVariant calls, not the GIAB variants, so
  its assessed variants differ slightly from the other two.

## Setup

**Machine:** VMware VM with 24 vCPUs (2 × Intel Xeon Gold 5418Y, 12 cores
each, no SMT, 2 NUMA nodes), 121 GiB RAM, ext4 on a virtual disk; Ubuntu
26.04.1, Linux 7.0. One benchmark ran at a time, but unrelated jobs (about one
core each) ran on the machine during parts of the queue, so timings carry some
noise.

**Software:** WhatsHap 2.8 (bioconda, Python 3.12.15), samtools/htslib 1.24,
braid 0.1.0 at commit `ae1f539` (then named whatshap-rs, as the raw results
call it), built with `cargo build --release` and rustc 1.99.1.

**Data (all HG002):**

| Data | Reads | Reference | Variants |
|---|---|---|---|
| PacBio Revio | Revio SPRQ, about 46×, aligned CRAM (local run, not a public download) | UCSC hg38 analysis set | DeepVariant small-variant calls, phasing removed with `whatshap unphase` |
| Illumina NovaSeq | [DeepVariant case-study data](https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam): NovaSeq PCR-free 35×, BAM | [GRCh38 no-alt analysis set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) (GCA_000001405.15) | GIAB v4.2.1 benchmark VCF, chr1–22 |
| ONT R10.4.1 | [ONT open data giab_2025.01](https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/wf-human-variation/sup/HG002/PAW70337/output/SAMPLE.haplotagged.cram), PAW70337, sup basecalls, CRAM | GRCh38 no-alt analysis set (GCA_000001405.15) | GIAB v4.2.1 benchmark VCF, chr1–22 |

The ONT CRAM's published `.crai` index has wrong container offsets (region
queries fail with a CRC32 error); the file itself is intact. Rebuild the index
with `samtools index` before use.

To decode the CRAMs without fetching reference sequences from EBI, point
`REF_PATH` and `REF_CACHE` at a local cache built with htslib's
`seq_cache_populate.pl` from the matching reference.

## Reproducing

Both tools ran with `--ignore-read-groups` (the `compare_wgs.sh` default) and
`-r REFERENCE`. For one data set:

```
tests/compare_wgs.sh ENV_BIN REFERENCE VARIANTS.vcf.gz READS.bam WORK_DIR 1 8 24
```

Here `phase` ran for all three data sets first and `haplotag` afterwards, which
`compare_wgs.sh` supports with `SKIP_HAPLOTAG=1` (phase only) followed by
`SKIP_PHASE=1` (haplotag only, appending to the same `WORK_DIR`).

WhatsHap's `haplotag` output is piped into `tests/record_md5.py` and never
written to disk. Each braid output BAM is hashed and then deleted, so
`WORK_DIR` needs room for one output BAM (about 40 GB for the Illumina data)
plus the haplotag lists (45 GB each for Illumina, under 0.5 GB for PacBio and
ONT).

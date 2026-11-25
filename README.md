# lrpipe

An all-in-one germline pipeline for long reads, from FASTQ to VCF.

From one sample's PacBio HiFi or Oxford Nanopore reads, aligned or unaligned, one command calls
phased SNVs and indels, structural variants, copy-number variants and loss of heterozygosity,
tandem repeats, read depth and CpG methylation.

**A whole human genome in 10 minutes, against 1 h 44 min for the original tools run one after
another** (HG002 HiFi, 16 threads, one consumer GPU; [Benchmarks](#benchmarks)).

lrpipe is built on three companion projects, also in Rust:

- [bonefish](https://github.com/DylanCosto/bonefish): the long-read aligner that maps unaligned
  reads, shipped in the release bundle
- [pasta](https://github.com/DylanCosto/pasta): reading and writing of SAM, BAM, CRAM, VCF, FASTA,
  BGZF and their indexes
- [zflate](https://github.com/DylanCosto/zflate): the DEFLATE and gzip codec behind BGZF
  compression

## Quick start

Download the Linux x86-64 bundle from the
[latest release](https://github.com/DylanCosto/lrpipe/releases/latest). It needs glibc 2.28 or
later, and nothing else installed.

```bash
tar -xzf lrpipe-0.3.0-x86_64-linux.tar.gz
export PATH=$PWD/lrpipe-0.3.0-x86_64-linux/bin:$PATH

# Once: reference, Clair3 models and catalogues; writes /ref/lrpipe/site.toml
lrpipe setup --genome GRCh38 --platform ont,hifi --dir /ref/lrpipe --accept-licences

# Per sample
lrpipe run --config /ref/lrpipe/site.toml --reads HG002.bam --sample HG002 \
  --platform hifi --output hg002
```

## Components

| Command | Reimplements | Agreement with the original on whole genomes |
|---|---|---|
| [`lrpipe snv`](crates/lrpipe-snv) | Clair3 2.0.3 | 99.98% of PASS calls (CPU); 99.78% with `--approximate` |
| [`lrpipe sv`](crates/lrpipe-sv) | Sniffles 2.8.1 | identical records, except for Sniffles2 bugs fixed on purpose |
| [`lrpipe cnv`](crates/lrpipe-cnv) | Spectre 0.2.1 | identical VCF and BED |
| [`lrpipe tr`](crates/lrpipe-tr) | none (own method; compared with TRGT 5.1.0) | see Accuracy |
| [`lrpipe depth`](crates/lrpipe-depth) | mosdepth 0.3.14 | identical output |
| [`lrpipe methyl`](crates/lrpipe-methyl) | modkit 0.6.4 `pileup`; pb-CpG-tools 3.0.0 | identical bedMethyl; identical CpG records |
| [`lrpipe phase`](crates/lrpipe-phase) | WhatsHap 2.8 | identical VCFs, read lists and haplotags |

Each command accepts the original tool's arguments and is also built as a standalone binary
(`lrpipe-snv`, `lrpipe-sv`, ...). Each directory's README describes its parity tests and the
deviations from the original, and why each was made.

## Benchmarks

**Setup.** All runs used an Intel Core i5-12600K (16 threads), 31 GB of RAM, an NVIDIA RTX 5060
Ti (16 GB) and inputs on NVMe. Wall time and peak memory are from `/usr/bin/time -v`.
[BENCHMARKS.md](BENCHMARKS.md) gives inputs, CPU time, every comparison, and the scripts to
reproduce them.

**Whole pipeline (`lrpipe run`, all analyses).** These were measured with 0.1.x and have not been
rerun since the 0.2 changes to BGZF compression and the SV and tandem repeat defaults.

| Sample | Input | Wall time | Without final phasing | Peak memory |
|---|---|---|---|---|
| HG002 PacBio Revio HiFi, 31x | 39 GB BAM | 10 min 2 s | 4 min 57 s ¹ | 13 GB |
| HG003 ONT R10.4.1 sup | 110 GB BAM | 11 min 5 s ¹ | 9–10 min ¹ | 14 GB |
| HG002 ONT R10.4.1 sup | 99 GB CRAM | 24 min 19 s | 14 min 53 s ¹ | 14 GB |

¹ Measured with 0.1.5, before haplotagging was added to the pipeline. Haplotagging adds 3–5
minutes ([BENCHMARKS.md](BENCHMARKS.md#haplotagging)).

**Components against the tools they replace**, on the same inputs and threads:

| Command | Job | Original | lrpipe | Speed-up |
|---|---|---|---|---|
| `snv` | HG003 ONT, CPU | Clair3: 83 min | 31 min (`--approximate`), 72 min exact | 2.7x |
| `snv` | HG002 HiFi, GPU | Clair3: 33 min 45 s | 4 min 18 s | 7.9x |
| `sv` | HG003 ONT, BAM | Sniffles2: 126 s | 68 s | 1.9x |
| `sv` | HG002 ONT, CRAM | Sniffles2: 573 s | 72 s | 8.0x |
| `cnv` | HG002 HiFi | Spectre: 148 s | 3.2 s | 46x |
| `tr` | HG002 HiFi, 937,204 loci | TRGT: 645 s | 78 s | 8.3x |
| `depth` | HG002 HiFi, BAM | mosdepth: 63 s | 24 s | 2.6x |
| `methyl pileup` | HG002 ONT, CRAM | modkit: 799 s, 13.5 GB | 251 s, 2.1 GB | 3.2x |
| `methyl cpg` | HG002 Revio | pb-CpG-tools: 1,932 s, 12.5 GB | 217 s, 1.3 GB | 8.9x |
| `phase` | HG002 HiFi, all contigs | WhatsHap (6 contigs in parallel): 1,326 s | 125 s | 10.6x |

On CRAM input, `lrpipe depth` is slower than mosdepth, because CRAM decoding dominates (HG003 ONT
chr20: 8.8 s against 3.3 s).

The total at the top sums the original tools' HG002 HiFi rows: Clair3 on the GPU (2,025 s),
pb-CpG-tools (1,932 s, on PacBio's own alignment of the same reads), WhatsHap (1,326 s), TRGT
(645 s), Spectre (148 s), Sniffles2 (73 s) and mosdepth (63 s): 6,212 s. Against it,
`lrpipe run` took 602 s on the whole genome with every analysis (measured with 0.1.x).

**Apple silicon.** On an M4 Pro Mac mini (12 cores, 24 GB, CPU only), the full HG002 ONT
pipeline took 32 min 35 s. Its outputs matched the original tools run on Linux aarch64, and its
GIAB scores were SNP F1 0.9984 and indel F1 0.9142. On the Apple GPU (`--device metal`),
`lrpipe snv` ran 16.7x (HiFi) and 23.9x (ONT) faster than Clair3 over a 2 Mb interval
([details](BENCHMARKS.md#lrpipe-snv-on-an-apple-gpu-metal)).

## Accuracy

**SNVs and indels.** Scored against GIAB v4.2.1 with RTG `vcfeval`, on PASS calls within the
benchmark regions of chr1–22:

| Sample | Caller | SNP F1 | Indel F1 |
|---|---|---|---|
| HG003 ONT | Clair3 2.0.3, CPU | 0.99704 | 0.89236 |
| | lrpipe, CPU | 0.99704 | 0.89231 |
| | lrpipe, GPU | 0.99703 | 0.89245 |
| HG002 HiFi | Clair3 2.0.3, GPU | 0.99857 | 0.96875 |
| | lrpipe, CPU | 0.99857 | 0.96879 |
| | lrpipe, GPU | 0.99858 | 0.96875 |

Every lrpipe mode is within 0.00015 of Clair3's overall F1, and within 0.0004 in every GIAB
stratification.

On the GPU, lrpipe uses the CPU path's network inputs. Clair3's GPU mode stores them in 8 bits,
which overflow at high depth: on *E. coli* K-12 HiFi at 295x, reproducing that mode gave SNP F1
0.43 against 0.96 on the CPU. As a result, lrpipe's GPU calls agree with its own CPU calls
(99.58% on HG002 HiFi) more than with Clair3's GPU run (99.34%), at the same F1
([details](BENCHMARKS.md#gpu-calling-in-022)).

**Other variant types:**

| Analysis | Truth set | Result |
|---|---|---|
| Phasing | GIAB HG002 Q100 v5.0q | switch error 0.31% (HiFi), 0.40% (ONT) |
| SVs (accuracy mode) | GIAB HG002 v5.0q | about 1,600 more true SVs found than Sniffles2 ([BENCHMARKS.md](BENCHMARKS.md#structural-variants-giab-hg002-v50q)) |
| Tandem repeats | GIAB HG002 TR v1.0.1, Tier 1 | both allele lengths exact at 90.6% of called regions (TRGT 5.1.0: 87.6%); measured with 0.1.x |
| Repeat expansions | 33 Coriell cell lines (PacBio PureTarget) | repeat counts matched Coriell in 24 of 38 samples (TRGT: 22) |
| Mendelian consistency | HG002/HG003/HG004 trio, ONT | [BENCHMARKS.md](BENCHMARKS.md#mendelian-consistency-the-hg002-hg003-hg004-trio-ont) |

## Install

The release bundle contains:

- `lrpipe` and each tool's standalone binary
- bonefish 0.5.2, the aligner used for unaligned reads
- the licences of every component compiled in

Clair3's models are either fetched by `lrpipe setup` or downloaded separately from the release,
as `lrpipe-VERSION-models.tar.gz` (167 MB). **ONT's models are licensed for research use with
Nanopore data only.**

`lrpipe setup` writes a site configuration for GRCh38, CHM13 or GRCh37. With `--offline`, it
prints the downloads for a machine without internet access instead
([`lrpipe/`](crates/lrpipe#setting-up)).

To build from source, you need Rust 1.91 or later and a C compiler:

```bash
cargo build --release                         # with NVIDIA GPU support
cargo build --release --no-default-features   # CPU only
cargo build --release --features metal        # with Apple GPU support (macOS)
```

## License

MIT ([LICENSE](LICENSE)), except lrpipe-snv, which is BSD-3-Clause like Clair3, the tool it ports.
Each tool's directory lists the licences of the projects it ports.

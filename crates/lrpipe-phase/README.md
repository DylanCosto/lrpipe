# lrpipe-phase

`lrpipe phase` (and `lrpipe haplotag`): an unofficial Rust port of [WhatsHap](https://github.com/whatshap/whatshap)
(read-based phasing) by Marcel Martin, Tobias Marschall and contributors. The
goal is **byte-identical output to WhatsHap 2.8**, just faster. It is not
affiliated with or endorsed by the WhatsHap project; WhatsHap's license is in
`LICENSE-WHATSHAP`.

It also builds as its own binary, `lrpipe-phase`, with WhatsHap's command line
(`lrpipe-phase phase`, `lrpipe-phase haplotag`, ...); in `lrpipe`, `lrpipe phase` implies
`phase` unless another of its commands is named, and `lrpipe haplotag` is `haplotag`.
lrpipe-snv phases with it in process.

Formerly named braid, and whatshap-rs before that, in its own repository. Commits quoted in
these documents are from its development history, which this repository does not include.

## Status

`lrpipe-phase phase` supports:

- BAM/CRAM input (one or more files, merged by position), read groups / `--sample`,
  `--ignore-read-groups`
- phased VCF/BCF files as additional phasing input (GT/PS or HP phasing)
- allele detection by re-alignment (`--reference`) or from the CIGAR (`--no-reference`)
- read selection (`--internal-downsampling`), the exact MEC/PedMEC dynamic program,
  `--merge-reads`, and `--algorithm heuristic` (with `--row-limit`) and
  `--algorithm hapchat`
- pedigree phasing: `--ped`, `--use-ped-samples`, `--recombrate`, `--genmap`,
  `--no-genetic-haplotyping`, `--recombination-list`
- `--distrust-genotypes` (GL or PL likelihoods), `--include-homozygous`,
  `--default-gq`, `--gl-regularizer`, `--changed-genotype-list`
- `--tag PS` or `--tag HP`, `--only-snvs`, `--mapq`, `--chromosome`,
  `--exclude-chromosome`, `--use-supplementary`, `--supplementary-distance`,
  `--output-read-list`
- `--threads N`: phase up to N/2 chromosomes at a time (fewer if they would
  exceed `--max-memory`), and use N threads for
  allele detection, for solving independent parts of a chromosome (single
  samples; parts that no read connects) and for decompressing the input. Output
  is identical for any N.
- `--dp-memory MB` (lrpipe-phase only, default 1024): memory the phasing DP may use
  to keep its backtrace tables instead of recomputing most columns from
  checkpoints, as WhatsHap does. Faster, but does not change the output.
- `--max-memory MB` (lrpipe-phase only, default: half of the system memory or of
  the cgroup limit, e.g. in a container or cluster job): approximate memory limit.
  With `--threads`, a chromosome is only started while the estimated memory of
  all chromosomes in progress fits next to `--dp-memory` and 96 MB per thread;
  one is always phased. The estimate comes from the memory per variant of the
  chromosomes read so far, since most memory goes to a chromosome's reads. Does
  not change the output.

`lrpipe-phase haplotag` tags reads by haplotype (HP/PC/PS tags):

- BAM/CRAM input and output (`-o x.cram` with `--reference`), indexed VCF.gz (tabix
  or CSI) or BCF input, `--regions`, `--sample`, `--ignore-read-groups`
- allele detection with `--reference` or `--no-reference`, `--ploidy`
- linked reads (BX tags, `--linked-read-distance-cutoff`, `--ignore-linked-read`)
- all `--tag-supplementary` strategies, `--supplementary-distance`,
  `--no-supplementary-strand-match`, `--skip-missing-contigs`
- `--output-haplotag-list` (gzipped if it ends in `.gz`)
- `--threads N`: prepare haplotype assignments for the next two chromosomes while
  one is written, detect alleles and decompress the input with N threads, and
  (unless `--output-threads` is given) compress the output with N threads. WhatsHap compresses with one thread by
  default. The records written are identical for any N.
- `--compression-level 0-9` (lrpipe-phase only): trade output size for speed. Only the
  compressed bytes change, not the records.
- `--list-only` (lrpipe-phase only): write the haplotag list alone, no alignments, which
  `lrpipe run` gives its SV phasing and CpG tables in place of a second BAM. The list is the
  one written beside the alignments.
- `--write-index` (lrpipe-phase only): index the output BAM as `samtools index` does
  (`.bai`, or `.csi` for contigs too long for a BAI).
- `--untag-missing-contigs` (lrpipe-phase only): write the reads on a contig the VCF does not
  have untagged, as on one without variants, where WhatsHap fails or, with
  `--skip-missing-contigs`, leaves them out of the output (lrpipe-snv's VCFs need not name
  every contig of the reads).
- `--regions-bed BED` (lrpipe-phase only), in place of `--regions`: the regions of a BED, sorted
  and merged, which `lrpipe run --regions` passes (as many as a panel has, which a command line
  would not hold). A read reaching two regions is written and listed once, in order, where
  WhatsHap's `--regions` writes it again for each region it reaches, after reads that start
  later (an unsorted BAM); haplotypes are WhatsHap's for the same regions.

The `@PG` header line names `lrpipe-phase` (ID/PN/VN/CL) instead of `whatshap`;
everything else in the output BAM matches WhatsHap, including BAM tag types.

`lrpipe-phase stats` prints phasing statistics of a VCF/BCF (blocks, block sizes
and lengths, NG50), with all options: `--tsv`, `--block-list`, `--gtf`,
`--sample`, `--chromosome` (looked up through the index if there is one),
`--chr-lengths` and `--only-snvs`. Output and files are identical to WhatsHap.

`lrpipe-phase compare` compares two or more phasings (switch errors, switch/flip
decomposition, Hamming distance, block intersections, multiway comparison),
for diploid and polyploid samples (`-p`), with `--names`, `--sample`,
`--ignore-sample-name`, `--only-snvs`, `--tsv-pairwise`, `--tsv-multiway`,
`--switch-error-bed` and `--longest-block-tsv`. Output and files are identical
to WhatsHap; the polyploid switch/flip calculator breaks ties exactly as
WhatsHap's C++ code does. Not supported: `--plot-blocksizes` and
`--plot-sum-of-blocksizes`, which make PDFs with matplotlib.

`lrpipe-phase split` splits reads by haplotype according to a haplotype list
(such as `haplotag --output-haplotag-list`): FASTA/FASTQ (plain or gzipped) or
BAM input, `--output-h1/--output-h2` or `-o` per haplotype, `--output-untagged`,
`--add-untagged`, `--only-largest-block`, `--discard-unknown-reads` and
`--read-lengths-histogram`. FASTQ outputs and lists can be compressed with
gzip, bzip2 or xz (by extension, as WhatsHap's xopen does; zstandard is not
supported). The records and the header of BAM outputs and the decompressed
contents of all other outputs are identical to WhatsHap's.

`lrpipe-phase genotype` genotypes variants like WhatsHap: priors from all reads,
then the forward-backward algorithm over phase-informative reads, for single
samples and pedigrees (`--ped`, `--genmap`, `--use-ped-samples`), with
`--prioroutput`, `--no-priors`, `--constant`, `--gt-qual-threshold`,
`--affine-gap` and `--use-kmerald` re-alignment and phased VCF blocks as
input. WhatsHap computes the DP in 80-bit `long double`, and the last bits
show in the GL of confident genotypes, so lrpipe-phase emulates x87 extended
precision exactly (`compat::f80`, checked against the hardware) and takes the
powers of ten from glibc's `powl`, which is not correctly rounded
(`core::powl_table`). Only two rare cases are approximations: families with
two or more trios and phred scores or recombination costs of 1024 and more
(glibc's general `powl` algorithm); the GL can then differ in the last digit.
`-t N` (lrpipe-phase only) genotypes N chromosomes at a time; the output does
not change. The software extended precision makes the DP slower than
WhatsHap's C++: on simulated two-sample data (two chromosomes), one thread took
200 s against WhatsHap's 133 s, two threads 118 s. A single chromosome is
about 1.5x slower than with WhatsHap; whole genomes are faster with `-t`.

`lrpipe-phase haplotagphase` phases variants from the HP/PS tags of haplotagged
reads (votes per phase set, `--gap-threshold`, `--cut-poly`, `--only-indels`,
multi-allelic sites unless `--no-mav`).

`lrpipe-phase unphase` removes phasing (PS, HP, PQ; sorted unphased genotypes),
also from standard input.

`lrpipe-phase hapcut2vcf` combines HapCUT/HapCUT2 output with the VCF. As in
WhatsHap, the output is uncompressed whatever the file name, and records after
the last chromosome in the HapCUT output are not written.

Not ported yet: the other subcommands (`polyphase`, `polyphasegenetic`,
`learn`, `find_snv_candidates`).

Like WhatsHap (pysam), lrpipe-phase rejects VCFs compressed with plain gzip
instead of bgzip.

Unlike WhatsHap, it phases a chromosome holding haploid and diploid calls, as chrX of an XY
sample does with haploid calls outside the pseudo-autosomal regions: the haploid calls are left
out of phasing and written unchanged. WhatsHap stops with "Inconsistent ploidy (2 and 1)".
They are left out sample by sample: in a VCF of several samples, a man's haploid call leaves
another sample's diploid call at the same site in its phasing.

### Known cases where WhatsHap itself is not deterministic

- `phase --algorithm hapchat` phases only the first block of overlapping reads,
  but writes super reads for the positions of all reads; for positions after
  the first block, WhatsHap reads past the end of its haplotype buffers. Usually
  the result is a homozygous "genotype change" and WhatsHap fails (as
  lrpipe-phase does), but with `--distrust-genotypes` the alleles written there
  depend on uninitialized memory; lrpipe-phase uses allele 0.

- `haplotag` with several samples whose reads share read names (in one BAM):
  WhatsHap lets the last processed sample win, and iterates samples via a Python
  `set` of strings, whose order depends on the per-process hash seed. Different
  WhatsHap runs then disagree with each other; lrpipe-phase processes samples in
  sorted order (matching one of WhatsHap's outcomes).
- If more than one INFO field is missing from the VCF header, WhatsHap adds the
  definitions in hash-seed-dependent order; lrpipe-phase uses first-seen order.
- `phase --use-ped-samples`: WhatsHap takes the samples from a Python `set`, so the
  order of family members (and with it, tie-breaking in the DP) depends on the
  hash seed; lrpipe-phase uses sorted order.
- `genotype` takes the samples from a Python `set`; with `--ped`, the order of
  family members (and with it the order of floating-point sums in the DP)
  depends on the hash seed. lrpipe-phase uses the order of the VCF.
- `compare --tsv-multiway --ignore-sample-name` with differently named samples:
  WhatsHap joins the names from a Python `set` (hash-seed-dependent order);
  lrpipe-phase uses the order of the VCFs.

WhatsHap overwrites `--recombination-list` and `--changed-genotype-list` for every
chromosome (and, for the recombination list, for every family), so the files only
contain the last one. lrpipe-phase does the same.

## Benchmarks

Whole-genome HG002 (PacBio Revio, Illumina NovaSeq 35×, ONT R10.4.1) on a
24-core machine; output identical to WhatsHap 2.8 in every run. Speed-up over
WhatsHap's wall time:

| Command | Threads | PacBio | Illumina | ONT |
|---|---:|---:|---:|---:|
| `phase` | 1 | 7.3× | 5.2× | 4.8× |
| `phase` | 24 | 101× | 51× | 79× |
| `haplotag` | 1 | 1.5× | 3.0× | 1.5× |
| `haplotag` | 24 | 27× | 20× | 25× |

Details, memory use, data sources and how to reproduce:
[`benchmarks/hg002/`](benchmarks/hg002/README.md).

The table predates three later speed-ups, each with the same output: read selection indexes
its priority queue and variants by position instead of hashing them; a single individual's DP
columns are computed with the partitionings in binary order, vectorized, with ties broken by Gray
code rank as WhatsHap's order breaks them; and BGZF went through libdeflate. On HG002 Revio chr1
(`phase --distrust-genotypes -t 2`, as Clair3 runs it) they took lrpipe-phase from 31.9 s to 16.8 s:
read selection 9.7 s to 2.5 s, the DP 10.5 s to 5.1 s, reading the BAM 9.8 s to 7.3 s. These times
are with libdeflate; pasta 0.9.0 replaced it with zflate, in Rust, and they have not been measured
since.

## How identical output is achieved

WhatsHap's results depend on a few implementation details beyond the algorithm:

- Reads starting at the same variant are ordered by libstdc++'s
  `std::hash<std::string>` (`ReadSet::sort`).
- Read selection breaks ties via the iteration order of CPython `set`s and a
  C++ `unordered_set`.

`src/compat/` emulates both bit-for-bit (CPython 3.12 `setobject.c`, libstdc++'s
`_Hash_bytes` and hashtable layout). The emulations are fuzz-tested against the
real implementations (`tests/compat/`).

All file I/O is on pasta (in Rust, with zflate for BGZF), a library begun from
[noodles](https://github.com/zaeleus/noodles) (BGZF, BAM, CRAM, CSI/BAI/tabix
indexes, FASTA). WhatsHap reads and writes through pysam, i.e. htslib 1.23, so
where htslib's behavior shows in the output it is ported on top of pasta:
VCF/BCF headers, record parsing and formatting (`src/hts/`, including htslib's
number formatting and its handling of undeclared tags), BAM record decoding and
encoding, and the MD/NM tags htslib adds when decoding CRAM (`src/alignment/`).

## Building

Part of the lrpipe workspace; requires Rust 1.91 or later (developed with 1.98), and a C compiler
for the mimalloc allocator, which the `lrpipe-phase` binary uses by default to keep its peak memory
down with many threads (glibc's allocator keeps freed tables reserved per thread). Without mimalloc
(`--no-default-features`) the output is the same; a program using lrpipe-phase as a library chooses
its own allocator, so it can depend on it with `default-features = false`.

From lrpipe's directory:

```
cargo build --release -p lrpipe-phase
./target/release/lrpipe-phase phase -r ref.fa -o phased.vcf variants.vcf reads.bam
```

## Testing against WhatsHap

You need WhatsHap 2.8 and samtools (e.g. `micromamba create -p env -c conda-forge -c bioconda whatshap=2.8 samtools minimap2`)
and a checkout of the WhatsHap repository at tag `v2.8`.

- `tests/compare.sh WHATSHAP SAMTOOLS WHATSHAP_SRC WORKDIR` runs both tools'
  `phase` on the WhatsHap test data (single samples, pedigrees, phased VCF inputs,
  distrusted genotypes, HP tags) and diffs the VCFs, read lists, recombination
  lists and changed-genotype lists. The `##commandline` header line is excluded
  because it names the tool.
- `tests/compare_haplotag.sh ENV_BIN WHATSHAP_SRC WORKDIR` does the same for
  `haplotag`: it compares the output BAMs record by record including BAM tag
  types (`tests/dump_bam.py`, needs pysam), the header except the tool's own
  `@PG` line, and the haplotag lists.
- `tests/compare_stats.sh ENV_BIN WHATSHAP_SRC WORKDIR [SEEDS [FIRST_SEED]]`
  compares `stats` (standard output, exit code, `--tsv`, `--block-list` and
  `--gtf` files) on every VCF of the WhatsHap test data and on random phased
  VCFs from `tests/stats/gen.py` (interleaved blocks, PS of every type, HP
  tags, missing values), plain and bgzipped with an index.
- `tests/compare_compare.sh ENV_BIN WHATSHAP_SRC WORKDIR [SEEDS [FIRST_SEED]]`
  compares `compare` (standard output, exit code and all TSV/BED files) on the
  WhatsHap test data and on random sets of two to four phasings of the same
  truth from `tests/compare/gen.py` (ploidy 2 to 4, switch and flip errors,
  changed, missing and unphased genotypes, multi-ALT sites, PS types).
- `tests/compare_split.sh ENV_BIN WHATSHAP_SRC WORKDIR [SEEDS [FIRST_SEED]]`
  compares `split` (exit code, decompressed outputs and the histogram) on the
  WhatsHap test data and on random FASTQ/BAM inputs and haplotype lists from
  `tests/split/gen.py` (multi-line, FASTA and CRLF records, unaligned BAMs with
  padded headers, compressed lists, every option).
- `tests/compare_genotype.sh ENV_BIN WHATSHAP_SRC WORKDIR` compares `genotype`
  (output and `--prioroutput` VCFs) on the WhatsHap test data: single samples,
  pedigrees, priors options, re-alignment modes, phased VCF inputs.
- `tests/compare_unphase_hapcut2vcf.sh ENV_BIN WHATSHAP_SRC WORKDIR [SEEDS]`
  compares `unphase` (every test VCF, also from standard input) and
  `hapcut2vcf` on random VCFs and HapCUT 1/2 files from `tests/hapcut2vcf/gen.py`.
- `tests/core/genotype_oracle.py SEED COUNT` + `examples/genotype_replay.rs`:
  random read sets and pedigrees through WhatsHap's compiled `GenotypeDPTable`
  and `compute_genotypes`; the likelihoods must be identical bit for bit.
- `tests/compat/f80_oracle.c` + `examples/f80_replay.rs`: x87 results of random
  extended precision operations against `compat::f80`.
  `tests/compat/powl_table.c` generates `src/core/powl_table.rs`.
- `tests/core/switchflip_oracle.py SEED COUNT` + `examples/switchflip_replay.rs`:
  random polyploid phasing pairs through WhatsHap's compiled
  `SwitchFlipCalculator` and through `core::switchflip`; totals and the
  per-position path must be identical.
- `tests/compare_wgs.sh ENV_BIN REF VCF BAM WORK_DIR [THREADS...]`: runs both
  tools on real data (e.g. a whole genome) and checks that phase (VCF, read list)
  and haplotag (BAM records hashed by `tests/record_md5.py`, haplotag list) are
  identical, with timings. `SKIP_HAPLOTAG=1` runs only the phase half and
  `SKIP_PHASE=1` only the haplotag half.
- `tests/vcf/fuzz.sh ENV_BIN WORK_DIR [FIRST_SEED [COUNT]]`: VCF I/O oracle.
  Random VCFs (all INFO/FORMAT types and numbers, missing values, undeclared
  tags, every integer width, floats in many spellings) and random phasing edits
  are replayed through pysam with WhatsHap's `PhasedVcfWriter`
  (`tests/vcf/replay.py`) and through lrpipe-phase (`examples/vcf_replay.rs`);
  the values read, the edited records, region queries and the output VCF must
  be identical. Build the replayer with
  `cargo build --release -p lrpipe-phase --example vcf_replay`.
- `tests/core/oracle.py SEED OUT` + `examples/core_replay.rs`: random
  differential tests of read sorting, read selection and the DP (single
  samples and trios).
- `tests/compat/`: fuzzers for the CPython/libstdc++ emulations.
- Profiling: `cargo build --profile profiling -p lrpipe-phase --features profile`, then run with
  `LRPIPE_PHASE_PROFILE=out` to get `out.folded` (stacks for `inferno-flamegraph`) and
  `out.txt` (top functions); `haplotag --debug` logs per-chromosome prepare/write times.
- `tests/sim/simulate.py`: simulates a diploid genome with ONT-like reads, optionally
  with chimeric reads, unmappable reads, segmental duplications, a trio with
  crossovers (`--trio`, writes a PED file and a genetic map), genotype-call errors
  and GL fields (align with minimap2) for larger end-to-end comparisons and
  benchmarks.

## License

WhatsHap is MIT licensed; see `LICENSE-WHATSHAP`. This port is derived from it.

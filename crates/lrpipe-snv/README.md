# lrpipe-snv

An unofficial Rust port of [Clair3](https://github.com/HKU-BAL/Clair3) v2.0.3, and
[lrpipe](..)'s SNV and indel caller (`lrpipe snv`, or the `lrpipe-snv` binary). It runs Clair3's
models with its own networks in Rust, on the CPU or an NVIDIA GPU, without libtorch or Python.
In comparisons using the same read order, its PASS calls agree with Clair3's on 99.8-100%
([Networks](#networks)). An opt-in [approximate mode](#approximate-mode) trades a little of that
agreement for speed. Not affiliated with the Clair3 authors; see
[Acknowledgements](#acknowledgements).

## Status

Clair3's calling pipelines are ported with their options: the C implementation that
`run_clair3.py` uses for ONT and PacBio HiFi, and the Python pipeline it uses for Illumina
(see [Illumina](#illumina)): `lrpipe-snv run` goes from the BAM to
`merge_output.vcf.gz` (and with `--gvcf`, `merge_output.gvcf.gz`) in one process, with no
Python, PyPy, GNU parallel, bgzip, tabix, lz4 or WhatsHap, and writes Clair3's output layout
(`pileup.vcf.gz`, `full_alignment.vcf.gz`, `merge_output.vcf.gz`, and `tmp/` with every
intermediate) so each step can be diffed with Clair3. `lrpipe-snv pileup` runs a single
`CallVariantsFromCffi --pileup` job (its `--threads` is accepted for Clair3's job command lines,
but the networks score on one thread) and `lrpipe-snv phase` a single phasing job.

| stage | Clair3 | crate |
|---|---|---|
| 1. pileup calling | `calculate_clair3_pileup`, `CreateTensorPileup`, `Clair3_P`, `CallVariants` | `snv-pileup`, `snv-infer`, `snv-call` |
| merge chunks | `SortVcf` | `snv-merge` |
| 2. het SNPs for phasing | `SelectQual --phase`, `SelectHetSnp` | `snv-select` |
| 3. phasing | `whatshap phase --distrust-genotypes --ignore-read-groups` | `snv-phase` (lrpipe-phase) |
| 5. full-alignment candidates | `SelectQual`, `SelectCandidates` | `snv-select` |
| 6. full-alignment calling | `calculate_clair3_full_alignment` (with haplotagging), `CreateTensorFullAlignment`, `Clair3_F`, `CallVariants` (indel-length heads) | `snv-fullalign`, `snv-infer`, `snv-call` |
| 7. merge | `SortVcf`, `MergeVcf`, `SortVcf` | `snv-merge` |
| gVCF (`--gvcf`) | `compute_PL`, `variantInfoCalculator`, `SortVcf`, `gvcfGenerator.mergeCalls` | `snv-pileup`, `snv-call`, `snv-merge` |

[VALIDATION.md](VALIDATION.md) has the GIAB validation, made while lrpipe-snv ran the networks
through libtorch (at commit `282d716`):

- **Parity:** with libtorch, every output file was the same as Clair3's on HG002 ONT, Revio HiFi
  and Illumina (chr20), and on the whole genome for HiFi and ONT on the CPU and ONT on the GPU.
- **Accuracy:** hap.py against GIAB v4.2.1. On the whole HG003 genome, ONT sup v5.0.0 scored
  SNP F1 0.9973 and INDEL F1 0.9042, and Illumina 35x scored SNP F1 0.9949 and INDEL F1 0.9942.

Since libtorch was removed, the probabilities come from lrpipe-snv's own networks and differ from
libtorch's in the last digits, so the output is no longer byte-identical to Clair3's. Its PASS
calls agree with Clair3's on 99.98% of the whole HG003 ONT genome on the CPU, on 100% of HG002
HiFi chr20 on the CPU and 99.47% on the GPU, and on 100% of every option case Clair3 was run
with on ONT, HiFi and Illumina regions, with the same GIAB results where measured
([Against Clair3](#against-clair3)).

### Options

`lrpipe-snv run` takes `run_clair3.py`'s options, written with dashes (`--bed-fn`,
`--snp-min-af`, ...) or as Clair3 writes them, so a `run_clair3.sh` command line works with
`lrpipe-snv run` in place of `run_clair3.sh`:
- **Spellings:** underscores (`--bam_fn=...`), the short options `-b -f -m -t -p -o`, and
  unique prefixes, as Python's argparse allows.
- **Boolean options:** `--flag`, `--flag=True` or `--flag False`.
- **`EMPTY`:** the value `EMPTY`, Clair3's "not given", for `--bed_fn`, `--vcf_fn`,
  `--ctg_name`, `--longphase`, `--device` and `--par_regions_bed`.
- **Tool paths:** `--samtools`, `--python`, `--pypy`, `--parallel` and `--whatshap` are
  accepted and ignored, since lrpipe-snv runs those steps in process, as is
  `--use_whatshap_for_intermediate_phasing`, which Clair3's pipeline does not read either.
- **A link named `run_clair3.sh`:** invoked through a link of that name (or `run_clair3.py`),
  lrpipe-snv runs `lrpipe-snv run` with `--read-order input` and `--ambiguous-bases n` by
  default: Clair3's order of reads that start at one position, and its `N` for an ambiguous
  reference base (which `lrpipe snv run` writes as the VCF specification has it, the first
  nucleotide the code stands for).

For the native ONT and PacBio HiFi pipelines, `lrpipe run` and `lrpipe snv run` use
`--read-order canonical`. Reads with the same alignment start are ordered by strand, read name
and the alignment fields used to build the tensors; more than 64 of them (an amplicon's or a
targeted capture's reads) by a hash of the read name first, so that the reads taken up to a
depth are of both strands. This gives pileup and full-alignment calling
the same read order when merging files or rewriting ties. Use `lrpipe snv run --read-order input`
to compare with a Clair3 run on the original file order.

Phasing and haplotagging still read tied reads in the file's order (they are WhatsHap's steps,
whose read selection follows it), so the phased VCF and what is phased from it are the same for
two files only where their tied reads are in the same order.

Every stored option case gives the same files when run with Clair3's own command line through
such a link. The options are resolved as `run_clair3.py` does:
- platform defaults for the AF and full-alignment proportions;
- `--min-mq >= 5` and `--min-coverage >= 2`;
- AF 0 when genotyping;
- the guppy5 phasing proportion;
- haplotagging implying final-output phasing.

`--bam`, `--ref`, `--model-path` and `--output` are required. `--threads` defaults to 8 and
`--platform` to `ont`. `--approximate` (lrpipe-snv only) gives up a little agreement for speed;
see [Approximate mode](#approximate-mode).

`tools/option_parity.sh NAME PLATFORM BAM REF MODEL CONTIGS [OPTION...]` runs both tools with
the same options and compares every file Clair3 writes (VCFs, gVCFs, BEDs, `CONTIGS`,
`CHUNK_LIST`, split BEDs, the haplotagged BAM's records), ignoring the header lines that name
the command. On 6 Mb fixtures (ONT HG003 chr20:30-36 Mb, Revio HiFi chr22:20-26 Mb, and ONT
chrX/chrY 2.5-4.5 Mb across the PAR1 boundary), with libtorch, every case below wrote the same
files as Clair3. With lrpipe-snv's own networks every file that depends on the networks' output
can differ, and `tools/concordance.py` compares the VCFs' PASS calls:

| fixture | options |
|---|---|
| ONT | defaults; `--bed_fn`; `--snp_min_af 0.1 --indel_min_af 0.2 --min_mq 10 --min_coverage 4 --qual 10 --var_pct_full 0.5 --ref_pct_full 0.2 --var_pct_phasing 0.6`; `--print_ref_calls --haploid_precise`; `--haploid_sensitive`; `--call_snp_only`; `--fast_mode`; `--enable_variant_calling_at_sequence_head_and_tail`; `--enable_long_indel`; `--no_phasing_for_fa`; `--pileup_only --print_ref_calls`; `--enable_phasing --use_whatshap_for_final_output_haplotagging`; `--chunk_size 1000000 --chunk_num 4`; `--vcf_fn` (1,053 sites, 774 not called); `--min_contig_size` above the contig; `--qual 0 --print_ref_calls`; `--gvcf`; `--use_longphase_for_intermediate_phasing`; `--use_longphase_for_final_output_phasing`; `--chunk_num -1`; `--chunk_num -1 --bed_fn --gvcf` |
| ONT with synthetic move tables (`tools/add_synthetic_moves.py`) | `--enable_dwell_time` with `r1041_e82_400bps_sup_v520_with_mv` (two move tables) and `r1041_e82_400bps_hac_v600_with_mv`, on the CPU; with `--use_gpu`, against Clair3's GPU run |
| ONT chrX/Y, no `--ctg_name` | `--gender male`; `--gender male --par_regions_bed`; `--gender female` |
| HiFi | defaults; `--print_ref_calls --snp_min_af 0.1 --no_phasing_for_fa`; `--longphase_for_phasing` |
| ONT, `--disable_c_impl` | defaults; `--fast_mode`; `--use_longphase_for_intermediate_phasing --use_longphase_for_final_output_phasing`; `--gvcf --print_ref_calls --enable_long_indel` |
| HiFi, `--disable_c_impl` | defaults; `--no_phasing_for_fa --snp_min_af 0.1 --enable_variant_calling_at_sequence_head_and_tail` |

LongPhase is the one external tool: like Clair3, lrpipe-snv runs the `longphase` binary (2.1.1
here; `--longphase PATH`).

Unlike Clair3, lrpipe-snv picks the full-alignment candidates before phasing (they depend on the
pileup calls only) and runs phasing alongside full-alignment calling: a contig's jobs start
once it is phased, so the two stages overlap and the contig's reads are read again while still
cached. The files written are the same.

Not ported: several GPUs (Clair3 can spread its GPU processes over several).
Contigs outside chr1-22/X/Y keep their reference order, where Clair3's order comes from a
Python set and changes between runs.

### Parity

The measurements in this section were made while lrpipe-snv ran the networks through libtorch,
the same build as Clair3's. What comes before the networks (tensors) and after them (calling,
selection, phasing, merging) is unchanged since, but the probabilities now come from
lrpipe-snv's own [networks](#networks), so a few calls near a cut-off change and the files are no
longer byte-identical.

End to end, `lrpipe-snv run --ctg-name chr20 --threads 16` on HG003 chr20 (ONT R10.4.1 sup
v5.0.0, model `r1041_e82_400bps_sup_v500`) wrote the same files as `run_clair3.sh` with the
same libtorch, apart from the `##cmdline` header lines that name the command:

- the 13 pileup chunk VCFs and `pileup.vcf.gz` (332,673 records);
- `phase_qual`, the 57,146 het SNPs to phase, and the phased VCF;
- the full-alignment `qual` file and 13 candidate BED chunks (120,785 regions);
- the 13 full-alignment chunk VCFs and `full_alignment.vcf.gz` (120,785 records);
- `merge_chr20.vcf` and `merge_output.vcf.gz` (147,739 records); tabix queries agree.

On 16 threads it took 98.3 s against 341 s for Clair3 (3.5x), on a quiet machine. Clair3 runs
one process per chunk, so each stage waits for its biggest chunk (one pileup chunk holds 124k
of chr20's 333k candidates); lrpipe-snv schedules every 200-candidate batch as its own task
(`crates/lrpipe-snv/src/parallel.rs`), which keeps each batch's arithmetic unchanged. Phasing is
about 8x faster (braid). What remained was libtorch inference, about 370 CPU-s for pileup
and 1,000 CPU-s for full alignment, the same kernels as Clair3: the full-alignment ResNet is about
420 MFLOP per candidate and ran near the CPU's FMA throughput, so with libtorch's arithmetic this
CPU could not go much below about 70-80 s on chr20. lrpipe-snv's own CPU networks take about 1.3x
libtorch's CPU time ([Networks](#networks)).

On the whole HG003 genome (ONT, 24 contigs, 631 pileup chunks, 458 full-alignment jobs, 5.74M
calls; 16 threads), CPU mode took 59 min 7 s against 1 h 22 min 55 s for `run_clair3.py`
(1.4x), and wrote the same 1,629 files apart from `##cmdline`: every pileup chunk VCF, the
phased VCFs, every candidate BED and full-alignment job VCF, every merge VCF, and the three
final VCFs. Clair3 spent 1,308 s on pileup calling, 759 s on phasing and 2,761 s on
full-alignment calling; lrpipe-snv 963 s on pileup calling, then phasing and full-alignment
calling together in 2,515 s. With 631 chunks Clair3's processes stay busy, so the gap was smaller
than on chr20: both were bound by the same libtorch kernels.

### GPU

`lrpipe-snv run --device cuda:0` runs both networks on an NVIDIA GPU (requirements under
[On the GPU](#on-the-gpu)), reproducing Clair3's `--use_gpu --device cuda:0` calling: CPU threads
build each job's tensors and one GPU thread scores them as they arrive, in batches of 1000
(`predictBatchSize * 5`), one VCF per stage, and `SortVcf --use_gpu` merging. The pileup
network takes the CPU path's inputs, with its high-depth scaling: Clair3's GPU path wraps the
counts to `int8` as it saves them and skips the scaling, which lrpipe-snv did too until it made
calls collapse at high depth ([Against Clair3](#against-clair3)).

With libtorch, on HG003 chr20 (RTX 5060 Ti, 16 threads) it took 32.6 s against 4 min 52 s for
Clair3's GPU mode, with the GPU busy in every `nvidia-smi` sample, and wrote the same files as
Clair3's GPU run apart from `##cmdline`: `pileup.vcf.gz`, `phase_qual`, the het SNPs, the phased
VCF, the candidate BEDs, `full_alignment.vcf.gz`, `merge_chr20.vcf` and `merge_output.vcf.gz`
(148,832 records); both networks' GPU probabilities were bit-identical to PyTorch's on the same
GPU. `crates/snv-infer/examples/score_npy.rs` and `tools/score_npy_torch.py` score the same saved
inputs with lrpipe-snv's networks and with PyTorch's, for comparing them.

While a thread waits for the GPU, the CUDA driver by default keeps it spinning on a core.
lrpipe-snv lets it sleep instead (`--gpu-wait block`, the default; lrpipe-snv only), which leaves
the core to phasing and to building tensors: on the whole HG002 HiFi genome (RTX 5060 Ti, 16
threads) 4 min 55 s against 5 min 10 s spinning (2,630 user-seconds against 2,807), on HG003
ONT 7 min 22 s against 7 min 25 s, with the same output. (With libtorch, sleeping had made pileup
calling 8% slower, so it was used only for ONT CRAMs.) `--gpu-wait spin` spins.

On the whole HG002 HiFi genome (RTX 5060 Ti, 16 threads) the run takes 4 min 7 s (3 min 57 s
to 4 min 16 s in three runs), 8.2 times as fast as Clair3's GPU run. In the fastest: pileup
calling 103 s, selecting het SNPs and full-alignment candidates 8 s, then phasing (lrpipe-phase,
on the CPU, 8 contigs at a time) until 207 s, which the full-alignment jobs wait for, ending 24 s
later, and the final merge 4 s. Memory peaks while phasing and full alignment overlap, where it
swings by 2-3 GB between the 5-second returns of freed heap, so the peak varies from run to run:
10.7 to 12.5 GB in four runs. The merges read each contig's records
through the intermediate VCFs' tabix indexes instead of scanning whole files once per contig, and
compress the output's BGZF blocks in parallel (cut where bgzip cuts them, so the files decompress
to the same text). Before the tensor cores, the faster LSTM step, sleeping GPU waits and these
merges it took 5 min 24 s, and 4 min 25 s before lrpipe-phase's faster read selection and DP (chr1's
phasing on 2 threads: 31.9 s, now 16.8 s, with WhatsHap's output). On the whole HG003 ONT
genome the run takes 7 min 1 s (7 min 22 s before lrpipe-phase's speed-up, 3,720 CPU-seconds against
4,281); there full-alignment calling, which reads the BAM again, ends 108 s after phasing.

With libtorch, on the whole HG003 genome (ONT, 24 contigs, 631 pileup chunks, 459 full-alignment
jobs, 5.74M calls; 16 threads, RTX 5060 Ti; the BAM on an NTFS disk that reads at about
400 MB/s), GPU mode took 11 min 22 s to 12 min 22 s (682 s and 742 s in two runs): pileup calling
260 s, bound by libtorch's `float32` LSTMs, then phasing overlapped with full-alignment calling,
bound by the disk. A contig's reads leave the page cache before its full-alignment jobs read them
again, so that stage reads the BAM about twice. Most of the process's memory there was heap that
glibc keeps after it is freed (10-17 GB, next to 2-5 GB in use), so from phasing on it is handed
back every 5 s (`malloc_trim`): on this 31 GB machine the peak fell from 16.7 to 11.7 GB and
swapping by about 70%, the output unchanged and the run time within its run-to-run spread. It
took 13 min 22 s before the per-contig selections ran in parallel, the LSTM weights were packed
for cuDNN and the pileup counts were built read by read; the first whole-genome run, with phasing
after pileup and full alignment after phasing, took 23 min 25 s. All wrote the same 563 files.
Against Clair3's own whole-genome GPU run on HG002 (L40S, 20 threads: 3,328 s for Clair3, 672 s
for lrpipe-snv), the final VCFs and every phasing and candidate file were identical. The pileup
and full-alignment intermediates held the same records, but Clair3 splits them across its GPU
workers (see [VALIDATION.md](VALIDATION.md#gpu)). lrpipe-snv's own GPU networks take the whole
HG002 ONT genome in 602 s against 771 s with libtorch ([Networks](#networks)).

GPU and CPU results differ, in Clair3 as here: GPU kernels change the probabilities in the last
bits (up to 7e-3 on chr20 in Clair3), the GPU path's inputs differ as above, and on chr20 Clair3's
GPU mode makes 148,832 calls against 147,739 on the CPU. So a GPU run is compared with Clair3's
GPU run, and a CPU run with its CPU run; CPU mode agrees the more closely (100% of PASS calls on
HG002 ONT chr20).

With `--gvcf` (`--base_err`, `--gq_bin_size` and `--output_all_contigs_in_gvcf_header` as in
`run_clair3.sh`) and libtorch, `lrpipe-snv run --ctg-name chr20 --threads 4` on a 6 Mb ONT slice
of chr20 (`ont.chr20_30-36m.bam`, 11 of its 13 chunks without reads) wrote the same files as
`run_clair3.py --gvcf`, apart from `##cmdline` (and since then the gVCF fixes in "Clair3 bugs
fixed on purpose", which change the blocks of chunks without reads and `merge_output.gvcf.gz`
wherever they apply): the chunk VCFs, `pileup.vcf.gz`,
`full_alignment.vcf.gz` and `merge_output.vcf.gz` (all with PL), the 13 chunk
`tmp/gvcf_tmp_output/*.tmp.gvcf` (byte-identical lz4 frames), `non_var.gvcf`,
`merge_chr20.gvcf` and `merge_output.gvcf.gz` (234,977 records; bgzip and tabix write the same
`.gz` and `.tbi` from lrpipe-snv's text). So did a run with `--base_err 0.002 --gq_bin_size 3
--output_all_contigs_in_gvcf_header`, and without `--gvcf` the slice matched Clair3's run
without it.

Below that, each piece is checked on its own against Clair3's code:

- pileup tensors and model inputs are bit-identical on all 13 chunks (`tools/parity_contig.sh`);
  with libtorch the probabilities were too, and chunk VCFs were byte-identical to
  `CallVariantsFromCffi --pileup`'s (`tools/parity_vcf.sh`), covering every allele class;
- full-alignment tensors and strings are bit-identical on all 13 jobs and 30 synthetic edge cases
  (`tools/parity_fullalign.sh`); with libtorch the network's probabilities were bit-identical on
  4 jobs (30,785 candidates); the indel-length VCF output rebuilds all 13 job VCFs from Clair3's
  probabilities and matches `output_with` on 1,471 synthetic candidates (both models);
- the CPU networks are held to 1e-4 of the probabilities libtorch computed for synthetic weights
  (`crates/snv-infer/tests/data`), and the GPU networks to the CPU ones;
- `SortVcf`, `MergeVcf`, `SelectQual`, `SelectHetSnp` and `SelectCandidates` match Clair3's
  scripts on synthetic cases: options (RefCalls, haploid, gender male/female with and without
  PAR, qual, ilmn merging), multiple contigs, empty and malformed inputs, ties, signed zeros and
  Python float repr;
- all 42 checkpoints Clair3 ships read identically to PyTorch;
- gVCF: PL on 170 random probability sets (`compute_PL`), the non-variant GQ and PL on 95,000
  depth combinations with five base error rates (2,255 kept as a fixture), non-variant blocks
  on synthetic sites (IUPAC and `N` references, depth swings, per-base mode) and `mergeCalls`
  on 13 synthetic cases; the lz4 frames match lz4 1.10's `lz4 -c` from empty input to
  multi-block files.

Whole chromosomes have been checked on ONT chr20 only; the option matrix above covers the
rest on smaller fixtures.


### Apple GPUs

`lrpipe-snv run --device metal`, in a build with `--features metal` (macOS only), runs both
networks on the Mac's GPU through Metal: the kernels (`crates/snv-infer/src/metal/kernels.metal`)
are compiled by the Metal framework when the models load, called from Rust through `objc2-metal`,
so building needs no Xcode. It schedules work as `--device cuda:0` does (one GPU thread, batches
of 1000), takes the same inputs and runs `float32` throughout, so its calls are the CPU's but where the order of its sums tips a close call; on HG002 chr20 they matched
record for record, QUALs within 0.01. The probabilities are within 6e-7 of the CPU's in the
crate's tests. Timings: [BENCHMARKS.md](../../BENCHMARKS.md#lrpipe-snv-on-an-apple-gpu-metal).
### Approximate mode

`lrpipe-snv run --approximate` (lrpipe-snv only) gives up a little agreement with Clair3 for
speed. Its calls agree with an exact run on more than 99%, counting PASS calls by position,
alleles and genotype (`tools/concordance.py BASE.vcf.gz OTHER.vcf.gz [TRUTH.vcf.gz TRUTH.bed]`),
and their accuracy against GIAB is unchanged. For scale, Clair3's own GPU and CPU runs agree on
98.4% of HG003 chr20.

- On a GPU `--approximate` changes nothing: the networks run in `float32`. Until `c56d769` it ran
  them on the tensor cores in `float16`, and scored the pileup candidates near a quality cut-off
  again in `float32`; with lrpipe-snv's own kernels that saved only 5% of a whole HG003 ONT
  lrpipe run from a BAM (643 s against 674 s, 99.996% of the PASS calls the same) and nothing
  from a CRAM, so it was removed.
- On a CPU with AVX-VNNI (Intel Alder Lake and later desktop CPUs) or AVX-512 VNNI (Xeons from
  Cascade Lake, AMD Zen 4 and later), or on aarch64 with the dot product instructions (Apple M1
  and later, Graviton 2 and later), the full-alignment convolutions after the first run in 8 bits
  (`crates/snv-infer/src/int8.rs`: `u8` activations scaled per candidate, `i8` weights per
  channel with the batch norm folded in, `vpdpbusd` sums, one 512-bit register per 16 channels
  with AVX-512, or `sdot` sums of the activations less 128 with the weights' sums added back),
  and candidates whose two likeliest GT21 or genotype classes are within 0.1 are scored again in
  `float32`. The convolutions read their inputs in place (a 3x3 window's inputs from one row are
  adjacent in the `[row][column][channel]` layout, so no im2col copy) and a batch's images share
  each pass over the weights; at 75% of the CPU's 8-bit peak they are 1.6 times as fast as they
  were, with the same bits. Half of the ONT variant calls go to full alignment instead of 70%
  unless `--var-pct-full` is given.
- On a CPU without these instructions (AMD Zen 2 and 3, Intel Skylake), `--approximate` changes
  nothing, with a warning: the networks run in `float32` and 70% of the ONT variant calls go to
  full alignment, as in an exact run.
- The pileup network's products run in 8 bits too (all but its first input product, on counts,
  and its heads): the recurrent ones and the second layer's input product (as one product per
  step, of `[x_t, h_prev]`, which keeps a step's gates in cache) and `L4`, with activations in
  (-1, 1) offset into `u8`. Candidates whose two likeliest classes are within 0.16 are scored
  again in `float32` by the network, and the pipeline then scores again in `float32` every
  candidate with a QUAL below 7.5 or within 2 of the `--qual` thresholds or the phasing cut-off,
  or within 0.75 of a full-alignment cut-off (`crates/lrpipe-snv/src/rescore.rs`): on HG002 ONT
  chr20, 18% of them, after which the pileup calls and the het SNPs picked for phasing were
  those of an exact run. How far 8 bits move a model's outputs depends on the model: ONT's moved
  probabilities by at most 0.053 and QUALs by 1.4, HiFi's by 0.96 and 17.6, through its
  recurrent products. So a run first scores about 1,000 of its own candidates (the middle
  sixty-fourth of four chunks spread over the run) both ways, and calls pileup in `float32` if
  8 bits moved a probability by more than 0.08 or a QUAL by more than 1.5; HiFi's model fails
  that check, ONT's passes it.

  On HG002 chr20 (16 threads for ONT, 12 for HiFi, i5-12600K):

  | | exact | `--approximate` before | `--approximate` | PASS calls agree with exact |
  |---|---|---|---|---|
  | ONT | 139 s | 83 s | 60 s | 99.21% (pileup calls 100%), as before (99.19%); GIAB SNP F1 0.73921 (before 0.73924, exact 0.74026) |
  | HiFi (pileup in `float32` by the check) | 98 s | | 55 s | 99.96% (pileup calls 100%) |

  On the whole HG003 ONT genome on the CPU (16 threads) `--approximate` took 31 min, against
  72 min exact and Clair3's 83 min (with libtorch, `--approximate` took 42.5 min): pileup
  inference 7,000 CPU-s instead of 15,900, full alignment 12,600 instead of 47,200, and
  2.37 million pileup candidates rescored in 3.6 min. Its pileup calls were exact's and Clair3's
  (100%); its PASS calls agreed with exact's on 99.80% and with Clair3's on 99.78%. Peak memory
  was 12.7 GB, against 11.9 GB exact.

  Under a full load the 8-bit pileup network scores twice as fast as the `float32` one; the
  rescoring takes back about half of that (ONT chr20: 221 CPU-s of inference and 7.5 s of
  rescoring on 16 threads, against 480 CPU-s in `float32`).
- Clair3's Python pipeline (Illumina, `--disable-c-impl`) takes the same 8-bit full-alignment
  convolutions on the CPU; its candidate selection is unchanged, and its pileup network stays
  `float32` (it has no rescoring near the cut-offs). On the Illumina fixture
  (HG003 chr20:30-36 Mb, 4 threads, with libtorch) full-alignment inference took 42.3 CPU-s
  instead of 70.2, and 99.95% of the PASS calls agreed with the exact run, with the same recall
  and precision.

Measured with libtorch, whose exact mode wrote Clair3's output:

| HG003, ONT sup v5.0.0 | exact | `--approximate` | agreement |
|---|---|---|---|
| chr20, CPU (16 threads) | 98.3 s | 66.8 s | 99.50% (against Clair3's CPU run) |
| chr20, GPU | 27.5 s | 18.0 s | 99.99% (against Clair3's GPU run) |
| whole genome, GPU | 682-742 s | 634 s | 99.99% (against the exact run) |
| whole genome, CPU (16 threads) | 59 min 7 s (before later exact-mode speedups) | 42 min 30 s | 99.5% on chr20 |

On the whole genome the GPU gained less: after pileup calling (185 s instead of 260 s), phasing
and full alignment were bound by reading the BAM again (see [GPU](#gpu)).

With lrpipe-snv's own networks, on HG002 ONT chr20 (16 threads, i5-12600K) `--approximate` takes
83 s on the CPU against 139 s exact, and its PASS calls agree with libtorch's `--approximate`
run on 99.994% and with the exact run on 99.185%. Against GIAB v4.2.1 (all of chr20, so
precision counts calls outside the benchmark regions and F1 is low, but comparable), recall is
the same (SNP 0.99787 and 0.99789, indel 0.81358) and precision a little lower: SNP F1 0.74026
exact and 0.73924 approximate, indel 0.62747 and 0.62676.

### Networks

lrpipe-snv runs Clair3's two networks with its own code (`crates/snv-infer`), reading the weights
from Clair3's `pileup.pt` and `full_alignment.pt`. Up to commit `55bad8a` it could also run them
through libtorch, as Clair3 does; libtorch has since been removed, with `--backend` and the
`tch` feature.

On the CPU (`crates/snv-infer/src/cpu*`): the LSTMs as one matrix product for all input steps
plus one per recurrent step, with a vectorized cell update; the convolutions with batch norm
folded in, stride 1 as three products over a zero-padded copy of the input and the strided ones
through im2col; products from the `gemm` crate, in `float32`, with its AVX-512 kernels on CPUs
that have AVX-512 and its AVX2/FMA ones elsewhere. Sums run in another order than
libtorch's (and between those two kernels), so probabilities differ in the sixth decimal (at most 3.3e-6 on 20,000 pileup and
2,000 full-alignment chr20 tensors, no top class changed). `--approximate` adds the 8-bit
convolutions above. The tests hold both networks to 1e-4 of the probabilities libtorch computed
for synthetic weights, saved before it was removed (`crates/snv-infer/tests/data`).

#### On the GPU

With `--device cuda:N` (or Clair3's `--use_gpu`, which means `cuda:0`), both networks run on an
NVIDIA GPU (`crates/snv-infer/src/gpu*`, the `gpu` feature, on by default): CUDA kernels written
for lrpipe-snv (`crates/snv-infer/src/gpu/kernels.cu`), compiled at start-up by NVRTC for the
GPU's compute capability and launched through the CUDA driver. Both libraries are loaded at run
time (`cudarc`), so building needs no CUDA; running on a GPU needs the NVIDIA driver, NVRTC
(`libnvrtc.so.12`, from CUDA or the `nvidia-cuda-nvrtc` wheel) and a GPU of compute capability
6.0 (Pascal) or later. NVRTC is looked for on the library path, then in `$CUDA_HOME`,
`$CUDA_PATH` and `/usr/local/cuda` (`lib64/`, `targets/x86_64-linux/lib/`); `--nvrtc PATH` names
the file instead. From compute capability 8.0 the TF32 products (the convolutions and the LSTMs'
input products) run on the tensor cores (`mma.sync` m16n8k8, double-buffered through shared
memory), as cuDNN runs them there; the other products are register-tiled on the CUDA cores
(convolutions as implicit im2col, long sums split across blocks). Each LSTM step is one kernel
computing the recurrent product and the cell update of both directions, each thread four rows
and two units with vector loads from shared memory.

libtorch's GPU "float32" is not plain `float32`: cuDNN, which libtorch lets use TF32 tensor
cores, multiplies TF32-rounded operands in the LSTMs' input projections and in the convolutions
(the convolution's own weights, the batch norm after). Plain `float32` probabilities differ from
libtorch's GPU ones by up to 2.6e-3, and through the quantile cut-offs that pick het SNPs for
phasing and candidates for full alignment that moved 0.85% of the PASS calls on chr20. So exact
mode rounds the same operands to TF32, with Ampere's `cvt.rna.tf32.f32` from compute capability
8.0 and the same rounding in integer arithmetic before it (products of TF32 values are exact in
`float32`, so only the order of the sums differs), and keeps the recurrence and the dense layers
in `float32`, as libtorch does: probabilities then differ from libtorch's GPU ones by at most
3.2e-4 (pileup) and 1.4e-3 (full alignment) on 70,000 chr20 candidates, no top class changed.
Compiled for compute capability 6.0 and 7.5 and run on an RTX 5060 Ti (12.0), the kernels give
probabilities within 1e-4 of the native build's, which uses the tensor cores. On the whole HG002
HiFi genome the tensor cores changed 0.04% of the PASS calls against the CUDA cores (calls near a
cut-off, and the phasing that follows), and agreed with Clair3's GPU run on 99.707% of them
against 99.710%.

#### Against Clair3

PASS calls (position, alleles, genotype; `tools/concordance.py`) of `merge_output.vcf.gz`
against Clair3 2.0.3's own runs, with the commit after lrpipe-snv joined lrpipe (i5-12600K, RTX
5060 Ti):

| run | Clair3 | lrpipe-snv | PASS calls agree |
|---|---|---|---|
| HG003 ONT, whole genome, CPU, 16 threads | 1 h 23 min, 5.8 GB | 1 h 12 min, 11.9 GB | 99.984% (5.7 million records; pileup 100%, full alignment 99.977%) |
| HG002 Revio HiFi, chr20, CPU, 12 threads | 349 s | 98 s | 100% |
| HG002 Revio HiFi, chr20, GPU | 219 s | 19 s | 99.469% (pileup 99.986%); GIAB SNP F1 0.75509 and 0.75526, indel 0.63180 and 0.63142 |
| HG002 Revio HiFi, whole genome, GPU, 16 threads | 33 min 45 s, 11.4 GB | 4 min 7 s (3:57-4:16), 10.7-12.5 GB | 99.707% (pileup 99.997%, full alignment 98.9%); against GIAB v4.2.1, SNP F1 0.99858 and Clair3's 0.99857, indel 0.96881 and 0.96875 ([VALIDATION.md](VALIDATION.md#hg002-hifi-whole-genome-chr1-22-gpu-against-clair3)) |
| `tools/option_agreement.sh`: every option case of `option_parity.sh` (ONT and HiFi, chr20 and chr22 regions of 6 Mb, 36 cases, one on the GPU) | | | 100% in all but `ont_snp_only`'s `full_alignment.vcf.gz` (99.987%; its `merge_output` 100%) |
| the same for Illumina (9 cases, against Clair3 with its realigner built; see [Illumina](#illumina)) | | | 100% in all |

On the GPU, Clair3 itself scores the pileup tensors as wrapping `int8` and skips the CPU's
high-depth scaling, so its GPU and CPU runs differ: 98.1% on that HiFi chr20. The GPU rows above
were measured when lrpipe-snv's GPU calling did as Clair3's GPU calling. At sites deep enough
that a count passes 127 it wraps around (-139 became 117 at chr20:29,505,541, depth 248) and the
network's output there is close to arbitrary (that site: 0/1 at 86% from libtorch's GPU, 96% from
lrpipe-snv's; 1/1 at 93% from both on the CPU; 0/0 with the CPU path's scaling). On a human
genome at 30-60x that touched few sites, but on a deep sample nearly every count wraps: K-12 HiFi
at 295x called SNPs at F1 42.5% and indels at 1.4% on CUDA, against 96.4% and 97.6% on the CPU
(#45). GPU calling now takes the CPU path's
inputs, as Metal does, and rounds the same operands to TF32 as libtorch's cuDNN, so its calls
differ from the CPU's only where TF32 tips a close call, and from Clair3's GPU run at the sites
its wrapping changes. On that HiFi chr20 the CPU path's inputs brought GPU pileup calls to
Clair3's CPU calls (99.988%); in plain `float32` instead of TF32 the GPU agreed with Clair3's GPU
run less (99.04%).

On the CPU, Rust's `float32` networks are slower than libtorch's (MKL) were: the whole HG003
ONT genome took 72 min against 59 min, and HG002 ONT chr20 131 s against 111 s. With
`--approximate` it took 31 min, against libtorch's 42.5 min ([Approximate mode](#approximate-mode)).

#### Against libtorch

HG002 ONT, 16 threads, i5-12600K and RTX 5060 Ti, against the former libtorch backend, whose exact
mode wrote Clair3's output:

| HG002 ONT | libtorch | Rust | PASS calls agree |
|---|---|---|---|
| chr20, CPU, exact | 111 s (1,312 CPU-s) | 131 s (1,696 CPU-s) | 100% (144,477); the same GIAB recall on chr20 |
| chr20, CPU, `--approximate` | 79 s | 83 s | 99.994% (against libtorch's `--approximate`) |
| chr20, GPU, exact | 44.2 s (pileup calling 12.8 s) | 37.6 s (7.8 s) | 99.945%; GIAB SNP and indel F1 the same to five decimals |
| whole genome, GPU, exact | 771 s (281 s) | 602 s (188 s) | 99.827% |

On the whole genome the Rust run's peak RSS was 14.6 GB against libtorch's 12.4 GB. The extra
comes from how phasing and full alignment overlapped in the two runs, not from the networks: on
chr1-chr3 alone the Rust run peaked at 9.0 GB against 9.45 GB.

The PASS calls that differ on the whole genome come mostly from phasing: one het SNP landing on
the other side of the phasing cut-off re-phases blocks and gives full alignment other haplotags.
In an earlier whole-genome comparison (99.83% agreement), only 36 of 10.46 million pileup records
differed in position, alleles, filter or genotype; per contig the calls agreed on 99.98% or more
on 9 contigs and on 99.33-99.96% on the others: on chr3 the het SNPs selected for phasing were
the same and so were all phased genotypes, while on chr7 one SNP of 118,299 fell on the other
side of the phasing cut-off (31,121 phased genotypes differed). chr20's GIAB F1 within those
whole-genome runs: SNP 0.99869 (libtorch) and 0.99863 (Rust), indel 0.85473 and 0.85440.

Scoring alone (batches of 1,000 chr20 candidates, RTX 5060 Ti), libtorch took 22 µs per pileup
candidate and 61 µs per full-alignment one in `float32`; Rust takes 9.2 and 50 µs.

### Illumina

For `--platform ilmn`, `run_clair3.py` runs Clair3's Python pipeline (`scripts/clair3_pipeline.py`)
instead of the C implementation, and `lrpipe-snv run --platform ilmn` follows it
(`crates/lrpipe-snv/src/python.rs`; ONT and HiFi take the same path with `--disable_c_impl`,
without realignment: `samtools mpileup` reads the haplotagged BAM within the candidate BED, and
the BED rows' centres are the candidates; on the ONT slice of chr20 with 8 threads it takes
31 s against 143 s for `run_clair3.py --disable_c_impl`, and on HiFi 24 s against 72 s):

| stage | Clair3 | crate |
|---|---|---|
| 1. pileup tensors | `samtools mpileup` text parsed by `CreateTensorPileup.py` | `snv-mpileup` |
| 4. haplotagging | `whatshap haplotag --regions CTG` into `tmp/phase_output/phase_bam/CTG.bam`, `samtools index` | `snv-phase` (lrpipe-phase) |
| 6. realignment | `RealignReads.py`: de Bruijn graph haplotypes (`debruijn_graph.cpp`) and realignment (`realigner.cpp`, SSW) | `snv-realign` |
| 6. full-alignment tensors | `samtools mpileup -` over the realigned reads, parsed by `CreateTensorFullAlignment.py` | `snv-realign` |

`samtools mpileup` is not called: its columns come from htslib's pileup engine with samtools
1.24's read filters (`mplp_func`: flags, MAPQ, anomalous pairs) and `pileup_seq`'s text, and
the Python that parses the text is ported as it is. Calling, phasing, candidate selection and
merging (`MergeVcf` in its Illumina, region-based mode) reuse the long-read crates.

Clair3 2.0.3 itself cannot run this path's full-alignment stage: the bioconda package ships no
compiled realigner (`RealignReads.py` fails to load `preprocess/realign/realigner`), and when it
is built from the 2.0.3 sources, `free_memory` is `static` (commit 583cf9a) so the library does
not export it and `RealignReads.py` fails at its first realigned window. Either way Clair3 writes
an empty `full_alignment.vcf.gz` and copies `pileup.vcf.gz` to `merge_output.vcf.gz`. lrpipe-snv
reproduces the pipeline as it runs with a working realigner: the oracle is Clair3 2.0.3 with
`realigner` and `debruijn_graph` built as its Dockerfile builds them and `free_memory` exported.

With libtorch, on the Illumina fixture (chr20:30-36 Mb, 4 threads), `lrpipe-snv run --platform
ilmn` wrote the same files as that Clair3 run, apart from `##cmdline` and the haplotagged BAM's
`@PG` line: the 2 pileup chunk VCFs and `pileup.vcf.gz` (38,529 records), the phased VCF,
`phase_bam/chr20.bam` (same records), the candidate BEDs, the 2 full-alignment job VCFs and
`full_alignment.vcf.gz` (31,673), `merge_chr20.vcf` and `merge_output.vcf.gz` (13,063). It took
89 s against about 8 min; realignment runs the windows of each 5 kb chunk in parallel and applies
their results in window order, as Clair3 applies them one after the other.
With `--bed_fn` (three regions), `--vcf_fn` (4,816 known sites) and `--no_phasing_for_fa
--print_ref_calls --enable_variant_calling_at_sequence_head_and_tail` the outputs were identical
too (`tools/option_parity.sh`, 21, 28 and 16 files). So were the gVCF runs: `--gvcf`, with
`--bed_fn`, with `--vcf_fn`, and with `--base_err 0.002 --gq_bin_size 3
--output_all_contigs_in_gvcf_header --print_ref_calls
--enable_variant_calling_at_sequence_head_and_tail` (34, 37, 44 and 34 files, among them every
chunk's `tmp/gvcf_tmp_output/*.tmp.gvcf`, `non_var.gvcf` and `merge_output.gvcf.gz`). With
`--gvcf` Clair3 runs `samtools mpileup -a`, whose columns without reads lrpipe-snv writes as
samtools does (checked on all of chr20, with and without `-l`); inside the chunk they become
non-variant sites, and in the 33 flanking positions their `*` counts as a deletion in the
tensors, so a candidate near a chunk edge can get a tensor it would not get without `--gvcf`.
Below that, the mpileup bases column matches `samtools mpileup` 1.24 on all 4.4M columns of a
chunk, and pileup tensors, realigned SAM lines (276,053 reads of one job) and full-alignment
tensors match Clair3's Python line for line (`tools/parity_ilmn_jobs.sh`); the SSW port matches
the C library on 100,000 random alignments, apart from two where nothing aligns and the C code
compares the byte before its reference, and the graph port the C++ on 20,000 random windows
(`tools/parity_realign.sh`).

The SSW kernels run on SSE2 with the C library's instructions in the same order, and a test checks
them against the lane-by-lane emulation on 16,000 random alignments. On HG003 chr20 (12 threads)
that took Illumina full alignment from 423 s to 242 s, and its CPU time from 2,055 to 1,385 CPU-s
(`perf`). Smith-Waterman went from 51% of that time to 27%; libtorch inference was then the
largest part (38%), and the 16-bit kernel accounts for most of SSW.

### Clair3 bugs fixed on purpose

lrpipe fixes real bugs of the tools it reproduces rather than copying them; output differs
from Clair3's only where one applies (review findings LRP-002, 003, 005, 012 and 020, and the
gVCF's, found reviewing `lrpipe run --gvcf`):

- **Leftovers of an earlier run.** A run starts by removing its `tmp/` and the final files an
  earlier run wrote in the output directory (`pileup.vcf.gz`, `merge_output.vcf.gz`,
  `phased_merge_output.vcf.gz`, their indexes and the gVCF and phased BAM). Clair3 lists its
  chunk VCFs by name, so an interrupted run's chunks joined the next run's calls, and its
  early exits kept a merged VCF from an earlier run while copying the new pileup index beside it.
  It removes them only after checking that every input can be read (the BAM or CRAM and its
  index, the FASTA and its `.fai`, `--bed_fn`/`--vcf_fn` and the model's checkpoints), so a
  mistyped path fails with the earlier results intact, and only what lrpipe-snv marked as its
  own with a `.lrpipe-snv` file (one in `tmp/`, one beside the results): any other `tmp/` (a
  project's own with `--output .`, even one made after an lrpipe-snv run there) or file with
  one of those names (an earlier Clair3 run, or a run of lrpipe-snv before 0.1.6) makes the
  run stop with them untouched. A rerun that fails after
  the checks still leaves none of the earlier run's results: keep them elsewhere (or use
  another output directory) if they matter.
- **PAR coordinates.** With `--gender male`, a PAR BED row `[start, end)` keeps VCF positions
  `start + 1` to `end` diploid. Clair3 compares the 1-based position with the 0-based BED, one
  base off at each end.
- **Phasing beside haploid calls.** With `--gender male` the final calls hold haploid genotypes
  on chrX and chrY outside the PARs and diploid ones inside them. lrpipe-phase leaves the
  haploid ones unphased; WhatsHap, so Clair3's `--enable_phasing`, stops with "Inconsistent
  ploidy (2 and 1)".
- **A BAI without read counts.** The counts in a BAI (its pseudo-bin) are optional, contig by
  contig. Clair3 takes a contig without them for one without reads and calls nothing there;
  lrpipe-snv calls it, as for a CRAM, unless the index has no bins for it, which means no reads.
  A count elsewhere in the index says nothing of this contig's.
- **The phased VCF of a run without calls.** With `--enable_phasing`, a run that ends early
  writes `phased_merge_output.vcf.gz` too: without records when there is nothing to call, and
  the pileup calls phased when no call goes to full alignment. Clair3 writes only the merged
  VCF there, so `lrpipe run` found the phased one it was promised missing.
- **`--pileup_only` per contig.** Each record is written as post-processing leaves it
  (haploid genotypes, `LowQual`), as Clair3 does on the GPU; on the CPU Clair3 tests the
  post-processed row but writes the original one, so only its RefCall filter applied.
  `--pileup_only` with `--gvcf` is refused: Clair3 stops before the gVCF, and lrpipe-snv
  did too, exiting 0 without it.
- **Targeted gVCFs.** With `--bed-fn`, reference blocks stay inside the original BED intervals.
  Clair3 can write reference blocks across gaps between targets.
- **gVCF: every call, and blocks that neither overlap nor go missing.** The gVCF holds every
  call of `merge_output.vcf.gz` and, of each non-variant block, the parts no call covers, with a
  part's new start given its reference base. Clair3's `gvcfGenerator.mergeCalls`
  (`preprocess/utils.py`) instead writes a second copy of the current block in place of a call
  that ends before it (a SNV inside an earlier deletion: the call is lost, and on the ONT
  fixture one of 467 is); writes the left part of a block an earlier call already cut with the
  block's original start, and the block holding a contig's last call whole after it, so blocks
  overlap calls and each other; writes the blocks after a contig's last call as read, over a
  deletion reaching into them; and drops the first block of a contig without calls, which joint
  genotyping takes for a no-call. A chunk without reads is one `N` block of its own sites
  (`start` to `end - 1`, and not past the contig's end): Clair3's `write_empty_pileup` ends it at
  the next chunk's first position, so on a whole genome every run of such chunks (centromeres,
  contigs without reads) overlaps by one base at each chunk boundary, and the last chunk's runs
  past the contig. `ci/parity.sh` checks the fixture's gVCF against Clair3's but for these
  (`ci/fixture/expected/merge_output.gvcf.fixes`).

### Behaviour reproduced on purpose

- Insertion alleles are written in khash bucket order: `khash.rs` emulates klib khash 0.2.8
  exactly (tested against Clair3's `khash.h` with `tools/khash_oracle.c`).
- An `N` (or other IUPAC) read base indexes feature -1, which adds 1 to the previous position's
  last feature.
- Python formats the region as `max(0, start - 1)` and `hts_parse_reg` subtracts one again.
- Outcome probabilities are numpy `float32` products, and QUAL is computed in `float32` up to
  the logarithm (NumPy 2 promotion rules).
- Python tie-breaking: `sorted(..., reverse=True)` keeps ties in order while
  `sorted(...)[::-1]` reverses them; Clair3 uses both.
- The decoding loop's `continue` re-tests `reference_base`/`alternate_base`, which persist
  across iterations, so some "give up" branches end the loop with the alleles already chosen.
- With indel-length heads, `hetero_DelDel_tuples_from` drops equal length pairs except
  (16, 16), and a proposed length of 16 ("16 or more") proposes nothing.
- Full-alignment read selection for deep candidates uses glibc's `rand()` (TYPE_3, seed 1),
  emulated per job because Clair3 runs every job in a fresh process.
- Illumina: realigned SAM lines carry PNEXT in their TLEN column (`self.TLEN = PNEXT`), which
  `samtools mpileup`'s overlap detection reads; read names become 5-character base-87 counters,
  so mates are never paired; soft clips count like insertions in the realignment pileup;
  `random.seed(0); random.sample(...)` picks the reads of deep candidates (MT19937 and
  `_randbelow` emulated); haplotype alignments are ordered by libstdc++'s unstable `std::sort`;
  SSW's lazy-F loop does not update `E`, so its 8- and 16-bit striped kernels are reproduced
  with the same SSE2 instructions (emulated lane by lane on other targets) rather than replaced
  by a textbook dynamic programme.
- `normalize_af` in the full-alignment C code expands without parentheses, and a missing int
  key reads as -1 through `size_t`; out-of-bounds table reads for `=` bases use libclair3's
  actual `.rodata` bytes.
- With `--enable_long_indel`, a long insertion's AD also counts insertions within 10% of its
  length, compared by the alt-info key's length (which includes the anchor base); the
  deletion version never counts any, because Clair3 calls it without `is_del=True`.
- `--chunk_num -1` (no chunking) leaves Clair3 one job, for a contig named `None` (`CHUNK_LIST`
  is `None 0 0`), so Clair3 2.0.3 calls nothing and writes header-only VCFs, on both pipelines.
  With `--bed_fn` that job does call the BED's regions, into `tmp/pileup_output/pileup_None_0.vcf`
  (and non-variant blocks with `--gvcf`), which the sort step then ignores. lrpipe-snv writes those
  files too, calling the merged regions in the order Python iterates the interval tree's `set`
  (CPython 3.11's table order, emulated); on the GPU, where Clair3 builds no tensors for that
  job, it skips them.
- `--no_phasing_for_fa` links `tmp/phase_output/phase_bam/<contig>.bam` to the input, and Clair3
  links only a `.bai` next to it; lrpipe-snv also links a `.csi` or `.crai`, so a CRAM input works.
- Dwell time (`--enable_dwell_time`, the `_with_mv` models): the full-alignment tensors get a
  ninth channel, each base's count of move-table blocks (`mv` tag) from its first `1` to the
  next, read in basecall order and reversed for reverse-strand reads; blocks before the first
  base are dropped. An insertion's blocks add to the base before it, also when a second
  insertion there replaced the first's bases. The count is written as `int8`, so stalls over
  127 blocks wrap. Reads without a usable tag get 0. As in `run_clair3.py`, move-table models
  need ONT and the option, and the option needs an `mv` tag with values in the first 50
  records. The test BAMs carry made-up move tables (no BAM with real ones here); a second
  table changes 2,528 of 2,542 full-alignment records, and both were identical to Clair3's
  output with libtorch.
- `--vcf_fn` chunks by known sites (`len // chunk_num`, rounded up) and keeps the extended BED
  filter of `split_extend_vcf` too; missed sites are added back with only their first REF base.
- gVCF: the C pileup's total depth adds only the alternate counts that raise the running
  maximum; chunks starting at 2-34 read their depths one position late; Clair3's
  `LOG_10 = 2.3025` and the out-of-bounds maximum in its `log10sumexp` (which behaves as 0, and
  matters because 2.3025 is not ln 10); `round(x, 6)` after each step; only a literal `N`
  reference splits blocks (an IUPAC base joins its neighbours). (Its merge of calls and blocks
  is fixed: see above.)

## Building

lrpipe-snv is part of the lrpipe workspace: `lrpipe snv` runs it (with `run` implied), and it also
builds as its own binary. Needs Rust 1.91 or later and a C compiler, for psm's build script (through
ndarray-npy); BGZF goes through pasta's zflate, in Rust. The lz4 frames of the gVCF intermediates
are written and read with `lz4_flex`, and lrpipe-phase is built without mimalloc, both in Rust.
Reading and writing alignments, VCFs and indexes goes through
[pasta](https://github.com/DylanCosto/pasta), in Rust; there is no htslib, and no libtorch or CUDA
is needed to build. From lrpipe's directory:

```
cargo build --release -p lrpipe                             # lrpipe, with `lrpipe snv`
cargo build --release -p lrpipe-snv                         # target/release/lrpipe-snv alone
cargo build --release -p lrpipe-snv --no-default-features   # CPU only
```

The binary links only libc, libm and libgcc_s. With the `gpu` feature (the default) it loads the
CUDA driver and NVRTC at run time for `--device cuda:N` ([On the GPU](#on-the-gpu)); without it,
`--device cuda:N` is an error.

pasta is a git dependency pinned, in lrpipe's `Cargo.toml`, to a commit of `DylanCosto/pasta`;
the phasing is lrpipe's own [lrpipe-phase](../lrpipe-phase) (WhatsHap in Rust). lrpipe-snv
began as its own repository; commits quoted in these documents are from its development
history, which this repository does not include.

On this machine `env.sh` puts the original's env (`TORCH_ENV`: Clair3 2.0.3 with PyTorch 2.13.0,
`micromamba create -p ... -c conda-forge -c bioconda clair3=2.0.3 pytorch=2.13.0`) on `PATH`, for
Clair3's models and for the tools that run Clair3 for parity (`tools/*.py`, the oracles); building
needs nothing from it.

### Releases

lrpipe's release bundle (`tools/release_build.sh` in lrpipe) carries `lrpipe` and
`lrpipe-snv` with the other tools' binaries, `LICENSE-CLAIR3`, and lrpipe-phase's and pasta's licenses;
lrpipe's `Dockerfile` puts it in a container image. The bundle needs glibc 2.28 or later and, on
a GPU, the NVIDIA driver and NVRTC at run time. Built in `manylinux_2_28`, its `lrpipe snv` and
`lrpipe-snv` pass `ci/parity.sh` in a clean environment, as does a run of the fixture on the GPU
([BENCHMARKING.md](BENCHMARKING.md#6-the-release-bundle)). Neither the bundle nor the image
includes Clair3's models. (lrpipe-snv's own 0.1.0 bundle was CPU-only and carried libtorch, MKL and OpenMP.)

## Parity testing

[BENCHMARKING.md](BENCHMARKING.md) walks through re-running VALIDATION.md's parity, accuracy and
speed checks on the current code, with `tools/bench_run.sh` timing and sampling each run.

- `tools/parity_vcf.sh BAM REF CONTIG CHUNKS OUTDIR [JOBS]`: runs Clair3's
  `CallVariantsFromCffi --pileup` and `lrpipe-snv pileup --call-fn` on every chunk and compares
  the VCFs byte for byte. They were identical with libtorch; with lrpipe-snv's own networks,
  records differ where the probabilities' last digits reach the output.
- `tools/parity_contig.sh BAM REF CONTIG CHUNKS OUTDIR [JOBS]`: runs `tools/oracle_pileup.py`
  (Clair3's own Python functions) and `lrpipe-snv pileup` on every chunk and compares tensors,
  inputs and probabilities bit for bit (`tools/compare_pileup.py`). The tensors and inputs still
  match; the probabilities only did with libtorch.
- `tools/oracle_fa_output.py`: runs one of Clair3's full-alignment jobs (from its
  `parallel_6_call_var_bam_full_alignment.log`) with `batch_output` wrapped to save each
  candidate's position, alt info and probabilities; `SNV_FA_ORACLE=DIR cargo test -p
  snv-call` then rebuilds every job's VCF from them and compares it with Clair3's.
- `tools/parity_fullalign.sh PARALLEL_LOG OUTDIR [PHASED_VCF] [JOBS]`: re-runs every
  full-alignment job of a Clair3 run (from `log/parallel_6_call_var_bam_full_alignment.log`)
  under `tools/oracle_fullalign.py`, which stops Clair3's own process once the tensors are
  built, and compares `snv-fullalign`'s tensors and strings bit for bit
  (`crates/snv-fullalign/tests/parity.rs`, opt-in via `SNV_FA_CASES`).
- `tools/option_agreement.sh OUTDIR [CASE...]`: reruns `option_parity.sh`'s cases (each with the
  options Clair3 recorded in its VCF header) and compares the VCFs by PASS calls, which is how
  they are compared since libtorch was removed ([Against Clair3](#against-clair3)).
- `ci/parity.sh LRPIPE_SNV MODEL_DIR`: the parity part of lrpipe's `tools/check.sh`
  (`ci/README.md`): `lrpipe snv` on a 40 kb ONT fixture, from BAM and from CRAM, against the VCFs Clair3 wrote
  for it. The CRAM must give exactly the BAM's VCFs, and each VCF must be the same as Clair3's
  or its PASS calls must agree with Clair3's on at least 99.9% (agreement and Jaccard).
- `cargo test`: unit tests, including khash iteration order against C fixtures in
  `crates/snv-pileup/tests/data` and the CPU networks against libtorch's probabilities in
  `crates/snv-infer/tests/data`; the GPU tests run when the CUDA driver, a GPU and NVRTC are
  found. With `SNV_MODELS` and `SNV_DIGESTS` (from `tools/checkpoint_digest.py`) it also checks
  every checkpoint against PyTorch.
- `tools/oracle_full_alignment.py`: runs one of Clair3's full-alignment jobs, with the arguments
  from its GNU parallel log, and records the model inputs, batches and probabilities. (The test
  that scored them with libtorch, bit for bit, was removed with libtorch.)
- `crates/snv-infer/examples/score_npy.rs` and `tools/score_npy_torch.py`: score saved model
  inputs (`.npy`) with lrpipe-snv's networks and with Clair3's PyTorch ones, on the CPU or a GPU,
  and write the probabilities for comparison.
- `tools/parity_ilmn_jobs.sh CLAIR3_BIN RUN_DIR BAM REF CONTIG CHUNK_ID CHUNK_NUM BED OUTDIR`:
  runs Clair3's `CreateTensorPileup`, `RealignReads` and `CreateTensorFullAlignment` on one
  Illumina chunk and one full-alignment job and compares them line by line with
  `examples/mpileup_dump.rs` and `examples/realign_dump.rs` (also `mpileup_dump text` against
  `samtools mpileup`).
- `tools/parity_realign.sh REALIGN_DIR OUTDIR [CASES]`: links `tools/ssw_harness.cpp` and
  `tools/dbg_harness.cpp` against Clair3's `preprocess/realign` sources and checks the SSW and
  de Bruijn graph ports on random inputs (`tools/make_realign_cases.py`; the graph needs Boost
  headers, `BOOST_INCLUDE`).
- `tools/option_parity.sh` takes `CLAIR3_RUN` to run another `run_clair3.py` (for Illumina, a
  Clair3 tree with a built realigner).
- `tools/concordance.py BASE.vcf.gz OTHER.vcf.gz [TRUTH.vcf.gz TRUTH.bed]`: the share of
  `BASE`'s PASS calls (position, alleles, genotype) that `OTHER` makes, and both runs' recall,
  precision and F1 against a truth set: lrpipe-snv against Clair3, `--approximate` against exact,
  and GPU runs.
- `tools/happy_summary.py happy.extended.csv`: hap.py's PASS SNP and INDEL metrics overall and
  in a few GIAB stratifications (difficult regions, homopolymers, low mappability, MHC), as in
  [VALIDATION.md](VALIDATION.md).

## Layout

| crate | what |
|---|---|
| `crates/snv-hts` | indexed BAM/CRAM reading, an htslib-compatible pileup and FASTA access, on pasta |
| `crates/snv-pileup` | pileup features, candidates and tensors (`hts.rs` holds all `unsafe`) |
| `crates/snv-infer` | model traits, checkpoint reader, the networks on the CPU (`cpu*`) and CUDA GPUs (`gpu*`, `gpu/kernels.cu`), and the 8-bit convolutions of `--approximate` (`int8.rs`) |
| `crates/snv-call` | genotype decoding and VCF output |
| `crates/snv-merge` | `SortVcf` and `MergeVcf`: merging, bgzip, tabix (`hts.rs` holds its `unsafe`) |
| `crates/snv-select` | `SelectQual`, `SelectHetSnp` and `SelectCandidates` |
| `crates/snv-phase` | in-process phasing through lrpipe-phase |
| `crates/snv-fullalign` | full-alignment tensors and haplotagging |
| `crates/snv-mpileup` | `samtools mpileup` columns and Illumina pileup tensors (`hts.rs` holds its `unsafe`) |
| `crates/snv-realign` | Illumina realignment (de Bruijn graph, SSW) and full-alignment tensors |
| `crates/lrpipe-snv` | the `lrpipe-snv` binary |

## Acknowledgements

lrpipe-snv is a Rust port of [Clair3](https://github.com/HKU-BAL/Clair3), the variant caller from
the University of Hong Kong's Department of Computer Science, and runs Clair3's trained models.
lrpipe-snv's own license, BSD 3-Clause like Clair3's, is in [`LICENSE`](LICENSE), Clair3's in
[`LICENSE-CLAIR3`](LICENSE-CLAIR3), and that of the Striped Smith-Waterman library its realignment
ports (MIT) in [`LICENSE-SSW`](LICENSE-SSW). lrpipe-snv is not
affiliated with or endorsed by the Clair3 authors.

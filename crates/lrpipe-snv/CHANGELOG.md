# Changelog

lrpipe-snv is released with lrpipe, under its version, since it joined lrpipe.

## 0.3.0

- With `--bed-fn` and `--gvcf`, reference blocks stay inside the BED's intervals, on the pileup
  and mpileup paths, and in chunks without reads. Clair3 writes blocks across the gaps between
  targets, from the pileup span of the whole BED (#53).

## 0.2.2

- `--pileup-only` with `--gvcf` is refused before the run starts. Clair3 stops after the
  pileup without a gVCF, and so did lrpipe-snv wherever the pileup had calls, exiting 0
  without `merge_output.gvcf.gz` (#54).
- CUDA calling takes the CPU path's pileup inputs, with its high-depth scaling, as Metal does.
  It reproduced Clair3's GPU path, which wraps the counts to `int8` and skips the scaling, and
  calls collapsed at high depth: K-12 HiFi at 295x gave SNP F1 42.5% and indel F1 1.4% against
  the CPU's 96.4% and 97.6% (#45).

## 0.2.1

- Apple GPUs: `--device metal`, with the `metal` feature (macOS), runs both networks as Metal
  kernels in `float32`, on the CPU path's inputs, so its calls are the CPU's but where the order of its sums tips a close call.

## 0.2.0

- A read's tags are read past a `d` (double) tag, which samtools writes to BAM: lrpipe-snv
  stopped at one, and so missed a move table (`mv`) after it, which dwell-time calling reads.
- BGZF goes through pasta's zflate, in Rust, in place of libdeflate (pasta `8e2d1e9`).
- Run by `lrpipe run`, lrpipe-snv reports how far it has got (its pileup chunks and
  full-alignment jobs) to the file `LRPIPE_PROGRESS` names, for the run's display; run on its
  own, it writes nothing.

## 0.1.9

- No changes to lrpipe-snv; released with lrpipe 0.1.9.

## 0.1.8

- No changes to lrpipe-snv; released with lrpipe 0.1.8.

## 0.1.7

- The CRAM queries of a run share each contig's reference sequence, which each pileup chunk
  and full-alignment job read again for itself. pasta `8041483` checked a CRAM slice without a
  reference checksum against `@SQ` M5, which hashed the whole contig again for each chunk too
  (0.2.0 hashes only a contig whose length differs from `LN`): with that pasta alone,
  chr17-chr22 of the HG002 ONT CRAM took 5% more CPU time (1,050 CPU-s, 0.1.6: 1,000); shared,
  1,007 CPU-s, with the same calls and a lower peak memory (7.1 GB, 0.1.6: 8.0 GB), as concurrent
  chunks no longer each hold a copy of the contig.

## 0.1.6

- **Clair3 bugs fixed: the gVCF** (`--gvcf`; README, "Clair3 bugs fixed on purpose"). Clair3's
  `mergeCalls` lost a call inside an earlier deletion (replacing it with a second copy of a
  block), left blocks overlapping calls and each other, dropped a contig's first block when it
  had no calls, and ran a readless chunk's `N` block one base into the next chunk (and past the
  contig's end). Every call is now written and blocks neither overlap nor go missing;
  `ci/parity.sh` checks the fixture's gVCF against Clair3's but for these differences.
- **No record joined to the `#CHROM` line** (GPU, with phasing): a contig without full-alignment
  calls took its header from a header-only VCF without the last newline, so the next contig's
  first record joined `#CHROM` in `merge_output.vcf.gz`. Header lines now end with a newline.
- **Earlier results are kept:** the BAM or CRAM and its index, the FASTA and its `.fai`,
  `--bed_fn`/`--vcf_fn` and the model's checkpoints are checked before an earlier run's results
  are removed, and `tmp/` and the results are removed only where lrpipe-snv marked them as its
  own (`.lrpipe-snv`), so `--output .` no longer deletes a project's `tmp/`. A standalone output
  directory of 0.1.5 or earlier is refused once.
- An empty or cut-off checkpoint is an error, not a crash; `-m DIR` without the networks names
  `DIR/pileup.pt`.
- `--approximate` on a CPU without 8-bit dot-product instructions changes nothing, with a
  warning (it sent fewer ONT calls to full alignment, for no speed-up).
- LongPhase ends with lrpipe-snv (`PR_SET_PDEATHSIG`).

## 0.1.5

- **Clair3 bugs fixed** (second review, LRP-005 and 020): a contig without a read count in a BAI
  that counts others is called, unless the index has no bins for it; with `--enable_phasing`,
  every early exit writes `phased_merge_output.vcf.gz` and its index, without records when there
  is nothing to call and as the pileup calls phased when none goes to full alignment.

## 0.1.4

- **Clair3 bugs fixed** (review findings LRP-002, 003, 005, 012; README, "Clair3 bugs fixed on
  purpose"): a run removes an earlier run's `tmp/` and final files first; a PAR BED row covers
  the VCF positions it names; final phasing leaves haploid calls unphased instead of stopping; a
  BAI without read counts no longer skips every contig; `--pileup_only` writes post-processed
  records on the CPU too.

- **AVX-512 on the CPU.** `--approximate`'s 8-bit products run on CPUs with AVX-512 VNNI
  (Xeons from Cascade Lake, AMD Zen 4 and later), with each block's 16 channels in one 512-bit
  register; they used to need AVX-VNNI, which most server CPUs lack, so there `--approximate`
  ran in `float32`. The sums are the same integers either way. The `float32` products use
  `gemm`'s AVX-512 kernels on CPUs with AVX-512, which sum in another order than its AVX2 ones:
  probabilities can differ in the sixth decimal between such CPUs. Not yet measured on an
  AVX-512 machine.
- **Apple silicon and Linux aarch64**: `--approximate`'s 8-bit products run with `sdot` on CPUs
  with the dot product instructions, where they fell back to `float32`: on an M4 Pro, HG002 ONT
  chr20 in 56 s instead of 126 s, with the same GIAB F1. The Smith-Waterman realignment has a
  NEON version. On macOS freed heap goes back with `malloc_zone_pressure_relief`.
- Phasing and full-alignment selection read the pileup VCF one contig at a time from disk
  instead of decompressing it whole into memory (HG002 ONT at 10x: 265 million records,
  18.6 GB). The calls are the same: on the whole HG003 ONT genome on the GPU and on chr20 on the
  CPU, exact and `--approximate`, every pileup, full-alignment and final record is as before.

## 0.1.3

- On the GPU with final phasing (`--enable_phasing`, which `lrpipe run` now uses), each contig
  is merged and phased as soon as its full-alignment jobs are scored, beside the other contigs'
  calling, instead of all of them at the end: the same files, 20-30 s sooner on whole genomes,
  with a higher peak memory on ONT (13.4-14.3 GiB against 9.8-11.7).

## 0.1.2

- **Phasing about twice as fast** (lrpipe-phase, formerly braid; the same output as WhatsHap): read selection
  indexes instead of hashing, the DP computes a column's partitionings in binary order,
  vectorized, with Gray code order's ties. With it the whole HG002 HiFi genome on the GPU takes
  3 min 57 s.
- **GPU calling about 18% faster**: the whole HG002 HiFi genome in 4 min 25 s against 5 min 24 s
  (RTX 5060 Ti, 16 threads; Clair3's GPU run takes 33 min 45 s), with 99.707% of the PASS calls
  Clair3's GPU run's (99.710% before):
  - From compute capability 8.0 the TF32 products (convolutions and the LSTMs' input products)
    run on the tensor cores, as cuDNN runs them; probabilities stay within 1e-4 of the CUDA
    cores'.
  - The LSTM step kernel reads its operands with vector loads, each thread computing four rows
    and two units, with the same sums.
  - Threads waiting for the GPU sleep by default (`--gpu-wait block`), leaving their cores to
    phasing; `--gpu-wait spin` keeps the old behaviour.
  - The merges read each contig's records through the intermediate VCFs' tabix indexes rather
    than scanning the whole files per contig, and compress BGZF blocks in parallel, cut where
    bgzip cuts them.

## 0.1.1

- **`--approximate` on the CPU about 1.4 times as fast as before** (HG002 ONT chr20: 60 s
  against 83 s; exact takes 139 s; the whole HG003 ONT genome in 31 min against 72 min exact and
  Clair3's 83 min, with 99.78% of the PASS calls Clair3's):
  - The 8-bit full-alignment convolutions read their inputs in place instead of through an
    im2col copy, run a batch's images through each pass over the weights, and quantize and
    finish in AVX2: 1.6 times as fast, with the same results bit for bit.
  - The pileup network's products run in 8 bits (all but its first input product and its
    heads), with close calls scored again in `float32` and, in the pipeline, every candidate
    with a QUAL below 7.5 or near the `--qual` thresholds, the phasing cut-off or a
    full-alignment cut-off: the pileup calls and the het SNPs picked for phasing are an exact
    run's on HG002 ONT chr20. A run first checks the 8-bit network on its own candidates and
    calls pileup in `float32` if it moves them too far, as it does with HiFi's model.
- The intermediate gVCFs' lz4 frames are written and read with `lz4_flex`, in Rust, instead of
  liblz4: the frames keep the `lz4` tool's layout (and the tool reads them), but compressed
  blocks are no longer its bytes. braid is built without mimalloc. libdeflate is the one C
  library left.
## In lrpipe 0.1.0

- **Part of lrpipe**: lrpipe-snv moved, with its history, into the lrpipe workspace (`lrpipe-snv/`
  there), on lrpipe's pasta and a braid on the same pasta. `lrpipe snv` runs it (`run` implied)
  and the pipeline uses that; the `lrpipe-snv` binary still builds on its own. Its checks are part
  of lrpipe's `tools/check.sh`, its releases part of lrpipe's bundle, and the minimum Rust is
  lrpipe's, 1.91.
- **libtorch removed**: both networks run in lrpipe-snv's own Rust code, on the CPU (the `gemm`
  crate) and, with the `gpu` feature (on by default), on NVIDIA GPUs (`--device cuda:N`): CUDA
  kernels compiled at start-up by NVRTC and launched through the CUDA driver, both loaded at run
  time, so building needs neither libtorch nor CUDA. The default binary links only libc, libm and
  libgcc_s; a GPU needs the NVIDIA driver and `libnvrtc.so.12` at run time (found on
  the library path, under `$CUDA_HOME`, `$CUDA_PATH` or `/usr/local/cuda`, or given with the new
  `--nvrtc PATH`) and compute capability 6.0 (Pascal) or later. `cargo build --no-default-features` builds a CPU-only binary. The
  `tch` feature and `--backend` are gone.
- **`--approximate` no longer changes GPU runs**: the networks run in `float32` there. The
  `float16` tensor-core kernels and the `float32` rescoring of pileup candidates near a quality
  cut-off are gone: with lrpipe-snv's own kernels they saved 5% of a whole HG003 ONT lrpipe run
  from a BAM and nothing from a CRAM. On the CPU `--approximate` is unchanged.
- **Output no longer byte-identical to Clair3's**: the networks' sums run in another order than
  libtorch's, and on the GPU exact mode rounds the operands libtorch's cuDNN runs in TF32 the same
  way. Against the former libtorch backend on HG002 ONT (16 threads, i5-12600K, RTX 5060 Ti), the
  PASS calls agree as follows (times: Rust against libtorch):
  - CPU, chr20: 100% (144,477), the same GIAB recall; 131 s against 111 s (1,696 against 1,312
    CPU-s); with `--approximate` 99.994%, 83 s against 79 s.
  - GPU, chr20: 99.945%, GIAB F1 the same to five decimals; 37.6 s against 44.2 s, pileup
    calling 7.8 s against 12.8 s.
  - GPU, whole genome: 99.827%, the differences mostly from phasing (one het SNP on the other side
    of the phasing cut-off re-phases blocks); 602 s against 771 s, pileup calling 188 s against
    281 s.
- `tools/option_agreement.sh` reruns every option case Clair3 was run with and compares the VCFs
  by PASS calls. Against Clair3 2.0.3, PASS calls agree on 99.984% of the whole HG003 ONT
  genome on the CPU, 100% of HG002 HiFi chr20 on the CPU and 99.469% on the GPU, and 100% of
  every ONT, HiFi and Illumina option case (one `full_alignment.vcf.gz` at 99.987%).
- The CPU networks' tests compare with libtorch's probabilities for synthetic weights, saved
  before its removal (`crates/snv-infer/tests/data`).
- `ci/parity.sh` runs the fixture from BAM and CRAM: the CRAM must give the BAM's VCFs exactly,
  and each VCF must be the same as Clair3's or agree on at least 99.9% of its PASS calls.
  The checks lint and test with and without the `gpu` feature, and run the GPU tests when NVRTC
  and a GPU are found (and then must run them).
- `env.sh` only puts the original's env on `PATH`, for the models and the parity tools;
  `ci/environment.yml` is removed.
- `lrpipe-snv pileup --threads` is accepted for Clair3's job command lines; the networks score on
  one thread.

## 0.1.0 (its own repository, DylanCosto/lrpipe-snv)

First release: a Rust port of Clair3 2.0.3 whose output is identical to Clair3's.

- `lrpipe-snv run` takes `run_clair3.sh`'s options, with dashes or exactly as Clair3 spells
  them, and writes Clair3's output layout. A link named `run_clair3.sh` runs existing Clair3
  command lines unchanged. It runs the whole pipeline in one process: pileup calling, het-SNP
  selection, WhatsHap phasing and haplotagging (braid), full-alignment calling, merging
  and gVCF. LongPhase is run as an external binary, as Clair3 does.
- **Platforms**: ONT, PacBio HiFi and Illumina (with the short-read realigner).
- **Devices**: CPU and CUDA (`--device cuda:0`).
- **Parity**: on GIAB HG002 the output files match Clair3's byte for byte: chr20 on all three
  platforms, and the whole genome for HiFi and ONT on the CPU. On the GPU the final VCFs match
  ([VALIDATION.md](VALIDATION.md)).
- **Speed**: about 10x faster than Clair3 on ONT and HiFi chr20, and 2.4-2.8x on a whole genome
  on the CPU and 5x on the GPU.
- `--approximate` (lrpipe-snv only) trades identical output for speed:
  - fp16 on the GPU and 8-bit full-alignment convolutions on the CPU;
  - calls near a quality cut-off are rescored in fp32;
  - PASS calls agree with exact mode at 99.5-99.99% on HG003.
- **CRAM** input, including CRAM 3.1, which Clair3 2.0.3's bundled htslib cannot decode.
- **Release files**: an x86_64 Linux CPU bundle with its libtorch and MKL (glibc 2.28 or later)
  and a container image at `ghcr.io/dylancosto/lrpipe-snv`. Neither includes Clair3's models:
  pass `--model-path`.

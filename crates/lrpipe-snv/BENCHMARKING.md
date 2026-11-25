# Benchmarking on the Linux VM

How to re-run the parity, accuracy and speed checks of [VALIDATION.md](VALIDATION.md) on the
current code, and what to add to them. The setup and data are VALIDATION.md's.

## What to measure

VALIDATION.md was run at `282d716`, when lrpipe-snv ran the networks through libtorch and wrote
the same files as Clair3. Since then:

| change | commit | what checks it |
|---|---|---|
| SSE2 Smith-Waterman kernels | `a39cec7` | Illumina parity, Illumina whole-genome time ([4](#4-speed)) |
| `--chunk_num -1 --bed_fn` writes Clair3's `None` job; `--no_phasing_for_fa` links `.crai`/`.csi` | `5afa5a3` | option parity ([2](#2-parity-on-chr20)), CRAM input ([7](#7-cram-input)) |
| `--approximate` uses 8-bit full alignment on the Python pipeline (Illumina) | `2ad7744` | accuracy of `--approximate` ([5](#5-accuracy-of---approximate)) |
| freed heap given back to the system from phasing on (`malloc_trim` every 5 s) | `ea66598` | every parity run, peak memory ([4](#4-speed)) |
| `run_clair3.sh` command lines accepted unchanged | `1e97b32` | [8](#8-a-clair3-command-line-unchanged) |
| libtorch removed: the networks run in Rust on the CPU and CUDA GPUs | `c56d769` | every step: output is compared by PASS calls ([2](#2-parity-on-chr20), [3](#3-parity-on-the-whole-genome)), speed ([4](#4-speed)), the release bundle ([6](#6-the-release-bundle)) |

Measured since (README, [Against Clair3](README.md#against-clair3)): agreement with Clair3
without libtorch on the whole HG003 ONT genome on the CPU, HiFi chr20 on the CPU and GPU, and
every option case on ONT, HiFi and Illumina regions (`tools/option_agreement.sh`); and
`--approximate` against GIAB on HG002 ONT chr20 with `tools/concordance.py`. Still to measure,
on the VM, which has the data and hap.py:
- **Agreement on the whole genome for HiFi and Illumina**, against Clair3's own runs.
- **`--approximate` accuracy:** hap.py on `--approximate` output; VALIDATION.md has none yet.
- **Illumina whole-genome time:** the 10,331 s figure predates the SSE2 kernels.

## 0. Rules

- **One run at a time,** with nothing else heavy on the VM. Timings are only comparable between
  runs on the same machine.
- **Same threads as VALIDATION.md:** 16 for chr20, 20 for the whole genome on the CPU, 20 on the
  GPU.
- **Two runs per timing:** report both. The page cache changes the second run, so note whether
  the BAM was cached. To start cold, drop the caches as root: `sync; echo 3 > /proc/sys/vm/drop_caches`.
- **Record the build:** `git rev-parse --short HEAD`, and `lscpu | grep -o 'avx_vnni\|avx512_vnni' | sort -u`
  (lrpipe-snv's 8-bit convolutions need one of them; with `avx512_vnni` they use 512-bit
  registers, and the `float32` products use `gemm`'s AVX-512 kernels on any CPU with `avx512f`).
- **Wrap every timed run** in `tools/bench_run.sh OUT COMMAND...`. It saves the output, `/usr/bin/time -v`,
  vmstat, RSS and GPU samples next to `OUT`, and prints the wall time, CPU time, peak RSS,
  swap-out, GPU use and lrpipe-snv's stage times. Its GPU numbers include anything else on the GPU.

## 1. Setup

`env.sh`'s paths are this repository's development machine's; set `TORCH_ENV` as on the VM. It
is the original's env (Clair3 2.0.3 with PyTorch 2.13.0), for Clair3's models and the runs of
Clair3; building lrpipe-snv needs nothing from it. One binary runs on the CPU and the GPU; on the
GPU it needs the NVIDIA driver and NVRTC (`libnvrtc.so.12`, from CUDA or the `nvidia-cuda-nvrtc`
wheel): on the library path, under `$CUDA_HOME` or `/usr/local/cuda`, or named with `--nvrtc PATH`. `tools/bench_run.sh` needs GNU time (`apt install time`) for CPU time
and peak RSS.

```bash
cd lrpipe && git pull && git rev-parse --short HEAD && cd lrpipe-snv
source env.sh                     # TORCH_ENV=<Clair3 env with pytorch 2.13.0>
cargo build --release -p lrpipe-snv     # ../target/release/lrpipe-snv, CPU and GPU
SNV_MODEL=$TORCH_ENV/bin/models/r1041_e82_400bps_sup_v500 ../tools/check.sh --no-msrv
```

`tools/check.sh` (lrpipe's) runs every test, the GPU ones where the driver, a GPU and NVRTC are
found, and `ci/parity.sh` on the fixture.

The variables below stand for VALIDATION.md's files:

```bash
P=$PWD/../target/release/lrpipe-snv
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta  UCSC_REF=hg38.analysisSet.fa
M=$TORCH_ENV/bin/models
HG002_ONT=...  HG002_HIFI=...  HG002_ILMN=...  HG003_ONT=...  HG003_ILMN=...
WORK=/path/with/space/bench  # the whole-genome outputs take about 2 GB each
CHR=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
```

## 2. Parity on chr20

`tools/option_parity.sh` runs Clair3 and lrpipe-snv with the same options and compares every
file. With `SKIP_CLAIR3=1` it reuses VALIDATION.md's Clair3 runs under `$WORK/<name>/clair3`.
For Illumina, `CLAIR3_RUN` points at the Clair3 tree with the realigner built
(VALIDATION.md, [Illumina oracle](VALIDATION.md#illumina-oracle)).

```bash
export WORK THREADS=16
SKIP_CLAIR3=1 tools/option_parity.sh ont_chr20  ont  $HG002_ONT  $REF      $M/r1041_e82_400bps_sup_v500 chr20
SKIP_CLAIR3=1 tools/option_parity.sh hifi_chr20 hifi $HG002_HIFI $UCSC_REF $M/hifi_revio chr20
SKIP_CLAIR3=1 CLAIR3_RUN=<clair3 with realigner>/run_clair3.py \
  tools/option_parity.sh ilmn_chr20 ilmn $HG002_ILMN $REF $M/ilmn chr20
```

Each ends in `compared N files, M differ` (48, 40 and 43 files in VALIDATION.md, where none
differed with libtorch). Without libtorch, every file that depends on the networks' output can
differ; the files written before them (`CONTIGS`, `CHUNK_LIST`, split BEDs) should not. Compare
the PASS calls of each top-level VCF (`merge_output`, `pileup`, `full_alignment`), or rerun every
case of `$WORK` that has a Clair3 run, by the options in its VCF header, with
`tools/option_agreement.sh OUTDIR` (`WORK=<the Illumina cases>` for those):

```bash
for name in ont_chr20 hifi_chr20 ilmn_chr20; do
  for vcf in merge_output pileup full_alignment; do
    echo "$name $vcf"
    python3 tools/concordance.py $WORK/$name/clair3/$vcf.vcf.gz $WORK/$name/lrpipe-snv/$vcf.vcf.gz
  done
done
```

`ci/parity.sh` holds the fixture to an agreement and Jaccard of at least 0.999; on HG002 ONT chr20
the CPU run's PASS calls agreed 100% with the libtorch build's. Then the option cases that
changed:

```bash
BED=<a BED with a few chr20 regions>
SKIP_CLAIR3=1 tools/option_parity.sh ont_chunk_neg_bed ont $HG002_ONT $REF $M/r1041_e82_400bps_sup_v500 chr20 \
  --chunk_num=-1 --bed_fn=$BED --gvcf
SKIP_CLAIR3=1 tools/option_parity.sh ont_no_phasing ont $HG002_ONT $REF $M/r1041_e82_400bps_sup_v500 chr20 \
  --no_phasing_for_fa
```

These need their Clair3 runs, so run them once without `SKIP_CLAIR3=1`.

## 3. Parity on the whole genome

As in VALIDATION.md, with `THREADS=20`:

```bash
THREADS=20 SKIP_CLAIR3=1 tools/option_parity.sh ont_wgs  ont  $HG002_ONT  $REF      $M/r1041_e82_400bps_sup_v500 $CHR
THREADS=20 SKIP_CLAIR3=1 tools/option_parity.sh hifi_wgs hifi $HG002_HIFI $UCSC_REF $M/hifi_revio $CHR
THREADS=20 SKIP_CLAIR3=1 tools/option_parity.sh ont_wgs_gpu ont $HG002_ONT $REF \
  $M/r1041_e82_400bps_sup_v500 $CHR --use_gpu --device=cuda:0
```

With libtorch no file differed on the CPU (of 1,633 and 1,236), and on the GPU only six
intermediates that held the same records. Now compare the PASS calls as in step 2. On the GPU
the whole HG002 ONT genome agreed 99.827% with the libtorch build's, the differences mostly from
phasing (one het SNP on the other side of the phasing cut-off re-phases blocks); the CPU runs
have not been compared on the whole genome yet.

## 4. Speed

Time lrpipe-snv alone on HG003, exact and `--approximate`, two runs each.

```bash
run() { # name, then lrpipe-snv's options
  local name=$1; shift
  tools/bench_run.sh $WORK/time/$name "$@" --output $WORK/time/$name.out
}
ONT="--bam $HG003_ONT --ref $REF --platform ont --model-path $M/r1041_e82_400bps_sup_v500 --ctg-name $CHR"
ILMN="--bam $HG003_ILMN --ref $REF --platform ilmn --model-path $M/ilmn --ctg-name $CHR"

run ont_cpu          $P  run $ONT  --threads 20
run ont_cpu_approx   $P  run $ONT  --threads 20 --approximate
run ont_gpu          $P  run $ONT  --threads 20 --device cuda:0
run ilmn_cpu         $P  run $ILMN --threads 20
run ilmn_cpu_approx  $P  run $ILMN --threads 20 --approximate
```

Compare with VALIDATION.md, measured with libtorch: ONT CPU 1,403 s on 22 threads, Illumina
10,331 s on 20 (8,953 s of it full alignment), and Clair3's own times there. On the development
machine (i5-12600K, RTX 5060 Ti, 16 threads) the Rust networks took 1.3x libtorch's CPU time on
HG002 ONT chr20, and the whole HG002 ONT genome on the GPU in 602 s against 771 s (README,
Networks). Report the wall time, CPU time, peak RSS (VALIDATION.md has none; on the development
machine the ONT GPU run went from 16.7 to 11.7 GB with `malloc_trim`, and was 14.6 GB with the
Rust networks against 12.4 GB with libtorch) and the stage times `bench_run.sh` prints.

## 5. Accuracy of `--approximate`

For each run of step 4, `merge_output.vcf.gz` against GIAB, with VALIDATION.md's command:

```bash
for name in ont_cpu ont_cpu_approx ont_gpu ilmn_cpu ilmn_cpu_approx; do
  hap.py HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz $WORK/time/$name.out/merge_output.vcf.gz \
    -f HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -r $REF \
    --engine vcfeval --engine-vcfeval-template GRCh38_no_alt.sdf \
    --stratification GRCh38@all/GRCh38-all-stratifications.tsv -o $WORK/happy/$name
  python3 tools/happy_summary.py $WORK/happy/$name.extended.csv
done
```

And how closely each approximate run agrees with its exact one (PASS calls, same position,
alleles and genotype):

```bash
python3 tools/concordance.py $WORK/time/ont_cpu.out/merge_output.vcf.gz  $WORK/time/ont_cpu_approx.out/merge_output.vcf.gz
python3 tools/concordance.py $WORK/time/ilmn_cpu.out/merge_output.vcf.gz $WORK/time/ilmn_cpu_approx.out/merge_output.vcf.gz
```

The README's figures, agreement above 99% (99.5% on the CPU on ONT chr20, 99.95% on the Illumina
fixture) with F1 unchanged, were measured with libtorch; these runs measure them with the Rust
networks. `--approximate` changes only CPU runs.

## 6. The release bundle

lrpipe's bundle carries lrpipe-snv. Build it in `manylinux_2_28` (needs Docker), then run
lrpipe-snv from it in a shell without the conda env:

```bash
cd ..   # lrpipe
docker run --rm -v $PWD:/src -w /src -e CARGO_TARGET_DIR=/src/target-release \
  quay.io/pypa/manylinux_2_28_x86_64 \
  bash -c "tools/release_build.sh 0.1.9; s=\$?; chown -R $(id -u):$(id -g) dist target-release; exit \$s"
env -i PATH=/usr/bin:/bin HOME=$HOME \
  lrpipe-snv/ci/parity.sh "dist/lrpipe-0.1.9-x86_64-linux/bin/lrpipe snv" $M/r1041_e82_400bps_sup_v500
```

`ci/parity.sh` needs `python3` on `PATH` for `tools/concordance.py`. For each VCF it prints either
"same as the expected output" or the PASS-call agreement and Jaccard, and fails below 0.999; note
which. The same bundle is meant for the GPU, so check that it runs there too, given NVRTC:

```bash
env -i PATH=/usr/bin:/bin HOME=$HOME \
  dist/lrpipe-0.1.9-x86_64-linux/bin/lrpipe snv --nvrtc <path to libnvrtc.so.12> \
  --bam lrpipe-snv/ci/fixture/ont.bam --ref lrpipe-snv/ci/fixture/chr20_window.fa.gz --platform ont \
  --model-path $M/r1041_e82_400bps_sup_v500 --ctg-name chr20 --device cuda:0 --output $WORK/bundle_gpu
```

## 7. CRAM input

The HG002 Revio input is CRAM 3.1 (VALIDATION.md converted it to BAM for Clair3). Run lrpipe-snv
on the CRAM and on the BAM, and compare:

```bash
run hifi_cram $P run --bam <HG002 Revio CRAM> --ref $UCSC_REF --platform hifi --model-path $M/hifi_revio --ctg-name chr20 --threads 16
run hifi_bam  $P run --bam <HG002 Revio BAM>  --ref $UCSC_REF --platform hifi --model-path $M/hifi_revio --ctg-name chr20 --threads 16
cmp <(zcat $WORK/time/hifi_cram.out/merge_output.vcf.gz | grep -v '^##') \
    <(zcat $WORK/time/hifi_bam.out/merge_output.vcf.gz  | grep -v '^##') && echo same
```

## 8. A Clair3 command line unchanged

Any `run_clair3.sh` command you use, run through a link named `run_clair3.sh`:

```bash
mkdir -p $WORK/link && ln -sf $P $WORK/link/run_clair3.sh
$WORK/link/run_clair3.sh --bam_fn=$HG002_ONT --ref_fn=$REF --threads=16 --platform=ont \
  --model_path=$M/r1041_e82_400bps_sup_v500 --output=$WORK/link/out --ctg_name=chr20
```

Its `merge_output.vcf.gz` should equal step 2's `ont_chr20` lrpipe-snv output, apart from
`##cmdline`.

## Report

For VALIDATION.md, per run:
- commit, CPU flags, threads, and whether the cache was warm;
- wall time, CPU time and peak RSS;
- parity: files compared and differing, and the PASS-call agreement and Jaccard of
  `merge_output`, `pileup` and `full_alignment` against Clair3;
- hap.py: SNP and INDEL F1, from `happy_summary.py`'s `all` row;
- for approximate runs: agreement and Jaccard against the exact run.

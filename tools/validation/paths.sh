# The inputs, tools and truth sets of the whole-genome validation (BENCHMARKS.md), as
# environment variables. Each defaults to where it was on the machine that ran BENCHMARKS.md;
# set the variable to run elsewhere. Sourced by every script here.

: "${REPO:=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
# Where the outputs, logs and done-markers go.
: "${VALIDATION_DIR:=/home/dyl/bench-data/validation}"
: "${LRPIPE:=$REPO/target/release/lrpipe}"
: "${LRPIPE_TR_BENCH:=$REPO/target/release/lrpipe-tr-bench}"

# Reference: GRCh38 without ALT contigs (GRCh38_no_alt_analysis_set), with .fai.
: "${REF:=/data/bioinformatics/projects/ont-aligner/hifi-work/ref/GRCh38_no_alt.fa}"
# HG003 ONT runs' SV reference (the same sequence, as Sniffles2's reference run used it).
: "${SNIFFLES_REF:=/data2/bioinformatics/sniffles-rs-bench/reference.fa}"

# Reads (BENCHMARKS.md, Inputs).
: "${HG003_ONT_BAM:=/data/bioinformatics/sniffles-rs-bench/hg003-wgs/hg003.wgs.bam}"
: "${HG002_HIFI_BAM:=/data2/bioinformatics/sniffles-rs-bench/hifi-wgs/HG002.hifi.bonefishY.bam}"
: "${HG002_ONT_CRAM:=/home/dyl/bench-data/hg002-ont-wgs/SAMPLE.haplotagged.cram}"
: "${HG002_REVIO_BAM:=/home/dyl/bench-data/hg002-revio/HG002.m84011_220902_175841_s1.GRCh38.bam}"

# Clair3's model directories (pileup.pt and full_alignment.pt; README.md lists their SHA-256).
: "${MODELS:=/data2/bioinformatics/envs/clair3-t213/bin/models}"
: "${NVRTC:=/data2/bioinformatics/envs/torch213-cuda/lib/libnvrtc.so.12}"
# Tandem repeat catalog: adotto + STRchive, hg38 (937,204 loci).
# Set TR_CATALOG to catalog.adotto_strchive.hg38.bed for the tandem repeat steps.
: "${TR_CATALOG:=}"
# Spectre's reference run: inputs/hg002_hifi_wgs/mosdepth.regions.bed.gz, grch38_reference.mdr,
# grch38_blacklist_spectre.bed, and Spectre 0.2.1's source in spectre-0.2.1/.
: "${SPECTRE_BENCH:=/data2/bioinformatics/spectre-rs-bench}"

# The original tools.
: "${MOSDEPTH:=/data/bioinformatics/core/env/wes-tools/bin/mosdepth}"
: "${SPECTRE_PYTHON:=/data/bioinformatics/micromamba/envs/spectre-0.2.1/bin/python}"
# Set TRGT to a TRGT 5.1.0 binary for the TRGT comparison.
: "${TRGT:=}"
: "${MODKIT:=/data2/bioinformatics/methyl-bench/tools/modkit/dist_modkit_v0.6.4_cd85862/modkit}"
: "${PB_CPG:=/data2/bioinformatics/methyl-bench/tools/pbcpg/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores}"
# WhatsHap 2.8, bgzip and tabix (bioconda).
: "${CONDA_BIN:=/data2/bioinformatics/envs/clair3-t213/bin}"

# The original tools' outputs kept from earlier runs on the same inputs and machine.
: "${SNIFFLES_HG003_ONT_VCF:=/data/bioinformatics/sniffles-rs-bench/hg003-wgs/oracle.vcf}"
: "${SNIFFLES_HG002_CRAM_VCF:=/data2/bioinformatics/sniffles-rs-bench/ont-wgs/sniffles2.vcf}"
: "${SNIFFLES_HG002_HIFI_VCF:=/data2/bioinformatics/sniffles-rs-bench/hifi-wgs/sniffles2.vcf}"
: "${CLAIR3_HG003_ONT_CPU:=/data/bioinformatics/pellucid-wgs/clair3-cpu/merge_output.vcf.gz}"
: "${CLAIR3_HG002_HIFI_GPU:=/home/dyl/bench-data/agreement/clair3-hifi-wg-gpu/merge_output.vcf.gz}"

# `lrpipe run`'s whole-genome outputs phased and scored here (pipeline-check runs of 2026-11-02).
: "${HG002_HIFI_RUN:=/home/dyl/bench-data/pipeline-check/hg002-hifi-early}"
: "${HG002_ONT_RUN:=/home/dyl/bench-data/pipeline-check/hg002-ont-cram-early}"

# Truth sets and scoring: GIAB v4.2.1 (truth/, truth-hg003/), GIAB stratifications v3.x
# (strat/GRCh38@all), GIAB HG002 Q100 v5.0q (truth-q100/), RTG Tools 3.13 (rtg/), and the
# scoring's working files (ref.sdf, split/, eval/).
: "${GIAB_DIR:=/home/dyl/bench-data/giab}"
: "${GIAB_TR_TIER1:=/data2/bioinformatics/trcall-bench/giab_tr_v1.0.1/HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz}"
: "${RTG:=$GIAB_DIR/rtg/rtg-tools-3.13/rtg}"

export REPO VALIDATION_DIR LRPIPE REF GIAB_DIR RTG

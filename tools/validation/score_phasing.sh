#!/usr/bin/env bash
# Phasing accuracy against GIAB's HG002 Q100 benchmark (v5.0q, GRCh38 small variants), whose
# genotypes are phased chromosome-wide from the diploid assembly. Each call set is reduced to
# PASS calls on chr1-22 inside the benchmark regions; whatshap compare then counts switch and
# flip errors over the variants heterozygous in both, and whatshap stats the phased blocks.
#
# usage: score_phasing.sh NAME PHASED.vcf[.gz] [NAME PHASED.vcf[.gz]]...
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/paths.sh"
V=$VALIDATION_DIR
Q=$GIAB_DIR/truth-q100
E=$CONDA_BIN
RTG="env RTG_MEM=16g $RTG"
FAI=$REF.fai
O=$V/out/phase-score
mkdir -p "$O"
autosomes='^chr([1-9]|1[0-9]|2[0-2])$'

# restrict IN OUT: PASS calls on chr1-22 inside the benchmark regions.
restrict() {
  { zcat -f "$1" | grep '^#'
    zcat -f "$1" | grep -v '^#' | awk -v re="$autosomes" '$1 ~ re && ($7 == "PASS" || $7 == ".")'
  } | $E/bgzip > "$2.all.vcf.gz"
  $E/tabix -f -p vcf "$2.all.vcf.gz"
  $RTG vcffilter -i "$2.all.vcf.gz" -o "$2" --include-bed "$Q/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
  rm "$2.all.vcf.gz" "$2.all.vcf.gz.tbi"
}

[[ -e $O/truth.vcf.gz ]] || restrict "$Q/HG002_GRCh38_v5.0q_smvar.vcf.gz" "$O/truth.vcf.gz"
awk -v re="$autosomes" '$1 ~ re' "$FAI" | cut -f1,2 > "$O/autosomes.len"

while (( $# >= 2 )); do
  name=$1 vcf=$2
  shift 2
  restrict "$vcf" "$O/$name.vcf.gz"
  $E/whatshap compare --names truth,"$name" --tsv-pairwise "$O/$name.compare.tsv" \
    "$O/truth.vcf.gz" "$O/$name.vcf.gz" > "$O/$name.compare.txt" 2>&1
  $E/whatshap stats --chr-lengths "$O/autosomes.len" --tsv "$O/$name.stats.tsv" \
    "$O/$name.vcf.gz" > "$O/$name.stats.txt" 2>&1
done

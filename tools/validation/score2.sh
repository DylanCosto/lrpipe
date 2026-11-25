#!/usr/bin/env bash
# Scores whole-genome calls against GIAB v4.2.1 with rtg vcfeval: SNPs and other variants apart,
# PASS calls only, within the benchmark's confident regions (chr1-22).
# usage: score2.sh SAMPLE(HG002|HG003) NAME CALLS.vcf.gz
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/paths.sh"
G=$GIAB_DIR
RTG="env RTG_MEM=16g $RTG"
sample=$1 name=$2 calls=$3
case $sample in
  HG002) T=$G/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz BED=$G/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ;;
  HG003) T=$G/truth-hg003/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz BED=$G/truth-hg003/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ;;
  *) echo "unknown sample $sample" >&2; exit 2 ;;
esac
[[ -d $G/ref.sdf ]] || $RTG format -o $G/ref.sdf $REF > /dev/null
mkdir -p $G/split $G/eval
for kind in snps non-snps; do
  [[ -e $G/split/truth-$sample.$kind.vcf.gz ]] || $RTG vcffilter -i $T -o $G/split/truth-$sample.$kind.vcf.gz --$kind-only
  $RTG vcffilter -i "$calls" -o $G/split/$name.$kind.vcf.gz --$kind-only > /dev/null 2>&1
  rm -rf "${G:?}/eval/$name.$kind"
  $RTG vcfeval -b $G/split/truth-$sample.$kind.vcf.gz -c $G/split/$name.$kind.vcf.gz -e $BED -t $G/ref.sdf \
    -o $G/eval/$name.$kind --threads 16 > /dev/null 2>&1
done
python3 -I "$(dirname "$0")"/giab_strat.py $G/eval/$name.snps $G/eval/$name.non-snps > $G/eval/$name.strat.tsv
cat $G/eval/$name.strat.tsv

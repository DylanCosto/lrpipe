#!/usr/bin/env bash
# score2.sh with touchstone vcfeval (RTG 3.13's port, in Rust) in place of RTG's: the same
# split (RTG vcffilter, streaming, 2 GB heap), the same regions and stratification script.
# usage: score_ts.sh SAMPLE NAME CALLS.vcf.gz EVAL_DIR
set -euo pipefail
. /home/dyl/Documents/lrpipe/tools/validation/paths.sh
TS=/data2/bioinformatics/lrpipe-build/target-touchstone/release/touchstone
G=$GIAB_DIR
sample=$1 name=$2 calls=$3 E=$4
case $sample in
  HG002) T=$G/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz BED=$G/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ;;
  HG003) T=$G/truth-hg003/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz BED=$G/truth-hg003/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ;;
esac
mkdir -p $E/split
for kind in snps non-snps; do
  split=$G/split/$name.$kind.vcf.gz
  [[ -e $split ]] || env RTG_MEM=2g $RTG vcffilter -i "$calls" -o $E/split/$name.$kind.vcf.gz --$kind-only > /dev/null 2>&1
  [[ -e $split ]] || split=$E/split/$name.$kind.vcf.gz
  rm -rf "$E/$name.$kind"
  /usr/bin/time -v -o $E/$name.$kind.time $TS vcfeval -b $G/split/truth-$sample.$kind.vcf.gz -c $split -e $BED \
    -t $REF -o $E/$name.$kind -T 8 > /dev/null 2>&1
done
python3 -I /home/dyl/Documents/lrpipe/tools/validation/giab_strat.py $E/$name.snps $E/$name.non-snps > $E/$name.strat.tsv
cat $E/$name.strat.tsv

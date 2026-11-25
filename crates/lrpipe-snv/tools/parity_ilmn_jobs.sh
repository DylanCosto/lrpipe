#!/usr/bin/env bash
# Runs Clair3's own Python on one Illumina chunk and one full-alignment job of a Clair3 run and
# compares snv-mpileup's and snv-realign's intermediate outputs line by line: the
# `samtools mpileup` bases column, CreateTensorPileup's tensors, RealignReads' SAM lines and
# CreateTensorFullAlignment's tensors (the sequence column compared at its centre base).
#
# usage: parity_ilmn_jobs.sh CLAIR3_BIN RUN_DIR BAM REF CONTIG CHUNK_ID CHUNK_NUM BED OUTDIR
#   CLAIR3_BIN: a Clair3 tree whose preprocess/realign holds built realigner/debruijn_graph
#   libraries; RUN_DIR: its run directory (for tmp/phase_output/phase_bam/CONTIG.bam).
# env: TORCH_ENV (samtools, pypy3).
set -eu
C=$1 RUN=$2 BAM=$3 REF=$4 CTG=$5 ID=$6 NUM=$7 BED=$8 OUT=$9
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
HERE=$(cd "$(dirname "$0")/.." && pwd)
(cd "$HERE" && cargo build --release -p snv-mpileup -p snv-realign --examples)
X=$HERE/../../target/release/examples
export PATH=$TORCH_ENV/bin:$PATH
mkdir -p "$OUT"
PBAM=$RUN/tmp/phase_output/phase_bam/$CTG.bam
pypy3 "$C/clair3.py" CreateTensorPileup --bam_fn "$BAM" --ref_fn "$REF" --ctgName "$CTG" --chunk_id "$ID" \
  --chunk_num "$NUM" --platform ilmn --snp_min_af 0.08 --indel_min_af 0.08 --minMQ 5 --minCoverage 2 \
  --samtools samtools --min_af 0.08 2>/dev/null > "$OUT/py_pileup.txt"
"$X/mpileup_dump" pileup "$BAM" "$REF" "$CTG" "$ID" "$NUM" > "$OUT/rs_pileup.txt"
pypy3 "$C/clair3.py" RealignReads --bam_fn "$PBAM" --ref_fn "$REF" --ctgName "$CTG" --samtools samtools \
  --full_aln_regions "$BED" 2>/dev/null > "$OUT/py_reads.sam"
"$X/realign_dump" reads "$PBAM" "$REF" "$CTG" "$BED" > "$OUT/rs_reads.sam"
pypy3 "$C/clair3.py" CreateTensorFullAlignment --bam_fn PIPE --ref_fn "$REF" --ctgName "$CTG" --min_af 0.08 \
  --platform ilmn --samtools samtools --minCoverage 2 --minMQ 5 --phasing_info_in_bam --full_aln_regions "$BED" \
  < "$OUT/py_reads.sam" 2>/dev/null > "$OUT/py_fa.txt"
"$X/realign_dump" tensors "$PBAM" "$REF" "$CTG" "$BED" > "$OUT/rs_fa.txt"
centre() { awk -F'\t' '{print $1"\t"$2"\t"substr($3,17,1)"\t"$4"\t"$5}' "$1"; }
bad=0
cmp -s <(centre "$OUT/py_pileup.txt") "$OUT/rs_pileup.txt" || { echo "DIFFERS: pileup tensors"; bad=1; }
cmp -s <(grep -v '^@' "$OUT/py_reads.sam") "$OUT/rs_reads.sam" || { echo "DIFFERS: realigned reads"; bad=1; }
cmp -s <(centre "$OUT/py_fa.txt") "$OUT/rs_fa.txt" || { echo "DIFFERS: full-alignment tensors"; bad=1; }
echo "pileup $(wc -l < "$OUT/rs_pileup.txt"), reads $(wc -l < "$OUT/rs_reads.sam"), full-alignment $(wc -l < "$OUT/rs_fa.txt")"
exit $bad

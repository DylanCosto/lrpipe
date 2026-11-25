#!/usr/bin/env bash
# Rebuilds ci/fixture: an ONT BAM of chr20:30,200,000-30,240,000 and a chr20 that keeps the
# reference only where those reads lie (30,140,001-30,310,000) and is N elsewhere, bgzipped, so
# the contig keeps its length and coordinates. Then runs Clair3 on it and stores its VCFs, and from
# a run with --gvcf its gVCF, without the header lines that name the command and the reference
# path, as ci/fixture/expected.
#
# usage: make_ci_fixture.sh ONT_BAM REF TORCH_ENV
#   ONT_BAM covers the region (HG003, R10.4.1 sup v5.0.0 here); REF is GRCh38 with chr20.
set -euo pipefail
BAM=$1 REF=$2 TORCH_ENV=$3
HERE=$(cd "$(dirname "$0")/.." && pwd)
OUT=$HERE/ci/fixture
ST=$TORCH_ENV/bin/samtools
WORK=$(mktemp -d)
mkdir -p "$OUT/expected"

"$ST" view --no-PG -b -o "$WORK/ont.bam" "$BAM" chr20:30200000-30240000
# Drop the @PG lines, which name local paths.
"$ST" view --no-PG -H "$WORK/ont.bam" | grep -v '^@PG' > "$WORK/header.sam"
"$ST" reheader -P "$WORK/header.sam" "$WORK/ont.bam" > "$OUT/ont.bam"
"$ST" index "$OUT/ont.bam"

length=$(awk '$1 == "chr20" { print $2 }' "$REF.fai")
"$ST" faidx "$REF" chr20:30140001-30310000 | grep -v '>' | tr -d '\n' > "$WORK/keep.txt"
python3 -I - "$length" 30140000 "$WORK/keep.txt" > "$WORK/chr20_window.fa" <<'PY'
import sys
length, start, keep = int(sys.argv[1]), int(sys.argv[2]), open(sys.argv[3]).read()
seq = "N" * start + keep + "N" * (length - start - len(keep))
sys.stdout.write(">chr20\n")
for i in range(0, len(seq), 60):
    sys.stdout.write(seq[i:i + 60] + "\n")
PY
"$ST" faidx "$WORK/chr20_window.fa"
"$TORCH_ENV/bin/bgzip" -c "$WORK/chr20_window.fa" > "$OUT/chr20_window.fa.gz"
"$ST" faidx "$OUT/chr20_window.fa.gz"
# The same reads as CRAM, which ci/parity.sh also calls from (the window holds chr20 alone, so
# samtools embeds the reference).
"$ST" view --no-PG -C -T "$OUT/chr20_window.fa.gz" -o "$OUT/ont.cram" "$OUT/ont.bam"
"$ST" index "$OUT/ont.cram"

PATH=$TORCH_ENV/bin:$PATH "$TORCH_ENV/bin/python3" "$TORCH_ENV/bin/run_clair3.py" \
  --bam_fn="$OUT/ont.bam" --ref_fn="$WORK/chr20_window.fa" --threads=8 --platform=ont \
  --model_path="$TORCH_ENV/bin/models/r1041_e82_400bps_sup_v500" --output="$WORK/clair3" --ctg_name=chr20 \
  > "$WORK/clair3.log" 2>&1
for name in merge_output pileup full_alignment; do
  gzip -dc "$WORK/clair3/$name.vcf.gz" | grep -vE '^##(cmdline|reference)=' > "$OUT/expected/$name.vcf"
done
# And its gVCF, from a run with --gvcf (which adds PL to the VCFs above, so it is a run of its own).
# expected/merge_output.gvcf.fixes, how lrpipe-snv's gVCF differs where it fixes Clair3's bugs, is
# ci/parity.sh's `diff` of the two, checked by hand.
PATH=$TORCH_ENV/bin:$PATH "$TORCH_ENV/bin/python3" "$TORCH_ENV/bin/run_clair3.py" \
  --bam_fn="$OUT/ont.bam" --ref_fn="$WORK/chr20_window.fa" --threads=8 --platform=ont \
  --model_path="$TORCH_ENV/bin/models/r1041_e82_400bps_sup_v500" --output="$WORK/clair3-gvcf" --ctg_name=chr20 \
  --gvcf > "$WORK/clair3-gvcf.log" 2>&1
gzip -dc "$WORK/clair3-gvcf/merge_output.gvcf.gz" | grep -vE '^##(cmdline|reference)=' | gzip -9n \
  > "$OUT/expected/merge_output.gvcf.gz"
rm -rf "$WORK"

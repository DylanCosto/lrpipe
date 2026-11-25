#!/usr/bin/env bash
# Runs `lrpipe-snv run` on ci/fixture, from the BAM and from the CRAM, and compares its VCFs with the
# ones the original (2.0.3) wrote for the same inputs (ci/fixture/expected), ignoring the header
# lines that name the command and the reference path. A third run, from the BAM with --gvcf,
# must write an indexed gVCF with Clair3's records but for the Clair3
# bugs lrpipe-snv fixes (expected/merge_output.gvcf.fixes). Use file order for tied read starts,
# as Clair3 did when writing those fixtures.
#
# lrpipe-snv's networks run in Rust, whose sums round differently from the original's libtorch, so
# the VCFs are held to near-perfect agreement of the PASS calls (position, alleles, genotype) rather
# than identical bytes; identical files are reported as such. The CRAM must give exactly what the
# BAM gives.
#
# usage: ci/parity.sh LRPIPE_SNV MODEL_DIR [OUTDIR]
#   LRPIPE_SNV is the lrpipe-snv binary, or "lrpipe snv" (words split at spaces).
#   MODEL_DIR is the original's r1041_e82_400bps_sup_v500 (pileup.pt and full_alignment.pt).
set -euo pipefail
read -r -a LRPIPE_SNV <<< "$1"
MODEL=$2
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${3:-$(mktemp -d)}
FIXTURE=$HERE/fixture
mkdir -p "$OUT"
for input in ont.bam ont.cram; do
  run=$OUT/run-${input##*.}
  rm -rf "$run"
  "${LRPIPE_SNV[@]}" run --bam "$FIXTURE/$input" --ref "$FIXTURE/chr20_window.fa.gz" --threads "${THREADS:-4}" \
    --platform ont --model-path "$MODEL" --output "$run" --ctg-name chr20 --read-order input > "$run.log" 2>&1 ||
    { tail -20 "$run.log"; exit 1; }
done
# And with a gVCF, as `lrpipe run --gvcf` asks for (its VCF records gain Clair3's PL).
run=$OUT/run-gvcf
rm -rf "$run"
"${LRPIPE_SNV[@]}" run --bam "$FIXTURE/ont.bam" --ref "$FIXTURE/chr20_window.fa.gz" --threads "${THREADS:-4}" \
  --platform ont --model-path "$MODEL" --output "$run" --ctg-name chr20 --read-order input --gvcf > "$run.log" 2>&1 ||
  { tail -20 "$run.log"; exit 1; }

# The gVCF is indexed and is the one Clair3 wrote for the same reads (expected/merge_output.gvcf.gz)
# but for the Clair3 bugs lrpipe-snv fixes, which expected/merge_output.gvcf.fixes lists as `diff`
# shows them: each record's position, alleles, filter, INFO and genotype (and a block's GQ,
# MIN_DP and PL, which no network computes), compared in order. QUAL and a call's GQ and PL are
# left out, as the networks' sums round differently.
gvcf=$OUT/run-gvcf/merge_output.gvcf.gz
for f in "$gvcf" "$gvcf.tbi"; do
  if [ ! -s "$f" ]; then echo "--gvcf did not write $f"; exit 1; fi
done
records() {
  gzip -dc "$1" | grep -v '^#' |
    awk -F'\t' '{ s = $5 == "<NON_REF>" ? $10 : substr($10, 1, index($10, ":") - 1); print $1, $2, $4, $5, $7, $8, s }'
}
if ! cmp -s "$FIXTURE/expected/merge_output.gvcf.fixes" <(diff <(records "$FIXTURE/expected/merge_output.gvcf.gz") <(records "$gvcf")); then
  echo "merge_output.gvcf.gz differs from Clair3's other than by the fixes in expected/merge_output.gvcf.fixes:"
  diff "$FIXTURE/expected/merge_output.gvcf.fixes" <(diff <(records "$FIXTURE/expected/merge_output.gvcf.gz") <(records "$gvcf")) | head -20
  exit 1
fi
echo "merge_output.gvcf.gz: indexed, Clair3's records but for the fixes listed ($(records "$gvcf" | wc -l) records)"

# The CRAM, which lrpipe-snv decodes only in part (no methylation tags, for one), gives the BAM's VCFs.
for name in merge_output pileup full_alignment; do
  if ! cmp -s <(gzip -dc "$OUT/run-bam/$name.vcf.gz" | grep -v '^##cmdline=') \
    <(gzip -dc "$OUT/run-cram/$name.vcf.gz" | grep -v '^##cmdline='); then
    echo "ont.cram: $name.vcf.gz differs from ont.bam's"
    exit 1
  fi
done

# The same fixture with no candidates selected for full alignment must still emit haploid calls
# and reference blocks. The CPU precision probe must not contribute reference blocks of its own.
# A tiny nonzero selection proportion avoids the preset default for zero.
run=$(mktemp -d "$OUT/haploid-pileup.XXXXXX")
RUST_LOG=info "${LRPIPE_SNV[@]}" run --bam "$FIXTURE/ont.bam" --ref "$FIXTURE/chr20_window.fa.gz" \
  --threads "${THREADS:-4}" --platform ont --model-path "$MODEL" --output "$run" --ctg-name chr20 \
  --read-order input --haploid-precise --no-phasing-for-fa --gvcf --approximate \
  --var-pct-full 0.000001 --ref-pct-full 0.000001 > "$run.log" 2>&1 ||
  { tail -20 "$run.log"; exit 1; }
if ! grep -q 'No Candidate found! Exit in selecting full-alignment candidates' "$run.log"; then
  echo "the haploid pileup fixture should skip full alignment"
  exit 1
fi
if ! grep -q '8-bit pileup network' "$run.log"; then
  echo "the haploid pileup fixture should probe the 8-bit network"
  exit 1
fi
for file in "$run/merge_output.vcf.gz" "$run/merge_output.gvcf.gz"; do
  if [ ! -s "$file.tbi" ]; then echo "missing index for $file"; exit 1; fi
  gzip -dc "$file" | awk -F '\t' '
    !/^#/ { split($10, sample, ":"); if (sample[1] !~ /^([0-9]+|\.)$/) bad = 1; ++records }
    END { exit (bad || !records) }
  ' || { echo "expected haploid calls in $file"; exit 1; }
done
gzip -dc "$run/merge_output.gvcf.gz" | awk -F '\t' '
  $5 == "<NON_REF>" {
    if ($1 != contig) { contig = $1; end = 0 }
    if ($2 <= end) bad = 1
    end = substr($8, 5) + 0
    ++blocks
  }
  END { exit (bad || !blocks) }
' || { echo "expected nonoverlapping reference blocks in the haploid gVCF"; exit 1; }
echo "haploid pileup-only fallback: indexed VCF and gVCF with nonoverlapping reference blocks"

for name in merge_output pileup full_alignment; do
  records=$(grep -vc '^#' "$FIXTURE/expected/$name.vcf")
  if cmp -s "$FIXTURE/expected/$name.vcf" <(gzip -dc "$OUT/run-bam/$name.vcf.gz" | grep -vE '^##(cmdline|reference)='); then
    echo "$name.vcf.gz: same as the expected output ($records records)"
    continue
  fi
  report=$(python3 -I "$HERE/../tools/concordance.py" "$FIXTURE/expected/$name.vcf" "$OUT/run-bam/$name.vcf.gz")
  read -r agreement jaccard < <(echo "$report" | awk '/^agreement/ { print $2, $4 }')
  if awk -v a="$agreement" -v j="$jaccard" 'BEGIN { exit !(a >= 0.999 && j >= 0.999) }'; then
    echo "$name.vcf.gz: PASS calls agree with the expected output $agreement (Jaccard $jaccard; $records records)"
  else
    echo "$name.vcf.gz: PASS calls agree with the expected output only $agreement (Jaccard $jaccard)"
    echo "$report"
    exit 1
  fi
done

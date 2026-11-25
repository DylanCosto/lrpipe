#!/usr/bin/env bash
# Rebuilds the test fixtures in tests/data.
#
# Real excerpts keep alignments, flags, mates and read groups but drop sequences, qualities and
# other tags, which depth does not read, so they stay small. Inputs are the clair3-rs-bench
# fixtures; override with HIFI, ONT, ILMN, ONTXY. Needs samtools (SAMTOOLS) and python3.
set -euo pipefail
cd "$(dirname "$0")/.."
SAMTOOLS=${SAMTOOLS:-samtools}
SRC=${SRC:-/data2/bioinformatics/clair3-rs-bench/fixtures}
HIFI=${HIFI:-$SRC/hifi.chr22_20-26m.bam}
ONT=${ONT:-$SRC/ont.chr20_30-36m.bam}
ILMN=${ILMN:-$SRC/ilmn.chr20_30-36m.bam}
ONTXY=${ONTXY:-$SRC/ont.chrXY.bam}
out=tests/data
mkdir -p "$out"

strip() { # in region out
    "$SAMTOOLS" view -h "$1" ${2:+"$2"} \
        | awk 'BEGIN{OFS="\t"} /^@/{if ($0 !~ /^@PG/) print; next}
               {rg=""; for (i=12;i<=NF;i++) if ($i ~ /^RG:Z:/) rg=$i; $10="*"; $11="*"; NF=11; if (rg!="") $12=rg; print}' \
        | "$SAMTOOLS" view -b --no-PG -o "$3"
    "$SAMTOOLS" index "$3"
}

strip "$HIFI" chr22:20000000-22000000 "$out/hifi.bam"
strip "$ONT" chr20:30000000-31500000 "$out/ont.bam"
strip "$ILMN" chr20:30000000-30400000 "$out/ilmn.bam"
strip "$ONTXY" "" "$out/ontxy.bam"

python3 tools/synth.py "$out"
"$SAMTOOLS" faidx "$out/synth.fa"
"$SAMTOOLS" sort --no-PG -o "$out/synth.bam" "$out/synth.sam"
"$SAMTOOLS" index "$out/synth.bam"
# The same records with a CSI index instead of BAI.
cp "$out/synth.bam" "$out/synth_csi.bam"
"$SAMTOOLS" index -c "$out/synth_csi.bam"
"$SAMTOOLS" view -C --no-PG -T "$out/synth.fa" -o "$out/synth.cram" "$out/synth.bam"
"$SAMTOOLS" index "$out/synth.cram"
rm "$out/synth.sam"
ls -la "$out"

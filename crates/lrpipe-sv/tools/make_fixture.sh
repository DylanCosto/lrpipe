#!/usr/bin/env bash
# Regenerates tests/data: a small slice of real reads plus the VCFs Sniffles2 2.8.1 writes for it,
# which tests/parity.rs compares lrpipe-sv against. Needs samtools, Python with pysam, and
# Sniffles 2.8.1; the result is checked in, so the tests themselves need none of them.
#
# usage: tools/make_fixture.sh ALIGNMENTS REFERENCE [REGION]
#   The fixture was cut from ONT's HG002 giab_2025.01 haplotagged CRAM (R10.4.1, sup) against
#   GRCh38 (no alt), REGION chr17:46200000-46400000.
#   SNIFFLES=/path/to/sniffles  SAMTOOLS=/path/to/samtools  PYTHON=/path/to/python (with pysam)
#   LRPIPE_SV=/path/to/lrpipe-sv  (for the --accuracy snapshot; default target/release)
set -euo pipefail

src=${1:?alignments}
reference=${2:?reference FASTA}
region=${3:-chr17:46200000-46400000}
env=/data/bioinformatics/micromamba/envs/vareval/bin
sniffles=${SNIFFLES:-$env/sniffles}
samtools=${SAMTOOLS:-$env/samtools}
python=${PYTHON:-$env/python}
here=$(cd "$(dirname "$0")" && pwd)
rs=${LRPIPE_SV:-$here/../../../target/release/lrpipe-sv}
data=$here/../tests/data
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

rm -rf "$data" && mkdir -p "$data/expected"
"$python" -I "$here/make_fixture.py" "$src" "$reference" "$region" "$data"

# The full-length contig the tests rebuild: N, then the window, then N.
"$python" -I - "$data/reference.window.fa" "$work/reference.fa" <<'EOF'
import sys
window, out = sys.argv[1:3]
with open(window) as f:
    name, ln = f.readline()[1:].split()
    seq = "".join(l.strip() for l in f)
contig, span = name.split(":")
start = int(span.split("-")[0]) - 1
length = int(ln.removeprefix("LN:"))
with open(out, "w") as f:
    f.write(f">{contig}\n")
    full = "N" * start + seq + "N" * (length - start - len(seq))
    for i in range(0, len(full), 60):
        f.write(full[i:i + 60] + "\n")
EOF
"$samtools" faidx "$work/reference.fa"
ref=$work/reference.fa

# CRAM against the external reference (only the fixture contig is in it, so no embedding).
for s in sample1 sample2; do
    "$samtools" view -C --output-fmt-option embed_ref=0 -T "$ref" -o "$data/$s.cram" "$data/$s.bam" 2> /dev/null
    "$samtools" index "$data/$s.cram"
done
# The second sample is only kept as CRAM, to keep the fixture small.
mv "$data"/sample2.bam* "$work/"

printf '%s\t%s\t%s\n' "${region%%:*}" 46250000 46300000 > "$data/regions.bed"
# lrpipe run --regions's (lrpipe/tests/run.rs): a DUP starting before the first region, an
# insertion in it and a deletion in the second.
printf '%s\t%s\t%s\n' "${region%%:*}" 46207000 46212000 "${region%%:*}" 46270000 46285000 > "$data/run_regions.bed"

# Inputs live in the work directory so that header lines naming files are the same as in tests.
cp "$data"/sample1.bam* "$data"/*.cram* "$data/regions.bed" "$data/run_regions.bed" "$work/"
cd "$work"
run() {
    local name=$1
    shift
    "$sniffles" --allow-overwrite -t 4 "$@" > "$name.log" 2>&1 || { cat "$name.log"; exit 1; }
}
run default --input sample1.bam --vcf default.vcf
run reference --input sample1.bam --vcf reference.vcf --reference "$ref"
run cram --input sample1.cram --vcf cram.vcf --reference "$ref"
run phase --input sample1.bam --vcf phase.vcf --phase
run mosaic --input sample1.bam --vcf mosaic.vcf --mosaic
run no_qc --input sample1.bam --vcf no_qc.vcf --no-qc
run symbolic --input sample1.bam --vcf symbolic.vcf --symbolic
run regions --input sample1.bam --vcf regions.vcf --regions regions.bed
run regions_run --input sample1.bam --vcf regions_run.vcf --regions run_regions.bed --reference "$ref" --sample-id S
run genotype --input sample2.cram --vcf genotype.vcf --reference "$ref" --genotype-vcf default.vcf
run snf1 --input sample1.bam --snf s1.snf --sample-id s1
run snf2 --input sample2.cram --snf s2.snf --sample-id s2 --reference "$ref"
run combine --input s1.snf s2.snf --vcf combine.vcf
for v in default reference cram phase mosaic no_qc symbolic regions regions_run genotype combine; do
    cp "$v.vcf" "$data/expected/"
done

# Not Sniffles2 output: a snapshot of lrpipe-sv --accuracy, so changes to it are deliberate.
"$rs" --allow-overwrite -t 4 --input sample1.bam --vcf accuracy.vcf --reference "$ref" --accuracy > /dev/null 2>&1
cp accuracy.vcf "$data/expected/"

ls -la "$data" "$data/expected"
du -sh "$data"

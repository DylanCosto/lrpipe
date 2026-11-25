#!/usr/bin/env bash
# Truvari (GIAB v5.0q, chr20) for lrpipe-sv on HG002 chr20 HiFi and ONT with extra options.
#
# usage: tools/eval_chr20.sh NAME [lrpipe-sv options...]
# Prints, per dataset: TP FP FN precision recall F1 for strict and --pick multi matching.
set -euo pipefail
name=$1
shift
E=/data/bioinformatics/micromamba/envs/vareval/bin
RS=$(dirname "$0")/../../../target/release/lrpipe-sv
REF=/data/bioinformatics/projects/ont-aligner/hifi-work/ref/GRCh38_no_alt.fa
SVT=/data/bioinformatics/projects/ont-aligner/hifi-work/truth_HG002_sv.vcf.gz
BED=/data/bioinformatics/projects/ont-aligner/hifi-data/truth/HG002_chr20_v5.0q_stvar.bed
OUT=/data2/bioinformatics/lrpipe-sv-bench/accuracy/$name
declare -A BAMS=(
    [hifi]=/data/bioinformatics/lrpipe-sv-bench/hifi/HG002.chr20.mm2Y.bam
    [ont]=/data/bioinformatics/projects/ont-aligner/ont-data/realign_eval/bonefish_dna/aln.bam
)
mkdir -p "$OUT"
for ds in hifi ont; do
    d=$OUT/$ds
    mkdir -p "$d"
    "$RS" -i "${BAMS[$ds]}" --vcf "$d/calls.vcf" --reference "$REF" -t 16 --allow-overwrite "$@" > "$d/rs.log" 2>&1
    $E/bcftools sort -Oz -o "$d/calls.vcf.gz" "$d/calls.vcf" 2> /dev/null
    $E/tabix -f -p vcf "$d/calls.vcf.gz"
    line="$name $ds"
    for mode in strict multi; do
        args=(--passonly --pctseq 0)
        [[ $mode == multi ]] && args+=(--pick multi)
        rm -rf "$d/tv_$mode"
        $E/truvari bench -b "$SVT" -c "$d/calls.vcf.gz" --includebed "$BED" "${args[@]}" -o "$d/tv_$mode" \
            > "$d/tv_$mode.log" 2>&1
        line+=$(python3 -c "import json; s=json.load(open('$d/tv_$mode/summary.json')); print(f\"  {'$mode'}: TP {s['TP-base']} FP {s['FP']} FN {s['FN']} P {s['precision']:.4f} R {s['recall']:.4f} F1 {s['f1']:.4f}\")")
    done
    echo "$line"
done

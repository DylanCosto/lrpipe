#!/usr/bin/env bash
# Runs tools/oracle_fullalign.py (direct mode) for every line of a parity cases.tsv, writing
# NAME.tensors.npy, NAME.info.tsv and NAME.seconds next to it.
#
# usage: tools/oracle_fullalign_cases.sh CASES_TSV [JOBS]
# env:   TORCH_ENV (Clair3 env)
set -euo pipefail
CASES=${1:?cases.tsv} JOBS=${2:-4}
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
HERE=$(cd "$(dirname "$0")/.." && pwd)
DIR=$(cd "$(dirname "$CASES")" && pwd)
export PATH="$TORCH_ENV/bin:$PATH" TORCH_ENV HERE DIR

one() {
    IFS=$'\t' read -r name bam ref contig bed vcf mq platform long <<< "$1"
    args=(--bam "$bam" --ref "$ref" --ctg-name "$contig" --full-aln-regions "$bed" --min-mq "$mq" --platform "$platform")
    if [ "$vcf" = "-" ]; then args+=(--no-phasing); else args+=(--phased-vcf "$vcf"); fi
    if [ "$long" = "1" ]; then args+=(--enable-long-indel); fi
    if "$TORCH_ENV/bin/python3" -I "$HERE/tools/oracle_fullalign.py" direct --clair3-dir "$TORCH_ENV/bin" \
        "${args[@]}" --out "$DIR/$name" > "$DIR/$name.log" 2>&1; then
        echo "$name: $(tail -1 "$DIR/$name.log")"
    else
        echo "$name: oracle FAILED (see $DIR/$name.log)"
    fi
}
export -f one
grep -v '^#' "$CASES" | grep . | xargs -d '\n' -P "$JOBS" -I{} bash -c 'one "$@"' _ {}

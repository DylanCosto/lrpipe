#!/usr/bin/env bash
# Reruns the Spectre reference runs in tools/reference_runs.tsv into $SPECTRE_BENCH/ref.
#   SPECTRE_BENCH=... SPECTRE_ENV=... SPECTRE_SRC=... REFERENCE=... tools/make_reference.sh [name...]
set -euo pipefail
cd "$(dirname "$0")/.."
: "${SPECTRE_BENCH:?}" "${SPECTRE_ENV:?}" "${SPECTRE_SRC:?}" "${REFERENCE:?}"
export SPECTRE_SRC PATH=$SPECTRE_ENV/bin:$PATH
tools=$PWD/tools
while IFS=$'\t' read -r name sample coverage snv only_chr extra; do
    [[ -z $name || $name == \#* ]] && continue
    if [[ $# -gt 0 && ! " $* " == *" $name "* ]]; then continue; fi
    seed=0
    [[ $name == *seed7 ]] && seed=7
    args=(CNVCaller -c "$SPECTRE_BENCH/$coverage" -s "$sample" -d out -r "$REFERENCE"
          -m "$SPECTRE_BENCH/inputs/grch38_reference.mdr" -l "$SPECTRE_BENCH/inputs/grch38_blacklist_spectre.bed")
    [[ $snv != - ]] && args+=(-v "$SPECTRE_BENCH/$snv")
    [[ $only_chr != - ]] && args+=(--only-chr "$only_chr")
    # shellcheck disable=SC2206
    [[ $extra != - ]] && args+=($extra)
    dir=$SPECTRE_BENCH/ref/$name
    rm -rf "$dir" && mkdir -p "$dir"
    (cd "$dir" && SPECTRE_NUMPY_SEED=$seed "$SPECTRE_ENV/bin/python" "$tools/run_spectre.py" dump -- "${args[@]}" > log.txt 2>&1) \
        && echo "ok   $name" || { echo "FAIL $name: $(tail -1 "$dir/log.txt")"; }
done < tools/reference_runs.tsv

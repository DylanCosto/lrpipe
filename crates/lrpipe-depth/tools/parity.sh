#!/usr/bin/env bash
# Runs mosdepth and lrpipe-depth on every case in tools/cases.tsv and compares their outputs
# (gzipped outputs after decompression). Prints a diff for each mismatch.
#
#   tools/parity.sh            compare
#   tools/parity.sh --record   also write mosdepth's output hashes to tests/data/expected.sha256,
#                              which `cargo test` checks lrpipe-depth against
#
# MOSDEPTH (default: mosdepth) must be mosdepth 0.3.14. With TABIX set, each regions.bed.gz
# lrpipe-depth writes is also queried through its CSI index and compared with a plain read.
set -euo pipefail
cd "$(dirname "$0")/.."
MOSDEPTH=${MOSDEPTH:-mosdepth}
record=0
[[ ${1:-} == --record ]] && record=1
"$MOSDEPTH" --version | grep -q 'mosdepth 0.3.14' || { echo "need mosdepth 0.3.14 (MOSDEPTH=...)" >&2; exit 2; }
cargo build --release -q
bin=../../target/release/lrpipe-depth
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
outputs=(mosdepth.summary.txt mosdepth.global.dist.txt mosdepth.region.dist.txt regions.bed.gz thresholds.bed.gz)
text() { if [[ $1 == *.gz ]]; then gzip -dc "$1"; else cat "$1"; fi; }
fail=0
n=0
expected=$work/expected.sha256
: > "$expected"
while IFS=$'\t' read -r name input prec opts; do
    [[ -z $name || $name == \#* ]] && continue
    n=$((n + 1))
    env=()
    [[ $prec != - ]] && env=(MOSDEPTH_PRECISION="$prec")
    # shellcheck disable=SC2086
    env "${env[@]}" "$MOSDEPTH" $opts "$work/m.$name" "tests/data/$input" 2> "$work/m.$name.err" \
        || { echo "FAIL $name: mosdepth failed: $(cat "$work/m.$name.err")"; fail=1; continue; }
    # shellcheck disable=SC2086
    env "${env[@]}" "$bin" -t 3 $opts "$work/d.$name" "tests/data/$input" 2> "$work/d.$name.err" \
        || { echo "FAIL $name: lrpipe-depth failed: $(cat "$work/d.$name.err")"; fail=1; continue; }
    ok=1
    for o in "${outputs[@]}"; do
        m=$work/m.$name.$o d=$work/d.$name.$o
        if [[ ! -e $m ]]; then
            [[ -e $d ]] && { echo "FAIL $name: lrpipe-depth wrote $o, mosdepth did not"; ok=0; }
            continue
        fi
        [[ -e $d ]] || { echo "FAIL $name: lrpipe-depth did not write $o"; ok=0; continue; }
        if ! cmp -s <(text "$m") <(text "$d"); then
            echo "FAIL $name: $o differs (< mosdepth, > lrpipe-depth)"
            diff <(text "$m") <(text "$d") | head -20 || true
            ok=0
        fi
        printf '%s\t%s\t%s\n' "$name" "$o" "$(text "$m" | sha256sum | cut -d' ' -f1)" >> "$expected"
        if [[ -n ${TABIX:-} && $o == *.bed.gz ]]; then
            for c in $(text "$d" | grep -v '^#' | cut -f1 | uniq | head -5); do
                cmp -s <("$TABIX" "$d" "$c") <(text "$d" | grep -v '^#' | awk -v c="$c" '$1 == c') \
                    || { echo "FAIL $name: $o index query of $c differs"; ok=0; }
            done
        fi
    done
    [[ $ok == 1 ]] && echo "ok   $name" || fail=1
done < tools/cases.tsv
if [[ $record == 1 ]]; then
    cp "$expected" tests/data/expected.sha256
    echo "recorded tests/data/expected.sha256"
fi
[[ $fail == 0 ]] && echo "all $n cases match" || { echo "mismatches"; exit 1; }

#!/usr/bin/env bash
# VCF I/O oracle: random VCFs and edits (tests/vcf/gen.py), replayed through
# pysam + WhatsHap (tests/vcf/replay.py) and through lrpipe-phase
# (examples/vcf_replay.rs). Compares what reading returns, the edited records,
# the region queries and the output VCF byte for byte.
#
# Usage: tests/vcf/fuzz.sh ENV_BIN WORK_DIR [FIRST_SEED [COUNT]]
#   ENV_BIN: bin directory of an environment with WhatsHap 2.8 (python, pysam)
# Build the replayer first: cargo build --release --example vcf_replay
set -u
ENV=$1; WORK=$2; FIRST=${3:-1}; COUNT=${4:-200}
HERE=$(cd "$(dirname "$0")" && pwd)
RS=${RS:-$HERE/../../../../target/release/examples/vcf_replay}
pass=0; fail=0; skip=0; failed=()
for seed in $(seq "$FIRST" $((FIRST + COUNT - 1))); do
    d=$WORK/$seed
    rm -rf "$d"; mkdir -p "$d"
    "$ENV/python" "$HERE/gen.py" "$seed" "$d"
    if ! "$ENV/python" "$HERE/replay.py" "$d" > "$d/expected.txt" 2> "$d/py.err"; then
        echo "SKIP $seed: pysam crashed"; skip=$((skip+1)); continue
    fi
    "$RS" "$d" > "$d/actual.txt" 2> "$d/rs.err"
    ok=1
    if ! cmp -s "$d/expected.txt" "$d/actual.txt"; then
        echo "FAIL $seed: reads differ"; diff "$d/expected.txt" "$d/actual.txt" | head -6 | sed 's/^/    /'
        ok=0
    fi
    if [ -f "$d/out.vcf" ] || [ -f "$d/rs.vcf" ]; then
        if ! cmp -s "$d/out.vcf" "$d/rs.vcf"; then
            echo "FAIL $seed: output VCF differs"; diff "$d/out.vcf" "$d/rs.vcf" | head -6 | sed 's/^/    /'
            ok=0
        fi
    fi
    if [ $ok -eq 1 ]; then
        pass=$((pass+1))
        [ -n "${KEEP:-}" ] || rm -rf "$d"
    else
        fail=$((fail+1)); failed+=("$seed")
    fi
done
echo "passed: $pass, failed: $fail, skipped: $skip ${failed[*]:-}"
[ $fail -eq 0 ]

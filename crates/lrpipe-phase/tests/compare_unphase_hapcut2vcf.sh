#!/usr/bin/env bash
# Compare `whatshap unphase` and `whatshap hapcut2vcf` (Python, the oracle)
# with lrpipe-phase: output (without the ##commandline line) and exit code,
# on every VCF of the WhatsHap test data (unphase, also from standard input)
# and on random inputs from tests/hapcut2vcf/gen.py.
#
# Usage: tests/compare_unphase_hapcut2vcf.sh ENV_BIN WHATSHAP_SRC_DIR WORK_DIR [SEEDS]
set -u
ENV=$1; SRC=$2; WORK=$3; SEEDS=${4:-200}
HERE=$(cd "$(dirname "$0")" && pwd)
RS=${RS:-$HERE/../../../target/release/lrpipe-phase}
pass=0; fail=0; bothfail=0; failed=()
mkdir -p "$WORK"
# run_case NAME STDIN SUBCOMMAND ARGS...
run_case() {
    local name=$1 stdin=$2; shift 2
    local d=$WORK/cases/$name
    rm -rf "$d"; mkdir -p "$d"
    "$ENV/whatshap" "$@" < "$stdin" > "$d/wh.out" 2> "$d/wh.err"; local e1=$?
    "$RS" "$@" < "$stdin" > "$d/rs.out" 2> "$d/rs.err"; local e2=$?
    if [ $e1 -ne 0 ] && [ $e2 -ne 0 ]; then
        echo "$name: $(tail -n1 "$d/wh.err")" >> "$WORK/bothfail.log"
        pass=$((pass+1)); bothfail=$((bothfail+1)); rm -rf "$d"; return
    fi
    if [ $e1 -ne $e2 ]; then
        echo "FAIL $name: exit codes whatshap=$e1 rs=$e2"; tail -n1 "$d/wh.err" "$d/rs.err" | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    if ! cmp -s <(grep -v '^##commandline' "$d/wh.out") <(grep -v '^##commandline' "$d/rs.out"); then
        echo "FAIL $name: output differs"; diff <(grep -v '^##commandline' "$d/wh.out") <(grep -v '^##commandline' "$d/rs.out") | head -6 | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    pass=$((pass+1)); rm -rf "$d"
}
D=$SRC/tests/data
for vcf in "$D"/*.vcf "$D"/*.vcf.gz "$D"/*.bcf; do
    b=$(basename "$vcf")
    run_case "unphase-$b" /dev/null unphase "$vcf"
done
for vcf in "$D"/phased1.vcf "$D"/pacbio/phased.vcf.gz; do
    run_case "unphase-stdin-$(basename "$vcf")" "$vcf" unphase -
done
run_case hapcut2vcf-pacbio /dev/null hapcut2vcf "$D/pacbio/variants.vcf" "$D/pacbio/hapcut.txt"
for seed in $(seq 1 "$SEEDS"); do
    g=$WORK/gen/$seed
    rm -rf "$g"; mkdir -p "$g"
    "$ENV/python" -I "$HERE/hapcut2vcf/gen.py" "$seed" "$g"
    run_case "seed$seed-unphase" /dev/null unphase "$g/in.vcf"
    run_case "seed$seed-hapcut2vcf" /dev/null hapcut2vcf "$g/in.vcf" "$g/hapcut.txt"
    [ -n "${KEEP:-}" ] || rm -rf "$g"
done
echo "unphase/hapcut2vcf: passed: $pass (both tools failed: $bothfail)  failed: $fail ${failed[*]:-}"
[ "$fail" -eq 0 ]

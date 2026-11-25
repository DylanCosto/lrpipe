#!/usr/bin/env bash
# Compare `whatshap compare` (Python, the oracle) with `lrpipe-phase compare`:
# standard output, exit code and the --tsv-pairwise, --tsv-multiway,
# --switch-error-bed and --longest-block-tsv files.
#
# Usage: tests/compare_compare.sh ENV_BIN WHATSHAP_SRC_DIR WORK_DIR [SEEDS [FIRST_SEED]]
#   ENV_BIN: bin directory of an environment with WhatsHap 2.8
#   WHATSHAP_SRC_DIR: a checkout of https://github.com/whatshap/whatshap (v2.8)
#   SEEDS: number of random VCF sets from tests/compare/gen.py (default 200)
#   FIRST_SEED: first random seed (default 1); NO_DATA=1 skips the test data
set -u
ENV=$1; SRC=$2; WORK=$3; SEEDS=${4:-200}; FIRST=${5:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
RS=${RS:-$HERE/../../../target/release/lrpipe-phase}
pass=0; fail=0; bothfail=0; failed=()
mkdir -p "$WORK"
# run_case NAME ARGS...: the words BED, LONGEST and MULTI are replaced by
# per-tool output files; --tsv-pairwise is always written
run_case() {
    local name=$1; shift
    local d=$WORK/cases/$name
    rm -rf "$d"; mkdir -p "$d"
    for t in wh rs; do
        local bin=("$ENV/whatshap"); [ $t = rs ] && bin=("$RS")
        local args=()
        for a in "$@"; do
            case $a in
                BED) args+=("$d/$t.bed") ;;
                LONGEST) args+=("$d/$t.longest") ;;
                MULTI) args+=("$d/$t.multi") ;;
                *) args+=("$a") ;;
            esac
        done
        "${bin[@]}" compare --tsv-pairwise "$d/$t.tsv" "${args[@]}" > "$d/$t.out" 2> "$d/$t.err"
        echo $? > "$d/$t.exit"
    done
    local e1 e2; e1=$(cat "$d/wh.exit"); e2=$(cat "$d/rs.exit")
    if [ "$e1" -ne 0 ] && [ "$e2" -ne 0 ]; then
        echo "$name: $(tail -n1 "$d/wh.err")" >> "$WORK/bothfail.log"
        pass=$((pass+1)); bothfail=$((bothfail+1)); rm -rf "$d"; return
    fi
    if [ "$e1" -ne "$e2" ]; then
        echo "FAIL $name: exit codes whatshap=$e1 rs=$e2"
        tail -n1 "$d/wh.err" "$d/rs.err" | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    for ext in out tsv bed longest multi; do
        if [ -e "$d/wh.$ext" ] || [ -e "$d/rs.$ext" ]; then
            if ! cmp -s "$d/wh.$ext" "$d/rs.$ext"; then
                echo "FAIL $name: $ext differs"; diff "$d/wh.$ext" "$d/rs.$ext" | head -6 | sed 's/^/    /'
                fail=$((fail+1)); failed+=("$name"); return
            fi
        fi
    done
    pass=$((pass+1)); rm -rf "$d"
}

# WhatsHap's test data
D=$SRC/tests/data
if [ -z "${NO_DATA:-}" ]; then
    run_case pair --sample sample1 --switch-error-bed BED --longest-block-tsv LONGEST "$D/phased1.vcf" "$D/phased2.vcf"
    run_case pair-names --sample sample1 --names A,B --only-snvs "$D/phased1.vcf" "$D/phased2.vcf"
    run_case three --sample sample1 --tsv-multiway MULTI --switch-error-bed BED "$D/phased1.vcf" "$D/phased2.vcf" "$D/phased3.vcf"
    run_case three-snvs --sample sample1 --tsv-multiway MULTI --only-snvs "$D/phased1.vcf" "$D/phased2.vcf" "$D/phased3.vcf"
    run_case poly -p 4 --sample sample1 "$D/phased.poly1.vcf" "$D/phased.poly2.vcf"
    run_case poly-snvs -p 4 --sample sample1 --only-snvs "$D/phased.poly1.vcf" "$D/phased.poly2.vcf"
    run_case overlapping "$D/phased_overlapping.vcf" "$D/phased_overlapping.vcf"
    run_case no-sample "$D/phased1.vcf" "$D/phased2.vcf"
    run_case ploidy-mismatch -p 3 --sample sample1 "$D/phased1.vcf" "$D/phased2.vcf"
    for vcf in "$D"/*.vcf; do
        b=$(basename "$vcf")
        sample=$(grep -m1 '^#CHROM' "$vcf" | cut -f10)
        [ -n "$sample" ] && run_case "self-$b" --sample "$sample" "$vcf" "$vcf"
    done
fi

# Random VCF sets
for seed in $(seq "$FIRST" $((FIRST + SEEDS - 1))); do
    g=$WORK/gen/$seed
    rm -rf "$g"; mkdir -p "$g"
    "$ENV/python" -I "$HERE/compare/gen.py" "$seed" "$g"
    vcfs=()
    for i in $(seq 0 $(($(cat "$g/nfiles") - 1))); do vcfs+=("$g/in$i.vcf"); done
    n=0
    while IFS= read -r line; do
        n=$((n+1))
        read -ra extra <<< "$line"
        run_case "seed$seed-$n" "${extra[@]}" "${vcfs[@]}"
    done < "$g/cases.txt"
    [ -n "${KEEP:-}" ] || rm -rf "$g"
done

echo "compare: passed: $pass (both tools failed: $bothfail)  failed: $fail ${failed[*]:-}"
[ "$fail" -eq 0 ]

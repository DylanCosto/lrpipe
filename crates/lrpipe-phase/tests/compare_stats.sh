#!/usr/bin/env bash
# Compare `whatshap stats` (Python, the oracle) with `lrpipe-phase stats`:
# standard output, exit code and the --tsv, --block-list and --gtf files.
#
# Usage: tests/compare_stats.sh ENV_BIN WHATSHAP_SRC_DIR WORK_DIR [SEEDS [FIRST_SEED]]
#   ENV_BIN: bin directory of an environment with WhatsHap 2.8, bgzip and tabix
#   WHATSHAP_SRC_DIR: a checkout of https://github.com/whatshap/whatshap (v2.8)
#   SEEDS: number of random VCFs from tests/stats/gen.py (default 200)
#   FIRST_SEED: first random seed (default 1); NO_DATA=1 skips the test data
# Every VCF of WhatsHap's test data is run with a few option sets; each random
# VCF is run as plain text and, with --chromosome, bgzipped and indexed.
set -u
ENV=$1; SRC=$2; WORK=$3; SEEDS=${4:-200}; FIRST=${5:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
RS=${RS:-$HERE/../../../target/release/lrpipe-phase}
pass=0; fail=0; bothfail=0; failed=()
mkdir -p "$WORK"
run_case() {
    local name=$1 vcf=$2; shift 2
    local d=$WORK/cases/$name
    rm -rf "$d"; mkdir -p "$d"
    for t in wh rs; do
        local bin=("$ENV/whatshap"); [ $t = rs ] && bin=("$RS")
        "${bin[@]}" stats --tsv "$d/$t.tsv" --block-list "$d/$t.blocks" --gtf "$d/$t.gtf" "$@" "$vcf" \
            > "$d/$t.out" 2> "$d/$t.err"
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
    for ext in out tsv blocks gtf; do
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
[ -n "${NO_DATA:-}" ] || for vcf in "$D"/*.vcf "$D"/*.vcf.gz "$D"/*.bcf; do
    [ -e "$vcf" ] || continue
    b=$(basename "$vcf")
    run_case "$b" "$vcf"
    run_case "$b-chrlen" "$vcf" --chr-lengths "$D/chr-lengths.txt" --only-snvs
    first=$("$ENV/bcftools" view -H "$vcf" 2>/dev/null | head -1 | cut -f1)
    [ -n "$first" ] && run_case "$b-chrom" "$vcf" --chromosome "$first"
done

# Random VCFs
for seed in $(seq "$FIRST" $((FIRST + SEEDS - 1))); do
    g=$WORK/gen/$seed
    rm -rf "$g"; mkdir -p "$g"
    "$ENV/python" -I "$HERE/stats/gen.py" "$seed" "$g"
    vcfs=("$g/in.vcf")
    if "$ENV/bgzip" -c "$g/in.vcf" > "$g/in.vcf.gz" 2>/dev/null && "$ENV/tabix" -p vcf "$g/in.vcf.gz" 2>/dev/null; then
        vcfs+=("$g/in.vcf.gz")
    fi
    n=0
    while IFS= read -r line; do
        n=$((n+1))
        read -ra extra <<< "${line//CHRLEN/$g/chr-lengths.txt}"
        for vcf in "${vcfs[@]}"; do
            # the indexed copy only matters for --chromosome
            [ "$vcf" != "$g/in.vcf" ] && [[ "$line" != *--chromosome* ]] && continue
            run_case "seed$seed-$n-$(basename "$vcf")" "$vcf" "${extra[@]}"
        done
    done < "$g/cases.txt"
    [ -n "${KEEP:-}" ] || rm -rf "$g"
done

echo "stats: passed: $pass (both tools failed: $bothfail)  failed: $fail ${failed[*]:-}"
[ "$fail" -eq 0 ]

#!/usr/bin/env bash
# Compare `whatshap split` (Python, the oracle) with `lrpipe-phase split`:
# exit code, the decompressed contents of every output file (BAM, FASTQ
# plain/gz/bz2/xz) and the read lengths histogram.
#
# Usage: tests/compare_split.sh ENV_BIN WHATSHAP_SRC_DIR WORK_DIR [SEEDS [FIRST_SEED]]
#   ENV_BIN: bin directory of an environment with WhatsHap 2.8 and bgzip
#   WHATSHAP_SRC_DIR: a checkout of https://github.com/whatshap/whatshap (v2.8)
#   SEEDS: number of random inputs from tests/split/gen.py (default 200)
#   FIRST_SEED: first random seed (default 1); NO_DATA=1 skips the test data
set -u
ENV=$1; SRC=$2; WORK=$3; SEEDS=${4:-200}; FIRST=${5:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
RS=${RS:-$HERE/../../../target/release/lrpipe-phase}
pass=0; fail=0; bothfail=0; failed=()
mkdir -p "$WORK"
decompress() {
    case $1 in
        *.bam|*.gz) gzip -dc "$1" ;;
        *.bz2) bzip2 -dc "$1" ;;
        *.xz) xz -dc "$1" ;;
        *) cat "$1" ;;
    esac
}
# run_case NAME ARGS...: words NAME:EXT (NAME one of H1..H9, UNTAGGED, HIST)
# are replaced by per-tool output files
run_case() {
    local name=$1; shift
    local d=$WORK/cases/$name
    rm -rf "$d"; mkdir -p "$d"
    local outs=()
    for t in wh rs; do
        local bin=("$ENV/whatshap"); [ $t = rs ] && bin=("$RS")
        local args=()
        for a in "$@"; do
            if [[ $a =~ ^(H[0-9]|UNTAGGED|HIST):(.*)$ ]]; then
                args+=("$d/$t.${BASH_REMATCH[1]}${BASH_REMATCH[2]}")
                [ $t = wh ] && outs+=("${BASH_REMATCH[1]}${BASH_REMATCH[2]}")
            else
                args+=("$a")
            fi
        done
        "${bin[@]}" split "${args[@]}" > "$d/$t.out" 2> "$d/$t.err"
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
    for o in "${outs[@]}"; do
        if ! cmp -s <(decompress "$d/wh.$o") <(decompress "$d/rs.$o"); then
            echo "FAIL $name: $o differs"
            cmp <(decompress "$d/wh.$o") <(decompress "$d/rs.$o") 2>&1 | head -2 | sed 's/^/    /'
            fail=$((fail+1)); failed+=("$name"); return
        fi
    done
    pass=$((pass+1)); rm -rf "$d"
}

# WhatsHap's test data
D=$SRC/tests/data
if [ -z "${NO_DATA:-}" ]; then
    run_case pacbio -o H1:.bam -o H2:.bam --output-untagged UNTAGGED:.bam --read-lengths-histogram HIST:.tsv \
        "$D/pacbio/pacbio.bam" "$D/pacbio/haplotags.txt"
    run_case pacbio-largest --only-largest-block -o H1:.bam -o H2:.bam --add-untagged \
        "$D/pacbio/pacbio.bam" "$D/pacbio/haplotags.txt"
    run_case no-sequence --output-h1 H1:.bam --output-h2 H2:.bam --output-untagged UNTAGGED:.bam \
        --read-lengths-histogram HIST:.tsv "$D/reads-no-sequence.bam" "$D/reads-no-sequence.haplotags.tsv"
    [ -e "$D/pacbio/pacbio.fastq.gz" ] && run_case fastq -o H1:.fastq.gz -o H2:.fastq --discard-unknown-reads \
        "$D/pacbio/pacbio.fastq.gz" "$D/pacbio/haplotags.txt"
fi

# Random inputs
for seed in $(seq "$FIRST" $((FIRST + SEEDS - 1))); do
    g=$WORK/gen/$seed
    rm -rf "$g"; mkdir -p "$g"
    "$ENV/python" -I "$HERE/split/gen.py" "$seed" "$g"
    [ -e "$g/reads.raw" ] && "$ENV/bgzip" -c "$g/reads.raw" > "$g/reads.bam"
    read -r reads list < "$g/files"
    n=0
    while IFS= read -r line; do
        n=$((n+1))
        read -ra extra <<< "$line"
        run_case "seed$seed-$n" "${extra[@]}" "$g/$reads" "$g/$list"
    done < "$g/cases.txt"
    [ -n "${KEEP:-}" ] || rm -rf "$g"
done

echo "split: passed: $pass (both tools failed: $bothfail)  failed: $fail ${failed[*]:-}"
[ "$fail" -eq 0 ]

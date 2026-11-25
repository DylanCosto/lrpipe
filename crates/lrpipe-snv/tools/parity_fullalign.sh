#!/usr/bin/env bash
# Full-alignment tensor parity for every job of a Clair3 run: re-runs each
# `CallVariantsFromCffi --full_aln_regions` command from the run's GNU parallel log under
# tools/oracle_fullalign.py (pipeline mode: Clair3's own process, stopped once the tensors are
# built), then compares snv-fullalign against them bit for bit (tests/parity.rs).
#
# usage: tools/parity_fullalign.sh PARALLEL_LOG OUTDIR [PHASED_VCF] [JOBS]
#   PARALLEL_LOG  log/parallel_6_call_var_bam_full_alignment.log of the Clair3 run
#   PHASED_VCF    replaces --phased_vcf_fn (when the run's tmp/ copy has been cleaned up)
# env: TORCH_ENV (Clair3 env), libtorch settings from env.sh
set -euo pipefail
LOG=${1:?parallel log} OUT=${2:?output dir} PHASED=${3:-} JOBS=${4:-4}
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)

# One line of arguments per job, after `CallVariantsFromCffi`.
tail -n +2 "$LOG" | awk -F'\t' '{print $NF}' | sed "s/'//g; s#^.*clair3.py CallVariantsFromCffi ##" > "$OUT/jobs.txt"
if [ -n "$PHASED" ]; then
    sed -i.bak "s#--phased_vcf_fn [^ ]*#--phased_vcf_fn $PHASED#" "$OUT/jobs.txt" && rm "$OUT/jobs.txt.bak"
fi

: > "$OUT/cases.tsv"
n=0
while read -r line; do
    arg() { sed -n "s#.*--$1 \([^ ]*\).*#\1#p" <<< "$line"; }
    bed=$(arg full_aln_regions)
    name=$(basename "$bed")
    vcf=$(arg phased_vcf_fn)
    if [ "$(arg no_phasing_for_fa)" = "True" ]; then vcf=-; fi
    long=0; if [ "$(arg enable_long_indel)" = "True" ]; then long=1; fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$(arg bam_fn)" "$(arg ref_fn)" "$(arg ctgName)" \
        "$bed" "$vcf" "$(arg minMQ)" "$(arg platform)" "$long" >> "$OUT/cases.tsv"
    printf '%s\t%s\n' "$name" "$line" >> "$OUT/jobs.named.txt.tmp"
    n=$((n + 1))
done < "$OUT/jobs.txt"
mv "$OUT/jobs.named.txt.tmp" "$OUT/jobs.named.txt"

export TORCH_ENV HERE OUT PATH="$TORCH_ENV/bin:$PATH"
oracle() {
    IFS=$'\t' read -r name line <<< "$1"
    # shellcheck disable=SC2086
    "$TORCH_ENV/bin/python3" -I "$HERE/tools/oracle_fullalign.py" pipeline --clair3-dir "$TORCH_ENV/bin" \
        --out "$OUT/$name" -- $line > "$OUT/$name.log" 2>&1 && echo "$name: $(tail -1 "$OUT/$name.log")"
}
export -f oracle
start=$SECONDS
xargs -d '\n' -P "$JOBS" -I{} bash -c 'oracle "$@"' _ {} < "$OUT/jobs.named.txt"
echo "oracle: $n jobs in $((SECONDS - start)) s wall ($JOBS parallel)"

cd "$HERE"
SNV_FA_CASES="$OUT/cases.tsv" cargo test --release -p snv-fullalign --test parity -- --nocapture

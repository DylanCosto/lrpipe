#!/usr/bin/env bash
# Pileup parity for every chunk of one contig: Clair3's Python (oracle) vs lrpipe-snv, both on one
# thread per chunk with Clair3's batch size, then a bit-for-bit comparison per chunk.
#
# usage: tools/parity_contig.sh BAM REF CONTIG CHUNKS OUTDIR [JOBS]
# env:   TORCH_ENV (Clair3 env with PyTorch 2.13.0), MODEL (pileup.pt)
set -euo pipefail

BAM=${1:?bam} REF=${2:?ref} CONTIG=${3:?contig} CHUNKS=${4:?chunk count} OUT=${5:?output dir} JOBS=${6:-8}
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
MODEL=${MODEL:-$TORCH_ENV/bin/models/r1041_e82_400bps_sup_v500/pileup.pt}
HERE=$(cd "$(dirname "$0")/.." && pwd)
PY=$TORCH_ENV/bin/python
mkdir -p "$OUT"

common=(--bam "$BAM" --ref "$REF" --ctg-name "$CONTIG" --chunk-num "$CHUNKS" --extend-bed /nonexistent)
export PY HERE OUT MODEL TORCH_ENV
export COMMON="${common[*]}"

oracle() { "$PY" -I "$HERE/tools/oracle_pileup.py" --clair3-dir "$TORCH_ENV/bin" $COMMON --chunk-id "$1" \
    --model "$MODEL" --out "$OUT/oracle.c$1" > "$OUT/oracle.c$1.log" 2>&1; }
rust() { "$HERE/../../target/release/lrpipe-snv" pileup $COMMON --chunk-id "$1" --model "$MODEL" \
    --out "$OUT/rust.c$1" > "$OUT/rust.c$1.log" 2>&1; }
export -f oracle rust

time_all() {
    local start=$SECONDS
    seq 1 "$CHUNKS" | xargs -P "$JOBS" -I{} bash -c "$1 {}"
    echo "$1: $((SECONDS - start)) s wall for $CHUNKS chunks ($JOBS parallel)"
}
time_all oracle
time_all rust

fail=0
for c in $(seq 1 "$CHUNKS"); do
    if "$PY" -I "$HERE/tools/compare_pileup.py" "$OUT/oracle.c$c" "$OUT/rust.c$c" > "$OUT/compare.c$c.txt"; then
        echo "chunk $c: $(grep -c . "$OUT/oracle.c$c.info.tsv") candidates, identical"
    else
        fail=1
        echo "chunk $c: DIFFERS"; sed 's/^/    /' "$OUT/compare.c$c.txt"
    fi
done
exit $fail

#!/usr/bin/env bash
# Pileup VCF parity for every chunk of one contig: Clair3's own `CallVariantsFromCffi --pileup`
# (exactly as run_clair3 invokes it, one thread per chunk) vs `lrpipe-snv pileup --call-fn`,
# compared byte for byte.
#
# usage: tools/parity_vcf.sh BAM REF CONTIG CHUNKS OUTDIR [JOBS]
# env:   TORCH_ENV (Clair3 env with PyTorch 2.13.0), MODEL_DIR (Clair3 model directory)
set -euo pipefail

BAM=${1:?bam} REF=${2:?ref} CONTIG=${3:?contig} CHUNKS=${4:?chunk count} OUT=${5:?output dir} JOBS=${6:-8}
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
MODEL_DIR=${MODEL_DIR:-$TORCH_ENV/bin/models/r1041_e82_400bps_sup_v500}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"
CMD_FN=$OUT/CMD
echo "lrpipe-snv parity run on $(basename "$BAM") $CONTIG" > "$CMD_FN"
export BAM REF CONTIG CHUNKS OUT TORCH_ENV MODEL_DIR HERE CMD_FN

clair3() {
    PATH=$TORCH_ENV/bin:$PATH "$TORCH_ENV/bin/python3" "$TORCH_ENV/bin/clair3.py" CallVariantsFromCffi \
        --chkpnt_fn "$MODEL_DIR/pileup.pt" --bam_fn "$BAM" --call_fn "$OUT/clair3.c$1.vcf" --sampleName SAMPLE \
        --ref_fn "$REF" --extend_bed "$OUT/no_such_bed" --bed_fn EMPTY --vcf_fn EMPTY --threads 1 \
        --ctgName "$CONTIG" --chunk_id "$1" --chunk_num "$CHUNKS" --platform ont --fast_mode False \
        --snp_min_af 0.08 --indel_min_af 0.15 --minMQ 5 --minCoverage 2 --call_snp_only False \
        --enable_variant_calling_at_sequence_head_and_tail False --gvcf False --base_err 0.001 --gq_bin_size 5 \
        --enable_long_indel False --samtools samtools --temp_file_dir "$OUT" --pileup --keep_iupac_bases False \
        --cmd_fn "$CMD_FN" --use_gpu False --cpu_threads 4 > "$OUT/clair3.c$1.log" 2>&1
}
rust() {
    "$HERE/../../target/release/lrpipe-snv" pileup --bam "$BAM" --ref "$REF" --ctg-name "$CONTIG" --chunk-id "$1" \
        --chunk-num "$CHUNKS" --extend-bed "$OUT/no_such_bed" --model "$MODEL_DIR/pileup.pt" \
        --call-fn "$OUT/rust.c$1.vcf" --cmd-fn "$CMD_FN" --out "$OUT/rust.c$1" > "$OUT/rust.c$1.log" 2>&1
}
export -f clair3 rust

time_all() {
    local start=$SECONDS
    seq 1 "$CHUNKS" | xargs -P "$JOBS" -I{} bash -c "$1 {}"
    echo "$1: $((SECONDS - start)) s wall for $CHUNKS chunks ($JOBS parallel)"
}
time_all clair3
time_all rust

fail=0
for c in $(seq 1 "$CHUNKS"); do
    a=$OUT/clair3.c$c.vcf b=$OUT/rust.c$c.vcf
    if cmp -s "$a" "$b"; then
        echo "chunk $c: $(grep -vc '^#' "$a") records, identical"
    else
        fail=1
        echo "chunk $c: DIFFERS ($(grep -vc '^#' "$a") vs $(grep -vc '^#' "$b") records)"
        diff "$a" "$b" | head -6 | sed 's/^/    /'
    fi
done
exit $fail

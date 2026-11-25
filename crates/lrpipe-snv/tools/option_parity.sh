#!/usr/bin/env bash
# Runs Clair3 (run_clair3.py) and `lrpipe-snv run --read-order input` on one BAM and compares
# every output file both write, ignoring the header lines that name the command. The merged gVCF
# differs where lrpipe-snv fixes Clair3's merge (ci/parity.sh checks it against the listed fixes).
#
# usage: option_parity.sh NAME PLATFORM BAM REF MODEL_DIR CONTIGS [OPTION...]
#   OPTIONs are Clair3's spelling (--gvcf, --print_ref_calls, --snp_min_af=0.1, ...); lrpipe-snv
#   takes the same names with dashes.
# env: WORK (default /data2/bioinformatics/clair3-rs-bench/options), THREADS (8),
#      LRPIPE_SNV (../../target/release/lrpipe-snv), TORCH_ENV, CLAIR3_RUN (run_clair3.py to use, default
#      the env's), SKIP_CLAIR3=1 to reuse an oracle run, ONLY_CLAIR3=1 to only run Clair3.
set -u
NAME=$1 PLATFORM=$2 BAM=$3 REF=$4 MODEL=$5 CONTIGS=$6; shift 6
WORK=${WORK:-/data2/bioinformatics/clair3-rs-bench/options}
THREADS=${THREADS:-8}
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
CLAIR3_RUN=${CLAIR3_RUN:-$TORCH_ENV/bin/run_clair3.py}
HERE=$(cd "$(dirname "$0")/.." && pwd)
LRPIPE_SNV=${LRPIPE_SNV:-$HERE/../../target/release/lrpipe-snv}
A=$WORK/$NAME/clair3 B=$WORK/$NAME/lrpipe-snv
mkdir -p "$WORK/$NAME"

rust_opts=()
for o in "$@"; do
  key=${o%%=*} val=; [[ $o == *=* ]] && val=${o#*=}
  key=--$(echo "${key#--}" | tr _ -)
  if [[ -n $val ]]; then rust_opts+=("$key=$val"); else rust_opts+=("$key"); fi
done

if [[ ${SKIP_CLAIR3:-0} != 1 ]]; then
  rm -rf "$A"
  start=$SECONDS
  PATH=$TORCH_ENV/bin:$PATH nice -n 19 "$TORCH_ENV/bin/python3" "$CLAIR3_RUN" \
    --bam_fn="$BAM" --ref_fn="$REF" --threads="$THREADS" --platform="$PLATFORM" --model_path="$MODEL" \
    --output="$A" --ctg_name="$CONTIGS" "$@" > "$WORK/$NAME/clair3.log" 2>&1 || { echo "clair3 failed: $WORK/$NAME/clair3.log"; exit 1; }
  echo "clair3: $((SECONDS - start)) s"
fi
[[ ${ONLY_CLAIR3:-0} == 1 ]] && exit 0
rm -rf "$B"
start=$SECONDS
nice -n 19 "$LRPIPE_SNV" run --read-order input --bam "$BAM" --ref "$REF" --threads "$THREADS" --platform "$PLATFORM" \
  --model-path "$MODEL" --output "$B" --ctg-name "$CONTIGS" ${rust_opts[@]+"${rust_opts[@]}"} > "$WORK/$NAME/lrpipe-snv.log" 2>&1 \
  || { echo "lrpipe-snv failed: $WORK/$NAME/lrpipe-snv.log"; tail -3 "$WORK/$NAME/lrpipe-snv.log"; exit 1; }
echo "lrpipe-snv: $((SECONDS - start)) s"

norm() { # file, its run dir; Clair3 writes its intermediate gVCFs with lz4
  if [[ $(head -c4 "$1" | od -An -tx1 | tr -d ' ') == 04224d18 ]]; then
    "$TORCH_ENV/bin/lz4" -dc "$1"
  else
    zcat -f "$1"
  fi | grep -avE '^##(cmdline|commandline)|^##source=|^##fileDate' | sed "s#$2/##g"
}
n=0 bad=0
for f in $(cd "$A" && find . -maxdepth 1 -type f \( -name '*.vcf.gz' -o -name '*.gvcf.gz' -o -name '*.bam' \) | sort) \
         $(cd "$A" && find tmp -type f \( -name '*.vcf' -o -name '*.vcf.gz' -o -name '*.gvcf' -o -name '*.bed' \
             -o -path '*candidate_bed/*' -o -name 'CONTIGS' -o -name 'CHUNK_LIST' -o -path '*split_beds/*' \) 2>/dev/null | sort); do
  n=$((n + 1))
  if [[ ! -e $B/$f ]]; then echo "MISSING: $f"; bad=$((bad + 1)); continue; fi
  case $f in
    # lrpipe-snv fixes Clair3's gVCF merge (its README's bugs fixed): differences here are those.
    ./merge_output.gvcf.gz)
      cmp -s <(norm "$A/$f" "$A") <(norm "$B/$f" "$B") || echo "EXPECTED TO DIFFER (Clair3's gVCF merge bugs fixed): $f"
      continue ;;
    *.bam) cmp -s <("$TORCH_ENV/bin/samtools" view "$A/$f") <("$TORCH_ENV/bin/samtools" view "$B/$f") ;;
    *_0.vcf|*gvcf_tmp_output*) cmp -s <(norm "$A/$f" "$A" | sort) <(norm "$B/$f" "$B" | sort) ;;
    *) cmp -s <(norm "$A/$f" "$A") <(norm "$B/$f" "$B") ;;
  esac || { echo "DIFFERS: $f"; bad=$((bad + 1)); }
done
echo "$NAME: compared $n files, $bad differ"

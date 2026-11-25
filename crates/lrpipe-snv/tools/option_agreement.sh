#!/usr/bin/env bash
# Reruns option_parity.sh's cases with lrpipe-snv and compares their final VCFs with Clair3's by
# PASS calls (tools/concordance.py). Each case is run with the options Clair3 recorded in its
# own VCF header (`##cmdline`), which lrpipe-snv takes in Clair3's spelling. It changes the output
# directory and preserves Clair3's file order for ties. A VCF passes when it is the same as
# Clair3's apart from header lines, or its PASS calls agree with Clair3's on at least MIN with
# a Jaccard index of at least MIN.
#
# usage: option_agreement.sh OUTDIR [CASE...]
#   OUTDIR gets one new directory per case (an existing one is an error); CASEs are directories
#   of WORK with a Clair3 run (default: all of them).
# env: WORK (default /data2/bioinformatics/clair3-rs-bench/options), LRPIPE_SNV (the binary, or
#      "lrpipe snv"; default ../../target/release/lrpipe-snv), NVRTC (for cases on the GPU),
#      TORCH_ENV (Clair3's env, for longphase and the like), MIN (default 0.999).
set -u
OUT=$1; shift
WORK=${WORK:-/data2/bioinformatics/clair3-rs-bench/options}
HERE=$(cd "$(dirname "$0")/.." && pwd)
read -r -a LRPIPE_SNV <<< "${LRPIPE_SNV:-$HERE/../../target/release/lrpipe-snv}"
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
MIN=${MIN:-0.999}
mkdir -p "$OUT"
cases=("$@")
if [[ ${#cases[@]} == 0 ]]; then
  for d in "$WORK"/*/clair3; do cases+=("$(basename "$(dirname "$d")")"); done
fi

body() { zcat -f "$1" | grep -v '^##'; }
pass=0 fail=0 failed=()
for name in "${cases[@]}"; do
  a=$WORK/$name/clair3 b=$OUT/$name
  vcfs=$(cd "$a" && ls *.vcf.gz *.gvcf.gz 2>/dev/null)
  if [[ -z $vcfs ]]; then echo "$name: no Clair3 VCF, skipped"; continue; fi
  if [[ -e $b ]]; then echo "$name: $b exists"; exit 1; fi
  header=$(for f in merge_output.vcf.gz pileup.vcf.gz $vcfs; do [[ -e $a/$f ]] && zcat "$a/$f" | grep -m1 '^##cmdline=' && break; done)
  read -r -a words <<< "${header#*=}"
  args=()
  for w in "${words[@]:1}"; do
    case $w in --output=*) args+=("--output=$b") ;; *) args+=("$w") ;; esac
  done
  if [[ " ${args[*]} " == *" --use_gpu"* || " ${args[*]} " == *"--device=cuda"* ]] && [[ -n ${NVRTC:-} ]]; then
    args+=("--nvrtc=$NVRTC")
  fi
  start=$SECONDS
  if ! PATH=$TORCH_ENV/bin:$PATH "${LRPIPE_SNV[@]}" run --read-order input "${args[@]}" > "$b.log" 2>&1; then
    echo "$name: lrpipe-snv failed ($b.log)"; tail -2 "$b.log"
    fail=$((fail + 1)); failed+=("$name"); continue
  fi
  ok=1 report=
  for f in $vcfs; do
    if [[ ! -e $b/$f ]]; then report+=" $f missing;"; ok=0; continue; fi
    if cmp -s <(body "$a/$f") <(body "$b/$f"); then report+=" $f same;"; continue; fi
    # Two VCFs without PASS calls agree.
    read -r agreement jaccard < <(python3 -I "$HERE/tools/concordance.py" "$a/$f" "$b/$f" |
      awk '/^calls:/ { none = $3 == 0 && $5 == 0 } /^agreement/ { print none ? "1.00000 1.00000" : $2 " " $4 }')
    report+=" $f ${agreement:-?} (Jaccard ${jaccard:-?});"
    awk -v a="${agreement:-0}" -v j="${jaccard:-0}" -v m="$MIN" 'BEGIN { exit !(a >= m && j >= m) }' || ok=0
  done
  echo "$name ($((SECONDS - start)) s):$report"
  if [[ $ok == 1 ]]; then pass=$((pass + 1)); else fail=$((fail + 1)); failed+=("$name"); fi
done
echo "passed: $pass  failed: $fail ${failed[*]}"
[[ $fail == 0 ]]

#!/bin/bash
# Runs Clair3's own MergeVcf (PyPy, Clair3 v2.0.3 env) on every case in $WORK/cases.tsv into
# py/<case>/: merge.vcf, the exit status and, for --bed_fn_prefix cases, a copy of the BED
# directory that MergeVcf writes full_aln_regions_<contig> into.
WORK=${WORK:-/data2/bioinformatics/clair3-rs-bench/m4-mergevcf-work}
ENV=${ENV:-/data2/bioinformatics/envs/clair3-t213/bin}
IN=$WORK/inputs
grep -v '^#' "$WORK/cases.tsv" | while IFS=$'\t' read -r name pileup fa extra; do
  out=$WORK/py/$name; rm -rf "$out"; mkdir -p "$out"
  eval "flags=($extra)"
  args=()
  for ((i = 0; i < ${#flags[@]}; i += 2)); do
    flag=${flags[i]}; value=${flags[i+1]}
    case $flag in
      --par_regions_bed) value=$IN/$value ;;
      --bed_fn_prefix) [ -d "$IN/$value" ] && cp -r "$IN/$value" "$out/$value"; value=$out/$value ;;
    esac
    args+=("$flag" "$value")
  done
  ( cd "$out" && PATH=$ENV:$PATH pypy3 $ENV/clair3.py MergeVcf --pileup_vcf_fn "$IN/$pileup" \
      --bed_fn_prefix "$out/no_beds" --full_alignment_vcf_fn "$IN/$fa" --output_fn "$out/merge.vcf" \
      --platform ont --print_ref_calls False --gvcf False --haploid_precise False --haploid_sensitive False \
      --gender unknown --par_regions_bed EMPTY --qual 2 "${args[@]}" ) > "$out/stdout" 2> "$out/stderr"
  echo $? > "$out/status"
  echo "$name status=$(cat "$out/status") files=$(ls "$out" | tr '\n' ' ')"
done

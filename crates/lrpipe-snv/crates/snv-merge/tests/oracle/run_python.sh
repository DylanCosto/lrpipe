#!/bin/bash
# Runs Clair3's own SortVcf (PyPy, Clair3 v2.0.3 env) on every case in cases.tsv into py/<case>/.
WORK=${WORK:-/data2/bioinformatics/clair3-rs-bench/m3-merge-work}
ENV=${ENV:-/data2/bioinformatics/envs/clair3-t213/bin}
REF=$(cat $WORK/inputs/REF)
grep -v '^#' $WORK/cases.tsv | while IFS=$'\t' read -r name input contigs cmd ref extra; do
  out=$WORK/py/$name; rm -rf "$out"; mkdir -p "$out"
  [ "$cmd" != None ] && cmd=$WORK/inputs/$cmd
  refarg=(); [ "$ref" = REF ] && refarg=(--ref_fn $REF)
  extra=${extra//par.bed/$WORK/inputs/par.bed}; extra=${extra//missing.bed/$WORK/inputs/missing.bed}
  ( cd $out && PATH=$ENV:$PATH pypy3 $ENV/clair3.py SortVcf --input_dir $WORK/inputs/$input --vcf_fn_prefix pileup \
      --output_fn $out/pileup.vcf --sampleName SAMPLE --cmd_fn "$cmd" "${refarg[@]}" --pileup_only False \
      --print_ref_calls False --haploid_precise False --haploid_sensitive False --gender unknown \
      --par_regions_bed EMPTY --qual 2 --use_gpu False --contigs_fn $WORK/inputs/$contigs $extra ) > $out/stdout 2> $out/stderr
  echo $? > $out/status
  echo "$name status=$(cat $out/status) files=$(ls $out | tr '\n' ' ')"
done

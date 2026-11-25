#!/usr/bin/env bash
# Runs lrpipe-cnv with each change to Spectre's method on the benchmark inputs and scores the
# calls (crates/lrpipe-cnv-bench). Usage: tools/bench_variants.sh [variant ...]
# Needs SPECTRE_BENCH (inputs/, as tools/make_reference.sh) and SPECTRE_REFERENCE; the GIAB
# HG002 truth defaults to $SPECTRE_BENCH/inputs/truth_HG002_sv.vcf.gz.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
B=${SPECTRE_BENCH:?} REF=${SPECTRE_REFERENCE:?}
TRUTH=${TRUTH:-$B/inputs/truth_HG002_sv.vcf.gz}
rs=$here/../../target/release/lrpipe-cnv bench=$here/../../target/release/lrpipe-cnv-bench
declare -A V=(
  [spectre]=""
  [sex]="--sex detect"
  [gc]="--gc-correction"
  [amp]="--keep-amplifications"
  [cn]="--refresh-copy-numbers"
  [mapq]="--all-mapq-coverage {Q0}"
  [improved]="--improved"
  [improved50k]="--improved --min-cnv-len 50000"
  [spectre50k]="--min-cnv-len 50000"
  [improved_mapq]="--improved --all-mapq-coverage {Q0}"
  [hmm]="--improved --caller hmm"
  [hmm_snv]="--improved --caller hmm -v {SNV}"
  [improved_sv]="--improved --sv-breakpoints {SV}"
  [hmm_sv]="--improved --caller hmm --sv-breakpoints {SV}"
)
variants=("$@"); [[ ${#variants[@]} -gt 0 ]] || variants=(spectre sex gc amp cn mapq improved improved_mapq hmm hmm_snv improved_sv hmm_sv)
# name, coverage directory, sample, truth, lrpipe-cnv-bench options, SNV VCF, SV VCF ("" for none;
# variants needing a missing one are skipped). Inputs: {dir}/mosdepth.regions.bed.gz, and
# {dir}/q0.regions.bed.gz for {Q0}.
SV_HIFI=${SV_HIFI:-$B/inputs/hg002_hifi_wgs/truffle-sv.vcf.gz}
SV_ONT=${SV_ONT:-$B/inputs/hg002_ont_wgs/truffle-sv.vcf.gz}
inputs=(
  "hifi|inputs/hg002_hifi_wgs|HG002|$TRUTH||$B/inputs/hg002_hifi_wgs/pellucid.vcf.gz|$SV_HIFI"
  "ont|inputs/hg002_ont_wgs|HG002|$TRUTH|||$SV_ONT"
  "synth|inputs/hg003_chr20_synth_only|HG003|$here/tools/synth_truth.bed|--contigs chr20||"
)
for v in "${variants[@]}"; do
  for spec in "${inputs[@]}"; do
    IFS='|' read -r name dir sample truth bopts snv sv <<< "$spec"
    q0=$B/$dir/q0.regions.bed.gz
    opts=${V[$v]}
    [[ $opts == *"{Q0}"* && ! -e $q0 ]] && continue
    [[ $opts == *"{SNV}"* && -z $snv ]] && continue
    [[ $opts == *"{SV}"* && ! -e $sv ]] && continue
    opts=${opts//\{Q0\}/$q0}; opts=${opts//\{SNV\}/$snv}; opts=${opts//\{SV\}/$sv}
    only=(); [[ $name == synth ]] && only=(-o chr20)
    out=$B/rs/variants/$v/$name; mkdir -p "$out"
    # shellcheck disable=SC2086
    "$rs" CNVCaller -c "$B/$dir/mosdepth.regions.bed.gz" -s "$sample" -d "$out" -r "$REF" \
      -m "$B/inputs/grch38_reference.mdr" -l "$B/inputs/grch38_blacklist_spectre.bed" "${only[@]}" $opts \
      2> "$out/log.txt"
    # shellcheck disable=SC2086
    "$bench" "$out/$sample.vcf.gz" "$truth" $bopts | awk -v v="$v" -v n="$name" '$2=="all"{print v"\t"n"\t"$0}'
  done
done

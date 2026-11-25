#!/usr/bin/env bash
# After run_all.sh: each lrpipe output against the original tool's, and the accuracy scores.
# Output-identity checks ignore only lines naming the program, its command or the date.
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/paths.sh"
V=$VALIDATION_DIR
O=$V/out
G=$GIAB_DIR
SNV=$REPO/crates/lrpipe-snv/tools/concordance.py
compare_dir=$(mktemp -d "$V/compare.XXXXXX")
trap 'rm -f -- "$compare_dir/left" "$compare_dir/right"; rmdir -- "$compare_dir"' EXIT
ok() {
  local rc
  if "$@"; then echo "  identical"; else
    rc=$?
    if [[ $rc -eq 1 ]]; then echo "  DIFFERENT"; else echo "  comparison failed" >&2; exit "$rc"; fi
  fi
}
body() {
  case $1 in *.gz) zcat "$1" ;; *) cat "$1" ;; esac | awk '!/^##(source|command|commandline|cmdline|fileDate|.*_?command)=/'
}
compare_body() {
  body "$1" > "$compare_dir/left" && body "$2" > "$compare_dir/right" || return 2
  cmp -s "$compare_dir/left" "$compare_dir/right"
}
records() {
  case $1 in *.gz) zcat "$1" ;; *) cat "$1" ;; esac | awk '!/^#/'
}
compare_records() {
  records "$1" > "$compare_dir/left" && records "$2" > "$compare_dir/right" || return 2
  cmp -s "$compare_dir/left" "$compare_dir/right"
}

echo "== lrpipe sv against Sniffles 2.8.1 (VCF records and header, but program, command and date)"
for pair in "sv-hg003-ont.vcf $SNIFFLES_HG003_ONT_VCF" \
            "sv-hg002-cram.vcf $SNIFFLES_HG002_CRAM_VCF" \
            "sv-hg002-hifi.vcf $SNIFFLES_HG002_HIFI_VCF"; do
  set -- $pair; echo " $1: $(grep -vc '^#' $O/$1) records"
  ok compare_body "$O/$1" "$2"
done

echo "== lrpipe depth against mosdepth 0.3.14 (every output file)"
for f in mosdepth.summary.txt mosdepth.global.dist.txt mosdepth.region.dist.txt regions.bed.gz; do
  echo " $f"; ok compare_body "$O/depth-lrpipe.$f" "$O/depth-mosdepth.$f"
done

echo "== lrpipe cnv against Spectre 0.2.1 (VCF and BED)"
for f in HG002.vcf.gz HG002_cnv.bed.gz; do
  echo " $f"; ok compare_body "$O/cnv-lrpipe/$f" "$O/cnv-spectre/$f"
done

echo "== lrpipe methyl pileup against modkit 0.6.4 (bedMethyl)"
ok cmp -s $O/methyl-lrpipe.bed $O/methyl-modkit.bed
echo "== lrpipe methyl cpg against pb-CpG-tools 3.0.0, count mode"
for t in combined hap1 hap2; do
  echo " $t (records)"; ok compare_records "$O/cpg-lrpipe.$t.bed.gz" "$O/cpg-pb.$t.bed.gz"
done

echo "== lrpipe phase against WhatsHap 2.8 (per contig, records without ##commandline)"
same=0 diff=0
for f in $O/phase/whatshap/*.vcf; do
  c=$(basename $f .vcf)
  records "$f" > "$compare_dir/left"
  records "$O/phase/lrpipe.vcf" | awk -F'\t' -v c="$c" '$1==c' > "$compare_dir/right"
  if cmp -s "$compare_dir/left" "$compare_dir/right"; then same=$((same+1)); else diff=$((diff+1)); echo "  $c differs"; fi
done
echo "  $same contigs identical, $diff differ"
echo " lrpipe phase against lrpipe-snv's own final phasing (phased_merge_output)"
ok compare_records "$O/phase/lrpipe.vcf" "$HG002_HIFI_RUN/snv/phased_merge_output.vcf.gz"

echo "== lrpipe tr and TRGT 5.1.0 against GIAB HG002 TR v1.0.1 Tier 1"
TRB=$LRPIPE_TR_BENCH
TRT=$GIAB_TR_TIER1
for t in tr-lrpipe.vcf.gz tr-trgt.vcf.gz; do
  echo " $t, every Tier 1 region"; $TRB $O/$t $TRT
done
echo " on the regions both score"
echo "  lrpipe"; $TRB --common $O/tr-trgt.vcf.gz $O/tr-lrpipe.vcf.gz $TRT
echo "  TRGT"; $TRB --common $O/tr-lrpipe.vcf.gz $O/tr-trgt.vcf.gz $TRT

echo "== lrpipe snv against Clair3 2.0.3 (PASS calls: position, alleles, genotype)"
for run in snv-ont-cpu-exact snv-ont-cpu-approx; do
  echo " $run against Clair3 CPU"; python3 -I $SNV $CLAIR3_HG003_ONT_CPU $O/$run/merge_output.vcf.gz | tail -n 1
done
echo " snv-hifi-gpu against Clair3 GPU"; python3 -I $SNV $CLAIR3_HG002_HIFI_GPU $O/snv-hifi-gpu/merge_output.vcf.gz | tail -n 1
echo " exact against approximate, CPU"
python3 -I $SNV $O/snv-ont-cpu-exact/merge_output.vcf.gz $O/snv-ont-cpu-approx/merge_output.vcf.gz | tail -n 1
python3 -I $SNV $O/snv-hifi-cpu-exact/merge_output.vcf.gz $O/snv-hifi-cpu-approx/merge_output.vcf.gz | tail -n 1

echo "== lrpipe snv and Clair3 against GIAB v4.2.1 (rtg vcfeval, stratified)"
for run in snv-ont-cpu-exact snv-ont-cpu-approx snv-ont-gpu; do echo " $run"; "$(dirname "$0")"/score2.sh HG003 $run $O/$run/merge_output.vcf.gz; done
echo " clair3-ont-cpu"; "$(dirname "$0")"/score2.sh HG003 clair3-ont-cpu $CLAIR3_HG003_ONT_CPU
for run in snv-hifi-cpu-exact snv-hifi-cpu-approx snv-hifi-gpu; do echo " $run"; "$(dirname "$0")"/score2.sh HG002 $run $O/$run/merge_output.vcf.gz; done
echo " clair3-hifi-gpu"; "$(dirname "$0")"/score2.sh HG002 clair3-hifi-gpu $CLAIR3_HG002_HIFI_GPU

echo "== phasing against GIAB HG002 Q100 v5.0q"
"$(dirname "$0")"/score_phasing.sh hifi-pipeline $HG002_HIFI_RUN/snv/phased_merge_output.vcf.gz \
  ont-pipeline $HG002_ONT_RUN/snv/phased_merge_output.vcf.gz hifi-phase $O/phase/lrpipe.vcf

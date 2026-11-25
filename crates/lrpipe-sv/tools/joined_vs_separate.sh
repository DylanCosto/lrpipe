#!/usr/bin/env bash
# One pass (SVs, depth, tandem repeats, CNVs) against the four tools run separately on the same
# reads: times each and checks every output is the same.
#   tools/joined_vs_separate.sh OUT_DIR READS REFERENCE CATALOG MDR BLACKLIST [KARYOTYPE] [PLATFORM]
# Needs lrpipe-sv (this repo), lrpipe-depth, lrpipe-tr and lrpipe-cnv release builds (paths below).
set -uo pipefail
OUT=$1 READS=$2 REF=$3 CAT=$4 MDR=$5 BL=$6 KARYO=${7:-xx} PLAT=${8:-hifi}
bin=$(cd "$(dirname "$0")/../../../target/release" && pwd)
RS=${RS:-$bin/lrpipe-sv}
DB=${DB:-$bin/lrpipe-depth}
TR=${TR:-$bin/lrpipe-tr}
SP=${SP:-$bin/lrpipe-cnv}
T=16
mkdir -p "$OUT/sep" "$OUT/joint"
t() { local n=$1; shift; /usr/bin/time -f "%e %M" -o "$OUT/$n.time" "$@" > "$OUT/$n.log" 2>&1 || { echo "$n FAILED: $(tail -3 "$OUT/$n.log")"; exit 1; }; echo "$n: $(awk '{printf "%.1f s, %.1f GB", $1, $2/1048576}' "$OUT/$n.time")"; }
cd "$OUT"
t sep_sv "$RS" --input "$READS" --vcf sep/sv.vcf --reference "$REF" --threads $T --allow-overwrite
t sep_depth "$DB" -t $T -n -x -Q 20 -b 1000 -f "$REF" sep/depth "$READS"
t sep_repeats "$TR" --reads "$READS" --genome "$REF" --repeats "$CAT" --output sep/repeats.vcf.gz --karyotype "$KARYO" --platform "$PLAT" -t $T
t sep_cnv "$SP" CNVCaller -c sep/depth.regions.bed.gz -s SAMPLE -d sep/cnv -r "$REF" -m "$MDR" -l "$BL" \
    --improved --sv-breakpoints sep/sv.vcf --vcf-style standard
t joint "$RS" --input "$READS" --vcf joint/sv.vcf --reference "$REF" --threads $T --allow-overwrite \
    --depth joint/depth --depth-by 1000 --depth-fast --depth-mapq 20 \
    --repeats "$CAT" --repeats-vcf joint/repeats.vcf.gz --karyotype "$KARYO" --repeats-platform "$PLAT" \
    --cnv joint/cnv --cnv-metadata "$MDR" --cnv-blacklist "$BL"
sep=$(cat sep_sv.time sep_depth.time sep_repeats.time sep_cnv.time | awk '{s+=$1} END {printf "%.1f", s}')
echo "separate total: $sep s; one pass: $(awk '{printf "%.1f", $1}' joint.time) s"
strip() { grep -v -e '^##source=' -e '^##command=' -e '^##fileDate=' "$1"; }
same() { cmp -s <(eval "$2") <(eval "$3") && echo "same: $1" || echo "DIFFERENT: $1"; }
same "SV VCF" "strip sep/sv.vcf" "strip joint/sv.vcf"
same "depth windows" "gzip -dc sep/depth.regions.bed.gz" "gzip -dc joint/depth.regions.bed.gz"
same "depth summary" "cat sep/depth.mosdepth.summary.txt" "cat joint/depth.mosdepth.summary.txt"
same "repeat VCF" "gzip -dc sep/repeats.vcf.gz" "gzip -dc joint/repeats.vcf.gz"
same "CNV VCF" "gzip -dc sep/cnv/SAMPLE.vcf.gz | grep -v '^##fileDate'" "gzip -dc joint/cnv/SAMPLE.vcf.gz | grep -v '^##fileDate'"
same "CNV BED" "gzip -dc sep/cnv/SAMPLE_cnv.bed.gz" "gzip -dc joint/cnv/SAMPLE_cnv.bed.gz"

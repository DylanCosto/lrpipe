#!/usr/bin/env bash
# Compares lrpipe-sv multi-sample calling with a Sniffles2 2.8.1 combine output on the same
# .snf inputs. Records are compared sorted: above 20 samples Sniffles2 writes calls that arrive out
# of position order to a separate result-*-unsorted.part.vcf in its temp directory instead of the
# output, so pass that file as UNSORTED to add those calls back.
#
# usage: tools/combine_parity.sh OUTDIR ORACLE.vcf [UNSORTED.part.vcf] -- INPUT.tsv|INPUT.snf... [options]
#   LRPIPE_SV=/path/to/lrpipe-sv  (default: target/release/lrpipe-sv)
set -euo pipefail

out=${1:?output directory}
oracle=${2:?Sniffles2 combine VCF}
shift 2
unsorted=/dev/null
if [[ ${1:-} != -- ]]; then
    unsorted=$1
    shift
fi
[[ ${1:-} == -- ]] && shift
rs=${LRPIPE_SV:-$(dirname "$0")/../../../target/release/lrpipe-sv}

mkdir -p "$out"
/usr/bin/time -f "%e s %M KB" -o "$out/rs.time" \
    "$rs" --input "$@" --vcf "$out/rs.vcf" --allow-overwrite > "$out/rs.log" 2>&1

headers() { grep '^#' "$1" | grep -v -e '^##source=' -e '^##command=' -e '^##fileDate='; }
records() { grep -hv '^#' "$@" | sort; }
if diff <(headers "$oracle") <(headers "$out/rs.vcf") > "$out/diff.txt" &&
    diff <(records "$oracle" "$unsorted") <(records "$out/rs.vcf") >> "$out/diff.txt"; then
    echo "IDENTICAL  $(grep -vc '^#' "$out/rs.vcf") records  rs $(cat "$out/rs.time")"
else
    echo "DIFFERENT  $(grep -c '^[<>]' "$out/diff.txt") differing lines (see $out/diff.txt)"
    exit 1
fi

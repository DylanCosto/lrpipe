#!/usr/bin/env bash
# Runs Sniffles2 (the oracle) and lrpipe-sv with the same options and compares the VCFs.
# Header lines that name the program, command or date are ignored.
#
# usage: tools/parity.sh OUTDIR INPUT.bam [extra options for both tools...]
#   SNIFFLES=/path/to/sniffles  (default: the vareval env's Sniffles 2.8.1)
#   LRPIPE_SV=/path/to/lrpipe-sv  (default: target/release/lrpipe-sv)
#   KEEP_ORACLE=1  reuse OUTDIR/oracle.vcf if present
set -euo pipefail

out=${1:?output directory}
bam=${2:?input BAM}
shift 2
sniffles=${SNIFFLES:-/data/bioinformatics/micromamba/envs/vareval/bin/sniffles}
rs=${LRPIPE_SV:-$(dirname "$0")/../../../target/release/lrpipe-sv}

mkdir -p "$out"
if [[ -z ${KEEP_ORACLE:-} || ! -s $out/oracle.vcf ]]; then
    /usr/bin/time -f "%e s %M KB" -o "$out/oracle.time" \
        "$sniffles" --input "$bam" --vcf "$out/oracle.vcf" --allow-overwrite "$@" > "$out/oracle.log" 2>&1
fi
/usr/bin/time -f "%e s %M KB" -o "$out/rs.time" \
    "$rs" --input "$bam" --vcf "$out/rs.vcf" --allow-overwrite "$@" > "$out/rs.log" 2>&1

strip() { grep -v -e '^##source=' -e '^##command=' -e '^##fileDate=' -e '^##genotypeSource=' -e '^##genotypeCommand=' -e '^##genotypeFileDate=' "$1"; }
if diff <(strip "$out/oracle.vcf") <(strip "$out/rs.vcf") > "$out/diff.txt"; then
    echo "IDENTICAL  $(grep -vc '^#' "$out/rs.vcf") records  oracle $(cat "$out/oracle.time")  rs $(cat "$out/rs.time")"
else
    echo "DIFFERENT  $(grep -c '^[<>]' "$out/diff.txt") differing lines (see $out/diff.txt)"
    exit 1
fi

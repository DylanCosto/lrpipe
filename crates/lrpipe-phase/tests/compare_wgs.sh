#!/usr/bin/env bash
# Compare WhatsHap 2.8 and lrpipe-phase on real data (e.g. a whole genome):
# phase, then haplotag with WhatsHap's phased VCF, timing every run.
#
# Usage: tests/compare_wgs.sh ENV_BIN_DIR REFERENCE VCF BAM WORK_DIR [THREADS...]
#   ENV_BIN_DIR: bin/ of an environment with whatshap 2.8, bgzip, tabix and python3
#   THREADS:     lrpipe-phase thread counts to run (default: 8)
# Environment:
#   RS          lrpipe-phase binary (default: target/release/lrpipe-phase)
#   EXTRA_ARGS  arguments for both tools and commands (default: --ignore-read-groups)
#   SKIP_HAPLOTAG=1  only compare phase
#   SKIP_PHASE=1     only compare haplotag, using WhatsHap's phased VCF from an
#                    earlier run in WORK_DIR (VCF is not used); results are
#                    appended to that run's compare.txt and timings.txt
#
# WhatsHap's haplotag output is hashed while it is written (it never touches
# the disk); each lrpipe-phase BAM is hashed and then deleted, so WORK_DIR needs
# room for one output BAM. Results go to WORK_DIR/compare.txt and timings to
# WORK_DIR/timings.txt; the exit status is 1 if any output differs.
set -u
if [ $# -lt 5 ]; then
    sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi
E=$1; REF=$2; VCF=$3; BAM=$4; WORK=$5; shift 5
THREADS=("${@:-8}")
TESTS=$(cd "$(dirname "$0")" && pwd)
RS=${RS:-$TESTS/../../../target/release/lrpipe-phase}
read -r -a ARGS <<< "${EXTRA_ARGS---ignore-read-groups}"
mkdir -p "$WORK"
cd "$WORK" || exit 2
T="%e s wall, %U s user, %P CPU, %M KB max RSS"
timed() { local name=$1; shift; /usr/bin/time -f "$name: $T" -a -o timings.txt "$@"; }
same() { cmp -s "$1" "$2" && echo identical || echo DIFFERENT; }
step() { echo "$(date +%T) $*" | tee -a progress.txt; }
record_md5() { "$E/bgzip" -dc -@4 | "$E/python3" -I "$TESTS/record_md5.py"; }
if [ "${SKIP_PHASE:-0}" = 1 ]; then
    [ -s wh.phased.vcf ] || { echo "SKIP_PHASE=1 needs $WORK/wh.phased.vcf from an earlier run"; exit 2; }
else
    : > timings.txt; : > compare.txt; : > progress.txt
    "$RS" --version >> compare.txt 2>&1
    "$E/whatshap" --version 2>&1 | sed 's/^/whatshap /' >> compare.txt

    step "phase: whatshap"
    timed "phase whatshap" "$E/whatshap" phase "${ARGS[@]}" -r "$REF" -o wh.phased.vcf --output-read-list wh.reads "$VCF" "$BAM" 2> wh.phase.log \
        || { echo "whatshap phase failed (see $WORK/wh.phase.log)"; exit 2; }
    grep -v '^##commandline' wh.phased.vcf > wh.phased.body
    for t in "${THREADS[@]}"; do
        step "phase: lrpipe-phase -t $t"
        timed "phase lrpipe-phase -t $t" "$RS" phase -t "$t" "${ARGS[@]}" -r "$REF" -o rs$t.phased.vcf --output-read-list rs$t.reads "$VCF" "$BAM" 2> rs$t.phase.log
        grep -v '^##commandline' rs$t.phased.vcf > rs$t.phased.body
        echo "phase -t $t: vcf $(same wh.phased.body rs$t.phased.body), read list $(same wh.reads rs$t.reads)" | tee -a compare.txt
    done
fi

if [ "${SKIP_HAPLOTAG:-0}" != 1 ]; then
    "$E/bgzip" -f -c wh.phased.vcf > phased.vcf.gz && "$E/tabix" -f -p vcf phased.vcf.gz
    step "haplotag: whatshap"
    timed "haplotag whatshap" "$E/whatshap" haplotag "${ARGS[@]}" -r "$REF" -o /dev/stdout --output-haplotag-list wh.ht.list phased.vcf.gz "$BAM" 2> wh.ht.log \
        | record_md5 > wh.ht.md5
    for t in "${THREADS[@]}"; do
        step "haplotag: lrpipe-phase -t $t"
        timed "haplotag lrpipe-phase -t $t" "$RS" haplotag -t "$t" "${ARGS[@]}" -r "$REF" -o rs$t.ht.bam --output-haplotag-list rs$t.ht.list phased.vcf.gz "$BAM" 2> rs$t.ht.log
        record_md5 < rs$t.ht.bam > rs$t.ht.md5
        rm -f rs$t.ht.bam
        echo "haplotag -t $t: records $(same wh.ht.md5 rs$t.ht.md5) ($(cut -d' ' -f1 rs$t.ht.md5)), list $(same wh.ht.list rs$t.ht.list)" | tee -a compare.txt
    done
fi
step done
cat timings.txt
! grep -q DIFFERENT compare.txt

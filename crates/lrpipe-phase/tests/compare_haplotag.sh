#!/usr/bin/env bash
# Compare `whatshap haplotag` (Python, the oracle) with `lrpipe-phase haplotag`.
#
# Usage: tests/compare_haplotag.sh ENV_BIN_DIR WHATSHAP_SRC_DIR WORK_DIR
#   ENV_BIN_DIR: bin/ of an environment with whatshap 2.8, samtools, python+pysam
# Compares the output BAM records (incl. BAM tag types), the header except for
# the tool's own @PG line, and the haplotag list.
set -u
E=$1; SRC=$2; WORK=$3
RS=${RS:-$(dirname "$0")/../../../target/release/lrpipe-phase}
DUMP=$(dirname "$0")/dump_bam.py
D=$SRC/tests/data
mkdir -p "$WORK"
pass=0; fail=0; failed=()
run_case() {
    local name=$1; shift
    local args=("$@")
    local o_wh=$WORK/$name.wh.bam o_rs=$WORK/$name.rs.bam l_wh=$WORK/$name.wh.list l_rs=$WORK/$name.rs.list
    rm -f "$o_wh" "$o_rs" "$l_wh" "$l_rs"
    "$E/whatshap" haplotag -o "$o_wh" --output-haplotag-list "$l_wh" "${args[@]}" > "$WORK/$name.wh.log" 2>&1; local e1=$?
    "$RS" haplotag ${RS_EXTRA:-} -o "$o_rs" --output-haplotag-list "$l_rs" "${args[@]}" > "$WORK/$name.rs.log" 2>&1; local e2=$?
    if [ $e1 -ne 0 ] || [ $e2 -ne 0 ]; then
        if [ $e1 -ne 0 ] && [ $e2 -ne 0 ]; then
            echo "BOTH-FAIL $name (whatshap: $(tail -1 "$WORK/$name.wh.log" | cut -c1-110)) (rs: $(tail -1 "$WORK/$name.rs.log" | cut -c1-110))"
            pass=$((pass+1)); return
        fi
        echo "FAIL $name: exit codes whatshap=$e1 rs=$e2"; tail -2 "$WORK/$name.rs.log" "$WORK/$name.wh.log" | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    "$E/python" "$DUMP" "$o_wh" > "$WORK/$name.wh.dump"; "$E/python" "$DUMP" "$o_rs" > "$WORK/$name.rs.dump"
    if ! cmp -s "$WORK/$name.wh.dump" "$WORK/$name.rs.dump"; then
        echo "FAIL $name: BAM differs"; diff "$WORK/$name.wh.dump" "$WORK/$name.rs.dump" | head -6 | cut -c1-250 | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    if ! cmp -s <(gzip -cdf "$l_wh") <(gzip -cdf "$l_rs"); then
        echo "FAIL $name: haplotag list differs"; diff <(gzip -cdf "$l_wh") <(gzip -cdf "$l_rs") | head -6 | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    echo "ok   $name ($(grep -c 'HP:C' "$WORK/$name.rs.dump") tagged records)"
    pass=$((pass+1))
}
NR=--no-reference
run_case basic1            $NR $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case basic2            $NR $D/haplotag_2.vcf.gz $D/haplotag.bam
run_case basic-bcf         $NR $D/haplotag_2.bcf $D/haplotag.bam
run_case basic-csi         $NR $D/haplotag_with_csi_index.vcf.gz $D/haplotag.bam
run_case ignore-rg         $NR --ignore-read-groups $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case noSM              $NR --ignore-read-groups $D/haplotag_1.vcf.gz $D/haplotag_noSM.bam
run_case noSM-norg         $NR $D/haplotag_1.vcf.gz $D/haplotag_noSM.bam
run_case missing-chr       $NR $D/haplotag.missing_chr.vcf.gz $D/haplotag.large.bam
run_case without-chr2      $NR $D/haplotag.without_chr2.vcf.gz $D/haplotag.large.bam
run_case without-chr2-skip $NR --skip-missing-contigs $D/haplotag.without_chr2.vcf.gz $D/haplotag.large.bam
run_case noRG              $NR $D/haplotag_1.vcf.gz $D/haplotag_noRG.bam
run_case noRG-ignore       $NR --ignore-read-groups $D/haplotag_1.vcf.gz $D/haplotag_noRG.bam
run_case sample-mother     $NR --sample mother $D/haplotag_sample.vcf.gz $D/haplotag_sample.bam
run_case sample-all        $NR $D/haplotag_sample.vcf.gz $D/haplotag_sample.bam
run_case sample-missing    $NR --sample nobody $D/haplotag_sample.vcf.gz $D/haplotag_sample.bam
run_case large             $NR $D/haplotag.large.vcf.gz $D/haplotag.large.bam
run_case large2            $NR $D/haplotag.large.2.vcf.gz $D/haplotag.large.bam
run_case 10x               $NR $D/haplotag.10X.vcf.gz $D/haplotag.10X.bam
run_case 10x-2             $NR $D/haplotag.10X_2.vcf.gz $D/haplotag.10X.bam
run_case 10x-3             $NR $D/haplotag.10X.vcf.gz $D/haplotag.10X_3.bam
run_case 10x-3-ignore      $NR --ignore-linked-read $D/haplotag.10X.vcf.gz $D/haplotag.10X_3.bam
run_case 10x-cutoff        $NR -d 100 $D/haplotag.10X.vcf.gz $D/haplotag.10X.bam
run_case supp              $NR --ignore-read-groups $D/haplotag.supplementary.vcf.gz $D/haplotag.supplementary.bam
run_case supp-copy         $NR --tag-supplementary --ignore-read-groups $D/haplotag.supplementary.vcf.gz $D/haplotag.supplementary.bam
SS=$D/supplementary_strategy_test.grch38
for strat in skip copy-primary independent-or-skip independent-or-copy-primary; do
    run_case strat-$strat  $NR --ignore-read-groups --tag-supplementary $strat $SS.vcf.gz $SS.bam
    run_case strat-$strat-nostrand $NR --ignore-read-groups --tag-supplementary $strat --no-supplementary-strand-match $SS.vcf.gz $SS.bam
    run_case strat-$strat-dist0 $NR --ignore-read-groups --tag-supplementary $strat --supplementary-distance 0 $SS.vcf.gz $SS.bam
done
run_case region-chr1       $NR --regions chr1 $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case region-chr2       $NR --regions chr2 $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case region-multi      $NR --regions chr1:1000-50000 --regions chr1:60000 $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case region-overlap    $NR --regions chr1:1000-80000 --regions chr1:60000 $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case region-reversed   $NR --regions chr1:60000 --regions chr1:1000-50000 $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case region-large      $NR --regions chr2:10000-30000000 --regions chr1 $D/haplotag.large.vcf.gz $D/haplotag.large.bam
run_case pacbio-ref        -r $D/pacbio/reference.fasta $D/pacbio/phased.vcf.gz $D/pacbio/pacbio.bam
run_case pacbio-noref      $NR $D/pacbio/phased.vcf.gz $D/pacbio/pacbio.bam
run_case pacbio-cram       -r $D/pacbio/reference.fasta $D/pacbio/phased.vcf.gz $D/pacbio/pacbio.cram
run_case pacbio-retag      $NR $D/pacbio/phased.vcf.gz $D/pacbio/haplotagged.bam
run_case unmapped          $NR $D/haplotag.10X.vcf.gz $D/unmapped.bam
run_case triploid          $NR --ploidy 3 $D/haplotag_triploid.vcf.gz $D/haplotag_triploid.bam
run_case triploid-p2       $NR $D/haplotag_triploid.vcf.gz $D/haplotag_triploid.bam
run_case tetraploid        $NR --ploidy 4 $D/haplotag_poly.vcf.gz $D/haplotag_poly.bam
run_case reads-no-seq      $NR --ignore-read-groups $D/haplotag_1.vcf.gz $D/reads-no-sequence.bam
echo
echo "passed: $pass  failed: $fail ${failed[*]:-}"
[[ $fail -eq 0 ]]

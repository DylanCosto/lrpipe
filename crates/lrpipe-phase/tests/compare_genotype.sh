#!/usr/bin/env bash
# Compare `whatshap genotype` (Python, the oracle) with `lrpipe-phase genotype`
# on the WhatsHap test data: output VCF (without the ##commandline line),
# the --prioroutput VCF and the exit code.
#
# Usage: tests/compare_genotype.sh ENV_BIN WHATSHAP_SRC_DIR WORK_DIR
#   ENV_BIN: bin directory of an environment with WhatsHap 2.8 and samtools
#   WHATSHAP_SRC_DIR: a checkout of https://github.com/whatshap/whatshap (v2.8)
set -u
ENV=$1; SRC=$2; WORK=$3
RS=${RS:-$(dirname "$0")/../../../target/release/lrpipe-phase}
D=$SRC/tests/data
P=$D/pacbio
mkdir -p "$WORK/bam"
for sam in $D/*.sam $D/short-genome/*.sam; do
    b=$WORK/bam/$(basename "${sam%.sam}").bam
    [ -f "$b" ] || { "$ENV/samtools" sort -o "$b" "$sam" 2>/dev/null && "$ENV/samtools" index "$b"; }
done
B=$WORK/bam
pass=0; bothfail=0; fail=0; failed=()
norm='/^##commandline=/d'
run_case() {
    local name=$1; shift
    for t in wh rs; do
        local bin=("$ENV/whatshap"); [ $t = rs ] && bin=("$RS")
        local args=()
        for a in "$@"; do [ "$a" = PRIORS ] && args+=("$WORK/$name.$t.priors.vcf") || args+=("$a"); done
        "${bin[@]}" genotype -o "$WORK/$name.$t.vcf" "${args[@]}" > "$WORK/$name.$t.log" 2>&1
        echo $? > "$WORK/$name.$t.exit"
    done
    local e1 e2; e1=$(cat "$WORK/$name.wh.exit"); e2=$(cat "$WORK/$name.rs.exit")
    if [ "$e1" -ne 0 ] && [ "$e2" -ne 0 ]; then
        echo "BOTH-FAIL $name ($(tail -n1 "$WORK/$name.wh.log" | cut -c1-100))"
        pass=$((pass+1)); bothfail=$((bothfail+1)); return
    fi
    if [ "$e1" -ne "$e2" ]; then
        echo "FAIL $name: exit codes whatshap=$e1 rs=$e2"; tail -n1 "$WORK/$name.wh.log" "$WORK/$name.rs.log" | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    for f in vcf priors.vcf; do
        [ -e "$WORK/$name.wh.$f" ] || continue
        if ! cmp -s <(gzip -cdf "$WORK/$name.wh.$f" | sed "$norm") <(gzip -cdf "$WORK/$name.rs.$f" | sed "$norm"); then
            echo "FAIL $name: $f differs"; diff <(gzip -cdf "$WORK/$name.wh.$f" | sed "$norm") <(gzip -cdf "$WORK/$name.rs.$f" | sed "$norm") | head -6 | sed 's/^/    /'
            fail=$((fail+1)); failed+=("$name"); return
        fi
    done
    echo "ok   $name"
    pass=$((pass+1))
}
TP=$B/trio.pacbio.bam
run_case pacbio            -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam
run_case pacbio-noref      $P/variants.vcf $P/pacbio.bam
run_case pacbio-priors     -r $P/reference.fasta --only-snvs -p PRIORS $P/variants.vcf $P/pacbio.bam
run_case pacbio-nopriors   -r $P/reference.fasta --no-priors $P/variants.vcf $P/pacbio.bam
run_case pacbio-constant   -r $P/reference.fasta --constant 0.1 $P/variants.vcf $P/pacbio.bam
run_case pacbio-gq10       -r $P/reference.fasta --gt-qual-threshold 10 $P/variants.vcf $P/pacbio.bam
run_case pacbio-gq40       --gt-qual-threshold 40 $P/variants.vcf $P/pacbio.bam
run_case pacbio-H5         -r $P/reference.fasta -H 5 $P/variants.vcf $P/pacbio.bam
run_case pacbio-H1         -H 1 $P/variants.vcf $P/pacbio.bam
run_case pacbio-overhang   -r $P/reference.fasta --overhang 3 $P/variants.vcf $P/pacbio.bam
run_case pacbio-affine     -r $P/reference.fasta --affine-gap $P/variants.vcf $P/pacbio.bam
run_case pacbio-affine-mm  -r $P/reference.fasta --affine-gap --mismatch 4 $P/variants.vcf $P/pacbio.bam
run_case pacbio-mapq60     -r $P/reference.fasta --mapq 60 $P/variants.vcf $P/pacbio.bam
run_case pacbio-gz         -r $P/reference.fasta $P/phased.vcf.gz $P/pacbio.bam
run_case pacbio-gz-out     -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam
run_case pacbio-cram       -r $P/reference.fasta $P/variants.vcf $P/pacbio.cram
run_case pacbio-twobams    -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $P/haplotagged.bam
run_case pacbio-phasedvcf  -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $P/phased.vcf
run_case oneread           $D/onevariant.vcf $D/oneread.bam
run_case noreadgroup       --ignore-read-groups $D/onevariant.vcf $D/no-readgroup.bam
run_case trio-individual   $D/trio.vcf $TP
run_case trio-one-sample   --sample HG003 $D/trio.vcf $TP
run_case trio-two-chrom    $D/trio-two-chromosomes.vcf $TP
run_case trio-chrom-select --chromosome 1 -p PRIORS $D/trio-two-chromosomes.vcf $TP
run_case trio-chrom-excl   --exclude-chromosome 1 $D/trio-two-chromosomes.vcf $TP
run_case trio-symbolic     $D/trio-symbolic-alt.vcf $TP
run_case paired-end        $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
run_case short             $D/short-genome/short.vcf $B/short.bam
run_case short-ignore-rg   --ignore-read-groups $D/short-genome/short.vcf $B/short.bam
run_case short-wrongchrom  $D/short-genome/wrongchromosome.vcf $B/short.bam
run_case indels            -r $D/random0.fasta $D/indels.vcf $B/indels.bam
run_case indels-noref      $D/indels.vcf $B/indels.bam
run_case dist-geno         $D/test_dist_geno.vcf $B/test_dist_geno.bam
run_case haplotag-large    $D/haplotag.large.vcf.gz $D/haplotag.large.bam
run_case haplotag-10x      $D/haplotag.10X.vcf.gz $D/haplotag.10X.bam
run_case multisample-s1    -r $D/alleledetection.fasta --ignore-read-groups --sample sample1 $D/multisample.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.multiallelic.01.bam
run_case multiallelic      -r $P/reference.fasta $D/multiallelic.vcf $P/pacbio.bam
run_case missing-headers   $D/missing-headers.vcf $D/haplotag.bam
run_case duplicate-pos     -r $P/reference.fasta --ignore-read-groups $D/duplicate-positions.vcf $P/pacbio.bam
run_case phased-blocks     $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf
run_case alleledet-bi      -r $D/alleledetection.fasta $D/alleledetection.biallelic.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.biallelic.02.bam $D/alleledetection.biallelic.03.bam
run_case issue586          -r $D/issue-586/sba_vgsc_mapping_ref_oneline.fa $D/issue-586/MUT011_S351_splitfromjoint.vcf.gz $D/issue-586/MUT011_S351.bam
# --- pedigree genotyping
run_case ped-trio          --ped $D/trio.ped $D/trio.vcf $TP
run_case ped-trio-ref      -r $P/reference.fasta --ped $D/trio.ped $D/trio.vcf $TP
run_case ped-trio-genmap   --ped $D/trio.ped --genmap $D/trio.map --chromosome 1 $D/trio.vcf $TP
run_case ped-trio-recomb   --ped $D/trio.ped --recombrate 100 $D/trio.vcf $TP
run_case ped-trio-2chrom   --ped $D/trio.ped -p PRIORS $D/trio-two-chromosomes.vcf $TP
run_case ped-trio-nopriors --ped $D/trio.ped --no-priors $D/trio.vcf $TP
run_case ped-samples       --ped $D/pedigree.ped $D/ped_samples.vcf $D/ped_samples.bam
run_case ped-use-samples   --ped $D/pedigree.ped --use-ped-samples $D/ped_samples.vcf $D/ped_samples.bam
run_case ped-paired-end    --ped $D/trio_paired_end.ped $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
run_case ped-recomb-breaks --ped $D/recombination_breaks.ped $D/quartet.vcf.gz $B/recombination_breaks.sorted.bam
run_case ped-quartet2      --ped $D/quartet2.ped $D/quartet2.vcf $B/quartet2.bam
echo "genotype: passed: $pass (both tools failed: $bothfail)  failed: $fail ${failed[*]:-}"
[ "$fail" -eq 0 ]

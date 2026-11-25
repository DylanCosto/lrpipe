#!/usr/bin/env bash
# Compare `whatshap phase` (Python, the oracle) with `lrpipe-phase phase`.
#
# Usage: tests/compare.sh WHATSHAP_BIN SAMTOOLS_BIN WHATSHAP_SRC_DIR WORK_DIR
#   WHATSHAP_SRC_DIR: a checkout of https://github.com/whatshap/whatshap (v2.8)
#   WHATSHAP_BIN should be a Linux build (bioconda's, x86_64 or aarch64): bioconda's macOS
#   arm64 build breaks ties differently in 20 cases here (and in haplotag's 10x-cutoff), where
#   lrpipe-phase, on either system, gives what WhatsHap on Linux gives.
# Each case runs both tools with the same arguments and diffs the VCF output
# (ignoring the ##commandline line written by the tool itself) and, when
# requested, the --output-read-list file.
set -u
WH=$1; SAMTOOLS=$2; SRC=$3; WORK=$4
RS=${RS:-$(dirname "$0")/../../../target/release/lrpipe-phase}
D=$SRC/tests/data
mkdir -p "$WORK/bam"
# Convert test SAMs to sorted, indexed BAMs
for sam in $D/*.sam $D/short-genome/*.sam; do
    b=$WORK/bam/$(basename "${sam%.sam}").bam
    [ -f "$b" ] || { "$SAMTOOLS" sort -o "$b" "$sam" 2>/dev/null && "$SAMTOOLS" index "$b"; }
done
B=$WORK/bam
pass=0; fail=0; failed=()
run_case() {
    local name=$1; shift
    local args=("$@")
    local out_wh=$WORK/$name.wh.vcf out_rs=$WORK/$name.rs.vcf
    local rl_wh=$WORK/$name.wh.reads rl_rs=$WORK/$name.rs.reads
    rm -f "$WORK/$name".{wh,rs}.{recomb,gtchange}
    # NO_RECOMB=1: without --recombination-list (which HapChat cannot do)
    local rec_wh=(--recombination-list "$WORK/$name.wh.recomb") rec_rs=(--recombination-list "$WORK/$name.rs.recomb")
    if [ -n "${NO_RECOMB:-}" ]; then rec_wh=(); rec_rs=(); fi
    "$WH" phase -o "$out_wh" --output-read-list "$rl_wh" ${rec_wh[@]+"${rec_wh[@]}"} --changed-genotype-list "$WORK/$name.wh.gtchange" "${args[@]}" > "$WORK/$name.wh.log" 2>&1; local e1=$?
    "$RS" phase ${RS_EXTRA:-} -o "$out_rs" --output-read-list "$rl_rs" ${rec_rs[@]+"${rec_rs[@]}"} --changed-genotype-list "$WORK/$name.rs.gtchange" "${args[@]}" > "$WORK/$name.rs.log" 2>&1; local e2=$?
    if [ $e1 -ne 0 ] || [ $e2 -ne 0 ]; then
        if [ $e1 -ne 0 ] && [ $e2 -ne 0 ]; then
            echo "BOTH-FAIL $name (whatshap: $(tail -1 "$WORK/$name.wh.log" | cut -c1-100)) (rs: $(tail -1 "$WORK/$name.rs.log" | cut -c1-100))"
            pass=$((pass+1)); return
        fi
        echo "FAIL $name: exit codes whatshap=$e1 rs=$e2"; tail -2 "$WORK/$name.rs.log" "$WORK/$name.wh.log" | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    local norm='/^##commandline="\((whatshap|lrpipe-phase)/d'
    if ! cmp -s <(gzip -cdf "$out_wh" | sed -E "$norm") <(gzip -cdf "$out_rs" | sed -E "$norm"); then
        echo "FAIL $name: VCF differs"; diff <(gzip -cdf "$out_wh" | sed -E "$norm") <(gzip -cdf "$out_rs" | sed -E "$norm") | head -6 | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    if ! cmp -s <(sort "$rl_wh") <(sort "$rl_rs"); then
        echo "FAIL $name: read list differs"; diff <(sort "$rl_wh") <(sort "$rl_rs") | head -6 | sed 's/^/    /'
        fail=$((fail+1)); failed+=("$name"); return
    fi
    for ext in recomb gtchange; do
        if [ -e "$WORK/$name.wh.$ext" ] || [ -e "$WORK/$name.rs.$ext" ]; then
            if ! cmp -s "$WORK/$name.wh.$ext" "$WORK/$name.rs.$ext"; then
                echo "FAIL $name: $ext list differs"; diff "$WORK/$name.wh.$ext" "$WORK/$name.rs.$ext" | head -6 | sed 's/^/    /'
                fail=$((fail+1)); failed+=("$name"); return
            fi
        fi
    done
    echo "ok   $name"
    pass=$((pass+1))
}
P=$D/pacbio
run_case pacbio            -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam
run_case pacbio-noref      --no-reference $P/variants.vcf $P/pacbio.bam
run_case pacbio-cram       -r $P/reference.fasta $P/variants.vcf $P/pacbio.cram
run_case pacbio-snvs       -r $P/reference.fasta --only-snvs $P/variants.vcf $P/pacbio.bam
run_case pacbio-ds5        -r $P/reference.fasta --internal-downsampling 5 $P/variants.vcf $P/pacbio.bam
run_case pacbio-ds1        -r $P/reference.fasta --internal-downsampling 1 $P/variants.vcf $P/pacbio.bam
run_case pacbio-mapq60     -r $P/reference.fasta --mapq 60 $P/variants.vcf $P/pacbio.bam
run_case pacbio-gz         -r $P/reference.fasta $P/phased.vcf.gz $P/pacbio.bam
run_case pacbio-twobams    -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $P/haplotagged.bam
run_case oneread           --no-reference $D/onevariant.vcf $D/oneread.bam
run_case oneread-cram      -r $D/oneread-ref.fasta $D/onevariant.vcf $D/oneread.cram
run_case noreadgroup-ignore --ignore-read-groups --no-reference $D/onevariant.vcf $D/no-readgroup.bam
run_case trio-individual   --no-reference $D/trio.vcf $B/trio.pacbio.bam
run_case trio-one-sample   --no-reference --sample HG003 $D/trio.vcf $B/trio.pacbio.bam
run_case trio-merged       --no-reference $D/trio-merged-blocks.vcf $B/trio-merged-blocks.bam
run_case trio-missing-gt   --no-reference $D/trio-missing-genotypes.vcf $B/trio.pacbio.bam
run_case trio-two-chrom    --no-reference $D/trio-two-chromosomes.vcf $B/trio.pacbio.bam
run_case trio-chrom-select --no-reference --chromosome 1 $D/trio-two-chromosomes.vcf $B/trio.pacbio.bam
run_case trio-chrom-excl   --no-reference --exclude-chromosome 1 $D/trio-two-chromosomes.vcf $B/trio.pacbio.bam
run_case trio-symbolic     --no-reference $D/trio-symbolic-alt.vcf $B/trio.pacbio.bam
run_case trio-mendel       --no-reference $D/trio-mendelian-conflict.vcf $B/trio.pacbio.bam
run_case paired-end        --no-reference $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
run_case recomb-breaks     --no-reference $D/quartet.vcf.gz $B/recombination_breaks.sorted.bam
run_case quartet2          --no-reference $D/quartet2.vcf $B/quartet2.bam
run_case short             --no-reference $D/short-genome/short.vcf $B/short.bam
cp -n $D/short-genome/short-ref.fasta $WORK/short-ref.fasta; [ -f $WORK/short-ref.fasta.fai ] || "$SAMTOOLS" faidx $WORK/short-ref.fasta
run_case short-ref         -r $WORK/short-ref.fasta --ignore-read-groups $D/short-genome/short.vcf $B/short.bam
run_case short-ignore-rg   --no-reference --ignore-read-groups $D/short-genome/short.vcf $B/short.bam
run_case short-dup         --no-reference $D/short-genome/short.vcf $B/short-one-read-duplicate.bam
run_case short-wrongchrom  --no-reference $D/short-genome/wrongchromosome.vcf $B/short.bam
run_case indels            -r $D/random0.fasta $D/indels.vcf $B/indels.bam
run_case indels-noref      --no-reference $D/indels.vcf $B/indels.bam
run_case dist-geno         --no-reference $D/test_dist_geno.vcf $B/test_dist_geno.bam
run_case haplotag          --no-reference $D/haplotag_1.vcf.gz $D/haplotag.bam
run_case haplotag-large    --no-reference $D/haplotag.large.vcf.gz $D/haplotag.large.bam
run_case haplotag-10x      --no-reference $D/haplotag.10X.vcf.gz $D/haplotag.10X.bam
run_case haplotag-noRG     --no-reference --ignore-read-groups $D/haplotag_noRG.vcf.gz $D/haplotag_noRG.bam
run_case haplotag-sample   --no-reference $D/haplotag_sample.vcf.gz $D/haplotag_sample.bam
run_case haplotag-bcf      --no-reference $D/haplotag_2.bcf $D/haplotag.bam
run_case supplementary     --no-reference $D/haplotag.supplementary.vcf.gz $D/haplotag.supplementary.bam
run_case supplementary-use --no-reference --use-supplementary $D/haplotag.supplementary.vcf.gz $D/haplotag.supplementary.bam
run_case supp-strategy     --no-reference --use-supplementary $D/supplementary_strategy_test.grch38.vcf.gz $D/supplementary_strategy_test.grch38.bam
run_case ped-samples       --no-reference $D/ped_samples.vcf $D/ped_samples.bam
run_case multisample-s1    -r $D/alleledetection.fasta --ignore-read-groups --sample sample1 $D/multisample.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.multiallelic.01.bam
run_case multisample-s2    -r $D/alleledetection.fasta --ignore-read-groups --sample sample2 $D/multisample.vcf $D/alleledetection.biallelic.02.bam $D/alleledetection.multiallelic.01.bam
run_case multisample-both  -r $D/alleledetection.fasta --ignore-read-groups --sample sample1 --sample sample2 $D/multisample.vcf $D/alleledetection.biallelic.01.bam
run_case issue586          -r $D/issue-586/sba_vgsc_mapping_ref_oneline.fa $D/issue-586/MUT011_S351_splitfromjoint.vcf.gz $D/issue-586/MUT011_S351.bam
run_case issue586-noref    --no-reference $D/issue-586/MUT011_S351_splitfromjoint.vcf.gz $D/issue-586/MUT011_S351.bam
run_case missing-headers   --no-reference $D/missing-headers.vcf $D/haplotag.bam
run_case duplicate-pos     -r $P/reference.fasta --ignore-read-groups $D/duplicate-positions.vcf $P/pacbio.bam
run_case unknown-gt        -r $D/alleledetection.fasta --ignore-read-groups --sample sample1 $D/unknown-genotype.vcf $D/alleledetection.biallelic.01.bam
run_case unknown-gt-s2     -r $D/alleledetection.fasta --ignore-read-groups --sample sample2 $D/unknown-genotype.vcf $D/alleledetection.biallelic.03.bam
run_case multiallelic      -r $P/reference.fasta $D/multiallelic.vcf $P/pacbio.bam
run_case phased-blocks     --no-reference $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam
run_case polyploid         --no-reference $D/polyploid.chr22.42M.12k.vcf $D/polyploid.chr22.42M.12k.bam
run_case alleledet-ignore  -r $D/alleledetection.fasta --ignore-read-groups $D/alleledetection.biallelic.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.biallelic.02.bam
run_case alleledet-bi      -r $D/alleledetection.fasta $D/alleledetection.biallelic.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.biallelic.02.bam $D/alleledetection.biallelic.03.bam $D/alleledetection.biallelic.04.bam $D/alleledetection.biallelic.05.bam
run_case alleledet-bi-noref --no-reference $D/alleledetection.biallelic.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.biallelic.02.bam
run_case alleledet-multi   -r $D/alleledetection.fasta $D/alleledetection.multiallelic.vcf $D/alleledetection.multiallelic.01.bam
# --- pedigree phasing
TP=$B/trio.pacbio.bam
run_case ped-trio          --no-reference --ped $D/trio.ped $D/trio.vcf $TP
run_case ped-trio-ref      -r $D/pacbio/reference.fasta --ped $D/trio.ped $D/trio.vcf $TP
run_case ped-trio-genmap   --no-reference --ped $D/trio.ped --genmap $D/trio.map --chromosome 1 $D/trio.vcf $TP
run_case ped-trio-zerogd   --no-reference --ped $D/trio.ped --genmap $D/zero-genetic-distance.map --chromosome 1 $D/trio.vcf $TP
run_case ped-trio-recomb   --no-reference --ped $D/trio.ped --recombrate 100 $D/trio.vcf $TP
run_case ped-trio-nogh     --no-reference --ped $D/trio.ped --no-genetic-haplotyping $D/trio.vcf $TP
run_case ped-trio-ds5      --no-reference --ped $D/trio.ped --internal-downsampling 5 $D/trio.vcf $TP
RS_EXTRA="-t 3" run_case ped-trio-threads  --no-reference --ped $D/trio.ped $D/trio-two-chromosomes.vcf $TP
run_case ped-trio-2chrom   --no-reference --ped $D/trio.ped $D/trio-two-chromosomes.vcf $TP
run_case ped-trio-merged   --no-reference --ped $D/trio.ped $D/trio-merged-blocks.vcf $B/trio-merged-blocks.bam
run_case ped-trio-missing  --no-reference --ped $D/trio.ped $D/trio-missing-genotypes.vcf $TP
run_case ped-trio-mendel   --no-reference --ped $D/trio.ped $D/trio-mendelian-conflict.vcf $TP
run_case ped-symbolic-only --no-reference --ped $D/trio.ped $D/trio-symbolic-alt.vcf
run_case ped-gh-only       --no-reference --ped $D/genetic-haplotyping.ped $D/genetic-haplotyping.vcf
run_case ped-gh-only-nogh  --no-reference --ped $D/genetic-haplotyping.ped --no-genetic-haplotyping $D/genetic-haplotyping.vcf
run_case ped-samples       --no-reference --ped $D/pedigree.ped $D/ped_samples.vcf $D/ped_samples.bam
run_case ped-samples-sel   --no-reference --ped $D/pedigree.ped --sample HG002 --sample HG003 --sample HG004 $D/ped_samples.vcf $D/ped_samples.bam
run_case ped-paired-end    --no-reference --ped $D/trio_paired_end.ped $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
run_case ped-recomb-breaks --no-reference --ped $D/recombination_breaks.ped $D/quartet.vcf.gz $B/recombination_breaks.sorted.bam
run_case ped-recomb-genmap --no-reference --ped $D/recombination_breaks.ped --genmap $D/recombination_breaks.map --chromosome 1 $D/quartet.vcf.gz $B/recombination_breaks.sorted.bam
run_case ped-quartet2      --no-reference --ped $D/quartet2.ped $D/quartet2.vcf $B/quartet2.bam
run_case ped-quartet2-gl   --no-reference --ped $D/quartet2.ped --distrust-genotypes $D/quartet2.vcf $B/quartet2.bam
run_case ped-trio-gl       --no-reference --ped $D/trio.ped --distrust-genotypes $D/trio_genotype_likelihoods.vcf $TP
run_case ped-trio-loggl    --no-reference --ped $D/trio.ped --distrust-genotypes $D/trio_genotype_log_likelihoods.vcf $TP
run_case ped-trio-gl-reg   --no-reference --ped $D/trio.ped --distrust-genotypes --gl-regularizer 0.01 $D/trio_genotype_likelihoods.vcf $TP
run_case ped-trio-hp       --no-reference --ped $D/trio.ped --tag HP $D/trio.vcf $TP
# --- --merge-reads
run_case merge-pacbio      -r $P/reference.fasta --merge-reads $P/variants.vcf $P/pacbio.bam
run_case merge-pacbio-ds5  -r $P/reference.fasta --merge-reads --internal-downsampling 5 $P/variants.vcf $P/pacbio.bam
run_case merge-paired-end  --no-reference --merge-reads $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
run_case merge-short       --no-reference --merge-reads $D/short-genome/short.vcf $B/short.bam
run_case merge-hap-large   --no-reference --merge-reads $D/haplotag.large.vcf.gz $D/haplotag.large.bam
run_case merge-multisample -r $D/alleledetection.fasta --ignore-read-groups --merge-reads --sample sample1 --sample sample2 $D/multisample.vcf $D/alleledetection.biallelic.01.bam
run_case merge-trio        --no-reference --ped $D/trio.ped --merge-reads $D/trio.vcf $TP
run_case merge-vcf-blocks  --no-reference --merge-reads $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf
# --- --algorithm heuristic
run_case heur-pacbio       -r $P/reference.fasta --algorithm heuristic $P/variants.vcf $P/pacbio.bam
run_case heur-pacbio-L16   -r $P/reference.fasta --algorithm heuristic --row-limit 16 $P/variants.vcf $P/pacbio.bam
run_case heur-pacbio-L1    -r $P/reference.fasta --algorithm heuristic -L 1 $P/variants.vcf $P/pacbio.bam
run_case heur-pacbio-merge -r $P/reference.fasta --algorithm heuristic --merge-reads $P/variants.vcf $P/pacbio.bam
run_case heur-paired-end   --no-reference --algorithm heuristic $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
run_case heur-short        --no-reference --algorithm heuristic $D/short-genome/short.vcf $B/short.bam
run_case heur-hap-large    --no-reference --algorithm heuristic $D/haplotag.large.vcf.gz $D/haplotag.large.bam
run_case heur-dist-geno    --no-reference --algorithm heuristic --distrust-genotypes $D/test_dist_geno.vcf $B/test_dist_geno.bam
run_case heur-multisample  -r $D/alleledetection.fasta --ignore-read-groups --algorithm heuristic --sample sample1 --sample sample2 $D/multisample.vcf $D/alleledetection.biallelic.01.bam
run_case heur-trio         --no-reference --algorithm heuristic --ped $D/trio.ped $D/trio.vcf $TP
run_case heur-trio-L4      --no-reference --algorithm heuristic -L 4 --ped $D/trio.ped $D/trio.vcf $TP
run_case heur-trio-gl      --no-reference --algorithm heuristic --ped $D/trio.ped --distrust-genotypes $D/trio_genotype_likelihoods.vcf $TP
run_case heur-trio-mendel  --no-reference --algorithm heuristic --ped $D/trio.ped $D/trio-mendelian-conflict.vcf $TP
run_case heur-trio-recomb  --no-reference --algorithm heuristic --ped $D/trio.ped --recombrate 100 $D/trio.vcf $TP
run_case heur-quartet2     --no-reference --algorithm heuristic --ped $D/quartet2.ped $D/quartet2.vcf $B/quartet2.bam
run_case heur-recomb-breaks --no-reference --algorithm heuristic --ped $D/recombination_breaks.ped $D/quartet.vcf.gz $B/recombination_breaks.sorted.bam
run_case heur-vcf-blocks   --no-reference --algorithm heuristic $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf
run_case heur-L-ignored    --no-reference -L 8 $D/short-genome/short.vcf $B/short.bam
# --- --algorithm hapchat
NO_RECOMB=1 run_case hapchat-pacbio    -r $P/reference.fasta --algorithm hapchat $P/variants.vcf $P/pacbio.bam
NO_RECOMB=1 run_case hapchat-paired    --no-reference --algorithm hapchat $D/paired_end.sorted.vcf $B/paired_end.sorted.bam
NO_RECOMB=1 run_case hapchat-short     --no-reference --algorithm hapchat $D/short-genome/short.vcf $B/short.bam
NO_RECOMB=1 run_case hapchat-short-gl  --no-reference --algorithm hapchat --distrust-genotypes $D/short-genome/short.vcf $B/short.bam
NO_RECOMB=1 run_case hapchat-hap-large --no-reference --algorithm hapchat $D/haplotag.large.vcf.gz $D/haplotag.large.bam
NO_RECOMB=1 run_case hapchat-dist-geno --no-reference --algorithm hapchat --distrust-genotypes $D/test_dist_geno.vcf $B/test_dist_geno.bam
NO_RECOMB=1 run_case hapchat-indels    -r $D/random0.fasta --algorithm hapchat $D/indels.vcf $B/indels.bam
NO_RECOMB=1 run_case hapchat-merge     -r $P/reference.fasta --algorithm hapchat --merge-reads $P/variants.vcf $P/pacbio.bam
NO_RECOMB=1 run_case hapchat-ds5       -r $P/reference.fasta --algorithm hapchat --internal-downsampling 5 $P/variants.vcf $P/pacbio.bam
run_case hapchat-trio      --no-reference --algorithm hapchat --ped $D/trio.ped $D/trio.vcf $TP
run_case hapchat-recomb    --no-reference --algorithm hapchat $D/short-genome/short.vcf $B/short.bam
NO_RECOMB=1 run_case hapchat-vcf-blocks --no-reference --algorithm hapchat $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf
NO_RECOMB=1 run_case hapchat-alleledet -r $D/alleledetection.fasta --algorithm hapchat $D/alleledetection.biallelic.vcf $D/alleledetection.biallelic.01.bam $D/alleledetection.biallelic.02.bam
# --- phased VCFs as phasing input
run_case vcf-blocks        --no-reference $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf
run_case vcf-blocks-only   --no-reference $D/phased-blocks.variants.vcf $D/phased-blocks.blocks.vcf
run_case vcf-blocks-ds1    --no-reference --internal-downsampling 1 $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf
run_case vcf-via-ps        -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $D/phased-via-PS.vcf
run_case vcf-via-hp        -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $D/phased-via-HP.vcf
run_case vcf-via-mixed     -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $D/phased-via-mixed-HP-PS.vcf
run_case vcf-pacbio-phased -r $P/reference.fasta $P/variants.vcf $P/pacbio.bam $P/phased.vcf
run_case vcf-two           --no-reference $D/phased-blocks.variants.vcf $D/phased-blocks.reads.bam $D/phased-blocks.blocks.vcf $D/phased-blocks.blocks.vcf
run_case vcf-ped           --no-reference --ped $D/trio.ped $D/trio.vcf $TP $D/phased1.vcf
# --- distrust genotypes
run_case distrust-pacbio   -r $P/reference.fasta --distrust-genotypes $P/variants.vcf $P/pacbio.bam
run_case distrust-homo     -r $P/reference.fasta --distrust-genotypes --include-homozygous $P/variants.vcf $P/pacbio.bam
run_case distrust-gq50     -r $P/reference.fasta --distrust-genotypes --default-gq 50 $P/variants.vcf $P/pacbio.bam
run_case distrust-dist     --no-reference --distrust-genotypes $D/test_dist_geno.vcf $B/test_dist_geno.bam
run_case distrust-dist-hom --no-reference --distrust-genotypes --include-homozygous $D/test_dist_geno.vcf $B/test_dist_geno.bam
run_case distrust-gl       --no-reference --distrust-genotypes $D/trio_genotype_likelihoods.vcf $TP
run_case distrust-gl-reg   --no-reference --distrust-genotypes --gl-regularizer 0.1 --include-homozygous $D/trio_genotype_likelihoods.vcf $TP
# --- HP tag
run_case hp-pacbio         -r $P/reference.fasta --tag HP $P/variants.vcf $P/pacbio.bam
run_case hp-phased-input   -r $P/reference.fasta --tag HP $P/phased.vcf $P/pacbio.bam
run_case hp-trio           --no-reference --tag HP $D/trio.vcf $TP
run_case hp-retag          --no-reference --tag HP $D/phased-via-HP.vcf $P/pacbio.bam
echo
echo "passed: $pass  failed: $fail ${failed[*]:-}"
[[ $fail -eq 0 ]]

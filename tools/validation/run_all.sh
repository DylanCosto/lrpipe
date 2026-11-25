#!/usr/bin/env bash
# Re-runs every lrpipe tool (lrpipe 0.1.3, commit 47a70e0) on whole genomes, and the original tools
# where no timing of theirs is kept, each step timed with /usr/bin/time -v into logs/STEP.time.
# A step whose done/STEP marker exists is skipped, so the script can be run again after a failure.
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/paths.sh"
V=$VALIDATION_DIR
L=$LRPIPE
HG003=$HG003_ONT_BAM
HIFI=$HG002_HIFI_BAM
CRAM=$HG002_ONT_CRAM
REVIO=$HG002_REVIO_BAM
CATALOG=$TR_CATALOG
SB=$SPECTRE_BENCH
mkdir -p "$V/logs" "$V/done" "$V/out" || exit 1
cd "$V/out" || exit 1
failed=0

# step NAME COMMAND...: runs COMMAND timed, in out/, unless done/NAME exists.
step() {
  local name=$1; shift
  if [[ -e $V/done/$name ]]; then echo "skip $name"; return 0; fi
  echo "[$(date -u +%FT%TZ)] start $name"
  if /usr/bin/time -v -o $V/logs/$name.time "$@" > $V/logs/$name.log 2>&1; then
    touch $V/done/$name
    echo "[$(date -u +%FT%TZ)] done $name: $(grep Elapsed $V/logs/$name.time | awk '{print $NF}')"
  else
    echo "[$(date -u +%FT%TZ)] FAILED $name (see logs/$name.log)"
    failed=1
  fi
}

snv() { # snv NAME BAM PLATFORM MODEL SAMPLE [extra...]
  local name=$1 bam=$2 platform=$3 model=$4 sample=$5; shift 5
  step $name $L snv run --bam $bam --ref $REF --model-path $MODELS/$model --platform $platform \
    --sample-name $sample --threads 16 --output $V/out/$name --remove-intermediate-dir=true "$@"
}

# lrpipe-snv: exact CPU, approximate CPU and GPU, HG003 ONT and HG002 HiFi.
snv snv-ont-cpu-exact $HG003 ont r1041_e82_400bps_sup_v500 HG003
snv snv-ont-cpu-approx $HG003 ont r1041_e82_400bps_sup_v500 HG003 --approximate
snv snv-ont-gpu $HG003 ont r1041_e82_400bps_sup_v500 HG003 --device cuda:0 --nvrtc $NVRTC
snv snv-hifi-cpu-exact $HIFI hifi hifi HG002
snv snv-hifi-cpu-approx $HIFI hifi hifi HG002 --approximate
snv snv-hifi-gpu $HIFI hifi hifi HG002 --device cuda:0 --nvrtc $NVRTC

# lrpipe sv on the three whole genomes Sniffles' outputs are kept for.
step sv-hg003-ont $L sv --input $HG003 --vcf $V/out/sv-hg003-ont.vcf --reference $SNIFFLES_REF --threads 16 --allow-overwrite
step sv-hg002-cram $L sv --input $CRAM --vcf $V/out/sv-hg002-cram.vcf --reference $REF -t 16 --allow-overwrite
step sv-hg002-hifi $L sv --input $HIFI --vcf $V/out/sv-hg002-hifi.vcf --reference $REF --threads 16 --allow-overwrite

# depth: lrpipe and mosdepth 0.3.14, HG002 HiFi.
step depth-lrpipe $L depth -n -x -Q 20 -b 1000 -t 16 $V/out/depth-lrpipe $HIFI
step depth-mosdepth $MOSDEPTH -n -x -Q 20 -b 1000 -t 16 $V/out/depth-mosdepth $HIFI

# CNVs: lrpipe and Spectre 0.2.1 on the HG002 HiFi coverage of Spectre's reference run.
CNV_ARGS=(CNVCaller -c $SB/inputs/hg002_hifi_wgs/mosdepth.regions.bed.gz -s HG002 -r $REF
  -m $SB/inputs/grch38_reference.mdr -l $SB/inputs/grch38_blacklist_spectre.bed)
step cnv-lrpipe $L cnv "${CNV_ARGS[@]}" -d $V/out/cnv-lrpipe
mkdir -p $V/out/cnv-spectre
# Spectre's own entry point, numpy seeded as the reference runs are (its dip test draws from it).
step cnv-spectre env SPECTRE_SRC=$SB/spectre-0.2.1 PATH=$(dirname $SPECTRE_PYTHON):/usr/bin:/bin \
  $SPECTRE_PYTHON -c \
  'import os, sys; sys.path.insert(0, os.environ["SPECTRE_SRC"]); import numpy as np; np.random.seed(0); from spectre import spectre as m; sys.argv = ["spectre"] + sys.argv[1:]; m.main()' \
  "${CNV_ARGS[@]}" -d $V/out/cnv-spectre

# tandem repeats: lrpipe and TRGT 5.1.0, HG002 HiFi, 937,204 loci.
step tr-lrpipe $L tr --reads $HIFI --genome $REF --repeats $CATALOG --output $V/out/tr-lrpipe.vcf.gz --sample HG002 --karyotype xy --platform hifi --threads 16
step tr-trgt $TRGT genotype --genome $REF --repeats $CATALOG --reads $HIFI --output-prefix $V/out/tr-trgt --karyotype XY --threads 16

# methylation: lrpipe and modkit 0.6.4 (HG002 ONT CRAM), lrpipe and pb-CpG-tools 3.0.0 (Revio).
step methyl-pileup-lrpipe $L methyl pileup $CRAM $V/out/methyl-lrpipe.bed --ref $REF -t 16 --modified-bases 5mC 5hmC --cpg --combine-strands
step methyl-pileup-modkit $MODKIT pileup $CRAM $V/out/methyl-modkit.bed --ref $REF -t 16 --modified-bases 5mC 5hmC --cpg --combine-strands
step methyl-cpg-lrpipe $L methyl cpg --bam $REVIO --ref $REF --output-prefix $V/out/cpg-lrpipe --modsites-mode reference --pileup-mode count --threads 16
step methyl-cpg-pb $PB_CPG --bam $REVIO --ref $REF --output-prefix $V/out/cpg-pb --modsites-mode reference --pileup-mode count --threads 16

# phasing: the HiFi pipeline's final calls, phased by lrpipe-phase (all contigs, 16 threads) and
# by WhatsHap 2.8 (one process per contig, 6 at a time).
PH=$V/out/phase; mkdir -p $PH/contigs $PH/whatshap
CALLS=$HG002_HIFI_RUN/snv/merge_output.vcf.gz
if [[ ! -e $V/done/phase-split ]]; then
  if zcat "$CALLS" | awk -v d="$PH/contigs" 'BEGIN{FS="\t"} /^#/{h=h $0 "\n"; next} {f=d "/" $1 ".vcf"; if(!(f in seen)){seen[f]=1; printf "%s", h > f} print > f} END {if (!NR || !length(seen)) exit 1}'; then
    touch "$V/done/phase-split" || exit 1
  else
    echo "[$(date -u +%FT%TZ)] FAILED phase-split (could not read $CALLS)"
    failed=1
  fi
fi
step phase-lrpipe $L phase --ignore-read-groups --reference $REF -t 16 -o $PH/lrpipe.vcf $CALLS $HIFI
phase_whatshap() {
  shopt -s nullglob
  local contigs=("$PH"/contigs/*.vcf)
  ((${#contigs[@]})) || return 1
  printf '%s\0' "${contigs[@]}" | xargs -0 -P 6 -I{} bash -c '
    contig=$(basename "$1" .vcf)
    "$CONDA_BIN/whatshap" phase --ignore-read-groups --reference "$REF" --chromosome "$contig" \
      -o "$PH/whatshap/$contig.vcf" "$1" "$HIFI" > "$PH/whatshap/$contig.log" 2>&1
  ' _ {}
}
export -f phase_whatshap; export PH REF HIFI CONDA_BIN
step phase-whatshap bash -c phase_whatshap
echo "[$(date -u +%FT%TZ)] all steps finished"
exit "$failed"

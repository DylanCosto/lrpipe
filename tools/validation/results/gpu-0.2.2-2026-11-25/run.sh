#!/usr/bin/env bash
# lrpipe 0.2.2 (the published bundle) on the GPU, both whole genomes, then GIAB scoring and the
# Clair3 GPU parity check, as tools/validation/run_all.sh and compare_all.sh did.
set -uo pipefail
R=/home/dyl/Documents/lrpipe/tools/validation
V=/data2/bioinformatics/lrpipe-build/val022
mkdir -p $V/bundle $V/logs $V/out
cd $V/bundle
[[ -x lrpipe-0.2.2-x86_64-linux/bin/lrpipe ]] || { gh release download v0.2.2 -R DylanCosto/lrpipe -p 'lrpipe-0.2.2-x86_64-linux.tar.gz*' --clobber && sha256sum -c *.sha256 && tar xzf lrpipe-0.2.2-x86_64-linux.tar.gz; }
export VALIDATION_DIR=$V LRPIPE=$V/bundle/lrpipe-0.2.2-x86_64-linux/bin/lrpipe
. $R/paths.sh
$LRPIPE --version
run() { # NAME BAM PLATFORM MODEL SAMPLE
  local name=$1 bam=$2 platform=$3 model=$4 sample=$5
  [[ -e $V/out/$name/merge_output.vcf.gz ]] && { echo "skip $name"; return; }
  /usr/bin/time -v -o $V/logs/$name.time $LRPIPE snv run --bam $bam --ref $REF --model-path $MODELS/$model \
    --platform $platform --sample-name $sample --threads 16 --output $V/out/$name --remove-intermediate-dir=true \
    --device cuda:0 --nvrtc $NVRTC > $V/logs/$name.log 2>&1 || { echo "FAILED $name"; tail -5 $V/logs/$name.log; exit 1; }
  echo "$name: $(grep -E 'Elapsed' $V/logs/$name.time | awk '{print $NF}'), $(awk -F': ' '/User time/{u=$2}/System time/{s=$2}END{printf "%.0f CPU-s", u+s}' $V/logs/$name.time), $(awk -F': ' '/Maximum resident/{printf "%.1f GB", $2/1048576}' $V/logs/$name.time)"
}
run snv-ont-gpu-0.2.2 $HG003_ONT_BAM ont r1041_e82_400bps_sup_v500 HG003
run snv-hifi-gpu-0.2.2 $HG002_HIFI_BAM hifi hifi HG002
SNV=/home/dyl/Documents/lrpipe/lrpipe-snv/tools/concordance.py
echo "== parity"
echo " hifi gpu 0.2.2 vs Clair3 GPU: $(python3 -I $SNV $CLAIR3_HG002_HIFI_GPU $V/out/snv-hifi-gpu-0.2.2/merge_output.vcf.gz | tail -n 1)"
echo " ont gpu 0.2.2 vs Clair3 CPU: $(python3 -I $SNV $CLAIR3_HG003_ONT_CPU $V/out/snv-ont-gpu-0.2.2/merge_output.vcf.gz | tail -n 1)"
O=/home/dyl/bench-data/validation/out
echo " ont gpu 0.2.2 vs its CPU exact run: $(python3 -I $SNV $O/snv-ont-cpu-exact/merge_output.vcf.gz $V/out/snv-ont-gpu-0.2.2/merge_output.vcf.gz | tail -n 1)"
echo " hifi gpu 0.2.2 vs its CPU exact run: $(python3 -I $SNV $O/snv-hifi-cpu-exact/merge_output.vcf.gz $V/out/snv-hifi-gpu-0.2.2/merge_output.vcf.gz | tail -n 1)"
echo " ont gpu 0.2.2 vs ont gpu before: $(python3 -I $SNV $O/snv-ont-gpu/merge_output.vcf.gz $V/out/snv-ont-gpu-0.2.2/merge_output.vcf.gz | tail -n 1)"
echo " hifi gpu 0.2.2 vs hifi gpu before: $(python3 -I $SNV $O/snv-hifi-gpu/merge_output.vcf.gz $V/out/snv-hifi-gpu-0.2.2/merge_output.vcf.gz | tail -n 1)"
echo "== GIAB"
$R/score2.sh HG003 snv-ont-gpu-0.2.2 $V/out/snv-ont-gpu-0.2.2/merge_output.vcf.gz
$R/score2.sh HG002 snv-hifi-gpu-0.2.2 $V/out/snv-hifi-gpu-0.2.2/merge_output.vcf.gz

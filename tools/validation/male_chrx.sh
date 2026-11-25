#!/usr/bin/env bash
# LRP-003: HG002 HiFi chrX and chrY called as the pipeline now calls an XY sample, then phased.
set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/paths.sh"
D=$VALIDATION_DIR/male-chrx
L=$LRPIPE
HIFI=$HG002_HIFI_BAM
$L snv run --bam $HIFI --ref $REF --model-path $MODELS/hifi \
  --platform hifi --sample-name HG002 --threads 16 --device cuda:0 \
  --nvrtc $NVRTC --output $D/snv \
  --ctg-name chrX,chrY --gender male --par-regions-bed $D/par.bed --approximate --enable-phasing=true \
  > $D/snv.log 2>&1
$L phase --ignore-read-groups --reference $REF -t 16 -o $D/phased.vcf $D/snv/merge_output.vcf.gz $HIFI > $D/phase.log 2>&1
echo ok

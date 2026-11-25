#!/usr/bin/env bash
# The tools whose output the review fixes can change, rerun with the fixed binary into out-fixed/.
set -uo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/paths.sh"
V=$VALIDATION_DIR
L=$LRPIPE
HG003=$HG003_ONT_BAM
HIFI=$HG002_HIFI_BAM
CRAM=$HG002_ONT_CRAM
CATALOG=$TR_CATALOG
SB=$SPECTRE_BENCH
F=$V/out-fixed; mkdir -p $F/cnv-lrpipe
cd $F
$L sv --input $HG003 --vcf $F/sv-hg003-ont.vcf --reference $SNIFFLES_REF --threads 16 --allow-overwrite > sv1.log 2>&1
$L sv --input $CRAM --vcf $F/sv-hg002-cram.vcf --reference $REF -t 16 --allow-overwrite > sv2.log 2>&1
$L sv --input $HIFI --vcf $F/sv-hg002-hifi.vcf --reference $REF --threads 16 --allow-overwrite > sv3.log 2>&1
$L cnv CNVCaller -c $SB/inputs/hg002_hifi_wgs/mosdepth.regions.bed.gz -s HG002 -r $REF \
  -m $SB/inputs/grch38_reference.mdr -l $SB/inputs/grch38_blacklist_spectre.bed -d $F/cnv-lrpipe > cnv.log 2>&1
$L tr --reads $HIFI --genome $REF --repeats $CATALOG --output $F/tr-lrpipe.vcf.gz --sample HG002 --karyotype xy --platform hifi --threads 16 > tr.log 2>&1
$L methyl pileup $CRAM $F/methyl-lrpipe.bed --ref $REF -t 16 --modified-bases 5mC 5hmC --cpg --combine-strands > methyl.log 2>&1
echo "exit codes done"

#!/bin/bash
# Builds the MergeVcf parity inputs under $WORK/inputs from the m4 oracle (real pileup and
# full-alignment VCFs of HG003 chr20) plus hand-written edge cases, and copies
# merge_vcf_cases.tsv to $WORK/cases.tsv. Then run merge_vcf_python.sh to make Clair3's
# outputs in $WORK/py, and the parity test with SNV_MERGEVCF_CASES=$WORK.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
ORACLE=${ORACLE:-/data2/bioinformatics/clair3-rs-bench/m4-oracle}
WORK=${WORK:-/data2/bioinformatics/clair3-rs-bench/m4-mergevcf-work}
IN=$WORK/inputs
rm -rf "$IN"; mkdir -p "$IN"
cp "$HERE/merge_vcf_cases.tsv" "$WORK/cases.tsv"
cp "$ORACLE/pileup.vcf.gz" "$IN/pileup.vcf.gz"
cp "$ORACLE/full_alignment.vcf.gz" "$IN/fa.vcf.gz"
# The same calls on chrX (plain text), for the gender and PAR cases.
zcat "$ORACLE/pileup.vcf.gz" | sed 's/^chr20\t/chrX\t/' > "$IN/pileup_x.vcf"
zcat "$ORACLE/full_alignment.vcf.gz" | sed 's/^chr20\t/chrX\t/' > "$IN/fa_x.vcf"
# A contig with no full-alignment calls: header only, and records of another contig only.
zcat "$ORACLE/full_alignment.vcf.gz" | grep '^#' > "$IN/fa_header.vcf"
zcat "$ORACLE/full_alignment.vcf.gz" | sed 's/^chr20\t/chr21\t/' | gzip > "$IN/fa_chr21.vcf.gz"
# PAR BED over parts of chrX: a comment, a zero-length row and a row for chrY.
P1=$(grep -v '^#' "$IN/pileup_x.vcf" | sed -n 20000p | cut -f2)
P2=$(grep -v '^#' "$IN/pileup_x.vcf" | sed -n 90000p | cut -f2)
P3=$(grep -v '^#' "$IN/pileup_x.vcf" | sed -n 150000p | cut -f2)
P4=$(grep -v '^#' "$IN/fa_x.vcf" | sed -n 60000p | cut -f2)
printf '#PAR\nchrX\t0\t%s\nchrX\t%s\t%s\nchrX\t%s\t%s\nchrY\t10\t20\n' "$P1" "$P2" "$P3" "$P4" "$P4" | gzip > "$IN/par.bed.gz"

H=$IN/fa_header.vcf
row() { # contig pos ref alt qual filter gt [info]
  printf '%s\t%s\t.\t%s\t%s\t%s\t%s\t%s\tGT:GQ:DP:AD:AF\t%s:9:10:4,6:0.6\n' "$1" "$2" "$3" "$4" "$5" "$6" "${8:-P}" "$7"
}
# Overlaps and duplicates: full-alignment calls (unsorted, duplicated, RefCall, low QUAL,
# multi-allelic, phased) replacing pileup calls at the same sites; pileup calls with no
# full-alignment record; identical rows in both; another contig.
{
  cat "$H"
  row chr20 100 A C 12.5 PASS 0/1 F
  row chr20 300 A . 3.1 RefCall 0/0 F
  row chr20 200 G T 1.5 PASS 1/1 F
  row chr20 200 G T 1.5 PASS 1/1 F
  row chr20 400 C A,T 30 PASS 1/2 F
  row chr20 500 T G 2 PASS '1|1' F
  row chr20 600 T . 20 RefCall 0/0 F
  row chr20 700 T G 8 PASS 0/1
  row chr21 100 A C 12.5 PASS 0/1 F
  row chr20 50 A C 0 PASS 0/1 F
} > "$IN/fa_edge.vcf"
{
  cat "$H"
  row chr20 50 A C 7 PASS 0/1
  row chr20 100 A G 30 PASS 1/1
  row chr20 150 A G 2.0 PASS 0/1
  row chr20 150 A G 2.0 PASS 0/1
  row chr20 160 A . 2.0 RefCall 0/0
  row chr20 250 A G 2.01 PASS 1/1
  row chr20 300 A G 40 PASS 0/1
  row chr20 350 A AG,G 40 PASS 1/2
  row chr20 600 T G 1 PASS 0/1
  row chr20 700 T G 8 PASS 0/1
  row chr20 650 T TA nan PASS 0/1
  row chr20 800 C C 9 PASS 0/0
  row chr21 120 A C 12.5 PASS 0/1
  row chr20 900 T G 5 PASS 1/0
} | gzip > "$IN/pileup_edge.vcf.gz"
# Text quirks: CRLF and CR line ends, a last line without a newline, a 7-column row whose
# newline LowQual drops, a blank pileup line (skipped under --ctgName) and spaces between
# full-alignment columns.
{
  printf '##fileformat=VCFv4.2\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\r\n'
  printf 'chr20\t10\t.\tA\tC\t1\tPASS\r\n'
  printf 'chr20 20 . A C 30 PASS P GT:GQ 0/1:3\n'
  printf 'chr20\t30\t.\tA\tC\t30\tPASS\tP\tGT:GQ\t1/1:3  \r'
  printf 'chr20\t40\t.\tA\tC\t30\tPASS\tP\tGT:GQ\t1/1:3'
} > "$IN/fa_text.vcf"
{
  printf '#header ignored\n'
  printf 'chr20\t5\t.\tA\tC\t1\tPASS\n'
  printf '\n'
  printf 'chr20\t15\t.\tA\tC\t30\tPASS\tP\tGT:GQ\t1/1:3 \t\n'
  printf 'chr20\t25\t.\tA\tC\t30\tPASS\tP\tGT:GQ\t0/1:3\n'
  printf 'chr20\t30\t.\tA\tC\t30\tPASS\tP\tGT:GQ\t0/1:3\n'
  printf 'chr20\t35\t.\tA\t.\t1\tRefCall\tP\tGT:GQ\t0/0:3'
} > "$IN/pileup_text.vcf"
# Error cases: a FA row whose columns are space-separated cannot be marked LowQual; a blank
# full-alignment line; a short row under a haploid mode; a bad position.
{ cat "$H"; printf 'chr20 20 . A C 1 PASS P GT:GQ 0/1:3\n'; } > "$IN/fa_spaces_lowqual.vcf"
{ cat "$H"; printf '\n'; } > "$IN/fa_blank.vcf"
{ cat "$H"; printf 'chr20\t20\t.\tA\t.\t1\tRefCall\n'; } > "$IN/fa_short.vcf"
{ cat "$H"; printf 'chr20\t2x\t.\tA\tC\t1\tPASS\n'; } > "$IN/pileup_badpos.vcf"

# Candidate BEDs for --platform ilmn. The oracle's (34-base windows) are all emptied by the
# 33-base padding, so Clair3 raises on them; beds_wide has windows wide enough to keep.
cp -r "$ORACLE/tmp/full_alignment_output/candidate_bed" "$IN/beds_oracle"
mkdir -p "$IN/beds_wide" "$IN/beds_none"
zcat "$ORACLE/full_alignment.vcf.gz" | grep -v '^#' | awk 'NR % 5 == 1 {print $1"\t"$2-40"\t"$2+60}' > "$IN/beds_wide/chr20.0_2"
printf 'chr20\t1000000\t1000066\nchr20\t2000000\t2000100\nchr21\t0\t10\n#c\nchr20\t3000000\t3500000\n' > "$IN/beds_wide/chr20.1_2"
# A file without a final newline fuses with the next one: either order gives the same regions.
mkdir -p "$IN/beds_fuse"
zcat "$ORACLE/full_alignment.vcf.gz" | grep -v '^#' | awk 'NR % 3 == 0 {print $1"\t"$2-50"\t"$2+50}' | head -c -1 > "$IN/beds_fuse/chr20.0_2"
printf '\tjunk\nchr20\t40000000\t50000000\n' > "$IN/beds_fuse/chr20.1_2"
printf 'chr20\t0\t64444167\n' > "$IN/beds_wide/chr2.0_1"
printf 'chr20\t0\t64444167\n' > "$IN/beds_wide/FULL_ALN_FILE_chr20"
printf 'chr20\t0\t64444167\n' > "$IN/beds_none/chr21.0_1"
mkdir -p "$IN/beds_bad"
printf 'chr20\t500\t400\n' > "$IN/beds_bad/chr20.0_1"
ls "$IN"

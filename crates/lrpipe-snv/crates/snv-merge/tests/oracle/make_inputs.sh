#!/bin/bash
# Builds the SortVcf parity cases under $WORK from the m3 oracle chunk VCFs: inputs/ and a copy
# of cases.tsv. Then run run_python.sh to make Clair3's outputs in $WORK/py, and the parity
# test with SNV_MERGE_CASES=$WORK.
set -eu
ORACLE=${ORACLE:-/data2/bioinformatics/clair3-rs-bench/m3-oracle}
WORK=${WORK:-/data2/bioinformatics/clair3-rs-bench/m3-merge-work}
REFERENCE=${REFERENCE:-/data2/bioinformatics/sniffles-rs-bench/reference.fa}
IN=$WORK/inputs
rm -rf "$IN"; mkdir -p "$IN"
# single contig: the oracle inputs
cp -r "$ORACLE/pileup_output" "$IN/single"
cp "$ORACLE/CONTIGS" "$IN/single.CONTIGS"
# multi contig: chunks 1-6 chr20, 7-10 chrX, 11-13 chr1_KI270706v1_random
mkdir -p "$IN/multi"
for i in $(seq 1 13); do
  src="$ORACLE/pileup_output/pileup_chr20_$i.vcf"
  if [ "$i" -le 6 ]; then cp "$src" "$IN/multi/pileup_chr20_$i.vcf";
  elif [ "$i" -le 10 ]; then sed 's/^chr20\t/chrX\t/' "$src" > "$IN/multi/pileup_chrX_$i.vcf";
  else sed 's/^chr20\t/chr1_KI270706v1_random\t/' "$src" > "$IN/multi/pileup_chr1_KI270706v1_random_$i.vcf"; fi
done
# a non-matching file and a full_alignment file that the prefix filter must skip
echo "junk" > "$IN/multi/other.vcf"
cp "$ORACLE/pileup_output/pileup_chr20_1.vcf" "$IN/multi/full_alignment_chr20_1.vcf"
printf 'chr1_KI270706v1_random\nchrX\nchr2\nchr20\nchr1\n' > "$IN/multi.CONTIGS"
# PAR bed for chrX (plain and gzipped), with a comment and a zero-length row
X0=$(grep -v '^#' "$IN/multi/pileup_chrX_7.vcf" | head -1 | cut -f2)
X1=$(grep -v '^#' "$IN/multi/pileup_chrX_8.vcf" | head -1 | cut -f2)
X2=$(grep -v '^#' "$IN/multi/pileup_chrX_9.vcf" | sed -n 100p | cut -f2)
printf '#PAR\nchrX\t0\t%s\nchrX\t%s\t%s\nchrY\t10\t20\n' "$X1" "$X2" "$X2" > "$IN/par.bed"
gzip -c "$IN/par.bed" > "$IN/par.bed.gz"
# empty cases
mkdir -p "$IN/noprefix"; echo x > "$IN/noprefix/foo.vcf"
mkdir -p "$IN/headers"; for i in 1 2; do grep '^#' "$ORACLE/pileup_output/pileup_chr20_$i.vcf" > "$IN/headers/pileup_chr20_$i.vcf"; done
mkdir -p "$IN/emptyfiles"; : > "$IN/emptyfiles/pileup_chr20_1.vcf"; : > "$IN/emptyfiles/pileup_chr20_2.vcf"
# malformed records
mkdir -p "$IN/malformed"
F=$ORACLE/pileup_output/pileup_chr20_1.vcf
{
  grep '^#' "$F" | grep -v '^#CHROM'
  # header line fused with the first record (get_header writes no newline after #CHROM)
  printf '%s' "$(grep '^#CHROM' "$F")"; grep -v '^#' "$F" | head -1
  # two records fused on one line
  grep -v '^#' "$F" | sed -n 2,3p | tr -d '\n'; echo
  # a record with >19 tab fields (three fused)
  grep -v '^#' "$F" | sed -n 4,6p | tr -d '\n'; echo
  # a short record (9 fields)
  grep -v '^#' "$F" | sed -n 7p | cut -f1-9
  grep -v '^#' "$F" | sed -n '8,200p'
} > "$IN/malformed/pileup_chr20_1.vcf"
# CRLF file, last line without newline
grep -v '^#' "$ORACLE/pileup_output/pileup_chr20_2.vcf" | head -50 > /tmp/m3crlf.$$
{ grep '^#' "$ORACLE/pileup_output/pileup_chr20_2.vcf"; cat /tmp/m3crlf.$$; } | sed 's/$/\r/' | head -c -2 > "$IN/malformed/pileup_chr20_2.vcf"
rm -f /tmp/m3crlf.$$
# a file whose header has a fused #CHROM line with a record of another contig
{ grep '^#' "$F" | grep -v '^#CHROM'; printf '%s' "$(grep '^#CHROM' "$F")"; grep -v '^#' "$ORACLE/pileup_output/pileup_chr20_3.vcf" | head -1 | sed 's/^chr20/chr21/'; grep -v '^#' "$ORACLE/pileup_output/pileup_chr20_3.vcf" | sed -n '2,40p'; } > "$IN/malformed/pileup_chr20_3.vcf"
printf 'chr20\n' > "$IN/malformed.CONTIGS"
# headers without ##cmdline (for the cmd_fn insertion)
mkdir -p "$IN/nocmd"; for i in 4 5; do grep -v '^##cmdline' "$ORACLE/pileup_output/pileup_chr20_$i.vcf" > "$IN/nocmd/pileup_chr20_$i.vcf"; done
# malformed records Clair3 survives (no fused header of another contig)
mkdir -p "$IN/malformed2"; cp "$IN/malformed/pileup_chr20_1.vcf" "$IN/malformed/pileup_chr20_2.vcf" "$IN/malformed2/"
cp "$ORACLE/CMD" "$IN/CMD"
echo "$REFERENCE" > "$IN/REF"
cp "$(dirname "$0")/cases.tsv" "$WORK/cases.tsv"
echo "inputs ready in $IN"

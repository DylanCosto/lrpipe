# Full E. coli K-12 small-variant benchmark

Both callers processed the complete ONT and HiFi archives from *E. coli* K-12 MG1655.
The runs used the lrpipe bacterial preset from PR #28 and Clair3 2.0.3 with the same model
files and the same Bonefish alignment for each sample and reference. Archive sizes and
MD5 hashes matched ENA, and both FASTQs passed gzip validation. `main` at `91b8ac0` contains
the tested pipeline source after the merge of PR #28.
[ONT archive](https://www.ebi.ac.uk/ena/browser/view/SRR30993121),
[HiFi archive](https://www.ncbi.nlm.nih.gov/sra/SRX7636524),
lrpipe source revision

| Platform | Full archive | Reads | Read bases | Mapped primary reads | Current reference covered | Mean depth |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ONT | `SRR30993121` | 50,430 | 508,242,655 | 50,430 | 99.25% | 105.04× |
| HiFi | `SRR10971019` | 95,514 | 1,389,500,381 | 95,514 | >99.99% | 295.28× |

## Precision and recall against assembly truth

The independently curated current MG1655 assembly `NC_000913.3` supplied the expected
alleles. Both full FASTQs were aligned once to each of two scoring references. The historical
MG1655 reference `U00096.1` supplies many known sequence corrections across both variant
classes. The BW25113 reference `CP009273.1` supplies differences between two laboratory
K-12 genomes. [MG1655 reference](https://www.ncbi.nlm.nih.gov/bioproject/57779),
[MG1655 resequencing study](https://pmc.ncbi.nlm.nih.gov/articles/PMC1681481/),
[BW25113 finished assembly](https://journals.asm.org/doi/10.1128/genomea.01038-14)

| Scoring reference | Scored bases | Truth SNVs | Truth small indels |
| --- | ---: | ---: | ---: |
| `U00096.1` | 4,638,164 | 144 | 137 |
| `CP009273.1` | 4,627,594 | 18 | 2 |

| Scoring reference | Reads | Caller | SNV precision | SNV recall | Indel precision | Indel recall |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| MG1655 historical | ONT | lrpipe | 94.74% | 95.14% | 95.86% | 96.35% |
| MG1655 historical | ONT | Clair3 | 94.74% | 95.14% | 95.86% | 96.35% |
| MG1655 historical | HiFi | lrpipe | 99.22% | 93.75% | 95.95% | 99.27% |
| MG1655 historical | HiFi | Clair3 | 99.22% | 93.75% | 95.95% | 99.27% |
| BW25113 | ONT | lrpipe | 78.95% | 83.33% | — | 0.00% |
| BW25113 | ONT | Clair3 | 78.95% | 83.33% | — | 0.00% |
| BW25113 | HiFi | lrpipe | 100.00% | 83.33% | 100.00% | 100.00% |
| BW25113 | HiFi | Clair3 | 100.00% | 83.33% | 100.00% | 100.00% |

The callers reported identical normalized PASS alleles within both scoring masks on both
platforms. On BW25113 with ONT, the mask excluded 174 lrpipe calls and 162 Clair3 calls.

Precision is undefined when a caller reports zero eligible variants in a class. The natural
BW25113 comparison has two small-indel truth events. Both callers missed those events on ONT,
where at least ten primary alignments covered each variant anchor at base quality 10 and mapping
quality 20 or higher. The scores measure agreement with the published MG1655 assembly. The ONT
dataset has a distinct BioSample and one public sequencing run. Stock-specific alleles can count
as disagreements.
[ONT BioProject](https://www.ncbi.nlm.nih.gov/bioproject/1172918)

## CPU runtime and current-reference agreement

Wall times are seconds for one fresh run per caller and input, using eight threads on an Intel
Xeon Gold 5418Y VM. The Bonefish time is shared. The lrpipe calling time is its SNV stage, and the
Clair3 calling time starts from the same BAM and ends with its indexed VCF. lrpipe used its default
CPU inference settings for `run --bacterial`. Both used haploid-precise calling with no phasing or gVCF.

| Reference | Reads | Shared alignment | lrpipe calling | Clair3 calling | lrpipe FASTQ to QC |
| --- | --- | ---: | ---: | ---: | ---: |
| MG1655 current | ONT | 8.18 | 7.74 | 20.50 | 16.79 |
| MG1655 current | HiFi | 19.24 | 24.27 | 86.28 | 43.77 |
| MG1655 historical | ONT | 9.02 | 8.69 | 25.41 | 18.65 |
| MG1655 historical | HiFi | 29.55 | 25.06 | 81.67 | 54.94 |
| BW25113 | ONT | 7.74 | 8.05 | 21.24 | 16.75 |
| BW25113 | HiFi | 26.36 | 22.66 | 75.51 | 49.37 |

On `NC_000913.3`, the complete reads also provided a comparison of normalized PASS calls
between the two tools. The numerator is shared alleles and the denominator is their union.

| Reads | Shared / union SNVs | Shared / union small indels |
| --- | ---: | ---: |
| ONT | 5 / 5 | 5 / 5 |
| HiFi | 1 / 1 | 3 / 3 |

## Scoring and reproduction

1. [Build the truth VCFs](build_truth.py) from the versioned FASTAs with minimap2 and paftools.
   The mask keeps reciprocal primary alignments of at least 5 kb with MAPQ 60. It trims 100 bp
   from alignment boundaries and excludes 100 bp on either side of assembly gaps longer than
   49 bp. The mask depends only on the assemblies. Truth and query alleles must fit entirely
   in the mask after normalization. Both truth VCFs scored perfectly against themselves.
2. [Run and time the callers](run_callers.py) with Bonefish 0.5.0, Clair3 2.0.3 and the lrpipe
   binary at the recorded revision. One CPU run of `lrpipe run --bacterial` produces the BAM
   and the lrpipe VCF. Clair3 takes that BAM with `--include_all_ctgs`, `--haploid_precise` and
   `--no_phasing_for_fa`. The ONT model is `r1041_e82_400bps_sup_v500`, and the HiFi model is
   `hifi_sequel2`. [Clair3 source](https://github.com/HKU-BAL/Clair3)
3. [Score the outputs](score_callers.py) with bcftools 1.24 and RTG Tools 3.13. Normalize and
   atomize PASS variants and require haploid alternate genotypes. Report SNVs and 1–49 bp
   indels. RTG uses haplotype matching at ploidy 1 with `--decompose`.
   Precision divides the call-side true positives by the eligible calls. Recall divides the
   truth-side true positives by the truth variants. Equivalent haplotypes can contain different
   numbers or classes of VCF records.

The [machine-readable metrics](metrics.json) contain archive and model hashes, scored territory,
all TP, FP and FN counts, calls outside the scoring masks, exact caller agreement, coverage and
runtime receipts. Prepare a scratch directory with `inputs/<accession>.fastq.gz`, the ENA `read_run`
metadata in `inputs/<accession>.tsv` and a checksum in `inputs/<accession>.sha256`. The TSV needs
`read_count` and `base_count`. Name the indexed
references `inputs/mg1655.fa`, `inputs/U00096.1.fa` and `inputs/CP009273.1.fa`. Run each script with
`--help` for the path options. Compress and index each `truth.vcf` and create `<output>/reference.sdf`
with `rtg format` before scoring its reference. The report uses `mg1655`, `U00096.1` and `CP009273.1`
as the reference labels and `ont` and `hifi` as the platform names.

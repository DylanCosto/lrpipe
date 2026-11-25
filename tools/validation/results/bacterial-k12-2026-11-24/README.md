# Bacterial preset smoke validation

`lrpipe run --bacterial` ran from FASTQ through Bonefish alignment, haploid small-variant
calling, final merge and QC on *Escherichia coli* K-12 MG1655. The preset uses the `snv`
analysis for SNVs and small indels. Both inputs were the first complete reads in archive order
through 125 million sequence bases, using the matching Clair3 models and the same 4,641,652 bp
reference. [MG1655 RefSeq record](https://www.ncbi.nlm.nih.gov/bioproject/57779)

| FASTQ | Caller | Primary reads mapped | Reference covered | PASS SNVs | PASS small indels | gVCF blocks | Wall time |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ONT SRR30993121 | CPU | 12,095 | 25.46% | 2 | 0 | 7,041 | 7.0 s |
| HiFi SRR10971019 | CPU | 8,580 | >99.99% | 1 | 2 | 78,547 | 8.6 s |
| ONT SRR30993121 | CUDA | 12,095 | 25.46% | 2 | 0 | 7,041 | 9.9 s |

ONT CPU and CUDA agreed on both PASS positions, alleles and genotypes. Wall times are
single runs from FASTQ through QC. Primary alignments covered 1,181,632 of
4,641,652 reference bases on ONT and 4,641,651 bases on HiFi. These are workflow
and output-integrity results. Precision and recall require an independent bacterial truth set.
[ONT archive record](https://www.ebi.ac.uk/ena/browser/view/SRR30993121), [HiFi archive record](https://www.ncbi.nlm.nih.gov/sra/SRX7636524)

## Output checks

- Alignment passed `samtools quickcheck`. The manifest recorded ploidy 1 and sex `not-applicable`.
  The calling stage skipped human sex inference and phasing.
- The sample's merged VCF, the caller VCF and the gVCF agreed on PASS positions, alleles and
  genotypes. `bcftools norm --check-ref e` verified REF alleles against the FASTA. Indexed reads
  with pysam checked the sample name, record bounds, one allele per genotype, one likelihood per
  allele whenever PL was present, and nonoverlapping gVCF reference blocks.
- Repeating each command skipped all completed stages and preserved all three compressed output
  hashes. Switching an existing ONT run from human contig selection to the bacterial preset
  reran calling, merge and QC and recovered both whole-reference PASS calls.
- A 400 bp ONT BED target had no full-alignment candidates. Both the VCF-only and gVCF runs kept
  the expected haploid SNV in an indexed VCF. The gVCF run also emitted indexed reference blocks.
- A second ONT run used two 20,001 bp fragments of the same K-12 reference as separate
  contigs. Each fragment had its expected haploid PASS SNV in both the indexed VCF and gVCF.
  This checks short-contig calling with the real aligner and caller.

## Regression checks

All 9 checks in `bash tools/check.sh --release` passed locally, including the Clair3
fixture, the pileup fallback and CPU gVCF overlap regression, WhatsHap comparison and
Bonefish alignment. The workspace passed 721 tests on Rust 1.98 and 721 on
the Rust 1.91 toolchain. The CPU-only inference tests and mixed-sample alias regression
also passed.

## Reproducing the runs

The site configuration pointed to the indexed MG1655 FASTA, Bonefish 0.5.0 and the model directory.
`[snv] device` was `cpu` for the CPU runs and `cuda:0` for the ONT CUDA run. The VM used an Intel
Xeon Gold 5418Y CPU and an NVIDIA L40S-24Q. Every command used eight threads and a 32 GiB memory budget.

```bash
lrpipe run --bacterial --config k12.toml --reads SRR30993121.smoke.fastq.gz \
  --sample K12-ont --platform ont --model r1041_e82_400bps_sup_v500 \
  --output k12-ont --threads 8 --memory 32 --gvcf

lrpipe run --bacterial --config k12.toml --reads SRR10971019.smoke.fastq \
  --sample K12-hifi --platform hifi --model hifi_sequel2 \
  --output k12-hifi --threads 8 --memory 32 --gvcf
```

| Model | `pileup.pt` SHA-256 | `full_alignment.pt` SHA-256 |
| --- | --- | --- |
| `r1041_e82_400bps_sup_v500` | `d84d9f1a1da0a68204404cac1817724ada4d15363d8a62d87703b6ed56c4a4ff` | `320f7083ae1b8eeba390471fbbeb6fbefb93dab079300ce2605a7222d7d65388` |
| `hifi_sequel2` | `71a762cb98e17f163496174bcc35b2ff5606db760383d91c9a885d6637fdc6a0` | `9ae2f5c75184f40649dde5c5bc82d39a65c3bfa5738f516fe20353008f0b4c26` |

Raw reads are available from the public archives. The [machine-readable results](metrics.json)
include FASTQ and reference hashes, stage timings, per-output counts and the tested commit.

The [full-read benchmark](../bacterial-k12-full-2026-11-24/README.md) compares lrpipe with
Clair3 on both complete archives against two assembly-based truth panels.

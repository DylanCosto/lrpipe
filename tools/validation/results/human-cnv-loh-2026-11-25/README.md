# Human CNV and copy-neutral LOH validation

25 November 2026 · lrpipe 0.2.0 · GRCh38 autosomes

Both independently supported HG002 one-copy deletions were recovered on HiFi and ONT. Exact tumor CNV dosage covered 25.1% of confirmed changed bases in HG008-T and 19.4% in HCC1395-T. Copy-neutral LOH performance varied sharply with the sample. These results support keeping somatic CNV and copy-neutral LOH experimental.

[Per-event results](events.tsv) contain 241 scored tumor event groups. [Runtimes](runtimes.tsv) include the three whole-genome runs, six CNV reruns and two HCC1395 diagnostics. [Metrics](metrics.json) record the denominators, call IDs, control results and commands.
In the event table, `truth_minor` is `-1` for HCC1395 dosage groups with an unassigned minor copy count. `NA` marks an unmatched call or an inapplicable field.

## Results

| Dataset | Endpoint | Assessed span | Event hits / truth | Event precision | Base precision | Base recall | Base F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| HG008-T ONT | Exact CNV dosage | 2285.60 Mb | 13/44 | 1.0% | 98.2% | 25.1% | 40.0% |
| HG008-T ONT | CN2 LOH | 2281.95 Mb | 2/3 | 16.7% | 87.4% | 82.4% | 84.9% |
| HCC1395-T HiFi | Exact CNV dosage | 1847.33 Mb | 21/120 | 14.3% | 69.7% | 19.4% | 30.4% |
| HCC1395-T HiFi | CN2 LOH | 1847.33 Mb | 0/74 | 0.0% | 49.0% | 27.7% | 35.4% |

Precision and recall apply only to the assessed masks. Event hits require a single call with reciprocal overlap. Base scores allow recovery by multiple fragments and expose partial matches. In HG008-T, 178 CN1 calls overlapped the independently confirmed chromosome 4 deletion `CNA_8`, covering 45.5 of 189.9 assessed Mb. The largest individual overlap was 1.59 Mb. The strict event score penalizes the fragmentation as well as the missed bases.

### Copy-neutral LOH

| Dataset | Confirmed CN2 LOH (Mb) | Any LOH over truth (Mb) | LOH labelled CN2 over truth (Mb) | CN2 LOH without a conflicting CNV (Mb) |
|---|---:|---:|---:|---:|
| HG008-T | 74.63 | 61.52 | 61.52 | 60.84 |
| HCC1395-T | 863.89 | 836.20 | 239.68 | 237.99 |
| HCC1395BL control | 863.89 | 1.00 | 1.00 | 1.00 |

The matched-normal row counts calls overlapping the tumor truth coordinates.

| HG008-T event | Chromosome | Independently assessed (Mb) | LOH labelled CN2 (Mb) | CN2 LOH without a conflicting CNV (Mb) |
|---|---|---:|---:|---:|
| CNA_1 | chr1 | 23.26 | 21.80 | 21.59 |
| CNA_21 | chr17 | 43.35 | 31.71 | 31.24 |
| CNA_32 | chr19 | 8.02 | 8.01 | 8.01 |

HCC1395-T detected LOH across 96.8% of confirmed copy-neutral bases and assigned CN2 to 27.7% of confirmed bases. It reported LOH with CN1 over the independently confirmed CN2 LOH in distal chr4, distal chr13 and proximal chr17q. The chr17 LOH call spanned roughly 22–83 Mb. Independent SNP arrays place a CN2 LOH block at about 27–67 Mb and a gain farther down the arm. The LOH caller estimates one copy state for the whole run of homozygous alternate calls. [SEQC2 orthogonal assays](https://zenodo.org/records/10081574), [SEQC2 clonal analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/).

| HCC1395 locus | Matched-normal heterozygous sites with depth | Normal median minor-allele fraction | Tumor median minor-allele fraction |
|---|---:|---:|---:|
| Distal chr4 | 192 | 43.4% | 0.0% |
| Distal chr13 | 203 | 43.5% | 0.0% |
| Proximal chr17q | 160 | 44.3% | 0.0% |

This read audit counted reference and alternate bases in both HiFi BAMs at dispersed, high-quality SNPs called heterozygous in the normal. It retained sites with at least 10 informative normal reads, at least 22% normal minor allele fraction and at least 8 informative tumor reads. The independent array calls define the truth in the table above.

### Dosage and normalization

| Dataset | CN0 agreement | CN1 agreement | CN3 agreement | CN4 agreement | Gain/loss direction recall | Estimated CN2 depth |
|---|---:|---:|---:|---:|---:|---:|
| HG008-T ONT | — | 22.7% | 96.3% | 79.3% | 25.2% | 48× |
| HCC1395-T HiFi | 89.8% | 84.6% | 10.2% | 47.1% | 26.6% | 31× |

The HCC1395 run normalized to 31×. The median 1 kb MAPQ20 depth over array-confirmed CN2 LOH was 20.99×. Adding `--cancer` gave identical copy-neutral LOH base scores and changed exact CNV dosage recall from 19.4% to 4.7%. A separate, truth-informed dosage-only run set normalization to 20.99×. Gain/loss direction recall moved from 26.6% to 66.2%, while exact CNV dosage precision moved from 69.7% to 44.8%. The override applies in the depth-only path. This diagnostic uses a known-copy anchor from the same sample and has no LOH output.

### Germline controls

| HG002 interval | Independent Illumina target/flank median depth | HiFi predicted CN | ONT predicted CN |
|---|---:|---:|---:|
| chr4:68,516,025–68,627,294 | 50.8% | 1 | 1 |
| chr22:22,759,209–22,904,986 | 54.1% | 1 | 1 |

HG002 had 8 HiFi and 9 ONT autosomal CN2 LOH calls. Seven of those matched across platforms at 50% reciprocal overlap. The Q100 benchmark corroborates long homozygous stretches with remaining heterozygous sites in several calls. These calls illustrate that the single-sample label can also describe regions of homozygosity. The germline dosage results cover two deletions. The other GIAB samples provide reproducibility controls.

### Runtime and storage

| Dataset | GPU SNV | CPU read pass | CPU CNV | Pipeline wall time |
|---|---:|---:|---:|---:|
| HCC1395-T HiFi | 318.9 s | 105.0 s | 73.7 s | 397.4 s |
| HCC1395BL HiFi | 228.0 s | 100.5 s | 57.8 s | 288.7 s |
| HG008-T ONT | 574.2 s | 359.0 s | 59.2 s | 636.5 s |

The CPU read pass overlapped GPU SNV calling. Pipeline wall time includes QC. The standalone CNV rows in [runtimes.tsv](runtimes.tsv) start from existing coverage and, where applicable, SNV files. `NA` marks a stage outside the timing scope. The three new BAMs total 216.2 GiB, with 539 GiB free on the staging volume at completion.

## Inputs and scoring

The test used lrpipe 0.2.0, SHA-256 `624957c41d8a259de4f1e286c09f4764ff72e6e709b77f4f6690a9359b4229b0`. The executable was built from `72a248a1c4f069d01f08670731b8282cfea23224`. The inspected main revision, `e171820c5f6a7930a331f5f58901ac2e1e9335b1`, has the same caller and pipeline code. The three new runs started from indexed public whole-genome alignments, with `lrpipe run --analyses cnv --haplotag off --merged-vcf off`, GPU SNV inference and the pipeline's CPU CNV stage with `--improved`, the GRCh38 blacklist and SV breakpoint refinement. The six GIAB reruns reused saved depth, SNV and SV files and reproduced the previous autosomal records. Each new and rerun CNV VCF passed strict REF and index checks.

- **HG008-T.** Northeastern ONT R10.4.1 ultra-long, vendor-labelled 54×, batch 0823p23. The local reference matches the provider's GRCh38-GIABv3 BAM dictionary. Truth starts with the NIST GRCh38 v0.6 PASS copy-number BED and its two haplotype copy counts. Each tested base lies in a concordant segment from an independent BCM Illumina tumor/normal analysis of the same batch. The Illumina segment needed at least 10 depth windows and total copy within 0.35 of the curated integer. Allele scoring also required at least 20 informative BAF sites and minor allele copy within 0.25. An event needed corroboration over at least half of its curated bases and enough corroborated sequence to meet the caller's size floor. [NIST v0.6 draft](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.6-20261006/), [BCM Illumina analysis](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Hartwig_Oncoanalyzer_BCM-ILMN_20260204/), [ONT alignment input](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern-ONT-UL-20241216/).
- **HCC1395 and HCC1395BL.** Public DeepSomatic PacBio Revio HiFi 30× titrations, run individually on GRCh38. The tumor mask intersects the published SEQC2 CNV benchmark with calls of the same copy state in at least two independent tumor SNP-array replicates. It excludes changes called in two or more matched-normal arrays. A copy-neutral LOH base requires array CN2 with LOH and a compatible LOH interval in the published benchmark without a benchmark gain or loss at that base. The array copy counts are independent estimates. The benchmark VCF's median copy numbers were excluded from exact-copy scoring. The same tumor intervals were used to measure overlap by any CN2 LOH calls in the separately called matched normal. Genome-wide precision cannot be inferred from this selected mask. [SEQC2 benchmark paper](https://doi.org/10.1186/s13059-024-03294-8), [published intervals](https://zenodo.org/records/10081398), [orthogonal assay calls](https://zenodo.org/records/10081574), [SNP array sample mapping](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212051), [HiFi data study](https://doi.org/10.1038/s41587-025-02839-x).
- **GIAB HG002, HG003 and HG004.** Reused whole-genome HiFi and ONT runs. For the HG002 dosage check, start from simple autosomal Q100 deletions at least 100 kb long. The independent local PCR-free Illumina alignments needed MAPQ20 target depth between 35% and 65% of the mean of the two 100 kb flank medians, similar flank depths and at least 80% of the MAPQ0 target median retained at MAPQ20. The three Q100 candidates had varying official BED coverage, so the two corroborated deletions are individual spot checks. The remaining chromosome 16 candidate is unscored in the two-event result because independent short-read depth was mapping-sensitive. [GIAB benchmark reference](https://www.nist.gov/programs-projects/faqs-genome-bottle), [HG002 Illumina dataset](https://opendata.nist.gov/pdrsrv/mds2-2336/).

Primary scores use GRCh38 autosomes. Coordinates in the tables are zero-based and end-exclusive. The event filters were at least 100 kb of supported dosage-changing sequence and more than 1 Mb of supported copy-neutral LOH. A matched event needed at least 50% overlap in both directions within the assessed mask and at least half of the prediction itself inside that mask. One prediction and one truth event can each participate in one match. The base-level measures distinguish any LOH call, an LOH call labelled CN2 and a CN2 LOH call with concordant dosage output. Absence of a dosage call was scored as CN2. Excluded bases were omitted from both precision and recall. An event score can be lower than the base-level score when fragmented truth or overextended calls fail reciprocal overlap.

The BAM and index downloads were verified against provider digests. Read depth and SV results were generated from those same alignments in lrpipe. CNV stage duration and wall time are measured from local inputs. The [metrics file](metrics.json) includes the public alignment and index URLs. Its command paths use `$LRPIPE` for the executable, `$RUN_ROOT` for the run root, `$GRCH38` for the indexed NCBI no-alt reference, `$CLAIR3_MODELS` for the model directory, and `$CNV_BLACKLIST` for the GRCh38 Spectre blacklist. Each study configuration supplies its reference, the GPU SNV device, models and blacklist.

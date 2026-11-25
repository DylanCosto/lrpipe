# lrpipe 0.2.1 option and resume validation

> **Status.** This records lrpipe 0.2.1 as released. Since then: the stale repeat index
> (#42, fixed in #50), damaged outputs accepted on resume (#46, #49), the standalone phasing
> `.fai` check (#43, #47), late XY repeat failures with unknown PARs (#44, #48), CUDA calls at
> high depth (#45, #51) and bonefish's MM/ML tags ([bonefish#46](https://github.com/DylanCosto/bonefish/issues/46)) are fixed on `main`, and so are the SV `.vcf.bgz` index
> (#52) and `--pileup-only --gvcf`
> (#54), both in 0.2.2, and gVCF blocks between
> disjoint targets (#53), after 0.2.2. Code links point at the 0.2.1 commit.

25 November 2026 · published Linux x86-64 release · `ce4486ff8857f9d03edb86ceeaad3eafa7a2f2af`

Resuming a repeat run can leave an earlier tabix index on the current VCF. In the panel-to-whole-reference reproduction, indexed site queries saw 1 of 10 repeat records. Default resume also reported success with unreadable results in deliberately damaged test directories. A forced rerun recovered all 108 checkpoint failures in the tested cases.

The follow-up covers run modes, GRCh38 genotype ploidy, scheduling, standalone SNV and SV output options, and gVCF coverage across disjoint targets. [Follow-up validation](run-modes.md).

## Coverage

The [coverage table](coverage.tsv) records the outcomes for each matrix. The [checkpoint trials](checkpoint-recovery.tsv) record the affected artifact, injected fault and recovery outcome. The [repeat index checks](repeat-index.tsv) record how many sites were visible before and after rebuilding the stale index. Per-input results cover the [output options](output-options-by-input.tsv), [resume transitions](resume-by-input.tsv) and [CUDA comparisons](cuda-parity-by-input.tsv).

The earlier read-layout matrix exercised every analysis subset with 12 layouts and four region states. All 7,776 valid runs completed and all 4,608 invalid requests were refused. Across equivalent inputs, 56,704 record stream comparisons agreed.

The expanded matrix checks each of the 127 nonempty analysis subsets, the default selection and explicit `all`. It crosses those selections with HiFi and ONT, 12 read layouts, four region states, and all combinations of `--merged-vcf`, `--gvcf`, `--haplotag` and `--haplotagged-reads`. All 198,144 invocations met their expected outcome. There were 78,336 valid plans and 119,808 expected refusals. The refusals were checked as actual invocations and left no run directory.

The 78,336 valid plans produced 18,624 distinct child command sets after substituting a placeholder for the output path. Cases with the same command set agreed on dependencies and required outputs. Of these sets, 5,772 match previously verified runs. All 12,852 new executions passed their stage and output checks. Raw reads were aligned once for each fixed fixture. The planning and runtime matrices reused the stamped alignment while checking every distinct downstream command set.

The test reference uses shifted GRCh38 chr20 fragments from HG002 and HG003 with small synthetic supporting contigs. The bacterial controls use K-12 alignments with about 45× coverage in a 15 kb interval. The 328 CUDA runs completed on an NVIDIA L40S. Across 1,736 comparisons of CUDA outputs with equivalent layouts or gVCF settings, the only differences were genotype-likelihood fields and formatting introduced by gVCF emission. The 96 controlled interruption trials used SIGINT and SIGTERM on a terminal or pipe. An independent invocation was refused while the interrupted run held its lock, and the resumed record streams matched clean runs. PR #39 records the macOS and Metal smoke results separately.

The 310 sample-label trials included 118 accepted runs and 192 expected refusals. The split-raw matrix included 256 accepted runs and 32 expected refusals across file order, aligner index reuse, gVCF and analysis settings. It compared alignment records independent of coordinate ties and haplotag assignments by read.

## Resume findings

### Stale tabix index on a current repeat VCF

With haplotagging on, the separate repeat stage writes a VCF and its tabix index. A later run in the same directory can use the shared read pass, which writes the same VCF path and leaves the prior index in place. On the tested panel-to-whole transition, the compressed VCF contained 10 correct rows on each of HiFi and ONT, and its index was byte-identical to the earlier one-record panel index. Individual locus lookups found 1 row. A whole-contig lookup found all 10. A newly built index on a copy restored all 10 individual lookups. Separate VCF writer, shared read-pass writer.

The full resume matrix exercised 9,856 transitions between 22 profiles across both technologies and all 12 read layouts. It recorded 7,648 clean passes, 1,488 expected refusals and 720 transitions with a stale repeat index. In 192 of the affected transitions, indexed queries missed nine records each. The resumed VCF record streams matched their clean controls. The [breakdown](resume-by-input.tsv) shows the incidence for each read layout and platform.

Adding `--rerun all` to the whole-reference invocation still left the earlier index on both platforms. `bcftools index --tbi --force "$repeat_vcf"` restored all 10 indexed site lookups in each test.

To exercise the transition, use a reference and repeat catalog with loci inside and outside the panel and an indexed BAM from the same reference. The output path is an isolated scratch directory.

```bash
run=("$LRPIPE" run --config "$SITE" --platform "$PLATFORM" --sample "$SAMPLE" \
    --reads "$BAM" --output "$RUN_ROOT" --analyses phase,tr)
"${run[@]}" --regions "$PANEL_BED"
repeat_vcf="$RUN_ROOT/repeats/$SAMPLE.repeats.vcf.gz"
cp "$repeat_vcf.tbi" "$RUN_ROOT/panel.tbi"
"${run[@]}" --haplotag off
cmp "$RUN_ROOT/panel.tbi" "$repeat_vcf.tbi"
bcftools query -f '%CHROM\t%POS\n' "$repeat_vcf" |
  while read -r contig pos; do
    printf '%s\t%s\t' "$contig" "$pos"
    bcftools view -H -r "$contig:$pos-$pos" "$repeat_vcf" | wc -l
  done
```

### Damaged output accepted or carried into a failed resume

In 424 recovery trials across HiFi and ONT BAM, CRAM, tagged FASTQ and BAM plus CRAM, the default rerun repaired 316 test copies. It completed with a damaged output in 56 and exited with an error in 52. All 108 cases recovered with `--rerun all`.

For example, toggling one bit in the stored CRC32 of the first BGZF block preserved the file length and left the paired index intact. The default rerun skipped the stage and reported a complete run for damaged final VCFs, gVCFs and haplotagged BAMs. It also accepted a truncated or invalid QC JSON and a missing or damaged repeat tabix index in the tested configurations. The checkpoint for the separate repeat stage lists the VCF itself. The generic checkpoint uses declared index contents and paired data length for indexed outputs. Completion-stamp checks, repeat-stage declaration.

To reproduce the final-VCF check on a completed test run, corrupt only the stored checksum and rerun the same command with default resume and then with `--rerun all`. Run `gzip -t "$MERGED_VCF"` after each invocation.

```bash
python3 - "$MERGED_VCF" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
data = bytearray(path.read_bytes())
assert data[:4] == b"\x1f\x8b\x08\x04" and data[12:14] == b"BC"
block_size = int.from_bytes(data[16:18], "little") + 1
data[block_size - 8] ^= 1
path.write_bytes(data)
PY
```

## Other release findings

### CPU and CUDA variant call differences

The default-analysis CUDA matrix covered all four output options across 12 read layouts, both platforms and whole-reference or panel calling. All 576 valid runs completed on CUDA and 192 requests produced the expected refusal. Matched CPU and CUDA runs differed in 480 of 1,440 SNV, phased SNV and gVCF record streams, across 192 of 576 pairs. The differences arose at three loci. The comparison checked coordinates, alleles, genotypes and phasing, including LowQual records.

- **HG002 HiFi insertion.** On the same BAM, CPU emitted `C>CTTTATTAA` at fixture `chr20:521964` as a `0/1` PASS call with QUAL 6.39 in three repeats. CUDA returned `0/0` at the full-alignment stage and the final variant VCF omitted the site. A run on each device with `--approximate` removed reproduced the difference. The corresponding GRCh38 site is `chr20:30191964`. The Q100 v5.0q VCF has a filtered record for this allele with a missing genotype, outside its benchmark BED. [Benchmark data sources](../../README.md).
- **HG003 ONT SNP.** On the same BAM, the default approximate CPU run retained a `0/1` `C>T` pileup call with QUAL 18.51 at fixture `chr20:545726`. CUDA rescored the site as `0/0` in full alignment. Both backends emitted the same pileup call with `--var-pct-full 0.5` and both omitted it with `--var-pct-full 0.7`. The published settings use 0.5 on approximate CPU and 0.7 on CUDA for this model. The GRCh38 site `chr20:30215726` lies outside HG003's GIAB v4.2.1 benchmark BED. ONT candidate selection settings.
- **HG003 ONT insertion.** On the same BAM produced from FASTQ, CUDA emitted a `0/1` LowQual insertion with QUAL 0.98 at fixture `chr20:588807`. CPU called `0/0` in full alignment. The difference persisted with equal `--var-pct-full` values and with `--approximate` removed. The GRCh38 site `chr20:30258807` lies outside HG003's GIAB v4.2.1 benchmark BED.

### Findings reproduced in earlier release checks

- **Bonefish modification tags.** Bundled Bonefish produced 31 malformed nonprimary tag sets per HiFi fixture format and 11 per ONT fixture format. Tagged FASTQ and unaligned BAM reproduced the defect. A local aligner patch cleared the invalid tags in all four tests and preserved the alignment records and primary modification positions and probabilities.
- **Standalone phasing FASTA index.** With a mismatched `.fai`, the published standalone phasing command exited successfully and changed phasing on 22 of 66 records. A local validator patch refused the bad index and preserved the 66 valid-index records. The separate pipeline and CNV preflight work is tracked in PR #39.
- **High-depth K-12 HiFi CUDA.** At 295× on the same alignment, CUDA SNP and indel F1 were 42.51% and 1.44%. CPU SNP and indel F1 were 96.41% and 97.58%. At a fixed 59× downsample, CUDA F1 was 97.86% for SNPs and 98.61% for indels. The compact 45× layout smoke tests measure execution and output consistency.
- **Region behavior.** A custom XY reference with unknown pseudoautosomal regions passed the dry run and failed at tandem-repeat calling. In a separate panel check, all 64 human and 80 K-12 gVCFs included reference blocks extending 33 bases beyond the BED edges. The matched ONT Clair3 run returned the same four boundary blocks.

The executable SHA-256 was `dca4f4d049296ec2f2cb7bfe00d1269cb0c8b84718bdc52498d4ab4c827afe98`. [Published release](https://github.com/DylanCosto/lrpipe/releases/tag/v0.2.1)

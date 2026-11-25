# The pipeline: `lrpipe run`

One command from a sample's long reads (PacBio HiFi or Oxford Nanopore), aligned (one sorted,
indexed BAM or CRAM, or [several](#several-read-files)) or not (dorado's or PacBio's unaligned
BAM, or FASTQ, which [bonefish](#unaligned-reads) aligns first), to:

| Result | Tool |
|---|---|
| SNVs and indels, phased | `lrpipe snv` ([lrpipe-snv](../lrpipe-snv): Clair3 in Rust), phased by `lrpipe phase` ([lrpipe-phase](../lrpipe-phase): WhatsHap in Rust) |
| The reads' haplotypes from the phased SNVs (a haplotag list; the haplotagged BAM on request) | `lrpipe haplotag` ([lrpipe-phase](../lrpipe-phase): WhatsHap's `haplotag`) |
| SVs, phased like the SNVs, and read depth from one pass over the reads | `lrpipe sv` ([lrpipe-sv](../lrpipe-sv): Sniffles2 in Rust, with [lrpipe-depth](../lrpipe-depth) built in) |
| Tandem repeat genotypes, their alleles split by the run's haplotypes (without phasing, in the read pass) | `lrpipe tr` ([lrpipe-tr](../lrpipe-tr)) |
| CNVs and loss of heterozygosity | `lrpipe cnv` ([lrpipe-cnv](../lrpipe-cnv): Spectre in Rust) |
| CpG methylation, combined and per haplotype | `lrpipe methyl` ([lrpipe-methyl](../lrpipe-methyl): modkit's and pb-CpG-tools' tables) |

and the variant calls (SNVs and indels, SVs, CNVs and tandem repeats) together as one sorted,
indexed VCF, `SAMPLE.vcf.gz`, with a manifest recording every input, tool binary (version and
SHA-256), command, outcome, time and peak memory, and a QC summary of the run ([QC
summary](#qc-summary)).

Each stage runs its tool as a program: one of lrpipe's own subcommands, with the binary running
the pipeline unless the site configuration names another.

It calls one germline sample: somatic, trio-aware calling, annotation, pharmacogenomics, HLA
typing and paralogous genes need other tools ([Not covered](../../README.md#not-covered)).

## Running

```bash
lrpipe run --config site.toml --reads HG002.bam --sample HG002 --platform hifi --output hg002
```

* `--reads` takes a sorted, indexed BAM or CRAM aligned to the site's reference, or unaligned
  reads ([Unaligned reads](#unaligned-reads)); given once for each file, several aligned files
  are merged first ([Several read files](#several-read-files)). The site's reference may be
  plain or bgzipped (`.fa.gz` with its `.fai` and `.gzi`).
* Before planning, lrpipe checks what a stage would otherwise meet first, in its log, or not at
  all, for every `--reads` (each of several; unaligned reads, whose format is checked at once,
  have their aligned BAM's contigs and index checked once the align stage has written it): the
  sample name (up to 200 letters, digits, `.`, `-` and `_`, not starting with `.` or `-`: it
  names the outputs' files and fills the VCFs' sample column); the reads (that each file exists,
  is BAM or CRAM as its name says, not SAM, is not sorted by read name (`SO:unsorted` with an
  index is warned of and taken as sorted), and has its index beside it
  as `READS.bai`, `READS.csi` or `READS.crai`, of as many contigs as its header); the reference
  and its `.fai` (and `.gzi` when bgzipped), naming the `samtools` command that makes what is
  missing; that the reference has the contigs the reads have reads on, at the header's length
  (contigs without reads are not checked, so a chr20 reference serves reads of chr20 whose
  header names a whole genome; a reference naming its contigs `20` for reads on `chr20` is
  refused, where it gave an empty VCF; a BAM with up to 1% of its mapped reads, by its index's
  counts, on contigs the reference lacks, as HLA, decoy or EBV contigs against the no-alt
  analysis set, is warned of with each contig's reads, and more is an error, as is any for
  CRAM, whose reads there cannot be decoded); `[snv] device`; and that lrpipe-snv's model is there,
  with its `pileup.pt` and `full_alignment.pt`. These hold for `--dry-run` too. It warns when
  `--sample` is not the sample the reads' read groups name (`@RG SM`, of every file, unaligned
  BAMs' included), and, at the end, of
  outputs in the output directory that an earlier run with other analyses, input or options
  wrote and this one does not, the aligned, merged and haplotagged reads' BAMs and the QC files
  among them (they
  are left, and are not in the manifest; haplotagged reads older than this run's haplotag list
  are said to be an earlier phasing's).
* When a stage fails, lrpipe prints the last 10 lines of its log that are not blank, after the
  error naming the log.
* `--analyses LIST` (default: all) runs only some analyses, comma-separated: `all`, `snv` (SNVs and
  indels), `phase` (the same phased), `sv`, `cnv`, `tr` (tandem repeats, when the configuration
  has a catalogue), `depth` and `methyl`. What an analysis needs runs too: `phase` runs `snv`'s
  calling, and `cnv` the SNVs and the read pass's depth and SVs. `--analyses snv,sv` runs
  lrpipe-snv without phasing and a read pass that calls SVs only; phasing the final calls adds
  1.5 to 3 min to a whole genome (HG002 HiFi 6 min 32 s against 4 min 57 s).
* `--merged-vcf off` leaves out `SAMPLE.vcf.gz`; each tool's own files are written either way.
* `--regions BED` restricts every analysis to the BED's regions (a panel's targets, adaptive
  sampling's, one chromosome for a test run), each by its tool's own region option, and leaves
  out the CNVs ([Regions](#regions)).
* `--gvcf` also writes lrpipe-snv's gVCF (Clair3's `--gvcf`), `snv/merge_output.gvcf.gz`, for
  joint genotyping. Turning it on or off reruns lrpipe-snv's stage and those reading its calls
  (cnv and merge), whose records gain or lose Clair3's PL, but not the read pass.
* `--sex auto` (the default) reads the sex chromosomes from the index in a few milliseconds:
  chrX's and chrY's reads per base (BAI or CSI read counts, or CRAM slice bytes) against the
  median autosome's. chrX below 0.75 with chrY from 0.2 is XY; chrX from 0.75 with chrY below
  0.2 is XX (HG004, a woman: chrX 0.986, chrY 0.083 from reads mapped to its repeats; HG002,
  her son: chrX 0.48, chrY 0.49 by CRAM and 1.4 by BAI); anything else (45,X, 47,XXY, an index without autosomes) stops the run and asks
  for `--sex female` or `--sex male`. `lrpipe sex READS` prints the estimate alone. Of several
  read files, the counts are added up over all of them (BAMs' reads and CRAMs' bytes each
  measured on their own, and the two averaged). The median is of the 22 autosomes, those the
  header lacks counting as without reads, so reads on fewer than 11 of them (a chr20 BAM) give no
  sex, saying so; for unaligned reads, a reference whose `.fai` has fewer than 11 autosomes, or
  no chrX, is refused before they are aligned. A `--sex` the index contradicts is a warning
  ([QC summary](#qc-summary)).
* `--memory GIB` (`32`, `32G` or `32GiB`, at least 1; `M` and `T` suffixes too; default: memory available at the start, or what is left under the limit of
  lrpipe's cgroup if that is less: a SLURM job's `--mem`, `docker run -m`, a Kubernetes pod's
  limit) caps what stages running together may use: a stage waits while its expected peak would
  not fit beside the running ones'. A stage's expected peak is an estimate from whole genomes
  with 16 threads, or the peak measured when it last succeeded in the same output directory plus
  15% if that is higher; a stage that was stopped or failed did not reach its peak, so what it
  measured counts only where it is higher, without the 15%. The
  estimates are lrpipe-snv's 13 GiB on HiFi and 14.5 on ONT, the read pass's 5 on HiFi and 7 with
  ONT's bedMethyl, `lrpipe methyl`'s CpG tables' 2.5 on HiFi, and merge-reads' 0.5 and 0.25 for
  each CRAM input (it streams: 0.13 GiB measured on two BAMs of HG003 ONT chr20), with room for lrpipe-snv's peak,
  which moves by 2 GiB from run to run. Fewer threads, or a smaller genome, need less, but a
  measured peak below the estimate does not lower it.
  On a machine with too little for both, lrpipe-snv and the read pass run one after the other.
  Less than a stage's estimate is warned of, naming the stages: each runs with no other beside
  it, and may run out of memory.
* `--threads N` (default: all cores; at least 1, and more than lrpipe may use is warned of) is the cores the run uses in all. lrpipe-snv and the read
  pass run at the same time, each with N threads on the same N cores: lrpipe-snv keeps fewer busy
  while it waits on the GPU, and has them all once the shorter read pass is done. On Linux lrpipe
  confines itself and its stages to the first N cores it may use; elsewhere lrpipe-snv takes the
  larger half of the N threads and the read pass the rest.
* lrpipe-snv runs with `--approximate`, which changes only CPU calling: there its networks run
  in 8 bits where a check on the run's own candidates allows it, with close calls scored again in
  `float32`. The whole HG003 ONT genome took 31 min against 72 min exact, and HG002 ONT chr20
  60 s against 139 s (HiFi 55 s against 98 s), with 99.2%
  of the PASS calls the same as exact (99.96% on HiFi), the same recall against GIAB and slightly
  lower precision. On a CPU without AVX-VNNI, AVX-512 VNNI or aarch64's dot product
  instructions it changes nothing, with a warning. On a GPU the networks run in `float32` either way (a `float16` mode saved 5% of a
  whole HG003 ONT run from a BAM, and nothing from a CRAM).
* `--display` (`auto`, the default; `tui`; `plain`) is how the run shows its stages. On a
  terminal it takes the screen: every stage with its tool and a progress bar, the run's overall
  progress and an estimate of the time left, the memory the running stages use against
  `--memory`, each stage's time and peak memory, and the end of a stage's log (↑ and ↓ choose
  which). Only the run's stages are listed, and the header names its analyses: a stage of an
  analysis not asked for has no row, and the read pass's row names only what it computes
  (`lrpipe-sv: SVs, depth`). From unaligned reads, the stages after alignment are planned once
  the reads are aligned; until then a line says which analyses follow. bonefish's bar follows
  how far into the reads it is, the reads aligned and its own estimate of the time left (from
  bonefish 0.5.1, which writes them to the file `$LRPIPE_PROGRESS` names), lrpipe-snv's its
  pileup chunks and full-alignment jobs, and the read pass's the contigs it has read; another
  stage's follows its time against what it took when it last succeeded in the output directory,
  marked `~`, or only shows it runs. `q` leaves the display
  and the run goes on printing its lines; Ctrl-C stops the run as SIGINT does. When the display
  ends, the lines a run without it prints are printed, so the terminal keeps the same record.
  Elsewhere (a pipe, a file, a batch job), or with `plain`, the run prints a line for each
  stage's start and end, as before.
* `--dry-run` prints lrpipe-snv's model and why it is the one (below), and each stage's command,
  quoted as a shell takes it, marking with `: done, would be skipped` each stage a run would skip
  (it succeeded with the same command and inputs, as did every stage it follows). It checks what
  the run checks before its first stage, and exits with 2 where the run would stop: the files
  the stages read (a missing `[sv] tandem_repeats`, `[repeats] catalog` or `[cnv] blacklist`),
  and an output directory that holds another sample's run or is in use by another run.
* `--model NAME` (a model in the site's `[snv] models`, as `r1041_e82_400bps_sup_v500`) or
  `--model DIR` (a model directory: a path with a `/`) calls the SNVs with that Clair3 model,
  in place of the one the site configuration gives or lrpipe chooses from the reads' header (see
  [Clair3's model](#clair3s-model)). `--model-mismatch warn` runs with a model that the header
  contradicts, which is otherwise an error.
* Haplotypes: when the run phases the SNVs and indels (`phase`, by default), it haplotags the
  reads from them (the **haplotag** stage), and the SVs' phases, the tandem repeats' alleles
  and the CpG tables per haplotype follow that phasing: their haplotypes and phase sets are the
  phased SNVs'. Reads straight from minimap2, pbmm2 or dorado, without `HP` tags, get phased
  SVs, repeats split by haplotype and `hap1`/`hap2` tables, and reads haplotagged by another
  pipeline get the run's haplotypes everywhere, not their own tags'. The manifest's `sample.haplotypes` is `lrpipe run's phasing of the SNVs (its haplotag
  list)`. No second BAM is written unless asked for: `--haplotagged-reads` writes
  `phase/SAMPLE.haplotagged.bam` and its index too, as WhatsHap's `haplotag` tags it, for IGV
  and other tools. `--haplotag off` keeps the SNVs phased but leaves the haplotagging out: the
  SVs, the repeats (genotyped in the read pass again) and the CpG tables run beside lrpipe-snv,
  as in 0.1.5, a few minutes less on a whole genome (BENCHMARKS.md has the measured cost), and
  follow the reads' own tags, with a warning when those are another pipeline's. Without phasing,
  or with `--haplotag off`, the SVs' phases, the repeats' alleles and the CpG tables per haplotype
  are read from the reads' integer `HP` tags. Before planning, the run looks for one in 200 mapped reads
  from the middle of each of the first 8 contigs with reads (of 1 Mb or more first), so a file
  tagged only on some contigs, or with few phased sites where a contig starts, counts as tagged:
  0.08 s for a 39 GB HiFi BAM, and of a CRAM only the flags, positions and the tag are decoded.
  Reads without one then get unphased SVs, repeats told apart by length alone and no
  `hap1`/`hap2` CpG tables (rather than empty ones), with a warning at the start and in the summary that adding `phase` to `--analyses`
  (or leaving out `--haplotag off`) phases them, and the manifest's `sample.haplotypes` is `none (input not haplotagged)`. `HP`
  tags of another type (`HP:Z`), which the tools read as none, count as none, and the warning
  says so.
* A run that stopped resumes where it stopped: a stage is skipped when its command, program and
  inputs are unchanged since it last succeeded and its outputs, indexes included, all exist.
  Thread counts and how paths are spelt do not count as changes: a run resumed with another
  `--threads` (or on a node with other cores), with `--config ./site.toml` for
  `--config /data/site.toml`, or from another directory, skips what is done. So does a run in an
  output directory moved or copied with its files' times (`mv`, `cp -a`, `rsync -a`): its paths
  count as paths under it. No stage reads a path out of an earlier stage's output: of the
  outputs, only the aligned BAM's `@PG` line (bonefish's command) names the directory it was
  written in, as a record of how it was made, which the haplotagged reads' header carries on.
  Inputs up to 64 MiB, such as the networks' weights (`pileup.pt`, `full_alignment.pt`),
  catalogues and indexes, are compared by contents, larger ones (the reads, the reference) by
  size and time, so a model replaced in the same directory reruns calling. A change reruns the
  stage and every stage after it; `--rerun all` reruns everything.
* `--keep-intermediates` keeps lrpipe-snv's `tmp/`.
* Ctrl-C (SIGINT), SIGTERM, SIGQUIT or SIGHUP (a terminal closed, an ssh session dropped)
  stops the running stages, records the run as `interrupted` and their stages as `cancelled`,
  and exits with 130; running it again resumes. This holds while unaligned reads are being
  aligned as after it (a signal as the align stage ends starts no other stage). Under `nohup`,
  which ignores SIGHUP, a hangup
  leaves the run going. No stage starts once the run is asked to stop. Each stage runs in a
  process group of its own, so a terminal's signals reach lrpipe alone, which kills the stages'
  groups, with whatever the stages started; a group is killed too when its stage ends, so
  nothing a stage left running outlives it. A stage that dies of SIGINT or SIGTERM itself
  (SLURM's `scancel` signals every process of the job) is stopped with the run too, not failed;
  one that crashes (SIGSEGV, SIGABRT) is failed. Another signal changes nothing: the run ends
  once its stages have. A stage still running 10 s after it was killed is recorded as stopped
  and left, and lrpipe prints its process ID (the manifest records it): do not run again into
  the directory until it has ended. Ctrl-Z suspends lrpipe alone, not its stages, which go on
  running. A run that cannot write its own files (a full disk) kills its stages and records
  itself as failed where it can. Stages die with lrpipe however it ends, even by SIGKILL (on Linux;
  elsewhere only SIGINT and SIGTERM stop them), and one output directory takes one run at a
  time.
* One output directory holds one sample: `snv/`, `logs/` and the manifest are not named after
  it, so a run of another sample into a directory that holds a sample's run (by
  `.lrpipe/sample`, or an older run's manifest) is refused before it writes anything there.

### Bacterial isolates

Use `--bacterial` for one haploid bacterial isolate. Set the site configuration to the isolate's
indexed reference and configure the caller model for the sequencing platform. The normal
Bonefish alignment stage handles FASTQ, FASTQ.gz, and unaligned BAM inputs.

```bash
lrpipe run --bacterial --config k12.toml --reads k12.fastq.gz --sample K12 \
  --platform ont --model r1041_e82_400bps_sup_v500 --output k12 --gvcf
```

The preset calls SNVs and small indels with a single allele per genotype and skips human sex
inference. It considers all reference contigs, including short replicons, or the requested BED
targets with `--regions`. `--gvcf` adds haploid reference blocks and variant calls to
`snv/merge_output.gvcf.gz`. The run records ploidy 1 in `manifest.json`. The supported analysis
is `snv`, including when `--analyses` is omitted or set to `all`.

The [K-12 smoke validation](../../tools/validation/results/bacterial-k12-2026-11-24/README.md)
records end-to-end runs on ONT and HiFi FASTQs.

### Several read files

A sample sequenced on several flow cells or SMRT cells, and aligned already, can be given as
it is, `--reads` once for each file, BAM and CRAM alike:

```bash
lrpipe run --config site.toml --reads cell1.bam --reads cell2.bam --reads cell3.cram \
  --sample HG002 --platform hifi --output hg002
```

A first stage, `merge-reads` (`lrpipe merge-reads`, in Rust, without samtools), merges them into
one sorted, indexed BAM, `OUT/reads/SAMPLE.bam`, which every other stage reads. It streams the
files, merging their records by position as `samtools merge` does (at one position forward strand
first, then in the order the files were given), with each file read by a thread of its own and
the output compressed by all the run's threads: two BAMs of HG003 ONT chr20 (2.5 GB in all) took
7.5 s into memory against `samtools merge -@16`'s 12.5 s, with a peak of 0.13 GiB (writing to this
machine's disk took about a minute, for both). The merged BAM needs disk space of about the
inputs' size (more where they are CRAM, which is about 60% of a BAM's), and is kept, so a resumed
run does not merge again: each file and its index are the stage's inputs, and one replaced
merges again and reruns what reads the merge. It is a BAM rather than a CRAM because the stages
decode BAM faster, which they all do, and it needs no reference to read. Unplaced unmapped reads
(without a position, at a file's end) are left out: no stage reads them. One `--reads` is read
as it is, with no merge. If a single aligned BAM or CRAM has several distinct `@RG SM`
labels, `lrpipe run` stops before starting any stage. For reads verified to come from one
donor, repeat `--sample-alias` for every label other than the output name given by `--sample`.

```bash
lrpipe run --config site.toml --reads hg002.bam --sample HG002 \
  --sample-alias HG002-cell2 --sample-alias HG002-cell3 --platform hifi --output hg002
```

Before anything runs, the files are checked to be one sample's, aligned to one reference:

* their `@SQ` lines must be the same (names, lengths and order), or the run stops naming the
  first that differs;
* their read groups are combined: an `@RG` ID in two files must have the same fields in both
  (the same read group given twice is one), or the run stops naming the ID and how to mend it
  (`samtools addreplacerg` or `samtools reheader`); and they must name one sample (`SM`), or it
  stops naming the samples and their files; every file's `@PG`
  chain is kept whole, as `samtools merge` keeps it (a line of an ID another file has with other
  fields is renamed `ID-N`, for the Nth file, and the lines after it follow it), and the merge
  adds its own `@PG` line after each; comments are kept, each once;
* a read in two of them is warned of (from 500 reads' names sampled at the same places of the
  same contigs in each file, the first long ones with reads in two files, a few BGZF blocks or
  CRAM containers of each): it would count twice. Two files sharing more than 90% of the smaller
  sample hold the same reads (a BAM and its CRAM, a file and its copy) and are refused; a third
  or more is warned of as what a file and a merge holding it share.

Clair3's model is chosen from all their read groups together (several of one basecall model are
fine; [different ones](#clair3s-model) need `--model`); without phasing, each file is sampled
for `HP` tags, with a warning when some are haplotagged and others not (with phasing, the run
haplotags the merged BAM itself); and the sex is read from all their indexes.
The manifest's `sample.reads` lists the files (a list for one file too).

merge-reads takes aligned files only: several unaligned ones are aligned into one BAM by the
align stage ([Unaligned reads](#unaligned-reads)), and aligned and unaligned files together are
refused.

### Unaligned reads

```bash
lrpipe run --config site.toml --reads flowcell1.bam --reads flowcell2.bam --sample HG002 --platform ont --output hg002
```

`--reads` also takes unaligned reads: unaligned BAMs, as dorado and PacBio's instruments write
them, or FASTQ or FASTQ.gz, given once for each file. A file's format is told from its contents:
a BAM without `@SQ` lines is unaligned, and so is one with them whose first 1,000 reads are all
unmapped. An align stage then maps every file, one after another,
into one sorted, indexed BAM, `OUT/align/SAMPLE.bam`, which every other stage reads as it would
aligned input:

```
bonefish -t N -x dna -y -o OUT/align/SAMPLE.bam --sort-tmp OUT/align/tmp/SAMPLE INDEX_OR_REFERENCE READS...
```

* [bonefish](https://github.com/DylanCosto/bonefish) (MIT) is the aligner: `-x dna` (minimap2's
  `lr:hq`) for `--platform ont` and `-x hifi` (`map-hifi`) for `hifi`, validated with Clair3 and
  Sniffles2 on whole HG002 and HG005 genomes at minimap2's SNP, indel and SV F1. `-y` keeps each
  read's tags (`MM`/`ML`/`MN` of base modifications, `RG`, `qs`, ...) and the input's `@RG`
  lines, which name the basecall model. It is `[tools] bonefish`, else a `bonefish` beside the
  `lrpipe` binary, else on `PATH`, checked with `bonefish --version` before anything runs and
  recorded in the manifest's `tools` (version and SHA-256).
* It maps to `[align] index`, bonefish's index of the reference (`lrpipe setup --bonefish`
  builds it; 8.5 GB for GRCh38, loaded in under a second), else to the reference, which bonefish
  indexes on every run (about a minute for GRCh38). Before aligning, lrpipe checks that the
  index's contigs (as `bonefish -a INDEX -` prints them) are the reference's `.fai`'s, names,
  lengths and order, and after it the aligned BAM's; either differing stops the run, saying
  where.
* Read groups: unaligned BAMs' `@RG` lines are kept, and must be of one sample (one `SM`), with an
  ID given different lines in two files refused; dorado's files of one run, with the same `@RG`
  line, are fine. FASTQ has no `@RG` lines, and dorado's carries `RG:Z` tags naming read groups
  no header holds, so its reads are given one read group named after the sample (bonefish's
  `-R '@RG\tID:SAMPLE\tSM:SAMPLE'`), in place of those tags. Unaligned BAMs and FASTQ together
  are refused, and so is a file given twice.
* The align stage runs first and alone, with every core. The sample's sex (`--sex auto`) and the
  reads' `HP` tags (which matter only without phasing) are read from the BAM it writes, so the other stages are planned once it is
  done; a `--dry-run` before then prints the align stage alone. What would stop the run after
  aligning stops it before, in seconds, in a dry run too: the Clair3 model is chosen from the
  unaligned BAMs' `@RG` lines, which bonefish keeps ([Clair3's model](#clair3s-model)), and the
  rest is planned (with the sex given, or XX) and its files checked (`--analyses tr` without a
  catalogue, a `.fai` or catalogue that is not there). A failure after the align stage and
  before the others, or Ctrl-C then, finishes the manifest as `failed` (with the run's `error`)
  or `interrupted`, with no other stage started. FASTQ has no header: the model is chosen from
  the basecall model dorado's FASTQ names in each read's `RG:Z` tag (`RG:Z:RUNID_MODEL`, of the
  first 100 reads), and an ONT run of FASTQ without them needs `--model` or
  `snv.ont_model`, and a HiFi run takes `hifi_revio`, with a warning.
* Memory: bonefish's documented peaks for a human genome with 16 threads, with `-y` carrying
  methylation tags, are about 23 GB with its defaults and 19 GB with `--sort-mem 1G`; with
  `-t 4 -K 10M` it maps in under 9 GB, to which the sorter's two 1 GB batches add 2. Each thread
  more or fewer moves a peak by 0.15 GiB (bonefish's thread-scaling measurements), and without an
  index the build before mapping peaks at 10.5 GB, which the stage's peak is at least. The align
  stage takes the first setting that fits `--memory`: with 16 threads, 23 GiB or more runs
  bonefish's defaults, 19 to 23 `--sort-mem 1G`, and less `-t 4 -K 10M --sort-mem 1G` (11 GiB).
* Resuming never realigns: the align stage is fingerprinted like any other, by the reads, the
  reference, the index and the bonefish binary. Its threads and batch sizes (`-t`, `-K`,
  `--sort-mem`) are left out, as none changes what bonefish writes (its output is the same at
  every thread count and batch size, and its sorter's batches are merged stably); so a resume
  with other threads or memory skips it. New reads, another reference (with an index or
  without), another index or another bonefish realign, and rerun everything after. The sorter's
  scratch directory, `OUT/align/tmp/`, is emptied before the stage starts and removed once it
  succeeds.
* `OUT/align/SAMPLE.bam` is the size of a BAM: bonefish writes no CRAM. Its index is a `.bai`,
  or a `.csi` where a contig is longer than 2^29 bases. The manifest records it as
  `sample.reads` and `outputs.aligned_bam`, and the unaligned files as `sample.unaligned_reads`.
* Methylation needs the reads' `MM`/`ML` tags, both. An unaligned BAM carries them on its
  reads, when basecalled with modified bases; a FASTQ only in its records' comments (dorado's
  FASTQ, or `samtools fastq -T MM,ML`), its bases on one line or wrapped. A file without them on
  its first 100 reads gives no methylation. When no `--reads` has them, `methyl` asked for in `--analyses`
  is an error, and the default analyses leave methylation out with a warning, as they leave out
  tandem repeats without a catalogue, rather than write empty tables; when only some do, the
  run warns that the others' reads give no methylation. Aligned `--reads` are held to the same
  rule, by the reads sampled for `HP` tags (200 from the middle of each of the first 8 long
  contigs): reads aligned from FASTQ without them, or basecalled without modified bases, give
  no methylation stage.
* Aligned and unaligned reads together are refused (several aligned files are merged:
  [Several read files](#several-read-files)), and so is an unaligned CRAM, which bonefish cannot
  read (`samtools view -b` makes a BAM of it). Aligned BAM or CRAM is never sent
  through bonefish.

### Regions

`--regions BED` restricts every analysis to the BED's regions: a targeted panel, the targets of
ONT adaptive sampling, or one chromosome for a test run.

```bash
lrpipe run --config site.toml --reads HG002.bam --sample HG002 --platform ont --output panel \
  --regions targets.bed
```

The BED is checked before planning, for `--dry-run` too: each line a contig, a start and an end
(`0 <= start <= end`; a zero-length row, an insertion point, is the base at its start, as Clair3
reads it; more columns are ignored, `#`, `track` and `browser` lines skipped, a byte order mark
dropped, a gzipped BED read too), on a contig of the reference's `.fai` and of every aligned
`--reads`' header, ending within it, and at least one region. A contig without reads is
accepted, as Clair3 takes it, with a warning that its regions get no calls. Regions on contigs a
run of the whole genome calls nothing on (neither lrpipe-snv's chr1-22, X and Y nor lrpipe-sv's
contigs of 1 Mbp or more: the alt, decoy, HLA, EBV and mitochondrial contigs) are left out with
a warning, and regions that cover every contig such a run calls, whole (a BED made from the
`.fai`), are no restriction: the run is the whole genome's, CNVs included. The rest is sorted in
the reference's order and merged where regions overlap or touch, into `OUT/.lrpipe/regions.bed`,
which the stages read; its contents are part of their fingerprints, so other regions rerun them
and the same regions in another order do not. The manifest records the BED given as
`resources.regions`.

Each analysis is restricted by its tool's own region option, with the original's edges:

| Analysis | Restricted by | At the regions' edges |
|---|---|---|
| snv, phase | lrpipe-snv's `--bed-fn` (Clair3's `--bed_fn`) | as Clair3: calls whose candidates overlap a region, from the regions widened by 33 bases. Its quality cut-offs (the share of pileup calls called again with full alignment, and phased from) are drawn from the regions' calls, so calls near them can differ from a whole genome's run; with the same BED, lrpipe-snv's PASS calls on the ONT fixture are Clair3 2.0.3's (`lrpipe-snv/ci/fixture/expected/window`). The gVCF covers the widened regions |
| sv | lrpipe-sv's `--regions` (Sniffles2's) | as Sniffles2: the reads overlapping a region are read, and an SV is called from signals that start in it, so one that reaches a region from before it is not (checked against Sniffles 2.8.1's `--regions` calls) |
| depth | mosdepth's `--by` with a BED of the whole genome's 1 kb windows, each cut to the regions (`OUT/.lrpipe/regions.windows.bed`), counted from the reads of the regions alone, in the read pass (`lrpipe sv --depth-in-regions`) | `depth/SAMPLE.regions.bed.gz` has the windows the regions reach, each with mosdepth's value from all the reads that reach it, so a window wholly in a region has the whole genome's; the summary has the regions' rows alone (`chr20_region`, `total_region`): the contigs' depth is not counted, and no `mosdepth.global.dist.txt` is written |
| tr | the catalogue's loci that a region reaches (`lrpipe tr --regions`, the read pass's `--repeats-regions`) | a locus is kept whole, genotyped from all its reads |
| methyl | ONT's bedMethyl: the read pass's `--bedmethyl-regions`, as `lrpipe methyl pileup --include-bed`; HiFi's CpG tables: `lrpipe methyl cpg --include-bed` | a CpG is counted when its C is in a region, both strands' calls, and the bedMethyl's pass thresholds are estimated from the calls of the regions' reads. modkit's `--include-bed` differs at the edges: its docs say it keeps the "positions overlapping intervals", and modkit 0.6.4 counts each strand's position on its own (a CpG split by an edge gets the inside strand's calls alone) and estimates the thresholds from the calls at the regions' positions alone |
| haplotagging | `lrpipe haplotag --regions-bed` | the reads overlapping a region, each once (WhatsHap's `--regions` writes one twice when it reaches two), with WhatsHap's haplotypes in those regions; `--haplotagged-reads` holds those reads alone |
| cnv | not called | lrpipe-cnv normalises each window's depth against the whole genome's, which the regions alone do not give, so its calls would mislead: the default analyses leave cnv out, saying so, and `--analyses cnv` with `--regions` is an error. Call CNVs in a run of the whole genome |

The SVs' phasing and the merge read what those stages wrote. When every region is a whole
contig (`chr20 0 64444167`), the run is restricted by contig instead: lrpipe-snv's `--ctg-name`
and the read pass's `--contig`, whose depth (in the whole genome's 1 kb windows, with mosdepth's
whole-contig rows), repeats and bedMethyl are counted in the pass itself. A run of chr20 then
gives a whole-genome run's chr20 records for the depth, the SVs and the tandem repeats' loci.
Not for every SNV: lrpipe-snv, as Clair3, draws its quality cut-offs from every contig it calls,
so a few calls near them differ (on the ONT fixture made two contigs, `lrpipe/tests/run.rs`, the
chr20 run shares 496 of the whole run's 497 PASS calls and has 3 others), and so may what
follows from the SNVs: their phasing, the haplotag list, and the SVs' and repeats' haplotypes.
The bedMethyl's thresholds are estimated from chr20's reads, as `lrpipe methyl pileup --region
chr20` estimates them.

The sex is still read from the index, which covers the whole genome, unless `--sex` gives it.
The align stage still aligns all the unaligned reads, and merge-reads merges all the files: the
regions restrict what is called from the reads, not the BAM.

## Comparing runs

```bash
lrpipe compare old-run new-run
```

checks that a change (to a tool, or a library under all of them) keeps every result: each
output's records, in order, and how many records only one run has; the first record where they
part, shown around the first differing character; header lines other than dates and command
lines (`##fileDate`, `##command`, `##cmdline`, and the merged VCF's `##snv_cmdline`,
`##sv_fileDate`, ...); tools whose binary changed; and each stage's wall time and peak memory,
old to new. A stage one run skipped (already done) shows the measurements of the run that did
it, marked so, with no change in per cent; so does the run's wall time when one of them ran no
stage, and a value that shows as 0 (`0 -> 7` is no +2404%).
It exits with 0 when every output's records match, 1 when any differ, and 2 on error. Outputs
are found under the given directories even if the runs were moved. An `.snf`'s records are its
blocks of SV candidates, decompressed. A BAM's (the aligned reads, the merged reads, the
haplotagged reads) are its alignments, each with every field and tag, those at one place as a
multiset (their order there is not the sort's to keep), and its header lines but for its
programs' command lines (`@PG CL:`); the two BAMs are read once, side by side, and a record
where they part is shown by its read name, place, flags and MAPQ. A BAM without its BGZF
end-of-file block (cut short) is an error, as is an output named as text that is not; any other
output that is not text is compared byte for byte.

## Stages

```
                                ┌─ snv (lrpipe-snv) ── haplotag (lrpipe haplotag) ─┬─ methylation (lrpipe methyl)
                                │                                                  ├─ repeats (lrpipe tr) ──────┐
reads ──[align | merge-reads]───┤                                                  │                            │
                                └─ reads (lrpipe sv) ─────────────────── sv-phase (lrpipe sv phase) ── cnv ── merge ── qc
```

The align stage runs only for unaligned reads ([Unaligned reads](#unaligned-reads)), and
merge-reads only for several aligned files ([Several read files](#several-read-files)); either
runs alone, before the rest, which reads the BAM it writes (the haplotag and repeats stages and
the CpG tables too).

(Without phasing there is no haplotag, repeats or sv-phase stage: the read pass writes the SVs
and the repeats, and methylation starts with the rest. With phasing but only CNVs among the
variants, the read pass writes the SVs for their breakpoints, unphased.) snv, haplotag, repeats,
sv-phase and cnv are the run's longest path, so they run at normal CPU priority; the read pass
and methylation run alongside at a lower one (`nice` 10), taking the CPU lrpipe-snv leaves idle
while it waits on the GPU without slowing it. Stages start together as far as their expected
memory fits.

* **snv**: `lrpipe snv`, with the Clair3 model for the reads ([below](#clair3s-model)) and the sample's sex: no chrY
  calls in XX samples, and in XY ones haploid chrX and chrY calls outside the pseudo-autosomal
  regions, which lrpipe writes to `.lrpipe/par.bed` from the reference (GRCh38, GRCh37 or
  T2T-CHM13, by chrX's length; on another reference XY samples' chrX and chrY are called diploid,
  with a warning); with `phase`, Clair3's `--enable_phasing`, which phases the final calls with
  lrpipe-phase as WhatsHap would, leaving the haploid ones unphased (HG002 HiFi chrX and chrY:
  111,068 and 13,027 haploid calls, all 13,113 diploid ones inside chrX's PARs, 8,643 of them
  phased). On the GPU each contig is merged and phased as soon as its full-alignment calls are done, beside the
  other contigs' calling (on HG002 HiFi 29 s sooner than phasing every contig at the end, with
  the same files).
* **reads**: `lrpipe sv`: SVs, with `--accuracy` unless `--sv-mode sniffles2` (a repeat
  expansion's two alleles called apart, a read's indels joined within 50 bp and long insertions
  genotyped from their reads, rather than Sniffles2's calls; [lrpipe-sv](../lrpipe-sv#accuracy-improvements-opt-in)), or none when not
  chosen (the pass then reads no SV evidence); depth in 1 kb windows (mosdepth `-x -Q 20 -b 1000`, what lrpipe-cnv needs);
  the reads' statistics for the QC summary (`--read-stats`), counted from the records it reads;
  without phasing, tandem repeats with lrpipe-tr's platform noise model and the sample's
  karyotype, from the reads' own `HP` tags. With phasing, the pass calls the SVs as if the reads had no
  haplotypes (`--phasing-evidence`) into `sv/SAMPLE.unphased.sv.vcf` and `.unphased.snf`, and
  writes what phasing them needs, `sv/SAMPLE.phasing-evidence.gz`: each call's reads by name.
* **haplotag**: `lrpipe haplotag --list-only` on `snv/phased_merge_output.vcf.gz` and the reads,
  after lrpipe-snv, as Clair3's WhatsHap runs it (`--ignore-read-groups`, the reference): the
  haplotype and phase set of each read, `phase/SAMPLE.haplotags.tsv.gz` (WhatsHap's
  `--output-haplotag-list`: read name, `H1`/`H2`/`none`, phase set, chromosome), the same as
  WhatsHap 2.8 writes. Supplementary alignments are not tagged (WhatsHap's default), and reads on
  contigs the phased VCF does not name are untagged (`--untag-missing-contigs`). With
  `--haplotagged-reads`, also the tagged BAM and its `.bai` (or `.csi` for contigs too long for
  it). Only when something reads haplotypes: the SVs, the tandem repeats, the CpG tables or
  `--haplotagged-reads`.
  On whole-genome HG002 (16 threads) it took 96 s and 2.8 GiB from the HiFi BAM and 253 s and
  5.6 GiB from the ONT CRAM ([BENCHMARKS.md](../../BENCHMARKS.md#haplotagging)).
* **repeats**: `lrpipe tr --haplotags`, after haplotag, with phasing: the tandem repeats from the
  reads, as the read pass genotypes them (the same records for the same haplotypes), with the
  haplotag list's haplotypes. A locus whose reads mostly carry one is split into alleles by
  haplotype, which on HG002 ONT chr20 gets 84.4% of the called regions' alleles exactly right
  against GIAB's TR benchmark, against 69.2% from lengths alone
  ([lrpipe-tr](../lrpipe-tr)). It reads the reads again: on whole-genome HG002 (16 threads, the
  machine loaded) 88 s and 3.9 GiB from the HiFi BAM, 231 s and 6.7 GiB from the ONT CRAM; on
  HiFi it runs beside the CpG tables, on ONT it lengthens the run by about that
  ([BENCHMARKS.md](../../BENCHMARKS.md#haplotagging)).
* **sv-phase**: `lrpipe sv phase`, after the read pass and haplotag: the SVs phased from the list,
  as Sniffles2 phases them from reads tagged as it says (each call's PHASE, genotype phase and
  PS, and its phasing rescue of MOSAIC_VAF calls), into `sv/SAMPLE.sv.vcf` and `sv/SAMPLE.snf`:
  the VCF and SNF the read pass would write from the haplotagged BAM
  ([lrpipe-sv](../lrpipe-sv#phasing-from-a-haplotag-list)). The CNVs and the merge read these.
* **cnv**: `lrpipe cnv` (lrpipe-cnv) with `--improved` (GC correction, amplifications, copy numbers from the
  final calls), the sample's sex, lrpipe-snv's heterozygous SNVs, and breakpoints from the SVs.
* **methylation**: `lrpipe methyl` (lrpipe-methyl), from the reads' `MM`/`ML` tags. On HG002 ONT
  (whole genome, CRAM, 16 threads) it took 404 s and 2.4 GiB,
  with output byte-identical to modkit 0.6.4's, which took 816 s and 14 GB. By default each platform gets its own tool's tables: for Oxford Nanopore,
  modkit's bedMethyl (`pileup --modified-bases 5mC 5hmC --cpg --combine-strands`), from the read
  pass, for all reads (not per haplotype); for PacBio HiFi, pb-CpG-tools' CpG tables for all
  reads and for haplotype 1 and haplotype 2 (count mode, as pb-CpG-tools writes it). With
  phasing, the tables wait for the haplotag list and take the reads' haplotypes from it
  (`--haplotags`), the tables of the alignments `haplotag` would write; without, they start
  with the run and per haplotype only for haplotagged reads. lrpipe-methyl warns in its log when
  the reads carry no tags.

* **merge**: `lrpipe merge`: the chosen variant analyses' VCFs (the SNVs phased with `phase`)
  as one, sorted in the reference's contig order, bgzipped and indexed. Each record keeps its
  own columns and gains `INFO/SOURCE` (`snv`, `sv`, `cnv` or `tr`); the header holds each input's
  definitions once, renaming an INFO or FORMAT ID that two inputs define differently for the
  later one (lrpipe-cnv's Float `GQ` becomes `CNV_GQ`), and keeps the inputs' other header lines
  prefixed with their source (`##snv_cmdline=`). On HG002 HiFi (5.4 million SNVs and indels,
  937,000 repeats) it takes 7 s. `lrpipe merge --reference REF --sample S --output OUT.vcf.gz
  --snv A --sv B ...` runs it alone. It writes under hidden names of its own beside `OUT.vcf.gz`
  and renames them into place when done, so no other file there is touched (`OUT.vcf`
  included) and a failed merge keeps an earlier `OUT.vcf.gz`; an `--output` that is one of the
  inputs is refused.
* **qc**: `lrpipe qc`, last, after every other stage: the run's QC summary
  ([below](#qc-summary)). It reruns when a file it reads changes, or the sample's sex. Its
  failure does not fail the run: the manifest records it (the stage `failed`, with why) and the
  run warns that there is no summary, but completes, and exits as its analyses make it.

If a stage fails, the stages still running are killed, none start, and the manifest and the
stage's log in `logs/` say why.

## Clair3's model

A Clair3 model is trained on one platform's reads, and for Oxford Nanopore on one chemistry and
basecall model; the wrong one costs accuracy with no error. The reads' header says what made
them: dorado writes `@RG ... PL:ONT DS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, and
PacBio `@RG ... PL:PACBIO PM:REVIO DS:READTYPE=CCS;...`. lrpipe-snv's model is, first found:

1. `--model`;
2. the site's `[snv] ont_model` or `hifi_model`, for the run's platform;
3. the model for the header from the site's `[snv] models` directory, which holds Clair3's
   models under their standard names: Clair3's own, as its Bioconda package and Docker image
   install them in `models/`, and any of ONT's Rerio models converted for Clair3 2
   (<https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/>), such as those for
   dorado v4.0.0, v4.2.0 and v4.3.0, for R10.4.1 at 260 bps, and for Guppy 6's R10.4.1 models.
   dorado's `basecall_model` gives the ONT model (`..._sup@v5.0.0` takes
   `r1041_e82_400bps_sup_v500`, `..._hac@v5.2.0` takes `r1041_e82_400bps_hac_v520`), and
   PacBio's instrument the HiFi one (`PM:REVIO` takes `hifi_revio`, `PM:SEQUELII`
   `hifi_sequel2`, the first Sequel's `PM:SEQUEL` `hifi_sequel2` too).

   Where the directory does not hold the model for the reads, or Clair3 has none (a version
   newer than lrpipe's table, or between two), the nearest it holds of the same chemistry is
   used: the same sampling rate (dorado's models from v4.2.0 are 5 kHz, earlier ones and
   Guppy's 4 kHz, unless the name says `5khz`), then the same accuracy (else the nearest, `fast`
   reads taking `hac`), then the same basecaller (Guppy's for a name without a version) at the
   nearest version, the newer of two as near, or the newest where the version is not known; a
   HiFi instrument Clair3 has no model for, or none named, takes `hifi_revio`. With Clair3's own
   models, dorado v4.3.0 SUP reads take `r1041_e82_400bps_sup_v500`, and Guppy's R10.4.1 SUP reads
   `r1041_e82_400bps_sup_v410`. Each is used with a warning that names it, the model for the
   reads and where that comes from. A chemistry lrpipe has no model for (R10.4, R10.4.1 at 300
   bps, RNA), an ONT header without a `basecall_model`, read groups basecalled with different
   models or on several instruments, and a directory holding no model of the chemistry are
   errors, until a model is given; a read group that does not say (`minimap2 -R
   '@RG\tID:x\tPL:ONT'`) beside ones that do is taken to be made as they were, with a warning.
   The fine-tuned models are used only when given. A `[snv] models` that does not exist, or that
   is one model's own directory (it holds `pileup.pt`: give the directory above it, or the model
   as `[snv] ont_model`), is an error saying so. The warning of a model other than the reads'
   gives why; the model's line (`snv model ...`) points to it rather than saying it again.

A model given (1 or 2) is checked against the header. One for another platform than the
reads' (a HiFi model for ONT reads), or another chemistry (an R9.4.1 model for R10.4.1 reads),
is an error, and so are Oxford Nanopore reads called `--platform hifi` or PacBio reads called
`ont`, and PacBio reads that are clearly not HiFi (`READTYPE` `SUBREAD`, `SCRAP`, `ZMW` or
`HQREGION`); `--model-mismatch warn` runs anyway, with a warning. A `PL` of another platform
(`PL:ILLUMINA`, stamped by an aligner's `-R`) and a `READTYPE` other than `CCS` that may be HiFi
(pbmm2 writes `UNKNOWN` for reads aligned from FASTQ) are taken as `--platform` says, with a
warning. A model for another basecall model or instrument of the same chemistry runs with a
warning that names Clair3's model for the reads; one whose name is not a standard Clair3
model's (a custom or fine-tuned model), or with nothing in the header to check it against,
runs with a warning that it is not checked. A move-table (`_with_mv`) model is refused: lrpipe
snv runs one only with `--enable_dwell_time`, which `lrpipe run` does not pass. `--model` on a
run without SNVs (`--analyses sv`) is not used, with a warning.

The manifest's `snv_model` records the model's directory, how it was found (`--model`,
`snv.ont_model`, `snv.hifi_model` or `snv.models`), why it is the one and the warnings; its
weights are in `resources` (`ont_model`, `hifi_model`). The model directory is in the snv
stage's command and its weights are the stage's inputs, so a run resumed with another model
calls again.

## Output

```
OUT/
  manifest.json
  SAMPLE.vcf.gz, .tbi                        the variant calls together (INFO/SOURCE)
  logs/{align,merge-reads,snv,reads,haplotag,repeats,sv-phase,cnv,methylation,merge,qc}.log
  align/SAMPLE.bam, .bai                     unaligned reads, aligned by bonefish (a BAM: no CRAM)
  reads/SAMPLE.bam, .bai                     several aligned --reads, merged
  snv/merge_output.vcf.gz                    SNVs and indels
  snv/phased_merge_output.vcf.gz             the same, phased
  snv/merge_output.gvcf.gz, .tbi             the same with reference blocks (gVCF), with --gvcf
  phase/SAMPLE.haplotags.tsv.gz              each read's haplotype and phase set, with phasing
  phase/SAMPLE.haplotagged.bam, .bai         the reads haplotagged, with --haplotagged-reads
  sv/SAMPLE.sv.vcf                           SVs (phased like the SNVs, with phasing)
  sv/SAMPLE.snf                              SV candidates, for multi-sample calling
  sv/SAMPLE.unphased.sv.vcf, .unphased.snf, .phasing-evidence.gz   the read pass's, with phasing
  depth/SAMPLE.regions.bed.gz, .mosdepth.summary.txt
  repeats/SAMPLE.repeats.vcf.gz              tandem repeats
  cnv/SAMPLE.vcf.gz, SAMPLE_cnv.bed.gz       CNVs
  methylation/SAMPLE.bedmethyl.gz            CpG methylation (bedMethyl), or
  methylation/SAMPLE.cpg.{combined,hap1,hap2}.bed.gz    CpG tables (hap1, hap2: with phasing, or haplotagged reads)
  qc/SAMPLE.qc.json, SAMPLE.qc.tsv           the QC summary
  qc/SAMPLE.read-stats.json                  the read pass's counts of the reads, which it summarises
```

The manifest's `outputs` lists each of these the run wrote (`qc_json`, `qc_tsv` and
`read_stats` for the QC's).

For small variants, tandem repeats and standard-mode CNVs, the VCF represents a
non-`N` IUPAC reference base with its first possible nucleotide in alphabetical order.
The input FASTA retains its original sequence. Reference blocks in a male gVCF are
split at pseudoautosomal boundaries. Non-PAR GT and PL use haploid cardinality.
With a target BED from `--regions`, gVCF reference blocks stay inside the target intervals.
A variant call keeps its GQ, the same with and without `--gvcf`; a reference block's is
taken from its projected likelihoods. If the selected genotype disagrees with the projected
likelihoods, PL is missing (and a reference block's GQ).

The read pass writes `SAMPLE.snf` whenever it calls SVs: it costs next to nothing (HG003 ONT
chr20: no measurable time, 20 MB of memory and a 0.9 MB file). For a cohort, run each sample,
then call the SVs of all of them together from their `.snf` files, and genotype the SNVs and
indels jointly from their gVCFs with GLnexus. GLnexus has no preset for Clair3's gVCFs (which
write `<NON_REF>`). Clair3 provides `clair3.yml`. Use the same reference and ordered contig
dictionary for every sample, unique SNF sample IDs and distinct gVCF basenames. The GLnexus
1.4.1 example covers human autosomes. Evaluate sex-chromosome joint calls with a ploidy-aware
joint caller. For targeted runs, give GLnexus a BED whose contigs are present in all input gVCF
headers.

```bash
lrpipe sv -i a/sv/A.snf b/sv/B.snf c/sv/C.snf -v cohort.sv.vcf --reference ref.fa
curl --fail --location -O https://www.bio8.cs.hku.hk/clair3_trio/config/clair3.yml
for dir in a b c; do
  gvcf="$PWD/$dir/snv/merge_output.gvcf.gz"
  ln -s "$gvcf" "$dir.g.vcf.gz"
  ln -s "$gvcf.tbi" "$dir.g.vcf.gz.tbi"
done
awk -F '\t' '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ {print $1 "\t0\t" $2}' ref.fa.fai > autosomes.bed
glnexus_cli --config clair3.yml --bed autosomes.bed a.g.vcf.gz b.g.vcf.gz c.g.vcf.gz > cohort.bcf
```

## QC summary

The qc stage writes `qc/SAMPLE.qc.json` and `qc/SAMPLE.qc.tsv` (one value a line, `name<TAB>value`,
`NA` for one not measured), and the run ends by printing its main values, as here for HG002 ONT
(the whole-genome CRAM, which has no supplementary alignments):

```
lrpipe: QC: 8,155,127 mapped reads (153.26 Gb), read N50 29.4 kb, mean length 18.8 kb, 97.71% mapped, of the alignments 0.00% supplementary, 0.00% secondary
lrpipe: QC: coverage 46.9x over the autosomes, chrX 24.2x, chrY 9.35x; sex male, as the index says
lrpipe: QC: PASS: 4,779,264 SNVs (Ti/Tv 1.70, het/hom 1.65), 889,493 indels (het/hom 1.82), 27,274 SVs (BND 191, DEL 11,479, DUP 222, INS 15,289, INV 93)
lrpipe: QC: 89.1% of 2,999,347 het SNVs (any FILTER) phased, block N50 1.88 Mb; 29,172,880 CpGs
```

Nothing is read again for it but the reads' index:

| Group | Values | From |
|---|---|---|
| Reads | mapped reads and their bases (the mapped yield), read N50, mean length and mean base quality; % mapped; secondary and supplementary alignments, and each as % of all alignments (primary, secondary and supplementary) | the read pass, which counts every record it reads (`lrpipe sv --read-stats`): the reads are the primary alignments of mapped reads, as cramino counts them, duplicates and QC failures included (reads sequenced, which long-read pipelines rarely mark). The unmapped reads' lengths are not read: the unplaced unmapped reads, which nothing reads, are counted from the index (a BAI's or CSI's count; a CRAM's containers of them), for % mapped |
| Coverage | mean depth over the autosomes (their bases over their length), of each autosome, chrX and chrY | the depth summary (MAPQ 20 or more, each read's span) |
| Sex | the sex the run called the sample with (`--sex`, or the index's), the index's, and whether they agree | the indexes, as the run estimated it as it started (over each of several `--reads`), as `lrpipe sex READS` reads it |
| SNVs and indels | PASS records, het (0/1, 1/2), hom ALT and haploid calls, het/hom; for the SNVs (REF and every ALT one base) transitions, transversions and Ti/Tv, by ALT allele. A record that is not an SNV is an indel, an MNP too (bcftools counts those apart; Clair3 calls none) | `snv/merge_output.vcf.gz` |
| SVs | PASS records by SVTYPE | `sv/SAMPLE.sv.vcf` |
| Phasing | het SNVs, those phased (in blocks of two or more variants), blocks and their N50 | `snv/phased_merge_output.vcf.gz`, by WhatsHap 2.8's `stats` (lrpipe-phase's), which counts the records of every FILTER, not only PASS, and leaves out multi-allelic sites: its het SNVs (`het_snvs_any_filter`) are not the PASS SNVs' `het`. Its `block_n50` is an NG50, the span of the shortest of the longest blocks that together span half the chromosomes with blocks |
| Methylation | CpGs reported: of all reads, and of hap1's and hap2's | the bedMethyl or the CpG tables |

The mean base quality is the arithmetic mean of the bases' Phred qualities, as samtools stats'
"average quality" (HG003 ONT chr20: 43.6). NanoPlot and cramino average the bases' error
probabilities instead (a read's Q is -10 log10 of its mean error), which the few low-quality
bases dominate, and report lower values for the same reads.

On whole-genome HG002 ONT it agreed with `bcftools stats` (records, Ti/Tv and genotypes;
bcftools also counts the 7,597 records of an SNV ALT beside an indel ALT among the SNVs, which
lrpipe counts as indels), `whatshap stats` (het and phased SNVs, blocks, NG50), `samtools
idxstats` (mapped and unmapped reads) and `zcat | cut -f1,2 | uniq` (CpGs), and took 14 s and
160 MB. The read pass's counting costs no measurable time (12 runs each, 16 threads, the machine
shared): on HG003 ONT chr20 (a 2.5 GB BAM) a median of 3.04 s against 3.01 s, with 2% more CPU
(summing the qualities) and 5 MB more memory; on HG002 ONT chr20 from the CRAM 2.41 s against
2.42 s, 0.5% more CPU.

What the run did not write is not measured, and the summary says why (`not_measured`): the
reads without the read pass (`--analyses snv` or `phase` alone), the coverage without depth, the
phasing without `phase`, the mean base quality of reads without qualities. From CRAM there is no
mean base quality: the read pass does not decode CRAM's qualities, which would make it 65%
slower (HG002 ONT chr20: 3.8 s against 2.3 s). An index that cannot be read for the unmapped
reads is a warning, with why.

In a run in regions (`--regions`), every value but the sex is of the regions, which the summary
names (`regions`: how many, merged, and their bases) and the run prints first: the read pass
counts each read that reaches them (its statistics say `restricted`), the coverage is of the
depth summary's regions' rows (`CONTIG_region`, the windows in them), and the variants, SVs and
CpGs are those called there. The unmapped reads are not measured: the unplaced ones are counted
only in the whole file, so a % mapped would compare the regions' reads with the file's.

A sex given with `--sex` that the index contradicts is a warning as the run starts, in the
summary (`warnings`), at the end of the run and in the manifest, once (the summary's estimate is
the run's): the index's chrX and chrY reads per base against the autosomes' say male, say, for a
run called as female (`lrpipe sex READS` prints them). A run of some contigs (a chr20 BAM)
cannot tell, and says so (`sex.note`).

`lrpipe qc --sample S --sex female|male --reads READS [--sex-reads F ...] [--read-stats F]
[--depth-summary F] [--snv-vcf F] [--phased-vcf F] [--sv-vcf F] [--bedmethyl F] [--cpg-combined
F --cpg-hap1 F --cpg-hap2 F] --json OUT.json --tsv OUT.tsv` writes it from any run's files.

## Setting up

```bash
lrpipe setup --genome GRCh38 --platform ont,hifi --dir /ref/lrpipe --accept-licences
lrpipe run --config /ref/lrpipe/site.toml --reads HG002.bam --sample HG002 --platform ont --output hg002
```

fetches a genome's files from a manifest built into lrpipe (`lrpipe/resources/*.toml`, each
file's URL, size, checksums and licence), checks each against its checksums, builds what is
missing and writes `DIR/site.toml` with absolute paths. A file with no published checksum
(Clair3's models, for now) is checked by its size alone (see **Unpinned files** below). A
rerun for another `--platform` adds its models and keeps the earlier platform's.

| File | GRCh38 | CHM13 | GRCh37 |
|---|---|---|---|
| reference (gunzipped; `.fai` built) | GRC's no-alt analysis set (NCBI) | T2T's `chm13v2.0_maskedY_rCRS` | hs37d5 (1000 Genomes) |
| `[sv] tandem_repeats` | Sniffles' TRF BED | pbsv's TRF BED | Sniffles' TRF BED |
| `[repeats] catalog` | adotto (Zenodo) and STRchive's disease loci, sorted together: `catalog.adotto_strchive.hg38.bed`, 937,204 loci | STRchive's disease loci | none |
| `[cnv] blacklist` | Spectre's | none (Spectre has none) | Spectre's refined one |
| `[cnv] metadata` | made from the reference by lrpipe-cnv | made | made |

and for each `--platform`, Clair3's PyTorch models from HKU into `DIR/models/NAME/`
(`pileup.pt`, `full_alignment.pt`), under Clair3's names: ONT's R10.4.1 E8.2 400 bps models
(`r1041_e82_400bps_{sup_v500,sup_v520,sup_v410,hac_v600,hac_v520,hac_v500,hac_v410}`) and
HiFi's `hifi_revio` and `hifi_sequel2`. `site.toml` sets `[snv] models = "DIR/models"`, so
each run takes the model for its reads' header (see [Clair3's model](#clair3s-model)); an
`ont_model` or `hifi_model` line, commented out, names a model to use whatever the header says.

* **Checked and resumable.** Each file is fetched to `NAME.part`, with progress every 10 s, and
  checked against the manifest's SHA-256 and the publisher's MD5 before it is renamed into
  place; one that does not match is fetched again from the start, three times at most. A
  transfer cut off, or a setup run again, resumes with an HTTP Range request (or from where a
  `file://` copy stopped), with `If-Range` and the version it began with (its `ETag` or
  `Last-Modified`, kept in `NAME.part.validator`), so a file that changed on the server is
  fetched whole; a part of an unknown version is resumed only where a checksum will check it.
  A transfer is resumed until five tries in a row bring no bytes; a name that does not resolve,
  a TLS or certificate failure, or a refusal (404, 403, 429) stops it at once. HTTPS is rustls
  with Mozilla's root certificates built in; `HTTPS_PROXY` is honoured.
* **Idempotent.** `DIR/.lrpipe-setup/state.json` records each file's size, time and SHA-256 as
  it was checked, and for a fetched file the size and checksums of what was fetched; a rerun
  skips what is present, unchanged and still as the manifest says (a manifest pinned since
  keeps it), checks a changed file again and fetches it again if it is wrong. One setup holds
  the directory at a time (`DIR/.lrpipe-setup/lock`). `site.toml` is written to a temporary
  file, loaded as a site configuration and then renamed into place, an earlier, different one
  kept as `site.toml.previous` (one the same is left alone). The last line setup prints is the
  file it wrote, or that held its configuration already. A `site.toml` changed since a setup wrote it (by its SHA-256 in the
  state, or for a setup from before lrpipe recorded it, of the one its recorded manifest and
  platforms make) is yours and is kept:
  the rerun writes its own as `site.toml.new` and says so, to merge your changes into and move
  into place. Replacing it would lose a hand-made change (a tool binary, `[snv] device`) but
  for the `.previous` copy, which the next rerun would replace too.
* **Verify.** `lrpipe setup --verify --dir DIR` reads every file again (the manifest and
  platforms are the recorded ones unless given), prints each one's state and exits with 1 if
  any is missing, wrong or not known to be good. A file the manifest has no checksum of passes
  if its size is the manifest's, and is listed and counted as `size only`, not as verified.
  `site.toml`'s comments give the command for its directory. A `site.toml` naming another
  reference than the setup's (one copied with its directory names the original's files) is a
  mismatch: setup run again in the copy writes the copy's own.
* **Unpinned files.** No checksum is published for Clair3's models (all of them, from HKU),
  hs37d5, pbsv's CHM13 TRF BED and STRchive's CHM13 catalogue: the manifest says so
  (`unpinned = "why"`). Such a file is fetched, resumed and fetched again like any other, but
  checked by its size alone, and is reported as checked by size only in the output, in
  `site.toml` and by `--verify`, never as verified. A maintainer records their SHA-256s with `lrpipe setup --pin --manifest
  lrpipe/resources/FILE.toml --dir SCRATCH`, which fetches every file of that manifest without
  one and writes it in; a directory set up before keeps the files that match.
* **Licences.** Each file's licence is in the manifest and in `site.toml`'s comments. ONT's
  models are under the Oxford Nanopore Technologies PLC. Public License 2.0: research use only,
  with Nanopore data only. That licence lets them be redistributed for research, so setup fetches
  them without asking. The
  TRF BEDs, all three pbsv's (GRCh38's and GRCh37's by way of Sniffles, whose
  `annotations/source.txt` names pbsv as their source) and so under PacBio's Software License
  Agreement, the only licence their origin gives, are fetched only with `--accept-licences`;
  without it, setup prints their terms and stops, and `--offline` prints them and marks each.
* **Offline.** `--offline` prints each file's URL, size, checksums and terms, the `curl -C -` (or
  `wget -c`) commands that fetch it to where setup looks for it, and the setup to run once the
  directory is on the machine without the internet: it then checks the files, builds the rest
  and fetches nothing. That setup names `--manifest FILE` when one was given, and `--bonefish
  PATH` (the bonefish there, which builds its index) when given.
* `--manifest FILE` sets up from a manifest of your own (with its models) instead.
* **bonefish.** `--bonefish PATH` builds bonefish's index of the reference,
  `DIR/REFERENCE.bonefish.raw` (`bonefish -d ... --idx-raw`; 8.5 GB for GRCh38), for runs from
  unaligned reads, and writes `[tools] bonefish` and `[align] index` into `site.toml`. PATH must
  be bonefish (its `--version` printing `bonefish VERSION`), checked before anything is fetched,
  as `lrpipe run` checks `[tools] bonefish`. It is
  built again only for another reference or another bonefish version; a later setup without
  `--bonefish` keeps both, and `--verify` checks the index.

## Site configuration

One TOML file per site; relative paths are taken from its directory. The stages run with the
binary running the pipeline unless `[tools] lrpipe` names another. Every section is optional,
but a run calling SNVs needs a model from `[snv]` (`models`, or the platform's `ont_model` or
`hifi_model`): without `[repeats]` no repeats are genotyped.

```toml
reference = "/ref/GRCh38_no_alt.fa"        # with .fai; also decodes CRAM

[tools]
lrpipe = "/opt/bin/lrpipe"                 # optional: the binary the stages run with
bonefish = "/opt/bin/bonefish"             # optional: aligns unaligned reads (default: beside
                                           # lrpipe, else on PATH)

[align]
index = "/ref/GRCh38_no_alt.fa.bonefish.raw"   # optional: bonefish's index of the reference

[snv]
models = "/opt/clair3/bin/models"          # Clair3's models by their standard names: each run
                                           # takes the one for its reads' header
# hifi_model = "/models/hifi_revio"        # optional: a model for every HiFi run instead
# ont_model = "/models/r1041_e82_400bps_sup_v500"   # and for every ONT run
device = "cuda:0"                          # or "cpu" (default), or "metal" (an Apple GPU;
                                           # a build with `--features metal`)
nvrtc = "/opt/cuda/lib64/libnvrtc.so.12"   # optional: for GPU calling, if not on the library path

[sv]
tandem_repeats = "/ref/human_GRCh38_no_alt_analysis_set.trf.bed"

[repeats]
catalog = "/ref/catalog.adotto_strchive.hg38.bed"

[cnv]
blacklist = "/ref/grch38_blacklist_spectre.bed"
metadata = "/ref/grch38_reference.mdr"     # optional: made here from the reference if absent
                                           # (name it after its reference: it is not remade when
                                           # the reference changes); without it, each run makes
                                           # one for its reference

[methylation]
report = "auto"                            # or "bedmethyl", "cpg", "none"
estimate = "count"                         # CpG tables: or "model" (experimental: lrpipe-methyl's model, fitted on one Revio run against EM-seq)
```

`[snv]` needs `models`, or `ont_model` or `hifi_model` for the platforms the site runs (or each
run's `--model`); with both, `ont_model` and `hifi_model` take precedence, checked against each
run's header ([Clair3's model](#clair3s-model)).

Only the files the chosen `--analyses` read need to exist: a run of `--analyses snv` does not
need the repeat catalogue, an ONT run does not need the HiFi model, and a run from aligned reads
needs neither bonefish nor its index.

## Building

```bash
cargo build --release
```

Linux and macOS (child processes' CPU time and peak memory come from `wait4`). On macOS lrpipe
cannot confine a run to `--threads` cores, so lrpipe-snv and the read pass split them, and the
memory available is the share macOS's memory pressure counts as free (`kern.memorystatus_level`).

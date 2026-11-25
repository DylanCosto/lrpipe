# Changelog

## 0.3.0

The first public release.

- Release bundles for Linux x86-64, Linux aarch64 and macOS on Apple silicon (with Apple GPU
  calling), built by the release workflow, each with bonefish 0.5.2. Clair3's models are a
  separate download, `lrpipe-0.3.0-models.tar.gz`.
- The crates live under `crates/`; the README is shorter, with the benchmarks.
- With `--regions` and `--gvcf`, gVCF reference blocks stay inside the BED targets. Between
  disjoint targets they covered the gaps, more than 100 Mb of confident reference calls on chrX
  with two 30 kb targets, as Clair3 writes them (#53). Variant calls, and runs without a BED, are
  unchanged.

## 0.2.2

- The bundle carries bonefish 0.5.2, which keeps a read's `MM`/`ML` modification tags only on
  records whose SEQ is the whole read. With tagged FASTQ or unaligned BAM input, bonefish 0.5.1
  copied them onto hard-clipped supplementary and SEQ-less secondary records, where they pointed
  past or into the wrong bases.
- The release also has `lrpipe-0.2.2-models.tar.gz`: Clair3's nine models, for `[snv] models`
  without running `lrpipe setup` (ONT's models for research use with Nanopore data only, as its
  `NOTICE.md` says).
- `lrpipe sv` names a bgzipped VCF's index after the VCF, `calls.vcf.bgz.tbi` (or `.csi`).
  It wrote `calls.vcf.gz.tbi` for a `.vcf.bgz`, as Sniffles2 does: bcftools did not find it,
  and it replaced the index of a `calls.vcf.gz` beside it, whose indexed lookups then found
  nothing (#52). Names ending in `.gz` are unchanged.
- `lrpipe snv --pileup-only --gvcf` is refused before the run starts. It exited 0 without the
  gVCF wherever the pileup had calls, as Clair3 does (#54).
- `lrpipe snv --device cuda` no longer collapses at high depth: the pileup network takes the CPU
  path's inputs, with their high-depth scaling, instead of Clair3's GPU-mode inputs, whose counts
  wrap past 127. On K-12 HiFi at 295x CUDA called SNPs at F1 42.5% and indels at 1.4%, against
  the CPU's 96.4% and 97.6% (#45).
- `lrpipe phase` (and `lrpipe-phase`'s `haplotag`, `genotype` and `haplotagphase`) refuses a
  `--reference` whose `.fai` is another FASTA's index, before reading or writing anything, as
  `lrpipe run` and `lrpipe cnv` already did. It used to exit 0 and phase 22 of 66 records
  differently (#43).
- `lrpipe run` and `lrpipe cnv` check that the reference's `.fai` is this FASTA's before they
  write anything: a rerun with another FASTA's index used to replace a completed run's manifest
  with a failed one, and `lrpipe cnv CNVCaller` accepted it. The gVCF and VCF chunk sorter skips
  macOS AppleDouble files (`._` names with that file signature), which on exFAT stopped it with
  an invalid UTF-8 error.
- A resume reruns a stage whose outputs were damaged since it finished (#46). Its stamp now
  holds a SHA-256 of each output up to 4 GB (about two seconds to check for a whole genome's VCFs,
  gVCF, BEDs and QC JSON on an NVMe disk), and for larger ones (the BAMs) their size,
  modification time and 33 evenly spaced 1 MiB pieces, first and last among them, and a BAM
  must end with its BGZF end-of-file block. Before, only the outputs' presence and indexes were
  checked, so a flipped bit or a truncated QC JSON was taken as done. A stamp of 0.2.1 or before
  is taken once and written anew with its outputs' state.
- A resumed run no longer leaves an old index beside a rewritten repeat VCF
  (`repeats/SAMPLE.repeats.vcf.gz`): the read pass now indexes the repeat VCF it writes, both
  writers remove an index of the other kind (`.tbi`/`.csi`), and the index is one of the stage's
  outputs, so a missing one reruns it. Before, a panel run then a whole-genome run in the same
  directory kept the panel's index, and position lookups found 1 of 10 records (#42).
- `lrpipe setup` fetches ONT's Clair3 models without `--accept-licences`. Their licence, now
  the Oxford Nanopore Technologies PLC. Public License 2.0, allows redistributing them for
  research use with Nanopore data, and the v0.2.1 release carries them as
  `lrpipe-0.2.1-models.tar.gz`. The manifest names version 2.0 and its current URL
  (`LICENSE.txt`; `LICENCE.txt` is gone). The TRF BEDs, under PacBio's licence, still need
  `--accept-licences`.
- A male sample on a reference whose chrX is not GRCh38's, GRCh37's or CHM13's (so whose
  pseudoautosomal regions lrpipe does not know) no longer fails `lrpipe run` at its tandem
  repeats after the dry run passed (#44). `lrpipe tr` and `lrpipe sv --repeats` count all of
  chrX as one copy there, as lrpipe-cnv already did, and warn.

## 0.2.1

- `lrpipe snv` on Apple GPUs: `--device metal` (`[snv] device = "metal"`), in a build with
  `--features metal` (macOS). The networks run as Metal kernels compiled at start-up by the
  Metal framework, called from Rust (`objc2-metal`), in `float32` with the CPU path's inputs,
  so the calls match the CPU's but where a sum's order tips a close one (HG002 chr20 at 10x and 46x: every record the same,
  QUALs within 0.01). On an M4 Pro, low-coverage (10x) HiFi chr20 took 52 s against 84 s with
  `--approximate` on 12 CPU threads (pileup calling 21 s against 54 s, its 8-bit network
  turned down by the run's check) and 19 CPU-s against 902. CPU and CUDA calling are unchanged.

- `lrpipe run` rejects an aligned BAM or CRAM whose read groups declare multiple samples. For
  one aligned file with verified same-donor labels, give `--sample-alias` once for each label other
  than `--sample`.
- `lrpipe tr` (and the repeats of `lrpipe run`): three faults in how reads' haplotypes and
  alignments were used.
  - Haplotypes of different phase sets were put together. A read's `PS` (its tag, or the
    haplotag list's) was dropped, so where two phase blocks meet, one block's H1 and the other's
    H2 were one group and both groups held both alleles. Only the haplotypes of the locus's
    most common phase set now count; the other reads go to the allele of the nearest length.
  - Haplotags that did not tell two alleles apart gave a homozygous call without looking for
    an expansion: with tags that mix the alleles (a switch error, another run's tags), a
    heterozygous expansion's reads were in both haplotypes and the normal allele was called
    twice. The expansion's much longer reads are now looked for, as without tags.
  - A read counted twice when its primary and a supplementary alignment both gave the locus an
    allele, as the two sides of an expansion the aligner split do where supplementary
    alignments keep their bases (`-Y`): two molecules were enough for an expansion that needs
    three. A read now gives a locus one observation, its primary alignment's.
  On HG002 ONT chr1 and chr20 (the haplotagged CRAM, GIAB TR v1.0.1 Tier 1), 147 of 99,178
  records change and 85.00% of scored regions are exactly right where 84.99% were.
- `lrpipe run --bacterial` calls haploid SNVs and small indels on the supplied reference contigs
  and skips human sex inference. `--gvcf` also writes haploid reference blocks. gVCF headers
  preserve the sample name on references with versioned contig names such as `NC_000913.3`.
  The caller finalizes indexed VCF and gVCF outputs when every call comes from the pileup.
- `lrpipe snv` writes gVCF reference blocks from the calling chunks only. The CPU precision probe
  previously added overlapping blocks for its sampled subregions.

- The release bundle holds bonefish (0.5.1), so `lrpipe run` goes from FASTQ or an unaligned BAM
  to VCF with nothing else installed; `lrpipe setup` builds bonefish's index of the reference
  with it (`--no-bonefish-index` leaves it out).
- `lrpipe run`'s display lists only the run's stages, with its analyses in the header and the
  read pass named for what it computes (`lrpipe-sv: SVs, depth`); a stage of an analysis not
  asked for has no row. The align stage's bar follows bonefish's own progress (reads aligned,
  time left), from bonefish 0.5.1's `--progress`, which the bundle's bonefish has.
- Found by running every tool on edge cases:
  - `lrpipe snv` checks the model's checkpoints, the BED, the sample name and the CUDA device
    before it removes a previous run's results; a corrupt `pileup.pt`, a BED with start after
    end, an empty `--sample-name` or `--device cuda:3` on a one-GPU machine deleted them first.
  - `lrpipe run` refuses an output directory where files no lrpipe run wrote are in its way
    (`SAMPLE.vcf.gz`, `manifest.json`, `qc/`, `logs/`, a stage's directory), as it refused one in
    `snv/` alone; any `manifest.json` no longer counts as an earlier run's.
  - `lrpipe qc` refuses `--json` or `--tsv` naming an input (or a link to one) or each other, and
    VCFs of another sample than `--sample`.
  - A `.fai` that is not its FASTA's index (copied from another FASTA of the assembly) is
    refused by every tool that opens a reference; `lrpipe snv` called 48,824 variants in 200 kb
    from one.
  - `lrpipe snv` and `lrpipe run` warn when the reads' index is older than the reads.
  - `lrpipe sv` refuses a reference whose contigs differ in length from the reads' (it wrote
    empty REFs), warns of contigs the reference lacks, and refuses `--regions`, `--region` or
    `--contig` that select nothing, where it called the whole genome or nothing.
  - `lrpipe tr` refuses a reference whose contigs differ in length from the reads', catalog
    loci past their contig's end, zero-length loci, motifs that are not bases or IUPAC codes and
    duplicate locus IDs, and reads a locus given twice once.
  - `lrpipe methyl pileup` estimates the same thresholds for a BAM and a CRAM of the same reads
    (and whatever `--window`), counts a call tied between the canonical base and a modification as
    modkit does with `--modified-bases`, and refuses repeated or conflicting filter options and `--modified-bases C`
    without `--combine-mods`, as modkit does.
  - The sex from the index is refused for reads concentrated on one autosome (a panel aligned
    genome-wide read as male from chrY at 10.5 times the autosomes), and chrX or chrY above 1.6
    times the autosomes fits no sex.
  - BED files with UCSC `track` and `browser` lines are read by `lrpipe snv`, `lrpipe tr` and
    `lrpipe cnv`'s blacklist.
- `lrpipe snv`, found in review:
  - With `--gender male` (as `lrpipe run` gives a male sample) and no full-alignment calls, chrX
    and chrY outside the PARs kept the pileup's diploid genotypes; they are now haploid, as with
    full-alignment calls.
  - An indel anchored on an ambiguous reference base (an IUPAC code such as R) had REF and ALT
    first bases that disagreed (`REF=ATT ALT=N`); the ALT's anchor is now written as the REF's.
  - With `--pileup_only --gvcf` and no pileup records, the gVCF's reference blocks were filtered
    as calls: all marked `LowQual`, and haploid ones dropped.
  - `--remove_intermediate_dir` now also removes `tmp/` when no full-alignment calls are made.
- `lrpipe sv` combines SNF files whose runs read different contigs of one reference (a panel's
  and a whole genome's), calling SVs on every contig any of them has; files that give a contig
  different lengths are refused. Since 0.2.0 any two files with different contig lists were.
- Dependencies: `cargo audit` reports no vulnerabilities. ndarray-npy no longer pulls in `zip`
  (lrpipe reads and writes no `.npz`), and pprof's flame graph feature, which pulled in `inferno`
  and `quick-xml`, is dropped. Builds with `--features profile` write folded stacks
  (`PREFIX.folded`, for `inferno-flamegraph`) in place of an SVG.
- pasta moves to `3050cd1`: CRAM codec streams that are corrupt fail instead of panicking, faster
  rANS and fqzcomp decoding, and htscodecs 1.6.7's read name tokens and fqzcomp parameters.
- The release bundle carries the license files of every crate compiled into its binaries, in
  `licenses/third-party/` with an index of their terms.
- `tools/check.sh` compares lrpipe-cnv with Spectre 0.2.1's reference runs where they are found
  (`$SPECTRE_BENCH`). These tests passed without running when the variable was unset.

## 0.2.0

- `lrpipe run` on a terminal shows its stages full-screen (`--display`, `auto` by default): every
  stage with its tool and a progress bar, the overall progress and time left, the memory in use
  against the budget, and the end of a stage's log. lrpipe-snv reports its pileup chunks and
  full-alignment jobs, and the read pass its contigs, to a file the run names in
  `LRPIPE_PROGRESS`; other stages' bars follow their last run's time. Off a terminal, or with
  `--display plain`, the run prints its lines as before, and the display prints them too when it
  ends.
- pasta `8e2d1e9` (0.9.1 and later), from `fe18544`. BGZF goes through zflate, in Rust, which
  replaced pasta's libdeflate, so no C library is compiled for it; BGZF speed with zflate is not
  yet measured. pasta also brings a multithreaded CRAM reader (`multithreaded_records`, not yet
  used here), HTSlib's `d` tag type, and errors instead of panics or aborts on corrupt BAM
  records and indexes and on reads without a sequence in the pileup.
- Reads with a `d` (double) tag, which samtools writes to BAM. `lrpipe snv` stopped reading a
  read's tags at one, so it missed a move table (`mv`) after it; `lrpipe sv` failed with
  "unknown tag type d". Both now read it. `lrpipe sv` also fails with an error,
  instead of a crash, on a read whose tags are cut short or hold an array of unknown type.
- Semver-compatible dependency updates (clap, flate2, log, lz4_flex, toml, ureq and others).

## 0.1.9

- SV genotypes of long insertions. Sniffles2 genotypes an insertion of 2.5 kb or more from its
  support scaled up by its length (about 2x at 2.7 kb), for reads that end inside it, though it
  already counts those reads, so heterozygous long insertions were called 1/1. lrpipe-sv's
  `--accuracy`, and so `lrpipe run`, now genotypes them from their supporting reads
  (`--read-support-genotypes`); VAF and every filter keep the scaled support, so the same calls
  pass. GIAB v5.0q genotype concordance on HG002 rises from 0.895 to 0.900 (HiFi) and 0.888 to
  0.896 (ONT), with precision and recall unchanged; on the HG002, HG003, HG004 ONT trio,
  Mendelian errors among insertions of 1 kb or more fall from 14.7% to 6.3% (all SVs: 5.05% to
  4.54%). `--sv-mode sniffles2` is unchanged.
- `lrpipe run` says when a CRAM's `.crai` is the index of another version of the file (its
  first, middle and last entries must each point where a container starts), as ONT's published
  indexes of their GIAB 2025.01 HG002 and HG004 CRAMs are (63 and 68 bytes off). Such an index
  failed partway into the run as "container header checksum mismatch", as if the CRAM were
  corrupt.

## 0.1.8

- `lrpipe tr` calls large repeat expansions. Aligners write the reads of an expansion of more
  than about 1 kb as alignments that stop in the repeat, the rest soft-clipped, and lrpipe-tr
  used only reads whose alignment spanned the locus, so it missed every such expansion (DMPK,
  FMR1 full mutations, FXN, RFC1) and called the normal allele twice. A clipped read now counts
  when its clipped bases start as the repeat and hold the far flank; a few much longer reads
  are an expanded allele (3 reads and 5%: long molecules are sampled less); and at most 200
  reads per locus are used spread evenly over them, not the first 200. On PacBio's PureTarget
  data of 33 Coriell cell lines, lrpipe-tr matches Coriell's repeat counts on 24 of 38 samples
  (0.1.7: 11 of 34, no call at 4; TRGT 22) and TRGT's calls on 41 of 45. Whole-genome HG002
  accuracy and speed are unchanged (HiFi 94.24%, ONT 85.79% exactly right against GIAB), and
  two large VNTR expansions are now called as GIAB has them.

## 0.1.7

- `lrpipe run` calls SVs with lrpipe-sv's `--accuracy` (a repeat expansion's two alleles called
  apart, a read's indels joined within 50 bp), which finds about 1,600 more of GIAB v5.0q's
  SVs on HG002 and genotypes them better: through the pipeline's phasing, strict F1 0.835 where
  Sniffles2's calls score 0.807 (HiFi) and 0.836 where they score 0.810 (ONT), genotype
  concordance 0.895 and 0.888 where 0.839 and 0.834, at the same speed. `--sv-mode sniffles2`
  keeps Sniffles2's calls, record for record; `lrpipe sv` on its own still writes them unless
  given `--accuracy`. The read pass's command changes, so it reruns once in an output
  directory from an earlier version.
- `lrpipe tr --platform ont` (and `lrpipe run` on Nanopore reads): an allele whose longest
  single-base run has 10 bases or more takes the length its reads make most likely under a
  model of Nanopore's homopolymer errors (bases missed more often than added, more so the
  longer the run), not the median read's, and that run is resized to it. On HG002 ONT
  against the GIAB tandem repeat benchmark, 85.9% of regions are exactly right where 85.0%
  were; on the chromosomes the model was not fitted on (all but chr1), 86.0% where 85.1%.
  HiFi output is unchanged.
- `lrpipe tr`: an insertion at a boundary two catalog loci share counts in one of them, split
  where it best continues each one's repeat; it counted in both. Catalogs that split a locus
  into adjacent repeats, as ExpansionHunter's do, got both parts' expansions wrong: on its
  pathogenic catalog, HG002 ATXN7's +9 bp read as +9 in both parts (now +6 and +3, as GIAB's
  total and ExpansionHunter have it), FXN's GAA expansion was also counted in the A tail before
  it, and NOP56's parts summed to GIAB's change only now. lrpipe's own catalog has no shared
  boundaries, so `lrpipe run` writes the same repeats as before (identical on HG002 HiFi and
  ONT, whole genome).
- pasta is updated from `1f68e79` to its first public release, 0.2.0 (`fe18544`), so the
  workspace builds from source without access to a private repository, and GitHub's check needs
  no deploy key. Every output of whole-genome runs of HG002 HiFi (BAM) and ONT (CRAM) is
  identical to 0.1.6's, record for record, and the release gate passes; on the ONT CRAM the run
  takes 3% less CPU time than with the previous pasta. From pasta: CRAM slices without a
  reference checksum are checked against the header's `@SQ` LN (and M5 where the length differs,
  or where the reference has a `.dict`), so a wrong reference is an error instead of wrong
  bases; headers that htslib reads (duplicate `@PG` IDs, empty values) are read; and the merge
  of several `--reads` adds its `@PG` line as samtools 1.24 does (`lrpipe-merge-reads`, then
  `lrpipe-merge-reads.1` after the end of each further chain, or one line after the last when
  no line follows another), where it named the others `lrpipe-merge-reads-LAST`.
- haplotag finds the contigs with reads of an indexed CRAM from its index instead of decoding
  each contig's first slice. With a pasta before 0.2.0 that decode also hashed every contig of
  the reference (the check above), 5-6 s on one thread; the haplotag list of chr17-chr22 from
  the HG002 ONT CRAM takes 29 s, where it took 36 s with 0.1.6 and 43 s with that pasta alone.
- lrpipe-snv's CRAM queries share each contig's reference sequence instead of reading (and,
  with a pasta before 0.2.0, hashing) it again for each chunk (lrpipe-snv's CHANGELOG): 1,007
  CPU-s on chr17-chr22 of the ONT CRAM, against 1,000 with 0.1.6 and 1,050 with that pasta
  alone.

## 0.1.6

- `lrpipe run --haplotag off` phases the SNVs and indels but does not haplotag the reads from
  them: the SVs, tandem repeats and CpG tables run beside lrpipe-snv and follow the reads' own
  `HP` tags, as in 0.1.5, saving the haplotag stage's minutes (3 on HG002 HiFi, 4 to 5 on the
  ONT CRAM). Reads tagged by another pipeline are warned of, as their haplotypes need not be the
  phased SNVs'. `--haplotagged-reads` needs it on.

- `lrpipe run` leaves methylation out (with a warning) when the aligned reads carry no base
  modification tags (`MM`/`ML`) on the reads it samples, and refuses `--analyses ...methyl`
  for them: such reads gave empty CpG tables and bedMethyl with an exit status of 0, as the
  HG002 HiFi benchmark BAM did. Found by the whole-genome rerun (BENCHMARKS.md).

- **`lrpipe compare ... | head` ends quietly**: a reader that stopped reading made it panic
  ("failed printing to stdout: Broken pipe") and abort; it now ends as other Unix tools do, as
  do `lrpipe setup --verify` and `--offline`.
- **`lrpipe setup` fixes**: `--offline` with `--manifest` printed a setup to run with `--genome
  NAME` (the manifest's genome, which is no `--genome`); it now names `--manifest FILE`, and
  takes `--bonefish`, which it refused, naming it too. The last line setup prints is the file it
  wrote: `site.toml.new` when it kept a `site.toml` changed by hand (it printed `site.toml`), and
  a rerun that changes nothing says so rather than "wrote". `--bonefish /bin/true` was taken
  (its `--version` ran) and failed after fetching, renaming an index it had not written: a
  bonefish whose `--version` is not `bonefish VERSION` is refused before anything is fetched,
  and by `lrpipe run` as `[tools] bonefish`. `--verify` of a directory copied with its
  `site.toml`, which names the original's files, called it valid: a `site.toml` naming another
  reference than the setup's is a mismatch.
- **Memory below a stage's peak is warned of**: with `--memory 4` (or that little available)
  every stage ran alone without a word, though lrpipe-snv expects 14.5 GiB. The run now warns,
  naming the stages expected to use more and that each runs with no other stage beside it and
  may run out of memory.
- **A failed or stopped stage leaves no empty outputs**: Ctrl-C while bonefish or the
  haplotagging had opened their outputs left 0-byte BAMs (`align/SAMPLE.bam`,
  `phase/SAMPLE.haplotagged.bam`, the haplotag list), and the align stage's scratch directory.
  The 0-byte outputs and the scratch directory of a stage that failed or was stopped are
  removed; resuming is as before.
- **Clearer model messages**: a `[snv] models` that does not exist is said not to (it was "holds
  no Clair3 model for ...: add r1041_e82_400bps_sup_v500"), and one model's own directory given
  as `[snv] models` is said to be one, naming the directory to give instead; a site
  configuration without `[snv]` parses (it stopped at the raw "missing field `snv`"), a run
  calling SNVs then saying it needs `[snv] models` or `[snv] ont_model`. A given model the header
  contradicts (a SUP model for HAC reads) is no longer said to be "checked against the header"
  after the warning, but "which the header says the reads are not (see the warning)"; HiFi's
  model for reads of another platform let through with `--model-mismatch warn` is not called
  "the nearest"; and why a model other than the reads' is used is given once, in the warning,
  the model's line pointing to it. A FASTQ run's model is chosen from the basecall model in
  dorado's `RG:Z` tags (`RG:Z:RUNID_dna_r10.4.1_e8.2_400bps_sup@v5.0.0`), and a model given for
  FASTQ without them is said to be unchecked as "FASTQ has no header saying what made the
  reads", not for want of `@RG` fields. A network file cut short (a zip archive without its end
  record, as an interrupted copy leaves `pileup.pt`) is refused at planning, not in lrpipe-snv.
- **Smaller wording fixes**: an empty FASTQ is "empty, with no reads", not "not a format lrpipe
  reads"; reads without HP tags are "no HP tag in any of 152 mapped reads sampled".
- **`--sex auto` of unaligned reads is refused before aligning when it cannot work**: the sex
  is estimated from the aligned reads' index against the autosomes, so reads aligned to a
  reference of some contigs (the fixture's chr20 window) cannot give it, which the run found
  only after aligning, saying "the index of align/SAMPLE.bam has no autosomes". A reference
  whose `.fai` has fewer than 11 of the 22 autosomes, or no chrX, is now refused before bonefish
  runs (`the reference REF has 1 of the 22 autosomes (chr20) and no chrX: ...; give --sex`).
  Whether the reads cover the autosomes is known only once they are aligned, and is said as for
  aligned reads: a header of chr20 alone and a whole genome's header with reads on chr20 alone
  both give "has reads on 1 of the 22 autosomes (chr20)" (a header without chrX was taken for
  one naming its contigs `1`, `2`, ...).
- **The leftovers warning names earlier runs' BAMs**: `align/SAMPLE.bam`, `reads/SAMPLE.bam` and
  `phase/SAMPLE.haplotagged.bam` left by a run of other input (unaligned reads, then aligned) or
  options (`--haplotagged-reads`, then a run without it) were kept silently. They are named with
  the other outputs this run does not write, and kept (another tool may read them); haplotagged
  reads older than this run's haplotag list are said to be an earlier phasing's.
- **A moved or copied output directory resumes**: the stages' stamps held the output
  directory's absolute path, so a directory moved (or copied with `cp -a`) reran every stage,
  the alignment included. Paths in the output directory, in a stage's arguments and inputs, are
  now stamped as paths under it; everything else in the stamp is as before. This changes every
  stamp once: those written by lrpipe 0.1.5 and before are still taken while the directory is
  where they were written, and are written anew as the stage is skipped, so an upgraded run
  resumes without rerunning; a directory moved before its first run with this lrpipe reruns
  once.
- **The same reads in two `--reads` are refused**: two files sharing more than 90% of the reads
  sampled from the smaller (a BAM and its CRAM, a file and its copy, a file and a merge holding
  it) were only warned of, and merged into a BAM of twice the depth. They are now an error; a
  smaller overlap is still a warning. Either names the file of the smaller sample, whose share
  it is, and a warning of a third or more shared asks whether one file is a merge holding the
  other (on a whole genome a file and such a merge share about half their samples).
- **`lrpipe run` records the align stage as this run's**: a run that aligned its reads recorded
  the align stage as `skipped` in its manifest, and printed `align: already done, skipped` after
  `align: done`, as the second phase (the stages planned from the aligned reads) judged it again.
  The record of the alignment the run did is kept as it is, `succeeded`, with its measurements.
- **`--dry-run` refuses what the run refuses**: a dry run of aligned reads did not check the
  site configuration's files (a missing `[sv] tandem_repeats`, `[cnv] blacklist` or `[repeats]
  catalog`) or the stages' programs, and no dry run checked the output directory, so it exited 0
  where the run exited 2. It now checks them, and refuses a directory holding another sample's
  run or in use by another run, as the run does. A dry run marks each stage a run would skip
  (`: done, would be skipped`), prints the snv model's line for unaligned reads too, and quotes
  the commands as a shell takes them (bonefish's `-R '@RG\tID:...'`). A missing file of the site
  configuration is named by its key (`the site configuration's [repeats] catalog ... does not
  exist`, not `repeat_catalog`), and a missing `[align] index` is refused before bonefish runs,
  not reported as bonefish's "cannot open reference". A configured `[cnv] metadata` not made yet
  no longer stops a run of unaligned reads before aligning: the run makes it, as for aligned
  reads.
- **`lrpipe compare` compares the BAMs**: it stopped with "stream did not contain valid UTF-8"
  on a run of unaligned reads (`align/SAMPLE.bam`) or several aligned files
  (`reads/SAMPLE.bam`). A BAM is now compared by its alignments, each with every field and tag,
  and its header lines but for the `@PG` command lines, those at one place as a multiset, both
  BAMs read once side by side; one without its BGZF end-of-file block (cut short) is an error.
  The haplotagged reads, which were not compared, are compared so too; an output named as text
  that is not text is an error, and any other that is not text is compared byte for byte.
- **A QC summary for each run** (F9): a new last stage, **qc** (`lrpipe qc`), writes
  `qc/SAMPLE.qc.json` and `qc/SAMPLE.qc.tsv`, which the manifest lists, and the run ends by
  printing their main values. Reads: mapped reads, mapped yield, read N50, mean length and base
  quality, % mapped, and the secondary and supplementary alignments, counted by the read pass from
  the records it reads anyway (`lrpipe sv --read-stats`, a new option: no measurable time on HG003
  ONT chr20, 2% more CPU; from CRAM without the mean quality, whose decoding would make the pass
  65% slower), with the unplaced unmapped reads from the index (a BAI's count, or a CRAM's
  containers; an index that cannot be read is a warning). Coverage over the autosomes, each
  autosome, chrX and chrY from the depth summary; the sex the run called the sample with against
  the index's; PASS SNVs and indels (het, hom, het/hom, Ti/Tv) and SVs by type; the phasing (het
  SNVs phased, blocks, block NG50) by WhatsHap's `stats` (lrpipe-phase's, now also
  `lrpipe_phase::stats::totals`); CpGs reported, combined and per haplotype. On whole-genome HG002
  ONT each value agreed with bcftools, WhatsHap and samtools, in 14 s. What a run does not write
  is said to be not measured. A `--sex` the index contradicts is now a warning as the run starts
  and at its end, in the summary and in the manifest, once (with several `--reads`, the summary's
  estimate is the run's, over each of them). The qc stage runs last, after every other stage, and
  reruns when a file it reads changes; its failure does not fail the run (the manifest records it,
  and the run warns that there is no summary). The mean base quality is the arithmetic mean of the
  bases' Phred qualities (samtools stats' "average quality"; NanoPlot and cramino report lower
  values). A run resumed from an earlier version reruns the read pass, whose command gains
  `--read-stats`, and the stages after it. In a run in regions (`--regions`) every value but the
  sex is of the regions, which the summary names: the read pass counts each read reaching them
  (`restricted`), the coverage is of the regions' rows of the depth summary, and the unmapped
  reads, counted in the whole file, are not measured.
- **`lrpipe run --regions BED`** (F7): every analysis restricted to a BED's regions (a targeted
  panel, adaptive sampling's targets, one chromosome for a test run), each by its tool's own
  region option: lrpipe-snv's `--bed-fn` (Clair3's `--bed_fn`; on the ONT fixture its PASS calls
  are Clair3 2.0.3's with the same BED), the read pass's `--regions` (Sniffles2's; its calls are
  Sniffles 2.8.1's with the same regions) for the SVs, the whole genome's 1 kb windows cut to the
  regions for the depth, counted from the regions' reads alone in the pass (`lrpipe sv
  --depth-in-regions`: the windows' values are mosdepth's, the summary has the regions' rows
  alone), the catalogue's loci a region reaches for the tandem repeats, `--include-bed` for the
  bedMethyl and the CpG tables (a CpG counted when its C is in a region; modkit's
  `--include-bed` counts each strand's position on its own), and the reads in the regions for
  the haplotagging. Regions that are whole contigs are restricted by contig (`--ctg-name`,
  `--contig`), so a chr20 run gives a whole-genome run's chr20 records for the depth, SVs and
  repeats' loci; a few SNVs near lrpipe-snv's quality cut-offs, which it draws from every contig
  it calls (as Clair3), differ, and so may their phasing and the haplotypes that follow. CNVs
  are left out of the default analyses, saying why (lrpipe-cnv normalises over the whole
  genome), and `--analyses cnv` with `--regions` is an error. The BED is checked before planning
  (a contig, start and end on each line, the reference's and the reads' contigs, ends within
  them; a byte order mark dropped, a zero-length row taken as one base, a contig without reads
  warned of), regions on contigs a whole-genome run calls nothing on (alt, decoy, chrM, EBV) are
  left out with a warning, and regions covering every contig such a run calls are the whole
  genome. The rest is sorted in the reference's order and merged into
  `OUT/.lrpipe/regions.bed`, whose contents are part of the restricted stages' fingerprints:
  other regions rerun them, the same in another order do not. The sex is still read from the
  index; unaligned reads are all aligned and several files all merged. New tool options for it:
  `lrpipe tr --regions`, `lrpipe sv --repeats-regions`, `--bedmethyl-regions` and
  `--depth-in-regions`, `lrpipe methyl pileup --include-bed` and `cpg --include-bed`, and
  `lrpipe haplotag --regions-bed`, which writes a read reaching two regions once, where
  WhatsHap's `--regions` writes it again for each, out of order.
- **lrpipe-snv: no record joined to the `#CHROM` line** (HG002 ONT chr20 from a CRAM, on the GPU
  with phasing): merging each contig as its calls finish, a contig without full-alignment calls
  took its header from SortVcf's header-only VCF, which ends without a newline after `#CHROM`,
  so the first pileup record of the next contig joined that line in `merge_output.vcf.gz` (one
  record lost, and lrpipe-cnv stopped with "expected exactly one sample"). The merge now ends
  each header line with a newline. `ci/parity.sh` still agrees on every record.
- **`lrpipe run` haplotags the reads from its own phasing** (F1b): when it phases the SNVs and
  indels, the SVs are phased like them, the tandem repeats' alleles are split by the same
  haplotypes and the HiFi CpG tables are written per haplotype, from reads without `HP` tags
  too, and reads haplotagged by another pipeline get the run's haplotypes rather than their own
  tags' (F1a's warning about those goes). A new **haplotag**
  stage after lrpipe-snv writes `phase/SAMPLE.haplotags.tsv.gz` (`lrpipe haplotag --list-only`:
  WhatsHap's haplotag list, the same as WhatsHap 2.8 writes) instead of a second BAM;
  `--haplotagged-reads` also writes `phase/SAMPLE.haplotagged.bam` with its index. The read pass
  still runs beside lrpipe-snv, calling the SVs unphased with what phasing them needs, and a new
  **sv-phase** stage (`lrpipe sv phase`) phases them from the list into `sv/SAMPLE.sv.vcf` and
  `.snf`, exactly as Sniffles2 phases them from reads tagged as the list says (phase, PS and the
  phasing rescue of MOSAIC_VAF calls); the CNVs and the merge read those. A new **repeats** stage
  (`lrpipe tr --haplotags`) genotypes the tandem repeats after the list instead of the read
  pass, giving the read pass's records for the same haplotypes (on HG002 ONT chr20 84.4% of the
  called regions' alleles exactly right with haplotags, 69.2% from lengths alone). The CpG
  tables wait for the list and take it. ONT's bedMethyl and depth are unchanged, and a run of
  CNVs alone with phasing does not haplotag. Without phasing, the run is as before but for the
  warning's advice, which now names `phase` (and the repeats). A new phased VCF reruns the
  haplotagging and the stages reading the list. The manifest's
  `sample.haplotypes` says `lrpipe run's phasing of the SNVs (its haplotag list)`.
- **Haplotag lists in place of haplotagged BAMs**: `lrpipe methyl cpg|pileup --haplotags LIST`
  and `lrpipe sv --haplotags LIST` take the reads' haplotypes from a haplotag list (WhatsHap's
  `--output-haplotag-list`), as `haplotag` would tag them (a listed read's primary alignment, no
  other), overriding the input's tags; `lrpipe sv --phasing-evidence FILE` and `lrpipe sv phase`
  phase the SVs after calling. On lrpipe-sv's parity fixture phased from its own tags as a list,
  the result is Sniffles 2.8.1's VCF from the tags and lrpipe-sv's SNF, as it is on reads split
  into primary and supplementary alignments, with secondary ones, tagged otherwise by another
  pipeline (`lrpipe tr --haplotags` likewise). A read listed more than once on a contig with
  different haplotypes (a list written with `--tag-supplementary`) is left untagged with a
  warning, and each tool says how many of its reads the list names, warning when fewer than half
  (a list of other reads). The list takes about 24 bytes and the name per listed read in
  memory.
- **`lrpipe haplotag`**: `--list-only` (the haplotag list without the alignments),
  `--write-index` (index the output BAM) and `--untag-missing-contigs` (write reads on contigs the
  VCF lacks untagged, where WhatsHap fails or drops them), lrpipe-phase's own options.

- **`lrpipe run` from its user's side** (a command-line audit, LUX):
  - Ctrl-C ends a run as `interrupted` with exit code 130 and its running stages `cancelled`.
    The terminal's SIGINT reached the stages as well as lrpipe, and a stage that died of it
    first was recorded as failed (`snv: FAILED: killed by signal 2`, exit code 2, manifest
    `failed`), every time. Each stage now runs in a process group of its own, which lrpipe kills
    whole, and a stage that dies of SIGINT or SIGTERM (SLURM's `scancel` signals every process
    of the job) is stopped with the run, not failed. Two signals in quick succession
    (`timeout -s INT`, which signals lrpipe and its group) ended lrpipe at once with the
    manifest left `running`; a second signal now changes nothing, and a killed stage that does
    not end within 10 s is recorded as stopped and left, its process ID printed and recorded
    (LUX-01). SIGHUP (a terminal closed, an ssh session dropped, which no longer reached the
    stages) and SIGQUIT stop the run as SIGINT does, but SIGHUP is left alone under `nohup`. A
    stage's process group (bonefish's, merge-reads' too) is killed when the stage ends, so nothing it left running outlives
    it; no stage starts once the run is asked to stop; and a run that cannot write its own
    files kills its stages and records itself as failed where it can. lrpipe-snv's LongPhase
    ends with lrpipe-snv (`PR_SET_PDEATHSIG`); it outlived an lrpipe-snv that was killed.
  - A reference without a contig the reads have reads on, or with it at another length than
    the reads' header gives, is an error before any stage, saying which and, for `20` against
    `chr20`, that the two name contigs differently. Before, lrpipe-snv wrote an empty VCF and
    the run exited 0 (with every analysis, the read pass then failed). Contigs without reads are
    not checked, as lrpipe-methyl checks (LPH-004), so a reference of chr20 serves chr20 reads
    whose header names a whole genome. A BAM with up to 1% of its mapped reads (by its
    index's counts) on contigs the reference lacks, as HLA, decoy or EBV contigs against the
    no-alt analysis set, gets a warning naming them and their reads; more, or every read, is
    an error, as is any for CRAM (LUX-02).
  - A failed stage's last 10 lines that are not blank are printed after the error naming its
    log; the terminal said only `stage reads failed: exit code 1`. `[snv] device = "gpu"`, which
    passed `--dry-run` and failed in the snv stage, is an error at planning (LUX-03).
  - `--sample` is up to 200 letters, digits, `.`, `-` and `_`, not starting with `.` or `-`. `../evil`
    wrote `evil.vcf.gz` outside `--output`, `a/b` made directories, an empty name hidden files,
    and a tab or space broke the merged VCF's sample column (LUX-04).
  - `--analyses all` is accepted, as the help's default; it was refused (LUX-07).
  - The reads are checked before anything reads their index (each of several `--reads`, and for
    unaligned reads the align stage's BAM before the rest is planned), each problem with its own message
    and the `samtools` command that fixes it: a missing file, SAM, BAM named `.cram` or CRAM
    named `.bam`, reads sorted by name (`SO:unsorted` with an index is a warning), no index under the name the tools look for
    (`READS.bai`, `READS.csi`, `READS.crai`; `reads.bai` beside `reads.bam` is named but is not
    one), and an index that cannot be read or indexes another number of contigs. Before, a
    missing file or SAM was said to have no index, and a wrong index gave `invalid BAI header`
    (LUX-08).
  - A missing reference is said to be missing (it was its `.fai`), and a missing `.fai` or
    `.gzi` names `samtools faidx` (LUX-09).
  - A `--sample` other than the reads' read groups' `SM` is warned of (LUX-11).
  - lrpipe-snv's model is checked at planning, `--dry-run` included: a `--model` name not in
    `[snv] models` (whose models the error lists), or a directory without `pileup.pt` or
    `full_alignment.pt` (or `pileup` and `full_alignment`, as lrpipe-snv also takes them), is
    an error; it was a warning that the model was "not checked" (LUX-13).
  - `--memory` takes `32`, `32G` or `32GiB` (and `M`, `T`), and refuses less than 1 GiB; `8G`
    was "invalid float literal", and `0` was taken. `--threads 0` is an error (it ran with 1),
    and more threads than lrpipe may use is warned of (LUX-14).
  - A run into a directory holding outputs of an earlier run with other analyses (a full run,
    then `--analyses snv`) warns of those it did not write, which are left (LUX-15).
  - `lrpipe methyl pileup READS -` (or `stdout`) writes the bedMethyl to the standard output,
    as modkit does; it wrote a file named `-` (LUX-05).
  - Every subcommand takes `--version` and `-V`, the tools' own subcommands too; `lrpipe snv`,
    `cnv`, `phase`, `haplotag`, `run` and `setup` refused it, though `lrpipe phase help` listed
    it (LUX-06).
  - `lrpipe tr` takes `--karyotype XY` as its help writes it (and `--platform ONT`), takes
    `--reference` for `--genome`, says what `-t` is, and indexes a `.vcf.gz` output (`.tbi`),
    which it left unindexed (LUX-16).
  - `lrpipe help TOOL [COMMAND]` is the tool's (or its command's) help (it was "unrecognized subcommand"); a mistyped command
    is told the nearest (`lrpipe snvv`: `snv`); and the help of `lrpipe snv` and `lrpipe phase`,
    whose `run` and `phase` are implied, names their other commands (LUX-17).
  - `lrpipe snv run -m DIR` without the networks names `DIR/pileup.pt`, not `DIR/pileup`.
  - `lrpipe sex` on reads of some autosomes (a chr20 BAM) says it has reads on 1 of its 22
    autosomes (chr20) and that the sex is estimated from a whole genome's reads; it said there
    were no reads on the autosomes (LUX-10).
  - `lrpipe compare` ignores the merged VCF's `##snv_cmdline`, `##sv_command` and
    `##sv_fileDate` as the README says (two identical runs' headers differed), gives no change
    in per cent for a stage one run skipped, a run that ran no stage, or a value that shows as
    0 (`run 0 -> 7 +2404%`), and says "1 record" (LUX-12). Sizes under 1 kB are bytes, not
    `0.0 kB`.
  - `lrpipe setup` keeps a `site.toml` changed by hand since a setup wrote it, writing the
    rerun's as `site.toml.new` with a message; it replaced it, keeping the edits only as
    `site.toml.previous` until the next rerun. A setup from before lrpipe recorded the
    checksum counts as unchanged when its `site.toml` is the one its recorded manifest and
    platforms make (LUX-18).
  - lrpipe-cnv `--vcf-style standard` (which `lrpipe run` passes) ends a call that reaches a
    contig's end at its last base: END was the last window's start plus the window size,
    past the contig (`END=64445000` on chr20, of 64,444,167 bases), which strict readers
    refuse. Spectre's own style, and Spectre, write the last window's start, so this was
    lrpipe-cnv's, not a Spectre bug.
  - Wording: "1 thread" and "1 core"; every warning starts with `lrpipe:`; the model's reason no
    longer repeats where it was given (`snv.ont_model: given by snv.ont_model: ...`); reads of
    another platform let through with `--model-mismatch warn` are no longer said to lack `@RG PM`
    or a basecall model; an `--output` that is a file is said to be one, not `File exists (os
    error 17)`.
- **Several aligned read files per sample** (F5): `lrpipe run --reads` may be given more than
  once, BAM and CRAM alike. A first stage, `merge-reads` (the new `lrpipe merge-reads`, in Rust,
  without samtools), merges them into one sorted, indexed BAM, `OUT/reads/SAMPLE.bam`, which
  every other stage reads: it streams the files, merging their records by position as `samtools
  merge` orders them (forward strand first at one position, then by file), reads each file with a
  thread of its own and compresses with all the run's threads, indexing as it writes (two BAMs of
  HG003 ONT chr20, 2.5 GB: 7.5 s into memory against `samtools merge -@16`'s 12.5 s, 0.13 GiB).
  The merged BAM takes about the inputs' size on disk and is kept: resumed, the run does not
  merge again, and an input or its index replaced merges again. Before anything runs, the
  files' `@SQ` lines must be the same (names, lengths, order), a read group ID in two files
  must have the same fields, and their read groups must name one sample (`SM`), or the run stops
  naming the first difference. Every file's `@PG` chain is kept whole (a clashing ID renamed
  `ID-N`, its followers' `PP` changed with it, as `samtools merge` does), and the merge adds its
  own `@PG` line; reads found in two
  files (from a sample of each one's names) are warned of. Clair3's model is chosen from all
  the files' read groups, each file is sampled for `HP` tags (one haplotagged beside one not is
  warned of), and `--sex auto` adds up all their index counts. The manifest's `sample.reads` is
  now a list of the files given, for one file too (it was that file's record, an object). One
  `--reads` runs as before. Unplaced unmapped reads are left out of the merge, as no stage reads
  them. (Several unaligned files are aligned into one BAM by the align stage, below; aligned and
  unaligned files together are refused.)
- **lrpipe-sv reads a bgzipped reference** (`ref.fa.gz` with its `.fai` and `.gzi`, as references
  are distributed and as the repository's own fixture is), told apart by contents as htslib and
  pysam do: the SV calls, CRAM decoding, depth, copy numbers, repeats and bedMethyl, with the
  same output as from the plain FASTA. It refused one before ("compressed reference FASTA is not
  supported yet"), so `lrpipe run`'s read pass failed on such a reference. A FASTA compressed
  with gzip alone is the error lrpipe-reads gives (recompress it with `bgzip`). lrpipe-reads'
  new `FastaFile` reads either at the `.fai`'s offsets, from several threads at once, and
  `indexed_fasta` opens either as pasta's reader. lrpipe-tr's reference windows and lrpipe-cnv's
  GC fractions and `.mdr` (in their own tools and in lrpipe-sv's pass) tell a bgzipped reference
  by contents too: they went by its name, and read a bgzipped `ref.fa` as text.
- **Unaligned reads, aligned by bonefish** (F6): `lrpipe run --reads` also takes unaligned
  reads, unaligned BAMs (no `@SQ` lines, or only unmapped reads, as dorado and PacBio write
  them) or FASTQ or FASTQ.gz, and is given once for each file. A new align stage, first and
  alone, maps them with bonefish (`-x dna` or `-x hifi` by `--platform`, and `-y`, which keeps
  `MM`/`ML`/`MN`, `RG` and the `@RG` lines; FASTQ, which has none, gets one read group named
  after the sample) into `OUT/align/SAMPLE.bam`, which every other stage reads; the sex and the
  `HP` tags are read from it once it is written. The Clair3 model is chosen before, from the
  unaligned BAMs' header, and the rest is planned and its files checked before aligning, so
  what would fail after it fails in seconds; a failure or Ctrl-C between the align stage and
  the others finishes the manifest as `failed` (with its new `error`) or `interrupted`. bonefish
  is `[tools] bonefish`, else beside `lrpipe`, else on `PATH`, checked with `--version` and
  recorded in the manifest; `[align] index` names its index of the reference, which `lrpipe
  setup --bonefish PATH` builds and writes into `site.toml`, and whose contigs, like the aligned
  BAM's, must be the reference's. The stage takes bonefish's documented settings that fit
  `--memory`, for its threads and with the index build where there is no index (`--sort-mem 1G`
  below 23 GiB with 16 threads, and `-t 4 -K 10M` below 19), which change neither its output nor
  its fingerprint, so a resumed run never realigns, and keeps every stage's measurements; the
  reads, reference, index and bonefish binary are its inputs. Its scratch directory,
  `OUT/align/tmp`, is removed once it succeeds. Reads without `MM`/`ML` tags give no
  methylation: `methyl` asked for is refused, and the default analyses leave it out with a
  warning, rather than writing empty tables. Refused: aligned and unaligned `--reads` together
  (several aligned ones are merged, above), unaligned BAMs with FASTQ, a file given twice, unaligned BAMs of several
  samples or with one read group ID of two `@RG` lines, and an unaligned CRAM, which bonefish
  cannot read. Aligned BAM or CRAM runs as before, never through bonefish. The manifest's
  `sample.reads` (a list) is empty until the aligned BAM it names exists, and
  `sample.unaligned_reads` lists unaligned input.
- **The Clair3 model from the reads' header** (F3): `lrpipe run` takes lrpipe-snv's model
  from the site's new `[snv] models` directory (Clair3's models under their standard names), by
  the reads' `@RG` header: dorado's `basecall_model` (`dna_r10.4.1_e8.2_400bps_sup@v5.0.0`
  takes `r1041_e82_400bps_sup_v500`) or PacBio's instrument (`PM:REVIO` takes `hifi_revio`,
  `PM:SEQUELII` `hifi_sequel2`). A basecall model Clair3 has none for gets the nearest of the
  same chemistry, with a warning naming it and why; an unknown chemistry, or read groups
  basecalled differently, are errors until a model is given. A model given (the new `--model`,
  or `snv.ont_model`/`snv.hifi_model`, which keep working as before and take precedence) is
  checked against the header: one for another platform or chemistry (a HiFi model for ONT
  reads, R9.4.1 for R10.4.1) is an error unless `--model-mismatch warn`, one for another
  basecall model or instrument a warning, and one lrpipe does not know by name (a custom or
  fine-tuned model) a warning that it is not checked. Before, a wrong model cost accuracy with
  no word. The manifest's `snv_model` and `--dry-run` say which model and why. From
  `[snv] models` the nearest model the directory holds is used where it lacks the one for the
  reads (Clair3's own models have no dorado v4.2.0 or v4.3.0 model, which are ONT's Rerio
  models converted for Clair3 2), with a warning naming the one to add and where it comes from;
  the table has Rerio's models for dorado v4.0.0 to v4.3.0, R10.4.1 at 260 bps and Guppy 6, and
  tells 4 kHz reads (dorado before v4.2.0, Guppy) from 5 kHz ones. Runs `lrpipe run` accepted
  before and now refuses (`--model-mismatch warn` runs the first three): a configured
  `ont_model` or `hifi_model` for another platform or chemistry than the reads' header says (an
  R9.4.1 `ont` model for R10.4.1 reads, a HiFi model for ONT reads); `--platform hifi` for
  reads the header says are Oxford Nanopore, or the reverse; PacBio reads whose `READTYPE` is
  `SUBREAD`, `SCRAP`, `ZMW` or `HQREGION`; a move-table (`_with_mv`) model, which failed in the
  snv stage after the other stages had started; and reads whose header cannot be read as BAM or
  CRAM, which failed in the stages.
- **`lrpipe setup`**: `lrpipe setup --genome GRCh38|CHM13|GRCh37 --platform ont,hifi --dir DIR`
  fetches the reference, Clair3's models for the platforms (into `DIR/models/NAME/`, under
  Clair3's names), Sniffles' TRF BED, the repeat catalogue and Spectre's blacklist from
  manifests checked into the repository (`lrpipe/resources/`: URL, size, SHA-256 or the
  publisher's MD5, and licence of each), builds the `.fai` (and a `.gzi` for a bgzipped
  reference), lrpipe-cnv's `.mdr` and GRCh38's `catalog.adotto_strchive.hg38.bed`, and writes
  `DIR/site.toml` with absolute paths and each file's licence. Downloads resume (HTTP Range),
  are checked before use and fetched again when wrong; a rerun skips what is present and
  verified. `--offline` prints the URLs, checksums and `curl`/`wget` commands, `--verify` checks
  every file again, `--accept-licences` accepts the terms of ONT's models and the TRF BEDs
  (all three pbsv's, under PacBio's licence), and files with no published checksum (for now
  all of Clair3's models) are fetched like any other but checked by their size alone, and
  reported as checked by size only (never as verified) until a maintainer records their
  SHA-256 with `--pin`; a directory set up before keeps what matches. Resumes send `If-Range`, a transfer
  is resumed until five tries in a row bring nothing (TLS, name and 4xx failures are not
  retried), progress is printed, one setup holds a directory at a time, and a rerun for
  another platform keeps the earlier one's models. GRCh37 has no
  repeat catalogue and CHM13 no blacklist (none is published for them).
- **No more lost files** (pre-release bug report, LRUN-01, LRUN-09, LSNV-001, LRUN-04):
  - `lrpipe merge` writes, compresses and indexes under hidden names of its own beside the
    output, and renames them into place: `--snv calls.vcf --output calls.vcf.gz` no longer
    overwrites and removes `calls.vcf` (leaving a truncated output, with exit 0), and `lrpipe
    run` no longer destroys a file `OUTPUT/SAMPLE.vcf`. An `--output` that is one of the
    inputs is refused, a failed merge keeps the earlier output, and one stopped between its
    renames leaves no index rather than the earlier one; `lrpipe run` removes a killed merge's
    hidden files.
  - lrpipe-snv checks that the BAM or CRAM and its index, the FASTA and its `.fai`, a
    `--bed_fn`/`--vcf_fn` and the model's checkpoints can be read before it removes an earlier
    run's results, so a mistyped path no longer costs them.
  - lrpipe-snv removes an earlier run's `tmp/` and final files only where it marked them as its
    own (a `.lrpipe-snv` file in `tmp/`, and one beside the results); elsewhere it stops with
    them untouched, so `--output .` or `~` no longer deletes a user's `tmp/` or files named like
    Clair3's. A standalone lrpipe-snv output directory from 0.1.5 or earlier is refused once:
    remove its `tmp/` and results, or use another directory. `lrpipe run` marks its own `snv/`,
    so its directories from 0.1.5 resume as before. It marks `snv/` only in a directory that
    held an lrpipe run (a `.lrpipe/sample` or manifest): in any other, an `snv/` with files in
    it (an earlier Clair3 run's, say) is refused untouched, where it used to be marked and
    then deleted by the snv stage. An empty `tmp/` left by a run killed before marking it no
    longer turns later runs away.
  - An empty or cut-off model checkpoint is an error ("not a zip archive (is the model file
    complete?)") instead of a crash of the snv stage.
  - `lrpipe run` refuses an output directory that holds another sample's run, before writing
    anything there: two samples in one directory replaced each other's `snv/`, logs and
    manifest. The sample is recorded in `.lrpipe/sample`; the error names the file to remove
    (`.lrpipe/sample`, or the manifest of a directory from 0.1.5) if it holds no such run.
- **No silent wrong output** (pre-release bug hunt, batch 3):
  - A bgzipped reference (`.fa.gz` with `.fai` and `.gzi`, as references are distributed) is
    read through its `.gzi` by every tool that reads the reference through lrpipe-reads
    (lrpipe-methyl, lrpipe-depth, lrpipe-tr's CRAM decoding). It was read as plain text, so
    `lrpipe methyl cpg` and `lrpipe run --analyses snv,methyl` wrote empty CpG tables and
    exited 0. A FASTA compressed with gzip but not bgzip is an error, and `lrpipe run` checks
    a bgzipped reference's `.gzi` before any stage runs (LRUN-06).
  - lrpipe-snv `--approximate` on a CPU without AVX-VNNI, AVX-512 VNNI or aarch64's dot
    product instructions changes nothing, with a warning. It sent 50% of the ONT variant calls
    to full alignment instead of 70% there, for no 8-bit speed-up; `lrpipe run` always passes
    it (LSNV-004).
  - lrpipe-methyl fails when the reference lacks a contig that has reads, or has it of another
    length (a reference naming contigs `1` against reads on `chr1`, or another assembly's),
    where it wrote no rows for them and exited 0 (LPH-004).
  - lrpipe-methyl `pileup` rejects a `--filter-threshold` other than modkit's forms (`0.8`,
    `C:0.8`) and a `--filter-percentile` outside [0, 1); `5mC:0.99` and `1.5` filtered
    nothing. A base's own threshold now wins over one for every base given after it, as in
    modkit (LPH-005).
- **`lrpipe run`: scheduling, memory and resume** (pre-release bug report, LRUN-02, 03, 05, 07,
  08 and 12):
  - A resumed run never expects a stage to need less than the plan's estimate: a peak measured
    when it succeeded counts, scaled by 1.15, where that is higher, and one of a run that was
    stopped or failed only where it is higher (how far it got). A stage stopped or failed
    early had recorded a tiny peak, so lrpipe-snv and the read pass started together beyond
    `--memory`.
  - The default `--memory` is what is left under lrpipe's cgroup limit (cgroup v1 or v2, up
    the hierarchy: SLURM `--mem`, `docker run -m`, Kubernetes) when that is less than
    `MemAvailable`.
  - Progress lines that cannot be written (stderr a pipe whose reader exited, or a closed
    terminal) are dropped; lrpipe aborted there, killing every stage, with no manifest.
  - Resume ignores `--threads` and how paths are spelt: stamps leave out thread counts and
    take paths with their directories resolved, and the configuration's relative paths are
    taken from its resolved directory. Stamps written by earlier versions do not match, so the
    first run after updating reruns its stages once.
  - Without `[cnv] metadata`, lrpipe-cnv's `.mdr` of N runs is named for the reference it is
    made from (`cnv/reference-KEY.mdr`), so a run with another reference in the same output
    directory no longer reuses the earlier reference's gaps. lrpipe-cnv writes it whole or not
    at all (to a temporary file, synced, then renamed).
  - A configured `[cnv] metadata` that does not exist yet is made from the reference, as the
    README said (its directory too), where `lrpipe run` refused to start; it is made after the
    output directory is locked and the reads are checked. Only the files the chosen `--analyses`
    read need to exist (a run without `tr` no longer needs the repeat catalogue, an ONT run the
    HiFi model).
- **lrpipe-sv: documented options** (pre-release bug report, LSV-002 to LSV-004):
  - lrpipe-sv: `--contig` with `--depth` and `--cnv` calls copy number on the `--contig` contigs
    (those longer than 1 Mbp, each once and in the reference's order, whatever the order of the
    options; none is an error before the pass). It failed after the whole pass
    with "no coverage for --only-chr contig" whenever the reference had another contig over
    1 Mbp.
  - lrpipe-sv: multi-sample calling reads any number of `.snf` files. `--combine-close-handles`
    opens each file again for every block read, as in Sniffles2, and is turned on when the inputs
    would not all stay open under the open-file limit (`ulimit -n`); combining 1,022 or more files
    under the usual limit of 1024 failed with "Too many open files", flag or not.
  - lrpipe-sv: `--genotype-vcf` with a sites-only VCF was fixed in 0.1.5 (LRP-008), for every
    SVTYPE ending the INFO, not only BND; a parity test now covers records with SVTYPE last.
- **lrpipe-snv: a correct gVCF** (`--gvcf`, and `lrpipe run --gvcf`; README, "Clair3 bugs fixed
  on purpose"). It had Clair3's `mergeCalls` bugs: a call inside an earlier deletion was
  replaced by a second copy of a block (1 of 467 calls lost on the ONT fixture), blocks cut by
  a call kept their original start and the block holding a contig's last call was written
  whole, so blocks overlapped calls and each other, and a contig without calls lost its first
  block. A chunk without reads ended its `N` block on the next chunk's first position, and the
  last past the contig's end. Now every call is written, blocks neither overlap nor go missing,
  and `ci/parity.sh` checks the fixture's gVCF against Clair3's but for these differences.
- **`lrpipe run`: cohort outputs.** `--gvcf` passes `--gvcf` to lrpipe-snv and lists
  `snv/merge_output.gvcf.gz` and its index among the snv stage's outputs and the manifest's
  (`gvcf`), for joint genotyping with GLnexus; turning it on or off reruns only that stage and
  what reads its calls. The read pass writes `sv/SAMPLE.snf` whenever it calls SVs (manifest
  `snf`), for multi-sample calling with `lrpipe sv -i a.snf b.snf ...`: on HG003 ONT chr20 it
  added no measurable time, 20 MB of memory and a 0.9 MB file, and the SV VCF is the same.
  As the read pass's command gains `--snf`, resuming a run made before this reruns the read
  pass, and with it cnv and merge.
  `lrpipe compare` compares `.snf` files block by block, decompressed, ignoring the command line.
- **`lrpipe run`: no empty haplotype tables from reads that are not haplotagged** (F1a; with
  phasing, the run now haplotags the reads itself, above, so this holds for runs without
  `phase`). Without phasing, the SVs' phases and the CpG tables per haplotype come from the
  reads' integer `HP` tags; on reads
  without them (minimap2's, pbmm2's, dorado's) a run wrote empty `hap1`/`hap2` tables and exited
  0, while the README promised methylation per haplotype. Before planning, `lrpipe run` now
  looks for an integer `HP` tag in 200 mapped reads from the middle of each of the first 8
  contigs with reads (lrpipe-reads' `IndexedReads::sample_tag`; 0.08 s on a 39 GB HiFi BAM,
  decoding only flags, positions and the tag of a CRAM). Without one it plans the combined CpG
  table alone (`lrpipe methyl cpg --combined-only`, a new option), removes an earlier run's
  `hap1`/`hap2` tables, warns at the start and in the summary (naming `phase` among
  `--analyses` as the remedy, and `HP` tags of another type, which the tools read as none), and
  records `sample.haplotypes: "none (input not haplotagged)"` and the warning (`warnings`) in
  the manifest. Reads haplotagged by another pipeline keep their tables. Resuming an earlier run
  on HiFi reads without `HP` tags reruns the CpG stage (its command gains `--combined-only`).
- **The README says what lrpipe does not cover** ("Not covered"): somatic and trio-aware
  calling, annotation, pharmacogenomics, HLA typing and paralogous genes, and AMD GPUs, each with
  the tools to add.

## 0.1.5

- **Fixes from a second external review** (findings LRP-005, 008 and 015 to 020; each tool's
  README lists the ones that change what the original tool writes, and the parity tests expect
  exactly that):
  - lrpipe-sv: `--genotype-vcf` gives a breakend only the support of the same breakend (mate
    contig, mate within `--cluster-merge-bnd`, orientation on both sides), where Sniffles2 took
    any breakend nearby with a mate on that contig; combining `.snf` files keeps breakends of
    opposite orientation apart. The genotyped VCF declares PS and names one sample, the
    genotyped one, whatever the input's columns (sites-only and multi-sample inputs gave
    headers that did not describe the records), and a sites-only record ending in
    `SVTYPE=BND` is genotyped.
  - lrpipe-snv: a contig without a count in a BAI that counts others is not taken for one
    without reads (unless the index has no bins for it); with `--enable_phasing`, a run that
    ends early writes `phased_merge_output.vcf.gz` too, and phases the pileup calls when none
    goes to full alignment.
  - lrpipe-phase: a sample's haploid call leaves another sample's diploid call at the same site
    in its phasing (only records where no sample has a diploid call stay out).
  - lrpipe-cnv: `--sex` takes chrX's PARs from the shared table, which has CHM13's.
  - `lrpipe run` takes its output directory's lock before writing the PAR BED, so a run turned
    away from a busy directory no longer replaces the running one's.
- The whole-genome validation scripts, with their inputs, versions, model hashes and raw results,
  are in [tools/validation](tools/validation). BENCHMARKS.md's F1 bound against Clair3 is 0.0003
  (it said 0.0001).
- Two tests no longer fail in containers: CPU confinement under a CPU quota, and niceness where
  `ps` fails.
- GitHub runs `tools/check.sh` on every push to main and every pull request; the GPU, parity and
  WhatsHap checks skip there, so a release is still checked locally with `--release`.

## 0.1.4

- **Fixes from an external review** (findings LRP-001 to LRP-014). Where a fix changes what an
  original tool would write, the tool's README lists it and its parity tests expect exactly that
  change:
  - Resume: a model replaced in the same directory reruns calling (inputs up to 64 MiB are
    compared by contents), a deleted index reruns the stage that writes it, and lrpipe-snv
    starts by removing an earlier run's `tmp/` and final files, so none can mix into its results.
  - Sex chromosomes: `lrpipe run` gives lrpipe-snv the sample's sex and the reference's
    pseudo-autosomal regions (GRCh38, GRCh37 or CHM13), so XX samples get no chrY calls and XY
    samples haploid chrX and chrY calls outside the PARs; a PAR BED row now covers the VCF
    positions it names (Clair3 is one base off); lrpipe-phase leaves haploid calls unphased
    instead of stopping. lrpipe-tr takes the PARs from the reference too (they were GRCh38's).
  - lrpipe-sv: breakend mates keep their names with a reference (Sniffles2 turned `chrY` into
    `chrN` and `chrUn_KI270746v1` into `chrUn_NI270746v1`: 23-33 records per genome), combining
    `.snf` files keeps breakends, `--genotype-vcf` declares PS in FORMAT, and a CSI index is
    read.
  - lrpipe-snv: a BAI without its optional read counts no longer means no reads, and
    `--pileup_only` on the CPU writes post-processed records.
  - lrpipe-methyl: `--modified-bases 5mC` keeps opposite-strand (`G-m`) calls.
  - lrpipe-cnv: `1/0` and `1|0` count as heterozygous; CNV depth and TR catalogs may be plain
    gzip.
- **Apple silicon and Linux aarch64**: `lrpipe snv --approximate` runs its 8-bit networks with
  `sdot` on CPUs with the dot product instructions (Apple M1 and later, Graviton 2 and later),
  where until now it fell back to `float32` and was no faster than exact; on an M4 Pro, HG002 ONT
  chr20 in 56 s instead of 126 s, with the same GIAB F1. The Smith-Waterman realignment has a
  NEON version. A whole HG002 ONT genome runs through `lrpipe run` in 32 min 35 s there; see the
  README's [Apple silicon and ARM](README.md#apple-silicon-and-arm).
- **lrpipe-snv uses AVX-512 on the CPU.** `--approximate`'s 8-bit products run on CPUs with AVX-512 VNNI
  (Xeons from Cascade Lake, AMD Zen 4 and later), with each block's 16 channels in one 512-bit
  register; they used to need AVX-VNNI, which most server CPUs lack, so there `--approximate`
  ran in `float32`. The sums are the same integers either way. The `float32` products use
  `gemm`'s AVX-512 kernels on CPUs with AVX-512, which sum in another order than its AVX2 ones:
  probabilities can differ in the sixth decimal between such CPUs. Not yet measured on an
  AVX-512 machine.
- **`lrpipe methyl cpg` and `lrpipe methyl pileup` are faster through `lrpipe`**: the mmap threshold
  `lrpipe` fixed for `lrpipe sv`'s read pass applied to every tool, so they mapped and unmapped a
  buffer per read. Fixed for `lrpipe sv` only, the whole HG002 Revio genome's CpG tables take
  217 s instead of 389 s, and the HG002 ONT CRAM's bedMethyl 251 s instead of 303 s, the same
  files. Found by [BENCHMARKS.md](BENCHMARKS.md), every tool rerun on whole genomes.
- **[BENCHMARKS.md](BENCHMARKS.md)**: every tool rerun on whole genomes against the tool it
  replaces, and phasing scored against GIAB's phased HG002 Q100 benchmark: `lrpipe run`'s
  phased calls have 0.31% switch errors on HiFi (block NG50 606 kb) and 0.40% on ONT (1.99 Mb).
- On macOS: peak memory in the manifest was 1024 times too large (`ru_maxrss` is in bytes
  there), `--memory` is no longer needed (the memory available comes from
  `kern.memorystatus_level`), lrpipe-snv gives freed heap back with
  `malloc_zone_pressure_relief`, and the scripts run with macOS's sed, gzip and bash 3.2.
- lrpipe-snv reads its pileup VCF one contig at a time from disk in phasing and full-alignment
  selection, instead of decompressing it whole into memory: on noisy low-depth reads it holds
  hundreds of millions of records (HG002 ONT at 10x: 18.6 GB), more than many machines have.
- `tools/check.sh` reports passed, failed and skipped checks; `--release` fails on a skip.
- `lrpipe-tr-bench` reports the call rate and accuracy over every region beside accuracy over
  the scored ones, and `--common` scores two call sets on the regions both score.

## 0.1.3

- **`lrpipe run` phases the SNVs and indels and writes one VCF**: `SAMPLE.vcf.gz` (indexed)
  holds the SNVs and indels, phased, the SVs, CNVs and tandem repeats, each record tagged with
  its analysis (`INFO/SOURCE`), from the new `lrpipe merge` (7 s on HG002 HiFi). lrpipe-snv runs
  with Clair3's `--enable_phasing`, so `snv/phased_merge_output.vcf.gz` is written too: on chr20
  of HG002 HiFi its records are WhatsHap 2.8's exactly. Each tool's own files stay as they were.
  Phasing the final calls makes a run 1.5-3 min longer: HG002 HiFi 6 min 32 s (4 min 57 s
  without), HG003 ONT 11 min 5 s (9 min 1 s), the HG002 ONT CRAM 17 min 39 s (14 min 53 s). On
  the GPU, lrpipe-snv merges and phases each contig as soon as its full-alignment calls are
  scored, beside the other contigs' calling, instead of all of them at the end (the same files,
  20-30 s sooner); lrpipe run expects it to peak at 14.5 GiB on ONT for that (13.4-14.3).
- **`lrpipe run --analyses`** runs only some analyses (`snv`, `phase`, `sv`, `cnv`, `tr`,
  `depth`, `methyl`), with what they need: `--analyses snv,sv` is lrpipe-snv without phasing and
  a read pass that calls SVs only. `--merged-vcf off` leaves out `SAMPLE.vcf.gz`.
- `lrpipe sv` counts depth, repeats or methylation without calling SVs when given no SV output.
- libdeflate (through pasta) is on for every tool in the workspace's own dependency on pasta, so
  a tool built alone (`cargo build -p lrpipe-depth`, ...) gets it too; until now only builds that
  included lrpipe-snv or lrpipe-phase did, as the release bundle and `cargo build` of the whole
  workspace do. `lrpipe-depth` built alone: the whole HG002 HiFi genome in 28.5 s instead of
  38 s, the same output. lrpipe-phase's `libdeflate` feature is gone.
- `lrpipe run`'s memory estimates follow the whole-genome peaks: lrpipe-snv 13 GiB on HiFi (14
  before; it peaks at 8.9-11.7) and 14.5 on ONT (16; 13.4-14.3 with phasing), the ONT read pass
  7 GiB with bedMethyl (6.5; it peaked at 6.3-6.9). On ONT the two run together with 21.5 GiB
  free.

## 0.1.2

- **`lrpipe phase` and `lrpipe haplotag`**: braid, the WhatsHap 2.8 port lrpipe-snv phases with,
  joined lrpipe with its history as lrpipe-phase (`git log -- lrpipe-phase/`). `lrpipe phase`
  takes WhatsHap's commands (`phase` implied, `haplotag`, `genotype`, `stats`, `compare`,
  `split`, ...), `lrpipe haplotag` is its haplotag, and `lrpipe-phase` is its own binary in the
  bundle, with WhatsHap's command line. Its output is still WhatsHap's.
- **Phasing about twice as fast** (still WhatsHap's output): `lrpipe snv` on the GPU takes the
  whole HG002 HiFi genome in 3 min 57 s, and `lrpipe run` 4 min 57 s. See
  `lrpipe-phase/README.md`.
- **`lrpipe snv` on the GPU about 18% faster** (the whole HG002 HiFi genome in 4 min 25 s against
  5 min 24 s, same calls): TF32 tensor cores from compute capability 8.0, a faster LSTM step,
  sleeping GPU waits by default (the pipeline no longer sets them only for ONT CRAMs), and
  indexed, parallel-compressed merges. See `lrpipe-snv/CHANGELOG.md`.

## 0.1.1

- **`lrpipe snv --approximate` on the CPU (what the pipeline runs there) 2.3 times as fast as
  exact calling**: the whole HG003 ONT genome in 31 min against 72 min exact and Clair3's 83 min,
  with 99.78% of the PASS calls Clair3's, and its pileup calls exactly an exact run's. The 8-bit
  full-alignment convolutions are 1.6 times as fast with the same results, and the pileup
  network runs in 8 bits too, with close calls scored again in `float32`, where a check on the
  run's own candidates allows it (ONT's model passes it; HiFi's keeps `float32`). See
  `lrpipe-snv/CHANGELOG.md`.
- One C library is left in the build, libdeflate (BGZF, through pasta). lrpipe-snv writes and
  reads its intermediate lz4 frames with `lz4_flex` instead of liblz4 (the `lz4` tool still
  reads them), and braid comes without mimalloc, which only braid's own binary used.
- lrpipe-cnv's `norm.sf` test allows 2 units in the last place: it goes through the C library's
  `exp`, which glibc 2.28 (the release build's) and macOS round differently.

## 0.1.0

First release: long-read calling in Rust, as one binary.

- `lrpipe run`: the pipeline on one sample, with `lrpipe snv` for SNVs and indels (with
  `--approximate`, which speeds up CPU calling): stages scheduled by their dependencies and memory,
  resumed after an interruption, and a manifest of every input, tool, command, time and peak
  memory. `lrpipe sex` reads a sample's
  sex chromosomes from its index; `lrpipe compare` compares two runs record by record.
- `lrpipe snv` (lrpipe-snv, merged from its own repository with its history): SNVs and indels
  from Clair3 2.0.3's models, its networks in Rust on the CPU or an NVIDIA GPU (the driver and
  NVRTC loaded at run time); PASS calls agree with Clair3's on 99.8-100%. See
  `lrpipe-snv/CHANGELOG.md`.
- `lrpipe sv` (lrpipe-sv): structural variants with VCF records identical to Sniffles
  2.8.1's, and from the same pass over the reads, read depth (mosdepth's files), tandem repeat
  genotypes (lrpipe-tr), copy number (lrpipe-cnv) and CpG bedMethyl (`--bedmethyl`, identical to
  `lrpipe methyl pileup`'s). The pipeline counts ONT methylation there, so the reads are decoded
  once for it and the SVs.
- `lrpipe cnv` (lrpipe-cnv): Spectre 0.2.1's `CNVCaller`, with identical output.
- `lrpipe tr` (lrpipe-tr): tandem repeat genotypes at every locus of a repeat catalog.
- `lrpipe depth` (lrpipe-depth): binned depth with output identical to mosdepth 0.3.14.
- `lrpipe methyl` (lrpipe-methyl): bedMethyl identical to modkit 0.6.4 `pileup` (`pileup`), and
  pb-CpG-tools' CpG tables (`cpg`): count mode byte-identical, and a model mode fitted against
  EM-seq that agrees with it better than pb-CpG-tools' model does.
- Each tool also builds as its own binary: `lrpipe-snv`, `lrpipe-sv`, `lrpipe-cnv`, `lrpipe-tr`,
  `lrpipe-depth` and `lrpipe-methyl`.
- **Release files**: a Linux x86_64 bundle (glibc 2.28 or later) with `lrpipe`, every tool's own
  binary and their licenses. GPU calling needs the NVIDIA driver and NVRTC at run time; Clair3's
  models are not included.

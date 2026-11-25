# lrpipe-methyl

Base modification pileups from aligned long reads, ONT and PacBio alike, in pure Rust: 5mC,
5hmC, 6mA or any other code in the reads' `MM`/`ML` tags, per reference position or per CpG
site, combined and per haplotype.

Two commands reproduce the two tools long-read pipelines use:

* `lrpipe-methyl pileup` writes modkit's bedMethyl (18 columns) with modkit `pileup`'s options.
* `lrpipe-methyl cpg` writes pb-CpG-tools' bed files (combined, hap1, hap2) with
  `aligned_bam_to_cpg_scores`' options.

Written from the SAMtags specification and both tools' documentation, and checked against
their output as black boxes; no code of either was read or used (modkit is licensed for
research use only, pb-CpG-tools for PacBio data only).

## Parity

| Comparison (HG002) | Rows | Result |
|---|---|---|
| modkit 0.6.4 `pileup`, ONT chr20:30-31 Mb, every mode below at fixed thresholds | 9k to 1.3M per mode | byte-identical |
| modkit, ONT chr20:29-32 Mb, every position, 5mC 5hmC 6mA | 4,733,489 | byte-identical |
| modkit, ONT CRAM chr20:30-40 Mb, every position | 15,509,870 | byte-identical |
| modkit, ONT CRAM whole chr20, `--modified-bases 5mC 5hmC --cpg --combine-strands` | 1,535,800 | byte-identical |
| modkit, ONT CRAM whole genome, the same, thresholds estimated by each tool | 58,345,760 | byte-identical |
| modkit's threshold estimate (all reads of the region) | | identical (A 0.8027344, C 0.7441406, from the same 17.8M and 12.1M calls) |
| pb-CpG-tools 3.0.0 count mode, reference sites, PacBio Revio chr20, chrX and chr11/chr15 imprinted regions | 2,111,861 / 860,416 / 820,591 (combined / hap1 / hap2) | byte-identical |
| pb-CpG-tools count mode, denovo sites, same data | 2,125,920 / 874,627 / 834,684 | identical values; 99 site choices (0.005%) differ |
| pb-CpG-tools count mode, reference sites, PacBio Revio whole genome (45 GB BAM, about 30x) | 28,934,010 / 25,314,245 / 25,118,292 | byte-identical |

modkit modes checked: default (every position, both strands, 5mC 5hmC 6mA), `--combine-mods`,
`--modified-bases 5mC 5hmC 6mA`, `--modified-bases 5mC 5hmC --cpg`, the same with
`--combine-strands`, `--modified-bases C --cpg --combine-strands --combine-mods`, and `--phased`
(combined, hp1 and hp2).

pb-CpG-tools' model mode is a neural network trained by PacBio and cannot be reproduced.
**Experimental:** `lrpipe-methyl cpg --pileup-mode model` instead takes each site's posterior mean level under a fixed
model: how methylated and unmethylated reads' values are distributed (in sixteen bins, plus no
call, which leans unmethylated), and how levels are distributed across sites. The three were
fitted once on HG002 Revio reads at chr20's CpGs, each site's level taken from EM-seq (EpiQC,
GSE186383, two replicates merged), and are judged on the other chromosomes only, CpGs with at
least 10 reads in both (whole genome, the same Revio BAM):

| PacBio estimate | vs EM-seq, 26,351,605 CpGs: r, MAE, bias | vs ONT 5mC (modkit), 27,451,916 CpGs: r, MAE, bias |
|---|---|---|
| count (pb-CpG-tools and lrpipe-methyl) | 0.882, 17.9, -12.2 | 0.881, 17.7, -10.4 |
| pb-CpG-tools model | 0.917, 11.5, -7.0 | 0.921, 11.1, -5.1 |
| lrpipe-methyl model | 0.923, 9.3, -0.7 | 0.925, 9.3, +1.2 |

The model is fixed rather than learned from each input. Learned per input (as lrpipe-methyl first
did), the two read distributions are not identifiable, since an "unmethylated" component heavy
in no calls fits as well as the real one, and the levels moved with the sites sampled: a
calibration that fitted one input read 11 points high on another. All of the above is one Revio
run (movie m84011_220902_175841_s1), so it is untested on other runs, and like PacBio's network
the model describes the chemistry it was fitted on.

Pooling neighbouring CpGs' evidence did not help (held-out r up 0.003, MAE worse at every
distance tried), nor did smoothing along reads.

### A modkit bug fixed on purpose

With `--modified-bases`, a group of the `MM` tag is kept by its modified base on the strand it
names, so `G-m` (5mC on the opposite strand, as in duplex reads) counts for `5mC`. modkit
compares the group's own base (`G`) and drops it, halving the coverage of such reads (review
finding LRP-006). Reads with only same-strand groups, as most simplex ONT and PacBio reads
are, give modkit's output.

## Speed

16 threads (Core i5-12600K). Wall times varied by a second or two between runs; CPU time is
steadier.

pb-CpG-tools on the PacBio extract (2.2 GB BAM: chr20, chrX, 5 Mb of chr11/chr15), all three
tables:

| Mode | pb-CpG-tools | lrpipe-methyl |
|---|---|---|
| count, reference sites | 105 s | 12 s |
| count, denovo sites | 111 s | 25 s |
| model, reference sites | 54 s | 15 s |

The whole Revio BAM (HG002, 45 GB, about 30x), all three tables (another job used about one
core throughout):

| Mode | pb-CpG-tools | lrpipe-methyl |
|---|---|---|
| count, reference sites | 1,868 s, 26,051 CPU-s, 13.0 GB | 237 s, 2,480 CPU-s, 1.3 GB |
| model, reference sites | 1,114 s, 15,625 CPU-s, 12.9 GB | 254 s, 2,748 CPU-s, 1.8 GB |

Rerun on 2026-11-03 through `lrpipe methyl cpg` (count mode, nothing else running): 217 s,
2,317 CPU-s, 1.3 GB, against pb-CpG-tools' 1,932 s, 26,746 CPU-s and 12.5 GB, with the same
records ([BENCHMARKS.md](../../BENCHMARKS.md)).

modkit on the ONT CRAM (HG002, about 30x):

| Job | modkit | lrpipe-methyl |
|---|---|---|
| whole chr20, 5mC and 5hmC at CpGs, strands combined | 14-17 s, 180 CPU-s, 5.8 GB | 12-14 s, 127 CPU-s, 2.8 GB |
| chr20:30-40 Mb, every position | 31.7 s, 178 CPU-s, 6.8 GB | 8.9 s, 42 CPU-s, 3.2 GB |
| whole genome, CpGs, thresholds estimated | 816 s, 8,251 CPU-s, 14.0 GB | 404 s, 4,986 CPU-s, 2.4 GB |

Decoding CRAM takes most of lrpipe-methyl's time on CpG jobs. A CRAM slice spans 100 kb
or more and is decoded again by every window it reaches, so on CRAM the default window is
1 Mb (250 kb otherwise); `--window` trades memory for that repeated decoding. Reads are
counted as they are decoded rather than gathered first (only denovo CpG sites, which depend
on every read, gather a window's reads), so a window's memory does not grow with its depth:
the unplaced contigs and centromeres where tens of thousands of reads (over a billion bases)
pile into one window once took the whole-genome peak to 7 GB on ONT and 10.7 GB on Revio,
now 2.4 and 1.3 GB. The windows in progress together may still decode at most 1.5 median
windows' data per thread, by the index's byte counts, and only a contig in use is kept of
the reference.

## What the tools do, as lrpipe-methyl reproduces it

modkit `pileup` (without `--modified-bases`):

* A call's canonical probability is one less its modifications' (`(ML + 0.5) / 256` each); a
  call passes when its most likely state's probability reaches the base's threshold.
* An unlisted base of an implicit (`.`) group is canonical, except one base that is not
  counted at all: the first aligned base (insertions skipped) after the read's last listed
  call in SEQ order, whichever base or group that call belongs to.
* `Ndiff` counts reads with another base that have no valid call of that other base there.
* The percentage modified is the share in single precision, times 100, printed to two
  places (94 of 320 prints 29.37).
* A group whose code is its base (`C+C`, a modification of unknown kind) says nothing about
  its bases: they are counted nowhere, but its listed bases still mark where the read's calls
  end, and keep the read out of other bases' `Ndiff` as a valid call would.
* Primary alignments only.
* `--region chr:start-end` takes a 0-based start (one base later than samtools reads it) and
  the same end.
* `--include-bed BED` (modkit's option name, in place of `--region`) writes only the rows of
  positions in the BED's regions (sorted and merged, in the header's order), with a CpG's row
  kept when its C is in one, both strands' calls combined, and estimates the pass thresholds
  from the calls of the regions' reads. modkit's own differs at the edges: its docs say it keeps
  the "positions overlapping intervals" and estimates the thresholds from them, and modkit 0.6.4
  counts each strand's position on its own (a CpG split by an edge gets the inside strand's
  calls alone) and samples only the calls at the regions' positions, so in small regions its
  thresholds, and so some counts, differ. `cpg --include-bed BED` (lrpipe-methyl's own option)
  scores only the CpG sites in the regions; `lrpipe run --regions` passes both.

With `--modified-bases`, rows need the reference base, the canonical probability is the middle
of the `ML` bin the modifications leave, every base is counted, and every read with another
base counts in `Ndiff`.

The threshold estimate is the 10th percentile of every call's most likely state on the `ML`
scale (unlisted implicit bases counting as certain), over whole reads (soft clips and
insertions too). modkit samples reads at random; lrpipe-methyl samples evenly and deterministically
(the first reads of 256 windows of 250 kb spread over the regions, the same for a BAM and a CRAM
of the same reads and whatever `--window`). On the HG002 genome both chose the same thresholds,
but with `--region` they can differ by several `ML` bins: HG002 ONT chr20 gave 0.7598 against
modkit's 0.7949, and chr20:29-32 Mb 0.7598 against 0.8223, which changes some counts. Give
`--filter-threshold` (modkit's estimate, say) for modkit's counts exactly: with the same
thresholds every mode compared was byte-identical, down to thresholds below one half, where a
call's tie between the canonical base and a modification goes to the modification as modkit's
does with `--modified-bases` (modkit's default mode gave no such tie to compare in 3 Mb of HG002
ONT). A tie between two modifications (5mC and 5hmC, seen only with `--no-filtering` or a
threshold below one half) can still go the other way: 6 of 62,514 rows in 100 kb of HG002 ONT
with `--modified-bases`, and 2,552 of 5,045,186 rows in chr20:29-32 Mb in the default mode.

pb-CpG-tools count mode:

* A read counts at a CpG when a base aligns to the C (either strand) and there are at least 20
  aligned bases, insertions and deletions aside, between that base and both ends of its
  alignment. Supplementary alignments count; MAPQ at least 1.
* Its value is the 5mC call it places at the site: a forward read's call where its C aligns, a
  reverse read's where the base after its C (as sequenced) aligns. A read without a call there
  (another base, a deletion at the G, an unlisted C) counts as unmodified with value 0.
* Modified means a value over 0.5. The share is rounded half up; the mean values are summed in
  single precision and rounded half to even.
* Denovo sites: where, in any of all reads, haplotype 1 and haplotype 2, the counted reads'
  most common base at the C is C and at the next position G (deletions aside; ties are no
  consensus). Every table lists every site.

## Running

```bash
lrpipe-methyl pileup reads.bam out.bed --ref ref.fa --modified-bases 5mC 5hmC --cpg --combine-strands -t 16
lrpipe-methyl pileup reads.cram phased/ --ref ref.fa --phased -t 16
lrpipe-methyl cpg --bam reads.bam --ref ref.fa --output-prefix sample --modsites-mode reference --threads 16
```

`pileup`'s output is plain text, or BGZF when its name ends in `.gz`; `-` or `stdout` writes
it to the standard output, as modkit does (but not with `--phased`, whose output is a
directory).

`cpg --combined-only` (lrpipe-methyl's own option) writes only `PREFIX.combined.bed.gz`, for
reads without `HP` tags, whose `hap1` and `hap2` tables would be empty; `lrpipe run` passes it
for such reads when it does not phase them.

`--haplotags LIST` (`cpg` and `pileup`, lrpipe-methyl's own option) takes the reads'
haplotypes from a haplotag list (WhatsHap's `haplotag --output-haplotag-list`, plain or
gzipped: read name, `H1`/`H2`/`none`, phase set, chromosome) in place of their `HP` tags, as
`haplotag` would have tagged them: a listed read's primary alignment gets its haplotype, its
supplementary and secondary alignments and every unlisted read none, whatever tags the input
carries. The tables are then those of the alignments `haplotag` writes. A read listed more than
once on a contig with different haplotypes (a list written with `--tag-supplementary`) is left
untagged, with a warning, and the log says how many of the primary alignments the list names,
warning when it is fewer than half. `lrpipe run` passes the list it makes from its own phased
SNVs. The list takes some 60 bytes per listed read in memory (24 and the name), and up to
twice that while it is read.

`--ref` is a FASTA with its `.fai`, plain or bgzipped (with its `.gzi` too, as `samtools faidx`
makes it). A table that reads the reference (`cpg`, `pileup --cpg` or `--modified-bases`) fails
when the reference lacks a contig that has reads, or has it of another length than the reads'
header does: such a reference is not the one the reads were aligned to (contigs without reads are not
checked, so a chr20 reference serves chr20 reads with a whole genome's header).
`--filter-threshold` takes modkit's forms, a probability for every base (`0.8`) or for one
primary base (`C:0.8`, with `A`, `C`, `G` or `T`), and `--filter-percentile` a fraction in
[0, 1); anything else is an error.

`tools/compare_bedmethyl.py` and `tools/compare_cpg.py` compare outputs row by row.

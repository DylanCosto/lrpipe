# lrpipe-cnv

A pure-Rust port of the [Spectre](https://github.com/fritzsedlazeck/Spectre) long-read CNV
caller (0.2.1, commit 90580943, MIT), aiming first at output identical to Spectre's, then at
improvements behind options. It runs end to end, with or without SNVs: on every reference run
its VCF and BED are byte-identical to Spectre's and its `.spc` identical apart from the order
of Python sets. Not ported: `--snfj` breakpoint support, population calling, and the plots.
Changes to the method are opt-in (see [Changes to Spectre's method](#changes-to-spectres-method)).

```bash
mosdepth -t 8 -x -b 1000 -Q 20 -n sample sample.bam      # or lrpipe-depth, or lrpipe-sv --depth
lrpipe-cnv CNVCaller -c sample.regions.bed.gz -s sample -d out -r ref.fa -l blacklist.bed \
    -v sample.lrpipe-snv.vcf.gz     # optional: Clair3 or lrpipe-snv SNVs, for copy-neutral depth and LoH
# or straight from the alignments (lrpipe-depth, identical windows; CRAM decoded with
# --cram-reference, default -r):
lrpipe-cnv CNVCaller --bam sample.bam -t 16 -s sample -d out -r ref.fa -l blacklist.bed
```

`--bam` gives output identical to the mosdepth route: HG002 HiFi in 31 s on 16 threads (28 s of
it counting), HG002 ONT CRAM in 199 s; each counting thread holds a 4-byte counter per base of
its contig (10-13 GB at 16 threads).

It takes Spectre's `CNVCaller` options, including its multi-letter short ones (`-mcl`, `-pc`,
...). On HG002 HiFi (31x whole genome) it takes 3.4 s and 0.35 GB where Spectre takes 184 s
and 1.1 GB; with lrpipe-snv's 5.4 M SNVs 22 s where Spectre takes 324 s; HG003 chr20 with SNVs
0.7 s where Spectre takes 7 s; making the `.mdr` of GRCh38 about 20 s where Spectre takes
400 s. It also reads the `AF` of current Clair3 and lrpipe-snv VCFs (`Number=A`), on which Spectre
fails.

## Changes to Spectre's method

Off by default, so a plain run stays Spectre-exact; `--improved` turns on all but the
mappability mask, which needs a second depth file.

| Option | Change |
|---|---|
| `--sex detect\|male\|female` | chrX outside the PARs and chrY are one copy in a male sample: their depth is doubled for calling and their copy numbers reported in copies of one. Spectre calls a male chrX as deletions and LoH end to end, and its loose deletion threshold partly came from those windows. A female sample's chrY is left out. |
| `--gc-correction`, `--gc-span` | Scale each window by the autosome median over the median of windows of its GC percentage, GC measured over 20 kb centred on the window (HiFi drops out with GC over a read's length, not a window's). |
| `--keep-amplifications` | Mask only short runs of windows above Spectre's outlier cutoff (fewer than 25), so amplifications are called; kept windows stay out of the threshold sample. |
| `--refresh-copy-numbers` | `CN`, `GT` and `CD` from each call's final windows, not the first raw run (Spectre's bug B3). |
| `--all-mapq-coverage`, `--min-mapq-ratio` | Mask windows whose MAPQ ≥ 20 depth is below half their all-MAPQ depth (`-Q 0`): ambiguously mapped, not deleted. |
| `--caller hmm` | Segment with a copy-number HMM (Viterbi over CN 0-6, depth emissions with a robust spread, heterozygous SNV allele fractions with `-v`) instead of Spectre's threshold runs; a segment is called when its median is beyond Spectre's strict thresholds. |
| `--sv-breakpoints` | A deletion or duplication from an SV caller's VCF (lrpipe-sv, Sniffles2) overlapping a call reciprocally by half gives it base-pair POS and END, and `SVID`. |
| `--vcf-style standard` | VCF 4.2 for strict readers. REF is taken from the configured FASTA and the floating-point score is `CNV_GQ`. `HO` and `CD` have one value each, END is the event's last base, and deletions have negative SVLEN. Existing `.mdr` metadata must match the FASTA's N runs. |

`crates/lrpipe-cnv-bench` scores calls against a truth set (GIAB's sequence-resolved SVs, where
duplications appear as insertions, or a BED), and `tools/bench_variants.sh` runs each change on
HG002 HiFi and ONT whole genomes against GIAB v5.0q (autosomes, events ≥ 50 kb) and on HG003
chr20 with injected events (`tools/synth_truth.bed`). Found / truth and true / calls:

| Run | HiFi DEL | HiFi DUP | ONT DEL | ONT DUP | synthetic DEL | synthetic DUP |
|---|---|---|---|---|---|---|
| Spectre | 0/21, 0/14 | 1/46, 1/2 | 3/21, 3/8 | 1/46, 1/5 | 3/6, 3/4 | 2/4, 2/2 |
| `--improved` | 3/21, 3/9 | 0/46, 0/2 | 3/21, 3/10 | 1/46, 1/6 | 3/6, 3/3 | 4/4, 4/4 |
| `--improved --caller hmm` | 2/21, 2/10 | 0/46, 0/1 | 3/21, 3/13 | 1/46, 1/11 | 3/6, 3/4 | 3/4, 3/4 |
| ... and `-v` (lrpipe-snv) | 1/21, 1/9 | 0/46, 0/1 | | | | |

The HMM does not beat the thresholds here; its breakpoints are closer (median error 3.6 kb
against 7.2 kb on ONT, 0.5 kb against 1.5 kb on the synthetic deletions), and SV breakpoints
halve the error of `--improved` on HiFi (3.2 kb against 6.1 kb), but only a few calls have a
matching SV.

Most GIAB deletions this size lie in segmental duplications (chr7:73 Mb, chr15:84 Mb,
chr16:18-33 Mb), where MAPQ ≥ 20 depth cannot see them, and 17 of the 21 are below the
100 kb minimum call length. The HiFi false deletions left are GC-rich subtelomeric stretches
at half depth that ONT covers normally.

### Spectre bugs fixed on purpose

- **Heterozygotes in either order.** A SNV is heterozygous when its genotype is `0/1`, `1/0`,
  `0|1` or `1|0`. Spectre takes only `0/1` and `0|1`, so phased calls with the ALT first counted
  for nothing (review finding LRP-011). lrpipe's own runs give it lrpipe-snv's unphased calls,
  which are always `0/1`.
- **Plain gzip depth.** `--coverage` reads a gzip file as well as a BGZF one (finding LRP-014).
- **CHM13's PARs.** `--sex` takes chrX's pseudo-autosomal regions from the reference as every
  lrpipe tool does (GRCh38, GRCh37 or T2T-CHM13, by chrX's length). lrpipe-cnv had its own table
  without CHM13, so on CHM13 a male's PARs counted as one copy and their depth was doubled
  (finding LRP-018). This is an extension of lrpipe's, not a Spectre bug: Spectre has no `--sex`.

## Layout

- `crates/pynum`: the Python, numpy and scipy numerics Spectre's results depend on, bit for
  bit: numpy's pairwise sums and the means, variances, medians and percentiles built on them;
  numpy rounding; Python float text; CPython's `random` (seeding, `sample`, `choice`); numpy's
  legacy RNG; scipy's `norm.sf`. Checked against the reference environment by
  `crates/pynum/tests/vectors.rs`. `norm.sf` calls the C library's `exp`, as scipy's does, so it
  is bit for bit only on the same `exp`: glibc 2.28 and macOS round a few values a unit or two
  in the last place away (the test allows 2 there).
- `crates/lrpipe-cnv`: the calling stages, quality scores and output files; `crates/lrpipe-cnv-cli`
  the `lrpipe-cnv` command. `tests/outputs.rs` compares whole runs' output files with
  Spectre's.
- `crates/lrpipe-cnv`'s stages: `tests/reference.rs` compares each with Spectre's
  state in the reference runs (`tools/reference_runs.tsv`; set `SPECTRE_BENCH` and
  `SPECTRE_REFERENCE` to run it): all of them match exactly so far, from the masked and
  normalized depth of every window, the random window sample and the thresholds to every
  raw candidate (16,275 on the HG002 HiFi whole genome) and every call, including the raw
  candidates' state after refinement changes them in place.
- `docs/SPECTRE_SPEC.md`: what Spectre's `CNVCaller` does, with file and line references,
  including the numerics, the random streams, the output formats and its bugs.
- `crates/lrpipe-cnv/src/extensions.rs` (with `sex.rs`, `gc.rs`, `hmm.rs`, `breakpoints.rs`): the
  opt-in changes;
  `crates/lrpipe-cnv-bench` the accuracy benchmark.
- `tools/`: the reference setup (below).

## Reference runs

Spectre's results depend on its Python stack, so expected outputs come from a pinned
environment (`tools/spectre-env.txt`: Python 3.10, numpy 1.24.3, pandas 2.0.1, scipy 1.10.1,
pysam 0.22.0, matplotlib 3.7.1, no bottleneck). numpy 2 crashes Spectre and pandas 3 changes
its sums.

```bash
micromamba create -p ENV -c conda-forge -c bioconda --override-channels --file tools/spectre-env.txt
tools/fetch_spectre.sh target/spectre-0.2.1   # not vendored: parts of Spectre are GPL-3.0
SPECTRE_SRC=target/spectre-0.2.1 ENV/bin/python tools/run_spectre.py DUMP_DIR -- CNVCaller ...
```

`tools/run_spectre.py` seeds numpy's global RNG (Spectre's dip test draws from it unseeded)
and writes Spectre's state after each stage, with float arrays summarized by a SHA-256 of
their bits, so the port can be compared stage by stage. `tools/inject_cnvs.py` adds known
deletions, duplications and a loss of heterozygosity to real mosdepth and SNV inputs.
`tools/pynum_vectors.py` regenerates the numerics test vectors and
`tools/htslib_float_vectors.py` the VCF float text vectors, `tools/dip_vectors.py` the dip-test
vectors (Spectre's VCF goes through
`bcftools sort`, so its floats read as htslib prints them). `tools/spc_canonical.py` sorts the
set-valued lists of a `.spc` for comparison.

## Notes on Spectre 0.2.1

- It skips SNV records whose `AF` is declared `Number=A`, as current Clair3 (and lrpipe-snv)
  declare it, and then fails. Reference runs rewrite the header to `Number=1`; lrpipe-cnv reads
  single-value `Number=A` directly, and on the original Clair3 VCF its output is identical to
  Spectre's on the rewritten one.
- Its dip test draws from numpy's global generator without seeding it, so with SNVs its
  output can change between runs. lrpipe-cnv seeds the same generator (`--dip-seed`, default
  0, as the reference runs). With heterozygous depths of a germline sample the test finds one
  mode by a wide margin whatever the draws; the seed matters only for multimodal depths.
- With `--snv` it normalizes by the depth of heterozygous SNVs; without, by the genome-wide
  median depth, so a file with reads on one contig only (and zero depth elsewhere) calls
  nothing unless the windows are restricted to that contig.
- Its outlier cutoff (mean + 5 SD) comes from every contig in the mosdepth file. With reads
  on one contig only, as in the chr20 test BAM, the zero-depth contigs pull the cutoff down
  to about 25x and most chr20 windows are masked; such runs are still parity tests, but
  their calls say little about Spectre.
- The same outlier cutoff removes high-level amplifications before calling: on chr20 with
  an injected 2.5x (CN 5) amplification, the cutoff is about 2.2x the median and most of its
  windows are masked, so it is never a candidate. Single-copy and two-copy duplications in
  well-covered regions are called.
- The dip-test code Spectre vendors (`spectre/libs/unidip`) is GPL-3.0; lrpipe-cnv
  implements the dip test and UniDip from their papers (Hartigan and Hartigan 1985, with
  Hartigan's AS 217; Maurus and Plant 2016) and checks them against that code run as a black
  box (`tools/dip_vectors.py`, `crates/lrpipe-cnv/tests/dip.rs`). Every p-value matches, dips
  agree to within an ulp, and UniDip agrees on every input with one mode. On multimodal inputs
  the modal interval, and so UniDip's recursion, follow the papers rather than that code's
  bookkeeping and can differ (25 of 33 such test inputs agree).

## License

MIT (`LICENSE`). Spectre's license is in `LICENSE-SPECTRE`.

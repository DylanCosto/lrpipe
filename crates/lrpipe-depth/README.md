# lrpipe-depth

Binned read depth from BAM and CRAM in pure Rust, with output identical to
[mosdepth](https://github.com/brentp/mosdepth) 0.3.14.

It takes mosdepth's options and writes mosdepth's files, so tools that read mosdepth output
(Spectre, for example) can use it unchanged. It is also a library: a tool that already reads the
alignments, such as an SV caller, can count depth in the same pass.

```bash
lrpipe-depth -t 8 -n -x -Q 20 -b 1000 sample sample.bam
```

## Identical to mosdepth

Every output lrpipe-depth writes is byte-identical to mosdepth 0.3.14's (gzipped files after
decompression), for these options:

| Option | Meaning |
|---|---|
| `-b/--by` | integer window size, or a BED file (plain or gzipped) |
| `-n/--no-per-base` | required: per-base output is not implemented |
| `-x/--fast-mode` | count each read's whole span |
| `-a/--fragment-mode` | count whole fragments of proper pairs |
| `-Q/--mapq` | minimum mapping quality |
| `-F/--flag`, `-i/--include-flag` | exclude / require flag bits (default `-F 1796`) |
| `-l/--min-frag-len`, `-u/--max-frag-len` | insert-size limits |
| `-T/--thresholds` | bases at or above each depth, per region |
| `-m/--use-median` | median instead of mean per region |
| `-R/--read-groups` | only these read groups |
| `-c/--chrom` | one chromosome |
| `-f/--fasta`, `$REF_PATH` | reference for CRAM |
| `MOSDEPTH_PRECISION` | decimal places |

Outputs: `{prefix}.mosdepth.summary.txt`, `{prefix}.mosdepth.global.dist.txt`, and with `--by`
`{prefix}.regions.bed.gz` (+ `.csi`), `{prefix}.mosdepth.region.dist.txt`, and with `-T`
`{prefix}.thresholds.bed.gz` (+ `.csi`).

The counting follows mosdepth's source line by line (cited in the code), including its quirks:

- In the default mode, deletions and skips are not counted, insertions don't split a read, and
  where the mates of a proper pair overlap, the overlap is counted once (same algorithm, same
  ordering of the start/end events).
- A contig whose reads were all filtered out is reported (zeros); a contig with no reads at
  all is not, though its region lines are still written.
- Region means add `depth / length` base by base, so the floating-point result matches bit for
  bit; medians are exact up to 65535. Numbers are rounded as C's `%#.*f` does (ties to even).
- With a BED file and no `-T`, contigs without BED regions are skipped entirely (no summary
  line, not in the totals). A BED file ends at its first empty line.

### How it is checked

- `tools/parity.sh` runs mosdepth and lrpipe-depth on 39 option combinations over the fixtures in
  `tests/data` and diffs every output. The fixtures are excerpts of real HiFi (chr22), ONT
  (chr20, chrX/Y) and Illumina (chr20) data with sequences stripped, plus synthetic BAM, CRAM
  and CSI-indexed BAM files of edge cases (`tools/synth.py`): overlapping pairs with split
  CIGARs, reads that are only clips or insertions, placed unmapped reads, every excluded flag,
  read groups, contig edges, a contig without reads.
- `cargo test` checks lrpipe-depth against the recorded hashes of mosdepth's outputs
  (`tests/data/expected.sha256`), so it needs no mosdepth.
- Index queries of `regions.bed.gz` and `thresholds.bed.gz` with `tabix` return the same lines
  as reading the file.
- Whole files: outputs are identical for HG003 chr20 ONT as BAM and as CRAM (fast and default
  modes, 1 to 16 threads) and for HG002 HiFi whole-genome BAM.

### Differences

- Per-base output (`*.per-base.bed.gz`, `--d4`) and `--quantize` are not implemented.
- `-t` is the total thread count: contigs with reads are counted in parallel, and the threads
  left over decompress BAM blocks or decode CRAM containers in parallel. In mosdepth it is the
  number of decompression threads.
- For BED regions that extend past the end of a contig, mosdepth's `-T` counts read memory past
  the end of its depth array, so its output there is undefined; lrpipe-depth counts only the bases
  on the contig. (Region means and medians are not affected and match.)
- Reads that extend past the end of a contig are counted up to the end instead of crashing.
- The `.csi` indexes are equivalent but not byte-identical to htslib's.

## Speed and memory

Wall time and peak memory on a 16-core machine (31 GB RAM), against the official mosdepth
0.3.14 binary with the same `-t`. Outputs were identical in every run.

| Input | Options | `-t` | mosdepth | lrpipe-depth |
|---|---|---|---|---|
| HG002 HiFi whole genome, BAM (39 GB, 31×) | `-n -x -Q 20 -b 1000` | 8 | 55.0 s, 2.0 GB | 39.0 s, 6.2 GB |
| | | 16 | 44.3 s, 2.0 GB | 28.5 s, 10.1 GB |
| HG003 ONT chr20, BAM (2.5 GB) | `-n -x -Q 20 -b 1000` | 1 | 14.4 s, 0.5 GB | 13.7 s, 0.4 GB |
| | | 4 | 5.6 s, 0.5 GB | 4.6 s, 0.4 GB |
| | | 8 | 4.4 s, 0.5 GB | 3.3 s, 0.4 GB |
| | `-n -b 500` | 4 | 7.7 s, 0.5 GB | 6.1 s, 0.5 GB |
| | | 8 | 6.5 s, 0.5 GB | 4.9 s, 0.5 GB |
| HG003 ONT chr20, CRAM (1.5 GB) | `-n -x -Q 20 -b 1000` | 8 | 3.3 s, 0.6 GB | 12.5 s, 0.4 GB |
| | | 16 | 3.3 s, 0.6 GB | 8.8 s, 0.5 GB |

(mosdepth's 1-thread figure is its default `-t 0`.)

Each contig being counted holds 4 bytes per base (about 1 GB for chr1), so memory grows with
the number of contigs counted at once: up to `-t`, and never more than the number of contigs
with reads. Contigs without reads cost nothing.

CRAM is slower than in mosdepth: pasta (like noodles, which it began from) decodes complete records (bases are rebuilt from the
reference, qualities and tags are decoded), where htslib can skip the fields depth does not
need. Containers are decoded in parallel to make up some of the difference.

## As a library

Count depth from records you already read: fill an `Alignment` per record and push it.

```rust
use lrpipe_depth::{Alignment, Cigar, ContigCounter, ContigReport, Filter, Mode, ReportOptions, Scratch, Writer};

let mut counter = ContigCounter::new(contig_len, Mode::Fast, Filter { min_mapq: 20, ..Filter::default() });
// for each record of the contig, in coordinate order:
counter.push(&Alignment {
    tid, pos, flags, mapq, mate_tid, mate_pos, template_len,
    name,
    cigar: Cigar::Packed(raw_bam_cigar), // or Cigar::Ops(&[(op, len)])
    read_group: None,
});
let coverage = counter.finish();          // per-base depth in coverage.depth
let report = ContigReport::new(contig_name, contig_len as u32, &coverage, &opts, &mut Scratch::new(&opts));
writer.write_contig(&report)?;            // Writer::new(prefix, &opts, max_contig_len), in contig order
```

`Alignment` borrows plain fields rather than a pasta record, so callers on other versions of
pasta, or with their own decoders, can use it without conversion.

To count a contig in parallel pieces (lrpipe-sv reads each contig in 2 Mbp pieces of reads, by
start position), give each piece a `ContigCounter::chunk`, which stores only the bases its reads
touch, and combine the pieces in order with `append`. The result is identical to one counter fed
every record, including overlapping mates split between pieces (`tests/chunks.rs`).
`ContigCounter::merge` adds counters in any order, which is exact when no mate overlaps need
correcting (`Mode::Fast`, or reads that are not proper pairs, which covers long reads).

`lrpipe_depth::run::run` is the whole-file driver the command line uses, and
`lrpipe_depth::run::count_contigs` counts a chosen set of contigs the same way.

## Development

```bash
cargo test                        # unit tests and parity against recorded mosdepth output
MOSDEPTH=/path/to/mosdepth-0.3.14 TABIX=tabix tools/parity.sh            # compare live
MOSDEPTH=/path/to/mosdepth-0.3.14 tools/parity.sh --record               # after changing cases
SAMTOOLS=samtools tools/make_fixtures.sh                                  # rebuild tests/data
```

Parity cases are listed in `tools/cases.tsv`.

## License

MIT (`LICENSE`). The counting logic is ported from mosdepth, MIT-licensed by Brent Pedersen
(`LICENSE-MOSDEPTH`).

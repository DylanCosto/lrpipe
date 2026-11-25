#!/usr/bin/env python3
"""Injects known copy-number changes and a loss of heterozygosity into real inputs.

Usage: inject_cnvs.py REGIONS.bed.gz SNV.vcf.gz OUT_DIR

Scales the mean depth of mosdepth windows inside each event (writing
OUT_DIR/mosdepth.regions.bed.gz with a CSI index), and turns heterozygous SNVs inside the
LoH interval into homozygous ALT calls (OUT_DIR/snv.vcf.gz with a tabix index). The events
are chosen to exercise Spectre's calling, merging, size filter and contig-start handling on
chr20, in well-covered regions (pericentromeric and contig-end windows have too little depth
for an injected change to show); they are not meant to look like real biology. Runs in the Spectre reference
environment (needs pysam).
"""

import gzip
import os
import sys

import pysam

# (contig, start, end, depth factor): 0-based, half-open.
EVENTS = [
    ("chr20", 100_000, 500_000, 0.5),  # deletion at the contig start
    ("chr20", 10_000_000, 10_400_000, 0.5),  # heterozygous deletion
    ("chr20", 20_000_000, 20_150_000, 0.0),  # homozygous deletion
    ("chr20", 32_000_000, 32_600_000, 1.5),  # duplication (CN 3)
    ("chr20", 40_000_000, 40_250_000, 2.5),  # amplification (CN 5): Spectre masks it as outliers
    ("chr20", 55_000_000, 55_300_000, 2.0),  # duplication (CN 4)
    ("chr20", 45_000_000, 45_080_000, 0.5),  # below the 100 kb minimum
    ("chr20", 50_000_000, 50_200_000, 0.5),  # two deletions 30 kb apart, to merge
    ("chr20", 50_230_000, 50_450_000, 0.5),
    ("chr20", 60_000_000, 60_400_000, 1.5),  # duplication (CN 3)
]
LOH = ("chr20", 25_000_000, 28_000_000)


def factor(contig, start, end):
    for c, s, e, f in EVENTS:
        if c == contig and start < e and end > s:
            return f
    return 1.0


def write_regions(src, dst):
    with gzip.open(src, "rt") as fin, open(dst, "w") as fout:
        for line in fin:
            contig, start, end, mean = line.rstrip("\n").split("\t")
            f = factor(contig, int(start), int(end))
            if f != 1.0:
                mean = f"{float(mean) * f:.2f}"
            fout.write(f"{contig}\t{start}\t{end}\t{mean}\n")
    pysam.tabix_index(dst, preset="bed", csi=True, force=True)


def write_snvs(src, dst):
    contig, s, e = LOH
    with pysam.VariantFile(src) as fin, pysam.VariantFile(dst, "w", header=fin.header) as fout:
        for rec in fin:
            if rec.chrom == contig and s <= rec.pos - 1 < e:
                sample = rec.samples[0]
                if sample["GT"] == (0, 1):
                    sample["GT"] = (1, 1)
                    sample["AF"] = 0.97
            fout.write(rec)
    pysam.tabix_index(dst, preset="vcf", force=True)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    regions, snvs, out = sys.argv[1:]
    os.makedirs(out, exist_ok=True)
    write_regions(regions, os.path.join(out, "mosdepth.regions.bed"))
    write_snvs(snvs, os.path.join(out, "snv.vcf"))


if __name__ == "__main__":
    main()

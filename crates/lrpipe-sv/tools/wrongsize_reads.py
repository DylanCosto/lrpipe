#!/usr/bin/env python3
"""For false negatives with a wrong-size PASS call nearby, compare per-read indel sizes to truth.

usage: wrongsize_reads.py TRUVARI_DIR CALLS.vcf BAM

For each such locus, every primary read (MAPQ >= 20) spanning the locus gets three sizes
(insertions for an INS truth, deletions for a DEL):
  largest  - its largest single indel >= 45 bp (what Sniffles2 mostly sees)
  joined45 - the sum of its indels >= 45 bp in the window (full per-read joining)
  joined10 - the sum of its indels >= 10 bp in the window (including pieces Sniffles2 ignores)
and reports which measure's per-read median is closest to the truth size.
"""
import bisect
import gzip
import statistics
import sys
from collections import Counter, defaultdict

import pysam

WINDOW = 1000
FLANK = 300


def info(i, k):
    for kv in i.split(";"):
        if kv.startswith(k + "="):
            return kv[len(k) + 1:]
    return None


def main():
    d, calls_path, bam_path = sys.argv[1:4]
    calls = defaultdict(list)
    with open(calls_path) as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.split("\t")
            ty = info(t[7], "SVTYPE")
            if t[6] == "PASS" and ty in ("INS", "DEL"):
                calls[(t[0], ty)].append((int(t[1]), abs(int(info(t[7], "SVLEN")))))
    for v in calls.values():
        v.sort()
    bam = pysam.AlignmentFile(bam_path)
    best = Counter()
    rows = []
    with gzip.open(f"{d}/fn.vcf.gz", "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.split("\t")
            chrom, pos = t[0], int(t[1])
            diff = len(t[4]) - len(t[3])
            ty, size = ("INS" if diff > 0 else "DEL"), abs(diff)
            lst = calls[(chrom, ty)]
            k = bisect.bisect_left(lst, (pos - WINDOW,))
            near = []
            while k < len(lst) and lst[k][0] <= pos + WINDOW:
                near.append(lst[k])
                k += 1
            ratio = lambda s: min(s, size) / max(s, size)
            if not near or any(ratio(c[1]) >= 0.7 for c in near):
                continue
            call = min(near, key=lambda c: abs(c[0] - pos))
            lo = min(pos, call[0]) - FLANK
            hi = max(pos + (size if ty == "DEL" else 0), call[0] + (call[1] if ty == "DEL" else 0)) + FLANK
            op = 1 if ty == "INS" else 2
            largest, j45, j10 = [], [], []
            for r in bam.fetch(chrom, lo, hi):
                if r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
                    continue
                if r.reference_start > lo or r.reference_end < hi:
                    continue
                rp = r.reference_start
                sizes = []
                for o, n in r.cigartuples:
                    if o == op and lo <= rp <= hi:
                        sizes.append(n)
                    if o in (0, 2, 3, 7, 8):
                        rp += n
                largest.append(max([s for s in sizes if s >= 45], default=0))
                j45.append(sum(s for s in sizes if s >= 45))
                j10.append(sum(s for s in sizes if s >= 10))
            if not largest:
                continue
            # Per-read medians over reads that show the SV at all.
            med = {}
            for name, v in (("largest", largest), ("joined45", j45), ("joined10", j10)):
                nz = [x for x in v if x > 0]
                med[name] = statistics.median(nz) if nz else 0
            win = max(med, key=lambda m: ratio(med[m]) if med[m] else 0)
            ok = ratio(med[win]) >= 0.7 if med[win] else False
            best[(ty, win if ok else "none within 0.7")] += 1
            rows.append((chrom, pos, ty, size, call[1], med["largest"], med["joined45"], med["joined10"]))
    print("closest per-read measure to the truth size (within 0.7x):")
    for k, v in sorted(best.items()):
        print(f"  {v:5d}  {k[0]} {k[1]}")
    print("examples: chrom pos type truth_size call_size | median largest joined45 joined10")
    for r in rows[:25]:
        print("  ", *r)


if __name__ == "__main__":
    main()

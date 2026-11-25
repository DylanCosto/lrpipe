#!/usr/bin/env python3
"""Break down Truvari results by SV type, size and tandem-repeat context.

usage: sv_breakdown.py TRUVARI_DIR TR.bed.gz TR_GE101.bed.gz [OTHER_TRUVARI_DIR]

With OTHER_TRUVARI_DIR (another caller benchmarked on the same truth), also counts how many
of this caller's false negatives the other caller found.
"""
import bisect
import gzip
import sys
from collections import defaultdict

SIZE_BINS = [(50, 100), (100, 300), (300, 1000), (1000, 5000), (5000, 10**9)]


def size_bin(n):
    for lo, hi in SIZE_BINS:
        if lo <= n < hi:
            return f"{lo}-{hi if hi < 10**9 else 'inf'}"
    return "<50"


def read_bed(path):
    by_chrom = defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            by_chrom[c].append((int(s), int(e)))
    merged = {}
    for c, iv in by_chrom.items():
        iv.sort()
        out = []
        for s, e in iv:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[c] = ([s for s, _ in out], [e for _, e in out])
    return merged


def overlaps(bed, chrom, start, end):
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    k = bisect.bisect_right(starts, end) - 1
    return k >= 0 and ends[k] > start


def records(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, info = t[0], int(t[1]), t[3], t[4], t[7]
            svtype = next((kv[7:] for kv in info.split(";") if kv.startswith("SVTYPE=")), None)
            svlen = next((kv[6:] for kv in info.split(";") if kv.startswith("SVLEN=")), None)
            if not alt.startswith("<") and "[" not in alt and "]" not in alt:
                d = len(alt) - len(ref)
                svtype = svtype or ("INS" if d > 0 else "DEL")
                size = abs(d)
            else:
                size = abs(int(svlen.split(",")[0])) if svlen and svlen != "." else 0
            if svtype not in ("INS", "DEL"):
                svtype = svtype or "OTHER"
            start = pos
            end = pos + (size if svtype == "DEL" else 1)
            yield (chrom, pos, ref, alt), chrom, start, end, svtype, size


def context(tr, tr101, chrom, s, e):
    if overlaps(tr101, chrom, s, e):
        return "TR>=101bp"
    if overlaps(tr, chrom, s, e):
        return "TR<101bp"
    return "non-TR"


def main():
    d = sys.argv[1]
    tr, tr101 = read_bed(sys.argv[2]), read_bed(sys.argv[3])
    other_tp = set()
    if len(sys.argv) > 4:
        other_tp = {k for k, *_ in records(f"{sys.argv[4]}/tp-base.vcf.gz")}

    counts = defaultdict(lambda: defaultdict(int))
    for name, path in (("TP", "tp-base"), ("FN", "fn"), ("FP", "fp")):
        for key, chrom, s, e, svtype, size in records(f"{d}/{path}.vcf.gz"):
            for group in (
                ("type", svtype),
                ("size", f"{svtype} {size_bin(size)}"),
                ("context", context(tr, tr101, chrom, s, e)),
                ("type x context", f"{svtype} {context(tr, tr101, chrom, s, e)}"),
            ):
                counts[group][name] += 1
                if name == "FN" and key in other_tp:
                    counts[group]["FN_other"] += 1

    print(f"{'group':<16} {'class':<22} {'truth':>7} {'TP':>7} {'FN':>7} {'recall':>7} {'FP':>6}"
          + (f" {'FN found by other':>18}" if other_tp else ""))
    for kind in ("type", "size", "context", "type x context"):
        for (k, cls), c in sorted((g, c) for g, c in counts.items() if g[0] == kind):
            truth = c["TP"] + c["FN"]
            recall = c["TP"] / truth if truth else float("nan")
            line = f"{k:<16} {cls:<22} {truth:>7} {c['TP']:>7} {c['FN']:>7} {recall:>7.3f} {c['FP']:>6}"
            if other_tp:
                line += f" {c['FN_other']:>18}"
            print(line)
        print()


if __name__ == "__main__":
    main()

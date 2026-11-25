#!/usr/bin/env python3
"""Classify Truvari false negatives by the nearest candidate in an all-candidates VCF.

usage: fn_diagnose.py TRUVARI_DIR NOQC.vcf TR_GE101.bed.gz

For each false negative, looks for candidates of the same type within WINDOW bp in the
`--no-qc` VCF and reports: no candidate, a candidate removed by a filter (which one), or a
PASS candidate that was not matched (with its size ratio and distance).
"""
import bisect
import gzip
import sys
from collections import Counter, defaultdict

WINDOW = 1000
SIZE_OK = 0.7


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def info_get(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


def truth_records(path):
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t")
            d = len(t[4]) - len(t[3])
            yield t[0], int(t[1]), ("INS" if d > 0 else "DEL"), abs(d)


def candidates(path):
    by = defaultdict(list)
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t")
            svtype = info_get(t[7], "SVTYPE")
            if svtype not in ("INS", "DEL"):
                continue
            svlen = abs(int(info_get(t[7], "SVLEN")))
            support = int(info_get(t[7], "SUPPORT") or 0)
            by[(t[0], svtype)].append((int(t[1]), svlen, t[6], support))
    for v in by.values():
        v.sort()
    return by


def read_bed(path):
    by = defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            c, s, e = line.split("\t")[:3]
            by[c].append((int(s), int(e)))
    for v in by.values():
        v.sort()
    return {c: ([s for s, _ in v], [e for _, e in v]) for c, v in by.items()}


def in_bed(bed, c, p):
    if c not in bed:
        return False
    s, e = bed[c]
    k = bisect.bisect_right(s, p) - 1
    # Intervals may overlap; check a few back.
    return any(s[j] <= p < e[j] for j in range(max(0, k - 3), k + 1))


def main():
    d, noqc, tr = sys.argv[1], sys.argv[2], read_bed(sys.argv[3])
    cands = candidates(noqc)
    classes = Counter()
    filters = Counter()
    by_ctx = defaultdict(Counter)
    for chrom, pos, svtype, size in truth_records(f"{d}/fn.vcf.gz"):
        near = [c for c in cands.get((chrom, svtype), [])
                if abs(c[0] - pos) <= WINDOW] if (chrom, svtype) in cands else []
        lst = cands.get((chrom, svtype), [])
        k = bisect.bisect_left(lst, (pos - WINDOW,))
        near = []
        while k < len(lst) and lst[k][0] <= pos + WINDOW:
            near.append(lst[k])
            k += 1
        ratio = lambda c: min(c[1], size) / max(c[1], size) if max(c[1], size) else 0
        passing = [c for c in near if c[2] == "PASS"]
        sized = [c for c in near if ratio(c) >= SIZE_OK]
        if not near:
            cls = "no candidate within 1 kb"
        elif any(c[2] == "PASS" and ratio(c) >= SIZE_OK for c in near):
            best = min((c for c in passing if ratio(c) >= SIZE_OK), key=lambda c: abs(c[0] - pos))
            cls = "PASS, size ok, >500 bp away" if abs(best[0] - pos) > 500 else "PASS, size ok (matched elsewhere)"
        elif passing:
            cls = "PASS nearby, size off (<0.7)"
        elif sized:
            best = max(sized, key=lambda c: c[3])
            cls = "filtered, size ok"
            filters[best[2]] += 1
        else:
            cls = "filtered, size off"
        classes[cls] += 1
        ctx = "TR>=101bp" if in_bed(tr, chrom, pos) else "other"
        by_ctx[ctx][cls] += 1
    total = sum(classes.values())
    print(f"{total} false negatives")
    for cls, n in classes.most_common():
        print(f"  {n:6d} {n / total:6.1%}  {cls}   (TR>=101bp {by_ctx['TR>=101bp'][cls]}, other {by_ctx['other'][cls]})")
    print("filters on size-matching candidates:")
    for f, n in filters.most_common():
        print(f"  {n:6d}  {f}")


if __name__ == "__main__":
    main()

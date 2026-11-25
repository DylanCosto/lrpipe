#!/usr/bin/env python3
"""Compares two pb-CpG-tools-style bed files (only contigs both have): rows only in one, and
per column how many shared rows differ, with examples. Optional third argument: contig."""
import gzip, sys
from collections import Counter
COLS = ["chrom","begin","end","mod_score","type","cov","c7","c8","c9","c10"]
def rows(path, contig=None):
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            c = line.rstrip("\n").split("\t")
            if contig and c[0] != contig: continue
            out[(c[0], int(c[1]))] = c
    return out
contig = sys.argv[3] if len(sys.argv) > 3 else None
a, b = rows(sys.argv[1], contig), rows(sys.argv[2], contig)
oa = [k for k in a if k not in b]; ob = [k for k in b if k not in a]
print(f"{len(a)} vs {len(b)} rows; {len(oa)} only in first, {len(ob)} only in second")
for k in oa[:3]: print("  only first: ", "\t".join(a[k]))
for k in ob[:3]: print("  only second:", "\t".join(b[k]))
d = Counter(); ex = {}
for k in a.keys() & b.keys():
    for i, (x, y) in enumerate(zip(a[k], b[k])):
        if x != y:
            d[COLS[i]] += 1; ex.setdefault(COLS[i], (a[k], b[k]))
for col, n in d.most_common():
    print(f"  {col}: {n} rows differ"); x, y = ex[col]
    print("    first: ", "\t".join(x)); print("    second:", "\t".join(y))

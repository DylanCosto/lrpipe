#!/usr/bin/env python3
"""Compares two bedMethyl files row by row (keyed by contig, start, code, strand):
rows only in one, and per column how many shared rows differ, with examples."""
import gzip, sys
from collections import Counter

COLS = ["chrom","start","end","code","score","strand","s2","e2","color","valid","pct","mod","canonical","other","delete","fail","diff","nocall"]

def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            out[(c[0], int(c[1]), c[3], c[5])] = c
    return out

a, b = rows(sys.argv[1]), rows(sys.argv[2])
only_a = [k for k in a if k not in b]
only_b = [k for k in b if k not in a]
print(f"{len(a)} vs {len(b)} rows; {len(only_a)} only in first, {len(only_b)} only in second")
for k in only_a[:3]: print("  only first: ", "\t".join(a[k]))
for k in only_b[:3]: print("  only second:", "\t".join(b[k]))
diff = Counter(); ex = {}
for k in a.keys() & b.keys():
    for i, (x, y) in enumerate(zip(a[k], b[k])):
        if x != y:
            diff[COLS[i]] += 1
            ex.setdefault(COLS[i], (a[k], b[k]))
for col, n in diff.most_common():
    print(f"  {col}: {n} rows differ")
    x, y = ex[col]
    print("    first: ", "\t".join(x)); print("    second:", "\t".join(y))
order_a = list(a.keys()); order_b = [k for k in b.keys() if k in a]
oa = [k for k in order_a if k in b]
print("same row order" if oa == order_b else "row order differs")

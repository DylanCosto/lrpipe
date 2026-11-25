#!/usr/bin/env python3
"""Prints a .spc.gz with every merged_sample_references list sorted, so two runs can be
compared: Spectre writes those Python sets in hash-randomized order.

Usage: spc_canonical.py FILE.spc.gz > canonical.txt
"""

import gzip
import re
import sys


def main():
    lines = gzip.open(sys.argv[1], "rt").read().split("\n")
    out, i = [], 0
    while i < len(lines):
        out.append(lines[i])
        if lines[i].endswith('"merged_sample_references": ['):
            j = i + 1
            while not re.match(r"^\s*\]", lines[j]):
                j += 1
            items = sorted(line.rstrip(",") for line in lines[i + 1 : j])
            out.extend(item + ("," if k < len(items) - 1 else "") for k, item in enumerate(items))
            i = j
            continue
        i += 1
    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main()

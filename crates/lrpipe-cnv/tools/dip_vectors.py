#!/usr/bin/env python3
"""Writes crates/lrpipe-cnv/tests/data/dip.json.gz: inputs and outputs of the dip test and UniDip
that Spectre vendors (spectre/libs/unidip, GPL-3.0), run as a black box.

Usage: SPECTRE_SRC=... dip_vectors.py OUT.json.gz

lrpipe-cnv implements the dip test and UniDip from their papers and checks them against these
outputs; the GPL code is only executed here, never copied. Results are recorded for:
- dip_fn(data): the dip and the modal interval indices diptst reports;
- diptst(data) with numpy's global RNG seeded first: dip, p-value, modal indices;
- UniDip(data).run() with numpy seeded first: the intervals.
"""

import gzip
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.environ["SPECTRE_SRC"])
from spectre.libs.unidip.dip import diptst, dip_fn  # noqa: E402
from spectre.libs.unidip.unidip import UniDip  # noqa: E402


def bits(x):
    return None if x is None else struct.pack("<d", float(x)).hex()


def datasets(rng):
    out = []
    for n in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]:
        for _ in range(4):
            out.append(rng.integers(0, 12, n).tolist())
    out.append([5] * 30)
    out.append([1, 2, 3, 4])
    out.append([1, 1, 2, 2, 3, 3, 4, 4, 5])
    for n in [30, 100, 300, 1000, 5000, 20000, 60000]:
        out.append(rng.poisson(33, n).tolist())
        out.append(np.concatenate([rng.poisson(16, n // 3), rng.poisson(33, n - n // 3)]).tolist())
    for n in [50, 200, 1000]:
        out.append(np.concatenate([rng.poisson(10, n // 2), rng.poisson(60, n // 2)]).tolist())
        out.append(np.round(rng.normal(30, 6, n), 3).tolist())
        out.append(np.round(np.concatenate([rng.normal(10, 2, n // 2), rng.normal(30, 2, n // 2)]), 3).tolist())
        out.append(rng.integers(0, 100, n).tolist())
    return out


def main():
    rng = np.random.default_rng(11)
    records = []
    for k, data in enumerate(datasets(rng)):
        arr = np.array(data, dtype=float)
        rec = {"data": data}
        d, (_, idxs, left, _, right, _) = dip_fn(arr)
        rec["dip"] = {"d": bits(d), "idxs": len(idxs), "left": len(left), "right": len(right)}
        tests = []
        for seed in [0, 1, 7]:
            np.random.seed(seed)
            d, p, mod = diptst(arr)
            tests.append({"seed": seed, "d": bits(d), "p": bits(p), "mod": [int(mod[0]), int(mod[1])]})
        rec["diptst"] = tests
        runs = []
        for seed in [0, 7]:
            np.random.seed(seed)
            runs.append({"seed": seed, "intervals": [[int(a), int(b)] for a, b in UniDip(data).run()]})
        rec["unidip"] = runs
        records.append(rec)
    with gzip.open(sys.argv[1], "wt") as f:
        json.dump(records, f)


if __name__ == "__main__":
    main()

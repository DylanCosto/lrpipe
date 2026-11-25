#!/usr/bin/env python3
"""Writes crates/pynum/tests/data/vectors.json.gz: inputs and results of the Python, numpy and
pandas operations pynum reproduces, from the Spectre reference environment.

Usage: pynum_vectors.py OUT.json.gz

Float results are stored as their IEEE-754 bit patterns (hex), so the Rust tests compare
bits, not approximately.
"""

import gzip
import hashlib
import json
import math
import random
import struct
import sys

import numpy as np
import pandas as pd
import scipy.stats


def bits(x):
    return struct.pack("<d", float(x)).hex()


def arrays(rng):
    """Arrays whose sum depends on the summation order, of lengths around numpy's blocks."""
    out = []
    for n in list(range(0, 20)) + [63, 64, 65, 127, 128, 129, 130, 255, 256, 257, 1000, 1023, 8191, 8192, 8193, 20000]:
        a = [rng.uniform(-1, 1) * 10 ** rng.randint(-8, 8) for _ in range(n)]
        out.append(a)
    out.append([-0.0])
    out.append([-0.0, -0.0])
    out.append([0.1] * 1000)
    out.append([1e16, 1.0, -1e16] * 50)
    return out


def with_nans(rng, a):
    return [math.nan if rng.random() < 0.2 else x for x in a]


def main():
    rng = random.Random(12345)
    data = arrays(rng)
    data += [with_nans(rng, a) for a in data]
    cases = {"sum": [], "nan": [], "percentile": []}
    for i, a in enumerate(data):
        x = np.array(a, dtype=np.float64)
        if not np.isnan(x).any():
            cases["sum"].append({"a": i, "sum": bits(np.sum(x)), "mean": bits(np.mean(x)) if len(a) else None,
                                 "std": bits(np.std(x)) if len(a) else None, "var": bits(np.var(x)) if len(a) else None})
            if len(a):
                for q in [0, 10, 10.0000001, 25, 50, 90, 98, 99.9, 100]:
                    cases["percentile"].append({"a": i, "q": bits(q), "r": bits(np.percentile(x, q))})
        if len(a) and not np.isnan(x).all():
            s = pd.Series(x)
            cases["nan"].append({"a": i, "nanmean": bits(np.nanmean(x)), "nanstd": bits(np.nanstd(x)),
                                 "nanmedian": bits(np.nanmedian(x)), "series_mean": bits(s.mean()),
                                 "series_std": bits(s.std()) if s.count() > 1 else None})

    rounds = []
    for _ in range(3000):
        x = rng.uniform(-100, 100)
        for d in [0, 1, 2, 3]:
            rounds.append({"x": bits(x), "d": d, "r": bits(round(np.float64(x), d))})
    for x in [0.075, 0.125, 0.375, 2.675, 1.005, 0.5, 1.5, 2.5, -0.5, -2.5, 1e300, -0.0]:
        for d in [0, 1, 2]:
            rounds.append({"x": bits(x), "d": d, "r": bits(round(np.float64(x), d))})

    reprs = []
    specials = [0.0, -0.0, 1.0, 2.0, 0.1, 1e-4, 9.99e-5, 1e-5, 1.5e-7, 1e16, 9999999999999998.0, 1.2345678e20, 123456.789,
                math.inf, -math.inf, math.nan, 5e-324, 1.7976931348623157e308, 100.0, 0.30000000000000004]
    for x in specials + [rng.uniform(-1, 1) * 10 ** rng.randint(-12, 20) for _ in range(3000)]:
        reprs.append({"x": bits(x), "repr": repr(float(x)), "np_str": str(np.float64(x))})

    # CPython's random: seeding from an md5 hex string, sample() in both algorithms, choice().
    py_random = []
    for sample_id in ["HG003", "HG002", "sample-1", ""]:
        seed = hashlib.md5(sample_id.encode()).hexdigest()
        r = random.Random()
        r.seed(seed)
        draws = {"seed": seed, "getrandbits32": [r.getrandbits(32) for _ in range(5)]}
        for n, k in [(10, 3), (100, 10), (1000, 500), (3_000_000, 3_000), (2_900_000, 10), (50, 50), (70, 8)]:
            r.seed(seed)
            draws[f"sample_{n}_{k}"] = r.sample(range(n), k)
        r.seed(seed)
        chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        draws["choice"] = "".join(r.choice(chars) for _ in range(64))
        py_random.append(draws)

    np_random = []
    for seed in [0, 7, 12345, 2**32 - 1]:
        np.random.seed(seed)
        np_random.append({"seed": seed, "uniform": [bits(v) for v in np.random.uniform(size=20)]})

    sf = []
    for x in [0.0, 0.1, 0.5, 0.7071, 0.70710678, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 20.0, 37.0, 38.5, -0.3, -1.0, -5.0, -40.0]:
        sf.append({"x": bits(x), "r": bits(scipy.stats.norm.sf(x))})
    for _ in range(2000):
        x = rng.uniform(-12, 40)
        sf.append({"x": bits(x), "r": bits(scipy.stats.norm.sf(x))})

    out = {"arrays": [[bits(v) for v in a] for a in data], "cases": cases, "round": rounds, "repr": reprs, "py_random": py_random, "np_random": np_random, "norm_sf": sf}
    with gzip.open(sys.argv[1], "wt", compresslevel=9) as f:
        json.dump(out, f)


if __name__ == "__main__":
    main()

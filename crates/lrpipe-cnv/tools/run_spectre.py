#!/usr/bin/env python3
"""Runs Spectre's CNVCaller with a seeded numpy RNG and dumps its state after each stage.

Usage: run_spectre.py DUMP_DIR -- <spectre CNVCaller arguments>

Needs the reference environment (tools/spectre-env.txt) and Spectre's source in
$SPECTRE_SRC (tools/fetch_spectre.sh). Spectre's dip test draws from numpy's global RNG
without seeding it; this runner seeds it with $SPECTRE_NUMPY_SEED (default 0) so runs
repeat, and lrpipe-cnv replays the same stream.

DUMP_DIR/NN_<stage>.json holds the analysis objects' fields after that stage. Float
arrays are summarized as {"len", "nan", "sha256"} over their float64 little-endian bytes
with NaNs normalized, so the port can compare whole arrays without storing them.
"""

import hashlib
import json
import math
import os
import sys

import numpy as np
import pandas as pd

SKIP_FIELDS = {"logger", "spectre_args", "args", "genome_info", "metadata", "mosdepth_data"}


def finite(value):
    """A float for JSON: NaN as null, infinities as "inf" / "-inf" (JSON has no non-finite numbers)."""
    v = float(value)
    if math.isnan(v):
        return None
    if math.isinf(v):
        return "inf" if v > 0 else "-inf"
    return v


def array_summary(values):
    a = np.asarray(values, dtype=np.float64).copy()
    nan = np.isnan(a)
    a[nan] = np.nan
    return {
        "len": int(a.size),
        "nan": int(nan.sum()),
        "sha256": hashlib.sha256(a.astype("<f8").tobytes()).hexdigest(),
        "head": [finite(v) for v in a[:5]],
    }


def plain(value, seen, depth=0):
    """A JSON-ready copy of `value`; objects seen before become {"@ref": type#n}."""
    if depth > 8:
        return "<deep>"
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return finite(value)
    if isinstance(value, np.ndarray):
        if value.dtype.kind in "fiub":
            return {"array": array_summary(value)}
        return [plain(v, seen, depth + 1) for v in value.tolist()]
    if isinstance(value, pd.Series):
        if value.dtype.kind in "fiub":
            return {"series": array_summary(value.to_numpy(dtype=np.float64, na_value=np.nan))}
        return {"series": [plain(v, seen, depth + 1) for v in value.tolist()]}
    if isinstance(value, pd.DataFrame):
        return {"frame": {str(c): plain(value[c], seen, depth + 1) for c in value.columns}, "rows": len(value)}
    if isinstance(value, dict):
        return {str(k): plain(v, seen, depth + 1) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        items = sorted(value, key=repr) if isinstance(value, set) else value
        if items and all(isinstance(v, (float, np.floating, int, np.integer)) and not isinstance(v, bool) for v in items):
            return {"array": array_summary(items)}
        return [plain(v, seen, depth + 1) for v in items]
    if hasattr(value, "__dict__"):
        key = id(value)
        if key in seen:
            return {"@ref": seen[key]}
        seen[key] = f"{type(value).__name__}#{len(seen)}"
        fields = {k: plain(v, seen, depth + 1) for k, v in vars(value).items() if k not in SKIP_FIELDS}
        return {"@type": type(value).__name__, "@id": seen[key], **fields}
    return repr(value)


class Dumper:
    def __init__(self, directory):
        self.directory = directory
        self.n = 0
        os.makedirs(directory, exist_ok=True)

    def dump(self, stage, obj):
        self.n += 1
        path = os.path.join(self.directory, f"{self.n:02d}_{stage}.json")
        with open(path, "w") as f:
            json.dump(plain(obj, {}), f, indent=1, sort_keys=True)


def wrap(cls, method, dumper, stage):
    original = getattr(cls, method)

    def wrapped(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        dumper.dump(stage, self)
        return result

    setattr(cls, method, wrapped)


def main():
    if "--" not in sys.argv:
        sys.exit(__doc__)
    split = sys.argv.index("--")
    dump_dir, spectre_args = sys.argv[1], sys.argv[split + 1 :]
    sys.path.insert(0, os.environ["SPECTRE_SRC"])
    np.random.seed(int(os.environ.get("SPECTRE_NUMPY_SEED", "0")))

    from spectre.analysis import analysis, snv_analysis
    from spectre import spectre as spectre_main

    dumper = Dumper(dump_dir)
    a = analysis.CNVAnalysis
    wrap(a, "data_normalization", dumper, "normalization")
    wrap(a, "call_cnv_coverage", dumper, "raw_calls")
    wrap(a, "refine_cnv_calls", dumper, "refined")
    wrap(a, "cnv_loh_candidate_merge", dumper, "final")
    original_metrics = a.get_cnv_metrics

    def metrics(self, refined_cnvs=False):
        result = original_metrics(self, refined_cnvs)
        dumper.dump("metrics_refined" if refined_cnvs else "metrics_raw", self)
        return result

    a.get_cnv_metrics = metrics
    s = snv_analysis.SNVAnalysis
    wrap(s, "snv_copy_number_state", dumper, "snv_cn_neutral")
    wrap(s, "loh", dumper, "loh")

    sys.argv = ["spectre"] + spectre_args
    spectre_main.main()


if __name__ == "__main__":
    main()

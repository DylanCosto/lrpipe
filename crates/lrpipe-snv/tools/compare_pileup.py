"""Compare `lrpipe-snv pileup` output with `oracle_pileup.py` output, bit for bit.

usage: compare_pileup.py ORACLE_PREFIX RUST_PREFIX
Exits non-zero on any difference and prints where the first ones are.
"""
import sys

import numpy as np


def main():
    want, got = sys.argv[1], sys.argv[2]
    ok = True

    w_info = open(want + ".info.tsv").read().splitlines()
    g_info = open(got + ".info.tsv").read().splitlines()
    if w_info != g_info:
        ok = False
        print(f"info: {len(w_info)} vs {len(g_info)} lines")
        w_set, g_set = set(w_info), set(g_info)
        for line in [l for l in w_info if l not in g_set][:5]:
            print(f"  only oracle: {line}")
        for line in [l for l in g_info if l not in w_set][:5]:
            print(f"  only rust:   {line}")
        for i, (a, b) in enumerate(zip(w_info, g_info)):
            if a != b:
                print(f"  first differing line {i}:\n    oracle {a}\n    rust   {b}")
                break

    for name in ("tensors", "inputs", "probs"):
        try:
            w = np.load(f"{want}.{name}.npy")
        except FileNotFoundError:
            continue
        g = np.load(f"{got}.{name}.npy")
        if w.shape != g.shape or w.dtype != g.dtype:
            ok = False
            print(f"{name}: shape/dtype {w.shape} {w.dtype} vs {g.shape} {g.dtype}")
            continue
        same = w.view(np.uint8).tobytes() == g.view(np.uint8).tobytes()
        if same:
            print(f"{name}: identical {w.shape} {w.dtype}")
            continue
        ok = False
        diff = np.argwhere(w != g)
        print(f"{name}: {len(np.unique(diff[:, 0]))} of {len(w)} candidates differ")
        if name == "probs":
            print(f"  max abs diff {np.abs(w.astype(np.float64) - g).max():.3e}")
        for idx in diff[:5]:
            idx = tuple(idx)
            print(f"  at {idx}: oracle {w[idx]} rust {g[idx]}  ({w_info[idx[0]] if idx[0] < len(w_info) else ''})")
    print("PARITY OK" if ok else "PARITY FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

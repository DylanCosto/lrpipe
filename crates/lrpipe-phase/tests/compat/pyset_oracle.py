"""Executes a set-operation script with real CPython sets (see examples/compat_replay.rs)."""
import sys
sets = {}
for line in sys.stdin:
    f = line.split()
    if not f:
        continue
    c = f[0]
    if c == "new": sets[f[1]] = set()
    elif c == "range": sets[f[1]] = set(range(int(f[2])))
    elif c == "add": sets[f[1]].add(int(f[2]))
    elif c == "discard": sets[f[1]].discard(int(f[2]))
    elif c == "isub":
        a = sets[f[1]]; a -= sets[f[2]]
    elif c == "copy": sets[f[1]] = set(sets[f[2]])
    elif c == "diff": sets[f[1]] = sets[f[2]].difference(sets[f[3]])
    elif c == "updlist": sets[f[1]].update([int(x) for x in f[2:]])
    elif c == "merge": sets[f[1]].update(sets[f[2]])
    elif c == "print": print(" ".join(str(x) for x in sets[f[1]]))
    else: raise ValueError(c)

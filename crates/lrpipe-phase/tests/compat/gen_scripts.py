"""Generate random operation scripts exercising set patterns used by readselect.pyx."""
import random, sys
seed = int(sys.argv[1]); kind = sys.argv[2]
rng = random.Random(seed)
out = []
if kind == "py":
    names = ["A", "B", "C", "D"]
    n = rng.choice([5, 50, 300, 2000, 20000])
    out.append(f"range A {n}")
    for nm in names[1:]:
        out.append(f"new {nm}")
    for _ in range(rng.randint(20, 300)):
        op = rng.random()
        a, b, c = rng.sample(names, 3)
        hi = n * rng.choice([1, 2, 10])
        if op < 0.25: out.append(f"add {a} {rng.randrange(hi)}")
        elif op < 0.4: out.append(f"discard {a} {rng.randrange(hi)}")
        elif op < 0.5: out.append(f"isub {a} {b}")
        elif op < 0.6: out.append(f"copy {a} {b}")
        elif op < 0.7: out.append(f"diff {a} {b} {c}")
        elif op < 0.85: out.append(f"updlist {a} " + " ".join(str(rng.randrange(hi)) for _ in range(rng.randint(0, 60))))
        elif op < 0.9: out.append(f"merge {a} {b}")
        else: out.append(f"print {a}")
    for nm in names:
        out.append(f"print {nm}")
else:
    out.append("N")
    for _ in range(rng.randint(5, 40)):
        base = rng.randrange(10**8)
        for _ in range(rng.choice([1, 3, 12, 13, 14, 40, 200, 3000])):
            out.append(f"I {base + rng.randrange(10**6) - 5000}")
        out.append("P")
        if rng.random() < 0.7:
            out.append("C")
        if rng.random() < 0.1:
            out.append("N")
    for _ in range(50):
        s = "".join(rng.choice("ACGT_/:0123456789abcdefXYZ") for _ in range(rng.randint(1, 60)))
        out.append(f"H {s}")
print("\n".join(out))

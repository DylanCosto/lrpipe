"""Oracle for core::switchflip: random polyploid phasing pairs run through
WhatsHap's compiled SwitchFlipCalculator.

Usage: python switchflip_oracle.py SEED COUNT > expected.txt
Each problem is printed as a line "ploidy switch_cost flip_cost P0 P1" (P0
and P1: haplotypes as digit strings, comma-separated), followed by the
result line. examples/switchflip_replay.rs prints the same for the problem
lines, for diffing.
"""

import random
import sys

from whatshap.polyphase.solver import SwitchFlipCalculator


def result_line(r):
    switches, flips, sw_cols, flip_cols, perm_cols = r
    return "{!r} {!r} {} {} {}".format(
        switches,
        flips,
        ",".join(map(str, sw_cols)),
        ";".join(" ".join(map(str, f)) for f in flip_cols),
        ";".join(" ".join(map(str, p)) for p in perm_cols),
    )


def main():
    seed, count = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    for _ in range(count):
        ploidy = rng.choice([2, 3, 3, 4, 4, 5, 6])
        n = rng.choice([1, 2, 3, 5, 10, 30])
        alleles = rng.choice([[0, 1], [0, 1], [0, 1, 2]])
        p0 = [[rng.choice(alleles) for _ in range(n)] for _ in range(ploidy)]
        # p1: p0 with haplotypes permuted in a few segments and some alleles changed
        order = list(range(ploidy))
        p1 = [[0] * n for _ in range(ploidy)]
        for i in range(n):
            if rng.random() < 0.2:
                rng.shuffle(order)
            for k in range(ploidy):
                p1[k][i] = p0[order[k]][i]
            if rng.random() < 0.15:
                p1[rng.randrange(ploidy)][i] = rng.choice(alleles)
        h0 = ["".join(map(str, h)) for h in p0]
        h1 = ["".join(map(str, h)) for h in p1]
        costs = rng.choice([(1, 1), (1, 1), (1, 2 * n * ploidy + 1), (2, 1)])
        print(ploidy, costs[0], costs[1], ",".join(h0), ",".join(h1))
        calc = SwitchFlipCalculator(ploidy, costs[0], costs[1])
        print(result_line(calc.compute_switch_flips_poly(h0, h1)))


if __name__ == "__main__":
    main()

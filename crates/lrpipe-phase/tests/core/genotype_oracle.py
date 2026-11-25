"""Oracle for core::genotypedp: random read sets and pedigrees through
WhatsHap's compiled GenotypeDPTable and compute_genotypes.

Usage: python genotype_oracle.py SEED COUNT > expected.txt
For each problem a block "P" line (problem) follows: positions, reads,
individuals with prior likelihoods, trios and recombination costs, then
the likelihoods WhatsHap computes ("L" lines, float.hex) and the prior
genotypes ("G" lines). examples/genotype_replay.rs prints the same for
the problem lines, for diffing.
"""

import random
import sys

from whatshap.core import (
    Genotype,
    GenotypeDPTable,
    NumericSampleIds,
    Pedigree,
    PhredGenotypeLikelihoods,
    Read,
    ReadSet,
    compute_genotypes,
)


def main():
    seed, count = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    for _ in range(count):
        n_trios = rng.choice([0, 0, 0, 1, 1, 2])
        n_ind = 1 if n_trios == 0 else (3 if n_trios == 1 else 4)
        n_pos = rng.randint(2, 40)
        positions = sorted(rng.sample(range(1, 10000), n_pos))
        haps = [[[rng.randint(0, 1) for _ in positions] for _ in range(2)] for _ in range(n_ind)]
        reads = []
        # coverage stays low enough for the exponential DP (as after read selection)
        for r in range(rng.randint(0, 4 * n_pos // (1 + n_trios))):
            ind = rng.randrange(n_ind)
            start = rng.randrange(n_pos)
            idxs = [i for i in range(start, min(n_pos, start + rng.randint(2, 6))) if rng.random() < 0.9]
            if len(idxs) < 2:
                continue
            h = rng.randint(0, 1)
            variants = []
            for i in idxs:
                allele = haps[ind][h][i] if rng.random() > 0.08 else 1 - haps[ind][h][i]
                quality = rng.choice([0, 3, 10, 20, 30, 30, 30, 40, 60, rng.randint(0, 300)])
                variants.append((positions[i], allele, quality))
            reads.append((f"r{r}", rng.choice([20, 60]), ind, variants))
        priors = []
        for ind in range(n_ind):
            p = []
            for _ in positions:
                if rng.random() < 0.3:
                    p.append((1 / 3, 1 / 3, 1 / 3))
                else:
                    a, b, c = (rng.random() + 1e-6 for _ in range(3))
                    s = a + b + c
                    p.append((a / s, b / s, c / s))
            priors.append(p)
        trios = []
        if n_trios >= 1:
            trios.append((0, 1, 2))
        if n_trios == 2:
            trios.append((0, 1, 3))
        recomb = [rng.choice([0, 1, 10, 20, 50, 80, rng.randint(0, 1200)]) for _ in positions]

        print("P", n_ind, " ".join(map(str, positions)))
        for name, mapq, ind, variants in reads:
            print("R", name, mapq, ind, " ".join(f"{p}:{a}:{q}" for p, a, q in variants))
        for ind in range(n_ind):
            print("I", ind, " ".join(",".join(x.hex() for x in t) for t in priors[ind]))
        for f, m, c in trios:
            print("T", f, m, c)
        print("C", " ".join(map(str, recomb)))

        ids = NumericSampleIds()
        pedigree = Pedigree(ids)
        for ind in range(n_ind):
            pedigree.add_individual(
                f"s{ind}", [Genotype([]) for _ in positions], [PhredGenotypeLikelihoods(list(t)) for t in priors[ind]]
            )
        for f, m, c in trios:
            pedigree.add_relationship(father_id=f"s{f}", mother_id=f"s{m}", child_id=f"s{c}")
        rs = ReadSet()
        for name, mapq, ind, variants in reads:
            read = Read(name, mapq, 0, ids[f"s{ind}"])
            for p, a, q in variants:
                read.add_variant(p, a, q)
            rs.add(read)
        rs.sort()
        accessible = sorted(rs.get_positions())
        if accessible:
            table = GenotypeDPTable(ids, rs, [recomb[positions.index(p)] for p in accessible], pedigree, accessible)
            for ind in range(n_ind):
                for k in range(len(accessible)):
                    print("L", ind, k, " ".join(x.hex() for x in table.get_genotype_likelihoods(f"s{ind}", k)))
        genotypes, gls = compute_genotypes(rs, positions)
        for g, l in zip(genotypes, gls):
            print("G", "/".join(map(str, g.as_vector())) or "none", " ".join(x.hex() for x in l))
        print("E")


if __name__ == "__main__":
    main()

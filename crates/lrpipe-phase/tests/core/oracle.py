"""
Differential test oracle for read selection and the PedMEC DP.

Usage: oracle.py SEED INPUT_OUT  > expected_output
Generates a random problem, writes it to INPUT_OUT (replayed by
examples/core_replay.rs) and prints WhatsHap's results to stdout.
"""
import random
import sys

from whatshap.core import Read, ReadSet, PedigreeDPTable, Pedigree, NumericSampleIds, Genotype
from whatshap.readselect import readselection
from whatshap.pedigree import UniformRecombinationCostComputer

seed = int(sys.argv[1])
rng = random.Random(seed)
trio = rng.random() < 0.3
n_samples = 3 if trio else 1
n_pos = rng.randint(5, 400)
positions = sorted(rng.sample(range(1, 200000), n_pos))
n_reads = rng.randint(1, 250)
const_q = rng.random() < 0.5
max_cov = rng.choice([1, 2, 3, 5, 8, 12, 15]) if not trio else rng.choice([1, 2, 3, 4, 5])
haps = [[[rng.randint(0, 1) for _ in positions] for _ in range(2)] for _ in range(n_samples)]
if trio:
    # child inherits one haplotype from each parent (with occasional recombination)
    haps[2] = [list(haps[0][rng.randint(0, 1)]), list(haps[1][rng.randint(0, 1)])]

lines = []
reads = []
names = set()
for r in range(n_reads):
    sample = rng.randrange(n_samples)
    start = rng.randrange(n_pos)
    length = rng.randint(2, 30)
    idxs = [i for i in range(start, min(n_pos, start + length)) if rng.random() < 0.85]
    if len(idxs) < 2:
        continue
    name = "r" + str(rng.randrange(10**6))
    source = rng.choice([0, 0, 0, 1])
    if (name, source) in names:
        continue
    names.add((name, source))
    hap = rng.randint(0, 1)
    read = Read(name, 60, source, sample)
    for i in idxs:
        allele = haps[sample][hap][i]
        if rng.random() < 0.08:
            allele = 1 - allele
        q = 30 if const_q else rng.randint(1, 60)
        read.add_variant(positions[i], allele, q)
    reads.append(read)

rs = ReadSet()
for read in reads:
    rs.add(read)
rs.sort()

def dump(rs):
    out = [f"READSET {len(rs)}"]
    for read in rs:
        out.append(" ".join([read.name, str(read.source_id), str(read.sample_id)] + [f"{v.position}:{v.allele}:{v.quality}" for v in read]))
    return out

lines += dump(rs)
lines.append(f"MAXCOV {max_cov}")
pref = None
if rng.random() < 0.3:
    pref = {1}
    lines.append("PREFERRED 1")
lines.append(f"SAMPLES {n_samples}")
with open(sys.argv[2], "w") as f:
    print("\n".join(lines), file=f)

print("SORTED", " ".join(r.name for r in rs))
if len(rs) == 0:
    print("EMPTY")
    sys.exit(0)

# --- read selection, per sample as in phase.py
selected_sets = []
for s in range(n_samples):
    sub = rs.subset([i for i, r in enumerate(rs) if r.sample_id == s])
    if len(sub) == 0:
        selected_sets.append(sub)
        print(f"SELECTED {s}")
        continue
    sel = readselection(sub, max_cov, pref)
    print(f"SELECTED {s}", " ".join(str(i) for i in sorted(sel)))
    selected_sets.append(sub.subset(sel))

all_reads = ReadSet()
for sub in selected_sets:
    for r in sub:
        all_reads.add(r)
all_reads.sort()
accessible = sorted(all_reads.get_positions())
ids = NumericSampleIds()
ped = Pedigree(ids)
gts = {}
pos_index = {p: i for i, p in enumerate(positions)}
with open(sys.argv[2] + ".gt", "w") as gtf:
    for s in range(n_samples):
        alleles = [sorted([haps[s][0][pos_index[p]], haps[s][1][pos_index[p]]]) for p in accessible]
        print(" ".join(f"{a}/{b}" for a, b in alleles), file=gtf)
        ped.add_individual(str(s), [Genotype(a) for a in alleles])
if trio:
    ped.add_relationship("0", "1", "2")
assert [ids[str(s)] for s in range(n_samples)] == list(range(n_samples))
recomb = UniformRecombinationCostComputer(1.26).compute(accessible)
try:
    dp = PedigreeDPTable(all_reads, recomb, ped, False, accessible)
    superreads, tv = dp.get_super_reads()
except RuntimeError as e:
    print("ERROR", e)
    sys.exit(0)
print("COST", dp.get_optimal_cost())
for s, sr in enumerate(superreads):
    for h in range(2):
        print(f"SUPER {s} {h}", " ".join(f"{v.position}:{v.allele}:{v.quality}" for v in sr[h]))
print("TV", " ".join(map(str, tv)))
print("PART", " ".join(map(str, dp.get_optimal_partitioning())))

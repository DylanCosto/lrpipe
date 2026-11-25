"""
Simulate a diploid genome with heterozygous variants and ONT-like long reads.

Usage: simulate.py OUTDIR [--seed N] [--contigs 3000000,2000000] [--coverage 30] [--samples 1]
Writes OUTDIR/ref.fa, OUTDIR/variants.vcf (unphased GT, plus truth phase in
OUTDIR/truth.vcf) and OUTDIR/reads_<sample>.fq.
"""
import argparse
import os
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("outdir")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--contigs", default="3000000,2000000")
p.add_argument("--coverage", type=float, default=30)
p.add_argument("--samples", type=int, default=1)
p.add_argument("--snv-rate", type=float, default=1 / 800)
p.add_argument("--indel-rate", type=float, default=1 / 5000)
p.add_argument("--chimeric-rate", type=float, default=0.0, help="fraction of reads joined from two loci")
p.add_argument("--junk-rate", type=float, default=0.0, help="fraction of random (unmappable) reads")
p.add_argument("--segdups", type=int, default=0, help="number of 5-20 kb segmental duplications per contig")
p.add_argument("--trio", action="store_true", help="simulate father, mother and child (with crossovers); writes trio.ped and trio.map")
p.add_argument("--crossovers", type=float, default=2.0, help="mean number of crossovers per parent and contig (--trio)")
p.add_argument("--gt-error-rate", type=float, default=0.0, help="fraction of genotype calls written wrongly to variants.vcf")
p.add_argument("--gl", action="store_true", help="write GL (log10 genotype likelihoods) to variants.vcf")
args = p.parse_args()
rng = np.random.default_rng(args.seed)
os.makedirs(args.outdir, exist_ok=True)
B = np.frombuffer(b"ACGT", dtype=np.uint8)
COMP = np.zeros(256, dtype=np.uint8)
for a, b in zip(b"ACGTN", b"TGCAN"):
    COMP[a] = b

def random_seq(n):
    s = B[rng.integers(0, 4, n)]
    # sprinkle homopolymers
    for start in rng.integers(0, n - 20, n // 2000):
        s[start:start + rng.integers(5, 12)] = s[start]
    return s

contig_lengths = [int(x) for x in args.contigs.split(",")]
names = [f"chr{i + 1}" for i in range(len(contig_lengths))]
refs = {n: random_seq(l) for n, l in zip(names, contig_lengths)}
# segmental duplications: near-identical copies (0.5% divergence) of segments
for n in names:
    L = len(refs[n])
    for _ in range(args.segdups):
        k = int(rng.integers(5000, 20000))
        src = int(rng.integers(0, L - k))
        dst = int(rng.integers(0, L - k))
        seg = refs[n][src:src + k].copy()
        mut = rng.random(k) < 0.005
        seg[mut] = B[rng.integers(0, 4, mut.sum())]
        refs[n][dst:dst + k] = seg
with open(f"{args.outdir}/ref.fa", "w") as f:
    for n in names:
        f.write(f">{n}\n")
        s = refs[n].tobytes().decode()
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + "\n")

samples = ["father", "mother", "child"] if args.trio else [f"sample{i + 1}" for i in range(args.samples)]
# crossover positions for the child's inherited haplotypes: (contig, parent) -> sorted positions
crossovers = {}
# variants: (pos, ref, alt, genotypes per sample as (h0, h1) or None for missing)
variants = {}
for n in names:
    L = len(refs[n])
    vs = []
    n_snv = rng.binomial(L, args.snv_rate)
    n_indel = rng.binomial(L, args.indel_rate)
    pos_snv = rng.choice(np.arange(100, L - 100), n_snv, replace=False)
    pos_ind = rng.choice(np.arange(100, L - 100), n_indel, replace=False)
    for pos in sorted(set(pos_snv.tolist()) | set(pos_ind.tolist())):
        r = chr(refs[n][pos])
        if pos in set(pos_ind.tolist()) and pos not in set(pos_snv.tolist()):
            k = int(rng.integers(1, 6))
            if rng.random() < 0.5:  # insertion
                alt = r + "".join(chr(c) for c in B[rng.integers(0, 4, k)])
                ref_allele = r
            else:  # deletion
                ref_allele = refs[n][pos:pos + k + 1].tobytes().decode()
                alt = r
        else:
            ref_allele = r
            alt = chr(B[(np.where(B == ord(r))[0][0] + rng.integers(1, 4)) % 4])
        gts = []
        for s in samples[:2] if args.trio else samples:
            u = rng.random()
            if u < 0.75:
                gts.append((0, 1) if rng.random() < 0.5 else (1, 0))
            elif u < 0.90:
                gts.append((1, 1))
            elif u < 0.97:
                gts.append((0, 0))
            else:
                gts.append(None)  # missing in VCF, but still present as het in truth
        if args.trio:
            # child: one haplotype from each parent, switching at crossovers
            child = []
            for parent in (0, 1):
                g = gts[parent] if gts[parent] is not None else (0, 1)
                key = (n, parent)
                if key not in crossovers:
                    k = rng.poisson(args.crossovers)
                    crossovers[key] = (int(rng.integers(0, 2)), sorted(rng.integers(0, L, k).tolist()))
                start_hap, xs = crossovers[key]
                hap = (start_hap + sum(1 for x in xs if x <= pos)) % 2
                child.append(g[hap])
            gts.append(tuple(child))
        vs.append((pos, ref_allele, alt, gts))
    # drop overlapping variants
    clean = []
    last_end = -1
    for v in vs:
        if v[0] > last_end:
            clean.append(v)
            last_end = v[0] + len(v[1]) - 1
    variants[n] = clean

def gt_str(g, phased):
    if g is None:
        return "./."
    return f"{g[0]}{'|' if phased else '/'}{g[1]}"

def gl_str(g):
    # log10 likelihoods for 0/0, 0/1, 1/1; the called genotype is most likely
    idx = sum(g)
    vals = [-(rng.random() * 4 + 0.5) for _ in range(3)]
    vals[idx] = -(rng.random() * 0.1)
    return ",".join(f"{v:.3f}" for v in vals)

for fname, phased in (("variants.vcf", False), ("truth.vcf", True)):
    with open(f"{args.outdir}/{fname}", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for n, l in zip(names, contig_lengths):
            f.write(f"##contig=<ID={n},length={l}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        use_gl = args.gl and not phased
        if use_gl:
            f.write('##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for n in names:
            for pos, r, a, gts in variants[n]:
                if not phased:
                    calls = []
                    for g in gts:
                        g = None if g is None else tuple(sorted(g))
                        if g is not None and rng.random() < args.gt_error_rate:
                            g = [(0, 0), (0, 1), (1, 1)][int(rng.integers(0, 3))]
                        c = gt_str(g, False) + f":{int(rng.integers(5, 99))}"
                        if use_gl:
                            c += ":" + (gl_str(g) if g is not None else ".")
                        calls.append(c)
                    fmt = "GT:GQ:GL" if use_gl else "GT:GQ"
                else:
                    calls = [gt_str(g, True) + ":50" for g in gts]
                    fmt = "GT:GQ"
                qual = f"{rng.random() * 100:.3f}"
                f.write(f"{n}\t{pos + 1}\t.\t{r}\t{a}\t{qual}\tPASS\t.\t{fmt}\t" + "\t".join(calls) + "\n")

if args.trio:
    with open(f"{args.outdir}/trio.ped", "w") as f:
        f.write("fam father 0 0 1 0\nfam mother 0 0 2 0\nfam child father mother 1 0\n")
    with open(f"{args.outdir}/trio.map", "w") as f:
        # genetic map with varying rates for the first contig
        f.write("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n")
        cum = 0.0
        for p in range(0, contig_lengths[0], 100000):
            f.write(f"{p} 0 {cum:.6f}\n")
            cum += rng.random() * 0.4

def haplotype(n, si, h):
    ref = refs[n]
    parts = []
    last = 0
    for pos, r, a, gts in variants[n]:
        g = gts[si] if gts[si] is not None else (0, 1)
        if g[h] == 1:
            parts.append(ref[last:pos])
            parts.append(np.frombuffer(a.encode(), dtype=np.uint8))
            last = pos + len(r)
    parts.append(ref[last:])
    return np.concatenate(parts)

def add_errors(s):
    n = len(s)
    u = rng.random(n)
    sub = u < 0.02
    ins = (u >= 0.02) & (u < 0.04)
    dele = (u >= 0.04) & (u < 0.065)
    s = s.copy()
    s[sub] = B[(np.searchsorted(B, s[sub]) + rng.integers(1, 4, sub.sum())) % 4]
    counts = np.ones(n, dtype=np.int64)
    counts[dele] = 0
    counts[ins] = 2
    out = np.repeat(s, counts)
    ins_idx = np.cumsum(counts)[ins] - 1
    out[ins_idx] = B[rng.integers(0, 4, len(ins_idx))]
    return out

for si, sample in enumerate(samples):
    with open(f"{args.outdir}/reads_{sample}.fq", "w") as f:
        rid = 0
        for n in names:
            for h in (0, 1):
                hap = haplotype(n, si, h)
                total = 0
                target = args.coverage / 2 * len(hap)
                while total < target:
                    L = int(min(len(hap) - 1, max(500, rng.lognormal(np.log(12000), 0.6))))
                    start = int(rng.integers(0, len(hap) - L))
                    s = add_errors(hap[start:start + L])
                    if rng.random() < 0.5:
                        s = COMP[s][::-1]
                    u = rng.random()
                    if u < args.chimeric_rate:
                        # join with a segment from elsewhere (either strand)
                        L2 = int(rng.integers(2000, 15000))
                        st2 = int(rng.integers(0, len(hap) - L2))
                        s2 = add_errors(hap[st2:st2 + L2])
                        if rng.random() < 0.5:
                            s2 = COMP[s2][::-1]
                        s = np.concatenate([s, s2])
                    elif u < args.chimeric_rate + args.junk_rate:
                        s = B[rng.integers(0, 4, len(s))]
                    q = (rng.integers(3, 40, len(s)) + 33).astype(np.uint8)
                    f.write(f"@{sample}_read{rid}_{n}_h{h}\n{s.tobytes().decode()}\n+\n{q.tobytes().decode()}\n")
                    rid += 1
                    total += L
print("done", {n: len(variants[n]) for n in names})

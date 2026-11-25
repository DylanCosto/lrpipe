"""Random inputs for comparing `whatshap hapcut2vcf` and `whatshap unphase`
with lrpipe-phase.

Usage: python gen.py SEED DIR

Writes DIR/in.vcf (variants with PS/HP/PQ tags and phased/unphased
genotypes, one to three samples, several chromosomes) and DIR/hapcut.txt
(HapCUT 1 or 2 blocks over a subset of the variants, with "-" alleles,
FV marks and occasional malformed lines).
"""

import random
import sys


def main():
    seed, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(seed)
    chroms = [f"chr{i + 1}" for i in range(rng.choice([1, 2, 3]))]
    samples = [f"S{i}" for i in range(rng.choice([1] * 8 + [2, 3]))]
    # haploid calls and GT not first in FORMAT make WhatsHap fail
    haploid = rng.random() < 0.1
    gt_not_first = rng.random() < 0.1
    lines = ["##fileformat=VCFv4.2"]
    if rng.random() < 0.3:
        lines.append("##phasing=partial")
    for c in chroms:
        lines.append(f"##contig=<ID={c},length=100000>")
    lines += [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        '##FORMAT=<ID=PQ,Number=1,Type=Float,Description="Phase quality">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">',
    ]
    if rng.random() < 0.5:
        lines.append('##FORMAT=<ID=HP,Number=.,Type=String,Description="Haplotype">')
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples))
    variants = []
    for c in chroms:
        pos = rng.randint(1, 100)
        for _ in range(rng.choice([0, 3, 20, 60])):
            pos += rng.choice([0, 1, 5, 50, 300]) if rng.random() < 0.05 else rng.randint(1, 300)
            ref = rng.choice("ACGT")
            alt = rng.choice([b for b in "ACGT" if b != ref]) if rng.random() < 0.8 else ref + "T"
            if rng.random() < 0.05:
                alt += "," + rng.choice([b for b in "ACGT" if b != ref])
            calls = []
            for _ in samples:
                a, b = rng.choice([(0, 1), (1, 0), (0, 0), (1, 1), (0, 1)])
                sep = "|" if rng.random() < 0.6 else "/"
                gt = f"{a}{sep}{b}"
                if rng.random() < 0.03:
                    gt = "./."
                elif haploid and rng.random() < 0.05:
                    gt = "1"
                fields = [gt, str(rng.randint(1, 99)) if rng.random() < 0.8 else ".", rng.choice(["30", "."]), "50"]
                calls.append(":".join(fields))
            fmt = "GQ:GT:PS:PQ" if gt_not_first and rng.random() < 0.3 else "GT:PS:PQ:GQ"
            if fmt.startswith("GQ"):
                calls = [":".join([c.split(":")[3]] + c.split(":")[:3]) for c in calls]
            variants.append((c, pos, ref, alt))
            lines.append("\t".join([c, str(pos), ".", ref, alt, "50", "PASS", ".", fmt] + calls))
    with open(f"{out}/in.vcf", "w") as f:
        f.write("\n".join(lines) + "\n")

    # HapCUT output for a subset of the variants, in VCF order
    hapcut2 = rng.random() < 0.5
    blocks = []
    block = []
    chosen = [v for v in variants if rng.random() < 0.8]
    if rng.random() < 0.1 and chosen:
        rng.shuffle(chosen)
    for k, v in enumerate(chosen):
        if block and (block[-1][0] != v[0] or rng.random() < 0.1):
            blocks.append(block)
            block = []
        block.append(v)
    if block:
        blocks.append(block)
    out_lines = []
    vid = 1
    for block in blocks:
        mec = "" if hapcut2 else "MECscore 1.50 "
        out_lines.append(f"BLOCK: offset: {vid} len: {len(block)} phased: {len(block)} SPAN: 100 {mec}fragments 3")
        for c, pos, ref, alt in block:
            h1 = rng.choice(["0", "1", "0", "1", "-", "2"]) if rng.random() < 0.1 else rng.choice(["0", "1"])
            h2 = "1" if h1 == "0" else "0"
            if hapcut2:
                rest = ["0/1", "0", "0.000000", "-0.000000"]
            else:
                rest = ["0/1", rng.choice(["1,0:-0.0,-0.0,-1.3:0.0:0.0", "2,0:-0.1,-0.1,-2.7:0.0:0.0:FV"])]
            out_lines.append("\t".join([str(vid), h1, h2, c, str(pos), ref, alt.split(",")[0]] + rest))
            vid += 1
        out_lines.append("********")
    if rng.random() < 0.05 and out_lines:
        out_lines.insert(rng.randrange(len(out_lines)), "garbage line")
    with open(f"{out}/hapcut.txt", "w") as f:
        f.write("\n".join(out_lines) + ("\n" if out_lines else ""))


if __name__ == "__main__":
    main()

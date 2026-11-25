"""Random sets of phased VCFs for comparing `whatshap compare` with
`lrpipe-phase compare`.

Usage: python gen.py SEED DIR

Writes DIR/in0.vcf, DIR/in1.vcf, ... (two to four phasings of the same
truth, with switch and flip errors, changed, missing and unphased genotypes
and dropped variants) and DIR/cases.txt, one set of extra command-line
arguments per line (the comparison script adds the output files and the
VCFs). Ploidy 2, 3 or 4; multi-ALT sites and PS of different types occur.
"""

import random
import sys


def main():
    seed, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(seed)
    ploidy = rng.choice([2, 2, 2, 3, 4])
    # a multiway comparison of polyploid phasings fails in WhatsHap
    n_files = rng.choice([2, 2, 2, 3, 4]) if ploidy == 2 or rng.random() < 0.05 else 2
    chroms = [f"chr{i + 1}" for i in range(rng.choice([1, 1, 2, 3]))]
    ps_type = rng.choice(["Integer"] * 6 + ["Float", "String"])
    two_samples = rng.random() < 0.1
    same_sample = two_samples or rng.random() < 0.85

    # truth: per chromosome a list of (pos, ref, alts, haplotypes)
    truth = {}
    for c in chroms:
        variants = []
        pos = rng.randint(1, 100)
        for _ in range(rng.choice([1, 2, 5, 20, 60, 150] * 4 + [0])):
            pos += rng.choice([1, 2, 5, 30, 100, 500])
            ref = rng.choice("ACGT")
            kind = rng.random()
            if kind < 0.75:
                alts = [rng.choice([b for b in "ACGT" if b != ref])]
            elif kind < 0.88:
                alts = [ref + rng.choice("ACGT")]
            else:
                alts = rng.sample([b for b in "ACGT" if b != ref], 2)
            n_alleles = len(alts) + 1
            haps = [rng.randrange(n_alleles) for _ in range(ploidy)]
            if rng.random() < 0.3:
                haps[rng.randrange(ploidy)] = 0
            variants.append((pos, ref, alts, haps))
        truth[c] = variants

    for f in range(n_files):
        sample = "S0" if same_sample else f"S{f}"
        samples = [sample, "T"] if two_samples else [sample]
        lines = [
            "##fileformat=VCFv4.2",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            f'##FORMAT=<ID=PS,Number=1,Type={ps_type},Description="Phase set">',
        ]
        file_chroms = [c for c in chroms if rng.random() < 0.9 or c == chroms[0]]
        for c in file_chroms:
            lines.append(f"##contig=<ID={c}>")
        lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples))
        switch_rate = rng.choice([0.0, 0.02, 0.1])
        flip_rate = rng.choice([0.0, 0.02, 0.1])
        for c in file_chroms:
            order = list(range(ploidy))
            block = rng.randint(1, 1000)
            for pos, ref, alts, haps in truth[c]:
                if rng.random() < 0.05:
                    continue
                if rng.random() < 0.08:
                    block = rng.randint(1, 1000)
                if rng.random() < switch_rate:
                    rng.shuffle(order)
                alleles = [haps[order[k]] for k in range(ploidy)]
                if rng.random() < flip_rate:
                    i, j = rng.sample(range(ploidy), 2)
                    alleles[i], alleles[j] = alleles[j], alleles[i]
                if rng.random() < 0.03:
                    alleles[rng.randrange(ploidy)] = rng.randrange(len(alts) + 1)
                gt = "|".join(map(str, alleles))
                r = rng.random()
                if r < 0.03:
                    gt = "/".join(map(str, alleles))
                elif r < 0.05:
                    gt = "/".join(["."] * ploidy)
                elif r < 0.055:
                    gt = "|".join(map(str, alleles[:-1] + [None])).replace("None", ".")
                if ps_type == "Float":
                    ps = rng.choice([str(block), f"{block}.5"])
                elif ps_type == "String":
                    ps = rng.choice([str(block), f"b{block}"])
                else:
                    ps = str(block)
                if rng.random() < 0.02:
                    ps = "."
                calls = [f"{gt}:{ps}"] + (["|".join(["0"] * (ploidy - 1) + ["1"]) + ":7"] if two_samples else [])
                lines.append("\t".join([c, str(pos), ".", ref, ",".join(alts), ".", "PASS", ".", "GT:PS"] + calls))
        with open(f"{out}/in{f}.vcf", "w") as fh:
            fh.write("\n".join(lines) + "\n")

    p = f"-p {ploidy}" if ploidy != 2 or rng.random() < 0.1 else ""
    if rng.random() < 0.03:
        p = f"-p {rng.choice([2, 3, 4])}"
    base = p + (" --sample S0" if two_samples else "") + ("" if same_sample else " --ignore-sample-name")
    cases = [base, base + " --only-snvs"]
    names = ",".join(f"set{i}" for i in range(n_files))
    cases.append(base + f" --names {names}")
    if ploidy == 2 and "-p" not in p:
        cases.append(base + " --switch-error-bed BED --longest-block-tsv LONGEST")
        # WhatsHap orders differing sample names in the multiway TSV by a
        # (randomly seeded) Python set
        if n_files > 2 and same_sample:
            cases.append(base + " --tsv-multiway MULTI")
    with open(f"{out}/cases.txt", "w") as fh:
        fh.write("\n".join(cases) + "\n")
    with open(f"{out}/nfiles", "w") as fh:
        fh.write(f"{n_files}\n")


if __name__ == "__main__":
    main()

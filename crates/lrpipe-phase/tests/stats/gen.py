"""Random phased VCFs for comparing `whatshap stats` with `lrpipe-phase stats`.

Usage: python gen.py SEED DIR

Writes DIR/in.vcf, DIR/chr-lengths.txt and DIR/cases.txt, one set of extra
command-line arguments per line (the comparison script adds the VCF and the
output files). Covers interleaved and nested blocks, PS of type Integer,
Float and String or undefined, missing PS values, HP tags, unphased and
missing genotypes, indels, multi-ALT and duplicate sites, contigs with and
without lengths, and occasionally unsorted records and inconsistent ploidy.
"""

import random
import sys


def main():
    seed, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(seed)
    n_chrom = rng.choice([1, 1, 2, 3, 4])
    chroms = [f"chr{i + 1}" for i in range(n_chrom)]
    if rng.random() < 0.2:
        rng.shuffle(chroms)
    lengths = {c: rng.choice([None, rng.randint(1, 5000), rng.randint(1000, 100000)]) for c in chroms}
    ps_type = rng.choice(["Integer"] * 5 + ["Float", "String", None])
    use_hp = rng.random() < 0.1
    mix_hp = use_hp and rng.random() < 0.3
    samples = [f"S{i}" for i in range(rng.choice([0, 1, 1, 1, 2, 3]) if rng.random() < 0.95 else 0)]
    if not samples and rng.random() < 0.9:
        samples = ["S0"]
    haploid_rate = 0.02 if rng.random() < 0.2 else 0.0
    unsorted = rng.random() < 0.03
    # missing PS values, and records without PS next to String PS values,
    # make WhatsHap's --block-list fail (ids that cannot be sorted)
    missing_ps = rng.random() < 0.1
    ps_always = ps_type == "String" and rng.random() < 0.9

    lines = ["##fileformat=VCFv4.2"]
    for c in chroms:
        if lengths[c] is not None:
            lines.append(f"##contig=<ID={c},length={lengths[c]}>")
        elif rng.random() < 0.7:
            lines.append(f"##contig=<ID={c}>")
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if ps_type is not None:
        lines.append(f'##FORMAT=<ID=PS,Number=1,Type={ps_type},Description="Phase set">')
    if use_hp:
        lines.append('##FORMAT=<ID=HP,Number=.,Type=String,Description="Phasing haplotype">')
    lines.append('##FORMAT=<ID=PQ,Number=1,Type=Float,Description="Phasing quality">')
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                           + (["FORMAT"] + samples if samples else [])))

    bases = "ACGT"
    for c in chroms:
        n = rng.choice([0, 1, 3, 10, 30, 60, 120])
        pos = rng.randint(1, 50)
        block_ids = [rng.randint(1, 3000) for _ in range(rng.randint(1, 6))]
        records = []
        for _ in range(n):
            step = rng.choice([0] + [1] * 3 + [rng.randint(1, 200) for _ in range(6)])
            pos += step
            ref = rng.choice(bases)
            kind = rng.random()
            if kind < 0.65:
                alts = [rng.choice([b for b in bases if b != ref])]
            elif kind < 0.75:
                alts = [ref + "".join(rng.choice(bases) for _ in range(rng.randint(1, 4)))]
            elif kind < 0.82:
                ref = ref + "".join(rng.choice(bases) for _ in range(rng.randint(1, 4)))
                alts = [ref[0]]
            elif kind < 0.9:
                alts = rng.sample([b for b in bases if b != ref], 2)
            elif kind < 0.93:
                alts = ["."]
            elif kind < 0.95:
                alts = [ref]
            else:
                alts = ["<DEL>"]
            records.append([pos, ref, alts])
        if unsorted and len(records) > 2:
            i = rng.randrange(len(records) - 1)
            records[i][0], records[i + 1][0] = records[i + 1][0] + 1, records[i][0]
        for pos, ref, alts in records:
            fields = [c, str(pos), ".", ref, ",".join(alts), ".", "PASS", "."]
            if samples:
                fmt = ["GT"]
                with_ps = (
                    ps_always
                    or ps_type is not None and rng.random() < 0.95
                    or ps_type is None and rng.random() < 0.3
                )
                with_hp = use_hp and rng.random() < 0.9
                if with_ps and (not use_hp or mix_hp):
                    fmt.append("PS")
                if with_hp:
                    fmt.append("HP")
                if rng.random() < 0.2:
                    fmt.append("PQ")
                calls = []
                for _ in samples:
                    calls.append(
                        ":".join(call(rng, f, ps_type, block_ids, haploid_rate, len(alts), missing_ps) for f in fmt)
                    )
                fields += [":".join(fmt)] + calls
            lines.append("\t".join(fields))

    with open(f"{out}/in.vcf", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(f"{out}/chr-lengths.txt", "w") as f:
        for c in chroms:
            if rng.random() < 0.85:
                f.write(f"{c}\t{rng.randint(1, 200000)}\n")
    cases = ["", "--only-snvs", "--chr-lengths CHRLEN"]
    if chroms:
        pick = rng.sample(chroms, rng.randint(1, len(chroms)))
        cases.append(" ".join(f"--chromosome {c}" for c in pick))
        cases.append(f"--chromosome {','.join(pick)},chrX --only-snvs")
    if len(samples) > 1:
        cases.append(f"--sample {samples[-1]} --chr-lengths CHRLEN")
    if rng.random() < 0.2:
        cases.append("--sample nosuchsample")
    with open(f"{out}/cases.txt", "w") as f:
        f.write("\n".join(cases) + "\n")


def call(rng, field, ps_type, block_ids, haploid_rate, n_alts, missing_ps):
    if field == "GT":
        if rng.random() < haploid_rate:
            return rng.choice(["0", "1", "."])
        a = rng.choice(["0", "1", "1", "."] + (["2"] if n_alts > 1 else []))
        b = rng.choice(["0", "1", "1"] + (["2"] if n_alts > 1 else []))
        if rng.random() < 0.03:
            b = "."
        sep = "|" if rng.random() < 0.75 else "/"
        return f"{a}{sep}{b}"
    if field == "PS":
        if missing_ps and rng.random() < 0.05:
            return "."
        b = rng.choice(block_ids)
        if ps_type == "Float":
            return rng.choice([str(b), f"{b}.5", f"{b}.0"])
        if ps_type in ("String", None):
            return rng.choice([str(b), f"b{b}"])
        return str(b)
    if field == "HP":
        if rng.random() < 0.1:
            return "."
        b = rng.choice(block_ids)
        return rng.choice([f"{b}-1,{b}-2", f"{b}-2,{b}-1"])
    if field == "PQ":
        return rng.choice(["30", "12.5", "."])
    raise ValueError(field)


if __name__ == "__main__":
    main()

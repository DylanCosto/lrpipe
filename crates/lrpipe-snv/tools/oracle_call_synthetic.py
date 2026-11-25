"""Synthetic `output_with` cases from Clair3's own code, for `crates/snv-call/tests`.

Draws random candidates (reference base, alt info, float32 probability rows) designed to reach
the rare branches of `output_from`: every outcome type winning, retries after an outcome is
zeroed out (no matching allele in the alt info, proposed indel lengths without reads), exact
ties between outcome probabilities, non-ACGT reference bases and inputs Clair3 rejects. Each case
runs through `output_with`, both with indel-length heads (full-alignment job settings) and
without (pileup job settings), and a sample of every (mode, outcome, retried, tie) bucket is
written as

    mode <TAB> position <TAB> alt_info <TAB> probabilities <TAB> expected

where mode is `fa` (90 probabilities, INFO F) or `p` (24, INFO P), probabilities are
comma-separated shortest float32 reprs and expected is the VCF line without its newline,
`None` when Clair3 writes nothing, or `ERROR` when it raises. Run with the Python of a Clair3
env:

    python tools/oracle_call_synthetic.py --clair3-dir $ENV/bin \
        --out crates/snv-call/tests/data/call_synthetic.tsv
"""
import argparse
import collections
import random
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--clair3-dir", required=True, help="directory holding clair3.py")
    p.add_argument("--out", required=True)
    p.add_argument("--draws", type=int, default=60000)
    p.add_argument("--per-bucket", type=int, default=12)
    p.add_argument("--seed", type=int, default=20261026)
    a = p.parse_args()

    sys.path.insert(0, a.clair3_dir)
    import numpy as np
    from clair3 import CallVariants as cv

    def config(add_indel_length):
        return cv.OutputConfig(
            is_show_reference=True, is_debug=False, is_haploid_precise_mode_enabled=False,
            is_haploid_sensitive_mode_enabled=False, is_output_for_ensemble=False, quality_score_for_pass=2,
            tensor_fn="PIPE", input_probabilities=False, add_indel_length=add_indel_length, gvcf=False,
            pileup=not add_indel_length, enable_long_indel=False, maximum_variant_length_that_need_infer=50,
            keep_iupac_bases=False,
        )

    rng = random.Random(a.seed)
    bases = "ACGT"

    def head(n, peaks, quantized):
        """A probability vector with mass on a few peaks; dyadic values make exact ties likely."""
        v = [0.0] * n
        for i in peaks:
            v[i] += rng.choice([1, 1, 2, 3, 4]) if quantized else rng.random()
        if rng.random() < 0.5:
            for i in rng.sample(range(n), rng.randint(0, min(n, 4))):
                v[i] += rng.choice([0.125, 0.25, 0.5]) if quantized else rng.random() * 0.3
        total = sum(v) or 1.0
        if quantized:
            # Round to multiples of 1/64 so that products of a few values are exact in float32.
            return [round(x / total * 64) / 64 for x in v]
        return [x / total for x in v]

    def alt_info(ref):
        items = []
        for b in rng.sample(bases, rng.randint(0, 4)):
            if b != ref:
                items.append(("X" + b, rng.choice([1, 2, 3, 5, 5, 8, 12, 20])))
        for _ in range(rng.choice([0, 0, 1, 1, 2, 3, 4])):
            length = rng.choice([1, 1, 2, 2, 3, 4, 5, 7, 10, 15, 16, 17, 30, 50, 60])
            items.append(("I" + ref + "".join(rng.choice(bases) for _ in range(length)), rng.choice([1, 2, 3, 3, 6, 9])))
        for _ in range(rng.choice([0, 0, 1, 1, 2, 3, 4])):
            length = rng.choice([1, 1, 2, 2, 3, 4, 5, 7, 10, 15, 16, 17, 30, 50, 51])
            items.append(("D" + "".join(rng.choice(bases) for _ in range(length)), rng.choice([1, 2, 3, 3, 6, 9])))
        rng.shuffle(items)
        if rng.random() < 0.8:
            items.append(("R" + ref, rng.randint(0, 30)))
        depth = sum(c for _, c in items) + rng.choice([0, 0, 1, 5])
        if rng.random() < 0.03:
            depth = 0
        return "{}-{}".format(depth, "".join(f"{k} {c} " for k, c in items))

    def draw():
        ref = rng.choice(bases * 30 + "NRYKMSW" + "X")
        quantized = rng.random() < 0.6
        gt21 = head(21, rng.sample(range(21), rng.choice([1, 1, 2, 3])), quantized)
        genotype = head(3, rng.sample(range(3), rng.choice([1, 2])), quantized)
        zero = cv.VariantLength.index_offset
        lengths = []
        for _ in range(2):
            peaks = rng.sample(range(33), rng.choice([1, 2]))
            if rng.random() < 0.4:
                peaks[0] = zero
            lengths.append(head(33, peaks, quantized))
        row = np.array(gt21 + genotype + lengths[0] + lengths[1], dtype=np.float32)
        position = f"chr{rng.choice(['1', '20', 'Un:alt'])}:{rng.randint(1, 10 ** 8)}:{ref}"
        return position, alt_info(ref.upper() if ref in bases else "A"), row

    def run(add_indel_length, position, alt, row):
        cfg = config(add_indel_length)
        gt21, genotype = row[:21], row[21:24]
        l1, l2 = (row[24:57], row[57:90]) if add_indel_length else (0, 0)
        try:
            line = cv.output_with(position, alt, gt21, genotype, l1, l2, cfg, None)
        except Exception:  # noqa: BLE001 (any exception is Clair3 rejecting the input)
            return "ERROR", ("error", False, False)
        # Classify: the winning outcome, whether it needed retries, whether the first maximum tied.
        ref = position.rsplit(":", 1)[1]
        alt_dict = dict(zip(alt.rstrip().split("-")[1].split(" ")[::2],
                            [int(x) for x in alt.rstrip().split("-")[1].split(" ")[1::2]]))
        lists = cv.possible_outcome_probabilites_from(gt21, genotype, l1, l2, cv.BASE2ACGT[ref], alt_dict,
                                                      add_indel_length=add_indel_length)
        if len(lists) == 1:
            first_max, tie = lists[0], False
        else:
            probs = [lists[0]] + [x for i in (1, 2, 4, 6, 9, 11, 13, 16, 18) for x in lists[i]]
            first_max = max(probs)
            tie = sum(1 for x in probs if x == first_max) > 1
        flags, _, final = cv.output_from(ref, "chr", 1, 0, gt21, genotype, l1, l2, cfg, None, alt_dict)
        names = ["ref", "homo_snp", "hetero_snp", "homo_ins", "hetero_acgt_ins", "hetero_ins_ins", "homo_del",
                 "hetero_acgt_del", "hetero_del_del", "ins_del"]
        kind = "+".join(n for n, f in zip(names, flags) if f)
        return ("None" if line is None else line.rstrip("\n")), (kind, bool(final != first_max), tie)

    buckets = collections.defaultdict(list)
    for _ in range(a.draws):
        position, alt, row = draw()
        for mode, add in (("fa", True), ("p", False)):
            expected, key = run(add, position, alt, row)
            bucket = buckets[(mode,) + key]
            if len(bucket) < a.per_bucket:
                values = row if add else row[:24]
                probs = ",".join(np.format_float_positional(x, unique=True, trim="-") for x in values)
                bucket.append(f"{mode}\t{position}\t{alt}\t{probs}\t{expected}\n")

    with open(a.out, "w") as fh:
        for key in sorted(buckets):
            fh.writelines(buckets[key])
    for key in sorted(buckets):
        print(key, len(buckets[key]), file=sys.stderr)
    print(f"{sum(len(b) for b in buckets.values())} cases in {len(buckets)} buckets", file=sys.stderr)


if __name__ == "__main__":
    main()

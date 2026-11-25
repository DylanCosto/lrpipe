#!/usr/bin/env python3
"""Writes crates/lrpipe-cnv/tests/data/htslib_floats.json: how htslib 1.18 (pysam 0.22) prints
float FORMAT values, which is how Spectre's VCF floats end up after `bcftools sort`.

Usage: htslib_float_vectors.py OUT.json

Each value is written as text into a FORMAT field of type Float and read back as htslib
formats it; the stored input is the text written (htslib parses it to float32).
"""

import json
import os
import random
import sys
import tempfile

import pysam
import pysam.bcftools


def main():
    rng = random.Random(3)
    texts = ["0", "0.0", "-0.0", "1", "60", "0.5", "0.88", "1.02", "3.47", "0.14", "nan", "1e-05", "0.0001",
             "0.00012345678", "999999", "999999.4", "999999.6", "1000000", "123456.7", "12345.678", "1.5e+20",
             "0.33333333", "2.0000001", "9.9999999", "0.099999999", "100000", "-2.5", "7.000001"]
    for _ in range(4000):
        x = rng.uniform(-1, 1) * 10 ** rng.uniform(-7, 8)
        texts.append(repr(round(x, rng.randint(0, 9))))
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "t.vcf")
        with open(path, "w") as f:
            f.write("##fileformat=VCFv4.2\n##contig=<ID=c,length=100000>\n")
            f.write('##FORMAT=<ID=X,Number=1,Type=Float,Description="x">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n")
            for i, t in enumerate(texts):
                f.write(f"c\t{i + 1}\t.\tN\t<DEL>\t.\t.\t.\tX\t{t}\n")
        out = os.path.join(d, "o.vcf")
        pysam.bcftools.sort("-o", out, path, catch_stdout=False)
        printed = [line.rstrip("\n").split("\t")[9] for line in open(out) if not line.startswith("#")]
    with open(sys.argv[1], "w") as f:
        json.dump([[t, p] for t, p in zip(texts, printed)], f)


if __name__ == "__main__":
    main()

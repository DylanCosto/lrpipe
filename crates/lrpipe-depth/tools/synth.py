#!/usr/bin/env python3
"""Writes a small reference and a SAM of edge cases for depth counting.

Usage: synth.py OUT_DIR. Writes OUT_DIR/synth.fa and OUT_DIR/synth.sam (unsorted).
Covers what the real excerpts rarely hit: overlapping proper pairs with split CIGARs, reads of
only soft clips or insertions, unmapped reads placed on a contig, every flag the default
filter excludes, low MAPQ, insert-size limits, read groups, single-operation =/X CIGARs,
reads at contig starts and ends, and a contig without reads.
"""
import random
import re
import sys

out = sys.argv[1]
rng = random.Random(7)
contigs = [("c1", 5000), ("c2", 3000), ("empty", 1000), ("c3", 800)]
seqs = {n: "".join(rng.choice("ACGT") for _ in range(L)) for n, L in contigs}

with open(f"{out}/synth.fa", "w") as f:
    for n, _ in contigs:
        f.write(f">{n}\n")
        s = seqs[n]
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")


def qlen(cigar):
    return sum(int(n) for n, o in re.findall(r"(\d+)([MIDNSHP=XB])", cigar) if o in "MIS=X")


def rlen(cigar):
    return sum(int(n) for n, o in re.findall(r"(\d+)([MIDNSHP=XB])", cigar) if o in "MDN=X")


def read_seq(contig, pos, cigar):
    """Bases matching the reference where the CIGAR aligns, random elsewhere."""
    ref, p, s = seqs[contig], pos, []
    for n, o in re.findall(r"(\d+)([MIDNSHP=XB])", cigar):
        n = int(n)
        if o in "M=":
            s.append(ref[p:p + n].ljust(n, "A"))
            p += n
        elif o == "X":
            s.append("".join("C" if b != "C" else "G" for b in ref[p:p + n].ljust(n, "A")))
            p += n
        elif o in "IS":
            s.append("".join(rng.choice("ACGT") for _ in range(n)))
        elif o in "DN":
            p += n
    return "".join(s)


lines = []


def rec(name, flag, contig, pos, mapq, cigar, mrnm="*", mpos=0, tlen=0, rg=None):
    """pos is 0-based; written 1-based."""
    seq = read_seq(contig, pos, cigar) if contig != "*" and cigar != "*" else "ACGT"
    if cigar != "*" and qlen(cigar) != len(seq):
        raise ValueError(f"{name}: {cigar} {len(seq)}")
    tags = f"\tRG:Z:{rg}" if rg else ""
    lines.append(f"{name}\t{flag}\t{contig}\t{pos + 1}\t{mapq}\t{cigar}\t{mrnm}\t{mpos + 1 if mrnm != '*' else 0}"
                 f"\t{tlen}\t{seq}\t{'I' * len(seq)}{tags}")


def pair(name, contig, p1, c1, p2, c2, rg="rgA", extra=0, mapq=60):
    end = max(p1 + rlen(c1), p2 + rlen(c2))
    t = end - p1
    rec(name, 99 | extra, contig, p1, mapq, c1, "=", p2, t, rg)
    rec(name, 147 | extra, contig, p2, mapq, c2, "=", p1, -t, rg)


# Overlapping proper pairs: single-op, split CIGARs, full containment, identical starts.
pair("ov1", "c1", 100, "100M", 150, "100M")
pair("ov2", "c1", 300, "40M5D30M2I28M", 330, "10S60M3D27M")
pair("ov3", "c1", 500, "150M", 520, "50M")
pair("ov4", "c1", 700, "80M", 700, "80M")
pair("ov5", "c1", 900, "30M100N30M", 920, "60M")
pair("ov6", "c1", 1100, "50=1X49=", 1120, "100M")
pair("adj", "c1", 1300, "50M", 1350, "50M")
pair("far", "c1", 1500, "50M", 2500, "50M")
# Supplementary and duplicate pairs, low MAPQ pair.
pair("sup", "c1", 1700, "60M", 1720, "60M", extra=0x800)
pair("dup", "c1", 1800, "60M", 1820, "60M", extra=0x400)
pair("lowq", "c1", 1900, "60M", 1920, "60M", mapq=5)
# Pairs in another read group and with large inserts.
pair("rgb", "c1", 2000, "60M", 2030, "60M", rg="rgB")
pair("big", "c2", 100, "50M", 1600, "50M")
# Mate on another contig, flagged proper anyway.
rec("xmate", 99, "c1", 2200, 60, "100M", "c2", 50, 0, "rgA")
rec("xmate", 147, "c2", 50, 60, "100M", "c1", 2200, 0, "rgA")
# Fragment-mode cases: read 2 first, and a pair whose read 1 is reverse.
pair("frag", "c2", 400, "70M", 380, "70M")
# Long-read-like singles.
for i in range(40):
    p = rng.randrange(0, 4600)
    ops = [f"{rng.randrange(1, 30)}S"]
    left = 5000 - p
    while left > 50 and len(ops) < 12:
        m = rng.randrange(10, min(200, left - 20))
        ops.append(f"{m}M")
        left -= m
        k = rng.choice("IDN=X")
        n = rng.randrange(1, 15)
        if k in "DN":
            if left - n < 30:
                break
            left -= n
        elif k in "=X":
            if left - n < 30:
                break
            left -= n
        ops.append(f"{n}{k}")
    ops.append("5M")
    # No zero-length operations: htslib writes them to CRAM but cannot decode them.
    ops.append(f"{rng.randrange(1, 30)}S")
    rec(f"lr{i}", rng.choice([0, 16, 0, 256, 2048, 1024, 512, 16]), "c1", p, rng.choice([0, 1, 20, 60, 255]),
        "".join(ops), rg=rng.choice(["rgA", "rgB", None]))
# Single-operation = and X CIGARs, and reads of only clips or insertions.
rec("eq", 0, "c2", 1000, 60, "100=", rg="rgA")
rec("x1", 0, "c2", 1050, 60, "40X", rg="rgA")
rec("clip", 0, "c2", 1200, 60, "30S", rg="rgA")
rec("ins", 0, "c2", 1210, 60, "5S10I5S", rg="rgA")
rec("clip0", 0, "c3", 0, 60, "20S", rg="rgA")
# Unmapped reads placed on a contig.
rec("unm", 4, "c2", 1500, 0, "*", rg="rgA")
rec("unm2", 4, "c3", 400, 0, "*", rg="rgA")
# Reads at contig edges.
rec("start", 0, "c3", 0, 60, "100M", rg="rgA")
rec("end", 16, "c3", 700, 60, "100M", rg="rgA")
rec("qcfail", 512, "c3", 300, 60, "100M", rg="rgA")
rec("second", 256, "c3", 310, 60, "100M", rg="rgA")
# Fully unmapped reads (not placed).
rec("lost", 4, "*", -1, 0, "*")

with open(f"{out}/synth.sam", "w") as f:
    f.write("@HD\tVN:1.6\tSO:unsorted\n")
    for n, L in contigs:
        f.write(f"@SQ\tSN:{n}\tLN:{L}\n")
    f.write("@RG\tID:rgA\tSM:s\n@RG\tID:rgB\tSM:s\n")
    for l in lines:
        f.write(l + "\n")

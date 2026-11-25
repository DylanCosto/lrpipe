"""Synthetic edge cases for full-alignment tensor parity (run with the Clair3 env's Python).

Writes into OUT a small reference and BAM full of awkward alignments (every CIGAR operation,
zero-length operations, long indels, insertions at the read start, N and '=' bases, missing
qualities, duplicate names, flags Clair3 keeps or drops, reads at both contig ends, a pile deep
enough to force downsampling), phased VCFs, candidate BEDs, and OUT/cases.tsv in the format
`crates/snv-fullalign/tests/parity.rs` reads. With --cigar-fixture DIR it also adds cases on
the 360-read CIGAR fixture (`DIR/cigar-cases.bam`, `.fa`) of the earlier Clair3 study.

    python tools/make_fullalign_cases.py --out WORK/synth [--cigar-fixture DIR] \
        [--chr20 BAM REF BED PHASED_VCF]
"""
import argparse
import os
import random
from array import array

import pysam

M, I, D, N, S, H, P, EQ, X = 0, 1, 2, 3, 4, 5, 6, 7, 8


def write_fasta(path, contigs):
    with open(path, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(path)


def random_read(rng, ref, start, name, deep=False):
    """A read starting at `start` with a random CIGAR; returns an AlignedSegment."""
    cigar, seq = [], []
    pos = start
    if rng.random() < 0.3:
        cigar.append((H, rng.randrange(1, 20)))
    if rng.random() < 0.5:
        n = rng.randrange(1, 30)
        cigar.append((S, n))
        seq.extend(rng.choices("ACGT", k=n))
    if rng.random() < 0.15:  # insertion or deletion right at the start
        op = rng.choice([I, D])
        n = rng.randrange(1, 8)
        cigar.append((op, n))
        if op == I:
            seq.extend(rng.choices("ACGT", k=n))
        else:
            pos += n
    target = rng.randrange(60, 900) if not deep else rng.randrange(200, 700)
    while pos - start < target and pos < len(ref) - 1:
        r = rng.random()
        if r < 0.55:
            op = rng.choice([M, M, M, EQ, X])
            n = min(rng.randrange(1, 60), len(ref) - pos)
            bases = list(ref[pos:pos + n].upper())
            for j in range(n):
                if op == X or rng.random() < 0.04:
                    bases[j] = rng.choice([b for b in "ACGT" if b != bases[j]])
                if rng.random() < 0.01:
                    bases[j] = rng.choice("N=")
            seq.extend(bases)
            cigar.append((op, n))
            pos += n
            continue
        elif r < 0.7:
            n = rng.choice([0, 1, 1, 2, 3, 5, 12, 40, 55, 70])
            cigar.append((I, n))
            seq.extend(rng.choices("ACGTN", weights=[5, 5, 5, 5, 1], k=n))
            continue
        elif r < 0.85:
            n = rng.choice([0, 1, 1, 2, 4, 9, 30, 51, 80])
            n = min(n, len(ref) - pos - 1)
            op = D
        elif r < 0.88:
            n = rng.randrange(5, 120)
            n = min(n, len(ref) - pos - 1)
            op = N
        elif r < 0.91:
            cigar.append((P, rng.randrange(0, 3)))
            continue
        else:
            continue
        if op in (D, N):
            pos += n
        cigar.append((op, n))
    if not cigar or cigar[-1][0] not in (M, EQ, X):
        n = min(5, len(ref) - pos)
        if n > 0:
            cigar.append((M, n))
            seq.extend(ref[pos:pos + n].upper())
            pos += n
    if rng.random() < 0.5:
        n = rng.randrange(1, 30)
        cigar.append((S, n))
        seq.extend(rng.choices("ACGT", k=n))
    if rng.random() < 0.2:
        cigar.append((H, rng.randrange(1, 9)))
    # Adjacent operations of one kind are kept on purpose.
    r = pysam.AlignedSegment()
    r.query_name = name
    r.flag = 16 if rng.random() < 0.5 else 0
    r.reference_start = start
    r.mapping_quality = rng.choice([0, 3, 5, 10, 19, 20, 21, 40, 59, 60, 60, 60, 255])
    r.cigartuples = cigar
    r.query_sequence = "".join(seq)
    if rng.random() < 0.05:
        r.query_qualities = None
    else:
        r.query_qualities = array("B", (rng.choice([0, 1, 7, 13, 20, 39, 40, 41, 60, 93]) for _ in seq))
    return r


def edge_inputs(out, rng):
    ref1 = list("".join(rng.choices("ACGT", k=4000)))
    for i in range(1000, 1100):
        ref1[i] = ref1[i].lower()
    for i in range(1500, 1510):
        ref1[i] = "N"
    ref1[1600] = "R"
    ref1[1601] = "y"
    ref1 = "".join(ref1)
    ref2 = "".join(rng.choices("ACGT", k=500))
    fasta = os.path.join(out, "edge.fa")
    write_fasta(fasta, [("ctg", ref1), ("ctg2", ref2)])
    refs = {"ctg": ref1, "ctg2": ref2}

    records = []
    k = 0
    for tid, (name, ref) in enumerate([("ctg", ref1), ("ctg2", ref2)]):
        starts = []
        if name == "ctg":
            starts += [0] * 5 + list(range(0, 40, 3))
            starts += [rng.randrange(150, 260) for _ in range(260)]  # deep pile
            starts += [rng.randrange(0, 3900) for _ in range(500)]
            starts += [len(ref) - rng.randrange(5, 80) for _ in range(30)]
        else:
            starts += [0, 0, 1, 2, 5] + [rng.randrange(0, 480) for _ in range(120)]
        for s in starts:
            r = random_read(rng, ref, s, f"r{k:05}", deep=150 <= s < 260 and name == "ctg")
            r.reference_id = tid
            u = rng.random()
            if u < 0.03:
                r.flag |= 256
            elif u < 0.06:
                r.flag |= 2048
            elif u < 0.08:
                r.flag |= 1 | 8
            elif u < 0.10:
                r.flag |= 512  # kept: not in 2316
            elif u < 0.12:
                r.flag |= 1024  # kept
            elif u < 0.14 and k > 10:
                r.query_name = f"r{k - 7:05}"  # duplicate name, skipped if the first was kept
            records.append(r)
            k += 1
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "ctg", "LN": len(ref1)}, {"SN": "ctg2", "LN": len(ref2)}]}
    unsorted = os.path.join(out, "edge.unsorted.bam")
    with pysam.AlignmentFile(unsorted, "wb", header=header) as fh:
        for r in records:
            fh.write(r)
    bam = os.path.join(out, "edge.bam")
    pysam.sort("-o", bam, unsorted)
    pysam.index(bam)
    os.remove(unsorted)

    # Phased SNPs, with unphased and odd genotypes and several phase sets.
    rows = []
    for name, ref in refs.items():
        for pos in range(0, len(ref), 7):
            if rng.random() < 0.4:
                continue
            base = ref[pos].upper()
            alt = rng.choice([b for b in "ACGT" if b != base])
            gt = rng.choice(["0|1", "1|0", "0|1", "1|0", "0/1", "1|1"])
            ps = 1 + pos // 600 * 600
            rows.append((name, pos + 1, ref[pos], alt, gt, ps))
    vcf = os.path.join(out, "edge.phased.vcf")
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for name, ref in refs.items():
            fh.write(f"##contig=<ID={name},length={len(ref)}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for name, pos, r, a, gt, ps in rows:
            fh.write(f"{name}\t{pos}\t.\t{r}\t{a}\t30\tPASS\t.\tGT:GQ:PS\t{gt}:30:{ps}\n")
    pysam.tabix_compress(vcf, vcf + ".gz", force=True)
    return bam, fasta, vcf + ".gz", refs


def bed_rows(contig, centers):
    """Rows whose centre (as CreateTensorFullAlignment computes it) is each 0-based centre."""
    for c in centers:
        if c - 16 > 0:
            yield f"{contig}\t{c - 16}\t{c + 18}\n"
        else:
            yield f"{contig}\t0\t{c + 19}\n"


def write_bed(path, contig, centers, snp_rows=()):
    with open(path, "w") as fh:
        for row in bed_rows(contig, sorted(centers)):
            fh.write(row)
        for row in snp_rows:
            fh.write(row)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", required=True)
    p.add_argument("--cigar-fixture")
    p.add_argument("--chr20", nargs=4, metavar=("BAM", "REF", "BED", "PHASED_VCF"))
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(20261026)
    cases = []

    def case(name, bam, ref, contig, bed, vcf="-", mq=5, platform="ont", long_indel=0):
        cases.append("\t".join(map(str, [name, bam, ref, contig, bed, vcf, mq, platform, long_indel])))

    bam, fasta, vcf, refs = edge_inputs(a.out, rng)
    L1, L2 = len(refs["ctg"]), len(refs["ctg2"])
    dense = write_bed(os.path.join(a.out, "edge.dense.bed"), "ctg", range(L1))
    sparse = write_bed(os.path.join(a.out, "edge.sparse.bed"), "ctg", rng.sample(range(L1), 300),
                       snp_rows=["ctg\t3990\t4005\tA-C-1-1\n", "ctg2\t1\t30\n"])
    ends = write_bed(os.path.join(a.out, "edge.ctg2.bed"), "ctg2", range(L2))
    het_only = os.path.join(a.out, "edge.hetonly.bed")
    with open(het_only, "w") as fh:
        fh.write("ctg\t100\t121\tA-C-1-1\nctg\t500\t521\tG-T-2-1\n")
    for mq in (0, 5, 20):
        case(f"edge_dense_mq{mq}", bam, fasta, "ctg", dense, vcf, mq)
    case("edge_dense_nophase", bam, fasta, "ctg", dense, "-", 5)
    case("edge_dense_missingvcf", bam, fasta, "ctg", dense, os.path.join(a.out, "missing.vcf.gz"), 5)
    case("edge_dense_hifi", bam, fasta, "ctg", dense, vcf, 0, "hifi")
    case("edge_dense_long", bam, fasta, "ctg", dense, vcf, 0, "ont", 1)
    case("edge_sparse", bam, fasta, "ctg", sparse, vcf, 5)
    case("edge_ctg2", bam, fasta, "ctg2", ends, vcf, 0)
    case("edge_ctg2_hifi_long", bam, fasta, "ctg2", ends, vcf, 0, "hifi", 1)
    case("edge_hetonly", bam, fasta, "ctg", het_only, vcf, 5)

    if a.cigar_fixture:
        cbam = os.path.join(a.cigar_fixture, "cigar-cases.bam")
        cfa = os.path.join(a.cigar_fixture, "cigar-cases.fa")
        seq = pysam.FastaFile(cfa).fetch("test")
        # The fixture's phased SNPs (verify.py): position, ref, C or A, genotype, phase set.
        cvcf = os.path.join(a.out, "cigar.phased.vcf")
        with open(cvcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            for pos in range(520, 8600, 47):
                alt = "C" if seq[pos] != "C" else "A"
                gt = "0|1" if pos % 2 else "1|0"
                fh.write(f"test\t{pos + 1}\t.\t{seq[pos]}\t{alt}\t30\tPASS\t.\tGT:PS\t{gt}:{500 + pos // 3500 * 3500}\n")
        pysam.tabix_compress(cvcf, cvcf + ".gz", force=True)
        cvcf += ".gz"
        edges = set()
        with pysam.AlignmentFile(cbam) as fh:
            for r in fh.fetch("test"):
                pos = r.reference_start
                for op, n in r.cigartuples:
                    edges.add(pos)
                    if op in (M, D, N, EQ, X):
                        pos += n
                edges.add(pos)
        crng = random.Random(907)
        sets = {
            "no_reads": [12000],
            "sparse": list(range(450, 8800, 137)),
            "overlapping_windows": list(range(580, 1100, 7)),
            "cigar_boundaries": sorted({p + d for p in edges for d in (-1, 0, 1)}),
            "long_gap": [600, 4100, 7650, 12000],
        }
        for i in range(3):
            sets[f"random_{i}"] = sorted(crng.sample(range(450, 8900), 100 + i * 73))
        for name, centers in sets.items():
            bed = write_bed(os.path.join(a.out, f"cigar.{name}.bed"), "test", centers)
            case(f"cigar_{name}", cbam, cfa, "test", bed, cvcf, 5)
            case(f"cigar_{name}_mq0_hifi", cbam, cfa, "test", bed, cvcf, 0, "hifi")
        case("cigar_cigar_boundaries_nophase", cbam, cfa, "test",
             os.path.join(a.out, "cigar.cigar_boundaries.bed"), "-", 0)

    if a.chr20:
        bam20, ref20, bed20, vcf20 = a.chr20
        case("chr20_nophase", bam20, ref20, "chr20", bed20, "-", 5)
        case("chr20_hifi_long", bam20, ref20, "chr20", bed20, vcf20, 5, "hifi", 1)

    with open(os.path.join(a.out, "cases.tsv"), "w") as fh:
        fh.write("\n".join(cases) + "\n")
    print(f"{len(cases)} cases in {a.out}/cases.tsv")


if __name__ == "__main__":
    main()

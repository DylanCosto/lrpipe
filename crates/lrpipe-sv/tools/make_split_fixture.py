"""Makes tests/data/split*.bam from sample1.bam: reads with split and secondary alignments, for
the phasing tests (tests/phasing.rs). Run by tools/make_split_fixture.sh.

The reads around chr17:46,237,501 (a 724 bp deletion on haplotype 2) are kept. Every third
read carrying the deletion is split at it into a primary alignment (the bases before it) and a
supplementary one (the bases after it), each with the other's SA tag, as an aligner writes a
split read; three reads also get a secondary alignment. Then:

- split.bam: tagged as `whatshap haplotag` (`--tag-supplementary skip`) tags them: each primary
  alignment keeps sample1.bam's HP and PS, supplementary and secondary alignments have none;
- split.haplotags.tsv: the haplotag list of those tags (primary alignments only);
- split.foreign.bam: the same alignments tagged by "another pipeline": every alignment HP 1 and
  another phase set, supplementary and secondary ones included.

usage: make_split_fixture.py SAMPLE1_BAM OUT_DIR
"""

import sys

import pysam

DEL_POS = 46237501
WINDOW = (46227000, 46248000)
CONTIG = "chr17"


def cigar_string(ops):
    return "".join(f"{n}{'MIDNSHP=X'[op]}" for op, n in ops)


def split(r):
    """The primary and supplementary halves of `r` at its long deletion near DEL_POS, or None."""
    ops = r.cigartuples
    ref, query = r.reference_start, 0
    for i, (op, n) in enumerate(ops):
        if op == 2 and n >= 700 and abs(ref - DEL_POS) < 300:
            qlen = r.query_length
            first = ops[:i] + [(4, qlen - query)]
            second = [(4, query)] + ops[i + 1 :]
            # Clips next to clips merge.
            def merge(cig):
                out = []
                for op, n in cig:
                    if out and out[-1][0] == op == 4:
                        out[-1] = (4, out[-1][1] + n)
                    else:
                        out.append((op, n))
                return out
            return merge(first), ref + n, merge(second), n
        if op in (0, 2, 3, 7, 8):
            ref += n
        if op in (0, 1, 4, 7, 8):
            query += n
    return None


def aligned(ops):
    return sum(n for op, n in ops if op in (0, 7, 8))


def main():
    src, out = sys.argv[1:3]
    bam = pysam.AlignmentFile(src)
    reads = [r for r in bam.fetch(CONTIG, *WINDOW)]
    records = []
    carriers = 0
    secondaries = 0
    for r in reads:
        hp = r.get_tag("HP") if r.has_tag("HP") else None
        ps = r.get_tag("PS") if r.has_tag("PS") else None
        halves = split(r)
        if halves is not None:
            carriers += 1
        if halves is not None and carriers % 3 == 1:
            first, second_start, second, dellen = halves
            nm = r.get_tag("NM") - dellen if r.has_tag("NM") else 0
            a1, a2 = aligned(first), aligned(second)
            nm1 = round(nm * a1 / (a1 + a2))
            nm2 = nm - nm1
            strand = "-" if r.is_reverse else "+"
            p = pysam.AlignedSegment(bam.header)
            p.query_name, p.flag, p.reference_id = r.query_name, r.flag, r.reference_id
            p.reference_start, p.mapping_quality, p.cigartuples = r.reference_start, r.mapping_quality, first
            p.query_sequence = r.query_sequence
            s = pysam.AlignedSegment(bam.header)
            s.query_name, s.flag, s.reference_id = r.query_name, r.flag | 0x800, r.reference_id
            s.reference_start, s.mapping_quality, s.cigartuples = second_start, r.mapping_quality, second
            s.query_sequence = r.query_sequence
            p.set_tag("NM", nm1)
            s.set_tag("NM", nm2)
            p.set_tag("SA", f"{CONTIG},{second_start + 1},{strand},{cigar_string(second)},{r.mapping_quality},{nm2};")
            s.set_tag("SA", f"{CONTIG},{r.reference_start + 1},{strand},{cigar_string(first)},{r.mapping_quality},{nm1};")
            records += [(p, hp, ps, "primary"), (s, hp, ps, "supplementary")]
            continue
        records.append((r, hp, ps, "primary"))
        if halves is None and secondaries < 3 and r.reference_start > WINDOW[0] + 2000:
            secondaries += 1
            sec = pysam.AlignedSegment(bam.header)
            sec.query_name, sec.flag, sec.reference_id = r.query_name, r.flag | 0x100, r.reference_id
            sec.reference_start, sec.mapping_quality, sec.cigartuples = r.reference_start + 500, 0, r.cigartuples
            sec.query_sequence = r.query_sequence
            records.append((sec, hp, ps, "secondary"))

    def write(path, tag):
        tmp = path + ".unsorted.bam"
        with pysam.AlignmentFile(tmp, "wb", header=bam.header) as o:
            for rec, hp, ps, kind in records:
                rec = pysam.AlignedSegment.fromstring(rec.to_string(), bam.header)
                for t in ("HP", "PS", "PC"):
                    rec.set_tag(t, None)
                for t, v in tag(hp, ps, kind):
                    rec.set_tag(t, v, value_type="i")
                o.write(rec)
        pysam.sort("-o", path, tmp)
        pysam.index(path)
        import os

        os.remove(tmp)

    def haplotagged(hp, ps, kind):
        return [("HP", hp), ("PS", ps)] if kind == "primary" and hp is not None else []

    def foreign(hp, ps, kind):
        return [("HP", 1), ("PS", 46000000)]

    write(f"{out}/split.bam", haplotagged)
    write(f"{out}/split.foreign.bam", foreign)
    with open(f"{out}/split.haplotags.tsv", "w") as f:
        f.write("#readname\thaplotype\tphaseset\tchromosome\n")
        for rec, hp, ps, kind in sorted(records, key=lambda x: x[0].reference_start):
            if kind == "primary":
                assigned = f"H{hp}\t{ps}" if hp is not None else "none\tnone"
                f.write(f"{rec.query_name}\t{assigned}\t{CONTIG}\n")
    print(f"{len(records)} alignments, {carriers} carrying the deletion, {secondaries} secondary", file=sys.stderr)


main()

#!/usr/bin/env python3
"""Cuts the test fixture in tests/data from an indexed BAM/CRAM.

usage: make_fixture.py ALIGNMENTS REFERENCE REGION OUTDIR

Writes, for the reads overlapping REGION (primary and supplementary alignments):
  sample1.bam(.bai)      the reads without base qualities, keeping only SA, NM, HP and PS
  sample2.bam            the reads not on haplotype 1 (a second sample for multi-sample calling)
  reference.window.fa    the reference bases the reads cover, named CONTIG:START-END (1-based),
                         with the contig length as LN:LENGTH
The tests rebuild the full-length contig from the window, padded with N.
"""
import hashlib
import sys

import pysam

KEEP_TAGS = {"SA", "NM", "HP", "PS"}
PAD = 1000


def main():
    src, reference, region, out = sys.argv[1:5]
    contig = region.split(":")[0]
    with pysam.AlignmentFile(src, reference_filename=reference) as f:
        header = f.header.to_dict()
        records = [r for r in f.fetch(region=region) if not r.is_secondary]
        contig_len = f.get_reference_length(contig)

    start = max(0, min(r.reference_start for r in records) - PAD)
    end = min(contig_len, max(r.reference_end for r in records) + PAD)

    with pysam.FastaFile(reference) as fa:
        seq = fa.fetch(contig, start, end).upper()
    # Keep the reference sequences (for SA and BND mates) and nothing that names the source. CRAM
    # needs each sequence's MD5; the fixture contig's is that of the N-padded contig.
    padded = hashlib.md5()
    padded.update(b"N" * start + seq.encode() + b"N" * (contig_len - end))
    header = {
        "HD": header["HD"],
        "SQ": [
            {"SN": sq["SN"], "LN": sq["LN"], "M5": padded.hexdigest() if sq["SN"] == contig else sq["M5"]}
            for sq in header["SQ"]
        ],
    }

    with open(f"{out}/reference.window.fa", "w") as f:
        f.write(f">{contig}:{start + 1}-{end} LN:{contig_len}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")

    def write(path, keep):
        with pysam.AlignmentFile(path, "wb", header=header) as w:
            for r in records:
                if not keep(r):
                    continue
                tags = [(k, v, t) for k, v, t in r.get_tags(with_value_type=True) if k in KEEP_TAGS]
                r.query_qualities = None
                r.set_tags(tags)
                w.write(r)
        pysam.index(path)

    write(f"{out}/sample1.bam", lambda r: True)
    write(f"{out}/sample2.bam", lambda r: not (r.has_tag("HP") and r.get_tag("HP") == 1))
    print(f"{len(records)} records; reference window {contig}:{start + 1}-{end}", file=sys.stderr)


if __name__ == "__main__":
    main()

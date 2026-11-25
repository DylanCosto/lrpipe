"""Dump a BAM in a canonical text form including BAM-level tag type codes,
so that outputs can be compared exactly (compression may differ)."""
import sys
import pysam

with pysam.AlignmentFile(sys.argv[1], check_sq=False) as f:
    # header without the @PG line(s) written by whatshap / lrpipe-phase
    for line in str(f.header).rstrip("\n").split("\n"):
        if line.startswith(("@PG\tPN:whatshap", "@PG\tPN:lrpipe-phase")):
            print("@PG(whatshap) <normalized>")
        else:
            print(line)
    for r in f.fetch(until_eof=True):
        tags = "\t".join(f"{t}:{c}:{v}" for t, v, c in r.get_tags(with_value_type=True))
        print(r.to_string().split("\t")[:11], tags, sep="\t")

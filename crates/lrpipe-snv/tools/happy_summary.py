"""Prints hap.py's PASS SNP and INDEL metrics for a few GIAB v3.3 stratifications.

usage: python3 happy_summary.py happy.extended.csv

Reads the overall rows (Subtype `*`) of the whole benchmark (`*`), all difficult regions and
their complement, homopolymers >= 7 bp, low mappability and segmental duplications, and MHC.
"""
import csv
import sys

SUBSETS = ["*", "alldifficultregions", "notinalldifficultregions",
           "AllHomopolymers_ge7bp_imperfectge11bp_slop5", "alllowmapandsegdupregions", "segdups", "MHC"]


def fmt(x):
    return f"{float(x):.4f}" if x else "-"


def main():
    with open(sys.argv[1]) as f:
        rows = {(r["Type"], r["Subset"]): r for r in csv.DictReader(f)
                if r["Filter"] == "PASS" and r["Subtype"] == "*"}
    for kind in ("SNP", "INDEL"):
        for subset in SUBSETS:
            r = rows.get((kind, subset))
            if r:
                print(f'{kind:6} {subset:45} TP {r["TRUTH.TP"]:>9} FN {r["TRUTH.FN"]:>7} '
                      f'FP {r["QUERY.FP"]:>7} R {fmt(r["METRIC.Recall"])} '
                      f'P {fmt(r["METRIC.Precision"])} F1 {fmt(r["METRIC.F1_Score"])}')


if __name__ == "__main__":
    main()

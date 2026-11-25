"""hap.py-style metrics per GIAB stratification from rtg vcfeval's output directories.

Recall counts truth variants (tp-baseline, fn) and precision call variants (tp, fp) whose
position falls in a stratification's regions; F1 combines the two, as hap.py does.

usage: giab_strat.py EVAL_SNPS_DIR EVAL_INDELS_DIR
"""
import bisect
import gzip
import os
import sys

STRAT_DIR = os.path.join(os.environ["GIAB_DIR"], "strat/GRCh38@all")
STRATS = [
    ("all", None),
    ("all difficult regions", "Union/GRCh38_alldifficultregions.bed.gz"),
    ("not in difficult regions", "Union/GRCh38_notinalldifficultregions.bed.gz"),
    ("homopolymers >=7 bp", "LowComplexity/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz"),
    ("low mappability and segdups", "Union/GRCh38_alllowmapandsegdupregions.bed.gz"),
    ("segdups", "SegmentalDuplications/GRCh38_segdups.bed.gz"),
    ("MHC", "OtherDifficult/GRCh38_MHC.bed.gz"),
]


def load_bed(path):
    """Per contig, sorted merged [start, end) intervals."""
    by_contig = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            by_contig.setdefault(c, []).append((int(s), int(e)))
    merged = {}
    for c, ivs in by_contig.items():
        ivs.sort()
        out = []
        for s, e in ivs:
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[c] = ([s for s, _ in out], [e for _, e in out])
    return merged


def inside(bed, contig, pos0):
    starts, ends = bed.get(contig, ((), ()))
    i = bisect.bisect_right(starts, pos0) - 1
    return i >= 0 and pos0 < ends[i]


def positions(path):
    """(contig, 0-based position) of each record of a vcfeval output VCF."""
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, p = line.split("\t", 2)[:2]
            out.append((c, int(p) - 1))
    return out


def metrics(eval_dir, beds):
    sets = {k: positions(os.path.join(eval_dir, f"{k}.vcf.gz")) for k in ("tp-baseline", "fn", "tp", "fp")}
    rows = []
    for name, bed in beds:
        n = {k: sum(1 for c, p in v if bed is None or inside(bed, c, p)) for k, v in sets.items()}
        recall = n["tp-baseline"] / max(1, n["tp-baseline"] + n["fn"])
        precision = n["tp"] / max(1, n["tp"] + n["fp"])
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append((name, n["tp-baseline"], n["fn"], n["fp"], recall, precision, f1))
    return rows


def main():
    beds = [(name, None if rel is None else load_bed(os.path.join(STRAT_DIR, rel))) for name, rel in STRATS]
    print("type\tstratification\tTP\tFN\tFP\trecall\tprecision\tF1")
    for kind, d in (("SNP", sys.argv[1]), ("INDEL", sys.argv[2])):
        for name, tp, fn, fp, r, p, f in metrics(d, beds):
            print(f"{kind}\t{name}\t{tp}\t{fn}\t{fp}\t{r:.5f}\t{p:.5f}\t{f:.5f}")


if __name__ == "__main__":
    main()

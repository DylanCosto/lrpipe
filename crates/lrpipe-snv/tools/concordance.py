"""Compares the PASS calls of two VCFs, and optionally scores each against a truth set.

usage: python3 concordance.py BASE.vcf.gz OTHER.vcf.gz [TRUTH.vcf.gz TRUTH.bed]

A call is (contig, position, REF, ALT, unphased GT) of a PASS record with a non-reference
genotype. Agreement is the share of BASE's calls OTHER makes exactly; Jaccard is shared calls
over the union. Against a truth set (calls inside the BED only) matching is exact on alleles and
genotype, without haplotype-aware normalisation, so F1 is lower than hap.py's but comparable
between runs.
"""
import bisect
import gzip
import sys


def calls(path, regions=None):
    out = set()
    with gzip.open(path, 'rt') if path.endswith('.gz') else open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\t', 10)
            if fields[6] not in ('PASS', '.') or (fields[6] == '.' and regions is None):
                continue
            keys = fields[8].split(':')
            values = fields[9].rstrip('\n').split(':')
            gt = values[keys.index('GT')].replace('|', '/')
            alleles = sorted(gt.split('/'))
            if all(a in ('0', '.') for a in alleles):
                continue
            contig, pos = fields[0], int(fields[1])
            if regions is not None and not inside(regions, contig, pos):
                continue
            out.add((contig, pos, fields[3], fields[4], '/'.join(alleles)))
    return out


def read_bed(path):
    regions = {}
    with open(path) as f:
        for line in f:
            contig, start, end = line.split('\t')[:3]
            regions.setdefault(contig, []).append((int(start), int(end)))
    return {c: (sorted(r), [s for s, _ in sorted(r)]) for c, r in regions.items()}


def inside(regions, contig, pos):
    if contig not in regions:
        return False
    spans, starts = regions[contig]
    i = bisect.bisect_right(starts, pos - 1) - 1
    return i >= 0 and spans[i][0] <= pos - 1 < spans[i][1]


def kind(call):
    ref, alts = call[2], call[3].split(',')
    return 'snp' if len(ref) == 1 and all(len(a) == 1 for a in alts) else 'indel'


def score(name, truth, called):
    for k in ('snp', 'indel'):
        t = {c for c in truth if kind(c) == k}
        p = {c for c in called if kind(c) == k}
        tp = len(t & p)
        recall, precision = tp / max(len(t), 1), tp / max(len(p), 1)
        f1 = 2 * recall * precision / max(recall + precision, 1e-12)
        print(f'{name} {k}: recall {recall:.5f} precision {precision:.5f} F1 {f1:.5f}')


def main():
    base, other = calls(sys.argv[1]), calls(sys.argv[2])
    shared = len(base & other)
    print(f'calls: base {len(base)} other {len(other)} shared {shared} '
          f'only-base {len(base - other)} only-other {len(other - base)}')
    print(f'agreement {shared / max(len(base), 1):.5f} jaccard {shared / max(len(base | other), 1):.5f}')
    if len(sys.argv) > 4:
        regions = read_bed(sys.argv[4])
        truth = calls(sys.argv[3], regions)
        inside_regions = lambda s: {c for c in s if inside(regions, c[0], c[1])}
        score('base ', truth, inside_regions(base))
        score('other', truth, inside_regions(other))


main()

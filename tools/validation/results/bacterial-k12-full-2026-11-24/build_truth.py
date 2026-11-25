"""Project a published haploid assembly onto another reference and make a small-variant truth VCF."""
import argparse
from collections import Counter
from pathlib import Path
import subprocess


MAPQ = 60
MIN_ALIGNMENT = 5_000
BOUNDARY_PADDING = 100
SMALL_INDEL_LIMIT = 49


def assembly_alignments(target, query):
    fields = subprocess.check_output([
        'minimap2', '-t', '4', '-cx', 'asm5', '--eqx', '--cs=short', '--secondary=yes', target, query,
    ], text=True).splitlines()
    for line in fields:
        row = line.split('\t')
        if int(row[11]) < MAPQ or 'tp:A:P' not in row[12:]:
            continue
        if min(int(row[3]) - int(row[2]), int(row[8]) - int(row[7])) < MIN_ALIGNMENT:
            continue
        yield row


def endpoints(row):
    return tuple(row[i] for i in (0, 2, 3)), row[4], tuple(row[i] for i in (5, 7, 8))


def masks(ranges, excluded):
    for beg, end in ranges:
        for lo, hi in excluded:
            if hi <= beg:
                continue
            if lo >= end:
                break
            if beg < lo:
                yield beg, lo
            beg = max(beg, hi)
            if beg >= end:
                break
        if beg < end:
            yield beg, end


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--truth-assembly', type=Path, required=True)
    parser.add_argument('--paftools', type=Path, required=True)
    parser.add_argument('--k8', type=Path, default=Path('k8'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    reference, truth = args.reference, args.truth_assembly
    contig, length, *_ = Path(str(reference) + '.fai').read_text().split()
    forward = sorted(assembly_alignments(reference, truth), key=lambda row: (row[5], int(row[7])))
    reciprocal = {(tgt, strand, qry) for qry, strand, tgt in map(endpoints, assembly_alignments(truth, reference))}
    paired = [row for row in forward if endpoints(row) in reciprocal]
    if not paired or any(row[5] != contig for row in paired):
        raise ValueError('Expected reciprocal alignments of a single reference chromosome')
    intervals = [(int(row[7]), int(row[8])) for row in paired]
    if any(right > left for (_, right), (left, _) in zip(intervals, intervals[1:])):
        raise ValueError('Reciprocal alignments overlap on the reference')
    paf = out / 'assembly.paf'
    paf.write_text(''.join('\t'.join(row) + '\n' for row in paired))
    variants = subprocess.check_output([
        args.k8, args.paftools, 'call', '-q', str(MAPQ), '-l', str(MIN_ALIGNMENT), '-L', str(MIN_ALIGNMENT),
        '-f', reference, paf,
    ], text=True).splitlines()
    records = [line.split('\t') for line in variants if not line.startswith('#')]
    long_events = []
    for row in records:
        if abs(len(row[3]) - len(row[4])) > SMALL_INDEL_LIMIT:
            start = int(row[1]) - 1
            long_events.append((start - BOUNDARY_PADDING, start + len(row[3]) + BOUNDARY_PADDING))
    trimmed = ((beg + BOUNDARY_PADDING, end - BOUNDARY_PADDING) for beg, end in intervals)
    intervals = list(masks(trimmed, sorted(long_events)))
    (out / 'confident.bed').write_text(''.join(f'{contig}\t{beg}\t{end}\n' for beg, end in intervals))
    plain = out / 'small-variants.vcf'
    with plain.open('w') as vcf:
        print('##fileformat=VCFv4.2', f'##contig=<ID={contig},length={length}>',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype">',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMG1655', sep='\n', file=vcf)
        for row in records:
            ref, alt = row[3:5]
            if not set(ref + alt) <= set('ACGT'):
                continue
            small_indel = 0 < abs(len(ref) - len(alt)) <= SMALL_INDEL_LIMIT
            snv = len(ref) == len(alt) == 1 and ref != alt
            if small_indel or snv:
                print(*row[:5], '.', 'PASS', '.', 'GT', '1', sep='\t', file=vcf)
    normalized = out / 'normalized.vcf'
    subprocess.run(['bcftools', 'norm', '-f', reference, '-c', 'e', '-Ov', '-o', normalized, plain], check=True)
    counts = Counter()
    with (out / 'truth.vcf').open('w') as vcf:
        for line in normalized.read_text().splitlines():
            if line.startswith('#'):
                if not line.startswith('##bcftools_'):
                    print(line, file=vcf)
                continue
            row = line.split('\t')
            start = int(row[1]) - 1
            if any(lo <= start and start + len(row[3]) <= hi for lo, hi in intervals):
                counts['snv' if len(row[3]) == len(row[4]) == 1 else 'indel'] += 1
                print(line, file=vcf)
    print(f'{contig}: {sum(hi - lo for lo, hi in intervals):,} scored bases, {dict(counts)}')


if __name__ == '__main__':
    main()

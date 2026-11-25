"""Normalize PASS haploid calls and score them against the independent assembly truth."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess

import pysam


TOOLS = ('lrpipe', 'clair3')
CLASSES = ('snv', 'indel')


def variant_class(ref, alt):
    if ref == alt or not set(ref + alt) <= set('ACGT'):
        return 'other'
    if len(ref) == len(alt) == 1:
        return 'snv'
    if 0 < abs(len(ref) - len(alt)) <= 49:
        return 'indel'
    return 'other'


def normalize(path, out, reference, expected_sample, mask=None):
    normalized = out.with_suffix('.normalized.vcf')
    subprocess.run(['bcftools', 'norm', '--check-ref', 'e', '--fasta-ref', reference,
                    '--atomize', '--atom-overlaps', '.', '-m', '-any', '-Ov', '-o', normalized, path], check=True)
    calls, all_pass, excluded = {}, Counter(), Counter()
    with pysam.VariantFile(normalized) as src, pysam.VariantFile(out, 'wz', header=src.header) as dst:
        if tuple(src.header.samples) != (expected_sample,):
            raise ValueError(f'Unexpected samples in {path}: {tuple(src.header.samples)}')
        for record in src:
            genotype = record.samples[expected_sample].get('GT')
            if genotype is None or len(genotype) != 1:
                raise ValueError(f'Nonhaploid call in {path}: {record.contig}:{record.pos}: {genotype}')
            if genotype != (1,) or tuple(record.filter) != ('PASS',):
                continue
            alternate, = record.alts
            category = variant_class(record.ref, alternate)
            all_pass[category] += 1
            if category == 'other':
                continue
            if mask is not None and not any(start <= record.start and record.stop <= stop
                                            for start, stop in mask.get(record.contig, ())):
                excluded[category] += 1
                continue
            key = record.contig, record.pos, record.ref, alternate
            if key in calls:
                raise ValueError(f'Duplicate normalized call in {path}: {key}')
            calls[key] = category
            dst.write(record)
    pysam.tabix_index(str(out), preset='vcf', force=True)
    return calls, {'pass_before_mask': {key: all_pass[key] for key in (*CLASSES, 'other')},
                   'outside_mask': {key: excluded[key] for key in CLASSES}}


def counts_of_vcf(path):
    with pysam.VariantFile(path) as records:
        return Counter(variant_class(record.ref, record.alts[0]) for record in records)


def rate(top, bottom):
    return top / bottom if bottom else None


def scores(command, scores_dir, truth, evaluated):
    with scores_dir.with_suffix('.log').open('w') as log:
        subprocess.run([*command, '-o', scores_dir], check=True, stdout=log, stderr=subprocess.STDOUT)
    files = {'recall_tp': 'tp-baseline', 'precision_tp': 'tp', 'false_negative': 'fn', 'false_positive': 'fp'}
    buckets = {name: counts_of_vcf(scores_dir / f'{prefix}.vcf.gz') for name, prefix in files.items()}
    result = {}
    for category in CLASSES:
        values = {key: counts[category] for key, counts in buckets.items()}
        tp_call, tp_base, fn, fp = (values[key] for key in ('precision_tp', 'recall_tp', 'false_negative', 'false_positive'))
        if tp_base + fn != truth[category] or tp_call + fp != evaluated[category]:
            raise ValueError(f'Unexpected {category} denominators {values} in {scores_dir}')
        result[category] = {**values, 'truth_variants': tp_base + fn, 'called_variants': tp_call + fp,
                            'precision': rate(tp_call, tp_call + fp), 'recall': rate(tp_base, tp_base + fn)}
    return result


def read_mask(path):
    masks = {}
    for row in path.read_text().splitlines():
        contig, start, stop = row.split('\t')
        masks.setdefault(contig, []).append((int(start), int(stop)))
    return masks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bench', type=Path, required=True)
    parser.add_argument('--rtg', type=Path, required=True)
    parser.add_argument('--platform', choices=('ont', 'hifi'), required=True)
    parser.add_argument('--reference-label', required=True)
    args = parser.parse_args()
    bench, label = args.bench.resolve(strict=True), args.reference_label
    reference = bench / 'inputs' / f'{label}.fa'
    run = bench / 'calls' / label / args.platform
    out = run / 'scores'
    out.mkdir()
    sample = f'K12-{args.platform}'
    use_truth = label != 'mg1655'
    truth_dir = bench / 'truth' / label
    mask = read_mask(truth_dir / 'confident.bed') if use_truth else None
    truth_vcf = truth_dir / 'truth.vcf.gz'
    comparison_dir = {'lrpipe': run / 'lrpipe' / f'{sample}.vcf.gz',
                      'clair3': run / 'clair3/merge_output.vcf.gz'}
    normalized = {}
    all_counts = {}
    for tool, path in comparison_dir.items():
        vcf = out / f'{tool}.vcf.gz'
        normalized[tool], counts = normalize(path, vcf, reference, sample, mask)
        all_counts[tool] = counts
        if use_truth:
            cmd = [args.rtg, 'RTG_MEM=4g', 'vcfeval', '-b', truth_vcf, '-c', vcf,
                   '-t', truth_dir / 'reference.sdf', '--sample-ploidy=1', '--sample', f'MG1655,{sample}',
                   '--threads=2', '--decompose', '--no-roc']
            counts['accuracy'] = scores(cmd, out / tool, counts_of_vcf(truth_vcf), Counter(normalized[tool].values()))
    lrpipe, clair3 = (normalized[name] for name in TOOLS)
    agreement = {}
    for category in CLASSES:
        lhs = {key for key, variant in lrpipe.items() if variant == category}
        rhs = {key for key, variant in clair3.items() if variant == category}
        agreement[category] = {'shared': len(lhs & rhs), 'union': len(lhs | rhs),
                               'lrpipe_only': len(lhs - rhs), 'clair3_only': len(rhs - lhs)}
    report = {'calls': all_counts, 'exact_pass_agreement': agreement}
    (run / 'score-metrics.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'reference': label, 'platform': args.platform, **report}, separators=(',', ':')))


if __name__ == '__main__':
    main()

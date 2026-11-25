"""Run both haploid callers serially on the same alignment from an archive FASTQ."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

PLATFORMS = {
    'ont': ('SRR30993121', 'r1041_e82_400bps_sup_v500'),
    'hifi': ('SRR10971019', 'hifi_sequel2'),
}


def sha256(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def timed(args, run, name, env):
    format_string = '%e\t%U\t%S\t%M'
    metrics = run / f'{name}.time.tsv'
    with (run / f'{name}.log').open('w') as log:
        subprocess.run(['/usr/bin/time', '-o', metrics, '-f', format_string, *args], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    wall, user, system, rss = metrics.read_text().strip().split('\t')
    return {'wall_seconds': float(wall), 'user_seconds': float(user), 'system_seconds': float(system),
            'max_rss_kib': int(rss)}


def check_alignment(alignment, reference_name, expected_reads):
    subprocess.run(['samtools', 'quickcheck', alignment], check=True)
    primary = int(subprocess.check_output(['samtools', 'view', '-c', '-F', '0x900', alignment], text=True))
    if primary != expected_reads:
        raise ValueError(f'Alignment contains {primary} primary records for {expected_reads} FASTQ reads')
    coverage = subprocess.check_output(['samtools', 'coverage', '--ff', '0x904', alignment], text=True)
    record = next(line.split('\t') for line in coverage.splitlines() if line.startswith(reference_name + '\t'))
    return {'primary_records': primary, 'mapped_primary_reads': int(record[3]),
            'reference_bases_covered': int(record[4]), 'mean_reference_depth': float(record[6])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bench', type=Path, required=True)
    parser.add_argument('--lrpipe', type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--clair3-bin', type=Path, required=True)
    parser.add_argument('--bonefish', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--memory', type=int, default=32, help='lrpipe memory budget in GiB')
    parser.add_argument('--platform', choices=PLATFORMS, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--reference-label', required=True)
    args = parser.parse_args()
    platform = args.platform
    accession, model_name = PLATFORMS[platform]
    bench = args.bench.resolve(strict=True)
    reference = args.reference.resolve(strict=True)
    lrpipe = args.lrpipe.resolve(strict=True)
    clair3 = args.clair3_bin.resolve(strict=True)
    bonefish = args.bonefish.resolve(strict=True)
    model = clair3 / 'models' / model_name
    reads = (bench / 'inputs' / f'{accession}.fastq.gz').resolve(strict=True)
    release_sha = (reads.parent / f'{accession}.sha256').read_text().split()[0]
    if sha256(reads) != release_sha:
        raise ValueError(f'Checksum mismatch for {reads}')
    accession_metadata = (reads.parent / f'{accession}.tsv').read_text().splitlines()
    metadata = dict(zip(accession_metadata[0].split('\t'), accession_metadata[1].split('\t'), strict=True))
    ref_contig, ref_bases, *_ = Path(str(reference) + '.fai').read_text().split()
    run = bench / 'calls' / args.reference_label / platform
    run.mkdir(parents=True)
    sample = f'K12-{platform}'
    env = dict(os.environ, LC_ALL='C', TMPDIR=str(bench / 'tmp'), RUST_LOG='info', CUDA_VISIBLE_DEVICES='')
    site = run / 'site.toml'
    site.write_text(f'reference = "{reference}"\n[tools]\nbonefish = "{bonefish}"\n'
                    f'[snv]\nmodels = "{model.parent}"\ndevice = "cpu"\n')
    lr_out = run / 'lrpipe'
    clair_out = run / 'clair3'
    print(f'{args.reference_label} {platform} lrpipe starts', flush=True)
    lr_elapsed = timed([lrpipe, 'run', '--bacterial', '--config', site, '--reads', reads, '--sample', sample,
                        '--platform', platform, '--model', model_name, '--output', lr_out, '--threads', str(args.threads),
                        '--memory', str(args.memory), '--display', 'plain'], run, 'lrpipe', env)
    manifest = json.loads((lr_out / 'manifest.json').read_text())
    if manifest['status'] != 'complete' or manifest['sample'].get('ploidy') != 1:
        raise ValueError(f'lrpipe did not complete as a haploid run in {lr_out}')
    stages = {stage['name']: {field: stage[field] for field in ('status', 'wall_seconds')} for stage in manifest['stages']}
    if list(stages) != ['align', 'snv', 'merge', 'qc'] or any(stage['status'] != 'succeeded' for stage in stages.values()):
        raise ValueError(f'Unexpected lrpipe stages {stages}')
    alignment = lr_out / 'align' / f'{sample}.bam'
    print(f'{args.reference_label} {platform} lrpipe {lr_elapsed["wall_seconds"]:.2f}s, Clair3 starts', flush=True)
    clair_env = dict(env, PATH=f'{clair3}:{env["PATH"]}')
    clair_elapsed = timed([
        clair3 / 'python3', clair3 / 'run_clair3.py', '--bam_fn', alignment, '--ref_fn', reference,
        '--threads', str(args.threads), '--platform', platform, '--model_path', model, '--output', clair_out,
        '--sample_name', sample, '--include_all_ctgs=True', '--haploid_precise=True', '--no_phasing_for_fa=True',
        '--use_gpu=False', '--device=cpu', '--remove_intermediate_dir=True',
    ], run, 'clair3', clair_env)
    expected = [lr_out / f'{sample}.vcf.gz', lr_out / 'snv/merge_output.vcf.gz', clair_out / 'merge_output.vcf.gz']
    for vcf in expected:
        subprocess.run(['bcftools', 'index', '-n', vcf], check=True, stdout=subprocess.DEVNULL)
    accession_reads = int(metadata['read_count'])
    metrics = {'input': {'accession': accession, 'records': accession_reads, 'bases': int(metadata['base_count']),
                         'sha256': release_sha}, 'reference': {'id': ref_contig, 'bases': int(ref_bases),
                                                              'sha256': sha256(reference)},
               'lrpipe': {'code_revision': args.revision, 'binary_sha256': sha256(lrpipe), 'stages': stages,
                          'fastq_to_qc': lr_elapsed}, 'clair3': {'bam_to_vcf': clair_elapsed},
               'alignment': {**check_alignment(alignment, ref_contig, accession_reads), 'sha256': sha256(alignment)},
               'model': {'name': model_name, 'files': {name: sha256(model / name)
                                                       for name in ('pileup.pt', 'full_alignment.pt')}},
               'threads': args.threads, 'memory_gib': args.memory, 'device': 'cpu'}
    (run / 'run-metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    print(f'{args.reference_label} {platform} finished, Clair3 {clair_elapsed["wall_seconds"]:.2f}s', flush=True)


if __name__ == '__main__':
    main()

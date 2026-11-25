#!/usr/bin/env bash
# Regenerates tests/data/split*.bam, split.haplotags.tsv and expected/split.vcf (Sniffles2 2.8.1
# on split.bam) from tests/data/sample1.bam: see tools/make_split_fixture.py. Needs Python with
# pysam and Sniffles 2.8.1; the result is checked in.
#   SNIFFLES=/path/to/sniffles  PYTHON=/path/to/python (with pysam)
set -euo pipefail
env=/data/bioinformatics/micromamba/envs/vareval/bin
sniffles=${SNIFFLES:-$env/sniffles}
python=${PYTHON:-$env/python}
here=$(cd "$(dirname "$0")" && pwd)
data=$here/../tests/data
"$python" -I "$here/make_split_fixture.py" "$data/sample1.bam" "$data"
cd "$data"
"$sniffles" --input split.bam --vcf expected/split.vcf --phase --allow-overwrite

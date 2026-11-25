#!/usr/bin/env bash
# Downloads Spectre 0.2.1 (commit 90580943) into DIR (default target/spectre-0.2.1) for
# tools/run_spectre.py. It is not vendored: parts of it (spectre/libs/unidip) are GPL-3.0.
set -euo pipefail
dir=${1:-target/spectre-0.2.1}
mkdir -p "$dir"
curl -sSL https://github.com/fritzsedlazeck/Spectre/archive/905809435a6bb00aa413a4cf4d4424365112a0af.tar.gz \
    | tar -xz --strip-components=1 -C "$dir"
echo "$dir"

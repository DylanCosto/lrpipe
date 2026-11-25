#!/usr/bin/env bash
# The checks to pass before pushing, run on this machine. GitHub runs them too on every push to
# main and every pull request (.github/workflows/check.yml), where the GPU, parity and WhatsHap
# checks skip for want of a GPU, model and WhatsHap: --release here stays the release gate.
#
#   1. rustfmt, and clippy with -D warnings, on the workspace, and on lrpipe-snv without its GPU
#      backend;
#   2. every crate's tests in release mode, including the parity tests against Sniffles 2.8.1's
#      and mosdepth's output on the bundled fixtures, and snv-infer's without the GPU backend.
#      lrpipe-snv's GPU tests run when a CUDA driver, a GPU and NVRTC are found (NVRTC from
#      $NVRTC_ENV, see below), and must then run rather than skip;
#   3. crates/lrpipe-snv/ci/parity.sh: `lrpipe snv` on its ONT fixture against the VCFs Clair3 2.0.3
#      wrote, with the model from $SNV_MODEL (skipped, saying so, where it is missing);
#   4. lrpipe-phase's `phase` and `haplotag` against WhatsHap 2.8 on WhatsHap's test data
#      (crates/lrpipe-phase/tests/compare.sh and compare_haplotag.sh), with WhatsHap 2.8, samtools and
#      pysam from $WHATSHAP_ENV and WhatsHap's v2.8 source from $WHATSHAP_SRC (skipped, saying
#      so, where they are missing);
#   5. `lrpipe setup --bonefish` and `lrpipe run` from unaligned reads with the real bonefish
#      (crates/lrpipe/tests/align.rs), with bonefish from $BONEFISH, which must be at least
#      $BONEFISH_MIN (skipped, saying so, where it is missing);
#   6. lrpipe-cnv's stages and output files against Spectre 0.2.1's reference runs (its tests
#      reference.rs and outputs.rs), with the runs from $SPECTRE_BENCH and the GRCh38 FASTA from
#      $SPECTRE_REFERENCE (skipped, saying so, where they are missing);
#   7. the tests with the oldest supported Rust (`rust-version` in Cargo.toml), unless --no-msrv.
#
# Every check runs, even after one fails, and the summary lists which passed, failed and were
# skipped. It exits non-zero when one failed; with --release, also when one was skipped, so a
# release is only checked by a run in which every check ran.
#
# usage: tools/check.sh [--no-msrv] [--release]
set -uo pipefail
cd "$(dirname "$0")/.."
msrv=$(sed -n '/^\[workspace.package\]/,/^\[/s/^rust-version = "\(.*\)"/\1/p' Cargo.toml)
no_msrv=false release=false
for arg in "$@"; do
  case $arg in
    --no-msrv) no_msrv=true ;;
    --release) release=true ;;
    *) echo "usage: tools/check.sh [--no-msrv] [--release]" >&2; exit 2 ;;
  esac
done

passed=() failed=() skipped=()
# `${a[@]+"${a[@]}"}` below: macOS's bash 3.2 calls an empty array unbound under `set -u`.
# check NAME COMMAND...: runs COMMAND and records whether it passed.
check() {
  local name=$1
  shift
  echo "== tools/check.sh: $name"
  if "$@"; then passed+=("$name"); else failed+=("$name"); echo "== tools/check.sh: $name FAILED"; fi
}
# skip NAME WHY: records a check that could not run.
skip() {
  echo "== tools/check.sh: $1 skipped: $2"
  skipped+=("$1 ($2)")
}

check rustfmt cargo fmt --all --check
check clippy cargo clippy --locked --workspace --all-targets --all-features -- -D warnings
check "clippy without GPU" cargo clippy --locked -p lrpipe-snv -p snv-infer --all-targets --no-default-features -- -D warnings
# lrpipe-snv's GPU backend compiles its kernels with NVRTC at run time. On this machine NVRTC is in
# a CUDA env; its whole lib directory on the library path would shadow other libraries, so only
# NVRTC is linked into a directory of its own.
NVRTC_ENV=${NVRTC_ENV:-/data2/bioinformatics/envs/torch213-cuda}
if [[ -e $NVRTC_ENV/lib/libnvrtc.so.12 ]]; then
  nvrtc=$(mktemp -d)
  trap 'rm -rf "$nvrtc"' EXIT
  ln -s "$NVRTC_ENV"/lib/libnvrtc.so.12 "$NVRTC_ENV"/lib/libnvrtc-builtins.so.12.* "$nvrtc"/
  export LD_LIBRARY_PATH=$nvrtc${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  if nvidia-smi -L > /dev/null 2>&1; then
    export SNV_REQUIRE_GPU=1
  else
    skip "GPU tests" "no GPU found"
  fi
else
  skip "GPU tests" "no NVRTC at $NVRTC_ENV (set NVRTC_ENV)"
fi
# Exported for lrpipe's test of several --reads too, which compares their SNVs where it is there.
export SNV_MODEL=${SNV_MODEL:-/data2/bioinformatics/envs/clair3-t213/bin/models/r1041_e82_400bps_sup_v500}
check tests cargo test --locked --workspace --release
check "snv-infer tests without GPU" cargo test --locked --release -p snv-infer --no-default-features
snv_parity() {
  cargo build --locked --release -p lrpipe && crates/lrpipe-snv/ci/parity.sh "target/release/lrpipe snv" "$SNV_MODEL"
}
if [[ -e $SNV_MODEL/pileup.pt ]]; then
  check "lrpipe-snv parity" snv_parity
else
  skip "lrpipe-snv parity" "no model at $SNV_MODEL (set SNV_MODEL)"
fi
WHATSHAP_ENV=${WHATSHAP_ENV:-/data2/bioinformatics/envs/clair3-t213}
WHATSHAP_SRC=${WHATSHAP_SRC:-/data2/bioinformatics/github-work/whatshap-audit-2026-10-25/upstream}
whatshap_comparison() {
  cargo build --locked --release -p lrpipe-phase || return 1
  rm -rf target/whatshap-compare
  crates/lrpipe-phase/tests/compare.sh "$WHATSHAP_ENV/bin/whatshap" "$WHATSHAP_ENV/bin/samtools" "$WHATSHAP_SRC" \
    target/whatshap-compare/phase | tail -n 1 || return 1
  crates/lrpipe-phase/tests/compare_haplotag.sh "$WHATSHAP_ENV/bin" "$WHATSHAP_SRC" target/whatshap-compare/haplotag | tail -n 1
}
if [[ -x $WHATSHAP_ENV/bin/whatshap && -d $WHATSHAP_SRC/tests/data ]]; then
  check "lrpipe-phase against WhatsHap" whatshap_comparison
else
  skip "lrpipe-phase against WhatsHap" "no WhatsHap 2.8 at $WHATSHAP_ENV or source at $WHATSHAP_SRC"
fi
SPECTRE_BENCH=${SPECTRE_BENCH:-/data2/bioinformatics/spectre-rs-bench}
SPECTRE_REFERENCE=${SPECTRE_REFERENCE:-/data2/bioinformatics/sniffles-rs-bench/reference.fa}
spectre_parity() {
  SPECTRE_BENCH=$SPECTRE_BENCH SPECTRE_REFERENCE=$SPECTRE_REFERENCE LRPIPE_REQUIRE_SPECTRE=1 \
    cargo test --locked --release -p lrpipe-cnv --test reference --test outputs
}
if [[ -d $SPECTRE_BENCH/inputs && -e $SPECTRE_REFERENCE.fai ]]; then
  check "lrpipe-cnv against Spectre" spectre_parity
else
  skip "lrpipe-cnv against Spectre" "no Spectre reference runs at $SPECTRE_BENCH or no $SPECTRE_REFERENCE (set SPECTRE_BENCH and SPECTRE_REFERENCE)"
fi
BONEFISH=${BONEFISH:-$(command -v bonefish || echo bonefish)}
# The oldest bonefish a release is checked with: the one the release bundle ships
# (tools/bonefish.version), so a stale local build cannot stand in for it.
BONEFISH_MIN=$(< tools/bonefish.version)
bonefish_alignment() {
  local version
  version=$("$BONEFISH" --version | awk '{print $2}')
  if [[ $(printf '%s\n%s\n' "$BONEFISH_MIN" "$version" | sort -V | head -n 1) != "$BONEFISH_MIN" ]]; then
    echo "bonefish $version at $BONEFISH is older than $BONEFISH_MIN: build bonefish's latest release (set BONEFISH)"
    return 1
  fi
  echo "bonefish $version"
  BONEFISH=$BONEFISH LRPIPE_REQUIRE_BONEFISH=1 cargo test --locked --release -p lrpipe --test align
}
if [[ -x $BONEFISH ]]; then
  check "lrpipe run with bonefish" bonefish_alignment
else
  skip "lrpipe run with bonefish" "no bonefish at $BONEFISH (set BONEFISH)"
fi
msrv_tests() {
  { rustup toolchain list | grep -q "^$msrv" || rustup toolchain install --profile minimal "$msrv"; } &&
    cargo "+$msrv" test --locked --workspace --release
}
if $no_msrv; then
  skip "tests with Rust $msrv" "--no-msrv"
else
  check "tests with Rust $msrv" msrv_tests
fi

echo
echo "tools/check.sh: ${#passed[@]} passed, ${#failed[@]} failed, ${#skipped[@]} skipped"
for f in ${failed[@]+"${failed[@]}"}; do echo "  failed: $f"; done
for s in ${skipped[@]+"${skipped[@]}"}; do echo "  skipped: $s"; done
if (( ${#failed[@]} > 0 )); then
  exit 1
elif $release && (( ${#skipped[@]} > 0 )); then
  echo "tools/check.sh: --release needs every check to run"
  exit 1
fi
if (( ${#skipped[@]} > 0 )); then
  echo "tools/check.sh: every check that ran passed, but not every check ran"
else
  echo "tools/check.sh: every check passed"
fi

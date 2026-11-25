#!/usr/bin/env bash
# Builds and tests the workspace, then packs the release bundle: lrpipe, the tools' own binaries
# (lrpipe-snv, lrpipe-sv, lrpipe-cnv, lrpipe-tr, lrpipe-depth, lrpipe-methyl, lrpipe-phase),
# bonefish (the aligner `lrpipe run` takes FASTQ and unaligned BAMs through, at the version in
# tools/bonefish.version, built from its tag), the README, the changelog and every license, into
# dist/lrpipe-VERSION-ARCH-OS.tar.gz with its .sha256.
#
# Run in quay.io/pypa/manylinux_2_28_*, the Linux binaries need only glibc 2.28 (RHEL 8 and
# later, Ubuntu 20.04 and later); run anywhere else, they need the build host's glibc, which the
# script prints. lrpipe-snv's GPU calling loads the NVIDIA driver and NVRTC at run time; the
# bundle includes neither.
#
# usage: tools/release_build.sh VERSION
# env: RUST_TOOLCHAIN (default 1.98.0; installed with rustup when cargo is missing),
#      CARGO_TARGET_DIR (default target),
#      LRPIPE_BUILD_FLAGS (extra cargo flags for the build and the tests, e.g.
#      --no-default-features for a bundle without GPU calling, as on macOS),
#      BONEFISH_SRC (a checkout of bonefish at that version, to build without cloning it)
set -euo pipefail
VERSION=$1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
TOOLCHAIN=${RUST_TOOLCHAIN:-1.98.0}
BINS=(lrpipe lrpipe-snv lrpipe-sv lrpipe-cnv lrpipe-tr lrpipe-depth lrpipe-methyl lrpipe-phase)

if ! command -v cargo > /dev/null; then
  curl -fsSL https://sh.rustup.rs | sh -s -- -y -q --profile minimal --default-toolchain "$TOOLCHAIN"
  export PATH=$HOME/.cargo/bin:$PATH
fi
cd "$ROOT"
# shellcheck disable=SC2086  # the flags are split on purpose
cargo build --release --locked --workspace --bins ${LRPIPE_BUILD_FLAGS:-}
# shellcheck disable=SC2086
cargo test --release --locked --workspace ${LRPIPE_BUILD_FLAGS:-}

# bonefish, at its release tag, with its own lockfile.
BONEFISH_VERSION=$(< tools/bonefish.version)
bonefish_src=${BONEFISH_SRC:-$ROOT/dist/bonefish-$BONEFISH_VERSION}
if [[ -z ${BONEFISH_SRC:-} ]]; then
  rm -rf "$bonefish_src"
  git clone --quiet --depth 1 --branch "v$BONEFISH_VERSION" https://github.com/DylanCosto/bonefish "$bonefish_src"
fi
bonefish_target=$ROOT/dist/bonefish-target
# As bonefish's own release builds it (scripts/package-release.sh): jemalloc, stripped.
CARGO_PROFILE_RELEASE_STRIP=true CARGO_TARGET_DIR=$bonefish_target \
  cargo build --release --locked --manifest-path "$bonefish_src/Cargo.toml" --bin bonefish --features jemalloc

arch=$(uname -m)
[[ $arch == arm64 ]] && arch=aarch64
os=$(uname -s | tr '[:upper:]' '[:lower:]')
[[ $os == darwin ]] && os=macos
name=lrpipe-$VERSION-$arch-$os
stage=$ROOT/dist/$name
rm -rf "$stage" && mkdir -p "$stage/bin" "$stage/licenses"
for b in "${BINS[@]}"; do cp "${CARGO_TARGET_DIR:-target}/release/$b" "$stage/bin/"; done
cp "$bonefish_target/release/bonefish" "$stage/bin/"
cp README.md CHANGELOG.md LICENSE "$stage/"
# Each tool's license and those of the projects it ports, by directory.
for f in crates/*/LICENSE*; do
  tool=${f#crates/}
  mkdir -p "$stage/licenses/${tool%%/*}" && cp "$f" "$stage/licenses/$tool"
done
# pasta, the reader library, with its NOTICE (it began from noodles and ports HTSlib's rules).
pasta=$(cargo metadata --format-version 1 --locked | grep -o '"manifest_path":"[^"]*/pasta/Cargo.toml"' | head -1 | cut -d'"' -f4)
mkdir -p "$stage/licenses/pasta" && cp "$(dirname "$(dirname "$pasta")")"/{LICENSE,NOTICE} "$stage/licenses/pasta/"
# Every other crate compiled into the binaries: its license files, and an index of their terms.
python3 tools/third_party_licenses.py "$stage/licenses/third-party" ${LRPIPE_BUILD_FLAGS:-}
# bonefish's license, and those of the crates compiled into it.
mkdir -p "$stage/licenses/bonefish" && cp "$bonefish_src/LICENSE" "$stage/licenses/bonefish/"
python3 tools/third_party_licenses.py "$stage/licenses/bonefish/third-party" --manifest-path "$bonefish_src/Cargo.toml" \
  --features jemalloc

# The bundle on its own, with nothing from the build environment.
env -i PATH=/usr/bin:/bin "$stage/bin/lrpipe" --version
bundled=$(env -i PATH=/usr/bin:/bin "$stage/bin/bonefish" --version)
[[ $bundled == "bonefish $BONEFISH_VERSION" ]] || { echo "bundled $bundled, not bonefish $BONEFISH_VERSION"; exit 1; }
echo "$bundled"
fixture=crates/lrpipe-sv/tests/data
env -i PATH=/usr/bin:/bin "$stage/bin/lrpipe" sv --input "$fixture/sample1.bam" --vcf "$ROOT/dist/smoke.vcf" \
  --threads 2 --allow-overwrite > /dev/null
cmp <(grep -v '^##' "$ROOT/dist/smoke.vcf") <(grep -v '^##' "$fixture/expected/default.vcf")
rm "$ROOT/dist/smoke.vcf"
if [[ $os == linux ]]; then
  glibc=$(objdump -T "$stage"/bin/* | grep -o 'GLIBC_[0-9.]*' | sort -uV | tail -1)
  echo "needs ${glibc/_/ }"
fi

tar -C "$ROOT/dist" -czf "$ROOT/dist/$name.tar.gz" "$name"
(cd "$ROOT/dist" && { sha256sum "$name.tar.gz" 2> /dev/null || shasum -a 256 "$name.tar.gz"; } > "$name.tar.gz.sha256")
echo "wrote dist/$name.tar.gz"

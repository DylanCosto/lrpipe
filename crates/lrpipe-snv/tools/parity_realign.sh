#!/usr/bin/env bash
# Compares snv-realign's SSW and de Bruijn graph ports with Clair3's C/C++ code on random
# inputs: builds tools/ssw_harness.cpp and tools/dbg_harness.cpp against REALIGN_DIR (Clair3's
# preprocess/realign; the graph needs Boost headers in BOOST_INCLUDE), then runs the opt-in
# tests of crates/snv-realign.
#
# usage: parity_realign.sh REALIGN_DIR OUTDIR [CASES]
set -eu
R=$1 OUT=$2 N=${3:-20000}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"
g++ -O1 -std=c++14 -I"$R" "$HERE/tools/ssw_harness.cpp" "$R/ssw_cpp.cpp" "$R/ssw.c" -o "$OUT/ssw_harness"
g++ -O1 -std=c++11 -I"$R" ${BOOST_INCLUDE:+-I"$BOOST_INCLUDE"} "$HERE/tools/dbg_harness.cpp" "$R/debruijn_graph.cpp" \
  -o "$OUT/dbg_harness"
python3 "$HERE/tools/make_realign_cases.py" ssw "$N" 7 > "$OUT/ssw_cases.tsv"
python3 "$HERE/tools/make_realign_cases.py" dbg "$((N / 5))" 11 > "$OUT/dbg_cases.tsv"
"$OUT/ssw_harness" < "$OUT/ssw_cases.tsv" > "$OUT/ssw_expected.tsv" 2>/dev/null
"$OUT/dbg_harness" < "$OUT/dbg_cases.tsv" > "$OUT/dbg_expected.tsv"
cd "$HERE"
SNV_SSW_CASES=$OUT/ssw_cases.tsv SNV_SSW_EXPECTED=$OUT/ssw_expected.tsv \
SNV_DBG_CASES=$OUT/dbg_cases.tsv SNV_DBG_EXPECTED=$OUT/dbg_expected.tsv \
  cargo test --release -p snv-realign -- matches_the

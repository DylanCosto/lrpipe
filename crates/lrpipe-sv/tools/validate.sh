#!/usr/bin/env bash
# Validates lrpipe-sv on this machine (a VM, an ARM server, ...) and writes a report to share.
#
# usage: tools/validate.sh [OUTDIR]          (default OUTDIR: validation/)
#
# Always: records the platform, builds the release binary, runs the unit tests and the parity
# tests (tests/parity.rs: output must match Sniffles 2.8.1's on the bundled fixture), and times
# the fixture.
#
# Optionally, a full-size comparison with your own data:
#   INPUT=sample.bam|.cram  [REFERENCE=ref.fa]  [ORACLE=sniffles2.vcf]  [SNIFFLES=/path/to/sniffles]
#   [OPTIONS="--minsupport auto ..."]  [THREADS=N]
# With ORACLE, the lrpipe-sv VCF is compared with it. Without ORACLE but with SNIFFLES (2.8.1),
# Sniffles2 is run here first; comparing with a Sniffles2 run on this machine rules out
# platform floating-point differences in the reference output itself.
set -uo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
out=$(mkdir -p "${1:-validation}" && cd "${1:-validation}" && pwd)
host=$(hostname -s 2>/dev/null || hostname)
report=$out/validation-$host-$(uname -m).txt
threads=${THREADS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}
failed=0

log() { echo "$*" | tee -a "$report"; }
step() { log ""; log "== $*"; }
check() { # name, command...
    local name=$1
    shift
    if "$@" >> "$report" 2>&1; then
        log "PASS  $name"
    else
        log "FAIL  $name (details above in $report)"
        failed=1
    fi
}
seconds() { # command... -> elapsed seconds on stdout, command output to the report
    local TIMEFORMAT=%R t
    if t=$({ time "$@" >> "$report" 2>&1; } 2>&1); then
        echo "$t"
    else
        echo "failed after $t"
        return 1
    fi
}

: > "$report"
step "Platform"
log "date:     $(date -u +%Y-%m-%dT%H:%M:%SZ)"
log "host:     $host"
log "kernel:   $(uname -srvm)"
if [[ -r /etc/os-release ]]; then
    log "os:       $(. /etc/os-release && echo "$PRETTY_NAME")"
elif command -v sw_vers > /dev/null; then
    log "os:       $(sw_vers -productName) $(sw_vers -productVersion) ($(sw_vers -buildVersion))"
fi
if command -v lscpu > /dev/null; then
    log "cpu:      $(lscpu | sed -n 's/^Model name: *//p' | head -1) ($(nproc) threads)"
    log "virt:     $(lscpu | sed -n 's/^Hypervisor vendor: *//p' | head -1 | grep . || echo none)"
elif command -v sysctl > /dev/null; then
    log "cpu:      $(sysctl -n machdep.cpu.brand_string 2>/dev/null) ($(sysctl -n hw.ncpu) threads)"
fi
if [[ -r /proc/meminfo ]]; then
    log "memory:   $(awk '/MemTotal/ {printf "%.1f GB", $2 / 1048576}' /proc/meminfo)"
elif command -v sysctl > /dev/null; then
    log "memory:   $(sysctl -n hw.memsize | awk '{printf "%.1f GB", $1 / 1073741824}')"
fi
if command -v ldd > /dev/null; then log "libc:     $(ldd --version 2>&1 | head -1)"; fi
log "rustc:    $(rustc -V)"
log "cargo:    $(cargo -V)"
log "commit:   $(git -C "$root" rev-parse --short HEAD 2>/dev/null)$(git -C "$root" diff --quiet 2>/dev/null || echo ' (modified)')"

step "Build and tests"
cd "$root"
check "release build" cargo build --locked --release
check "unit and parity tests (Sniffles 2.8.1 fixture)" cargo test --locked --release
rs=$root/../../target/release/lrpipe-sv

step "Fixture timing (tests/data/sample1.bam, 5 runs)"
for i in 1 2 3 4 5; do
    s=$(seconds "$rs" --input tests/data/sample1.bam --vcf "$out/fixture.vcf" --threads "$threads" --allow-overwrite)
    log "run $i: $s s"
done

if [[ -n ${INPUT:-} ]]; then
    step "Full-size comparison: $INPUT"
    ref_args=()
    if [[ -n ${REFERENCE:-} ]]; then ref_args=(--reference "$REFERENCE"); fi
    extra=()
    if [[ -n ${OPTIONS:-} ]]; then read -r -a extra <<< "$OPTIONS"; fi
    oracle=${ORACLE:-}
    if [[ -z $oracle && -n ${SNIFFLES:-} ]]; then
        oracle=$out/sniffles2.vcf
        s=$(seconds "$SNIFFLES" --input "$INPUT" --vcf "$oracle" ${ref_args[@]+"${ref_args[@]}"} --threads "$threads" \
            --allow-overwrite ${extra[@]+"${extra[@]}"})
        log "Sniffles2: $s s ($("$SNIFFLES" --version 2>&1 | tail -1))"
    fi
    s=$(seconds "$rs" --input "$INPUT" --vcf "$out/lrpipe-sv.vcf" ${ref_args[@]+"${ref_args[@]}"} --threads "$threads" \
        --allow-overwrite ${extra[@]+"${extra[@]}"})
    log "lrpipe-sv: $s s, $(grep -vc '^#' "$out/lrpipe-sv.vcf" 2>/dev/null) records"
    if [[ -n $oracle ]]; then
        strip() { grep -v -e '^##source=' -e '^##command=' -e '^##fileDate=' "$1"; }
        if diff <(strip "$oracle") <(strip "$out/lrpipe-sv.vcf") > "$out/full.diff"; then
            log "PASS  full-size output identical to $oracle"
        else
            log "FAIL  full-size output differs: $(grep -c '^[<>]' "$out/full.diff") lines ($out/full.diff)"
            failed=1
        fi
    fi
fi

step "Result"
if ((failed)); then log "VALIDATION FAILED on $(uname -m)"; else log "VALIDATION PASSED on $(uname -m)"; fi
log "report: $report"
exit "$failed"

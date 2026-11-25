#!/usr/bin/env bash
# Runs one command (usually `lrpipe-snv run ...` or `run_clair3.py ...`) as a benchmark and
# records next to OUT:
#   OUT.log     its output
#   OUT.time    /usr/bin/time -v (wall time, CPU time, peak RSS; on macOS -l, which adds the
#               peak memory footprint and the instructions retired and cycles)
#   OUT.vmstat  vmstat every 5 s (CPU, I/O wait, swap; on macOS vm_stat, for swap only)
#   OUT.rss     the RSS of the command and its children in KB, every 5 s
#   OUT.gpu     GPU utilization (%) and memory (MiB) every 5 s, when nvidia-smi is there
# then prints a summary: wall time, CPU time, peak RSS, swap-out, GPU use and lrpipe-snv's stage
# times. Nothing else heavy should run on the machine meanwhile.
#
# usage: tools/bench_run.sh OUT COMMAND [ARG...]
set -u
OUT=$1; shift
mkdir -p "$(dirname "$OUT")"
rm -f "$OUT".{log,time,vmstat,rss,gpu}
macos=$([[ $(uname -s) == Darwin ]] && echo 1)

# Seconds since the epoch; BSD date has no %N.
now() {
  if [[ $macos ]]; then perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; else date +%s.%N; fi
}

samplers=()
if [[ $macos ]]; then
  vm_stat 5 > "$OUT.vmstat" & samplers+=($!)
else
  vmstat -n 5 > "$OUT.vmstat" & samplers+=($!)
fi
if command -v nvidia-smi > /dev/null && nvidia-smi > /dev/null 2>&1; then
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits -l 5 > "$OUT.gpu" &
  samplers+=($!)
fi

start=$(now)
if [[ $macos ]]; then
  /usr/bin/time -l -o "$OUT.time" "$@" > "$OUT.log" 2>&1 &
elif [[ -x /usr/bin/time ]]; then
  /usr/bin/time -v -o "$OUT.time" "$@" > "$OUT.log" 2>&1 &
else
  "$@" > "$OUT.log" 2>&1 &
fi
pid=$!
# The command (under time) and every process below it.
descendants() {
  echo "$1"
  local child
  for child in $(pgrep -P "$1"); do descendants "$child"; done
}
( while kill -0 "$pid" 2> /dev/null; do
    ps -o rss= -p "$(descendants "$pid" | paste -sd, -)" | awk '{s += $1} END {print s + 0}'
    sleep 5
  done > "$OUT.rss" ) &
sampler=$!
wait "$pid"
status=$?
end=$(now)
kill "${samplers[@]}" 2> /dev/null
wait "$sampler" 2> /dev/null

printf 'command: %s\n' "$*"
printf 'exit status: %d\n' "$status"
printf 'wall: %.1f s\n' "$(echo "$end - $start" | bc)"
if [[ -s $OUT.time && $macos ]]; then
  # `time -l`: "W real U user S sys", then "VALUE  name" lines, memory in bytes.
  awk 'NR == 1 {u = $3; s = $5} / maximum resident set size$/ {m = $1} / peak memory footprint$/ {f = $1}
    / instructions retired$/ {i = $1} / cycles elapsed$/ {c = $1}
    END {printf "CPU: %.0f s (user %.0f, system %.0f)\npeak RSS (largest process): %.2f GB\n", u + s, u, s, m / 1073741824
      if (f) printf "peak memory footprint: %.2f GB\n", f / 1073741824
      if (c) printf "instructions: %.0f billion, cycles: %.0f billion (time -l)\n", i / 1e9, c / 1e9}' "$OUT.time"
elif [[ -s $OUT.time ]]; then
  awk -F': ' '/User time/ {u = $2} /System time/ {s = $2} /Maximum resident/ {m = $2}
    END {printf "CPU: %.0f s (user %.0f, system %.0f)\npeak RSS (largest process): %.2f GB\n", u + s, u, s, m / 1048576}' \
    "$OUT.time"
fi
[[ -s $OUT.rss ]] && sort -n "$OUT.rss" | tail -1 | awk '{printf "peak RSS (all processes, 5 s samples): %.2f GB\n", $1 / 1048576}'
if [[ $macos ]]; then
  # vm_stat repeats its header, and the first row of numbers after each counts since boot; the
  # others count per 5 s. The last column is pages swapped out (a K or M suffix is thousands or
  # millions).
  awk '/page size of/ {page = $8; total = 1} $1 ~ /^[0-9]/ {if (total) total = 0; else {v = $NF; so += v * (v ~ /K$/ ? 1e3 : v ~ /M$/ ? 1e6 : 1); rows++}}
    END {if (rows) printf "swapped out: %.2f GB\n", so * page / 1073741824}' "$OUT.vmstat"
else
  # vmstat's `so` column is KB/s swapped out; over 5 s samples this is the total swapped out.
  awk 'NR > 2 {so += $8 * 5; wa += $16; n++} END {if (n) printf "swapped out: %.2f GB, mean I/O wait: %.0f%%\n", so / 1048576, wa / n}' \
    "$OUT.vmstat"
fi
if [[ -s $OUT.gpu ]]; then
  awk -F', ' '{n++; u += $1; if ($1 > 0) busy++; if ($2 > m) m = $2}
    END {printf "GPU: mean %.0f%%, busy in %.0f%% of samples, peak memory %.1f GB\n", u / n, 100 * busy / n, m / 1024}' \
    "$OUT.gpu"
fi
grep -E 'done in|calls in' "$OUT.log" | sed 's/.*\] /  /'
exit "$status"

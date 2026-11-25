//! `lrpipe-snv`: the command-line driver of the variant caller.

fn main() -> anyhow::Result<()> {
    #[cfg(feature = "profile")]
    let profiler = profiling::start();
    let result = lrpipe_snv::cli::main(std::env::args_os().collect());
    #[cfg(feature = "profile")]
    profiling::finish(profiler);
    result
}

/// A sampling CPU profiler: with `LRPIPE_SNV_PROFILE=prefix`, a run writes folded stacks for a flame
/// graph and the functions by inclusive and self time to `prefix.PID.{folded,txt}`.
#[cfg(feature = "profile")]
mod profiling {
    use std::collections::{HashMap, HashSet};
    use std::io::Write;

    pub fn start() -> Option<(String, pprof::ProfilerGuard<'static>)> {
        let prefix = format!("{}.{}", std::env::var("LRPIPE_SNV_PROFILE").ok()?, std::process::id());
        Some((prefix, pprof::ProfilerGuardBuilder::default().frequency(997).build().expect("start profiler")))
    }

    pub fn finish(profiler: Option<(String, pprof::ProfilerGuard<'static>)>) {
        let Some((prefix, guard)) = profiler else { return };
        let report = guard.report().build().expect("build profile");
        // Folded stacks, one line per stack as pprof's flame graph reads them (thread, then the
        // frames from the root, inlined ones expanded): `inferno-flamegraph` draws them.
        let mut folded =
            std::io::BufWriter::new(std::fs::File::create(format!("{prefix}.folded")).expect("create folded stacks"));
        for (frames, count) in &report.data {
            let mut line = frames.thread_name_or_id();
            for symbol in frames.frames.iter().rev().flat_map(|frame| frame.iter().rev()) {
                line.push(';');
                line.push_str(&symbol.name());
            }
            writeln!(folded, "{line} {count}").expect("write folded stacks");
        }
        folded.flush().expect("write folded stacks");
        let (mut inclusive, mut own) = (HashMap::<String, isize>::new(), HashMap::<String, isize>::new());
        let mut total = 0;
        for (frames, &count) in &report.data {
            total += count;
            let names: Vec<String> =
                frames.frames.iter().map(|f| f.first().map(|s| s.name()).unwrap_or_default()).collect();
            if let Some(leaf) = names.first() {
                *own.entry(leaf.clone()).or_default() += count;
            }
            let mut seen = HashSet::new();
            for n in names {
                if seen.insert(n.clone()) {
                    *inclusive.entry(n).or_default() += count;
                }
            }
        }
        let mut out = std::fs::File::create(format!("{prefix}.txt")).expect("create profile report");
        for (title, map) in [("inclusive", &inclusive), ("self", &own)] {
            let mut v: Vec<_> = map.iter().collect();
            v.sort_by(|a, b| b.1.cmp(a.1));
            writeln!(out, "== {title} ({total} samples)").expect("write profile report");
            for (name, count) in v.into_iter().take(80) {
                writeln!(out, "{:6.2}%  {name}", 100.0 * *count as f64 / total as f64).expect("write profile report");
            }
        }
    }
}

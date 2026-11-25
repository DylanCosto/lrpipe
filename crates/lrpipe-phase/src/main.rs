//! `lrpipe-phase`: WhatsHap 2.8 in Rust ([`lrpipe_phase::cli`]).

// glibc's allocator keeps freed DP tables reserved per thread, which with many
// threads multiplies peak memory (the `mimalloc` feature, on by default).
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    #[cfg(feature = "profile")]
    let profiler = profiling::start();
    let result = lrpipe_phase::cli::main(std::env::args_os().collect());
    #[cfg(feature = "profile")]
    profiling::finish(profiler);
    if let Err(e) = result {
        eprintln!("lrpipe-phase error: {e:#}");
        std::process::exit(1);
    }
}

#[cfg(feature = "profile")]
mod profiling {
    use std::collections::HashMap;
    use std::io::Write;

    pub fn start() -> Option<(String, pprof::ProfilerGuard<'static>)> {
        let prefix = std::env::var("LRPIPE_PHASE_PROFILE").ok()?;
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
        // inclusive and self sample counts per function
        let (mut inclusive, mut own) = (HashMap::<String, isize>::new(), HashMap::<String, isize>::new());
        let mut total = 0;
        for (frames, &count) in &report.data {
            total += count;
            let names: Vec<String> =
                frames.frames.iter().map(|f| f.first().map(|s| s.name()).unwrap_or_default()).collect();
            if let Some(leaf) = names.first() {
                *own.entry(leaf.clone()).or_default() += count;
            }
            let mut seen = std::collections::HashSet::new();
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
            writeln!(out, "== {title} ({total} samples)").unwrap();
            for (name, count) in v.into_iter().take(60) {
                writeln!(out, "{:6.2}%  {name}", 100.0 * *count as f64 / total as f64).unwrap();
            }
        }
    }
}

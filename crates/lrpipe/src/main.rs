//! lrpipe: long-read calling in one binary. `lrpipe run` runs the pipeline on a sample; `snv`,
//! `sv`, `cnv`, `tr`, `depth`, `methyl` and `phase` run its tools, each taking exactly the
//! arguments of the tool's own binary, and `haplotag` is `lrpipe phase haplotag`.

use std::ffi::OsString;
use std::process::ExitCode;

/// The tools' subcommands.
#[derive(Clone, Copy, Debug)]
enum Tool {
    Snv,
    Sv,
    Cnv,
    Tr,
    Depth,
    Methyl,
    Phase,
    Haplotag,
}

impl Tool {
    fn named(name: &str) -> Option<Self> {
        Some(match name {
            "snv" => Tool::Snv,
            "sv" => Tool::Sv,
            "cnv" => Tool::Cnv,
            "tr" => Tool::Tr,
            "depth" => Tool::Depth,
            "methyl" => Tool::Methyl,
            "phase" => Tool::Phase,
            "haplotag" => Tool::Haplotag,
            _ => return None,
        })
    }
}

fn main() -> ExitCode {
    let args: Vec<OsString> = std::env::args_os().collect();
    if args.get(1).is_some_and(|a| a == "sv") {
        fix_mmap_threshold();
    }
    #[cfg(feature = "profile")]
    let profiler = profiling::start(&args.get(1).map_or_else(String::new, |a| a.to_string_lossy().into_owned()));
    let code = dispatch(args);
    #[cfg(feature = "profile")]
    profiling::finish(profiler);
    code
}

/// For `lrpipe sv`, the pipeline's read pass: keeps glibc's size from which an allocation gets a
/// mapping of its own at its default, 128 KiB. glibc raises it each time such a block is freed,
/// after which the read pass's large buffers (a long read's bases and calls run to hundreds of
/// KiB) come from the threads' heaps, which keep the memory once they are freed: on HG002 ONT
/// chr1-chr3 the read pass with methylation peaked at 7.9 GB, and at 5.6 GB with the threshold
/// fixed, as fast. Not for the other tools, which then map and unmap a buffer per read: `lrpipe
/// methyl cpg` on the whole Revio genome took 375 s with it (900 s of it in the kernel) and 212 s
/// without, at the same peak.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
#[expect(unsafe_code, reason = "mallopt has no safe binding")]
fn fix_mmap_threshold() {
    // SAFETY: mallopt only sets an allocator parameter, before any other thread starts.
    unsafe {
        libc::mallopt(libc::M_MMAP_THRESHOLD, 128 * 1024);
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
fn fix_mmap_threshold() {}

/// The tools' subcommands, as [`Tool::named`] knows them.
const TOOLS: [&str; 8] = ["snv", "sv", "cnv", "tr", "depth", "methyl", "phase", "haplotag"];

fn dispatch(mut args: Vec<OsString>) -> ExitCode {
    // `lrpipe help TOOL [COMMAND]` is `lrpipe TOOL [COMMAND] --help`, as `lrpipe help run` is clap's.
    if args.get(1).is_some_and(|a| a == "help")
        && args.get(2).and_then(|a| a.to_str()).is_some_and(|a| Tool::named(a).is_some())
    {
        // The tool's own command too: `lrpipe help snv pileup` is `lrpipe snv pileup --help`.
        let mut rewritten = vec![args[0].clone()];
        rewritten.extend(args.drain(2..));
        rewritten.push("--help".into());
        args = rewritten;
    }
    if let Some(code) = unknown_command(&args) {
        return code;
    }
    let tool = args.get(1).and_then(|a| a.to_str()).and_then(|a| Tool::named(a).map(|t| (t, a.to_owned())));
    let Some((tool, name)) = tool else { return lrpipe::cli::main(args) };
    // The tool sees `lrpipe <TOOL>` as its program name, in help and in VCF headers.
    let mut argv: Vec<OsString> = vec![format!("lrpipe {name}").into()];
    let mut rest = args.into_iter().skip(2).peekable();
    let result = match tool {
        Tool::Snv => {
            // lrpipe-snv's `run` is implied unless another of its commands is named.
            if rest.peek().is_none_or(|a| !matches!(a.to_str(), Some("run" | "pileup" | "phase" | "help"))) {
                argv.push("run".into());
            }
            lrpipe_snv::cli::main(with(argv, rest))
        }
        Tool::Sv => lrpipe_sv::cli::main(with(argv, rest)),
        Tool::Cnv => {
            // lrpipe-cnv's one subcommand is implied.
            if rest.peek().is_none_or(|a| a != "CNVCaller") {
                argv.push("CNVCaller".into());
            }
            lrpipe_cnv_cli::cli::main(with(argv, rest))
        }
        Tool::Tr => lrpipe_tr::cli::main(with(argv, rest)),
        Tool::Depth => lrpipe_depth::cli::main(with(argv, rest)),
        Tool::Methyl => lrpipe_methyl::cli::main(with(argv, rest)),
        Tool::Phase => {
            // WhatsHap's `phase` is implied unless another of its commands is named.
            const COMMANDS: [&str; 10] = [
                "phase",
                "haplotag",
                "genotype",
                "stats",
                "compare",
                "split",
                "unphase",
                "hapcut2vcf",
                "haplotagphase",
                "help",
            ];
            // Then the program is `lrpipe` and the command `phase`, so help reads `lrpipe phase`.
            if rest.peek().is_none_or(|a| a.to_str().is_none_or(|a| !COMMANDS.contains(&a))) {
                argv = vec!["lrpipe".into(), "phase".into()];
            }
            lrpipe_phase::cli::main(with(argv, rest))
        }
        Tool::Haplotag => lrpipe_phase::cli::main(with(vec!["lrpipe".into(), "haplotag".into()], rest)),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {e:?}");
            ExitCode::FAILURE
        }
    }
}

/// For a first argument that is no command (nor an option), clap's error with the nearest of
/// the tools' and lrpipe's own commands, which clap alone does not know the tools of; none
/// when it is a command, or nothing is near it (clap then says it is unknown).
fn unknown_command(args: &[OsString]) -> Option<ExitCode> {
    let given = args.get(1)?.to_str()?;
    let own = lrpipe::cli::commands();
    if given.starts_with('-') || Tool::named(given).is_some() || own.iter().any(|c| c == given) {
        return None;
    }
    let all = TOOLS.iter().map(|t| (*t).to_owned()).chain(own);
    let (near, distance) = all.map(|c| (edit_distance(given, &c), c)).min().map(|(d, c)| (c, d))?;
    if distance > 2 || distance >= given.len() {
        return None;
    }
    eprintln!(
        "error: unrecognized subcommand '{given}'\n\n  tip: a similar subcommand exists: '{near}'\n\nUsage: lrpipe \
         <COMMAND>\n\nFor more information, try '--help'."
    );
    Some(ExitCode::from(2))
}

/// Edit distance between `a` and `b`, by characters, with a swap of neighbours one edit
/// (`rnu` is one from `run`).
fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    // d[i][j]: the distance between a's first i characters and b's first j.
    let mut d: Vec<Vec<usize>> = (0..=a.len()).map(|i| vec![i; b.len() + 1]).collect();
    d[0] = (0..=b.len()).collect();
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let substitution = d[i - 1][j - 1] + usize::from(a[i - 1] != b[j - 1]);
            d[i][j] = (d[i - 1][j] + 1).min(d[i][j - 1] + 1).min(substitution);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                d[i][j] = d[i][j].min(d[i - 2][j - 2] + 1);
            }
        }
    }
    d[a.len()][b.len()]
}

/// The program name followed by the remaining arguments.
fn with(mut argv: Vec<OsString>, rest: impl IntoIterator<Item = OsString>) -> Vec<OsString> {
    argv.extend(rest);
    argv
}

#[cfg(feature = "profile")]
mod profiling {
    use std::collections::{HashMap, HashSet};
    use std::io::Write;

    /// Starts sampling when `LRPIPE_PROFILE` is set; the report goes to
    /// `$LRPIPE_PROFILE.COMMAND.PID.{folded,txt}`, so a pipeline's stages each write their own.
    pub fn start(command: &str) -> Option<(String, pprof::ProfilerGuard<'static>)> {
        let prefix = format!("{}.{command}.{}", std::env::var("LRPIPE_PROFILE").ok()?, std::process::id());
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

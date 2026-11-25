//! The binary's subcommands reach the pipeline and the tools.

use std::process::Command;

fn lrpipe(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_lrpipe")).args(args).output().expect("run lrpipe")
}

#[test]
fn every_tool_subcommand_reaches_its_tool() {
    for (command, expect) in [
        ("sv", "Usage: lrpipe sv"),
        ("cnv", "Usage: lrpipe cnv CNVCaller"),
        ("tr", "Usage: lrpipe tr"),
        ("depth", "Usage: lrpipe depth"),
        ("methyl", "Usage: lrpipe methyl"),
        ("phase", "Usage: lrpipe phase [OPTIONS]"),
        ("haplotag", "Usage: lrpipe haplotag [OPTIONS]"),
    ] {
        let out = lrpipe(&[command, "--help"]);
        assert!(out.status.success(), "{command}");
        let text = String::from_utf8_lossy(&out.stdout);
        assert!(text.contains(expect), "{command}: {text}");
    }
}

#[test]
fn the_pipeline_commands_run() {
    for (command, expect) in
        [("run", "Usage: lrpipe run"), ("sex", "Usage: lrpipe sex"), ("compare", "Usage: lrpipe compare")]
    {
        let out = lrpipe(&[command, "--help"]);
        assert!(out.status.success(), "{command}");
        let text = String::from_utf8_lossy(&out.stdout);
        assert!(text.contains(expect), "{command}: {text}");
    }
    let out = lrpipe(&["sex", "/nonexistent.bam"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).starts_with("lrpipe: error:"));
}

#[test]
fn help_lists_the_tools_and_unknown_commands_fail() {
    let help = lrpipe(&["--help"]);
    assert!(help.status.success());
    let text = String::from_utf8_lossy(&help.stdout);
    for tool in ["sv", "cnv", "tr", "depth", "methyl", "phase", "haplotag"] {
        assert!(text.contains(&format!("\n  {tool} ")), "{tool}: {text}");
    }
    let out = lrpipe(&["bogus"]);
    assert_eq!(out.status.code(), Some(2));
    assert_eq!(lrpipe(&[]).status.code(), Some(2));
}

#[test]
fn phase_implies_whatshap_phase_unless_another_command_is_named() {
    let out = lrpipe(&["phase", "stats", "--help"]);
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).contains("Usage: lrpipe phase stats"));
}

#[test]
fn tool_errors_fail() {
    let out = lrpipe(&["depth", "out", "/nonexistent.bam", "-n"]);
    assert_eq!(out.status.code(), Some(1));
}

/// LUX-06: every subcommand takes --version and -V, the tools' own subcommands too.
#[test]
fn every_subcommand_prints_its_version() {
    let commands: [&[&str]; 18] = [
        &[],
        &["run"],
        &["setup"],
        &["sex"],
        &["merge"],
        &["compare"],
        &["snv"],
        &["snv", "pileup"],
        &["sv"],
        &["cnv"],
        &["tr"],
        &["depth"],
        &["methyl"],
        &["methyl", "pileup"],
        &["methyl", "cpg"],
        &["phase"],
        &["phase", "stats"],
        &["haplotag"],
    ];
    for command in commands {
        for flag in ["--version", "-V"] {
            let out = lrpipe(&[command, &[flag]].concat());
            let text = String::from_utf8_lossy(&out.stdout);
            assert!(out.status.success(), "{command:?} {flag}: {}", String::from_utf8_lossy(&out.stderr));
            assert!(text.trim_end().ends_with(env!("CARGO_PKG_VERSION")), "{command:?} {flag}: {text}");
        }
    }
}

/// LUX-17: `lrpipe help TOOL` is the tool's help; a mistyped command is told the nearest; and
/// a tool whose command is implied lists its others in that command's help.
#[test]
fn help_names_the_tools_and_their_commands() {
    let out = lrpipe(&["help", "snv"]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stdout).contains("Usage: lrpipe snv run"));
    // The tool's own command too.
    for (command, usage) in [
        (&["help", "snv", "pileup"][..], "Usage: lrpipe snv pileup"),
        (&["help", "phase", "stats"], "Usage: lrpipe phase stats"),
    ] {
        let out = lrpipe(command);
        assert!(out.status.success(), "{command:?}");
        assert!(
            String::from_utf8_lossy(&out.stdout).contains(usage),
            "{command:?}: {}",
            String::from_utf8_lossy(&out.stdout)
        );
    }
    for (typo, near) in [("snvv", "'snv'"), ("rnu", "'run'"), ("methly", "'methyl'")] {
        let out = lrpipe(&[typo]);
        assert_eq!(out.status.code(), Some(2), "{typo}");
        let text = String::from_utf8_lossy(&out.stderr);
        assert!(text.contains(&format!("a similar subcommand exists: {near}")), "{typo}: {text}");
    }
    for (tool, others) in [("snv", "pileup and phase"), ("phase", "haplotag, genotype, stats")] {
        let out = lrpipe(&[tool, "--help"]);
        let text = String::from_utf8_lossy(&out.stdout);
        assert!(text.contains(&format!("Other commands: {others}")), "{tool}: {text}");
    }
}

/// LUX-16: `lrpipe tr` takes --karyotype XY as its help writes it, --reference as the other
/// tools name the reference, and indexes a bgzipped VCF.
#[test]
fn tr_takes_its_helps_spellings_and_indexes_its_vcf() {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("tr-spellings");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch directory");
    std::fs::write(dir.join("catalog.bed"), "chr20\t30200000\t30200030\tAT\n").expect("catalog");
    let path = |p: &std::path::Path| p.to_str().expect("UTF-8").to_owned();
    let (reads, reference) = (path(&fixture.join("ont.bam")), path(&fixture.join("chr20_window.fa.gz")));
    let (catalog, vcf) = (path(&dir.join("catalog.bed")), path(&dir.join("repeats.vcf.gz")));
    let out = lrpipe(&[
        "tr",
        "--reads",
        &reads,
        "--reference",
        &reference,
        "--repeats",
        &catalog,
        "--output",
        &vcf,
        "--karyotype",
        "XY",
        "--platform",
        "ont",
        "-t",
        "2",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let index = dir.join("repeats.vcf.gz.tbi");
    let tabix = std::fs::read(&index).expect("a tabix index");
    let mut magic = Vec::new();
    std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(&tabix[..]), &mut magic).expect("BGZF");
    assert!(magic.starts_with(b"TBI\x01"), "{}", index.display());
    let help = String::from_utf8_lossy(&lrpipe(&["tr", "--help"]).stdout).into_owned();
    assert!(help.lines().any(|l| l.contains("-t, --threads") && l.contains("Worker threads")), "{help}");
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn compare_into_a_closed_pipe_ends_quietly() {
    // `lrpipe compare a b | head` panicked ("failed printing to stdout: Broken pipe") and aborted.
    let dir = std::env::temp_dir().join(format!("lrpipe-compare-pipe-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    for run in ["a", "b"] {
        let root = dir.join(run);
        std::fs::create_dir_all(&root).expect("directory");
        let mut outputs = serde_json::Map::new();
        for i in 0..200 {
            let path = root.join(format!("out{i}.txt"));
            std::fs::write(&path, format!("#h\n{run}\n")).expect("written");
            outputs.insert(format!("kind{i}"), path.to_string_lossy().into_owned().into());
        }
        let manifest = serde_json::json!({"output": root, "stages": [], "tools": {}, "outputs": outputs});
        std::fs::write(root.join("manifest.json"), manifest.to_string()).expect("written");
    }
    let mut child = std::process::Command::new(env!("CARGO_BIN_EXE_lrpipe"))
        .arg("compare")
        .args([dir.join("a"), dir.join("b")])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("lrpipe runs");
    // The reader gone before lrpipe writes, as `head` once it has its lines.
    drop(child.stdout.take());
    let out = child.wait_with_output().expect("lrpipe ends");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(1), "the outputs differ: {stderr}");
    assert!(stderr.is_empty(), "{stderr}");
    let _ = std::fs::remove_dir_all(dir);
}

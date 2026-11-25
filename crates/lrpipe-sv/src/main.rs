//! lrpipe-sv: structural variant calling from long reads, ported from Sniffles2 2.8.1.

fn main() -> anyhow::Result<()> {
    lrpipe_sv::cli::main(std::env::args_os().collect())
}

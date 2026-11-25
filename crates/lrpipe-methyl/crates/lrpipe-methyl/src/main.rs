//! lrpipe-methyl: base modification pileups from long reads.

fn main() -> anyhow::Result<()> {
    lrpipe_methyl::cli::main(std::env::args_os().collect())
}

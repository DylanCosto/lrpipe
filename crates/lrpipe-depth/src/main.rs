//! lrpipe-depth: binned read depth identical to mosdepth.

fn main() -> anyhow::Result<()> {
    lrpipe_depth::cli::main(std::env::args_os().collect())
}

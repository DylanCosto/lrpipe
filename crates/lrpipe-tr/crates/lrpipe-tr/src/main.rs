//! lrpipe-tr: tandem repeat genotyping from long reads.

fn main() -> anyhow::Result<()> {
    lrpipe_tr::cli::main(std::env::args_os().collect())
}

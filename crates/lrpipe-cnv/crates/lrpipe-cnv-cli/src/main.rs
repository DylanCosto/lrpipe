//! lrpipe-cnv: long-read copy number calling, a Spectre 0.2.1 port.

fn main() -> anyhow::Result<()> {
    lrpipe_cnv_cli::cli::main(std::env::args_os().collect())
}

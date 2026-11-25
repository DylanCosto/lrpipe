//! BGZF files with tabix indexes, as `pysam.tabix_compress` and `tabix_index` make them.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use pasta::bgzf;
use pasta::core::Position;
use pasta::csi::binning_index::index::header::Builder;
use pasta::csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::Error;

/// A data line and the 1-based inclusive interval it covers.
#[derive(Debug)]
pub(crate) struct Line {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub text: String,
}

/// Which tabix preset the index records.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Preset {
    Bed,
    Vcf,
}

/// Writes `header` (meta lines, not indexed) and `lines` to a BGZF file at `path`, and a
/// tabix index at `path.tbi`.
pub(crate) fn write(path: &Path, header: &str, lines: &[Line], preset: Preset) -> Result<(), Error> {
    let file = File::create(path).map_err(Error::io(path))?;
    let mut writer = bgzf::io::Writer::new(BufWriter::new(file));
    writer.write_all(header.as_bytes()).map_err(Error::io(path))?;
    let mut indexer = pasta::tabix::index::Indexer::default();
    indexer.set_header(match preset {
        Preset::Bed => Builder::bed().build(),
        Preset::Vcf => Builder::vcf().build(),
    });
    for l in lines {
        let start = writer.virtual_position();
        writer.write_all(l.text.as_bytes()).and_then(|()| writer.write_all(b"\n")).map_err(Error::io(path))?;
        let end = writer.virtual_position();
        // htslib clamps intervals that start before 1 (a VCF record at POS 0).
        let first = Position::try_from(l.start.max(1) as usize).expect("at least 1");
        let last = Position::try_from(l.end.max(l.start).max(1) as usize).expect("at least 1");
        indexer.add_record(&l.contig, first, last, Chunk::new(start, end)).map_err(Error::io(path))?;
    }
    writer.finish().and_then(|mut w| w.flush()).map_err(Error::io(path))?;
    let index = indexer.build();
    let tbi = path.with_added_extension("tbi");
    pasta::tabix::fs::write(&tbi, &index).map_err(Error::io(&tbi))?;
    Ok(())
}

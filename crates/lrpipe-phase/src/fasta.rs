//! Reference FASTA access (faidx-indexed, plain or bgzip-compressed).

use std::fs::File;
use std::io::{Seek, SeekFrom};

use anyhow::{Context, Result};
use pasta::core::{Position, Region};
use pasta::fasta;

/// Refuses a reference whose `.fai` is another FASTA's index (it would read wrong bases without
/// an error); a missing `.fai` is reported where the reference is opened.
pub fn check_reference(path: &str) -> Result<()> {
    if !std::path::Path::new(&format!("{path}.fai")).exists() {
        return Ok(());
    }
    lrpipe_reads::check_fai(std::path::Path::new(path))
        .map_err(lrpipe_reads::OpenError::into_io_error)
        .with_context(|| format!("checking the reference {path} (and its .fai)"))
}

pub struct Fasta {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
}

impl Fasta {
    /// Opens `path` with the index `path.fai` (and `path.gzi` if it ends in
    /// .gz or .bgz).
    pub fn open(path: &str) -> Result<Fasta> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("could not open {path} with its FASTA index"))?;
        Ok(Fasta { reader })
    }

    pub fn names(&self) -> Vec<String> {
        self.reader.index().as_ref().iter().map(|r| r.name().to_string()).collect()
    }

    /// The whole sequence of `name` (as stored, without changing case).
    pub fn fetch(&mut self, name: &str) -> Result<Vec<u8>> {
        let index = self.reader.index();
        let Some(record) = index.as_ref().iter().find(|r| r.name() == name) else {
            return Ok(Vec::new());
        };
        let length = usize::try_from(record.length())?;
        let Some(end) = Position::new(length) else {
            return Ok(Vec::new());
        };
        // Read into one buffer of the final size (IndexedReader::query would
        // return the sequence in a buffer that has to be copied).
        let start = index.query(&Region::new(name, Position::MIN..=end))?;
        let mut seq = Vec::with_capacity(length);
        self.reader.get_mut().seek(SeekFrom::Start(start))?;
        self.reader
            .read_sequence(&mut seq)
            .with_context(|| format!("could not fetch sequence {name} from the reference"))?;
        seq.truncate(length);
        Ok(seq)
    }
}

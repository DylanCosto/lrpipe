//! `{sample}_cnv.bed.gz` (outputWriter.py:13-39).

use std::path::Path;

use pynum::{repr::py_repr, round::np_round};

use super::tabix::{self, Line, Preset};
use crate::Error;
use crate::candidate::Arena;
use crate::merged::Call;

/// One line per call, per contig sorted by start (stable): contig, start and end (first and
/// last window start for a CNV, first and last SNV position for LoH), type, size, and
/// `median_cov_norm` (a LoH call's copy number) rounded to 2 decimals the numpy way.
///
/// # Errors
///
/// If the file or its index cannot be written.
pub fn write(path: &Path, arena: &Arena, calls: &[(String, Vec<Call>)]) -> Result<(), Error> {
    let mut lines = Vec::new();
    for (_, list) in calls {
        let mut sorted: Vec<&Call> = list.iter().collect();
        sorted.sort_by_key(|c| c.start(arena));
        for call in sorted {
            let (chrom, start, end, kind, size, median) = match call {
                Call::Cnv(id) => {
                    let c = arena.get(*id);
                    (&c.chromosome, c.start, c.end, c.cnv_type.as_str(), c.size, c.median_cov_norm)
                }
                Call::Loh(l) => (&l.chromosome, l.start, l.end, "LOH", l.size, l.median_cov_norm()),
            };
            let text = format!("{chrom}\t{start}\t{end}\t{kind}\t{size}\t{}", py_repr(np_round(median, 2)));
            lines.push(Line { contig: chrom.clone(), start: start + 1, end, text });
        }
    }
    tabix::write(path, "", &lines, Preset::Bed)
}

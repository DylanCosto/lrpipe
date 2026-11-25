//! `check_malformed_records` from `SortVcf.py`: repairs lines where two VCF lines were fused.
//!
//! The repairs are heuristics that can themselves produce odd rows (for example a bare contig
//! name); they are reproduced exactly, since the merged output keeps whatever they return.

use std::path::Path;

use crate::RowError;
use crate::pytext::{char_len, find, head, split_whitespace, strip, tab_fields, tail};

/// What a header line is split at when it holds a fused record.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Anchor<'a> {
    /// The contig being merged (`contig=`; the per-contig path).
    Contig(&'a str),
    /// The sample name (`sample_name=`; the all-files path, where `contig` is `None`).
    Sample(&'a str),
}

/// The `(row, extra_row)` pair `check_malformed_records` returns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Checked {
    /// The line to keep in place of the input, if any.
    pub(crate) row: Option<String>,
    /// A second record split off the line, if any.
    pub(crate) extra: Option<String>,
}

/// `check_malformed_records(row, ..., check_header=True)` for a header line.
pub(crate) fn check_header(row: &str, anchor: Anchor<'_>, file: &Path) -> Result<Checked, RowError> {
    if !(row.starts_with("#CHROM") && tab_fields(row) != 10) {
        return Ok(Checked { row: Some(row.to_owned()), extra: None });
    }
    log::warn!("[WARNING] Malformed header line:\n {row} in {}, split into two lines", file.display());
    let split_at = match anchor {
        Anchor::Sample(sample) => find(row, sample) + char_len(sample),
        Anchor::Contig(contig) => find(row, contig),
    };
    let header_row = format!("{}\n", head(row, split_at));
    let record = tail(row, split_at);
    if tab_fields(record) > 10 {
        return Ok(Checked { row: None, extra: None });
    }
    let columns = split_whitespace(strip(record), Some(3));
    let ctg_name = columns.first().ok_or(RowError::MissingColumn { index: 0 })?;
    columns.get(1).ok_or(RowError::MissingColumn { index: 1 })?;
    // With a sample anchor `contig` is None, which no contig name equals: the header line and
    // its record are both dropped.
    let same_contig = matches!(anchor, Anchor::Contig(contig) if *ctg_name == contig);
    let row = same_contig.then_some(header_row);
    Ok(Checked { row, extra: Some(record.to_owned()) })
}

/// `check_malformed_records(row, contig=contig, ...)` for a record line.
pub(crate) fn check_record(row: &str, contig: &str, file: &Path) -> Checked {
    let fields = tab_fields(row);
    if fields == 10 {
        return Checked { row: Some(row.to_owned()), extra: None };
    }
    if fields > 19 {
        log::warn!("[WARNING] Malformed header line:\n {row} in {}, skip!", file.display());
        return Checked { row: None, extra: None };
    }
    log::warn!("[WARNING] Malformed VCF line:\n {row} in {}, split into two lines", file.display());
    // The second occurrence of the contig name; if there is none, find() returns -1 and the
    // split lands one character before the end of the first name.
    let len = char_len(contig);
    let split_at = find(tail(row, len), contig) + len;
    Checked { row: Some(format!("{}\n", head(row, split_at))), extra: Some(tail(row, split_at).to_owned()) }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILE: &str = "pileup_chr20_1.vcf";

    fn record(pos: u32) -> String {
        format!("chr20\t{pos}\t.\tA\tC\t5.55\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:5:3:1,2:0.6667\n")
    }

    #[test]
    fn well_formed_lines_pass_through() {
        let header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        let got = check_header(header, Anchor::Contig("chr20"), Path::new(FILE)).expect("valid");
        assert_eq!(got, Checked { row: Some(header.to_owned()), extra: None });
        let got = check_record(&record(1), "chr20", Path::new(FILE));
        assert_eq!(got, Checked { row: Some(record(1)), extra: None });
    }

    #[test]
    fn fused_header_is_split() {
        let fused = format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE{}", record(7));
        let header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        let got = check_header(&fused, Anchor::Contig("chr20"), Path::new(FILE)).expect("valid");
        assert_eq!(got, Checked { row: Some(header.to_owned()), extra: Some(record(7)) });
        // The all-files path drops both halves.
        let got = check_header(&fused, Anchor::Sample("SAMPLE"), Path::new(FILE)).expect("valid");
        assert_eq!(got, Checked { row: None, extra: Some(record(7)) });
        // A record of another contig: find() misses and the remainder is a lone newline.
        let other = fused.replace("SAMPLEchr20", "SAMPLEchr21");
        let got = check_header(&other, Anchor::Contig("chr20"), Path::new(FILE));
        assert_eq!(got, Err(RowError::MissingColumn { index: 0 }));
    }

    #[test]
    fn fused_records_are_split() {
        let two = format!("{}{}", record(1).trim_end(), record(2));
        let got = check_record(&two, "chr20", Path::new(FILE));
        assert_eq!(got, Checked { row: Some(record(1)), extra: Some(record(2)) });
        let three = format!("{}{}{}", record(1).trim_end(), record(2).trim_end(), record(3));
        assert_eq!(check_record(&three, "chr20", Path::new(FILE)), Checked { row: None, extra: None });
        // A short record has no second contig name: it is cut after "chr2".
        let short = "chr20\t9\t.\tA\tC\t5.55\tPASS\tP\tGT\n";
        let got = check_record(short, "chr20", Path::new(FILE));
        assert_eq!(got, Checked { row: Some("chr2\n".to_owned()), extra: Some(short[4..].to_owned()) });
    }
}

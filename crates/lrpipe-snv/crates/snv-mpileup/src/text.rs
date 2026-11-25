//! The read-bases column of `samtools mpileup --reverse-del` without a reference (`pileup_seq`
//! in samtools 1.24's `bam_plcmd.c`), which the original's Python parses back into alleles.

use std::fmt::Write as _;

use crate::hts::{Column, Entry};

/// Forward-strand base letters by `seq_nt16` code; `=` (code 0) becomes `.`.
const UPPER: &[u8; 16] = b".ACMGRSVTWYHKDBN";
/// Reverse-strand base letters; `=` becomes `,`.
const LOWER: &[u8; 16] = b",acmgrsvtwyhkdbn";
/// `seq_nt16_str`.
const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// The bases column of one pileup column: every read's entry, in pileup order (the original
/// keeps every base, as it passes `--min-BQ 0`), or `*` for a column without reads.
pub fn read_bases(column: &Column<'_>) -> String {
    if column.entries().len() == 0 {
        return "*".to_owned();
    }
    let mut out = String::new();
    for entry in column.entries() {
        push_entry(&mut out, &entry);
    }
    out
}

/// `pileup_seq` for one read.
pub fn push_entry(out: &mut String, entry: &Entry<'_>) {
    let record = entry.record();
    let reverse = record.is_reverse();
    if entry.is_head() {
        out.push('^');
        let mapq = record.mapq();
        out.push(char::from(if mapq > 93 { 126 } else { mapq + 33 }));
    }
    if entry.is_del() {
        out.push(match (entry.is_refskip(), reverse) {
            (true, true) => '<',
            (true, false) => '>',
            (false, true) => '#',
            (false, false) => '*',
        });
    } else {
        let code = if entry.qpos() < record.len() { record.base_code(entry.qpos()) } else { 15 };
        let letters = if reverse { LOWER } else { UPPER };
        out.push(char::from(letters[usize::from(code)]));
    }
    let mut del_len = -entry.indel();
    if entry.indel() > 0 {
        let (bases, deletion) = insertion(entry);
        del_len = deletion;
        write!(out, "+{}", bases.len()).expect("writing to a String cannot fail");
        for b in bases {
            out.push(match (reverse, b) {
                (true, b'*') => '#',
                (true, b) => char::from(b.to_ascii_lowercase()),
                (false, b) => char::from(b.to_ascii_uppercase()),
            });
        }
    }
    if del_len > 0 {
        write!(out, "{}", -del_len).expect("writing to a String cannot fail");
        // No reference: every deleted base is written as N.
        for _ in 0..del_len {
            out.push(if reverse { 'n' } else { 'N' });
        }
    }
    if entry.is_tail() {
        out.push('$');
    }
}

/// CIGAR operation codes (htslib's `BAM_C*`).
const BAM_CINS: u32 = 1;
const BAM_CDEL: u32 = 2;
const BAM_CPAD: u32 = 6;

/// `bam_plp_insertion_mod`: the inserted bases after this column (`*` for padding) and the
/// length of a deletion that directly follows them.
fn insertion(entry: &Entry<'_>) -> (Vec<u8>, i32) {
    let record = entry.record();
    let ops: Vec<(u32, u32)> = record.cigar().collect();
    let mut bases = Vec::new();
    let mut deletion = 0;
    let mut j = 1;
    let is_del = i32::from(entry.is_del());
    for &(op, len) in ops.iter().skip(entry.cigar_index() + 1) {
        match op {
            BAM_CPAD => bases.extend(std::iter::repeat_n(b'*', len as usize)),
            BAM_CINS => {
                for _ in 0..len {
                    let qpos = entry.qpos() + j - is_del;
                    bases.push(if qpos < record.len() { NT16[usize::from(record.base_code(qpos))] } else { b'N' });
                    j += 1;
                }
            }
            BAM_CDEL => {
                deletion = len as i32;
                break;
            }
            _ => break,
        }
    }
    (bases, deletion)
}

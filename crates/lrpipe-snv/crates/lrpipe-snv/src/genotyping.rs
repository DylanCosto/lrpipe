//! `postprocess/AddBackMissingVariantsInGenotyping.py`: with `--vcf_fn`, every site of the input
//! VCF ends up in `merge_output.vcf.gz`, with the original's call where there is one and `./.`
//! where there is not.

use std::collections::HashMap;
use std::fs;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};
use snv_merge::compress_and_index;

use crate::contigs::open_text;

/// `major_contigs_order`.
fn major_contigs() -> Vec<String> {
    let names: Vec<String> = (1..=22).map(|i| i.to_string()).chain(["X".to_owned(), "Y".to_owned()]).collect();
    names.iter().map(|n| format!("chr{n}")).chain(names.iter().cloned()).collect()
}

/// `VcfReader(keep_row_str=True, save_header=True).read_vcf()`: the header lines, and each
/// `(contig, position)`'s last row (with its newline) in first-seen order.
#[derive(Debug, Default)]
struct Rows {
    header: String,
    keys: Vec<(String, i64)>,
    rows: HashMap<(String, i64), String>,
}

impl Rows {
    fn read(path: &Path) -> Result<Self> {
        let mut out = Self::default();
        let mut reader = open_text(path)?;
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line).with_context(|| format!("reading {}", path.display()))? == 0 {
                break;
            }
            let mut columns = line.split_whitespace();
            let Some(contig) = columns.next() else { continue };
            if contig.starts_with('#') {
                out.header.push_str(&line);
                continue;
            }
            let pos: i64 = columns
                .next()
                .and_then(|p| p.parse().ok())
                .with_context(|| format!("bad VCF row in {}: {line}", path.display()))?;
            let key = (contig.to_owned(), pos);
            if !out.rows.contains_key(&key) {
                out.keys.push(key.clone());
            }
            out.rows.insert(key, line.clone());
        }
        Ok(out)
    }
}

/// A site the calls missed, as `--switch_genotype` writes it: the first REF base, no ALT, QUAL,
/// FILTER or INFO, and `GT` `./.`.
fn missing_row(row: &str) -> String {
    let mut columns: Vec<String> = row.trim_end().split('\t').map(str::to_owned).collect();
    if columns.len() < 10 {
        columns.resize(10, ".".to_owned());
    }
    columns[3] = columns[3].chars().next().map_or_else(|| ".".to_owned(), String::from);
    for (i, value) in [(4, "."), (5, "."), (6, "."), (7, "."), (8, "GT"), (9, "./.")] {
        columns[i] = value.to_owned();
    }
    columns.join("\t") + "\n"
}

/// Rewrites `output` (`merge_output.vcf`, then bgzipped and indexed) from the sites of `vcf`:
/// the row of `calls` at each site, or a `./.` row. Returns the original's summary message.
///
/// # Errors
///
/// Fails if a file cannot be read or written.
pub(crate) fn add_back_missing_variants(vcf: &Path, calls: &Path, output: &Path) -> Result<String> {
    let known = Rows::read(vcf)?;
    let called = Rows::read(calls)?;
    let major = major_contigs();
    let mut contigs: Vec<&str> = Vec::new();
    for (contig, _) in &known.keys {
        if !contigs.contains(&contig.as_str()) {
            contigs.push(contig);
        }
    }
    // sorted(contigs, key=(major_contigs_order + contigs).index); the other contigs keep their
    // first-seen order here, where Python's set order is arbitrary.
    let rank = |c: &str| major.iter().position(|m| m == c).unwrap_or(major.len());
    contigs.sort_by_key(|c| rank(c));

    let mut text = called.header.clone();
    let mut added = 0usize;
    for contig in contigs {
        let mut sites: Vec<(i64, String)> = Vec::new();
        for key in known.keys.iter().filter(|k| k.0 == contig) {
            let row = match called.rows.get(key) {
                Some(row) => row.clone(),
                None => {
                    added += 1;
                    missing_row(&known.rows[key])
                }
            };
            sites.push((key.1, row));
        }
        sites.sort_by_key(|(pos, _)| *pos);
        for (_, row) in sites {
            text.push_str(&row);
        }
    }
    fs::write(output, text).with_context(|| format!("writing {}", output.display()))?;
    compress_and_index(output)?;
    Ok(format!(
        "[INFO] Total variants for genotyping: {}, total variant output: {}, added {added} variants into output VCF",
        known.rows.len(),
        called.rows.len()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_sites_become_no_calls() {
        assert_eq!(missing_row("chr1\t10\trs1\tACG\tA\t50\tPASS\tDP=3\n"), "chr1\t10\trs1\tA\t.\t.\t.\t.\tGT\t./.\n");
        assert_eq!(
            missing_row("chr1\t10\t.\tA\tC\t5\tPASS\t.\tGT:GQ\t0/1:3\t1/1:4\n"),
            "chr1\t10\t.\tA\t.\t.\t.\t.\tGT\t./.\t1/1:4\n"
        );
    }
}

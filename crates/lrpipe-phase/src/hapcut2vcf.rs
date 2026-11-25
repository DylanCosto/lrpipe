//! Port of `cli/hapcut2vcf.py`: combine HapCUT/HapCUT2 output with the
//! original VCF and write a phased VCF.

use std::collections::HashMap;
use std::io::BufReader;

use anyhow::{Context, Result, bail};

use crate::compat::pystr::{py_int, py_isspace};
use crate::core::read::{Read, ReadSet};
use crate::vcf::PhasedVcfWriter;
use crate::xopen::TextLines;

/// A phased variant of a HapCUT block (`HapCutVariant`)
struct HapCutVariant {
    chromosome: String,
    position: i32,
    haplotype1: i32,
    haplotype2: i32,
    component_id: i32,
}

/// Python's `int(s)` converted to a C int (as `Read.add_variant` takes it)
fn c_int(s: &str) -> Result<i32> {
    let v = py_int(s).with_context(|| format!("ValueError: invalid literal for int() with base 10: '{s}'"))?;
    i32::try_from(v).map_err(|_| anyhow::anyhow!("OverflowError: value too large to convert to int"))
}

/// `BLOCK: offset: N len: N phased: N SPAN: N [MECscore N.N ]fragments N`
fn block_line_matches(line: &str) -> bool {
    fn number(s: &str, decimal: bool) -> Option<&str> {
        let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        if end == 0 {
            return None;
        }
        let rest = &s[end..];
        if !decimal {
            return Some(rest);
        }
        let rest = rest.strip_prefix('.')?;
        let end = rest.find(|c: char| !c.is_ascii_digit()).unwrap_or(rest.len());
        if end == 0 { None } else { Some(&rest[end..]) }
    }
    let mut s = line;
    for label in ["BLOCK: offset: ", " len: ", " phased: ", " SPAN: "] {
        let Some(rest) = s.strip_prefix(label) else { return false };
        let Some(rest) = number(rest, false) else { return false };
        s = rest;
    }
    let Some(mut s) = s.strip_prefix(' ') else { return false };
    if let Some(rest) = s.strip_prefix("MECscore ")
        && let Some(rest) = number(rest, true).and_then(|r| r.strip_prefix(' '))
    {
        s = rest;
    }
    s.strip_prefix("fragments ").and_then(|r| number(r, false)).is_some()
}

/// `HapCutParser.parse_blocks`: the connected components
fn parse_blocks(path: &str) -> Result<Vec<Vec<HapCutVariant>>> {
    let file = std::fs::File::open(path).with_context(|| format!("No such file or directory: '{path}'"))?;
    let mut blocks = Vec::new();
    let mut expect_block = true;
    let mut block: Vec<HapCutVariant> = Vec::new();
    for line in TextLines::new(BufReader::new(file), path) {
        let line = line?;
        if expect_block {
            expect_block = false;
            if !line.starts_with("BLOCK:") {
                bail!("ParseError: Expected a new block (line starting with \"BLOCK:\")");
            }
            if !block_line_matches(&line) {
                bail!("ParseError: BLOCK line malformed");
            }
        } else if line.starts_with("********") {
            if !block.is_empty() {
                blocks.push(std::mem::take(&mut block));
            }
            expect_block = true;
        } else {
            let fields: Vec<&str> = line.split(py_isspace).filter(|f| !f.is_empty()).collect();
            if fields.len() != 9 && fields.len() != 11 {
                bail!(
                    "ParseError: Expected nine fields (for hapCUT 1) or eleven fields (for hapCUT 2) in variant line"
                );
            }
            if fields.len() == 9 {
                // only checked for the right format
                let mut rest: Vec<&str> = fields[8].split(':').collect();
                if rest.len() == 5 {
                    if rest[4] != "FV" {
                        bail!("ParseError: Expected \"FV\" after last colon");
                    }
                    rest.pop();
                }
                if rest.len() != 4 {
                    bail!("ParseError: Too few elements in last (colon-separated) field");
                }
            }
            let (variant_id, haplotype1, haplotype2, chromosome, position) =
                (fields[0], fields[1], fields[2], fields[3], fields[4]);
            if haplotype1 == "-" || haplotype2 == "-" {
                // This happens in hapCUT 2 sometimes
                continue;
            }
            py_int(variant_id).with_context(|| format!("ValueError: invalid literal for int(): '{variant_id}'"))?;
            let haplotype1 = c_int(haplotype1)?;
            let haplotype2 = c_int(haplotype2)?;
            let position = c_int(position)?.checked_sub(1).context("OverflowError: position")?;
            let component_id = block.first().map_or(position, |v| v.position);
            block.push(HapCutVariant {
                chromosome: chromosome.to_string(),
                position,
                haplotype1,
                haplotype2,
                component_id,
            });
        }
    }
    if !block.is_empty() {
        blocks.push(block);
    }
    Ok(blocks)
}

/// `run_hapcut2vcf`; `output` "-" is standard output.
pub fn run_hapcut2vcf(hapcut: &str, vcf: &str, output: &str, command_line: &str) -> Result<()> {
    let mut writer = PhasedVcfWriter::new_plain_output(vcf, Some(command_line), output)?;
    if writer.samples.len() > 1 {
        bail!("There is more than one sample in this VCF");
    }
    let sample = writer.samples.first().context("IndexError: list index out of range")?.clone();
    let blocks = parse_blocks(hapcut)?;
    // blocks grouped by consecutive chromosome (itertools.groupby)
    let mut i = 0;
    while i < blocks.len() {
        let chromosome = blocks[i][0].chromosome.clone();
        let mut j = i;
        while j < blocks.len() && blocks[j][0].chromosome == chromosome {
            j += 1;
        }
        log::info!("Read {} phased blocks for chromosome {}", j - i, chromosome);
        // one read for each haplotype, and the connected components
        let mut haplotypes = [Read::new("1", 0, 0, 0), Read::new("2", 0, 0, 0)];
        let mut components: HashMap<i32, i32> = HashMap::new();
        for block in &blocks[i..j] {
            for v in block {
                haplotypes[0].add_variant(v.position, v.haplotype1, 0);
                haplotypes[1].add_variant(v.position, v.haplotype2, 0);
                components.insert(v.position, v.component_id);
            }
        }
        let mut superreads = ReadSet::new();
        let [h1, h2] = haplotypes;
        superreads.add(h1)?;
        superreads.add(h2)?;
        let sample_components = HashMap::from([(sample.clone(), components)]);
        writer.write(&chromosome, &[(sample.clone(), superreads)], &sample_components)?;
        i = j;
    }
    // WhatsHap never writes the records after the last phased chromosome
    writer.close()
}

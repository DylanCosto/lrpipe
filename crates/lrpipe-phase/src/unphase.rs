//! Port of `cli/unphase.py`: remove all phasing information from a VCF.
//!
//! The HP, PS and PQ tags are removed, and genotypes are written unphased
//! with the alleles sorted (`1|0` becomes `0/1`).

use anyhow::{Context, Result, bail};

use crate::hts::{BCF_HL_FMT, Header, Record, VcfFile, VcfWriter};

const TAGS_TO_REMOVE: [&str; 3] = ["HP", "PQ", "PS"];

/// `unphase_header`
fn unphase_header(header: &mut Header) {
    header.remove_first_record_with_key("phasing");
    for tag in TAGS_TO_REMOVE {
        if header.is_defined(BCF_HL_FMT, tag) {
            header.remove(BCF_HL_FMT, tag);
        }
    }
}

/// `run_unphase`: reads `vcf_path` ("-": standard input) and writes the
/// unphased VCF to standard output.
pub fn run_unphase(vcf_path: &str) -> Result<()> {
    let mut reader = VcfFile::open(vcf_path)?;
    unphase_header(&mut reader.header);
    let mut writer = VcfWriter::create("-", &reader.header)?;
    let mut record = Record::new();
    while reader.read(&mut record)? {
        let header = &reader.header;
        for tag in TAGS_TO_REMOVE {
            if record.has_format(header, tag)? {
                record.format_delete(header, tag)?;
            }
        }
        for sample in 0..header.nsamples() {
            // call["GT"] raises KeyError without a GT field
            let Some(gt) = record.gt(header, sample)? else { bail!("KeyError: 'invalid FORMAT: GT'") };
            // call["GT"][0] is not None and call["GT"][1] is not None
            if gt.first().context("IndexError: tuple index out of range")?.is_some()
                && gt.get(1).context("IndexError: tuple index out of range")?.is_some()
            {
                // sorted() of ints and None raises TypeError
                let mut alleles: Vec<u32> = gt
                    .iter()
                    .map(|a| a.context("TypeError: '<' not supported between instances of 'NoneType' and 'int'"))
                    .collect::<Result<_>>()?;
                alleles.sort_unstable();
                record.set_gt(header, sample, &alleles)?;
            }
            record.set_phased(header, sample, false)?;
        }
        writer.write(&mut record, &reader.header)?;
    }
    writer.close()
}

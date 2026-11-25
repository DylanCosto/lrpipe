//! Dumps what `snv-realign` builds, for comparison with the original's Python.
//!
//! - `realign_dump reads BAM REF CONTIG BED`: the SAM lines `RealignReads.py
//!   --full_aln_regions BED` writes after the header.
//! - `realign_dump tensors BAM REF CONTIG BED [ont|hifi|ilmn]`: `CreateTensorFullAlignment.py`'s
//!   rows for the same job (the sequence column holds only the centre base); Illumina by default.

use std::io::{BufWriter, Write};

use snv_mpileup::{AlignmentFile, Regions};
use snv_pileup::Fasta;
use snv_realign::{FullAlignmentJob, Platform, RealignInput, RealignedReads, full_alignment_candidates};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out = BufWriter::new(std::io::stdout().lock());
    match args.iter().map(String::as_str).collect::<Vec<_>>()[..] {
        ["reads", bam, reference, contig, bed] => {
            let mut view_file = AlignmentFile::open(bam, reference)?;
            let file = AlignmentFile::open(bam, reference)?;
            let fasta = Fasta::open(reference)?;
            let text = std::fs::read_to_string(bed)?;
            let (mut start, mut end) = (i64::MAX, 0);
            for row in text.lines() {
                let cols: Vec<&str> = row.split('\t').collect();
                if cols.first() == Some(&contig) {
                    start = start.min(cols[1].parse::<i64>()? + 1);
                    end = end.max(cols[2].parse::<i64>()? + 1);
                }
            }
            let input = RealignInput { contig, span: (start, end) };
            for line in RealignedReads::new(&mut view_file, &file, &fasta, input, Regions::read(bed.as_ref())?)? {
                writeln!(out, "{}", line?)?;
            }
        }
        ["tensors", bam, reference, contig, bed, ref platform @ ..] => {
            let platform = match platform.first().copied() {
                None | Some("ilmn") => Platform::Ilmn,
                Some("ont") => Platform::Ont,
                Some("hifi") => Platform::Hifi,
                Some(other) => return Err(format!("unknown platform {other}").into()),
            };
            let job = FullAlignmentJob {
                platform,
                bam: bam.as_ref(),
                reference: reference.as_ref(),
                contig,
                regions: bed.as_ref(),
                min_mapq: 5,
                min_coverage: 2,
                known_sites: None,
            };
            for c in full_alignment_candidates(&job)? {
                let mut fields = c.position.rsplitn(3, ':');
                let (center, pos, contig) = (fields.next(), fields.next(), fields.next());
                let tensor: Vec<String> = c.tensor.iter().map(i8::to_string).collect();
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    contig.unwrap_or_default(),
                    pos.unwrap_or_default(),
                    center.unwrap_or_default(),
                    tensor.join(" "),
                    c.alt_info
                )?;
            }
        }
        _ => {
            return Err(
                "usage: realign_dump reads BAM REF CONTIG BED | tensors BAM REF CONTIG BED [ont|hifi|ilmn]".into()
            );
        }
    }
    out.flush()?;
    Ok(())
}

//! Dumps what `snv-mpileup` builds, for comparison with samtools and the original's Python.
//!
//! - `mpileup_dump text BAM REF REGION [BED|- [-a]]`: columns 1-5 of `samtools mpileup BAM -r
//!   REGION --reverse-del --min-MQ 5 --min-BQ 0 --excl-flags 2316 [-l BED [-a]]`.
//! - `mpileup_dump pileup BAM REF CONTIG CHUNK_ID CHUNK_NUM`: `CreateTensorPileup`'s output
//!   rows (the sequence column holds only the centre base).

use std::io::{BufWriter, Write};

use snv_mpileup::{
    AlignmentFile, EXCLUDED_FLAGS, FastMode, PileupJob, Positions, ReadFilter, Regions, pileup_candidates,
    text::read_bases,
};
use snv_pileup::{CandidateAlleles, ChunkRegion, Fasta, FlankRequirement};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out = BufWriter::new(std::io::stdout().lock());
    match args.iter().map(String::as_str).collect::<Vec<_>>()[..] {
        ["text", bam, reference, region, ref rest @ ..] => {
            let mut bam = AlignmentFile::open(bam, reference)?;
            let filter = ReadFilter {
                excluded_flags: EXCLUDED_FLAGS,
                min_mapq: 5,
                bed: rest
                    .first()
                    .filter(|bed| **bed != "-")
                    .map(|bed| Regions::read(std::path::Path::new(bed)))
                    .transpose()?,
                positions: if rest.get(1) == Some(&"-a") { Positions::All } else { Positions::Covered },
            };
            let mut pileup = bam.pileup(region, filter)?;
            while let Some(column) = pileup.next_column()? {
                let depth = column.entries().len();
                writeln!(out, "{}\t{}\tN\t{depth}\t{}", column.contig, column.pos + 1, read_bases(&column))?;
            }
        }
        ["pileup", bam, reference, contig, id, count] => {
            let mut bam = AlignmentFile::open(bam, reference)?;
            let fasta = Fasta::open(reference)?;
            let job = PileupJob {
                contig,
                chunk: ChunkRegion::FaiChunk { id: id.parse::<i64>()? - 1, count: count.parse()? },
                extend_bed: None,
                confident_bed: None,
                min_snp_af: 0.08,
                min_indel_af: 0.08,
                min_coverage: 2,
                min_mapq: 5,
                alleles: CandidateAlleles::SnpsAndIndels,
                flanks: FlankRequirement::Contiguous,
                fast_mode: FastMode::Off,
            };
            for c in pileup_candidates(&mut bam, &fasta, &job)? {
                let mut fields = c.position.rsplitn(3, ':');
                let (center, pos, contig) = (fields.next(), fields.next(), fields.next());
                let tensor: Vec<String> = c.tensor.iter().map(i32::to_string).collect();
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
                "usage: mpileup_dump text BAM REF REGION [BED [-a]] | pileup BAM REF CONTIG CHUNK_ID CHUNK_NUM".into(),
            );
        }
    }
    out.flush()?;
    Ok(())
}

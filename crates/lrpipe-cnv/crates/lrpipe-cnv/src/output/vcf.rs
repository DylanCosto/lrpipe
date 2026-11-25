//! `{sample}.vcf.gz` (outputWriter.py:69-242), as it reads after Spectre's `bcftools sort`:
//! records in reference order, INFO without the trailing `;`, floats as htslib prints them.

use std::fmt::Write as _;
use std::path::Path;

use pynum::{repr::py_repr, round::np_round};

use super::htsfloat::vcf_float;
use super::tabix::{self, Line, Preset};
use crate::Error;
use crate::candidate::{Arena, CnvType};
use crate::genome::Genome;
use crate::merged::Call;

/// The header Spectre writes (outputWriter.py:96-137), with a contig line per kept contig.
pub fn header(genome: &Genome, sample: &str) -> String {
    let mut h =
        String::from("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##source=Spectre\n");
    for c in &genome.contigs {
        let _ = writeln!(h, "##contig=<ID={},length={}>", c.name, c.len);
    }
    h.push_str(concat!(
        "##ALT=<ID=DEL,Description=\"Deletion\">\n",
        "##ALT=<ID=DUP,Description=\"Duplications\">\n",
        "##ALT=<ID=LOH,Description=\"Loss of Heterozygosity,\">\n",
        "##FILTER=<ID=UNRESOLVED,Description=\"An insertion that is longer than the read and thus we cannot predict the full size.\">\n",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of copy number variant\">\n",
        "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Estimated copy number status\">\n",
        "##INFO=<ID=SVSUPPORT,Number=1,Type=String,Description=\"Indicator if a SV support was found in a provided SNFJ file\">\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##FORMAT=<ID=HO,Number=2,Type=Float,Description=\"Homozygosity proportion\">\n",
        "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">\n",
        "##FORMAT=<ID=CD,Number=2,Type=Float,Description=\"Coverage read depth\">\n",
        "##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Population ID of supporting CNV calls\">\n",
    ));
    let _ = writeln!(h, "##Spectre_sample={sample}");
    let _ = writeln!(h, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}");
    h
}

/// How the VCF is written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Style {
    /// As Spectre writes it, byte for byte: END is a CNV's last window start, and the header
    /// declares two values for `HO` and `CD`, of which records hold one.
    #[default]
    Spectre,
    /// Valid VCF 4.2 that strict readers accept: `HO` and `CD` declared with one value, END
    /// the event's last base, at most the contig's last (a contig's last window is cut short
    /// by its end; POS stays the base before it, as for symbolic alleles), SVLEN
    /// negative for deletions, an LoH call's POS the base before its first SNV, and
    /// the actual reference anchor, the CNV score as `CNV_GQ`, and `##source=lrpipe-cnv`.
    Standard,
}

/// Writes the calls. For a CNV, POS is the first window's 0-based start and END the last
/// window's start, as Spectre writes them (or, with SV breakpoints, the SV's POS and END and
/// an `SVID`); GQ is the z-test score, CD `median_raw_cov` rounded
/// to 2 decimals. For LoH, POS and END are the first and last SNV positions. [`Style::Standard`]
/// changes this as it says.
///
/// # Errors
///
/// If the file or its index cannot be written, or a CNV call has no quality score.
pub fn write(
    path: &Path,
    reference: &Path,
    genome: &Genome,
    sample: &str,
    arena: &Arena,
    calls: &[(String, Vec<Call>)],
    style: Style,
) -> Result<(), Error> {
    let standard = style == Style::Standard;
    let mut fasta = if standard {
        Some(lrpipe_reads::indexed_fasta(reference).map_err(|e| Error::Input(format!("opening reference: {e}")))?)
    } else {
        None
    };
    let mut anchor = |contig: &str, pos: i64| -> Result<char, Error> {
        let Some(fasta) = &mut fasta else { return Ok('N') };
        use pasta::core::{Position, Region};
        let pos = usize::try_from(pos)
            .ok()
            .and_then(Position::new)
            .ok_or_else(|| Error::Input(format!("invalid CNV anchor {contig}:{pos}")))?;
        let region = Region::new(contig.as_bytes(), pos..=pos);
        let record = fasta.query(&region).map_err(Error::io(reference))?;
        let bases: &[u8] = record.sequence().as_ref();
        let base = bases.first().ok_or_else(|| Error::Input(format!("reference has no base at {region}")))?;
        Ok(char::from(lrpipe_reads::vcf_reference_base(*base)))
    };
    let score_key = if standard { "CNV_GQ" } else { "GQ" };
    let mut records = Vec::new();
    for (_, list) in calls {
        for call in list {
            let (chrom, start, end, text, alt) = match call {
                Call::Cnv(id) => {
                    let c = arena.get(*id);
                    let score =
                        c.statistics.as_ref().ok_or_else(|| Error::Input(format!("{} has no quality score", c.id)))?;
                    let alt = format!("<{}>", c.cnv_type.as_str());
                    let format = [
                        c.gt.to_owned(),
                        vcf_float(&py_repr(0.0)),
                        vcf_float(&score.sample_score.to_string()),
                        vcf_float(&py_repr(np_round(c.median_raw_cov, 2))),
                    ]
                    .join(":");
                    // With SV breakpoints: the SV's POS and END, and its ID.
                    let (start, end, size, svid) = match &c.breakpoints {
                        Some(b) => (b.start, b.end, b.end - b.start, format!(";SVID={}", b.sv_id)),
                        None if standard => {
                            // The last window ends at the contig's end, if that comes first.
                            let contig = genome.contigs.iter().find(|g| g.name == c.chromosome);
                            let last = contig.map_or(i64::MAX, |g| i64::try_from(g.len).unwrap_or(i64::MAX));
                            let end = (c.end + c.bin_size).min(last);
                            (c.start.max(1), end, end - c.start.max(1), String::new())
                        }
                        None => (c.start, c.end, c.size, String::new()),
                    };
                    let size = if standard && c.cnv_type == CnvType::Del { -size } else { size };
                    let reference_base = anchor(&c.chromosome, start)?;
                    let text = format!(
                        "{}\t{start}\t{}\t{reference_base}\t{alt}\t.\t.\tEND={end};SVLEN={size};SVTYPE={};CN={}{svid}\tGT:HO:{score_key}:CD\t{format}",
                        c.chromosome,
                        c.id,
                        c.cnv_type.as_str(),
                        c.cn_status
                    );
                    (&c.chromosome, start, end, text, alt)
                }
                Call::Loh(l) => {
                    // A passing LoH call: QUAL the mean SNV QUAL, HO the mean AF, GQ the size score.
                    let qual = l.quality.map_or_else(|| ".".to_owned(), |q| vcf_float(&py_repr(q)));
                    let gq = l.score.map_or(0, |s| s as i64);
                    let format = [
                        "1/1".to_owned(),
                        vcf_float(&py_repr(np_round(l.het_score, 2))),
                        vcf_float(&gq.to_string()),
                        if standard { ".".to_owned() } else { vcf_float("nan") },
                    ]
                    .join(":");
                    let (start, size) =
                        if standard { ((l.start - 1).max(1), l.end - l.start + 1) } else { (l.start, l.size) };
                    let reference_base = anchor(&l.chromosome, start)?;
                    let text = format!(
                        "{}\t{start}\t{}\t{reference_base}\t<LOH>\t{qual}\t{}\tEND={};SVLEN={size};SVTYPE=LOH;CN={}\tGT:HO:{score_key}:CD\t{format}",
                        l.chromosome, l.id, l.filter, l.end, l.cn_status
                    );
                    (&l.chromosome, start, l.end, text, "<LOH>".to_owned())
                }
            };
            let rid = genome.contigs.iter().position(|g| &g.name == chrom).unwrap_or(usize::MAX);
            records.push((rid, start, alt, Line { contig: chrom.clone(), start, end, text }));
        }
    }
    // bcftools sort: by contig, position, then alleles; REF is the same at a position.
    records.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)).then_with(|| a.2.cmp(&b.2)));
    let lines: Vec<Line> = records.into_iter().map(|r| r.3).collect();
    let mut header = header(genome, sample);
    let snapped =
        calls.iter().flat_map(|(_, l)| l).any(|c| matches!(c, Call::Cnv(id) if arena.get(*id).breakpoints.is_some()));
    if snapped {
        header = header.replacen(
            "##FORMAT=<ID=GT,",
            "##INFO=<ID=SVID,Number=1,Type=String,Description=\"ID of the SV call giving the breakpoints\">\n##FORMAT=<ID=GT,",
            1,
        );
    }
    if standard {
        header = header
            .replacen("##source=Spectre", "##source=lrpipe-cnv", 1)
            .replacen(
                "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype quality\">",
                "##FORMAT=<ID=CNV_GQ,Number=1,Type=Float,Description=\"Copy-number or LOH score\">",
                1,
            )
            .replacen("##FORMAT=<ID=HO,Number=2", "##FORMAT=<ID=HO,Number=1", 1)
            .replacen("##FORMAT=<ID=CD,Number=2", "##FORMAT=<ID=CD,Number=1", 1);
    }
    tabix::write(path, &header, &lines, Preset::Vcf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::Candidate;
    use crate::genome::Contig;
    use crate::score::ZScore;

    fn records(path: &Path) -> Vec<String> {
        let mut text = String::new();
        std::io::Read::read_to_string(
            &mut flate2::read::MultiGzDecoder::new(std::fs::File::open(path).expect("written")),
            &mut text,
        )
        .expect("gzip");
        text.lines().filter(|l| !l.starts_with('#')).map(str::to_owned).collect()
    }

    #[test]
    fn a_standard_end_stops_at_the_contigs_end() {
        // lrpipe run on HG003 chr20 wrote END=64445000 for a deletion to chr20's end, whose
        // length is 64444167: its last window, at 64444000, is cut short by the contig's end.
        let dir = std::env::temp_dir().join(format!("lrpipe-cnv-vcf-end-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let genome =
            Genome { contigs: vec![Contig { name: "chr20".into(), len: 4167, ploidy_multiplier: 2 }], ploidy: 2 };
        let reference = dir.join("ref.fa");
        std::fs::write(&reference, format!(">chr20\n{}\n", "G".repeat(4167))).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "chr20\t4167\t7\t4167\t4168\n").unwrap();
        let mut arena = Arena::default();
        let mut c = Candidate::new("chr20", CnvType::Del, 1000, vec![1000, 4000], vec![0.0, 0.0]);
        (c.start, c.end, c.size, c.id) = (1000, 4000, 3000, "Spectre.DEL.X".into());
        c.statistics = Some(ZScore::unrefined());
        let id = arena.push(c);
        let calls = vec![("chr20".to_owned(), vec![Call::Cnv(id)])];
        let path = dir.join("s.vcf.gz");
        write(&path, &reference, &genome, "S", &arena, &calls, Style::Standard).expect("written");
        let record = records(&path).remove(0);
        assert!(record.contains("\tG\t<DEL>\t"), "{record}");
        assert!(record.contains("\tEND=4167;SVLEN=-3167;"), "{record}");
        // Spectre's own style writes the last window's start, as Spectre does.
        write(&path, &reference, &genome, "S", &arena, &calls, Style::Spectre).expect("written");
        assert!(records(&path)[0].contains("\tEND=4000;SVLEN=3000;"), "{}", records(&path)[0]);
        let _ = std::fs::remove_dir_all(dir);
    }
}

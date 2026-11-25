//! VCF output: one record per locus. REF is the reference repeat with the base before it;
//! ALT the distinct non-reference alleles likewise. FORMAT has the genotype, allele lengths
//! (`AL`, repeat bases only), the range of supporting reads' lengths (`ALLR`), supporting
//! reads (`SD`) and motif counts (`MC`, one count per motif joined by `_`).

use std::fmt::Write as _;

use crate::catalog::Locus;
use crate::genotype::{Allele, motif_counts};

/// The header, with a `##contig` line per reference contig.
pub fn header(contigs: &[(String, u64)], sample: &str) -> String {
    let mut h = String::from("##fileformat=VCFv4.2\n##source=lrpipe-tr\n");
    for (name, len) in contigs {
        let _ = writeln!(h, "##contig=<ID={name},length={len}>");
    }
    h.push_str(concat!(
        "##INFO=<ID=TRID,Number=1,Type=String,Description=\"Locus ID\">\n",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the repeat\">\n",
        "##INFO=<ID=MOTIFS,Number=.,Type=String,Description=\"Repeat units\">\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##FORMAT=<ID=AL,Number=.,Type=Integer,Description=\"Allele lengths (repeat bases)\">\n",
        "##FORMAT=<ID=ALLR,Number=.,Type=String,Description=\"Shortest-longest supporting read per allele\">\n",
        "##FORMAT=<ID=SD,Number=.,Type=Integer,Description=\"Reads supporting each allele\">\n",
        "##FORMAT=<ID=MC,Number=.,Type=String,Description=\"Copies of each motif per allele, joined by _\">\n",
    ));
    let _ = writeln!(h, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}");
    h
}

/// One locus's record. `anchor` is the reference base before the repeat, `reference` the
/// repeat's reference bases; `copies` the expected allele count (for a no-call's GT).
pub fn record(locus: &Locus, anchor: u8, reference: &[u8], alleles: &[Allele], copies: usize) -> String {
    let anchor = lrpipe_reads::vcf_reference_base(anchor);
    let reference: Vec<u8> = reference.iter().copied().map(lrpipe_reads::vcf_reference_base).collect();
    let mut alts = Vec::new();
    let mut gt = Vec::new();
    for a in alleles {
        let sequence: Vec<u8> = a
            .seq
            .iter()
            .map(|base| match base.to_ascii_uppercase() {
                base @ (b'A' | b'C' | b'G' | b'T') => base,
                _ => b'N',
            })
            .collect();
        if sequence == reference {
            gt.push("0".to_owned());
        } else {
            let i = alts.iter().position(|s| *s == sequence).unwrap_or_else(|| {
                alts.push(sequence);
                alts.len() - 1
            });
            gt.push((i + 1).to_string());
        }
    }
    let with_anchor = |s: &[u8]| {
        let mut v = String::with_capacity(s.len() + 1);
        v.push(char::from(anchor));
        v.push_str(&String::from_utf8_lossy(s));
        v
    };
    let alt = if alts.is_empty() {
        ".".to_owned()
    } else {
        alts.iter().map(|s| with_anchor(s)).collect::<Vec<_>>().join(",")
    };
    let motifs: Vec<String> = locus.motifs.iter().map(|m| String::from_utf8_lossy(m).into_owned()).collect();
    let format = if alleles.is_empty() {
        vec![
            ".".repeat(copies).chars().map(String::from).collect::<Vec<_>>().join("/"),
            ".".into(),
            ".".into(),
            ".".into(),
            ".".into(),
        ]
    } else {
        let join = |f: &dyn Fn(&Allele) -> String| alleles.iter().map(f).collect::<Vec<_>>().join(",");
        vec![
            gt.join("/"),
            join(&|a| a.seq.len().to_string()),
            join(&|a| format!("{}-{}", a.min_len, a.max_len)),
            join(&|a| a.support.to_string()),
            join(&|a| motif_counts(&a.seq, &locus.motifs).iter().map(usize::to_string).collect::<Vec<_>>().join("_")),
        ]
    };
    format!(
        "{}\t{}\t.\t{}\t{alt}\t.\t.\tTRID={};END={};MOTIFS={}\tGT:AL:ALLR:SD:MC\t{}\n",
        locus.contig,
        locus.start,
        with_anchor(&reference),
        locus.id,
        locus.end,
        motifs.join(","),
        format.join(":")
    )
}

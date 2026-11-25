//! CNV and LoH calls together, as the BED and VCF list them (`MergeCNVLoH`, spec §3.14).

use crate::candidate::{Arena, CandidateId};
use crate::loh::{LohCalls, LohCandidate};

/// One output record.
#[derive(Clone, Debug)]
pub enum Call {
    Cnv(CandidateId),
    Loh(Box<LohCandidate>),
}

impl Call {
    /// The start the BED is sorted by: a 0-based window start or a 1-based SNV position.
    pub fn start(&self, arena: &Arena) -> i64 {
        match self {
            Call::Cnv(id) => arena.get(*id).start,
            Call::Loh(l) => l.start,
        }
    }
}

/// Per analysed contig, CNV calls and passing LoH calls interleaved by start (a CNV first on
/// ties). Without SNV data, the CNV calls alone.
pub fn merge(arena: &Arena, cnv: &[(String, Vec<CandidateId>)], loh: Option<&LohCalls>) -> Vec<(String, Vec<Call>)> {
    // With no analysed contig but SNV data, Spectre outputs the LoH calls alone.
    if let (true, Some(l)) = (cnv.is_empty(), loh) {
        return l
            .passing
            .iter()
            .map(|(c, v)| (c.clone(), v.iter().map(|l| Call::Loh(Box::new(l.clone()))).collect()))
            .collect();
    }
    cnv.iter()
        .map(|(contig, cnvs)| {
            let lohs: &[LohCandidate] =
                loh.and_then(|l| l.passing.iter().find(|(c, _)| c == contig)).map_or(&[], |(_, v)| v.as_slice());
            let (mut i, mut j) = (0, 0);
            let mut out = Vec::with_capacity(cnvs.len() + lohs.len());
            while i < cnvs.len() || j < lohs.len() {
                let take_cnv = j >= lohs.len() || (i < cnvs.len() && lohs[j].start >= arena.get(cnvs[i]).start);
                if take_cnv {
                    out.push(Call::Cnv(cnvs[i]));
                    i += 1;
                } else {
                    out.push(Call::Loh(Box::new(lohs[j].clone())));
                    j += 1;
                }
            }
            (contig.clone(), out)
        })
        .collect()
}

//! Port of `possible_outcome_probabilites_from`, `possible_outcome_probabilites_with_indel_length_from`
//! and `output_from` (`CallVariants.py`), with and without indel-length heads.
//!
//! The probabilities are numpy `float32` scalars in the original, so every product here is computed
//! in `f32`, and equality tests (`maximum_probability in list`) compare those exact values.

use crate::alt_info::AltInfo;
use crate::variant_length;

/// Index of each GT21 label (`GT21_Type`).
mod gt21 {
    pub(super) const AA: usize = 0;
    pub(super) const AC: usize = 1;
    pub(super) const AG: usize = 2;
    pub(super) const AT: usize = 3;
    pub(super) const CC: usize = 4;
    pub(super) const CG: usize = 5;
    pub(super) const CT: usize = 6;
    pub(super) const GG: usize = 7;
    pub(super) const GT: usize = 8;
    pub(super) const TT: usize = 9;
    pub(super) const DEL_DEL: usize = 10;
    pub(super) const A_DEL: usize = 11;
    pub(super) const INS_INS: usize = 15;
    pub(super) const A_INS: usize = 16;
    pub(super) const INS_DEL: usize = 20;
}

/// Number of GT21 classes.
pub const GT21_CLASSES: usize = 21;
/// Number of zygosity classes (0/0, 1/1, 0/1).
pub const GENOTYPE_CLASSES: usize = 3;
/// Number of classes of each of the two indel-length heads (lengths -16 to 16).
pub const VARIANT_LENGTH_CLASSES: usize = variant_length::CLASSES;

const HOMO_SNP_GT21: [usize; 4] = [gt21::AA, gt21::CC, gt21::GG, gt21::TT];
const HOMO_SNP_LABELS: [[u8; 2]; 4] = [*b"AA", *b"CC", *b"GG", *b"TT"];
const HETERO_SNP_GT21: [usize; 6] = [gt21::AC, gt21::AG, gt21::AT, gt21::CG, gt21::CT, gt21::GT];
const HETERO_SNP_LABELS: [[u8; 2]; 6] = [*b"AC", *b"AG", *b"AT", *b"CG", *b"CT", *b"GT"];
const ACGT: &[u8; 4] = b"ACGT";

/// The network outputs of one candidate.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Heads<'a> {
    /// GT21 softmax.
    pub(crate) gt21: &'a [f32],
    /// Zygosity softmax (0/0, 1/1, 0/1).
    pub(crate) genotype: &'a [f32],
    /// The two indel-length softmaxes (`add_indel_length=True`), or `None`.
    pub(crate) variant_lengths: Option<[&'a [f32]; 2]>,
}

/// Which outcome lists contain the winning probability (all of them can be true at once).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[expect(clippy::struct_excessive_bools, reason = "mirrors the ten flags the original returns and tests independently")]
pub(crate) struct Outcome {
    pub(crate) reference: bool,
    pub(crate) homo_snp: bool,
    pub(crate) hetero_snp: bool,
    pub(crate) homo_ins: bool,
    pub(crate) hetero_acgt_ins: bool,
    pub(crate) hetero_ins_ins: bool,
    pub(crate) homo_del: bool,
    pub(crate) hetero_acgt_del: bool,
    pub(crate) hetero_del_del: bool,
    pub(crate) ins_del: bool,
}

/// The decoded call: outcome flags, REF and ALT, and the winning probability.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Decoded {
    pub(crate) outcome: Outcome,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) probability: f32,
}

/// The tuple `possible_outcome_probabilites_from` returns: one probability list per outcome
/// type, with the indel lengths (and bases) each entry stands for. Without indel-length heads
/// the length lists are empty, as in the original, and the ACGT bases are `ACGT[idx]`.
#[derive(Debug, Clone, Default, PartialEq)]
struct Outcomes {
    homo_ref: f32,
    homo_snp: Vec<f32>,
    hetero_snp: Vec<f32>,
    homo_ins: Vec<f32>,
    homo_ins_lengths: Vec<usize>,
    hetero_ins_ins: Vec<f32>,
    hetero_ins_ins_lengths: Vec<(usize, usize)>,
    hetero_acgt_ins: Vec<f32>,
    hetero_acgt_ins_bases: Vec<u8>,
    hetero_acgt_ins_lengths: Vec<usize>,
    homo_del: Vec<f32>,
    homo_del_lengths: Vec<usize>,
    hetero_del_del: Vec<f32>,
    /// Ascending, as the original stores them.
    hetero_del_del_lengths: Vec<(usize, usize)>,
    hetero_acgt_del: Vec<f32>,
    hetero_acgt_del_bases: Vec<u8>,
    hetero_acgt_del_lengths: Vec<usize>,
    hetero_ins_del: Vec<f32>,
    /// `(deletion length, insertion length)`.
    hetero_ins_del_lengths: Vec<(usize, usize)>,
}

/// What `possible_outcome_probabilites_from` returns.
#[derive(Debug, Clone, PartialEq)]
enum Possible {
    /// `[homo_Ref_probability]`: a confident reference call, decided without the other outcomes.
    Reference(f32),
    /// All outcome lists.
    Outcomes(Box<Outcomes>),
}

/// `reference_base_ACGT`'s GT21 class (`gt21_enum_from_label(ref + ref)`).
fn reference_gt21(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(gt21::AA),
        b'C' => Some(gt21::CC),
        b'G' => Some(gt21::GG),
        b'T' => Some(gt21::TT),
        _ => None,
    }
}

/// `IUPAC_base_to_ACGT_base_dict`.
pub(crate) fn iupac_to_acgt(base: u8) -> Option<u8> {
    let i = b"ACGTURYSWKMBDHVN".iter().position(|&b| b == base)?;
    Some(b"ACGTTACCAGACAAAA"[i])
}

/// `max(lst)` of a non-empty list.
fn list_max(v: &[f32]) -> f32 {
    v.iter().copied().fold(f32::NEG_INFINITY, |m, x| if x > m { x } else { m })
}

/// `lst.index(x)`.
fn index_of(v: &[f32], x: f32) -> usize {
    v.iter().position(|&p| p == x).expect("the maximum was found in this list")
}

/// `np.argmax(lst)`.
fn argmax(v: &[f32]) -> usize {
    index_of(v, list_max(v))
}

/// The `[(GT21_Type.AIns, "A"), ...]` loop: `p * gt21[base] * hetero` for every length and base.
fn hetero_acgt(by_length: &[(usize, f32)], gt21: &[f32], first: usize, hetero: f32) -> (Vec<u8>, Vec<usize>, Vec<f32>) {
    let mut bases = Vec::with_capacity(4 * by_length.len());
    let mut lengths = Vec::with_capacity(4 * by_length.len());
    let mut probabilities = Vec::with_capacity(4 * by_length.len());
    for &(length, p) in by_length {
        for (b, &base) in ACGT.iter().enumerate() {
            bases.push(base);
            lengths.push(length);
            probabilities.push(p * gt21[first + b] * hetero);
        }
    }
    (bases, lengths, probabilities)
}

/// `possible_outcome_probabilites_with_indel_length_from`: every outcome list, weighted by the
/// two indel-length heads `l1` and `l2`.
fn outcomes_with_indel_length(gt21: &[f32], genotype: &[f32], l1: &[f32], l2: &[f32], ref_gt21: usize) -> Outcomes {
    let (homo_ref_p, homo_var_p, hetero_p) = (genotype[0], genotype[1], genotype[2]);
    let length_0 = l1[variant_length::INDEX_OFFSET] * l2[variant_length::INDEX_OFFSET];
    let (homo_ins_lengths, homo_ins) =
        variant_length::homo_ins(l1, l2, homo_var_p * gt21[gt21::INS_INS]).into_iter().unzip();
    let (hetero_ins_ins_lengths, hetero_ins_ins) =
        variant_length::hetero_ins_ins(l1, l2, hetero_p * gt21[gt21::INS_INS]).into_iter().unzip();
    let (hetero_acgt_ins_bases, hetero_acgt_ins_lengths, hetero_acgt_ins) =
        hetero_acgt(&variant_length::hetero_ins(l1, l2), gt21, gt21::A_INS, hetero_p);
    let (homo_del_lengths, homo_del) =
        variant_length::homo_del(l1, l2, homo_var_p * gt21[gt21::DEL_DEL]).into_iter().unzip();
    let (hetero_del_del_lengths, hetero_del_del) =
        variant_length::hetero_del_del(l1, l2, hetero_p * gt21[gt21::DEL_DEL]).into_iter().unzip();
    let (hetero_acgt_del_bases, hetero_acgt_del_lengths, hetero_acgt_del) =
        hetero_acgt(&variant_length::hetero_del(l1, l2), gt21, gt21::A_DEL, hetero_p);
    let (hetero_ins_del_lengths, hetero_ins_del) =
        variant_length::hetero_ins_del(l1, l2, hetero_p * gt21[gt21::INS_DEL]).into_iter().unzip();
    Outcomes {
        homo_ref: length_0 * homo_ref_p * gt21[ref_gt21],
        homo_snp: HOMO_SNP_GT21.iter().map(|&g| length_0 * homo_var_p * gt21[g]).collect(),
        hetero_snp: HETERO_SNP_GT21.iter().map(|&g| length_0 * hetero_p * gt21[g]).collect(),
        homo_ins,
        homo_ins_lengths,
        hetero_ins_ins,
        hetero_ins_ins_lengths,
        hetero_acgt_ins,
        hetero_acgt_ins_bases,
        hetero_acgt_ins_lengths,
        homo_del,
        homo_del_lengths,
        hetero_del_del,
        hetero_del_del_lengths,
        hetero_acgt_del,
        hetero_acgt_del_bases,
        hetero_acgt_del_lengths,
        hetero_ins_del,
        hetero_ins_del_lengths,
    }
}

/// `possible_outcome_probabilites_from`: with indel-length heads it returns early only when
/// both heads also give length 0 at least 0.5.
fn possible_outcomes(heads: &Heads<'_>, ref_gt21: usize) -> Possible {
    let Heads { gt21, genotype, variant_lengths } = *heads;
    let (homo_ref_p, homo_var_p, hetero_p) = (genotype[0], genotype[1], genotype[2]);
    if let Some([l1, l2]) = variant_lengths {
        let (zero_1, zero_2) = (l1[variant_length::INDEX_OFFSET], l2[variant_length::INDEX_OFFSET]);
        if zero_1 >= 0.5 && zero_2 >= 0.5 && homo_ref_p >= 0.5 && gt21[ref_gt21] >= 0.5 {
            return Possible::Reference(zero_1 * zero_2 * homo_ref_p * gt21[ref_gt21]);
        }
        return Possible::Outcomes(Box::new(outcomes_with_indel_length(gt21, genotype, l1, l2, ref_gt21)));
    }
    let homo_ref = homo_ref_p * gt21[ref_gt21];
    if homo_ref_p >= 0.5 && gt21[ref_gt21] >= 0.5 {
        return Possible::Reference(homo_ref);
    }
    Possible::Outcomes(Box::new(Outcomes {
        homo_ref,
        homo_snp: HOMO_SNP_GT21.iter().map(|&g| homo_var_p * gt21[g]).collect(),
        hetero_snp: HETERO_SNP_GT21.iter().map(|&g| hetero_p * gt21[g]).collect(),
        homo_ins: vec![homo_var_p * gt21[gt21::INS_INS]],
        hetero_ins_ins: vec![hetero_p * gt21[gt21::INS_INS]],
        hetero_acgt_ins: (0..4).map(|i| gt21[gt21::A_INS + i] * hetero_p).collect(),
        hetero_acgt_ins_bases: ACGT.to_vec(),
        homo_del: vec![homo_var_p * gt21[gt21::DEL_DEL]],
        hetero_del_del: vec![hetero_p * gt21[gt21::DEL_DEL]],
        hetero_acgt_del: (0..4).map(|i| gt21[gt21::A_DEL + i] * hetero_p).collect(),
        hetero_acgt_del_bases: ACGT.to_vec(),
        hetero_ins_del: vec![hetero_p * gt21[gt21::INS_DEL]],
        ..Outcomes::default()
    }))
}

/// Decodes one candidate (`output_from`), or `None` where the original would raise (a non-IUPAC
/// reference base).
///
/// `max_indel_len` is `maximum_variant_length_that_need_infer`.
pub(crate) fn decode(
    reference_sequence: &str,
    heads: &Heads<'_>,
    alt: &AltInfo,
    max_indel_len: usize,
) -> Option<Decoded> {
    // `tensor_position_center` is 0 for the one-base reference the C code writes.
    let center = *reference_sequence.as_bytes().first()?;
    let ref_acgt = iupac_to_acgt(center)?;
    let reference_call = |probability| {
        let base = char::from(ref_acgt).to_string();
        Decoded {
            outcome: Outcome { reference: true, ..Outcome::default() },
            reference: base.clone(),
            alternate: base,
            probability,
        }
    };
    let mut o = match possible_outcomes(heads, reference_gt21(ref_acgt)?) {
        Possible::Reference(p) => return Some(reference_call(p)),
        Possible::Outcomes(o) => o,
    };

    let ref_char = char::from(center);
    let ref_str = ref_char.to_string();
    let insertion = |length: Option<usize>, ignore: &str| {
        alt.insertion_bases(variant_length::proposed(length), ignore, max_indel_len)
    };
    let deletion = |length: Option<usize>, ignore: &str| {
        alt.deletion_bases(variant_length::proposed(length), ignore, max_indel_len)
    };
    // `reference_base` and `alternate_base` live across iterations, and every `continue` goes
    // back through `while reference_base is None or alternate_base is None`, so a branch that
    // gives up after setting both still ends the loop with them.
    let mut reference: Option<String> = None;
    let mut alternate: Option<String> = None;
    loop {
        let max = [
            o.homo_ref,
            list_max(&o.homo_snp),
            list_max(&o.hetero_snp),
            list_max(&o.homo_ins),
            list_max(&o.homo_del),
            list_max(&o.hetero_acgt_ins),
            list_max(&o.hetero_ins_ins),
            list_max(&o.hetero_acgt_del),
            list_max(&o.hetero_del_del),
            list_max(&o.hetero_ins_del),
        ]
        .into_iter()
        .fold(f32::NEG_INFINITY, |m, x| if x > m { x } else { m });

        if max == o.homo_ref {
            return Some(reference_call(max));
        }
        let outcome = Outcome {
            reference: false,
            homo_snp: o.homo_snp.contains(&max),
            hetero_snp: o.hetero_snp.contains(&max),
            homo_ins: o.homo_ins.contains(&max),
            hetero_acgt_ins: o.hetero_acgt_ins.contains(&max),
            hetero_ins_ins: o.hetero_ins_ins.contains(&max),
            homo_del: o.homo_del.contains(&max),
            hetero_acgt_del: o.hetero_acgt_del.contains(&max),
            hetero_del_del: o.hetero_del_del.contains(&max),
            ins_del: o.hetero_ins_del.contains(&max),
        };

        if outcome.homo_snp {
            reference = Some(ref_str.clone());
            let idx = index_of(&o.homo_snp, max);
            let [b1, b2] = HOMO_SNP_LABELS[argmax(&o.homo_snp)].map(char::from);
            let proposed = if b1 != ref_char { b1 } else { b2 }.to_string();
            alternate = alt.find_alt_base(Some(&proposed)).1;
            if alternate.as_deref().is_none_or(|a| a == ref_str) {
                o.homo_snp[idx] = 0.0;
            }
        } else if outcome.hetero_snp {
            let [b1, b2] = HETERO_SNP_LABELS[argmax(&o.hetero_snp)].map(char::from);
            let idx = index_of(&o.hetero_snp, max);
            reference = Some(ref_str.clone());
            if b1 != ref_char && b2 != ref_char {
                let (sorted, _) = alt.find_alt_base(None);
                if sorted.len() < 2 {
                    o.hetero_snp[idx] = 0.0;
                } else {
                    alternate = Some(sorted[..2].join(","));
                }
            } else {
                let proposed = if b1 != ref_char { b1 } else { b2 }.to_string();
                alternate = alt.find_alt_base(Some(&proposed)).1;
                if alternate.as_deref().is_none_or(|a| a == ref_str) {
                    o.hetero_snp[idx] = 0.0;
                }
            }
        } else if outcome.homo_ins {
            let idx = index_of(&o.homo_ins, max);
            match insertion(o.homo_ins_lengths.get(idx).copied(), "") {
                Some(ins) => {
                    reference = Some(ref_str.clone());
                    alternate = Some(ins.to_owned());
                }
                None => o.homo_ins[idx] = 0.0,
            }
        } else if outcome.hetero_acgt_ins {
            let idx = index_of(&o.hetero_acgt_ins, max);
            let hetero_base = char::from(o.hetero_acgt_ins_bases[idx]);
            match insertion(o.hetero_acgt_ins_lengths.get(idx).copied(), "") {
                Some(ins) => {
                    reference = Some(ref_str.clone());
                    alternate = Some(ins.to_owned());
                    if hetero_base != ref_char {
                        match alt.find_alt_base(None).0.first() {
                            Some(snp) => alternate = Some(format!("{snp},{ins}")),
                            None => o.hetero_acgt_ins[idx] = 0.0,
                        }
                    }
                }
                None => o.hetero_acgt_ins[idx] = 0.0,
            }
        } else if outcome.hetero_ins_ins {
            let idx = index_of(&o.hetero_ins_ins, max);
            // With lengths, the proposed pair first; the two best insertions otherwise.
            let proposed = o.hetero_ins_ins_lengths.get(idx).and_then(|&(l1, l2)| {
                let first = insertion(Some(l1), "")?;
                Some(vec![first, insertion(Some(l2), first)?])
            });
            match proposed.unwrap_or_else(|| alt.insertion_bases_multi(max_indel_len))[..] {
                [first, second] => {
                    reference = Some(ref_str.clone());
                    alternate = Some(first.to_owned());
                    if second != first {
                        alternate = Some(format!("{second},{first}"));
                    } else {
                        o.hetero_ins_ins[idx] = 0.0;
                    }
                }
                _ => o.hetero_ins_ins[idx] = 0.0,
            }
        } else if outcome.homo_del {
            let idx = index_of(&o.homo_del, max);
            match deletion(o.homo_del_lengths.get(idx).copied(), "") {
                Some(del) => {
                    reference = Some(format!("{ref_char}{del}"));
                    alternate = Some(ref_str.clone());
                }
                None => o.homo_del[idx] = 0.0,
            }
        } else if outcome.hetero_acgt_del {
            let idx = index_of(&o.hetero_acgt_del, max);
            let hetero_base = char::from(o.hetero_acgt_del_bases[idx]);
            match deletion(o.hetero_acgt_del_lengths.get(idx).copied(), "") {
                Some(del) => {
                    reference = Some(format!("{ref_char}{del}"));
                    alternate = Some(if hetero_base != ref_char {
                        format!("{ref_char},{hetero_base}{del}")
                    } else {
                        ref_str.clone()
                    });
                }
                None => o.hetero_acgt_del[idx] = 0.0,
            }
        } else if outcome.hetero_del_del {
            let idx = index_of(&o.hetero_del_del, max);
            // With lengths, the proposed pair (longer length first); the two best deletions
            // otherwise. Either way the longer deletion comes first.
            let proposed = o.hetero_del_del_lengths.get(idx).and_then(|&(short, long)| {
                let first = deletion(Some(long), "")?;
                let second = deletion(Some(short), first)?;
                Some(if first.chars().count() > second.chars().count() { (first, second) } else { (second, first) })
            });
            match proposed.or_else(|| alt.deletion_bases_multi(max_indel_len)) {
                Some((longer, shorter)) => {
                    let ref_allele = format!("{ref_char}{longer}");
                    let alt1 = ref_str.clone();
                    let alt2: String =
                        std::iter::once(ref_char).chain(ref_allele.chars().skip(shorter.chars().count() + 1)).collect();
                    alternate = Some(alt1.clone());
                    if alt1 != alt2 && ref_allele != alt1 && ref_allele != alt2 {
                        alternate = Some(format!("{alt1},{alt2}"));
                    } else {
                        o.hetero_del_del[idx] = 0.0;
                    }
                    reference = Some(ref_allele);
                }
                None => o.hetero_del_del[idx] = 0.0,
            }
        } else if outcome.ins_del {
            let idx = index_of(&o.hetero_ins_del, max);
            let (del_length, ins_length) = o.hetero_ins_del_lengths.get(idx).copied().unzip();
            match (insertion(ins_length, ""), deletion(del_length, "")) {
                (Some(ins), Some(del)) => {
                    reference = Some(format!("{ref_char}{del}"));
                    alternate = Some(format!("{ref_char},{ins}{del}"));
                }
                _ => o.hetero_ins_del[idx] = 0.0,
            }
        } else {
            unreachable!("the maximum {max} comes from one of the outcome lists");
        }

        if let (Some(reference), Some(alternate)) = (&reference, &alternate) {
            return Some(Decoded {
                outcome,
                reference: reference.clone(),
                alternate: alternate.clone(),
                probability: max,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn heads<'a>(gt21: &'a [f32], genotype: &'a [f32], lengths: Option<[&'a [f32]; 2]>) -> Heads<'a> {
        Heads { gt21, genotype, variant_lengths: lengths }
    }

    #[test]
    fn indel_length_reference_needs_both_heads_at_zero() {
        let mut gt21 = [0.0f32; GT21_CLASSES];
        gt21[gt21::AA] = 0.9;
        let genotype = [0.9f32, 0.05, 0.05];
        let mut zero = [0.0f32; VARIANT_LENGTH_CLASSES];
        zero[variant_length::INDEX_OFFSET] = 0.8;
        let mut del = [0.0f32; VARIANT_LENGTH_CLASSES];
        del[variant_length::INDEX_OFFSET - 2] = 0.6;
        del[variant_length::INDEX_OFFSET] = 0.4;
        let early = possible_outcomes(&heads(&gt21, &genotype, Some([&zero, &zero])), gt21::AA);
        assert_eq!(early, Possible::Reference(0.8 * 0.8 * 0.9 * 0.9));
        let full = possible_outcomes(&heads(&gt21, &genotype, Some([&zero, &del])), gt21::AA);
        let Possible::Outcomes(o) = full else { panic!("no early return when one head proposes a deletion") };
        assert_eq!(o.homo_ref, 0.8f32 * 0.4 * 0.9 * 0.9);
        assert_eq!((o.homo_ins.len(), o.hetero_ins_ins.len(), o.hetero_acgt_ins.len()), (16, 136, 64));
        assert_eq!((o.homo_del.len(), o.hetero_del_del.len(), o.hetero_acgt_del.len()), (16, 241, 64));
        assert_eq!(o.hetero_ins_del.len(), 256);
        assert_eq!(&o.hetero_acgt_del_bases[..5], b"ACGTA");
        assert_eq!(&o.hetero_acgt_del_lengths[..5], &[1, 1, 1, 1, 2]);
    }

    #[test]
    fn proposed_deletion_length_picks_the_allele() {
        // A homo deletion of length 2, although the 1-base deletion has more reads.
        let mut gt21 = [0.0f32; GT21_CLASSES];
        gt21[gt21::DEL_DEL] = 0.9;
        let genotype = [0.05f32, 0.9, 0.05];
        let mut del2 = [0.0f32; VARIANT_LENGTH_CLASSES];
        del2[variant_length::INDEX_OFFSET - 2] = 0.9;
        let alt = AltInfo::parse("30-DT 12 DTT 8 RA 10").expect("valid alt info");
        let call = decode("A", &heads(&gt21, &genotype, Some([&del2, &del2])), &alt, 50).expect("decodable");
        assert!(call.outcome.homo_del);
        assert_eq!((call.reference.as_str(), call.alternate.as_str()), ("ATT", "A"));
        // Without the heads, the most supported deletion wins.
        let call = decode("A", &heads(&gt21, &genotype, None), &alt, 50).expect("decodable");
        assert_eq!((call.reference.as_str(), call.alternate.as_str()), ("AT", "A"));
    }
}

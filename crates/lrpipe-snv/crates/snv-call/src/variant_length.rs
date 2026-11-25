//! `VariantLength` (`task/variant_length.py`) and the `*_tuples_from` helpers of
//! `CallVariants.py` that pair indel lengths with their `float32` probabilities.
//!
//! Each indel-length head has [`CLASSES`] outputs for lengths `-16..=16`; index
//! [`INDEX_OFFSET`]` + n` is length `n` (negative for deletions). Products keep Python's
//! left-to-right order, all in `f32` as numpy scalar arithmetic does.

/// `VariantLength.index_offset`: the index of length 0.
pub(crate) const INDEX_OFFSET: usize = 16;
/// `VariantLength.max`: the longest length a head distinguishes ("16 or more").
pub(crate) const MAX: usize = INDEX_OFFSET;
/// `VariantLength.output_label_count`.
pub(crate) const CLASSES: usize = 2 * INDEX_OFFSET + 1;

/// Probability of insertion length `n` (`p[n + index_offset]`).
fn ins(p: &[f32], n: usize) -> f32 {
    p[INDEX_OFFSET + n]
}

/// Probability of deletion length `n` (`p[-n + index_offset]`).
fn del(p: &[f32], n: usize) -> f32 {
    p[INDEX_OFFSET - n]
}

/// `homo_Ins_tuples_from`.
pub(crate) fn homo_ins(l1: &[f32], l2: &[f32], extra: f32) -> Vec<(usize, f32)> {
    (1..=MAX).map(|i| (i, ins(l1, i) * ins(l2, i) * extra)).collect()
}

/// `hetero_Ins_tuples_from`.
pub(crate) fn hetero_ins(l1: &[f32], l2: &[f32]) -> Vec<(usize, f32)> {
    (1..=MAX).map(|i| (i, ins(l1, 0) * ins(l2, i))).collect()
}

/// `hetero_InsIns_tuples_from`: `i <= j`.
pub(crate) fn hetero_ins_ins(l1: &[f32], l2: &[f32], extra: f32) -> Vec<((usize, usize), f32)> {
    (1..=MAX).flat_map(|i| (i..=MAX).map(move |j| ((i, j), ins(l1, i) * ins(l2, j) * extra))).collect()
}

/// `homo_Del_tuples_from`.
pub(crate) fn homo_del(l1: &[f32], l2: &[f32], extra: f32) -> Vec<(usize, f32)> {
    (1..=MAX).map(|i| (i, del(l1, i) * del(l2, i) * extra)).collect()
}

/// `hetero_Del_tuples_from`.
pub(crate) fn hetero_del(l1: &[f32], l2: &[f32]) -> Vec<(usize, f32)> {
    (1..=MAX).map(|i| (i, del(l1, i) * del(l2, 0))).collect()
}

/// `hetero_DelDel_tuples_from`: the lengths sorted ascending.
///
/// The original skips `i == j` unless `i` equals `index_offset`, so `(16, 16)` is kept.
pub(crate) fn hetero_del_del(l1: &[f32], l2: &[f32], extra: f32) -> Vec<((usize, usize), f32)> {
    (1..=MAX)
        .flat_map(|i| (1..=MAX).map(move |j| (i, j)))
        .filter(|&(i, j)| !(i == j && i != INDEX_OFFSET && j != INDEX_OFFSET))
        .map(|(i, j)| ((i.min(j), i.max(j)), del(l1, i) * del(l2, j) * extra))
        .collect()
}

/// `hetero_InsDel_tuples_from`: `(deletion length, insertion length)`.
pub(crate) fn hetero_ins_del(l1: &[f32], l2: &[f32], extra: f32) -> Vec<((usize, usize), f32)> {
    (1..=MAX).flat_map(|i| (1..=MAX).map(move |j| ((i, j), del(l1, i) * ins(l2, j) * extra))).collect()
}

/// `variant_length if variant_length and variant_length < VariantLength.max else None`: the
/// length to look for in the alt info, if the head proposes a specific one.
pub(crate) fn proposed(length: Option<usize>) -> Option<usize> {
    length.filter(|&n| n != 0 && n < MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_counts_match_the_original() {
        let p = [1.0f32; CLASSES];
        assert_eq!(homo_ins(&p, &p, 1.0).len(), 16);
        assert_eq!(hetero_ins_ins(&p, &p, 1.0).len(), 136);
        // 256 pairs minus the 15 equal ones below 16.
        let del_del = hetero_del_del(&p, &p, 1.0);
        assert_eq!(del_del.len(), 241);
        assert!(del_del.iter().any(|&(t, _)| t == (16, 16)));
        assert!(del_del.iter().all(|&((a, b), _)| a <= b));
        assert_eq!(hetero_ins_del(&p, &p, 1.0).len(), 256);
    }

    #[test]
    fn lengths_index_both_sides_of_the_offset() {
        let p: Vec<f32> = (0..CLASSES).map(|i| i as f32).collect();
        // Deletion of 3 is index 13, insertion of 3 is index 19.
        assert_eq!(homo_del(&p, &p, 1.0)[2], (3, 169.0));
        assert_eq!(homo_ins(&p, &p, 1.0)[2], (3, 361.0));
        assert_eq!(hetero_ins_del(&p, &p, 1.0)[0], ((1, 1), 15.0 * 17.0));
    }

    #[test]
    fn only_short_lengths_are_proposed() {
        assert_eq!(proposed(None), None);
        assert_eq!(proposed(Some(0)), None);
        assert_eq!(proposed(Some(15)), Some(15));
        assert_eq!(proposed(Some(16)), None);
    }
}

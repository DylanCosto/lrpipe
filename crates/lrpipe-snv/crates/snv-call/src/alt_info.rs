//! The candidate "alt info" as `CallVariants.py` reads it: a Python dict of allele keys
//! (`XC`, `IACG`, `DTT`, `RA`, ...) to read counts, kept in insertion order because the
//! original's tie-breaking depends on it.

/// Read depth and allele counts of one candidate (`read_depth` and `alt_info_dict`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AltInfo {
    pub(crate) depth: i64,
    /// `(key, count)` in Python dict order.
    entries: Vec<(String, i64)>,
}

/// Python's `max(d, key=d.get)`: the first key with the largest value.
fn first_max<'a>(items: impl Iterator<Item = (&'a str, i64)>) -> Option<&'a str> {
    let mut best: Option<(&str, i64)> = None;
    for (key, count) in items {
        if best.is_none_or(|(_, c)| count > c) {
            best = Some((key, count));
        }
    }
    best.map(|(k, _)| k)
}

/// `[k for k, _ in sorted(items, key=count)][::-1]`: ascending stable sort, then reversed, so
/// equal counts come out in reverse insertion order.
fn sorted_desc_reversed(mut items: Vec<(&str, i64)>) -> Vec<&str> {
    items.sort_by_key(|&(_, c)| c);
    items.into_iter().rev().map(|(k, _)| k).collect()
}

/// The body of `insertion_bases_using_alt_info_from` and `deletion_bases_using_alt_info_from`
/// without `return_multi`: keys of length `proposed` go to one dict, the others of length
/// `1..=max_len` to another (both skip `ignore`), and the first dict wins when non-empty.
fn indel_bases<'a>(
    items: impl Iterator<Item = (&'a str, i64)>,
    proposed: Option<usize>,
    ignore: &str,
    max_len: usize,
) -> Option<&'a str> {
    let mut proposed_items = Vec::new();
    let mut other_items = Vec::new();
    for (key, count) in items.filter(|&(k, _)| k != ignore) {
        let len = key.chars().count();
        if proposed == Some(len) {
            proposed_items.push((key, count));
        } else if (1..=max_len).contains(&len) {
            other_items.push((key, count));
        }
    }
    let items = if proposed_items.is_empty() { other_items } else { proposed_items };
    first_max(items.into_iter()).filter(|k| !k.is_empty())
}

impl AltInfo {
    /// Parses `"{depth}-{key} {count} {key} {count} ..."` like `output_with` does.
    pub(crate) fn parse(alt_info: &str) -> Option<Self> {
        let fields: Vec<&str> = alt_info.trim_end().split('-').collect();
        let depth = fields.first()?.parse().ok()?;
        let indel_str = fields.get(1).copied().unwrap_or("");
        let seqs: Vec<&str> = indel_str.split(' ').collect();
        let mut entries: Vec<(String, i64)> = Vec::new();
        // dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])): every count is converted
        // first, so one bad count fails the whole parse, as `int()` would raise.
        let counts: Vec<i64> = seqs.iter().skip(1).step_by(2).map(|s| s.parse()).collect::<Result<_, _>>().ok()?;
        for (key, count) in seqs.iter().step_by(2).zip(counts) {
            match entries.iter_mut().find(|(k, _)| k == key) {
                Some(entry) => entry.1 = count,
                None => entries.push(((*key).to_owned(), count)),
            }
        }
        Some(Self { depth, entries })
    }

    fn with_prefix(&self, prefix: char) -> impl Iterator<Item = (&str, i64)> {
        self.entries.iter().filter_map(move |(k, c)| k.strip_prefix(prefix).map(|rest| (rest, *c)))
    }

    /// `find_alt_base`: SNP bases by count, highest first (ties keep dict order), and the
    /// alternate base to use: `proposed`, unless it is absent or at least 9 reads behind the top.
    pub(crate) fn find_alt_base(&self, proposed: Option<&str>) -> (Vec<String>, Option<String>) {
        let mut sorted: Vec<(String, i64)> = self
            .entries
            .iter()
            .filter(|(k, _)| k.starts_with('X'))
            .map(|(k, c)| (k.chars().nth(1).map(String::from).unwrap_or_default(), *c))
            .collect();
        // sorted(..., reverse=True) is stable: ties keep their order.
        sorted.sort_by_key(|&(_, c)| std::cmp::Reverse(c));
        if sorted.is_empty() {
            return (Vec::new(), None);
        }
        let alt_count = sorted.iter().find(|(b, _)| Some(b.as_str()) == proposed).map(|(_, c)| *c);
        let alternate = match alt_count {
            Some(c) if sorted[0].1 - c < 9 => proposed.map(str::to_owned),
            _ => Some(sorted[0].0.clone()),
        };
        (sorted.into_iter().map(|(b, _)| b).collect(), alternate)
    }

    /// `insertion_bases_using_alt_info_from(propose_insertion_length=proposed,
    /// maximum_insertion_length=max_len, insertion_bases_to_ignore=ignore)`, `None` for `""`.
    /// Keys include the reference base.
    ///
    /// With a proposed length, the most supported insertion of exactly that length wins; if there
    /// is none, the most supported one of any other length up to `max_len`.
    pub(crate) fn insertion_bases(&self, proposed: Option<usize>, ignore: &str, max_len: usize) -> Option<&str> {
        // `propose_insertion_length += 1`: the key includes the reference base.
        indel_bases(self.with_prefix('I'), proposed.map(|n| n + 1), ignore, max_len)
    }

    /// The same with `return_multi=True`: the two best insertions.
    pub(crate) fn insertion_bases_multi(&self, max_len: usize) -> Vec<&str> {
        let items: Vec<_> = self.with_prefix('I').filter(|(k, _)| (1..=max_len).contains(&k.chars().count())).collect();
        sorted_desc_reversed(items).into_iter().take(2).collect()
    }

    /// `deletion_bases_using_alt_info_from(propose_deletion_length=proposed,
    /// maximum_deletion_length=max_len, deletion_bases_to_ignore=ignore)`: the deleted bases
    /// (without the reference base), `None` for `""`.
    pub(crate) fn deletion_bases(&self, proposed: Option<usize>, ignore: &str, max_len: usize) -> Option<&str> {
        indel_bases(self.with_prefix('D'), proposed, ignore, max_len)
    }

    /// The same with `return_multi=True`: the two best deletions, longer first, or nothing if
    /// there are fewer than two.
    pub(crate) fn deletion_bases_multi(&self, max_len: usize) -> Option<(&str, &str)> {
        let items: Vec<_> = self.with_prefix('D').filter(|(k, _)| (1..=max_len).contains(&k.chars().count())).collect();
        let sorted = sorted_desc_reversed(items);
        let [a, b, ..] = sorted[..] else { return None };
        Some(if a.chars().count() > b.chars().count() { (a, b) } else { (b, a) })
    }

    /// `decode_alt_info`: per-type counts and the reference count (at least 0).
    pub(crate) fn decode(&self) -> DecodedAlleles {
        let mut decoded = DecodedAlleles::default();
        for (key, count) in &self.entries {
            let mut chars = key.chars();
            match chars.next() {
                Some('X') => {
                    let base = chars.next().map(String::from).unwrap_or_default();
                    set(&mut decoded.snp, base, *count);
                }
                Some('I') => set(&mut decoded.ins, chars.as_str().to_owned(), *count),
                Some('D') => set(&mut decoded.del, chars.as_str().to_owned(), *count),
                Some('R') => decoded.ref_count = *count,
                _ => {}
            }
        }
        decoded.ref_count = decoded.ref_count.max(0);
        decoded
    }
}

/// `d[key] = value` on an insertion-ordered list.
fn set(items: &mut Vec<(String, i64)>, key: String, value: i64) {
    match items.iter_mut().find(|(k, _)| *k == key) {
        Some(entry) => entry.1 = value,
        None => items.push((key, value)),
    }
}

/// `alt_type_list` and `ref_count` from `decode_alt_info`.
#[derive(Debug, Default)]
pub(crate) struct DecodedAlleles {
    pub(crate) snp: Vec<(String, i64)>,
    pub(crate) ins: Vec<(String, i64)>,
    pub(crate) del: Vec<(String, i64)>,
    pub(crate) ref_count: i64,
}

/// `d[key] if key in d else 0`.
pub(crate) fn count_of(items: &[(String, i64)], key: &str) -> i64 {
    items.iter().find(|(k, _)| k == key).map_or(0, |(_, c)| *c)
}

/// `get_long_indel_read_count` for a proposed insertion (`--enable_long_indel`): the reads of
/// the other insertions whose length is within 10% of it, for insertions longer than 50 bases.
/// The length compared is the alt-info key's, which includes the anchor base, against a range
/// around the inserted length, as in the original.
///
/// The original calls it for deletions without `is_del=True`, which makes the range empty, so long
/// deletions never gain reads and have no counterpart here.
pub(crate) fn long_insertion_reads(items: &[(String, i64)], proposed: &str) -> i64 {
    // `param.maximum_variant_length_that_need_infer`, whatever the long-indel setting.
    const MAX_INFERRED: f64 = 50.0;
    const DISTANCE: f64 = 0.1;
    let len = proposed.chars().count();
    if len as f64 <= MAX_INFERRED {
        return 0;
    }
    let inserted = (len - 1) as f64;
    let (min, max) = ((inserted * (1.0 - DISTANCE)).max(MAX_INFERRED), inserted * (1.0 + DISTANCE));
    items
        .iter()
        .filter(|(k, _)| k != proposed)
        .filter(|(k, _)| (min..=max).contains(&(k.chars().count() as f64)))
        .map(|(_, c)| c)
        .sum()
}

/// `[d[k] for k in d if len(k) == n]`, first element or 0.
pub(crate) fn first_count_with_len(items: &[(String, i64)], n: i64) -> i64 {
    items.iter().find(|(k, _)| k.chars().count() as i64 == n).map_or(0, |(_, c)| *c)
}

#[cfg(test)]
mod long_indel_tests {
    use super::*;

    #[test]
    fn nearby_long_insertions_add_their_reads() {
        let ins = |n: usize| "A".repeat(n);
        let items = vec![(ins(101), 5), (ins(95), 3), (ins(111), 2), (ins(112), 7), (ins(89), 1), (ins(40), 9)];
        // Inserted length 100: keys of 90..=110 characters other than the proposed one.
        assert_eq!(long_insertion_reads(&items, &ins(101)), 3);
        assert_eq!(long_insertion_reads(&items, &ins(50)), 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_like_python_dict() {
        let a = AltInfo::parse("30-XC 5 XG 5 DTT 3 IAC 2 RA 15 ").unwrap();
        assert_eq!(a.depth, 30);
        assert_eq!(a.find_alt_base(Some("G")), (vec!["C".into(), "G".into()], Some("G".into())));
        assert_eq!(a.find_alt_base(None).1, Some("C".into()));
        assert_eq!(a.insertion_bases(None, "", 50), Some("AC"));
        assert_eq!(a.deletion_bases(None, "", 50), Some("TT"));
        assert_eq!(AltInfo::parse("12-").unwrap().find_alt_base(None), (vec![], None));
    }

    #[test]
    fn far_behind_proposal_is_replaced() {
        let a = AltInfo::parse("40-XC 20 XG 11").unwrap();
        assert_eq!(a.find_alt_base(Some("G")).1, Some("C".into()));
        assert_eq!(a.find_alt_base(Some("T")).1, Some("C".into()));
    }

    #[test]
    fn multi_ties_come_out_reversed() {
        let a = AltInfo::parse("30-IAC 4 IAG 4 IAT 1 DT 2 DTT 2 DTTT 1").unwrap();
        assert_eq!(a.insertion_bases_multi(50), vec!["AG", "AC"]);
        // Reversed order gives TT then T; longer first.
        assert_eq!(a.deletion_bases_multi(50), Some(("TT", "T")));
    }

    #[test]
    fn proposed_lengths_win_then_fall_back() {
        let a = AltInfo::parse("30-IAC 9 IACG 4 IAGG 5 DT 9 DTT 2 DTTT 3").unwrap();
        // Insertion keys include the reference base: a proposed length of 2 matches ACG and AGG.
        assert_eq!(a.insertion_bases(Some(2), "", 50), Some("AGG"));
        assert_eq!(a.insertion_bases(Some(2), "AGG", 50), Some("ACG"));
        assert_eq!(a.insertion_bases(Some(5), "", 50), Some("AC"));
        assert_eq!(a.insertion_bases(Some(5), "AC", 50), Some("AGG"));
        assert_eq!(a.insertion_bases(None, "", 2), Some("AC"));
        assert_eq!(a.deletion_bases(Some(2), "", 50), Some("TT"));
        assert_eq!(a.deletion_bases(Some(2), "TT", 50), Some("T"));
        assert_eq!(a.deletion_bases(Some(4), "T", 50), Some("TTT"));
        assert_eq!(a.deletion_bases(None, "", 0), None);
    }
}

//! Counters that iterate in the same order as the original's `KH_COUNTER` and `KH_INT_COUNTER`.
//!
//! The alt-info string lists insertions and deletions in khash bucket order, so matching the
//! original means reproducing klib khash 0.2.8: its hash functions (X31 for strings, the identity
//! for `khint32_t`), quadratic probing, the 0.77 load factor and the in-place "kick-out" rehash on
//! resize. Keys are only ever inserted, never deleted, which is all the original does with
//! these counters. (`snv-pileup` has the string-only twin of this module.)

/// `__ac_HASH_UPPER` in `khash.h`.
const HASH_UPPER: f64 = 0.77;

/// A key type with khash's hash function for it.
pub(crate) trait KhKey: Clone + PartialEq {
    /// The `khint_t` hash khash computes for this key.
    fn kh_hash(&self) -> u32;
}

/// `KHASH_MAP_INIT_INT`: `kh_int_hash_func` is the key itself.
impl KhKey for i32 {
    fn kh_hash(&self) -> u32 {
        *self as u32
    }
}

/// `KHASH_MAP_INIT_STR`: `__ac_X31_hash_string`; `char` is signed on the platforms the
/// original builds for.
impl KhKey for Box<[u8]> {
    fn kh_hash(&self) -> u32 {
        let mut bytes = self.iter().map(|&b| b as i8 as u32);
        let Some(first) = bytes.next() else { return 0 };
        bytes.fold(first, |h, c| (h << 5).wrapping_sub(h).wrapping_add(c))
    }
}

/// `kroundup32`: the next power of two at or above `x` (0 stays 0).
fn roundup32(x: u32) -> u32 {
    if x == 0 { 0 } else { x.next_power_of_two() }
}

fn upper_bound_for(n_buckets: u32) -> u32 {
    (f64::from(n_buckets) * HASH_UPPER + 0.5) as u32
}

/// Insert-only key-to-`int` map with khash's bucket layout (`kh_counter_add`,
/// `kh_int_counter_add`).
#[derive(Debug, Clone)]
pub(crate) struct KhCounter<K> {
    slots: Vec<Option<(K, i32)>>,
    size: u32,
    n_occupied: u32,
    upper_bound: u32,
}

impl<K: KhKey> Default for KhCounter<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: KhKey> KhCounter<K> {
    pub(crate) fn new() -> Self {
        Self { slots: Vec::new(), size: 0, n_occupied: 0, upper_bound: 0 }
    }

    fn n_buckets(&self) -> u32 {
        self.slots.len() as u32
    }

    /// Adds `val` to the count of `key`, inserting it with `val` if it is new.
    pub(crate) fn add(&mut self, key: K, val: i32) {
        let x = self.put(&key);
        match &mut self.slots[x] {
            Some((_, cur)) => *cur = cur.wrapping_add(val),
            slot @ None => *slot = Some((key, val)),
        }
    }

    /// The count of `key`, if present (`kh_get`).
    pub(crate) fn get(&self, key: &K) -> Option<i32> {
        let n = self.n_buckets();
        if n == 0 {
            return None;
        }
        let mask = n - 1;
        let mut i = key.kh_hash() & mask;
        let last = i;
        let mut step = 0;
        loop {
            match &self.slots[i as usize] {
                None => return None,
                Some((k, v)) if k == key => return Some(*v),
                Some(_) => {}
            }
            step += 1;
            i = (i + step) & mask;
            if i == last {
                return None;
            }
        }
    }

    /// Entries in bucket order, as `for (k = kh_begin(h); k != kh_end(h); ++k)` visits them.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&K, i32)> {
        self.slots.iter().filter_map(|slot| slot.as_ref().map(|(k, v)| (k, *v)))
    }

    /// `kh_put`: returns the bucket holding `key`, leaving it empty if the key is new. The
    /// caller fills a new bucket, so `size` counts it already.
    fn put(&mut self, key: &K) -> usize {
        if self.n_occupied >= self.upper_bound {
            if self.n_buckets() > self.size << 1 {
                self.resize(self.n_buckets() - 1);
            } else {
                self.resize(self.n_buckets() + 1);
            }
        }
        let n = self.n_buckets();
        let mask = n - 1;
        let mut i = key.kh_hash() & mask;
        // No bucket is ever deleted, so the probe stops at the key or at the first empty bucket.
        let mut step = 0;
        let last = i;
        loop {
            match &self.slots[i as usize] {
                None => break,
                Some((k, _)) if k == key => return i as usize,
                Some(_) => {}
            }
            step += 1;
            i = (i + step) & mask;
            assert!(i != last, "khash table full with {} of {n} buckets used", self.size);
        }
        self.size += 1;
        self.n_occupied += 1;
        i as usize
    }

    /// `kh_resize`, including its in-place kick-out rehash.
    fn resize(&mut self, requested: u32) {
        let new_n = roundup32(requested).max(4);
        if self.size >= upper_bound_for(new_n) {
            return;
        }
        let old_n = self.n_buckets() as usize;
        let new_mask = new_n - 1;
        if old_n < new_n as usize {
            self.slots.resize(new_n as usize, None);
        }
        // `pending[j]` is true while old bucket `j` still holds an entry that has not moved.
        let mut pending: Vec<bool> = self.slots[..old_n].iter().map(Option::is_some).collect();
        let mut taken = vec![false; new_n as usize];
        for j in 0..old_n {
            if !pending[j] {
                continue;
            }
            pending[j] = false;
            let mut entry = self.slots[j].take();
            loop {
                let Some((key, _)) = &entry else { unreachable!("only live entries are rehashed") };
                let mut i = key.kh_hash() & new_mask;
                let mut step = 0;
                while taken[i as usize] {
                    step += 1;
                    i = (i + step) & new_mask;
                }
                taken[i as usize] = true;
                let i = i as usize;
                if i < old_n && pending[i] {
                    // Kick out the old entry sitting in our new bucket and place it next.
                    pending[i] = false;
                    entry = std::mem::replace(&mut self.slots[i], entry);
                } else {
                    self.slots[i] = entry;
                    break;
                }
            }
        }
        self.slots.truncate(new_n as usize);
        self.n_occupied = self.size;
        self.upper_bound = upper_bound_for(new_n);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(s: &str) -> Box<[u8]> {
        s.as_bytes().into()
    }

    #[test]
    fn x31_matches_reference_values() {
        assert_eq!(key("").kh_hash(), 0);
        assert_eq!(key("A").kh_hash(), 65);
        assert_eq!(key("ACG").kh_hash(), (65 * 31 + 67) * 31 + 71);
    }

    #[test]
    fn add_and_get() {
        let mut c = KhCounter::new();
        for (k, v) in [(3, 1), (7, 1), (3, 1), (-1, 5)] {
            c.add(k, v);
        }
        assert_eq!(c.get(&3), Some(2));
        assert_eq!(c.get(&-1), Some(5));
        assert_eq!(c.get(&4), None);
        assert_eq!(KhCounter::<i32>::new().get(&0), None);
    }

    /// Iteration order against the original's own khash (fixtures from
    /// `tools/khash_oracle_fullalign.c`).
    #[test]
    fn iteration_order_matches_c_khash() {
        let cases = include_str!("../tests/data/khash_cases.txt").lines();
        let expected = include_str!("../tests/data/khash_expected.txt").lines();
        for (i, (case, want)) in cases.zip(expected).enumerate() {
            let (kind, keys) = case.split_once(' ').expect("case has a kind and keys");
            let got: Vec<String> = match kind {
                "int" => {
                    let mut c = KhCounter::new();
                    for k in keys.split(' ') {
                        c.add(k.parse::<i32>().expect("int key"), 1);
                    }
                    c.iter().map(|(k, v)| format!("{k}={v}")).collect()
                }
                _ => {
                    let mut c = KhCounter::new();
                    for k in keys.split(' ') {
                        c.add(key(k), 1);
                    }
                    c.iter().map(|(k, v)| format!("{}={v}", String::from_utf8_lossy(k))).collect()
                }
            };
            assert_eq!(got.join(" "), want, "case {i}: {case}");
        }
    }
}

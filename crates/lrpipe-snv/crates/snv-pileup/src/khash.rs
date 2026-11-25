//! A string counter that iterates in the same order as the original's `khash_t(KH_COUNTER)`.
//!
//! The original writes insertion alleles into the candidate string in khash bucket order, so
//! matching its output means reproducing klib khash 0.2.8 exactly: the X31 string hash, quadratic
//! probing, the 0.77 load factor and the in-place "kick-out" rehash on resize. Keys are only ever
//! inserted, never deleted, which is all `kh_counter_increment` does.

/// `__ac_HASH_UPPER` in `khash.h`.
const HASH_UPPER: f64 = 0.77;

#[derive(Debug, Clone, PartialEq, Eq)]
enum Slot {
    Empty,
    Live(Vec<u8>, i32),
}

/// Insert-only string-to-count map with khash's bucket layout.
#[derive(Debug, Clone, Default)]
pub(crate) struct StrCounter {
    slots: Vec<Slot>,
    size: u32,
    n_occupied: u32,
    upper_bound: u32,
}

/// `__ac_X31_hash_string`; `char` is signed on the platforms the original builds for.
fn x31_hash(key: &[u8]) -> u32 {
    let mut bytes = key.iter().map(|&b| b as i8 as u32);
    let Some(first) = bytes.next() else { return 0 };
    bytes.fold(first, |h, c| (h << 5).wrapping_sub(h).wrapping_add(c))
}

/// `kroundup32`: the next power of two at or above `x` (0 stays 0).
fn roundup32(x: u32) -> u32 {
    if x == 0 { 0 } else { x.next_power_of_two() }
}

fn upper_bound_for(n_buckets: u32) -> u32 {
    (f64::from(n_buckets) * HASH_UPPER + 0.5) as u32
}

impl StrCounter {
    fn n_buckets(&self) -> u32 {
        self.slots.len() as u32
    }

    /// `kh_counter_increment`.
    pub(crate) fn increment(&mut self, key: &[u8]) {
        let x = self.put(key);
        match &mut self.slots[x] {
            Slot::Live(_, val) => *val += 1,
            slot @ Slot::Empty => *slot = Slot::Live(key.to_vec(), 1),
        }
    }

    /// Entries in bucket order, as `for (k = kh_begin(h); k != kh_end(h); ++k)` visits them.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&[u8], i32)> {
        self.slots.iter().filter_map(|slot| match slot {
            Slot::Live(key, val) => Some((key.as_slice(), *val)),
            Slot::Empty => None,
        })
    }

    /// Sum and maximum of the counts (`kh_counter_stats`).
    pub(crate) fn stats(&self) -> (usize, usize) {
        self.iter().fold((0, 0), |(sum, max), (_, val)| (sum + val as usize, max.max(val as usize)))
    }

    /// `kh_put`: returns the bucket holding `key`, leaving it `Empty` if the key is new. The
    /// caller fills a new bucket, so `size` counts it already.
    fn put(&mut self, key: &[u8]) -> usize {
        if self.n_occupied >= self.upper_bound {
            if self.n_buckets() > self.size << 1 {
                self.resize(self.n_buckets() - 1);
            } else {
                self.resize(self.n_buckets() + 1);
            }
        }
        let n = self.n_buckets();
        let mask = n - 1;
        let mut i = x31_hash(key) & mask;
        // No bucket is ever deleted, so the probe stops at the key or at the first empty bucket.
        let mut step = 0;
        let last = i;
        loop {
            match &self.slots[i as usize] {
                Slot::Empty => break,
                Slot::Live(k, _) if k == key => return i as usize,
                Slot::Live(..) => {}
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
            self.slots.resize(new_n as usize, Slot::Empty);
        }
        // `pending[j]` is true while old bucket `j` still holds an entry that has not moved.
        let mut pending: Vec<bool> = self.slots[..old_n].iter().map(|s| matches!(s, Slot::Live(..))).collect();
        let mut taken = vec![false; new_n as usize];
        for j in 0..old_n {
            if !pending[j] {
                continue;
            }
            pending[j] = false;
            let mut entry = std::mem::replace(&mut self.slots[j], Slot::Empty);
            loop {
                let Slot::Live(key, _) = &entry else { unreachable!("only live entries are rehashed") };
                let mut i = x31_hash(key) & new_mask;
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

    #[test]
    fn x31_matches_reference_values() {
        // h = s[0]; h = h*31 + c for the rest.
        assert_eq!(x31_hash(b""), 0);
        assert_eq!(x31_hash(b"A"), 65);
        assert_eq!(x31_hash(b"AC"), 65 * 31 + 67);
        assert_eq!(x31_hash(b"ACG"), (65 * 31 + 67) * 31 + 71);
    }

    #[test]
    fn counts_and_stats() {
        let mut c = StrCounter::default();
        for key in [&b"A"[..], b"AC", b"A", b"TTT", b"A"] {
            c.increment(key);
        }
        let mut got: Vec<_> = c.iter().map(|(k, v)| (k.to_vec(), v)).collect();
        got.sort();
        assert_eq!(got, vec![(b"A".to_vec(), 3), (b"AC".to_vec(), 1), (b"TTT".to_vec(), 1)]);
        assert_eq!(c.stats(), (5, 3));
    }

    /// Iteration order against the original's own khash (fixtures from `tools/khash_oracle.c`).
    #[test]
    fn iteration_order_matches_c_khash() {
        let cases = include_str!("../tests/data/khash_cases.txt").lines();
        let expected = include_str!("../tests/data/khash_expected.txt").lines();
        for (i, (case, want)) in cases.zip(expected).enumerate() {
            let mut c = StrCounter::default();
            for key in case.split(' ') {
                c.increment(key.as_bytes());
            }
            let got: Vec<String> = c.iter().map(|(k, v)| format!("{}={v}", String::from_utf8_lossy(k))).collect();
            assert_eq!(got.join(" "), want, "case {i}: {case}");
        }
    }

    #[test]
    fn grows_through_several_resizes() {
        let mut c = StrCounter::default();
        let keys: Vec<Vec<u8>> = (0..500).map(|i| format!("K{i}").into_bytes()).collect();
        for key in &keys {
            c.increment(key);
        }
        assert_eq!(c.iter().count(), keys.len());
        assert!(c.n_buckets().is_power_of_two());
        assert!(c.size <= c.upper_bound, "size {} bound {}", c.size, c.upper_bound);
    }
}

//! Bit-exact emulation of the libstdc++ behaviour that WhatsHap's C++ core
//! (implicitly) depends on for tie-breaking.
//!
//! * `std::hash<std::string>` is used by `ReadSet::sort()` to order reads that
//!   start at the same variant position.
//! * Iteration order of `std::unordered_set<int>` leaks into read selection
//!   (`readselect.pyx`), where it determines insertion order into Python sets.
//! * Iteration order of `std::unordered_map<Permutation, ...>` breaks ties in
//!   the polyploid switch/flip calculator used by `compare`.

/// libstdc++ `std::_Hash_bytes` (the 64-bit Murmur-derived variant).
pub fn hash_bytes(data: &[u8], seed: u64) -> u64 {
    const MUL: u64 = (0xc6a4_a793u64 << 32) + 0x5bd1_e995;
    #[inline]
    fn shift_mix(v: u64) -> u64 {
        v ^ (v >> 47)
    }
    let len = data.len();
    let len_aligned = len & !7usize;
    let mut hash = seed ^ (len as u64).wrapping_mul(MUL);
    let (chunks, _) = data[..len_aligned].as_chunks::<8>();
    for chunk in chunks {
        let d = u64::from_le_bytes(*chunk);
        let d = shift_mix(d.wrapping_mul(MUL)).wrapping_mul(MUL);
        hash ^= d;
        hash = hash.wrapping_mul(MUL);
    }
    if len & 7 != 0 {
        let mut d: u64 = 0;
        for &b in data[len_aligned..].iter().rev() {
            d = (d << 8).wrapping_add(b as u64);
        }
        hash ^= d;
        hash = hash.wrapping_mul(MUL);
    }
    hash = shift_mix(hash).wrapping_mul(MUL);
    shift_mix(hash)
}

/// `std::hash<std::string>()(s)`
pub fn std_hash_string(s: &[u8]) -> u64 {
    hash_bytes(s, 0xc70f_6907)
}

/// `std::hash<int>()(i)`: identity after conversion to `size_t` (sign-extending).
pub fn std_hash_int(i: i32) -> u64 {
    i as i64 as u64
}

/// libstdc++ `__prime_list` entries reachable from the default growth policy
/// (starting at 13 buckets and at least doubling each time). Verified against
/// libstdc++ by `tests/compat/cxx_oracle.cpp`.
const GROWTH_BUCKETS: &[u64] = &[
    13, 29, 59, 127, 257, 541, 1109, 2357, 5087, 10273, 20753, 42043, 85229, 172933, 351061, 712697, 1447153, 2938679,
    5967347, 12117689, 24607243, 49969847,
];

/// A key hashed as libstdc++'s `std::hash` does for integers (identity).
pub trait CxxHash: Copy + Eq {
    fn cxx_hash(self) -> u64;
}

impl CxxHash for i32 {
    fn cxx_hash(self) -> u64 {
        std_hash_int(self)
    }
}

impl CxxHash for u64 {
    fn cxx_hash(self) -> u64 {
        self
    }
}

const NONE: usize = usize::MAX;
const BEFORE_BEGIN: usize = usize::MAX - 1;

struct Node<K, V> {
    key: K,
    value: V,
    next: usize,
}

/// `std::unordered_map<K, V>` (unique keys, max_load_factor 1.0) with
/// libstdc++'s node layout, so that iteration order matches. Only the
/// operations used by WhatsHap are provided.
pub struct CxxHashMap<K, V> {
    /// Nodes; erased nodes stay in the vector but are unlinked. The list in
    /// iteration order starts at `head` and follows `next`.
    nodes: Vec<Node<K, V>>,
    len: usize,
    /// For each bucket, the node *before* the first node of that bucket
    /// (`BEFORE_BEGIN` for the bucket containing the list head), or `NONE`.
    buckets: Vec<usize>,
    head: usize,
    bucket_count: u64,
    next_resize: u64,
}

impl<K: CxxHash, V> Default for CxxHashMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: CxxHash, V> CxxHashMap<K, V> {
    pub fn new() -> Self {
        CxxHashMap { nodes: Vec::new(), len: 0, buckets: vec![NONE; 1], head: NONE, bucket_count: 1, next_resize: 0 }
    }

    fn bucket_of(&self, k: K) -> usize {
        (k.cxx_hash() % self.bucket_count) as usize
    }

    fn next_of(&self, before: usize) -> usize {
        if before == BEFORE_BEGIN { self.head } else { self.nodes[before].next }
    }

    fn set_next(&mut self, before: usize, n: usize) {
        if before == BEFORE_BEGIN {
            self.head = n;
        } else {
            self.nodes[before].next = n;
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// `_M_find_before_node`: the node before the node with key `k`
    fn find_before(&self, k: K) -> Option<usize> {
        let b = self.bucket_of(k);
        let mut prev = self.buckets[b];
        if prev == NONE {
            return None;
        }
        let mut n = self.next_of(prev);
        loop {
            if self.nodes[n].key == k {
                return Some(prev);
            }
            let next = self.nodes[n].next;
            if next == NONE || self.bucket_of(self.nodes[next].key) != b {
                return None;
            }
            prev = n;
            n = next;
        }
    }

    fn find(&self, k: K) -> Option<usize> {
        self.find_before(k).map(|prev| self.next_of(prev))
    }

    pub fn contains_key(&self, k: K) -> bool {
        self.find(k).is_some()
    }

    pub fn get(&self, k: K) -> Option<&V> {
        self.find(k).map(|n| &self.nodes[n].value)
    }

    /// `_Prime_rehash_policy::_M_need_rehash` with max_load_factor 1.0.
    fn need_rehash(&mut self, n_elt: u64, n_ins: u64) -> Option<u64> {
        if n_elt + n_ins > self.next_resize {
            let min_bkts = std::cmp::max(n_elt + n_ins, if self.next_resize != 0 { 0 } else { 11 }) as f64;
            if min_bkts >= self.bucket_count as f64 {
                let want = std::cmp::max(min_bkts.floor() as u64 + 1, self.bucket_count * 2);
                return Some(self.next_bkt(want));
            }
            self.next_resize = self.bucket_count;
        }
        None
    }

    /// `_Prime_rehash_policy::_M_next_bkt`, restricted to the sizes reachable here.
    fn next_bkt(&mut self, n: u64) -> u64 {
        const FAST_BKT: [u64; 14] = [2, 2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13];
        if (n as usize) < FAST_BKT.len() {
            if n == 0 {
                return 1;
            }
            self.next_resize = FAST_BKT[n as usize];
            return FAST_BKT[n as usize];
        }
        let p = *GROWTH_BUCKETS.iter().find(|&&p| p >= n).expect("unordered_map grew beyond emulated prime table");
        self.next_resize = p;
        p
    }

    fn rehash(&mut self, new_count: u64) {
        let mut new_buckets = vec![NONE; new_count as usize];
        let mut p = self.head;
        self.head = NONE;
        let mut bbegin_bkt = 0usize;
        while p != NONE {
            let next = self.nodes[p].next;
            let bkt = (self.nodes[p].key.cxx_hash() % new_count) as usize;
            if new_buckets[bkt] == NONE {
                self.nodes[p].next = self.head;
                self.head = p;
                new_buckets[bkt] = BEFORE_BEGIN;
                if self.nodes[p].next != NONE {
                    new_buckets[bbegin_bkt] = p;
                }
                bbegin_bkt = bkt;
            } else {
                let before = new_buckets[bkt];
                let after = self.next_of(before);
                self.nodes[p].next = after;
                self.set_next(before, p);
            }
            p = next;
        }
        self.buckets = new_buckets;
        self.bucket_count = new_count;
    }

    /// `_M_insert_unique_node` for a key that is not in the map
    fn insert_new(&mut self, k: K, v: V) -> usize {
        if let Some(n) = self.need_rehash(self.len as u64, 1) {
            self.rehash(n);
        }
        let node = self.nodes.len();
        self.nodes.push(Node { key: k, value: v, next: NONE });
        self.len += 1;
        let bkt = self.bucket_of(k);
        let before = self.buckets[bkt];
        if before != NONE {
            self.nodes[node].next = self.next_of(before);
            self.set_next(before, node);
        } else {
            self.nodes[node].next = self.head;
            self.head = node;
            let next = self.nodes[node].next;
            if next != NONE {
                let b = self.bucket_of(self.nodes[next].key);
                self.buckets[b] = node;
            }
            self.buckets[bkt] = BEFORE_BEGIN;
        }
        node
    }

    /// `insert(value)`: false if the key was already present (value unchanged)
    pub fn insert(&mut self, k: K, v: V) -> bool {
        if self.contains_key(k) {
            return false;
        }
        self.insert_new(k, v);
        true
    }

    /// `operator[]`: the value for `k`, inserting `V::default()` if absent
    pub fn entry_or_default(&mut self, k: K) -> &mut V
    where
        V: Default,
    {
        let n = match self.find(k) {
            Some(n) => n,
            None => self.insert_new(k, V::default()),
        };
        &mut self.nodes[n].value
    }

    /// `erase(key)`
    pub fn remove(&mut self, k: K) -> bool {
        let Some(prev) = self.find_before(k) else { return false };
        let n = self.next_of(prev);
        let bkt = self.bucket_of(k);
        let next = self.nodes[n].next;
        let next_bkt = if next != NONE { self.bucket_of(self.nodes[next].key) } else { 0 };
        if prev == self.buckets[bkt] {
            // _M_remove_bucket_begin
            if next == NONE || next_bkt != bkt {
                if next != NONE {
                    self.buckets[next_bkt] = self.buckets[bkt];
                }
                if self.buckets[bkt] == BEFORE_BEGIN {
                    self.head = next;
                }
                self.buckets[bkt] = NONE;
            }
        } else if next != NONE && next_bkt != bkt {
            self.buckets[next_bkt] = prev;
        }
        self.set_next(prev, next);
        self.len -= 1;
        true
    }

    /// `clear()` keeps the bucket array and the rehash policy state.
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.len = 0;
        self.head = NONE;
        for b in self.buckets.iter_mut() {
            *b = NONE;
        }
    }

    /// Entries in iteration order
    pub fn iter(&self) -> impl Iterator<Item = (K, &V)> + '_ {
        let mut n = self.head;
        std::iter::from_fn(move || {
            if n == NONE {
                None
            } else {
                let node = &self.nodes[n];
                n = node.next;
                Some((node.key, &node.value))
            }
        })
    }

    pub fn bucket_count(&self) -> u64 {
        self.bucket_count
    }
}

impl<K: CxxHash, V: Clone> CxxHashMap<K, V> {
    /// The range constructor `unordered_map(other.begin(), other.end())`:
    /// for unique keys, an empty map into which the entries are inserted in
    /// iteration order.
    pub fn copied(&self) -> Self {
        let mut m = CxxHashMap::new();
        for (k, v) in self.iter() {
            m.insert(k, v.clone());
        }
        m
    }
}

/// `std::unordered_set<int>` with libstdc++'s node layout, so that iteration
/// order matches.
#[derive(Default)]
pub struct CxxUnorderedSetI32(CxxHashMap<i32, ()>);

impl CxxUnorderedSetI32 {
    pub fn new() -> Self {
        CxxUnorderedSetI32(CxxHashMap::new())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn contains(&self, v: i32) -> bool {
        self.0.contains_key(v)
    }

    pub fn insert(&mut self, v: i32) -> bool {
        self.0.insert(v, ())
    }

    pub fn clear(&mut self) {
        self.0.clear()
    }

    pub fn iter(&self) -> impl Iterator<Item = i32> + '_ {
        self.0.iter().map(|(k, _)| k)
    }

    pub fn bucket_count(&self) -> u64 {
        self.0.bucket_count()
    }
}

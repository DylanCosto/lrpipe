//! Bit-exact emulation of CPython 3.12's `set` for non-negative integers.
//!
//! WhatsHap's read selection iterates over Python sets of read indices, and
//! the iteration order decides how ties in the priority queue are broken. To
//! select the same reads, we reproduce `Objects/setobject.c`: the open
//! addressing scheme (linear probing + perturbation), dummy entries left by
//! deletions, and the resize policy. Integers hash to themselves.

const MINSIZE: usize = 8;
const LINEAR_PROBES: usize = 9;
const PERTURB_SHIFT: u32 = 5;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Slot {
    Empty,
    Dummy,
    Active(u64),
}

#[derive(Clone, Debug)]
pub struct PySet {
    table: Vec<Slot>,
    mask: usize,
    fill: usize,
    used: usize,
}

impl Default for PySet {
    fn default() -> Self {
        Self::new()
    }
}

impl PySet {
    pub fn new() -> Self {
        PySet { table: vec![Slot::Empty; MINSIZE], mask: MINSIZE - 1, fill: 0, used: 0 }
    }

    /// `set(iterable)` for an iterable that is not a set or dict.
    pub fn from_iter_adds(iter: impl IntoIterator<Item = usize>) -> Self {
        let mut s = PySet::new();
        for v in iter {
            s.add(v);
        }
        s
    }

    /// `set(other_set)` / `other_set.copy()`
    pub fn copy_of(other: &PySet) -> Self {
        let mut s = PySet::new();
        s.merge(other);
        s
    }

    pub fn len(&self) -> usize {
        self.used
    }

    pub fn is_empty(&self) -> bool {
        self.used == 0
    }

    /// `set_insert_clean`
    fn insert_clean(table: &mut [Slot], mask: usize, key: u64) {
        let mut perturb = key;
        let mut i = (key as usize) & mask;
        loop {
            if table[i] == Slot::Empty {
                table[i] = Slot::Active(key);
                return;
            }
            if i + LINEAR_PROBES <= mask {
                for j in 1..=LINEAR_PROBES {
                    if table[i + j] == Slot::Empty {
                        table[i + j] = Slot::Active(key);
                        return;
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }

    /// `set_table_resize`
    fn resize(&mut self, minused: usize) {
        let mut newsize = MINSIZE;
        while newsize <= minused {
            newsize <<= 1;
        }
        // Staying in the inline "smalltable" without dummies is a no-op.
        if newsize == MINSIZE && self.mask == MINSIZE - 1 && self.fill == self.used {
            return;
        }
        let old = std::mem::replace(&mut self.table, vec![Slot::Empty; newsize]);
        self.mask = newsize - 1;
        // Dummies are not copied, so fill becomes used (the "fill == used" branch
        // of CPython is equivalent since then no dummies exist).
        self.fill = self.used;
        for slot in old {
            if let Slot::Active(k) = slot {
                Self::insert_clean(&mut self.table, self.mask, k);
            }
        }
    }

    /// `set_add_entry`
    pub fn add(&mut self, key: usize) {
        let key = key as u64;
        let mask = self.mask;
        let mut i = (key as usize) & mask;
        let mut freeslot: Option<usize> = None;
        let mut perturb = key;
        loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            for j in 0..=probes {
                let idx = i + j;
                match self.table[idx] {
                    Slot::Empty => {
                        if let Some(f) = freeslot {
                            self.used += 1;
                            self.table[f] = Slot::Active(key);
                            return;
                        }
                        self.fill += 1;
                        self.used += 1;
                        self.table[idx] = Slot::Active(key);
                        if self.fill * 5 < mask * 3 {
                            return;
                        }
                        let m = if self.used > 50000 { self.used * 2 } else { self.used * 4 };
                        self.resize(m);
                        return;
                    }
                    Slot::Active(k) => {
                        if k == key {
                            return;
                        }
                    }
                    // CPython 3.12 remembers the *last* dummy seen, not the first.
                    Slot::Dummy => freeslot = Some(idx),
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }

    /// `set_lookkey`: index of the slot holding `key`, if present.
    fn lookup(&self, key: u64) -> Option<usize> {
        let mask = self.mask;
        let mut i = (key as usize) & mask;
        let mut perturb = key;
        loop {
            let probes = if i + LINEAR_PROBES <= mask { LINEAR_PROBES } else { 0 };
            for j in 0..=probes {
                match self.table[i + j] {
                    Slot::Empty => return None,
                    Slot::Active(k) if k == key => return Some(i + j),
                    _ => {}
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }

    pub fn contains(&self, key: usize) -> bool {
        self.lookup(key as u64).is_some()
    }

    /// `set_discard_entry`; returns whether the key was present.
    pub fn discard(&mut self, key: usize) -> bool {
        match self.lookup(key as u64) {
            Some(idx) => {
                self.table[idx] = Slot::Dummy;
                self.used -= 1;
                true
            }
            None => false,
        }
    }

    /// `set.remove(key)`; panics like a KeyError if missing.
    pub fn remove(&mut self, key: usize) {
        assert!(self.discard(key), "KeyError: {key}");
    }

    /// `set.update(iterable)` for a non-set iterable (e.g. a list).
    pub fn update_iter(&mut self, iter: impl IntoIterator<Item = usize>) {
        for v in iter {
            self.add(v);
        }
    }

    /// `set_merge`, i.e. `self.update(other_set)`.
    pub fn merge(&mut self, other: &PySet) {
        if other.used == 0 {
            return;
        }
        if (self.fill + other.used) * 5 >= self.mask * 3 {
            self.resize((self.used + other.used) * 2);
        }
        if self.fill == 0 && self.mask == other.mask && other.fill == other.used {
            self.table.copy_from_slice(&other.table);
            self.fill = other.fill;
            self.used = other.used;
            return;
        }
        if self.fill == 0 {
            self.fill = other.used;
            self.used = other.used;
            for slot in &other.table {
                if let Slot::Active(k) = *slot {
                    Self::insert_clean(&mut self.table, self.mask, k);
                }
            }
            return;
        }
        for slot in &other.table {
            if let Slot::Active(k) = *slot {
                self.add(k as usize);
            }
        }
    }

    /// `set_intersection` (only used internally by difference_update).
    fn intersection(&self, other: &PySet) -> PySet {
        let (a, b) = if other.len() > self.len() { (other, self) } else { (self, other) };
        // CPython: if len(other) > len(so): swap, then iterate over the smaller one
        // (`b` here) and keep elements contained in the larger one.
        let mut result = PySet::new();
        for k in b.iter() {
            if a.contains(k) {
                result.add(k);
            }
        }
        result
    }

    /// `set_difference_update_internal`, i.e. `self -= other`.
    pub fn difference_update(&mut self, other: &PySet) {
        let reduced;
        let other = if (other.len() >> 3) > self.len() {
            reduced = self.intersection(other);
            &reduced
        } else {
            other
        };
        let keys: Vec<usize> = other.iter().collect();
        for k in keys {
            self.discard(k);
        }
        if self.fill - self.used <= self.mask / 4 {
            return;
        }
        let m = if self.used > 50000 { self.used * 2 } else { self.used * 4 };
        self.resize(m);
    }

    /// `set.difference(other_set)`
    pub fn difference(&self, other: &PySet) -> PySet {
        if (self.len() >> 2) > other.len() {
            let mut result = PySet::copy_of(self);
            result.difference_update(other);
            return result;
        }
        let mut result = PySet::new();
        for k in self.iter() {
            if !other.contains(k) {
                result.add(k);
            }
        }
        result
    }

    /// Iterate in CPython's order (table slot order).
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.table.iter().filter_map(|s| match *s {
            Slot::Active(k) => Some(k as usize),
            _ => None,
        })
    }
}

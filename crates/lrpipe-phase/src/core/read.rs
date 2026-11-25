//! Port of `read.{h,cpp}`, `readset.{h,cpp}` and `entry.{h,cpp}`.

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::hash::{BuildHasher, RandomState};

use anyhow::{Result, bail};
use hashbrown::HashTable;

use crate::compat::cxx::{std_hash_int, std_hash_string};

/// Allele codes as used by `Entry::allele_t`.
pub const REF_ALLELE: i32 = 0;
pub const ALT_ALLELE: i32 = 1;
pub const BLANK: i32 = 2;
pub const EQUAL_SCORES: i32 = 3;

/// One cell of the read/variant matrix (`Entry` in C++).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Entry {
    pub read_id: u32,
    pub allele: i32,
    pub quality: u32,
}

/// A variant observed on a read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadVariant {
    pub position: i32,
    pub allele: i32,
    pub quality: u32,
}

#[derive(Clone, Debug)]
pub struct Read {
    pub name: String,
    pub mapqs: Vec<i32>,
    pub source_id: i32,
    pub sample_id: i32,
    pub id: i32,
    pub reference_start: i32,
    pub chromosome: String,
    pub sub_alignment_id: String,
    pub is_supplementary: bool,
    pub is_reverse: bool,
    pub reference_end: i32,
    pub bx_tag: String,
    pub hp_tag: i32,
    pub ps_tag: i32,
    pub variants: Vec<ReadVariant>,
}

impl Read {
    pub fn new(name: &str, mapq: i32, source_id: i32, sample_id: i32) -> Self {
        Read {
            name: name.to_string(),
            mapqs: vec![mapq],
            source_id,
            sample_id,
            id: -1,
            reference_start: -1,
            chromosome: String::new(),
            sub_alignment_id: String::new(),
            is_supplementary: false,
            is_reverse: false,
            reference_end: -1,
            bx_tag: String::new(),
            hp_tag: -1,
            ps_tag: -1,
            variants: Vec::new(),
        }
    }

    pub fn add_variant(&mut self, position: i32, allele: i32, quality: i32) {
        // Entry stores the quality as unsigned int.
        self.variants.push(ReadVariant { position, allele, quality: quality as u32 });
    }

    pub fn len(&self) -> usize {
        self.variants.len()
    }

    pub fn is_empty(&self) -> bool {
        self.variants.is_empty()
    }

    pub fn sort(&mut self) -> Result<()> {
        self.variants.sort_by_key(|v| v.position);
        for w in self.variants.windows(2) {
            if w[0].position == w[1].position {
                bail!("Duplicate variant in read {} at position {}", self.name, w[1].position);
            }
        }
        Ok(())
    }

    pub fn is_sorted(&self) -> bool {
        self.variants.windows(2).all(|w| w[0].position < w[1].position)
    }

    pub fn first_position(&self) -> i32 {
        self.variants.first().expect("No variants present").position
    }

    pub fn last_position(&self) -> i32 {
        self.variants.last().expect("No variants present").position
    }

    pub fn entry(&self, i: usize) -> Entry {
        let v = &self.variants[i];
        Entry { read_id: self.id as u32, allele: v.allele, quality: v.quality }
    }

    pub fn set_id(&mut self, id: i32) {
        self.id = id;
    }

    /// Approximate memory used by the read, including its heap allocations.
    pub fn memory_size(&self) -> usize {
        std::mem::size_of::<Read>()
            + self.name.capacity()
            + self.mapqs.capacity() * std::mem::size_of::<i32>()
            + self.chromosome.capacity()
            + self.sub_alignment_id.capacity()
            + self.bx_tag.capacity()
            + self.variants.capacity() * std::mem::size_of::<ReadVariant>()
    }
}

/// `ReadSet::read_comparator_t`
fn read_cmp(r1: &Read, r2: &Read) -> Ordering {
    if !r1.is_empty() || !r2.is_empty() {
        if r1.is_empty() {
            return Ordering::Less;
        }
        if r2.is_empty() {
            return Ordering::Greater;
        }
        let (p1, p2) = (r1.first_position(), r2.first_position());
        if p1 != p2 {
            return p1.cmp(&p2);
        }
    }
    let h1 = std_hash_string(r1.name.as_bytes()) ^ std_hash_int(r1.source_id);
    let h2 = std_hash_string(r2.name.as_bytes()) ^ std_hash_int(r2.source_id);
    if h1 != h2 {
        return h1.cmp(&h2);
    }
    match r1.name.as_bytes().cmp(r2.name.as_bytes()) {
        Ordering::Equal => r1.source_id.cmp(&r2.source_id),
        o => o,
    }
}

#[derive(Clone, Debug, Default)]
pub struct ReadSet {
    reads: Vec<Read>,
    /// Indices into `reads`, by (name, source_id); unlike a map keyed by the
    /// pair it needs no copy of the names.
    by_name: HashTable<usize>,
    hasher: RandomState,
}

impl ReadSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(capacity: usize) -> Self {
        ReadSet { reads: Vec::with_capacity(capacity), by_name: HashTable::with_capacity(capacity), ..Self::default() }
    }

    fn name_hash(&self, name: &str, source_id: i32) -> u64 {
        self.hasher.hash_one((name, source_id))
    }

    fn find(&self, name: &str, source_id: i32) -> Option<usize> {
        let reads = &self.reads;
        self.by_name
            .find(self.name_hash(name, source_id), |&i| reads[i].source_id == source_id && reads[i].name == name)
            .copied()
    }

    /// Adds `reads[i]` to the name index.
    fn index(&mut self, i: usize) {
        let hash = self.name_hash(&self.reads[i].name, self.reads[i].source_id);
        let (reads, hasher) = (&self.reads, &self.hasher);
        self.by_name.insert_unique(hash, i, |&j| hasher.hash_one((reads[j].name.as_str(), reads[j].source_id)));
    }

    pub fn add(&mut self, read: Read) -> Result<()> {
        if self.find(&read.name, read.source_id).is_some() {
            bail!("ReadSet::add: duplicate read name.");
        }
        self.reads.push(read);
        self.index(self.reads.len() - 1);
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.reads.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reads.is_empty()
    }

    pub fn get(&self, i: usize) -> &Read {
        &self.reads[i]
    }

    pub fn get_mut(&mut self, i: usize) -> &mut Read {
        &mut self.reads[i]
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Read> {
        self.reads.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Read> {
        self.reads.iter_mut()
    }

    pub fn get_by_name(&self, name: &str, source_id: i32) -> Option<&Read> {
        self.find(name, source_id).map(|i| &self.reads[i])
    }

    /// Sort reads by position of first variant; ties broken as in C++.
    pub fn sort(&mut self) {
        // The comparator is a strict total order on distinct (name, source_id)
        // pairs, so a stable sort yields the same result as std::sort.
        self.reads.sort_by(read_cmp);
        self.by_name.clear();
        for i in 0..self.reads.len() {
            self.index(i);
        }
    }

    /// Sorted list of all variant positions.
    pub fn get_positions(&self) -> Vec<i32> {
        let mut s: Vec<i32> = self.reads.iter().flat_map(|r| r.variants.iter().map(|v| v.position)).collect();
        s.sort_unstable();
        s.dedup();
        s
    }

    /// `ReadSet::subset`: copies of the reads at the given indices, in
    /// ascending index order (the C++ code collects indices in a std::set).
    pub fn subset(&self, indices: impl IntoIterator<Item = usize>) -> ReadSet {
        let idx: BTreeSet<usize> = indices.into_iter().collect();
        let mut result = ReadSet::new();
        for i in idx {
            result.add(self.reads[i].clone()).expect("subset of a valid ReadSet");
        }
        result
    }

    /// Like `subset`, but moves the reads instead of copying them.
    pub fn into_subset(self, indices: impl IntoIterator<Item = usize>) -> ReadSet {
        let idx: BTreeSet<usize> = indices.into_iter().collect();
        let mut reads: Vec<Option<Read>> = self.reads.into_iter().map(Some).collect();
        let mut result = ReadSet::with_capacity(idx.len());
        for i in idx {
            result.add(reads[i].take().expect("index listed once")).expect("subset of a valid ReadSet");
        }
        result
    }

    pub fn reassign_read_ids(&mut self) {
        for (i, r) in self.reads.iter_mut().enumerate() {
            r.set_id(i as i32);
        }
    }

    pub fn into_reads(self) -> Vec<Read> {
        self.reads
    }
}

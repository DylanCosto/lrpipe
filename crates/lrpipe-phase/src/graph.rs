//! Port of `graph.py` (`ComponentFinder`).

use std::collections::HashMap;
use std::hash::Hash;

/// Union-find where the representative of a set is always its minimum value.
pub struct ComponentFinder<T: Clone + Ord + Hash> {
    index: HashMap<T, usize>,
    values: Vec<T>,
    parent: Vec<usize>,
}

impl<T: Clone + Ord + Hash> ComponentFinder<T> {
    pub fn new(values: impl IntoIterator<Item = T>) -> Self {
        let mut cf = ComponentFinder { index: HashMap::new(), values: Vec::new(), parent: Vec::new() };
        for v in values {
            if !cf.index.contains_key(&v) {
                let i = cf.values.len();
                cf.index.insert(v.clone(), i);
                cf.values.push(v);
                cf.parent.push(i);
            }
        }
        cf
    }

    fn find_root(&mut self, i: usize) -> usize {
        let mut root = i;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut node = i;
        while self.parent[node] != node {
            let next = self.parent[node];
            self.parent[node] = root;
            node = next;
        }
        root
    }

    pub fn merge(&mut self, x: T, y: T) {
        assert!(x != y);
        let xr = self.find_root(self.index[&x]);
        let yr = self.find_root(self.index[&y]);
        if xr == yr {
            return;
        }
        if self.values[xr] < self.values[yr] {
            self.parent[yr] = xr;
        } else {
            self.parent[xr] = yr;
        }
    }

    pub fn find(&mut self, x: T) -> T {
        let r = self.find_root(self.index[&x]);
        self.values[r].clone()
    }
}

//! Grouping leads into SV clusters, as in Sniffles2's `cluster.py`.

use std::collections::HashMap;
use std::sync::Arc;

use crate::config::Config;
use crate::pymath;
use crate::reads::{Lead, SvType};

/// All leads of one SV type on a contig, binned by `ref_start`.
///
/// Clusters refer to leads by index, so that merging a read's leads updates the shared lead
/// (Sniffles2 mutates the lead objects, and DUP genotyping later reads the INS leads).
#[derive(Debug, Default)]
pub struct Leadtab {
    pub leads: Vec<Lead>,
    /// (bin, lead indices in read order), sorted by bin.
    pub bins: Vec<(i64, Vec<usize>)>,
    /// Leads per haplotype (HP 0/1/2) in each bin, parallel to `bins`.
    pub hap_counts: Vec<[u16; 3]>,
}

impl Leadtab {
    /// Builds the table from (bin, lead) pairs in read order.
    pub fn new(entries: Vec<(i64, Lead)>, consensus_max_reads_bin: usize) -> Self {
        let mut order: Vec<usize> = (0..entries.len()).collect();
        order.sort_by_key(|&i| entries[i].0);
        let mut bins: Vec<(i64, Vec<usize>)> = Vec::new();
        for &i in &order {
            let bin = entries[i].0;
            match bins.last_mut() {
                Some((b, v)) if *b == bin => v.push(i),
                _ => bins.push((bin, vec![i])),
            }
        }
        let mut leads: Vec<Lead> = entries.into_iter().map(|(_, l)| l).collect();
        let mut hap_counts = Vec::with_capacity(bins.len());
        for (_, idx) in &bins {
            let mut hc = [0u16; 3];
            for (k, &i) in idx.iter().enumerate() {
                // Sniffles2 drops the sequence of every lead beyond the first few in a bin.
                if k >= consensus_max_reads_bin {
                    leads[i].seq = None;
                }
                // array('H') raises OverflowError (caught) instead of wrapping.
                let h = &mut hc[leads[i].hap as usize];
                *h = h.saturating_add(1);
            }
            hap_counts.push(hc);
        }
        Leadtab { leads, bins, hap_counts }
    }

    pub fn bin(&self, bin: i64) -> Option<&[usize]> {
        self.bins.binary_search_by_key(&bin, |(b, _)| *b).ok().map(|k| self.bins[k].1.as_slice())
    }
}

/// Reference (all reads) haplotype counts per bin, from read spans.
#[derive(Debug, Default)]
pub struct HapRef {
    /// Per HP: sorted first bins and sorted end bins of reads.
    starts: [Vec<i64>; 3],
    ends: [Vec<i64>; 3],
}

impl HapRef {
    pub fn new(hap_bins: &[(i64, i64, u8)]) -> Self {
        let mut h = HapRef::default();
        for &(s, e, hp) in hap_bins {
            if e > s {
                h.starts[hp as usize].push(s);
                h.ends[hp as usize].push(e);
            }
        }
        for v in h.starts.iter_mut().chain(h.ends.iter_mut()) {
            v.sort_unstable();
        }
        h
    }

    /// Reads per HP whose bins [first, end) contain `bin` (`leadhapcount["REF"][bin]`).
    pub fn counts(&self, bin: i64) -> [i64; 3] {
        std::array::from_fn(|hp| {
            let started = self.starts[hp].partition_point(|&s| s <= bin);
            let ended = self.ends[hp].partition_point(|&e| e <= bin);
            // A bin entry is a u16 array in Sniffles2 that stops counting at 65535.
            ((started - ended) as i64).min(u16::MAX as i64)
        })
    }
}

#[derive(Debug, Clone)]
pub struct Cluster {
    pub svtype: SvType,
    pub start: i64,
    pub end: i64,
    pub leads: Vec<usize>,
    pub repeat: bool,
    /// Soft-clip based insertion leads (INS only).
    pub leads_long: Option<Arc<Vec<usize>>>,
    /// (SV HP0, HP1, HP2, ref HP0, HP1, HP2)
    pub hap_counts: [i64; 6],
    pub sa_counts: (usize, f64),
    mean_svlen: f64,
    stdev_start: f64,
}

impl Cluster {
    fn compute_metrics(&mut self, tab: &Leadtab) {
        let len = self.leads.len();
        let n = len.min(100);
        if n == 0 {
            self.mean_svlen = 0.0;
            self.stdev_start = 0.0;
            return;
        }
        let step = len / n;
        if n > 1 {
            // Sums every step-th lead (possibly more than n of them) but divides by n, as Sniffles2.
            let sampled: Vec<&Lead> = (0..len).step_by(step).map(|i| &tab.leads[self.leads[i]]).collect();
            let sum: i64 = sampled.iter().map(|l| l.svlen.expect("cluster leads have a length")).sum();
            self.mean_svlen = sum as f64 / n as f64;
            let starts: Vec<i64> = sampled.iter().map(|l| l.ref_start).collect();
            self.stdev_start = pymath::stdev(&starts).expect("at least two sampled leads");
        } else {
            self.mean_svlen = tab.leads[self.leads[0]].svlen.expect("cluster leads have a length") as f64;
            self.stdev_start = 0.0;
        }
    }

    /// A cluster like this one but with other leads (without copying this one's leads).
    pub(crate) fn with_leads(&self, leads: Vec<usize>, leads_long: Option<Arc<Vec<usize>>>) -> Cluster {
        Cluster {
            svtype: self.svtype,
            start: self.start,
            end: self.end,
            leads,
            repeat: self.repeat,
            leads_long,
            hap_counts: self.hap_counts,
            sa_counts: self.sa_counts,
            mean_svlen: self.mean_svlen,
            stdev_start: self.stdev_start,
        }
    }

    /// A cluster of the given leads, for tests.
    #[cfg(test)]
    pub(crate) fn for_test(svtype: SvType, leads: Vec<usize>) -> Cluster {
        Cluster {
            svtype,
            start: 0,
            end: 100,
            leads,
            repeat: false,
            leads_long: None,
            hap_counts: [0; 6],
            sa_counts: (0, 0.0),
            mean_svlen: 0.0,
            stdev_start: 0.0,
        }
    }

    /// `Cluster.get_sa_count`: leads from supplementary alignments, and their share of all leads.
    pub fn compute_sa_counts(&mut self, tab: &Leadtab) {
        let long = self.leads_long.as_deref().map_or(&[][..], Vec::as_slice);
        let all = self.leads.iter().chain(long);
        let total = self.leads.len() + long.len();
        let sa = all.filter(|&&i| tab.leads[i].is_sa).count();
        self.sa_counts = (sa, sa as f64 / total as f64);
    }
}

/// Clusters one SV type's leads (`cluster.resolve`). Mutates leads when merging a read's leads.
pub fn resolve(
    svtype: SvType,
    tab: &mut Leadtab,
    hap_ref: &HapRef,
    tr: Option<&[(i64, i64)]>,
    cfg: &Config,
) -> Vec<Cluster> {
    if tab.bins.is_empty() {
        return Vec::new();
    }
    let mut clusters = seed_clusters(svtype, tab, hap_ref, tr, cfg);
    merge_neighbours(&mut clusters, svtype, tab, cfg);

    let mut out = Vec::new();
    for mut cluster in clusters {
        if cluster.leads.is_empty() {
            continue;
        }
        if svtype == SvType::Bnd {
            if cfg.args.dev_no_resplit {
                out.push(cluster);
            } else {
                out.extend(resplit_bnd(cluster, cfg.args.cluster_merge_bnd, tab));
            }
            continue;
        }
        if matches!(svtype, SvType::Ins | SvType::Del) {
            let threshold = if cluster.repeat { -1 } else { cfg.cluster_merge_pos };
            merge_inner(&mut cluster, threshold, tab);
        }
        if !cfg.args.dev_no_resplit_repeat && !cfg.args.dev_no_resplit {
            out.extend(resplit(cluster, cfg.args.cluster_resplit_binsize, cfg.minsvlen, cfg.cluster_merge_len, tab));
        } else {
            out.push(cluster);
        }
    }
    out
}

/// One cluster per leadtab bin with enough leads, flagged when inside a tandem repeat.
fn seed_clusters(
    svtype: SvType,
    tab: &Leadtab,
    hap_ref: &HapRef,
    tr: Option<&[(i64, i64)]>,
    cfg: &Config,
) -> Vec<Cluster> {
    let tr = tr.filter(|t| !t.is_empty());
    let mut tr_index = 0;
    let mut clusters: Vec<Cluster> = Vec::new();
    for (k, (seed, idx)) in tab.bins.iter().enumerate() {
        let seed = *seed;
        let mut within_tr = false;
        if let Some(tr) = tr
            && tr_index < tr.len()
        {
            while tr[tr_index].1 < seed && tr_index + 1 < tr.len() {
                tr_index += 1;
            }
            let (tr_start, tr_end) = tr[tr_index];
            within_tr = tr_start < seed && seed < tr_end;
        }
        let (leads, leads_long) = if svtype == SvType::Ins {
            let (with_len, long): (Vec<usize>, Vec<usize>) = idx.iter().partition(|&&i| tab.leads[i].svlen.is_some());
            (with_len, Some(Arc::new(long)))
        } else {
            (idx.clone(), None)
        };
        let sv = tab.hap_counts[k];
        let r = hap_ref.counts(seed);
        let hap_counts = [sv[0] as i64, sv[1] as i64, sv[2] as i64, r[0], r[1], r[2]];
        if leads.len() >= cfg.dev_min_leads_cluster {
            let mut c = Cluster {
                svtype,
                start: seed,
                end: seed + cfg.args.cluster_binsize,
                leads,
                repeat: within_tr || cfg.args.repeat,
                leads_long,
                hap_counts,
                sa_counts: (0, 0.0),
                mean_svlen: 0.0,
                stdev_start: 0.0,
            };
            c.compute_metrics(tab);
            clusters.push(c);
        }
    }
    clusters
}

/// Merges neighbouring clusters that are close relative to their spread (or, in repeats, their
/// lengths). Steps back after each merge exactly as Sniffles2's loop does.
fn merge_neighbours(clusters: &mut Vec<Cluster>, svtype: SvType, tab: &Leadtab, cfg: &Config) {
    let mut i = 0;
    while i + 1 < clusters.len() {
        let (curr, next) = (&clusters[i], &clusters[i + 1]);
        let inner_dist = (next.start - curr.end) as f64;
        let outer_dist = (next.end - curr.start) as f64;
        let repeat_reach = cfg
            .args
            .cluster_repeat_h_max
            .min((curr.mean_svlen.abs() + next.mean_svlen.abs()) * cfg.args.cluster_repeat_h);
        let merge = inner_dist <= curr.stdev_start.min(next.stdev_start) * cfg.args.cluster_r
            || ((cfg.args.repeat || curr.repeat || next.repeat) && outer_dist <= repeat_reach)
            || (svtype == SvType::Bnd && inner_dist <= cfg.args.cluster_merge_bnd as f64);
        if merge {
            let next = clusters.remove(i + 1);
            let curr = &mut clusters[i];
            curr.leads.extend(next.leads);
            if svtype == SvType::Ins {
                // Unshared until resplit, so this extends in place.
                let long = Arc::make_mut(curr.leads_long.get_or_insert_with(Default::default));
                long.extend(next.leads_long.as_deref().into_iter().flatten());
            }
            curr.end = next.end;
            curr.repeat = curr.repeat || next.repeat;
            curr.compute_metrics(tab);
            i = i.saturating_sub(2);
        }
        i += 1;
    }
}

/// `merge_inner`: joins consecutive leads of the same read (summing lengths).
fn merge_inner(cluster: &mut Cluster, threshold: i64, tab: &mut Leadtab) {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut by_name: HashMap<&[u8], usize> = HashMap::new();
    for &i in &cluster.leads {
        let g = *by_name.entry(&tab.leads[i].qname).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[g].push(i);
    }
    let mut merged = Vec::with_capacity(cluster.leads.len());
    for mut group in groups {
        group.sort_by_key(|&i| tab.leads[i].ref_start);
        let mut curr = group[0];
        let first = &tab.leads[curr];
        let (mut last_ref_end, mut last_qry_end) = (first.ref_end, first.qry_end);
        let (mut last_ref_start, mut last_qry_start) = (first.ref_start, first.qry_start);
        for &to in &group[1..] {
            let t = &tab.leads[to];
            let c = &tab.leads[curr];
            // Python precedence: `t == -1 or (near_ref and near_qry and same_strand)`.
            let merge = threshold == -1
                || (((t.ref_start - last_ref_end).abs() < threshold
                    || (t.ref_start - last_ref_start).abs() < threshold)
                    && ((t.qry_start - last_qry_end).abs() < threshold
                        || (t.qry_start - last_qry_start).abs() < threshold)
                    && c.reverse == t.reverse);
            let (t_ref_start, t_ref_end, t_qry_start, t_qry_end) = (t.ref_start, t.ref_end, t.qry_start, t.qry_end);
            if merge {
                let t_svlen = t.svlen.expect("INS/DEL leads have a length");
                let t_seq = t.seq.clone();
                let c = &mut tab.leads[curr];
                c.svlen = Some(c.svlen.expect("INS/DEL leads have a length") + t_svlen);
                match (&mut c.seq, t_seq) {
                    (Some(s), Some(ts)) => s.extend_from_slice(&ts),
                    (s, _) => *s = None,
                }
            } else {
                merged.push(curr);
                curr = to;
            }
            last_ref_end = t_ref_end;
            last_qry_end = t_qry_end;
            last_ref_start = t_ref_start;
            last_qry_start = t_qry_start;
        }
        merged.push(curr);
    }
    cluster.leads = merged;
}

/// `resplit`: splits a cluster by SV length.
fn resplit(
    cluster: Cluster,
    binsize: i64,
    merge_threshold_min: i64,
    merge_threshold_frac: f64,
    tab: &Leadtab,
) -> Vec<Cluster> {
    let mut keys: Vec<i64> = Vec::new();
    let mut bins: HashMap<i64, Vec<usize>> = HashMap::new();
    for &i in &cluster.leads {
        let svlen = tab.leads[i].svlen.expect("resplit leads have a length");
        let bin = svlen.abs() / binsize * binsize;
        bins.entry(bin)
            .or_insert_with(|| {
                keys.push(bin);
                Vec::new()
            })
            .push(i);
    }
    let mut nc = keys;
    nc.sort_unstable();
    // Python indexing: i - 1 may be -1, meaning the last element.
    let mut i: usize = 1;
    while nc.len() > 1 && i < nc.len() {
        let prev = if i == 0 { nc.len() - 1 } else { i - 1 };
        let (last, curr) = (nc[prev], nc[i]);
        let threshold = (merge_threshold_min as f64).max(curr.min(last) as f64 * merge_threshold_frac);
        if ((curr - last).abs() as f64) <= threshold {
            // The merged-away bin is never used again, so move its leads instead of copying.
            let moved = bins.remove(&last).unwrap_or_default();
            bins.get_mut(&curr).expect("bin exists").extend(moved);
            nc.remove(prev);
            i = i.saturating_sub(2);
        } else {
            i += 1;
        }
    }
    let leads_long = cluster.leads_long.clone();
    nc.into_iter().map(|k| cluster.with_leads(bins.remove(&k).expect("bin exists"), leads_long.clone())).collect()
}

/// `resplit_bnd`: splits a BND cluster by mate contig, orientation and mate position.
fn resplit_bnd(cluster: Cluster, merge_threshold: i64, tab: &Leadtab) -> Vec<Cluster> {
    if cluster.leads.len() <= 1 {
        return vec![cluster];
    }
    type Ident = (Arc<str>, bool);
    let mut idents: Vec<Ident> = Vec::new();
    let mut groups: HashMap<Ident, HashMap<i64, Vec<usize>>> = HashMap::new();
    for &i in &cluster.leads {
        let b = tab.leads[i].bnd.as_ref().expect("BND leads have mate info");
        let ident = (b.mate_contig.clone(), b.is_first);
        let pos_bin = if merge_threshold > 0 { b.mate_ref_start / merge_threshold * merge_threshold } else { 0 };
        groups
            .entry(ident.clone())
            .or_insert_with(|| {
                idents.push(ident);
                HashMap::new()
            })
            .entry(pos_bin)
            .or_default()
            .push(i);
    }
    let mut out = Vec::new();
    let new_cluster = |leads: Vec<usize>| cluster.with_leads(leads, None);
    for ident in idents {
        let pos_bins = &groups[&ident];
        let mut bins: Vec<i64> = pos_bins.keys().copied().collect();
        bins.sort_unstable();
        let mut curr: Vec<usize> = pos_bins[&bins[0]].clone();
        let mut last_bin = bins[0];
        for &b in &bins[1..] {
            if b - last_bin <= merge_threshold {
                curr.extend_from_slice(&pos_bins[&b]);
            } else {
                if !curr.is_empty() {
                    out.push(new_cluster(std::mem::take(&mut curr)));
                }
                curr = pos_bins[&b].clone();
            }
            last_bin = b;
        }
        if !curr.is_empty() {
            out.push(new_cluster(curr));
        }
    }
    out
}

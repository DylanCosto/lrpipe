//! Splitting a cluster into two alleles (lrpipe-sv improvement, off by default).
//!
//! Sniffles2 reports one call per cluster. In tandem repeats the two haplotypes often carry
//! insertions or deletions of different lengths (for example +300 and +330 bp), which
//! Sniffles2's length re-split merges unless they differ by more than max(50 bp, 22%). This
//! module looks for two length modes among a cluster's leads and, when both are well
//! supported and separated by a clear gap, returns two clusters.

use std::collections::HashSet;

use crate::cluster::{Cluster, Leadtab};
use crate::config::Config;

/// Splits an INS/DEL cluster into two by lead length when the lengths are clearly bimodal.
pub fn split(cluster: Cluster, tab: &Leadtab, cfg: &Config) -> Vec<Cluster> {
    let mut lens: Vec<(i64, usize)> =
        cluster.leads.iter().map(|&i| (tab.leads[i].svlen.expect("INS/DEL leads have a length").abs(), i)).collect();
    let n = lens.len();
    let min_group = (cfg.args.split_min_support as usize).max((cfg.args.split_min_frac * n as f64).ceil() as usize);
    if n < 2 * min_group || min_group == 0 {
        return vec![cluster];
    }
    lens.sort_unstable();
    // Best two-way split by within-group sum of squares, via prefix sums.
    let mut s = vec![0f64; n + 1];
    let mut ss = vec![0f64; n + 1];
    for (k, &(l, _)) in lens.iter().enumerate() {
        s[k + 1] = s[k] + l as f64;
        ss[k + 1] = ss[k] + (l as f64) * (l as f64);
    }
    let sse = |a: usize, b: usize| {
        let m = (b - a) as f64;
        (ss[b] - ss[a]) - (s[b] - s[a]) * (s[b] - s[a]) / m
    };
    let k = (min_group..=n - min_group)
        .min_by(|&x, &y| (sse(0, x) + sse(x, n)).total_cmp(&(sse(0, y) + sse(y, n))))
        .expect("non-empty range");
    let mean = |a: usize, b: usize| (s[b] - s[a]) / (b - a) as f64;
    let sd = |a: usize, b: usize| (sse(a, b) / (b - a) as f64).max(0.0).sqrt();
    // Two alleles leave an empty stretch between the groups' lengths; a single noisy allele
    // split in two does not, even though the two halves' means are well apart.
    let boundary_gap = (lens[k].0 - lens[k - 1].0) as f64;
    let noise = sd(0, k).max(sd(k, n)).max(cfg.args.split_noise * mean(k, n)).max(1.0);
    if boundary_gap < cfg.args.split_min_gap_sd * noise {
        return vec![cluster];
    }
    let long = cluster.leads_long.clone();
    let shorter: HashSet<usize> = lens[..k].iter().map(|x| x.1).collect();
    let (a, b): (Vec<usize>, Vec<usize>) = cluster.leads.iter().partition(|i| shorter.contains(i));
    vec![cluster.with_leads(a, long.clone()), cluster.with_leads(b, long)]
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use clap::Parser;

    use super::*;
    use crate::config::Args;
    use crate::reads::{Lead, Source, SvType};

    fn table(lens: &[i64]) -> Leadtab {
        let leads = lens
            .iter()
            .enumerate()
            .map(|(k, &l)| {
                let lead = Lead {
                    read_id: k as u64,
                    qname: Arc::from(format!("r{k}").as_bytes()),
                    ref_start: 1000,
                    ref_end: 1000,
                    qry_start: 0,
                    qry_end: l,
                    reverse: k % 2 == 0,
                    mapq: 60,
                    nm: 0.0,
                    source: Source::Inline,
                    svlen: Some(l),
                    seq: None,
                    bnd: None,
                    hap: 0,
                    phase_set: None,
                    is_sa: false,
                    read_len: 10000,
                };
                (1000, lead)
            })
            .collect();
        Leadtab::new(leads, 10)
    }

    fn cfg() -> Config {
        let args = Args::parse_from(["lrpipe-sv", "-i", "x.bam", "--split-alleles"]);
        Config::from_args(args, String::new()).expect("valid config")
    }

    fn split_sizes(lens: &[i64]) -> Vec<usize> {
        let tab = table(lens);
        let c = Cluster::for_test(SvType::Ins, (0..lens.len()).collect());
        split(c, &tab, &cfg()).iter().map(|c| c.leads.len()).collect()
    }

    #[test]
    fn splits_two_clear_alleles() {
        let lens = [300, 302, 299, 301, 300, 330, 331, 329, 332, 330];
        assert_eq!(split_sizes(&lens), vec![5, 5]);
    }

    #[test]
    fn keeps_one_noisy_allele() {
        let lens = [295, 300, 305, 298, 302, 310, 290, 300, 304, 296];
        assert_eq!(split_sizes(&lens), vec![10]);
    }

    #[test]
    fn keeps_a_poorly_supported_minor_allele() {
        // Two reads are below both the 3-read and 20% minimums.
        let lens = [300, 301, 299, 300, 302, 300, 301, 299, 400, 401];
        assert_eq!(split_sizes(&lens), vec![10]);
    }
}

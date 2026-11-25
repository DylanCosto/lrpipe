//! Replays a problem written by tests/core/oracle.py and prints results in the
//! same format, for diffing against WhatsHap.
use std::collections::HashSet;
use std::io::{self, Read as _};

use lrpipe_phase::core::dp;
use lrpipe_phase::core::genotype::Genotype;
use lrpipe_phase::core::pedigree::Pedigree;
use lrpipe_phase::core::read::{Read, ReadSet};
use lrpipe_phase::pedigree::uniform_recombination_map;
use lrpipe_phase::readselect::readselection;

fn main() {
    let path = std::env::args().nth(1).unwrap();
    let text = std::fs::read_to_string(path).unwrap();
    let _ = io::stdin().read_to_string(&mut String::new());
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let n: usize = header[1].parse().unwrap();
    let mut rs = ReadSet::new();
    let mut haps_info: Vec<Vec<(i32, i32)>> = Vec::new();
    for _ in 0..n {
        let f: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
        let mut r = Read::new(f[0], 60, f[1].parse().unwrap(), f[2].parse().unwrap());
        for v in &f[3..] {
            let p: Vec<i32> = v.split(':').map(|x| x.parse().unwrap()).collect();
            r.add_variant(p[0], p[1], p[2]);
        }
        haps_info.push(vec![]);
        rs.add(r).unwrap();
    }
    let mut max_cov = 0;
    let mut pref: Option<HashSet<i32>> = None;
    let mut n_samples = 1;
    for line in lines {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "MAXCOV" => max_cov = f[1].parse().unwrap(),
            "PREFERRED" => pref = Some(f[1..].iter().map(|x| x.parse().unwrap()).collect()),
            "SAMPLES" => n_samples = f[1].parse().unwrap(),
            _ => {}
        }
    }
    // The oracle input is already sorted by WhatsHap; sort again to test our sort.
    let order_before: Vec<String> = rs.iter().map(|r| r.name.clone()).collect();
    rs.sort();
    let order_after: Vec<String> = rs.iter().map(|r| r.name.clone()).collect();
    assert_eq!(order_before, order_after, "sort order differs");
    println!("SORTED {}", order_after.join(" "));
    if rs.is_empty() {
        println!("EMPTY");
        return;
    }
    let mut selected_sets = Vec::new();
    for s in 0..n_samples {
        let sub = rs.subset((0..rs.len()).filter(|&i| rs.get(i).sample_id == s));
        if sub.is_empty() {
            println!("SELECTED {s}");
            selected_sets.push(sub);
            continue;
        }
        let sel = readselection(&sub, max_cov, pref.as_ref(), true).unwrap();
        let strs: Vec<String> = sel.iter().map(|x| x.to_string()).collect();
        println!("SELECTED {s} {}", strs.join(" "));
        selected_sets.push(sub.subset(sel));
    }
    let mut all_reads = ReadSet::new();
    for sub in selected_sets {
        for r in sub.into_reads() {
            all_reads.add(r).unwrap();
        }
    }
    all_reads.sort();
    let accessible = all_reads.get_positions();
    // Genotypes are not part of the input file: reconstruct them from the
    // GENOTYPES line the oracle appends.
    let gt_text = std::fs::read_to_string(std::env::args().nth(2).unwrap()).unwrap();
    let mut ped = Pedigree::new();
    for (s, line) in gt_text.lines().enumerate() {
        let gts: Vec<Genotype> = line
            .split_whitespace()
            .map(|g| {
                let a: Vec<u32> = g.split('/').map(|x| x.parse().unwrap()).collect();
                Genotype::new(&a).unwrap()
            })
            .collect();
        assert_eq!(gts.len(), accessible.len());
        let n = gts.len();
        ped.add_individual(s as u32, gts, vec![None; n]);
    }
    if n_samples == 3 {
        ped.add_relationship(0, 1, 2).unwrap();
    }
    let recomb = uniform_recombination_map(1.26, &accessible);
    all_reads.reassign_read_ids();
    let solution = match dp::solve(&all_reads, &recomb, &ped, false, Some(&accessible)) {
        Ok(s) => s,
        Err(e) => {
            println!("ERROR {e}");
            return;
        }
    };
    println!("COST {}", solution.optimal_score);
    for (s, sr) in solution.superreads.iter().enumerate() {
        for h in 0..2 {
            let v: Vec<String> =
                sr.get(h).variants.iter().map(|v| format!("{}:{}:{}", v.position, v.allele, v.quality)).collect();
            println!("SUPER {s} {h} {}", v.join(" "));
        }
    }
    let tvs: Vec<String> = solution.transmission_vector.iter().map(|x| x.to_string()).collect();
    println!("TV {}", tvs.join(" "));
    let p: Vec<String> = solution.partitioning.iter().map(|x| x.to_string()).collect();
    println!("PART {}", p.join(" "));
}

//! Replays the problems written by tests/core/genotype_oracle.py and prints
//! the results in the same format, for diffing against WhatsHap.
use std::io::{self, BufRead};

use lrpipe_phase::core::genotype::{Genotype, PhredGenotypeLikelihoods};
use lrpipe_phase::core::genotypedp::{GenotypeDPTable, compute_genotypes};
use lrpipe_phase::core::pedigree::Pedigree;
use lrpipe_phase::core::read::{Read, ReadSet};

/// Python's float.hex()
fn hex(x: f64) -> String {
    if x == 0.0 {
        return if x.is_sign_negative() { "-0x0.0p+0".into() } else { "0x0.0p+0".into() };
    }
    let bits = x.to_bits();
    let sign = if bits >> 63 == 1 { "-" } else { "" };
    let e = ((bits >> 52) & 0x7ff) as i64;
    let m = bits & ((1 << 52) - 1);
    let (lead, exp) = if e == 0 { (0, -1022) } else { (1, e - 1023) };
    let frac = format!("{m:013x}");
    format!("{sign}0x{lead}.{frac}p{}{}", if exp >= 0 { "+" } else { "-" }, exp.abs())
}

fn parse_hex(s: &str) -> f64 {
    // 0x1.<hex>p<exp> as written by float.hex()
    let (neg, s) = s.strip_prefix('-').map_or((false, s), |r| (true, r));
    let s = s.strip_prefix("0x").unwrap();
    let (mant, exp) = s.split_once('p').unwrap();
    let (int, frac) = mant.split_once('.').unwrap();
    let mut v = int.parse::<u64>().unwrap() as f64;
    let mut scale = 1.0 / 16.0;
    for c in frac.chars() {
        v += c.to_digit(16).unwrap() as f64 * scale;
        scale /= 16.0;
    }
    let v = v * 2f64.powi(exp.parse::<i32>().unwrap());
    if neg { -v } else { v }
}

fn main() {
    let lines: Vec<String> = io::stdin().lock().lines().map(|l| l.unwrap()).collect();
    let mut i = 0;
    while i < lines.len() {
        let header: Vec<&str> = lines[i].split(' ').collect();
        assert_eq!(header[0], "P");
        let n_ind: usize = header[1].parse().unwrap();
        let positions: Vec<i32> = header[2..].iter().map(|p| p.parse().unwrap()).collect();
        println!("{}", lines[i]);
        i += 1;
        let mut rs = ReadSet::new();
        let mut pedigree = Pedigree::new();
        let mut recomb: Vec<u32> = Vec::new();
        while !lines[i].starts_with('L') && !lines[i].starts_with('G') && !lines[i].starts_with('E') {
            let f: Vec<&str> = lines[i].split(' ').collect();
            match f[0] {
                "R" => {
                    let mut read = Read::new(f[1], f[2].parse().unwrap(), 0, f[3].parse().unwrap());
                    for v in &f[4..] {
                        let p: Vec<i32> = v.split(':').map(|x| x.parse().unwrap()).collect();
                        read.add_variant(p[0], p[1], p[2]);
                    }
                    rs.add(read).unwrap();
                }
                "I" => {
                    let gls = f[2..]
                        .iter()
                        .map(|t| {
                            let g: Vec<f64> = t.split(',').map(parse_hex).collect();
                            Some(PhredGenotypeLikelihoods::new(g, 2, 2).unwrap())
                        })
                        .collect();
                    pedigree.add_individual(f[1].parse().unwrap(), vec![Genotype::none(); positions.len()], gls);
                }
                "T" => pedigree
                    .add_relationship(f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap())
                    .unwrap(),
                "C" => recomb = f[1..].iter().map(|c| c.parse().unwrap()).collect(),
                _ => panic!("unexpected line"),
            }
            println!("{}", lines[i]);
            i += 1;
        }
        rs.sort();
        let mut accessible = rs.get_positions();
        accessible.sort_unstable();
        if !accessible.is_empty() {
            let costs: Vec<u32> =
                accessible.iter().map(|p| recomb[positions.iter().position(|q| q == p).unwrap()]).collect();
            let table = GenotypeDPTable::new(&mut rs, &costs, &pedigree, &accessible).unwrap();
            for ind in 0..n_ind {
                for k in 0..accessible.len() {
                    let l = table.genotype_likelihoods(ind as u32, k);
                    println!("L {ind} {k} {}", l.iter().map(|x| hex(*x)).collect::<Vec<_>>().join(" "));
                }
            }
        }
        let (genotypes, gls) = compute_genotypes(&rs, &positions).unwrap();
        for (g, l) in genotypes.iter().zip(&gls) {
            let v = g.as_vector();
            let gs = if v.is_empty() {
                "none".to_string()
            } else {
                v.iter().map(|a| a.to_string()).collect::<Vec<_>>().join("/")
            };
            println!("G {gs} {}", l.iter().map(|x| hex(*x)).collect::<Vec<_>>().join(" "));
        }
        println!("E");
        while lines[i] != "E" {
            i += 1;
        }
        i += 1;
    }
}

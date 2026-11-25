//! Replays the problems written by tests/core/switchflip_oracle.py and prints
//! the results in the same format, for diffing against WhatsHap.
use std::io::{self, BufRead};

use lrpipe_phase::compat::pyfloat;
use lrpipe_phase::core::switchflip::SwitchFlipCalculator;

fn main() {
    let stdin = io::stdin();
    let lines: Vec<String> = stdin.lock().lines().map(|l| l.unwrap()).collect();
    for problem in lines.iter().step_by(2) {
        let f: Vec<&str> = problem.split(' ').collect();
        let ploidy: u32 = f[0].parse().unwrap();
        let (switch_cost, flip_cost): (f64, f64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        let haps = |s: &str| -> Vec<Vec<u32>> {
            s.split(',').map(|h| h.bytes().map(|c| u32::from(c - b'0')).collect()).collect()
        };
        let (h0, h1) = (haps(f[3]), haps(f[4]));
        // position-wise, as the Cython wrapper converts them
        let n = h0[0].len();
        let pos = |h: &[Vec<u32>]| -> Vec<Vec<u32>> { (0..n).map(|i| h.iter().map(|hap| hap[i]).collect()).collect() };
        let r = SwitchFlipCalculator::new(ploidy, switch_cost, flip_cost).compare(&pos(&h0), &pos(&h1));
        println!("{problem}");
        let join = |v: &[Vec<u32>]| {
            v.iter()
                .map(|x| x.iter().map(|a| a.to_string()).collect::<Vec<_>>().join(" "))
                .collect::<Vec<_>>()
                .join(";")
        };
        println!(
            "{} {} {} {} {}",
            pyfloat::repr(r.switches),
            pyfloat::repr(r.flips),
            r.switches_in_column.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(","),
            join(&r.flipped_haps_in_column),
            join(&r.perm_in_column)
        );
    }
}

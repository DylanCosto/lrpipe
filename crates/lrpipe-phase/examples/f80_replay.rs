//! Checks compat::f80 against cases from tests/compat/f80_oracle.c.
use std::io::{self, BufRead};

use lrpipe_phase::compat::f80::{Divisor, F80};

fn parse(h: &str) -> F80 {
    F80::from_bits(u16::from_str_radix(&h[..4], 16).unwrap(), u64::from_str_radix(&h[4..], 16).unwrap())
}

fn main() {
    let (mut n, mut bad) = (0, 0);
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let f: Vec<F80> = line.split_whitespace().map(parse).collect();
        let (a, b) = (f[0], f[1]);
        let got = [a + b, a - b, a * b, a / b, F80::from_f64(a.to_f64())];
        if a.div_by(&Divisor::new(b)) != f[5] {
            bad += 1;
            eprintln!("div_by differs: {line}");
        }
        for (k, g) in got.iter().enumerate() {
            if *g != f[2 + k] {
                bad += 1;
                if bad <= 10 {
                    eprintln!("op {k} differs: {line}");
                }
            }
        }
        n += 1;
    }
    println!("{n} cases, {bad} differences");
}

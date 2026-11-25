//! Replays a script of set operations (see tests/compat/*.py) against the
//! Rust emulations and prints the results, for diffing against the real
//! CPython / libstdc++ implementations.
use lrpipe_phase::compat::cxx::{CxxUnorderedSetI32, std_hash_string};
use lrpipe_phase::compat::pyset::PySet;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};

fn main() {
    let stdin = io::stdin();
    let mut out = io::BufWriter::new(io::stdout());
    let mut sets: HashMap<String, PySet> = HashMap::new();
    let mut cxx = CxxUnorderedSetI32::new();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.is_empty() {
            continue;
        }
        let num = |i: usize| f[i].parse::<usize>().unwrap();
        match f[0] {
            "new" => {
                sets.insert(f[1].into(), PySet::new());
            }
            "range" => {
                sets.insert(f[1].into(), PySet::from_iter_adds(0..num(2)));
            }
            "add" => sets.get_mut(f[1]).unwrap().add(num(2)),
            "discard" => {
                sets.get_mut(f[1]).unwrap().discard(num(2));
            }
            "isub" => {
                let b = sets[f[2]].clone();
                sets.get_mut(f[1]).unwrap().difference_update(&b);
            }
            "copy" => {
                let c = PySet::copy_of(&sets[f[2]]);
                sets.insert(f[1].into(), c);
            }
            "diff" => {
                let d = sets[f[2]].difference(&sets[f[3]]);
                sets.insert(f[1].into(), d);
            }
            "updlist" => {
                let v: Vec<usize> = (2..f.len()).map(num).collect();
                sets.get_mut(f[1]).unwrap().update_iter(v);
            }
            "merge" => {
                let b = sets[f[2]].clone();
                sets.get_mut(f[1]).unwrap().merge(&b);
            }
            "print" => {
                let v: Vec<String> = sets[f[1]].iter().map(|x| x.to_string()).collect();
                writeln!(out, "{}", v.join(" ")).unwrap();
            }
            // libstdc++ commands (same syntax as cxx_oracle.cpp)
            "H" => writeln!(out, "{}", std_hash_string(f[1].as_bytes())).unwrap(),
            "N" => cxx = CxxUnorderedSetI32::new(),
            "I" => {
                cxx.insert(f[1].parse::<i32>().unwrap());
            }
            "C" => cxx.clear(),
            "P" => {
                let v: Vec<String> = cxx.iter().map(|x| format!("{x} ")).collect();
                writeln!(out, "{}| {}", v.join(""), cxx.bucket_count()).unwrap();
            }
            other => panic!("unknown command {other}"),
        }
    }
}

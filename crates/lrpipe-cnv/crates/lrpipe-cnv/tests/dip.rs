//! The dip test and UniDip against the results of Spectre's (GPL) implementation, run as a
//! black box by `tools/dip_vectors.py`.
//!
//! Exact: every p-value (same numpy draws), and UniDip on every input whose first test finds
//! one mode, which is what heterozygous SNV depths of a germline sample give. Within a few ulps:
//! the dip (arithmetic order differs). Not reproduced: the modal interval and so UniDip's
//! recursion on multimodal inputs, whose bookkeeping is specific to that implementation; the
//! agreement there is printed, not asserted.

use std::io::Read;

use lrpipe_cnv::dip::{dip, dip_test, unidip};
use pynum::numpy_random::RandomState;
use serde_json::Value;

fn f(v: &Value) -> Option<f64> {
    let hex = v.as_str()?;
    let bytes: Vec<u8> = (0..16).step_by(2).map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap()).collect();
    Some(f64::from_le_bytes(bytes.try_into().unwrap()))
}

fn vectors() -> Vec<Value> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/dip.json.gz");
    let mut s = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap()).read_to_string(&mut s).unwrap();
    serde_json::from_str(&s).unwrap()
}

#[test]
fn matches_spectres_dip_and_unidip() {
    let (mut dips, mut modal, mut pvals, mut runs, mut total) = (0, 0, 0, 0, 0);
    let (mut unimodal, mut unimodal_runs) = (0, 0);
    let mut far_dips = Vec::new();
    let mut notes = Vec::new();
    for (i, rec) in vectors().iter().enumerate() {
        let data: Vec<f64> = rec["data"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        total += 1;
        let d = dip(&data);
        let want = f(&rec["dip"]["d"]).unwrap();
        let (wl, wr, wm) = (
            rec["dip"]["left"].as_i64().unwrap(),
            rec["dip"]["right"].as_i64().unwrap(),
            rec["dip"]["idxs"].as_i64().unwrap(),
        );
        let want_modal = ((wl - 1) as isize, (wm - wr) as isize);
        if d.dip.to_bits() == want.to_bits() {
            dips += 1;
        }
        if (d.dip - want).abs() > 1e-15 {
            far_dips.push(format!("#{i}: dip {} want {want}", d.dip));
        }
        if d.modal == want_modal {
            modal += 1;
        }
        let t = &rec["diptst"][0];
        let mut rng = RandomState::new(t["seed"].as_u64().unwrap() as u32);
        let (_, p) = dip_test(&data, 1000, &mut rng);
        if p.map(f64::to_bits) == f(&t["p"]).map(f64::to_bits) {
            pvals += 1;
        } else if notes.len() < 25 {
            notes.push(format!("#{i}: p {p:?} want {:?}", f(&t["p"])));
        }
        if std::env::var_os("DIP_ONLY").is_some() {
            continue;
        }
        let u = &rec["unidip"][0];
        let mut rng = RandomState::new(u["seed"].as_u64().unwrap() as u32);
        let got: Vec<(usize, usize)> = unidip(&data, &mut rng);
        let first = dip_test(&data, 100, &mut RandomState::new(u["seed"].as_u64().unwrap() as u32)).1;
        let want: Vec<(usize, usize)> = u["intervals"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| (p[0].as_u64().unwrap() as usize, p[1].as_u64().unwrap() as usize))
            .collect();
        if first.is_some_and(|p| p > 0.05) {
            unimodal += 1;
            unimodal_runs += usize::from(got == want);
        }
        if got == want {
            runs += 1;
        } else {
            notes.push(format!("#{i} n={} distinct={}: unidip {got:?} want {want:?}", data.len(), d.distinct));
        }
    }
    eprintln!(
        "dip bit-exact {dips}/{total}, modal interval {modal}/{total}, p-value {pvals}/{total}, \
         UniDip {runs}/{total} (unimodal inputs {unimodal_runs}/{unimodal})\n{}",
        notes.join("\n")
    );
    assert_eq!(pvals, total, "p-values differ:\n{}", notes.join("\n"));
    assert!(far_dips.is_empty(), "dips differ by more than 1e-15:\n{}", far_dips.join("\n"));
    assert_eq!(unimodal_runs, unimodal, "UniDip differs on unimodal inputs:\n{}", notes.join("\n"));
}

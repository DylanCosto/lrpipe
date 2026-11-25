//! Every function against results from the reference environment
//! (`tools/pynum_vectors.py`), compared bit for bit.

use std::io::Read;
use std::path::Path;
use std::sync::OnceLock;

use md5::{Digest, Md5};
use pynum::{numpy_random, random, repr, round, special, sum};
use serde_json::Value;

fn vectors() -> &'static Value {
    static V: OnceLock<Value> = OnceLock::new();
    V.get_or_init(|| {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vectors.json.gz");
        let mut text = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(path).expect("opening the vectors"))
            .read_to_string(&mut text)
            .expect("reading the vectors");
        serde_json::from_str(&text).expect("parsing the vectors")
    })
}

fn f(v: &Value) -> f64 {
    let hex = v.as_str().expect("a float's bits as hex");
    let bytes: Vec<u8> = (0..16).step_by(2).map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap()).collect();
    f64::from_le_bytes(bytes.try_into().unwrap())
}

fn array(i: &Value) -> Vec<f64> {
    vectors()["arrays"][i.as_u64().unwrap() as usize].as_array().unwrap().iter().map(f).collect()
}

/// Equal bits, treating every NaN as equal.
fn same(got: f64, want: f64) -> bool {
    got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan())
}

fn check(failures: &mut Vec<String>, what: &str, got: f64, want: f64) {
    if !same(got, want) {
        failures.push(format!("{what}: got {got:e} ({:#x}), want {want:e} ({:#x})", got.to_bits(), want.to_bits()));
    }
}

/// As [`check`], within `ulps` units in the last place.
fn check_ulps(failures: &mut Vec<String>, what: &str, got: f64, want: f64, ulps: u64) {
    if !same(got, want)
        && (got.is_sign_negative() != want.is_sign_negative() || got.to_bits().abs_diff(want.to_bits()) > ulps)
    {
        failures.push(format!("{what}: got {got:e} ({:#x}), want {want:e} ({:#x})", got.to_bits(), want.to_bits()));
    }
}

fn report(failures: Vec<String>) {
    let n = failures.len();
    assert!(failures.is_empty(), "{n} mismatches, first: {:#?}", &failures[..n.min(10)]);
}

#[test]
fn sums_means_and_variances() {
    let mut failures = Vec::new();
    for c in vectors()["cases"]["sum"].as_array().unwrap() {
        let a = array(&c["a"]);
        let n = a.len();
        check(&mut failures, &format!("sum n={n}"), sum::sum(&a), f(&c["sum"]));
        if n > 0 {
            check(&mut failures, &format!("mean n={n}"), sum::mean(&a), f(&c["mean"]));
            check(&mut failures, &format!("var n={n}"), sum::var(&a), f(&c["var"]));
            check(&mut failures, &format!("std n={n}"), sum::std(&a), f(&c["std"]));
        }
    }
    report(failures);
}

#[test]
fn nan_aware_reductions() {
    let mut failures = Vec::new();
    for c in vectors()["cases"]["nan"].as_array().unwrap() {
        let a = array(&c["a"]);
        let n = a.len();
        check(&mut failures, &format!("nanmean n={n}"), sum::nanmean(&a), f(&c["nanmean"]));
        check(&mut failures, &format!("nanstd n={n}"), sum::nanstd(&a), f(&c["nanstd"]));
        check(&mut failures, &format!("nanmedian n={n}"), sum::nanmedian(&a), f(&c["nanmedian"]));
        check(&mut failures, &format!("Series.mean n={n}"), sum::series_mean(&a), f(&c["series_mean"]));
        if !c["series_std"].is_null() {
            check(&mut failures, &format!("Series.std n={n}"), sum::series_std(&a), f(&c["series_std"]));
        }
    }
    report(failures);
}

#[test]
fn percentiles() {
    let mut failures = Vec::new();
    for c in vectors()["cases"]["percentile"].as_array().unwrap() {
        let a = array(&c["a"]);
        let q = f(&c["q"]);
        check(&mut failures, &format!("percentile n={} q={q}", a.len()), sum::percentile(&a, q), f(&c["r"]));
    }
    report(failures);
}

#[test]
fn numpy_rounding() {
    let mut failures = Vec::new();
    for c in vectors()["round"].as_array().unwrap() {
        let (x, d) = (f(&c["x"]), c["d"].as_u64().unwrap() as u32);
        check(&mut failures, &format!("round({x:e}, {d})"), round::np_round(x, d), f(&c["r"]));
    }
    report(failures);
}

#[test]
fn float_text() {
    let mut failures = Vec::new();
    for c in vectors()["repr"].as_array().unwrap() {
        let x = f(&c["x"]);
        let got = repr::py_repr(x);
        for key in ["repr", "np_str"] {
            if got != c[key].as_str().unwrap() {
                failures.push(format!("{key}({x:e}): got {got}, want {}", c[key]));
            }
        }
    }
    report(failures);
}

#[test]
fn cpython_random() {
    let chars: Vec<char> = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789".chars().collect();
    for c in vectors()["py_random"].as_array().unwrap() {
        let seed = c["seed"].as_str().unwrap();
        let sample_id_digest = seed;
        let mut r = random::Random::from_str_seed(sample_id_digest);
        let bits: Vec<u64> = (0..5).map(|_| u64::from(r.getrandbits(32))).collect();
        let want: Vec<u64> = c["getrandbits32"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap()).collect();
        assert_eq!(bits, want, "getrandbits(32), seed {seed}");
        for (key, value) in c.as_object().unwrap() {
            let Some(nk) = key.strip_prefix("sample_") else { continue };
            let (n, k) = nk.split_once('_').unwrap();
            let mut r = random::Random::from_str_seed(seed);
            let got = r.sample_range(n.parse().unwrap(), k.parse().unwrap());
            let want: Vec<u64> = value.as_array().unwrap().iter().map(|v| v.as_u64().unwrap()).collect();
            assert!(got == want, "sample(range({n}), {k}), seed {seed}");
        }
        let mut r = random::Random::from_str_seed(seed);
        let got: String = (0..64).map(|_| *r.choice(&chars)).collect();
        assert_eq!(got, c["choice"].as_str().unwrap(), "choice, seed {seed}");
    }
}

#[test]
fn seeds_are_md5_hex_of_the_sample_id() {
    let hex = |s: &str| Md5::digest(s.as_bytes()).iter().map(|b| format!("{b:02x}")).collect::<String>();
    let seeds: Vec<&str> =
        vectors()["py_random"].as_array().unwrap().iter().map(|c| c["seed"].as_str().unwrap()).collect();
    assert_eq!(seeds, [hex("HG003"), hex("HG002"), hex("sample-1"), hex("")]);
}

#[test]
fn numpy_legacy_uniform() {
    let mut failures = Vec::new();
    for c in vectors()["np_random"].as_array().unwrap() {
        let seed = c["seed"].as_u64().unwrap() as u32;
        let mut rs = numpy_random::RandomState::new(seed);
        for (i, want) in c["uniform"].as_array().unwrap().iter().enumerate() {
            check(&mut failures, &format!("uniform seed {seed} #{i}"), rs.uniform(), f(want));
        }
    }
    report(failures);
}

/// `norm.sf` goes through the C library's `exp`, as SciPy's does, and `exp` is not correctly
/// rounded: glibc 2.28 (the release build's) and macOS round a few values a unit or two in the
/// last place away from the glibc the vectors were made with. Bit for bit on the same `exp`.
#[test]
fn normal_tail() {
    let mut failures = Vec::new();
    for c in vectors()["norm_sf"].as_array().unwrap() {
        let x = f(&c["x"]);
        check_ulps(&mut failures, &format!("norm.sf({x:e})"), special::norm_sf(x), f(&c["r"]), 2);
    }
    report(failures);
}

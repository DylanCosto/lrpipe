//! Reads every checkpoint the original ships and compares it with PyTorch's own reading.
//!
//! Opt-in: set `SNV_MODELS` to the original's `models/` directory and `SNV_DIGESTS` to
//! the output of `tools/checkpoint_digest.py` for that directory. Each line is
//! `<model>/<file> <tensor count> <sha256 prefix>` over names, shape reprs and raw bytes.

use std::path::Path;

use sha2::{Digest, Sha256};
use snv_infer::{TensorValues, read_state_dict};

/// Python's `str(tuple(shape))`.
fn tuple_repr(shape: &[usize]) -> String {
    match shape {
        [d] => format!("({d},)"),
        _ => format!("({})", shape.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")),
    }
}

#[test]
fn all_checkpoints_match_pytorch() {
    let (Ok(models), Ok(digests)) = (std::env::var("SNV_MODELS"), std::env::var("SNV_DIGESTS")) else {
        eprintln!("skipped: set SNV_MODELS and SNV_DIGESTS");
        return;
    };
    let digests = std::fs::read_to_string(&digests).expect("SNV_DIGESTS is readable");
    let mut checked = 0;
    for line in digests.lines() {
        let [file, count, want] = line.split(' ').collect::<Vec<_>>()[..] else { panic!("bad digest line {line:?}") };
        let tensors = read_state_dict(&Path::new(&models).join(file)).unwrap_or_else(|e| panic!("{file}: {e}"));
        let mut h = Sha256::new();
        for (name, t) in &tensors {
            h.update(name.as_bytes());
            h.update(tuple_repr(&t.shape).as_bytes());
            match &t.values {
                TensorValues::F32(v) => v.iter().for_each(|x| h.update(x.to_le_bytes())),
                TensorValues::I64(v) => v.iter().for_each(|x| h.update(x.to_le_bytes())),
            }
        }
        let got: String = h.finalize().iter().take(8).map(|b| format!("{b:02x}")).collect();
        assert_eq!((tensors.len().to_string(), got.as_str()), (count.to_owned(), want), "{file}");
        checked += 1;
    }
    assert!(checked > 0, "no digests in SNV_DIGESTS");
}

//! Python numbers that may be an `int` or a `float`.

/// A value Spectre leaves as a Python `int` or `float` depending on the path taken, which
/// shows when it is written out (`1` versus `1.0`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Number {
    Int(i64),
    Float(f64),
}

impl Number {
    pub fn value(self) -> f64 {
        match self {
            Number::Int(i) => i as f64,
            Number::Float(f) => f,
        }
    }
}

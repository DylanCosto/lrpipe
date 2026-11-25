//! Python's `repr()`/`str()` of a float, and numbers that are either a
//! Python `int` or `float`.

/// A Python number: `int` or `float` (they print differently).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PyNum {
    Int(i64),
    Float(f64),
}

impl PyNum {
    pub fn as_f64(self) -> f64 {
        match self {
            PyNum::Int(v) => v as f64,
            PyNum::Float(v) => v,
        }
    }

    pub fn is_nan(self) -> bool {
        matches!(self, PyNum::Float(v) if v.is_nan())
    }

    /// `str(value)`
    pub fn py_str(self) -> String {
        match self {
            PyNum::Int(v) => v.to_string(),
            PyNum::Float(v) => repr(v),
        }
    }
}

impl std::ops::Add for PyNum {
    type Output = PyNum;

    /// `a + b`: an int only if both are ints
    fn add(self, other: PyNum) -> PyNum {
        match (self, other) {
            (PyNum::Int(a), PyNum::Int(b)) => PyNum::Int(a + b),
            (a, b) => PyNum::Float(a.as_f64() + b.as_f64()),
        }
    }
}

impl std::ops::AddAssign for PyNum {
    fn add_assign(&mut self, other: PyNum) {
        *self = *self + other;
    }
}

impl std::fmt::Display for PyNum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.py_str())
    }
}

/// `repr(x)` for a Python float: the shortest digits that round-trip, in
/// fixed notation for decimal exponents -4..16 and scientific otherwise
/// (`1e+16`, `1.5e-05`).
pub fn repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    // Rust's `{:e}` also gives the shortest round-trip digits: "-1.2345e-5"
    let s = format!("{x:e}");
    let (mantissa, exp) = s.split_once('e').expect("LowerExp output has an exponent");
    let exp: i32 = exp.parse().expect("LowerExp exponent is an integer");
    let (sign, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => ("-", m),
        None => ("", mantissa),
    };
    let digits: String = mantissa.chars().filter(|&c| c != '.').collect();
    let n = digits.len() as i32;
    if (-4..16).contains(&exp) {
        if exp < 0 {
            format!("{sign}0.{}{digits}", "0".repeat((-exp - 1) as usize))
        } else if exp + 1 >= n {
            format!("{sign}{digits}{}.0", "0".repeat((exp + 1 - n) as usize))
        } else {
            let (int, frac) = digits.split_at((exp + 1) as usize);
            format!("{sign}{int}.{frac}")
        }
    } else {
        let exp_sign = if exp < 0 { '-' } else { '+' };
        let (first, rest) = digits.split_at(1);
        let rest = if rest.is_empty() { String::new() } else { format!(".{rest}") };
        format!("{sign}{first}{rest}e{exp_sign}{:02}", exp.abs())
    }
}

#[cfg(test)]
mod tests {
    use super::repr;

    #[test]
    fn matches_python() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (5.0, "5.0"),
            (0.5, "0.5"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (1e15, "1000000000000000.0"),
            (123456.789, "123456.789"),
            (0.0001, "0.0001"),
            (1e-5, "1e-05"),
            (-2.5e-7, "-2.5e-07"),
            (1.0 / 3.0, "0.3333333333333333"),
            (1e300, "1e+300"),
            (f64::NAN, "nan"),
            (f64::NEG_INFINITY, "-inf"),
        ];
        for (x, want) in cases {
            assert_eq!(repr(x), want, "repr({x})");
        }
    }
}

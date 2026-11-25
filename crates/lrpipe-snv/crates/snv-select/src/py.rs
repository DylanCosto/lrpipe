//! The few Python built-ins whose exact behaviour reaches the original's output bytes:
//! `str.split()`, `float()`, `int()`, `str(float)`, `round()` and slicing with
//! `int(pct * len(xs))`.

/// `str.isspace()` for one character, the set `str.split()` and `str.rstrip()` use: Unicode
/// `White_Space` plus the ASCII separators `\x1c`..`\x1f`, which Python also treats as space.
pub(crate) fn is_space(c: char) -> bool {
    c.is_whitespace() || ('\x1c'..='\x1f').contains(&c)
}

/// `s.split()`: fields separated by runs of whitespace, ignoring leading and trailing space.
pub(crate) fn split(s: &str) -> Vec<&str> {
    s.split(is_space).filter(|f| !f.is_empty()).collect()
}

/// `s.rstrip()`.
pub(crate) fn rstrip(s: &str) -> &str {
    s.trim_end_matches(is_space)
}

/// `s.strip()`.
pub(crate) fn strip(s: &str) -> &str {
    s.trim_matches(is_space)
}

/// Removes the underscores Python accepts in numeric literals passed to `float()` and `int()`:
/// each must sit between two digits. `None` if one does not.
fn without_underscores(s: &str) -> Option<String> {
    let bytes = s.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'_' {
            let before = i.checked_sub(1).and_then(|j| bytes.get(j));
            let after = bytes.get(i + 1);
            if !before.is_some_and(u8::is_ascii_digit) || !after.is_some_and(u8::is_ascii_digit) {
                return None;
            }
        }
    }
    Some(s.replace('_', ""))
}

/// `float(s)` for a string: surrounding whitespace, a sign, `inf`/`infinity`/`nan` in any case,
/// and underscores between digits are accepted. (Non-ASCII digits, which Python also accepts,
/// are not.)
pub(crate) fn parse_float(s: &str) -> Option<f64> {
    let s = strip(s);
    let s = if s.contains('_') { without_underscores(s)? } else { s.to_owned() };
    // Rust's `f64` grammar is Python's for ASCII input: optional sign, `inf`/`infinity`/`nan`
    // in any case, digits with an optional `.` and exponent.
    s.parse().ok()
}

/// `int(s)` for a string in base 10: surrounding whitespace, a sign, leading zeros and
/// underscores between digits. `None` if invalid or outside `i64`.
pub(crate) fn parse_int(s: &str) -> Option<i64> {
    let s = strip(s);
    let s = if s.contains('_') { without_underscores(s)? } else { s.to_owned() };
    let digits = s.strip_prefix(['+', '-']).unwrap_or(&s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse().ok()
}

/// `str(x)` (`repr(x)`) for a float: the shortest digits that round-trip, in fixed notation with
/// at least one fractional digit (`16.39`, `0.0`, `0.0001`), or in exponent notation when the
/// decimal exponent is below -4 or at least 16 (`1e-05`, `1.5e+16`).
pub(crate) fn float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_owned();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_owned();
    }
    let sci = shortest_closest_sci(x.abs());
    let (mantissa, exp) = sci.split_once('e').expect("`{:e}` output always contains an exponent");
    let exp: i32 = exp.parse().expect("`{:e}` exponent is an integer");
    let digits: String = mantissa.chars().filter(|&c| c != '.').collect();
    let sign = if x.is_sign_negative() { "-" } else { "" };
    // Python's `decpt`: the value is 0.<digits> * 10^decpt.
    let decpt = exp + 1;
    let n = i32::try_from(digits.len()).expect("a float has at most 17 significant digits");
    let body = if decpt <= -4 || decpt > 16 {
        let (first, rest) = digits.split_at(1);
        let dot = if rest.is_empty() { String::new() } else { format!(".{rest}") };
        let exp_sign = if exp < 0 { '-' } else { '+' };
        format!("{first}{dot}e{exp_sign}{:02}", exp.unsigned_abs())
    } else if decpt <= 0 {
        format!("0.{}{digits}", "0".repeat(decpt.unsigned_abs() as usize))
    } else if decpt >= n {
        format!("{digits}{}.0", "0".repeat((decpt - n) as usize))
    } else {
        let (int, frac) = digits.split_at(decpt as usize);
        format!("{int}.{frac}")
    };
    format!("{sign}{body}")
}

/// The shortest digits that round-trip to `x` (finite, non-negative) as `d.ddde<exp>`, choosing
/// the one closest to `x` when several have that length, as Python's `repr` does.
fn shortest_closest_sci(x: f64) -> String {
    // Rust's `{:e}` finds the shortest length but not always the closest digits of that length
    // (1005369574750092.2 comes out as ...092.3); the correctly rounded string of the same length
    // is the closest one, so use it whenever it round-trips.
    let shortest = format!("{x:e}");
    let mantissa = shortest.split_once('e').map_or(shortest.as_str(), |(m, _)| m);
    let len = mantissa.bytes().filter(u8::is_ascii_digit).count();
    let closest = format!("{x:.prec$e}", prec = len - 1);
    if closest.parse::<f64>() == Ok(x) { closest } else { shortest }
}

/// `str(round(x))` with no digits: rounds half to even and prints the exact integer. `None` for
/// infinity or NaN, where Python raises.
pub(crate) fn round_to_int_repr(x: f64) -> Option<String> {
    if !x.is_finite() {
        return None;
    }
    let rounded = x.round_ties_even();
    // `{:.0}` prints the exact integer value of the float; `int()` has no negative zero.
    Some(if rounded == 0.0 { "0".to_owned() } else { format!("{rounded:.0}") })
}

/// `str(round(x, 0))`: a float rounded half to even, keeping the sign of zero.
pub(crate) fn round_to_float_repr(x: f64) -> String {
    float_repr(if x.is_finite() { x.round_ties_even() } else { x })
}

/// `len(xs[:int(pct * len(xs))])`: how many of the lowest values the original keeps. `int()`
/// truncates toward zero, and a negative end counts from the back, so a negative proportion drops
/// the highest values instead. `None` when the product is infinite or NaN, where `int()` raises.
pub(crate) fn prefix_len(pct: f64, len: usize) -> Option<usize> {
    // Python converts the int to a float for the product; lengths here are far below 2^53.
    #[expect(clippy::cast_precision_loss, reason = "Python's int * float does the same conversion")]
    let product = pct * len as f64;
    if !product.is_finite() {
        return None;
    }
    let end = product.trunc();
    #[expect(clippy::cast_precision_loss, reason = "compared against an exact float of the length")]
    let len_f = len as f64;
    Some(if end >= 0.0 {
        if end >= len_f { len } else { end as usize }
    } else if -end >= len_f {
        0
    } else {
        len - (-end) as usize
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_repr_matches_python() {
        // Expected strings from CPython/PyPy `repr`.
        let cases: &[(f64, &str)] = &[
            (16.39, "16.39"),
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (23.04, "23.04"),
            (1.0 - 0.7, "0.30000000000000004"),
            // Rust's shortest `{:e}` gives ...092.3 here; Python picks the closer ...092.2.
            (f64::from_bits(0x430c_9306_c2f7_bc62), "1005369574750092.2"),
            (0.1, "0.1"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (9999999999999998.0, "9999999999999998.0"),
            (123456789012345680.0, "1.2345678901234568e+17"),
            (1e-5, "1e-05"),
            (1.5e-5, "1.5e-05"),
            (0.0001, "0.0001"),
            (0.00012, "0.00012"),
            (1e100, "1e+100"),
            (5e-324, "5e-324"),
            (f64::MAX, "1.7976931348623157e+308"),
            (-2.5, "-2.5"),
            (100.0, "100.0"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (f64::NAN, "nan"),
        ];
        for &(x, want) in cases {
            assert_eq!(float_repr(x), want, "repr({x:e})");
        }
    }

    #[test]
    fn round_matches_python() {
        assert_eq!(round_to_int_repr(16.39).as_deref(), Some("16"));
        assert_eq!(round_to_int_repr(2.5).as_deref(), Some("2"));
        assert_eq!(round_to_int_repr(3.5).as_deref(), Some("4"));
        assert_eq!(round_to_int_repr(-0.4).as_deref(), Some("0"));
        assert_eq!(round_to_int_repr(-1.5).as_deref(), Some("-2"));
        assert_eq!(
            round_to_int_repr(1e300).as_deref(),
            Some(
                "1000000000000000052504760255204420248704468581108159154915854115511802457988908195786371375080447864\
                 0437044438328838781769425232353604305756447921847867069828483872009265758037378302337947880900593689\
                 53234970799945081119038967640880074652742780142494579258788820056842838115669472196386865459400540160"
            )
        );
        assert_eq!(round_to_int_repr(f64::INFINITY), None);
        assert_eq!(round_to_float_repr(16.39), "16.0");
        assert_eq!(round_to_float_repr(22.5), "22.0");
        assert_eq!(round_to_float_repr(-0.4), "-0.0");
        assert_eq!(round_to_float_repr(f64::INFINITY), "inf");
    }

    #[test]
    fn parse_float_matches_python() {
        assert_eq!(parse_float("16.39"), Some(16.39));
        assert_eq!(parse_float(" 16.39\n"), Some(16.39));
        assert_eq!(parse_float("1_0.5"), Some(10.5));
        assert_eq!(parse_float("1e5"), Some(1e5));
        assert_eq!(parse_float("1."), Some(1.0));
        assert_eq!(parse_float(".5"), Some(0.5));
        assert_eq!(parse_float("-Infinity"), Some(f64::NEG_INFINITY));
        assert_eq!(parse_float("iNf"), Some(f64::INFINITY));
        assert!(parse_float("NaN").is_some_and(f64::is_nan));
        for bad in ["", ".", "_1", "1_", "1__0", "1e", "e5", "16.39x", "."] {
            assert_eq!(parse_float(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn parse_int_matches_python() {
        assert_eq!(parse_int("60230"), Some(60230));
        assert_eq!(parse_int("+007"), Some(7));
        assert_eq!(parse_int("1_000"), Some(1000));
        assert_eq!(parse_int(" -3 "), Some(-3));
        for bad in ["", "+", "1.0", "1e3", "_1", "0x10"] {
            assert_eq!(parse_int(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn split_uses_python_whitespace() {
        assert_eq!(split(" a\tb  c\x1fd\u{a0}e\r\n"), ["a", "b", "c", "d", "e"]);
        assert_eq!(rstrip("a b \x1c\r\n"), "a b");
        assert!(split("").is_empty());
    }

    #[test]
    fn prefix_len_follows_int_and_slicing() {
        assert_eq!(prefix_len(1.0 - 0.7, 10), Some(3));
        assert_eq!(prefix_len(0.3, 10), Some(3));
        // 0.7 * 10 == 7.000000000000001 truncates to 7; 0.29 * 100 == 28.999999999999996.
        assert_eq!(prefix_len(0.7, 10), Some(7));
        assert_eq!(prefix_len(0.29, 100), Some(28));
        assert_eq!(prefix_len(2.0, 10), Some(10));
        assert_eq!(prefix_len(-0.25, 10), Some(8));
        assert_eq!(prefix_len(-5.0, 10), Some(0));
        // int(-0.5) is 0, not -1.
        assert_eq!(prefix_len(-0.05, 10), Some(0));
        assert_eq!(prefix_len(0.5, 0), Some(0));
        assert_eq!(prefix_len(f64::NAN, 3), None);
        assert_eq!(prefix_len(f64::INFINITY, 3), None);
        // Python: int(inf * 0) raises too, since inf * 0 is nan.
        assert_eq!(prefix_len(f64::INFINITY, 0), None);
    }
}

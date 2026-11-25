//! Python's float text: `repr(float)` (and `str()`, f-strings without a format spec).

/// `repr(x)`: the shortest digits that round-trip, positional when `1e-4 <= |x| < 1e16`
/// (with at least one fractional digit, `2.0`), otherwise scientific with a signed exponent
/// of at least two digits (`1e-05`, `1.5e+16`); `nan`, `inf`, `-inf`.
pub fn py_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0".into() } else { "0.0".into() };
    }
    // Rust's `{:e}` gives the shortest round-trip length. Where two strings of that length
    // round-trip, Rust takes the larger and Python the correctly rounded one (ties to even),
    // which `{:.*e}` gives.
    let shortest = format!("{:e}", x.abs());
    let len = shortest.split_once('e').expect("{:e} has an exponent").0.replace('.', "").len();
    let sci = format!("{:.*e}", len - 1, x.abs());
    let (mantissa, exp) = sci.split_once('e').expect("{:e} has an exponent");
    let exp: i32 = exp.parse().expect("{:e} exponent is an integer");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    // Position of the decimal point relative to the digits, as in CPython's `decpt`.
    let decpt = exp + 1;
    let mut out = String::with_capacity(digits.len() + 8);
    if x < 0.0 {
        out.push('-');
    }
    if -4 < decpt && decpt <= 16 {
        if decpt <= 0 {
            out.push_str("0.");
            out.extend(std::iter::repeat_n('0', (-decpt) as usize));
            out.push_str(&digits);
        } else if decpt as usize >= digits.len() {
            out.push_str(&digits);
            out.extend(std::iter::repeat_n('0', decpt as usize - digits.len()));
            out.push_str(".0");
        } else {
            out.push_str(&digits[..decpt as usize]);
            out.push('.');
            out.push_str(&digits[decpt as usize..]);
        }
    } else {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        out.push_str(&format!("{:02}", exp.abs()));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn examples() {
        let cases = [
            (0.1, "0.1"),
            (2.0, "2.0"),
            (1e-4, "0.0001"),
            (1e-5, "1e-05"),
            (1.5e16, "1.5e+16"),
            (1e16, "1e+16"),
            (9999999999999998.0, "9999999999999998.0"),
            (-0.0, "-0.0"),
            (123456.789, "123456.789"),
        ];
        for (x, want) in cases {
            assert_eq!(py_repr(x), want, "{x:e}");
        }
    }
}

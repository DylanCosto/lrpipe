//! Python's string whitespace handling and `int()` parsing.

/// Python's `str.isspace()` for one character
pub fn py_isspace(c: char) -> bool {
    c.is_whitespace() || ('\x1c'..='\x1f').contains(&c)
}

/// Python's `str.strip()`
pub fn py_strip(s: &str) -> &str {
    s.trim_matches(py_isspace)
}

/// Python's `int(s)` for a decimal literal
pub fn py_int(s: &str) -> Option<i128> {
    let t = py_strip(s);
    let (negative, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if digits.is_empty() || digits.starts_with('_') || digits.ends_with('_') || digits.contains("__") {
        return None;
    }
    let mut value: i128 = 0;
    for c in digits.chars().filter(|&c| c != '_') {
        value = value.checked_mul(10)?.checked_add(i128::from(c.to_digit(10)?))?;
    }
    Some(if negative { -value } else { value })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_strings() {
        assert_eq!(py_strip("\u{1f} x\t\u{a0}"), "x");
        assert_eq!(py_int("100\n"), Some(100));
        assert_eq!(py_int(" -1_000 "), Some(-1000));
        assert_eq!(py_int("1__0"), None);
        assert_eq!(py_int(""), None);
        assert_eq!(py_int("1.5"), None);
    }
}

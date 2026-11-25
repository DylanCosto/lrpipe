//! The Python `str` and text-file semantics `SortVcf.py` relies on, as PyPy 3.11 implements them.
//!
//! The original runs `SortVcf` under PyPy, whose whitespace sets differ slightly from CPython's,
//! so these helpers follow PyPy.

use std::collections::VecDeque;
use std::io::{self, BufRead};

/// Whether `str.split()` (no separator) splits on `c` in PyPy: Unicode `White_Space` plus U+180E.
fn is_split_space(c: char) -> bool {
    c.is_whitespace() || c == '\u{180e}'
}

/// Whether `str.strip()` removes `c` in PyPy: the split set plus U+001C..=U+001F, which PyPy's
/// `split()` does not split on.
fn is_strip_space(c: char) -> bool {
    is_split_space(c) || ('\u{1c}'..='\u{1f}').contains(&c)
}

/// `s.strip()`.
pub(crate) fn strip(s: &str) -> &str {
    s.trim_matches(is_strip_space)
}

/// `s.rstrip()`.
pub(crate) fn rstrip(s: &str) -> &str {
    s.trim_end_matches(is_strip_space)
}

/// `s.split()`, or `s.split(maxsplit=n)` when `max_split` is given.
pub(crate) fn split_whitespace(s: &str, max_split: Option<usize>) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut rest = s.trim_start_matches(is_split_space);
    while !rest.is_empty() {
        if max_split.is_some_and(|n| fields.len() == n) {
            // After the last split Python keeps the remainder, trailing whitespace included.
            fields.push(rest);
            break;
        }
        let end = rest.find(is_split_space).unwrap_or(rest.len());
        fields.push(&rest[..end]);
        rest = rest[end..].trim_start_matches(is_split_space);
    }
    fields
}

/// `len(s.split('\t'))`.
pub(crate) fn tab_fields(s: &str) -> usize {
    s.bytes().filter(|&b| b == b'\t').count() + 1
}

/// Strips the underscores PEP 515 allows between digits, or `None` if one is misplaced.
fn without_digit_underscores(s: &str) -> Option<String> {
    let bytes = s.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'_' {
            let digit_before = i > 0 && bytes[i - 1].is_ascii_digit();
            let digit_after = bytes.get(i + 1).is_some_and(u8::is_ascii_digit);
            if !(digit_before && digit_after) {
                return None;
            }
        }
    }
    Some(s.replace('_', ""))
}

/// `int(s)` for ASCII decimal input, or `None` where Python raises `ValueError`.
///
/// Values outside `i64` (which Python accepts) are also `None`.
pub(crate) fn parse_int(s: &str) -> Option<i64> {
    let s = without_digit_underscores(strip(s))?;
    let digits = s.strip_prefix(['+', '-']).unwrap_or(&s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse().ok()
}

/// `float(s)`, or `None` where Python raises `ValueError`.
///
/// Rust's float grammar matches Python's (`inf`, `infinity` and `nan` in any case, optional
/// sign and exponent) once PEP 515 underscores are removed.
pub(crate) fn parse_float(s: &str) -> Option<f64> {
    without_digit_underscores(strip(s))?.parse().ok()
}

/// Converts a Python character index (negative counts from the end) to a clamped byte offset.
fn byte_offset(s: &str, index: isize) -> usize {
    let chars = s.chars().count();
    let index = if index < 0 { chars.saturating_sub(index.unsigned_abs()) } else { index.unsigned_abs().min(chars) };
    s.char_indices().nth(index).map_or(s.len(), |(i, _)| i)
}

/// `s[:index]`.
pub(crate) fn head(s: &str, index: isize) -> &str {
    &s[..byte_offset(s, index)]
}

/// `s[index:]`.
pub(crate) fn tail(s: &str, index: isize) -> &str {
    &s[byte_offset(s, index)..]
}

/// `len(s)` in characters.
pub(crate) fn char_len(s: &str) -> isize {
    isize::try_from(s.chars().count()).expect("string length fits in isize")
}

/// `haystack.find(needle)`: the character index of the first match, or -1.
pub(crate) fn find(haystack: &str, needle: &str) -> isize {
    haystack.find(needle).map_or(-1, |byte| char_len(&haystack[..byte]))
}

/// `str` read in text mode: universal newlines (`\r\n` and `\r` become `\n`).
pub(crate) fn translate_newlines(s: &str) -> String {
    s.replace("\r\n", "\n").replace('\r', "\n")
}

/// Why [`TextLines`] could not produce a line.
#[derive(Debug)]
pub(crate) enum TextError {
    /// Reading failed.
    Io(io::Error),
    /// The line is not UTF-8.
    Utf8,
}

/// The lines of a file opened in text mode (`open(path)`), as Python iterates them: split at
/// `\n`, `\r\n` or `\r`, each translated to a trailing `\n`; the last line keeps no terminator
/// if the file has none.
///
/// Lines are decoded as UTF-8 one at a time, so an invalid byte only fails the line holding it.
#[derive(Debug)]
pub(crate) struct TextLines<R> {
    reader: R,
    chunk: Vec<u8>,
    pending: VecDeque<Vec<u8>>,
}

impl<R: BufRead> TextLines<R> {
    pub(crate) fn new(reader: R) -> Self {
        Self { reader, chunk: Vec::new(), pending: VecDeque::new() }
    }

    /// Splits one `read_until('\n')` chunk at its bare `\r`s. A `\r` can only end a chunk at
    /// end of file, so a `\r` directly before the final `\n` is the only CRLF.
    fn queue_chunk(&mut self) {
        let mut rest: &[u8] = &self.chunk;
        while let Some(i) = rest.iter().position(|&b| b == b'\r') {
            let mut line = rest[..i].to_vec();
            line.push(b'\n');
            self.pending.push_back(line);
            rest = &rest[i + 1..];
            if rest == b"\n" {
                rest = &[];
            }
        }
        if !rest.is_empty() {
            self.pending.push_back(rest.to_vec());
        }
    }
}

impl<R: BufRead> Iterator for TextLines<R> {
    type Item = Result<String, TextError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pending.is_empty() {
            self.chunk.clear();
            match self.reader.read_until(b'\n', &mut self.chunk) {
                Ok(0) => return None,
                Ok(_) => self.queue_chunk(),
                Err(e) => return Some(Err(TextError::Io(e))),
            }
        }
        let line = self.pending.pop_front()?;
        Some(String::from_utf8(line).map_err(|_| TextError::Utf8))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_matches_python() {
        assert_eq!(split_whitespace("  a \tb\n", None), ["a", "b"]);
        assert_eq!(split_whitespace("a b c d e ", Some(3)), ["a", "b", "c", "d e "]);
        assert_eq!(split_whitespace("a b", Some(3)), ["a", "b"]);
        assert!(split_whitespace(" \n", None).is_empty());
        // PyPy splits on U+180E but not on U+001C, which strip() removes.
        assert_eq!(split_whitespace("a\u{180e}b\u{1c}c", None), ["a", "b\u{1c}c"]);
        assert_eq!(strip("\u{1c} a\n"), "a");
    }

    #[test]
    fn numbers_match_python() {
        assert_eq!(parse_int("1_000"), Some(1000));
        assert_eq!(parse_int("+07"), Some(7));
        assert_eq!(parse_int("-3"), Some(-3));
        for bad in ["", "_1", "1_", "1__0", "1.0", "0x10", "+", "1e3"] {
            assert_eq!(parse_int(bad), None, "{bad:?}");
        }
        assert_eq!(parse_float("13.58"), Some(13.58));
        assert_eq!(parse_float("1_0.5"), Some(10.5));
        assert_eq!(parse_float("-Infinity"), Some(f64::NEG_INFINITY));
        assert!(parse_float("nan").is_some_and(f64::is_nan));
        assert_eq!(parse_float("."), None);
        assert_eq!(parse_float("1._5"), None);
    }

    #[test]
    fn slicing_matches_python() {
        assert_eq!(head("abc", -1), "ab");
        assert_eq!(tail("abc", -1), "c");
        assert_eq!(head("abc", 10), "abc");
        assert_eq!(tail("abc", 10), "");
        assert_eq!(head("abc", -10), "");
        assert_eq!(tail("é1é2", 2), "é2");
        assert_eq!(find("éxchr", "chr"), 2);
        assert_eq!(find("abc", "z"), -1);
        assert_eq!(find("abc", ""), 0);
    }

    #[test]
    fn text_lines_use_universal_newlines() {
        let lines: Vec<String> =
            TextLines::new(&b"a\r\nb\rc\n\nd\r"[..]).collect::<Result<_, _>>().expect("valid UTF-8");
        assert_eq!(lines, ["a\n", "b\n", "c\n", "\n", "d\n"]);
        let lines: Vec<String> = TextLines::new(&b"x\ny"[..]).collect::<Result<_, _>>().expect("valid UTF-8");
        assert_eq!(lines, ["x\n", "y"]);
        let mut bad = TextLines::new(&b"ok\n\xff\n"[..]);
        assert_eq!(bad.next().expect("a line").expect("valid"), "ok\n");
        assert!(matches!(bad.next(), Some(Err(TextError::Utf8))));
    }
}

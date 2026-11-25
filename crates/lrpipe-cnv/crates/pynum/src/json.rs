//! Python's `json.dump` text (`ensure_ascii=True`, `allow_nan=True`), with an `indent`.

use std::fmt::Write as _;

use crate::repr::py_repr;

/// A Python value as `json.dump` sees it. Dicts keep insertion order.
#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    List(Vec<Json>),
    Dict(Vec<(String, Json)>),
}

impl Json {
    /// `json.dump(self, f, indent=indent)`: nested containers on their own lines, `", "` and
    /// `": "` separators, empty containers as `[]` and `{}`, floats as `repr`, and NaN and
    /// infinities as `NaN`, `Infinity` and `-Infinity`. No trailing newline.
    pub fn write_indented(&self, indent: &str, out: &mut String) {
        self.write(indent, 0, out);
    }

    fn write(&self, indent: &str, level: usize, out: &mut String) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Int(i) => {
                let _ = write!(out, "{i}");
            }
            Json::Float(f) => out.push_str(&float(*f)),
            Json::Str(s) => string(s, out),
            Json::List(items) if items.is_empty() => out.push_str("[]"),
            Json::Dict(items) if items.is_empty() => out.push_str("{}"),
            Json::List(items) => {
                out.push('[');
                for (k, item) in items.iter().enumerate() {
                    out.push_str(if k == 0 { "\n" } else { ",\n" });
                    push_indent(indent, level + 1, out);
                    item.write(indent, level + 1, out);
                }
                out.push('\n');
                push_indent(indent, level, out);
                out.push(']');
            }
            Json::Dict(items) => {
                out.push('{');
                for (k, (key, value)) in items.iter().enumerate() {
                    out.push_str(if k == 0 { "\n" } else { ",\n" });
                    push_indent(indent, level + 1, out);
                    string(key, out);
                    out.push_str(": ");
                    value.write(indent, level + 1, out);
                }
                out.push('\n');
                push_indent(indent, level, out);
                out.push('}');
            }
        }
    }
}

fn push_indent(indent: &str, level: usize, out: &mut String) {
    for _ in 0..level {
        out.push_str(indent);
    }
}

/// `json`'s `floatstr`: `repr`, except for the non-finite values.
fn float(f: f64) -> String {
    if f.is_nan() {
        "NaN".into()
    } else if f == f64::INFINITY {
        "Infinity".into()
    } else if f == f64::NEG_INFINITY {
        "-Infinity".into()
    } else {
        py_repr(f)
    }
}

/// `json`'s `py_encode_basestring_ascii`.
fn string(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (' '..='~').contains(&c) => out.push(c),
            c => {
                let mut buf = [0u16; 2];
                for unit in c.encode_utf16(&mut buf) {
                    let _ = write!(out, "\\u{unit:04x}");
                }
            }
        }
    }
    out.push('"');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn like_python() {
        let v = Json::Dict(vec![
            ("a".into(), Json::List(vec![Json::Int(1), Json::Float(f64::NAN), Json::Float(0.1)])),
            ("b".into(), Json::Dict(vec![])),
            ("c".into(), Json::List(vec![])),
            ("d\u{e9}".into(), Json::Str("x\"y\u{1}".into())),
        ]);
        let mut s = String::new();
        v.write_indented("\t", &mut s);
        // json.dumps({"a": [1, float("nan"), 0.1], "b": {}, "c": [], "dé": 'x"y\x01'}, indent="\t")
        assert_eq!(
            s,
            "{\n\t\"a\": [\n\t\t1,\n\t\tNaN,\n\t\t0.1\n\t],\n\t\"b\": {},\n\t\"c\": [],\n\t\"d\\u00e9\": \"x\\\"y\\u0001\"\n}"
        );
    }
}

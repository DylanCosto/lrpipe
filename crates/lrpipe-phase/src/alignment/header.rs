//! The header of a BAM/CRAM file: its text and its reference sequences.

use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct Header {
    pub text: String,
    /// (name, length)
    pub refs: Vec<(String, u32)>,
}

impl Header {
    pub fn tid(&self, name: &str) -> Option<usize> {
        self.refs.iter().position(|(n, _)| n == name)
    }

    pub fn target_names(&self) -> impl Iterator<Item = &str> {
        self.refs.iter().map(|(n, _)| n.as_str())
    }

    /// Header records by type (`CO` lines skipped), each as (tag, value)
    /// pairs in order; a later tag replaces an earlier one with the same name.
    pub fn records(&self) -> HashMap<String, Vec<Vec<(String, String)>>> {
        let mut map: HashMap<String, Vec<Vec<(String, String)>>> = HashMap::new();
        for line in self.text.split('\n').filter(|l| !l.is_empty()) {
            let mut parts = line.split('\t').filter(|p| !p.is_empty());
            let Some(first) = parts.next() else { continue };
            let Some(typ) = first.strip_prefix('@').filter(|t| t.len() >= 2).map(|t| &t[..2]) else { continue };
            if typ == "CO" {
                continue;
            }
            let mut fields: Vec<(String, String)> = Vec::new();
            for part in parts {
                let Some((tag, value)) = part.split_once(':') else { continue };
                match fields.iter_mut().find(|(t, _)| t == tag) {
                    Some(f) => f.1 = value.to_string(),
                    None => fields.push((tag.to_string(), value.to_string())),
                }
            }
            map.entry(typ.to_string()).or_default().push(fields);
        }
        map
    }
}

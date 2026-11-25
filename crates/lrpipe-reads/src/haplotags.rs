//! Reads' haplotypes by name, from a haplotag list: WhatsHap's `haplotag
//! --output-haplotag-list` (and `lrpipe haplotag`'s), with a read name, haplotype (`H1`, `H2`,
//! or `none`), phase set and chromosome per line, tab-separated, plain or gzipped.
//!
//! A tool given a list tags each read as `haplotag` would have tagged it in the alignments: the
//! primary alignment of a read the list assigns gets its haplotype and phase set, and every
//! other alignment (supplementary, secondary, or of a read the list does not assign) gets none,
//! whatever tags the input carries. That is `haplotag`'s default (`--tag-supplementary skip`),
//! whose list names primary alignments only.
//!
//! A read is looked up by its name on the contig its primary alignment is on, as `haplotag`
//! assigns them: two reads of one name on different contigs keep their own haplotypes. A list
//! written with `--tag-supplementary` also has lines for supplementary alignments, which do not
//! say which line is the primary's: a read whose lines on a contig disagree is left untagged
//! and counted ([`Haplotags::ambiguous`]), rather than given what may be a supplementary
//! alignment's haplotype. Reads the list names without a haplotype (`none`) are kept too, so
//! that a tool can say how many of its reads the list names at all ([`Matches`]): a list
//! written for other reads names none of them.

use std::sync::atomic::{AtomicU64, Ordering};

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

/// A read's haplotype and phase set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Assignment {
    /// Its haplotype, from 1 (`H1`).
    pub haplotype: u8,
    /// The phase set (`PS`): the position that names its phase block.
    pub phase_set: i64,
}

/// One listed read of a contig.
#[derive(Clone, Copy, Debug)]
struct Entry {
    hash: u64,
    /// Where its name starts in the contig's names.
    name_at: u32,
    name_len: u16,
    /// From 1; 0 for `none`.
    haplotype: u8,
    phase_set: i64,
}

/// The assigned reads of one contig, sorted by the hash of their names.
#[derive(Debug, Default)]
struct Contig {
    names: Vec<u8>,
    entries: Vec<Entry>,
}

/// Why a haplotag list could not be read.
#[derive(Debug, thiserror::Error)]
pub enum HaplotagsError {
    /// The file could not be opened or read.
    #[error("reading {}", path.display())]
    Io {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A line is not a read name, haplotype, phase set and chromosome.
    #[error("{}, line {line}: {reason}", path.display())]
    Format {
        /// The file.
        path: PathBuf,
        /// Its line, from 1.
        line: usize,
        /// What is wrong with it.
        reason: String,
    },
}

/// The haplotypes a haplotag list assigns, by contig and read name.
///
/// Its memory is about 24 bytes and a read name per listed read (some 60 bytes for ONT's
/// 36-character names); while it is read, up to about twice that.
#[derive(Default)]
pub struct Haplotags {
    contigs: HashMap<String, Contig>,
    assigned: usize,
    listed: usize,
    ambiguous: usize,
    matches: Matches,
}

/// The assigned reads of one contig of a [`Haplotags`].
#[derive(Clone, Copy)]
pub struct ContigHaplotags<'a>(&'a Contig);

impl std::fmt::Debug for Haplotags {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Haplotags").field("contigs", &self.contigs.len()).field("assigned", &self.assigned).finish()
    }
}

impl std::fmt::Debug for ContigHaplotags<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ContigHaplotags").field("assigned", &self.len()).finish()
    }
}

impl Haplotags {
    /// Reads a haplotag list, gzipped (told apart by its first bytes) or not. Of a read listed
    /// twice on the same contig, the first line is kept.
    ///
    /// # Errors
    ///
    /// If the file cannot be read, or a line that is neither empty nor a `#` comment does not
    /// have four tab-separated columns, a haplotype `H<n>` (n from 1) or `none`, and an integer
    /// phase set for an assigned read.
    pub fn read(path: &Path) -> Result<Self, HaplotagsError> {
        let io_error = |source| HaplotagsError::Io { path: path.into(), source };
        let mut file = BufReader::new(File::open(path).map_err(io_error)?);
        let gzipped = file.fill_buf().map_err(io_error)?.starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> =
            if gzipped { Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))) } else { Box::new(file) };
        Self::from_reader(reader).map_err(|e| match e {
            ParseError::Io(source) => HaplotagsError::Io { path: path.into(), source },
            ParseError::Format(line, reason) => HaplotagsError::Format { path: path.into(), line, reason },
        })
    }

    fn from_reader(mut reader: impl BufRead) -> Result<Self, ParseError> {
        let mut list = Haplotags::default();
        let mut line = Vec::new();
        let mut number = 0;
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line).map_err(ParseError::Io)? == 0 {
                break;
            }
            number += 1;
            let text = line.strip_suffix(b"\n").unwrap_or(&line);
            let text = text.strip_suffix(b"\r").unwrap_or(text);
            if text.is_empty() || text.starts_with(b"#") {
                continue;
            }
            let fail = |reason: String| ParseError::Format(number, reason);
            let fields: Vec<&[u8]> = text.split(|&b| b == b'\t').collect();
            let [name, haplotype, phase_set, contig, ..] = fields[..] else {
                return Err(fail(format!(
                    "{} columns, not the read name, haplotype, phase set and chromosome",
                    fields.len()
                )));
            };
            let text_of = |b: &[u8]| String::from_utf8_lossy(b).into_owned();
            let (haplotype, phase_set) = if haplotype == b"none" {
                (0, 0)
            } else {
                let haplotype = haplotype
                    .strip_prefix(b"H")
                    .and_then(|n| std::str::from_utf8(n).ok()?.parse::<u8>().ok())
                    .filter(|&n| n >= 1)
                    .ok_or_else(|| fail(format!("haplotype {:?} is not H1, H2, ... or none", text_of(haplotype))))?;
                let phase_set = std::str::from_utf8(phase_set)
                    .ok()
                    .and_then(|p| p.parse::<i64>().ok())
                    .ok_or_else(|| fail(format!("phase set {:?} is not an integer", text_of(phase_set))))?;
                (haplotype, phase_set)
            };
            let name_len = u16::try_from(name.len()).map_err(|_| fail("read name too long".into()))?;
            let contig = list.contigs.entry(text_of(contig)).or_default();
            let name_at = u32::try_from(contig.names.len()).map_err(|_| fail("too many read names".into()))?;
            contig.names.extend_from_slice(name);
            contig.entries.push(Entry { hash: hash(name), name_at, name_len, haplotype, phase_set });
        }
        for contig in list.contigs.values_mut() {
            // By hash, then name: a name's lines come together.
            let names = &contig.names;
            let name = |e: &Entry| &names[e.name_at as usize..e.name_at as usize + usize::from(e.name_len)];
            contig.entries.sort_by(|a, b| a.hash.cmp(&b.hash).then_with(|| name(a).cmp(name(b))));
            // In place: one entry per name, none for a name whose lines disagree.
            let entries = &mut contig.entries;
            let (mut kept, mut i) = (0, 0);
            while i < entries.len() {
                let same = |e: &Entry| e.hash == entries[i].hash && name(e) == name(&entries[i]);
                let end = i + entries[i..].iter().take_while(|e| same(e)).count();
                let first = entries[i];
                if entries[i..end].iter().all(|e| (e.haplotype, e.phase_set) == (first.haplotype, first.phase_set)) {
                    entries[kept] = first;
                    kept += 1;
                    list.assigned += usize::from(first.haplotype > 0);
                } else {
                    list.ambiguous += 1;
                }
                i = end;
            }
            entries.truncate(kept);
            entries.shrink_to_fit();
            contig.names.shrink_to_fit();
            list.listed += kept;
        }
        Ok(list)
    }

    /// The assigned reads of `contig`, if it has any.
    #[must_use]
    pub fn contig(&self, contig: &str) -> Option<ContigHaplotags<'_>> {
        self.contigs.get(contig).map(ContigHaplotags)
    }

    /// The haplotype and phase set of the read `name` whose primary alignment is on `contig`.
    #[must_use]
    pub fn get(&self, contig: &str, name: &[u8]) -> Option<Assignment> {
        self.contig(contig)?.get(name)
    }

    /// Reads assigned a haplotype.
    #[must_use]
    pub fn len(&self) -> usize {
        self.assigned
    }

    /// Reads the list names, with a haplotype or `none` (ambiguous ones aside).
    #[must_use]
    pub fn listed(&self) -> usize {
        self.listed
    }

    /// The count of the alignments looked up in the list, which its users keep.
    #[must_use]
    pub fn matches(&self) -> &Matches {
        &self.matches
    }

    /// Reads whose lines on a contig give different haplotypes or phase sets (a list written
    /// with `--tag-supplementary`), which are left untagged.
    #[must_use]
    pub fn ambiguous(&self) -> usize {
        self.ambiguous
    }

    /// Whether no read is assigned one.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.assigned == 0
    }

    /// The bytes the assignments take in memory, about.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        self.contigs
            .iter()
            .map(|(name, c)| name.len() + c.names.capacity() + c.entries.capacity() * size_of::<Entry>())
            .sum()
    }

    /// What a tool tagging an alignment as `haplotag` would tag it with: `name`'s assignment on
    /// `contig` for a primary alignment (`flags` without the secondary 0x100 and supplementary
    /// 0x800 bits), and none for any other. Counts the primary ones in [`Self::matches`].
    #[must_use]
    pub fn of_alignment(&self, contig: &str, name: &[u8], flags: u16) -> Option<Assignment> {
        if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            return None;
        }
        let listed = self.contig(contig);
        self.matches.count(listed.is_some_and(|c| c.lists(name)));
        listed?.get(name)
    }

    /// [`Self::of_alignment`] on a contig already looked up ([`Self::contig`]).
    #[must_use]
    pub fn of_contig_alignment(
        &self,
        listed: Option<ContigHaplotags<'_>>,
        name: &[u8],
        flags: u16,
    ) -> Option<Assignment> {
        if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            return None;
        }
        self.matches.count(listed.is_some_and(|c| c.lists(name)));
        listed?.get(name)
    }
}

/// SAM's secondary alignment flag.
const FLAG_SECONDARY: u16 = 0x100;
/// SAM's supplementary alignment flag.
const FLAG_SUPPLEMENTARY: u16 = 0x800;

impl ContigHaplotags<'_> {
    /// The haplotype and phase set of the read `name`.
    #[must_use]
    pub fn get(&self, name: &[u8]) -> Option<Assignment> {
        self.entry(name)
            .filter(|e| e.haplotype > 0)
            .map(|e| Assignment { haplotype: e.haplotype, phase_set: e.phase_set })
    }

    /// Whether the list names the read `name`, with a haplotype or `none`.
    #[must_use]
    pub fn lists(&self, name: &[u8]) -> bool {
        self.entry(name).is_some()
    }

    fn entry(&self, name: &[u8]) -> Option<&Entry> {
        let c = self.0;
        let h = hash(name);
        let first = c.entries.partition_point(|e| e.hash < h);
        c.entries[first..].iter().take_while(|e| e.hash == h).find(|e| {
            let at = e.name_at as usize;
            &c.names[at..at + usize::from(e.name_len)] == name
        })
    }

    /// Reads of the contig the list names.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.entries.len()
    }

    /// Whether none is.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.entries.is_empty()
    }
}

/// How many of a tool's primary alignments a haplotag list names, counted as they are looked
/// up (from any thread), for a warning when a list is not the reads': for another sample, or
/// another alignment of the reads.
#[derive(Debug, Default)]
pub struct Matches {
    looked_up: AtomicU64,
    listed: AtomicU64,
}

/// Below this share of the reads named by the list, [`Matches::warning`] warns.
pub const FEW_LISTED: f64 = 0.5;

impl Matches {
    /// Counts one primary alignment looked up, and whether the list names it.
    pub fn count(&self, listed: bool) {
        self.looked_up.fetch_add(1, Ordering::Relaxed);
        if listed {
            self.listed.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// The alignments looked up, and how many of them the list names.
    #[must_use]
    pub fn counts(&self) -> (u64, u64) {
        (self.looked_up.load(Ordering::Relaxed), self.listed.load(Ordering::Relaxed))
    }

    /// What to say of `list` once the reads are looked up: how many it named, with a warning
    /// when that is under [`FEW_LISTED`] of them (`haplotag` lists every primary alignment);
    /// nothing when none was looked up.
    #[must_use]
    pub fn report(&self, list: &Path) -> Option<String> {
        let (looked_up, listed) = self.counts();
        if looked_up == 0 {
            return None;
        }
        let share = listed as f64 / looked_up as f64;
        let mut text = format!("{} names {listed} of {looked_up} primary alignments", list.display());
        if share < FEW_LISTED {
            text = format!(
                "warning: {text}, which haplotag would all have listed: is it the list of these reads, \
                 aligned to this reference? The reads it does not name have no haplotype"
            );
        }
        Some(text)
    }
}

enum ParseError {
    Io(io::Error),
    Format(usize, String),
}

/// FNV-1a: read names are short, and looked up once per alignment.
fn hash(name: &[u8]) -> u64 {
    name.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3))
}

/// A list from its text.
#[cfg(test)]
fn parse(text: &str) -> Result<Haplotags, ParseError> {
    Haplotags::from_reader(io::Cursor::new(text.as_bytes().to_vec()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const LIST: &str = "#readname\thaplotype\tphaseset\tchromosome\n\
        r1\tH1\t100\tchr1\n\
        r2\tH2\t100\tchr1\n\
        r3\tnone\tnone\tchr1\n\
        r1\tH2\t500\tchr2\n\
        r2\tH2\t100\tchr1\n\
        r4\tnone\tnone\tchr1\n\
        r4\tH1\t100\tchr1\n";

    #[test]
    fn reads_are_found_by_contig_and_name_and_disagreeing_lines_tag_none() {
        let Ok(list) = parse(LIST) else { panic!("a valid list") };
        assert_eq!(list.get("chr1", b"r1"), Some(Assignment { haplotype: 1, phase_set: 100 }));
        assert_eq!(list.get("chr2", b"r1"), Some(Assignment { haplotype: 2, phase_set: 500 }));
        // Listed twice on chr1 alike: as once.
        assert_eq!(list.get("chr1", b"r2"), Some(Assignment { haplotype: 2, phase_set: 100 }));
        assert_eq!(list.get("chr1", b"r3"), None, "none is no haplotype");
        assert!(list.contig("chr1").is_some_and(|c| c.lists(b"r3") && !c.lists(b"r")));
        // A primary alignment's `none` and a supplementary one's H1 (`--tag-supplementary
        // independent-or-skip`), which the first line would have taken for the primary's,
        // whichever came first: the read is untagged.
        assert_eq!(list.get("chr1", b"r4"), None);
        assert!(!list.contig("chr1").is_some_and(|c| c.lists(b"r4")));
        assert_eq!(list.ambiguous(), 1);
        assert_eq!(list.get("chr1", b"r"), None);
        assert_eq!(list.get("chr3", b"r1"), None);
        assert_eq!((list.len(), list.listed()), (3, 4));
        assert!(list.heap_bytes() > 0);
    }

    #[test]
    fn a_list_naming_few_of_the_reads_is_warned_of() {
        let matches = Matches::default();
        for listed in [true, false, false] {
            matches.count(listed);
        }
        assert_eq!(matches.counts(), (3, 1));
        assert_eq!(Matches::default().report(Path::new("l.tsv")), None);
        let report = matches.report(Path::new("l.tsv")).expect("a report");
        assert!(report.starts_with("warning: l.tsv names 1 of 3 primary alignments"), "{report}");
        matches.count(true);
        matches.count(true);
        assert_eq!(matches.report(Path::new("l.tsv")).as_deref(), Some("l.tsv names 3 of 5 primary alignments"));
    }

    #[test]
    fn only_primary_alignments_are_tagged() {
        let Ok(list) = parse(LIST) else { panic!("a valid list") };
        assert_eq!(list.of_alignment("chr1", b"r1", 0x10).map(|a| a.haplotype), Some(1));
        assert_eq!(list.of_alignment("chr1", b"r1", 0x800), None);
        assert_eq!(list.of_alignment("chr1", b"r1", 0x100), None);
    }

    #[test]
    fn malformed_lines_are_refused() {
        for (text, reason) in [
            ("r1\tH1\t100\n", "3 columns"),
            ("r1\tHX\t100\tchr1\n", "haplotype \"HX\""),
            ("r1\tH0\t100\tchr1\n", "haplotype \"H0\""),
            ("r1\tH1\tnone\tchr1\n", "phase set \"none\""),
        ] {
            match parse(text) {
                Err(ParseError::Format(1, r)) => assert!(r.contains(reason), "{r}"),
                _ => panic!("{text:?} is refused"),
            }
        }
    }

    #[test]
    fn gzipped_lists_are_read() {
        let path = std::env::temp_dir().join(format!("lrpipe-reads-haplotags-{}.tsv.gz", std::process::id()));
        let mut gz = flate2::write::GzEncoder::new(File::create(&path).expect("file"), flate2::Compression::default());
        gz.write_all(LIST.as_bytes()).expect("written");
        gz.finish().expect("finished");
        let list = Haplotags::read(&path).expect("read");
        assert_eq!(list.get("chr2", b"r1").map(|a| a.phase_set), Some(500));
        let _ = std::fs::remove_file(path);
        assert!(matches!(Haplotags::read(Path::new("/nonexistent/list")), Err(HaplotagsError::Io { .. })));
    }
}

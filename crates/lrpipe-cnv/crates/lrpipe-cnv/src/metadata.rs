//! Regions excluded from analysis: reference N runs (`.mdr`) and a blacklist (spec §2.3).

use std::fs;
use std::path::Path;

use crate::Error;

/// Where an excluded region comes from. The two differ in more than origin: Spectre keeps
/// `.mdr` coordinates as integers snapped to the bin grid and blacklist coordinates as the
/// strings read (`metadataCollector.py:101-127, 194-214`), and later steps compare them as such.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Source {
    /// An N run from the `.mdr`, snapped outwards to the bin grid.
    NRun,
    /// A blacklist line's start and end text (unsnapped, possibly with surrounding blanks).
    Blacklist { start: String, end: String },
}

/// An excluded region: positions `start..=end` (both ends inclusive when masking bins).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    pub start: i64,
    pub end: i64,
    pub source: Source,
}

/// Excluded regions per contig, in Spectre's dict order: `.mdr` contigs in file order, then
/// blacklist contigs not in the `.mdr`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub contigs: Vec<(String, Vec<Region>)>,
}

impl Metadata {
    /// `FastaRef.load_metadata`: the `.mdr` N runs snapped to `bin_size`, merged with the
    /// blacklist, if given.
    ///
    /// # Errors
    ///
    /// If a file cannot be read or a line cannot be parsed.
    pub fn load(mdr: &Path, blacklist: Option<&Path>, bin_size: i64) -> Result<Self, Error> {
        let mut meta = Metadata::from_mdr(mdr, bin_size)?;
        if let Some(path) = blacklist {
            meta.merge_blacklist(read_blacklist(path)?);
        }
        Ok(meta)
    }

    /// `extact_metadata_from_mdr`: `NSPOS` lines, start snapped down and end up to the grid.
    fn from_mdr(path: &Path, bin_size: i64) -> Result<Self, Error> {
        let text = fs::read_to_string(path).map_err(Error::io(path))?;
        let mut meta = Metadata::default();
        for (i, line) in text.split_inclusive('\n').enumerate() {
            if !line.starts_with("NSPOS") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let int = |k: usize| -> Result<i64, Error> {
                fields.get(k).and_then(|f| py_int(f)).ok_or_else(|| Error::format(path, i + 1, "bad NSPOS line"))
            };
            let (mut start, mut end) = (int(2)?, int(3)?);
            if start % bin_size != 0 {
                start -= start % bin_size;
            }
            if end % bin_size != 0 {
                end += bin_size - end % bin_size;
            }
            meta.regions_mut(fields[1]).push(Region { start, end, source: Source::NRun });
        }
        Ok(meta)
    }

    /// `merge_metadata`: blacklist regions are appended to the contig's list unless an equal
    /// one is there (only exact duplicate blacklist lines are), new contigs go last.
    fn merge_blacklist(&mut self, blacklist: Vec<(String, Vec<Region>)>) {
        for (contig, regions) in blacklist {
            match self.contigs.iter_mut().find(|(c, _)| *c == contig) {
                Some((_, existing)) => {
                    for r in regions {
                        if !existing.contains(&r) {
                            existing.push(r);
                        }
                    }
                }
                None => self.contigs.push((contig, regions)),
            }
        }
    }

    fn regions_mut(&mut self, contig: &str) -> &mut Vec<Region> {
        let i = match self.contigs.iter().position(|(c, _)| c == contig) {
            Some(i) => i,
            None => {
                self.contigs.push((contig.to_owned(), Vec::new()));
                self.contigs.len() - 1
            }
        };
        &mut self.contigs[i].1
    }

    /// The excluded regions on `contig`.
    pub fn regions(&self, contig: &str) -> &[Region] {
        self.contigs.iter().find(|(c, _)| c == contig).map_or(&[], |(_, r)| r.as_slice())
    }
}

/// `extract_blacklisted_regions`: each line loses its last character (meant to be the
/// newline, so a last line without one loses a real character), and its first three
/// tab-separated fields are kept as text. Contigs in file order, regions in line order. `#`
/// comments and UCSC `track` and `browser` lines are skipped.
fn read_blacklist(path: &Path) -> Result<Vec<(String, Vec<Region>)>, Error> {
    let bytes = fs::read(path).map_err(Error::io(path))?;
    let text = String::from_utf8(bytes).map_err(|_| Error::format(path, 0, "not UTF-8"))?;
    let mut out: Vec<(String, Vec<Region>)> = Vec::new();
    for (i, line) in text.split_inclusive('\n').enumerate() {
        // UCSC track and browser lines and # comments carry no regions (Spectre fails on them).
        let first = line.split_whitespace().next().unwrap_or_default();
        if first.starts_with('#') || first == "track" || first == "browser" {
            continue;
        }
        let mut term = line.to_owned();
        term.pop();
        let fields: Vec<&str> = term.split('\t').collect();
        let &[chrom, start, end, ..] = &fields[..] else {
            return Err(Error::format(path, i + 1, "fewer than 3 fields"));
        };
        // Masking converts with `int()`, which crashes on anything that is not an integer.
        let parse = |s: &str| py_int(s).ok_or_else(|| Error::format(path, i + 1, format!("not an integer: {s:?}")));
        let region = Region {
            start: parse(start)?,
            end: parse(end)?,
            source: Source::Blacklist { start: start.into(), end: end.into() },
        };
        match out.iter_mut().find(|(c, _)| c == chrom) {
            Some((_, regions)) => regions.push(region),
            None => out.push((chrom.to_owned(), vec![region])),
        }
    }
    Ok(out)
}

/// `FastaRef.get_n_regions` (metadataCollector.py:25-99, 148-169): writes the `.mdr` report of
/// base counts and N runs of a FASTA (gzip-compressed if its path contains "gz").
///
/// Its quirks are kept (Q-MDR): a run of `n_size - 1` Ns qualifies, a run at a contig's
/// start is reported from 1, a run at its end is not reported, GC% divides by the last
/// contig's length only, and headers are split on spaces only.
///
/// The report is written beside `out` and renamed to it, so a run killed while writing leaves
/// no partial `.mdr` for the next to take as complete (it is made only when absent), and runs
/// making a shared one at once each put a whole one in place.
///
/// # Errors
///
/// If the FASTA cannot be read or the report written.
pub fn write_mdr(reference: &Path, n_size: i64, out: &Path) -> Result<(), Error> {
    let text = mdr_text(reference, n_size)?;
    // Each writer publishes a complete file, including when the path is shared over NFS.
    let host = fs::read_to_string("/proc/sys/kernel/hostname").unwrap_or_default();
    let since = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap_or_default();
    let mut name = out.file_name().unwrap_or_default().to_owned();
    name.push(format!(".{}-{}-{}.tmp", host.trim(), std::process::id(), since.as_nanos()));
    let tmp = out.with_file_name(name);
    let written = fs::File::create(&tmp)
        .and_then(|mut file| {
            std::io::Write::write_all(&mut file, text.as_bytes())?;
            file.sync_all()
        })
        .and_then(|()| fs::rename(&tmp, out));
    if written.is_err() {
        let _ = fs::remove_file(&tmp);
    }
    written.map_err(Error::io(out))
}

/// Checks that the N runs in an existing metadata file agree with the current FASTA and N-run size.
///
/// # Errors
/// Fails if either file cannot be read or the mask does not match the reference.
pub fn check_mdr(reference: &Path, n_size: i64, path: &Path) -> Result<(), Error> {
    fn n_runs(text: &str) -> Vec<Vec<&str>> {
        let mut rows: Vec<_> = text
            .lines()
            .filter(|line| line.starts_with("NSPOS\t"))
            .map(|line| line.split_whitespace().collect())
            .collect();
        rows.sort_unstable();
        rows.dedup();
        rows
    }
    let expected = mdr_text(reference, n_size)?;
    let metadata = fs::read_to_string(path).map_err(Error::io(path))?;
    if n_runs(&expected) != n_runs(&metadata) {
        return Err(Error::Input(format!(
            "CNV metadata {} does not match the N runs of reference {}; use metadata for this reference and N-run size",
            path.display(),
            reference.display()
        )));
    }
    Ok(())
}

fn mdr_text(reference: &Path, n_size: i64) -> Result<String, Error> {
    use std::io::BufRead;
    let mut file = fs::File::open(reference).map_err(Error::io(reference))?;
    // Compressed or not by its contents (gzip's magic bytes, bgzip's too), not its name, which
    // Spectre went by ("gz" in it): a bgzipped `ref.fa` was read as text.
    let mut magic = Vec::with_capacity(2);
    std::io::Read::read_to_end(&mut std::io::Read::take(&mut file, 2), &mut magic).map_err(Error::io(reference))?;
    std::io::Seek::rewind(&mut file).map_err(Error::io(reference))?;
    let mut reader: Box<dyn BufRead> = if magic == [0x1f, 0x8b] {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(std::io::BufReader::new(file))))
    } else {
        Box::new(std::io::BufReader::with_capacity(1 << 20, file))
    };
    let mut runs: Vec<(String, Vec<(i64, i64)>)> = Vec::new();
    let (mut i, mut j) = (0i64, 0i64);
    let mut counts = [0u64; 4]; // A, T, C, G
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if reader.read_until(b'\n', &mut buf).map_err(Error::io(reference))? == 0 {
            break;
        }
        // Python reads the file in text mode: "\r\n" ends a line as "\n" does.
        let line = buf.strip_suffix(b"\n").unwrap_or(&buf);
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.first().is_some_and(|&b| b == b'>' || b == b';') {
            let name = String::from_utf8_lossy(&line[1..]).split(' ').next().unwrap_or("").to_owned();
            match runs.iter_mut().find(|(c, _)| *c == name) {
                Some((_, r)) => r.clear(),
                None => runs.push((name, Vec::new())),
            }
            i = 0;
            j = 0;
            continue;
        }
        for &b in line {
            match b.to_ascii_uppercase() {
                b'N' => j += 1,
                c => {
                    if j - i >= n_size {
                        let (_, r) = runs
                            .last_mut()
                            .ok_or_else(|| Error::format(reference, 0, "sequence before the first header"))?;
                        r.push((i + 1, j - 1));
                    }
                    i = j;
                    j += 1;
                    match c {
                        b'A' => counts[0] += 1,
                        b'T' => counts[1] += 1,
                        b'C' => counts[2] += 1,
                        b'G' => counts[3] += 1,
                        _ => {}
                    }
                }
            }
        }
    }
    let [a, t, c, g] = counts;
    let gc = (g + c) as f64 / j as f64 * 100.0;
    let mut text = format!(
        "#####\tStatistic report\n#####\tBasepair statistic\nBPDEF\tA\tT\tC\tG\tGC%\tTotal bp\nBPSTA\t{a}\t{t}\t{c}\t{g}\t{}\t{j}\n",
        pynum::repr::py_repr(gc)
    );
    text.push_str("#####\tN-Sequence positions in given reference file\nNSDEF\tFrom\tTo\n");
    for (contig, r) in &runs {
        for (s, e) in r {
            text.push_str(&format!("NSPOS\t{contig}\t{s}\t{e}\n"));
        }
    }
    Ok(text)
}

/// Python's `int(str)` for plain decimal integers: surrounding whitespace allowed.
pub(crate) fn py_int(s: &str) -> Option<i64> {
    let t = s.trim();
    let digits = t.strip_prefix(['+', '-']).unwrap_or(t);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit() || b == b'_') {
        return None;
    }
    t.replace('_', "").parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str, text: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("spectre-metadata-{}-{name}", std::process::id()));
        fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn mdr_snaps_outwards_and_blacklist_keeps_text() {
        let mdr = tmp("a.mdr", "#####\tStatistic report\nNSPOS\tchr1\t1\t9999\nNSPOS\tchr1\t20000\t20500\n");
        let bl = tmp("a.bed", "chr1\t5\t100\tx\nchr1\t5\t100\tx\nchr2\t10\t2000\n");
        let m = Metadata::load(&mdr, Some(&bl), 1000).unwrap();
        let c1: Vec<(i64, i64)> = m.regions("chr1").iter().map(|r| (r.start, r.end)).collect();
        assert_eq!(c1, [(0, 10000), (20000, 21000), (5, 100)]);
        assert_eq!(m.regions("chr2")[0].source, Source::Blacklist { start: "10".into(), end: "2000".into() });
        // The last line had no newline, so its last character was cut.
        let bl2 = tmp("b.bed", "chr3\t10\t2000");
        let m = Metadata::load(&mdr, Some(&bl2), 1000).unwrap();
        assert_eq!(m.regions("chr3")[0].end, 200);
    }
}

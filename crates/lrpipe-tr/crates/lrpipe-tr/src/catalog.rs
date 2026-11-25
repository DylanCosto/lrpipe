//! Repeat catalogs: BED lines of `contig start end` and the locus, either TRGT-style
//! (`ID=..;MOTIFS=A,B;...`) or a bare motif list (`CAG` or `CAG,CAA`).

use std::io::BufRead;
use std::path::Path;

use crate::Error;

/// One repeat locus: `[start, end)` on `contig`, 0-based.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Locus {
    /// The catalog's `ID`, else `contig_start_end`.
    pub id: String,
    /// Reference contig name.
    pub contig: String,
    /// 0-based start of the repeat.
    pub start: i64,
    /// Exclusive end of the repeat.
    pub end: i64,
    /// Repeat units, uppercase.
    pub motifs: Vec<Vec<u8>>,
    /// The catalog's fourth column, as written (repeated in the VCF).
    pub info: String,
    /// Where another catalog locus ends at this one's start, the boundary they share.
    pub shared_start: Option<SharedBoundary>,
    /// Where another catalog locus starts at this one's end, the boundary they share.
    pub shared_end: Option<SharedBoundary>,
}

/// A position where catalog loci meet: some end there and others start there (as
/// ExpansionHunter's catalogs split a locus such as FXN or ATXN7 into adjacent repeats). A
/// read's insertion at it is split once between the two sides ([`crate::extract`]), not counted
/// in both.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SharedBoundary {
    /// Motif lengths of the loci ending here.
    pub left: Vec<usize>,
    /// Motif lengths of the loci starting here.
    pub right: Vec<usize>,
}

/// Reads a catalog (plain or gzip).
///
/// # Errors
///
/// If the file cannot be read, or a line lacks the four columns, valid coordinates or motifs.
pub fn read(path: &Path) -> Result<Vec<Locus>, Error> {
    let file = std::fs::File::open(path).map_err(Error::io(path))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        gzip_text(file).map_err(Error::io(path))?
    } else {
        Box::new(std::io::BufReader::new(file))
    };
    let mut out = Vec::new();
    // Each locus's line, for the duplicate check after reading.
    let mut lines: Vec<usize> = Vec::new();
    for (n, line) in reader.lines().enumerate() {
        let line = line.map_err(Error::io(path))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let bad = |m: &str| Error::Format { path: path.to_owned(), line: n + 1, message: m.to_owned() };
        let f: Vec<&str> = line.split('\t').collect();
        let [contig, start, end, info, ..] = f[..] else { return Err(bad("expected contig, start, end and motifs")) };
        let start: i64 = start.parse().map_err(|_| bad("bad start"))?;
        let end: i64 = end.parse().map_err(|_| bad("bad end"))?;
        if start < 0 || end <= start {
            return Err(bad("bad interval: the end must be after the start"));
        }
        let field = |key: &str| info.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
        let (id, motifs) = match field("MOTIFS") {
            Some(m) => (field("ID").map_or_else(|| format!("{contig}_{start}_{end}"), str::to_owned), m),
            None => (format!("{contig}_{start}_{end}"), info),
        };
        let motifs: Vec<Vec<u8>> =
            motifs.split(',').filter(|m| !m.is_empty()).map(|m| m.as_bytes().to_ascii_uppercase()).collect();
        // Bases or IUPAC codes (STRchive's GCN repeats).
        if motifs.is_empty() || motifs.iter().any(|m| !m.iter().all(|b| b"ACGTURYSWKMBDHVN".contains(b))) {
            return Err(bad("no valid motif: motifs are bases or IUPAC codes, such as CAG or GCN"));
        }
        lines.push(n + 1);
        out.push(Locus {
            id,
            contig: contig.to_owned(),
            start,
            end,
            motifs,
            info: info.to_owned(),
            shared_start: None,
            shared_end: None,
        });
    }
    let mut out = without_duplicates(out, &lines, path)?;
    mark_shared_boundaries(&mut out);
    Ok(out)
}

/// `loci` with each ID once: a locus given twice at the same place (two catalogs concatenated)
/// is kept once, at its first line; an ID at two places is an error naming the later line. By
/// sorting indices rather than hashing the IDs, a whole catalog's check copies no ID.
fn without_duplicates(loci: Vec<Locus>, lines: &[usize], path: &Path) -> Result<Vec<Locus>, Error> {
    let mut order: Vec<u32> = (0..loci.len() as u32).collect();
    order.sort_unstable_by(|&a, &b| loci[a as usize].id.cmp(&loci[b as usize].id).then(a.cmp(&b)));
    let mut drop = vec![false; loci.len()];
    for pair in order.windows(2) {
        let (first, later) = (&loci[pair[0] as usize], &loci[pair[1] as usize]);
        if first.id != later.id {
            continue;
        }
        if (&first.contig, first.start, first.end) != (&later.contig, later.start, later.end) {
            return Err(Error::Format {
                path: path.to_owned(),
                line: lines[pair[1] as usize],
                message: format!(
                    "locus ID {} is already {}:{}-{}: each locus needs its own ID",
                    later.id, first.contig, first.start, first.end
                ),
            });
        }
        drop[pair[1] as usize] = true;
    }
    if !drop.contains(&true) {
        return Ok(loci);
    }
    Ok(loci.into_iter().zip(drop).filter_map(|(locus, dropped)| (!dropped).then_some(locus)).collect())
}

/// Marks the boundaries where one locus ends and another starts. Marked on the whole catalog,
/// before any `--regions`, so a locus's record does not depend on which of its neighbours run.
fn mark_shared_boundaries(loci: &mut [Locus]) {
    use std::collections::HashMap;
    let mut sides: HashMap<(&str, i64), SharedBoundary> = HashMap::new();
    for l in loci.iter().filter(|l| l.end > l.start) {
        let lengths = l.motifs.iter().map(Vec::len);
        sides.entry((&l.contig, l.end)).or_default().left.extend(lengths.clone());
        sides.entry((&l.contig, l.start)).or_default().right.extend(lengths);
    }
    let shared: HashMap<(String, i64), SharedBoundary> = sides
        .into_iter()
        .filter(|(_, b)| !b.left.is_empty() && !b.right.is_empty())
        .map(|((contig, p), mut b)| {
            for v in [&mut b.left, &mut b.right] {
                v.sort_unstable();
                v.dedup();
            }
            ((contig.to_owned(), p), b)
        })
        .collect();
    for l in loci.iter_mut().filter(|l| l.end > l.start) {
        l.shared_start = shared.get(&(l.contig.clone(), l.start)).cloned();
        l.shared_end = shared.get(&(l.contig.clone(), l.end)).cloned();
    }
}

/// Keeps the loci that overlap `regions` (share a base with one), in their order. A locus is
/// kept whole, and genotyped from all the reads that reach it, however little of it a region
/// covers.
pub fn retain_in(loci: &mut Vec<Locus>, regions: &lrpipe_reads::regions::Regions) {
    loci.retain(|l| regions.overlaps(&l.contig, l.start, l.end));
}

/// The text of a gzip file: BGZF through pasta, other gzip (one member or several) through
/// flate2. A BGZF-only reader stops at a plain gzip file's first block (LRP-014).
fn gzip_text(file: std::fs::File) -> std::io::Result<Box<dyn BufRead>> {
    use std::io::BufReader;
    let mut reader = BufReader::new(file);
    let head = reader.fill_buf()?;
    let bgzf = head.len() >= 14 && head[..4] == [0x1f, 0x8b, 8, 4] && head[12..14] == *b"BC";
    Ok(if bgzf {
        Box::new(BufReader::new(pasta::bgzf::io::Reader::new(reader)))
    } else {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(reader)))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_trgt_and_bare_motif_lines() {
        let dir = std::env::temp_dir().join(format!("lrpipe-tr-catalog-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("c.bed");
        std::fs::write(&p, "chr4\t3074876\t3074933\tID=HTT;MOTIFS=CAG,CCG;STRUC=<TR>\nchr1\t10\t20\tca\n").unwrap();
        let loci = read(&p).unwrap();
        assert_eq!(loci[0].id, "HTT");
        assert_eq!(loci[0].motifs, vec![b"CAG".to_vec(), b"CCG".to_vec()]);
        assert_eq!((loci[1].id.as_str(), loci[1].motifs[0].as_slice()), ("chr1_10_20", &b"CA"[..]));
    }

    #[test]
    fn catalog_lines_are_checked() {
        // Edge-case testing: lines a catalog should not have were called.
        let dir = std::env::temp_dir().join(format!("lrpipe-tr-checks-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let read_text = |text: &str| {
            let p = dir.join("c.bed");
            std::fs::write(&p, text).unwrap();
            read(&p).map_err(|e| e.to_string())
        };
        let line = "chr1\t10\t20\tID=A;MOTIFS=CA\n";
        // UCSC browser and track lines, CRLF endings, IUPAC motifs and the same locus twice.
        let loci =
            read_text(&format!("browser position chr1\r\ntrack name=x\n{line}{line}chr1\t30\t40\tID=B;MOTIFS=GCN\r\n"))
                .expect("a valid catalog");
        assert_eq!(loci.iter().map(|l| l.id.as_str()).collect::<Vec<_>>(), ["A", "B"]);
        assert!(read_text("chr1\t10\t20\tID=A;MOTIFS=NNX\n").unwrap_err().contains("no valid motif"));
        assert!(read_text("chr1\t10\t10\tID=A;MOTIFS=CA\n").unwrap_err().contains("bad interval"));
        let other = format!("{line}chr2\t10\t20\tID=A;MOTIFS=CA\n");
        assert!(read_text(&other).unwrap_err().contains("locus ID A is already chr1:10-20"));
    }

    #[test]
    fn gzip_and_bgzf_catalogs_read_alike() {
        // LRP-014: a catalog compressed with gzip, not bgzip, reads as the plain one.
        use std::io::Write as _;
        let dir = std::env::temp_dir().join(format!("lrpipe-tr-gzip-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let text = "chr4\t3074876\t3074933\tID=HTT;MOTIFS=CAG,CCG;STRUC=<TR>\nchr1\t10\t20\tca\n";
        let plain = dir.join("c.bed");
        std::fs::write(&plain, text).unwrap();
        let gzip = dir.join("c.gzip.bed.gz");
        let mut w =
            flate2::write::GzEncoder::new(std::fs::File::create(&gzip).unwrap(), flate2::Compression::default());
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap();
        let bgzf = dir.join("c.bgzf.bed.gz");
        let mut w = pasta::bgzf::io::Writer::new(std::fs::File::create(&bgzf).unwrap());
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap();
        let want = read(&plain).unwrap();
        assert_eq!(read(&gzip).unwrap(), want);
        assert_eq!(read(&bgzf).unwrap(), want);
    }

    #[test]
    fn loci_overlapping_the_regions_are_kept_whole() {
        let locus = |contig: &str, start, end| Locus {
            id: format!("{contig}_{start}"),
            contig: contig.into(),
            start,
            end,
            motifs: vec![b"CA".to_vec()],
            info: "CA".into(),
            shared_start: None,
            shared_end: None,
        };
        let mut loci =
            vec![locus("chr1", 0, 10), locus("chr1", 95, 120), locus("chr1", 120, 130), locus("chr2", 0, 50)];
        let regions = lrpipe_reads::regions::Regions::from_intervals([("chr1".into(), 100, 120)]);
        retain_in(&mut loci, &regions);
        assert_eq!(loci, vec![locus("chr1", 95, 120)]);
    }

    #[test]
    fn adjacent_loci_share_their_boundary() {
        let dir = std::env::temp_dir().join(format!("lrpipe-tr-shared-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("c.bed");
        std::fs::write(
            &p,
            "chr3\t100\t130\tID=A;MOTIFS=GCA\nchr3\t130\t142\tID=B;MOTIFS=GCC,GCT\nchr3\t90\t150\tID=W;MOTIFS=GCA\nchr4\t130\t142\tID=C;MOTIFS=GCC\n",
        )
        .unwrap();
        let loci = read(&p).unwrap();
        let shared = SharedBoundary { left: vec![3], right: vec![3] };
        assert_eq!((loci[0].shared_start.as_ref(), loci[0].shared_end.as_ref()), (None, Some(&shared)));
        assert_eq!((loci[1].shared_start.as_ref(), loci[1].shared_end.as_ref()), (Some(&shared), None));
        // An overlapping locus, and one on another contig, share nothing.
        assert!(loci[2..].iter().all(|l| l.shared_start.is_none() && l.shared_end.is_none()));
    }
}

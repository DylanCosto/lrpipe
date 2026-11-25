//! mosdepth window depths (spec §2.1).

use std::io::BufRead;
use std::path::Path;

use crate::Error;

/// One contig's windows from a mosdepth `regions.bed.gz`, in file order.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContigDepth {
    pub name: String,
    /// 0-based window starts.
    pub positions: Vec<i64>,
    /// Mean depth per window, as printed by mosdepth.
    pub depth: Vec<f64>,
}

/// Reads every line of a mosdepth `regions.bed.gz` (BGZF or plain gzip), grouped by contig in
/// order of first appearance. Spectre reads it through a tabix index, one contig at a time,
/// which gives the same lines in the same order.
///
/// # Errors
///
/// If the file cannot be read, or a line does not have exactly four tab-separated fields
/// (Spectre's `[chrom, start, _, depth] = line.split("\t")`).
pub fn read_regions(path: &Path) -> Result<Vec<ContigDepth>, Error> {
    let file = std::fs::File::open(path).map_err(Error::io(path))?;
    parse_regions(gzip_text(file).map_err(Error::io(path))?, path)
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

/// Parses mosdepth `regions.bed` lines, as [`read_regions`]; `source` names them in errors.
///
/// # Errors
///
/// As [`read_regions`].
pub fn parse_regions(reader: impl BufRead, source: &Path) -> Result<Vec<ContigDepth>, Error> {
    let path = source;
    let mut out: Vec<ContigDepth> = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(Error::io(path))?;
        let bad = |m: &str| Error::format(path, i + 1, m);
        let fields: Vec<&str> = line.split('\t').collect();
        let &[chrom, start, _, depth] = &fields[..] else { return Err(bad("expected 4 tab-separated fields")) };
        let start = crate::metadata::py_int(start).ok_or_else(|| bad("bad window start"))?;
        let depth: f64 = depth.trim().parse().map_err(|_| bad("bad depth"))?;
        if out.last().is_none_or(|c| c.name != chrom) {
            if out.iter().any(|c| c.name == chrom) {
                return Err(bad("contig lines are not contiguous (is the file sorted?)"));
            }
            out.push(ContigDepth { name: chrom.to_owned(), ..Default::default() });
        }
        let c = out.last_mut().expect("pushed above");
        c.positions.push(start);
        c.depth.push(depth);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    #[test]
    fn gzip_and_bgzf_depths_read_alike() {
        // LRP-014: mosdepth writes BGZF, but a regions.bed.gz recompressed with gzip reads too.
        let dir = std::env::temp_dir().join(format!("lrpipe-cnv-gzip-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let text = "chr1\t0\t1000\t30.5\nchr1\t1000\t2000\t31.25\nchr2\t0\t1000\t12\n";
        let gzip = dir.join("gzip.regions.bed.gz");
        let mut w =
            flate2::write::GzEncoder::new(std::fs::File::create(&gzip).expect("file"), flate2::Compression::default());
        w.write_all(text.as_bytes()).expect("written");
        w.finish().expect("finished");
        let bgzf = dir.join("bgzf.regions.bed.gz");
        let mut w = pasta::bgzf::io::Writer::new(std::fs::File::create(&bgzf).expect("file"));
        w.write_all(text.as_bytes()).expect("written");
        w.finish().expect("finished");
        let (a, b) = (read_regions(&gzip).expect("gzip reads"), read_regions(&bgzf).expect("BGZF reads"));
        assert_eq!(a.len(), 2);
        for (x, y) in a.iter().zip(&b) {
            assert_eq!((&x.name, &x.positions, &x.depth), (&y.name, &y.positions, &y.depth));
        }
    }
}

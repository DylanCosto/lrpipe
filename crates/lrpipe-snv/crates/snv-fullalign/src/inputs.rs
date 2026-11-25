//! The two text inputs of a full-alignment job, parsed as `CreateTensorFullAlignmentFromCffi.py`
//! parses them: the `--full_aln_regions` BED and the `--phased_vcf_fn` VCF.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

use crate::{Error, FLANKING_BASES};

/// Python's `int(text)` for decimal strings: surrounding whitespace, an optional sign and
/// single underscores between digits are accepted.
pub(crate) fn py_int(text: &str) -> Option<i64> {
    let s = text.trim_matches(py_isspace);
    let (negative, digits) = match s.as_bytes().first()? {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || digits.starts_with('_') || digits.ends_with('_') || digits.contains("__") {
        return None;
    }
    let mut value: i64 = 0;
    for c in digits.chars().filter(|&c| c != '_') {
        let d = c.to_digit(10)?;
        value = value.checked_mul(10)?.checked_add(i64::from(d))?;
    }
    Some(if negative { -value } else { value })
}

/// `str.isspace` for one character (Python counts the ASCII separators 0x1c-0x1f as space).
pub(crate) fn py_isspace(c: char) -> bool {
    c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c)
}

/// Splits like Python's `str.rstrip()`.
pub(crate) fn py_rstrip(s: &str) -> &str {
    s.trim_end_matches(py_isspace)
}

/// The candidates of one full-alignment job, from its `--full_aln_regions` BED.
///
/// Each row `contig start end` on the job's contig names one candidate at the centre of the
/// interval; rows with a fourth column (`ref-alt-genotype-phase_set`, heterozygous SNPs kept for
/// haplotagging) only widen the region. The region spans every row on the contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CandidateRegions {
    contig: String,
    start: i64,
    end: i64,
    centers: Vec<u64>,
}

impl CandidateRegions {
    /// Reads the rows of `path` on `contig`.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read, a row on `contig` is malformed (where Python raises),
    /// or no row is on `contig` (the original's C code then fails to read the BAM and
    /// exits).
    pub fn from_bed(path: impl AsRef<Path>, contig: &str) -> Result<Self, Error> {
        let path = path.as_ref();
        let io_err = |source| Error::Io { path: path.to_owned(), source };
        let file = File::open(path).map_err(io_err)?;
        Self::from_reader(BufReader::new(file), contig, path)
    }

    fn from_reader(reader: impl BufRead, contig: &str, path: &Path) -> Result<Self, Error> {
        let mut span: Option<(i64, i64)> = None;
        let mut centers = Vec::new();
        for (n, line) in reader.lines().enumerate() {
            let line = line.map_err(|source| Error::Io { path: path.to_owned(), source })?;
            let bad = || Error::Bed { path: path.to_owned(), line: n + 1 };
            let row: Vec<&str> = py_rstrip(&line).split('\t').collect();
            if row[0] != contig {
                continue;
            }
            let field = |i: usize| row.get(i).and_then(|s| py_int(s)).ok_or_else(bad);
            let position = field(1)? + 1;
            let end = field(2)? + 1;
            span = Some(span.map_or((position, end), |(s, e)| (s.min(position), e.max(end))));
            if row.len() > 3 {
                // `ref_base, alt_base, genotype, phase_set = row[3].split('-')`.
                // Its `center_pos = position + extend_bp + 1` is never used.
                if row[3].split('-').count() != 4 {
                    return Err(bad());
                }
            } else if position == 1 {
                centers.push(end - FLANKING_BASES as i64 - 2);
            } else {
                centers.push(position + (end - position).div_euclid(2) - 1);
            }
        }
        let (start, end) =
            span.ok_or_else(|| Error::EmptyRegions { path: path.to_owned(), contig: contig.to_owned() })?;
        // `sorted(set(item - 1 for item in candidates_set if ctg_start <= item <= ctg_end))`.
        let mut zero_based = Vec::with_capacity(centers.len());
        for center in centers.into_iter().filter(|c| (start..=end).contains(c)) {
            // The list becomes a `size_t[]`, where cffi rejects negative values.
            let c = u64::try_from(center - 1).map_err(|_| Error::Position { position: center - 1 })?;
            zero_based.push(c);
        }
        zero_based.sort_unstable();
        zero_based.dedup();
        Ok(Self { contig: contig.to_owned(), start, end, centers: zero_based })
    }

    /// Builds the candidates directly: `centers` are 0-based positions and `start..=end` the
    /// 1-based region (`ctg_start`, `ctg_end`). Centers outside the region are dropped, the rest
    /// sorted and deduplicated, as Python does.
    pub fn new(contig: &str, start: i64, end: i64, centers: impl IntoIterator<Item = u64>) -> Self {
        let mut centers: Vec<u64> =
            centers.into_iter().filter(|&c| i64::try_from(c).is_ok_and(|c| (start..=end).contains(&(c + 1)))).collect();
        centers.sort_unstable();
        centers.dedup();
        Self { contig: contig.to_owned(), start, end, centers }
    }

    /// The contig.
    pub fn contig(&self) -> &str {
        &self.contig
    }

    /// 1-based first position of the region (`ctg_start`).
    pub fn start(&self) -> i64 {
        self.start
    }

    /// 1-based last position of the region (`ctg_end`).
    pub fn end(&self) -> i64 {
        self.end
    }

    /// 0-based candidate positions, ascending.
    pub fn centers(&self) -> &[u64] {
        &self.centers
    }
}

/// Which haplotype of a phased heterozygous SNP carries the alternative allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PhasedGenotype {
    /// `0|1`.
    AltOnSecond,
    /// Any other phased genotype (`1|0` and anything else with a `|`).
    AltOnFirst,
}

impl PhasedGenotype {
    /// The C `genotype` field: 1 for `0|1`, else 2.
    pub(crate) fn code(self) -> i32 {
        match self {
            Self::AltOnSecond => 1,
            Self::AltOnFirst => 2,
        }
    }
}

/// One phased heterozygous SNP (`struct Variant`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PhasedSnp {
    /// 0-based position.
    pub position: i32,
    /// Reference base.
    pub ref_base: u8,
    /// Alternative base.
    pub alt_base: u8,
    /// Phased genotype.
    pub genotype: PhasedGenotype,
    /// Phase set (`PS`, the last FORMAT field).
    pub phase_set: i32,
}

/// The phased heterozygous SNPs of one contig, in file order, used to haplotag reads.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PhasedSnps {
    snps: Vec<PhasedSnp>,
}

impl PhasedSnps {
    /// Reads the phased records on `contig` from a VCF, bgzipped or plain (`pigz -fdc`).
    /// Unphased records are skipped.
    ///
    /// # Errors
    ///
    /// Fails if the file cannot be read, a record is malformed where Python or cffi would raise
    /// (a missing column, a non-integer `POS` or `PS`, an allele that is not one character, a
    /// value outside a C `int`), or the positions decrease, which the original's C code does
    /// not handle.
    pub fn from_vcf(path: impl AsRef<Path>, contig: &str) -> Result<Self, Error> {
        let path = path.as_ref();
        let io_err = |source| Error::Io { path: path.to_owned(), source };
        let mut file = BufReader::new(File::open(path).map_err(io_err)?);
        let gzipped = file.fill_buf().map_err(io_err)?.starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> =
            if gzipped { Box::new(BufReader::new(MultiGzDecoder::new(file))) } else { Box::new(file) };
        Self::from_reader(reader, contig, path)
    }

    fn from_reader(mut reader: impl BufRead, contig: &str, path: &Path) -> Result<Self, Error> {
        let mut snps: Vec<PhasedSnp> = Vec::new();
        let mut line = Vec::new();
        let mut n = 0;
        loop {
            line.clear();
            n += 1;
            let read =
                reader.read_until(b'\n', &mut line).map_err(|source| Error::Io { path: path.to_owned(), source })?;
            if read == 0 {
                break;
            }
            let bad = || Error::Vcf { path: path.to_owned(), line: n };
            let text = String::from_utf8_lossy(&line);
            let row = py_rstrip(&text);
            match row.chars().next() {
                None => return Err(bad()),
                Some('#') => continue,
                Some(_) => {}
            }
            let columns: Vec<&str> = row.trim_matches(py_isspace).split('\t').collect();
            if !contig.is_empty() && columns[0] != contig {
                continue;
            }
            let column = |i: usize| columns.get(i).copied().ok_or_else(bad);
            let pos = py_int(column(1)?).ok_or_else(bad)?;
            let (ref_base, alt_base) = (column(3)?, column(4)?);
            let genotype_info: Vec<&str> = column(9)?.split(':').collect();
            let (genotype, phase_set) = (genotype_info[0], genotype_info[genotype_info.len() - 1]);
            if !genotype.contains('|') {
                continue;
            }
            let genotype = if genotype == "0|1" { PhasedGenotype::AltOnSecond } else { PhasedGenotype::AltOnFirst };
            let c_int = |v: i64| i32::try_from(v).map_err(|_| bad());
            let one_char = |s: &str| match s.as_bytes() {
                [b] => Ok(*b),
                _ => Err(bad()),
            };
            let snp = PhasedSnp {
                position: c_int(pos - 1)?,
                ref_base: one_char(ref_base)?,
                alt_base: one_char(alt_base)?,
                genotype,
                phase_set: c_int(py_int(phase_set).ok_or_else(bad)?)?,
            };
            if snps.last().is_some_and(|last| last.position > snp.position) {
                return Err(Error::UnsortedVcf { path: path.to_owned(), line: n });
            }
            snps.push(snp);
        }
        Ok(Self { snps })
    }

    /// Wraps SNPs given directly.
    ///
    /// # Errors
    ///
    /// Fails if the positions decrease, which the original's C code does not handle.
    pub fn new(snps: Vec<PhasedSnp>) -> Result<Self, Error> {
        if let Some(i) = snps.windows(2).position(|w| w[0].position > w[1].position) {
            return Err(Error::UnsortedVcf { path: PathBuf::new(), line: i + 2 });
        }
        Ok(Self { snps })
    }

    /// The SNPs, in file order.
    pub fn as_slice(&self) -> &[PhasedSnp] {
        &self.snps
    }

    /// Number of SNPs.
    pub fn len(&self) -> usize {
        self.snps.len()
    }

    /// True if there are no SNPs.
    pub fn is_empty(&self) -> bool {
        self.snps.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bed(text: &str) -> Result<CandidateRegions, Error> {
        CandidateRegions::from_reader(text.as_bytes(), "chr20", Path::new("test.bed"))
    }

    fn vcf(text: &str) -> Result<PhasedSnps, Error> {
        PhasedSnps::from_reader(text.as_bytes(), "chr20", Path::new("test.vcf"))
    }

    #[test]
    fn python_int() {
        assert_eq!(py_int(" 12\n"), Some(12));
        assert_eq!(py_int("+1_000"), Some(1000));
        assert_eq!(py_int("-7"), Some(-7));
        for bad in ["", "1__0", "_1", "1_", "1.5", "0x10", "-"] {
            assert_eq!(py_int(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn centers_and_region() {
        let r = bed("chr20\t100\t134\nchr1\t5\t39\nchr20\t0\t34\nchr20\t300\t334\tA-C-1-200\nchr20\t100\t134\n")
            .expect("valid");
        // position 101 -> 101 + 17 - 1 = 117 (0-based 116); position 1 -> 35 - 16 - 2 = 17.
        assert_eq!(r.centers(), [16, 116]);
        assert_eq!((r.start(), r.end()), (1, 335));
        assert!(matches!(bed("chr20\t1\t35\tA-C\n"), Err(Error::Bed { line: 1, .. })));
        assert!(matches!(bed("chr1\t1\t35\n"), Err(Error::EmptyRegions { .. })));
        assert!(matches!(bed("chr20\tx\t35\n"), Err(Error::Bed { .. })));
    }

    #[test]
    fn odd_widths_use_floor_division() {
        let r = bed("chr20\t10\t13\n").expect("valid");
        // position 11, end 14: 11 + 3 // 2 - 1 = 11.
        assert_eq!(r.centers(), [10]);
    }

    #[test]
    fn phased_records_only() {
        let text = "##x\n#CHROM\nchr20\t10\t.\tA\tG\t1\tPASS\t.\tGT:PS\t0|1:5\n\
                    chr20\t11\t.\tA\tG\t1\tPASS\t.\tGT:PS\t0/1:5\n\
                    chr1\t12\t.\tA\tG\t1\tPASS\t.\tGT:PS\t1|0:5\n\
                    chr20\t13\t.\tC\tT\t1\tPASS\t.\tGT:PS\t1|0:9\n";
        let snps = vcf(text).expect("valid");
        let got: Vec<_> = snps.as_slice().iter().map(|s| (s.position, s.genotype.code(), s.phase_set)).collect();
        assert_eq!(got, [(9, 1, 5), (12, 2, 9)]);
        assert!(vcf("chr20\t10\t.\tAC\tG\t1\tPASS\t.\tGT:PS\t0|1:5\n").is_err());
        assert!(vcf("chr20\t10\t.\tA\tG\t1\tPASS\t.\tGT:PS\t0|1:.\n").is_err());
        assert!(vcf("\n").is_err());
        let unsorted = "chr20\t10\t.\tA\tG\t1\tPASS\t.\tGT:PS\t0|1:1\nchr20\t9\t.\tA\tG\t1\tPASS\t.\tGT:PS\t0|1:1\n";
        assert!(matches!(vcf(unsorted), Err(Error::UnsortedVcf { line: 2, .. })));
    }
}

//! The SNV VCF and the copy-number-neutral depth estimated from heterozygous sites (spec
//! §2.4, §3.1; `SNVAnalysis`, snv_analysis.py:22-237).

use std::io::BufRead;
use std::path::Path;

use pynum::sum;

use crate::Error;
use crate::dip::unidip;
use crate::genome::Genome;

/// A genotype as pysam gives it: allele indexes, `None` for `.`.
pub type Genotype = Vec<Option<u32>>;

/// The FORMAT allele frequency of a record.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Af {
    /// One value. pysam gives a scalar only for `Number=1`; lrpipe-cnv also reads a single
    /// `Number=A` value (current Clair3 and lrpipe-snv) this way, where Spectre 0.2.1 skips the
    /// record and then fails.
    Value(f64),
    /// Several values (multi-allelic): Spectre skips the record.
    Several,
}

/// One VCF record, with the fields Spectre reads.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    /// 1-based POS.
    pub pos: i64,
    /// Reference span: `len(REF)`, or from `INFO/END`.
    pub rlen: i64,
    /// QUAL as htslib stores it (float32).
    pub qual: f64,
    pub gt: Genotype,
    pub dp: i64,
    pub af: Af,
}

impl Record {
    /// A REF/ALT heterozygote, in either order and phased or not (`0/1`, `1/0`, `0|1`, `1|0`).
    /// Spectre takes only `0/1` and `0|1`, so phased calls with ALT first count for nothing
    /// (LRP-011).
    pub fn is_het(&self) -> bool {
        self.gt == [Some(0), Some(1)] || self.gt == [Some(1), Some(0)]
    }

    pub fn is_hom_alt(&self) -> bool {
        self.gt == [Some(1), Some(1)]
    }
}

/// A single-sample SNV VCF, held in memory.
#[derive(Clone, Debug, Default)]
pub struct SnvVcf {
    /// Header contigs (`##contig` lines) in order, with their lengths.
    pub header_contigs: Vec<(String, Option<u64>)>,
    /// Records per contig, in file order; contigs in order of first record.
    pub contigs: Vec<(String, Vec<Record>)>,
    /// Per contig of [`contigs`](Self::contigs): the longest record span when its records are
    /// sorted by position, so [`region`](Self::region) can search; `None` when unsorted (or not
    /// yet indexed).
    spans: Vec<Option<i64>>,
}

impl SnvVcf {
    /// Reads a VCF (plain, gzip or BGZF). It must declare FORMAT `GT`, `DP` and `AF` (or
    /// DeepVariant's `VAF`) and have one sample.
    ///
    /// # Errors
    ///
    /// If the file cannot be read, lacks those tags, has several samples, or a record misses
    /// QUAL, DP or the allele frequency (all failures in Spectre).
    pub fn read(path: &Path) -> Result<Self, Error> {
        let file = std::fs::File::open(path).map_err(Error::io(path))?;
        let mut head = [0u8; 2];
        let gz = {
            use std::io::Read;
            let mut f = std::fs::File::open(path).map_err(Error::io(path))?;
            f.read(&mut head).map_err(Error::io(path))? == 2 && head == [0x1f, 0x8b]
        };
        let reader: Box<dyn BufRead> = if gz {
            Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(std::io::BufReader::new(file))))
        } else {
            Box::new(std::io::BufReader::new(file))
        };
        let mut vcf = SnvVcf::default();
        let mut formats: Vec<String> = Vec::new();
        let mut af_tag = None;
        for (n, line) in reader.lines().enumerate() {
            let line = line.map_err(Error::io(path))?;
            let bad = |m: &str| Error::format(path, n + 1, m);
            if let Some(meta) = line.strip_prefix("##") {
                if let Some(c) = meta.strip_prefix("contig=<") {
                    let id = field(c, "ID").ok_or_else(|| bad("##contig without ID"))?;
                    vcf.header_contigs.push((id.to_owned(), field(c, "length").and_then(|l| l.parse().ok())));
                } else if let Some(f) = meta.strip_prefix("FORMAT=<") {
                    formats.push(field(f, "ID").unwrap_or_default().to_owned());
                }
                continue;
            }
            if line.starts_with('#') {
                if line.split('\t').count() != 10 {
                    return Err(bad("expected exactly one sample"));
                }
                if !(formats.iter().any(|f| f == "GT") && formats.iter().any(|f| f == "DP")) {
                    return Err(bad("FORMAT GT and DP must be declared"));
                }
                af_tag = Some(if formats.iter().any(|f| f == "AF") {
                    "AF"
                } else if formats.iter().any(|f| f == "VAF") {
                    "VAF"
                } else {
                    return Err(bad("FORMAT AF (or VAF) must be declared"));
                });
                continue;
            }
            let af_tag = af_tag.ok_or_else(|| bad("record before the #CHROM line"))?;
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 10 {
                return Err(bad("expected 10 columns"));
            }
            let pos: i64 = f[1].parse().map_err(|_| bad("bad POS"))?;
            let end = f[7].split(';').find_map(|kv| kv.strip_prefix("END=")).and_then(|e| e.parse::<i64>().ok());
            let rlen = match end {
                Some(e) if e >= pos => e - pos + 1,
                _ => f[3].len() as i64,
            };
            // htslib stores QUAL and FORMAT floats as float32.
            let qual = f[5].parse::<f64>().map_err(|_| bad("QUAL is missing (Spectre fails on it)"))? as f32;
            let keys: Vec<&str> = f[8].split(':').collect();
            let values: Vec<&str> = f[9].split(':').collect();
            let get = |k: &str| keys.iter().position(|x| *x == k).and_then(|i| values.get(i).copied());
            let gt = get("GT").map_or_else(Vec::new, parse_gt);
            let dp =
                get("DP").and_then(|v| v.parse().ok()).ok_or_else(|| bad("DP is missing (Spectre fails on it)"))?;
            let af = match get(af_tag) {
                Some(v) if v.contains(',') => Af::Several,
                Some(v) => Af::Value(f64::from(
                    v.parse::<f64>().map_err(|_| bad("AF is missing (Spectre fails on it)"))? as f32,
                )),
                None => return Err(bad("AF is missing (Spectre fails on it)")),
            };
            let rec = Record { pos, rlen, qual: f64::from(qual), gt, dp, af };
            match vcf.contigs.last_mut() {
                Some((c, recs)) if c == f[0] => recs.push(rec),
                _ => vcf.contigs.push((f[0].to_owned(), vec![rec])),
            }
        }
        vcf.index();
        Ok(vcf)
    }

    /// The records of `contig`, in file order (pysam's `fetch(contig)`).
    pub fn records(&self, contig: &str) -> &[Record] {
        self.contigs.iter().find(|(c, _)| c == contig).map_or(&[], |(_, r)| r.as_slice())
    }

    /// The records overlapping the 1-based region `start..=end` of `contig`, in file order
    /// (pysam's `fetch("contig:start-end")`): a record covers `POS..POS + rlen - 1`. On sorted
    /// records only the positions that can overlap are visited.
    pub fn region(&self, contig: &str, start: i64, end: i64) -> impl Iterator<Item = &Record> {
        let at = self.contigs.iter().position(|(c, _)| c == contig);
        let records = at.map_or(&[][..], |i| self.contigs[i].1.as_slice());
        let (from, to) = match at.and_then(|i| self.spans.get(i).copied().flatten()) {
            // A record overlapping `start` begins after `start - span`; none past `end` does.
            Some(span) => {
                (records.partition_point(|r| r.pos <= start - span), records.partition_point(|r| r.pos <= end))
            }
            None => (0, records.len()),
        };
        records[from..to.max(from)].iter().filter(move |r| r.pos <= end && r.pos + r.rlen.max(1) > start)
    }

    /// Notes for each contig whether its records are sorted, and their longest span.
    fn index(&mut self) {
        self.spans = self
            .contigs
            .iter()
            .map(|(_, records)| {
                let sorted = records.windows(2).all(|w| w[0].pos <= w[1].pos);
                sorted.then(|| records.iter().map(|r| r.rlen.max(1)).max().unwrap_or(1))
            })
            .collect();
    }
}

/// The value of `key` in a structured header line's `<...>` body.
fn field<'a>(body: &'a str, key: &str) -> Option<&'a str> {
    body.trim_end_matches('>').split(',').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

/// `0/1`, `1|1`, `./.`, `1` as allele index lists.
fn parse_gt(s: &str) -> Genotype {
    s.split(['/', '|']).map(|a| a.parse().ok()).collect()
}

/// The copy-number-neutral depth (`cn_neutral_coverage_med`) and the depths it came from.
#[derive(Clone, Debug, PartialEq)]
pub struct CnNeutral {
    pub median: f64,
    pub mean: f64,
    pub stdev: f64,
    /// The heterozygous depths after cleaning (`perfect_het_depth`).
    pub het_depths: Vec<i64>,
}

const PERFECT_HET_AF: f64 = 0.5;
const MAX_DEVIATION_PERFECT_HET_AF: f64 = 0.12;
const MIN_SNV_QUAL: f64 = 5.0;
const MIN_COVERAGE_SITE: i64 = 3;
const MIN_HOM_ALT_AF: f64 = 0.75;
const MIN_HET_CHR_SITES: f64 = 0.4;

fn near_half(af: f64) -> bool {
    (PERFECT_HET_AF - af).abs() < MAX_DEVIATION_PERFECT_HET_AF
}

/// `snv_copy_number_state`: depths of "perfect" heterozygous SNVs (AF within 0.12 of 0.5,
/// QUAL > 5) on contigs where most passing sites are heterozygous, cleaned of extreme values,
/// then the median of the first unimodal interval UniDip finds (spec §3.1). `rng` is numpy's
/// global generator, which Spectre leaves unseeded.
///
/// # Errors
///
/// If there is no such SNV (Spectre fails on it).
pub fn cn_neutral(
    vcf: &SnvVcf,
    genome: &Genome,
    rng: &mut pynum::numpy_random::RandomState,
) -> Result<CnNeutral, Error> {
    // snv_nhet_chromosome: per contig, heterozygous and homozygous ALT sites with DP > 3.
    let mut het = std::collections::HashMap::<&str, u64>::new();
    let mut alt = std::collections::HashMap::<&str, u64>::new();
    for (contig, recs) in &vcf.contigs {
        for r in recs {
            let Af::Value(af) = r.af else { continue };
            if r.qual > MIN_SNV_QUAL && r.dp > MIN_COVERAGE_SITE {
                if r.is_het() {
                    if near_half(af) {
                        *het.entry(contig).or_default() += 1;
                    }
                } else if af >= MIN_HOM_ALT_AF {
                    *alt.entry(contig).or_default() += 1;
                }
            }
        }
    }
    let used: Vec<&str> = genome
        .contigs
        .iter()
        .map(|c| c.name.as_str())
        .filter(|c| {
            let (h, a) = (het.get(c).copied().unwrap_or(0), alt.get(c).copied().unwrap_or(0));
            h + a > 0 && h as f64 / (h + a) as f64 > MIN_HET_CHR_SITES
        })
        .collect();
    let mut depths: Vec<i64> = Vec::new();
    for (contig, recs) in &vcf.contigs {
        if !used.contains(&contig.as_str()) {
            continue;
        }
        for r in recs {
            if let Af::Value(af) = r.af
                && r.qual > MIN_SNV_QUAL
                && r.is_het()
                && near_half(af)
            {
                depths.push(r.dp);
            }
        }
    }
    if depths.is_empty() {
        return Err(Error::Input("no heterozygous SNV with an allele frequency near 0.5 (Spectre fails here)".into()));
    }
    // het_clean: drop depths above the 98th percentile when some depth is extreme.
    let as_f64 = |d: &[i64]| d.iter().map(|&x| x as f64).collect::<Vec<f64>>();
    let max = *depths.iter().max().expect("not empty") as f64;
    if max > 10.0 * sum::nanmean(&as_f64(&depths)) {
        let p98 = sum::percentile(&as_f64(&depths), 98.0);
        depths.retain(|&d| d as f64 <= p98);
    }
    let intervals = unidip(&as_f64(&depths), rng);
    // `depths[start:end]` with UniDip's inclusive end: one depth short (quirk Q-UNIDIP).
    let part: Vec<f64> = match intervals.first() {
        Some(&(s, e)) => as_f64(&depths[s.min(depths.len())..e.min(depths.len()).max(s.min(depths.len()))]),
        None => as_f64(&depths),
    };
    Ok(CnNeutral {
        median: sum::nanmedian(&part),
        mean: sum::nanmean(&part),
        stdev: sum::nanstd(&part),
        het_depths: depths,
    })
}

#[cfg(test)]
mod het_tests {
    use super::*;

    fn record(gt: &[Option<u32>]) -> Record {
        Record { pos: 1, rlen: 1, qual: 30.0, gt: gt.to_vec(), dp: 20, af: Af::Value(0.5) }
    }

    #[test]
    fn heterozygotes_do_not_depend_on_allele_order() {
        // LRP-011: 0/1, 1/0, 0|1 and 1|0 are the same evidence (the separator is not kept).
        assert!(record(&[Some(0), Some(1)]).is_het());
        assert!(record(&[Some(1), Some(0)]).is_het());
        for gt in [&[Some(1), Some(1)][..], &[Some(0), Some(0)], &[None, Some(1)], &[Some(1)], &[Some(1), Some(2)]] {
            assert!(!record(gt).is_het(), "{gt:?}");
        }
    }
}

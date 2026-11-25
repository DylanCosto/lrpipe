//! Changes to Spectre's method. Each is off by default, so a default run stays
//! Spectre-exact; [`Extensions::all`] turns on every one that needs no extra input.
//!
//! They act at three points: on the window depths before normalization ([`prepare`]: masking
//! ambiguously mapped windows, sex chromosome copy numbers, GC correction), on the outlier
//! masking during normalization ([`crate::normalize::OutlierRemoval::IsolatedWindows`]), and
//! on the calls' copy numbers after scoring ([`finish_copy_numbers`]).

use std::path::{Path, PathBuf};

use pynum::round::np_round;
use pynum::sum;

use crate::Error;
use crate::candidate::{Arena, CandidateId, gt};
use crate::coverage::{ContigDepth, read_regions};
use crate::gc;
use crate::genome::Genome;
use crate::hmm::HmmOptions;
use crate::sex::{Sex, SexDepth, sex_contigs, x_pars};

/// The changes to make.
#[derive(Clone, Debug, Default)]
pub struct Extensions {
    pub sex_chromosomes: SexChromosomes,
    pub mappability: Option<Mappability>,
    pub gc: GcCorrection,
    pub copy_numbers: CopyNumbers,
    pub caller: Caller,
    /// An SV caller's VCF (lrpipe-sv, Sniffles2) whose deletions and duplications give calls
    /// they overlap their base-pair breakpoints.
    pub sv_breakpoints: Option<PathBuf>,
}

/// How windows become calls.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum Caller {
    /// Spectre's runs beyond thresholds, merged and refined.
    #[default]
    Thresholds,
    /// The copy-number HMM, with allele fractions when SNVs are given ([`crate::hmm`]).
    Hmm(HmmOptions),
}

impl Extensions {
    /// Every change that needs no extra input: sex detection, GC correction and copy
    /// numbers from the final calls.
    pub fn all() -> Self {
        Extensions {
            sex_chromosomes: SexChromosomes::Detect,
            mappability: None,
            gc: GcCorrection::Median { span: DEFAULT_GC_SPAN },
            copy_numbers: CopyNumbers::FromCall,
            caller: Caller::Thresholds,
            sv_breakpoints: None,
        }
    }
}

/// How chrX and chrY are treated.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SexChromosomes {
    /// As autosomes, with the global ploidy (Spectre; `--ploidy-chr` still applies). A male
    /// sample's chrX then comes out as one long deletion.
    #[default]
    AsAutosomes,
    /// Read the sex from chrX and chrY depth; if it is not clear, as autosomes.
    Detect,
    /// This sex.
    Given(Sex),
}

/// Masking of windows where reads map ambiguously (segmental duplications and the like),
/// whose MAPQ ≥ 20 depth drops without any change in copy number.
#[derive(Clone, Debug)]
pub struct Mappability {
    /// Window depths over all reads, whatever their MAPQ (`mosdepth -x -b 1000 -Q 0`, the same
    /// windows as the main coverage).
    pub all_mapq_coverage: PathBuf,
    /// Windows whose depth is below this share of their all-MAPQ depth are masked.
    pub min_ratio: f64,
}

/// Correction of depth for GC content.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GcCorrection {
    /// None (Spectre).
    #[default]
    Off,
    /// Scale each window by the autosome median over the median of windows of its GC
    /// percentage, GC measured over `span` bases centred on the window. Long reads drop out
    /// with the GC content of a read's length, so the span should be near it.
    Median { span: i64 },
}

/// [`GcCorrection::Median`]'s span in [`Extensions::all`].
pub const DEFAULT_GC_SPAN: i64 = 20_000;

/// Where a call's copy number (VCF `CN`, `GT` and `CD`) comes from.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CopyNumbers {
    /// The median of the raw run of windows first called, before merging and refinement
    /// changed the call (Spectre's bug B3).
    #[default]
    FromRawRun,
    /// The median of the call's final windows.
    FromCall,
}

/// What [`prepare`] did, which later steps and the caller need.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Prepared {
    /// chrX and chrY depth relative to the autosomes, when sex chromosomes are handled.
    pub sex_depth: Option<SexDepth>,
    /// The sex used: given, or detected.
    pub sex: Option<Sex>,
    /// Contigs with one copy, and the regions of them with two (the PARs). Their depth was
    /// doubled (by the ploidy) and their calls' copy numbers are halved back.
    pub haploid: Vec<(String, Vec<(i64, i64)>)>,
    /// Contigs left out: chrY in a female sample.
    pub dropped: Vec<String>,
    /// Windows masked for ambiguous mapping.
    pub masked_windows: usize,
}

impl Prepared {
    /// Copies expected at window `pos` of `contig`, for a sample of ploidy `ploidy`.
    pub fn expected_copies(&self, contig: &str, pos: i64, ploidy: i64) -> i64 {
        match self.haploid.iter().find(|(c, _)| c == contig) {
            Some((_, diploid)) if !diploid.iter().any(|&(s, e)| pos >= s && pos < e) => 1,
            _ => ploidy,
        }
    }
}

/// Applies the depth changes, in order: mappability masking, sex chromosome copy numbers,
/// GC correction.
///
/// # Errors
///
/// If the all-MAPQ coverage or the reference cannot be read, or the all-MAPQ windows differ
/// from the coverage's.
pub fn prepare(
    depths: &mut Vec<ContigDepth>,
    ext: &Extensions,
    genome: &Genome,
    reference: &Path,
) -> Result<Prepared, Error> {
    let mut out = Prepared::default();
    let (x, y) = sex_contigs(genome.contigs.iter().map(|c| c.name.as_str()));
    let autosomes: Vec<String> = genome
        .contigs
        .iter()
        .map(|c| c.name.clone())
        .filter(|c| Some(c) != x.as_ref() && Some(c) != y.as_ref())
        .collect();

    if let Some(m) = &ext.mappability {
        out.masked_windows = mask_ambiguous(depths, m)?;
    }

    if ext.sex_chromosomes != SexChromosomes::AsAutosomes {
        let measured = SexDepth::measure(depths, x.as_deref(), y.as_deref(), &autosomes);
        out.sex_depth = Some(measured);
        out.sex = match ext.sex_chromosomes {
            SexChromosomes::Given(sex) => Some(sex),
            _ => measured.sex(),
        };
    }
    match out.sex {
        Some(Sex::Male) => {
            if let Some(x) = &x {
                let len = genome.contig(x).map_or(0, |c| c.len);
                out.haploid.push((x.clone(), x_pars(len)));
            }
            if let Some(y) = &y {
                out.haploid.push((y.clone(), Vec::new()));
            }
            for d in depths.iter_mut() {
                for (p, v) in d.positions.iter().zip(&mut d.depth) {
                    let copies = out.expected_copies(&d.name, *p, genome.ploidy);
                    *v *= (genome.ploidy / copies) as f64;
                }
            }
        }
        Some(Sex::Female) => {
            if let Some(y) = &y {
                depths.retain(|d| &d.name != y);
                out.dropped.push(y.clone());
            }
        }
        None => {}
    }

    if let GcCorrection::Median { span } = ext.gc {
        let bin = depths.iter().find_map(|d| d.positions.get(1).zip(d.positions.first()).map(|(b, a)| b - a));
        if let Some(bin) = bin {
            let fractions = gc::fractions(reference, depths, bin, span)?;
            gc::correct(depths, &fractions, &autosomes);
        }
    }
    Ok(out)
}

/// Masks (sets to NaN) windows whose depth is below `min_ratio` of their all-MAPQ depth.
/// Returns how many it masked.
fn mask_ambiguous(depths: &mut [ContigDepth], m: &Mappability) -> Result<usize, Error> {
    let all = read_regions(&m.all_mapq_coverage)?;
    let mut masked = 0;
    for d in depths.iter_mut() {
        let Some(a) = all.iter().find(|a| a.name == d.name) else { continue };
        if a.positions != d.positions {
            return Err(Error::Input(format!(
                "{}: the windows of {} differ from the coverage's",
                m.all_mapq_coverage.display(),
                d.name
            )));
        }
        for (v, &total) in d.depth.iter_mut().zip(&a.depth) {
            if total > 0.0 && *v < m.min_ratio * total {
                *v = f64::NAN;
                masked += 1;
            }
        }
    }
    Ok(masked)
}

/// Sets each call's copy number, genotype and depth (`CN`, `GT`, `CD`) as the extensions
/// say: from its final windows ([`CopyNumbers::FromCall`]) and in copies of a haploid
/// region's chromosome (one copy normal). Calls whose copy number then is the expected one
/// are dropped. Without either change the calls are untouched.
pub fn finish_copy_numbers(
    arena: &mut Arena,
    calls: &mut [(String, Vec<CandidateId>)],
    ext: &Extensions,
    prepared: &Prepared,
    ploidy: i64,
) {
    for &id in calls.iter().flat_map(|(_, ids)| ids) {
        let c = arena.get_mut(id);
        let copies = prepared.expected_copies(&c.chromosome, c.start + (c.end - c.start) / 2, ploidy);
        if ext.copy_numbers == CopyNumbers::FromCall {
            let median = sum::nanmedian(&c.cov);
            if median.is_finite() {
                c.median_cov_norm = median;
                c.median_raw_cov = median;
                c.cn_status = np_round(median, 0) as i64;
            }
        } else if copies == ploidy {
            continue;
        }
        if copies == ploidy {
            c.gt = gt(c.cn_status);
        } else {
            // Depth was scaled by ploidy / copies; scale the copy number back.
            c.median_cov_norm *= copies as f64 / ploidy as f64;
            c.median_raw_cov = c.median_cov_norm;
            c.cn_status = np_round(c.median_cov_norm, 0) as i64;
            c.gt = if c.cn_status == copies { "0" } else { "1" };
        }
    }
    for (_, ids) in calls.iter_mut() {
        ids.retain(|&id| {
            let c = arena.get(id);
            let copies = prepared.expected_copies(&c.chromosome, c.start + (c.end - c.start) / 2, ploidy);
            let recomputed = ext.copy_numbers == CopyNumbers::FromCall || copies != ploidy;
            !recomputed || c.cn_status != copies
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::{Candidate, CnvType};

    fn male() -> Prepared {
        Prepared { haploid: vec![("chrX".into(), vec![(0, 3000)])], ..Prepared::default() }
    }

    /// A male at 30x, chrX outside its PARs at 15x: after `prepare` every window reads 30x, in
    /// GRCh38, GRCh37 and CHM13 alike, with the PARs' first and last bases diploid (LRP-018).
    #[test]
    fn male_pars_of_every_build_are_diploid() {
        use crate::genome::Contig;
        use lrpipe_reads::par::Build;
        for build in [Build::Grch38, Build::Grch37, Build::Chm13] {
            let [par1, par2] = build.x();
            let x_len = u64::try_from(par2.1).expect("positive") + if build == Build::Chm13 { 0 } else { 10_000 };
            let contig = |name: &str, len: u64| Contig { name: name.into(), len, ploidy_multiplier: 2 };
            let genome = Genome { contigs: vec![contig("chr1", 1_000_000), contig("chrX", x_len)], ploidy: 2 };
            let par = |p: i64| [par1, par2].iter().any(|&(s, e)| p >= s && p < e);
            let sites = [par1.0, par1.1 - 1, par1.1, par2.0 - 1, par2.0, par2.1 - 1, 10_000_000];
            let mut depths = vec![
                ContigDepth { name: "chr1".into(), positions: vec![0, 1000], depth: vec![30.0, 30.0] },
                ContigDepth {
                    name: "chrX".into(),
                    positions: sites.to_vec(),
                    depth: sites.iter().map(|&p| if par(p) { 30.0 } else { 15.0 }).collect(),
                },
            ];
            let ext = Extensions { sex_chromosomes: SexChromosomes::Given(Sex::Male), ..Extensions::default() };
            let prepared = prepare(&mut depths, &ext, &genome, Path::new("unused.fa")).expect("prepared");
            for (&p, &d) in sites.iter().zip(&depths[1].depth) {
                assert_eq!(prepared.expected_copies("chrX", p, 2), if par(p) { 2 } else { 1 }, "{build:?} chrX:{p}");
                assert_eq!(d, 30.0, "{build:?} chrX:{p} normalized to the autosomes' depth");
            }
        }
    }

    #[test]
    fn haploid_outside_pars_only() {
        let p = male();
        assert_eq!(p.expected_copies("chrX", 1000, 2), 2);
        assert_eq!(p.expected_copies("chrX", 5000, 2), 1);
        assert_eq!(p.expected_copies("chr1", 5000, 2), 2);
    }

    #[test]
    fn copy_numbers_from_final_windows_and_haploid_scale() {
        let mut arena = Arena::default();
        // A chrX duplication at doubled depth 4.0 (2 copies of one), and an autosome deletion
        // whose raw run read 0.4 but whose final windows read 1.1.
        let dup = arena.push(Candidate::new("chrX", CnvType::Dup, 1000, vec![10_000, 11_000], vec![4.0, 4.0]));
        let del = arena.push(Candidate::new("chr1", CnvType::Del, 1000, vec![0, 1000], vec![0.4, 0.4]));
        arena.get_mut(del).cov = vec![1.1, 1.1, f64::NAN];
        let mut calls = vec![("chrX".to_owned(), vec![dup]), ("chr1".to_owned(), vec![del])];
        let ext = Extensions { copy_numbers: CopyNumbers::FromCall, ..Extensions::default() };
        finish_copy_numbers(&mut arena, &mut calls, &ext, &male(), 2);
        assert_eq!((arena.get(dup).cn_status, arena.get(dup).gt), (2, "1"));
        assert_eq!((arena.get(del).cn_status, arena.get(del).gt), (1, "0/1"));
    }

    #[test]
    fn spectre_default_leaves_calls_alone() {
        let mut arena = Arena::default();
        let del = arena.push(Candidate::new("chr1", CnvType::Del, 1000, vec![0, 1000], vec![0.4, 0.4]));
        arena.get_mut(del).cov = vec![1.1, 1.1];
        let before = arena.get(del).clone();
        let mut calls = vec![("chr1".to_owned(), vec![del])];
        finish_copy_numbers(&mut arena, &mut calls, &Extensions::default(), &Prepared::default(), 2);
        assert_eq!(calls[0].1, vec![del]);
        assert_eq!(*arena.get(del), before);
    }
}

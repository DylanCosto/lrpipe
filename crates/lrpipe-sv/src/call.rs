//! SV calls from clusters, QC, genotyping and phasing (`sv.py`, `postprocessing.py`,
//! `genotyping.py`).

use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;

use crate::cluster::{Cluster, Leadtab};
use crate::config::{Args, Config, MinSupport};
use crate::consensus;
use crate::fasta::NRuns;
use crate::pymath::{self, PyNum};
use crate::reads::{BndInfo, Lead, Source, SvType};

/// Reads above which a call shorter than `--minsvlen` (with a `~` prefix) is still reported.
const SVLEN_MIN_SUPPORT_OVERRIDE: usize = 10;
/// Reads from supplementary alignments needed before INLINE_SA can filter a call.
const INLINE_SA_MIN_READS: usize = 5;
/// Ratio below which coverage on one side of a long DEL/DUP counts as sloped.
const COV_SLOPE_MAX_RATIO: f64 = 0.7;
/// Ratio between upstream and downstream coverage below which a long DEL/DUP is filtered.
const COV_UP_DOWN_MIN_RATIO: f64 = 0.5;
/// Read counts above this are scaled down to it before computing genotype likelihoods.
const GENOTYPE_NORMALIZATION_TARGET: i64 = 250;
/// Maximum genotype quality (and z score).
const MAX_GQ: i64 = 60;
/// Phasing rescue: minimum share of a haplotype's reads that support the call.
const RESCUE_MIN_HAP_FRACTION: f64 = 0.75;
/// Phasing rescue: the call needs more supporting leads than this.
const RESCUE_MIN_LEADS: usize = 3;
/// Phasing rescue: minimum support as a fraction of `--dev-minreads-extra`.
const RESCUE_SUPPORT_FRACTION: f64 = 0.60;

/// An INFO value, formatted as Sniffles2's `format_info` would.
#[derive(Debug, Clone, PartialEq)]
pub enum Info {
    Int(i64),
    Float(f64),
    Str(String),
    /// A flag (`True` is written as the bare key).
    Bool(bool),
    /// A key whose value is `None` (omitted from the output).
    None,
}

impl From<PyNum> for Info {
    fn from(v: PyNum) -> Self {
        match v {
            PyNum::Int(i) => Info::Int(i),
            PyNum::Float(f) => Info::Float(f),
        }
    }
}

/// HP and PS of a genotype's phase (`None` where Sniffles2 has `None`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Phase {
    pub hp: Option<String>,
    pub ps: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    pub a: i64,
    pub b: i64,
    pub gq: i64,
    pub dr: i64,
    pub dv: i64,
    pub phase: Phase,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Alt {
    /// `<INS>`, `<DEL>`, ... for the call's type.
    Symbolic,
    Seq(Vec<u8>),
}

#[derive(Debug)]
pub struct SvCall {
    pub svtype: SvType,
    pub pos: i64,
    pub end: i64,
    pub id: String,
    pub alt: Alt,
    pub qual: i64,
    pub filter: &'static str,
    pub info: BTreeMap<&'static str, Info>,
    pub svlen: i64,
    pub svlens: Option<Vec<i64>>,
    pub genotype: Option<Genotype>,
    pub precise: bool,
    pub support: usize,
    pub rnames: Vec<Arc<[u8]>>,
    pub qc: bool,
    /// Whether the phasing rescue (`rescue_phasing`) looked at the call: a MOSAIC_VAF call it
    /// re-passes when its phase says so, as [`rescued_by_phase`] decides.
    pub rescuable: bool,
    pub nm: f64,
    pub fwd: usize,
    pub rev: usize,
    /// Coverage upstream, at start, center, end and downstream.
    pub cov_up: i64,
    pub cov_start: i64,
    pub cov_center: i64,
    pub cov_end: i64,
    pub cov_down: i64,
    pub bnd: Option<BndInfo>,
    pub cluster: Cluster,
}

impl SvCall {
    fn info_num(&self, k: &str) -> Option<f64> {
        match self.info.get(k)? {
            Info::Int(i) => Some(*i as f64),
            Info::Float(f) => Some(*f),
            _ => None,
        }
    }

    fn info_int(&self, k: &str) -> Option<i64> {
        match self.info.get(k)? {
            Info::Int(i) => Some(*i),
            _ => None,
        }
    }
}

/// Coverage of a contig, from the read spans (Sniffles2 keeps a uint16 array per base).
#[derive(Debug)]
pub struct Coverage<'a> {
    starts: Vec<i64>,
    ends: Vec<i64>,
    len: i64,
    /// `N` bases of the contig's reference, if given (coverage there reads as 0).
    n_mask: Option<&'a NRuns>,
}

impl<'a> Coverage<'a> {
    pub fn new(spans: &[(i64, i64)], len: i64, n_mask: Option<&'a NRuns>) -> Self {
        let mut starts: Vec<i64> = spans.iter().map(|s| s.0.min(len)).collect();
        let mut ends: Vec<i64> = spans.iter().map(|s| s.1.min(len)).collect();
        starts.sort_unstable();
        ends.sort_unstable();
        Coverage { starts, ends, len, n_mask }
    }

    /// `cv[i]` with numpy indexing: negative indexes count from the end, out of range is `None`.
    pub fn get(&self, i: i64) -> Option<i64> {
        let i = if i < 0 { i + self.len } else { i };
        if i < 0 || i >= self.len {
            return None;
        }
        if self.n_mask.is_some_and(|m| m.is_n(i)) {
            return Some(0);
        }
        let started = self.starts.partition_point(|&s| s <= i);
        let ended = self.ends.partition_point(|&e| e <= i);
        Some(((started - ended) as i64) % 65536)
    }

    /// `coverage.mean()` over the whole contig (N bases count as 0).
    pub fn mean(&self) -> f64 {
        // The sum is an integer below 2^53, so numpy's float sum is exact in any order.
        let total: i64 = self.ends.iter().sum::<i64>() - self.starts.iter().sum::<i64>();
        let masked: i64 = self
            .n_mask
            .iter()
            .flat_map(|m| m.iter())
            .map(|(s, e)| (s..e.min(self.len)).map(|p| self.get_unmasked(p)).sum::<i64>())
            .sum();
        (total - masked) as f64 / self.len as f64
    }

    /// The coverage summed over consecutive windows of `width` bases from position 0, for every
    /// window that starts inside the contig (as Sniffles2's zero-padded numpy reshape).
    pub fn window_sums(&self, width: i64) -> Vec<i64> {
        let n = (self.len + width - 1) / width;
        let mut below = CoveredBelow::new(self);
        let mut prev = 0;
        let mut sums: Vec<i64> = (1..=n)
            .map(|w| {
                let next = below.at((w * width).min(self.len));
                let sum = next - prev;
                prev = next;
                sum
            })
            .collect();
        let mut below = CoveredBelow::new(self);
        for (s, e) in self.n_mask.iter().flat_map(|m| m.iter()) {
            let (mut p, e) = (s.max(0), e.min(self.len));
            while p < e {
                let w = p / width;
                let to = e.min((w + 1) * width);
                sums[w as usize] -= below.at(to) - below.at(p);
                p = to;
            }
        }
        sums
    }

    fn get_unmasked(&self, i: i64) -> i64 {
        let started = self.starts.partition_point(|&s| s <= i);
        let ended = self.ends.partition_point(|&e| e <= i);
        (started - ended) as i64
    }
}

/// Per-contig calling state shared by all calls.
#[derive(Debug)]
pub struct ContigCtx<'a> {
    pub cfg: &'a Config,
    pub tabs: &'a [Leadtab; 5],
}

/// `sv.call_from`: an SV call from a cluster, or `None` if it is too short.
///
/// The call's `id` is left empty; the caller numbers calls in Sniffles2's order.
pub fn call_from(mut cluster: Cluster, tab: &Leadtab, cfg: &Config) -> Option<SvCall> {
    let svtype = cluster.svtype;
    let leads: Vec<&Lead> = cluster.leads.iter().map(|&i| &tab.leads[i]).collect();
    let svlens: Vec<i64> = leads.iter().map(|l| l.svlen.expect("clustered leads have a length")).collect();
    let svlen = pymath::center(&svlens);
    if svtype != SvType::Bnd && svlen.abs() < cfg.minsvlen_screen {
        return None;
    }
    let mut support_set: HashSet<&[u8]> = HashSet::new();
    let mut rnames: Vec<Arc<[u8]>> = Vec::new();
    for l in &leads {
        if support_set.insert(&l.qname) {
            rnames.push(l.qname.clone());
        }
    }
    let mut support_long = 0;
    if svtype == SvType::Ins && svlen >= cfg.args.long_ins_length {
        let long = cluster.leads_long.as_deref().map_or(&[][..], Vec::as_slice);
        let mut long_set: HashSet<&[u8]> = HashSet::new();
        for &i in long {
            let q = &tab.leads[i].qname;
            long_set.insert(q);
            if support_set.insert(q) {
                rnames.push(q.clone());
            }
        }
        support_long = long_set.len();
    }
    let support = support_set.len();
    let starts: Vec<i64> = leads.iter().map(|l| l.ref_start).collect();
    let ref_start = pymath::center(&starts);
    let stdev_pos = pymath::util_stdev(&pymath::trim(&starts));
    let (stdev_len, precise) = if svtype != SvType::Bnd {
        let sl = pymath::util_stdev(&pymath::trim(&svlens));
        (Some(sl), stdev_pos.as_f64() + sl.as_f64() < cfg.precise)
    } else {
        (None, stdev_pos.as_f64() < cfg.precise)
    };
    let (svstart, svend) = match svtype {
        SvType::Ins => (ref_start, ref_start),
        SvType::Del => (ref_start + svlen, ref_start),
        _ => (ref_start, ref_start + svlen.abs()),
    };
    let mapqs: Vec<i64> = leads.iter().map(|l| i64::from(l.mapq)).collect();
    let qual = pymath::py_int(pymath::mean_int(&mapqs));
    let fwd = leads.iter().filter(|l| !l.reverse).count();
    let rev = leads.len() - fwd;
    let nm = if cfg.qc_nm_measure { leads.iter().fold(0.0, |s, l| s + l.nm) / leads.len() as f64 } else { -1.0 };
    let support_sa = leads.iter().filter(|l| l.source != Source::Inline).count();
    drop(leads);

    let mut info = BTreeMap::new();
    let mut alt = Alt::Symbolic;
    let mut bnd = None;
    let mut support = support;
    match svtype {
        SvType::Bnd => {
            let r = resolve_bnd(&mut cluster, tab);
            alt = Alt::Seq(r.alt.into_bytes());
            support = r.support;
            info.insert("CHR2", Info::Str(r.info.mate_contig.to_string()));
            bnd = Some(r.info);
        }
        SvType::Ins => {
            info.insert("SUPPORT_LONG", Info::Int(support_long as i64));
        }
        SvType::Del => {
            info.insert("SUPPORT_SA", Info::Int(support_sa as i64));
        }
        _ => {}
    }
    info.insert("STDEV_POS", stdev_pos.into());
    if let Some(sl) = stdev_len {
        info.insert("STDEV_LEN", sl.into());
    }
    let call = SvCall {
        svtype,
        pos: svstart,
        end: svend,
        id: String::new(),
        alt,
        qual,
        filter: "PASS",
        info,
        svlen,
        svlens: cfg.args.dev_emit_sv_lengths.then(|| svlens.clone()),
        genotype: None,
        precise,
        support,
        rnames,
        qc: true,
        rescuable: false,
        nm,
        fwd,
        rev,
        cov_up: 0,
        cov_start: 0,
        cov_center: 0,
        cov_end: 0,
        cov_down: 0,
        bnd,
        cluster,
    };
    Some(call)
}

/// The mate and ALT of a BND call (`sv.resolve_bnd`).
struct ResolvedBnd {
    alt: String,
    support: usize,
    info: BndInfo,
}

/// Keeps the cluster's leads to the most common mate contig and describes the breakend.
fn resolve_bnd(cluster: &mut Cluster, tab: &Leadtab) -> ResolvedBnd {
    let bnd = |i: usize| tab.leads[i].bnd.as_deref().expect("BND leads have mate info");
    let mate_contig = pymath::most_common_top(cluster.leads.iter().map(|&i| bnd(i).mate_contig.clone()));
    let selected: Vec<usize> = cluster.leads.iter().copied().filter(|&i| bnd(i).mate_contig == mate_contig).collect();
    let starts: Vec<i64> = selected.iter().map(|&i| bnd(i).mate_ref_start).collect();
    let mate_ref_start = pymath::center(&starts);
    let is_first = pymath::most_common_top(selected.iter().map(|&i| bnd(i).is_first));
    let is_reverse = pymath::most_common_top(selected.iter().map(|&i| bnd(i).is_reverse));
    let bracket = if is_reverse { "]" } else { "[" };
    let alt = format!(
        "{}{bracket}{mate_contig}:{mate_ref_start}{bracket}{}",
        if is_first { "N" } else { "" },
        if is_first { "" } else { "N" }
    );
    let support = selected.iter().map(|&i| &tab.leads[i].qname).collect::<HashSet<_>>().len();
    cluster.leads = selected;
    ResolvedBnd { alt, support, info: BndInfo { mate_contig, mate_ref_start, is_first, is_reverse } }
}

/// The summed coverage of positions [0, x) for non-decreasing x: reads started before x
/// times their covered length up to x, minus the same for reads ended.
struct CoveredBelow<'c> {
    starts: &'c [i64],
    ends: &'c [i64],
    n_started: usize,
    n_ended: usize,
    sum_started: i64,
    sum_ended: i64,
}

impl<'c> CoveredBelow<'c> {
    fn new(cv: &'c Coverage<'_>) -> Self {
        CoveredBelow { starts: &cv.starts, ends: &cv.ends, n_started: 0, n_ended: 0, sum_started: 0, sum_ended: 0 }
    }

    fn at(&mut self, x: i64) -> i64 {
        while self.n_started < self.starts.len() && self.starts[self.n_started] < x {
            self.sum_started += self.starts[self.n_started];
            self.n_started += 1;
        }
        while self.n_ended < self.ends.len() && self.ends[self.n_ended] < x {
            self.sum_ended += self.ends[self.n_ended];
            self.n_ended += 1;
        }
        (self.n_started as i64 * x - self.sum_started) - (self.n_ended as i64 * x - self.sum_ended)
    }
}

/// `postprocessing.coverage`: annotates calls with coverage around them.
pub fn annotate_coverage(calls: &mut [SvCall], cv: &Coverage<'_>, cfg: &Config) {
    let bs = cfg.coverage_binsize;
    let ud = bs * cfg.args.coverage_updown_bins;
    for c in calls {
        let mut start = c.pos;
        if c.svtype == SvType::Bnd && c.bnd.as_ref().is_some_and(|b| b.is_first) {
            start -= 1;
        }
        let end = match c.svtype {
            SvType::Ins | SvType::Bnd => start + 1,
            _ => c.pos + c.svlen.abs(),
        };
        if matches!(c.svtype, SvType::Ins | SvType::Bnd) {
            if let Some(v) = cv.get(start - bs) {
                c.cov_start = v;
            }
            if let Some(v) = cv.get(start) {
                c.cov_center = v;
            }
            if let Some(v) = cv.get(end + bs) {
                c.cov_end = v;
            }
        } else {
            if let Some(v) = cv.get(start) {
                c.cov_start = v;
            }
            if let Some(v) = cv.get((start + end) / 2) {
                c.cov_center = v;
            }
            if let Some(v) = cv.get(end - bs) {
                c.cov_end = v;
            }
        }
        if let Some(v) = cv.get(start - ud) {
            c.cov_up = v;
        }
        if let Some(v) = cv.get(end + ud) {
            c.cov_down = v;
        }
    }
}

fn rescale_support(c: &SvCall, cfg: &Config) -> i64 {
    if c.svtype != SvType::Ins || c.svlen < cfg.args.long_ins_length {
        c.support as i64
    } else {
        let scale = cfg.long_ins_rescale_mult * (c.svlen as f64 / cfg.args.long_ins_length as f64);
        pymath::py_round(c.support as f64 * (cfg.long_ins_rescale_base + scale))
    }
}

fn qc_sv_support(c: &mut SvCall, coverage_global: f64, cfg: &Config) -> bool {
    let ok = match cfg.minsupport {
        MinSupport::Auto => {
            let support = rescale_support(c, cfg);
            let mut covs: Vec<i64> = [c.cov_up, c.cov_down].into_iter().filter(|&v| v != 0).collect();
            if covs.is_empty() {
                covs = [c.cov_start, c.cov_center, c.cov_end].into_iter().filter(|&v| v != 0).collect();
            }
            let regional = if covs.is_empty() {
                coverage_global
            } else {
                let r = pymath::py_round(covs.iter().sum::<i64>() as f64 / covs.len() as f64);
                if r == 0 { coverage_global } else { r as f64 }
            };
            let w = cfg.minsupport_auto_regional_coverage_weight;
            let cov = regional * w + coverage_global * (1.0 - w);
            let min_support = pymath::py_round(cfg.minsupport_auto_base + cfg.minsupport_auto_mult * cov);
            support >= min_support
        }
        MinSupport::Const(m) => c.support >= m,
    };
    if !ok {
        c.filter = "SUPPORT_MIN";
    }
    ok
}

fn strands(c: &SvCall, tab: &Leadtab) -> usize {
    let mut fwd = false;
    let mut rev = false;
    for &i in &c.cluster.leads {
        if tab.leads[i].reverse {
            rev = true;
        } else {
            fwd = true;
        }
    }
    usize::from(fwd) + usize::from(rev)
}

/// `postprocessing.qc_sv`. Sets the filter and returns false on failure.
fn qc_sv(c: &mut SvCall, tab: &Leadtab, cfg: &Config) -> bool {
    macro_rules! fail {
        ($f:expr) => {{
            c.filter = $f;
            return false;
        }};
    }
    let svlen_abs = c.svlen.abs() as f64;
    if cfg.args.qc_stdev {
        let stdev_pos = c.info_num("STDEV_POS").expect("STDEV_POS is always set");
        if stdev_pos > cfg.args.qc_stdev_abs_max {
            fail!("STDEV_POS");
        }
        if c.svtype != SvType::Bnd && stdev_pos / svlen_abs > 2.0 {
            fail!("STDEV_POS");
        }
        if let Some(stdev_len) = c.info_num("STDEV_LEN")
            && stdev_len != 0.0
        {
            if c.svtype != SvType::Bnd && stdev_len / svlen_abs > 1.0 {
                fail!("STDEV_LEN");
            }
            if stdev_len > cfg.args.qc_stdev_abs_max {
                fail!("STDEV_LEN");
            }
        }
    }
    if c.svlen.abs() < cfg.minsvlen
        && c.svtype != SvType::Bnd
        && (c.support < SVLEN_MIN_SUPPORT_OVERRIDE || cfg.minsvlen_hard_cap)
    {
        fail!("SVLEN_MIN");
    }
    if c.svtype == SvType::Bnd && cfg.args.qc_bnd_filter_strand && strands(c, tab) < 2 {
        fail!("STRAND_BND");
    }
    let max_cov_diff = COV_SLOPE_MAX_RATIO;
    let ud_diff = COV_UP_DOWN_MIN_RATIO;
    let (up, center, down) = (c.cov_up, c.cov_center, c.cov_down);
    let ratio = |a: i64, b: i64| a as f64 / b as f64;
    if c.svtype == SvType::Del
        && cfg.args.long_del_length != -1
        && c.svlen.abs() >= cfg.args.long_del_length
        && !cfg.args.mosaic
        && c.svlen.abs() <= cfg.dev_longer_del
    {
        if center as f64 > (up + down) as f64 * (cfg.args.long_del_coverage / 2.0) {
            if up > center && center > down {
                if ratio(down, up) < max_cov_diff {
                    fail!("COV_CHANGE_DEL");
                }
            } else if up < center && center < down && ratio(up, down) < max_cov_diff {
                fail!("COV_CHANGE_DEL");
            }
        }
        if up > down {
            if ud_diff > ratio(down, up) || center > down {
                fail!("COV_CHANGE_DEL");
            }
        } else if up < down && (ud_diff > ratio(up, down) || up < center) {
            fail!("COV_CHANGE_DEL");
        }
    } else if c.svtype == SvType::Dup
        && cfg.args.long_dup_length != -1
        && c.svlen.abs() >= cfg.args.long_dup_length
        && !cfg.args.mosaic
        && c.svlen.abs() <= cfg.dev_longer_dup
    {
        if (center as f64) < (up + down) as f64 * (cfg.args.long_dup_coverage / 2.0) {
            if up > center && center > down {
                if ratio(down, up) < max_cov_diff {
                    fail!("COV_CHANGE_DUP");
                }
            } else if up < center && center < down && ratio(up, down) < max_cov_diff {
                fail!("COV_CHANGE_DUP");
            }
            // Unlike for DEL, Sniffles2 nests this check inside the coverage condition.
            if up > down {
                if ud_diff > ratio(down, up) || center < down {
                    fail!("COV_CHANGE_DUP");
                }
            } else if up < down && (ud_diff > ratio(up, down) || up > center) {
                fail!("COV_CHANGE_DUP");
            }
        }
    } else if c.svtype == SvType::Ins && (up < cfg.args.qc_coverage || down < cfg.args.qc_coverage) {
        fail!("COV_CHANGE_INS");
    }
    if matches!(c.svtype, SvType::Ins | SvType::Del) {
        let (sa_inline, sap_inline) = c.cluster.sa_counts;
        let sa_split = c.info_int("SUPPORT_SA");
        let no_split_sa = matches!(sa_split, None | Some(0));
        if sap_inline > cfg.args.dev_inline_sa_support_max && sa_inline > INLINE_SA_MIN_READS && no_split_sa {
            fail!("INLINE_SA");
        }
    }
    // The coverage variance sampler is never fed in Sniffles2 2.8.1, so COV_VAR never fails.
    c.info.insert("COVERAGE_VAR", Info::None);
    let frac = cfg.args.qc_coverage_max_change_frac;
    if frac != -1.0 {
        let f = |v: i64| if v != 0 { v as f64 } else { 1.0 };
        let (u, s, ce, e, d) = (f(c.cov_up), f(c.cov_start), f(c.cov_center), f(c.cov_end), f(c.cov_down));
        if (u - s).abs() / u.max(s) > frac {
            fail!("COV_CHANGE_FRAC_US");
        }
        if (s - ce).abs() / s.max(ce) > frac {
            fail!("COV_CHANGE_FRAC_SC");
        }
        if (ce - e).abs() / ce.max(e) > frac {
            fail!("COV_CHANGE_FRAC_CE");
        }
        if (e - d).abs() / e.max(d) > frac {
            fail!("COV_CHANGE_FRAC_ED");
        }
    }
    true
}

/// `postprocessing.qc_sv_post_annotate`.
fn qc_sv_post_annotate(c: &mut SvCall, tab: &Leadtab, cfg: &Config, coverage_total: f64, nm_threshold: f64) -> bool {
    macro_rules! fail {
        ($f:expr) => {{
            c.filter = $f;
            return false;
        }};
    }
    let args = &cfg.args;
    let af = c.info_num("VAF").unwrap_or(0.0);
    let is_mosaic = af <= args.mosaic_af_max;
    let gt_low = c.genotype.as_ref().is_none_or(|g| g.a + g.b < 2);
    if c.cov_center < args.qc_coverage && gt_low && c.svtype != SvType::Del && c.svlen.abs() > args.long_del_length {
        fail!("COV_MIN_GT");
    }
    if args.mosaic && !is_mosaic && !qc_sv_support(c, coverage_total, cfg) {
        return false;
    }
    let qc_nm = if args.mosaic && is_mosaic { true } else { args.qc_nm };
    let qc_nm_threshold = nm_threshold * args.qc_nm_mult;
    if qc_nm && c.nm > qc_nm_threshold && c.genotype.as_ref().is_none_or(|g| g.b == 0) {
        fail!("ALN_NM");
    }
    if !args.mosaic && is_mosaic {
        let skip_this_dup = c.svtype == SvType::Dup && af >= args.dev_min_dup_vaf;
        if !skip_this_dup {
            fail!("MOSAIC_VAF");
        }
    }
    if args.mosaic && is_mosaic {
        let stdev_pos = c.info_num("STDEV_POS");
        let stdev_len = c.info_num("STDEV_LEN");
        let mut min_support = args.mosaic_min_reads;
        if let (Some(sp), Some(sl)) = (stdev_pos, stdev_len) {
            // `svlen` here is the SVLEN info key, which is never set during calling (so 1).
            let filter_low_supp =
                (!c.precise || sl / c.svlen.abs() as f64 > 0.1 || sp > 5.0) && 1 <= args.max_svlen_mosaic;
            min_support = if matches!(c.svtype, SvType::Bnd | SvType::Inv) || filter_low_supp {
                args.mosaic_min_reads
            } else {
                args.mosaic_min_reads - 1
            };
        }
        if c.support < min_support {
            fail!("SUPPORT_MIN");
        }
        if c.svtype != SvType::Bnd && c.svlen.abs() > args.max_svlen_mosaic {
            fail!("SVLEN_MAX_MOSAIC");
        }
    }
    if c.svtype != SvType::Bnd {
        let is_long_ins = c.svtype == SvType::Ins && c.svlen >= args.long_ins_length;
        if !(args.mosaic && is_mosaic) && args.qc_strand {
            if !is_long_ins && strands(c, tab) < 2 {
                fail!("STRAND");
            }
        } else if args.mosaic
            && is_mosaic
            && args.mosaic_qc_strand
            && !is_long_ins
            && strands(c, tab) < 2
            && c.support >= args.mosaic_use_strand_thresholds
        {
            fail!("STRAND_MOSAIC");
        }
    }
    if args.mosaic
        && is_mosaic
        && matches!(c.svtype, SvType::Inv | SvType::Dup)
        && c.svlen < args.mosaic_qc_invdup_min_length
    {
        fail!("SVLEN_MIN_MOSAIC");
    }
    if c.cov_center < args.qc_coverage && !matches!(c.svtype, SvType::Del | SvType::Ins) {
        // Python: `(svtype == "INV" and svlen) > long_inv_length`
        let lhs = if c.svtype == SvType::Inv { c.svlen } else { 0 };
        if !(lhs > args.long_inv_length && !(args.mosaic && is_mosaic)) {
            fail!("COV_MIN");
        }
    }
    if args.mosaic {
        if is_mosaic && (af < args.mosaic_af_min || af > args.mosaic_af_max) {
            fail!("MOSAIC_VAF");
        } else if !is_mosaic && !args.mosaic_include_germline {
            fail!("NOT_MOSAIC_VAF");
        }
        if is_mosaic {
            let close = c
                .cluster
                .leads
                .iter()
                .map(|&i| &tab.leads[i])
                .filter(|l| {
                    l.qry_start <= args.dev_min_close_edge_dist
                        || (l.read_len - l.qry_start).abs() <= args.dev_min_close_edge_dist
                })
                .count();
            if close as f64 / c.support as f64 >= args.dev_min_read_close_edge_prop {
                fail!("MOSAIC_SV_CLOSE_EDGE");
            }
        }
    }
    true
}

/// `postprocessing.phase_sv`: sets PHASE and returns the (HP, PS) for the genotype.
fn phase_sv(c: &mut SvCall, tab: &Leadtab, cfg: &Config) -> Phase {
    // One entry per alignment (read_id), keeping the last lead's values.
    let mut ids: Vec<u64> = Vec::new();
    let mut vals: std::collections::HashMap<u64, ReadPhase> = std::collections::HashMap::new();
    for &i in &c.cluster.leads {
        let l = &tab.leads[i];
        if vals.insert(l.read_id, (l.hap.to_string(), l.phase_set.as_deref().map(str::to_string))).is_none() {
            ids.push(l.read_id);
        }
    }
    let reads: Vec<ReadPhase> = ids.iter().map(|id| vals[id].clone()).collect();
    let (info, phase) = phase_of(&reads, cfg.args.phase_conflict_threshold);
    c.info.insert("PHASE", Info::Str(info));
    phase
}

/// One supporting alignment's haplotype (its HP tag, `"0"` without one) and phase set (its PS
/// tag, `"NULL"` without one; `None` for a BND lead), as `phase_sv` counts them.
pub type ReadPhase = (String, Option<String>);

/// `phase_sv` from the haplotypes and phase sets of a call's supporting alignments, one each:
/// the PHASE INFO value and the genotype's phase.
pub fn phase_of(reads: &[ReadPhase], threshold: f64) -> (String, Phase) {
    let hp_list = pymath::most_common(reads.iter().map(|r| r.0.clone()));
    // PS values are strings, or None for BND leads (never mixed within a cluster).
    let ps_list = pymath::most_common(reads.iter().map(|r| r.1.clone()));
    let (hp_support, hp) = hp_list[0].clone();
    let (ps_support, ps) = ps_list[0].clone();
    let ps = ps.unwrap_or_else(|| "NULL".into());
    let other_hp: usize = hp_list.iter().filter(|(_, h)| *h != hp && h != "NULL").map(|(s, _)| s).sum();
    let other_ps: usize = ps_list
        .iter()
        .filter(|(_, p)| p.as_deref() != Some(ps.as_str()) && p.as_deref() != Some("NULL"))
        .map(|(s, _)| s)
        .sum();
    let t = threshold;
    let hp_filter = if (other_hp as f64 / (hp_support + other_hp) as f64) < t && hp != "NULL" && hp_support > 0 {
        "PASS"
    } else {
        "FAIL"
    };
    let ps_filter = if (other_ps as f64 / (ps_support + other_ps) as f64) < t && ps != "NULL" && ps_support > 0 {
        "PASS"
    } else {
        "FAIL"
    };
    let info = format!("{hp},{ps},{hp_support},{ps_support},{hp_filter},{ps_filter}");
    let phase = Phase {
        hp: ((hp == "1" || hp == "2") && hp_filter == "PASS").then_some(hp),
        ps: (ps_filter == "PASS").then_some(ps),
    };
    (info, phase)
}

/// `annotate_sv`'s phase of a homozygous call: a 1/1 genotype takes the supporting reads' main
/// haplotype and phase set from PHASE (`info`) unless that haplotype is `0`, which also passes
/// PHASE's haplotype filter. Returns the PHASE value then.
pub fn phase_homozygous(g: &mut Genotype, info: &str) -> Option<String> {
    if g.a != 1 || g.b != 1 {
        return None;
    }
    let f: Vec<&str> = info.split(',').collect();
    if f[0] == "0" {
        return None;
    }
    g.phase = Phase { hp: Some(f[0].to_string()), ps: Some(f[1].to_string()) };
    Some(format!("{},{},{},{},PASS,{}", f[0], f[1], f[2], f[3], f[5]))
}

/// Genotype likelihoods (`Genotyper.calculate`), with the per-type coverage and support rules.
fn genotype_sv(c: &mut SvCall, cfg: &Config, phase: Phase, ins_tab: &Leadtab) {
    let args = &cfg.args;
    // Support
    let support: i64 = match c.svtype {
        SvType::Ins => {
            c.info.insert("SUPPORT_UNSCALED", Info::Int(c.support as i64));
            rescale_support(c, cfg)
        }
        SvType::Dup => {
            let mut ins_support = 0;
            if c.support >= cfg.minsupport_const {
                let lower = c.svlen as f64 * 0.9;
                let upper = c.svlen as f64 * 1.1;
                let mut i = c.pos.div_euclid(100) * 100;
                while i < c.end + 100 {
                    if let Some(idx) = ins_tab.bin(i) {
                        for &k in idx {
                            if let Some(l) = ins_tab.leads[k].svlen
                                && lower < l as f64
                                && (l as f64) < upper
                            {
                                ins_support += 1;
                            }
                        }
                    }
                    i += 100;
                }
                if ins_support > 0 {
                    c.info.insert("SUPPORT_INLINE", Info::Int(ins_support));
                    c.support += ins_support as usize;
                }
            }
            c.support as i64
        }
        _ => c.support as i64,
    };
    // Coverage
    let from_list = |list: &[i64]| -> Option<i64> {
        let v: Vec<i64> = list.iter().copied().filter(|&x| x != 0).collect();
        (!v.is_empty()).then(|| pymath::py_round(v.iter().sum::<i64>() as f64 / v.len() as f64))
    };
    let coverage = match c.svtype {
        SvType::Ins => from_list(&[c.cov_center]),
        SvType::Dup => from_list(&[c.cov_up, c.cov_down]),
        SvType::Inv => from_list(&[c.cov_up, c.cov_down]).map(|v| v + pymath::py_round(support as f64 * 0.5)),
        SvType::Del => match c.info_int("SUPPORT_SA") {
            Some(sa) if sa != 0 => from_list(&[c.cov_start + sa, c.cov_center + sa, c.cov_end + sa]),
            _ => from_list(&[c.cov_start, c.cov_center, c.cov_end]),
        },
        SvType::Bnd => from_list(&[c.cov_start, c.cov_center, c.cov_end]),
    };
    let Some(read_coverage) = coverage else {
        c.filter = "GT_FAILED";
        c.qc = false;
        return;
    };
    let coverage = read_coverage.max(support);
    let af = support as f64 / coverage as f64;
    let (gt1, gq, z) = genotype_likelihoods(support, coverage, args);
    // --read-support-genotypes: only GT, GQ, DR and DV come from the reads; VAF and every filter
    // keep the scaled support.
    let (gt1, gq, support, coverage) = if args.read_support_genotypes && c.svtype == SvType::Ins {
        let reads = c.support as i64;
        let reads_coverage = read_coverage.max(reads);
        let (gt, gq, _) = genotype_likelihoods(reads, reads_coverage, args);
        (gt, gq, reads, reads_coverage)
    } else {
        (gt1, gq, support, coverage)
    };
    let update_this_dup = c.svtype == SvType::Dup && af >= args.dev_min_dup_vaf;
    let mut filter_by_z = z < cfg.genotype_min_z_score && !args.mosaic;
    if c.svtype == SvType::Ins && filter_by_z && c.svlen >= args.long_ins_length && args.detect_large_ins {
        filter_by_z = false;
    }
    if c.filter == "PASS" && filter_by_z {
        c.filter = if update_this_dup { "PASS" } else { "GT" };
        c.qc = !args.pass_only;
    }
    let (mut ga, mut gb) = gt1;
    if update_this_dup && gt1 == (0, 0) {
        (ga, gb) = (0, 1);
    }
    c.genotype = Some(Genotype { a: ga, b: gb, gq, dr: coverage - support, dv: support, phase });
    c.info.insert("VAF", Info::Float(af));
}

/// The likeliest genotype of `support` variant reads out of `coverage`, its GQ, and the Z score
/// against hom-ref (`Genotyper.calculate`).
fn genotype_likelihoods(support: i64, coverage: i64, args: &Args) -> ((i64, i64), i64, i64) {
    let (hom_ref_p, het_p, hom_var_p) =
        (args.genotype_error, 1.0 / args.genotype_ploidy as f64, 1.0 - args.genotype_error);
    let max_lead = support.max(coverage);
    let (ns, nc) = if max_lead > GENOTYPE_NORMALIZATION_TARGET {
        let norm = GENOTYPE_NORMALIZATION_TARGET as f64 / max_lead as f64;
        (pymath::py_round(support as f64 * norm), pymath::py_round(coverage as f64 * norm))
    } else {
        (support, coverage)
    };
    let binom = |p: f64| p.powf(ns as f64) * (1.0 - p).powf((nc - ns) as f64);
    let mut lik = [((0, 0), binom(hom_ref_p)), ((0, 1), binom(het_p)), ((1, 1), binom(hom_var_p))];
    lik.sort_by(|x, y| y.1.partial_cmp(&x.1).expect("likelihoods are not NaN"));
    let sum: f64 = lik.iter().fold(0.0, |s, x| s + x.1);
    let norm: Vec<((i64, i64), f64)> = lik.iter().map(|&(g, q)| (g, q / sum)).collect();
    let (gt1, q1) = norm[0];
    let q2 = norm[1].1;
    let qz = norm.iter().find(|x| x.0 == (0, 0)).expect("hom-ref present").1;
    let lr = |a: f64, b: f64| if a / b > 0.0 { pymath::log10_py(a / b) } else { 0.0 };
    let z = pymath::py_int(-10.0 * lr(qz, q1)).min(MAX_GQ);
    let gq = pymath::py_int(-10.0 * lr(q2, q1)).min(MAX_GQ);
    (gt1, gq, z)
}

/// `annotate_sv`: phase, genotype, and the INS sequence.
fn annotate_sv(c: &mut SvCall, ctx: &ContigCtx<'_>) {
    let cfg = ctx.cfg;
    let tab = &ctx.tabs[c.svtype as usize];
    let phase = phase_sv(c, tab, cfg);
    genotype_sv(c, cfg, phase, &ctx.tabs[SvType::Ins as usize]);
    if let Some(g) = &mut c.genotype
        && let Some(Info::Str(p)) = c.info.get("PHASE")
        && let Some(s) = phase_homozygous(g, p)
    {
        c.info.insert("PHASE", Info::Str(s));
    }
    if c.svtype == SvType::Ins && !cfg.args.symbolic {
        let mut merged: Vec<&Lead> =
            c.cluster.leads.iter().map(|&i| &tab.leads[i]).filter(|l| l.seq.is_some()).collect();
        if merged.is_empty() {
            return;
        }
        let diff = |l: &Lead| {
            (l.seq.as_ref().expect("filtered").len() as i64 - c.svlen).abs() as f64
                + (l.ref_start - c.pos).abs() as f64 * 1.5
        };
        let mut best = 0;
        let mut best_diff = diff(merged[0]);
        for (i, l) in merged.iter().enumerate().skip(1) {
            let d = diff(l);
            if d < best_diff {
                best = i;
                best_diff = d;
            }
        }
        let best_lead = merged.remove(best);
        let best_seq = best_lead.seq.as_deref().expect("filtered");
        c.alt = if merged.len() >= cfg.consensus_min_reads && !cfg.args.no_consensus {
            let skip =
                cfg.consensus_kmer_skip_base + (best_seq.len() as f64 * cfg.consensus_kmer_skip_seqlen_mult) as usize;
            let others: Vec<&[u8]> = merged.iter().map(|l| l.seq.as_deref().expect("filtered")).collect();
            Alt::Seq(consensus::novel_from_reads(best_seq, &others, cfg.consensus_kmer_len, skip, skip))
        } else {
            Alt::Seq(best_seq.to_vec())
        };
    }
}

/// Whether `Task.rescue_phasing` can re-pass a call that failed QC: a MOSAIC_VAF call with
/// more than three leads whose mean mismatch rate is within `--genotype-error`. Whether it does
/// depends on its phase alone ([`rescued_by_phase`]).
fn rescue_candidate(c: &SvCall, tab: &Leadtab, cfg: &Config) -> bool {
    let nm: Vec<f64> = c.cluster.leads.iter().map(|&i| tab.leads[i].nm).collect();
    let sv_nm = pymath::numpy_nanmean(&nm);
    !(sv_nm > cfg.args.genotype_error || nm.len() <= RESCUE_MIN_LEADS) && c.filter == "MOSAIC_VAF"
}

/// `Task.rescue_phasing`'s test of a [`rescue_candidate`]'s phase: its PHASE value (`info`)
/// passes the haplotype filter for haplotype 1 or 2, and at least three quarters of that
/// haplotype's reads at the call's first bin (`hap_counts`: of the SV type's leads there, then of
/// the reads over it, for HP 0, 1 and 2) support it.
pub fn rescued_by_phase(info: &str, hap_counts: [i64; 6]) -> bool {
    let f: Vec<&str> = info.split(',').collect();
    if f[4] != "PASS" {
        return false;
    }
    // int(hp) raises for non-numeric HP in Sniffles2 (caught, no rescue).
    let Ok(hp) = f[0].parse::<i64>() else { return false };
    let [_, sv1, sv2, _, hap1, hap2] = hap_counts;
    let (all, sv) = match hp {
        1 => (hap1, sv1),
        2 => (hap2, sv2),
        _ => return false,
    };
    all != 0 && sv as f64 / all as f64 >= RESCUE_MIN_HAP_FRACTION
}

/// What re-passing a call does (`Task.rescue_phasing`).
pub fn rescue(c: &mut SvCall) {
    c.filter = "PASS";
    if let Some(g) = &mut c.genotype {
        g.b = 1;
    }
    c.qc = true;
}

/// `Task.rescue_phasing`: re-passes MOSAIC_VAF calls that are well phased.
fn rescue_phasing(c: &mut SvCall, tab: &Leadtab, cfg: &Config) {
    c.rescuable = rescue_candidate(c, tab, cfg);
    if !c.rescuable {
        return;
    }
    let Some(Info::Str(p)) = c.info.get("PHASE") else { return };
    if rescued_by_phase(p, c.cluster.hap_counts) {
        rescue(c);
    }
}

/// `Task.finalize_candidates` for one call: QC, annotation, phasing rescue.
pub fn finalize(c: &mut SvCall, ctx: &ContigCtx<'_>, coverage_total: f64, nm_threshold: f64) {
    let cfg = ctx.cfg;
    let tab = &ctx.tabs[c.svtype as usize];
    c.qc = c.qc && qc_sv(c, tab, cfg);
    if !cfg.args.mosaic && c.qc {
        c.qc = qc_sv_support(c, coverage_total, cfg);
    }
    annotate_sv(c, ctx);
    c.qc = c.qc && qc_sv_post_annotate(c, tab, cfg, coverage_total, nm_threshold);
    let phasing_rescue = c.svtype != SvType::Bnd
        && c.svlen.abs() <= cfg.args.dev_maxsvlen_extra
        && c.support as f64 >= (cfg.args.dev_minreads_extra as f64 * RESCUE_SUPPORT_FRACTION).trunc();
    if !c.qc && phasing_rescue {
        rescue_phasing(c, tab, cfg);
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::*;

    #[test]
    fn long_insertion_support_scaled_by_length_reads_as_homozygous() {
        // HG002 ONT, chr1:38,829,924, a 2,748 bp insertion GIAB has as 0/1: 18 reads of 33.
        let args = Args::parse_from(["lrpipe-sv", "-i", "x.bam"]);
        let reads = 18;
        let scaled = pymath::py_round(reads as f64 * (1.66 + 0.33 * 2748.0 / 2500.0));
        assert_eq!(scaled, 36);
        assert_eq!(genotype_likelihoods(scaled, scaled.max(33), &args).0, (1, 1));
        assert_eq!(genotype_likelihoods(reads, 33, &args).0, (0, 1));
    }
}

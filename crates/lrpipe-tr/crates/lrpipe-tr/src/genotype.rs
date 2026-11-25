//! Alleles from the repeat sequences of a locus's reads.
//!
//! Reads are split into haplotypes by their `HP` tags when most carry one, otherwise into the
//! two length groups that best explain the lengths (when the split is clear: both groups
//! substantial and their medians apart). Each allele takes the most common sequence among its
//! reads of the median length, or the one closest to the others when all differ.

use crate::align::edit_distance;
use crate::consensus::consensus;

/// Copies of the locus expected.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ploidy {
    /// Haploid: chrX and chrY of an XY sample, chrM.
    One,
    /// Diploid.
    Two,
}

/// One read's observation: its repeat sequence and haplotag.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Observation {
    /// The read's bases between the locus boundaries.
    pub seq: Vec<u8>,
    /// The read's `HP` tag.
    pub haplotype: Option<u8>,
    /// The read's `PS` tag: haplotypes of different phase sets are unrelated, so only those of
    /// the locus's most common phase set count ([`call`]).
    pub phase_set: Option<i64>,
    /// A hash of the read's name, when known: a read gives a locus one observation
    /// ([`one_per_read`]).
    pub read: Option<u64>,
    /// Whether the observation is from a supplementary alignment.
    pub supplementary: bool,
}

/// `obs` with one observation per read, in their order: a read seen twice (its primary and a
/// supplementary alignment both reach the locus, as the two sides of an expansion the aligner
/// split do when supplementary alignments keep their bases, minimap2's `-Y`) keeps its primary
/// alignment's observation, else its first. Observations of unknown reads are all kept. Counted
/// twice, two molecules were enough for an expansion that needs three.
#[must_use]
pub fn one_per_read(obs: Vec<Observation>) -> Vec<Observation> {
    let mut kept: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    let mut out: Vec<Observation> = Vec::with_capacity(obs.len());
    for o in obs {
        let Some(read) = o.read else {
            out.push(o);
            continue;
        };
        match kept.get(&read) {
            None => {
                kept.insert(read, out.len());
                out.push(o);
            }
            Some(&at) if out[at].supplementary && !o.supplementary => out[at] = o,
            Some(_) => {}
        }
    }
    out
}

/// An allele and the reads behind it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Allele {
    /// The allele's repeat sequence.
    pub seq: Vec<u8>,
    /// Reads counted for it (homozygous calls split them between the copies).
    pub support: usize,
    /// Shortest supporting read sequence.
    pub min_len: usize,
    /// Longest supporting read sequence.
    pub max_len: usize,
}

/// The sequencing platform, which sets how much read repeat lengths vary.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Platform {
    /// PacBio HiFi.
    #[default]
    Hifi,
    /// Oxford Nanopore: noisier lengths, mostly in homopolymers.
    Ont,
}

/// Settings of genotyping.
#[derive(Clone, Debug, PartialEq)]
pub struct GenotypeOptions {
    /// Each allele of a heterozygous call needs this share of the reads (and two reads).
    pub min_allele_share: f64,
    /// Variance of a read's repeat length around its allele's, in squared bases, at length 0:
    /// the floor that keeps exact agreement within an allele from looking infinitely
    /// significant.
    pub noise_floor: f64,
    /// Standard deviation of a read's repeat length per base of allele length, added to the
    /// floor: long alleles' reads vary more.
    pub noise_per_base: f64,
    /// A split into two alleles must improve the fit by more than this many times `ln(reads)`
    /// (a BIC-style penalty for the extra allele).
    pub split_penalty: f64,
    /// Haplotags' two alleles stand only when this share of the tagged reads, net, prefer
    /// their own haplotype's allele (see [`call`]); otherwise their difference is read noise.
    pub min_haplotag_support: f64,
    /// Two exact read lengths that each hold this share of the reads (the second at least
    /// half the first's count) are two alleles even when the length test does not split them:
    /// single-base differences in homopolymers, where read noise is symmetric around one
    /// length but two alleles make two peaks.
    pub min_peak_share: f64,
    /// Fewer reads than this give no call.
    pub min_reads: usize,
    /// An expanded allele, at least twice as long as the other and 100 bases longer, stands on
    /// this many reads (see [`split_expansion`])...
    pub min_expansion_reads: usize,
    /// ... and this share of them.
    pub min_expansion_share: f64,
    /// How read lengths err in long homopolymers, when an allele's length is taken from it
    /// rather than from its consensus alone (ONT).
    pub length_model: Option<LengthModel>,
}

/// A read's repeat length error `d` (read minus allele) where the allele's longest homopolymer
/// has `h` bases: `P(d)` proportional to `exp(d / (c + f h))` for `d < 0` and
/// `exp(-d / (c + g h))` for `d >= 0`, mixed with a uniform `floor` over `|d| <= `
/// [`LengthModel::RANGE`]. Nanopore reads miss bases in homopolymers more often than they add
/// them, and more so the longer the run (fitted on HG002 chr1 against the GIAB tandem repeat
/// benchmark). An allele whose longest run has at least `min_run` bases, and at least
/// `min_run_share` of the allele's, takes the length most likely under it, and that run is
/// lengthened or shortened to match. Where the run is a small part of the allele (a long
/// expansion, a minisatellite with a short run in its unit), the rest of the allele's noise
/// outweighs it and the consensus stands.
#[derive(Clone, Debug, PartialEq)]
pub struct LengthModel {
    /// Width of both tails at a run of length 0.
    pub c: f64,
    /// Growth of the short tail (missed bases) per base of the run.
    pub f: f64,
    /// Growth of the long tail (added bases) per base of the run.
    pub g: f64,
    /// Share of reads whose error follows no pattern.
    pub floor: f64,
    /// Shortest run the model applies to.
    pub min_run: usize,
    /// Smallest share of the allele the run must make up.
    pub min_run_share: f64,
}

impl LengthModel {
    /// Errors beyond this many bases count as this many.
    pub const RANGE: i64 = 15;

    /// `ln P(d)` for `d` in `-RANGE..=RANGE`, at a run of `h` bases.
    fn log_probs(&self, h: usize) -> Vec<f64> {
        let (short, long) = (self.c + self.f * h as f64, self.c + self.g * h as f64);
        let w: Vec<f64> = (-Self::RANGE..=Self::RANGE)
            .map(|d| if d < 0 { (d as f64 / short).exp() } else { (-d as f64 / long).exp() })
            .collect();
        let z: f64 = w.iter().sum();
        let uniform = 1.0 / (2 * Self::RANGE + 1) as f64;
        w.iter().map(|x| ((1.0 - self.floor) * x / z + self.floor * uniform).ln()).collect()
    }

    /// The allele length most likely to give `lens`, at a longest run of `h` bases; ties go
    /// to the length nearest `prior`.
    fn best_length(&self, lens: &[usize], h: usize, prior: usize) -> usize {
        let lp = self.log_probs(h);
        let (lo, hi) = (lens.iter().min().copied().unwrap_or(prior), lens.iter().max().copied().unwrap_or(prior));
        let score = |len: usize| -> f64 {
            lens.iter()
                .map(|&l| lp[((l as i64 - len as i64).clamp(-Self::RANGE, Self::RANGE) + Self::RANGE) as usize])
                .sum()
        };
        (lo.saturating_sub(3)..=hi + 3)
            .map(|len| (score(len), std::cmp::Reverse(len.abs_diff(prior)), len))
            .max_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)))
            .map_or(prior, |(_, _, len)| len)
    }
}

/// The longest run of one base in `seq`: its start and length.
fn longest_run(seq: &[u8]) -> (usize, usize) {
    let (mut best, mut start) = ((0, 0), 0);
    for i in 0..=seq.len() {
        if i == seq.len() || seq[i] != seq[start] {
            if i - start > best.1 {
                best = (start, i - start);
            }
            start = i;
        }
    }
    best
}

/// `seq` with its longest run lengthened or shortened to make it `len` long, when the run
/// can absorb the change (it keeps at least one base).
fn resize_longest_run(seq: &[u8], len: usize) -> Option<Vec<u8>> {
    let (start, run) = longest_run(seq);
    let new_run = (run + len).checked_sub(seq.len()).filter(|&r| r >= 1)?;
    let base = seq[start];
    Some([&seq[..start], &vec![base; new_run][..], &seq[start + run..]].concat())
}

impl GenotypeOptions {
    /// The settings for a platform, tuned on HG002 chr1 against the GIAB tandem repeat
    /// benchmark.
    pub fn for_platform(platform: Platform) -> Self {
        match platform {
            Platform::Hifi => GenotypeOptions {
                min_allele_share: 0.2,
                noise_floor: 0.5,
                noise_per_base: 0.006,
                split_penalty: 5.0,
                min_haplotag_support: 0.3,
                min_peak_share: 0.3,
                min_reads: 2,
                min_expansion_reads: 3,
                min_expansion_share: 0.05,
                length_model: None,
            },
            Platform::Ont => GenotypeOptions {
                min_allele_share: 0.2,
                noise_floor: 8.0,
                noise_per_base: 0.03,
                split_penalty: 5.0,
                min_haplotag_support: 0.3,
                min_peak_share: 0.5,
                min_reads: 2,
                min_expansion_reads: 3,
                min_expansion_share: 0.05,
                length_model: Some(LengthModel {
                    c: 0.5,
                    f: 0.04,
                    g: 0.02,
                    floor: 0.05,
                    min_run: 10,
                    min_run_share: 0.1,
                }),
            },
        }
    }
}

impl Default for GenotypeOptions {
    fn default() -> Self {
        GenotypeOptions::for_platform(Platform::Hifi)
    }
}

fn median_of_sorted(v: &[usize]) -> usize {
    v[(v.len() - 1) / 2]
}

fn median_len(obs: &[Observation], group: &[usize]) -> usize {
    let mut l: Vec<usize> = group.iter().map(|&i| obs[i].seq.len()).collect();
    l.sort_unstable();
    median_of_sorted(&l)
}

/// The reads an allele's sequence is chosen from (`core`), and all the reads counted for it.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Group {
    core: Vec<usize>,
    reads: Vec<usize>,
}

/// `n ln((SS1 + n v) / (SS2 + n v))`: how much better two groups fit the reads' lengths than
/// one, with per-read length variance `v` from the noise settings at the reads' mean length.
fn split_gain(n: usize, one: f64, two: f64, mean_len: f64, opts: &GenotypeOptions) -> f64 {
    let noise = opts.noise_floor + (opts.noise_per_base * mean_len).powi(2);
    let floor = n as f64 * noise;
    n as f64 * ((one + floor) / (two + floor)).ln()
}

/// Two alleles by read length: the split of the sorted lengths into two groups with the
/// least within-group squared deviation (the optimal two-group split in one dimension), kept
/// when it fits the lengths enough better than one group. With per-read length variance `v`
/// (from [`GenotypeOptions::noise_floor`] and [`GenotypeOptions::noise_per_base`]) added to
/// both, the split is kept when `n * ln((SS1 + n v) / (SS2 + n v))` exceeds
/// [`GenotypeOptions::split_penalty`] `* ln(n)`, both groups hold
/// [`GenotypeOptions::min_allele_share`] of the reads, and their lengths differ.
fn split_lengths(obs: &[Observation], opts: &GenotypeOptions) -> Option<[Vec<usize>; 2]> {
    let n = obs.len();
    let min_group = ((n as f64 * opts.min_allele_share).ceil() as usize).max(2);
    if n < 2 * min_group {
        return None;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| obs[i].seq.len());
    let lens: Vec<f64> = order.iter().map(|&i| obs[i].seq.len() as f64).collect();
    // Prefix sums of lengths and squares give any range's squared deviation in O(1).
    let mut sum = vec![0.0; n + 1];
    let mut sq = vec![0.0; n + 1];
    for (i, &l) in lens.iter().enumerate() {
        sum[i + 1] = sum[i] + l;
        sq[i + 1] = sq[i] + l * l;
    }
    let ss = |a: usize, b: usize| {
        let (s, q, m) = (sum[b] - sum[a], sq[b] - sq[a], (b - a) as f64);
        (q - s * s / m).max(0.0)
    };
    let (two, k) = (min_group..=n - min_group).map(|k| (ss(0, k) + ss(k, n), k)).min_by(|a, b| a.0.total_cmp(&b.0))?;
    if lens[k - 1] == lens[k] {
        return None;
    }
    let gain = split_gain(n, ss(0, n), two, sum[n] / n as f64, opts);
    (gain > opts.split_penalty * (n as f64).ln()).then(|| [order[..k].to_vec(), order[k..].to_vec()])
}

/// Two alleles from the two most common exact lengths (see
/// [`GenotypeOptions::min_peak_share`]): the reads at each length.
fn split_peaks(obs: &[Observation], opts: &GenotypeOptions) -> Option<[Vec<usize>; 2]> {
    let mut counts: Vec<(usize, usize)> = Vec::new();
    for o in obs {
        let l = o.seq.len();
        match counts.iter_mut().find(|(len, _)| *len == l) {
            Some((_, c)) => *c += 1,
            None => counts.push((l, 1)),
        }
    }
    counts.sort_by_key(|&(l, c)| (std::cmp::Reverse(c), l));
    let (&(l1, c1), &(l2, c2)) = (counts.first()?, counts.get(1)?);
    let needed = ((obs.len() as f64 * opts.min_peak_share).ceil() as usize).max(3);
    if c2 < needed || c2 * 2 < c1 {
        return None;
    }
    let at = |l: usize| (0..obs.len()).filter(|&i| obs[i].seq.len() == l).collect::<Vec<_>>();
    Some([at(l1), at(l2)])
}

/// Two alleles when some reads are much longer than the rest, however few (at least
/// [`GenotypeOptions::min_expansion_reads`] and [`GenotypeOptions::min_expansion_share`]): an
/// expansion, whose long molecules are sampled less than the other allele's (in targeted
/// sequencing and whole genomes alike) and whose lengths vary from cell to cell. The reads split
/// at the widest gap between consecutive lengths, which must at least double the length below
/// it and add 100 bases.
fn split_expansion(obs: &[Observation], opts: &GenotypeOptions) -> Option<[Vec<usize>; 2]> {
    let n = obs.len();
    let min = opts.min_expansion_reads.max((n as f64 * opts.min_expansion_share).ceil() as usize);
    if n < 2 * min.max(1) {
        return None;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| obs[i].seq.len());
    let len = |k: usize| obs[order[k]].seq.len();
    let k = (min..=n - min).max_by_key(|&k| (len(k) - len(k - 1), std::cmp::Reverse(k)))?;
    let (below, above) = (len(k - 1), len(k));
    (above >= 2 * below && above >= below + 100).then(|| [order[..k].to_vec(), order[k..].to_vec()])
}

/// Where a locus's two groups of reads came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Split {
    Haplotags,
    Lengths,
}

/// What [`groups`] may split a locus's reads by.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum By {
    /// Their haplotags when usable, else their lengths.
    Haplotags,
    /// An expansion's much longer reads alone ([`split_expansion`]): what is left when the
    /// haplotags turned out not to tell the alleles apart.
    Expansion,
}

/// Splits observations into the reads of each copy: one group, or two (equal when
/// homozygous).
fn groups(obs: &[Observation], ploidy: Ploidy, opts: &GenotypeOptions, by: By) -> (Vec<Group>, Split) {
    let all: Vec<usize> = (0..obs.len()).collect();
    let n = obs.len();
    // Haplotags, when most reads have one and both haplotypes have reads. Those of the most
    // common phase set only (the lowest of equally common ones): H1 of one phase block is as
    // likely the other's H2, so where two blocks meet their reads' tags cannot be put together.
    let mut sets: Vec<(Option<i64>, usize)> = Vec::new();
    for o in obs.iter().filter(|o| o.haplotype.is_some()) {
        match sets.iter_mut().find(|(set, _)| *set == o.phase_set) {
            Some((_, count)) => *count += 1,
            None => sets.push((o.phase_set, 1)),
        }
    }
    let set = sets.iter().max_by_key(|&&(set, count)| (count, std::cmp::Reverse(set))).map(|&(set, _)| set);
    let tagged: Vec<(usize, u8)> = all
        .iter()
        .filter(|&&i| Some(obs[i].phase_set) == set)
        .filter_map(|&i| obs[i].haplotype.map(|h| (i, h)))
        .collect();
    let by_tag = |h: u8| tagged.iter().filter(|(_, t)| *t == h).map(|(i, _)| *i).collect::<Vec<_>>();
    let (h1, h2) = (by_tag(1), by_tag(2));
    let tags_usable = by == By::Haplotags && tagged.len() * 2 >= n && h1.len() >= 2 && h2.len() >= 2;
    let (cores, split) = match ploidy {
        Ploidy::Two if tags_usable => (vec![h1, h2], Split::Haplotags),
        Ploidy::Two if by == By::Expansion => {
            (split_expansion(obs, opts).map_or_else(|| vec![all.clone()], Vec::from), Split::Lengths)
        }
        Ploidy::Two => (
            split_expansion(obs, opts)
                .or_else(|| split_lengths(obs, opts))
                .or_else(|| split_peaks(obs, opts))
                .map_or_else(|| vec![all.clone()], Vec::from),
            Split::Lengths,
        ),
        Ploidy::One => (vec![all.clone()], Split::Lengths),
    };
    // Every read counts for the allele whose core length is nearest.
    let medians: Vec<usize> = cores.iter().map(|c| median_len(obs, c)).collect();
    let mut reads: Vec<Vec<usize>> = vec![Vec::new(); cores.len()];
    for i in all {
        let in_core = cores.iter().position(|c| c.contains(&i));
        let k = in_core.unwrap_or_else(|| {
            let l = obs[i].seq.len();
            (0..medians.len()).min_by_key(|&k| medians[k].abs_diff(l)).expect("one core at least")
        });
        reads[k].push(i);
    }
    let mut out: Vec<Group> = cores.into_iter().zip(reads).map(|(core, reads)| Group { core, reads }).collect();
    if ploidy == Ploidy::Two && out.len() == 1 {
        out.push(out[0].clone());
    }
    (out, split)
}

/// Whether the reads back their haplotags' two alleles `a` (of `h1`'s reads) and `b` (of
/// `h2`'s): each read prefers the one at the smaller edit distance; the tags stand when both
/// haplotypes' reads lean to their own allele and, net, at least
/// [`GenotypeOptions::min_haplotag_support`] of them do. Alleles that differ only by read
/// noise draw preferences unrelated to the tags.
fn tags_supported(obs: &[Observation], h1: &[usize], h2: &[usize], a: &[u8], b: &[u8], opts: &GenotypeOptions) -> bool {
    let lean = |reads: &[usize], own: &[u8], other: &[u8]| -> i64 {
        reads
            .iter()
            .map(|&i| {
                let seq = &obs[i].seq;
                match edit_distance(seq, own).cmp(&edit_distance(seq, other)) {
                    std::cmp::Ordering::Less => 1,
                    std::cmp::Ordering::Greater => -1,
                    std::cmp::Ordering::Equal => 0,
                }
            })
            .sum()
    };
    let (s1, s2) = (lean(h1, a, b), lean(h2, b, a));
    s1 > 0 && s2 > 0 && (s1 + s2) as f64 >= opts.min_haplotag_support * (h1.len() + h2.len()) as f64
}

/// Reads an allele's consensus is voted from, at most: more add time, not accuracy.
const MAX_CONSENSUS_READS: usize = 30;

/// The sequence an allele's reads agree on: their consensus, voted from a backbone that is the
/// most common sequence of the reads' median length; with a [`LengthModel`] and a long enough
/// run, resized to the length the reads' lengths make most likely.
fn allele_sequence(obs: &[Observation], group: &[usize], opts: &GenotypeOptions) -> Vec<u8> {
    let seq = voted_sequence(obs, group);
    let Some(model) = &opts.length_model else { return seq };
    let (_, run) = longest_run(&seq);
    if run < model.min_run || (run as f64) < model.min_run_share * seq.len() as f64 {
        return seq;
    }
    let lens: Vec<usize> = group.iter().map(|&i| obs[i].seq.len()).collect();
    let len = model.best_length(&lens, run, seq.len());
    if len == seq.len() {
        return seq;
    }
    resize_longest_run(&seq, len).unwrap_or(seq)
}

fn voted_sequence(obs: &[Observation], group: &[usize]) -> Vec<u8> {
    let target = median_len(obs, group);
    let mut counts: Vec<(&[u8], usize)> = Vec::new();
    for s in group.iter().map(|&i| obs[i].seq.as_slice()).filter(|s| s.len() == target) {
        match counts.iter_mut().find(|(c, _)| *c == s) {
            Some((_, n)) => *n += 1,
            None => counts.push((s, 1)),
        }
    }
    let (backbone, _) = counts.iter().copied().max_by_key(|&(_, n)| n).expect("the median length occurs");
    // The reads nearest the median length, in read order among equals: an expanded allele's
    // reads vary by hundreds of bases (somatic instability), and a consensus voted from reads
    // of very different lengths has none of them.
    let mut nearest: Vec<usize> = group.to_vec();
    nearest.sort_by_key(|&i| obs[i].seq.len().abs_diff(target));
    let seqs: Vec<&[u8]> = nearest.iter().take(MAX_CONSENSUS_READS).map(|&i| obs[i].seq.as_slice()).collect();
    consensus(&seqs, backbone)
}

/// Calls the alleles of a locus: none with too few reads, else one per copy (a homozygous
/// diploid call repeats the allele).
pub fn call(obs: &[Observation], ploidy: Ploidy, opts: &GenotypeOptions) -> Vec<Allele> {
    if obs.len() < opts.min_reads {
        return Vec::new();
    }
    let allele = |g: &Group| {
        let lens = g.reads.iter().map(|&i| obs[i].seq.len());
        Allele {
            seq: allele_sequence(obs, &g.core, opts),
            support: g.reads.len(),
            min_len: lens.clone().min().unwrap_or(0),
            max_len: lens.max().unwrap_or(0),
        }
    };
    let (mut gs, split) = groups(obs, ploidy, opts, By::Haplotags);
    if split == Split::Haplotags && gs[0] != gs[1] {
        let (a, b) = (allele(&gs[0]), allele(&gs[1]));
        if a.seq == b.seq || !tags_supported(obs, &gs[0].core, &gs[1].core, &a.seq, &b.seq, opts) {
            // The tags do not tell two alleles apart: one allele, unless some reads are an
            // expansion's, which tags that mix the alleles (a switch error, another run's tags)
            // hide in both haplotypes' medians. The other length splits are not tried: lengths
            // the tags call one allele's noise would split too often.
            gs = groups(obs, ploidy, opts, By::Expansion).0;
        } else {
            let mut alleles = vec![a, b];
            alleles.sort_by_key(|a| a.seq.len());
            return alleles;
        }
    }
    if let [a, b] = &gs[..]
        && a == b
    {
        // Homozygous: the reads are shared between the two copies.
        let mut first = allele(a);
        let mut second = first.clone();
        (first.support, second.support) = (a.reads.len().div_ceil(2), a.reads.len() / 2);
        return vec![first, second];
    }
    let mut alleles: Vec<Allele> = gs.iter().map(allele).collect();
    alleles.sort_by_key(|a| a.seq.len());
    alleles
}

/// Copies of each motif in `seq`, scanning left to right and taking at each position the
/// longest motif that matches exactly there.
pub fn motif_counts(seq: &[u8], motifs: &[Vec<u8>]) -> Vec<usize> {
    let mut counts = vec![0; motifs.len()];
    let mut order: Vec<usize> = (0..motifs.len()).collect();
    order.sort_by_key(|&m| std::cmp::Reverse(motifs[m].len()));
    let mut i = 0;
    while i < seq.len() {
        match order.iter().find(|&&m| !motifs[m].is_empty() && seq[i..].starts_with(&motifs[m])) {
            Some(&m) => {
                counts[m] += 1;
                i += motifs[m].len();
            }
            None => i += 1,
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(lens: &[(usize, Option<u8>)]) -> Vec<Observation> {
        lens.iter()
            .map(|&(l, h)| Observation { seq: b"CAG".repeat(l / 3), haplotype: h, ..Observation::default() })
            .collect()
    }

    #[test]
    fn splits_by_length_when_clear() {
        let o = obs(&[(30, None), (30, None), (30, None), (33, None), (60, None), (60, None), (60, None)]);
        let a = call(&o, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| (a.seq.len(), a.support)).collect::<Vec<_>>(), vec![(30, 4), (60, 3)]);
        let hom = obs(&[(30, None), (30, None), (33, None), (30, None)]);
        let a = call(&hom, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| a.seq.len()).collect::<Vec<_>>(), vec![30, 30]);
    }

    #[test]
    fn one_base_apart_needs_two_clear_groups() {
        let poly = |lens: &[usize]| -> Vec<Observation> {
            lens.iter().map(|&l| Observation { seq: vec![b'A'; l], ..Observation::default() }).collect()
        };
        let het = poly(&[[20; 12].as_slice(), &[21; 11], &[19; 2]].concat());
        let a = call(&het, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| a.seq.len()).collect::<Vec<_>>(), vec![20, 21]);
        let noise = poly(&[[20; 19].as_slice(), &[21; 3], &[19; 3]].concat());
        let a = call(&noise, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| a.seq.len()).collect::<Vec<_>>(), vec![20, 20]);
    }

    #[test]
    fn haplotags_need_their_reads_support() {
        // Tags split reads of one allele whose lengths vary by noise: no read prefers its own
        // haplotype's consensus, so the call is homozygous.
        let lens = [(20, 1), (20, 1), (21, 1), (21, 1), (20, 2), (21, 2), (20, 2), (21, 2), (19, 2)];
        let o: Vec<Observation> = lens
            .iter()
            .map(|&(l, h)| Observation { seq: vec![b'A'; l], haplotype: Some(h), ..Observation::default() })
            .collect();
        let a = call(&o, Ploidy::Two, &GenotypeOptions::for_platform(Platform::Ont));
        assert_eq!(a[0].seq, a[1].seq);
    }

    #[test]
    fn haplotags_decide_when_present() {
        let o = obs(&[(30, Some(1)), (30, Some(1)), (36, Some(2)), (36, Some(2)), (36, None)]);
        let a = call(&o, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| (a.seq.len(), a.support)).collect::<Vec<_>>(), vec![(30, 2), (36, 3)]);
    }

    #[test]
    fn counts_motifs() {
        assert_eq!(motif_counts(b"CAGCAGCAACAGTT", &[b"CAG".to_vec(), b"CAA".to_vec()]), vec![3, 1]);
    }

    #[test]
    fn ont_homopolymers_take_their_most_likely_length() {
        // Nanopore reads of a 20-base A run miss bases more often than they add them: the
        // median read is one short, the length model reads 20.
        let lens = [17, 18, 19, 19, 19, 20, 20, 20, 20, 21, 16, 18];
        let o: Vec<Observation> = lens
            .iter()
            .map(|&l| Observation { seq: [b"CG", &vec![b'A'; l][..], b"TC"].concat(), ..Observation::default() })
            .collect();
        let a = call(&o, Ploidy::Two, &GenotypeOptions::for_platform(Platform::Ont));
        assert_eq!(a.iter().map(|a| a.seq.len()).collect::<Vec<_>>(), vec![24, 24]);
        assert_eq!(a[0].seq, [b"CG", &[b'A'; 20][..], b"TC"].concat());
        // HiFi has no length model: its consensus stands.
        let a = call(&o, Ploidy::Two, &GenotypeOptions::for_platform(Platform::Hifi));
        assert_eq!(a[0].seq.len(), 23);
    }

    #[test]
    fn short_runs_and_long_alleles_keep_their_consensus() {
        let model = GenotypeOptions::for_platform(Platform::Ont).length_model.expect("ONT has one");
        assert_eq!(longest_run(b"CAGAAAAAG"), (3, 5));
        assert_eq!(resize_longest_run(b"CAGAAAAAG", 11).as_deref(), Some(&b"CAGAAAAAAAG"[..]));
        assert_eq!(resize_longest_run(b"CAGAAAAAG", 4), None);
        // A run under min_run, or under min_run_share of the allele, is left to the consensus.
        assert!(model.min_run > 5);
        let long = [b"CAG".repeat(40), vec![b'A'; 10]].concat();
        assert!((10.0) < model.min_run_share * long.len() as f64);
    }

    #[test]
    fn a_few_much_longer_reads_are_an_expansion() {
        // 20 reads of a normal allele and 4 of an expansion (long molecules are sampled less):
        // too few for the length split, enough for an expansion.
        let lens = [vec![39; 20], vec![2100; 4]].concat();
        let o: Vec<Observation> =
            lens.iter().map(|&l| Observation { seq: b"CAG".repeat(l / 3), ..Observation::default() }).collect();
        let a = call(&o, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| (a.seq.len(), a.support)).collect::<Vec<_>>(), vec![(39, 20), (2100, 4)]);
        // Two reads are not enough.
        let lens = [vec![39; 20], vec![2100; 2]].concat();
        let o: Vec<Observation> =
            lens.iter().map(|&l| Observation { seq: b"CAG".repeat(l / 3), ..Observation::default() }).collect();
        assert_eq!(call(&o, Ploidy::Two, &GenotypeOptions::default())[1].seq.len(), 39);
    }

    #[test]
    fn mixed_haplotags_do_not_hide_an_expansion() {
        // Tags that put both alleles' reads in both haplotypes: each haplotype's median read is
        // the normal allele's, so the tags give one allele; the expansion's reads still stand.
        let lens = [vec![39; 20], vec![2100; 4]].concat();
        let o: Vec<Observation> = lens
            .iter()
            .enumerate()
            .map(|(i, &l)| Observation {
                seq: b"CAG".repeat(l / 3),
                haplotype: Some(1 + (i % 2) as u8),
                phase_set: Some(1),
                ..Observation::default()
            })
            .collect();
        let a = call(&o, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| (a.seq.len(), a.support)).collect::<Vec<_>>(), vec![(39, 20), (2100, 4)]);
    }

    #[test]
    fn haplotags_of_two_phase_sets_are_not_put_together() {
        // Two phase blocks meet at the locus: in the second, H1 is the first's H2. Taken
        // together the tags mix the alleles; only the more common block's count, and the other
        // reads go to the allele of the nearest length.
        let tagged = |l: usize, h: u8, ps: i64| Observation {
            seq: b"CAG".repeat(l / 3),
            haplotype: Some(h),
            phase_set: Some(ps),
            ..Observation::default()
        };
        let mut o = Vec::new();
        for _ in 0..4 {
            o.extend([tagged(30, 1, 100), tagged(60, 2, 100)]);
        }
        for _ in 0..3 {
            o.extend([tagged(60, 1, 900), tagged(30, 2, 900)]);
        }
        let a = call(&o, Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| (a.seq.len(), a.support)).collect::<Vec<_>>(), vec![(30, 7), (60, 7)]);
    }

    #[test]
    fn a_read_counts_once() {
        let of = |read: u64, supplementary: bool, l: usize| Observation {
            seq: b"CAG".repeat(l),
            read: Some(read),
            supplementary,
            ..Observation::default()
        };
        // Read 7's supplementary alignment comes first; its primary one is kept, in its place.
        let o = one_per_read(vec![of(7, true, 700), of(8, false, 13), of(7, false, 701), of(8, true, 13)]);
        assert_eq!(o, vec![of(7, false, 701), of(8, false, 13)]);
        // Unknown reads are never one another.
        let unknown = vec![Observation { seq: b"CAG".to_vec(), ..Observation::default() }; 3];
        assert_eq!(one_per_read(unknown).len(), 3);
        // Two molecules of an expansion, each seen from both sides, are two reads: not the three
        // an expansion needs.
        let mut o: Vec<Observation> = (0..20).map(|i| of(i, false, 13)).collect();
        o.extend([of(100, false, 700), of(100, true, 700), of(101, false, 700), of(101, true, 700)]);
        let a = call(&one_per_read(o), Ploidy::Two, &GenotypeOptions::default());
        assert_eq!(a.iter().map(|a| a.seq.len()).collect::<Vec<_>>(), vec![39, 39]);
    }
}

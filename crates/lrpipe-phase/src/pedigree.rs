//! Port of the parts of `pedigree.py` needed for phasing.

use crate::core::genotype::Genotype;

pub fn centimorgen_to_phred(distance: f64) -> f64 {
    assert!(distance >= 0.0);
    if distance == 0.0 {
        panic!("Cannot convert genetic distance of zero to phred.");
    } else if distance < 1e-10 {
        -10.0 * (distance.log10() - 2.0)
    } else {
        let p = (1.0 - (-(2.0 * distance) / 100.0).exp()) / 2.0;
        -10.0 * p.log10()
    }
}

/// `UniformRecombinationCostComputer.compute`
pub fn uniform_recombination_map(recombrate: f64, positions: &[i32]) -> Vec<u32> {
    // Python: `[0] + [...]`, so even an empty position list yields [0]
    let mut result = Vec::with_capacity(positions.len().max(1));
    result.push(0);
    for i in 1..positions.len() {
        let d = (positions[i] - positions[i - 1]) as f64 * 1e-6 * recombrate;
        // Python's round() rounds half to even
        result.push(centimorgen_to_phred(d).round_ties_even() as u32);
    }
    result
}

/// `pedigree.mendelian_conflict`
pub fn mendelian_conflict(mother: &Genotype, father: &Genotype, child: &Genotype) -> bool {
    let m = mother.as_vector();
    let f = father.as_vector();
    let c = child.as_vector();
    if m.contains(&c[0]) && f.contains(&c[1]) {
        return false;
    }
    !(m.contains(&c[1]) && f.contains(&c[0]))
}

/// A trio from a PED file (`pedigree.Trio`); parents are None if unknown ("0").
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Trio {
    pub child: String,
    pub father: Option<String>,
    pub mother: Option<String>,
}

/// `PedReader`: parses a PED/FAM file (columns 2-4: individual, father, mother).
pub fn read_ped(path: &str) -> anyhow::Result<Vec<Trio>> {
    let text = std::fs::read_to_string(path).map_err(|e| anyhow::anyhow!("{e}: '{path}'"))?;
    let mut trios = Vec::new();
    // Python iterates lines including their newline; "\n" lines and "#" lines are skipped
    for line in text.split_inclusive('\n') {
        if line.starts_with('#') || line == "\n" {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            anyhow::bail!("Less than six fields found in PED/FAM file");
        }
        let parent = |s: &str| if s != "0" { Some(s.to_string()) } else { None };
        trios.push(Trio { child: fields[1].to_string(), father: parent(fields[2]), mother: parent(fields[3]) });
    }
    let mut counts: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for t in &trios {
        *counts.entry(t.child.as_str()).or_default() += 1;
    }
    // Counter.most_common()[0]: the first individual (in file order) with the maximum count
    if let Some(&max) = counts.values().max()
        && max > 1
    {
        let id = trios.iter().find(|t| counts[t.child.as_str()] == max).expect("an individual has the maximum count");
        anyhow::bail!("Individual '{}' occurs more than once in PED file", id.child);
    }
    Ok(trios)
}

/// `PedReader.samples()`. WhatsHap returns `list(set(...))`, whose order
/// depends on Python's hash seed; we return the samples sorted.
pub fn ped_samples(trios: &[Trio]) -> Vec<String> {
    let mut samples = std::collections::BTreeSet::new();
    for t in trios {
        if let (Some(f), Some(m)) = (&t.father, &t.mother) {
            samples.insert(f.clone());
            samples.insert(m.clone());
            samples.insert(t.child.clone());
        }
    }
    samples.into_iter().collect()
}

#[derive(Clone, Copy, Debug)]
pub struct RecombinationMapEntry {
    pub position: i64,
    pub cum_distance: f64,
}

/// `GeneticMapRecombinationCostComputer.load_genetic_map`
pub fn load_genetic_map(path: &str) -> anyhow::Result<Vec<RecombinationMapEntry>> {
    let text = std::fs::read_to_string(path).map_err(|e| anyhow::anyhow!("{e}: '{path}'"))?;
    let mut genetic_map: Vec<RecombinationMapEntry> = Vec::new();
    let mut warned_zero_distance = false;
    for (i, line) in text.lines().enumerate() {
        let line_number = i + 1;
        if line_number == 1 {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if fields.len() != 3 {
            anyhow::bail!(
                "Error at line {line_number} of genetic map file '{path}': Found {} fields instead of 3",
                fields.len()
            );
        }
        let position: i64 = fields[0].parse().map_err(|_| {
            anyhow::anyhow!("Error at line {line_number} of genetic map file '{path}': invalid literal for int() with base 10: '{}'", fields[0])
        })?;
        let cum_distance: f64 = fields[2].parse().map_err(|_| {
            anyhow::anyhow!(
                "Error at line {line_number} of genetic map file '{path}': could not convert string to float: '{}'",
                fields[2]
            )
        })?;
        genetic_map.push(RecombinationMapEntry { position, cum_distance });
        let n = genetic_map.len();
        if n >= 2 && !warned_zero_distance && genetic_map[n - 2].cum_distance == genetic_map[n - 1].cum_distance {
            log::warn!("Zero genetic distances encountered in {path}");
            warned_zero_distance = true;
        }
    }
    Ok(genetic_map)
}

fn interpolate(point: i64, start_pos: i64, end_pos: i64, start_value: f64, end_value: f64) -> f64 {
    assert!(start_pos <= point && point <= end_pos);
    if start_pos == point && point == end_pos {
        assert_eq!(start_value, end_value);
        return start_value;
    }
    start_value + ((point - start_pos) as f64 * (end_value - start_value) / (end_pos - start_pos) as f64)
}

const MINIMUM_GENETIC_DISTANCE: f64 = 1e-10;

/// `recombination_cost_map`
pub fn recombination_cost_map(genetic_map: &[RecombinationMapEntry], positions: &[i32]) -> Vec<u32> {
    assert!(!genetic_map.is_empty());
    let mut cumulative = Vec::with_capacity(positions.len());
    let mut i: Option<usize> = None;
    let mut j: Option<usize> = Some(0);
    for &position in positions {
        let position = position as i64;
        if i.is_none() && genetic_map[0].position <= position {
            i = Some(0);
        }
        while let Some(ii) = i {
            if ii + 1 < genetic_map.len() && genetic_map[ii + 1].position <= position {
                i = Some(ii + 1);
            } else {
                break;
            }
        }
        while let Some(jj) = j {
            if genetic_map[jj].position < position {
                j = if jj + 1 < genetic_map.len() { Some(jj + 1) } else { None };
            } else {
                break;
            }
        }
        let d = match (i, j) {
            (None, None) => panic!("AssertionError: position {position} outside of genetic map"),
            (None, Some(j)) => interpolate(position, 0, genetic_map[j].position, 0.0, genetic_map[j].cum_distance),
            (_, None) => {
                let last = genetic_map[genetic_map.len() - 1];
                let avg_rate = last.cum_distance / last.position as f64;
                last.cum_distance + (position - last.position) as f64 * avg_rate
            }
            (Some(i), Some(j)) => interpolate(
                position,
                genetic_map[i].position,
                genetic_map[j].position,
                genetic_map[i].cum_distance,
                genetic_map[j].cum_distance,
            ),
        };
        cumulative.push(d);
    }
    let mut result = Vec::with_capacity(positions.len().max(1));
    result.push(0);
    for k in 1..cumulative.len() {
        let d = (cumulative[k] - cumulative[k - 1]).max(MINIMUM_GENETIC_DISTANCE);
        result.push(centimorgen_to_phred(d).round_ties_even() as u32);
    }
    result
}

/// `RecombinationEvent` (field order defines the sort order, as in Python)
#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct RecombinationEvent {
    pub position1: i32,
    pub position2: i32,
    pub transmitted_hap_father1: u32,
    pub transmitted_hap_father2: u32,
    pub transmitted_hap_mother1: u32,
    pub transmitted_hap_mother2: u32,
    pub recombination_cost: u32,
}

/// `find_recombination`
pub fn find_recombination(
    transmission_vector: &[u32],
    components: &std::collections::HashMap<i32, i32>,
    positions: &[i32],
    recombcost: &[u32],
) -> anyhow::Result<Vec<RecombinationEvent>> {
    if !(transmission_vector.len() == positions.len() && positions.len() == recombcost.len()) {
        anyhow::bail!("AssertionError: lengths of transmission vector, positions and recombination costs differ");
    }
    let position_to_index: std::collections::HashMap<i32, usize> =
        positions.iter().enumerate().map(|(i, &p)| (p, i)).collect();
    let mut blocks: std::collections::HashMap<i32, Vec<i32>> = std::collections::HashMap::new();
    for (&position, &block_id) in components {
        blocks.entry(block_id).or_default().push(position);
    }
    let mut events = Vec::new();
    let mut cum_recomb_cost: u64 = 0;
    for block in blocks.values_mut() {
        block.sort_unstable();
        if block.len() <= 2 {
            continue;
        }
        let tv: Vec<u32> = block.iter().map(|p| transmission_vector[position_to_index[p]]).collect();
        let cost: Vec<u32> = block.iter().map(|p| recombcost[position_to_index[p]]).collect();
        for i in 2..block.len() {
            if tv[i - 1] != tv[i] {
                events.push(RecombinationEvent {
                    position1: block[i - 1],
                    position2: block[i],
                    transmitted_hap_father1: tv[i - 1] % 2,
                    transmitted_hap_father2: tv[i] % 2,
                    transmitted_hap_mother1: tv[i - 1] / 2,
                    transmitted_hap_mother2: tv[i] / 2,
                    recombination_cost: cost[i],
                });
                cum_recomb_cost += cost[i] as u64;
            }
        }
    }
    log::info!("Cost accounted for by recombination events: {cum_recomb_cost}");
    events.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Ok(events)
}

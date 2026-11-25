//! Port of `pedigree.{h,cpp}` and `pedigreepartitions.{h,cpp}`.

use std::collections::HashMap;

use anyhow::{Result, anyhow};

use super::genotype::{Genotype, PhredGenotypeLikelihoods};

/// (father index, mother index, child index)
pub type Triple = [usize; 3];

#[derive(Default)]
pub struct Pedigree {
    variant_count: Option<usize>,
    triples: Vec<Triple>,
    individual_ids: Vec<u32>,
    id_to_index_map: HashMap<u32, usize>,
    genotypes: Vec<Vec<Genotype>>,
    genotype_likelihoods: Vec<Vec<Option<PhredGenotypeLikelihoods>>>,
}

impl Pedigree {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_individual(
        &mut self,
        id: u32,
        genotypes: Vec<Genotype>,
        genotype_likelihoods: Vec<Option<PhredGenotypeLikelihoods>>,
    ) {
        if self.variant_count.is_none() {
            self.variant_count = Some(genotypes.len());
        }
        assert_eq!(Some(genotypes.len()), self.variant_count);
        assert_eq!(genotype_likelihoods.len(), genotypes.len());
        self.genotypes.push(genotypes);
        self.genotype_likelihoods.push(genotype_likelihoods);
        self.individual_ids.push(id);
        self.id_to_index_map.insert(id, self.individual_ids.len() - 1);
    }

    pub fn add_relationship(&mut self, father_id: u32, mother_id: u32, child_id: u32) -> Result<()> {
        let t = [self.id_to_index(father_id)?, self.id_to_index(mother_id)?, self.id_to_index(child_id)?];
        self.triples.push(t);
        Ok(())
    }

    pub fn id_to_index(&self, id: u32) -> Result<usize> {
        self.id_to_index_map
            .get(&id)
            .copied()
            .ok_or_else(|| anyhow!("Individual with ID {id} not present in pedigree."))
    }

    pub fn index_to_id(&self, index: usize) -> u32 {
        self.individual_ids[index]
    }

    pub fn genotype(&self, individual_index: usize, variant_index: usize) -> &Genotype {
        &self.genotypes[individual_index][variant_index]
    }

    pub fn genotype_likelihoods(
        &self,
        individual_index: usize,
        variant_index: usize,
    ) -> Option<&PhredGenotypeLikelihoods> {
        self.genotype_likelihoods[individual_index][variant_index].as_ref()
    }

    pub fn size(&self) -> usize {
        self.individual_ids.len()
    }

    pub fn triple_count(&self) -> usize {
        self.triples.len()
    }

    pub fn triples(&self) -> &[Triple] {
        &self.triples
    }
}

/// Partitions of haplotypes that are identical by descent for a fixed
/// transmission vector.
pub struct PedigreePartitions {
    partition_count: usize,
    haplotype_to_partition_map: Vec<[i32; 2]>,
}

impl PedigreePartitions {
    pub fn new(pedigree: &Pedigree, transmission_vector: u32) -> Self {
        let n = pedigree.size();
        let mut pp = PedigreePartitions {
            partition_count: 2 * (n - pedigree.triple_count()),
            haplotype_to_partition_map: vec![[-1, -1]; n],
        };
        let mut triple_indices = vec![-1i32; n];
        for (i, t) in pedigree.triples().iter().enumerate() {
            triple_indices[t[2]] = i as i32;
        }
        let mut p = 0;
        for (i, &ti) in triple_indices.iter().enumerate() {
            if ti == -1 {
                pp.haplotype_to_partition_map[i] = [p, p + 1];
                p += 2;
            }
        }
        for i in 0..n {
            pp.compute_rec(i, &triple_indices, pedigree, transmission_vector);
        }
        pp
    }

    fn compute_rec(&mut self, i: usize, triple_indices: &[i32], pedigree: &Pedigree, tv: u32) {
        if self.haplotype_to_partition_map[i][0] != -1 {
            return;
        }
        let ti = triple_indices[i];
        assert!(ti >= 0);
        let ti = ti as usize;
        let parent0 = pedigree.triples()[ti][0];
        let parent1 = pedigree.triples()[ti][1];
        self.compute_rec(parent0, triple_indices, pedigree, tv);
        self.compute_rec(parent1, triple_indices, pedigree, tv);
        let h0 = ((tv >> (2 * ti)) & 1 == 0) as usize;
        let h1 = ((tv >> (2 * ti + 1)) & 1 == 0) as usize;
        self.haplotype_to_partition_map[i] =
            [self.haplotype_to_partition_map[parent0][h0], self.haplotype_to_partition_map[parent1][h1]];
    }

    pub fn count(&self) -> usize {
        self.partition_count
    }

    pub fn haplotype_to_partition(&self, individual_index: usize, haplotype: usize) -> usize {
        self.haplotype_to_partition_map[individual_index][haplotype] as usize
    }
}

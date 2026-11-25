//! Port of `genotype.{h,cpp}`, `phredgenotypelikelihoods.{h,cpp}` and `binomial.cpp`.

use anyhow::{Result, bail};

pub const MAX_ALLELES: u32 = 16;
pub const MAX_PLOIDY: u32 = 15;

pub fn binomial_coefficient(n: i32, k: i32) -> i32 {
    if k < 0 || n < 0 || n < k {
        return 0;
    }
    let mut k = k;
    if k > n - k {
        k = n - k;
    }
    let mut result: i32 = 1;
    for i in 0..k {
        result = result.wrapping_mul(n - i);
        result /= i + 1;
    }
    result
}

/// Unordered genotype packed into 4-bit chunks; same encoding as the C++ class
/// so that equality and `get_index()` behave identically.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Genotype {
    gt: u64,
}

impl std::fmt::Debug for Genotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Genotype({})", self)
    }
}

impl std::fmt::Display for Genotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_none() {
            return write!(f, ".");
        }
        let p = self.ploidy();
        write!(f, "{}", self.get_position(p - 1))?;
        for i in 1..p {
            write!(f, "/{}", self.get_position(p - i - 1))?;
        }
        Ok(())
    }
}

impl Genotype {
    pub fn none() -> Self {
        Genotype { gt: 0 }
    }

    pub fn new(alleles: &[u32]) -> Result<Self> {
        let ploidy = alleles.len() as u32;
        if ploidy >= MAX_PLOIDY {
            bail!("Error: Maximum ploidy for genotype exceeded!");
        }
        let mut sorted = alleles.to_vec();
        sorted.sort_unstable();
        let mut g = Genotype { gt: 0 };
        for (i, &a) in sorted.iter().enumerate() {
            if a >= MAX_ALLELES {
                bail!("Error: Maximum alleles for genotype exceeded!");
            }
            g.set_position(ploidy - i as u32 - 1, a);
        }
        g.set_position(MAX_PLOIDY, ploidy);
        Ok(g)
    }

    /// Construct from canonical index (see C++ docs).
    pub fn from_index(index: u64, ploidy: u32) -> Result<Self> {
        let alleles = convert_index_to_alleles(index, ploidy);
        Genotype::new(&alleles)
    }

    fn get_position(&self, pos: u32) -> u32 {
        ((self.gt >> (pos * 4)) & 15) as u32
    }

    fn set_position(&mut self, pos: u32, allele: u32) {
        let set_mask = (allele as u64) << (pos * 4);
        let delete_mask = !(15u64 << (pos * 4));
        self.gt &= delete_mask;
        self.gt |= set_mask;
    }

    pub fn ploidy(&self) -> u32 {
        self.get_position(MAX_PLOIDY)
    }

    pub fn is_none(&self) -> bool {
        self.ploidy() == 0
    }

    /// Alleles in the internal order, i.e. *descending* (as `as_vector()` in C++).
    pub fn as_vector(&self) -> Vec<u32> {
        (0..self.ploidy()).map(|i| self.get_position(i)).collect()
    }

    pub fn get_index(&self) -> u64 {
        let ploidy = self.ploidy();
        let mut index: u32 = 0;
        for (i, k) in (0..ploidy).zip(1u32..) {
            let allele = self.get_position(ploidy - i - 1);
            index = index.wrapping_add(binomial_coefficient((k + allele) as i32 - 1, allele as i32 - 1) as u32);
        }
        index as u64
    }

    pub fn is_homozygous(&self) -> bool {
        if self.is_none() {
            return false;
        }
        let a = self.get_position(0);
        (1..self.ploidy()).all(|i| self.get_position(i) == a)
    }

    pub fn is_diploid_and_biallelic(&self) -> bool {
        self.ploidy() == 2 && (0..2).all(|i| self.get_position(i) <= 1)
    }

    pub fn code(&self) -> u64 {
        self.gt
    }
}

pub fn convert_index_to_alleles(index: u64, ploidy: u32) -> Vec<u32> {
    let mut genotype = vec![0u32; ploidy as usize];
    let mut pth = ploidy;
    let mut max_allele_index = index as u32;
    let mut leftover = index as u32;
    while pth > 0 {
        let mut allele_index: u32 = 0;
        while allele_index <= max_allele_index {
            let i = binomial_coefficient((pth + allele_index) as i32 - 1, pth as i32) as u32;
            if i >= leftover || allele_index == max_allele_index {
                if i > leftover {
                    allele_index -= 1;
                }
                leftover -= binomial_coefficient((pth + allele_index) as i32 - 1, pth as i32) as u32;
                pth -= 1;
                max_allele_index = allele_index;
                genotype[pth as usize] = allele_index;
                break;
            }
            allele_index += 1;
        }
    }
    genotype
}

#[derive(Clone, Debug, PartialEq)]
pub struct PhredGenotypeLikelihoods {
    pub gl: Vec<f64>,
    pub ploidy: u32,
    pub nr_alleles: u32,
}

impl PhredGenotypeLikelihoods {
    pub fn new(gl: Vec<f64>, ploidy: u32, nr_alleles: u32) -> Result<Self> {
        let expected = binomial_coefficient((ploidy + nr_alleles) as i32 - 1, nr_alleles as i32 - 1) as usize;
        if expected != gl.len() {
            bail!("Error: wrong number of given genotype likelihoods given.");
        }
        Ok(PhredGenotypeLikelihoods { gl, ploidy, nr_alleles })
    }

    pub fn get(&self, genotype: &Genotype) -> f64 {
        self.gl[genotype.get_index() as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diploid_indices() {
        assert_eq!(Genotype::new(&[0, 0]).unwrap().get_index(), 0);
        assert_eq!(Genotype::new(&[1, 0]).unwrap().get_index(), 1);
        assert_eq!(Genotype::new(&[1, 1]).unwrap().get_index(), 2);
        assert_eq!(Genotype::new(&[0, 2]).unwrap().get_index(), 3);
        assert_eq!(Genotype::new(&[0, 1]).unwrap(), Genotype::new(&[1, 0]).unwrap());
        assert_eq!(Genotype::new(&[0, 1]).unwrap().as_vector(), vec![1, 0]);
        for i in 0..10 {
            assert_eq!(Genotype::from_index(i, 2).unwrap().get_index(), i);
            assert_eq!(Genotype::from_index(i, 4).unwrap().get_index(), i);
        }
    }
}

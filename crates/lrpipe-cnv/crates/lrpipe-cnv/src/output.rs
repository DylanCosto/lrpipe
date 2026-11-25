//! Spectre's output files: `{sample}_cnv.bed.gz`, `{sample}.vcf.gz` (both with `.tbi`
//! indexes) and the `{sample}.spc.gz` intermediate (spec §4).

pub mod bed;
pub mod htsfloat;
pub mod spc;
pub mod vcf;

mod tabix;

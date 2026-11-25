//! The site configuration: reference, models, catalogues and tool binaries, shared by every
//! sample run at a site.
//!
//! Relative paths in the file are taken from the file's directory, with its symbolic links,
//! `.` and `..` resolved: the same paths however `--config` was spelt (LRUN-07). A tool without a
//! path is looked up on `PATH` under its binary's usual name.

use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;

/// Why the site configuration could not be loaded.
#[derive(Debug, thiserror::Error)]
pub enum ConfigError {
    /// The file could not be read.
    #[error("reading {path}")]
    Read {
        /// The configuration file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: std::io::Error,
    },
    /// The file is not valid TOML for a site configuration.
    #[error("parsing {path}")]
    Parse {
        /// The configuration file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: Box<toml::de::Error>,
    },
}

/// Everything the pipeline needs besides the sample.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Site {
    /// Reference FASTA with a `.fai` index; also decodes CRAM input.
    pub reference: PathBuf,
    /// Tool binaries.
    #[serde(default)]
    pub tools: Tools,
    /// Small variant calling (lrpipe-snv); a run calling SNVs needs a model from it.
    #[serde(default)]
    pub snv: Snv,
    /// Structural variant calling (`lrpipe sv`).
    #[serde(default)]
    pub sv: Sv,
    /// Tandem repeat genotyping (lrpipe-tr, inside `lrpipe sv`'s read pass); none without it.
    pub repeats: Option<Repeats>,
    /// Copy number calling (lrpipe-cnv).
    #[serde(default)]
    pub cnv: Cnv,
    /// Methylation from the reads' `MM`/`ML` tags (lrpipe-methyl).
    #[serde(default)]
    pub methylation: Methylation,
    /// Aligning unaligned reads (bonefish).
    #[serde(default)]
    pub align: Align,
}

/// Paths to the tool binaries.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Tools {
    /// The lrpipe binary the stages run with: SNVs (`lrpipe snv`), the read pass (`lrpipe sv`),
    /// CNVs (`lrpipe cnv`) and methylation (`lrpipe methyl`) (default: the one running the
    /// pipeline).
    pub lrpipe: Option<PathBuf>,
    /// The bonefish binary the align stage maps unaligned reads with (default: a `bonefish` beside
    /// the lrpipe binary, else on `PATH`).
    pub bonefish: Option<PathBuf>,
}

/// How lrpipe-snv runs.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snv {
    /// A directory of Clair3's models under their standard names (`r1041_e82_400bps_sup_v500`,
    /// `hifi_revio`, ...), of which each run takes the one for its reads' header
    /// ([`crate::snv_model`]).
    pub models: Option<PathBuf>,
    /// Clair3 model directory for PacBio HiFi reads, used for every HiFi run in place of one
    /// from `models`.
    pub hifi_model: Option<PathBuf>,
    /// Clair3 model directory for Oxford Nanopore reads, used for every ONT run in place of one
    /// from `models`.
    pub ont_model: Option<PathBuf>,
    /// Where the networks run: `cpu`, `cuda:N` or `metal`.
    #[serde(default = "cpu")]
    pub device: String,
    /// NVRTC (`libnvrtc.so.12`), which lrpipe-snv's GPU calling compiles its kernels with, when
    /// it is not on the library path, `$CUDA_HOME`, `$CUDA_PATH` or `/usr/local/cuda`.
    pub nvrtc: Option<PathBuf>,
}

fn cpu() -> String {
    "cpu".into()
}

impl Default for Snv {
    fn default() -> Self {
        Self { models: None, hifi_model: None, ont_model: None, device: cpu(), nvrtc: None }
    }
}

/// How `lrpipe sv` calls SVs.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sv {
    /// Tandem repeat BED for better SV calls in repeats (Sniffles2's `--tandem-repeats`).
    pub tandem_repeats: Option<PathBuf>,
}

/// The tandem repeat catalogue.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Repeats {
    /// Loci to genotype: BED with TRGT-style `ID=..;MOTIFS=..`, or a motif in column 4.
    pub catalog: PathBuf,
}

/// lrpipe-cnv's inputs.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Cnv {
    /// Blacklist BED (Spectre's `grch38_blacklist_spectre.bed`).
    pub blacklist: Option<PathBuf>,
    /// The reference's `.mdr` of N runs, shared by the runs that name it; made from the
    /// reference if it does not exist. Without it, each run makes one in its CNV directory for
    /// the reference it uses.
    pub metadata: Option<PathBuf>,
}

/// How the align stage maps unaligned reads.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Align {
    /// bonefish's index of the reference (`bonefish -d ref.raw --idx-raw ref.fa`, as `lrpipe
    /// setup --bonefish` builds it), which loads in a second; without it, bonefish indexes the
    /// reference itself on every run. It must be built from `reference`.
    pub index: Option<PathBuf>,
}

/// How lrpipe-methyl reports methylation.
#[derive(Clone, Copy, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Methylation {
    /// Which tables.
    #[serde(default)]
    pub report: MethylationReport,
    /// How the CpG tables estimate a site's level.
    #[serde(default)]
    pub estimate: Estimate,
}

/// Which methylation tables a run writes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum MethylationReport {
    /// The platform's own tool's: bedMethyl for Oxford Nanopore, CpG tables for PacBio HiFi.
    #[default]
    Auto,
    /// modkit's bedMethyl: 5mC and 5hmC at CpGs, both strands together.
    Bedmethyl,
    /// pb-CpG-tools' CpG tables: all reads, haplotype 1 and haplotype 2.
    Cpg,
    /// No methylation stage.
    None,
}

/// How a CpG table estimates a site's methylation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Estimate {
    /// The share of reads called methylated (pb-CpG-tools' count mode, reproduced exactly).
    #[default]
    Count,
    /// lrpipe-methyl's fixed model of reads' values, fitted against EM-seq (in place of pb-CpG-tools'
    /// trained model, which it beats held out).
    Model,
}

impl Site {
    /// Reads a site configuration from a TOML file.
    ///
    /// # Errors
    ///
    /// If the file cannot be read or is not a valid site configuration.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, ConfigError> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|source| ConfigError::Read { path: path.into(), source })?;
        let dir = path.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
        let base = fs::canonicalize(dir).or_else(|_| std::path::absolute(dir)).unwrap_or_else(|_| dir.to_path_buf());
        Self::parse(&text, &base).map_err(|source| ConfigError::Parse { path: path.into(), source: Box::new(source) })
    }

    /// Parses a site configuration, taking relative paths from `base`.
    ///
    /// # Errors
    ///
    /// If `text` is not a valid site configuration.
    pub fn parse(text: &str, base: &Path) -> Result<Self, toml::de::Error> {
        let mut site: Self = toml::from_str(text)?;
        site.resolve(base);
        Ok(site)
    }

    fn resolve(&mut self, base: &Path) {
        let fix = |p: &mut PathBuf| {
            if p.is_relative() {
                *p = base.join(&*p);
            }
        };
        let tools = &mut self.tools;
        let optional = [
            &mut tools.lrpipe,
            &mut tools.bonefish,
            &mut self.align.index,
            &mut self.snv.models,
            &mut self.snv.hifi_model,
            &mut self.snv.ont_model,
            &mut self.snv.nvrtc,
            &mut self.sv.tandem_repeats,
            &mut self.cnv.blacklist,
            &mut self.cnv.metadata,
        ];
        let catalog = self.repeats.as_mut().map(|r| &mut r.catalog);
        let paths = std::iter::once(&mut self.reference).chain(optional.into_iter().flatten()).chain(catalog);
        for p in paths {
            fix(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_paths_are_taken_from_the_file() {
        let site = Site::parse(
            r#"
            reference = "ref/genome.fa"
            [tools]
            lrpipe = "bin/lrpipe"
            bonefish = "bin/bonefish"
            [snv]
            models = "models"
            hifi_model = "models/hifi"
            nvrtc = "cuda/libnvrtc.so.12"
            [repeats]
            catalog = "catalog.bed"
            [align]
            index = "ref/genome.raw"
            "#,
            Path::new("/site"),
        )
        .expect("valid configuration");
        assert_eq!(site.reference, Path::new("/site/ref/genome.fa"));
        assert_eq!(site.tools.lrpipe.as_deref(), Some(Path::new("/site/bin/lrpipe")));
        assert_eq!(site.tools.bonefish.as_deref(), Some(Path::new("/site/bin/bonefish")));
        assert_eq!(site.align.index.as_deref(), Some(Path::new("/site/ref/genome.raw")));
        assert_eq!(site.snv.models.as_deref(), Some(Path::new("/site/models")));
        assert_eq!(site.snv.hifi_model.as_deref(), Some(Path::new("/site/models/hifi")));
        assert_eq!(site.snv.device, "cpu");
        assert_eq!(site.snv.nvrtc.as_deref(), Some(Path::new("/site/cuda/libnvrtc.so.12")));
        assert_eq!(site.repeats.map(|r| r.catalog), Some(PathBuf::from("/site/catalog.bed")));
        assert_eq!(site.methylation.report, MethylationReport::Auto);
    }

    #[test]
    fn paths_are_the_same_however_the_file_is_named() {
        // LRUN-07: `--config site.toml`, `./site.toml` and its absolute path give the same paths,
        // which the stages' fingerprints hold.
        let dir = std::env::temp_dir().join(format!("lrpipe-config-{}", std::process::id()));
        fs::create_dir_all(dir.join("sub")).expect("scratch directory");
        let text = "reference = \"ref/genome.fa\"\n[tools]\nlrpipe = \"lrpipe\"\n[snv]\n";
        fs::write(dir.join("site.toml"), text).expect("configuration");
        let real = fs::canonicalize(&dir).expect("canonical");
        for spelling in [dir.join("site.toml"), dir.join("sub/../site.toml"), dir.join("./site.toml")] {
            let site = Site::load(&spelling).expect("loads");
            assert_eq!(site.reference, real.join("ref/genome.fa"), "{}", spelling.display());
            assert_eq!(site.tools.lrpipe, Some(real.join("lrpipe")), "a tool beside it is not looked up on PATH");
        }
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn methylation_choices_parse() {
        let site = Site::parse(
            "reference = \"r.fa\"\n[snv]\n[methylation]\nreport = \"cpg\"\nestimate = \"model\"\n",
            Path::new("/"),
        )
        .expect("valid configuration");
        assert_eq!(site.methylation.report, MethylationReport::Cpg);
        assert_eq!(site.methylation.estimate, Estimate::Model);
    }

    #[test]
    fn a_configuration_without_snv_parses() {
        // Without [snv], parsing stopped at "missing field `snv`"; a run calling SNVs now says
        // which key gives it a model.
        let site = Site::parse("reference = \"r.fa\"\n", Path::new("/")).expect("valid configuration");
        assert_eq!((site.snv.models, site.snv.device.as_str()), (None, "cpu"));
    }

    #[test]
    fn unknown_keys_are_rejected() {
        let err = Site::parse("reference = \"r.fa\"\n[snv]\ndevise = \"cuda:0\"\n", Path::new("/"));
        assert!(err.is_err(), "a misspelt key must not be ignored");
    }
}

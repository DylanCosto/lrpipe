//! The run's manifest: inputs, tools, each stage's command, outcome, time and memory, and the
//! output files; rewritten after every stage so a stopped run leaves an accurate one.

use std::collections::BTreeMap;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::plan::{Stage, StageName};
use crate::provenance::{FileRecord, ToolRecord};

/// Where the whole run stands.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    /// Still running (or stopped without finishing the manifest).
    Running,
    /// Every stage succeeded or was already done.
    Complete,
    /// A stage failed.
    Failed,
    /// Stopped by a signal; the stages that were running were killed.
    Interrupted,
}

/// Where one stage stands.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StageStatus {
    /// Waiting for its turn.
    Pending,
    /// Running now.
    Running,
    /// Done by an earlier run with the same command and inputs.
    Skipped,
    /// Ran and succeeded.
    Succeeded,
    /// Ran and failed.
    Failed,
    /// Killed because another stage failed.
    Cancelled,
    /// Never started because another stage failed.
    NotRun,
}

/// One stage in the manifest.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StageRecord {
    /// Which stage.
    pub name: StageName,
    /// Where it stands.
    pub status: StageStatus,
    /// Program and arguments.
    pub command: Vec<String>,
    /// Its log.
    pub log: PathBuf,
    /// When it started (UTC).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub started: Option<String>,
    /// Wall-clock seconds.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wall_seconds: Option<f64>,
    /// User CPU seconds.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub user_seconds: Option<f64>,
    /// Kernel CPU seconds.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub system_seconds: Option<f64>,
    /// Peak resident memory in GiB.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_rss_gib: Option<f64>,
    /// How it ended, when it did not succeed.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

impl StageRecord {
    /// A pending record of `stage`, logging to `log`.
    #[must_use]
    pub fn pending(stage: &Stage, log: PathBuf) -> Self {
        let command = std::iter::once(stage.program.as_os_str())
            .chain(stage.args.iter().map(|a| a.as_os_str()))
            .map(|a| a.to_string_lossy().into_owned())
            .collect();
        Self {
            name: stage.name,
            status: StageStatus::Pending,
            command,
            log,
            started: None,
            wall_seconds: None,
            user_seconds: None,
            system_seconds: None,
            max_rss_gib: None,
            error: None,
        }
    }

    /// Takes the start, times and memory of `earlier`, the record of the run that did this
    /// stage, for a stage skipped as already done.
    pub fn keep_measurements(&mut self, earlier: &Self) {
        self.started.clone_from(&earlier.started);
        self.wall_seconds = earlier.wall_seconds;
        self.user_seconds = earlier.user_seconds;
        self.system_seconds = earlier.system_seconds;
        self.max_rss_gib = earlier.max_rss_gib;
    }
}

/// The stage records of the manifest at `path`, or none if it cannot be read.
#[must_use]
pub fn earlier_stages(path: &Path) -> Vec<StageRecord> {
    #[derive(Deserialize)]
    struct Stages {
        stages: Vec<StageRecord>,
    }
    fs::read_to_string(path)
        .ok()
        .and_then(|text| serde_json::from_str::<Stages>(&text).ok())
        .map_or_else(Vec::new, |m| m.stages)
}

/// Whether `path` is a manifest an lrpipe run wrote: a JSON object with lrpipe's `sample` or
/// `stages` fields. Another tool's `manifest.json` is not.
#[must_use]
pub fn is_lrpipe(path: &Path) -> bool {
    let Ok(text) = fs::read_to_string(path) else { return false };
    serde_json::from_str::<serde_json::Value>(&text)
        .ok()
        .and_then(|v| v.as_object().map(|o| o.contains_key("sample") || o.contains_key("stages")))
        .unwrap_or(false)
}

/// The sample named in the manifest at `path`, or none if it cannot be read.
#[must_use]
pub fn earlier_sample(path: &Path) -> Option<String> {
    #[derive(Deserialize)]
    struct Named {
        name: String,
    }
    #[derive(Deserialize)]
    struct Sampled {
        sample: Named,
    }
    let text = fs::read_to_string(path).ok()?;
    serde_json::from_str::<Sampled>(&text).ok().map(|m| m.sample.name)
}

/// The sample as given.
#[derive(Clone, Debug, Serialize)]
pub struct SampleRecord {
    /// Its name.
    pub name: String,
    /// Its platform.
    pub platform: String,
    /// Its sex chromosomes.
    pub sex: String,
    /// What the index said, when the sex was read from it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sex_from_index: Option<SexEvidence>,
    /// Genome-wide ploidy when the caller uses one value for every contig.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ploidy: Option<u8>,
    /// The aligned reads the run calls from: the BAMs or CRAMs given (several are merged for the
    /// stages), or for unaligned input the align stage's BAM once it is written (none before).
    /// Always a list.
    pub reads: Vec<FileRecord>,
    /// The unaligned input, which the align stage mapped.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub unaligned_reads: Vec<FileRecord>,
    /// Where the SVs' phases, the tandem repeats' alleles and the CpG tables per haplotype come
    /// from: the run's phasing ([`HAPLOTYPES_RUN`]), the reads' own `HP` tags, or none
    /// ([`HAPLOTYPES_NONE`]).
    pub haplotypes: String,
}

/// [`SampleRecord::haplotypes`] of reads without `HP` tags.
pub const HAPLOTYPES_NONE: &str = "none (input not haplotagged)";

/// [`SampleRecord::haplotypes`] of haplotagged reads.
pub const HAPLOTYPES_INPUT: &str = "the input's HP tags (not lrpipe run's phasing)";

/// [`SampleRecord::haplotypes`] when the run haplotags the reads from its phased SNVs.
pub const HAPLOTYPES_RUN: &str = "lrpipe run's phasing of the SNVs (its haplotag list)";

/// chrX's and chrY's reads per base against the autosomes', from the reads' index.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct SexEvidence {
    /// chrX's.
    pub chrx: f64,
    /// chrY's.
    pub chry: f64,
}

/// The Clair3 model lrpipe-snv called with, and why ([`crate::snv_model::Choice`]).
#[derive(Clone, Debug, Serialize)]
pub struct ModelRecord {
    /// Its directory.
    pub path: PathBuf,
    /// How it was found: `--model`, `snv.ont_model` or `snv.hifi_model`, or `snv.models` for one
    /// chosen from the reads' header.
    pub by: String,
    /// Why it is the model.
    pub reason: String,
    /// What the run was warned of about it.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

/// The manifest.
#[derive(Clone, Debug, Serialize)]
pub struct Manifest {
    /// lrpipe's version.
    pub lrpipe: String,
    /// lrpipe's command line.
    pub command: Vec<String>,
    /// The output directory, which every output path below starts with.
    pub output: PathBuf,
    /// GiB the stages running together were allowed.
    pub memory_budget_gib: f64,
    /// When the run started (UTC).
    pub started: String,
    /// When it ended (UTC).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub finished: Option<String>,
    /// Where it stands.
    pub status: RunStatus,
    /// Wall-clock seconds from start to end.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wall_seconds: Option<f64>,
    /// Why the run stopped outside any stage: one from unaligned reads that failed after its
    /// align stage, before the others started.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// The sample.
    pub sample: SampleRecord,
    /// lrpipe-snv's model, when it runs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snv_model: Option<ModelRecord>,
    /// Reference data and catalogues, by role.
    pub resources: BTreeMap<String, FileRecord>,
    /// Tool binaries, by name.
    pub tools: BTreeMap<String, ToolRecord>,
    /// Each stage.
    pub stages: Vec<StageRecord>,
    /// The results, by kind.
    pub outputs: BTreeMap<String, PathBuf>,
    /// What the run warned of before its stages ran.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

impl Manifest {
    /// The record of `name`.
    ///
    /// # Panics
    ///
    /// If the manifest has no such stage.
    pub fn stage_mut(&mut self, name: StageName) -> &mut StageRecord {
        self.stages.iter_mut().find(|s| s.name == name).unwrap_or_else(|| panic!("no stage {name} in the manifest"))
    }

    /// Writes the manifest as JSON to `path`, replacing it atomically.
    ///
    /// # Errors
    ///
    /// If the file cannot be written.
    pub fn write(&self, path: &Path) -> io::Result<()> {
        let json = serde_json::to_string_pretty(self).map_err(io::Error::other)?;
        let tmp = path.with_extension("json.tmp");
        fs::write(&tmp, json + "\n")?;
        fs::rename(tmp, path)
    }
}

impl Serialize for StageName {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for StageName {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let name = String::deserialize(d)?;
        StageName::named(&name).ok_or_else(|| serde::de::Error::custom(format!("unknown stage {name}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every stage's record is read back from an earlier manifest, merge-reads's too: one name
    /// not read made the whole list unread, and every stage lost its times and measured peak.
    #[test]
    fn earlier_stages_of_every_name_are_read_back() {
        let path = std::env::temp_dir().join(format!("lrpipe-manifest-{}.json", std::process::id()));
        let records: Vec<String> = StageName::ALL
            .iter()
            .map(|s| format!(r#"{{"name":"{s}","status":"succeeded","command":[],"log":"l","wall_seconds":2.0}}"#))
            .collect();
        fs::write(&path, format!(r#"{{"stages":[{}]}}"#, records.join(","))).expect("written");
        let stages = earlier_stages(&path);
        assert_eq!(stages.iter().map(|s| s.name).collect::<Vec<_>>(), StageName::ALL);
        assert!(stages.iter().all(|s| s.wall_seconds == Some(2.0)));
        fs::remove_file(path).expect("removed");
    }
}

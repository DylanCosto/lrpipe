//! Which Clair3 model lrpipe-snv calls a sample's reads with, from the reads' header.
//!
//! A Clair3 model is trained on one platform's reads, and for Oxford Nanopore on one chemistry
//! and basecall model; the wrong one costs accuracy without an error. The header says what made
//! the reads: dorado writes `@RG ... PL:ONT DS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0`,
//! and PacBio `@RG ... PL:PACBIO PM:REVIO DS:READTYPE=CCS;...`. [`choose`] picks the model for
//! them from the Clair3 models ([`table`]) in the site's `[snv] models` directory, or checks a
//! model given explicitly against them:
//!
//! * a model for another platform or chemistry than the reads' is an error (a [`Mismatch::Warn`]
//!   run goes on with a warning), as are Oxford Nanopore reads called PacBio or the reverse;
//! * the nearest model the directory holds, where it holds none for the reads' basecall model or
//!   Clair3 has none, is used with a warning that names it and why; so is a given model for
//!   another basecall model, and a given model lrpipe does not know (a custom or fine-tuned one),
//!   which it cannot check;
//! * read groups basecalled differently leave no one model to choose: an error, unless one is
//!   given;
//! * a move-table (`_with_mv`) model, which lrpipe-snv runs only with dwell time, is refused.

pub mod table;

use std::fmt;
use std::io;
use std::path::{Path, PathBuf};

use pasta::sam;

use crate::plan::Platform;
use table::{Basecaller, Model, Pick, Reads, Source, Tier, Version};

/// What a run does when the reads' header contradicts the model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mismatch {
    /// Stops before running anything.
    #[default]
    Error,
    /// Runs with the model, with a warning.
    Warn,
}

/// A model given explicitly, and where.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Given {
    /// The model's directory.
    pub path: PathBuf,
    /// Where it was given: `--model`, `snv.ont_model`, ...
    pub by: String,
}

/// What [`choose`] chooses from.
#[derive(Clone, Debug)]
pub struct Request<'a> {
    /// The platform the run was given.
    pub platform: Platform,
    /// A model given explicitly, which is checked rather than chosen.
    pub given: Option<Given>,
    /// The site's directory of Clair3's models under their standard names (`[snv] models`).
    pub models: Option<&'a Path>,
    /// What a contradiction does.
    pub mismatch: Mismatch,
    /// Whether the reads are FASTQ without read group tags, which has no header: where the header
    /// would say nothing, that is said instead.
    pub fastq: bool,
}

/// The model a run calls with, and why.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Choice {
    /// Its directory.
    pub path: PathBuf,
    /// How it was found: `--model`, `snv.ont_model`, or `snv.models` for one chosen from the
    /// header.
    pub by: String,
    /// Why it is the model, in a sentence.
    pub reason: String,
    /// What the run should be warned of.
    pub warnings: Vec<String>,
}

/// Why no model could be chosen.
#[derive(Debug, thiserror::Error)]
pub enum ModelError {
    /// The header could not be read.
    #[error("reading the header of {}", path.display())]
    Header {
        /// The reads.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// No model is configured or given for the platform.
    #[error(
        "no Clair3 model for {0} reads: the site configuration has neither [snv] {0}_model nor [snv] models, and \
         no --model was given"
    )]
    NoModel(&'static str),
    /// The header contradicts the model.
    #[error("{0}; give the model and --platform the reads were made for, or --model-mismatch warn to run anyway")]
    Mismatch(String),
    /// A move-table model, which lrpipe-snv runs only with dwell time, which `lrpipe run` does
    /// not ask for.
    #[error(
        "{} is a move-table (signal-aware, _with_mv) model, which lrpipe snv runs only with --enable_dwell_time, \
         which lrpipe run does not pass: give another model",
        .0.display()
    )]
    MoveTable(PathBuf),
    /// Read groups were basecalled differently, so no one model fits them all.
    #[error(
        "the reads' read groups were made differently ({}), so no one model fits them all: give the model with \
         --model",
        .0.join("; ")
    )]
    Several(Vec<String>),
    /// The header does not say enough to choose a model.
    #[error("{reason}: give the model with --model or snv.{platform}_model")]
    Unknown {
        /// The platform.
        platform: &'static str,
        /// What the header lacks.
        reason: String,
    },
    /// The models directory is not there.
    #[error("the site configuration's [snv] models {} does not exist: give the directory of Clair3's models", .0.display())]
    NoModelsDir(PathBuf),
    /// The models directory is one model's own.
    #[error(
        "the site configuration's [snv] models {} is a model's own directory (it holds pileup.pt): give the \
         directory that holds the models ({}), or this one as [snv] {platform}_model",
        dir.display(),
        dir.parent().map_or_else(|| "its parent".to_owned(), |p| p.display().to_string())
    )]
    OneModel {
        /// The directory.
        dir: PathBuf,
        /// The run's platform.
        platform: &'static str,
    },
    /// The models directory holds no model for the reads.
    #[error(
        "{} holds no Clair3 model for {reads}: add {name} ({origin}), or give a model with --model",
        dir.display()
    )]
    NotInstalled {
        /// The models directory.
        dir: PathBuf,
        /// The reads.
        reads: String,
        /// The model for them.
        name: &'static str,
        /// Where it comes from.
        origin: Source,
    },
}

/// The network `name` (`pileup`, `full_alignment`) of the model directory `model`, as lrpipe-snv
/// finds it: `NAME.pt`, or `NAME` when only that is there (as Clair3's `Path.with_suffix`).
#[must_use]
pub fn network(model: &Path, name: &str) -> PathBuf {
    let bare = model.join(name);
    let pt = bare.with_extension("pt");
    if !pt.exists() && bare.is_file() { bare } else { pt }
}

/// The SAM header text of the BAM or CRAM at `reads`.
///
/// # Errors
///
/// If the file cannot be opened or its header read.
pub fn header_text(reads: &Path) -> Result<String, ModelError> {
    let err = |source| ModelError::Header { path: reads.to_path_buf(), source };
    let header = crate::input::read_header(reads).map_err(err)?;
    text_of(&header).map_err(err)
}

/// `header` as SAM text, which [`choose`] reads.
///
/// # Errors
///
/// If a field of the header cannot be written.
pub fn text_of(header: &sam::Header) -> io::Result<String> {
    let mut text = Vec::new();
    sam::io::Writer::new(&mut text).write_header(header)?;
    Ok(String::from_utf8_lossy(&text).into_owned())
}

/// PacBio read types that are not HiFi. pbmm2 writes `READTYPE=UNKNOWN` for reads aligned from
/// FASTQ or FASTA (pbbam's default read type), which may well be HiFi.
const NOT_HIFI: [&str; 4] = ["SUBREAD", "SCRAP", "ZMW", "HQREGION"];

/// What made a read group's reads, by its `@RG` line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Sequenced {
    /// Oxford Nanopore, with the basecall model in `DS:basecall_model=`, if any.
    Ont(Option<BasecallModel>),
    /// PacBio, from the instrument `PM`, with `DS`'s `READTYPE` (`CCS` for HiFi).
    PacBio {
        /// The instrument (`REVIO`, `SEQUELII`).
        instrument: Option<String>,
        /// The read type (`CCS`, `SUBREAD`, `UNKNOWN`).
        read_type: Option<String>,
    },
    /// Another platform, as `PL` names it.
    Other(String),
}

impl Sequenced {
    /// What a read group says, if it says anything about how its reads were made. A
    /// `basecall_model` (or PacBio's `READTYPE`) says more than a `PL` of another platform,
    /// which aligners stamp out of habit (`minimap2 -R '@RG\tID:x\tPL:ILLUMINA'`).
    #[must_use]
    pub fn of(read_group: &str) -> Option<Self> {
        let field =
            |tag: &str| read_group.split('\t').skip(1).find_map(|f| f.strip_prefix(tag)).filter(|v| !v.is_empty());
        let description = field("DS:").unwrap_or("");
        let basecall =
            description.split([' ', ';']).find_map(|kv| kv.strip_prefix("basecall_model=")).map(BasecallModel::parse);
        let read_type = description.split([' ', ';']).find_map(|kv| kv.strip_prefix("READTYPE=")).map(str::to_owned);
        let pacbio = || Self::PacBio { instrument: field("PM:").map(str::to_uppercase), read_type: read_type.clone() };
        match field("PL:").map(str::to_uppercase).as_deref() {
            Some("ONT" | "ONTUNKNOWN" | "NANOPORE" | "OXFORD_NANOPORE") => Some(Self::Ont(basecall)),
            Some("PACBIO" | "PACBIO_SMRT") => Some(pacbio()),
            _ if basecall.is_some() => Some(Self::Ont(basecall)),
            _ if read_type.is_some() => Some(pacbio()),
            Some(other) => Some(Self::Other(other.to_owned())),
            None => None,
        }
    }

    /// The platform `--platform` names for these reads, if lrpipe calls them.
    fn platform(&self) -> Option<Platform> {
        match self {
            Self::Ont(_) => Some(Platform::Ont),
            Self::PacBio { .. } => Some(Platform::Hifi),
            Self::Other(_) => None,
        }
    }

    /// What the model depends on: the basecall model, or the instrument; none if the read group
    /// does not say.
    fn key(&self) -> Option<&str> {
        match self {
            Self::Ont(Some(b)) => Some(&b.name),
            Self::PacBio { instrument: Some(i), .. } => Some(i),
            _ => None,
        }
    }
}

impl fmt::Display for Sequenced {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ont(Some(b)) => write!(f, "Oxford Nanopore reads basecalled with {}", b.name),
            Self::Ont(None) => f.write_str("Oxford Nanopore reads (@RG PL:ONT, no basecall_model)"),
            Self::PacBio { instrument, read_type } => {
                write!(f, "PacBio reads")?;
                match (instrument, read_type) {
                    (Some(i), Some(t)) => write!(f, " (@RG PM:{i}, READTYPE={t})"),
                    (Some(i), None) => write!(f, " (@RG PM:{i})"),
                    (None, Some(t)) => write!(f, " (READTYPE={t})"),
                    (None, None) => Ok(()),
                }
            }
            Self::Other(platform) => write!(f, "{platform} reads (@RG PL:{platform})"),
        }
    }
}

/// A dorado (or Guppy) basecall model, as `DS:basecall_model=` names it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BasecallModel {
    /// The name as written (`dna_r10.4.1_e8.2_400bps_sup@v5.0.0`).
    pub name: String,
    /// The pore (`r10.4.1`), if the name is a DNA model's.
    pub pore: Option<String>,
    /// Bases per second (`400`), if the name says.
    pub bps: Option<u16>,
    /// The signal's sampling rate in kHz: as the name says (`5khz`), else 5 for dorado's models
    /// from v4.2.0 and 4 for earlier ones and Guppy's.
    pub khz: Option<u8>,
    /// The accuracy, if the name says.
    pub tier: Option<Tier>,
    /// dorado's version (`@v5.0.0`), else not known (Guppy's names have none).
    pub basecaller: Basecaller,
}

impl BasecallModel {
    /// Reads `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, Guppy's `dna_r10.4.1_e8.2_400bps_5khz_sup`, or
    /// a path ending in one; what it cannot read (an RNA model) is left unknown.
    #[must_use]
    pub fn parse(text: &str) -> Self {
        let name = text.rsplit('/').next().unwrap_or(text);
        let (stem, version) = name.split_once('@').map_or((name, None), |(s, v)| (s, Some(v)));
        let mut model = Self {
            name: name.to_owned(),
            pore: None,
            bps: None,
            khz: None,
            tier: None,
            basecaller: Basecaller::Unknown,
        };
        let Some(fields) = stem.strip_prefix("dna_") else { return model };
        for field in fields.split('_') {
            if let Some(pore) = field.strip_prefix('r').filter(|p| Version::parse(p).is_some()) {
                model.pore = Some(format!("r{pore}"));
            } else if let Some(bps) = field.strip_suffix("bps").and_then(|b| b.parse().ok()) {
                model.bps = Some(bps);
            } else if let Some(khz) = field.strip_suffix("khz").and_then(|k| k.parse().ok()) {
                model.khz = Some(khz);
            } else if let Some(tier) = Tier::parse(field) {
                model.tier = Some(tier);
            }
        }
        if let Some(v) = version.and_then(Version::parse) {
            model.basecaller = Basecaller::Dorado(v);
        }
        if model.pore.is_some() && model.khz.is_none() {
            let five = matches!(model.basecaller, Basecaller::Dorado(v) if v >= Version([4, 2, 0]));
            model.khz = Some(if five { 5 } else { 4 });
        }
        model
    }
}

/// The `@RG` lines of a SAM header.
fn read_groups(header: &str) -> impl Iterator<Item = &str> {
    header.lines().filter(|l| l.starts_with("@RG\t"))
}

/// What the read groups of `header` say made the reads, each once, in order.
#[must_use]
pub fn sequenced(header: &str) -> Vec<Sequenced> {
    let mut all: Vec<Sequenced> = Vec::new();
    for s in read_groups(header).filter_map(Sequenced::of) {
        if !all.contains(&s) {
            all.push(s);
        }
    }
    all
}

/// Whether the model directory `path` is a move-table model's, by its name, as lrpipe-snv tells.
fn is_move_table(path: &Path) -> bool {
    let name = path.file_name().map(|n| n.to_string_lossy().to_lowercase()).unwrap_or_default();
    name.contains("with_mv") || name.contains("with_move")
}

/// The model for the reads `header` describes, or the given one checked against them.
///
/// # Errors
///
/// If the header contradicts a given model ([`ModelError::Mismatch`], unless
/// [`Mismatch::Warn`]) or it is a move-table model; or, with none given, if there is no models
/// directory, the header does not say enough or names several basecall models, or the directory
/// holds no model for the reads.
pub fn choose(header: &str, request: &Request<'_>) -> Result<Choice, ModelError> {
    let platform = request.platform;
    let reads = sequenced(header);
    let mut warnings = Vec::new();
    let mut problems = Vec::new();
    for s in &reads {
        match s.platform() {
            Some(p) if p == platform => {}
            Some(_) => problems
                .push(format!("the header says these are {s}, not {} reads as --platform says", platform.as_str())),
            None => warnings.push(format!(
                "the header says these are {s}, which lrpipe does not call; taken as {} reads, as --platform says",
                platform.as_str()
            )),
        }
        if let Sequenced::PacBio { read_type: Some(t), .. } = s
            && t != "CCS"
        {
            if NOT_HIFI.contains(&t.as_str()) {
                problems.push(format!("the header says these are {s}, not HiFi (READTYPE=CCS) reads"));
            } else {
                warnings.push(format!("the header says these are {s}, not READTYPE=CCS; taken as HiFi reads"));
            }
        }
    }
    let ours = described(&reads, platform, &mut warnings);

    if let Some(given) = &request.given {
        if is_move_table(&given.path) {
            return Err(ModelError::MoveTable(given.path.clone()));
        }
        let name = given.path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        let reason = match table::find(&name) {
            Some(model) => {
                let before = (problems.len(), warnings.len());
                check(model, platform, &reads, &mut problems, &mut warnings);
                let differs = (problems.len(), warnings.len()) != before;
                // Read groups of another platform (a contradiction, an error unless let through)
                // say what made the reads too.
                let telling = if ours.is_empty() { reads.iter().collect() } else { ours.clone() };
                if telling.iter().all(|s| s.key().is_none()) {
                    let unsaid = if request.fastq {
                        "FASTQ has no header saying what made the reads"
                    } else {
                        "the header does not say what made the reads (no @RG basecall_model or PM)"
                    };
                    warnings.push(format!("{unsaid}, so {name} (given by {}) is not checked against them", given.by));
                    let unsaid = unsaid.split(" (").next().unwrap_or(unsaid);
                    format!("Clair3's model for {}; {unsaid}", model.reads)
                } else if differs {
                    format!(
                        "Clair3's model for {}, which the header says the reads are not (see the warning)",
                        model.reads
                    )
                } else {
                    format!("Clair3's model for {}, checked against the header", model.reads)
                }
            }
            None => {
                warnings.push(format!(
                    "{} (given by {}) is not one of Clair3's models by name, so lrpipe cannot check it against the reads",
                    given.path.display(),
                    given.by
                ));
                "not one of Clair3's models by name, so not checked".to_owned()
            }
        };
        settle(problems, request.mismatch, &mut warnings)?;
        return Ok(Choice { path: given.path.clone(), by: given.by.clone(), reason, warnings });
    }

    let Some(dir) = request.models else { return Err(ModelError::NoModel(platform.as_str())) };
    if !dir.is_dir() {
        return Err(ModelError::NoModelsDir(dir.to_path_buf()));
    }
    if network(dir, "pileup").is_file() {
        return Err(ModelError::OneModel { dir: dir.to_path_buf(), platform: platform.as_str() });
    }
    settle(problems, request.mismatch, &mut warnings)?;
    let one = match ours[..] {
        // None of the platform's (a contradiction let through): what the header says, if anything.
        [] => reads.first(),
        [one] => Some(one),
        _ => return Err(ModelError::Several(ours.iter().map(ToString::to_string).collect())),
    };
    let installed = |m: &Model| dir.join(m.name).is_dir();
    let found = pick(one, platform, installed)?;
    let Some(chosen) = found.chosen else {
        let (name, origin) = (found.ideal.name, found.ideal.source);
        return Err(ModelError::NotInstalled { dir: dir.to_path_buf(), reads: found.reads, name, origin });
    };
    let why = if chosen != found.ideal {
        let ideal = found.ideal;
        let theirs = if found.exact { "Clair3's model for them" } else { "the nearest Clair3 model for them" };
        Some(format!(
            "{} does not hold {}, {theirs} ({}); {} ({}) is the nearest it holds",
            dir.display(),
            ideal.name,
            ideal.source,
            chosen.name,
            chosen.reads
        ))
    } else {
        found.not_exact.map(|why| {
            let what =
                if found.nearest { "the nearest" } else { "the HiFi model lrpipe takes for no known instrument" };
            format!("{why}; {} ({}) is {what}", chosen.name, chosen.reads)
        })
    };
    // The why once, in the warning; the model's reason points to it.
    let reason = match &why {
        None => format!("chosen from the header: Clair3's model for {}", found.reads),
        Some(why) => {
            warnings.push(format!("using {} for {}: {why}", chosen.name, found.reads));
            format!("chosen from the header for {}, not Clair3's model for them (see the warning)", found.reads)
        }
    };
    Ok(Choice { path: dir.join(chosen.name), by: "snv.models".into(), reason, warnings })
}

/// The read groups of the platform that say what made them, each once, with a warning for those
/// that do not beside them (a read group added by `minimap2 -R`); all of them if none says.
fn described<'a>(reads: &'a [Sequenced], platform: Platform, warnings: &mut Vec<String>) -> Vec<&'a Sequenced> {
    let ours: Vec<&Sequenced> = reads.iter().filter(|s| s.platform() == Some(platform)).collect();
    let mut keyed: Vec<&Sequenced> = Vec::new();
    for s in ours.iter().filter(|s| s.key().is_some()) {
        if keyed.iter().all(|k| k.key() != s.key()) {
            keyed.push(s);
        }
    }
    if keyed.is_empty() {
        return ours;
    }
    for s in ours.iter().filter(|s| s.key().is_none()) {
        warnings.push(format!(
            "a read group of {s} does not say what made its reads; taken to be made as the others ({})",
            keyed.iter().map(ToString::to_string).collect::<Vec<_>>().join("; ")
        ));
    }
    keyed
}

/// Contradictions are errors, or with [`Mismatch::Warn`] warnings.
fn settle(problems: Vec<String>, mismatch: Mismatch, warnings: &mut Vec<String>) -> Result<(), ModelError> {
    match (problems.is_empty(), mismatch) {
        (true, _) => Ok(()),
        (false, Mismatch::Error) => Err(ModelError::Mismatch(problems.join("; "))),
        (false, Mismatch::Warn) => {
            warnings.extend(problems.into_iter().map(|p| format!("{p} (--model-mismatch warn)")));
            Ok(())
        }
    }
}

/// Checks a given Clair3 model against the run's platform and the reads: another platform or
/// chemistry is a problem, another basecall model or instrument a warning.
fn check(
    model: &Model,
    platform: Platform,
    reads: &[Sequenced],
    problems: &mut Vec<String>,
    warnings: &mut Vec<String>,
) {
    let name = model.name;
    let model_platform = match model.reads {
        Reads::Ont { .. } => Some(Platform::Ont),
        Reads::Hifi { .. } => Some(Platform::Hifi),
        Reads::Illumina => None,
    };
    if model_platform != Some(platform) {
        problems.push(format!("{name} is Clair3's model for {}, not {} reads", model.reads, platform.as_str()));
        return;
    }
    if let Pick::No(why) = model.pick {
        warnings.push(format!("{name} is {why}"));
    }
    for s in reads {
        match s.platform() {
            // Warned of already.
            None => continue,
            Some(p) if Some(p) != model_platform => {
                problems.push(format!("{name} is Clair3's model for {}, but these are {s}", model.reads));
                continue;
            }
            Some(_) => {}
        }
        match (&model.reads, s) {
            (Reads::Ont { pore, bps, khz, tier, basecaller, .. }, Sequenced::Ont(Some(b))) => {
                let Some(their_pore) = &b.pore else {
                    warnings
                        .push(format!("lrpipe cannot read the basecall model {}, so {name} is not checked", b.name));
                    continue;
                };
                if their_pore != pore || b.bps.is_some_and(|their| their != *bps) {
                    problems.push(format!("{name} is Clair3's model for {}, but these are {s}", model.reads));
                } else if b.tier != Some(*tier) || b.basecaller != *basecaller || b.khz != Some(*khz) {
                    let own = nearest_ont(b, their_pore, |_| true).map(|(m, _)| m.name).filter(|own| *own != name);
                    warnings.push(format!(
                        "{name} is Clair3's model for {}, but these are {s}{}",
                        model.reads,
                        own.map(|own| format!(" (Clair3's model for them: {own})")).unwrap_or_default()
                    ));
                }
            }
            (Reads::Hifi { instruments, .. }, Sequenced::PacBio { instrument: Some(i), .. })
                if !instruments.contains(&i.as_str()) =>
            {
                warnings.push(format!("{name} is Clair3's model for {}, but these are {s}", model.reads));
            }
            _ => {}
        }
    }
}

/// The model for the reads, and the one picked from those `installed`.
struct Found {
    /// The reads, described.
    reads: String,
    /// Clair3's model for them, or the nearest.
    ideal: &'static Model,
    /// Whether `ideal` is the model for them.
    exact: bool,
    /// Why it is not, if it is not.
    not_exact: Option<String>,
    /// Whether `ideal`, when it is not the model for them, is the nearest to theirs (not just the
    /// model taken for HiFi reads of an instrument lrpipe has no model for, or none named).
    nearest: bool,
    /// The model picked of those installed: `ideal`, or the nearest to it.
    chosen: Option<&'static Model>,
}

/// The model lrpipe picks for `reads` (none: the header does not say) on `platform`, of all and
/// of those `installed`.
fn pick(
    reads: Option<&Sequenced>,
    platform: Platform,
    installed: impl Fn(&Model) -> bool,
) -> Result<Found, ModelError> {
    let unknown = |reason: String| ModelError::Unknown { platform: platform.as_str(), reason };
    match (platform, reads) {
        (Platform::Hifi, reads) => {
            let instrument = match reads {
                Some(Sequenced::PacBio { instrument, .. }) => instrument.as_deref(),
                _ => None,
            };
            // Trained on the instrument, then the nearest for it, then the first in the table.
            let rank = |m: &Model| match m.reads {
                Reads::Hifi { instruments, .. } if instrument.is_some_and(|i| instruments.contains(&i)) => 0,
                Reads::Hifi { nearest_for, .. } if instrument.is_some_and(|i| nearest_for.contains(&i)) => 1,
                _ => 2,
            };
            let hifi = || table::MODELS.iter().filter(|m| m.pick == Pick::Yes && matches!(m.reads, Reads::Hifi { .. }));
            let ideal = hifi().min_by_key(|m| rank(m)).expect("a HiFi model");
            let chosen = hifi().filter(|m| installed(m)).min_by_key(|m| rank(m));
            let exact = rank(ideal) == 0;
            let not_exact = (!exact).then(|| match (instrument, reads) {
                (Some(i), _) => format!("Clair3 has no model for the instrument {i}"),
                // Another platform's reads, let through with --model-mismatch warn: their PM is
                // no PacBio instrument.
                (None, Some(Sequenced::Ont(_) | Sequenced::Other(_))) => {
                    "these are not PacBio reads, so the header names no PacBio instrument".to_owned()
                }
                (None, _) => "the header does not name the instrument (@RG PM)".to_owned(),
            });
            let nearest = rank(ideal) == 1;
            let reads =
                reads.map_or_else(|| "PacBio HiFi reads (the header does not say)".to_owned(), ToString::to_string);
            Ok(Found { reads, ideal, exact, not_exact, nearest, chosen })
        }
        (Platform::Ont, None) => {
            Err(unknown("the reads' header has no @RG line saying how they were basecalled".into()))
        }
        (Platform::Ont, Some(Sequenced::Ont(None))) => Err(unknown(
            "the reads' header does not name their basecall model (@RG DS:basecall_model=), which the model depends on"
                .into(),
        )),
        (Platform::Ont, Some(s @ Sequenced::Ont(Some(b)))) => {
            let Some(pore) = &b.pore else {
                return Err(unknown(format!("lrpipe cannot read the basecall model {}", b.name)));
            };
            let (ideal, exact) = nearest_ont(b, pore, |_| true)
                .ok_or_else(|| unknown(format!("Clair3 has no model lrpipe knows for {s}")))?;
            let chosen = nearest_ont(b, pore, installed).map(|(m, _)| m);
            let not_exact = (!exact).then(|| format!("Clair3 has no model for {}", b.name));
            Ok(Found { reads: s.to_string(), ideal, exact, not_exact, nearest: true, chosen })
        }
        (Platform::Ont, Some(other)) => Err(unknown(format!("the header says these are {other}"))),
    }
}

/// The model lrpipe picks for dorado's `basecall` on `pore` of those `usable`, and whether it
/// is the model for those reads: of the models for the same pore and speed, the same sampling
/// rate, then the same tier (else the nearest, the more accurate of two), then the same basecaller
/// (Guppy's for a name without a version) at the nearest version (the newer of two), or the
/// newest where the version is not known.
fn nearest_ont(
    basecall: &BasecallModel,
    pore: &str,
    usable: impl Fn(&Model) -> bool,
) -> Option<(&'static Model, bool)> {
    let tier_distance = |t: Tier| basecall.tier.map_or(0, |their| (t as i32 - their as i32).abs());
    let version = |b: Basecaller| match b {
        Basecaller::Dorado(v) | Basecaller::Guppy(v) => Some(v),
        Basecaller::Unknown => None,
    };
    let same_caller = |b: Basecaller| match (basecall.basecaller, b) {
        (Basecaller::Unknown, Basecaller::Guppy(_)) => true,
        (theirs, ours) => std::mem::discriminant(&theirs) == std::mem::discriminant(&ours),
    };
    table::MODELS
        .iter()
        .filter(|m| m.pick == Pick::Yes && usable(m))
        .filter_map(|m| match m.reads {
            Reads::Ont { pore: p, bps, khz, tier, basecaller, .. }
                if p == pore && basecall.bps.is_none_or(|their| their == bps) =>
            {
                Some((m, khz, tier, basecaller))
            }
            _ => None,
        })
        .min_by_key(|&(_, khz, tier, basecaller)| {
            let ours = version(basecaller).map_or(0, Version::ordinal);
            let distance = match version(basecall.basecaller) {
                Some(theirs) if same_caller(basecaller) => theirs.ordinal().abs_diff(ours),
                // Not known, or another basecaller's: the newest.
                _ => u64::MAX - ours,
            };
            let other_rate = basecall.khz.is_some_and(|their| their != khz);
            (
                other_rate,
                tier_distance(tier),
                std::cmp::Reverse(tier),
                !same_caller(basecaller),
                distance,
                std::cmp::Reverse(ours),
            )
        })
        .map(|(m, khz, tier, basecaller)| {
            let exact = basecall.tier == Some(tier)
                && basecall.basecaller == basecaller
                && basecall.bps.is_some()
                && basecall.khz == Some(khz);
            (m, exact)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ont(basecall_model: &str) -> String {
        format!("@HD\tVN:1.6\n@RG\tID:a\tPL:ONT\tPM:PC48B260\tDS:basecall_model={basecall_model} runid=x\tSM:S\n")
    }

    fn pacbio(instrument: &str, read_type: &str) -> String {
        format!("@RG\tID:a\tPL:PACBIO\tPM:{instrument}\tDS:READTYPE={read_type};BINDINGKIT=103-426-500\tSM:S\n")
    }

    /// A scratch models directory holding `names`, removed when dropped.
    struct ModelsDir(PathBuf);

    impl ModelsDir {
        fn new(test: &str, names: &[&str]) -> Self {
            let dir = std::env::temp_dir().join(format!("lrpipe-models-{test}-{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&dir);
            for name in names {
                std::fs::create_dir_all(dir.join(name)).expect("model directory");
            }
            std::fs::create_dir_all(&dir).expect("models directory");
            Self(dir)
        }

        /// Every model in the table.
        fn all(test: &str) -> Self {
            Self::new(test, &table::MODELS.iter().map(|m| m.name).collect::<Vec<_>>())
        }
    }

    impl Drop for ModelsDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn auto(platform: Platform, dir: &ModelsDir) -> Request<'_> {
        Request { platform, given: None, models: Some(&dir.0), mismatch: Mismatch::Error, fastq: false }
    }

    fn given<'a>(platform: Platform, path: &str, mismatch: Mismatch) -> Request<'a> {
        let given = Some(Given { path: PathBuf::from(path), by: "--model".into() });
        Request { platform, given, models: None, mismatch, fastq: false }
    }

    fn name(choice: &Choice) -> String {
        choice.path.file_name().expect("a name").to_string_lossy().into_owned()
    }

    /// A header naming the reads `model` was trained on.
    fn header_for(model: &Model) -> (Platform, String) {
        match model.reads {
            Reads::Ont { pore, motor, bps, tier, basecaller, .. } => {
                let header = match basecaller {
                    Basecaller::Dorado(v) => ont(&format!("dna_{pore}_{motor}_{bps}bps_{}@{v}", tier.as_str())),
                    // Guppy's names carry no version.
                    _ => ont(&format!("dna_{pore}_{motor}_{bps}bps_{}", tier.as_str())),
                };
                (Platform::Ont, header)
            }
            Reads::Hifi { instruments, .. } => (Platform::Hifi, pacbio(instruments[0], "CCS")),
            Reads::Illumina => unreachable!("not picked"),
        }
    }

    #[test]
    fn every_model_lrpipe_picks_is_picked_for_its_own_reads() {
        let all = ModelsDir::all("every");
        for model in table::MODELS.iter().filter(|m| m.pick == Pick::Yes) {
            let (platform, header) = header_for(model);
            let guppy = matches!(model.reads, Reads::Ont { basecaller: Basecaller::Guppy(_), .. });
            // With every model there, each of dorado's and PacBio's is picked for its reads.
            if !guppy {
                let choice = choose(&header, &auto(platform, &all)).unwrap_or_else(|e| panic!("{}: {e}", model.name));
                assert_eq!(name(&choice), model.name, "{header}");
                assert_eq!(choice.path, all.0.join(model.name));
                assert_eq!(choice.by, "snv.models");
                assert!(choice.warnings.is_empty(), "{}: {:?}", model.name, choice.warnings);
            }
            // Alone, each is picked for its reads; Guppy's, whose names have no version, warned.
            let alone = ModelsDir::new("alone", &[model.name]);
            let choice = choose(&header, &auto(platform, &alone)).unwrap_or_else(|e| panic!("{}: {e}", model.name));
            assert_eq!(name(&choice), model.name, "{header}");
            assert_eq!(choice.warnings.is_empty(), !guppy, "{}: {:?}", model.name, choice.warnings);
        }
    }

    #[test]
    fn the_fixtures_basecall_model_picks_sup_v500() {
        let dir = ModelsDir::new("fixture", table::tests::BUNDLED);
        let choice = choose(&ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), &auto(Platform::Ont, &dir)).expect("chosen");
        assert_eq!(name(&choice), "r1041_e82_400bps_sup_v500");
        assert!(choice.warnings.is_empty(), "{:?}", choice.warnings);
        assert!(choice.reason.contains("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), "{}", choice.reason);
    }

    #[test]
    fn an_unknown_version_gets_the_nearest_model_with_a_warning() {
        let dir = ModelsDir::all("nearest");
        for (basecall, expect) in [
            // Newer than any: the newest of the tier.
            ("dna_r10.4.1_e8.2_400bps_sup@v7.1.0", "r1041_e82_400bps_sup_v520"),
            ("dna_r10.4.1_e8.2_400bps_hac@v6.1.0", "r1041_e82_400bps_hac_v600"),
            // Between two: the nearer, the newer of two as near.
            ("dna_r10.4.1_e8.2_400bps_hac@v4.4.0", "r1041_e82_400bps_hac_v430"),
            ("dna_r10.4.1_e8.2_400bps_hac@v5.0.1", "r1041_e82_400bps_hac_v500"),
            ("dna_r10.4.1_e8.2_400bps_sup@v5.1.0", "r1041_e82_400bps_sup_v520"),
            // Older than any of its sampling rate: the oldest.
            ("dna_r10.4.1_e8.2_400bps_sup@v3.9.0", "r1041_e82_400bps_sup_v400"),
            // No SUP model for v6: the SUP one nearest it, before another tier.
            ("dna_r10.4.1_e8.2_400bps_sup@v6.0.0", "r1041_e82_400bps_sup_v520"),
            // Fast 5 kHz reads: HAC at 5 kHz, before Guppy's 4 kHz fast models.
            ("dna_r10.4.1_e8.2_400bps_fast@v5.0.0", "r1041_e82_400bps_hac_v500"),
            // No version (Guppy's names, 4 kHz): Guppy's models, the newest of the tier.
            ("dna_r10.4.1_e8.2_400bps_sup", "r1041_e82_400bps_sup_g615"),
            ("dna_r10.4.1_e8.2_400bps_hac", "r1041_e82_400bps_hac_g632"),
            // Guppy's 5 kHz: dorado's 5 kHz models.
            ("dna_r10.4.1_e8.2_400bps_5khz_sup", "r1041_e82_400bps_sup_v520"),
            // dorado's R9.4.1 models: Clair3's are Guppy's.
            ("dna_r9.4.1_e8_sup@v3.3", "r941_prom_sup_g5014"),
            ("dna_r9.4.1_e8_hac@v3.3", "r941_prom_hac_g360+g422"),
            // A path to the model.
            ("/opt/dorado/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0", "r1041_e82_400bps_sup_v520"),
        ] {
            let choice =
                choose(&ont(basecall), &auto(Platform::Ont, &dir)).unwrap_or_else(|e| panic!("{basecall}: {e}"));
            assert_eq!(name(&choice), expect, "{basecall}");
            if basecall.starts_with('/') {
                assert!(choice.warnings.is_empty(), "{basecall}: {:?}", choice.warnings);
                continue;
            }
            let [warning] = &choice.warnings[..] else { panic!("{basecall}: one warning, not {:?}", choice.warnings) };
            assert!(warning.starts_with(&format!("using {expect} for ")), "{warning}");
            assert!(warning.contains(basecall) && warning.contains("nearest"), "{warning}");
            // The why is the warning's, once.
            assert!(
                choice.reason.contains("(see the warning)") && !choice.reason.contains("nearest"),
                "{}",
                choice.reason
            );
        }
    }

    #[test]
    fn clair3s_own_models_serve_reads_whose_model_comes_from_rerio() {
        // The models Clair3 installs lack dorado's v4.2.0 and v4.3.0 models (Rerio's): those
        // reads get the nearest it holds, with a warning naming the one to add.
        let dir = ModelsDir::new("bundled", table::tests::BUNDLED);
        for (basecall, expect, missing) in [
            ("dna_r10.4.1_e8.2_400bps_sup@v4.3.0", "r1041_e82_400bps_sup_v500", Some("r1041_e82_400bps_sup_v430")),
            ("dna_r10.4.1_e8.2_400bps_hac@v4.3.0", "r1041_e82_400bps_hac_v500", Some("r1041_e82_400bps_hac_v430")),
            ("dna_r10.4.1_e8.2_400bps_sup@v4.2.0", "r1041_e82_400bps_sup_v500", Some("r1041_e82_400bps_sup_v420")),
            ("dna_r10.4.1_e8.2_400bps_sup@v4.0.0", "r1041_e82_400bps_sup_v410", Some("r1041_e82_400bps_sup_v400")),
            // Guppy's R10.4.1 reads (4 kHz): the 4 kHz dorado model, not 5 kHz v5.2.0.
            ("dna_r10.4.1_e8.2_400bps_sup", "r1041_e82_400bps_sup_v410", Some("r1041_e82_400bps_sup_g615")),
            ("dna_r10.4.1_e8.2_400bps_4khz_hac", "r1041_e82_400bps_hac_v410", Some("r1041_e82_400bps_hac_g632")),
            ("dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "r1041_e82_400bps_sup_v500", None),
        ] {
            let choice =
                choose(&ont(basecall), &auto(Platform::Ont, &dir)).unwrap_or_else(|e| panic!("{basecall}: {e}"));
            assert_eq!(name(&choice), expect, "{basecall}");
            match missing {
                None => assert!(choice.warnings.is_empty(), "{:?}", choice.warnings),
                Some(missing) => {
                    let [warning] = &choice.warnings[..] else { panic!("{basecall}: {:?}", choice.warnings) };
                    assert!(warning.contains(&format!("does not hold {missing}, ")), "{warning}");
                    assert!(warning.contains("ONT's Rerio models converted for Clair3 2"), "{warning}");
                    assert!(
                        warning.contains(&format!("{expect} (")) && warning.contains("nearest it holds"),
                        "{warning}"
                    );
                }
            }
        }
        // None of the chemistry there: the error names the model and where it comes from.
        let e = choose(&ont("dna_r10.4.1_e8.2_260bps_sup@v4.1.0"), &auto(Platform::Ont, &dir)).expect_err("none");
        assert!(matches!(e, ModelError::NotInstalled { name: "r1041_e82_260bps_sup_v410", .. }), "{e}");
        assert!(e.to_string().contains("add r1041_e82_260bps_sup_v410 (one of ONT's Rerio models"), "{e}");
    }

    #[test]
    fn an_unknown_chemistry_or_basecall_model_is_an_error() {
        let dir = ModelsDir::all("unknown");
        for (header, says) in [
            (ont("dna_r10.4_e8.1_sup@v3.5.2"), "no model lrpipe knows"),
            (ont("dna_r10.4.1_e8.2_300bps_sup@v4.0.0"), "no model lrpipe knows"),
            (ont("rna004_130bps_sup@v3.0.1"), "cannot read the basecall model rna004_130bps_sup@v3.0.1"),
            ("@RG\tID:a\tPL:ONT\tSM:S\n".into(), "does not name their basecall model"),
            ("@HD\tVN:1.6\n".into(), "no @RG line"),
        ] {
            let e = choose(&header, &auto(Platform::Ont, &dir)).expect_err("no model");
            assert!(matches!(e, ModelError::Unknown { .. }), "{e}");
            let e = e.to_string();
            assert!(e.contains(says) && e.ends_with("give the model with --model or snv.ont_model"), "{e}");
        }
    }

    #[test]
    fn pacbio_instruments_pick_their_models() {
        let dir = ModelsDir::all("pacbio");
        let hifi = |header: &str| choose(header, &auto(Platform::Hifi, &dir));
        assert_eq!(name(&hifi(&pacbio("REVIO", "CCS")).expect("chosen")), "hifi_revio");
        let sequel = hifi(&pacbio("SEQUELII", "CCS")).expect("chosen");
        assert_eq!(name(&sequel), "hifi_sequel2");
        assert!(sequel.warnings.is_empty());
        // The first Sequel: Sequel II's model, the nearest, with a warning.
        let first = hifi(&pacbio("SEQUEL", "CCS")).expect("chosen");
        assert_eq!(name(&first), "hifi_sequel2");
        assert_eq!(first.warnings.len(), 1, "{:?}", first.warnings);
        // An instrument Clair3 has no model for, or none named: the newest, with a warning.
        for header in [pacbio("VEGA", "CCS"), "@RG\tID:a\tPL:PACBIO\tSM:S\n".into(), String::new()] {
            let choice = hifi(&header).expect("chosen");
            assert_eq!(name(&choice), "hifi_revio");
            assert_eq!(choice.warnings.len(), 1, "{header}: {:?}", choice.warnings);
            assert!(choice.warnings[0].starts_with("using hifi_revio for "), "{:?}", choice.warnings);
        }
        // Revio's model not there: Sequel II's, warned.
        let sequel_only = ModelsDir::new("sequel-only", &["hifi_sequel2"]);
        let choice = choose(&pacbio("REVIO", "CCS"), &auto(Platform::Hifi, &sequel_only)).expect("chosen");
        assert_eq!(name(&choice), "hifi_sequel2");
        assert!(
            choice.warnings[0].contains("does not hold hifi_revio, Clair3's model for them"),
            "{:?}",
            choice.warnings
        );
        // Subreads are not HiFi; pbmm2's UNKNOWN (reads aligned from FASTQ) may be.
        let e = hifi(&pacbio("SEQUELII", "SUBREAD")).expect_err("not HiFi").to_string();
        assert!(e.contains("not HiFi (READTYPE=CCS) reads") && e.contains("--model-mismatch warn"), "{e}");
        let unknown = hifi(&pacbio("REVIO", "UNKNOWN")).expect("taken as HiFi");
        assert_eq!(name(&unknown), "hifi_revio");
        let [warning] = &unknown.warnings[..] else { panic!("{:?}", unknown.warnings) };
        assert!(warning.contains("READTYPE=UNKNOWN") && warning.ends_with("taken as HiFi reads"), "{warning}");
        assert!(hifi("@RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\n").is_ok());
    }

    #[test]
    fn reads_of_the_other_platform_are_refused_and_of_another_warned_of() {
        let dir = ModelsDir::all("platform");
        let e = choose(&pacbio("REVIO", "CCS"), &auto(Platform::Ont, &dir)).expect_err("PacBio reads").to_string();
        assert!(e.contains("PacBio reads (@RG PM:REVIO, READTYPE=CCS), not ont reads as --platform says"), "{e}");
        let e = choose(&ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), &auto(Platform::Hifi, &dir)).expect_err("ONT");
        assert!(matches!(e, ModelError::Mismatch(_)), "{e}");
        // A PL lrpipe does not call, stamped by an aligner's -R: taken as --platform says.
        let illumina = "@RG\tID:a\tPL:ILLUMINA\tSM:S\n";
        let choice = choose(illumina, &auto(Platform::Hifi, &dir)).expect("taken as HiFi");
        assert_eq!(name(&choice), "hifi_revio");
        assert!(
            choice.warnings[0].contains("ILLUMINA reads (@RG PL:ILLUMINA), which lrpipe does not call"),
            "{:?}",
            choice.warnings
        );
        let choice = choose(illumina, &given(Platform::Ont, "/m/r1041_e82_400bps_sup_v500", Mismatch::Error))
            .expect("a model given");
        assert!(choice.warnings.iter().any(|w| w.contains("taken as ont reads")), "{:?}", choice.warnings);
        // A basecall model says more than such a PL.
        let stamped = ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0").replace("PL:ONT", "PL:ILLUMINA");
        let choice = choose(&stamped, &auto(Platform::Ont, &dir)).expect("ONT by its basecall model");
        assert_eq!(name(&choice), "r1041_e82_400bps_sup_v500");
        assert!(choice.warnings.is_empty(), "{:?}", choice.warnings);
    }

    #[test]
    fn several_basecall_models_need_a_given_model() {
        let dir = ModelsDir::all("several");
        let second = |model: &str| ont(model).replace("ID:a", "ID:b").replace("@HD\tVN:1.6\n", "");
        let two =
            format!("{}{}", ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), second("dna_r10.4.1_e8.2_400bps_sup@v5.2.0"));
        let e = choose(&two, &auto(Platform::Ont, &dir)).expect_err("two basecall models");
        assert!(matches!(e, ModelError::Several(ref s) if s.len() == 2), "{e}");
        assert!(e.to_string().contains("sup@v5.0.0") && e.to_string().contains("sup@v5.2.0"), "{e}");
        // Given, it is checked against each.
        let choice =
            choose(&two, &given(Platform::Ont, "/m/r1041_e82_400bps_sup_v520", Mismatch::Error)).expect("given");
        assert_eq!(choice.warnings.len(), 1, "{:?}", choice.warnings);
        assert!(choice.warnings[0].contains("r1041_e82_400bps_sup_v500"), "{:?}", choice.warnings);
        // Read groups basecalled alike are one.
        let same =
            format!("{}{}", ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), second("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"));
        let choice = choose(&same, &auto(Platform::Ont, &dir)).expect("one model");
        assert_eq!(name(&choice), "r1041_e82_400bps_sup_v500");
        assert!(choice.warnings.is_empty(), "{:?}", choice.warnings);
        // A read group that does not say (`minimap2 -R '@RG\tID:x\tPL:ONT'`) is taken as the others.
        let bare = format!("{}@RG\tID:x\tPL:ONT\tSM:S\n", ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"));
        let choice = choose(&bare, &auto(Platform::Ont, &dir)).expect("one model");
        assert_eq!(name(&choice), "r1041_e82_400bps_sup_v500");
        let [warning] = &choice.warnings[..] else { panic!("{:?}", choice.warnings) };
        assert!(warning.contains("does not say what made its reads; taken to be made as the others"), "{warning}");
    }

    #[test]
    fn a_given_model_the_header_contradicts_is_not_said_to_be_checked() {
        // A SUP model for HAC reads: warned of, and its reason says so rather than "checked
        // against the header".
        let header = ont("dna_r10.4.1_e8.2_400bps_hac@v4.3.0");
        let choice =
            choose(&header, &given(Platform::Ont, "/m/r1041_e82_400bps_sup_v500", Mismatch::Error)).expect("warned");
        assert_eq!(choice.warnings.len(), 1, "{:?}", choice.warnings);
        assert!(
            choice.reason.ends_with("which the header says the reads are not (see the warning)"),
            "{}",
            choice.reason
        );
        let header = ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
        let choice = choose(&header, &given(Platform::Ont, "/m/r1041_e82_400bps_sup_v500", Mismatch::Error))
            .expect("the reads' own");
        assert!(choice.reason.ends_with(", checked against the header"), "{}", choice.reason);
    }

    #[test]
    fn a_models_directory_that_is_missing_or_one_model_is_said_to_be() {
        let header = ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
        let gone = Path::new("/nonexistent/models");
        let request = Request {
            platform: Platform::Ont,
            given: None,
            models: Some(gone),
            mismatch: Mismatch::Error,
            fastq: false,
        };
        let e = choose(&header, &request).expect_err("no directory").to_string();
        assert_eq!(
            e,
            "the site configuration's [snv] models /nonexistent/models does not exist: give the directory of Clair3's models"
        );
        let dir = ModelsDir::new("one-model", &["r1041_e82_400bps_sup_v500"]);
        let model = dir.0.join("r1041_e82_400bps_sup_v500");
        std::fs::write(model.join("pileup.pt"), "pileup").expect("written");
        let request = Request {
            platform: Platform::Ont,
            given: None,
            models: Some(&model),
            mismatch: Mismatch::Error,
            fastq: false,
        };
        let e = choose(&header, &request).expect_err("one model").to_string();
        assert!(
            e.contains("is a model's own directory (it holds pileup.pt): give the directory that holds the models"),
            "{e}"
        );
        assert!(e.contains(&format!("({})", dir.0.display())) && e.contains("[snv] ont_model"), "{e}");
        let request = Request {
            platform: Platform::Ont,
            given: None,
            models: Some(&dir.0),
            mismatch: Mismatch::Error,
            fastq: false,
        };
        assert_eq!(name(&choose(&header, &request).expect("its parent")), "r1041_e82_400bps_sup_v500");
    }

    #[test]
    fn a_given_model_for_other_reads_is_refused() {
        let header = ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
        for (model, says) in [
            ("/m/hifi_revio", "hifi_revio is Clair3's model for PacBio HiFi reads from REVIO, not ont reads"),
            ("/m/ilmn", "ilmn is Clair3's model for Illumina reads, not ont reads"),
            ("/m/r941_prom_sup_g5014", "r941_prom_sup_g5014 is Clair3's model for Oxford Nanopore R9.4.1 E8 450 bps"),
            ("/m/ont", "but these are Oxford Nanopore reads basecalled with dna_r10.4.1_e8.2_400bps_sup@v5.0.0"),
        ] {
            let e = choose(&header, &given(Platform::Ont, model, Mismatch::Error)).expect_err(model);
            assert!(matches!(e, ModelError::Mismatch(_)), "{e}");
            let e = e.to_string();
            assert!(e.contains(says) && e.ends_with("or --model-mismatch warn to run anyway"), "{e}");
            // With --model-mismatch warn, the run goes on with the model, warned.
            let choice = choose(&header, &given(Platform::Ont, model, Mismatch::Warn)).expect("warned");
            assert_eq!(choice.path, PathBuf::from(model));
            assert!(
                choice.warnings.iter().any(|w| w.contains(says) && w.ends_with("(--model-mismatch warn)")),
                "{:?}",
                choice.warnings
            );
        }
        // A HiFi model for a HiFi run of ONT reads: the reads are not what --platform says.
        let e = choose(&header, &given(Platform::Hifi, "/m/hifi_revio", Mismatch::Error)).expect_err("ONT reads");
        assert!(e.to_string().contains("not hifi reads as --platform says"), "{e}");
        // A model of the right chemistry for another basecall model, or instrument: a warning.
        let choice = choose(&header, &given(Platform::Ont, "/m/r1041_e82_400bps_hac_v520", Mismatch::Error))
            .expect("same chemistry");
        let [warning] = &choice.warnings[..] else { panic!("{:?}", choice.warnings) };
        assert!(warning.contains("(Clair3's model for them: r1041_e82_400bps_sup_v500)"), "{warning}");
        let choice =
            choose(&pacbio("REVIO", "CCS"), &given(Platform::Hifi, "/m/hifi_sequel2", Mismatch::Error)).expect("HiFi");
        assert_eq!(choice.warnings.len(), 1, "{:?}", choice.warnings);
        // The model for the reads: nothing to say.
        let choice =
            choose(&header, &given(Platform::Ont, "/m/r1041_e82_400bps_sup_v500", Mismatch::Error)).expect("right");
        assert!(choice.warnings.is_empty(), "{:?}", choice.warnings);
        // The reason follows `by` (`--model: Clair3's model for ...`), so it does not repeat it.
        assert!(choice.reason.starts_with("Clair3's model for"), "{}", choice.reason);
        // Nothing in the header to check it against: said.
        let choice = choose("@HD\tVN:1.6\n", &given(Platform::Ont, "/m/r1041_e82_400bps_sup_v500", Mismatch::Error))
            .expect("unchecked");
        let [warning] = &choice.warnings[..] else { panic!("{:?}", choice.warnings) };
        assert!(warning.contains("the header does not say what made the reads"), "{warning}");
    }

    #[test]
    fn reads_of_another_platform_let_through_are_not_said_to_lack_a_header() {
        // `--platform hifi` on ONT reads with --model-mismatch warn: the header says what made
        // them, so it is not said to lack the instrument or the basecall model.
        let dir = ModelsDir::all("let-through");
        let header = ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
        let request = Request {
            platform: Platform::Hifi,
            given: None,
            models: Some(&dir.0),
            mismatch: Mismatch::Warn,
            fastq: false,
        };
        let choice = choose(&header, &request).expect("let through");
        let using = choice.warnings.iter().find(|w| w.starts_with("using ")).expect("a warning of the model");
        assert!(!using.contains("does not name the instrument"), "{using}");
        assert!(using.contains("these are not PacBio reads, so the header names no PacBio instrument"), "{using}");
        // hifi_revio is no nearer to ONT reads than another HiFi model.
        assert!(!using.contains("nearest") && using.contains("the HiFi model lrpipe takes"), "{using}");
        let choice = choose(&header, &given(Platform::Hifi, "/m/hifi_revio", Mismatch::Warn)).expect("let through");
        assert!(
            choice.warnings.iter().all(|w| !w.contains("does not say what made the reads")),
            "{:?}",
            choice.warnings
        );
    }

    #[test]
    fn a_move_table_model_is_refused_whatever_the_mismatch_setting() {
        // lrpipe snv refuses one without --enable_dwell_time, after the other stages started.
        let header = ont("dna_r10.4.1_e8.2_400bps_sup@v5.2.0");
        for model in ["/m/r1041_e82_400bps_sup_v520_with_mv", "/m/my_model_WITH_MV"] {
            for mismatch in [Mismatch::Error, Mismatch::Warn] {
                let e = choose(&header, &given(Platform::Ont, model, mismatch)).expect_err(model);
                assert!(matches!(e, ModelError::MoveTable(_)), "{e}");
                assert!(e.to_string().contains("--enable_dwell_time, which lrpipe run does not pass"), "{e}");
            }
        }
    }

    #[test]
    fn a_given_model_lrpipe_does_not_know_is_warned_of_not_refused() {
        let header = ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
        let choice = choose(&header, &given(Platform::Ont, "/m/my_finetuned_model", Mismatch::Error)).expect("custom");
        assert_eq!(choice.path, PathBuf::from("/m/my_finetuned_model"));
        let [warning] = &choice.warnings[..] else { panic!("{:?}", choice.warnings) };
        assert!(warning.contains("not one of Clair3's models by name"), "{warning}");
        // Still the reads must be the platform's.
        let e =
            choose(&pacbio("REVIO", "CCS"), &given(Platform::Ont, "/m/custom", Mismatch::Error)).expect_err("PacBio");
        assert!(matches!(e, ModelError::Mismatch(_)), "{e}");
        // A model lrpipe never picks, given: a warning saying why.
        let finetuned = given(Platform::Ont, "/m/r1041_e82_400bps_sup_v430_bacteria_finetuned", Mismatch::Error);
        let choice = choose(&header, &finetuned).expect("given");
        assert!(choice.warnings.iter().any(|w| w.contains("fine-tuned on bacterial genomes")), "{:?}", choice.warnings);
    }

    #[test]
    fn a_chosen_model_must_be_installed() {
        let dir = ModelsDir::new("none", &[]);
        let e = choose(&ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), &auto(Platform::Ont, &dir)).expect_err("not there");
        assert!(matches!(e, ModelError::NotInstalled { name: "r1041_e82_400bps_sup_v500", .. }), "{e}");
        assert!(e.to_string().contains("(one of Clair3's own models"), "{e}");
        let request = Request { models: None, ..auto(Platform::Ont, &dir) };
        let e = choose(&ont("dna_r10.4.1_e8.2_400bps_sup@v5.0.0"), &request).expect_err("no directory");
        assert!(matches!(e, ModelError::NoModel("ont")), "{e}");
    }

    #[test]
    fn the_header_is_read_from_bam_and_cram() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../lrpipe-snv/ci/fixture");
        for reads in ["ont.bam", "ont.cram"] {
            let text = header_text(&fixture.join(reads)).unwrap_or_else(|e| panic!("{reads}: {e}"));
            let [Sequenced::Ont(Some(b))] = &sequenced(&text)[..] else { panic!("{reads}: {text}") };
            assert_eq!(b.name, "dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "{reads}");
        }
    }

    #[test]
    fn basecall_models_parse() {
        let b = BasecallModel::parse("dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
        assert_eq!(b.pore.as_deref(), Some("r10.4.1"));
        assert_eq!((b.bps, b.khz, b.tier), (Some(400), Some(5), Some(Tier::Sup)));
        assert_eq!(b.basecaller, Basecaller::Dorado(Version([5, 0, 0])));
        assert_eq!(BasecallModel::parse("dna_r10.4.1_e8.2_400bps_hac@v4.1.0").khz, Some(4));
        let guppy = BasecallModel::parse("dna_r9.4.1_450bps_hac");
        assert_eq!((guppy.pore.as_deref(), guppy.bps, guppy.khz), (Some("r9.4.1"), Some(450), Some(4)));
        assert_eq!(guppy.basecaller, Basecaller::Unknown);
        assert_eq!(BasecallModel::parse("dna_r10.4.1_e8.2_400bps_5khz_sup").khz, Some(5));
        let rna = BasecallModel::parse("rna004_130bps_sup@v3.0.1");
        assert_eq!((rna.pore, rna.tier, rna.khz), (None, None, None));
        // The fixture's read group, with its modified-base models after the basecall model.
        let rg = "@RG\tID:x\tPU:PAY87954\tPM:PC48B260\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
                  modbase_models=dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mC_5hmC@v2.0.1 runid=70eb\tLB:hg003\tSM:hg003";
        assert_eq!(Sequenced::of(rg), Some(Sequenced::Ont(Some(b))));
        assert_eq!(Sequenced::of("@RG\tID:x\tSM:s"), None);
    }
}

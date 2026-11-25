//! The versioned manifests of what `lrpipe setup` fetches: `lrpipe/resources/*.toml`, built into
//! the binary, one per genome, each read with `models.toml` (Clair3's models) after it.
//!
//! A [`FileEntry`] is fetched from its URL to its path under the setup directory and checked
//! against its SHA-256 and the publisher's MD5 where the manifest has them. One with neither says
//! why in `unpinned`: it is fetched and used all the same, checked by its size alone and never
//! counted as verified, and `lrpipe setup --pin` writes its SHA-256 into the manifest
//! ([`pin_text`]).

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Component, Path};

use anyhow::{Context, Result, bail, ensure};
use serde::Deserialize;

use crate::plan::Platform;

/// The format of the manifests this lrpipe reads.
pub const FORMAT: u32 = 1;

/// Clair3's models, read after each genome's manifest.
const MODELS: &str = include_str!("../../resources/models.toml");

/// A genome `lrpipe setup` has a manifest for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Genome {
    /// GRCh38's no-alt analysis set.
    Grch38,
    /// T2T-CHM13 v2.0's analysis set.
    Chm13,
    /// GRCh37 as hs37d5.
    Grch37,
}

impl Genome {
    /// Its name in manifests and on the command line.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Grch38 => "GRCh38",
            Self::Chm13 => "CHM13",
            Self::Grch37 => "GRCh37",
        }
    }

    /// The manifest built into lrpipe for it, with Clair3's models.
    #[must_use]
    pub fn manifest_text(self) -> String {
        let genome = match self {
            Self::Grch38 => include_str!("../../resources/grch38.toml"),
            Self::Chm13 => include_str!("../../resources/chm13.toml"),
            Self::Grch37 => include_str!("../../resources/grch37.toml"),
        };
        with_models(genome)
    }
}

/// A genome's manifest followed by Clair3's models, as one manifest.
#[must_use]
pub fn with_models(genome: &str) -> String {
    format!("{genome}\n{MODELS}")
}

/// A manifest: the files to fetch, what is built from them, and the site configuration they make.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    /// The manifest format ([`FORMAT`]).
    pub format: Option<u32>,
    /// The genome, as `--genome` names it.
    pub genome: Option<String>,
    /// Resources this genome has none of, said when it is set up.
    #[serde(default)]
    pub missing: Vec<String>,
    /// Which file or build each site configuration key takes.
    #[serde(default)]
    pub site: SiteKeys,
    /// The model written to the site configuration for each platform.
    #[serde(default)]
    pub models: Models,
    /// The files to fetch.
    #[serde(default, rename = "file")]
    pub files: Vec<FileEntry>,
    /// What is built from them.
    #[serde(default, rename = "build")]
    pub builds: Vec<Build>,
}

/// The site configuration's keys, each naming the [`FileEntry`] or [`Build`] (by id) that fills it.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SiteKeys {
    /// The reference FASTA.
    pub reference: Option<String>,
    /// `[sv] tandem_repeats`.
    pub tandem_repeats: Option<String>,
    /// `[repeats] catalog`.
    pub catalog: Option<String>,
    /// `[cnv] blacklist`.
    pub blacklist: Option<String>,
    /// `[cnv] metadata`: a path under the setup directory, made from the reference by
    /// lrpipe-cnv.
    pub metadata: Option<String>,
}

/// The model each platform's runs use.
#[derive(Clone, Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Models {
    /// `[snv] ont_model`.
    pub ont: Option<String>,
    /// `[snv] hifi_model`.
    pub hifi: Option<String>,
}

impl Models {
    /// The model for `platform`.
    #[must_use]
    pub fn of(&self, platform: Platform) -> Option<&str> {
        match platform {
            Platform::Ont => self.ont.as_deref(),
            Platform::Hifi => self.hifi.as_deref(),
        }
    }
}

/// How a fetched file is unpacked.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Unpack {
    /// gzip or BGZF, decompressed to the file's path.
    Gunzip,
}

/// One file to fetch.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FileEntry {
    /// Its name in the manifest.
    pub id: String,
    /// Where it is published: `https://`, `http://` or `file://`.
    pub url: String,
    /// Where it goes, relative to the setup directory.
    pub path: String,
    /// How it is unpacked to `path`, having been fetched to [`FileEntry::download`].
    pub unpack: Option<Unpack>,
    /// Bytes fetched.
    pub size: Option<u64>,
    /// SHA-256 of the bytes fetched (lowercase hex).
    pub sha256: Option<String>,
    /// The publisher's MD5 of the bytes fetched (lowercase hex).
    pub md5: Option<String>,
    /// Why it has no checksum, where it has none.
    pub unpinned: Option<String>,
    /// Its licence.
    pub licence: String,
    /// Where the licence is.
    pub licence_url: Option<String>,
    /// Whether its terms must be accepted (`--accept-licences`) before it is fetched.
    #[serde(default)]
    pub accept: bool,
    /// The terms to accept, as printed.
    pub terms: Option<String>,
    /// The only platform that needs it (a model's files); every platform's otherwise.
    pub platform: Option<Platform>,
    /// The Clair3 model whose directory it is in.
    pub model: Option<String>,
}

impl FileEntry {
    /// Where it is fetched to, relative to the setup directory: its path, or for a file
    /// unpacked there, `sources/` and the path's file name with `.gz`.
    #[must_use]
    pub fn download(&self) -> String {
        match self.unpack {
            None => self.path.clone(),
            Some(Unpack::Gunzip) => {
                let name = Path::new(&self.path).file_name().map_or_else(String::new, |n| n.to_string_lossy().into());
                format!("sources/{name}.gz")
            }
        }
    }

    /// Whether the manifest has a checksum of it.
    #[must_use]
    pub fn pinned(&self) -> bool {
        self.sha256.is_some() || self.md5.is_some()
    }

    /// Whether a setup for `platforms` needs it.
    #[must_use]
    pub fn needed_for(&self, platforms: &[Platform]) -> bool {
        self.platform.is_none_or(|p| platforms.contains(&p))
    }

    /// Its checksums, as `ALGORITHM:HEX`, in a fixed order: what a file was verified against.
    #[must_use]
    pub fn expected(&self) -> String {
        let sums: Vec<String> = [("sha256", &self.sha256), ("md5", &self.md5)]
            .into_iter()
            .filter_map(|(algorithm, hex)| Some(format!("{algorithm}:{}", hex.as_ref()?)))
            .collect();
        sums.join(" ")
    }
}

/// What is built from fetched files.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BuildKind {
    /// A repeat catalogue: the inputs' BED lines (gzip, BGZF or plain) together, sorted by
    /// contig name and start as `LC_ALL=C sort -k1,1 -k2,2n` sorts them.
    Catalog,
}

/// A file built from fetched ones.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Build {
    /// Its name in the manifest.
    pub id: String,
    /// How it is built.
    pub kind: BuildKind,
    /// The files (by id) it is built from.
    pub inputs: Vec<String>,
    /// Where it goes, relative to the setup directory.
    pub path: String,
    /// SHA-256 it must have, where known.
    pub sha256: Option<String>,
}

impl Manifest {
    /// Parses and checks a manifest.
    ///
    /// # Errors
    ///
    /// If `text` is not a manifest of format [`FORMAT`]: unknown keys, a file with neither a
    /// checksum nor a reason it has none (or with both), a path outside the setup directory, a
    /// malformed checksum, ids used twice, or a site key, build input or model naming nothing.
    pub fn parse(text: &str) -> Result<Self> {
        let manifest: Self = toml::from_str(text).context("parsing the manifest")?;
        manifest.check()?;
        Ok(manifest)
    }

    /// The file with this id.
    #[must_use]
    pub fn file(&self, id: &str) -> Option<&FileEntry> {
        self.files.iter().find(|f| f.id == id)
    }

    /// The build with this id.
    #[must_use]
    pub fn build(&self, id: &str) -> Option<&Build> {
        self.builds.iter().find(|b| b.id == id)
    }

    /// The path, relative to the setup directory, of the file or build with this id.
    #[must_use]
    pub fn path_of(&self, id: &str) -> Option<&str> {
        self.file(id).map(|f| f.path.as_str()).or_else(|| self.build(id).map(|b| b.path.as_str()))
    }

    /// The files of `model`.
    pub fn model_files<'a>(&'a self, model: &'a str) -> impl Iterator<Item = &'a FileEntry> {
        self.files.iter().filter(move |f| f.model.as_deref() == Some(model))
    }

    /// Checks what [`Manifest::parse`] promises, and that a genome's manifest names the genome
    /// and its reference.
    fn check(&self) -> Result<()> {
        if let Some(format) = self.format {
            ensure!(format == FORMAT, "manifest format {format}: this lrpipe reads format {FORMAT}");
        }
        let mut ids = BTreeSet::new();
        for f in &self.files {
            ensure!(ids.insert(f.id.as_str()), "the id {} is used twice", f.id);
            relative(&f.path).with_context(|| format!("file {}", f.id))?;
            ensure!(
                ["https://", "http://", "file://"].iter().any(|s| f.url.starts_with(s)),
                "file {}: {} is not an https://, http:// or file:// URL",
                f.id,
                f.url
            );
            for (name, hex, len) in [("sha256", &f.sha256, 64), ("md5", &f.md5, 32)] {
                if let Some(hex) = hex {
                    hex_digest(hex, len).with_context(|| format!("file {}: {name}", f.id))?;
                }
            }
            match (f.pinned(), &f.unpinned) {
                (false, None) => bail!("file {} has no sha256 or md5, and no `unpinned` saying why", f.id),
                (true, Some(_)) if f.sha256.is_some() => bail!("file {} has a sha256 and is marked unpinned", f.id),
                (true, Some(_)) => bail!("file {} has an md5 and is marked unpinned", f.id),
                _ => {}
            }
            ensure!(!f.accept || f.terms.is_some(), "file {}: terms to accept, but no `terms`", f.id);
            if let Some(model) = &f.model {
                ensure!(f.platform.is_some(), "file {} of model {model} has no platform", f.id);
                ensure!(
                    Path::new(&f.path).parent() == Some(&Path::new("models").join(model)),
                    "file {}: a file of model {model} goes in models/{model}/",
                    f.id
                );
            }
        }
        for b in &self.builds {
            ensure!(ids.insert(b.id.as_str()), "the id {} is used twice", b.id);
            relative(&b.path).with_context(|| format!("build {}", b.id))?;
            ensure!(!b.inputs.is_empty(), "build {} has no inputs", b.id);
            for input in &b.inputs {
                ensure!(self.file(input).is_some(), "build {}: no file {input}", b.id);
            }
            if let Some(hex) = &b.sha256 {
                hex_digest(hex, 64).with_context(|| format!("build {}: sha256", b.id))?;
            }
        }
        let keys = [
            ("reference", &self.site.reference),
            ("tandem_repeats", &self.site.tandem_repeats),
            ("catalog", &self.site.catalog),
            ("blacklist", &self.site.blacklist),
        ];
        for (key, id) in keys {
            if let Some(id) = id {
                ensure!(self.path_of(id).is_some(), "[site] {key}: no file or build {id}");
            }
        }
        if let Some(reference) = &self.site.reference {
            ensure!(self.file(reference).is_some(), "[site] reference: {reference} is not a file");
        }
        if let Some(mdr) = &self.site.metadata {
            relative(mdr).context("[site] metadata")?;
        }
        for platform in [Platform::Ont, Platform::Hifi] {
            if let Some(model) = self.models.of(platform) {
                let files: Vec<&FileEntry> = self.model_files(model).collect();
                ensure!(
                    !files.is_empty() && files.iter().all(|f| f.platform == Some(platform)),
                    "[models] {}: no {} model {model}",
                    platform.as_str(),
                    platform.as_str()
                );
            }
        }
        if self.genome.is_some() {
            ensure!(self.format.is_some(), "the manifest has no format");
            ensure!(self.site.reference.is_some(), "[site] has no reference");
        }
        self.check_paths()
    }

    /// Checks that no two things setup writes are at one path, and none is inside another or in
    /// `.lrpipe-setup/`: each file's path and download (and their `.part` files and their
    /// versions), each build's path, the reference's `.fai` and `.gzi`, the `.mdr`, and
    /// `site.toml` with its temporary and previous copies.
    fn check_paths(&self) -> Result<()> {
        let mut claims: BTreeMap<String, String> = BTreeMap::new();
        let mut claim = |path: String, what: String| -> Result<()> {
            for p in [path.clone(), format!("{path}.part"), format!("{path}.part.validator")] {
                if let Some(other) = claims.insert(p.clone(), what.clone()) {
                    bail!("{what} and {other} are both written at {p}");
                }
            }
            Ok(())
        };
        for name in ["site.toml", "site.toml.tmp", "site.toml.previous"] {
            claim(name.into(), "lrpipe setup's site configuration".into())?;
        }
        for f in &self.files {
            claim(f.path.clone(), format!("file {}", f.id))?;
            if f.unpack.is_some() {
                claim(f.download(), format!("file {}'s download", f.id))?;
            }
        }
        for b in &self.builds {
            claim(b.path.clone(), format!("build {}", b.id))?;
        }
        if let Some(reference) = self.site.reference.as_deref().and_then(|r| self.file(r)) {
            claim(format!("{}.fai", reference.path), "the reference's .fai".into())?;
            claim(format!("{}.gzi", reference.path), "the reference's .gzi".into())?;
        }
        if let Some(mdr) = &self.site.metadata {
            claim(mdr.clone(), "[site] metadata".into())?;
        }
        for (path, what) in &claims {
            ensure!(!path.starts_with(".lrpipe-setup"), "{what}: {path} is in lrpipe setup's own .lrpipe-setup/");
            for dir in Path::new(path).ancestors().skip(1).filter(|d| !d.as_os_str().is_empty()) {
                if let Some(other) = claims.get(&*dir.to_string_lossy()) {
                    bail!("{what} is written inside {}, where {other} is written", dir.display());
                }
            }
        }
        Ok(())
    }
}

/// Checks that `path` stays inside the setup directory.
fn relative(path: &str) -> Result<()> {
    let p = Path::new(path);
    ensure!(
        !path.is_empty() && p.components().all(|c| matches!(c, Component::Normal(_))),
        "{path} is not a path inside the setup directory"
    );
    Ok(())
}

/// Checks that `hex` is `len` lowercase hexadecimal digits.
fn hex_digest(hex: &str, len: usize) -> Result<()> {
    ensure!(
        hex.len() == len && hex.bytes().all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)),
        "{hex} is not {len} lowercase hexadecimal digits"
    );
    Ok(())
}

/// `text`, a manifest, with `sha256` set for each `(id, sha256)` given: written after the
/// file's `url`, replacing any it had, and with its `unpinned` line removed. Comments and
/// layout are kept.
///
/// # Errors
///
/// If a file is not in the manifest, or the result does not parse with those checksums.
pub fn pin_text(text: &str, sums: &[(String, String)]) -> Result<String> {
    let mut out: Vec<String> = Vec::new();
    // The block being read: its id and the checksum it gets, if it is pinned.
    let mut pinning: Option<(String, String)> = None;
    let mut done = BTreeSet::new();
    for line in text.lines() {
        let trimmed = line.trim_start();
        if trimmed.starts_with('[') {
            pinning = None;
        }
        if let Some(id) = value_of(trimmed, "id") {
            pinning = sums.iter().find(|(i, _)| *i == id).cloned();
        }
        if pinning.is_some() && (value_of(trimmed, "sha256").is_some() || value_of(trimmed, "unpinned").is_some()) {
            continue;
        }
        out.push(line.to_owned());
        if value_of(trimmed, "url").is_some()
            && let Some((id, sum)) = &pinning
            && done.insert(id.clone())
        {
            out.push(format!("sha256 = \"{sum}\""));
        }
    }
    for (id, _) in sums {
        ensure!(done.contains(id), "no file {id} with a url in the manifest");
    }
    let mut pinned = out.join("\n");
    if text.ends_with('\n') {
        pinned.push('\n');
    }
    let parsed: Manifest = toml::from_str(&pinned).context("the pinned manifest does not parse")?;
    parsed.check().context("the pinned manifest")?;
    for (id, sum) in sums {
        ensure!(parsed.file(id).and_then(|f| f.sha256.as_deref()) == Some(sum.as_str()), "pinning {id}");
    }
    Ok(pinned)
}

/// The string value of `key = "..."` on a manifest line, if that is what the line holds.
fn value_of(line: &str, key: &str) -> Option<String> {
    let rest = line.strip_prefix(key)?.trim_start().strip_prefix('=')?.trim();
    let value: toml::Value =
        toml::from_str(&format!("v = {rest}")).ok().and_then(|t: toml::Table| t.get("v").cloned())?;
    value.as_str().map(str::to_owned)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SHA: &str = "0000000000000000000000000000000000000000000000000000000000000000";

    fn file(id: &str, checksum: &str) -> String {
        format!(
            "[[file]]\nid = \"{id}\"\nurl = \"file:///x/{id}\"\npath = \"{id}.bed\"\n{checksum}\nlicence = \"MIT\"\n"
        )
    }

    #[test]
    fn the_builtin_manifests_parse() {
        for genome in [Genome::Grch38, Genome::Chm13, Genome::Grch37] {
            let manifest =
                Manifest::parse(&genome.manifest_text()).unwrap_or_else(|e| panic!("{}: {e:#}", genome.as_str()));
            assert_eq!(manifest.genome.as_deref(), Some(genome.as_str()));
            assert!(manifest.file("reference").is_some_and(|f| f.unpack == Some(Unpack::Gunzip)));
            for platform in [Platform::Ont, Platform::Hifi] {
                let model = manifest.models.of(platform).expect("a model for each platform");
                let files: Vec<&str> = manifest.model_files(model).map(|f| f.path.as_str()).collect();
                assert_eq!(files, [format!("models/{model}/pileup.pt"), format!("models/{model}/full_alignment.pt")]);
            }
        }
        let grch38 = Manifest::parse(&Genome::Grch38.manifest_text()).expect("GRCh38");
        assert_eq!(grch38.path_of("catalog"), Some("catalog.adotto_strchive.hg38.bed"));
        assert_eq!(
            grch38.file("reference").map(FileEntry::download).as_deref(),
            Some("sources/GRCh38_no_alt_analysis_set.fa.gz")
        );
        // ONT's models are licensed for research only, which setup records, but fetches them without
        // asking: the licence permits redistributing them for research.
        let ont = grch38.file("r1041_e82_400bps_sup_v500/pileup").expect("the ONT model");
        assert!(!ont.accept && ont.platform == Some(Platform::Ont) && !ont.pinned());
    }

    #[test]
    fn every_trf_bed_from_pbsv_needs_its_terms_accepted() {
        // pbsv's annotations/ is their origin, whose only licence is PacBio's (Sniffles'
        // annotations/source.txt says so for GRCh38's and hs37d5's): each is gated alike.
        let mut terms = BTreeSet::new();
        for genome in [Genome::Grch38, Genome::Chm13, Genome::Grch37] {
            let manifest =
                Manifest::parse(&genome.manifest_text()).unwrap_or_else(|e| panic!("{}: {e:#}", genome.as_str()));
            let trf = manifest.site.tandem_repeats.as_deref().and_then(|id| manifest.file(id)).expect("a TRF BED");
            assert!(trf.accept && trf.licence.starts_with("PacBio Software License Agreement"), "{}", genome.as_str());
            terms.insert(trf.terms.clone());
        }
        assert_eq!(terms.len(), 1, "the same terms for each");
    }

    #[test]
    fn the_readmes_setup_commands_have_the_flags_they_need() {
        // Each `lrpipe setup --genome G --platform P ...` in lrpipe's README runs as written:
        // with --accept-licences where a file it fetches has terms to accept, and no flag
        // lrpipe no longer has.
        let readme = include_str!("../../README.md");
        let mut seen = 0;
        for line in readme.lines().filter(|l| l.starts_with("lrpipe setup --genome")) {
            let words: Vec<&str> = line.split_whitespace().collect();
            let value = |flag: &str| words.iter().position(|w| *w == flag).map(|i| words[i + 1]);
            let genome = match value("--genome").expect("a genome") {
                "GRCh38" => Genome::Grch38,
                "CHM13" => Genome::Chm13,
                "GRCh37" => Genome::Grch37,
                other => panic!("{other}"),
            };
            let platforms: Vec<Platform> = value("--platform")
                .expect("platforms")
                .split(',')
                .map(|p| if p == "ont" { Platform::Ont } else { Platform::Hifi })
                .collect();
            let manifest = Manifest::parse(&genome.manifest_text()).expect("parses");
            let files: Vec<&FileEntry> = manifest.files.iter().filter(|f| f.needed_for(&platforms)).collect();
            assert!(!words.contains(&"--allow-unpinned"), "{line}");
            assert_eq!(words.contains(&"--accept-licences"), files.iter().any(|f| f.accept), "{line}");
            seen += 1;
        }
        assert!(seen > 0, "the README shows a setup command");
    }

    #[test]
    fn a_file_says_why_it_has_no_checksum() {
        assert!(Manifest::parse(&file("a", &format!("sha256 = \"{SHA}\""))).is_ok());
        assert!(Manifest::parse(&file("a", "md5 = \"00000000000000000000000000000000\"")).is_ok());
        assert!(Manifest::parse(&file("a", "unpinned = \"none published\"")).is_ok());
        let none = Manifest::parse(&file("a", "")).expect_err("neither");
        assert!(format!("{none:#}").contains("no sha256 or md5"), "{none:#}");
        assert!(Manifest::parse(&file("a", &format!("sha256 = \"{SHA}\"\nunpinned = \"x\""))).is_err());
        assert!(Manifest::parse(&file("a", "sha256 = \"ABC\"")).is_err(), "a malformed checksum");
    }

    #[test]
    fn paths_stay_in_the_setup_directory() {
        for path in ["../x", "/etc/x", "a/../../x", ""] {
            let text = format!(
                "[[file]]\nid = \"a\"\nurl = \"file:///a\"\npath = \"{path}\"\nunpinned = \"x\"\nlicence = \"x\"\n"
            );
            assert!(Manifest::parse(&text).is_err(), "{path}");
        }
    }

    #[test]
    fn nothing_setup_writes_collides() {
        let reference = "[site]\nreference = \"a\"\nmetadata = \"m.mdr\"\n";
        let other = |path: &str| {
            format!("[[file]]\nid = \"b\"\nurl = \"file:///b\"\npath = \"{path}\"\nunpinned = \"x\"\nlicence = \"x\"\n")
        };
        let base = format!("{reference}{}", file("a", "unpinned = \"x\""));
        assert!(Manifest::parse(&base).is_ok());
        for path in [
            "a.bed.fai",
            "a.bed.gzi",
            "m.mdr",
            "site.toml",
            "a.bed.part",
            ".lrpipe-setup/x",
            "a.bed/x",
            "site.toml.previous",
        ] {
            assert!(Manifest::parse(&format!("{base}{}", other(path))).is_err(), "{path}");
        }
        // An unpacked file's download, in sources/.
        let unpacked = format!(
            "{reference}[[file]]\nid = \"a\"\nurl = \"file:///a\"\npath = \"a.fa\"\nunpack = \"gunzip\"\nunpinned = \"x\"\nlicence = \"x\"\n"
        );
        assert!(Manifest::parse(&unpacked).is_ok());
        assert!(Manifest::parse(&format!("{unpacked}{}", other("sources/a.fa.gz"))).is_err());
    }

    #[test]
    fn names_must_name_something() {
        let base = file("a", "unpinned = \"x\"");
        assert!(
            Manifest::parse(&format!(
                "{base}[[build]]\nid = \"b\"\nkind = \"catalog\"\ninputs = [\"zz\"]\npath = \"b.bed\"\n"
            ))
            .is_err()
        );
        assert!(Manifest::parse(&format!("[site]\ncatalog = \"zz\"\n{base}")).is_err());
        assert!(Manifest::parse(&format!("[models]\nont = \"zz\"\n{base}")).is_err());
        assert!(Manifest::parse(&format!("{base}{}", file("a", "unpinned = \"x\""))).is_err(), "an id twice");
        assert!(Manifest::parse(&format!("format = 2\n{base}")).is_err(), "another format");
    }

    #[test]
    fn pinning_writes_the_checksum_and_keeps_the_rest() {
        let text =
            format!("# kept\n{}\n# also kept\n{}", file("a", "unpinned = \"none\""), file("b", "unpinned = \"none\""));
        let pinned = pin_text(&text, &[("a".into(), SHA.into())]).expect("pinned");
        let manifest = Manifest::parse(&pinned).expect("parses");
        assert_eq!(manifest.file("a").and_then(|f| f.sha256.as_deref()), Some(SHA));
        assert!(manifest.file("a").is_some_and(|f| f.unpinned.is_none()));
        assert!(manifest.file("b").is_some_and(|f| f.unpinned.is_some()), "only the file pinned changes");
        assert!(pinned.contains("# kept") && pinned.contains("# also kept"));
        // Pinning again replaces the checksum.
        let other = "1".repeat(64);
        let again = pin_text(&pinned, &[("a".into(), other.clone())]).expect("pinned again");
        assert_eq!(Manifest::parse(&again).expect("parses").file("a").and_then(|f| f.sha256.clone()), Some(other));
        assert!(pin_text(&text, &[("zz".into(), SHA.into())]).is_err());
    }
}

//! `lrpipe setup`: fetches a genome's reference, Clair3's models and the catalogues from a
//! versioned manifest ([`manifest`]), builds what is missing (the reference's `.fai`, a `.gzi`
//! for a bgzipped one, lrpipe-cnv's `.mdr`, a repeat catalogue made of several, and with
//! `--bonefish` bonefish's index of the reference, for aligning unaligned reads), and writes the
//! site configuration, `DIR/site.toml`, with absolute paths.
//!
//! Every file is checked against the manifest's checksums before it is used. What was checked
//! is recorded in `DIR/.lrpipe-setup/state.json` (each file's size, time and SHA-256, and for a
//! fetched one the size and checksums of what was fetched), so a rerun skips what is present,
//! unchanged and still as the manifest says, without reading it again; `--verify` reads every
//! file again. A file the manifest has no checksum of is fetched and used like any other, checked
//! by its size where the manifest gives one, and is reported as checked by size only, never as
//! verified. One setup at a time holds the directory (`DIR/.lrpipe-setup/lock`).

pub mod fetch;
pub mod manifest;

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::UNIX_EPOCH;

use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use sha2::{Digest as _, Sha256};

use self::fetch::{Digests, Expected, Fetcher};
use self::manifest::{BuildKind, FileEntry, Manifest, Unpack};
use crate::config::Site;
use crate::plan::Platform;
use crate::provenance::hex;

/// Where a setup directory's record of what was verified is.
const STATE: &str = ".lrpipe-setup/state.json";

/// The lock one setup of a directory at a time holds.
const LOCK: &str = ".lrpipe-setup/lock";

/// How a setup runs.
#[derive(Clone, Debug)]
pub struct Options {
    /// The setup directory.
    pub dir: PathBuf,
    /// The platforms whose models are fetched.
    pub platforms: Vec<Platform>,
    /// Whether the terms of files that need them accepted are accepted.
    pub accept_licences: bool,
    /// The bonefish to build the reference's index with, which `site.toml` then names with the
    /// index; without it, a setup keeps the bonefish and index of the one before.
    pub bonefish: Option<PathBuf>,
}

/// What a setup directory's files were checked against, by path under the directory.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
struct State {
    /// The manifest the directory was set up from.
    manifest: String,
    /// The platforms it was set up for.
    platforms: Vec<String>,
    /// Each file's record.
    files: BTreeMap<String, Record>,
    /// The SHA-256 of the `site.toml` the last setup made, which tells one edited by hand since.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    site_toml: Option<String>,
    /// The bonefish the reference's index was built with.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    bonefish: Option<PathBuf>,
}

/// A file as it was when it was checked or made.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
struct Record {
    /// Bytes.
    size: u64,
    /// When it was last changed: seconds and nanoseconds since 1970.
    modified: (u64, u32),
    /// Its SHA-256.
    sha256: String,
    /// What it was made from: for a file built here, the SHA-256s of its inputs.
    basis: String,
    /// For a fetched file, what was fetched (for one gunzipped here, the download), which the
    /// manifest's checksums are compared with: a manifest that gains a checksum (`--pin`) keeps
    /// the files it matches.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    fetched: Option<Digests>,
}

impl Record {
    /// Records `path`, whose SHA-256 is `sha256`.
    fn of(path: &Path, sha256: String, basis: String, fetched: Option<Digests>) -> Result<Self> {
        let (size, modified) = stamp(path)?;
        Ok(Self { size, modified, sha256, basis, fetched })
    }

    /// Whether `path` is still the file recorded, by size and time.
    fn unchanged(&self, path: &Path) -> bool {
        stamp(path).is_ok_and(|s| s == (self.size, self.modified))
    }
}

/// What a fetched `file` must be.
fn expected(file: &FileEntry) -> Expected<'_> {
    Expected { size: file.size, sha256: file.sha256.as_deref(), md5: file.md5.as_deref() }
}

/// A file's size and modification time.
fn stamp(path: &Path) -> Result<(u64, (u64, u32))> {
    let meta = fs::metadata(path).with_context(|| format!("reading {}", path.display()))?;
    let since = meta.modified()?.duration_since(UNIX_EPOCH).unwrap_or_default();
    Ok((meta.len(), (since.as_secs(), since.subsec_nanos())))
}

impl State {
    fn load(dir: &Path) -> Result<Self> {
        let path = dir.join(STATE);
        match fs::read_to_string(&path) {
            Ok(text) => serde_json::from_str(&text).with_context(|| format!("parsing {}", path.display())),
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(Self::default()),
            Err(e) => Err(e).with_context(|| format!("reading {}", path.display())),
        }
    }

    fn save(&self, dir: &Path) -> Result<()> {
        let path = dir.join(STATE);
        let parent = path.parent().expect("in .lrpipe-setup");
        fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
        let tmp = path.with_extension("json.tmp");
        fs::write(&tmp, serde_json::to_string_pretty(self)?).with_context(|| format!("writing {}", tmp.display()))?;
        fs::rename(&tmp, &path).with_context(|| format!("writing {}", path.display()))
    }

    /// The record of the built `path` if it still holds the file, made from `basis`.
    fn current(&self, dir: &Path, path: &str, basis: &str) -> Option<&Record> {
        self.files.get(path).filter(|r| r.basis == basis && r.unchanged(&dir.join(path)))
    }

    /// Whether the fetched `file` is there unchanged since it was fetched, and what was fetched
    /// is what the manifest says (by size alone for a file it has no checksum of).
    fn fetched_current(&self, dir: &Path, file: &FileEntry) -> bool {
        self.files.get(&file.path).is_some_and(|r| {
            r.unchanged(&dir.join(&file.path))
                && r.fetched.as_ref().is_some_and(|d| d.mismatch(&expected(file)).is_none())
        })
    }

    /// The platforms recorded.
    fn platforms(&self) -> Result<Vec<Platform>> {
        self.platforms
            .iter()
            .map(|p| match p.as_str() {
                "ont" => Ok(Platform::Ont),
                "hifi" => Ok(Platform::Hifi),
                other => bail!("{STATE}: platform {other}"),
            })
            .collect()
    }
}

/// Holds `dir` against other setups until it is dropped. A lock still held after
/// [`LOCK_WAIT`] is another setup's: a process that another thread of this one forks shares the
/// lock's file until it execs, so a setup started just after an earlier one of the same process
/// can find the lock held for a moment (the release gate's tests did).
fn lock(dir: &Path) -> Result<File> {
    let path = dir.join(LOCK);
    let parent = path.parent().expect("in .lrpipe-setup");
    fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
    let file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
    let deadline = std::time::Instant::now() + LOCK_WAIT;
    loop {
        match file.try_lock() {
            Ok(()) => return Ok(file),
            Err(fs::TryLockError::WouldBlock) if std::time::Instant::now() < deadline => {
                std::thread::sleep(std::time::Duration::from_millis(20));
            }
            Err(fs::TryLockError::WouldBlock) => bail!("another lrpipe setup is using {}", dir.display()),
            Err(fs::TryLockError::Error(e)) => return Err(e).with_context(|| format!("locking {}", path.display())),
        }
    }
}

/// How long [`lock`] waits for a held lock before taking it for another setup's.
const LOCK_WAIT: std::time::Duration = std::time::Duration::from_secs(2);

/// Writes a progress line to stderr, dropping it if stderr is gone.
fn say(line: &str) {
    let _ = writeln!(io::stderr(), "lrpipe setup: {line}");
}

/// The files of `manifest` a setup for `platforms` needs.
fn needed<'a>(manifest: &'a Manifest, platforms: &[Platform]) -> Vec<&'a FileEntry> {
    manifest.files.iter().filter(|f| f.needed_for(platforms)).collect()
}

/// The terms of the files in `files` that must be accepted, as printed: each licence's once,
/// after the files (or models) under it.
fn terms(files: &[&FileEntry]) -> String {
    // (licence, terms, licence URL, the files or models under them)
    type Terms<'a> = (&'a str, &'a str, Option<&'a str>, Vec<&'a str>);
    let mut groups: Vec<Terms<'_>> = Vec::new();
    for f in files.iter().filter(|f| f.accept) {
        let what = f.model.as_deref().unwrap_or(&f.id);
        let key = (f.licence.as_str(), f.terms.as_deref().unwrap_or_default(), f.licence_url.as_deref());
        match groups.iter_mut().find(|g| (g.0, g.1, g.2) == key) {
            Some(group) if group.3.contains(&what) => {}
            Some(group) => group.3.push(what),
            None => groups.push((key.0, key.1, key.2, vec![what])),
        }
    }
    let mut out = String::new();
    for (licence, terms, url, names) in groups {
        let _ = writeln!(out, "  {}: {licence}", names.join(", "));
        let _ = writeln!(out, "    {terms}");
        if let Some(url) = url {
            let _ = writeln!(out, "    {url}");
        }
    }
    out
}

/// Sets up `options.dir` from `manifest`: fetches and checks what is missing, builds what is
/// missing, and writes `site.toml`. Returns the path of the site configuration it wrote, or that
/// held it already: `site.toml`, or `site.toml.new` beside one changed by hand.
///
/// # Errors
///
/// If a file needs its terms accepted and they are not, cannot be fetched or does not match the
/// manifest, or a build fails.
pub fn setup(manifest_text: &str, options: &Options, fetcher: &Fetcher) -> Result<PathBuf> {
    let manifest = Manifest::parse(manifest_text)?;
    let genome = manifest.genome.clone().context("the manifest names no genome")?;
    fs::create_dir_all(&options.dir).with_context(|| format!("creating {}", options.dir.display()))?;
    let dir = fs::canonicalize(&options.dir).with_context(|| format!("resolving {}", options.dir.display()))?;
    let _lock = lock(&dir)?;
    let mut state = State::load(&dir)?;
    // The SHA-256 of the site.toml the last setup made: recorded, or for a setup from before
    // lrpipe recorded it, of the one its manifest and platforms make.
    let made_before = state.site_toml.clone().or_else(|| {
        let earlier = Manifest::parse(&state.manifest).ok()?;
        let reference = earlier.site.reference.as_deref().and_then(|r| earlier.file(r))?;
        let index = bonefish_index_path(&reference.path);
        let align = state.bonefish.as_deref().map(|b| (b, index.as_str()));
        let text = site_toml(&earlier, &dir, &state.platforms().ok()?, align);
        Some(hex(&Sha256::digest(text.as_bytes())))
    });
    // A rerun for other platforms keeps the earlier ones whose models are still there.
    let mut platforms = options.platforms.clone();
    let same_genome = Manifest::parse(&state.manifest).is_ok_and(|m| m.genome == manifest.genome);
    for earlier in if same_genome { state.platforms()? } else { Vec::new() } {
        if platforms.contains(&earlier) {
            continue;
        }
        if manifest.files.iter().filter(|f| f.platform == Some(earlier)).all(|f| state.fetched_current(&dir, f)) {
            say(&format!("keeping the {} models set up before", earlier.as_str()));
            platforms.push(earlier);
        } else {
            say(&format!("the {} models set up before are not all there: leaving them out", earlier.as_str()));
        }
    }
    let files = needed(&manifest, &platforms);
    let to_fetch: Vec<&FileEntry> = files.iter().copied().filter(|f| !state.fetched_current(&dir, f)).collect();
    let to_accept = terms(&to_fetch);
    if !to_accept.is_empty() {
        if !options.accept_licences {
            bail!(
                "these files' terms must be accepted before they are fetched:\n{to_accept}Give --accept-licences to accept them"
            );
        }
        say(&format!("accepting the terms of:\n{to_accept}").trim_end().replace('\n', "\n  "));
    }

    state.manifest = manifest_text.to_owned();
    state.platforms = platforms.iter().map(|p| p.as_str().to_owned()).collect();
    for file in &files {
        get_file(&dir, file, fetcher, &mut state)?;
        state.save(&dir)?;
    }
    for build in &manifest.builds {
        if build.inputs.iter().all(|i| files.iter().any(|f| f.id == *i)) {
            make_catalog(&dir, &manifest, build, &mut state)?;
            state.save(&dir)?;
        }
    }
    let reference = manifest.site.reference.as_deref().and_then(|r| manifest.file(r)).expect("checked by parse");
    index_reference(&dir, &reference.path, &mut state)?;
    state.save(&dir)?;
    if let Some(mdr) = &manifest.site.metadata {
        make_mdr(&dir, &reference.path, mdr, &mut state)?;
        state.save(&dir)?;
    }
    // The bonefish asked for, else the one of the setup before, while it is still there.
    let earlier = state.bonefish.clone().filter(|_| same_genome);
    let bonefish = match (&options.bonefish, earlier) {
        (Some(bonefish), _) => Some(bonefish.clone()),
        (None, Some(earlier)) if earlier.is_file() => {
            say(&format!("keeping the bonefish index of the setup before ({})", earlier.display()));
            Some(earlier)
        }
        (None, Some(earlier)) => {
            say(&format!("{} is not there any more: leaving out the bonefish index", earlier.display()));
            None
        }
        (None, None) => None,
    };
    state.bonefish = None;
    let align = match &bonefish {
        Some(bonefish) => {
            let index = bonefish_index(&dir, &reference.path, bonefish, &mut state)?;
            state.bonefish = Some(bonefish.clone());
            Some((bonefish.as_path(), index))
        }
        None => None,
    };
    state.save(&dir)?;
    for missing in &manifest.missing {
        say(&format!("not available for {genome}: {missing}"));
    }
    let site = dir.join("site.toml");
    let text = site_toml(&manifest, &dir, &platforms, align.as_ref().map(|(b, i)| (*b, i.as_str())));
    let written = write_if_changed(&site, &text, made_before.as_deref())?;
    state.site_toml = Some(hex(&Sha256::digest(text.as_bytes())));
    state.save(&dir)?;
    match written {
        Written::Unchanged => {
            say(&format!("{} holds this setup's configuration already", site.display()));
            Ok(site)
        }
        Written::Site => {
            say(&format!("wrote {}", site.display()));
            Ok(site)
        }
        Written::Beside(new) => {
            say(&format!(
                "kept {}, which was changed since lrpipe setup wrote it, and wrote this setup's as {}: merge what you \
                 changed into it and move it to site.toml, or keep yours",
                site.display(),
                new.display()
            ));
            Ok(new)
        }
    }
}

/// Where [`write_if_changed`] put the site configuration.
#[derive(Debug, PartialEq, Eq)]
enum Written {
    /// Nowhere: `site.toml` held it already.
    Unchanged,
    /// In `site.toml`.
    Site,
    /// Beside a `site.toml` changed by hand, which is kept: in this file (`site.toml.new`).
    Beside(PathBuf),
}

/// Makes `file` present and checked at its path under `dir`.
fn get_file(dir: &Path, file: &FileEntry, fetcher: &Fetcher, state: &mut State) -> Result<()> {
    if state.fetched_current(dir, file) {
        say(&format!("{}: present and {}", file.path, verdict(file)));
        return Ok(());
    }
    let dest = dir.join(&file.path);
    let expected = expected(file);
    let download = dir.join(file.download());
    // A file put there by hand (fetched with --offline's commands) or not recorded.
    let mut digests = None;
    if download.exists() {
        let found = Digests::of(&download).with_context(|| format!("reading {}", download.display()))?;
        match found.mismatch(&expected) {
            None => digests = Some(found),
            Some(why) => {
                say(&format!("{}: {why}; fetching it again", file.download()));
                fs::remove_file(&download).with_context(|| format!("removing {}", download.display()))?;
            }
        }
    }
    let digests = match digests {
        Some(d) => d,
        None => {
            let size = file.size.map_or_else(String::new, |s| format!(" ({})", fetch::bytes(s)));
            say(&format!("fetching {}{size} from {}", file.download(), file.url));
            fetcher.fetch(&file.url, &download, &expected, &say)?
        }
    };
    let verdict = verdict(file);
    match file.unpack {
        None => {
            let record = Record::of(&dest, digests.sha256.clone(), String::new(), Some(digests))?;
            state.files.insert(file.path.clone(), record);
            say(&format!("{}: {verdict}", file.path));
        }
        Some(Unpack::Gunzip) => {
            say(&format!("unpacking {} to {}", file.download(), file.path));
            let sha256 = gunzip(&download, &dest)?;
            state.files.insert(file.path.clone(), Record::of(&dest, sha256, String::new(), Some(digests))?);
            // Recorded before the download goes: without it, the file could not be told good.
            state.save(dir)?;
            say(&format!("{}: {verdict}", file.path));
            if let Err(e) = fs::remove_file(&download) {
                say(&format!("could not remove {} ({e}); it is not needed", download.display()));
            }
        }
    }
    Ok(())
}

/// How a fetched `file` that matches the manifest was checked: against its checksums, or, for one
/// the manifest has none of, by its size alone (or, with no size either, not at all).
fn verdict(file: &FileEntry) -> String {
    if file.pinned() {
        format!("verified ({})", file.expected())
    } else if let Some(size) = file.size {
        format!("checked by size only ({size} bytes; the manifest has no checksum of it)")
    } else {
        "not checked (the manifest has no checksum or size of it)".into()
    }
}

/// A writer that computes the SHA-256 of what it writes.
struct Hashing<W> {
    inner: W,
    hasher: Sha256,
}

impl<W: Write> Write for Hashing<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.hasher.update(&buf[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Decompresses `src` (gzip or BGZF) to `dest`, by way of `dest.part`. Returns the SHA-256 of
/// what it wrote.
fn gunzip(src: &Path, dest: &Path) -> Result<String> {
    let input = File::open(src).with_context(|| format!("opening {}", src.display()))?;
    let mut reader = flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, input));
    write_atomically(dest, None, |out| {
        io::copy(&mut reader, out).with_context(|| format!("decompressing {}", src.display()))?;
        Ok(())
    })
}

/// Writes `dest` through `write`, to `dest.part` renamed into place when it is done and, if
/// `sha256` is given, has it (otherwise `dest.part` is removed and `dest` left as it was).
/// Returns the SHA-256 of what was written.
fn write_atomically(
    dest: &Path,
    sha256: Option<&str>,
    write: impl FnOnce(&mut dyn Write) -> Result<()>,
) -> Result<String> {
    let part = fetch::part_path(dest);
    if let Some(parent) = dest.parent() {
        fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
    }
    let file = File::create(&part).with_context(|| format!("creating {}", part.display()))?;
    let mut out = Hashing { inner: BufWriter::with_capacity(1 << 20, file), hasher: Sha256::new() };
    write(&mut out)?;
    out.flush().with_context(|| format!("writing {}", part.display()))?;
    let Hashing { inner, hasher } = out;
    inner.into_inner().map_err(io::IntoInnerError::into_error)?.sync_all()?;
    let written = hex(&hasher.finalize());
    if let Some(want) = sha256.filter(|w| *w != written) {
        fs::remove_file(&part).with_context(|| format!("removing {}", part.display()))?;
        bail!("{} was made with sha256 {written}, not the manifest's {want}", dest.display());
    }
    fs::rename(&part, dest).with_context(|| format!("renaming {} to {}", part.display(), dest.display()))?;
    Ok(written)
}

/// Builds a repeat catalogue from its inputs, unless it is there and made from them as they are.
fn make_catalog(dir: &Path, manifest: &Manifest, build: &manifest::Build, state: &mut State) -> Result<()> {
    let BuildKind::Catalog = build.kind;
    let inputs: Vec<&str> = build.inputs.iter().filter_map(|i| manifest.file(i)).map(|f| f.path.as_str()).collect();
    let sums: Vec<&str> = inputs.iter().filter_map(|p| state.files.get(*p)).map(|r| r.sha256.as_str()).collect();
    let basis = format!("catalogue of {}", sums.join(" "));
    if state.current(dir, &build.path, &basis).is_some_and(|r| build.sha256.as_ref().is_none_or(|w| *w == r.sha256)) {
        say(&format!("{}: present", build.path));
        return Ok(());
    }
    say(&format!("making {} from {}", build.path, inputs.join(" and ")));
    let mut lines: Vec<Vec<u8>> = Vec::new();
    for input in &inputs {
        let path = dir.join(input);
        read_lines(&path, &mut lines).with_context(|| format!("reading {}", path.display()))?;
    }
    lines.sort_by(|a, b| sort_key(a).cmp(&sort_key(b)).then_with(|| a.cmp(b)));
    let dest = dir.join(&build.path);
    let sha256 = write_atomically(&dest, build.sha256.as_deref(), |out| {
        for line in &lines {
            out.write_all(line)?;
            out.write_all(b"\n")?;
        }
        Ok(())
    })?;
    state.files.insert(build.path.clone(), Record::of(&dest, sha256, basis, None)?);
    Ok(())
}

/// The non-empty lines of `path`, gzip, BGZF or plain.
fn read_lines(path: &Path, lines: &mut Vec<Vec<u8>>) -> Result<()> {
    let mut file = BufReader::with_capacity(1 << 20, File::open(path)?);
    let gzip = file.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    let reader: Box<dyn BufRead> =
        if gzip { Box::new(BufReader::new(flate2::bufread::MultiGzDecoder::new(file))) } else { Box::new(file) };
    for line in reader.split(b'\n') {
        let mut line = line?;
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if !line.is_empty() {
            lines.push(line);
        }
    }
    Ok(())
}

/// A BED line's place as `LC_ALL=C sort -k1,1 -k2,2n` orders lines: its first field's bytes,
/// then its second field's number (leading blanks skipped, the digits before anything else, 0
/// without any). Lines equal in both are then ordered by their bytes.
fn sort_key(line: &[u8]) -> (&[u8], u64) {
    let blank = |b: &u8| *b == b' ' || *b == b'\t';
    let end = line.iter().position(blank).unwrap_or(line.len());
    let rest = &line[end..];
    let rest = &rest[rest.iter().position(|b| !blank(b)).unwrap_or(rest.len())..];
    let digits = rest.iter().take_while(|b| b.is_ascii_digit()).count();
    let start = std::str::from_utf8(&rest[..digits]).ok().and_then(|s| s.parse().ok()).unwrap_or(0);
    (&line[..end], start)
}

/// Whether `path` starts as a BGZF file does (gzip with a `BC` extra subfield).
fn is_bgzf(path: &Path) -> Result<bool> {
    let mut head = Vec::with_capacity(18);
    File::open(path)?.take(18).read_to_end(&mut head)?;
    Ok(head.len() == 18 && head[..4] == [0x1f, 0x8b, 8, 4] && head[12..14] == *b"BC")
}

/// Builds the reference's `.fai`, and a bgzipped one's `.gzi`, unless they are there and made
/// from the reference as it is.
fn index_reference(dir: &Path, reference: &str, state: &mut State) -> Result<()> {
    let path = dir.join(reference);
    let sha = state.files.get(reference).map(|r| r.sha256.clone()).unwrap_or_default();
    let basis = format!("index of {sha}");
    let bgzf = is_bgzf(&path).with_context(|| format!("reading {}", path.display()))?;
    let fai = format!("{reference}.fai");
    if state.current(dir, &fai, &basis).is_none() {
        say(&format!("indexing {reference}"));
        let index = if bgzf {
            let reader = pasta::bgzf::io::Reader::new(BufReader::new(File::open(&path)?));
            let mut indexer = pasta::fasta::io::Indexer::new(reader);
            let mut records = Vec::new();
            while let Some(record) = indexer.index_record().with_context(|| format!("indexing {}", path.display()))? {
                records.push(record);
            }
            pasta::fasta::fai::Index::from(records)
        } else {
            pasta::fasta::fs::index(&path).with_context(|| format!("indexing {}", path.display()))?
        };
        let dest = dir.join(&fai);
        let sha256 = write_atomically(&dest, None, |out| {
            pasta::fasta::fai::io::Writer::new(out).write_index(&index)?;
            Ok(())
        })?;
        state.files.insert(fai.clone(), Record::of(&dest, sha256, basis.clone(), None)?);
    }
    let gzi = format!("{reference}.gzi");
    if bgzf && state.current(dir, &gzi, &basis).is_none() {
        say(&format!("indexing {reference}'s BGZF blocks"));
        let index = gzi_of(&path).with_context(|| format!("indexing {}", path.display()))?;
        let dest = dir.join(&gzi);
        let sha256 = write_atomically(&dest, None, |out| {
            pasta::bgzf::gzi::io::Writer::new(out).write_index(&index)?;
            Ok(())
        })?;
        state.files.insert(gzi, Record::of(&dest, sha256, basis, None)?);
    }
    Ok(())
}

/// The `.gzi` of a BGZF file, as `bgzip -r` makes it: each block with data after the first, by
/// its compressed offset and the uncompressed bytes before it.
fn gzi_of(path: &Path) -> Result<pasta::bgzf::gzi::Index> {
    let mut reader = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut frame = Vec::new();
    let (mut compressed, mut uncompressed, mut started) = (0u64, 0u64, false);
    let mut index = Vec::new();
    while pasta::bgzf::io::raw::read_frame(&mut reader, &mut frame)? {
        // ISIZE, the frame's last four bytes: its uncompressed size.
        let tail: [u8; 4] =
            frame.get(frame.len().saturating_sub(4)..).and_then(|t| t.try_into().ok()).context("a short BGZF block")?;
        let size = u64::from(u32::from_le_bytes(tail));
        if size > 0 {
            if started {
                index.push((compressed, uncompressed));
            }
            started = true;
            uncompressed += size;
        }
        compressed += frame.len() as u64;
    }
    Ok(pasta::bgzf::gzi::Index::from(index))
}

/// Makes lrpipe-cnv's `.mdr` of the reference's N runs, unless it is there and made from the
/// reference as it is.
fn make_mdr(dir: &Path, reference: &str, mdr: &str, state: &mut State) -> Result<()> {
    let sha = state.files.get(reference).map(|r| r.sha256.clone()).unwrap_or_default();
    let basis = format!("N runs of {sha}");
    if state.current(dir, mdr, &basis).is_some() {
        say(&format!("{mdr}: present"));
        return Ok(());
    }
    say(&format!("making {mdr} from {reference} (its N runs, for lrpipe-cnv)"));
    let dest = dir.join(mdr);
    let part = fetch::part_path(&dest);
    let defaults = lrpipe_cnv::pipeline::Options::new(PathBuf::new(), String::new(), PathBuf::new(), PathBuf::new());
    lrpipe_cnv::metadata::write_mdr(&dir.join(reference), defaults.n_size, &part)
        .with_context(|| format!("making {}", dest.display()))?;
    fs::rename(&part, &dest).with_context(|| format!("renaming {} to {}", part.display(), dest.display()))?;
    let sha256 = fetch::Digests::of(&dest)?.sha256;
    state.files.insert(mdr.to_owned(), Record::of(&dest, sha256, basis, None)?);
    Ok(())
}

/// The path under `dir` of bonefish's index of `reference`.
fn bonefish_index_path(reference: &str) -> String {
    format!("{reference}.bonefish.raw")
}

/// Builds bonefish's raw index of the reference (`bonefish -d INDEX --idx-raw REFERENCE`, which
/// lrpipe run's align stage maps to and which loads in a second), unless it is there and was made
/// from the reference as it is by the same version of bonefish. Returns its path under `dir`.
fn bonefish_index(dir: &Path, reference: &str, bonefish: &Path, state: &mut State) -> Result<String> {
    let tool = crate::checks::bonefish(bonefish)?;
    let sha = state.files.get(reference).map(|r| r.sha256.clone()).unwrap_or_default();
    let basis = format!("raw index of {sha} by {}", tool.version);
    let index = bonefish_index_path(reference);
    if state.current(dir, &index, &basis).is_some() {
        say(&format!("{index}: present"));
        return Ok(index);
    }
    say(&format!("building {index} with {} (bonefish's index of {reference})", tool.version));
    let dest = dir.join(&index);
    let part = fetch::part_path(&dest);
    let out = Command::new(&tool.path)
        .arg("-d")
        .arg(&part)
        .arg("--idx-raw")
        .arg(dir.join(reference))
        .output()
        .with_context(|| format!("running {}", tool.path.display()))?;
    if !out.status.success() {
        let _ = fs::remove_file(&part);
        let stderr = String::from_utf8_lossy(&out.stderr);
        bail!("{} -d failed ({}): {}", tool.path.display(), out.status, stderr.trim());
    }
    fs::rename(&part, &dest).with_context(|| format!("renaming {} to {}", part.display(), dest.display()))?;
    let sha256 = fetch::Digests::of(&dest)?.sha256;
    state.files.insert(index.clone(), Record::of(&dest, sha256, basis, None)?);
    Ok(index)
}

/// Writes the site configuration `text` to `path` unless it holds it already: to a temporary
/// file first, renamed over `path` once it loads as a site configuration. An earlier, different
/// one that a setup made (its SHA-256 is `made`) is copied to `site.toml.previous` and replaced;
/// one changed since, or not known to be a setup's (no setup is recorded), is someone's own and
/// is kept, `text` going to `site.toml.new` beside it. Replacing a configuration edited by hand (a tool binary, a GPU,
/// a blacklist added) lost the edits but for the `.previous` copy, which the next rerun replaced
/// too.
fn write_if_changed(path: &Path, text: &str, made: Option<&str>) -> Result<Written> {
    let old = fs::read_to_string(path).ok();
    if old.as_deref() == Some(text) {
        return Ok(Written::Unchanged);
    }
    let tmp = path.with_extension("toml.tmp");
    fs::write(&tmp, text).with_context(|| format!("writing {}", tmp.display()))?;
    if let Err(e) = Site::load(&tmp) {
        let _ = fs::remove_file(&tmp);
        return Err(e).context("the site configuration lrpipe setup made does not load");
    }
    if let Some(old) = &old {
        if made != Some(hex(&Sha256::digest(old.as_bytes())).as_str()) {
            let new = path.with_extension("toml.new");
            fs::rename(&tmp, &new).with_context(|| format!("writing {}", new.display()))?;
            return Ok(Written::Beside(new));
        }
        let previous = path.with_extension("toml.previous");
        fs::copy(path, &previous).with_context(|| format!("keeping {}", previous.display()))?;
        say(&format!("kept the earlier {} as {}", path.display(), previous.display()));
    }
    fs::rename(&tmp, path).with_context(|| format!("writing {}", path.display()))?;
    Ok(Written::Site)
}

/// A TOML string holding `path`.
fn quoted(path: &Path) -> String {
    toml::Value::String(path.to_string_lossy().into_owned()).to_string()
}

/// The site configuration for a setup of `dir`, with `align`'s bonefish and its index (under
/// `dir`) when it has one.
fn site_toml(manifest: &Manifest, dir: &Path, platforms: &[Platform], align: Option<(&Path, &str)>) -> String {
    let genome = manifest.genome.as_deref().unwrap_or_default();
    let names: Vec<&str> = platforms.iter().map(|p| p.as_str()).collect();
    let mut out = String::new();
    let _ = writeln!(
        out,
        "# lrpipe's site configuration for {genome} ({} reads), written by `lrpipe setup`",
        names.join(" and ")
    );
    let _ = writeln!(
        out,
        "# (lrpipe {}). Every path is absolute. To check the files against the manifest:",
        env!("CARGO_PKG_VERSION")
    );
    let _ = writeln!(out, "#   lrpipe setup --verify --dir {}", shell(&dir.to_string_lossy()));
    let _ = writeln!(out, "#\n# Licences:");
    let mut seen: Vec<&str> = Vec::new();
    for f in needed(manifest, platforms) {
        let what = f.model.as_deref().unwrap_or(&f.path);
        if seen.contains(&what) {
            continue;
        }
        seen.push(what);
        let accepted = if f.accept { " (its terms accepted with --accept-licences)" } else { "" };
        let _ = writeln!(out, "#   {what}: {}{accepted}", f.licence);
        if let Some(url) = &f.licence_url {
            let _ = writeln!(out, "#     {url}");
        }
    }
    let size_only: Vec<&str> =
        needed(manifest, platforms).into_iter().filter(|f| !f.pinned()).map(|f| f.path.as_str()).collect();
    if !size_only.is_empty() {
        let _ = writeln!(
            out,
            "#\n# Checked by size only (the manifest has no checksum of them; `lrpipe setup --pin` records one):"
        );
        for p in size_only {
            let _ = writeln!(out, "#   {p}");
        }
    }
    if !manifest.missing.is_empty() {
        let _ = writeln!(out, "#\n# Not available for {genome}:");
        for m in &manifest.missing {
            let _ = writeln!(out, "#   {m}");
        }
    }
    let path = |id: &str| manifest.path_of(id).map(|p| quoted(&dir.join(p)));
    let reference = manifest.site.reference.as_deref().and_then(path).unwrap_or_default();
    let _ = writeln!(out, "\nreference = {reference}");
    let _ = writeln!(out, "\n[snv]");
    // The directory, not one model per platform: lrpipe run takes each run's model from the
    // reads' header, which a fixed ont_model or hifi_model would override.
    let models = quoted(&dir.join("models"));
    let _ = writeln!(out, "# Every model fetched, under Clair3's names: each run takes the one for its reads' header.");
    let _ = writeln!(out, "models = {models}");
    let _ = writeln!(out, "# To use one model whatever the header says:");
    for platform in platforms {
        if let Some(model) = manifest.models.of(*platform) {
            let _ = writeln!(out, "# {}_model = {}", platform.as_str(), quoted(&dir.join("models").join(model)));
        }
    }
    if let Some(trf) = manifest.site.tandem_repeats.as_deref().and_then(path) {
        let _ = writeln!(out, "\n[sv]\ntandem_repeats = {trf}");
    }
    if let Some(catalog) = manifest.site.catalog.as_deref().and_then(path) {
        let _ = writeln!(out, "\n[repeats]\ncatalog = {catalog}");
    }
    let blacklist = manifest.site.blacklist.as_deref().and_then(path);
    let metadata = manifest.site.metadata.as_ref().map(|m| quoted(&dir.join(m)));
    if blacklist.is_some() || metadata.is_some() {
        let _ = writeln!(out, "\n[cnv]");
        if let Some(b) = blacklist {
            let _ = writeln!(out, "blacklist = {b}");
        }
        if let Some(m) = metadata {
            let _ = writeln!(out, "metadata = {m}");
        }
    }
    if let Some((bonefish, index)) = align {
        let _ = writeln!(out, "\n[tools]\n# Aligns unaligned --reads.\nbonefish = {}", quoted(bonefish));
        let _ = writeln!(out, "\n[align]\n# bonefish's index of the reference.\nindex = {}", quoted(&dir.join(index)));
    }
    out
}

/// One file's state, as `--verify` found it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Check {
    /// Present, and as the manifest (or the record of what it was made from) says.
    Verified(String),
    /// Present and of the size the manifest gives, which has no checksum to check it against
    /// (or, with no size either, present only).
    SizeOnly(String),
    /// Present, but with no record of how it was made, so not known good.
    Unrecorded(String),
    /// Not there.
    Missing,
    /// Present, and not as it should be.
    Mismatch(String),
}

/// What `--verify` found: each file's path under the setup directory and its [`Check`].
#[derive(Clone, Debug, Default)]
pub struct Verification {
    /// The files, in the manifest's order.
    pub files: Vec<(String, Check)>,
}

impl Verification {
    /// Whether every file is present and as the manifest says: verified, or, for one the
    /// manifest has no checksum of, of its size.
    #[must_use]
    pub fn passed(&self) -> bool {
        self.files.iter().all(|(_, c)| matches!(c, Check::Verified(_) | Check::SizeOnly(_)))
    }
}

impl std::fmt::Display for Verification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (path, check) in &self.files {
            let (word, why) = match check {
                Check::Verified(why) => ("ok", why.as_str()),
                Check::SizeOnly(why) => ("size only", why.as_str()),
                Check::Unrecorded(why) => ("UNRECORDED", why.as_str()),
                Check::Missing => ("MISSING", ""),
                Check::Mismatch(why) => ("MISMATCH", why.as_str()),
            };
            writeln!(f, "{word:<11}{path}{}{why}", if why.is_empty() { "" } else { "  " })?;
        }
        let size_only = self.files.iter().filter(|(_, c)| matches!(c, Check::SizeOnly(_))).count();
        if size_only > 0 {
            writeln!(
                f,
                "{size_only} checked by size only: the manifest has no checksum of them (`lrpipe setup --pin` \
                 records one)"
            )?;
        }
        Ok(())
    }
}

/// Checks a setup directory: every file a setup for `platforms` fetches, against the
/// manifest's checksums (a file gunzipped there, against the record of the download it was
/// made from); what was built, against the manifest or the record of what it was made from;
/// and `site.toml`, by loading it. Every file is read again.
///
/// # Errors
///
/// If the manifest does not parse, or the setup directory's record cannot be read.
pub fn verify(manifest_text: &str, dir: &Path, platforms: &[Platform]) -> Result<Verification> {
    let manifest = Manifest::parse(manifest_text)?;
    let state = State::load(dir)?;
    let mut out = Verification::default();
    let files = needed(&manifest, platforms);
    let recorded = |path: &str, sha256: &str| -> Check {
        match state.files.get(path) {
            Some(r) if r.sha256 == sha256 => Check::Verified(r.basis.clone()),
            Some(r) => Check::Mismatch(format!("sha256 {sha256}, not {} as recorded ({})", r.sha256, r.basis)),
            None => Check::Unrecorded("no record of how it was made".into()),
        }
    };
    let sha256_of = |path: &Path| Digests::of(path).map(|d| d.sha256);
    for f in &files {
        let path = dir.join(&f.path);
        let check = if !path.exists() {
            Check::Missing
        } else if f.unpack.is_some() {
            let sha256 = sha256_of(&path)?;
            let record = state.files.get(&f.path);
            match (record, record.and_then(|r| r.fetched.as_ref())) {
                (Some(r), _) if r.sha256 != sha256 => {
                    Check::Mismatch(format!("sha256 {sha256}, not {} as unpacked", r.sha256))
                }
                (Some(_), Some(download)) => match download.mismatch(&expected(f)) {
                    Some(why) => Check::Mismatch(format!("unpacked from a download of {why}")),
                    None if f.pinned() => Check::Verified(format!("gunzipped from {}", f.expected())),
                    None => Check::SizeOnly(size_only(f, "its download")),
                },
                _ => Check::Unrecorded("no record of the download it was unpacked from".into()),
            }
        } else {
            let d = Digests::of(&path).with_context(|| format!("reading {}", path.display()))?;
            match d.mismatch(&expected(f)) {
                Some(why) => Check::Mismatch(why),
                None if f.pinned() => Check::Verified(f.expected()),
                None => Check::SizeOnly(size_only(f, "it")),
            }
        };
        out.files.push((f.path.clone(), check));
    }
    for b in manifest.builds.iter().filter(|b| b.inputs.iter().all(|i| files.iter().any(|f| f.id == *i))) {
        let path = dir.join(&b.path);
        let check = if !path.exists() {
            Check::Missing
        } else {
            let sha256 = sha256_of(&path)?;
            match &b.sha256 {
                Some(want) if *want == sha256 => Check::Verified(format!("sha256:{want}")),
                Some(want) => Check::Mismatch(format!("sha256 {sha256}, not {want}")),
                None => recorded(&b.path, &sha256),
            }
        };
        out.files.push((b.path.clone(), check));
    }
    let reference = manifest.site.reference.as_deref().and_then(|r| manifest.file(r)).expect("checked by parse");
    let mut derived = vec![format!("{}.fai", reference.path)];
    if state.files.contains_key(&format!("{}.gzi", reference.path)) {
        derived.push(format!("{}.gzi", reference.path));
    }
    derived.extend(manifest.site.metadata.clone());
    let index = bonefish_index_path(&reference.path);
    if state.files.contains_key(&index) {
        derived.push(index);
    }
    for d in derived {
        let path = dir.join(&d);
        let check = if path.exists() { recorded(&d, &sha256_of(&path)?) } else { Check::Missing };
        out.files.push((d, check));
    }
    let site = dir.join("site.toml");
    // A site.toml copied with the directory names the original's files.
    let ours = fs::canonicalize(dir).map(|d| d.join(&reference.path)).ok();
    let check = if !site.exists() {
        Check::Missing
    } else {
        match Site::load(&site) {
            Ok(loaded) if ours.as_ref().is_some_and(|r| *r != loaded.reference) => Check::Mismatch(format!(
                "names the reference {}, not this setup's {} (a site.toml copied with its directory names the \
                 original's files: run lrpipe setup again here)",
                loaded.reference.display(),
                ours.as_ref().expect("checked").display()
            )),
            Ok(_) => Check::Verified("a valid site configuration".into()),
            Err(e) => Check::Mismatch(format!("{e:#}")),
        }
    };
    out.files.push(("site.toml".into(), check));
    Ok(out)
}

/// What `--verify` says of `file`, which the manifest has no checksum of: what was checked of
/// `what` (the file, or the download it was unpacked from).
fn size_only(file: &FileEntry, what: &str) -> String {
    match file.size {
        Some(size) => format!("{size} bytes as the manifest says; it has no checksum of {what}"),
        None => format!("present; the manifest has no checksum or size of {what}"),
    }
}

/// The manifest and platforms `dir` was set up with, for `--verify` without `--genome`.
///
/// # Errors
///
/// If the directory has no record of a setup.
pub fn recorded_setup(dir: &Path) -> Result<(String, Vec<Platform>)> {
    let state = State::load(dir)?;
    ensure!(
        !state.manifest.is_empty(),
        "{} has no record of a setup ({STATE}): give --genome and --platform",
        dir.display()
    );
    let platforms = state.platforms()?;
    Ok((state.manifest, platforms))
}

/// A shell word holding `s`: itself, or quoted when a shell would split or expand it.
pub(crate) fn shell(s: &str) -> String {
    if s.bytes().all(|b| b.is_ascii_alphanumeric() || b"/._-+:=@,%".contains(&b)) {
        s.to_owned()
    } else {
        format!("'{}'", s.replace('\'', r"'\''"))
    }
}

/// What the setup that [`offline`] prints to run after fetching names besides the genome's files.
#[derive(Clone, Copy, Debug, Default)]
pub struct Rerun<'a> {
    /// The manifest of one's own it sets up from (`--manifest`), in place of the genome's.
    pub manifest: Option<&'a Path>,
    /// The bonefish it builds the reference's index with (`--bonefish`).
    pub bonefish: Option<&'a Path>,
}

/// For a machine without the internet: the URLs, checksums and `curl` and `wget` commands
/// that fetch every file a setup of `dir` for `platforms` needs to where `lrpipe setup` looks
/// for it, and the setup to run after them (`rerun`'s), which then fetches nothing.
///
/// # Errors
///
/// If the manifest does not parse.
pub fn offline(manifest_text: &str, dir: &Path, platforms: &[Platform], rerun: &Rerun<'_>) -> Result<String> {
    let manifest = Manifest::parse(manifest_text)?;
    let genome = manifest.genome.as_deref().unwrap_or_default();
    let files = needed(&manifest, platforms);
    let dir_text = dir.to_string_lossy();
    let mut out = String::new();
    let _ = writeln!(
        out,
        "# lrpipe setup --offline: the files of {genome} for {}.",
        platforms.iter().map(|p| p.as_str()).collect::<Vec<_>>().join(" and ")
    );
    let _ = writeln!(out, "# Fetch them into {dir_text} (curl's -C - and wget's -c resume a cut-off transfer), copy");
    let _ = writeln!(out, "# the directory to the machine that runs lrpipe, and there run the setup at the end: it");
    let _ = writeln!(out, "# checks them, builds the rest and writes site.toml, fetching nothing.");
    let to_accept = terms(&files);
    if !to_accept.is_empty() {
        let _ = writeln!(out, "#\n# These files' terms must be accepted to use them:");
        for line in to_accept.lines() {
            let _ = writeln!(out, "#{line}");
        }
    }
    let _ = writeln!(out);
    for f in &files {
        let dest = format!("{dir_text}/{}", f.download());
        let mut facts = vec![f.size.map_or_else(|| "size unknown".to_owned(), |s| format!("{s} bytes"))];
        facts.extend(f.sha256.iter().map(|s| format!("sha256 {s}")));
        facts.extend(f.md5.iter().map(|s| format!("md5 {s}")));
        if let Some(why) = &f.unpinned {
            facts.push(format!("no checksum: {why}"));
        }
        if f.accept {
            facts.push("its terms must be accepted (--accept-licences), see above".into());
        }
        let _ = writeln!(out, "# {}: {} ({})", f.id, f.url, facts.join(", "));
        let parent = Path::new(&dest).parent().map(|p| p.to_string_lossy().into_owned()).unwrap_or_default();
        let _ = writeln!(out, "mkdir -p {}", shell(&parent));
        let _ = writeln!(out, "curl -fL -C - -o {} {}", shell(&dest), shell(&f.url));
        let _ = writeln!(out, "#   or: wget -c -O {} {}", shell(&dest), shell(&f.url));
    }
    let sums: Vec<(&str, String)> = files
        .iter()
        .flat_map(|f| {
            let path = f.download();
            [("sha256sum", f.sha256.clone()), ("md5sum", f.md5.clone())]
                .into_iter()
                .filter_map(move |(tool, sum)| Some((tool, format!("{}  {path}", sum?))))
        })
        .collect();
    for tool in ["sha256sum", "md5sum"] {
        let lines: Vec<&str> = sums.iter().filter(|(t, _)| *t == tool).map(|(_, l)| l.as_str()).collect();
        if !lines.is_empty() {
            let _ = writeln!(out, "\n# Checksums ({tool} -c, in {dir_text}):");
            for line in lines {
                let _ = writeln!(out, "# {line}");
            }
        }
    }
    let platform_list = platforms.iter().map(|p| p.as_str()).collect::<Vec<_>>().join(",");
    let source = match rerun.manifest {
        Some(manifest) => format!("--manifest {}", shell(&manifest.to_string_lossy())),
        None => format!("--genome {genome}"),
    };
    let mut setup = format!("lrpipe setup {source} --platform {platform_list} --dir {}", shell(&dir_text));
    if files.iter().any(|f| f.accept) {
        setup.push_str(" --accept-licences");
    }
    if let Some(bonefish) = rerun.bonefish {
        setup.push_str(&format!(" --bonefish {}", shell(&bonefish.to_string_lossy())));
        let _ = writeln!(out, "\n# It builds bonefish's index of the reference there, with that bonefish.");
    }
    let _ = writeln!(out, "\n# Then, where lrpipe runs:\n# {setup}");
    Ok(out)
}

/// Maintainer mode: fetches every file of the manifest at `manifest_path` that has no SHA-256
/// (into `dir`, checked against its size and MD5 where the manifest has them) and writes their
/// SHA-256s into the manifest. Returns the files pinned, by id, with their SHA-256s.
///
/// # Errors
///
/// If the manifest cannot be read, parsed or written, a file needs its terms accepted and
/// `accept_licences` is false, or a file cannot be fetched or does not match.
pub fn pin(
    manifest_path: &Path,
    dir: &Path,
    accept_licences: bool,
    fetcher: &Fetcher,
) -> Result<Vec<(String, String)>> {
    let text = fs::read_to_string(manifest_path).with_context(|| format!("reading {}", manifest_path.display()))?;
    let manifest = Manifest::parse(&text).with_context(|| format!("in {}", manifest_path.display()))?;
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    let _lock = lock(dir)?;
    let files: Vec<&FileEntry> = manifest.files.iter().filter(|f| f.sha256.is_none()).collect();
    let to_accept = terms(&files);
    if !to_accept.is_empty() && !accept_licences {
        bail!(
            "these files' terms must be accepted before they are fetched:\n{to_accept}Give --accept-licences to accept them"
        );
    }
    let mut sums = Vec::new();
    for f in files {
        let dest = dir.join(f.download());
        let expected = Expected { sha256: None, ..expected(f) };
        let digests = match Digests::of(&dest) {
            Ok(d) if d.mismatch(&expected).is_none() => d,
            _ => {
                say(&format!("fetching {} from {}", f.download(), f.url));
                fetcher.fetch(&f.url, &dest, &expected, &say)?
            }
        };
        say(&format!("{}: sha256 {}", f.id, digests.sha256));
        sums.push((f.id.clone(), digests.sha256));
    }
    if !sums.is_empty() {
        let pinned = manifest::pin_text(&text, &sums)?;
        fs::write(manifest_path, pinned).with_context(|| format!("writing {}", manifest_path.display()))?;
    }
    Ok(sums)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalogue_lines_sort_as_sort_does() {
        // `LC_ALL=C sort -k1,1 -k2,2n`: contigs by bytes (chr10 before chr2), starts by number,
        // then whole lines.
        let mut lines: Vec<&[u8]> =
            vec![b"chr2\t5\tb", b"chr10\t100\ta", b"chr2\t40\ta", b"chr1\t7\tz", b"chr1\t7\ta", b"chr1\t10\ta"];
        lines.sort_by(|a, b| sort_key(a).cmp(&sort_key(b)).then_with(|| a.cmp(b)));
        let want: Vec<&[u8]> =
            vec![b"chr1\t7\ta", b"chr1\t7\tz", b"chr1\t10\ta", b"chr10\t100\ta", b"chr2\t5\tb", b"chr2\t40\ta"];
        assert_eq!(lines, want);
    }

    #[test]
    fn a_site_configuration_that_does_not_load_leaves_the_old_one() {
        let dir = std::env::temp_dir().join(format!("lrpipe-site-write-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("site.toml");
        let good = "reference = \"/r.fa\"\n[snv]\n";
        let made = |text: &str| hex(&Sha256::digest(text.as_bytes()));
        assert_eq!(write_if_changed(&path, good, None).expect("written"), Written::Site);
        assert!(write_if_changed(&path, "reference = \"/r.fa\"\n[snv]\nunknown = 1\n", Some(&made(good))).is_err());
        assert_eq!(fs::read_to_string(&path).expect("there"), good, "the old one is untouched");
        assert!(!dir.join("site.toml.previous").exists() && !dir.join("site.toml.tmp").exists());
        let newer = "reference = \"/s.fa\"\n[snv]\n";
        write_if_changed(&path, newer, Some(&made(good))).expect("written");
        assert_eq!(fs::read_to_string(&path).expect("there"), newer);
        assert_eq!(fs::read_to_string(dir.join("site.toml.previous")).expect("kept"), good);
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn a_site_configuration_edited_by_hand_is_kept() {
        // LUX-18: a rerun of setup replaced site.toml whatever was in it.
        let dir = std::env::temp_dir().join(format!("lrpipe-site-edited-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("site.toml");
        let made = "reference = \"/r.fa\"\n[snv]\n";
        let edited = "reference = \"/r.fa\"\n[snv]\ndevice = \"cuda:0\"\n";
        fs::write(&path, edited).expect("written");
        let newer = "reference = \"/s.fa\"\n[snv]\n";
        let sha = hex(&Sha256::digest(made.as_bytes()));
        let written = write_if_changed(&path, newer, Some(&sha)).expect("written");
        assert_eq!(written, Written::Beside(dir.join("site.toml.new")));
        assert_eq!(fs::read_to_string(&path).expect("there"), edited, "the edited one is kept");
        assert_eq!(fs::read_to_string(dir.join("site.toml.new")).expect("beside it"), newer);
        // Not known to be a setup's (from before setups recorded theirs): kept too.
        assert!(matches!(write_if_changed(&path, newer, None).expect("written"), Written::Beside(_)));
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn shell_words_are_quoted() {
        assert_eq!(shell("/ref/a.fa"), "/ref/a.fa");
        assert_eq!(shell("https://x/y?download=1"), "'https://x/y?download=1'");
        assert_eq!(shell("it's"), r"'it'\''s'");
    }
}

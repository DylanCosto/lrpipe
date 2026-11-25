//! Fetching a file: over HTTP(S) with ureq (rustls, no OpenSSL), or from a `file://` URL by
//! copying. Bytes go to a `.part` file that a later fetch resumes from (an HTTP Range request,
//! or a seek), and the whole is checked against the manifest's size and checksums before it is
//! renamed into place; one that does not match is fetched again from the start.
//!
//! What identifies the version being fetched (the server's `ETag` or `Last-Modified`, or a local
//! file's size and time) is kept beside the `.part` file, in `.part.validator`. A resume sends it
//! as `If-Range`, so a server whose file changed sends the whole new one, and a local file that
//! changed is copied again from the start. A `.part` without one is resumed only for a file the
//! manifest has a checksum of, which catches a mixture; one without a checksum starts again.

use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant, UNIX_EPOCH};

use anyhow::{Context, Result, anyhow, bail};
use md5::Md5;
use serde::{Deserialize, Serialize};
use sha2::Sha256;

use crate::provenance::hex;

/// Times a file is fetched whole before a checksum mismatch is an error.
pub const ATTEMPTS: usize = 3;

/// Times in a row a transfer cut off by the network (or a server error) is resumed without a
/// byte arriving in between; one that brings bytes starts the count again.
pub const RESUMES: usize = 5;

/// What a fetched file must be.
#[derive(Clone, Copy, Debug, Default)]
pub struct Expected<'a> {
    /// Its size in bytes.
    pub size: Option<u64>,
    /// Its SHA-256 (lowercase hex).
    pub sha256: Option<&'a str>,
    /// Its MD5 (lowercase hex).
    pub md5: Option<&'a str>,
}

impl Expected<'_> {
    /// Whether there is a checksum to check against.
    #[must_use]
    pub fn pinned(&self) -> bool {
        self.sha256.is_some() || self.md5.is_some()
    }
}

/// A file's size and checksums.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Digests {
    /// Bytes.
    pub size: u64,
    /// SHA-256 (lowercase hex).
    pub sha256: String,
    /// MD5 (lowercase hex).
    pub md5: String,
}

impl Digests {
    /// Reads `path` and computes its size and checksums.
    ///
    /// # Errors
    ///
    /// If the file cannot be read.
    pub fn of(path: &Path) -> io::Result<Self> {
        use md5::Digest as _;
        use sha2::Digest as _;
        let mut file = File::open(path)?;
        let (mut sha256, mut md5) = (Sha256::new(), Md5::new());
        let mut buf = vec![0; 1 << 20];
        let mut size = 0;
        loop {
            let n = file.read(&mut buf)?;
            if n == 0 {
                break;
            }
            sha256.update(&buf[..n]);
            md5.update(&buf[..n]);
            size += n as u64;
        }
        Ok(Self { size, sha256: hex(&sha256.finalize()), md5: hex(&md5.finalize()) })
    }

    /// What differs from `expected`, if anything.
    #[must_use]
    pub fn mismatch(&self, expected: &Expected<'_>) -> Option<String> {
        if let Some(size) = expected.size.filter(|&s| s != self.size) {
            return Some(format!("{} bytes, not {size}", self.size));
        }
        for (name, want, got) in [("sha256", expected.sha256, &self.sha256), ("md5", expected.md5, &self.md5)] {
            if let Some(want) = want.filter(|w| w != got) {
                return Some(format!("{name} {got}, not {want}"));
            }
        }
        None
    }
}

/// Fetches files.
#[derive(Debug)]
pub struct Fetcher {
    agent: ureq::Agent,
    /// Pause before resuming a cut-off transfer, times the resumes in a row.
    pause: Duration,
    /// Time between progress lines.
    progress: Duration,
}

impl Default for Fetcher {
    fn default() -> Self {
        Self::new()
    }
}

/// Why a transfer stopped.
#[derive(Debug)]
enum Failure {
    /// Worth resuming: the network, a timeout, or a server error.
    Retry(anyhow::Error),
    /// Not worth it: the URL, the host's name or certificate, the server's refusal (a 4xx,
    /// rate limits and quotas included), or a local file that cannot be written.
    Stop(anyhow::Error),
}

/// Whether a request that failed with `e` is worth trying again: a dropped connection, a
/// timeout, a malformed response or a server error, not a name that does not resolve, a TLS or
/// certificate failure, a bad URL or proxy, or a refusal (4xx, including 429's rate limits).
#[must_use]
pub fn retryable(e: &ureq::Error) -> bool {
    match e {
        ureq::Error::StatusCode(code) => *code >= 500 || *code == 408,
        ureq::Error::Io(_) | ureq::Error::Timeout(_) | ureq::Error::ConnectionFailed | ureq::Error::Protocol(_) => true,
        _ => false,
    }
}

impl Fetcher {
    /// A fetcher with lrpipe's user agent, `$HTTPS_PROXY`/`$HTTP_PROXY` honoured, timeouts for
    /// connecting and for each read and none for a whole transfer (a reference takes minutes),
    /// 2 s between resumes and a progress line every 10 s.
    #[must_use]
    pub fn new() -> Self {
        let config = ureq::Agent::config_builder()
            .user_agent(concat!("lrpipe/", env!("CARGO_PKG_VERSION")))
            .timeout_connect(Some(Duration::from_secs(30)))
            .timeout_per_read(Some(Duration::from_secs(120)))
            .build();
        Self { agent: config.into(), pause: Duration::from_secs(2), progress: Duration::from_secs(10) }
    }

    /// The same, pausing `pause` (times the resumes in a row) before resuming a cut-off
    /// transfer.
    #[must_use]
    pub fn with_pause(self, pause: Duration) -> Self {
        Self { pause, ..self }
    }

    /// The same, with `every` between progress lines.
    #[must_use]
    pub fn with_progress(self, every: Duration) -> Self {
        Self { progress: every, ..self }
    }

    /// Fetches `url` to `dest`, by way of `dest.part` (resumed if it is there), and checks it
    /// against `expected`; a file that does not match is fetched again from the start, up to
    /// [`ATTEMPTS`] times in all. Progress and resumes are reported through `say`. Returns the
    /// file's digests.
    ///
    /// # Errors
    ///
    /// If it cannot be fetched, or does not match `expected` after [`ATTEMPTS`] tries.
    pub fn fetch(&self, url: &str, dest: &Path, expected: &Expected<'_>, say: &dyn Fn(&str)) -> Result<Digests> {
        let part = part_path(dest);
        if let Some(dir) = dest.parent() {
            fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
        }
        let label = dest.file_name().map_or_else(String::new, |n| n.to_string_lossy().into_owned());
        for attempt in 1..=ATTEMPTS {
            self.transfer(url, &part, expected, &label, say)?;
            let digests = Digests::of(&part).with_context(|| format!("reading {}", part.display()))?;
            match digests.mismatch(expected) {
                None => {
                    fs::rename(&part, dest)
                        .with_context(|| format!("renaming {} to {}", part.display(), dest.display()))?;
                    let _ = fs::remove_file(validator_path(&part));
                    return Ok(digests);
                }
                Some(why) => {
                    say(&format!("{url}: fetched {why} (attempt {attempt} of {ATTEMPTS}); fetching it again"));
                    start_again(&part)?;
                }
            }
        }
        bail!("{url} does not match the manifest after {ATTEMPTS} tries")
    }

    /// Fetches `url` into `part`, from where `part` ends, resuming a transfer the network cut
    /// off until [`RESUMES`] tries in a row bring nothing.
    fn transfer(&self, url: &str, part: &Path, expected: &Expected<'_>, label: &str, say: &dyn Fn(&str)) -> Result<()> {
        let mut stalls = 0;
        loop {
            let before = length(part)?;
            match expected.size {
                Some(size) if before == size => return Ok(()),
                // Longer than the whole: not a prefix of it.
                Some(size) if before > size => start_again(part)?,
                _ => {}
            }
            if before > 0 && !expected.pinned() && !validator_path(part).exists() {
                // Nothing tells what it is a part of, and no checksum would catch a mixture.
                start_again(part)?;
            }
            let result = if let Some(path) = url.strip_prefix("file://") {
                self.copy_file(Path::new(path), part, expected.size, label, say)
            } else {
                self.get(url, part, expected.size, label, say)
            };
            match result {
                Ok(()) => return Ok(()),
                Err(Failure::Retry(e)) => {
                    if length(part)? > before {
                        stalls = 0;
                    }
                    stalls += 1;
                    if stalls > RESUMES {
                        return Err(e.context(format!("fetching {url}: {RESUMES} resumes in a row brought nothing")));
                    }
                    say(&format!("{url}: {e:#}; resuming ({stalls} of {RESUMES} in a row)"));
                    std::thread::sleep(self.pause * u32::try_from(stalls).unwrap_or(u32::MAX));
                }
                Err(Failure::Stop(e)) => return Err(e.context(format!("fetching {url}"))),
            }
        }
    }

    /// One HTTP(S) request for what `part` lacks: a Range request from its end if it has any
    /// bytes, with `If-Range` and the version it is a part of. A server that answers with the
    /// whole file (200: it ignores ranges, or the file changed) has it written from the start,
    /// its version recorded; one that says the range starts past the end (416) leaves `part` as
    /// it is, for its size and checksum to judge; one that sends another range than asked has
    /// `part` started again.
    fn get(&self, url: &str, part: &Path, size: Option<u64>, label: &str, say: &dyn Fn(&str)) -> Result<(), Failure> {
        let offset = length(part).map_err(Failure::Stop)?;
        let mut request = self.agent.get(url);
        if offset > 0 {
            request = request.header("Range", format!("bytes={offset}-"));
            if let Some(validator) = read_validator(part) {
                request = request.header("If-Range", validator);
            }
        }
        let response = match request.call() {
            Ok(response) => response,
            Err(ureq::Error::StatusCode(416)) if offset > 0 => return Ok(()),
            Err(e) if retryable(&e) => return Err(Failure::Retry(e.into())),
            Err(e) => return Err(Failure::Stop(e.into())),
        };
        let header = |name: &str| response.headers().get(name).and_then(|v| v.to_str().ok()).map(str::to_owned);
        let validator = header("etag").or_else(|| header("last-modified"));
        let (append, start) = match response.status().as_u16() {
            206 => {
                let range = header("content-range").unwrap_or_default();
                if content_range_start(&range) != Some(offset) {
                    start_again(part).map_err(Failure::Stop)?;
                    return Err(Failure::Retry(anyhow!(
                        "asked for bytes from {offset}, sent `{range}`; starting again"
                    )));
                }
                // A server that ignored If-Range, sending part of another version.
                if let (Some(old), Some(new)) = (read_validator(part), &validator)
                    && old != *new
                {
                    start_again(part).map_err(Failure::Stop)?;
                    return Err(Failure::Retry(anyhow!("the file changed ({old}, now {new}); starting again")));
                }
                (true, offset)
            }
            200 => {
                if offset > 0 {
                    say(&format!("{label}: the server sent the whole file; starting again"));
                }
                write_validator(part, validator.as_deref()).map_err(Failure::Stop)?;
                (false, 0)
            }
            code => return Err(Failure::Stop(anyhow!("the server answered {code}"))),
        };
        let total = size.or_else(|| response.body().content_length().map(|n| n + start));
        let mut body = response.into_body().into_reader();
        self.write_body(&mut body, part, append, start, total, label, say)
    }

    /// Copies what `part` lacks of the local file `src` onto its end; a `part` longer than
    /// `src`, or of another version of it (by its size and time), is started again.
    fn copy_file(
        &self,
        src: &Path,
        part: &Path,
        size: Option<u64>,
        label: &str,
        say: &dyn Fn(&str),
    ) -> Result<(), Failure> {
        let stop = |e: io::Error, what: &str, path: &Path| {
            Failure::Stop(anyhow::Error::new(e).context(format!("{what} {}", path.display())))
        };
        let mut from = File::open(src).map_err(|e| stop(e, "opening", src))?;
        let meta = from.metadata().map_err(|e| stop(e, "reading", src))?;
        let since = meta.modified().ok().and_then(|t| t.duration_since(UNIX_EPOCH).ok()).unwrap_or_default();
        let version = format!("file {} {}.{:09}", meta.len(), since.as_secs(), since.subsec_nanos());
        let mut offset = length(part).map_err(Failure::Stop)?;
        if offset > meta.len() || (offset > 0 && read_validator(part).is_some_and(|v| v != version)) {
            start_again(part).map_err(Failure::Stop)?;
            offset = 0;
        }
        write_validator(part, Some(&version)).map_err(Failure::Stop)?;
        from.seek(SeekFrom::Start(offset)).map_err(|e| stop(e, "reading", src))?;
        self.write_body(&mut from, part, true, offset, size.or(Some(meta.len())), label, say)
    }

    /// Writes `body` to `part`, appended (at `start`) or from the start, saying how far it has
    /// got every so often. What arrived before an error stays in `part`, for the next transfer
    /// to resume from. An error reading `body` is worth retrying; one writing `part` (a full
    /// disk, a quota) is not.
    #[expect(clippy::too_many_arguments, reason = "one transfer's facts")]
    fn write_body(
        &self,
        body: &mut dyn Read,
        part: &Path,
        append: bool,
        start: u64,
        total: Option<u64>,
        label: &str,
        say: &dyn Fn(&str),
    ) -> Result<(), Failure> {
        let writing =
            |e: io::Error| Failure::Stop(anyhow::Error::new(e).context(format!("writing {}", part.display())));
        let file =
            OpenOptions::new().create(true).append(append).write(true).truncate(!append).open(part).map_err(writing)?;
        let mut out = BufWriter::with_capacity(1 << 20, file);
        let mut buf = vec![0; 1 << 16];
        let (mut done, mut said) = (start, Instant::now());
        let read = loop {
            let n = match body.read(&mut buf) {
                Ok(0) => break Ok(()),
                Ok(n) => n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => break Err(e),
            };
            out.write_all(&buf[..n]).map_err(writing)?;
            done += n as u64;
            if said.elapsed() >= self.progress {
                said = Instant::now();
                say(&match total {
                    Some(total) if total > 0 => {
                        format!("{label}: {} of {} ({}%)", bytes(done), bytes(total), done * 100 / total)
                    }
                    _ => format!("{label}: {}", bytes(done)),
                });
            }
        };
        out.flush().map_err(writing)?;
        read.map_err(|e| Failure::Retry(e.into()))
    }
}

/// `n` bytes for people.
#[must_use]
pub fn bytes(n: u64) -> String {
    let x = n as f64;
    match n {
        1 => "1 byte".to_owned(),
        0..1_000 => format!("{n} bytes"),
        1_000..1_000_000 => format!("{:.1} kB", x / 1e3),
        1_000_000..1_000_000_000 => format!("{:.1} MB", x / 1e6),
        _ => format!("{:.1} GB", x / 1e9),
    }
}

/// `dest` with `.part` added: where it is fetched to.
#[must_use]
pub fn part_path(dest: &Path) -> PathBuf {
    with_suffix(dest, ".part")
}

/// Where the version a `.part` file is a part of is kept.
#[must_use]
pub fn validator_path(part: &Path) -> PathBuf {
    with_suffix(part, ".validator")
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// The size of `part`, 0 if it is not there.
fn length(part: &Path) -> Result<u64> {
    match fs::metadata(part) {
        Ok(m) => Ok(m.len()),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(0),
        Err(e) => Err(e).with_context(|| format!("reading {}", part.display())),
    }
}

/// Removes `part` and the record of its version.
fn start_again(part: &Path) -> Result<()> {
    for path in [part.to_path_buf(), validator_path(part)] {
        match fs::remove_file(&path) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => {
                return Err(e).with_context(|| format!("removing {}", path.display()));
            }
            _ => {}
        }
    }
    Ok(())
}

fn read_validator(part: &Path) -> Option<String> {
    fs::read_to_string(validator_path(part)).ok().map(|v| v.trim().to_owned()).filter(|v| !v.is_empty())
}

/// Records the version `part` is a part of, or that it is not known.
fn write_validator(part: &Path, validator: Option<&str>) -> Result<()> {
    let path = validator_path(part);
    match validator {
        Some(v) => fs::write(&path, v).with_context(|| format!("writing {}", path.display())),
        None => match fs::remove_file(&path) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => {
                Err(e).with_context(|| format!("removing {}", path.display()))
            }
            _ => Ok(()),
        },
    }
}

/// The first byte of a `Content-Range: bytes START-END/TOTAL` header.
fn content_range_start(header: &str) -> Option<u64> {
    header.strip_prefix("bytes ")?.split('-').next()?.trim().parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_sizes_are_bytes() {
        // A small file was "0.0 kB".
        assert_eq!(bytes(0), "0 bytes");
        assert_eq!(bytes(1), "1 byte");
        assert_eq!(bytes(34), "34 bytes");
        assert_eq!(bytes(300_000), "300.0 kB");
        assert_eq!(bytes(2_500_000_000), "2.5 GB");
    }

    #[test]
    fn content_range_gives_its_start() {
        assert_eq!(content_range_start("bytes 100-199/200"), Some(100));
        assert_eq!(content_range_start("bytes */200"), None);
        assert_eq!(content_range_start(""), None);
    }

    #[test]
    fn names_certificates_and_refusals_are_not_retried() {
        for e in [
            ureq::Error::HostNotFound,
            ureq::Error::Tls("bad certificate"),
            ureq::Error::BadUri("x".into()),
            ureq::Error::StatusCode(404),
            ureq::Error::StatusCode(403),
            ureq::Error::StatusCode(429),
            ureq::Error::InvalidProxyUrl,
        ] {
            assert!(!retryable(&e), "{e:?}");
        }
        for e in [
            ureq::Error::Io(io::Error::from(io::ErrorKind::ConnectionReset)),
            ureq::Error::ConnectionFailed,
            ureq::Error::StatusCode(503),
            ureq::Error::StatusCode(408),
        ] {
            assert!(retryable(&e), "{e:?}");
        }
    }

    #[test]
    fn digests_find_what_differs() {
        let dir = std::env::temp_dir().join(format!("lrpipe-digests-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("abc");
        fs::write(&path, b"abc").expect("written");
        let d = Digests::of(&path).expect("digests");
        assert_eq!(d.sha256, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        assert_eq!(d.md5, "900150983cd24fb0d6963f7d28e17f72");
        assert_eq!(d.mismatch(&Expected { size: Some(3), sha256: Some(&d.sha256), md5: Some(&d.md5) }), None);
        assert!(d.mismatch(&Expected { size: Some(4), ..Expected::default() }).is_some());
        assert!(
            d.mismatch(&Expected { md5: Some("00000000000000000000000000000000"), ..Expected::default() }).is_some()
        );
        let _ = fs::remove_dir_all(dir);
    }
}

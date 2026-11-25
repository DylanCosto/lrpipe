//! What a run was made from: files by size, time and (when small enough) SHA-256, and tool
//! binaries by path, version and SHA-256.

use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::sync::{LazyLock, Mutex, PoisonError};
use std::time::{SystemTime, UNIX_EPOCH};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// Files larger than this are recorded by size and time only (a reference or the reads).
pub const HASH_LIMIT: u64 = 1 << 30;

/// A file as it was when the run used it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FileRecord {
    /// Where it is.
    pub path: PathBuf,
    /// Its size in bytes.
    pub bytes: u64,
    /// When it was last changed (UTC).
    pub modified: String,
    /// SHA-256, for files up to [`HASH_LIMIT`].
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

impl FileRecord {
    /// Records `path`.
    ///
    /// # Errors
    ///
    /// If the file cannot be read.
    pub fn of(path: &Path) -> io::Result<Self> {
        let meta = fs::metadata(path)?;
        let sha256 = if meta.len() <= HASH_LIMIT { Some(sha256(path)?) } else { None };
        Ok(Self { path: path.to_path_buf(), bytes: meta.len(), modified: utc(meta.modified()?), sha256 })
    }

    /// [`Self::of`], read once in this process while the file keeps its size and time: a run
    /// records the site's files several times (checked before aligning, recorded with the align
    /// stage, then with the others), and each up to [`HASH_LIMIT`] would be read whole each time.
    ///
    /// # Errors
    ///
    /// If the file cannot be read.
    pub fn cached(path: &Path) -> io::Result<Self> {
        type Key = (PathBuf, u64, SystemTime);
        static RECORDS: LazyLock<Mutex<HashMap<Key, FileRecord>>> = LazyLock::new(Mutex::default);
        let meta = fs::metadata(path)?;
        let key = (path.to_path_buf(), meta.len(), meta.modified()?);
        let records = || RECORDS.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(record) = records().get(&key) {
            return Ok(record.clone());
        }
        let record = Self::of(path)?;
        records().insert(key, record.clone());
        Ok(record)
    }
}

/// A tool binary as it was when the run used it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ToolRecord {
    /// The binary, with `PATH` resolved.
    pub path: PathBuf,
    /// The first line it prints for `--version`.
    pub version: String,
    /// SHA-256 of the binary.
    pub sha256: String,
}

impl ToolRecord {
    /// Records the tool at `program`, a path or a name on `PATH`.
    ///
    /// # Errors
    ///
    /// If it cannot be found or read, or does not run with `--version`.
    pub fn of(program: &Path) -> io::Result<Self> {
        let path = resolve(program)
            .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("{} not found", program.display())))?;
        let out = Command::new(&path).arg("--version").output()?;
        if !out.status.success() {
            return Err(io::Error::other(format!("{} --version failed: {}", path.display(), out.status)));
        }
        let version = String::from_utf8_lossy(&out.stdout).lines().next().unwrap_or_default().trim().to_owned();
        let sha256 = sha256(&path)?;
        Ok(Self { path, version, sha256 })
    }
}

/// The executable `program` names: itself if it has a directory part, else the first match
/// on `PATH`.
#[must_use]
pub fn resolve(program: &Path) -> Option<PathBuf> {
    if program.components().count() > 1 {
        return program.is_file().then(|| program.to_path_buf());
    }
    env::split_paths(&env::var_os("PATH")?).map(|dir| dir.join(program)).find(|p| p.is_file())
}

/// `path` spelled the same way however it was given (LRUN-07): its directory with symbolic links,
/// `.` and `..` resolved, as far as it exists (the rest taken as written), and its own name as
/// given. The name is not followed, so a link to a file keeps the files beside it (a reference's
/// `.fai`). A relative path, or one whose directory cannot be resolved, is returned as is.
#[must_use]
pub fn canonical(path: &Path) -> PathBuf {
    let (Some(dir), Some(name)) = (path.parent(), path.file_name()) else { return path.to_path_buf() };
    if !path.is_absolute() {
        return path.to_path_buf();
    }
    for existing in dir.ancestors() {
        let Ok(mut out) = fs::canonicalize(existing) else { continue };
        let rest = dir.strip_prefix(existing).expect("an ancestor of the directory");
        for part in rest.components() {
            match part {
                Component::ParentDir => {
                    out.pop();
                }
                Component::Normal(part) => out.push(part),
                Component::RootDir | Component::CurDir | Component::Prefix(_) => {}
            }
        }
        out.push(name);
        return out;
    }
    path.to_path_buf()
}

/// SHA-256 of a file, in hex.
///
/// # Errors
///
/// If the file cannot be read.
pub fn sha256(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = vec![0; 1 << 20];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex(&hasher.finalize()))
}

/// Lowercase hex of `bytes`.
#[must_use]
pub fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write as _;
    bytes.iter().fold(String::with_capacity(bytes.len() * 2), |mut s, b| {
        let _ = write!(s, "{b:02x}");
        s
    })
}

/// `time` as `YYYY-MM-DDTHH:MM:SSZ`.
#[must_use]
pub fn utc(time: SystemTime) -> String {
    let secs = time.duration_since(UNIX_EPOCH).map_or(0, |d| d.as_secs());
    let (days, rem) = (secs / 86_400, secs % 86_400);
    // Howard Hinnant's days-to-civil algorithm, for days since 1970-01-01.
    let z = i64::try_from(days).unwrap_or(0) + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z", rem / 3600, rem % 3600 / 60, rem % 60)
}

#[cfg(test)]
mod tests {

    #[test]
    fn a_file_is_read_once_while_it_keeps_its_size_and_time() {
        let path = std::env::temp_dir().join(format!("lrpipe-cached-{}", std::process::id()));
        fs::write(&path, "first").expect("written");
        let first = FileRecord::cached(&path).expect("recorded");
        let modified = fs::metadata(&path).and_then(|m| m.modified()).expect("its time");
        // Changed in place, its size and time kept: taken from the record, not read again.
        fs::write(&path, "other").expect("written");
        File::options().write(true).open(&path).and_then(|f| f.set_modified(modified)).expect("its time kept");
        assert_eq!(FileRecord::cached(&path).expect("recorded"), first);
        assert_ne!(FileRecord::of(&path).expect("read").sha256, first.sha256);
        // Another time: read again.
        fs::write(&path, "third").expect("written");
        File::options()
            .write(true)
            .open(&path)
            .and_then(|f| f.set_modified(modified + std::time::Duration::from_secs(1)))
            .expect("time");
        assert_ne!(FileRecord::cached(&path).expect("recorded").sha256, first.sha256);
        let _ = fs::remove_file(path);
    }

    use super::*;
    use std::time::Duration;

    #[test]
    fn utc_formats_known_instants() {
        assert_eq!(utc(UNIX_EPOCH), "1970-01-01T00:00:00Z");
        assert_eq!(utc(UNIX_EPOCH + Duration::from_secs(951_782_400)), "2000-02-29T00:00:00Z");
        assert_eq!(utc(UNIX_EPOCH + Duration::from_secs(1_793_232_245)), "2026-10-29T00:04:05Z");
    }

    #[test]
    fn resolves_names_on_path() {
        let sh = resolve(Path::new("sh")).expect("sh is on PATH");
        assert!(sh.is_absolute() && sh.ends_with("sh"));
        assert_eq!(resolve(Path::new("/no/such/tool")), None);
    }

    #[test]
    fn hashes_files() {
        let dir = env::temp_dir().join(format!("lrpipe-hash-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("temp dir");
        let f = dir.join("abc");
        fs::write(&f, "abc").expect("written");
        assert_eq!(sha256(&f).expect("hashed"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
        let _ = fs::remove_dir_all(dir);
    }
}

//! How far a tool has got, for `lrpipe run`'s display: where [`ENV`] names a file, a tool writes
//! its fraction done and what it is doing there as one line, `0.413\tpileup 261/631 chunks`,
//! replaced as it goes. Without it (a tool run on its own) nothing is written.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::{Duration, Instant};

/// The environment variable naming the file a tool's progress goes to.
pub const ENV: &str = "LRPIPE_PROGRESS";

/// How often the file is written at most: the display reads it a few times a second.
const INTERVAL: Duration = Duration::from_millis(500);

/// Where a tool's progress goes, if anywhere.
#[derive(Debug)]
pub struct Progress {
    path: Option<PathBuf>,
    /// When the file was last written, so that threads reporting at once write it only every
    /// [`INTERVAL`].
    last: Mutex<Option<Instant>>,
}

impl Progress {
    /// The file [`ENV`] names, or none.
    #[must_use]
    pub fn from_env() -> Self {
        Self { path: std::env::var_os(ENV).filter(|p| !p.is_empty()).map(PathBuf::from), last: Mutex::new(None) }
    }

    /// Records `fraction` (0 to 1) done, doing `what`, unless the file was written less than
    /// [`INTERVAL`] ago: for what is reported often, as each of many jobs ends. A fraction of 1
    /// is always written. A file that cannot be written is left: the display is not worth
    /// failing a tool for.
    pub fn report(&self, fraction: f64, what: &str) {
        self.write(fraction, what, fraction < 1.0);
    }

    /// Records `fraction` done, doing `what`, however recently the file was written: for a step
    /// that starts or ends, which a later report may not follow for a while.
    pub fn mark(&self, fraction: f64, what: &str) {
        self.write(fraction, what, false);
    }

    fn write(&self, fraction: f64, what: &str, throttled: bool) {
        let Some(path) = &self.path else { return };
        let mut last = self.last.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        if throttled && last.is_some_and(|t| t.elapsed() < INTERVAL) {
            return;
        }
        *last = Some(Instant::now());
        let line = format!("{:.4}\t{}\n", fraction.clamp(0.0, 1.0), what.replace(['\t', '\n'], " "));
        // Renamed into place, so the display never reads half a line.
        let tmp = path.with_extension("tmp");
        if fs::write(&tmp, line).is_ok() {
            let _ = fs::rename(&tmp, path);
        }
    }

    /// Whether the progress goes anywhere.
    #[must_use]
    pub fn is_on(&self) -> bool {
        self.path.is_some()
    }
}

/// The fraction done and what is being done, as [`Progress::report`] wrote them to `path`, if it
/// has.
#[must_use]
pub fn read(path: &Path) -> Option<(f64, String)> {
    let text = fs::read_to_string(path).ok()?;
    let (fraction, what) = text.trim_end_matches('\n').split_once('\t')?;
    let fraction: f64 = fraction.parse().ok()?;
    fraction.is_finite().then(|| (fraction.clamp(0.0, 1.0), what.to_owned()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_report_is_read_back_and_none_is_written_without_a_file() {
        let path = std::env::temp_dir().join(format!("lrpipe-progress-{}", std::process::id()));
        let progress = Progress { path: Some(path.clone()), last: Mutex::new(None) };
        progress.report(0.25, "pileup\t5/20 chunks");
        assert_eq!(read(&path), Some((0.25, "pileup 5/20 chunks".to_owned())));
        progress.report(0.5, "too soon");
        assert_eq!(read(&path).map(|r| r.0), Some(0.25), "written at most every {INTERVAL:?}");
        progress.mark(0.5, "a step");
        assert_eq!(read(&path), Some((0.5, "a step".to_owned())), "a mark is always written");
        progress.report(1.0, "done");
        assert_eq!(read(&path), Some((1.0, "done".to_owned())), "the end is always written");
        for (garbage, read_as) in [
            ("nan\tx\n", None),
            ("inf\tx\n", None),
            ("-1\tx\n", Some(0.0)),
            ("2.5\tx\n", Some(1.0)),
            ("0.5 no tab\n", None),
            ("", None),
            ("\u{fffd}\t\u{fffd}", None),
        ] {
            fs::write(&path, garbage).expect("written");
            assert_eq!(read(&path).map(|r| r.0), read_as, "{garbage:?}");
        }
        let _ = fs::remove_file(&path);
        let off = Progress { path: None, last: Mutex::new(None) };
        assert!(!off.is_on());
        off.report(0.5, "nowhere");
        assert_eq!(read(&path), None);
    }
}

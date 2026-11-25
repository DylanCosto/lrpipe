//! A child process whose CPU time and peak memory are read when it exits, and which another
//! thread may kill meanwhile without ever signalling a reused process ID.
//!
//! The standard library's `Child::wait` drops the resource usage the kernel reports for the
//! child, so the child is reaped here with `wait4`. It is first waited for with `WNOWAIT`,
//! which leaves it a zombie: its ID stays reserved until it is reaped under the same lock that
//! `kill` takes.
#![expect(unsafe_code, reason = "waitid and wait4 have no safe binding that reports resource usage")]

use std::io;
use std::mem::MaybeUninit;
use std::os::unix::process::{CommandExt, ExitStatusExt};
use std::process::{Child, Command, ExitStatus};
use std::sync::Mutex;
use std::time::Duration;

/// Resources a child used.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Usage {
    /// CPU time in user mode.
    pub user: Duration,
    /// CPU time in the kernel.
    pub system: Duration,
    /// Peak resident memory in KiB.
    pub max_rss_kib: u64,
}

/// How a child ended.
#[derive(Clone, Copy, Debug)]
pub struct Exit {
    /// Its exit status.
    pub status: ExitStatus,
    /// What it used.
    pub usage: Usage,
}

/// A running child.
#[derive(Debug)]
pub struct Running {
    pid: libc::pid_t,
    /// Whether the child leads a process group of its own, which [`Running::kill`] then kills
    /// whole.
    group: bool,
    /// `None` once reaped.
    child: Mutex<Option<Child>>,
}

impl Running {
    /// Spawns `command` as a child that the kernel kills when the thread spawning it ends,
    /// lrpipe being killed included (a stage outliving its run would hold the machine and
    /// race a rerun for its files). That is Linux's `PR_SET_PDEATHSIG`; elsewhere a stage
    /// outlives an lrpipe that is killed outright (stopping it with SIGINT or SIGTERM still
    /// stops its stages).
    ///
    /// Positive `nice_delta` lowers the child's CPU priority relative to lrpipe's. The child
    /// runs in its own process group: a terminal's Ctrl-C (SIGINT to the terminal's foreground process
    /// group) then reaches lrpipe alone, which stops its stages itself and records them as
    /// stopped, rather than each stage dying of the signal on its own and being taken for a
    /// failure. [`Running::kill`] kills the group, so whatever the child started goes too.
    ///
    /// # Errors
    ///
    /// If the program cannot be started, or lrpipe ended before the child could ask to end
    /// with it.
    pub fn spawn(command: &mut Command, nice_delta: i32) -> io::Result<Self> {
        // SAFETY: getpid is async-signal-safe.
        let parent = unsafe { libc::getpid() };
        let child_nice = (nice_delta > 0).then(|| {
            // SAFETY: PRIO_PROCESS with ID 0 reads the calling thread's inherited niceness.
            let inherited = unsafe { libc::getpriority(libc::PRIO_PROCESS, 0) };
            inherited.saturating_add(nice_delta).min(19)
        });
        let hook = move || {
            // SAFETY: prctl, getppid, setpgid and setpriority are async-signal-safe, and touch no
            // memory of ours.
            unsafe {
                #[cfg(target_os = "linux")]
                if libc::prctl(libc::PR_SET_PDEATHSIG, libc::SIGKILL) != 0 {
                    return Err(io::Error::last_os_error());
                }
                // The parent may have gone before the request took hold.
                if libc::getppid() != parent {
                    return Err(io::Error::other("lrpipe ended while starting a stage"));
                }
                // Lowering a process's own priority needs no privilege; a failure leaves it as is.
                if let Some(priority) = child_nice {
                    libc::setpriority(libc::PRIO_PROCESS, 0, priority);
                }
                if libc::setpgid(0, 0) != 0 {
                    return Err(io::Error::last_os_error());
                }
            }
            Ok(())
        };
        // SAFETY: the hook runs between fork and exec, calling only async-signal-safe
        // functions and allocating nothing (the error above is a constant message).
        let child = unsafe { command.pre_exec(hook) }.spawn()?;
        Ok(Self { group: true, ..Self::new(child) })
    }

    /// Takes over a spawned child, which nothing else may wait for.
    ///
    /// # Panics
    ///
    /// If the child's ID does not fit a `pid_t`, which no Unix allows.
    #[must_use]
    pub fn new(child: Child) -> Self {
        let pid = libc::pid_t::try_from(child.id()).expect("process IDs fit pid_t");
        Self { pid, group: false, child: Mutex::new(Some(child)) }
    }

    /// The child's process ID.
    #[must_use]
    pub fn id(&self) -> u32 {
        self.pid.unsigned_abs()
    }

    /// Blocks until the child exits, then kills what is left of its process group when it leads
    /// one (a background process it started and did not wait for, which would outlive it and
    /// hold the machine) and reaps it.
    ///
    /// # Errors
    ///
    /// If waiting fails other than by interruption.
    ///
    /// # Panics
    ///
    /// If the lock is poisoned, or the child was already reaped.
    pub fn wait(&self) -> io::Result<Exit> {
        loop {
            let mut info = MaybeUninit::<libc::siginfo_t>::zeroed();
            // SAFETY: `info` is a valid, writable siginfo_t, and WNOWAIT leaves the child for
            // wait4 below.
            let r = unsafe {
                libc::waitid(libc::P_PID, self.pid as libc::id_t, info.as_mut_ptr(), libc::WEXITED | libc::WNOWAIT)
            };
            if r == 0 {
                break;
            }
            let e = io::Error::last_os_error();
            if e.kind() != io::ErrorKind::Interrupted {
                return Err(e);
            }
        }
        let mut child = self.child.lock().expect("lock poisoned");
        assert!(child.is_some(), "process {} reaped twice", self.pid);
        if self.group {
            // SAFETY: kill has no memory effects. The child is a zombie not yet reaped, so its
            // ID, the group's, is reserved: it signals what is left of the group alone.
            unsafe {
                libc::kill(-self.pid, libc::SIGKILL);
            }
        }
        let mut status = 0;
        let mut usage = MaybeUninit::<libc::rusage>::zeroed();
        loop {
            // SAFETY: `status` and `usage` are valid and writable; the child has exited and
            // is reaped here only, while the lock keeps `kill` from signalling its ID.
            let r = unsafe { libc::wait4(self.pid, &mut status, 0, usage.as_mut_ptr()) };
            if r == self.pid {
                break;
            }
            let e = io::Error::last_os_error();
            if e.kind() != io::ErrorKind::Interrupted {
                return Err(e);
            }
        }
        *child = None;
        // SAFETY: wait4 succeeded, so it filled `usage`.
        let usage = unsafe { usage.assume_init() };
        Ok(Exit { status: ExitStatus::from_raw(status), usage: Usage::from(&usage) })
    }

    /// Kills the child if it has not been reaped, and the rest of its process group when it
    /// leads one ([`Running::spawn`]'s).
    ///
    /// # Panics
    ///
    /// If the lock is poisoned.
    pub fn kill(&self) {
        if let Some(child) = self.child.lock().expect("lock poisoned").as_mut() {
            if self.group {
                // SAFETY: kill has no memory effects. The child is not reaped while the lock is
                // held, so its ID, which is the group's, is not anyone else's.
                unsafe {
                    libc::kill(-self.pid, libc::SIGKILL);
                }
            }
            // An error means the child has already exited, which is what killing it is for.
            let _ = child.kill();
        }
    }
}

/// Whether `signal` is ignored, as `nohup` leaves SIGHUP for what it runs: a handler installed
/// for it would undo that.
#[must_use]
pub fn ignored(signal: libc::c_int) -> bool {
    let mut old = MaybeUninit::<libc::sigaction>::zeroed();
    // SAFETY: a null new action only reads the current one into `old`, which is valid and
    // writable.
    let r = unsafe { libc::sigaction(signal, std::ptr::null(), old.as_mut_ptr()) };
    // SAFETY: sigaction succeeded, so it filled `old`.
    r == 0 && unsafe { old.assume_init() }.sa_sigaction == libc::SIG_IGN
}

/// `ru_maxrss` per KiB: Linux reports KiB, macOS bytes.
#[cfg(target_os = "macos")]
const MAX_RSS_UNIT: u64 = 1024;
#[cfg(not(target_os = "macos"))]
const MAX_RSS_UNIT: u64 = 1;

fn duration(t: libc::timeval) -> Duration {
    let secs = u64::try_from(t.tv_sec).unwrap_or(0);
    let micros = u32::try_from(t.tv_usec).unwrap_or(0);
    Duration::from_secs(secs) + Duration::from_micros(micros.into())
}

impl From<&libc::rusage> for Usage {
    fn from(r: &libc::rusage) -> Self {
        Self {
            user: duration(r.ru_utime),
            system: duration(r.ru_stime),
            max_rss_kib: u64::try_from(r.ru_maxrss).unwrap_or(0) / MAX_RSS_UNIT,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::process::Command;
    use std::sync::Arc;
    use std::thread;

    #[test]
    fn reports_status_and_usage() {
        let child = Command::new("sh").args(["-c", "exit 3"]).spawn().expect("sh runs");
        let exit = Running::new(child).wait().expect("waits");
        assert_eq!(exit.status.code(), Some(3));
        assert!(exit.usage.max_rss_kib > 0, "{:?}", exit.usage);
    }

    #[test]
    fn kill_from_another_thread_ends_the_child() {
        let child = Command::new("sleep").arg("30").spawn().expect("sleep runs");
        let running = Arc::new(Running::new(child));
        let waiter = {
            let running = Arc::clone(&running);
            thread::spawn(move || running.wait())
        };
        running.kill();
        let exit = waiter.join().expect("waiter finishes").expect("waits");
        assert_eq!(exit.status.signal(), Some(libc::SIGKILL));
        running.kill(); // after reaping: does nothing
    }
}

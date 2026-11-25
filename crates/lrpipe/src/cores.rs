//! Keeping a run to `--threads` cores in all: lrpipe-snv and the read pass run at the same
//! time, each with that many threads, on the same cores.
//!
//! On Linux lrpipe confines itself to the first cores it may use, and every thread and program
//! it starts afterwards inherits that, so the stages share them: each has them all while the
//! other waits (lrpipe-snv on the GPU) or is done, as on a machine with that many cores.
//! Elsewhere the two split the threads between them.

use std::io;
use std::num::NonZero;

/// Where a run's stages run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cores {
    /// Every core lrpipe may use.
    All,
    /// This many of them, which the stages share.
    Confined(usize),
    /// Fewer than lrpipe may use, but it cannot confine itself to them here.
    Unconfinable,
}

/// Confines the calling thread, and the threads and programs it starts afterwards, to the first
/// `n` cores it may use (at least one), unless it may use no more than that.
///
/// # Errors
///
/// If the cores cannot be read or set.
#[cfg(target_os = "linux")]
#[expect(unsafe_code, reason = "sched_getaffinity and sched_setaffinity have no safe binding")]
pub fn confine(n: usize) -> io::Result<Cores> {
    let n = n.max(1);
    let size = size_of::<libc::cpu_set_t>();
    let allowed = allowed()?;
    if n >= allowed.len() {
        return Ok(Cores::All);
    }
    // SAFETY: as above.
    unsafe {
        let mut set: libc::cpu_set_t = std::mem::zeroed();
        for &cpu in &allowed[..n] {
            libc::CPU_SET(cpu, &mut set);
        }
        if libc::sched_setaffinity(0, size, &raw const set) != 0 {
            return Err(io::Error::last_os_error());
        }
    }
    Ok(Cores::Confined(n))
}

/// The cores the calling thread may run on (its affinity). [`available`] can be fewer: it also
/// counts a container's CPU quota.
///
/// # Errors
///
/// If they cannot be read.
#[cfg(target_os = "linux")]
#[expect(unsafe_code, reason = "sched_getaffinity has no safe binding")]
fn allowed() -> io::Result<Vec<usize>> {
    // SAFETY: a zeroed cpu_set_t is an empty set, and the call is given its size.
    unsafe {
        let mut set: libc::cpu_set_t = std::mem::zeroed();
        if libc::sched_getaffinity(0, size_of::<libc::cpu_set_t>(), &raw mut set) != 0 {
            return Err(io::Error::last_os_error());
        }
        Ok((0..libc::CPU_SETSIZE as usize).filter(|&cpu| libc::CPU_ISSET(cpu, &set)).collect())
    }
}

/// Whether `n` cores are fewer than lrpipe may use: it cannot confine itself here.
///
/// # Errors
///
/// None here; they are those of the Linux version.
#[cfg(not(target_os = "linux"))]
pub fn confine(n: usize) -> io::Result<Cores> {
    let all = std::thread::available_parallelism().map_or(usize::MAX, NonZero::get);
    Ok(if n.max(1) >= all { Cores::All } else { Cores::Unconfinable })
}

/// The cores lrpipe may use.
#[must_use]
pub fn available() -> usize {
    std::thread::available_parallelism().map_or(8, NonZero::get)
}

#[cfg(all(test, target_os = "linux"))]
mod tests {
    use super::*;

    #[test]
    fn confining_leaves_the_thread_that_many_cores() -> io::Result<()> {
        // Run on a thread of its own, which alone is confined.
        std::thread::spawn(|| -> io::Result<()> {
            // Under a CPU quota (a container) `available` counts fewer cores than the affinity
            // allows, and confining to them is right.
            let (all, allowed) = (available(), allowed()?.len());
            assert_eq!(confine(all)?, if all >= allowed { Cores::All } else { Cores::Confined(all) });
            if all > 1 {
                assert_eq!(confine(1)?, Cores::Confined(1));
                assert_eq!(available(), 1);
                let child = std::thread::spawn(available).join().expect("the thread ran");
                assert_eq!(child, 1, "a thread started afterwards inherits the cores");
            }
            Ok(())
        })
        .join()
        .expect("the thread ran")
    }
}

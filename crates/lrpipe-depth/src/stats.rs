//! Summary statistics and depth distributions, as mosdepth computes and prints them.

use std::io::{self, Write};

/// Nim's `formatFloat(x, ffDecimal, precision)`, which is C's `%#.*f`: correctly rounded,
/// ties to even, and a trailing `.` when `precision` is 0.
pub fn format_float(out: &mut String, x: f64, precision: usize) {
    use std::fmt::Write as _;
    if x.is_nan() {
        out.push_str("nan");
    } else if x.is_infinite() {
        out.push_str(if x > 0.0 { "inf" } else { "-inf" });
    } else {
        let _ = write!(out, "{x:.precision$}");
        if precision == 0 {
            out.push('.');
        }
    }
}

pub fn fmt_float(x: f64, precision: usize) -> String {
    let mut s = String::new();
    format_float(&mut s, x, precision);
    s
}

/// mosdepth's `depth_stat` (depthstat.nim:3-57).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DepthStat {
    pub cum_length: u64,
    pub cum_depth: u64,
    pub min_depth: u32,
    pub max_depth: u32,
}

impl Default for DepthStat {
    fn default() -> Self {
        DepthStat { cum_length: 0, cum_depth: 0, min_depth: u32::MAX, max_depth: 0 }
    }
}

impl DepthStat {
    /// `newDepthStat`. Depths are reinterpreted as unsigned, as Nim's `uint32` conversion does.
    pub fn of(depths: &[i32]) -> Self {
        let mut s = DepthStat { cum_length: depths.len() as u64, ..Default::default() };
        for &d in depths {
            s.cum_depth = s.cum_depth.wrapping_add(i64::from(d) as u64);
            s.min_depth = s.min_depth.min(d as u32);
            s.max_depth = s.max_depth.max(d as u32);
        }
        s
    }

    pub fn add(&mut self, o: &DepthStat) {
        self.cum_length += o.cum_length;
        self.cum_depth = self.cum_depth.wrapping_add(o.cum_depth);
        self.min_depth = self.min_depth.min(o.min_depth);
        self.max_depth = self.max_depth.max(o.max_depth);
    }

    /// One `*.mosdepth.summary.txt` line (mosdepth.nim:460-480).
    pub fn summary_line(&self, name: &str, precision: usize) -> String {
        let mean = if self.cum_length > 0 { self.cum_depth as f64 / self.cum_length as f64 } else { 0.0 };
        let min = if self.min_depth == u32::MAX { 0 } else { self.min_depth };
        format!(
            "{name}\t{}\t{}\t{}\t{min}\t{}\n",
            self.cum_length,
            self.cum_depth,
            fmt_float(mean, precision),
            self.max_depth
        )
    }
}

pub const SUMMARY_HEADER: &str = "chrom\tlength\tbases\tmean\tmin\tmax\n";

/// A histogram for medians (depthstat.nim:9-35): depths above `size - 1` count as `size - 1`.
pub struct CountStat {
    counts: Vec<u32>,
    n: usize,
}

impl CountStat {
    pub fn new(size: usize) -> Self {
        CountStat { counts: vec![0; size], n: 0 }
    }

    #[inline]
    pub fn add(&mut self, v: i32) {
        assert!(v >= 0, "error setting negative depth value:{v}");
        self.n += 1;
        let i = (v as usize).min(self.counts.len() - 1);
        self.counts[i] = self.counts[i].wrapping_add(1);
    }

    pub fn median(&self) -> i64 {
        let stop_n = (0.5 + self.n as f64 * 0.5) as u64;
        let mut cum = 0u64;
        for (i, &c) in self.counts.iter().enumerate() {
            cum += u64::from(c);
            if cum >= stop_n {
                return i as i64;
            }
        }
        -1
    }

    pub fn clear(&mut self) {
        if self.n == 0 {
            return;
        }
        self.n = 0;
        self.counts.fill(0);
    }
}

/// A depth histogram for `*.dist.txt`.
#[derive(Clone, Debug)]
pub struct Distribution {
    counts: Vec<i64>,
}

const MAX_COVERAGE: i32 = 400_000;

impl Default for Distribution {
    fn default() -> Self {
        Distribution { counts: vec![0; 512] }
    }
}

impl Distribution {
    /// Counts each depth in `depths` (mosdepth.nim:417-437).
    pub fn add_depths(&mut self, depths: &[i32]) {
        for &d in depths {
            let v = if d > MAX_COVERAGE { MAX_COVERAGE - 10 } else { d };
            if v < 0 {
                continue;
            }
            let v = v as usize;
            if v >= self.counts.len() - 1 {
                self.counts.resize(v + 10, 0);
            }
            self.counts[v] += 1;
        }
    }

    /// Counts one value (window means, rounded).
    pub fn add_capped(&mut self, i: i64) {
        let i = i.clamp(0, self.counts.len() as i64 - 1) as usize;
        self.counts[i] += 1;
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.counts.iter().all(|&c| c == 0)
    }

    pub fn merge(&mut self, o: &Distribution) {
        if o.counts.len() > self.counts.len() {
            self.counts.resize(o.counts.len(), 0);
        }
        for (a, b) in self.counts.iter_mut().zip(&o.counts) {
            *a += b;
        }
    }

    /// Lines of "name, depth, proportion of bases with at least that depth"
    /// (mosdepth.nim:439-458).
    pub fn write(&self, name: &str, precision: usize, out: &mut impl Write) -> io::Result<()> {
        let sum: i64 = self.counts.iter().sum();
        if sum < 1 {
            return Ok(());
        }
        let mut cum = 0.0f64;
        let mut line = String::new();
        for (irev, &v) in self.counts.iter().enumerate().rev() {
            if irev > 300 && v == 0 {
                continue;
            }
            cum += v as f64 / sum as f64;
            if cum < 8e-5 {
                continue;
            }
            line.clear();
            line.push_str(name);
            line.push('\t');
            line.push_str(&irev.to_string());
            line.push('\t');
            format_float(&mut line, cum, precision);
            line.push('\n');
            out.write_all(line.as_bytes())?;
        }
        Ok(())
    }
}

/// Nim's `toInt` for floats: round half away from zero, computed as Nim does.
pub fn nim_to_int(x: f64) -> i64 {
    if x >= 0.0 { (x + 0.5) as i64 } else { (x - 0.5) as i64 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_like_printf() {
        assert_eq!(fmt_float(0.125, 2), "0.12");
        assert_eq!(fmt_float(0.375, 2), "0.38");
        assert_eq!(fmt_float(2.5, 0), "2.");
        assert_eq!(fmt_float(1234.565, 2), "1234.57");
        assert_eq!(fmt_float(0.0, 2), "0.00");
    }

    #[test]
    fn median_like_mosdepth() {
        let mut c = CountStat::new(32768);
        for i in 0..=10 {
            c.add(i);
        }
        assert_eq!(c.median(), 5);
        for _ in 0..4 {
            c.add(6);
        }
        assert_eq!(c.median(), 6);
        c.clear();
        assert_eq!(c.median(), 0);
    }

    #[test]
    fn distribution_skips_tiny_tail() {
        let mut d = Distribution::default();
        d.add_depths(&[0, 1, 1, 2, 500_000]);
        let mut out = Vec::new();
        d.write("c", 2, &mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.starts_with("c\t399990\t0.20\n"), "{s}");
        assert!(s.ends_with("c\t2\t0.40\nc\t1\t0.80\nc\t0\t1.00\n"), "{s}");
    }
}

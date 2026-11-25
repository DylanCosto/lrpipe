//! Number parsing and formatting as htslib does it (`textutils_internal.h`,
//! `kstring.c`) and as glibc's `strtod`/`printf("%g")` do.

/// `hts_str2int(in, &end, bits, &failed)`: (value, end, failed). No leading
/// whitespace is skipped.
pub fn str2int(s: &[u8], bits: u32) -> (i64, usize, bool) {
    let mut n: u64 = 0;
    let mut limit: u64 = (1u64 << (bits - 1)) - 1;
    let mut fast = (bits - 1) * 1000 / 3322 + 1;
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let mut v = 0;
    let neg = match at(0) {
        b'-' => {
            limit += 1;
            v = 1;
            true
        }
        b'+' => {
            v = 1;
            false
        }
        _ => false,
    };
    loop {
        fast -= 1;
        if fast == 0 || !at(v).is_ascii_digit() {
            break;
        }
        n = n * 10 + (at(v) - b'0') as u64;
        v += 1;
    }
    let mut failed = false;
    if at(v) >= b'0' && fast == 0 {
        let (limit_d_10, limit_m_10) = (limit / 10, limit % 10);
        while at(v).is_ascii_digit() {
            let d = (at(v) - b'0') as u64;
            if n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10) {
                n = n * 10 + d;
                v += 1;
            } else {
                while at(v).is_ascii_digit() {
                    v += 1;
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    let value = if neg { (n as i64).wrapping_neg() } else { n as i64 };
    (value, v, failed)
}

/// `hts_str2uint(in, &end, bits, &failed)`
pub fn str2uint(s: &[u8], bits: u32) -> (u64, usize, bool) {
    let mut n: u64 = 0;
    let limit: u64 = if bits < 64 { (1u64 << bits) - 1 } else { u64::MAX };
    let mut fast = bits * 1000 / 3322 + 1;
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let mut v = 0;
    if at(0) == b'+' {
        v = 1;
    }
    loop {
        fast -= 1;
        if fast == 0 || !at(v).is_ascii_digit() {
            break;
        }
        // wraps like C for the huge `bits` vcf_parse passes for GT
        n = n.wrapping_mul(10).wrapping_add((at(v) - b'0') as u64);
        v += 1;
    }
    let mut failed = false;
    if at(v).is_ascii_digit() && fast == 0 {
        let (limit_d_10, limit_m_10) = (limit / 10, limit % 10);
        while at(v).is_ascii_digit() {
            let d = (at(v) - b'0') as u64;
            if n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10) {
                n = n * 10 + d;
                v += 1;
            } else {
                while at(v).is_ascii_digit() {
                    v += 1;
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    (n, v, failed)
}

fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// glibc `strtod`: (value, end); end 0 if nothing was converted.
pub fn strtod(s: &[u8]) -> (f64, usize) {
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let mut i = 0;
    while is_c_space(at(i)) {
        i += 1;
    }
    let mut neg = false;
    if at(i) == b'+' || at(i) == b'-' {
        neg = at(i) == b'-';
        i += 1;
    }
    let sign = |v: f64| if neg { -v } else { v };
    let rest = &s[i.min(s.len())..];
    let lower: Vec<u8> = rest.iter().take(8).map(|c| c.to_ascii_lowercase()).collect();
    if lower.starts_with(b"inf") {
        let n = if lower.starts_with(b"infinity") { 8 } else { 3 };
        return (sign(f64::INFINITY), i + n);
    }
    if lower.starts_with(b"nan") {
        let mut end = i + 3;
        if at(end) == b'(' {
            let mut j = end + 1;
            while at(j).is_ascii_alphanumeric() || at(j) == b'_' {
                j += 1;
            }
            if at(j) == b')' {
                end = j + 1;
            }
        }
        return (sign(f64::NAN), end);
    }
    if at(i) == b'0' && (at(i + 1) == b'x' || at(i + 1) == b'X') {
        if let Some((v, n)) = parse_hex_float(&s[i + 2..]) {
            return (sign(v), i + 2 + n);
        }
        // "0x" without hex digits converts the "0"
        return (sign(0.0), i + 1);
    }
    let start = i;
    let mut j = i;
    let mut digits = 0;
    while at(j).is_ascii_digit() {
        j += 1;
        digits += 1;
    }
    if at(j) == b'.' {
        j += 1;
        while at(j).is_ascii_digit() {
            j += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return (0.0, 0);
    }
    if at(j) == b'e' || at(j) == b'E' {
        let mut k = j + 1;
        if at(k) == b'+' || at(k) == b'-' {
            k += 1;
        }
        if at(k).is_ascii_digit() {
            while at(k).is_ascii_digit() {
                k += 1;
            }
            j = k;
        }
    }
    let text = std::str::from_utf8(&s[start..j]).unwrap_or("0");
    let mut v: f64 = text.parse().unwrap_or(0.0);
    if text.starts_with('.') && v == 0.0 {
        v = format!("0{text}").parse().unwrap_or(0.0);
    }
    (sign(v), j)
}

/// Hexadecimal floating point (after "0x"): (value, bytes used).
fn parse_hex_float(s: &[u8]) -> Option<(f64, usize)> {
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let mut i = 0;
    let mut mantissa: u64 = 0;
    let mut exp: i64 = 0;
    let mut any = false;
    let mut overflow_bits = false;
    while at(i).is_ascii_hexdigit() {
        let d = (at(i) as char).to_digit(16).unwrap() as u64;
        if mantissa >> 60 == 0 {
            mantissa = mantissa << 4 | d;
        } else {
            exp += 4;
            overflow_bits |= d != 0;
        }
        any = true;
        i += 1;
    }
    if at(i) == b'.' {
        i += 1;
        while at(i).is_ascii_hexdigit() {
            let d = (at(i) as char).to_digit(16).unwrap() as u64;
            if mantissa >> 60 == 0 {
                mantissa = mantissa << 4 | d;
                exp -= 4;
            }
            any = true;
            i += 1;
        }
    }
    if !any {
        return None;
    }
    if at(i) == b'p' || at(i) == b'P' {
        let mut k = i + 1;
        let mut eneg = false;
        if at(k) == b'+' || at(k) == b'-' {
            eneg = at(k) == b'-';
            k += 1;
        }
        if at(k).is_ascii_digit() {
            let mut e: i64 = 0;
            while at(k).is_ascii_digit() {
                e = (e * 10 + (at(k) - b'0') as i64).min(100_000);
                k += 1;
            }
            exp += if eneg { -e } else { e };
            i = k;
        }
    }
    let _ = overflow_bits;
    Some(((mantissa as f64) * 2f64.powi(exp.clamp(-20000, 20000) as i32), i))
}

/// C `atof`
pub fn atof(s: &[u8]) -> f64 {
    strtod(s).0
}

/// `hts_str2dbl(in, &end, &failed)`: (value, end, failed).
pub fn str2dbl(s: &[u8]) -> (f64, usize, bool) {
    const D: [f64; 22] = [
        1.0, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18,
        1e19, 1e20,
    ];
    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let fallback = || {
        let (d, end) = strtod(s);
        (d, end, end == 0)
    };
    let mut v = 0;
    while is_c_space(at(v)) {
        v += 1;
    }
    let mut neg = false;
    if at(v) == b'-' {
        neg = true;
        v += 1;
    } else if at(v) == b'+' {
        v += 1;
    }
    match at(v) {
        b'1'..=b'9' => {}
        b'0' if at(v + 1) != b'x' && at(v + 1) != b'X' => {}
        _ => return fallback(),
    }
    while at(v) == b'0' {
        v += 1;
    }
    let start = v;
    let mut n: u64 = 0;
    let mut max_len = 15;
    let mut point: Option<usize> = None;
    loop {
        max_len -= 1;
        if max_len == 0 || !at(v).is_ascii_digit() {
            break;
        }
        n = n * 10 + (at(v) - b'0') as u64;
        v += 1;
    }
    if max_len != 0 && at(v) == b'.' {
        point = Some(v - start);
        v += 1;
        loop {
            max_len -= 1;
            if max_len == 0 || !at(v).is_ascii_digit() {
                break;
            }
            n = n * 10 + (at(v) - b'0') as u64;
            v += 1;
        }
    }
    let point = point.unwrap_or(v - start);
    if max_len == 0 || at(v) == b'e' || at(v) == b'E' {
        return fallback();
    }
    let d = n as f64 / D[v - start - point];
    (if neg { -d } else { d }, v, false)
}

/// `kputw`
pub fn put_int(v: i64, out: &mut Vec<u8>) {
    out.extend_from_slice(v.to_string().as_bytes());
}

/// glibc `printf("%g", d)` for a finite, positive or negative non-zero
/// value, or inf/nan.
fn format_g(d: f64) -> String {
    if d.is_nan() {
        return if d.is_sign_negative() { "-nan".into() } else { "nan".into() };
    }
    if d.is_infinite() {
        return if d < 0.0 { "-inf".into() } else { "inf".into() };
    }
    if d == 0.0 {
        return if d.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    // Rust's exact formatting rounds correctly, like glibc
    let sci = format!("{:.5e}", d);
    let (mant, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent");
    let trim = |s: &str| -> String {
        if s.contains('.') {
            let t = s.trim_end_matches('0');
            t.trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if !(-4..6).contains(&exp) {
        let m = trim(mant);
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", exp.abs())
    } else {
        let decimals = (5 - exp).max(0) as usize;
        trim(&format!("{:.*}", decimals, d))
    }
}

/// htslib's `kputd` (a custom `%g`).
pub fn put_double(mut d: f64, out: &mut Vec<u8>) {
    if d == 0.0 {
        out.extend_from_slice(if d.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    if d < 0.0 {
        out.push(b'-');
        d = -d;
    }
    if !(0.0001..=999999.0).contains(&d) {
        out.extend_from_slice(format_g(d).as_bytes());
        return;
    }
    let (i, shift) = if d < 0.001 {
        ((d * 1_000_000_000.0).round_ties_even(), 1)
    } else if d < 0.01 {
        ((d * 100_000_000.0).round_ties_even(), 2)
    } else if d < 0.1 {
        ((d * 10_000_000.0).round_ties_even(), 3)
    } else if d < 1.0 {
        ((d * 1_000_000.0).round_ties_even(), 4)
    } else if d < 10.0 {
        ((d * 100_000.0).round_ties_even(), 5)
    } else if d < 100.0 {
        ((d * 10_000.0).round_ties_even(), 6)
    } else if d < 1000.0 {
        ((d * 1000.0).round_ties_even(), 7)
    } else if d < 10000.0 {
        ((d * 100.0).round_ties_even(), 8)
    } else if d < 100000.0 {
        ((d * 10.0).round_ties_even(), 9)
    } else {
        (d.round_ties_even(), 10)
    };
    let mut i = i as u32;
    // buf[0..21]; cp starts at 20 - shift; digits written backwards
    let mut buf = [0u8; 21];
    let mut cp: usize = 20 - shift;
    for k in 0..3 {
        let two = i % 100;
        cp -= 2;
        buf[cp] = b'0' + (two / 10) as u8;
        buf[cp + 1] = b'0' + (two % 10) as u8;
        if k < 2 {
            i /= 100;
        }
    }
    if i >= 100 {
        cp -= 1;
        buf[cp] = b'0' + (i / 100) as u8;
    }
    let mut p = 20 - cp;
    let ep: usize;
    if p <= 10 {
        ep = cp + 5;
        while p < 10 {
            cp -= 1;
            buf[cp] = b'0';
            p += 1;
        }
        cp -= 1;
        buf[cp] = b'.';
        cp -= 1;
        buf[cp] = b'0';
    } else {
        cp -= 1;
        let mut xp = cp;
        ep = cp + 6;
        while p > 10 {
            buf[xp] = buf[xp + 1];
            xp += 1;
            p -= 1;
        }
        buf[xp] = b'.';
    }
    let mut ep = ep;
    while buf[ep] == b'0' && ep > cp {
        ep -= 1;
    }
    if buf[ep] != 0 && buf[ep] != b'.' {
        ep += 1;
    }
    out.extend_from_slice(&buf[cp..ep]);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kputd(d: f64) -> String {
        let mut v = Vec::new();
        put_double(d, &mut v);
        String::from_utf8(v).unwrap()
    }

    #[test]
    fn kputd_matches_htslib() {
        assert_eq!(kputd(0.0), "0");
        assert_eq!(kputd(-0.0), "-0");
        assert_eq!(kputd(1.0), "1");
        assert_eq!(kputd(0.5), "0.5");
        assert_eq!(kputd(123.456), "123.456");
        assert_eq!(kputd(99.99999), "100");
        assert_eq!(kputd(0.009999999), "0.01");
        assert_eq!(kputd(0.1234567), "0.123457");
        assert_eq!(kputd(1e-5), "1e-05");
        assert_eq!(kputd(1234567.0), "1.23457e+06");
        assert_eq!(kputd(f32::MAX as f64), "3.40282e+38");
        assert_eq!(kputd(-2.5), "-2.5");
        assert_eq!(kputd(f64::INFINITY), "inf");
    }

    #[test]
    fn str2dbl_fast_path() {
        assert_eq!(str2dbl(b"1.5,"), (1.5, 3, false));
        assert_eq!(str2dbl(b"-0.25"), (-0.25, 5, false));
        assert_eq!(str2dbl(b"1e3").0, 1000.0);
        assert!(str2dbl(b"x").2);
    }
}

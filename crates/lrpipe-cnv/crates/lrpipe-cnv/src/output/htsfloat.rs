//! How htslib 1.18 prints a VCF float (`kputd`, kstring.c), which is how Spectre's VCF floats
//! read after `bcftools sort` rewrites the file.

/// The text htslib writes for a float FORMAT value given as `text`: parsed to the nearest
/// double, stored as float32, printed by [`kputd`].
pub fn vcf_float(text: &str) -> String {
    let v: f64 = text.parse().unwrap_or(f64::NAN);
    kputd(f64::from(v as f32))
}

/// htslib's `kputd`: up to 6 significant digits, trailing zeros dropped; C's `%g` outside
/// 0.0001..=999999.
pub fn kputd(d: f64) -> String {
    if d == 0.0 {
        return if d.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    let mut out = String::new();
    let mut d = d;
    if d < 0.0 {
        out.push('-');
        d = -d;
    }
    if !(0.0001..=999_999.0).contains(&d) {
        out.push_str(&c_g(d));
        return out;
    }
    // `uint64_t i = d*10000000000LL`, plus a rounding correction at the 7th digit.
    let mut i = (d * 10_000_000_000.0) as u64;
    i += match d {
        d if d < 0.001 => 5,
        d if d < 0.01 => 50,
        d if d < 0.1 => 500,
        d if d < 1.0 => 5000,
        d if d < 10.0 => 50_000,
        d if d < 100.0 => 500_000,
        d if d < 1000.0 => 5_000_000,
        d if d < 10_000.0 => 50_000_000,
        d if d < 100_000.0 => 500_000_000,
        _ => 5_000_000_000,
    };
    let digits = i.to_string().into_bytes();
    let p = digits.len();
    // `cp` is the string from its first character; `ep` the last character kept before
    // trailing zeros are culled.
    let (mut s, ep) = if p <= 10 {
        // d < 1: "0." + zeros + the first 6 digits.
        let mut s = b"0.".to_vec();
        s.extend(std::iter::repeat_n(b'0', 10 - p));
        s.extend_from_slice(&digits[..6]);
        let ep = s.len() - 1;
        (s, ep)
    } else {
        let int = p - 10;
        let mut s = digits[..int].to_vec();
        s.push(b'.');
        s.extend_from_slice(&digits[int..]);
        s.truncate(7);
        if s[6] == b'.' {
            s.truncate(6);
        }
        (s, 6)
    };
    // Cull trailing zeros, only when the number has a decimal point.
    let mut ep = ep.min(s.len());
    if ep < s.len() {
        while ep > 0 && s[ep] == b'0' {
            ep -= 1;
        }
        let z = ep + 1;
        let mut e = ep;
        while e > 0 {
            if s[e] == b'.' {
                if s[z - 1] == b'.' {
                    s.truncate(z - 1);
                } else {
                    s.truncate(z);
                }
                break;
            }
            e -= 1;
        }
    }
    out.push_str(std::str::from_utf8(&s).expect("ASCII digits"));
    out
}

/// C's `printf("%g", d)` for finite positive `d` (and `nan`/`inf`): 6 significant digits,
/// positional when the exponent is in -4..6, trailing zeros removed.
fn c_g(d: f64) -> String {
    if d.is_nan() {
        return "nan".into();
    }
    if d.is_infinite() {
        return "inf".into();
    }
    let sci = format!("{d:.5e}");
    let (mantissa, exp) = sci.split_once('e').expect("{:e} has an exponent");
    let exp: i32 = exp.parse().expect("integer exponent");
    if (-4..6).contains(&exp) {
        let fixed = format!("{d:.*}", (5 - exp) as usize);
        strip_zeros(&fixed).to_owned()
    } else {
        format!("{}e{}{:02}", strip_zeros(mantissa), if exp < 0 { '-' } else { '+' }, exp.abs())
    }
}

fn strip_zeros(s: &str) -> &str {
    if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.') } else { s }
}

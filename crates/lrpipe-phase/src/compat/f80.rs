//! The x87 80-bit extended precision `long double` of x86-64 C/C++ code:
//! a 64-bit significand, rounding to nearest (ties to even) after every
//! operation. WhatsHap's genotyping DP computes in `long double`, and the
//! last bits of its results are visible in the output (GL of near-certain
//! genotypes), so they are emulated exactly. The exponent range is assumed
//! to stay far from the limits (no subnormal or overflowing extended values).

use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, Copy, Debug, PartialEq)]
enum Kind {
    Finite,
    Inf,
    Nan,
}

/// An 80-bit extended precision value: `mant * 2^(exp - 63)`, with the top
/// bit of `mant` set (or `mant == 0` for zero).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F80 {
    kind: Kind,
    neg: bool,
    exp: i32,
    mant: u64,
}

impl F80 {
    pub const ZERO: F80 = F80 { kind: Kind::Finite, neg: false, exp: 0, mant: 0 };
    pub const ONE: F80 = F80 { kind: Kind::Finite, neg: false, exp: 0, mant: 1 << 63 };
    const NAN: F80 = F80 { kind: Kind::Nan, neg: false, exp: 0, mant: 0 };

    fn inf(neg: bool) -> F80 {
        F80 { kind: Kind::Inf, neg, exp: 0, mant: 0 }
    }

    /// From the in-memory representation (sign and biased exponent, significand)
    pub const fn from_bits(sign_exp: u16, mant: u64) -> F80 {
        if mant == 0 {
            return F80 { kind: Kind::Finite, neg: sign_exp & 0x8000 != 0, exp: 0, mant: 0 };
        }
        F80 { kind: Kind::Finite, neg: sign_exp & 0x8000 != 0, exp: (sign_exp & 0x7fff) as i32 - 16383, mant }
    }

    pub fn is_zero(self) -> bool {
        self.kind == Kind::Finite && self.mant == 0
    }

    /// Round `m * 2^e` (plus a sticky bit for discarded nonzero bits below
    /// `m`) to 64 bits, ties to even.
    fn round(neg: bool, e: i32, m: u128, sticky: bool) -> F80 {
        if m == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        let p = 127 - m.leading_zeros() as i32;
        if p <= 63 {
            debug_assert!(!sticky);
            return F80 { kind: Kind::Finite, neg, exp: e + p, mant: (m << (63 - p)) as u64 };
        }
        let shift = (p - 63) as u32;
        let mut keep = m >> shift;
        let rem = m & ((1u128 << shift) - 1);
        let half = 1u128 << (shift - 1);
        if rem > half || (rem == half && (sticky || keep & 1 == 1)) {
            keep += 1;
        }
        let mut exp = e + p;
        if keep == 1u128 << 64 {
            keep >>= 1;
            exp += 1;
        }
        F80 { kind: Kind::Finite, neg, exp, mant: keep as u64 }
    }

    /// Round `keep` (top bit set) with the 64 bits below it in `rem`
    #[inline]
    fn round_top(neg: bool, mut exp: i32, keep: u64, rem: u64) -> F80 {
        const HALF: u64 = 1 << 63;
        let mut mant = keep;
        if rem > HALF || (rem == HALF && keep & 1 == 1) {
            mant = keep.wrapping_add(1);
            if mant == 0 {
                mant = HALF;
                exp += 1;
            }
        }
        F80 { kind: Kind::Finite, neg, exp, mant }
    }

    pub fn from_f64(x: f64) -> F80 {
        if x.is_nan() {
            return F80::NAN;
        }
        if x.is_infinite() {
            return F80::inf(x < 0.0);
        }
        let bits = x.to_bits();
        let neg = bits >> 63 == 1;
        let e = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & ((1 << 52) - 1);
        if e == 0 {
            if frac == 0 {
                return F80 { neg, ..F80::ZERO };
            }
            // subnormal: frac * 2^-1074
            return F80::round(neg, -1074, frac as u128, false);
        }
        F80 { kind: Kind::Finite, neg, exp: e - 1023, mant: ((1u64 << 52) | frac) << 11 }
    }

    pub fn from_u64(v: u64) -> F80 {
        F80::round(false, 0, v as u128, false)
    }

    /// Conversion to `double`, rounding to nearest (ties to even)
    pub fn to_f64(self) -> f64 {
        match self.kind {
            Kind::Nan => return f64::NAN,
            Kind::Inf => return if self.neg { f64::NEG_INFINITY } else { f64::INFINITY },
            Kind::Finite => {}
        }
        if self.mant == 0 {
            return if self.neg { -0.0 } else { 0.0 };
        }
        let sign = (self.neg as u64) << 63;
        if self.exp > 1023 {
            return f64::from_bits(sign | 0x7ff0_0000_0000_0000);
        }
        // bits to drop below the f64 significand (normal: 11)
        let drop = if self.exp >= -1022 { 11 } else { 11 + (-1022 - self.exp) as u32 };
        let m = self.mant as u128;
        let (mut q, rem, half) = if drop >= 128 {
            (0u128, m, u128::MAX)
        } else {
            (m >> drop, m & ((1u128 << drop) - 1), 1u128 << (drop - 1))
        };
        if rem > half || (rem == half && q & 1 == 1) {
            q += 1;
        }
        if self.exp >= -1022 {
            let mut exp = self.exp;
            if q == 1u128 << 53 {
                q >>= 1;
                exp += 1;
            }
            if exp > 1023 {
                return f64::from_bits(sign | 0x7ff0_0000_0000_0000);
            }
            f64::from_bits(sign | (((exp + 1023) as u64) << 52) | (q as u64 & ((1 << 52) - 1)))
        } else {
            // subnormal (q == 2^52 is the smallest normal number)
            f64::from_bits(sign | q as u64)
        }
    }

    /// Same-sign sum of two finite nonzero values with |x| >= |y|
    #[inline]
    fn add_magnitudes(x: F80, y: F80) -> F80 {
        let d = (x.exp - y.exp) as u32;
        let (xm, ym) = (x.mant, y.mant);
        if d == 0 {
            let (s, _) = xm.overflowing_add(ym);
            // both top bits set: the sum has 65 bits
            return F80::round_top(x.neg, x.exp + 1, (s >> 1) | (1 << 63), s << 63);
        }
        if d < 64 {
            let (y_hi, y_lo) = (ym >> d, ym << (64 - d));
            let (s, carry) = xm.overflowing_add(y_hi);
            if !carry {
                return F80::round_top(x.neg, x.exp, s, y_lo);
            }
            let rem = (s << 63) | (y_lo >> 1) | (y_lo & 1);
            return F80::round_top(x.neg, x.exp + 1, (s >> 1) | (1 << 63), rem);
        }
        if d < 66 {
            let rem = if d == 64 { ym } else { (ym >> 1) | (ym & 1) };
            return F80::round_top(x.neg, x.exp, xm, rem);
        }
        // y is below half an ulp of x
        x
    }

    #[inline]
    fn add_signed(self, other: F80) -> F80 {
        if self.kind == Kind::Finite
            && other.kind == Kind::Finite
            && self.mant != 0
            && other.mant != 0
            && self.neg == other.neg
        {
            return if (self.exp, self.mant) >= (other.exp, other.mant) {
                F80::add_magnitudes(self, other)
            } else {
                F80::add_magnitudes(other, self)
            };
        }
        match (self.kind, other.kind) {
            (Kind::Nan, _) | (_, Kind::Nan) => return F80::NAN,
            (Kind::Inf, Kind::Inf) => return if self.neg == other.neg { self } else { F80::NAN },
            (Kind::Inf, _) => return self,
            (_, Kind::Inf) => return other,
            _ => {}
        }
        if other.mant == 0 {
            return if self.mant == 0 { F80 { neg: self.neg && other.neg, ..F80::ZERO } } else { self };
        }
        if self.mant == 0 {
            return other;
        }
        let (x, y) = if (self.exp, self.mant) >= (other.exp, other.mant) { (self, other) } else { (other, self) };
        let d = (x.exp - y.exp) as u32;
        if x.neg == y.neg {
            // y is below half an ulp of x: the sum rounds to x
            if d >= 66 {
                return x;
            }
            // significands with the top bit at 125, y shifted with the
            // discarded bits jammed into the lowest bit
            let xm = (x.mant as u128) << 62;
            let ym_full = (y.mant as u128) << 62;
            let ym = (ym_full >> d) | u128::from(ym_full & ((1u128 << d) - 1) != 0);
            let sum = xm + ym;
            let shift = if sum >> 126 == 1 { 63 } else { 62 };
            let keep = (sum >> shift) as u64;
            let rem = (sum << (64 - shift)) as u64;
            return F80::round_top(x.neg, x.exp + shift - 62, keep, rem);
        }
        let xm = (x.mant as u128) << 63;
        let ym_full = (y.mant as u128) << 63;
        let (mut ym, lost) = if d >= 127 { (0, true) } else { (ym_full >> d, ym_full & ((1u128 << d) - 1) != 0) };
        if lost {
            // jam the discarded bits into the lowest bit
            ym |= 1;
        }
        let e = x.exp - 126;
        if x.neg == y.neg {
            F80::round(x.neg, e, xm + ym, false)
        } else {
            let diff = xm - ym;
            if diff == 0 {
                return F80::ZERO;
            }
            F80::round(x.neg, e, diff, false)
        }
    }
}

impl Add for F80 {
    type Output = F80;
    #[inline]
    fn add(self, other: F80) -> F80 {
        self.add_signed(other)
    }
}

impl Sub for F80 {
    type Output = F80;
    #[inline]
    fn sub(self, other: F80) -> F80 {
        self.add_signed(F80 { neg: !other.neg, ..other })
    }
}

impl Mul for F80 {
    type Output = F80;
    #[inline]
    fn mul(self, other: F80) -> F80 {
        let neg = self.neg != other.neg;
        match (self.kind, other.kind) {
            (Kind::Nan, _) | (_, Kind::Nan) => return F80::NAN,
            (Kind::Inf, _) | (_, Kind::Inf) => {
                return if self.is_zero() || other.is_zero() { F80::NAN } else { F80::inf(neg) };
            }
            _ => {}
        }
        if self.mant == 0 || other.mant == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        // the product of two normalized significands has its top bit at
        // position 127 or 126
        let m = self.mant as u128 * other.mant as u128;
        let (hi, lo) = ((m >> 64) as u64, m as u64);
        let (keep, rem, exp) = if hi >> 63 == 1 {
            (hi, lo, self.exp + other.exp + 1)
        } else {
            ((hi << 1) | (lo >> 63), lo << 1, self.exp + other.exp)
        };
        F80::round_top(neg, exp, keep, rem)
    }
}

impl Div for F80 {
    type Output = F80;
    #[inline]
    fn div(self, other: F80) -> F80 {
        let neg = self.neg != other.neg;
        match (self.kind, other.kind) {
            (Kind::Nan, _) | (_, Kind::Nan) | (Kind::Inf, Kind::Inf) => return F80::NAN,
            (Kind::Inf, _) => return F80::inf(neg),
            (_, Kind::Inf) => return F80 { neg, ..F80::ZERO },
            _ => {}
        }
        if other.mant == 0 {
            return if self.mant == 0 { F80::NAN } else { F80::inf(neg) };
        }
        if self.mant == 0 {
            return F80 { neg, ..F80::ZERO };
        }
        // a 65-bit quotient, then two more bits by long division and the
        // sticky remainder
        let d = other.mant as u128;
        let n = (self.mant as u128) << 64;
        let (mut q, mut r) = (n / d, n % d);
        for _ in 0..2 {
            r <<= 1;
            q <<= 1;
            if r >= d {
                r -= d;
                q |= 1;
            }
        }
        F80::round(neg, self.exp - other.exp - 66, q, r != 0)
    }
}

/// A finite nonzero divisor with a precomputed reciprocal, for fast exact
/// division (Möller and Granlund, "Improved division by invariant
/// integers", 2011). `x.div_by(&Divisor::new(d)) == x / d`.
#[derive(Clone, Copy, Debug)]
pub struct Divisor {
    value: F80,
    /// floor((2^128 - 1) / mant) - 2^64
    reciprocal: u64,
}

impl Divisor {
    pub fn new(value: F80) -> Divisor {
        let reciprocal = if value.kind == Kind::Finite && value.mant != 0 {
            (u128::MAX / value.mant as u128 - (1u128 << 64)) as u64
        } else {
            0
        };
        Divisor { value, reciprocal }
    }

    /// (u1:u0) / d for u1 < d: (quotient, remainder)
    #[inline]
    fn div_2by1(&self, u1: u64, u0: u64) -> (u64, u64) {
        let d = self.value.mant;
        let q = (self.reciprocal as u128 * u1 as u128).wrapping_add(((u1 as u128 + 1) << 64) | u0 as u128);
        let (mut q1, q0) = ((q >> 64) as u64, q as u64);
        let mut r = u0.wrapping_sub(q1.wrapping_mul(d));
        if r > q0 {
            q1 = q1.wrapping_sub(1);
            r = r.wrapping_add(d);
        }
        if r >= d {
            q1 = q1.wrapping_add(1);
            r -= d;
        }
        (q1, r)
    }
}

impl F80 {
    /// `self / divisor`, identical to `/`
    #[inline]
    pub fn div_by(self, divisor: &Divisor) -> F80 {
        let other = divisor.value;
        if self.kind != Kind::Finite || other.kind != Kind::Finite || self.mant == 0 || other.mant == 0 {
            return self / other;
        }
        let neg = self.neg != other.neg;
        let d = other.mant;
        // 65-bit quotient of (mant << 64) / d
        let (hi, a) = if self.mant >= d { (1u128, self.mant - d) } else { (0, self.mant) };
        let (q_lo, mut r) = divisor.div_2by1(a, 0);
        let mut q = (hi << 64) | q_lo as u128;
        for _ in 0..2 {
            let carry = r >> 63 == 1;
            r <<= 1;
            q <<= 1;
            if carry || r >= d {
                r = r.wrapping_sub(d);
                q |= 1;
            }
        }
        // q has 66 or 67 bits: keep the top 64, the rest (and the sticky
        // remainder) decide the rounding
        let shift = if hi == 1 { 3 } else { 2 };
        let keep = (q >> shift) as u64;
        let rem = ((q << (64 - shift)) as u64) | u64::from(r != 0);
        F80::round_top(neg, self.exp - other.exp - 66 + 63 + shift, keep, rem)
    }
}

impl std::ops::AddAssign for F80 {
    #[inline]
    fn add_assign(&mut self, other: F80) {
        *self = *self + other;
    }
}

impl std::ops::MulAssign for F80 {
    #[inline]
    fn mul_assign(&mut self, other: F80) {
        *self = *self * other;
    }
}

impl std::ops::DivAssign for F80 {
    #[inline]
    fn div_assign(&mut self, other: F80) {
        *self = *self / other;
    }
}

#[cfg(test)]
mod tests {
    use super::F80;

    #[test]
    fn exact_operations() {
        let third = F80::ONE / F80::from_u64(3);
        // 1/3 in extended precision: 0xaaaaaaaaaaaaaaab * 2^-65
        assert_eq!(third, F80::from_bits(0x3ffd, 0xaaaa_aaaa_aaaa_aaab));
        let tenth = F80::ONE / F80::from_u64(10);
        assert_eq!(tenth, F80::from_bits(0x3ffb, 0xcccc_cccc_cccc_cccd));
        assert_eq!((F80::from_f64(0.5) + F80::from_f64(0.25)).to_f64(), 0.75);
        assert_eq!((F80::from_f64(1.0) - F80::from_f64(1.0)).to_f64(), 0.0);
        assert_eq!(F80::from_f64(0.1).to_f64(), 0.1);
        assert_eq!(tenth.to_f64(), 0.1);
        assert_eq!(F80::from_f64(5e-324).to_f64(), 5e-324);
        assert_eq!((F80::from_f64(3.0) * F80::from_f64(7.0)).to_f64(), 21.0);
    }
}

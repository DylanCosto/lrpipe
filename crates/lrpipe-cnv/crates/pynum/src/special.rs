//! scipy's normal distribution tail, from the Cephes `ndtr`, `erf` and `erfc` that scipy
//! 1.10 ships (`scipy/special/cephes/ndtr.c`, Cephes Math Library 2.8, Stephen L. Moshier).
//! The coefficients are copied from that file as written.

#![expect(clippy::excessive_precision, reason = "Cephes coefficients, copied as written")]

/// `scipy.stats.norm.sf(x)`: `ndtr(-x)`.
pub fn norm_sf(x: f64) -> f64 {
    ndtr(-x)
}

/// The standard normal CDF, Cephes `ndtr`.
pub fn ndtr(a: f64) -> f64 {
    if a.is_nan() {
        return f64::NAN;
    }
    let x = a * std::f64::consts::FRAC_1_SQRT_2;
    let z = x.abs();
    if z < std::f64::consts::FRAC_1_SQRT_2 {
        0.5 + 0.5 * erf(x)
    } else {
        let y = 0.5 * erfc(z);
        if x > 0.0 { 1.0 - y } else { y }
    }
}

/// Cephes `MAXLOG`: `ln(DBL_MAX)`.
const MAXLOG: f64 = 7.09782712893383996843E2;

const P: [f64; 9] = [
    2.46196981473530512524E-10,
    5.64189564831068821977E-1,
    7.46321056442269912687E0,
    4.86371970985681366614E1,
    1.96520832956077098242E2,
    5.26445194995477358631E2,
    9.34528527171957607540E2,
    1.02755188689515710272E3,
    5.57535335369399327526E2,
];
const Q: [f64; 8] = [
    1.32281951154744992508E1,
    8.67072140885989742329E1,
    3.54937778887819891062E2,
    9.75708501743205489753E2,
    1.82390916687909736289E3,
    2.24633760818710981792E3,
    1.65666309194161350182E3,
    5.57535340817727675546E2,
];
const R: [f64; 6] = [
    5.64189583547755073984E-1,
    1.27536670759978104416E0,
    5.01905042251180477414E0,
    6.16021097993053585195E0,
    7.40974269950448939160E0,
    2.97886665372100240670E0,
];
const S: [f64; 6] = [
    2.26052863220117276590E0,
    9.39603524938001434673E0,
    1.20489539808096656605E1,
    1.70814450747565897222E1,
    9.60896809063285878198E0,
    3.36907645100081516050E0,
];
const T: [f64; 5] = [
    9.60497373987051638749E0,
    9.00260197203842689217E1,
    2.23200534594684319226E3,
    7.00332514112805075473E3,
    5.55923013010394962768E4,
];
const U: [f64; 5] = [
    3.35617141647503099647E1,
    5.21357949780152679795E2,
    4.59432382970980127987E3,
    2.26290000613890934246E4,
    4.92673942608635921086E4,
];

/// Cephes `polevl`: `coef[0] x^N + ... + coef[N]`, by Horner's rule.
fn polevl(x: f64, coef: &[f64]) -> f64 {
    coef[1..].iter().fold(coef[0], |ans, &c| ans * x + c)
}

/// Cephes `p1evl`: like [`polevl`] with an implied leading coefficient of 1.
fn p1evl(x: f64, coef: &[f64]) -> f64 {
    coef[1..].iter().fold(x + coef[0], |ans, &c| ans * x + c)
}

fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return -erf(-x);
    }
    if x.abs() > 1.0 {
        return 1.0 - erfc(x);
    }
    let z = x * x;
    x * polevl(z, &T) / p1evl(z, &U)
}

fn erfc(a: f64) -> f64 {
    if a.is_nan() {
        return f64::NAN;
    }
    let x = a.abs();
    if x < 1.0 {
        return 1.0 - erf(a);
    }
    let z = -a * a;
    let under = if a < 0.0 { 2.0 } else { 0.0 };
    if z < -MAXLOG {
        return under;
    }
    let z = z.exp();
    let (p, q) = if x < 8.0 { (polevl(x, &P), p1evl(x, &Q)) } else { (polevl(x, &R), p1evl(x, &S)) };
    let mut y = (z * p) / q;
    if a < 0.0 {
        y = 2.0 - y;
    }
    if y == 0.0 { under } else { y }
}

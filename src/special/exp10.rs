/// Base 10 exponential function (Common antilogarithm)
///
/// Returns 10 raised to the x power.
///
/// Range reduction is accomplished by expressing the argument
/// as `10**x = 2**n 10**f`, with `|f| < 0.5 log10(2)`.
/// The Pade' form
///
/// ```text
/// 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
/// ```
///
/// is used to approximate `10**f`.
///
/// # Accuracy
///
/// |   domain    |  # trials |  peak    |    rms |
/// | - | - | - | - |
/// |  -307,+307  |   30000   | 2.2e-16  |  5.5e-17 |
///
/// # Errors
///
/// | condition | value returned |
/// | --- | --- |
/// |  x < -MAXL10 | 0.0 |
/// |  x > MAXL10  | MAXNUM |
///
/// MAXL10 = 308.2547155599167.
#[expect(clippy::cast_possible_truncation)]
#[must_use]
pub fn exp10(mut x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x > MAXL10 {
        return f64::INFINITY;
    }
    if x < -MAXL10 {
        return 0.0f64;
    }
    let mut px = (f64::consts::LOG2_10 * x + 0.5f64).floor();
    let n = px as i16;
    x -= px * LG102A;
    x -= px * LG102B;
    let xx = x * x;
    px = x * evaluate_polynomial(xx, &P);
    x = px / (evaluate_polynomial_1(xx, &Q) - px);
    x = 1.0f64 + mul_exp2(x, 1);
    mul_exp2(x, i32::from(n))
}

const P: [f64; 4] = [
    4.099_625_197_985_870_6E-2_f64,
    1.174_527_325_543_440_5E1_f64,
    4.067_172_899_368_727E2_f64,
    2.394_237_412_073_882_8E3_f64,
];
const Q: [f64; 3] = [
    8.509_361_608_493_066E1_f64,
    1.272_092_711_783_451_3E3_f64,
    2.079_608_192_860_019E3_f64,
];
const LG102A: f64 = 3.010_253_906_25E-1_f64;
const LG102B: f64 = 4.605_038_981_195_214E-6_f64;
const MAXL10: f64 = 308.254_715_559_916_7_f64;

#[test]
fn works() {
    assert_float_eq::assert_float_relative_eq!(exp10(43.3922), 2.467_175_251_851_615_7e43);
}

use super::evaluate_polynomial;
use super::evaluate_polynomial_1;
use super::mul_exp2;
use std::f64;

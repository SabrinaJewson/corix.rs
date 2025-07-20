/// Circular cotangent
///
/// Returns the circular cotangent of the radian argument x.
///
/// Range reduction is modulo pi/4.  A rational function
/// ```text
/// x + x**3 P(x**2)/Q(x**2)
/// ```
/// is employed in the basic interval [0, pi/4].
///
/// # Accuracy
///
/// Relative error:
///
/// | domain | # trials | peak | rms |
/// | - | - | - | - |
/// | +-1.07e9 | 30000 | 2.9e-16 | 8.2e-17 |
///
/// # Errors
///
/// | condition | value returned |
/// | --- | --- |
/// | x > 1.073741824e9 | 0.0 |
/// | x = 0 | INFINITY |
#[must_use]
pub fn cot(x: f64) -> f64 {
    if x == 0.0f64 {
        return f64::INFINITY;
    }
    let (x, negative) = if x < 0.0 { (-x, true) } else { (x, false) };
    if x > LOSSTH {
        return 0.0f64;
    }
    let mut y = (x / f64::consts::FRAC_PI_4).floor();
    let mut z = mul_exp2(y, -3);
    z = z.floor();
    z = y - mul_exp2(z, 3);
    #[expect(clippy::cast_possible_truncation)]
    let mut j = z as i32;
    if j & 1 != 0 {
        j += 1;
        y += 1.0f64;
    }
    z = x - y * DP1 - y * DP2 - y * DP3;
    let zz = z * z;
    if zz > 1.0e-14f64 {
        y = z + z * (zz * evaluate_polynomial(zz, &P) / evaluate_polynomial_1(zz, &Q));
    } else {
        y = z;
    }
    if j & 2 != 0 {
        y = -y;
    } else {
        y = 1.0f64 / y;
    }
    if negative {
        y = -y;
    }
    y
}

const P: [f64; 3] = [
    -1.309_369_391_813_837_9E4_f64,
    1.153_516_648_385_874_2E6_f64,
    -1.795_652_519_764_848_8E7_f64,
];
const Q: [f64; 4] = [
    1.368_129_634_706_929_6E4_f64,
    -1.320_892_344_402_109_7E6_f64,
    2.500_838_018_233_579E7_f64,
    -5.386_957_559_294_546_4E7_f64,
];
const DP1: f64 = 7.853_981_554_508_209E-1_f64;
const DP2: f64 = 7.946_627_356_147_928E-9_f64;
const DP3: f64 = 3.061_616_997_868_383E-17_f64;
const LOSSTH: f64 = 1.073_741_824e9_f64;

#[test]
fn works() {
    assert_float_eq::assert_float_relative_eq!(cot(2.5), -1.338_648_128_304_151_4);
}

use super::evaluate_polynomial;
use super::evaluate_polynomial_1;
use super::mul_exp2;
use std::f64;

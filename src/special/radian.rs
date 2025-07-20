/// Degrees, minutes, seconds to radians:
///
/// 1 arc second, in radians = 4.8481368110953599358991410e-5.
#[must_use]
pub fn radian(d: f64, m: f64, s: f64) -> f64 {
    ((d * 60.0f64 + m) * 60.0f64 + s) * 4.848_136_811_095_36e-5_f64
}

#[test]
fn works() {
    assert_float_eq::assert_float_relative_eq!(radian(12.0, 49.0, 2.0), 2.237_027_287_375_621);
}

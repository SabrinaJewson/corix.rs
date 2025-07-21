unsafe extern "C" {
    /// Evaluates `x * 2 ^ exp`.
    #[link_name = "ldexp"]
    pub safe fn mul_exp2(x: f64, exp: i32) -> f64;

    safe fn tgamma(n: f64) -> f64;

    safe fn lgamma_r(n: f64, s: &mut i32) -> f64;

    /// The error function.
    pub safe fn erf(x: f64) -> f64;

    /// The complementary error function (`1 - erf(x)`).
    pub safe fn erfc(x: f64) -> f64;
}

/// Gamma function.
#[must_use]
pub fn gamma(n: f64) -> f64 {
    tgamma(n)
}

/// Natural logarithm of the gamma function.
///
/// Returns the value and the sign, which is 1 if it’s nonnegative and -1 if it’s negative.
#[must_use]
pub fn ln_gamma(n: f64) -> (f64, i32) {
    let mut s = 0;
    let res = lgamma_r(n, &mut s);
    (res, s)
}

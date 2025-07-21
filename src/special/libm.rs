unsafe extern "C" {
    /// Evaluates `x * 2 ^ exp`.
    ///
    /// # Examples
    ///
    /// ```
    /// assert_float_relative_eq!(corix::mul_exp2(3.5, 3), 28.0);
    /// assert_float_relative_eq!(corix::mul_exp2(240.0, -9), 0.468_75);
    /// # use assert_float_eq::assert_float_relative_eq;
    /// ```
    #[link_name = "ldexp"]
    pub safe fn mul_exp2(x: f64, exp: i32) -> f64;

    safe fn tgamma(n: f64) -> f64;

    safe fn lgamma_r(n: f64, s: &mut i32) -> f64;

    /// The error function.
    ///
    /// See also [`erfc`], its complement (`1 - erf(x)`).
    ///
    /// # Examples
    ///
    /// ```
    /// assert_float_relative_eq!(corix::erf(0.2), 0.222_702_589_210_478_45);
    /// # use assert_float_eq::assert_float_relative_eq;
    /// ```
    pub safe fn erf(x: f64) -> f64;

    /// The complementary error function (`1 - `[`erf(x)`](erf)).
    ///
    /// # Examples
    ///
    /// ```
    /// assert_float_relative_eq!(corix::erfc(0.2), 0.777_297_410_789_521_5);
    /// assert_float_relative_eq!(corix::erfc(0.2) + corix::erf(0.2), 1.0);
    /// # use assert_float_eq::assert_float_relative_eq;
    /// ```
    pub safe fn erfc(x: f64) -> f64;
}

/// Gamma function.
///
/// See also [`ln_gamma`], its logarithm.
///
/// # Examples
///
/// ```
/// // 5! = 120
/// assert_float_relative_eq!(corix::gamma(6.0), 120.0);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn gamma(n: f64) -> f64 {
    tgamma(n)
}

/// Natural logarithm of the [gamma function](gamma).
///
/// Returns a 2-tuple of:
/// - `ln(|Γ(x)|)`, and
/// - 1 if `Γ(x)` is positive and -1 if it is negative.
///
/// # Examples
///
/// ```
/// let (val, sign) = corix::ln_gamma(6.0);
/// assert_float_relative_eq!(val, 4.787_491_742_782_046);
/// assert_eq!(sign, 1);
///
/// let (val, sign) = corix::ln_gamma(-2.5);
/// assert_float_relative_eq!(val, -0.056_243_716_497_674);
/// assert_eq!(sign, -1);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn ln_gamma(n: f64) -> (f64, i32) {
    let mut s = 0;
    let res = lgamma_r(n, &mut s);
    (res, s)
}

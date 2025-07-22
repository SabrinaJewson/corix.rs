/// Evaluate polynomial
///
/// Evaluates a polynomial, with coefficients stored in reverse order:
///
/// ```text
/// y = coef[0] x^(N - 1) + coef[1] x^(N - 2) + … + coef[N - 2] x + coef[N - 1]
/// ```
///
/// The function [`eval_polynomial_1`] adds an additional `x^N` term.
/// Its calling arguments are otherwise the same as `eval_polynomial`.
///
/// # Speed
///
/// In the interest of speed, there are no checks for out-of-bounds arithmetic.
///
/// # Examples
///
/// ```
/// // x² - 3x - 2.5 evaluated at x = 2
/// assert_float_relative_eq!(corix::eval_polynomial(2.0, &[1.0, -3.0, -2.5]), -4.5);
///
/// // The empty polynomial evaluates to zero.
/// assert_float_relative_eq!(corix::eval_polynomial(2.0, &[]), 0.0);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn eval_polynomial(x: f64, coef: &[f64]) -> f64 {
    let mut ans: f64 = 0.0;
    for &coef in coef {
        ans = ans.mul_add(x, coef);
    }
    ans
}

/// Evaluate polynomial when coefficient of x is 1.0.
/// Otherwise same as [`eval_polynomial`].
///
/// # Examples
///
/// ```
/// // x³ + x² - 3x - 2.5 evaluated at x = 2
/// assert_float_relative_eq!(corix::eval_polynomial_1(2.0, &[1.0, -3.0, -2.5]), 3.5);
///
/// // This polynomial `1x⁰` always evaluates to 1.0.
/// assert_float_relative_eq!(corix::eval_polynomial_1(2.0, &[]), 1.0);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn eval_polynomial_1(x: f64, coef: &[f64]) -> f64 {
    let mut ans: f64 = 1.0;
    for &coef in coef {
        ans = ans.mul_add(x, coef);
    }
    ans
}

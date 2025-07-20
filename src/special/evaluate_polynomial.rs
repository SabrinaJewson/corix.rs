/// Evaluate polynomial
///
/// Evaluates polynomial of degree N:
///
/// ```text
///                     2          N
/// y  =  C  + C x + C x  +...+ C x
///        0    1     2          N
/// ```
///
/// Coefficients are stored in reverse order:
///
/// ```text
/// coef[0] = C  , ..., coef[N] = C  .
///            N                   0
/// ```
///
/// The function [`evaluate_polynomial_1`] assumes that `coef[N] = 1.0` and is
/// omitted from the array.  Its calling arguments are
/// otherwise the same as `evaluate_polynomial`.
///
/// # Speed
///
/// In the interest of speed, there are no checks for out
/// of bounds arithmetic.  This routine is used by most of
/// the functions in the library.
#[must_use]
pub fn evaluate_polynomial(x: f64, coef: &[f64]) -> f64 {
    let mut ans: f64 = 0.0;
    for coef in coef {
        ans = ans * x + coef;
    }
    ans
}

/// Evaluate polynomial when coefficient of x is 1.0.
/// Otherwise same as [`evaluate_polynomial`].
#[must_use]
pub fn evaluate_polynomial_1(x: f64, coef: &[f64]) -> f64 {
    let mut ans = 1.0;
    for coef in coef {
        ans = ans * x + coef;
    }
    ans
}

#[test]
fn works() {
    assert_float_eq::assert_float_relative_eq!(evaluate_polynomial(2.0, &[1.0, -3.0, -2.5]), -4.5);
    assert_float_eq::assert_float_relative_eq!(evaluate_polynomial(2.0, &[]), 0.0);
    assert_float_eq::assert_float_relative_eq!(evaluate_polynomial_1(2.0, &[1.0, -3.0, -2.5]), 3.5);
    assert_float_eq::assert_float_relative_eq!(evaluate_polynomial_1(2.0, &[]), 1.0);
}

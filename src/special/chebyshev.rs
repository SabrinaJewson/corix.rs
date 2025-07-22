/// Evaluate Chebyshev series
///
/// Evaluates the series
///
/// ```text
///     N-1
/// y =  ∑ coef[i] Tᵢ (x/2)
///     i=0
/// ```
///
/// of Chebyshev polynomials `Tᵢ` at argument `x/2`.
///
/// Coefficients are stored in reverse order, i.e. the zero
/// order term is last in the array.
///
/// If coefficients are for the interval a to b, x must
/// have been transformed to x -> 2(2x - b - a)/(b-a) before
/// entering the routine.  This maps x from (a, b) to (-1, 1),
/// over which the Chebyshev polynomials are defined.
///
/// If the coefficients are for the inverted interval, in
/// which (a, b) is mapped to (1/b, 1/a), the transformation
/// required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
/// this becomes x -> 4a/x - 1.
///
/// # Speed
///
/// Taking advantage of the recurrence properties of the
/// Chebyshev polynomials, the routine requires one more
/// addition per loop than evaluating a nested polynomial of
/// the same degree.
#[must_use]
pub fn eval_chebyshev(x: f64, coef: &[f64]) -> f64 {
    let mut b0 = 0.0;
    let mut b1 = 0.0;
    let mut b2 = 0.0;
    for v in coef {
        b2 = b1;
        b1 = b0;
        b0 = x.mul_add(b1, -b2) + v;
    }
    0.5f64 * (b0 - b2)
}

#[test]
fn works() {
    assert_float_eq::assert_float_relative_eq!(eval_chebyshev(f64::INFINITY, &[]), 0.0);
    assert_float_eq::assert_float_relative_eq!(
        eval_chebyshev(2.0, &[-2.7, 8.3, 37.0, 12.0, -6.0]),
        51.6
    );
}

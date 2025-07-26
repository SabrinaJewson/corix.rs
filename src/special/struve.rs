//! [Struve functions](https://en.wikipedia.org/wiki/Struve_function).

/// Struve function `Hᵥ(x)`.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::struve(1.0, 5.0), 0.8078119457940645);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn struve(v: f64, x: f64) -> f64 {
    unsafe { crate::xsf::struve_h(v, x) }
}

/// Modified Struve function `Lᵥ(x)`.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::struve_mod(1.0, 5.0), 23.728215780408284);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn struve_mod(v: f64, x: f64) -> f64 {
    unsafe { crate::xsf::struve_l(v, x) }
}

/// Integrate the Struve function with `v = 0`.
///
/// Calculates `∫₀ˣ H₀(t) dt`.
///
/// For the function itself, see [`struve`].
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::integrate_struve0(5.0), 2.044243662660246);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn integrate_struve0(x: f64) -> f64 {
    unsafe { crate::xsf::itstruve0(x) }
}

/// Integrate an alternate form of the Struve function with `v = 0`.
///
/// Calculates `∫ₓ^∞ H₀(t) / t dt`.
///
/// For the function itself, see [`struve`].
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::integrate_struve0_alt(5.0), 0.079545750553913);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn integrate_struve0_alt(x: f64) -> f64 {
    unsafe { crate::xsf::it2struve0(x) }
}


/// Integrate the modified Struve function with `v = 0`.
///
/// Calculates `∫₀¹ L₀(t) / t dt`.
///
/// For the function itself, see [`struve_mod`].
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::integrate_struve_mod0(5.0), 30.03078811315891);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn integrate_struve_mod0(x: f64) -> f64 {
    unsafe { crate::xsf::itmodstruve0(x) }
}

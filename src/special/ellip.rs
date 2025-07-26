//! Elliptic functions and integrals.

/// The result of the Jacobi elliptic functions `sn`, `cn`, `dn` and `am`.
/// Returned by [`jacobi`].
#[derive(Debug, Default, Clone, Copy)]
pub struct Jacobi {
    /// The result of the `sn` function.
    pub sn: f64,
    /// The result of the `cn` function.
    pub cn: f64,
    /// The result of the `dn` function.
    pub dn: f64,
    /// The Jacobi amplitude, otherwise known as φ.
    pub am: f64,
}

/// Calculate `sn`, `cn`, `dn` and the Jacobi amplitude φ.
///
/// `m` should be between 0 and 1.
///
/// # Examples
///
/// ```
/// let j = corix::ellip::jacobi(2.0, 0.5);
/// assert_float_relative_eq!(j.sn, 0.9946623253580177);
/// assert_float_relative_eq!(j.cn, -0.10318361552776205);
/// assert_float_relative_eq!(j.dn, 0.710861047784087);
/// assert_float_relative_eq!(j.am, 1.6741639220482394);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn jacobi(u: f64, m: f64) -> Jacobi {
    let mut jacobi = Jacobi::default();
    unsafe {
        crate::xsf::ellipj(
            u,
            m,
            &raw mut jacobi.sn,
            &raw mut jacobi.cn,
            &raw mut jacobi.dn,
            &raw mut jacobi.am,
        );
    };
    jacobi
}

/// Complete elliptic integral of the first kind.
///
/// See also [`k_1m`] for greater accuracy when passing in `1 - m`.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::ellip::k(0.7), 2.075363135292469);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn k(m: f64) -> f64 {
    unsafe { crate::xsf::ellipk(m) }
}

/// Complete elliptic integral of the first kind, using `1 - m`
/// (i.e. like [`k`]`(1 - m)` but with greater accuracy).
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::ellip::k_1m(0.3), 2.075363135292469);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn k_1m(m: f64) -> f64 {
    unsafe { crate::xsf::ellipkm1(m) }
}

/// Incomplete elliptic integral of the first kind.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::ellip::k_incomplete(PI / 2.0, 0.7), 2.075363135292469);
/// # use assert_float_eq::assert_float_relative_eq;
/// # use std::f64::consts::PI;
/// ```
#[must_use]
pub fn k_incomplete(φ: f64, m: f64) -> f64 {
    unsafe { crate::xsf::ellipkinc(φ, m) }
}

/// Complete elliptic integral of the second kind.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::ellip::e(0.7), 1.2416705679458229);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn e(m: f64) -> f64 {
    unsafe { crate::xsf::ellipe(m) }
}

/// Incomplete elliptic integral of the second kind.
///
/// # Examples
///
/// ```
/// assert_float_relative_eq!(corix::ellip::e_incomplete(PI / 2.0, 0.7), 1.2416705679458229);
/// # use assert_float_eq::assert_float_relative_eq;
/// # use std::f64::consts::PI;
/// ```
#[must_use]
pub fn e_incomplete(φ: f64, m: f64) -> f64 {
    unsafe { crate::xsf::ellipeinc(φ, m) }
}

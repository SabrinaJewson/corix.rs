//! [Bessel functions](https://en.wikipedia.org/wiki/Bessel_function).
// NB. `kn` weirdness: https://github.com/scipy/scipy/issues/23386
// TODO: Zeros of Bessel functions. The problem is that the choice of functions offered by specfun
// is seemingly extremely arbitrary :/

/// Bessel function of the first kind, of order `v`.
///
/// See also [`j0`] and [`j1`] for special cases when `v` is 0 or 1 respectively.
/// See [`j_exp`] for an exponentially scaled version.
#[must_use]
pub fn j<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_j, xsf::cyl_bessel_j_complex);
    unsafe { f(v, x) }
}

/// Exponentially scaled Bessel function of the first kind, of order `v`.
///
/// See [`j`] for the basic version. Defined as:
///
/// ```text
/// j_exp(v, z) = j(v, z) × exp(-|Im(z)|)
/// ```
#[must_use]
pub fn j_exp<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_je, xsf::cyl_bessel_je_complex);
    unsafe { f(v, x) }
}

/// Special case of [`j`] where `v = 0`.
///
/// See also [`j1`] for when `v = 1`.
#[must_use]
pub fn j0(x: f64) -> f64 {
    unsafe { xsf::cyl_bessel_j0(x) }
}

/// Special case of [`j`] where `v = 1`.
///
/// See also [`j0`] for when `v = 0`.
#[must_use]
pub fn j1(x: f64) -> f64 {
    unsafe { xsf::cyl_bessel_j1(x) }
}

/// Bessel function of the second kind, of order `v`.
///
/// See also [`yi`] for the special case where `v` is an integer.
/// See [`y_exp`] for an exponentially scaled version.
#[must_use]
pub fn y<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_y, xsf::cyl_bessel_y_complex);
    unsafe { f(v, x) }
}

/// Exponentially scaled Bessel function of the second kind, of order `v`.
///
/// See [`y`] for the basic version. Defined as:
///
/// ```text
/// y_exp(v, z) = y(v, z) × exp(-|Im(z)|)
/// ```
#[must_use]
pub fn y_exp<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_ye, xsf::cyl_bessel_ye_complex);
    unsafe { f(v, x) }
}

/// Special case of [`y`] where `v` is an integer.
#[must_use]
pub fn yi(v: i32, x: f64) -> f64 {
    // This check is actually duplicated in `yn`, but if we put it here then the optimizer can see
    // through it.
    match v {
        -1 => -unsafe { xsf::cyl_bessel_y1(x) },
        0 => unsafe { xsf::cyl_bessel_y0(x) },
        1 => unsafe { xsf::cyl_bessel_y1(x) },
        _ => unsafe { xsf::cyl_bessel_yn(v, x) },
    }
}

/// Modified Bessel function of the first kind, of order `v`.
///
/// See also [`i0`] and [`i1`] for special cases when `v` is 0 or 1 respectively.
/// See [`i_exp`] for an exponentially scaled version.
#[must_use]
pub fn i<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_i, xsf::cyl_bessel_i_complex);
    unsafe { f(v, x) }
}

/// Exponentially scaled modified Bessel function of the first kind, of order `v`.
///
/// See [`i`] for the basic version. Defined as:
///
/// ```text
/// i_exp(v, z) = i(v, z) × exp(-|Re(z)|)
/// ```
/// See also [`i0_exp`] and [`i1_exp`] for special cases when `v` is 0 or 1 respectively.
#[must_use]
pub fn i_exp<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_ie, xsf::cyl_bessel_ie_complex);
    unsafe { f(v, x) }
}

/// Special case of [`i`] where `v = 0`.
///
/// See also [`i1`] for when `v = 1`.
#[must_use]
pub fn i0(x: f64) -> f64 {
    unsafe { xsf::cyl_bessel_i0(x) }
}

/// Special case of [`i_exp`] where `v = 0`.
///
/// See also [`i1_exp`] for when `v = 1`.
#[must_use]
pub fn i0_exp(x: f64) -> f64 {
    unsafe { xsf::cyl_bessel_i0e(x) }
}

/// Special case of [`i`] where `v = 1`.
///
/// See also [`i0`] for when `v = 0`.
#[must_use]
pub fn i1(x: f64) -> f64 {
    unsafe { xsf::cyl_bessel_i1(x) }
}

/// Special case of [`i_exp`] where `v = 1`.
///
/// See also [`i0_exp`] for when `v = 0`.
#[must_use]
pub fn i1_exp(x: f64) -> f64 {
    unsafe { xsf::cyl_bessel_i1e(x) }
}

/// Modified Bessel function of the second kind, of order `v`.
///
/// See also [`ki`] for the special case where `v` is an integer.
/// See [`k_exp`] for an exponentially scaled version.
#[must_use]
pub fn k<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_k, xsf::cyl_bessel_k_complex);
    unsafe { f(v, x) }
}

/// Exponentially scaled modified Bessel function of the second kind, of order `v`.
///
/// See [`k`] for the basic version. Defined as:
///
/// ```text
/// k_exp(v, z) = k(v, z) × exp(z)
/// ```
///
/// See also [`ki_exp`] for the special case where `v` is an integer.
#[must_use]
pub fn k_exp<K: ComplexFloat<Real = f64>>(v: f64, x: K) -> K {
    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(f64, K) -> K;
    }
    let f = util::elim_complex_float::<K, C>(xsf::cyl_bessel_ke, xsf::cyl_bessel_ke_complex);
    unsafe { f(v, x) }
}

/// Special case of [`k`] where `v` is an integer.
#[must_use]
pub fn ki(v: i32, x: f64) -> f64 {
    // This check does not exist in CEPHES `kn`.
    // TODO: Measure the accuracy/speed of `k0(x)` vs `kn(0, x)`.
    match v {
        0 => unsafe { xsf::cyl_bessel_k0(x) },
        -1 | 1 => unsafe { xsf::cyl_bessel_k1(x) },
        _ => unsafe { xsf::cyl_bessel_kn(v, x) },
    }
}

/// Special case of [`k_exp`] where `v` is an integer.
#[must_use]
pub fn ki_exp(v: i32, x: f64) -> f64 {
    match v {
        0 => unsafe { xsf::cyl_bessel_k0e(x) },
        -1 | 1 => unsafe { xsf::cyl_bessel_k1e(x) },
        // There is no `kne` in CEPHES, so we fall back to the float version here, judging this
        // option to be better than:
        // - not exposing `k0e` and `k1e`, or
        // - having `k0_exp` but not `k0`, or
        // - adding `k0`, which would be redundant given `ki` exists.
        _ => k_exp(f64::from(v), x),
    }
}

/// Integrate Bessel functions of the first and second kinds of order zero.
///
/// Calculates `∫₀ˣ J₀(t) dt` and `∫₀ˣ Y₀(t) dt`.
///
/// See [`integrate_j0_y0_alt`] for an alternate form.
/// See [`j0`] and [`yi`] for evaluating these functions.
#[must_use]
pub fn integrate_j0_y0(x: f64) -> (f64, f64) {
    let [mut j0, mut y0] = [0.0; 2];
    unsafe { xsf::it1j0y0(x, &raw mut j0, &raw mut y0) };
    (j0, y0)
}

/// Integrate Bessel functions of the first and second kinds of order zero, in an alternate form.
///
/// Calculates `∫₀ˣ (1 − J₀(t)) / t dt` and `∫ₓ^∞ Y₀(t)/t dt`.
///
/// See [`integrate_j0_y0`] for the basic integral.
/// See [`j0`] and [`yi`] for evaluating these functions.
#[must_use]
pub fn integrate_j0_y0_alt(x: f64) -> (f64, f64) {
    let [mut j0, mut y0] = [0.0; 2];
    unsafe { xsf::it2j0y0(x, &raw mut j0, &raw mut y0) };
    (j0, y0)
}

/// Integrate modified Bessel functions of the first and second kinds of order zero.
///
/// Calculates `∫₀ˣ I₀(t) dt` and `∫₀ˣ K₀(t) dt`.
///
/// See [`integrate_i0_k0_alt`] for an alternate form.
/// See [`i0`] and [`ki`] for evaluating these functions.
#[must_use]
pub fn integrate_i0_k0(x: f64) -> (f64, f64) {
    let [mut i0, mut k0] = [0.0; 2];
    unsafe { xsf::it1i0k0(x, &raw mut i0, &raw mut k0) };
    (i0, k0)
}

/// Integrate modified Bessel functions of the first and second kinds of order zero, in an alternate
/// form.
///
/// Calculates `∫₀ˣ (I₀(t) − 1) / t dt` and `∫ₓ^∞ K₀(t)/t dt`.
///
/// See [`integrate_i0_k0`] for the basic integral.
/// See [`i0`] and [`ki`] for evaluating these functions.
#[must_use]
pub fn integrate_i0_k0_alt(x: f64) -> (f64, f64) {
    let [mut i0, mut k0] = [0.0; 2];
    unsafe { xsf::it2i0k0(x, &raw mut i0, &raw mut k0) };
    (i0, k0)
}

/// Weighted integral of the Bessel function of the first kind.
///
/// Calculates `∫₀¹ x ^ λ × Jᵥ(2ax) dx`.
///
/// For `Jᵥ` itself, see [`j`].
#[must_use]
pub fn poly(a: f64, λ: f64, v: f64) -> f64 {
    unsafe { xsf::besselpoly(a, λ, v) }
}

/// Hankel function of the first kind.
///
/// See [`hankel1_exp`] for an exponentially scaled version.
#[must_use]
pub fn hankel1(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_1(v, z) }
}

/// Exponentially scaled Hankel function of the first kind.
///
/// See [`hankel1`] for the basic version.
#[must_use]
pub fn hankel1_exp(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_1e(v, z) }
}

/// Hankel function of the second kind.
///
/// See [`hankel2_exp`] for an exponentially scaled version.
#[must_use]
pub fn hankel2(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_2(v, z) }
}

/// Exponentially scaled Hankel function of the second kind.
///
/// See [`hankel2`] for the basic version.
#[must_use]
pub fn hankel2_exp(v: f64, z: Complex<f64>) -> Complex<f64> {
    unsafe { xsf::cyl_hankel_2e(v, z) }
}

/// Wright’s generalized Bessel function.
///
/// Only non-negative arguments are implemented.
///
/// See [`ln_wright`] for its natural logarithm.
#[must_use]
pub fn wright(a: f64, b: f64, x: f64) -> f64 {
    unsafe { xsf::wright_bessel(a, b, x) }
}

/// The natural logarithm of [`wright`].
///
/// Only non-negative arguments are implemented.
#[must_use]
pub fn ln_wright(a: f64, b: f64, x: f64) -> f64 {
    unsafe { xsf::log_wright_bessel(a, b, x) }
}

/// Compute the nth derivative of [`j`], the Bessel function of the first kind.
///
/// Computes the derivative with respect to `x`, with fixed order `v`.
#[must_use]
pub fn j_deriv<K: ComplexFloat<Real = f64>>(v: f64, x: K, n: u32) -> K {
    bessel_diff_formula(v, x, n, j, -1.0)
}

/// Compute the nth derivative of [`y`], the Bessel function of the second kind.
///
/// Computes the derivative with respect to `x`, with fixed order `v`.
#[must_use]
pub fn y_deriv<K: ComplexFloat<Real = f64>>(v: f64, x: K, n: u32) -> K {
    bessel_diff_formula(v, x, n, y, -1.0)
}

/// Compute the nth derivative of [`k`], the modified Bessel function of the second kind.
///
/// Computes the derivative with respect to `x`, with fixed order `v`.
#[must_use]
pub fn k_deriv<K: ComplexFloat<Real = f64>>(v: f64, x: K, n: u32) -> K {
    let d = bessel_diff_formula(v, x, n, k, 1.0);
    if n % 2 == 0 { d } else { -d }
}

/// Compute the nth derivative of [`i`], the modified Bessel function of the first kind.
///
/// Computes the derivative with respect to `x`, with fixed order `v`.
#[must_use]
pub fn i_deriv<K: ComplexFloat<Real = f64>>(v: f64, x: K, n: u32) -> K {
    bessel_diff_formula(v, x, n, k, 1.0)
}

/// Compute the nth derivative of [`hankel1`], the Hankel function of the first kind.
///
/// Computes the derivative with respect to `x`, with fixed order `v`.
#[must_use]
pub fn hankel1_deriv(v: f64, x: Complex<f64>, n: u32) -> Complex<f64> {
    bessel_diff_formula(v, x, n, hankel1, -1.0)
}

/// Compute the nth derivative of [`hankel2`], the Hankel function of the second kind.
///
/// Computes the derivative with respect to `x`, with fixed order `v`.
#[must_use]
pub fn hankel2_deriv(v: f64, x: Complex<f64>, n: u32) -> Complex<f64> {
    bessel_diff_formula(v, x, n, hankel2, -1.0)
}

#[expect(clippy::many_single_char_names)]
fn bessel_diff_formula<K: ComplexFloat<Real = f64>>(
    v: f64,
    x: K,
    n: u32,
    l: impl Fn(f64, K) -> K,
    phase: f64,
) -> K {
    // Ported from:
    // https://github.com/scipy/scipy/blob/v1.16.0/scipy/special/_basic.py#L814-L825
    let mut p = 1.0;
    let mut s = l(v - f64::from(n), x);
    for i in 1..=n {
        p = phase * (p * f64::from(n - i + 1)) / f64::from(i);
        s = s + K::from(p).unwrap() * l(v - f64::from(n) + f64::from(i) * 2.0, x);
    }
    s / K::from(f64::from(n).exp2()).unwrap()
}

/// Jahnke-Emden Lambda function.
///
/// This function will compute `Λ_v(x)` (and its first derivative) at all integer increments of `v`
/// starting from zero up to the `v` passed into the function.
///
/// The number of computed values – which may be lower than `v as usize + 1` – is returned.
///
/// # Examples
///
/// ```
/// let [mut Λ, mut Λ_deriv] = [[f64::NAN; 3]; 2];
/// let n = corix::bessel::lambda(2.5, 15.0, &mut Λ, &mut Λ_deriv);
/// let [Λ, Λ_deriv] = [&Λ[..n], &Λ_deriv[..n]];
///
/// // The values at v = 0.5.
/// assert_float_relative_eq!(Λ[0], 0.04335252267714108);
/// assert_float_relative_eq!(Λ_deriv[0], -0.053536029035730764);
///
/// // The values at v = 1.5.
/// assert_float_relative_eq!(Λ[1], 0.010707205807146153);
/// assert_float_relative_eq!(Λ_deriv[1], 0.0065290633739989844);
///
/// // The values at v = 2.5.
/// assert_float_relative_eq!(Λ[2], -0.0021763544579996613);
/// assert_float_relative_eq!(Λ_deriv[2], 0.004294520088381938);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
///
/// # Panics
///
/// Panics if:
/// - `Λ` and `Λ_deriv` are of different lengths, or
/// - `v` is zero or positive and `Λ`’s length is not at least `v as usize + 1`
///   (note that `v as usize` rounds toward zero and saturates infinities).
#[expect(non_snake_case)]
#[must_use]
pub fn lambda(v: f64, x: f64, Λ: &mut [f64], Λ_deriv: &mut [f64]) -> usize {
    // See also:
    // https://github.com/scipy/scipy/blob/v1.16.0/scipy/special/_basic.py#L2208-L2255

    assert_eq!(Λ.len(), Λ_deriv.len());

    // Filter out NaNs and negatives.
    if !(v >= 0.0) {
        return 0;
    }

    // Check the function’s preconditions.
    #[expect(clippy::cast_possible_truncation)]
    #[expect(clippy::cast_sign_loss)]
    let () = assert!((v as usize) < Λ.len());

    #[expect(clippy::cast_possible_truncation)]
    let n = v as i32;
    debug_assert!(0 <= n);

    let v_frac = v - f64::from(n);
    debug_assert!(v_frac < 1.0);

    let call = |Λ: &mut [f64], Λ_deriv: &mut [f64]| {
        if v_frac == 0.0 {
            let mut nm = 0;
            unsafe { xsf::lamn(n, x, &raw mut nm, Λ.as_mut_ptr(), Λ_deriv.as_mut_ptr()) };
            assert!(0 <= nm && nm <= n);
            usize::try_from(nm).unwrap() + 1
        } else {
            let mut vm = 0.0;
            unsafe { xsf::lamv(v, x, &raw mut vm, Λ.as_mut_ptr(), Λ_deriv.as_mut_ptr()) };
            assert!(0.0 <= vm && vm <= v);
            #[expect(clippy::cast_possible_truncation)]
            #[expect(clippy::cast_sign_loss)]
            let res = vm as usize + 1;
            res
        }
    };

    if Λ.len() == 1 {
        // The `lamn` and `lamv` functions _always_ write to the first two array elements. So we
        // handle this special case with a stack buffer.
        let mut Λ_buf = [0.0; 2];
        let mut Λ_deriv_buf = [0.0; 2];
        let ret = call(&mut Λ_buf, &mut Λ_deriv_buf);
        Λ[0] = Λ_buf[0];
        Λ_deriv[0] = Λ_deriv_buf[0];
        ret
    } else {
        call(Λ, Λ_deriv)
    }
}

pub mod spherical {
    //! Spherical Bessel functions.

    enum C {}
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(u64, K) -> K;
    }

    /// Spherical Bessel function of the first kind.
    #[must_use]
    pub fn j<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f = util::elim_complex_float::<K, C>(xsf::sph_bessel_j, xsf::sph_bessel_j_complex);
        unsafe { f(n, x) }
    }
    /// Spherical Bessel function of the second kind.
    #[must_use]
    pub fn y<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f = util::elim_complex_float::<K, C>(xsf::sph_bessel_y, xsf::sph_bessel_y_complex);
        unsafe { f(n, x) }
    }
    /// Modified spherical Bessel function of the first kind.
    #[must_use]
    pub fn i<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f = util::elim_complex_float::<K, C>(xsf::sph_bessel_i, xsf::sph_bessel_i_complex);
        unsafe { f(n, x) }
    }
    /// Modified spherical Bessel function of the second kind.
    #[must_use]
    pub fn k<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f = util::elim_complex_float::<K, C>(xsf::sph_bessel_k, xsf::sph_bessel_k_complex);
        unsafe { f(n, x) }
    }
    /// First derivative of the spherical Bessel function of the first kind.
    #[must_use]
    pub fn j_deriv<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f =
            util::elim_complex_float::<K, C>(xsf::sph_bessel_j_jac, xsf::sph_bessel_j_jac_complex);
        unsafe { f(n, x) }
    }
    /// First derivative of the spherical Bessel function of the second kind.
    #[must_use]
    pub fn y_deriv<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f =
            util::elim_complex_float::<K, C>(xsf::sph_bessel_y_jac, xsf::sph_bessel_y_jac_complex);
        unsafe { f(n, x) }
    }
    /// First derivative of the modified spherical Bessel function of the first kind.
    #[must_use]
    pub fn i_deriv<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f =
            util::elim_complex_float::<K, C>(xsf::sph_bessel_i_jac, xsf::sph_bessel_i_jac_complex);
        unsafe { f(n, x) }
    }
    /// First derivative of the modified spherical Bessel function of the second kind.
    #[must_use]
    pub fn k_deriv<K: ComplexFloat<Real = f64>>(n: u64, x: K) -> K {
        let f =
            util::elim_complex_float::<K, C>(xsf::sph_bessel_k_jac, xsf::sph_bessel_k_jac_complex);
        unsafe { f(n, x) }
    }

    use crate::util;
    use crate::util::ComplexFloatMotive;
    use crate::xsf;
    use num_complex::ComplexFloat;
}

/// Compute Riccati–Bessel functions of the first kind and their first derivatives.
///
/// Computes the functions of order `n = 0, 1, 2, …`.
/// Returns the number of orders computed, which may be less than the input slice length.
///
/// # Examples
///
/// ```
/// let [mut s, mut s_deriv] = [[0.0; 3]; 2];
/// let n = corix::bessel::riccati_s(0.5, &mut s, &mut s_deriv);
/// let [s, s_deriv] = [&s[..n], &s_deriv[..n]];
///
/// // The Riccati–Bessel function and its first derivative with order 0.
/// assert_float_relative_eq!(s[0], 0.479425538604203);
/// assert_float_relative_eq!(s_deriv[0], 0.8775825618903728);
///
/// // The Riccati–Bessel function and its first derivative with order 1.
/// assert_float_relative_eq!(s[1], 0.08126851531803328);
/// assert_float_relative_eq!(s_deriv[1], 0.31688850796813645);
///
/// // The Riccati–Bessel function and its first derivative with order 2.
/// assert_float_relative_eq!(s[2], 0.008185553303996706);
/// assert_float_relative_eq!(s_deriv[2], 0.048526302102046455);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
///
/// # Panics
///
/// Panics if
/// - `s` is not the same length as `s_deriv`, or
/// - `s` is too long.
pub fn riccati_s(x: f64, s: &mut [f64], s_deriv: &mut [f64]) -> usize {
    assert_eq!(s.len(), s_deriv.len());
    if s.is_empty() {
        return 0;
    }

    let len = s.len().try_into().unwrap();

    let call = |s: &mut [f64], s_deriv: &mut [f64]| -> usize {
        let mut nm = 0;
        unsafe { xsf::rctj(x, &raw mut nm, s.as_mut_ptr(), s_deriv.as_mut_ptr(), len) };
        (nm + 1).try_into().unwrap()
    };

    if len == 1 {
        let [mut s_buf, mut s_deriv_buf] = [[0.0; 2]; 2];
        let res = call(&mut s_buf, &mut s_deriv_buf);
        s[0] = s_buf[0];
        s_deriv[0] = s_deriv_buf[0];
        res
    } else {
        call(s, s_deriv)
    }
}

#[test]
fn riccati_s_short_slice() {
    let [mut s, mut s_deriv] = [[0.0; 1]; 2];
    assert_eq!(riccati_s(0.5, &mut s, &mut s_deriv), 1);
    assert_float_eq::assert_float_relative_eq!(s[0], 0.479_425_538_604_203);
    assert_float_eq::assert_float_relative_eq!(s_deriv[0], 0.877_582_561_890_372_8);
}

/// Compute Riccati–Bessel functions of the second kind and their first derivatives.
///
/// Computes the functions of order `n = 0, 1, 2, …`.
/// Returns the number of orders computed, which may be less than the input slice length.
///
/// # Examples
///
/// ```
/// let [mut c, mut c_deriv] = [[0.0; 3]; 2];
/// let n = corix::bessel::riccati_c(0.5, &mut c, &mut c_deriv);
/// let [c, c_deriv] = [&c[..n], &c_deriv[..n]];
/// dbg!(c, c_deriv);
///
/// // The Riccati–Bessel function and its first derivative with order 0.
/// assert_float_relative_eq!(c[0], -0.8775825618903728);
/// assert_float_relative_eq!(c_deriv[0], 0.479425538604203);
///
/// // The Riccati–Bessel function and its first derivative with order 1.
/// assert_float_relative_eq!(c[1], -2.2345906623849485);
/// assert_float_relative_eq!(c_deriv[1], 3.5915987628795243);
///
/// // The Riccati–Bessel function and its first derivative with order 2.
/// assert_float_relative_eq!(c[2], -12.529961412419318);
/// assert_float_relative_eq!(c_deriv[2], 47.88525498729233);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
///
/// # Panics
///
/// Panics if
/// - `c` is not the same length as `c_deriv`, or
/// - `c` is too long.
pub fn riccati_c(x: f64, c: &mut [f64], c_deriv: &mut [f64]) -> usize {
    assert_eq!(c.len(), c_deriv.len());
    if c.is_empty() {
        return 0;
    }

    let len = c.len().try_into().unwrap();

    let call = |c: &mut [f64], c_deriv: &mut [f64]| -> usize {
        let mut nm = 0;
        unsafe { xsf::rcty(x, &raw mut nm, c.as_mut_ptr(), c_deriv.as_mut_ptr(), len) };
        (nm + 1).try_into().unwrap()
    };

    if len == 1 {
        let [mut buf, mut deriv_buf] = [[0.0; 2]; 2];
        let res = call(&mut buf, &mut deriv_buf);
        c[0] = buf[0];
        c_deriv[0] = deriv_buf[0];
        res
    } else {
        call(c, c_deriv)
    }
}

use crate::util;
use crate::util::ComplexFloatMotive;
use crate::xsf;
use num_complex::Complex;
use num_complex::ComplexFloat;

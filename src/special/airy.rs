//! [Airy functions](https://en.wikipedia.org/wiki/Airy_function).

/// The result of [`airy`] and [`airy_exp`].
#[derive(Debug, Default, Clone, Copy)]
pub struct Airy<K = f64> {
    /// The result of `Ai(x)`.
    pub ai: K,
    /// The result of `Bi(x)`.
    pub bi: K,
    /// The derivative of `ai`, `Ai'(x)`.
    pub ai_deriv: K,
    /// The derivative of `bi`, `Bi'(x)`.
    pub bi_deriv: K,
}

/// Compute the airy functions `Ai` and `Bi`, and their first derivatives.
///
/// # Example
///
/// ```
/// let airy = corix::airy(0.5);
/// assert_float_relative_eq!(airy.ai, 0.231_693_606_480_833);
/// assert_float_relative_eq!(airy.bi, 0.854_277_043_103_155);
/// assert_float_relative_eq!(airy.ai_deriv, -0.224_910_532_664_684);
/// assert_float_relative_eq!(airy.bi_deriv, 0.544_572_564_140_592_4);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn airy<K: ComplexFloat<Real = f64>>(x: K) -> Airy<K> {
    struct C;
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(K, *mut K, *mut K, *mut K, *mut K);
    }
    let f = crate::util::elim_complex_float::<K, C>(crate::xsf::airy, crate::xsf::airy_complex);

    let [mut ai, mut ai_deriv, mut bi, mut bi_deriv] = [K::zero(); 4];
    unsafe {
        f(
            x,
            &raw mut ai,
            &raw mut ai_deriv,
            &raw mut bi,
            &raw mut bi_deriv,
        );
    };
    Airy {
        ai,
        bi,
        ai_deriv,
        bi_deriv,
    }
}

/// Exponentially scaled Airy functions and their first derivatives.
///
/// The scaling is as follows:
///
/// ```text
/// ai *= exp(2/3 × z × sqrt(z))
/// bi *= exp(-|2/3 × Re(z × sqrt(z))|)
/// ```
///
/// (with the derivates scaled accordingly.)
///
/// # Example
///
/// ```
/// let airy = corix::airy_exp(0.5);
/// assert_float_relative_eq!(airy.ai, 0.293_277_159_129_947);
/// assert_float_relative_eq!(airy.bi, 0.674_892_411_115_630);
/// assert_float_relative_eq!(airy.ai_deriv, -0.284_691_162_091_942);
/// assert_float_relative_eq!(airy.bi_deriv, 0.430_220_961_463_770);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn airy_exp<K: ComplexFloat<Real = f64>>(x: K) -> Airy<K> {
    struct C;
    impl ComplexFloatMotive for C {
        type Output<K: ComplexFloat> = unsafe extern "C" fn(K, *mut K, *mut K, *mut K, *mut K);
    }
    let f = crate::util::elim_complex_float::<K, C>(crate::xsf::airye, crate::xsf::airye_complex);

    let [mut ai, mut ai_deriv, mut bi, mut bi_deriv] = [K::zero(); 4];
    unsafe {
        f(
            x,
            &raw mut ai,
            &raw mut ai_deriv,
            &raw mut bi,
            &raw mut bi_deriv,
        );
    };
    Airy {
        ai,
        bi,
        ai_deriv,
        bi_deriv,
    }
}

/// Compute zeros of the Airy function `Ai`.
///
/// This function takes four slice out-parameters, required to be of equal length.
/// - `zeros` will be filled with the `x`-values of the zeros of `Ai(x)`.
/// - `deriv_zeros` will be filled with the `x`-values of the zeros of `Ai'(x)`.
/// - `ai` will be filled with `Ai(deriv_zeros)` (the y-values of `Ai`’s stationary points).
/// - `ai_deriv` will be filled with `Ai'(zeros)` (the gradients at `Ai`’s zeros).
///
/// See also [`bi_zeros`].
///
/// # Panics
///
/// Panics if the slices do not all have equal lengths.
///
/// # Examples
///
/// ```
/// let [mut zeros, mut deriv_zeros, mut ai, mut ai_deriv] = [[0.0; 2]; 4];
/// corix::ai_zeros(&mut zeros, &mut deriv_zeros, &mut ai, &mut ai_deriv);
///
/// // The first zero of the Ai function, and the derivative at that point.
/// assert_float_relative_eq!(zeros[0], -2.3381074104597674);
/// assert_float_relative_eq!(ai_deriv[0], 0.7012108227206915);
///
/// // The second zero, and the derivative at that point.
/// assert_float_relative_eq!(zeros[1], -4.08794944413097);
/// assert_float_relative_eq!(ai_deriv[1], -0.8031113696548628);
///
/// // The first stationary point, and the value at that point.
/// assert_float_relative_eq!(deriv_zeros[0], -1.0187929716474715);
/// assert_float_relative_eq!(ai[0], 0.5356566560156998);
///
/// // The second stationary point, and the value at that point.
/// assert_float_relative_eq!(deriv_zeros[1], -3.2481975821798375);
/// assert_float_relative_eq!(ai[1], -0.4190154780325641);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
pub fn ai_zeros(zeros: &mut [f64], deriv_zeros: &mut [f64], ai: &mut [f64], ai_deriv: &mut [f64]) {
    assert_eq!(zeros.len(), deriv_zeros.len());
    assert_eq!(zeros.len(), ai.len());
    assert_eq!(zeros.len(), ai_deriv.len());
    let nt = zeros.len().try_into().unwrap();
    unsafe {
        crate::xsf::airyzo(
            nt,
            1,
            zeros.as_mut_ptr(),
            deriv_zeros.as_mut_ptr(),
            ai.as_mut_ptr(),
            ai_deriv.as_mut_ptr(),
        );
    };
}

/// Compute zeros of the Airy function `Bi`.
///
/// This function takes four slice out-parameters, required to be of equal length.
/// - `zeros` will be filled with the `x`-values of the zeros of `Bi(x)`.
/// - `deriv_zeros` will be filled with the `x`-values of the zeros of `Bi'(x)`.
/// - `bi` will be filled with `Bi(deriv_zeros)` (the y-values of `Bi`’s stationary points).
/// - `bi_deriv` will be filled with `Bi'(zeros)` (the gradients at `Bi`’s zeros).
///
/// See also [`ai_zeros`].
///
/// # Panics
///
/// Panics if the slices do not all have equal lengths.
///
/// # Examples
///
/// ```
/// let [mut zeros, mut deriv_zeros, mut bi, mut bi_deriv] = [[0.0; 2]; 4];
/// corix::bi_zeros(&mut zeros, &mut deriv_zeros, &mut bi, &mut bi_deriv);
/// dbg!(zeros, deriv_zeros, bi, bi_deriv);
///
/// // The first zero of the Ai function, and the derivative at that point.
/// assert_float_relative_eq!(zeros[0], -1.173713222709127);
/// assert_float_relative_eq!(bi_deriv[0], 0.6019578879762394);
///
/// // The second zero, and the derivative at that point.
/// assert_float_relative_eq!(zeros[1], -3.271093302836352);
/// assert_float_relative_eq!(bi_deriv[1], -0.7603101414928012);
///
/// // The first stationary point, and the value at that point.
/// assert_float_relative_eq!(deriv_zeros[0], -2.2944396826141222);
/// assert_float_relative_eq!(bi[0], -0.45494438363965733);
///
/// // The second stationary point, and the value at that point.
/// assert_float_relative_eq!(deriv_zeros[1], -4.073155089071831);
/// assert_float_relative_eq!(bi[1], 0.3965228360944649);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
pub fn bi_zeros(zeros: &mut [f64], deriv_zeros: &mut [f64], bi: &mut [f64], bi_deriv: &mut [f64]) {
    assert_eq!(zeros.len(), deriv_zeros.len());
    assert_eq!(zeros.len(), bi.len());
    assert_eq!(zeros.len(), bi_deriv.len());
    let nt = zeros.len().try_into().unwrap();
    unsafe {
        crate::xsf::airyzo(
            nt,
            2,
            zeros.as_mut_ptr(),
            deriv_zeros.as_mut_ptr(),
            bi.as_mut_ptr(),
            bi_deriv.as_mut_ptr(),
        );
    };
}

/// The result of [`integrate_airy`].
#[derive(Debug, Default, Clone, Copy)]
pub struct AiryIntegral {
    /// The integral of `Ai(x)` from 0 to x.
    pub ai: f64,
    /// The integral of `Bi(x)` from 0 to x.
    pub bi: f64,
    /// The integral of `Ai(-x)` from 0 to x.
    pub ai_neg: f64,
    /// The integral of `Bi(-x)` from 0 to x.
    pub bi_neg: f64,
}

/// Calculate integrals of Airy functions.
///
/// # Examples
///
/// ```
/// let integral = corix::integrate_airy(0.5);
///
/// // ∫_0^0.5 Ai(x) dx
/// assert_float_relative_eq!(integral.ai, 0.14595330516334704);
/// // ∫_0^0.5 Bi(x) dx
/// assert_float_relative_eq!(integral.bi, 0.36533846506717255);
/// // ∫_(-0.5)^0 Ai(x) dx
/// assert_float_relative_eq!(integral.ai_neg, 0.20880954780848898);
/// // ∫_(-0.5)^0 Bi(x) dx
/// assert_float_relative_eq!(integral.bi_neg, 0.2500627549959552);
/// # use assert_float_eq::assert_float_relative_eq;
/// ```
#[must_use]
pub fn integrate_airy(x: f64) -> AiryIntegral {
    let mut res = AiryIntegral::default();
    unsafe {
        crate::xsf::itairy(
            x,
            &raw mut res.ai,
            &raw mut res.bi,
            &raw mut res.ai_neg,
            &raw mut res.bi_neg,
        );
    };
    res
}

use crate::util::ComplexFloatMotive;
use num_complex::ComplexFloat;

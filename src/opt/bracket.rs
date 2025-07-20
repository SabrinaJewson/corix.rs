//! [`Bracket`]: Bracketing the minimum or maximum of a function.

/// Bracketing the minimum or maximum of a function.
///
/// # Example
///
/// Finding a maximum bracket of a complicated curve:
///
/// ```
/// let max = corix::opt::Bracket::new(0.0, 0.1)
///     .maximum(|x| -x.powi(4) - x.powi(3) + 6.0 * x.powi(2) - x - 1.0)
///     .unwrap();
/// assert_float_eq::assert_float_relative_eq!(max.x[0], -2.641_640_786_499_874);
/// assert_float_eq::assert_float_relative_eq!(max.x[1], -2.123_994_917_242_184_6);
/// assert_float_eq::assert_float_relative_eq!(max.x[2], -1.532_623_792_124_926_6);
/// assert_float_eq::assert_float_relative_eq!(max.f[0], 13.249_112_265_752_327);
/// assert_float_eq::assert_float_relative_eq!(max.f[1], 17.421_896_602_708_316);
/// assert_float_eq::assert_float_relative_eq!(max.f[2], 12.708_773_775_896_23);
/// assert_eq!(max.iters, 5);
/// ```
#[derive(Debug)]
#[non_exhaustive]
pub struct Bracket {
    /// Where searching should start.
    pub start: f64,
    /// Where searching should end.
    pub end: f64,
    /// Maximum grow limit.
    /// Defaults to 110.
    pub grow_limit: f64,
    /// The maximum number of iterations before the algorithm fails to converge.
    /// Defaults to 1000.
    pub max_iters: u64,
}

impl Bracket {
    /// Construct a new `Bracket` with default parameters.
    ///
    /// A local minimum need not be contained between `start` and `end`.
    #[must_use]
    pub fn new(start: f64, end: f64) -> Self {
        Self {
            start,
            end,
            // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L2954
            grow_limit: 110.0,
            max_iters: 1000,
        }
    }

    crate::util::setter!(with_grow_limit => grow_limit: f64);
    crate::util::setter!(with_max_iters => max_iters: u64);

    /// Find a maximum bracket of the given function.
    pub fn maximum(&self, mut f: impl FnMut(f64) -> f64) -> Result<Output, Error> {
        self.try_maximum(|x| Ok(f(x)))
    }

    /// Find a minimum bracket of the given function.
    pub fn minimum(&self, mut f: impl FnMut(f64) -> f64) -> Result<Output, Error> {
        self.try_minimum(|x| Ok(f(x)))
    }

    /// Find a maximum bracket of the given fallible function.
    pub fn try_maximum<E>(
        &self,
        mut f: impl FnMut(f64) -> Result<f64, E>,
    ) -> Result<Output, Error<E>> {
        let mut output = self.try_minimum(|x| Ok(-f(x)?))?;
        output.f = output.f.map(|f| -f);
        Ok(output)
    }

    /// Find a minimum bracket of the given fallible function.
    #[expect(clippy::many_single_char_names)]
    pub fn try_minimum<E>(
        &self,
        mut f: impl FnMut(f64) -> Result<f64, E>,
    ) -> Result<Output, Error<E>> {
        // https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_optimize.py#L2954-L3110

        let mut f = |arg: f64| -> Result<f64, Error<E>> {
            let res = f(arg).map_err(Error::Function)?;
            if res.is_nan() {
                return Err(Error::Nan(NanError));
            }
            Ok(res)
        };

        let very_small = 1e-21;

        let [mut a, mut b] = [self.start, self.end];
        let mut fa = f(a)?;
        let mut fb = f(b)?;
        if fa < fb {
            [a, b] = [b, a];
            [fa, fb] = [fb, fa];
        }
        let mut c = b + GOLDEN_RATIO * (b - a);
        let mut fc = f(c)?;
        let mut iters = 0;
        while fc < fb {
            let tmp1 = (b - a) * (fb - fc);
            let tmp2 = (b - c) * (fb - fa);
            let val = tmp2 - tmp1;
            let denom = 2.0
                * if val.abs() < very_small {
                    very_small
                } else {
                    val
                };
            let mut w = b - ((b - c) * tmp2 - (b - a) * tmp1) / denom;
            let wlim = b + self.grow_limit * (c - b);
            if iters > self.max_iters {
                return Err(Error::Convergence(ConvergenceError {
                    calls: false,
                    val: iters,
                }));
            }
            iters += 1;
            let mut fw;
            if (w - c) * (b - w) > 0.0 {
                fw = f(w)?;
                if fw < fc {
                    [a, fa] = [b, fb];
                    [b, fb] = [w, fw];
                    break;
                } else if fw > fb {
                    [c, fc] = [w, fw];
                    break;
                }
                w = c + GOLDEN_RATIO * (c - b);
                fw = f(w)?;
            } else if (w - wlim) * (wlim - c) >= 0.0 {
                [w, fw] = [wlim, f(wlim)?];
            } else if (w - wlim) * (c - w) > 0.0 {
                fw = f(w)?;
                if fw < fc {
                    [b, fb] = [c, fc];
                    [c, fc] = [w, fw];
                    w = c + GOLDEN_RATIO * (c - b);
                    fw = f(w)?;
                }
            } else {
                w = c + GOLDEN_RATIO * (c - b);
                fw = f(w)?;
            }
            [a, fa] = [b, fb];
            [b, fb] = [c, fc];
            [c, fc] = [w, fw];
        }

        // three conditions for a valid bracket
        let cond1 = (fb < fc && fb <= fa) || (fb < fa && fb <= fc);
        let cond2 = (a < b && b < c) || (c < b && b < a);
        let cond3 = a.is_finite() && b.is_finite() && c.is_finite();

        if !(cond1 && cond2 && cond3) {
            return Err(Error::Terminated(Terminated));
        }

        if c < a {
            [a, c] = [c, a];
            [fa, fc] = [fc, fa];
        }

        Ok(Output {
            x: [a, b, c],
            f: [fa, fb, fc],
            iters,
        })
    }
}

/// The output of [`Bracket`].
#[derive(Debug)]
#[non_exhaustive]
pub struct Output {
    /// The x-values `[a, b, c]` that form the bracket.
    /// It holds that `a < b < c`. All three are finite.
    pub x: [f64; 3],
    /// The values `[f(a), f(b), f(c)]`.
    /// It holds that `f(b) ≤ f(a)` and `f(b) ≤ f(c)`,
    /// where one inequality holds strictly;
    /// if finding the maximum, it holds that
    /// `f(a) ≤ f(b)` and `f(c) ≤ f(b)`,
    /// again with one strict.
    pub f: [f64; 3],
    /// The number of iterations.
    pub iters: u64,
}

/// An error from [`Bracket`].
#[derive(Debug)]
pub enum Error<E = Infallible> {
    /// Failed to converge after a number of iterations.
    Convergence(ConvergenceError),

    /// The algorithm could not find a valid bracket.
    /// Try different initial points.
    Terminated(Terminated),

    /// The function returned `NaN`.
    Nan(NanError),

    /// An error from the function being backeted.
    Function(E),
}

impl<E> Display for Error<E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("finding extrema bracket failed")
    }
}

impl<E: 'static + std::error::Error> std::error::Error for Error<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Convergence(e) => Some(e),
            Self::Terminated(e) => Some(e),
            Self::Nan(e) => Some(e),
            Self::Function(e) => Some(e),
        }
    }
}

/// The algorithm could not find a valid bracket.
/// Try different initial points.
#[derive(Debug)]
#[non_exhaustive]
pub struct Terminated;

impl Display for Terminated {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("could not find a valid bracket; perhaps try different initial points")
    }
}

impl std::error::Error for Terminated {}

const GOLDEN_RATIO: f64 = 1.618_034;

use super::ConvergenceError;
use super::NanError;
use core::f64;
use std::convert::Infallible;
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;

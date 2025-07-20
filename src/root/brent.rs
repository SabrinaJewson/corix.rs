//! [`Brent`]: Root-finding using Brent’s method.

/// Root-finding using Brent’s method.
///
/// # Examples
///
/// Solving a non-solvable quintic numerically:
///
/// ```
/// // We know the root is in the range [-4.0, 4.0].
/// let root = corix::root::Brent::new(-4.0, 4.0).root(|x| x.powi(5) - x - 1.0).unwrap();
/// assert_float_eq::assert_float_relative_eq!(root, 1.167_303_978_261_801);
/// ```
#[derive(Debug)]
#[non_exhaustive]
pub struct Brent {
    /// The starting point of the interval which contains the root.
    /// `f(start)` must be of a different sign to `f(end)`.
    pub start: f64,

    /// The ending point of the interval which contains the root.
    /// `f(end)` must be of a different sign to `f(start)`.
    pub end: f64,

    /// Absolute tolerance for what is considered a root.
    /// Defaults to [`DEFAULT_ABS_TOL`].
    pub abs_tol: f64,

    /// Relative tolerance for what is considered a root.
    /// Defaults to [`DEFAULT_REL_TOL`].
    pub rel_tol: f64,

    /// The maximum number of iterations before the algorithm fails to converge.
    /// Defaults to [`DEFAULT_MAX_ITERS`].
    pub max_iters: u64,
}

impl Brent {
    /// Construct a [`Brent`] with default parameters.
    #[must_use]
    pub fn new(start: f64, end: f64) -> Self {
        Self {
            start,
            end,
            abs_tol: DEFAULT_ABS_TOL,
            rel_tol: DEFAULT_REL_TOL,
            max_iters: DEFAULT_MAX_ITERS,
        }
    }

    crate::util::setter!(with_abs_tol => abs_tol: f64);
    crate::util::setter!(with_rel_tol => rel_tol: f64);
    crate::util::setter!(with_max_iters => max_iters: u64);

    /// Run the solver.
    ///
    /// # Errors
    ///
    /// Fails if a tolerance parameter is invalid,
    /// or the signs of the starting points are identical,
    /// or convergence does not succeed.
    pub fn root(&self, mut f: impl FnMut(f64) -> f64) -> Result<f64, Error> {
        self.try_root(|x| Ok(f(x)))
    }

    /// Run the solver with a fallible function.
    ///
    /// # Errors
    ///
    /// Fails if a tolerance parameter is invalid,
    /// or the signs of the starting points are identical,
    /// or convergence does not succeed,
    /// or the function itself returns an error.
    pub fn try_root<E>(&self, mut f: impl FnMut(f64) -> Result<f64, E>) -> Result<f64, Error<E>> {
        // Ported from:
        // - https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L712-L847
        // - https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/Zeros/brentq.c

        InvalidTol::check_abs(false, self.abs_tol).map_err(Error::InvalidTol)?;
        InvalidTol::check_rel(false, self.rel_tol).map_err(Error::InvalidTol)?;

        let mut xpre = self.start;
        let mut xcur = self.end;

        let mut fpre = f(xpre).map_err(Error::Function)?;
        let mut fcur = f(xcur).map_err(Error::Function)?;
        if fpre == 0.0 {
            return Ok(xpre);
        }
        if fcur == 0.0 {
            return Ok(xcur);
        }
        if !(fpre * fcur < 0.0) {
            return Err(Error::SameSigns(SameSigns {
                start: xpre,
                f_start: fpre,
                end: xcur,
                f_end: fcur,
            }));
        }

        let mut xblk = 0.0;
        let mut fblk = 0.0;
        let mut spre = 0.0;
        let mut scur = 0.0;

        for _ in 0..self.max_iters {
            if fpre != 0.0 && fcur != 0.0 && fpre.signum() != fcur.signum() {
                (xblk, fblk) = (xpre, fpre);
                scur = xcur - xpre;
                spre = scur;
            }
            if fblk.abs() < fcur.abs() {
                (xpre, xcur, xblk) = (xcur, xblk, xcur);
                (fpre, fcur, fblk) = (fcur, fblk, fcur);
            }
            let delta = f64::midpoint(self.abs_tol, self.rel_tol * xcur.abs());
            let sbis = (xblk - xcur) / 2.0;
            if fcur == 0.0 || sbis.abs() < delta {
                return Ok(xcur);
            }

            if spre.abs() > delta && fcur.abs() < fpre.abs() {
                #[expect(clippy::float_cmp)]
                let stry = if xpre == xblk {
                    // interpolate
                    -fcur * (xcur - xpre) / (fcur - fpre)
                } else {
                    // extrapolate
                    let dpre = (fpre - fcur) / (xpre - xcur);
                    let dblk = (fblk - fcur) / (xblk - xcur);
                    -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
                };
                if 2.0 * stry.abs() < f64::min(spre.abs(), 3.0 * sbis.abs() - delta) {
                    // good short step
                    (spre, scur) = (scur, stry);
                } else {
                    // bisect
                    (spre, scur) = (sbis, sbis);
                }
            } else {
                // bisect
                (spre, scur) = (sbis, sbis);
            }

            (xpre, fpre) = (xcur, fcur);
            xcur += if scur.abs() > delta {
                scur
            } else {
                delta.copysign(sbis)
            };

            fcur = f(xcur).map_err(Error::Function)?;

            if fcur.is_nan() {
                return Err(Error::Nan(NanError));
            }
        }

        Err(Error::Convergence(ConvergenceError {
            calls: false,
            val: self.max_iters,
        }))
    }
}

/// Errors from [`Brent::root`].
#[derive(Debug)]
pub enum Error<E = Infallible> {
    /// A tolerance parameter was invalid.
    InvalidTol(InvalidTol),

    /// `f(start)` and `f(end)` are of the same sign.
    SameSigns(SameSigns),

    /// Failed to converge after a number of iterations.
    Convergence(ConvergenceError),

    /// The function returned `NaN`.
    Nan(NanError),

    /// An error from the function being optimized.
    Function(E),
}

impl<E> Display for Error<E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("Brent root-finding failed")
    }
}

impl<E: 'static + std::error::Error> std::error::Error for Error<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidTol(e) => Some(e),
            Self::SameSigns(e) => Some(e),
            Self::Convergence(e) => Some(e),
            Self::Nan(e) => Some(e),
            Self::Function(e) => Some(e),
        }
    }
}

/// `f(start)` and `f(end)` are of the same sign.
#[derive(Debug)]
#[non_exhaustive]
pub struct SameSigns {
    start: f64,
    f_start: f64,
    end: f64,
    f_end: f64,
}

impl Display for SameSigns {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "`f(start)` and `f(end)` are not of different signs: f({:?}) = {:?} and f({:?}) = {:?}",
            self.start, self.f_start, self.end, self.f_end
        )
    }
}

impl std::error::Error for SameSigns {}

use super::ConvergenceError;
use super::DEFAULT_ABS_TOL;
use super::DEFAULT_MAX_ITERS;
use super::DEFAULT_REL_TOL;
use super::InvalidTol;
use super::NanError;
use std::convert::Infallible;
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;

//! Optimization.

pub mod nelder_mead;
pub use nelder_mead::NelderMead;

/// The minimum relative tolerance.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L11
pub const MIN_REL_TOL: f64 = 4.0 * f64::EPSILON;

/// Error when the provided tolerance levels are invalid.
///
/// Absolute tolerances must be positive, and relative tolerances must be above [`MIN_REL_TOL`].
#[derive(Debug)]
pub struct InvalidTol {
    rel: bool,
    f: bool,
    tol: f64,
}

impl InvalidTol {
    // callable from the `root` module
    pub(crate) fn check_abs(f: bool, tol: f64) -> Result<(), InvalidTol> {
        if !(tol > 0.0) {
            return Err(Self { rel: false, f, tol });
        }
        Ok(())
    }
    pub(crate) fn check_rel(f: bool, tol: f64) -> Result<(), InvalidTol> {
        if !(tol >= MIN_REL_TOL) {
            return Err(Self { rel: true, f, tol });
        }
        Ok(())
    }
}

impl Display for InvalidTol {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self.rel {
            false => f.write_str("absolute ")?,
            true => f.write_str("relative ")?,
        }
        match self.f {
            false => f.write_str("x")?,
            true => f.write_str("f(…)")?,
        }
        write!(f, "-tolerance {:?} is ", self.tol)?;
        match self.rel {
            false => f.write_str("nonpositive")?,
            true => write!(f, "not above minimum ({MIN_REL_TOL:?})")?,
        }
        Ok(())
    }
}

impl Error for InvalidTol {}

/// Error when the lower bound is not less than the upper bound.
#[derive(Debug)]
pub struct InvalidBound {
    i: usize,
    lower: f64,
    upper: f64,
}

impl Display for InvalidBound {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "invalid bound on parameter {}: {:?} is not less than {:?}",
            self.i, self.lower, self.upper
        )
    }
}

impl Error for InvalidBound {}

/// Error when convergence failed after a certain number of iterations or function calls.
#[derive(Debug)]
#[non_exhaustive]
pub struct ConvergenceError {
    // constructible from the `root` module
    pub(crate) calls: bool,
    pub(crate) val: u64,
}

impl Display for ConvergenceError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self.calls {
            false => write!(f, "convergence failed after {} iterations", self.val),
            true => write!(f, "convergence failed after {} calls", self.val),
        }
    }
}

impl Error for ConvergenceError {}

/// An error caused by the function being optimized returning `NaN`.
#[derive(Debug)]
#[non_exhaustive]
pub struct NanError;

impl Display for NanError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("function returned NaN")
    }
}

impl Error for NanError {}

use core::fmt;
use std::error::Error;
use std::fmt::Display;
use std::fmt::Formatter;

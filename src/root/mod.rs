//! Root-finding.

pub mod brent;
pub use brent::Brent;

/// The default absolute tolerance.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L10
pub const DEFAULT_ABS_TOL: f64 = 2e-12;

/// The default relative tolerance.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L11
pub const DEFAULT_REL_TOL: f64 = 4.0 * f64::EPSILON;

/// The default maximum number of iterations before convergence fails.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L9
pub const DEFAULT_MAX_ITERS: u64 = 100;

/// Error when the provided tolerance levels are invalid.
#[derive(Debug)]
pub struct InvalidTol {
    // constructible from the `opt` module
    pub(crate) rel: bool,
    pub(crate) f: bool,
    pub(crate) tol: f64,
}

impl InvalidTol {
    pub(crate) fn check_abs(f: bool, tol: f64) -> Result<(), InvalidTol> {
        if !(tol > 0.0) {
            return Err(Self { rel: false, f, tol });
        }
        Ok(())
    }
    pub(crate) fn check_rel(f: bool, tol: f64) -> Result<(), InvalidTol> {
        if !(tol >= DEFAULT_REL_TOL) {
            return Err(Self { rel: false, f, tol });
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
            true => write!(f, "not above minimum ({DEFAULT_REL_TOL:?})")?,
        }
        Ok(())
    }
}

impl Error for InvalidTol {}

/// Error when convergence failed after a certain number of iterations.
#[derive(Debug)]
#[non_exhaustive]
// constructible from the `opt` module
pub struct ConvergenceError(pub(crate) u64);

impl Display for ConvergenceError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "convergence failed after {} iterations", self.0)
    }
}

impl Error for ConvergenceError {}

use std::error::Error;
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;

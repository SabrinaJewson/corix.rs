//! Root-finding.

pub mod brent;
pub use brent::Brent;

/// The default absolute tolerance.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L10
pub const DEFAULT_ABS_TOL: f64 = 2e-12;

/// The default relative tolerance.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L11
pub const DEFAULT_REL_TOL: f64 = MIN_REL_TOL;

/// The default maximum number of iterations before convergence fails.
// https://github.com/scipy/scipy/blob/v1.16.0/scipy/optimize/_zeros_py.py#L9
pub const DEFAULT_MAX_ITERS: u64 = 100;

pub use crate::opt::ConvergenceError;
pub use crate::opt::InvalidTol;
pub use crate::opt::MIN_REL_TOL;
pub use crate::opt::NanError;

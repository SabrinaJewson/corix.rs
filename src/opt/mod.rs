//! Optimization.

pub use crate::root::ConvergenceError;
pub use crate::root::InvalidTol;

pub mod nelder_mead;
pub use nelder_mead::NelderMead;

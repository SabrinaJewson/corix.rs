//! Special functions.

mod libm;
pub use libm::*;

mod chebyshev;
pub use chebyshev::*;

mod polynomial;
pub use polynomial::*;

pub mod airy;
pub use airy::*;

pub mod ellip;

pub mod bessel;

#![doc = include_str!("../README.md")]
#![warn(clippy::pedantic)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::match_bool)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::items_after_test_module)]

pub mod opt;
pub mod root;
pub mod util;

/// An error caused by a float value becoming `NaN`.
///
/// In general, this is only returned from “high-level” routines that are already fallible.
/// For more low-level/infallible routines, we just propagate `NaN`s.
// TODO: move to opt module
#[derive(Debug)]
#[non_exhaustive]
pub struct NanError;

impl Display for NanError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("found NaN")
    }
}

impl Error for NanError {}

use std::error::Error;
use std::fmt;
use std::fmt::Display;
use std::fmt::Formatter;

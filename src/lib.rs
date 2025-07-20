#![doc = include_str!("../README.md")]
#![warn(clippy::pedantic)]
#![warn(missing_docs)]
#![warn(missing_debug_implementations)]
#![allow(mixed_script_confusables)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::items_after_test_module)]
#![allow(clippy::match_bool)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
// TODO: libm

pub mod opt;
pub mod root;

mod util;

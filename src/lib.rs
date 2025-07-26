#![doc = include_str!("../README.md")]
#![warn(clippy::pedantic)]
#![warn(missing_docs)]
#![warn(missing_debug_implementations)]
#![warn(unused_qualifications)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::items_after_test_module)]
#![allow(clippy::match_bool)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(mixed_script_confusables)]
// TODO: libm

pub mod opt;
pub mod root;

pub mod special;
pub use special::*;

pub mod stats;
pub use stats::*;

mod xsf {
    use num_complex::Complex as std_complex;

    include!(concat!(env!("OUT_DIR"), "/xsf_bindings.rs"));
}

mod util;

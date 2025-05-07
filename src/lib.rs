//! A stub for a Rust crate for BED file manipulation
//! 

#![warn(rust_2021_compatibility)]
#![warn(rust_2018_idioms)]

pub mod extract;
pub mod merge;
pub mod structs;

pub use crate::extract::*;
pub use crate::merge::*;
pub use crate::structs::*;
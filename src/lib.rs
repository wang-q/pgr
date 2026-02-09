// Add these imports to use the stdsimd library
#![feature(portable_simd)]

pub mod libs;
pub use libs::io::{reader, writer};

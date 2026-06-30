// Add these imports to use the stdsimd library
#![feature(portable_simd)]

pub mod libs;
pub use libs::alignment::coords::reverse_range;
pub use libs::io::{is_bgzf, read_lines, read_sizes, reader, writer};

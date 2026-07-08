//! pbit: native population genome compression format (2bit record + delta layer).
//!
//! Borrows the LZ-diff V2 and segment-level reference compression algorithms
//! from AGC v3.2.3, but uses a native `.pbit` file format (fixed-size
//! little-endian integers, no varint/prefix coding). The reference layer
//! reuses standard 2bit records via `libs::fmt::twobit::read_2bit_record` /
//! `write_2bit_record`.

pub mod format;
pub mod lz_diff;

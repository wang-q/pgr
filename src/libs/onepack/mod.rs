//! ONEcode binary container and codecs used by FastGA's `.1aln` format.
//!
//! This module ports the ONElib container (header, binary line I/O, object
//! index, footer) and its integer / Huffman / DNA codecs. It is self-contained
//! and does not depend on the rest of pgr, so it can be unit-tested in
//! isolation against FastGA's golden `.1aln` files.

pub mod container;
pub mod expand;
pub mod ltf;
pub mod record;
pub mod schema;
pub mod vc;
pub mod write;

//! Overlap-layout-consensus assembly on unitig pseudo-reads.
//!
//! Design: `notes/design/olc.md`. The pipeline treats unitigs produced at
//! several k values as pseudo-reads, finds exact overlaps between them
//! (stage S1), chains them into layouts (S2), and stitches each layout into
//! a consensus contig (S3).

pub mod consensus;
pub mod layout;
pub mod overlap;

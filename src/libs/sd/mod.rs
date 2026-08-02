//! Segmental duplication (SD) detection and analysis.

pub mod cluster;
pub mod cover;
pub mod decompose;
pub mod search_lastz;
pub mod search_pgi;

use crate::libs::fmt::psl::Psl;

/// Alignment block length of a PSL record
/// (matches + mismatches + repeats + Ns + query/target insert bases).
pub fn psl_block_len(p: &Psl) -> u32 {
    p.match_count
        + p.mismatch_count
        + p.rep_match
        + p.n_count
        + p.q_base_insert.max(0) as u32
        + p.t_base_insert.max(0) as u32
}

/// Block identity of a PSL record: `(matches + repeats) / block_len`.
pub fn psl_identity(p: &Psl) -> f64 {
    let blk = psl_block_len(p);
    if blk == 0 {
        0.0
    } else {
        (p.match_count + p.rep_match) as f64 / blk as f64
    }
}

/// Whether a PSL record passes the SD filters (min block length + identity).
pub fn passes_sd_filters(p: &Psl, min_len: u32, min_identity: f64) -> bool {
    psl_block_len(p) >= min_len && psl_identity(p) >= min_identity
}

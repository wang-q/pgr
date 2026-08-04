//! Segmental duplication (SD) detection and analysis.

pub mod cluster;
pub mod cover;
pub mod decompose;
pub mod search_lastz;
pub mod search_pgi;

use crate::libs::fmt::psl::Psl;
use std::io::Read;

/// Whether `path` (a file) looks like a `.pgi` index, by magic or extension.
///
/// The SD search filters score the alignment blocks; a `.pgi` input aligns
/// without extension sequences and every block scores 0, so the search would
/// silently return nothing. Refuse it up front with a clear error instead.
pub fn is_pgi_input(path: &str) -> bool {
    if let Ok(mut f) = std::fs::File::open(path) {
        let mut magic = [0u8; 4];
        if f.read_exact(&mut magic).is_ok() {
            return &magic == crate::libs::pgi::PGI_MAGIC;
        }
    }
    std::path::Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        == Some("pgi")
}

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

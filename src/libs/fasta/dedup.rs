//! FASTA deduplication signature computation.
//!
//! Provides a single [`record_signature`] entry point that computes a 64-bit
//! rapidhash over the chosen record field (sequence / description / name),
//! with optional case-insensitivity and both-strand comparison.

/// Dedup signature options.
///
/// `is_seq` / `is_desc` / name-only selection follows the precedence:
/// `is_seq` → `is_desc` (with description) → name-only.
///
/// `is_both` (both-strand comparison) only applies to sequence mode and
/// implies case-insensitive comparison.
#[derive(Debug, Clone, Copy, Default)]
pub struct DedupOptions {
    /// Compare by sequence content.
    pub is_seq: bool,
    /// Compare by name + description.
    pub is_desc: bool,
    /// Compare both strands (forward and reverse complement); implies
    /// case-insensitive. Only meaningful in sequence mode.
    pub is_both: bool,
    /// Case-insensitive comparison.
    pub is_insensitive: bool,
}

/// Compute a 64-bit rapidhash signature for a FASTA record.
///
/// Returns the minimum of the forward and reverse-complement hashes when
/// `opts.is_both` is set (sequence mode only).
pub fn record_signature(
    name: &[u8],
    desc: Option<&[u8]>,
    seq: &[u8],
    opts: &DedupOptions,
) -> anyhow::Result<u64> {
    Ok(if opts.is_seq {
        if opts.is_both {
            let fwd = rapidhash::rapidhash(&seq.to_ascii_uppercase());
            let rc: Vec<u8> = crate::libs::nt::rev_comp(seq).collect();
            let rev = rapidhash::rapidhash(&rc.to_ascii_uppercase());
            fwd.min(rev)
        } else if opts.is_insensitive {
            rapidhash::rapidhash(&seq.to_ascii_uppercase())
        } else {
            rapidhash::rapidhash(seq)
        }
    } else if opts.is_desc && desc.is_some() {
        let full = [name, desc.unwrap()].concat();
        if opts.is_insensitive {
            rapidhash::rapidhash(&full.to_ascii_uppercase())
        } else {
            rapidhash::rapidhash(&full)
        }
    } else if opts.is_insensitive {
        rapidhash::rapidhash(&name.to_ascii_uppercase())
    } else {
        rapidhash::rapidhash(name)
    })
}

//! FASTA record filtering and formatting.
//!
//! Provides size / N-count / uniqueness filters and a per-base formatter
//! (uppercase / IUPAC collapsing / dash stripping) shared by `pgr fa filter`.

use std::collections::BTreeSet;

/// Sentinel value indicating "no limit" for a size/N-count filter.
pub const NO_LIMIT: usize = usize::MAX;

/// Check whether a sequence passes size, N-count, and uniqueness filters.
///
/// `seq_len` is the sequence length; `seq` is the raw bases (for N-count).
/// When `is_uniq` is true, `name` is inserted into `seen` and the record is
/// rejected if it was already present.
pub fn pass_filters(
    seq: &[u8],
    minsize: usize,
    maxsize: usize,
    maxn: usize,
    is_uniq: bool,
    seen: &mut BTreeSet<String>,
    name: &str,
) -> bool {
    let len = seq.len();
    if minsize != NO_LIMIT && len < minsize {
        return false;
    }
    if maxsize != NO_LIMIT && len > maxsize {
        return false;
    }
    if maxn != NO_LIMIT && crate::libs::nt::count_n(seq) > maxn {
        return false;
    }
    if is_uniq && !seen.insert(name.to_string()) {
        return false;
    }
    true
}

/// Format a sequence by optionally stripping dashes, collapsing IUPAC codes
/// to `N`, and upper-casing the result.
pub fn format_sequence(seq: &[u8], is_dash: bool, is_iupac: bool, is_upper: bool) -> String {
    let mut out = String::with_capacity(seq.len());
    for &nt in seq {
        if is_dash && nt == b'-' {
            continue;
        }
        if is_iupac {
            let c = char::from(crate::libs::nt::to_n(nt));
            out.push(if is_upper { c.to_ascii_uppercase() } else { c });
        } else if is_upper {
            out.push(char::from(nt).to_ascii_uppercase());
        } else {
            out.push(char::from(nt));
        }
    }
    out
}

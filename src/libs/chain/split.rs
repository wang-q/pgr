//! Helpers for chain-split naming.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Derive a 3-digit lump bucket name from a sequence name.
///
/// Scans `name` for the first ASCII-digit run and returns `val % lump` formatted
/// as 3 digits. If no digits are found, falls back to a stable hash of the name.
pub fn lump_name(name: &str, lump: usize) -> String {
    // Look for integer part of name
    let mut s = name;
    while let Some(idx) = s.find(|c: char| c.is_ascii_digit()) {
        s = &s[idx..];
        let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        let digits = &s[..end];
        if let Ok(val) = digits.parse::<usize>() {
            return format!("{:03}", val % lump);
        }
        s = &s[end..];
    }

    // If no digits found, hash it
    let mut hasher = DefaultHasher::new();
    name.hash(&mut hasher);
    let hash = hasher.finish();
    format!("{:03}", (hash as usize) % lump)
}

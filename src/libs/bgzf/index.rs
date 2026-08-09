//! `.gzi` index: (compressed offset, uncompressed offset) block entries.

use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

use anyhow::Context;

/// A BGZF `.gzi` index (bgzip-compatible).
///
/// Entries are `(compressed_offset, uncompressed_offset)` pairs for each
/// non-empty block except the first (offset 0, 0 is implicit), matching the
/// format written by `bgzip -i`.
#[derive(Debug, Default)]
pub struct GziIndex {
    entries: Vec<(u64, u64)>,
}

impl GziIndex {
    /// Builds an index from block entries.
    pub fn from_entries(entries: Vec<(u64, u64)>) -> Self {
        Self { entries }
    }

    /// Returns the number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns true when the index has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Maps an uncompressed offset to a virtual position `(compressed offset,
    /// offset within the block)`.
    pub fn query(&self, pos: u64) -> Option<(u64, u16)> {
        let i = self.entries.partition_point(|&(_, upos)| upos <= pos);
        let (cpos, upos) = if i == 0 { (0, 0) } else { self.entries[i - 1] };
        let block_pos = u16::try_from(pos - upos).ok()?;
        Some((cpos, block_pos))
    }

    /// Reads an index from a file.
    pub fn read(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        let mut file = File::open(path.as_ref())
            .with_context(|| format!("could not open {}", path.as_ref().display()))?;
        let mut buf = Vec::new();
        file.read_to_end(&mut buf)
            .context("could not read gzi index")?;
        Self::from_bytes(&buf).context("invalid gzi index")
    }

    /// Writes the index to a file.
    pub fn write(&self, path: impl AsRef<Path>) -> anyhow::Result<()> {
        let mut file = File::create(path.as_ref())
            .with_context(|| format!("could not create {}", path.as_ref().display()))?;
        file.write_all(&self.to_bytes())
            .context("could not write gzi index")?;
        Ok(())
    }

    /// Serializes the index (u64 little-endian entry count, then pairs).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(8 + self.entries.len() * 16);
        out.extend_from_slice(&(self.entries.len() as u64).to_le_bytes());
        for &(cpos, upos) in &self.entries {
            out.extend_from_slice(&cpos.to_le_bytes());
            out.extend_from_slice(&upos.to_le_bytes());
        }
        out
    }

    /// Parses serialized index bytes.
    pub fn from_bytes(bytes: &[u8]) -> io::Result<Self> {
        if bytes.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "gzi index too short",
            ));
        }
        let n = u64::from_le_bytes(bytes[..8].try_into().unwrap()) as usize;
        if bytes.len() != 8 + n * 16 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "gzi index size mismatch",
            ));
        }
        let mut entries = Vec::with_capacity(n);
        for i in 0..n {
            let off = 8 + i * 16;
            let cpos = u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap());
            let upos = u64::from_le_bytes(bytes[off + 8..off + 16].try_into().unwrap());
            entries.push((cpos, upos));
        }
        Ok(Self { entries })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let index = GziIndex::from_entries(vec![(100, 65536), (200, 131072)]);
        let bytes = index.to_bytes();
        let parsed = GziIndex::from_bytes(&bytes).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed.query(0), Some((0, 0)));
        assert_eq!(parsed.query(10_000), Some((0, 10_000)));
        assert_eq!(parsed.query(70_000), Some((100, 4464)));
        assert_eq!(parsed.query(140_000), Some((200, 8928)));
    }
}

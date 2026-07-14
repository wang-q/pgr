//! A simple bit-vector for tracking used coordinate ranges on a sequence.

/// A bitmap over a fixed `size` of bits, supporting range set/test operations.
pub struct BitMap {
    /// Number of bits tracked (logical size).
    pub size: u64,
    bits: Vec<u64>,
}

impl BitMap {
    /// Create a new bitmap of `size` bits, all cleared.
    pub fn new(size: u64) -> Self {
        let num_words = size.div_ceil(64);
        Self {
            size,
            bits: vec![0; num_words as usize],
        }
    }

    /// Set `len` bits starting at `start` (clamped to `size`).
    pub fn set_range(&mut self, start: u64, len: u64) {
        if len == 0 {
            return;
        }
        let end = start.saturating_add(len).min(self.size);
        let start = start.min(self.size);
        if start >= end {
            return;
        }

        let start_word = (start / 64) as usize;
        let end_word = ((end - 1) / 64) as usize;

        let start_bit = start % 64;
        let end_bit = (end - 1) % 64;

        if start_word == end_word {
            let mask = (!0u64 << start_bit) & (!0u64 >> (63 - end_bit));
            self.bits[start_word] |= mask;
        } else {
            // First word
            self.bits[start_word] |= !0u64 << start_bit;

            // Middle words
            for i in (start_word + 1)..end_word {
                self.bits[i] = !0u64;
            }

            // Last word
            self.bits[end_word] |= !0u64 >> (63 - end_bit);
        }
    }

    /// Return true iff all bits in `[start, start+len)` are set.
    pub fn is_fully_set(&self, start: u64, len: u64) -> bool {
        if len == 0 {
            return true;
        }
        let end = start.saturating_add(len).min(self.size);
        let start = start.min(self.size);
        if start >= end {
            return true;
        }

        let start_word = (start / 64) as usize;
        let end_word = ((end - 1) / 64) as usize;

        let start_bit = start % 64;
        let end_bit = (end - 1) % 64;

        if start_word == end_word {
            let mask = (!0u64 << start_bit) & (!0u64 >> (63 - end_bit));
            return (self.bits[start_word] & mask) == mask;
        } else {
            // First word
            let mask1 = !0u64 << start_bit;
            if (self.bits[start_word] & mask1) != mask1 {
                return false;
            }

            // Middle words
            for i in (start_word + 1)..end_word {
                if self.bits[i] != !0u64 {
                    return false;
                }
            }

            // Last word
            let mask2 = !0u64 >> (63 - end_bit);
            if (self.bits[end_word] & mask2) != mask2 {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_range_does_not_overflow() {
        let mut bm = BitMap::new(100);
        // start near u64::MAX with a length that would overflow if added directly.
        bm.set_range(u64::MAX - 5, 10);
        // Bits are clamped to size, so nothing beyond the logical size is set.
        assert!(!bm.is_fully_set(0, 100));
    }

    #[test]
    fn test_is_fully_set_does_not_overflow() {
        let bm = BitMap::new(100);
        // Query range whose end would overflow if computed with plain addition.
        assert!(bm.is_fully_set(u64::MAX - 5, 10));
    }
}

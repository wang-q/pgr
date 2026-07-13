//! Size-based chunking state machine for FASTA splitting.
//!
//! Tracks cumulative sequence length and record count to decide when to
//! rotate output files. Used by `pgr fa split about`.

/// Chunker that rotates output files after a size or record-count threshold.
///
/// `max_size` is the approximate byte threshold per file (used in size mode).
/// When `is_even` is true, files are rotated every 2 records (binary split).
/// `max_files` caps the number of output files; `max_files_exceeded()`
/// signals the caller to stop.
pub struct SizeChunker {
    max_size: usize,
    is_even: bool,
    max_files: usize,
    cur_size: usize,
    record_count: usize,
    file_index: usize,
}

impl SizeChunker {
    /// Create a new chunker. `max_size` is ignored when `is_even` is true.
    pub fn new(max_size: usize, is_even: bool, max_files: usize) -> Self {
        Self {
            max_size,
            is_even,
            max_files,
            cur_size: 0,
            record_count: 0,
            file_index: 0,
        }
    }

    /// Current output file index (0-based).
    pub fn file_index(&self) -> usize {
        self.file_index
    }

    /// Returns true once the file index has reached `max_files`.
    pub fn max_files_exceeded(&self) -> bool {
        self.file_index >= self.max_files
    }

    /// Account for a record of `seq_len` bytes that was just written.
    /// Rotates to the next file if the threshold is reached.
    pub fn advance(&mut self, seq_len: usize) {
        self.cur_size += seq_len;
        self.record_count += 1;
        let rotate = if self.is_even {
            self.record_count.is_multiple_of(2)
        } else {
            self.cur_size > self.max_size
        };
        if rotate {
            self.cur_size = 0;
            self.record_count = 0;
            self.file_index += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_mode_rotates_after_threshold() {
        let mut c = SizeChunker::new(100, false, 999);
        c.advance(60);
        assert_eq!(c.file_index(), 0);
        assert!(!c.max_files_exceeded());
        // 60 + 50 = 110 > 100 → rotate
        c.advance(50);
        assert_eq!(c.file_index(), 1);
    }

    #[test]
    fn test_even_mode_rotates_every_two_records() {
        let mut c = SizeChunker::new(0, true, 999);
        c.advance(10);
        assert_eq!(c.file_index(), 0);
        c.advance(20);
        assert_eq!(c.file_index(), 1);
        c.advance(30);
        assert_eq!(c.file_index(), 1);
        c.advance(40);
        assert_eq!(c.file_index(), 2);
    }

    #[test]
    fn test_max_files_exceeded() {
        let mut c = SizeChunker::new(10, false, 2);
        // file 0: 5 bytes, rotate (5 > 10? no) → stay
        c.advance(5);
        assert_eq!(c.file_index(), 0);
        assert!(!c.max_files_exceeded());
        // file 0: 5 + 6 = 11 > 10 → rotate to file 1
        c.advance(6);
        assert_eq!(c.file_index(), 1);
        assert!(!c.max_files_exceeded());
        // file 1: 11 > 10 → rotate to file 2 (reached max_files)
        c.advance(11);
        assert_eq!(c.file_index(), 2);
        assert!(c.max_files_exceeded());
    }
}

//! Size-based chunking state machine for FASTA splitting.
//!
//! Tracks cumulative sequence length and record count to decide when to
//! rotate output files. Used by `pgr fa split about`.

/// Chunker that rotates output files after a size or record-count threshold.
///
/// `max_size` is the approximate byte threshold per file (used in size mode).
/// When `is_even` is true, rotation is deferred until the current file holds
/// an even number of records (so paired reads stay together).
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
    /// Create a new chunker. Rotation is size-based; when `is_even` is true it
    /// is deferred until the current file holds an even number of records.
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
    /// Rotates to the next file once the size threshold is reached; when
    /// `is_even` is set, rotation is deferred until the record count is even.
    pub fn advance(&mut self, seq_len: usize) {
        self.cur_size += seq_len;
        self.record_count += 1;
        let size_reached = self.cur_size > self.max_size;
        let rotate = if self.is_even {
            size_reached && self.record_count.is_multiple_of(2)
        } else {
            size_reached
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
    fn test_even_mode_rotates_when_size_reached_and_even() {
        let mut c = SizeChunker::new(100, true, 999);
        // 60 < 100, count 1 (odd) → stay
        c.advance(60);
        assert_eq!(c.file_index(), 0);
        // 60 + 50 = 110 > 100, count 2 (even) → rotate
        c.advance(50);
        assert_eq!(c.file_index(), 1);
    }

    #[test]
    fn test_even_mode_holds_until_record_count_is_even() {
        let mut c = SizeChunker::new(100, true, 999);
        // 110 > 100 but count 1 (odd) → hold, do not rotate
        c.advance(110);
        assert_eq!(c.file_index(), 0);
        // 110 + 10 = 120 > 100, count 2 (even) → rotate
        c.advance(10);
        assert_eq!(c.file_index(), 1);
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

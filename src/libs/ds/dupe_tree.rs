//! A depth-tracking interval tree for signed 1D intervals.

/// A non-overlapping run of constant integer depth.
pub struct Segment {
    /// Start coordinate (0-based, inclusive).
    pub start: u64,
    /// End coordinate (0-based, exclusive).
    pub end: u64,
    /// Depth over this interval.
    pub depth: i32,
}

/// Interval tree tracking signed coverage depth over a 1D coordinate line.
///
/// Intervals are added with positive or negative weights. After all intervals
/// are recorded, call [`DupeTree::build`] to flatten them into non-overlapping
/// constant-depth segments. Then [`DupeTree::count_over`] can report how many
/// bases within a query range have depth above a threshold.
pub struct DupeTree {
    intervals: Vec<(u64, u64, i32)>,
    segments: Vec<Segment>,
}

impl Default for DupeTree {
    fn default() -> Self {
        Self::new()
    }
}

impl DupeTree {
    /// Creates an empty DupeTree.
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
            segments: Vec::new(),
        }
    }

    /// Records a positive depth contribution over `[start, end)`.
    pub fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, 1));
        }
    }

    /// Records a negative depth contribution over `[start, end)`.
    pub fn subtract(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, -1));
        }
    }

    /// Flattens recorded intervals into non-overlapping constant-depth segments.
    pub fn build(&mut self) {
        if self.intervals.is_empty() {
            return;
        }

        let mut events = Vec::new();
        for (s, e, d) in &self.intervals {
            events.push((*s, *d));
            events.push((*e, -*d));
        }
        // Sort by position, then by delta so that end events (negative delta)
        // are processed before start events (positive delta) at the same coordinate.
        events.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        let mut current_depth = 0;
        let mut segments = Vec::new();

        for i in 0..events.len() - 1 {
            let (pos, delta) = events[i];
            current_depth += delta;

            let next_pos = events[i + 1].0;
            if next_pos > pos {
                segments.push(Segment {
                    start: pos,
                    end: next_pos,
                    depth: current_depth,
                });
            }
        }

        self.segments = segments;
    }

    /// Returns total bases covered by segments with `depth >= threshold`.
    pub fn count_over(&self, start: u64, end: u64, threshold: i32) -> u64 {
        if self.segments.is_empty() {
            return 0;
        }

        // Binary search for first segment ending after start
        let idx = self.segments.partition_point(|seg| seg.end <= start);

        let mut count = 0;
        for seg in &self.segments[idx..] {
            if seg.start >= end {
                break;
            }

            if seg.depth >= threshold {
                let s = seg.start.max(start);
                let e = seg.end.min(end);
                if s < e {
                    count += e - s;
                }
            }
        }
        count
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dupe_tree_adjacent_intervals() {
        let mut dt = DupeTree::new();
        dt.add(0, 10);
        dt.add(10, 20);
        dt.build();

        // Each adjacent interval contributes depth 1 on its own range.
        assert_eq!(dt.count_over(0, 20, 1), 20);
        assert_eq!(dt.count_over(0, 20, 2), 0);
    }

    #[test]
    fn test_dupe_tree_overlapping_intervals() {
        let mut dt = DupeTree::new();
        dt.add(0, 15);
        dt.add(10, 25);
        dt.build();

        assert_eq!(dt.count_over(0, 25, 1), 25);
        assert_eq!(dt.count_over(0, 25, 2), 5);
        assert_eq!(dt.count_over(0, 10, 2), 0);
        assert_eq!(dt.count_over(20, 25, 2), 0);
    }

    #[test]
    fn test_dupe_tree_subtract() {
        let mut dt = DupeTree::new();
        dt.add(0, 20);
        dt.subtract(5, 15);
        dt.build();

        assert_eq!(dt.count_over(0, 20, 1), 10);
        assert_eq!(dt.count_over(0, 20, 2), 0);
        assert_eq!(dt.count_over(5, 15, 1), 0);
    }
}

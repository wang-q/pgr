//! Synteny classification helper: depth-tracking interval tree.
//!
//! `DupeTree` accumulates signed intervals (added by fills, subtracted by
//! nested gaps) and flattens them into non-overlapping `Segment`s of constant
//! depth, so a fill's query-overlap with duplications can be queried.

/// A non-overlapping run of constant duplication depth.
pub struct Segment {
    pub start: u64,
    pub end: u64,
    pub depth: i32,
}

/// Interval tree tracking query-side duplication depth for synteny classification.
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

    /// Records a +1 depth contribution over `[start, end)`.
    pub fn add(&mut self, start: u64, end: u64) {
        if start < end {
            self.intervals.push((start, end, 1));
        }
    }

    /// Records a -1 depth contribution over `[start, end)`.
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
        // Sort by position, then by delta (to process all updates at same pos)
        // Actually, order of processing at same position doesn't matter for final segments between positions.
        events.sort_by_key(|a| a.0);

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

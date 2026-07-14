//! Top-K category purity detector.
//!
//! Penalizes a score when the largest K categories account for more than a
//! configurable ratio of the total observations. This is useful for detecting
//! low-complexity or dominated distributions in sequences, samples, or any
//! categorical data.

/// Tracks category counts and reports whether the top K categories dominate.
#[derive(Clone, Debug)]
pub struct TopKPurity {
    counts: Vec<usize>,
    total: usize,
    k: usize,
    ok_ratio: f64,
}

impl TopKPurity {
    /// Creates a new detector.
    ///
    /// # Arguments
    ///
    /// * `num_classes` - Number of distinct categories.
    /// * `k` - Number of top categories to monitor.
    /// * `ok_ratio` - Maximum acceptable ratio for the top K categories.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if `k == 0` or `k > num_classes`.
    pub fn new(num_classes: usize, k: usize, ok_ratio: f64) -> Self {
        debug_assert!(k > 0, "k must be greater than 0");
        debug_assert!(k <= num_classes, "k must not exceed num_classes");
        Self {
            counts: vec![0; num_classes],
            total: 0,
            k,
            ok_ratio,
        }
    }

    /// Records one observation for `class`.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if `class` is out of range.
    pub fn increment(&mut self, class: usize) {
        debug_assert!(class < self.counts.len(), "class out of range");
        self.counts[class] += 1;
        self.total += 1;
    }

    /// Returns the total number of observations recorded.
    pub fn total(&self) -> usize {
        self.total
    }

    /// Returns the ratio of the top K categories to the total observations.
    ///
    /// Returns `0.0` when no observations have been recorded.
    pub fn observed_ratio(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let top_k = self.top_k_sum();
        top_k as f64 / self.total as f64
    }

    /// Returns `true` if the observed top-K ratio is within the acceptable bound.
    pub fn is_acceptable(&self) -> bool {
        self.observed_ratio() <= self.ok_ratio
    }

    /// Returns the penalty factor if the distribution is too dominated.
    ///
    /// Returns `None` when `observed_ratio <= ok_ratio` (no penalty).
    /// Otherwise returns `Some(1.01 - (observed - ok_ratio) / (1 - ok_ratio))`.
    pub fn penalty_factor(&self) -> Option<f64> {
        let observed = self.observed_ratio();
        if observed <= self.ok_ratio {
            return None;
        }
        let factor = 1.01 - (observed - self.ok_ratio) / (1.0 - self.ok_ratio);
        Some(factor)
    }

    fn top_k_sum(&self) -> usize {
        let mut sorted = self.counts.clone();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        sorted.iter().take(self.k).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_is_acceptable() {
        let detector = TopKPurity::new(4, 2, 0.80);
        assert!(detector.is_acceptable());
        assert_eq!(detector.penalty_factor(), None);
    }

    #[test]
    fn test_below_threshold() {
        let mut detector = TopKPurity::new(4, 2, 0.80);
        // 5/8 for top 2
        detector.increment(0);
        detector.increment(0);
        detector.increment(1);
        detector.increment(1);
        detector.increment(2);
        detector.increment(2);
        detector.increment(3);
        detector.increment(3);
        assert!(detector.is_acceptable());
        assert_eq!(detector.penalty_factor(), None);
    }

    #[test]
    fn test_above_threshold() {
        let mut detector = TopKPurity::new(4, 2, 0.80);
        // 9/10 for top 2
        for _ in 0..5 {
            detector.increment(0);
        }
        for _ in 0..4 {
            detector.increment(1);
        }
        detector.increment(2);
        assert!(!detector.is_acceptable());
        let factor = detector.penalty_factor().unwrap();
        assert!(factor < 1.01);
    }
}

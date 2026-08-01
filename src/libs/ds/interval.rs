//! General interval utilities.

/// Merge overlapping or adjacent intervals in a list (sorted internally),
/// returning a non-overlapping list.
pub fn merge_intervals(intervals: &mut [(u64, u64)]) -> Vec<(u64, u64)> {
    intervals.sort_unstable();

    let mut merged: Vec<(u64, u64)> = Vec::new();
    for &(s, e) in intervals.iter() {
        if let Some(last) = merged.last_mut() {
            if s <= last.1 {
                last.1 = last.1.max(e);
                continue;
            }
        }
        merged.push((s, e));
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_overlapping_and_unsorted() {
        let mut intervals = vec![(10, 20), (0, 5), (4, 8), (15, 25), (30, 35)];
        assert_eq!(
            merge_intervals(&mut intervals),
            vec![(0, 8), (10, 25), (30, 35)]
        );
    }

    #[test]
    fn test_adjacent_merged() {
        let mut intervals = vec![(0, 10), (10, 20)];
        assert_eq!(merge_intervals(&mut intervals), vec![(0, 20)]);
    }

    #[test]
    fn test_nested_intervals() {
        let mut intervals = vec![(0, 100), (20, 30), (40, 60)];
        assert_eq!(merge_intervals(&mut intervals), vec![(0, 100)]);
    }
}

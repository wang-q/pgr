//! Set of integer spans (`IntSpan`), vendored from the external `intspan`
//! crate (kept API-identical so the two can be swapped without touching
//! call sites).

use anyhow::anyhow;
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::vec::Vec;

/// `IntSpan` handles of sets containing integer spans.
///
/// # SYNOPSIS
///
/// ```
/// use pgr::libs::ds::IntSpan;
///
/// let mut ints = IntSpan::new();
/// for i in vec![1, 2, 3, 5, 7, 9] {
///     ints.add_n(i);
/// }
/// ints.add_pair(100, 10000);
/// ints.remove_n(1000);
///
/// let expected = "1-3,5,7,9,100-999,1001-10000";
/// assert_eq!(ints.to_string(), expected);
/// assert_eq!(ints.cardinality(), 9906);
/// assert_eq!(ints.is_empty(), false);
/// assert_eq!(ints.is_universal(), false);
/// assert_eq!(ints.is_infinite(), false);
/// assert_eq!(ints.is_finite(), true);
/// assert_eq!(ints.is_pos_inf(), false);
/// assert_eq!(ints.is_neg_inf(), false);
/// ```
///
/// ```
/// # use pgr::libs::ds::IntSpan;
/// let ints = IntSpan::from("1-3,5,7,9,100-999,1001-10000");
/// assert_eq!(ints.to_string(), "1-3,5,7,9,100-999,1001-10000");
/// assert_eq!(ints.cardinality(), 9906);
/// ```
///
/// # DESCRIPTION
///
/// `IntSpan` (ints for abbr.) represents sets of integers as a number of inclusive ranges, for example
/// `1-10,19-23,45-48`. Because many of its operations involve linear searches of the list of ranges its
/// overall performance tends to be proportional to the number of distinct ranges. This is fine for
/// small sets but suffers compared to other possible set representations (bit vectors, hash keys) when
/// the number of ranges grows large.
///
/// This module also represents sets as ranges of values but stores those ranges in order and uses a
/// binary search for many internal operations so that overall performance tends towards O log N where N
/// is the number of ranges.
///
/// The internal representation used by this module is extremely simple: a set is represented as a list
/// of integers. Integers in even numbered positions (0, 2, 4 etc) represent the start of a run of
/// numbers while those in odd numbered positions represent the ends of runs. As an example the set (1,
/// 3-7, 9, 11, 12) would be represented internally as (1, 2, 3, 8, 11, 13).
///
/// Sets may be infinite - assuming you're prepared to accept that infinity is actually no more than a
/// fairly large integer. Specifically the constants `neg_inf` and `pos_inf` are defined to be (-2^31+1)
/// and (2^31-2) respectively. To create an infinite set invert an empty one:
///
/// ```
/// # use pgr::libs::ds::IntSpan;
/// let mut ints = IntSpan::new();
/// ints.invert();
/// let expected = format!("{}-{}", ints.get_neg_inf(), ints.get_pos_inf());
/// assert_eq!(ints.to_string(), expected);
/// assert_eq!(ints.is_empty(), false);
/// assert_eq!(ints.is_universal(), true);
/// assert_eq!(ints.is_infinite(), true);
/// assert_eq!(ints.is_finite(), false);
/// assert_eq!(ints.is_pos_inf(), true);
/// assert_eq!(ints.is_neg_inf(), true);
/// ```
///
/// Sets need only be bounded in one direction - for example this is the set of all positive integers
/// (assuming you accept the slightly feeble definition of infinity we're using):
///
/// ```
/// # use pgr::libs::ds::IntSpan;
/// let mut ints = IntSpan::new();
/// ints.add_pair(1, ints.get_pos_inf());
/// let expected = format!("{}-{}", 1, ints.get_pos_inf());
/// assert_eq!(ints.to_string(), expected);
/// assert_eq!(ints.is_empty(), false);
/// assert_eq!(ints.is_universal(), false);
/// assert_eq!(ints.is_infinite(), true);
/// assert_eq!(ints.is_finite(), false);
/// assert_eq!(ints.is_pos_inf(), true);
/// assert_eq!(ints.is_neg_inf(), false);
/// ```
///
/// This Rust crate is ported from the Java class `jintspan` and the Perl module `AlignDB::IntSpan`,
/// which contains many codes from `Set::IntSpan`, `Set::IntSpan::Fast` and `Set::IntSpan::Island`.
///
#[derive(Debug, Default, Clone)]
pub struct IntSpan {
    edges: VecDeque<i32>,
}

const POS_INF: i32 = 2_147_483_647 - 1; // INT_MAX - 1, Real Largest int is POS_INF - 1
const NEG_INF: i32 = -2_147_483_648 + 1;
static EMPTY_STRING: std::sync::LazyLock<String> = std::sync::LazyLock::new(|| "-".to_string());

/// INTERFACE: Set creation and contents
///
/// ----
/// ----
impl IntSpan {
    #[inline]
    pub fn new() -> Self {
        IntSpan {
            edges: VecDeque::new(),
        }
    }

    pub fn from(runlist: &str) -> Self {
        Self::try_from(runlist).expect("invalid runlist string")
    }

    /// Parse `runlist` into a set, returning an error instead of panicking
    /// on invalid input.
    pub fn try_from(runlist: &str) -> anyhow::Result<Self> {
        let mut new = Self::new();
        if !runlist.is_empty() && runlist != *EMPTY_STRING {
            let ranges = new.runlist_to_ranges(runlist)?;
            new.add_ranges(&ranges);
        }
        Ok(new)
    }

    pub fn valid(runlist: &str) -> bool {
        let new = Self::new();
        new.runlist_to_ranges(runlist).is_ok()
    }

    pub fn from_pair(lower: i32, upper: i32) -> Self {
        let mut new = Self::new();
        new.add_pair(lower, upper);

        new
    }

    /// Build a set from inclusive `(lower, upper)` pairs, sorting and merging
    /// them in a single pass (O(n log n)). Pairs with `lower > upper` or an
    /// `upper` beyond the representable maximum (`POS_INF - 1`) are skipped.
    pub fn from_pairs(pairs: impl IntoIterator<Item = (i32, i32)>) -> Self {
        let mut pairs: Vec<(i32, i32)> = pairs
            .into_iter()
            .filter(|(l, u)| l <= u && *u < POS_INF)
            .collect();
        pairs.sort_unstable();
        let mut edges: VecDeque<i32> = VecDeque::new();
        let mut cur: Option<(i32, i32)> = None;
        for (lo, hi) in pairs {
            match cur {
                Some((clo, chi)) if lo <= chi + 1 => cur = Some((clo, chi.max(hi))),
                Some((clo, chi)) => {
                    edges.push_back(clo);
                    edges.push_back(chi + 1);
                    cur = Some((lo, hi));
                }
                None => cur = Some((lo, hi)),
            }
        }
        if let Some((lo, hi)) = cur {
            edges.push_back(lo);
            edges.push_back(hi + 1);
        }
        IntSpan { edges }
    }

    /// Returns the constant of POS_INF
    ///
    /// Typically used to construct infinite sets
    #[inline]
    pub fn get_pos_inf(&self) -> i32 {
        POS_INF - 1
    }

    /// Returns the constant of NEG_INF
    ///
    /// Typically used to construct infinite sets
    #[inline]
    pub fn get_neg_inf(&self) -> i32 {
        NEG_INF
    }

    /// Clears all contents of ints
    #[inline]
    pub fn clear(&mut self) {
        self.edges.clear();
    }

    #[inline]
    pub fn edge_size(&self) -> usize {
        self.edges.len()
    }

    #[inline]
    pub fn span_size(&self) -> usize {
        self.edge_size() / 2
    }

    pub fn to_vec(&self) -> Vec<i32> {
        self.spans()
            .into_iter()
            .flat_map(|(lower, upper)| (lower..=upper).collect::<Vec<_>>())
            .collect()
    }

    #[inline]
    pub fn contains(&self, n: i32) -> bool {
        let pos = self.find_pos(n + 1, 0);
        (pos & 1) == 1
    }

    /// Number of bases of the inclusive range `[start, end]` that are in the
    /// set. Binary-searches the sorted spans, so it is O(log n + k) in the
    /// number of spans overlapping the query.
    pub fn covered(&self, start: i32, end: i32) -> i32 {
        if self.is_empty() || end < start {
            return 0;
        }
        let n = self.span_size();
        // VecDeque keeps its buffer as one or two contiguous slices; reading
        // through them avoids per-access ring-offset arithmetic (the common
        // case is a single slice after append-built sets).
        let (s1, s2) = self.edges.as_slices();
        let n1 = s1.len();
        let edge_at = |i: usize| -> i32 {
            if i < n1 {
                s1[i]
            } else {
                s2[i - n1]
            }
        };
        // First span whose end >= start.
        let mut lo = 0usize;
        let mut hi = n;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if edge_at(2 * mid + 1) - 1 < start {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let first = lo;
        // One past the last span whose start <= end.
        let mut lo = 0usize;
        let mut hi = n;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if edge_at(2 * mid) <= end {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let last = lo;
        let mut total: i64 = 0;
        for i in first..last {
            let s = edge_at(2 * i);
            let e = edge_at(2 * i + 1) - 1;
            total += i64::from(end.min(e) - start.max(s) + 1);
        }
        total.min(i64::from(i32::MAX)) as i32
    }

    #[inline]
    pub fn min(&self) -> i32 {
        if self.is_empty() {
            panic!("Can't get extrema for empty IntSpan");
        }

        *self.edges.front().unwrap()
    }

    #[inline]
    pub fn max(&self) -> i32 {
        if self.is_empty() {
            panic!("Can't get extrema for empty IntSpan");
        }

        *self.edges.back().unwrap() - 1
    }
}

#[cfg(test)]
mod create {
    use super::*;

    #[test]
    fn test_create() {
        let tests = vec![
            ("", "-", vec![]),
            ("-", "-", vec![]),
            ("0", "0", vec![0]),
            ("1", "1", vec![1]),
            ("-1", "-1", vec![-1]),
            ("1-2", "1-2", vec![1, 2]),
            ("-2--1", "-2--1", vec![-2, -1]),
            ("-2-1", "-2-1", vec![-2, -1, 0, 1]),
            ("1,3-4", "1,3-4", vec![1, 3, 4]),
            ("1-1", "1", vec![1]),
            ("1,2-4", "1-4", vec![1, 2, 3, 4]),
            ("1-3,4", "1-4", vec![1, 2, 3, 4]),
            ("1-3,4,5-7", "1-7", vec![1, 2, 3, 4, 5, 6, 7]),
            ("1,2,3,4,5,6,7", "1-7", vec![1, 2, 3, 4, 5, 6, 7]),
        ];

        // create new
        for (runlist, exp_runlist, exp_elements) in &tests {
            let mut intspan = IntSpan::new();
            intspan.add_runlist(runlist);

            assert_eq!(intspan.cardinality(), exp_elements.len() as i32);
            assert_eq!(intspan.size(), exp_elements.len() as i32);
            assert_eq!(intspan.to_string(), *exp_runlist);
            assert_eq!(intspan.runlist(), *exp_runlist);
            assert_eq!(intspan.to_vec(), *exp_elements);
            assert_eq!(intspan.elements(), *exp_elements);
        }

        for (runlist, exp_runlist, exp_elements) in &tests {
            let intspan = IntSpan::from(runlist);

            assert_eq!(intspan.cardinality(), exp_elements.len() as i32);
            assert_eq!(intspan.to_string(), *exp_runlist);
            assert_eq!(intspan.to_vec(), *exp_elements);
        }

        for (_, exp_runlist, exp_elements) in &tests {
            let mut intspan = IntSpan::new();
            intspan.add_vec(exp_elements);

            assert_eq!(intspan.cardinality(), exp_elements.len() as i32);
            assert_eq!(intspan.to_string(), *exp_runlist);
            assert_eq!(intspan.to_vec(), *exp_elements);
        }
    }

    #[test]
    fn test_valid() {
        let tests = vec![
            ("", true),
            ("-", true),
            ("-2--1", true),
            ("1-3,4,5-7", true),
            ("abc", false),
            ("abc-def", false),
            ("abc,def", false),
            // Trailing dash, empty upper and reversed pairs must be invalid,
            // not panics (previously `1-` panicked on an out-of-bounds read
            // and `1-0` panicked in add_pair).
            ("1-", false),
            ("1-0", false),
            ("5-3", false),
            ("1--1", false),
            // Oversized digit runs are invalid instead of overflowing.
            ("99999999999", false),
            // Digit runs longer than the i64 accumulator used to panic in
            // `lower * radix`; they must be invalid instead.
            ("99999999999999999999", false),
            ("1-99999999999999999999", false),
            ("-99999999999999999999", false),
            ("2147483648", false),
            ("-2147483649", false),
            // Coordinates above POS_INF - 1 are unrepresentable (add_pair
            // would overflow when storing the upper + 1 edge).
            ("2147483647", false),
            ("-2147483648", true),
            // `add_pair` stores upper + 1 as an edge, so coordinates above
            // POS_INF - 1 (2147483645) are rejected instead of overflowing.
            ("1-2147483645", true),
            ("1-2147483646", false),
            ("2147483646", false),
        ];

        // create new
        for (runlist, exp) in &tests {
            assert_eq!(IntSpan::valid(runlist), *exp);
        }
    }

    #[test]
    fn try_from_parses_once_and_errors() {
        assert_eq!(IntSpan::try_from("1-3,5").unwrap().to_string(), "1-3,5");
        assert!(IntSpan::try_from("").unwrap().is_empty());
        assert!(IntSpan::try_from("-").unwrap().is_empty());
        assert!(IntSpan::try_from("abc").is_err());
        assert!(IntSpan::try_from("1-").is_err());
        assert!(IntSpan::try_from("5-3").is_err());
    }

    #[test]
    fn from_pairs_matches_add_pair_build() {
        let pairs = [
            (5, 5),
            (1, 3),
            (10, 12),
            (4, 5), // overlaps 1-3 (adjacent) and 5
            (20, 25),
            (3, 1), // reversed, skipped
        ];
        let mut expected = IntSpan::new();
        for (l, u) in pairs {
            if l <= u {
                expected.add_pair(l, u);
            }
        }
        assert_eq!(IntSpan::from_pairs(pairs).to_string(), expected.to_string());
        assert!(IntSpan::from_pairs([(1, 2147483647)]).is_empty());
    }

    #[test]
    #[should_panic(expected = "Bad order: 1,-1")]
    fn panic_pair() {
        let mut set = IntSpan::new();
        set.add_pair(1, -1);
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Bad order: 1,-1")]
    fn panic_runlist() {
        let mut set = IntSpan::new();
        set.add_runlist("1--1");
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Number format error: a at 0 of abc")]
    fn panic_runlist_2() {
        let mut set = IntSpan::new();
        set.add_runlist("abc");
        println!("{:?}", set.ranges());
    }

    // Read as 1-11
    //#[test]
    //#[should_panic(expected = "Bad order: 1,-1")]
    //fn panic_runlist_3() {
    //    let mut set = IntSpan::new();
    //    set.add_runlist("1-1--1");
    //    println!("{:?}", set.ranges());
    //}
}

/// INTERFACE: Span contents
///
/// ----
/// ----
impl IntSpan {
    /// Returns the runs in IntSpan, as a vector of Tuple(lower, upper)
    ///
    /// ```
    /// let ints = pgr::libs::ds::IntSpan::from("1-2,4-7");
    /// assert_eq!(ints.spans(), vec![(1, 2), (4, 7)]);
    /// ```
    pub fn spans(&self) -> Vec<(i32, i32)> {
        (0..self.span_size())
            .map(|i| {
                let lower = self.edges[i * 2];
                let upper = self.edges[i * 2 + 1] - 1;
                (lower, upper)
            })
            .collect()
    }

    /// Returns the runs in IntSpan, as a vector of lower, upper
    ///
    /// ```
    /// let ints = pgr::libs::ds::IntSpan::from("1-2,4-7");
    /// assert_eq!(ints.ranges(), vec![1, 2, 4, 7]);
    /// ```
    pub fn ranges(&self) -> Vec<i32> {
        self.spans()
            .into_iter()
            .flat_map(|(lower, upper)| vec![lower, upper])
            .collect()
    }

    /// Returns the runs in IntSpan, as a vector of String
    ///
    /// ```
    /// let ints = pgr::libs::ds::IntSpan::from("1-2,4-7");
    /// assert_eq!(ints.runs(), vec!["1-2".to_string(), "4-7".to_string()]);
    /// ```
    pub fn runs(&self) -> Vec<String> {
        self.spans()
            .into_iter()
            .map(|(lower, upper)| Self::from_pair(lower, upper).to_string())
            .collect()
    }

    /// Returns the runs in IntSpan, as a vector of IntSpan
    ///
    /// ```
    /// let ints = pgr::libs::ds::IntSpan::from("1-2,4-7");
    /// assert_eq!(ints.intses().iter().map(|e| e.to_string()).collect::<Vec<String>>(),
    ///     vec!["1-2".to_string(), "4-7".to_string()]);
    /// ```
    pub fn intses(&self) -> Vec<IntSpan> {
        self.spans()
            .into_iter()
            .map(|(lower, upper)| Self::from_pair(lower, upper))
            .collect()
    }
}

#[cfg(test)]
mod content {
    use super::*;

    #[test]
    fn test_content() {
        let tests = vec![
            ("-", "-", vec![]),
            ("0", "0", vec![(0, 0)]),
            ("1", "1", vec![(1, 1)]),
            ("-1", "-1", vec![(-1, -1)]),
            ("1-2", "1-2", vec![(1, 2)]),
            ("-2--1", "-2--1", vec![(-2, -1)]),
            ("-2-1", "-2-1", vec![(-2, 1)]),
            ("1,3-4", "1,3-4", vec![(1, 1), (3, 4)]),
            ("1-2,4-7", "1-2,4-7", vec![(1, 2), (4, 7)]),
        ];

        // spans
        for (runlist, _, exp_spans) in &tests {
            let mut ints = IntSpan::new();
            ints.add_runlist(runlist);

            let res = ints.spans();

            assert_eq!(res, *exp_spans);
        }
    }
}

/// INTERFACE: Set cardinality
///
/// ----
/// ----
impl IntSpan {
    pub fn cardinality(&self) -> i32 {
        if self.is_empty() {
            return 0;
        }

        // Spans can cover more than i32::MAX integers (e.g. a nearly full
        // i32 range); accumulate in i64 and saturate instead of panicking
        // (debug) or wrapping (release).
        let mut total: i64 = 0;
        for (lower, upper) in self.spans() {
            total += i64::from(upper) - i64::from(lower) + 1;
        }
        total.min(i64::from(i32::MAX)) as i32
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    #[inline]
    pub fn is_neg_inf(&self) -> bool {
        self.edges.front().is_some_and(|e| *e == NEG_INF)
    }

    #[inline]
    pub fn is_pos_inf(&self) -> bool {
        self.edges.back().is_some_and(|e| *e == POS_INF)
    }

    #[inline]
    pub fn is_infinite(&self) -> bool {
        self.is_neg_inf() || self.is_pos_inf()
    }

    #[inline]
    pub fn is_finite(&self) -> bool {
        !self.is_infinite()
    }

    #[inline]
    pub fn is_universal(&self) -> bool {
        self.edge_size() == 2 && self.is_pos_inf() && self.is_neg_inf()
    }
}

/// INTERFACE: Member operations (mutate original set)
///
/// ----
/// ----
impl IntSpan {
    pub fn add_pair(&mut self, mut lower: i32, mut upper: i32) {
        if lower > upper {
            panic!("Bad order: {},{}", lower, upper)
        }

        upper += 1;

        let mut lower_pos = self.find_pos(lower, 0);
        let mut upper_pos = self.find_pos(upper + 1, lower_pos);

        if lower_pos & 1 == 1 {
            lower_pos -= 1;
            lower = *self.edges.get(lower_pos).unwrap();
        }

        if upper_pos & 1 == 1 {
            upper = *self.edges.get(upper_pos).unwrap();
            upper_pos += 1;
        }

        for _i in lower_pos..upper_pos {
            self.edges.remove(lower_pos);
        }
        self.edges.insert(lower_pos, lower);
        self.edges.insert(lower_pos + 1, upper);
    }

    pub fn add_n(&mut self, n: i32) {
        self.add_pair(n, n);
    }

    pub fn add_ranges(&mut self, ranges: &[i32]) {
        if !ranges.len().is_multiple_of(2) {
            panic!("Number of ranges must be even")
        }

        for i in 0..(ranges.len() / 2) {
            let lower = *ranges.get(i * 2).unwrap();
            let upper = *ranges.get(i * 2 + 1).unwrap();

            self.add_pair(lower, upper);
        }

        // CAUTIONS: can't capture bad orders
        //        // When this IntSpan is empty, just convert ranges to edges
        //        if self.is_empty() {
        //            for i in 0..ranges.len() {
        //                // odd index means upper
        //                if (i & 1) == 1 {
        //                    self.edges.push(*ranges.get(i).unwrap() + 1);
        //                } else {
        //                    self.edges.push(*ranges.get(i).unwrap());
        //                }
        //            }
        //        } else {
        //            for i in 0..(ranges.len() / 2) {
        //                let lower = *ranges.get(i * 2).unwrap();
        //                let upper = *ranges.get(i * 2 + 1).unwrap();
        //
        //                self.add_pair(lower, upper);
        //            }
        //        }
    }

    pub fn merge(&mut self, other: &Self) {
        let ranges = other.ranges();

        self.add_ranges(&ranges);
    }

    pub fn add_vec(&mut self, ints: &[i32]) {
        let ranges = self.list_to_ranges(ints);

        self.add_ranges(&ranges);
    }

    // https://hermanradtke.com/2015/05/06/creating-a-rust-function-that-accepts-string-or-str.html
    pub fn add_runlist(&mut self, runlist: &str) {
        // skip empty runlist
        if !runlist.is_empty() && !runlist.eq(&*EMPTY_STRING) {
            let ranges = self.runlist_to_ranges(runlist).unwrap();
            self.add_ranges(&ranges);
        }
    }

    pub fn invert(&mut self) {
        if self.is_empty() {
            // Universal set
            self.edges.push_back(NEG_INF);
            self.edges.push_back(POS_INF);
        } else {
            // Either add or remove infinity from each end. The net effect is always an even number
            // of additions and deletions

            if self.is_neg_inf() {
                self.edges.pop_front(); // shift
            } else {
                self.edges.push_front(NEG_INF); // unshift
            }

            if self.is_pos_inf() {
                self.edges.pop_back(); // pop
            } else {
                self.edges.push_back(POS_INF); // push
            }
        }
    }

    pub fn remove_pair(&mut self, lower: i32, upper: i32) {
        self.invert();
        self.add_pair(lower, upper);
        self.invert();
    }

    pub fn remove_n(&mut self, n: i32) {
        self.remove_pair(n, n);
    }

    pub fn remove_ranges(&mut self, ranges: &[i32]) {
        if !ranges.len().is_multiple_of(2) {
            panic!("Number of ranges must be even");
        }

        self.invert();
        self.add_ranges(ranges);
        self.invert();
    }

    pub fn subtract(&mut self, other: &Self) {
        let ranges = other.ranges();

        self.remove_ranges(&ranges);
    }

    pub fn remove_vec(&mut self, array: &[i32]) {
        let ranges = self.list_to_ranges(array);

        self.remove_ranges(&ranges);
    }

    pub fn remove_runlist(&mut self, runlist: &str) {
        // skip empty runlist
        if !runlist.is_empty() && !runlist.eq(&*EMPTY_STRING) {
            let ranges = self.runlist_to_ranges(runlist).unwrap();
            self.remove_ranges(&ranges);
        }
    }
}

#[cfg(test)]
mod mutate {
    use super::*;

    #[test]
    fn test_mutate() {
        let sets = ["-", "1", "1-2", "1,3-5"];

        let contains = [
            [false, false, false, false],
            [true, false, false, false],
            [true, true, false, false],
            [true, false, true, true],
        ];

        let added = [
            ["1", "2", "3", "4"],
            ["1", "1-2", "1,3", "1,4"],
            ["1-2", "1-2", "1-3", "1-2,4"],
            ["1,3-5", "1-5", "1,3-5", "1,3-5"],
        ];

        let removed = [
            ["-", "-", "-", "-"],
            ["-", "1", "1", "1"],
            ["2", "1", "1-2", "1-2"],
            ["3-5", "1,3-5", "1,4-5", "1,3,5"],
        ];

        for i in 0..4 {
            for j in 0..4 {
                let n = j + 1;

                let set = IntSpan::from(sets[i]);
                let mut set_added = set.copy();
                set_added.add_n(n);

                let mut set_removed = set.copy();
                set_removed.remove_n(n);

                // contains
                assert_eq!(set.contains(n), contains[i][j as usize]);

                // added
                assert_eq!(set_added.to_string(), added[i][j as usize].to_string());

                // removed
                assert_eq!(set_removed.to_string(), removed[i][j as usize].to_string());
            }
        }
    }
}

/// INTERFACE: Set binary operations (create new set)
///
/// ----
/// ----
impl IntSpan {
    #[inline]
    pub fn copy(&self) -> Self {
        IntSpan {
            edges: self.edges.clone(),
        }
    }

    pub fn union(&self, other: &Self) -> Self {
        if self.is_empty() {
            return other.copy();
        }
        if other.is_empty() {
            return self.copy();
        }
        // Linear merge of the two sorted, disjoint span lists (O(n + m)),
        // coalescing overlapping or adjacent spans.
        let a = self.spans();
        let b = other.spans();
        let mut edges: VecDeque<i32> = VecDeque::new();
        let (mut i, mut j) = (0usize, 0usize);
        let mut cur: Option<(i32, i32)> = None;
        loop {
            let next = if i < a.len() && (j >= b.len() || a[i].0 <= b[j].0) {
                let s = a[i];
                i += 1;
                Some(s)
            } else if j < b.len() {
                let s = b[j];
                j += 1;
                Some(s)
            } else {
                None
            };
            let Some((lo, hi)) = next else { break };
            match cur {
                Some((clo, chi)) if lo <= chi + 1 => cur = Some((clo, chi.max(hi))),
                Some((clo, chi)) => {
                    edges.push_back(clo);
                    edges.push_back(chi + 1);
                    cur = Some((lo, hi));
                }
                None => cur = Some((lo, hi)),
            }
        }
        if let Some((lo, hi)) = cur {
            edges.push_back(lo);
            edges.push_back(hi + 1);
        }
        IntSpan { edges }
    }

    pub fn complement(&self) -> Self {
        let mut new = self.copy();
        new.invert();
        new
    }

    pub fn diff(&self, other: &Self) -> Self {
        if self.is_empty() {
            Self::new()
        } else if other.is_empty() {
            self.copy()
        } else {
            // Linear walk subtracting the sorted other spans from each self
            // span (O(n + m)); pieces are emitted in ascending order.
            let a = self.spans();
            let b = other.spans();
            let mut edges: VecDeque<i32> = VecDeque::new();
            let mut j = 0usize;
            for &(a_lo, a_hi) in &a {
                while j < b.len() && b[j].1 < a_lo {
                    j += 1;
                }
                let mut cur = a_lo;
                let mut k = j;
                while k < b.len() && b[k].0 <= a_hi {
                    let (b_lo, b_hi) = b[k];
                    if b_lo > cur {
                        edges.push_back(cur);
                        edges.push_back(b_lo);
                    }
                    cur = cur.max(b_hi + 1);
                    k += 1;
                }
                if cur <= a_hi {
                    edges.push_back(cur);
                    edges.push_back(a_hi + 1);
                }
            }
            IntSpan { edges }
        }
    }

    pub fn intersect(&self, other: &Self) -> Self {
        if self.is_empty() || other.is_empty() {
            Self::new()
        } else {
            // Linear two-pointer merge over the sorted, disjoint span lists
            // (O(n + m)). The previous complement+merge+invert chain called
            // `add_pair` on large sets, whose VecDeque shifts made intersect
            // O(n·m) in the worst case.
            let a = self.spans();
            let b = other.spans();
            let mut edges: VecDeque<i32> = VecDeque::new();
            let (mut i, mut j) = (0usize, 0usize);
            while i < a.len() && j < b.len() {
                let (a_lo, a_hi) = a[i];
                let (b_lo, b_hi) = b[j];
                let lo = a_lo.max(b_lo);
                let hi = a_hi.min(b_hi);
                if lo <= hi {
                    // Result spans arrive in ascending order; extend the
                    // last span when adjacent (defensive, normally impossible
                    // given the non-adjacent input span invariant).
                    if edges.back() == Some(&lo) {
                        *edges.back_mut().unwrap() = hi + 1;
                    } else {
                        edges.push_back(lo);
                        edges.push_back(hi + 1);
                    }
                }
                if a_hi < b_hi {
                    i += 1;
                } else {
                    j += 1;
                }
            }
            IntSpan { edges }
        }
    }

    pub fn xor(&self, other: &Self) -> Self {
        self.union(other).diff(&self.intersect(other))
    }
}

#[cfg(test)]
mod binary {
    use super::*;

    #[test]
    fn test_binary() {
        //   A    B    U    I    X    A-B  B-A
        let tests = vec![
            ("-", "-", "-", "-", "-", "-", "-"),
            ("1", "1", "1", "1", "-", "-", "-"),
            ("1", "2", "1-2", "-", "1-2", "1", "2"),
            ("3-9", "1-2", "1-9", "-", "1-9", "3-9", "1-2"),
            ("3-9", "1-5", "1-9", "3-5", "1-2,6-9", "6-9", "1-2"),
            ("3-9", "4-8", "3-9", "4-8", "3,9", "3,9", "-"),
            ("3-9", "5-12", "3-12", "5-9", "3-4,10-12", "3-4", "10-12"),
            ("3-9", "10-12", "3-12", "-", "3-12", "3-9", "10-12"),
            (
                "1-3,5,8-11",
                "1-6",
                "1-6,8-11",
                "1-3,5",
                "4,6,8-11",
                "8-11",
                "4,6",
            ),
        ];

        for (a, b, u, i, x, ab, ba) in tests {
            let ia = IntSpan::from(a);
            let ib = IntSpan::from(b);

            // union
            assert_eq!(ia.union(&ib).to_string(), u);

            // intersect
            assert_eq!(ia.intersect(&ib).to_string(), i);

            // xor
            assert_eq!(ia.xor(&ib).to_string(), x);

            // diff A-B
            assert_eq!(ia.diff(&ib).to_string(), ab);

            // diff B-A
            assert_eq!(ib.diff(&ia).to_string(), ba);
        }
    }

    #[test]
    fn covered_matches_intersect() {
        let set = IntSpan::from("1-10,20-30,50-100");
        for (start, end, expected) in [
            (1, 10, 10),   // full span
            (5, 25, 12),   // 5-10 + 20-25
            (11, 19, 0),   // gap
            (30, 35, 1),   // boundary
            (0, 5, 5),     // before the set
            (1, 100, 72),  // 10 + 11 + 51
            (120, 130, 0), // after the set
            (5, 5, 1),     // single point inside
            (15, 15, 0),   // single point in gap
        ] {
            let mut range = IntSpan::new();
            range.add_pair(start, end);
            assert_eq!(
                set.covered(start, end),
                set.intersect(&range).cardinality(),
                "covered({start},{end})"
            );
            assert_eq!(set.covered(start, end), expected);
        }
        assert_eq!(set.covered(5, 3), 0); // reversed
        assert_eq!(IntSpan::new().covered(1, 10), 0); // empty set
    }

    #[test]
    fn intersect_matches_slow_implementation() {
        // The old complement+merge+invert implementation, kept as an oracle.
        fn slow(a: &IntSpan, b: &IntSpan) -> IntSpan {
            if a.is_empty() || b.is_empty() {
                return IntSpan::new();
            }
            let mut new = a.complement();
            new.merge(&b.complement());
            new.invert();
            new
        }

        // Deterministic pseudo-random disjoint runlists over [-100, 100].
        let mut x = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) as i32
        };
        let mut runlists: Vec<String> = Vec::new();
        for _ in 0..200 {
            let mut runs: Vec<(i32, i32)> = Vec::new();
            let n = next().rem_euclid(8);
            let mut pos = -100i32;
            for _ in 0..n {
                pos += next().rem_euclid(20) + 1;
                let len = next().rem_euclid(15) + 1;
                runs.push((pos, pos + len - 1));
                pos += len;
            }
            runlists.push(if runs.is_empty() {
                "-".to_string()
            } else {
                runs.iter()
                    .map(|(l, u)| format!("{}-{}", l, u))
                    .collect::<Vec<_>>()
                    .join(",")
            });
        }
        for a in &runlists {
            for b in &runlists {
                let ia = IntSpan::from(a);
                let ib = IntSpan::from(b);
                assert_eq!(
                    ia.intersect(&ib).to_string(),
                    slow(&ia, &ib).to_string(),
                    "a={a} b={b}"
                );
            }
        }
    }

    #[test]
    fn union_diff_xor_match_slow_implementations() {
        fn union_slow(a: &IntSpan, b: &IntSpan) -> IntSpan {
            let mut new = a.copy();
            new.merge(b);
            new
        }
        fn diff_slow(a: &IntSpan, b: &IntSpan) -> IntSpan {
            if a.is_empty() {
                return IntSpan::new();
            }
            let mut new = a.copy();
            new.subtract(b);
            new
        }
        fn xor_slow(a: &IntSpan, b: &IntSpan) -> IntSpan {
            let mut new = union_slow(a, b);
            // The old `intersect` short-circuits on empty inputs.
            let isect = if a.is_empty() || b.is_empty() {
                IntSpan::new()
            } else {
                let mut isect = a.complement();
                isect.merge(&b.complement());
                isect.invert();
                isect
            };
            new.subtract(&isect);
            new
        }

        let runlists = ["-", "1", "1-3,7-10", "-2--1,5,100-200", "3-9", "5-12,20-25"];
        for a in runlists {
            for b in runlists {
                let ia = IntSpan::from(a);
                let ib = IntSpan::from(b);
                assert_eq!(
                    ia.union(&ib).to_string(),
                    union_slow(&ia, &ib).to_string(),
                    "union a={a} b={b}"
                );
                assert_eq!(
                    ia.diff(&ib).to_string(),
                    diff_slow(&ia, &ib).to_string(),
                    "diff a={a} b={b}"
                );
                assert_eq!(
                    ia.xor(&ib).to_string(),
                    xor_slow(&ia, &ib).to_string(),
                    "xor a={a} b={b}"
                );
            }
        }
    }

    #[test]
    fn distance_and_islands_match_slow_implementations() {
        fn distance_slow(a: &IntSpan, b: &IntSpan) -> i32 {
            if a.is_empty() || b.is_empty() {
                return 0;
            }
            let overlap = a.overlap(b);
            if overlap > 0 {
                return -overlap;
            }
            let mut min_d = 0;
            for (lower1, upper1) in a.spans() {
                for (lower2, upper2) in b.spans() {
                    let d1 = (lower1 - upper2).abs();
                    let d2 = (upper1 - lower2).abs();
                    let d = d1.min(d2);
                    if min_d == 0 || d < min_d {
                        min_d = d;
                    }
                }
            }
            min_d
        }
        fn islands_slow(a: &IntSpan, b: &IntSpan) -> IntSpan {
            let mut island = IntSpan::new();
            if !a.intersect(b).is_empty() {
                for (lower, upper) in a.spans() {
                    let subints = IntSpan::from_pair(lower, upper);
                    if !subints.intersect(b).is_empty() {
                        island.merge(&subints);
                    }
                }
            }
            island
        }

        let runlists = [
            "-",
            "1",
            "1-3,7-10",
            "-2--1,5,100-200",
            "3-9",
            "5-12,20-25",
            "1-5,6-10",
        ];
        for a in runlists {
            for b in runlists {
                let ia = IntSpan::from(a);
                let ib = IntSpan::from(b);
                assert_eq!(
                    ia.distance(&ib),
                    distance_slow(&ia, &ib),
                    "distance a={a} b={b}"
                );
                assert_eq!(
                    ia.find_islands_ints(&ib).to_string(),
                    islands_slow(&ia, &ib).to_string(),
                    "islands a={a} b={b}"
                );
            }
        }
    }
}

/// INTERFACE: Set relations
///
/// ----
/// ----
impl IntSpan {
    pub fn equals(&self, other: &Self) -> bool {
        let edges = &self.edges;
        let edges_other = &other.edges;

        if edges.len() != edges_other.len() {
            return false;
        }

        for i in 0..edges.len() {
            if edges.get(i) != edges_other.get(i) {
                return false;
            }
        }

        true
    }

    pub fn subset(&self, other: &Self) -> bool {
        self.diff(other).is_empty()
    }

    pub fn superset(&self, other: &Self) -> bool {
        other.diff(self).is_empty()
    }
}

#[cfg(test)]
mod relation {
    use super::*;

    #[test]
    fn test_relation() {
        let sets = ["-", "1", "5", "1-5", "3-7", "1-3,8,10-23"];

        let equals = [
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ];

        let subset = [
            [1, 1, 1, 1, 1, 1],
            [0, 1, 0, 1, 0, 1],
            [0, 0, 1, 1, 1, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ];

        let superset = [
            [1, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0],
            [1, 1, 1, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [1, 1, 0, 0, 0, 1],
        ];

        for i in 0..6 {
            for j in 0..6 {
                let a = IntSpan::from(sets[i]);
                let b = IntSpan::from(sets[j]);

                // equals
                assert_eq!(a.equals(&b), equals[i][j] != 0);

                // subset
                assert_eq!(a.subset(&b), subset[i][j] != 0);

                // superset
                assert_eq!(a.superset(&b), superset[i][j] != 0);
            }
        }
    }
}

/// INTERFACE: Indexing
///
/// ----
/// ----
impl IntSpan {
    fn at_pos(&self, index: i32) -> i32 {
        let mut element = self.min();
        let mut ele_before = 0;

        for i in 0..self.span_size() {
            let lower = *self.edges.get(i * 2).unwrap();
            let upper = *self.edges.get(i * 2 + 1).unwrap() - 1;

            let span_len = upper - lower + 1;

            if index > ele_before + span_len {
                ele_before += span_len;
            } else {
                element = index - ele_before - 1 + lower;
                break;
            }
        }

        element
    }

    fn at_neg(&self, index: i32) -> i32 {
        let mut element = self.max();
        let mut ele_after = 0;

        for i in (0..self.span_size()).rev() {
            let lower = *self.edges.get(i * 2).unwrap();
            let upper = *self.edges.get(i * 2 + 1).unwrap() - 1;

            let span_len = upper - lower + 1;

            if index > ele_after + span_len {
                ele_after += span_len;
            } else {
                element = upper - (index - ele_after) + 1;
                break;
            }
        }

        element
    }

    /// Returns the index-th element of set, indices start from `1`.
    ///
    /// Negative indices count backwards from the end of the set.
    pub fn at(&self, index: i32) -> i32 {
        if self.is_empty() {
            panic!("Indexing on an empty set");
        }
        if i32::abs(index) < 1 {
            panic!("Index can't be 0");
        }
        if i32::abs(index) > self.cardinality() {
            panic!("Out of max index");
        }

        if index > 0 {
            self.at_pos(index)
        } else {
            self.at_neg(-index)
        }
    }

    /// Returns the index of an element in the set, indices start from `1`
    pub fn index(&self, element: i32) -> i32 {
        if self.is_empty() {
            panic!("Indexing on an empty set");
        }
        if !self.contains(element) {
            panic!("Element doesn't exist");
        }

        let mut index = -1; // not valid
        let mut ele_before = 0;

        for i in 0..self.span_size() {
            let lower = *self.edges.get(i * 2).unwrap();
            let upper = *self.edges.get(i * 2 + 1).unwrap() - 1;
            let span_len = upper - lower + 1;

            if element >= lower && element <= upper {
                index = element - lower + 1 + ele_before;
            } else {
                ele_before += span_len;
            }
        }

        index
    }

    pub fn slice(&self, from: i32, to: i32) -> IntSpan {
        if self.is_empty() {
            panic!("Indexing on an empty set");
        }
        if from < 1 {
            panic!("Index can't be 0 or negative");
        }
        if to > self.cardinality() {
            panic!("Out of max index");
        }
        if from > to {
            panic!("Bad order: {},{}", from, to)
        }

        let lower = self.at(from);
        let upper = self.at(to);

        let new = IntSpan::from_pair(lower, upper);
        new.intersect(self)
    }
}

#[cfg(test)]
mod index {
    use super::*;

    #[test]
    fn test_index() {
        // runlist, n, exp_index, exp_element
        let tests = vec![
            // None
            ("-", 1, None, None),
            ("-", -1, None, None),
            ("1-10,21-30", 25, None, Some(15)),
            ("1-10,21-30", -25, None, None),
            // at_pos
            ("0-9", 1, Some(0), Some(2)),
            ("0-9", 6, Some(5), Some(7)),
            ("0-9", 10, Some(9), None),
            ("0-9", 11, None, None),
            // at_neg
            ("0-9", -1, Some(9), None),
            ("0-9", -5, Some(5), None),
            ("0-9", -10, Some(0), None),
            ("0-9", -11, None, None),
            // at_pos
            ("1-10,21-30,41-50", 6, Some(6), Some(6)),
            ("1-10,21-30,41-50", 16, Some(26), None),
            ("1-10,21-30,41-50", 26, Some(46), Some(16)),
            ("1-10,21-30,41-50", 31, None, None),
            // at_neg
            ("1-10,21-30,41-50", -1, Some(50), None),
            ("1-10,21-30,41-50", -11, Some(30), None),
            ("1-10,21-30,41-50", -21, Some(10), None),
            ("1-10,21-30,41-50", -30, Some(1), None),
            ("1-10,21-30,41-50", -31, None, None),
        ];

        for (runlist, n, exp_index, exp_element) in tests {
            let set = IntSpan::from(runlist);

            // at
            if let Some(exp_index) = exp_index {
                assert_eq!(set.at(n), exp_index);
            }

            // index
            if let Some(exp_element) = exp_element {
                assert_eq!(set.index(n), exp_element);
            }
        }
    }

    #[test]
    fn test_slice() {
        // runlist, from, to, exp
        let tests = vec![
            ("1-10,21-30,41-50", 1, 3, "1-3"),
            ("1-10,21-30,41-50", 6, 8, "6-8"),
            ("1-10,21-30,41-50", 8, 10, "8-10"),
            ("1-10,21-30,41-50", 10, 10, "10"),
        ];

        for (runlist, from, to, exp) in tests {
            let set = IntSpan::from(runlist);

            assert_eq!(set.slice(from, to).to_string(), exp);
        }
    }

    #[test]
    #[should_panic(expected = "Indexing on an empty set")]
    fn panic_at_1() {
        let set = IntSpan::new();
        set.at(1);
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Index can't be 0")]
    fn panic_at_2() {
        let set = IntSpan::from("0-9");
        set.at(0);
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Out of max index")]
    fn panic_at_3() {
        let set = IntSpan::from("0-9");
        set.at(15);
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Indexing on an empty set")]
    fn panic_index_1() {
        let set = IntSpan::new();
        set.index(1);
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Element doesn't exist")]
    fn panic_index_2() {
        let set = IntSpan::from("0-9");
        set.index(15);
        println!("{:?}", set.ranges());
    }

    #[test]
    #[should_panic(expected = "Indexing on an empty set")]
    fn panic_slice_1() {
        let set = IntSpan::new();
        set.slice(1, 2);
        println!("{:?}", set.ranges());
    }
}

/// INTERFACE: Spans Ops
///
/// ----
/// ----
impl IntSpan {
    pub fn cover(&self) -> Self {
        let mut new = IntSpan::new();
        if !self.is_empty() {
            new.add_pair(self.min(), self.max());
        }
        new
    }

    pub fn holes(&self) -> Self {
        let mut new = IntSpan::new();
        if self.is_empty() || self.is_universal() {
            // empty and universal set have no holes
            return new;
        }
        let complement = self.complement();
        let mut ranges = complement.ranges();

        // Remove infinite arms of complement set
        if complement.is_neg_inf() {
            ranges.remove(0);
            ranges.remove(0);
        }
        if complement.is_pos_inf() {
            ranges.pop();
            ranges.pop();
        }

        new.add_ranges(&ranges);

        new
    }

    pub fn inset(&self, n: i32) -> Self {
        let mut new = IntSpan::new();

        for i in 0..self.span_size() {
            let mut lower = *self.edges.get(i * 2).unwrap();
            let mut upper = *self.edges.get(i * 2 + 1).unwrap() - 1;

            if lower != self.get_neg_inf() {
                lower = lower.saturating_add(n).clamp(NEG_INF, POS_INF - 1);
            }
            if upper != self.get_pos_inf() {
                upper = upper.saturating_sub(n).clamp(NEG_INF, POS_INF - 1);
            }

            if lower <= upper {
                new.add_pair(lower, upper);
            }
        }

        new
    }

    pub fn trim(&self, n: i32) -> Self {
        self.inset(n)
    }

    pub fn pad(&self, n: i32) -> Self {
        self.inset(n.saturating_neg())
    }

    pub fn excise(&self, min_len: i32) -> Self {
        let mut new = IntSpan::new();

        for i in 0..self.span_size() {
            let lower = *self.edges.get(i * 2).unwrap();
            let upper = *self.edges.get(i * 2 + 1).unwrap() - 1;

            let span_len = i64::from(upper) - i64::from(lower) + 1;
            if span_len >= i64::from(min_len) {
                new.add_pair(lower, upper);
            }
        }

        new
    }

    pub fn fill(&self, max_len: i32) -> Self {
        let mut new = self.copy();
        let holes = self.holes();

        for i in 0..holes.span_size() {
            let lower = *holes.edges.get(i * 2).unwrap();
            let upper = *holes.edges.get(i * 2 + 1).unwrap() - 1;

            let span_len = i64::from(upper) - i64::from(lower) + 1;
            if span_len <= i64::from(max_len) {
                new.add_pair(lower, upper);
            }
        }

        new
    }

    /// Removes elements inside the range, and all elements greater than this range are shifted
    /// towards the negative direction
    pub fn banish(&self, start: i32, end: i32) -> Self {
        let mut new = IntSpan::new();
        if start > end {
            return self.copy(); // nothing to banish
        }
        // i64 keeps the length and the shifted coordinates from overflowing
        // i32 for extreme arguments.
        let remove_len = i64::from(end) - i64::from(start) + 1;

        // No elements in the tmp ints intersect with the range
        let ints = self.diff(&IntSpan::from_pair(start, end));
        for (lower, upper) in ints.spans().iter().rev() {
            if *upper < start {
                new.add_pair(*lower, *upper);
            } else if *lower > end {
                let shifted_lower = (i64::from(*lower) - remove_len)
                    .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                    as i32;
                let shifted_upper = (i64::from(*upper) - remove_len)
                    .clamp(i64::from(i32::MIN), i64::from(i32::MAX))
                    as i32;
                if shifted_lower <= shifted_upper {
                    new.add_pair(shifted_lower, shifted_upper);
                }
            } else {
                panic!("Something went wrong while banishing {}-{}", start, end);
            }
        }

        new
    }
}

#[cfg(test)]
mod span {
    use super::*;

    #[test]
    fn cover_holes() {
        // runlist expCover expHoles
        let tests = vec![
            ("-", "-", "-"),
            ("1", "1", "-"),
            ("5", "5", "-"),
            ("1,3,5", "1-5", "2,4"),
            ("1,3-5", "1-5", "2"),
            ("1-3,5,8-11", "1-11", "4,6-7"),
        ];

        for (runlist, exp_cover, exp_holes) in tests {
            let set = IntSpan::from(runlist);

            // cover
            assert_eq!(set.cover().to_string(), exp_cover);

            // holes
            assert_eq!(set.holes().to_string(), exp_holes);
        }
    }

    #[test]
    fn inset() {
        let neg = IntSpan::new().get_neg_inf();
        let pos = IntSpan::new().get_pos_inf();

        let uni = format!("{}-{}", neg, pos);

        // runlist n expected
        let tests = vec![
            ("-".to_string(), -2, "-".to_string()),
            ("-".to_string(), -1, "-".to_string()),
            ("-".to_string(), 0, "-".to_string()),
            ("-".to_string(), 1, "-".to_string()),
            ("-".to_string(), 2, "-".to_string()),
            (uni.clone(), -2, uni.clone()),
            (uni.clone(), 2, uni.clone()),
            (format!("{}-0", neg), -2, format!("{}-2", neg)),
            (format!("{}-0", neg), 2, format!("{}--2", neg)),
            (format!("0-{}", pos), -2, format!("-2-{}", pos)),
            (format!("0-{}", pos), 2, format!("2-{}", pos)),
            (
                "0,2-3,6-8,12-15,20-24,30-35".to_string(),
                -2,
                "-2-26,28-37".to_string(),
            ),
            (
                "0,2-3,6-8,12-15,20-24,30-35".to_string(),
                -1,
                "-1-9,11-16,19-25,29-36".to_string(),
            ),
            (
                "0,2-3,6-8,12-15,20-24,30-35".to_string(),
                0,
                "0,2-3,6-8,12-15,20-24,30-35".to_string(),
            ),
            (
                "0,2-3,6-8,12-15,20-24,30-35".to_string(),
                1,
                "7,13-14,21-23,31-34".to_string(),
            ),
            (
                "0,2-3,6-8,12-15,20-24,30-35".to_string(),
                2,
                "22,32-33".to_string(),
            ),
        ];

        // inset
        for (runlist, n, expected) in tests {
            let set = IntSpan::from(&runlist);
            assert_eq!(set.inset(n).to_string(), expected);
        }

        // trim and pad
        assert_eq!(IntSpan::from("1-3").pad(1).cardinality(), 5);
        assert_eq!(IntSpan::from("1-3").pad(2).cardinality(), 7);
        assert_eq!(IntSpan::from("1-3").trim(1).cardinality(), 1);
        assert_eq!(IntSpan::from("1-3").trim(2).cardinality(), 0);
    }

    #[test]
    fn excise_fill() {
        // runlist n expExcise expFill
        let tests = vec![
            ("1-5", 1, "1-5", "1-5"),
            ("1-5,7", 1, "1-5,7", "1-7"),
            ("1-5,7", 2, "1-5", "1-7"),
            ("1-5,7-8", 1, "1-5,7-8", "1-8"),
            ("1-5,7-8", 3, "1-5", "1-8"),
            ("1-5,7-8", 6, "-", "1-8"),
            ("1-5,7,9-10", 0, "1-5,7,9-10", "1-5,7,9-10"),
            ("1-5,9-10", 2, "1-5,9-10", "1-5,9-10"),
            ("1-5,9-10", 3, "1-5", "1-10"),
            ("1-5,9-10,12-13,15", 2, "1-5,9-10,12-13", "1-5,9-15"),
            ("1-5,9-10,12-13,15", 3, "1-5", "1-15"),
        ];

        for (runlist, n, exp_excise, exp_fill) in tests {
            let set = IntSpan::from(runlist);

            // excise
            assert_eq!(set.excise(n).to_string(), exp_excise);

            // fill
            assert_eq!(set.fill(n).to_string(), exp_fill);
        }
    }

    #[test]
    fn extreme_ops_do_not_overflow() {
        // trim/pad with huge n used to overflow i32 arithmetic; the result
        // is clamped to the representable coordinate range instead.
        assert_eq!(
            IntSpan::from("1-2").pad(2_147_483_647).to_string(),
            "-2147483646-2147483645"
        );
        assert_eq!(IntSpan::from("1-2").trim(2_147_483_647).to_string(), "-");

        // excise/fill/cardinality on a nearly full i32 range used to overflow
        // `upper - lower + 1`; they must not panic.
        let huge = IntSpan::from("-2147483647-2147483645");
        assert_eq!(huge.excise(1).to_string(), "-2147483647-2147483645");
        assert_eq!(huge.cardinality(), i32::MAX);
        assert_eq!(huge.fill(100).to_string(), "-2147483647-2147483645");
    }

    #[test]
    fn banish() {
        // runlist n expExcise expFill
        let tests = vec![
            ("-", 3, 3, "-"),
            ("1", 3, 3, "1"),
            ("5", 3, 3, "4"),
            ("1,3,5", 3, 3, "1,4"),
            ("1,3-5", 3, 3, "1,3-4"),
            ("1-3,5,8-11", 3, 3, "1-2,4,7-10"),
            ("1-3,5,8-11", 3, 5, "1-2,5-8"),
            ("1-3,5,8-11", -5, -3, "-2-0,2,5-8"),
        ];

        for (runlist, start, end, expected) in tests {
            let ints = IntSpan::from(runlist);

            assert_eq!(ints.banish(start, end).to_string(), expected);
        }
    }

    #[test]
    fn banish_extreme_args_do_not_overflow() {
        // start > end: nothing to banish, the set is returned unchanged.
        let set = IntSpan::from("1-10");
        assert_eq!(set.banish(5, 3).to_string(), "1-10");
        // Extreme span length used to overflow `end - start + 1`.
        let full = IntSpan::from("-2147483647-2147483645");
        assert_eq!(full.banish(-2147483647, 2147483645).to_string(), "-");
        // Shifting elements above a huge range clamps to the i32 range.
        let shifted = IntSpan::from("100-200").banish(-2147483647, 2147483645);
        assert_eq!(shifted.to_string(), "-");
        let partial = IntSpan::from("1-10,100-110").banish(1, 10);
        assert_eq!(partial.to_string(), "90-100");
    }
}

/// INTERFACE: Inter-set OPs
///
/// ----
/// ----
impl IntSpan {
    /// Returns the size of intersection of two sets.
    ///
    /// `set.overlap(&other)` equivalent to `set.intersect(&other).cardinality()`
    ///
    /// ```
    /// # use pgr::libs::ds::IntSpan;
    /// let set = IntSpan::from("1");
    /// let other = IntSpan::from("1");
    /// assert_eq!(set.overlap(&other), 1);
    /// let other = IntSpan::from("2");
    /// assert_eq!(set.overlap(&other), 0);
    /// let set = IntSpan::from("1-5");
    /// let other = IntSpan::from("1-10");
    /// assert_eq!(set.overlap(&other), 5);
    /// let set = IntSpan::from("1-5,6");
    /// let other = IntSpan::from("6-10");
    /// assert_eq!(set.overlap(&other), 1);
    /// ```
    pub fn overlap(&self, other: &Self) -> i32 {
        self.intersect(other).cardinality()
    }

    /// Returns the distance between sets, measured as follows.
    ///
    /// * If the sets overlap, then the distance is negative and given by `- set.overlap(&other)`
    ///
    /// * If the sets do not overlap, $d is positive and given by the distance on the integer line
    ///   between the two closest islands of the sets.
    ///
    /// ```
    /// # use pgr::libs::ds::IntSpan;
    /// let set = IntSpan::from("1");
    /// let other = IntSpan::from("1");
    /// assert_eq!(set.distance(&other), -1);
    /// let other = IntSpan::from("");
    /// assert_eq!(set.distance(&other), 0);
    /// let other = IntSpan::from("2");
    /// assert_eq!(set.distance(&other), 1);
    ///
    /// let set = IntSpan::from("1-5");
    /// let other = IntSpan::from("1-10");
    /// assert_eq!(set.distance(&other), -5);
    /// let other = IntSpan::from("10-15");
    /// assert_eq!(set.distance(&other), 5);
    /// let set = IntSpan::from("1-5,6");
    /// let other = IntSpan::from("6-10");
    /// assert_eq!(set.distance(&other), -1);
    ///
    /// let set = IntSpan::from("1-5,10-15");
    /// let other = IntSpan::from("5-9");
    /// assert_eq!(set.distance(&other), -1);
    /// let other = IntSpan::from("6");
    /// assert_eq!(set.distance(&other), 1);
    /// let other = IntSpan::from("7");
    /// assert_eq!(set.distance(&other), 2);
    /// let other = IntSpan::from("7-9");
    /// assert_eq!(set.distance(&other), 1);
    /// let other = IntSpan::from("16-20");
    /// assert_eq!(set.distance(&other), 1);
    /// let other = IntSpan::from("17-20");
    /// assert_eq!(set.distance(&other), 2);
    /// ```
    pub fn distance(&self, other: &Self) -> i32 {
        if self.is_empty() || other.is_empty() {
            0
        } else {
            let overlap = self.overlap(other);

            if overlap > 0 {
                -overlap
            } else {
                // Two-pointer walk over the sorted span lists (O(n + m));
                // with no overlap the minimum distance is the closest gap
                // between the boundaries of adjacent spans.
                let a = self.spans();
                let b = other.spans();
                let mut min_gap = i64::from(i32::MAX);
                let (mut i, mut j) = (0usize, 0usize);
                while i < a.len() && j < b.len() {
                    let (a_lo, a_hi) = a[i];
                    let (b_lo, b_hi) = b[j];
                    if a_hi < b_lo {
                        min_gap = min_gap.min(i64::from(b_lo) - i64::from(a_hi));
                        i += 1;
                    } else if b_hi < a_lo {
                        min_gap = min_gap.min(i64::from(a_lo) - i64::from(b_hi));
                        j += 1;
                    } else {
                        break; // unreachable: overlap is handled above
                    }
                }
                min_gap.min(i64::from(i32::MAX)) as i32
            }
        }
    }
}

/// INTERFACE: Islands
///
/// ----
/// ----
impl IntSpan {
    /// Returns an ints equals to the island containing the integer
    ///
    /// If the integer is not in the ints, an empty ints is returned
    ///
    /// ```
    /// # use pgr::libs::ds::IntSpan;
    /// let tests = vec![
    ///     ("1-5", 1, "1-5"),
    ///     ("1-5,7", 1, "1-5"),
    ///     ("1-5,7", 6, "-"),
    ///     ("1-5,7", 7, "7"),
    ///     ("1-5,7-8", 7, "7-8"),
    ///     ("1-5,8", 7, "-"),
    /// ];
    ///
    /// for (runlist, val, exp_ints) in &tests {
    ///     let ints = IntSpan::from(runlist);
    ///
    ///     let res = ints.find_islands_n(*val);
    ///
    ///     assert_eq!(res.to_string().as_str(), *exp_ints);
    /// }
    /// ```
    pub fn find_islands_n(&self, val: i32) -> IntSpan {
        let mut island = Self::new();

        // if pos & 1, i.e. pos is an odd number, val is in the ints
        // same as contains()
        let pos = self.find_pos(val + 1, 0);
        if (pos & 1) == 1 {
            let ranges = self.ranges();
            island.add_pair(ranges[pos - 1], ranges[pos]);
        }

        island
    }

    /// Returns an ints containing all islands intersecting `other`
    ///
    /// If `ints` and `other` don't intersect, an empty ints is returned
    ///
    /// ```
    /// # use pgr::libs::ds::IntSpan;
    /// let tests = vec![
    ///     ("1-8", "7-8", "1-8"),
    ///     ("1-5,7-8", "7-8", "7-8"),
    ///     ("1-5,8-9", "7-8", "8-9"),
    ///     ("1-5,8-9,11-15", "9-11", "8-9,11-15"),
    /// ];
    ///
    /// for (runlist, other, exp_ints) in &tests {
    ///     let ints = IntSpan::from(runlist);
    ///     let other = IntSpan::from(other);
    ///
    ///     let res = ints.find_islands_ints(&other);
    ///
    ///     assert_eq!(res.to_string().as_str(), *exp_ints);
    /// }
    /// ```
    pub fn find_islands_ints(&self, other: &Self) -> IntSpan {
        let mut island = Self::new();

        if !self.is_empty() && !other.is_empty() {
            // Two-pointer walk (O(n + m)): every self span overlapping any
            // other span is an island, emitted whole.
            let a = self.spans();
            let b = other.spans();
            let mut edges: VecDeque<i32> = VecDeque::new();
            let (mut i, mut j) = (0usize, 0usize);
            while i < a.len() && j < b.len() {
                let (a_lo, a_hi) = a[i];
                let (b_lo, b_hi) = b[j];
                if a_hi < b_lo {
                    i += 1;
                } else if b_hi < a_lo {
                    j += 1;
                } else {
                    edges.push_back(a_lo);
                    edges.push_back(a_hi + 1);
                    i += 1;
                }
            }
            island.edges = edges;
        }

        island
    }
}

/// INTERFACE: Aliases
///
/// ----
/// ----
impl IntSpan {
    #[inline]
    pub fn size(&self) -> i32 {
        self.cardinality()
    }

    #[inline]
    pub fn runlist(&self) -> String {
        self.to_string()
    }

    #[inline]
    pub fn elements(&self) -> Vec<i32> {
        self.to_vec()
    }
}

/// Private methods
///
/// ----
/// ----
impl IntSpan {
    #[inline]
    fn find_pos(&self, val: i32, mut low: usize) -> usize {
        let mut high = self.edge_size();

        while low < high {
            let mid = (low + high) / 2;
            let mid_edge = self.edges.get(mid).unwrap();
            match val.cmp(mid_edge) {
                Ordering::Less => high = mid,
                Ordering::Greater => low = mid + 1,
                Ordering::Equal => return mid,
            }
        }

        low
    }

    fn list_to_ranges(&self, array: &[i32]) -> Vec<i32> {
        let mut ranges: Vec<i32> = Vec::new();

        let mut vec = array.to_owned();
        vec.sort_unstable();
        vec.dedup();

        let len = vec.len();
        let mut pos: usize = 0;

        while pos < len {
            let mut end = pos + 1;
            while (end < len) && (vec[end] <= vec[end - 1] + 1) {
                end += 1;
            }
            ranges.push(vec[pos]);
            ranges.push(vec[end - 1]);
            pos = end;
        }

        ranges
    }

    fn runlist_to_ranges(&self, runlist: &str) -> anyhow::Result<Vec<i32>> {
        let mut ranges: Vec<i32> = Vec::new();

        let bytes = runlist.as_bytes();

        let radix = 10i64;
        let mut idx = 0; // index in runlist
        let len = bytes.len();

        let mut lower_is_neg = false;
        let mut upper_is_neg = false;
        let mut in_upper = false;

        while idx < len {
            let mut i = 0; // index in one run
            if *bytes.get(idx).unwrap() == b'-' {
                lower_is_neg = true;
                i += 1;
            }

            // Ported from Java Integer.parseInt(), accumulated negative so
            // that i32::MIN parses without overflow. i64 accumulation makes
            // oversized digit runs an error instead of a panic (debug) or a
            // silent wrap (release).
            let mut lower: i64 = 0;
            let mut upper: i64 = 0;

            while idx + i < len {
                let ch = bytes[idx + i];
                if ch.is_ascii_digit() {
                    if !in_upper {
                        lower = lower * radix - i64::from(ch - b'0');
                        // Digit runs longer than the i32 range would overflow
                        // the i64 accumulator; the accumulation is monotonic
                        // below i32::MIN, so an early exit is safe.
                        if lower < i64::from(i32::MIN) {
                            return Err(anyhow!(
                                "Number format error: out of range at {} of {}",
                                idx + i,
                                runlist
                            ));
                        }
                    } else {
                        upper = upper * radix - i64::from(ch - b'0');
                        if upper < i64::from(i32::MIN) {
                            return Err(anyhow!(
                                "Number format error: out of range at {} of {}",
                                idx + i,
                                runlist
                            ));
                        }
                    }
                } else if ch == b'-' {
                    if !in_upper {
                        in_upper = true;
                        if idx + i + 1 < len && bytes[idx + i + 1] == b'-' {
                            upper_is_neg = true;
                        }
                    }
                } else if ch == b',' {
                    i += 1;
                    break; // end of run
                } else {
                    return Err(anyhow!(
                        "Number format error: {} at {} of {}",
                        ch as char,
                        idx + i,
                        runlist
                    ));
                }

                i += 1;
            }

            let lower_val = if lower_is_neg { lower } else { -lower };
            let upper_val = if in_upper {
                if upper_is_neg {
                    upper
                } else {
                    -upper
                }
            } else {
                lower_val
            };
            let i32_range = i64::from(i32::MIN)..=i64::from(i32::MAX);
            if !i32_range.contains(&lower_val) || !i32_range.contains(&upper_val) {
                return Err(anyhow!(
                    "Number format error: out of range at {} of {}",
                    idx,
                    runlist
                ));
            }
            // `add_pair` stores upper + 1 as an edge, so the largest
            // representable coordinate is POS_INF - 1; larger values used
            // to panic in `find_pos(upper + 1)`.
            if upper_val > i64::from(POS_INF - 1) {
                return Err(anyhow!(
                    "Number format error: out of range at {} of {}",
                    idx,
                    runlist
                ));
            }
            if lower_val > upper_val {
                return Err(anyhow!("Bad order: {},{}", lower_val, upper_val));
            }
            ranges.push(lower_val as i32);
            ranges.push(upper_val as i32);

            // reset boolean flags
            lower_is_neg = false;
            upper_is_neg = false;
            in_upper = false;

            // start next run
            idx += i;
        }

        Ok(ranges)
    }
}

impl std::fmt::Display for IntSpan {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_empty() {
            return write!(f, "{}", *EMPTY_STRING);
        }

        let runlist = self
            .spans()
            .into_iter()
            .map(|(lower, upper)| {
                if lower == upper {
                    Cow::from(lower.to_string())
                } else {
                    Cow::from(format!("{}-{}", lower, upper))
                }
            })
            .collect::<Vec<_>>()
            .join(",");

        write!(f, "{}", runlist)
    }
}

// ── Runlist JSON helpers (migrated from the intspan crate's `utils`) ──────

/// Convert a per-chromosome `IntSpan` map into a runlist JSON map.
///
/// Each value is the runlist string of its `IntSpan` (e.g. `"1-3,5"`).
pub fn set2json(
    set: &std::collections::BTreeMap<String, IntSpan>,
) -> std::collections::BTreeMap<String, serde_json::Value> {
    let mut json: std::collections::BTreeMap<String, serde_json::Value> =
        std::collections::BTreeMap::new();
    for (chr, value) in set {
        let runlist = value.to_string();
        json.insert(chr.into(), serde_json::to_value(runlist).unwrap());
    }
    json
}

/// Convert a nested map (`name` -> chromosome -> `IntSpan`) into a runlist
/// JSON map with each name holding its own chromosome map.
pub fn set2json_m(
    set_of: &std::collections::BTreeMap<String, std::collections::BTreeMap<String, IntSpan>>,
) -> std::collections::BTreeMap<String, serde_json::Value> {
    let mut out_json: std::collections::BTreeMap<String, serde_json::Value> =
        std::collections::BTreeMap::new();
    for (name, set) in set_of {
        let json = set2json(set);
        out_json.insert(name.to_string(), serde_json::to_value(json).unwrap());
    }
    out_json
}

/// Pretty-print a runlist JSON map to `output` (`stdout` supported).
pub fn write_json(
    output: &str,
    json: &std::collections::BTreeMap<String, serde_json::Value>,
) -> anyhow::Result<()> {
    use std::io::Write;
    let mut writer = crate::writer(output)?;
    let mut s = serde_json::to_string_pretty(json)?;
    s.push('\n');
    writer.write_all(s.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod display {
    use super::*;

    #[test]
    fn test_display() {
        let ints = IntSpan::from("1-3,5,7,9,100-999,1001-10000");
        assert_eq!(ints.to_string(), "1-3,5,7,9,100-999,1001-10000");

        let empty = IntSpan::new();
        assert_eq!(empty.to_string(), "-");
    }

    #[test]
    fn infinity_predicates_on_empty_set() {
        // `is_neg_inf`/`is_pos_inf` used to panic on empty sets (unwrap of
        // `front()`/`back()`); the empty set is finite and neither infinity.
        let empty = IntSpan::new();
        assert!(!empty.is_neg_inf());
        assert!(!empty.is_pos_inf());
        assert!(!empty.is_infinite());
        assert!(empty.is_finite());
        assert!(!empty.is_universal());
    }
}

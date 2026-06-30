use std::ops::Sub;

use anyhow::anyhow;
use intspan::IntSpan;

/// Coordinate transforming - from chr to align
///
/// ```
/// use pgr::libs::alignment::{indel_intspan, seq_intspan};
/// let tests : Vec<(&str, i32, i32, &str, i32)> = vec![
///     // seq, pos, chr_start, strand, expected
///     ("AAAATTTTTG", 4, 1, "+", 4),
///     ("AAAATTTTTG", 4, 1, "-", 7),
///     ("-AA--TTTGG", 5, 1, "+", 8),
///     ("-AA--TTTGG", 5, 1, "-", 6),
///     ("-AA--TTTGG", 105, 101, "+", 8),
///     ("-AA--TTTGG", 105, 101, "-", 6),
/// ];
/// for (seq, pos, chr_start, strand, expected) in tests {
///     let ints = seq_intspan(seq.as_ref());
///     // eprintln!("ints.to_string() = {:#?}", ints.to_string());
///     let result = pgr::libs::alignment::chr_to_align(&ints, pos, chr_start, strand).unwrap();
///     assert_eq!(result, expected);
/// }
/// ```
pub fn chr_to_align(ints: &IntSpan, pos: i32, chr_start: i32, strand: &str) -> anyhow::Result<i32> {
    let chr_end = chr_start + ints.size() - 1;

    if pos < chr_start || pos > chr_end {
        return Err(anyhow!(
            "[{}] out of ranges [{}, {}]",
            pos,
            chr_start,
            chr_end
        ));
    }

    let aln_pos = match strand {
        "+" => ints.at(pos - chr_start + 1),
        "-" => ints.at(-(pos - chr_start + 1)),
        _ => {
            return Err(anyhow!("Unrecognized strand: {}", strand));
        }
    };

    Ok(aln_pos)
}

/// Coordinate transforming - from align to chr
///
/// ```
/// use pgr::libs::alignment::{indel_intspan, seq_intspan};
/// let data : Vec<(&str, i32, i32, &str, i32)> = vec![
///     // seq, pos, chr_start, strand, expected
///     ("AAAATTTTTG", 4, 1, "+", 4),
///     ("AAAATTTTTG", 4, 1, "-", 7),
///     ("AAAATTTTTG", 4, 101, "+", 104),
///     ("AAAATTTTTG", 4, 101, "-", 107),
///     ("-AA--TTTGG", 6, 1, "+", 3),
///     ("-AA--TTTGG", 6, 1, "-", 5),
///     ("-AA--TTTGG", 6, 101, "+", 103),
///     ("-AA--TTTGG", 6, 101, "-", 105),
///     ("-AA--TTTGG", 1, 1, "+", 1),
///     ("-AA--TTTGG", 1, 1, "-", 7),
///     ("-AA--TTTGG-", 10, 1, "+", 7),
///     ("-AA--TTTGG-", 10, 1, "-", 1),
///     ("-AA--TTTGG", 4, 101, "+", 102),
///     ("-AA--TTTGG", 4, 101, "-", 106),
/// ];
/// for (seq, pos, chr_start, strand, expected) in data {
///     let ints = seq_intspan(seq.as_ref());
///     // eprintln!("ints.to_string() = {:#?}", ints.to_string());
///     let result = pgr::libs::alignment::align_to_chr(&ints, pos, chr_start, strand).unwrap();
///     assert_eq!(result, expected);
/// }
/// ```
pub fn align_to_chr(ints: &IntSpan, pos: i32, chr_start: i32, strand: &str) -> anyhow::Result<i32> {
    let chr_end = chr_start + ints.size() - 1;

    if pos < 1 {
        return Err(anyhow!("align pos [{}] out of ranges", pos,));
    }

    let mut chr_pos = if ints.contains(pos) {
        ints.index(pos)
    } else if pos < ints.min() {
        1
    } else if pos > ints.max() {
        ints.size()
    } else {
        // pos is in the holes
        // pins to the left base
        let spans = ints.spans();
        let mut cursor = pos;
        for i in 0..spans.len() {
            if spans[i].1 < cursor {
                continue;
            } else {
                cursor = spans[i - 1].1;
                break;
            }
        }

        ints.index(cursor)
    };

    chr_pos = match strand {
        "+" => chr_pos + chr_start - 1,
        "-" => chr_end - chr_pos + 1,
        _ => {
            return Err(anyhow!("Unrecognized strand: {}", strand));
        }
    };

    Ok(chr_pos)
}

/// Reverse a `[start, end]` range against a total `size`, in place.
///
/// Both endpoints are mapped as `new = size - old`, so a forward range
/// becomes a reverse-strand range (and vice versa).
///
/// ```
/// let mut s = 10;
/// let mut e = 20;
/// pgr::reverse_range(&mut s, &mut e, 100);
/// assert_eq!(s, 80);
/// assert_eq!(e, 90);
/// ```
pub fn reverse_range<T: Copy + Sub<Output = T>>(start: &mut T, end: &mut T, size: T) {
    let s = *start;
    let e = *end;
    *start = size - e;
    *end = size - s;
}

/// Reverse a 0-based half-open `[start, end)` against `size`, returning the
/// reversed range as a tuple `(size - end, size - start)`. Non-mutating
/// counterpart of [`reverse_range`].
///
/// ```
/// let (s, e) = pgr::libs::alignment::coords::reverse_range_pair(10, 20, 100);
/// assert_eq!(s, 80);
/// assert_eq!(e, 90);
/// ```
pub fn reverse_range_pair<T: Copy + Sub<Output = T>>(start: T, end: T, size: T) -> (T, T) {
    (size - end, size - start)
}

/// Reverse a 1-based inclusive `[start, end]` range against a total `size`,
/// in place. The forward 1-based inclusive range becomes a reverse-strand
/// 1-based inclusive range: `start -> size - end + 1`, `end -> size - start + 1`.
///
/// ```
/// let mut s = 11;
/// let mut e = 20;
/// pgr::libs::alignment::coords::reverse_range_1based(&mut s, &mut e, 100);
/// assert_eq!(s, 81);
/// assert_eq!(e, 90);
/// ```
pub fn reverse_range_1based(start: &mut usize, end: &mut usize, size: usize) {
    let s = *start;
    let e = *end;
    *start = size - e + 1;
    *end = size - s + 1;
}

/// Reverse a 1-based inclusive `[start, end]` range against `size`, returning the
/// reversed range as a tuple `(size - end + 1, size - start + 1)`. Non-mutating
/// counterpart of [`reverse_range_1based`].
///
/// ```
/// let (s, e) = pgr::libs::alignment::coords::reverse_range_1based_pair(11, 20, 100);
/// assert_eq!(s, 81);
/// assert_eq!(e, 90);
/// ```
pub fn reverse_range_1based_pair(start: usize, end: usize, size: usize) -> (usize, usize) {
    (size - end + 1, size - start + 1)
}

use intspan::Range;

use crate::libs::loc::{fetch_range_seq, open_indexed};

/// ```
/// let seq = pgr::libs::alignment::get_seq_loc("tests/fas/NC_000932.fa", "NC_000932:1-10").unwrap();
/// assert_eq!(seq, "ATGGGCGAAC".to_string());
/// let seq = pgr::libs::alignment::get_seq_loc("tests/fas/NC_000932.fa", "NC_000932(-):1-10").unwrap();
/// assert_eq!(seq, "GTTCGCCCAT".to_string());
/// let res = pgr::libs::alignment::get_seq_loc("tests/fas/NC_000932.fa", "FAKE:1-10");
/// assert_eq!(res.unwrap(), "".to_string());
/// ```
// cargo test --doc alignment::get_seq_loc
pub fn get_seq_loc(file: &str, range: &str) -> anyhow::Result<String> {
    let range = Range::from_str(range);
    if !range.is_valid() {
        return Ok("".to_string());
    }

    let (mut reader, loc_of) = open_indexed(file, false)?;

    if !loc_of.contains_key(range.chr()) {
        return Ok("".to_string());
    }

    let seq = fetch_range_seq(&mut reader, &loc_of, &range)?;

    Ok(seq)
}

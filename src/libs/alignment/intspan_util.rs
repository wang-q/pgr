use intspan::IntSpan;

/// ```
/// use pgr::libs::alignment::{indel_intspan, seq_intspan};
/// let tests : Vec<(&str, &str)> = vec![
///     // seq, expected
///     ("ATAA", "-"),
///     ("CcGc", "-"),
///     ("TAGggATaaC", "-"),
///     ("C-Gc", "2"),
///     ("C--c", "2-3"),
///     ("---c", "1-3"),
///     ("C---", "2-4"),
///     ("GCaN--NN--NNNaC", "5-6,9-10"),
/// ];
/// for (seq, expected) in tests {
///     let result = indel_intspan(seq.as_ref());
///     assert_eq!(result.to_string(), expected.to_string());
/// }
/// ```
pub fn indel_intspan(seq: &[u8]) -> IntSpan {
    let mut positions = vec![];

    for (i, base) in seq.iter().enumerate() {
        if *base == b'-' {
            positions.push(i as i32 + 1);
        }
    }

    let mut ints = IntSpan::new();
    ints.add_vec(&positions);

    ints
}

pub fn seq_intspan(seq: &[u8]) -> IntSpan {
    IntSpan::from_pair(1, seq.len() as i32).diff(&indel_intspan(seq))
}

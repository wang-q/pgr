use super::condensed::{get_condensed_index, CondensedMatrix};
use super::named::NamedMatrix;
use super::scoring::ScoringMatrix;

#[test]
fn test_condensed_matrix_indexing() {
    // N=4
    // (0,1) -> 0
    // (0,2) -> 1
    // (0,3) -> 2
    // (1,2) -> 3
    // (1,3) -> 4
    // (2,3) -> 5
    let m = CondensedMatrix::new(4);
    assert_eq!(get_condensed_index(m.size(), 0, 1), 0);
    assert_eq!(get_condensed_index(m.size(), 0, 2), 1);
    assert_eq!(get_condensed_index(m.size(), 0, 3), 2);
    assert_eq!(get_condensed_index(m.size(), 1, 2), 3);
    assert_eq!(get_condensed_index(m.size(), 1, 3), 4);
    assert_eq!(get_condensed_index(m.size(), 2, 3), 5);
}

#[test]
fn test_condensed_matrix_rw() {
    let mut m = CondensedMatrix::new(3);
    m.set(0, 1, 1.0);
    m.set(2, 0, 2.0); // set (0,2) via swap
    m.set(1, 2, 3.0);

    assert_eq!(m.get(0, 1), 1.0);
    assert_eq!(m.get(1, 0), 1.0);
    assert_eq!(m.get(0, 2), 2.0);
    assert_eq!(m.get(2, 0), 2.0);
    assert_eq!(m.get(1, 2), 3.0);
    assert_eq!(m.get(0, 0), 0.0);

    // Test underlying data access
    let data = m.data();
    assert_eq!(data.len(), 3); // 3*2/2 = 3
                               // Order: (0,1), (0,2), (1,2) -> 1.0, 2.0, 3.0
    assert_eq!(data[0], 1.0);
    assert_eq!(data[1], 2.0);
    assert_eq!(data[2], 3.0);
}

#[test]
fn test_condensed_matrix_from_vec() {
    let data = vec![1.0, 2.0, 3.0];
    let m = CondensedMatrix::from_vec(3, data);
    assert_eq!(m.get(0, 1), 1.0);
    assert_eq!(m.get(0, 2), 2.0);
    assert_eq!(m.get(1, 2), 3.0);
}

#[test]
#[should_panic(expected = "Data length 2 does not match expected length 3 for size 3")]
fn test_condensed_matrix_from_vec_invalid_len() {
    CondensedMatrix::from_vec(3, vec![1.0, 2.0]);
}

#[test]
fn test_scoring_matrix_basic() {
    let mut m = ScoringMatrix::with_defaults(0.0, -1.0);
    m.set(0, 1, 5.0);
    m.set(2, 1, 10.0);

    // Check set values (symmetric)
    assert_eq!(m.get(0, 1), 5.0);
    assert_eq!(m.get(1, 0), 5.0);
    assert_eq!(m.get(1, 2), 10.0);

    // Check diagonal default
    assert_eq!(m.get(0, 0), 0.0);
    assert_eq!(m.get(3, 3), 0.0);

    // Check missing default
    assert_eq!(m.get(0, 2), -1.0);
    assert_eq!(m.get(3, 4), -1.0);
}

#[test]
fn test_named_matrix_basic() {
    let names = vec!["A".to_string(), "B".to_string()];
    let mut m = NamedMatrix::new(names);

    m.set(0, 1, 0.5);
    assert_eq!(m.get(0, 1), 0.5);
    assert_eq!(m.get(1, 0), 0.5);
    assert_eq!(m.get(0, 0), 0.0);

    assert_eq!(m.get_by_name("A", "B"), Some(0.5));
}

#[test]
fn test_named_matrix_indexing() {
    let names = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let m = NamedMatrix::new(names);

    // Size 3 -> len 3
    assert_eq!(m.values().len(), 3);

    // Index check
    // (0,1) -> 0
    // (0,2) -> 1
    // (1,2) -> 2
    assert_eq!(m.index(0, 1), 0);
    assert_eq!(m.index(0, 2), 1);
    assert_eq!(m.index(1, 2), 2);
}

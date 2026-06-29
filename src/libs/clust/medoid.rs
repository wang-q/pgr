//! Medoid selection for flat clustering output.
//!
//! A medoid is the member of a cluster whose sum of pairwise distances
//! (or similarities) to the other members is extremal. For distance
//! matrices the medoid minimizes the sum; for similarity matrices it
//! maximizes the sum. Ties are broken by the iteration order — callers
//! should pass `members` sorted by name so the alphabetically-first
//! member wins.

use crate::libs::pairmat::ScoringMatrix;

/// Find the medoid of a cluster.
///
/// Iterates `members` and returns the index of the member whose sum of
/// `matrix.get(candidate, member)` over all `members` is minimal
/// (`find_max = false`, distance matrix) or maximal (`find_max = true`,
/// similarity matrix). Returns `None` for an empty `members` slice.
///
/// Tie-breaking: the first member achieving the extremal sum wins, so
/// callers should sort `members` by name beforehand for deterministic
/// alphabetical tie-breaking.
pub fn find_medoid(
    matrix: &ScoringMatrix<f32>,
    members: &[usize],
    find_max: bool,
) -> Option<usize> {
    if members.is_empty() {
        return None;
    }
    let mut best_rep = members[0];
    let mut best_sum = if find_max {
        f32::NEG_INFINITY
    } else {
        f32::MAX
    };

    for &candidate in members {
        let mut current_sum = 0.0;
        for &member in members {
            current_sum += matrix.get(candidate, member);
        }
        if find_max {
            if current_sum > best_sum {
                best_sum = current_sum;
                best_rep = candidate;
            }
        } else if current_sum < best_sum {
            best_sum = current_sum;
            best_rep = candidate;
        }
    }

    Some(best_rep)
}

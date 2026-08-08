//! Profile -> repeat-run extraction (Profex `-z` equivalent).

use anyhow::Context;
use std::io::Write;

/// Write one `prof.<sn>.rg` file per chromosome from profile runs.
///
/// A run is a maximal stretch of equal profile values above zero (`0` is a
/// separator, matching FastK/Profex). Each run becomes a 1-based inclusive
/// `chr:start-end` line: `start = first k-mer position + 1`, `end = start0 +
/// len + k - 1`. The profile spans the whole chromosome, so the final run is
/// closed correctly (its last k-mer covers through the sequence end) without
/// needing the chromosome length — fixing the old Profex quirk that dropped
/// or guessed tail runs. `min_depth` skips runs below the threshold; `None`
/// keeps every run. `rg_files` receives the written file names.
pub fn write_rg(
    profiles: &[Vec<u16>],
    chrs: &[String],
    k: usize,
    min_depth: Option<u16>,
    rg_files: &mut Vec<String>,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        profiles.len() == chrs.len(),
        "profile/chromosome count mismatch ({} vs {})",
        profiles.len(),
        chrs.len()
    );
    for (i, profile) in profiles.iter().enumerate() {
        let rg_file = format!("prof.{}.rg", i + 1);
        let mut writer = crate::writer(&rg_file)?;
        for (start, end) in constant_runs(profile, min_depth) {
            writer
                .write_fmt(format_args!("{}:{}-{}\n", chrs[i], start + 1, end + k - 1))
                .with_context(|| format!("writing {rg_file}"))?;
        }
        rg_files.push(rg_file);
    }
    Ok(())
}

/// Runs of constant profile value above zero (and at or above `min_depth`),
/// as 0-based half-open `[start, end)` k-mer position ranges.
fn constant_runs(profile: &[u16], min_depth: Option<u16>) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut start = 0usize;
    while start < profile.len() {
        let value = profile[start];
        if value == 0 || min_depth.is_some_and(|m| value < m) {
            start += 1;
            continue;
        }
        let mut end = start + 1;
        while end < profile.len() && profile[end] == value {
            end += 1;
        }
        runs.push((start, end));
        start = end;
    }
    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_split_on_constant_value() {
        // Profile [2,2,0,1,1,1] with k=3: run1 start0=0 len=2 -> 1-4,
        // run2 start0=3 len=3 -> 4-8.
        assert_eq!(
            constant_runs(&[2, 2, 0, 1, 1, 1], None),
            vec![(0, 2), (3, 6)]
        );
    }

    #[test]
    fn tail_run_is_kept() {
        // The old Profex quirk dropped/guessed tail runs; the native extract
        // must close the final run from the profile alone. A 24 bp sequence
        // with k=4 has 21 windows; a constant profile yields one run ending
        // at 24 (start0 + len + k - 1 = 0 + 21 + 4 - 1).
        assert_eq!(constant_runs(&[2u16; 21], None), vec![(0, 21)]);
    }

    #[test]
    fn min_depth_filters_runs() {
        // Values 3, 1, 3 split into three runs; min_depth 2 keeps the 3s.
        assert_eq!(
            constant_runs(&[3, 3, 1, 1, 3], Some(2)),
            vec![(0, 2), (4, 5)]
        );
    }

    #[test]
    fn zero_separates_and_empty_profiles_emit_nothing() {
        assert_eq!(
            constant_runs(&[0u16, 0], None),
            Vec::<(usize, usize)>::new()
        );
        assert_eq!(constant_runs(&[], None), Vec::<(usize, usize)>::new());
    }
}

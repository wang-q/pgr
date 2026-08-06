use std::collections::BTreeMap;
use std::io::Write;

use crate::libs::ds::IntSpan;
use anyhow::anyhow;

use crate::libs::alignment::{align_to_chr, chr_to_align, indel_intspan, seq_intspan};
use crate::libs::fmt::fas::FasBlock;

/// Slice a FasBlock by a set of chromosome runlists, writing each subslice
/// to `writer` as one `>range\nseq\n` entry per species.
///
/// `name` is the reference species whose range determines the chr lookup in
/// `set`. Returns `Ok(())` if no slicing happened (e.g., name not found,
/// chr not in `set`, or empty intersection).
pub fn slice_block<W: Write>(
    block: &FasBlock,
    name: &str,
    set: &BTreeMap<String, IntSpan>,
    writer: &mut W,
) -> anyhow::Result<()> {
    let idx = match block.names.iter().position(|x| x == name) {
        Some(i) => i,
        None => return Ok(()),
    };
    let trange = block.entries[idx].range().clone();

    // chr present in the requested set
    let i_ints_chr = match set.get(trange.chr()) {
        Some(s) if !s.is_empty() => trange.intspan().intersect(s),
        _ => return Ok(()),
    };
    if i_ints_chr.is_empty() {
        return Ok(());
    }

    // target sequence intspan
    let t_ints_seq = seq_intspan(block.entries[idx].seq());

    // per-species align intspans + shared indel regions
    let mut ints_seq_of: BTreeMap<&str, IntSpan> = BTreeMap::new();
    let mut indel_ints = IntSpan::new();
    for (i, n) in block.names.iter().enumerate() {
        let seq = block.entries[i].seq();
        ints_seq_of.insert(n.as_str(), seq_intspan(seq));
        indel_ints.merge(&indel_intspan(seq));
    }

    // collect subslices (chr-position intersections)
    let mut sub_slices: Vec<IntSpan> = vec![];
    for (lower, upper) in i_ints_chr.spans() {
        // On the reverse strand, increasing chr positions map to decreasing
        // alignment columns, so the two endpoints come out reversed. Normalize
        // to the [min, max] column span so the subslice is well-formed.
        let (mut ss_start, mut ss_end) = match (
            chr_to_align(&t_ints_seq, lower, trange.start, trange.strand()),
            chr_to_align(&t_ints_seq, upper, trange.start, trange.strand()),
        ) {
            (Ok(l), Ok(u)) => (l, u),
            // A reference species with its own gaps covers fewer genomic
            // positions than its range length, so a requested position can
            // fall outside the reference's non-gap span. Skip that subspan
            // (like an out-of-range coordinate) instead of aborting the whole
            // slice run.
            (l, u) => {
                let msg = l
                    .err()
                    .or(u.err())
                    .map(|e| e.to_string())
                    .unwrap_or_default();
                log::warn!("skipping slice subspan {}-{}: {}", lower, upper, msg);
                continue;
            }
        };
        if ss_start > ss_end {
            std::mem::swap(&mut ss_start, &mut ss_end);
        }
        // A single-base request maps both endpoints to the same column, so the
        // subslice is one column wide (a valid, length-1 slice). An empty span
        // is handled below by the `ss_ints.is_empty()` check.
        let mut ss_ints = IntSpan::from_pair(ss_start, ss_end);

        // trim indel borders
        for n in [ss_start, ss_end] {
            if indel_ints.contains(n) {
                let island = indel_ints.find_islands_n(n);
                ss_ints.subtract(&island);
            }
        }
        if ss_ints.is_empty() {
            // The whole subslice fell inside an indel island (every column has
            // a gap in some species), so there is nothing to slice out.
            continue;
        }
        sub_slices.push(ss_ints);
    }

    // emit entries per subslice per species
    for ss in &sub_slices {
        let ss_start = ss.min();
        let ss_end = ss.max();

        for (i, n) in block.names.iter().enumerate() {
            let range = block.entries[i].range();
            let mut start = align_to_chr(
                ints_seq_of.get(n.as_str()).unwrap(),
                ss_start,
                range.start,
                range.strand(),
            )?;
            let mut end = align_to_chr(
                ints_seq_of.get(n.as_str()).unwrap(),
                ss_end,
                range.start,
                range.strand(),
            )?;
            // On the reverse strand, the leftmost alignment column maps to the
            // largest chr coordinate, so the endpoints come out backwards.
            // Report the genomic span as [min, max] with the strand preserved.
            if start > end {
                std::mem::swap(&mut start, &mut end);
            }
            let ss_range = crate::libs::ds::Range::from_full(
                range.name(),
                range.chr(),
                range.strand(),
                start,
                end,
            );

            let seq = block.entries[i].seq();
            let seq_len_i32 = i32::try_from(seq.len())
                .map_err(|_| anyhow!("sequence length {} exceeds i32 range", seq.len()))?;
            if ss_start < 1 || ss_end > seq_len_i32 {
                anyhow::bail!(
                    "slice range {}..{} out of sequence length {}",
                    ss_start,
                    ss_end,
                    seq.len()
                );
            }
            let start_idx = (ss_start - 1) as usize;
            let end_idx = ss_end as usize;
            let ss_seq = &seq[start_idx..end_idx];

            let seq_str = std::str::from_utf8(ss_seq)
                .map_err(|e| anyhow!("invalid UTF-8 in sliced sequence: {}", e))?;
            writer.write_all(format!(">{}\n{}\n", ss_range, seq_str).as_ref())?;
        }
    }

    // blank line separating blocks
    writer.write_all(b"\n")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::ds::Range;
    use crate::libs::fmt::fas::{FasBlock, FasEntry};

    fn block_with_all_gap_second_species() -> FasBlock {
        let target = Range::from_full("T", "chr1", "+", 1, 4);
        let other = Range::from_full("O", "chr2", "+", 1, 4);
        FasBlock {
            entries: vec![
                FasEntry::from(&target, b"ACGT"),
                FasEntry::from(&other, b"----"),
            ],
            names: vec!["T".to_string(), "O".to_string()],
            headers: vec![target.to_string(), other.to_string()],
        }
    }

    /// Regression: when the whole subslice is covered by an indel island (a
    /// second species is entirely gaps), trimming the indel borders emptied
    /// the subslice and `ss.min()`/`ss.max()` panicked on the empty IntSpan.
    #[test]
    fn slice_block_all_gap_second_species_no_panic() {
        let block = block_with_all_gap_second_species();
        let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
        set.insert("chr1".to_string(), IntSpan::from("1-4"));

        let mut out: Vec<u8> = vec![];
        // Target has no gaps and the subslice covers all four columns, while
        // the other species is all gaps -> the subslice is a pure indel island.
        let result = slice_block(&block, "T", &set, &mut out);
        assert!(
            result.is_ok(),
            "slicing must not panic or error: {:?}",
            result
        );
    }

    /// Regression: a reverse-strand reference must still produce a subslice
    /// (previously the reversed endpoints made `ss_start >= ss_end`, so every
    /// subslice was dropped), and a reverse-strand non-reference species must
    /// report its genomic span as `start-end` with `start <= end` (previously
    /// it emitted a backwards range like `9-4`).
    #[test]
    fn slice_block_reverse_strand_reports_valid_ranges() {
        let target = Range::from_full("Ref", "chr1", "-", 100, 112);
        let other = Range::from_full("Oth", "chr2", "-", 1, 12);
        let block = FasBlock {
            entries: vec![
                FasEntry::from(&target, b"ACGTACGTACGT"),
                FasEntry::from(&other, b"TGCATGCATGCA"),
            ],
            names: vec!["Ref".to_string(), "Oth".to_string()],
            headers: vec![target.to_string(), other.to_string()],
        };
        let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
        set.insert("chr1".to_string(), IntSpan::from("103-108"));

        let mut out: Vec<u8> = vec![];
        slice_block(&block, "Ref", &set, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        // Reverse-strand reference must be sliced (was previously empty).
        assert!(
            text.contains(">Ref.chr1(-):103-108\nTACGTA\n"),
            "got: {text}"
        );
        // Reverse-strand non-reference species range must be non-backwards.
        assert!(text.contains(">Oth.chr2(-):4-9\nATGCAT\n"), "got: {text}");
    }

    /// Regression: a reference species whose own sequence contains a gap
    /// (fewer non-gap bases than its genomic range length) makes
    /// `chr_to_align` reject a runlist position beyond
    /// `chr_start + non_gap - 1`. Slicing must skip such a subspan instead
    /// of aborting the whole command.
    #[test]
    fn slice_block_gapped_reference_no_abort() {
        let target = Range::from_full("Ref", "chr1", "+", 1, 5);
        let other = Range::from_full("Oth", "chr2", "+", 1, 5);
        let block = FasBlock {
            entries: vec![
                // Reference has a gap at column 2 -> only 4 non-gap bases over
                // a 5-base genomic range.
                FasEntry::from(&target, b"A-CGT"),
                FasEntry::from(&other, b"ATGCA"),
            ],
            names: vec!["Ref".to_string(), "Oth".to_string()],
            headers: vec![target.to_string(), other.to_string()],
        };
        let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
        set.insert("chr1".to_string(), IntSpan::from("1-5"));

        let mut out: Vec<u8> = vec![];
        let result = slice_block(&block, "Ref", &set, &mut out);
        assert!(
            result.is_ok(),
            "slicing a gapped reference must skip, not abort: {:?}",
            result
        );
    }

    #[test]
    fn slice_block_single_base_subslice() {
        let target = Range::from_full("Ref", "chr1", "+", 1, 4);
        let other = Range::from_full("Oth", "chr2", "+", 1, 4);
        let block = FasBlock {
            entries: vec![
                FasEntry::from(&target, b"ACGT"),
                FasEntry::from(&other, b"TGCA"),
            ],
            names: vec!["Ref".to_string(), "Oth".to_string()],
            headers: vec![target.to_string(), other.to_string()],
        };
        let mut set: BTreeMap<String, IntSpan> = BTreeMap::new();
        set.insert("chr1".to_string(), IntSpan::from("2"));

        let mut out: Vec<u8> = vec![];
        slice_block(&block, "Ref", &set, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains(">Ref.chr1(+):2\nC\n"), "got: {text}");
        assert!(text.contains(">Oth.chr2(+):2\nG\n"), "got: {text}");
    }
}

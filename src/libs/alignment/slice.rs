use std::collections::BTreeMap;
use std::io::Write;

use anyhow::anyhow;
use intspan::IntSpan;

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
        let ss_start = chr_to_align(&t_ints_seq, lower, trange.start, trange.strand())?;
        let ss_end = chr_to_align(&t_ints_seq, upper, trange.start, trange.strand())?;
        if ss_start >= ss_end {
            continue;
        }
        let mut ss_ints = IntSpan::from_pair(ss_start, ss_end);

        // trim indel borders
        for n in [ss_start, ss_end] {
            if indel_ints.contains(n) {
                let island = indel_ints.find_islands_n(n);
                ss_ints.subtract(&island);
            }
        }
        sub_slices.push(ss_ints);
    }

    // emit entries per subslice per species
    for ss in &sub_slices {
        let ss_start = ss.min();
        let ss_end = ss.max();

        for (i, n) in block.names.iter().enumerate() {
            let range = block.entries[i].range();
            let start = align_to_chr(
                ints_seq_of.get(n.as_str()).unwrap(),
                ss_start,
                range.start,
                range.strand(),
            )?;
            let end = align_to_chr(
                ints_seq_of.get(n.as_str()).unwrap(),
                ss_end,
                range.start,
                range.strand(),
            )?;
            let ss_range =
                intspan::Range::from_full(range.name(), range.chr(), range.strand(), start, end);

            let ss_seq = &block.entries[i].seq()[((ss_start - 1) as usize)..(ss_end as usize)];

            let seq_str = std::str::from_utf8(ss_seq)
                .map_err(|e| anyhow!("invalid UTF-8 in sliced sequence: {}", e))?;
            writer.write_all(format!(">{}\n{}\n", ss_range, seq_str).as_ref())?;
        }
    }

    // blank line separating blocks
    writer.write_all(b"\n")?;
    Ok(())
}

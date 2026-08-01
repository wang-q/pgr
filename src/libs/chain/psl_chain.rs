use crate::libs::chain::{
    calc_block_score, chain_blocks, clean_input_overlaps, Chain, ChainableBlock, GapCalc,
    ScoreContext, SubMatrix,
};
use crate::libs::fmt::psl::Psl;
use crate::libs::io::SequenceReader;
use indexmap::IndexMap;
use std::io::{BufRead, Write};
use std::str::FromStr;

/// PSL alignment blocks grouped by (target, query, strand) for chaining.
pub struct GroupData {
    /// Target sequence size.
    pub t_size: u32,
    /// Query sequence size.
    pub q_size: u32,
    /// Alignment blocks in this group.
    pub blocks: Vec<ChainableBlock>,
}

/// Group key: (target_name, query_name, query_strand).
pub type GroupKey = (String, String, char);

/// Read PSL records and group alignment blocks by (target, query, strand).
///
/// Returns the grouped blocks plus any `#` comment lines collected from the
/// input header (UCSC axtChain propagates these via lineFileSetMetaDataOutput).
pub fn group_psl_blocks<R: BufRead, S: SequenceReader>(
    reader: R,
    score_ctx: &mut Option<ScoreContext<S>>,
) -> anyhow::Result<(IndexMap<GroupKey, GroupData>, Vec<String>)> {
    // IndexMap preserves the first-seen (PSL input) order of the groups, which
    // UCSC axtChain mirrors by prepending new seqPairs with slAddHead.
    let mut groups: IndexMap<GroupKey, GroupData> = IndexMap::new();
    let mut comments: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            comments.push(line);
            continue;
        }

        let psl = match Psl::from_str(&line) {
            Ok(p) => p,
            Err(_) => continue, // Skip invalid lines (e.g. headers)
        };

        let t_name = psl.t_name.clone();
        let q_name = psl.q_name.clone();
        let q_strand = psl.strand.chars().next().unwrap_or('+');

        let key = (t_name.clone(), q_name.clone(), q_strand);
        let entry = groups.entry(key).or_insert_with(|| GroupData {
            t_size: psl.t_size,
            q_size: psl.q_size,
            blocks: Vec::new(),
        });

        if psl.strand.len() > 1 && psl.strand.chars().nth(1) == Some('-') {
            log::warn!(
                "Skipping PSL record with negative target strand: {} {} {}",
                psl.q_name,
                psl.strand,
                psl.t_name
            );
            continue;
        }

        for i in 0..psl.block_count as usize {
            let size = psl.block_sizes[i] as u64;
            let t_start = psl.t_starts[i] as u64;
            let t_end = t_start + size;

            let (q_start, q_end) = {
                let s = psl.q_starts[i] as u64;
                (s, s + size)
            };

            let mut block = ChainableBlock {
                t_start,
                t_end,
                q_start,
                q_end,
                score: size as f64 * 100.0,
            };

            if let Some(ctx) = score_ctx.as_mut() {
                if let Ok(exact) =
                    calc_block_score(&block, ctx, &q_name, &t_name, psl.q_size as u64, q_strand)
                {
                    block.score = exact;
                }
            }

            entry.blocks.push(block);
        }
    }

    Ok((groups, comments))
}

/// Chain PSL alignments and write chains filtered by `min_score`.
///
/// Reads PSL records, groups by (target, query, strand), chains each group
/// via dynamic programming, sorts chains by descending score, and writes
/// chains with score >= `min_score` to `writer`.
pub fn chain_psl<R: BufRead, W: Write, S: SequenceReader>(
    reader: R,
    writer: &mut W,
    gap_calc: &GapCalc,
    min_score: f64,
    score_context: &mut Option<ScoreContext<S>>,
    matrix: &SubMatrix,
) -> anyhow::Result<()> {
    let (mut groups, comments) = group_psl_blocks(reader, score_context)?;

    let mut all_chains: Vec<Chain> = Vec::new();
    let mut chain_id_counter = 1;

    // UCSC axtChain accumulates each pair's chains with slAddHead, so the
    // global list before the score sort is (first-seen pair first, chains in
    // reverse order): readPslBlocks prepends new pairs, so the first-read pair
    // is chained last and its chains end up at the head of the global list via
    // slAddHead.  slSort is stable, so equal-score chains keep that order.
    // Mirror it by iterating groups in first-seen order and reversing each
    // group's chains before the final stable score sort.
    let group_order: Vec<GroupKey> = groups.keys().cloned().collect();
    for (t_name, q_name, q_strand) in group_order {
        let mut data = groups
            .shift_remove(&(t_name.clone(), q_name.clone(), q_strand))
            .expect("group key present");
        if data.blocks.is_empty() {
            continue;
        }

        // UCSC axtChain cleans the per-pair block list (slReverse +
        // removeExactOverlaps, which sorts by qStart then tStart and folds
        // same-start blocks), then chainBlocks sorts the leaves by tStart.
        // kent's slSort is stable for equal keys, so replicating both steps
        // gives the same tie order (and therefore the same KD-tree
        // split/search order) for the DP.
        clean_input_overlaps(&mut data.blocks);
        data.blocks.sort_by_key(|a| a.t_start);

        log::debug!("Group: {} {} {}", t_name, q_name, q_strand);
        for b in &data.blocks {
            log::debug!(
                "Block: T {}-{} Q {}-{} Score {}",
                b.t_start,
                b.t_end,
                b.q_start,
                b.q_end,
                b.score
            );
        }

        let chains = chain_blocks(
            &data.blocks,
            gap_calc,
            score_context,
            &q_name,
            data.q_size as u64,
            q_strand,
            &t_name,
            data.t_size as u64,
            &mut chain_id_counter,
        )?;
        all_chains.extend(chains.into_iter().rev());
    }

    all_chains.sort_by(|a, b| b.header.score.total_cmp(&a.header.score));

    // UCSC axtChain sorts all chains globally by score, then assigns IDs
    // 1, 2, 3, ... in that order via chainWriteHead -> chainIdNext.
    // Renumber to match, so downstream tools see the same ID ordering.
    for (i, chain) in all_chains.iter_mut().enumerate() {
        chain.header.id = (i + 1) as u64;
    }

    // UCSC axtChain writes the scoring scheme header first (via
    // axtScoreSchemeDnaWrite), then propagates `##` metadata from the PSL
    // input (via lineFileSetMetaDataOutput).
    writeln!(writer, "{}", matrix.axt_chain_header())?;
    for comment in &comments {
        writeln!(writer, "{}", comment)?;
    }

    for chain in all_chains {
        if chain.header.score < min_score {
            continue;
        }
        chain.write(writer)?;
    }

    Ok(())
}

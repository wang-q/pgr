use crate::libs::chain::{
    calc_block_score, chain_blocks, Chain, ChainableBlock, GapCalc, ScoreContext,
};
use crate::libs::fmt::psl::Psl;
use crate::libs::io::SequenceReader;
use std::collections::HashMap;
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
pub fn group_psl_blocks<R: BufRead, S: SequenceReader>(
    reader: R,
    score_ctx: &mut Option<ScoreContext<S>>,
) -> anyhow::Result<HashMap<GroupKey, GroupData>> {
    let mut groups: HashMap<GroupKey, GroupData> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
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

    Ok(groups)
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
) -> anyhow::Result<()> {
    let groups = group_psl_blocks(reader, score_context)?;

    let mut all_chains: Vec<Chain> = Vec::new();
    let mut chain_id_counter = 1;

    for ((t_name, q_name, q_strand), mut data) in groups {
        if data.blocks.is_empty() {
            continue;
        }

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
        all_chains.extend(chains);
    }

    all_chains.sort_by(|a, b| b.header.score.total_cmp(&a.header.score));

    for chain in all_chains {
        if chain.header.score < min_score {
            continue;
        }
        chain.write(writer)?;
    }

    Ok(())
}

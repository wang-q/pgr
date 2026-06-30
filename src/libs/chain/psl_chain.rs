use crate::libs::chain::{calc_block_score, ChainableBlock, ScoreContext};
use crate::libs::fmt::psl::Psl;
use crate::libs::io::SequenceReader;
use std::collections::HashMap;
use std::io::BufRead;
use std::str::FromStr;

/// PSL alignment blocks grouped by (target, query, strand) for chaining.
pub struct GroupData {
    pub t_size: u32,
    pub q_size: u32,
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
                if let Some(exact) =
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

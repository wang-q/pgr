//! Write pairwise or MSA block FASTA from PAF query results.

use std::io::Write;

use anyhow::anyhow;

use super::fasta::FastaStore;
use super::index::PafIndex;
use super::msa_build::{build_msa_entries, build_pairwise_block, run_poa_msa};
use super::query::QueryGroup;
use crate::libs::poa::AlignmentParams;

/// Write pairwise FAS blocks (target + query per result) to `writer`.
pub fn write_pairwise_fas<W: Write>(
    idx: &PafIndex,
    all_results: &[QueryGroup],
    fasta_store: &mut FastaStore,
    writer: &mut W,
) -> anyhow::Result<()> {
    for (_, results) in all_results {
        for result in results {
            let blk = build_pairwise_block(idx, result, fasta_store)?;
            let t_start_display = blk
                .t_start
                .checked_add(1)
                .ok_or_else(|| anyhow!("t_start {} overflow on display", blk.t_start))?;
            writeln!(
                writer,
                ">{0}(+):{1}-{2}",
                blk.tname, t_start_display, blk.t_end
            )?;
            writeln!(writer, "{}", blk.t_aln)?;
            writeln!(
                writer,
                ">{0}({1}):{2}-{3}",
                blk.qname,
                blk.q_strand,
                blk.q_start_fwd.saturating_add(1),
                blk.q_end_fwd
            )?;
            writeln!(writer, "{}", blk.q_aln)?;
            writeln!(writer)?;
        }
    }
    Ok(())
}

/// Write multi-way MSA FAS blocks (one block per region) to `writer`.
pub fn write_msa_fas<W: Write>(
    idx: &PafIndex,
    all_results: &[QueryGroup],
    fasta_store: &mut FastaStore,
    params: AlignmentParams,
    writer: &mut W,
) -> anyhow::Result<()> {
    for ((tname_region, _, _), results) in all_results {
        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;
        if entries.is_empty() {
            continue;
        }

        let msa = run_poa_msa(&entries, params.clone());

        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = i32::try_from(aln.chars().filter(|c| *c != '-').count())
                .map_err(|_| anyhow!("alignment length exceeds i32 range"))?;
            writeln!(
                writer,
                ">{0}({3}):{1}-{2}",
                e.name,
                e.start.saturating_add(1),
                e.start.saturating_add(size),
                e.strand
            )?;
            writeln!(writer, "{}", aln)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

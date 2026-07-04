//! MAF output helpers for PAF query results.

use std::io::Write;

use crate::libs::paf::fasta::FastaStore;
use crate::libs::paf::index::{PafIndex, QueryResult};
use crate::libs::paf::msa_build::{build_msa_entries, build_pairwise_block, run_poa_msa};
use crate::libs::poa::AlignmentParams;

/// One region's target coords plus its query results, as produced by `run_query`.
pub type RegionResults = ((String, i32, i32), Vec<QueryResult>);

/// Write pairwise MAF blocks: one `a` block per `QueryResult`, two `s` lines each.
/// Sequences are fetched on demand via `FastaStore`; CIGAR is walked directly
/// (no POA refinement).
pub fn write_pairwise_maf(
    idx: &PafIndex,
    all_results: &[RegionResults],
    fasta_store: &mut FastaStore,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    writeln!(writer, "##maf version=1")?;
    for (_, results) in all_results {
        for result in results {
            let blk = build_pairwise_block(idx, result, fasta_store)?;

            // size = number of non-gap bases
            let q_size = blk.q_aln.chars().filter(|c| *c != '-').count();
            let t_size = blk.t_aln.chars().filter(|c| *c != '-').count();

            writeln!(writer, "a")?;
            writeln!(
                writer,
                "s\t{0}\t{1}\t{2}\t+\t{3}\t{4}",
                blk.tname, blk.t_start, t_size, blk.t_src_size, blk.t_aln
            )?;
            writeln!(
                writer,
                "s\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                blk.qname, blk.q_start_maf, q_size, blk.q_strand, blk.q_src_size, blk.q_aln
            )?;
            writeln!(writer)?;
        }
    }
    Ok(())
}

/// Write multi-way MAF blocks via POA: one `a` block per region with N `s` lines.
/// For each region, collect target + all query sequences (queries RC'd if '-' strand),
/// feed them into the POA engine, and emit the resulting MSA. CIGAR is ignored.
pub fn write_msa_maf(
    idx: &PafIndex,
    all_results: &[RegionResults],
    fasta_store: &mut FastaStore,
    params: AlignmentParams,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    writeln!(writer, "##maf version=1")?;
    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;
        let msa = run_poa_msa(&entries, params.clone());

        // Emit the MAF block.
        writeln!(writer, "a")?;
        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = aln.chars().filter(|c| *c != '-').count();
            writeln!(
                writer,
                "s\t{}\t{}\t{}\t{}\t{}\t{}",
                e.name, e.start, size, e.strand, e.src_size, aln
            )?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

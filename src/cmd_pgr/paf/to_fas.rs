use clap::*;
use pgr::libs::nt;
use pgr::libs::paf::fasta::FastaStore;
use pgr::libs::paf::index::PafIndex;
use pgr::libs::paf::msa::{build_maf_block, build_msa_entries};

use super::common::{self, QueryGroup};
use super::query;

fn output_fas_pairwise(
    idx: &PafIndex,
    all_results: &[QueryGroup],
    fasta_store: &mut FastaStore,
) -> anyhow::Result<()> {
    for (_, results) in all_results {
        for (query_id, q_iv, t_iv, cigar, rec_ts, rec_qs, strand) in results {
            let qname = idx.id_to_name(*query_id).unwrap_or("?");
            let tname = idx.id_to_name(t_iv.metadata).unwrap_or("?");

            let (qs, qe) = if q_iv.first <= q_iv.last {
                (q_iv.first, q_iv.last)
            } else {
                (q_iv.last, q_iv.first)
            };
            let (ts, te) = if t_iv.first <= t_iv.last {
                (t_iv.first, t_iv.last)
            } else {
                (t_iv.last, t_iv.first)
            };

            let (q_seq_fwd, q_src_size) = fasta_store.fetch_range(qname, qs, qe)?;
            let (t_seq, _t_src_size) = fasta_store.fetch_range(tname, ts, te)?;

            let (q_seq_for_aln, rec_qs_eff, qs_eff, q_strand, _) = if *strand == '-' {
                let rc = nt::rev_comp(&q_seq_fwd).collect::<Vec<u8>>();
                let aligned_q_len: i32 = cigar.iter().map(|op| op.query_delta() as i32).sum();
                let rec_qe = *rec_qs + aligned_q_len;
                let rc_sub_start = rec_qe - qe;
                (rc, 0, rc_sub_start, '-', q_src_size as i32 - qe)
            } else {
                (q_seq_fwd, *rec_qs, qs, '+', qs)
            };

            let (q_aln, t_aln) = build_maf_block(
                cigar,
                *rec_ts,
                rec_qs_eff,
                ts,
                te,
                qs_eff,
                &q_seq_for_aln,
                &t_seq,
            );

            println!(">{tname}(+):{}-{}", ts + 1, te);
            println!("{}", t_aln);
            println!(">{qname}({}):{}-{}", q_strand, qs + 1, qe);
            println!("{}", q_aln);
            println!();
        }
    }
    Ok(())
}

fn output_fas_msa(
    idx: &PafIndex,
    all_results: &[QueryGroup],
    fasta_store: &mut FastaStore,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> anyhow::Result<()> {
    for ((tname_region, _, _), results) in all_results {
        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;
        if entries.is_empty() {
            continue;
        }

        let params = pgr::libs::poa::AlignmentParams {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        };
        let mut poa = pgr::libs::poa::Poa::new(params, pgr::libs::poa::AlignmentType::Global);
        for e in &entries {
            poa.add_sequence(&e.seq);
        }
        let msa = poa.msa();

        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = aln.chars().filter(|c| *c != '-').count() as i32;
            println!(
                ">{0}({3}):{1}-{2}",
                e.name,
                e.start + 1,
                e.start + size,
                e.strand
            );
            println!("{}", aln);
        }
        println!();
    }
    Ok(())
}

pub fn make_subcommand() -> Command {
    query::add_poa_args(query::add_query_args(
        Command::new("to-fas")
            .arg(
                Arg::new("fasta_tsv")
                    .long("fasta-tsv")
                    .short('f')
                    .required(true)
                    .num_args(1)
                    .help("TSV file: genome_name <tab> bgzf_fasta_path"),
            )
            .arg(
                Arg::new("msa")
                    .long("msa")
                    .num_args(0)
                    .help("Merge results per region into a multi-way block FASTA via POA"),
            ),
    ))
    .about("Query PAF index and output pairwise or multi-way block FASTA")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs block FASTA records.

Default mode (pairwise): each query result becomes a block of two FASTA
records (target first, query second) restored directly from CIGAR.
Alignments are assumed to be already refined by chain/net — no POA
refinement is performed.

--msa mode (multi-way): merge all query results of each region into a
single multi-sequence block FASTA via POA. Sequences (target first, then
each query, '-' strand reverse-complemented) are fed into the POA
engine; CIGAR is ignored. Best used with --transitive to gather all
homologous fragments of a region.

Output format (per block):
  >seq_name(+):start-end
  ATGC--ATGC
  >seq_name(-):start-end
  ATGCAT--GC
  (blank line)

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  A FASTA file may be referenced by multiple genome_names (multi-chrom).
  All genome names in the PAF index must be present in the TSV (strict
  validation — missing entries cause a hard error).

Notes:
* Input PAF files should contain cg:Z: tags for accurate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'
* Output is compatible with `pgr fas to-vcf`

Examples:
1. Single region to pairwise FAS:
   pgr paf to-fas alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way MSA with transitive BFS:
   pgr paf to-fas alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv

3. Pipeline to VCF:
   pgr paf to-fas alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv | pgr fas to-vcf

4. Batch query from BED regions:
   pgr paf to-fas alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let (idx, all_results, mut fasta_store) = common::prepare_query(args)?;
    if args.get_flag("msa") {
        let params = common::get_poa_params(args);
        output_fas_msa(
            &idx,
            &all_results,
            &mut fasta_store,
            params.match_score,
            params.mismatch_score,
            params.gap_open,
            params.gap_extend,
        )
    } else {
        output_fas_pairwise(&idx, &all_results, &mut fasta_store)
    }
}

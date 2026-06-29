use clap::*;
use pgr::libs::nt;
use pgr::libs::paf::fasta::{load_fasta_tsv, FastaStore};
use pgr::libs::paf::index::{PafIndex, QueryResult};
use pgr::libs::paf::msa::{build_maf_block, build_msa_entries};

use super::query;

// Output pairwise MAF blocks. Each QueryResult becomes one `a` block with two
// `s` lines (target first, query second). Sequences are fetched on demand via
// FastaStore; CIGAR is walked directly (no POA refinement).
#[allow(clippy::type_complexity)]
fn output_maf(
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut FastaStore,
) -> anyhow::Result<()> {
    println!("##maf version=1");
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
            let (t_seq, t_src_size) = fasta_store.fetch_range(tname, ts, te)?;

            // For '-' strand records: PAF query coords are on the forward
            // strand, but CIGAR describes alignment columns against the
            // reverse-complemented query. RC the fetched forward sequence
            // and walk CIGAR from offset 0 so column order matches.
            //
            // `q_seq_for_aln` is RC(forward[qs..qe)) and covers RC offset
            // [rec_qe - qe, rec_qe - qs) where rec_qe = rec_qs + aligned_q_len.
            // `build_maf_block` indexes q_seq via `(cq + skip_t) - qs_eff`, so
            // qs_eff must be the RC offset of the sub-interval start
            // (rec_qe - qe), not 0 — otherwise sub-interval queries index
            // past the end of q_seq. Full-overlap queries have qs_eff = 0.
            //
            // MAF `start` for '-' strand = srcSize - qe (position on forward
            // strand of the first displayed base, per MAF spec).
            let (q_seq_for_aln, rec_qs_eff, qs_eff, q_strand, q_start_maf) = if *strand == '-' {
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

            // size = number of non-gap bases
            let q_size = q_aln.chars().filter(|c| *c != '-').count();
            let t_size = t_aln.chars().filter(|c| *c != '-').count();

            println!("a");
            println!("s\t{tname}\t{ts}\t{t_size}\t+\t{t_src_size}\t{t_aln}");
            println!("s\t{qname}\t{q_start_maf}\t{q_size}\t{q_strand}\t{q_src_size}\t{q_aln}");
            println!();
        }
    }
    Ok(())
}

// Output multi-way MAF blocks via POA. For each region, collect target +
// all query sequences (queries RC'd if '-' strand), feed them into the POA
// engine, and emit one `a` block with N `s` lines. CIGAR is ignored.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn output_maf_msa(
    idx: &PafIndex,
    all_results: &[((String, i32, i32), Vec<QueryResult>)],
    fasta_store: &mut FastaStore,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> anyhow::Result<()> {
    let params = pgr::libs::poa::AlignmentParams {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };

    println!("##maf version=1");
    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;

        // Run POA MSA.
        let mut poa =
            pgr::libs::poa::Poa::new(params.clone(), pgr::libs::poa::AlignmentType::Global);
        for e in &entries {
            poa.add_sequence(&e.seq);
        }
        let msa = poa.msa();

        // Emit the MAF block.
        println!("a");
        for (e, aln) in entries.iter().zip(msa.iter()) {
            let size = aln.chars().filter(|c| *c != '-').count();
            println!(
                "s\t{}\t{}\t{}\t{}\t{}\t{}",
                e.name, e.start, size, e.strand, e.src_size, aln
            );
        }
        println!();
    }
    Ok(())
}

pub fn make_subcommand() -> Command {
    query::add_poa_args(query::add_query_args(
        Command::new("to-maf")
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
                    .help("Merge results per region into a multi-way MAF block via POA"),
            ),
    ))
    .about("Query PAF index and output pairwise or multi-way MAF")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs MAF blocks.

Default mode (pairwise): each query result becomes one 2-sequence MAF
block restored directly from CIGAR. Alignments are assumed to be
already refined by chain/net — no POA refinement is performed.

--msa mode (multi-way): merge all query results of each region into a
single multi-sequence MAF block via POA. Sequences (target first, then
each query, '-' strand reverse-complemented) are fed into the POA
engine; CIGAR is ignored. Best used with --transitive to gather all
homologous fragments of a region.

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

Examples:
1. Single region to pairwise MAF:
   pgr paf to-maf alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way MSA with transitive BFS:
   pgr paf to-maf alignments.paf chr1:1000-5000 -t --msa -f genomes.tsv

3. Batch query from BED regions:
   pgr paf to-maf alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = query::run_query(args)?;

    pgr::libs::paf::fasta::validate_tsv_covers_index(&seq_to_file, &idx)?;

    let mut fasta_store = FastaStore::new(&seq_to_file)?;
    if args.get_flag("msa") {
        let match_score = *args.get_one::<i32>("match_score").unwrap();
        let mismatch_score = *args.get_one::<i32>("mismatch_score").unwrap();
        let gap_open = *args.get_one::<i32>("gap_open").unwrap();
        let gap_extend = *args.get_one::<i32>("gap_extend").unwrap();
        output_maf_msa(
            &idx,
            &all_results,
            &mut fasta_store,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        )?;
    } else {
        output_maf(&idx, &all_results, &mut fasta_store)?;
    }
    Ok(())
}

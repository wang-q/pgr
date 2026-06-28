use clap::*;
use std::io::Write;

use pgr::libs::paf::index::PafIndex;
use pgr::libs::paf::index::QueryResult;
use pgr::libs::poa::Poa;

use super::query;
use super::to_maf::{build_msa_entries, load_fasta_tsv, FastaStore};

pub fn make_subcommand() -> Command {
    query::add_query_args(
        Command::new("to-gfa")
            .arg(
                Arg::new("fasta_tsv")
                    .long("fasta-tsv")
                    .short('f')
                    .required(true)
                    .num_args(1)
                    .help("TSV file: genome_name <tab> bgzf_fasta_path"),
            )
            .arg(
                Arg::new("match_score")
                    .long("match")
                    .num_args(1)
                    .default_value("5")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA match score (default: 5)"),
            )
            .arg(
                Arg::new("mismatch_score")
                    .long("mismatch")
                    .num_args(1)
                    .default_value("-4")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA mismatch score (default: -4)"),
            )
            .arg(
                Arg::new("gap_open")
                    .long("gap-open")
                    .num_args(1)
                    .default_value("-8")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA gap open penalty (default: -8)"),
            )
            .arg(
                Arg::new("gap_extend")
                    .long("gap-extend")
                    .num_args(1)
                    .default_value("-6")
                    .value_parser(clap::value_parser!(i32))
                    .allow_negative_numbers(true)
                    .help("POA gap extend penalty (default: -6)"),
            )
            .arg(
                Arg::new("outfile")
                    .long("outfile")
                    .short('o')
                    .num_args(1)
                    .default_value("stdout")
                    .help("Output filename. [stdout] for screen"),
            ),
    )
    .about("Query PAF index and output local GFA via POA graph")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs a local GFA (v1.0) graph induced from a POA multiple sequence
alignment of each region.

For each region, all homologous fragments (target first, then each
query, '-' strand reverse-complemented) are fed into the POA engine.
The POA graph itself — nodes are bases, edges are adjacencies, paths
trace each input sequence — is exported directly as GFA S/L/P lines.
Bubbles (SNPs, indels) appear as graph branches, preserving the full
variant structure without lossy MSA-to-GFA conversion.

Each region produces an independent GFA block (node IDs restart at 1).
Multiple regions are separated by `# region: <name>` comment lines.

Recommended with --transitive to gather all homologous fragments of
each region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  All genome names in the PAF index must be present in the TSV.

Notes:
* Nodes are single bases (one base per S line); unchopping is not done
* Input PAF files should contain cg:Z: tags (used for query projection)
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to local GFA:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way GFA with transitive BFS:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -t -f genomes.tsv

3. Batch query from BED regions:
   pgr paf to-gfa alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = query::run_query(args)?;

    // Validate: every name in the PAF index must be present in the TSV.
    let mut missing: Vec<&str> = idx
        .names
        .keys()
        .filter(|n| !seq_to_file.contains_key(*n))
        .map(|n| n.as_str())
        .collect();
    missing.sort_unstable();
    if !missing.is_empty() {
        anyhow::bail!(
            "FASTA TSV is missing {} genome(s) present in PAF index: {}",
            missing.len(),
            missing.join(", ")
        );
    }

    let mut fasta_store = FastaStore::new(&seq_to_file)?;

    let match_score = *args.get_one::<i32>("match_score").unwrap();
    let mismatch_score = *args.get_one::<i32>("mismatch_score").unwrap();
    let gap_open = *args.get_one::<i32>("gap_open").unwrap();
    let gap_extend = *args.get_one::<i32>("gap_extend").unwrap();

    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    output_gfa(
        &mut writer,
        &idx,
        &all_results,
        &mut fasta_store,
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    )?;

    writer.flush()?;
    Ok(())
}

// Output local GFA per region from POA graph. The POA graph (nodes = bases,
// edges = adjacencies, paths = per-sequence traversals) maps directly to
// GFA S/L/P lines. Each region gets an independent GFA block with node IDs
// restarting at 1; multiple regions are separated by `# region:` comments.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn output_gfa<W: Write>(
    writer: &mut W,
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

    let multi_region = all_results.iter().filter(|(_, r)| !r.is_empty()).count() > 1;

    for ((tname_region, _, _), results) in all_results {
        if results.is_empty() {
            continue;
        }

        let entries = build_msa_entries(idx, tname_region, results, fasta_store)?;

        // Run POA to build the graph.
        let mut poa = Poa::new(params.clone(), pgr::libs::poa::AlignmentType::Global);
        for e in &entries {
            poa.add_sequence(&e.seq);
        }

        // Region separator comment for multi-region output.
        if multi_region {
            writeln!(writer, "# region: {tname_region}")?;
        }

        write_poa_as_gfa(writer, &poa, &entries)?;
    }

    Ok(())
}

// Convert a POA graph to GFA S/L/P lines. Node IDs are 1-based, assigned by
// petgraph NodeIndex order. Each node is a single base. Edges become L lines
// with `0M` overlap (adjacency). Each input sequence's path becomes a P line.
fn write_poa_as_gfa<W: Write>(
    writer: &mut W,
    poa: &Poa,
    entries: &[super::to_maf::MsaEntry],
) -> anyhow::Result<()> {
    let graph = poa.graph().graph.node_indices().collect::<Vec<_>>();

    // Build a NodeIndex -> 1-based GFA id map.
    let max_idx = graph.iter().map(|n| n.index()).max().unwrap_or(0);
    let mut node_id: Vec<u32> = vec![0; max_idx + 1];
    for (i, n) in graph.iter().enumerate() {
        node_id[n.index()] = (i + 1) as u32;
    }

    // S lines: one per node, sequence is the single base.
    for n in &graph {
        let id = node_id[n.index()];
        let base = poa.graph().graph.node_weight(*n).unwrap().base;
        let base_char = char::from(base);
        writeln!(writer, "S\t{id}\t{base_char}")?;
    }

    // L lines: one per edge. Orientation is always '+' (POA graph is a DAG
    // on the forward strand; '-' strand queries were RC'd before POA).
    for edge_ref in poa.graph().graph.edge_indices() {
        let (from, to) = poa.graph().graph.edge_endpoints(edge_ref).unwrap();
        let from_id = node_id[from.index()];
        let to_id = node_id[to.index()];
        writeln!(writer, "L\t{from_id}\t+\t{to_id}\t+\t0M")?;
    }

    // P lines: one per input sequence (entries), using the recorded path.
    let paths = poa.paths();
    for (i, entry) in entries.iter().enumerate() {
        let path = &paths[i];
        if path.is_empty() {
            continue;
        }
        let path_str: Vec<String> = path
            .iter()
            .map(|n| format!("{}+", node_id[n.index()]))
            .collect();
        // Overlaps: 0M between each consecutive pair (adjacent single-base nodes).
        let overlaps = vec!["0M"; path.len().saturating_sub(1)];
        writeln!(
            writer,
            "P\t{}\t{}\t{}",
            entry.name,
            path_str.join(","),
            overlaps.join(",")
        )?;
    }

    Ok(())
}

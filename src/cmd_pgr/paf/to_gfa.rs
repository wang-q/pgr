use clap::*;
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

use petgraph::graph::NodeIndex;
use petgraph::visit::NodeIndexable;
use petgraph::Direction;

use pgr::libs::paf::index::PafIndex;
use pgr::libs::paf::index::QueryResult;
use pgr::libs::poa::Poa;

use super::query;
use pgr::libs::paf::fasta::{load_fasta_tsv, FastaStore};
use pgr::libs::paf::msa::{build_msa_entries, MsaEntry};

pub fn make_subcommand() -> Command {
    query::add_poa_args(query::add_query_args(
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
                Arg::new("crush")
                    .long("crush")
                    .action(ArgAction::SetTrue)
                    .help("Compress SNP bubbles (impg `crush` style; loses base-level ALT info)"),
            )
            .arg(crate::cmd_pgr::args::outfile_arg()),
    ))
    .about("Query PAF index and output local GFA via POA graph")
    .after_help(
        r###"
Queries a PAF file or saved index (same logic as `pgr paf query`) and
outputs a local GFA (v1.0) graph induced from a POA multiple sequence
alignment of each region.

For each region, all homologous fragments (target first, then each
query, '-' strand reverse-complemented) are fed into the POA engine.
The POA graph — nodes are bases, edges are adjacencies, paths trace
each input sequence — is then compacted and exported as GFA S/L/P.

Compaction (unchopping): linear stretches of single-base nodes with
no branching are merged into multi-base segments, reducing node count
by ~1 order of magnitude. SNP/indel bubbles remain as graph branches,
preserving the full variant structure.

--crush: optional impg `crush` style post-processing. Compresses simple
bubbles (nodes sharing the same in/out neighbors) into a single node,
keeping the highest-weight allele. Useful for SV overview graphs but
LOSES base-level ALT information — paths through ALTs are rewritten to
the kept allele. Off by default.

Each region produces an independent GFA block (node IDs restart at 1).
Multiple regions are separated by `# region: <name>` comment lines.

Recommended with --transitive to gather all homologous fragments of
each region.

-f/--fasta-tsv (required):
  TSV with two columns: genome_name <tab> bgzf_fasta_path
  Each genome_name must match a query/target name in the PAF index.
  All genome names in the PAF index must be present in the TSV.

Notes:
* Input PAF files should contain cg:Z: tags (used for query projection)
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Single region to local GFA:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -f genomes.tsv

2. Multi-way GFA with transitive BFS:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -t -f genomes.tsv

3. SV overview graph with bubbles crushed:
   pgr paf to-gfa alignments.paf chr1:1000-5000 -t -f genomes.tsv --crush

4. Batch query from BED regions:
   pgr paf to-gfa alignments.paf.idx -b regions.bed -f genomes.tsv

"###,
    )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let tsv_path = args.get_one::<String>("fasta_tsv").unwrap();
    let seq_to_file = load_fasta_tsv(tsv_path)?;

    let (idx, all_results) = query::run_query(args)?;

    pgr::libs::paf::fasta::validate_tsv_covers_index(&seq_to_file, &idx)?;

    let mut fasta_store = FastaStore::new(&seq_to_file)?;

    let match_score = *args.get_one::<i32>("match_score").unwrap();
    let mismatch_score = *args.get_one::<i32>("mismatch_score").unwrap();
    let gap_open = *args.get_one::<i32>("gap_open").unwrap();
    let gap_extend = *args.get_one::<i32>("gap_extend").unwrap();
    let crush = args.get_flag("crush");

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
        crush,
    )?;

    writer.flush()?;
    Ok(())
}

// Output local GFA per region from POA graph. The POA graph (nodes = bases,
// edges = adjacencies, paths = per-sequence traversals) is compacted (linear
// stretches merged) and exported as GFA S/L/P. Each region gets an independent
// GFA block with node IDs restarting at 1; multiple regions are separated by
// `# region:` comments.
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
    crush: bool,
) -> anyhow::Result<()> {
    let params = pgr::libs::poa::AlignmentParams {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };

    let multi_region = all_results.iter().filter(|(_, r)| !r.is_empty()).count() > 1;

    // GFA header (once, at the top).
    writeln!(writer, "H\tVN:Z:1.0")?;

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

        let mut graph = build_compacted(&poa);
        if crush {
            graph = crush_bubbles(graph);
        }
        write_compacted_gfa(writer, &graph, &entries)?;
    }

    Ok(())
}

// Compacted graph: linear stretches of single-base POA nodes merged into
// multi-base segments. Segment ids are 1-based (GFA convention).
struct CompactedGraph {
    // Segment sequences, indexed by (id - 1).
    segments: Vec<String>,
    // Segment weights (sum of POA node weights in the segment). Used by
    // `crush_bubbles` to pick the kept allele.
    weights: Vec<u32>,
    // Deduplicated edges between segment ids (1-based).
    edges: BTreeSet<(u32, u32)>,
    // Per-path segment id lists (adjacent duplicates removed).
    paths: Vec<Vec<u32>>,
}

// Compact (unchop) the POA graph: merge linear stretches of single-base nodes
// into multi-base segments. A node v merges into its predecessor p's segment
// when p has out-degree 1 (sole successor v) and v has in-degree 1 (sole
// predecessor p). This preserves all branching structure (bubbles remain as
// graph branches) while collapsing non-branching runs.
fn build_compacted(poa: &Poa) -> CompactedGraph {
    let g = &poa.graph().graph;
    let topo = poa.graph().topological_sort();

    // 1. Compute segment head for each node (topological order guarantees
    //    head[p] is final when we process v).
    let mut head: Vec<NodeIndex> = vec![NodeIndex::new(0); g.node_bound()];
    for &v in &topo {
        let h = {
            let preds: Vec<_> = g.neighbors_directed(v, Direction::Incoming).collect();
            if preds.len() == 1 {
                let p = preds[0];
                let succs: Vec<_> = g.neighbors_directed(p, Direction::Outgoing).collect();
                if succs.len() == 1 && succs[0] == v {
                    head[p.index()]
                } else {
                    v
                }
            } else {
                v
            }
        };
        head[v.index()] = h;
    }

    // 2. Assign 1-based ids to segment heads in topological order of first
    //    appearance (the head itself is always first in its segment).
    let mut head_to_id: HashMap<NodeIndex, u32> = HashMap::new();
    let mut next_id: u32 = 1;
    for &v in &topo {
        let h = head[v.index()];
        head_to_id.entry(h).or_insert_with(|| {
            let id = next_id;
            next_id += 1;
            id
        });
    }

    let n_segs = (next_id - 1) as usize;
    let mut segments: Vec<String> = vec![String::new(); n_segs];
    let mut weights: Vec<u32> = vec![0; n_segs];

    // 3. Build segment sequences by appending bases in topological order.
    for &v in &topo {
        let h = head[v.index()];
        let id = (head_to_id[&h] - 1) as usize;
        let base = char::from(g.node_weight(v).unwrap().base);
        segments[id].push(base);
        weights[id] += g.node_weight(v).unwrap().weight;
    }

    // 4. Edges: map original edges to (head_from, head_to), dedup.
    let mut edges: BTreeSet<(u32, u32)> = BTreeSet::new();
    for e in g.edge_indices() {
        let (u, v) = g.edge_endpoints(e).unwrap();
        let hu = head_to_id[&head[u.index()]];
        let hv = head_to_id[&head[v.index()]];
        if hu != hv {
            edges.insert((hu, hv));
        }
    }

    // 5. Paths: map nodes to segment ids, remove adjacent dups.
    let mut paths: Vec<Vec<u32>> = Vec::with_capacity(poa.paths().len());
    for path in poa.paths() {
        let mut compact: Vec<u32> = Vec::with_capacity(path.len());
        let mut last: Option<u32> = None;
        for n in path {
            let id = head_to_id[&head[n.index()]];
            if Some(id) != last {
                compact.push(id);
                last = Some(id);
            }
        }
        paths.push(compact);
    }

    CompactedGraph {
        segments,
        weights,
        edges,
        paths,
    }
}

// Crush simple bubbles (impg `crush` style). A bubble is a set of nodes
// sharing the same in-neighbor set and out-neighbor set; all but the
// highest-weight member are removed and their path references rewritten to
// the survivor. This loses base-level ALT info — paths through ALTs are
// rewritten to the kept allele.
fn crush_bubbles(graph: CompactedGraph) -> CompactedGraph {
    let n = graph.segments.len();
    if n == 0 {
        return graph;
    }

    // Build adjacency sets (1-based ids).
    let mut in_adj: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); n];
    let mut out_adj: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); n];
    for &(u, v) in &graph.edges {
        out_adj[(u - 1) as usize].insert(v);
        in_adj[(v - 1) as usize].insert(u);
    }

    // Group nodes by (in_set, out_set) signature.
    let mut groups: HashMap<(BTreeSet<u32>, BTreeSet<u32>), Vec<u32>> = HashMap::new();
    for id in 1..=n as u32 {
        let sig = (
            in_adj[(id - 1) as usize].clone(),
            out_adj[(id - 1) as usize].clone(),
        );
        groups.entry(sig).or_default().push(id);
    }

    // remap[old_id - 1] = new_id (the survivor of its bubble, or itself).
    let mut remap: Vec<u32> = (1..=n as u32).collect();
    for members in groups.values() {
        if members.len() < 2 {
            continue;
        }
        // Keep max weight; tie-break: lowest id (deterministic).
        let keep = *members
            .iter()
            .max_by_key(|&&id| (graph.weights[(id - 1) as usize], i64::MIN + id as i64))
            .unwrap();
        for &id in members {
            if id != keep {
                remap[(id - 1) as usize] = keep;
            }
        }
    }

    // If no crushable bubbles, return as-is.
    if remap.iter().enumerate().all(|(i, &r)| r == (i + 1) as u32) {
        return graph;
    }

    // Assign new sequential ids to survivors, in original id order.
    let mut new_id: Vec<u32> = vec![0; n];
    let mut next: u32 = 1;
    // Follow remap chains to the survivor.
    let mut survivor: Vec<u32> = vec![0; n];
    for id in 1..=n as u32 {
        let mut cur = id;
        while remap[(cur - 1) as usize] != cur {
            cur = remap[(cur - 1) as usize];
        }
        survivor[(id - 1) as usize] = cur;
    }
    for id in 1..=n as u32 {
        if survivor[(id - 1) as usize] == id {
            new_id[(id - 1) as usize] = next;
            next += 1;
        }
    }

    // Build new segments/weights.
    let mut new_segments: Vec<String> = Vec::with_capacity(next as usize - 1);
    let mut new_weights: Vec<u32> = Vec::with_capacity(next as usize - 1);
    for id in 1..=n as u32 {
        if survivor[(id - 1) as usize] == id {
            new_segments.push(graph.segments[(id - 1) as usize].clone());
            new_weights.push(graph.weights[(id - 1) as usize]);
        }
    }

    // Build new edges (remap + dedup).
    let mut new_edges: BTreeSet<(u32, u32)> = BTreeSet::new();
    for &(u, v) in &graph.edges {
        let su = new_id[(survivor[(u - 1) as usize] as usize - 1) as usize];
        let sv = new_id[(survivor[(v - 1) as usize] as usize - 1) as usize];
        if su != sv {
            new_edges.insert((su, sv));
        }
    }

    // Build new paths (remap, then collapse adjacent dups).
    let mut new_paths: Vec<Vec<u32>> = Vec::with_capacity(graph.paths.len());
    for path in &graph.paths {
        let mut compact: Vec<u32> = Vec::with_capacity(path.len());
        let mut last: Option<u32> = None;
        for &id in path {
            let s = new_id[(survivor[(id - 1) as usize] as usize - 1) as usize];
            if Some(s) != last {
                compact.push(s);
                last = Some(s);
            }
        }
        new_paths.push(compact);
    }

    CompactedGraph {
        segments: new_segments,
        weights: new_weights,
        edges: new_edges,
        paths: new_paths,
    }
}

// Write a compacted graph as GFA S/L/P lines. S lines include the LN tag
// (segment length). L lines use `0M` overlap (segments are adjacent, no
// overlap). P lines use `0M` overlaps between consecutive segments.
fn write_compacted_gfa<W: Write>(
    writer: &mut W,
    graph: &CompactedGraph,
    entries: &[MsaEntry],
) -> anyhow::Result<()> {
    // S lines.
    for (i, seq) in graph.segments.iter().enumerate() {
        let id = (i + 1) as u32;
        writeln!(writer, "S\t{id}\t{seq}\tLN:i:{}", seq.len())?;
    }

    // L lines.
    for &(from, to) in &graph.edges {
        writeln!(writer, "L\t{from}\t+\t{to}\t+\t0M")?;
    }

    // P lines.
    for (i, entry) in entries.iter().enumerate() {
        if i >= graph.paths.len() {
            break;
        }
        let path = &graph.paths[i];
        if path.is_empty() {
            continue;
        }
        let path_str: Vec<String> = path.iter().map(|&id| format!("{id}+")).collect();
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

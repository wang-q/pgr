use clap::*;
use pgr::libs::paf::graph::PafGraph;
use std::collections::HashMap;

pub fn make_subcommand() -> Command {
    Command::new("graph")
        .about("Induces a coarse GFA graph from PAF alignments")
        .after_help(
            r###"
Builds a coarse pangenome graph (GFA v1.0) from pairwise PAF alignments
and a set of FASTA sequences, using a seqwish-style segment-level DSU
algorithm.

Algorithm:
* Splits each alignment at indels >= --min-var-len into match segments.
* Unions aligned segments via disjoint-set union (transitive closure).
* Derives graph nodes (DSU classes) + edges (path adjacencies) + novel
  segments (unaligned gaps), then emits S/L/P lines.

This builds a coarse global GFA: only large structural variations
(>= --min-var-len) split nodes; small indels stay within a node as
variations. For base-level regional graphs, see `pgr paf to-gfa`.

Notes:
* Input PAF files should contain cg:Z: tags for accurate splitting
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads PAF from stdin if input file is 'stdin'
* FASTA files (-f) are required to populate node sequences and lengths
* GFA node ids are 1-based; node 1 is the earliest segment by (seq, start)
* S lines carry rGFA tags: SN:Z (source seq), SO:i (0-based start), SR:i:0

Examples:
1. Build a coarse graph with default SV threshold (100bp):
   pgr paf graph alignments.paf -f refs.fa -o graph.gfa

2. Stricter threshold (only >= 500bp SVs split nodes):
   pgr paf graph alignments.paf -f refs.fa --min-var-len 500 -o graph.gfa

3. Read PAF from stdin:
   cat alignments.paf | pgr paf graph stdin -f refs.fa > graph.gfa

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PAF file (or 'stdin' for piped input)"),
        )
        .arg(
            Arg::new("fasta")
                .long("fasta")
                .short('f')
                .required(true)
                .num_args(1..)
                .help("FASTA file(s) providing sequence content and lengths"),
        )
        .arg(
            Arg::new("min_var_len")
                .long("min-var-len")
                .num_args(1)
                .default_value("100")
                .value_parser(clap::value_parser!(i32))
                .help("Minimum indel length to split at (default: 100)"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output GFA filename. [stdout] for screen"),
        )
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let infile = matches.get_one::<String>("infile").unwrap();
    let fasta_files: Vec<String> = matches
        .get_many::<String>("fasta")
        .unwrap()
        .cloned()
        .collect();
    let min_var_len = matches
        .get_one::<i32>("min_var_len")
        .copied()
        .unwrap_or(100);
    let outfile = matches.get_one::<String>("outfile").unwrap();

    // Load FASTA sequences into a name -> bytes map.
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    for fa_path in &fasta_files {
        let reader = pgr::reader(fa_path);
        let mut fa_in = noodles_fasta::io::Reader::new(reader);
        for result in fa_in.records() {
            let record = result?;
            let name = String::from_utf8(record.name().into())?;
            let seq_bytes: Vec<u8> = record.sequence()[..].to_vec();
            seqs.insert(name, seq_bytes);
        }
    }

    // Read PAF.
    let paf_reader = pgr::reader(infile);

    // Build the graph.
    let graph = PafGraph::build(paf_reader, &seqs, min_var_len)?;

    // Write GFA.
    let writer = pgr::writer(outfile);
    graph.write_gfa(writer)?;

    Ok(())
}

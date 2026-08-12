use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::asm::assemble::{assemble_unitigs, AssembleOptions};
use std::io::Write;

/// Build the clap subcommand for unitig.
pub fn make_subcommand() -> Command {
    Command::new("unitig")
        .about("Assembles reads into maximal unitigs (non-branching paths)")
        .after_help(
            r###"
This command assembles reads into maximal unitigs through the k-mer graph,
following the BCALM 2 compaction semantics (GATB `ograph.cpp` `graph3`):
every solid k-mer (count >= 3) extends in both directions only while it has
exactly one solid successor whose own predecessor is also unique, so the
assembly stops at branches, junctions, coverage gaps, and loops. Parallel
paths stay separate (no bubble popping), and the result is independent of
the k-mer scan order.

This is the strict graph-compression counterpart of `pgr asm contig`, whose
seeded contig mode keeps extending through weak branches (tadpole-compatible
behavior). Unitigs are best suited to high-coverage or error-corrected input,
such as the anchr `unitigs` step's `pe.cor.fa`.

Notes:
* Input is 1 interleaved file or 2 paired files; FASTA and FASTQ both work
* Unitigs are written longest-first with a `unitig_<id>` FASTA header
  carrying length, coverage, GC, and dimer composition fields
* Processing is ordered and deterministic, independent of scan order
* Output sequences are wrapped at 70 columns, like BBTools FASTA output
* Supports both plain text and gzipped (.gz) files

Examples:
1. Assemble unitigs from corrected reads (anchr unitigs step):
   pgr asm unitig pe.cor.fa -o unitigs_K31.fasta --kmer 31

2. Assemble from paired-end reads:
   pgr asm unitig R1.fq.gz R2.fq.gz -o unitigs.fasta

3. Raise the solid k-mer threshold (like bcalm `-abundance-min`):
   pgr asm unitig in.fq -o out.fasta --min-count-seed 5
"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input file(s): 1 interleaved or 2 paired (R1, R2)",
            1..=2,
        ))
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(clap::builder::RangedU64ValueParser::<usize>::new().range(1..))
                .help("K-mer length"),
        )
        .arg(
            Arg::new("min_contig_len")
                .long("min-contig-len")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum unitig length (default: max(124, 2*k))"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Solid k-mer count threshold (default 3, like bcalm -abundance-min)"),
        )
        .arg(
            Arg::new("links")
                .long("links")
                .action(ArgAction::SetTrue)
                .help("Append L: links to unitig FASTA headers (bcalm format)"),
        )
        .arg(
            Arg::new("gfa")
                .long("gfa")
                .action(ArgAction::SetTrue)
                .help("Emit a GFA graph instead of FASTA"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Accepted for compatibility; ignored (deterministic processing)"),
        )
}

/// Execute the unitig command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    // Reject `-o` that would overwrite an input file (the writer is opened
    // before the reads are consumed).
    crate::cmd_pgr::args::ensure_outfile_distinct(outfile, infiles.iter().map(|s| s.as_str()))?;
    crate::cmd_pgr::args::parse_parallel_auto(args.get_one::<String>("parallel").unwrap())?;
    let opts = AssembleOptions {
        k: *args.get_one::<usize>("kmer").unwrap(),
        min_contig_len: args
            .get_one::<usize>("min_contig_len")
            .copied()
            .unwrap_or(0),
        min_count_seed: args
            .get_one::<usize>("min_count_seed")
            .copied()
            .unwrap_or(3),
        emit_links: args.get_flag("links"),
        emit_gfa: args.get_flag("gfa"),
        ..AssembleOptions::default()
    };

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let stats = assemble_unitigs(&infiles, &mut out, &opts)?;
    out.flush()?;
    eprintln!(
        "Reads in: {}  Unitigs: {}  Bases: {}  Longest: {}",
        stats.reads_in, stats.contigs_built, stats.bases_built, stats.longest_contig
    );
    Ok(())
}

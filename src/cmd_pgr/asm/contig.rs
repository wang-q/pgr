use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::asm::assemble::{assemble, AssembleOptions};
use std::io::Write;

/// Build the clap subcommand for contig.
pub fn make_subcommand() -> Command {
    Command::new("contig")
        .about("Assembles reads into contigs via k-mer graph traversal (tadpole-compatible)")
        .after_help(
            r###"
This command assembles reads into contigs through the k-mer graph, reproducing
the BBTools `tadpole.sh` contig mode (the default mode when no `ecc`/`extend`
flag is set): k-mers are counted with a quality gate (`--min-prob`), contigs
are seeded from k-mers above a depth threshold and extended greedily in both
directions, stopping at branches and dead ends. This replaces the tadpole
assembly steps of the anchr `2_insert_size` and `unitigs` flows.

Notes:
* Input is 1 interleaved file or 2 paired files; FASTA and FASTQ both work
* Contigs are written longest-first with a `contig_<id>` FASTA header carrying
  length, coverage, GC, and dimer composition fields (BBTools SHORT_NAMES)
* Processing is ordered and deterministic (equivalent to `threads=1`)
* Bubble-popping resolutions may differ slightly from BBTools on some
  overlapping structures (its expand order depends on a memory-dependent
  hash layout); the contig set and total bases match, and the output is
  reproducible across runs
* Bubble popping is on by default (tadpole `popbubbles=t`); pass
  `--no-bubbles` to keep parallel-path contigs separate (tadpole
  `popbubbles=f`)
* Output sequences are wrapped at 70 columns, like BBTools FASTA output
* Supports both plain text and gzipped (.gz) files

Examples:
1. Assemble contigs from corrected reads (anchr unitigs step):
   pgr asm contig pe.cor.fa -o unitigs_K31.fasta --kmer 31

2. Assemble from paired-end reads (anchr 2_insert_size step):
   pgr asm contig R1.fq.gz R2.fq.gz -o contigs.fasta

3. Raise the minimum contig length:
   pgr asm contig in.fq -o out.fasta --min-contig-len 500

4. Raise the seeding depth threshold (tadpole `mincountseed`):
   pgr asm contig in.fq -o out.fasta --min-count-seed 5

5. Drop low-coverage contigs (tadpole `mincoverage`):
   pgr asm contig in.fq -o out.fasta --min-coverage 5
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
                .help("Minimum contig length (default: max(124, 2*k))"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum k-mer depth to seed a contig (tadpole mincountseed, default 3)"),
        )
        .arg(
            Arg::new("min_coverage")
                .long("min-coverage")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help(
                    "Minimum mean k-mer coverage for a contig (tadpole mincoverage, default 1.0)",
                ),
        )
        .arg(
            Arg::new("no_bubbles")
                .long("no-bubbles")
                .action(ArgAction::SetTrue)
                .help("Keep parallel-path contigs separate (disable bubble popping)"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Accepted for tadpole.sh compatibility; ignored (deterministic single-pass)"),
        )
}

/// Execute the contig command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    // Validate the thread-count value; processing stays deterministic
    // single-pass (see the design notes), so the result is not used.
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
        min_coverage: args.get_one::<f32>("min_coverage").copied().unwrap_or(1.0),
        pop_bubbles: !args.get_flag("no_bubbles"),
        ..AssembleOptions::default()
    };

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let stats = assemble(&infiles, &mut out, &opts)?;
    out.flush()?;
    eprintln!(
        "Reads in: {}  Contigs: {}  Bases: {}  Longest: {}",
        stats.reads_in, stats.contigs_built, stats.bases_built, stats.longest_contig
    );
    Ok(())
}

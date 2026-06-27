use clap::*;
use pgr::libs::paf::index::PafIndex;

pub fn make_subcommand() -> Command {
    Command::new("index")
        .about("Build interval-tree index from PAF files")
        .after_help(
            r###"
Builds a per-target interval-tree index from one or more PAF files.

The index can be used by `pgr paf query` for fast coordinate projection
and transitive closure traversal.

Notes:
* Input PAF files should contain `cg:Z:` tags for accurate coordinate projection
* Reads from stdin if input file is 'stdin'

Examples:
1. Index a single PAF file:
   pgr paf index alignments.paf

2. Index multiple PAF files:
   pgr paf index *.paf

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input PAF file(s) to index"),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    eprintln!("Building PAF index...");

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::reader(infile);
        let idx = PafIndex::build(reader)?;

        eprintln!("{}", infile);
        eprintln!("  sequences: {}", idx.names.len());
        eprintln!("  targets:   {}", idx.num_targets());
    }

    Ok(())
}

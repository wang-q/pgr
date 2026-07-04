use clap::{ArgMatches, Command};
use pgr::libs::paf::index::PafIndex;
/// Build the clap subcommand for index.
pub fn make_subcommand() -> Command {
    Command::new("index")
        .about("Builds interval-tree index from PAF files")
        .after_help(
            r###"
Builds a per-target interval-tree index from one or more PAF files.

Multiple input files are merged into a single unified index:
sequences with the same name across files share the same internal ID.

Use -o to persist the index to a .paf.idx file for reusable queries.

Notes:
* Input PAF files should contain `cg:Z:` tags for accurate coordinate projection
* Supports both plain text and gzipped (.gz) files (including BGZF)
* Reads from stdin if input file is 'stdin'

Examples:
1. Index a single PAF file and print summary:
   pgr paf index alignments.paf

2. Index and save for later queries:
   pgr paf index alignments.paf -o alignments.paf.idx

3. Merge multiple PAF files into one index:
   pgr paf index a.paf b.paf -o merged.paf.idx

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("PAF"))
        .arg(crate::cmd_pgr::args::outfile_arg_optional())
}
/// Execute the index command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let infiles: Vec<&String> = args.get_many::<String>("infiles").unwrap().collect();
    let count = infiles.len();

    log::info!("Building PAF index from {count} file(s)...");

    // For single-file input: use build_from_path (enables lazy CIGAR for BGZF).
    // For multi-file input: use build_multi (in-memory merge).
    let idx = if count == 1 {
        PafIndex::build_from_path(infiles[0])?
    } else {
        let readers: Vec<_> = infiles
            .iter()
            .map(|f| pgr::reader(f))
            .collect::<Result<Vec<_>, _>>()?;
        PafIndex::build_multi(readers)?
    };

    let lazy = idx.is_lazy();
    log::info!("  sequences: {}", idx.names.len());
    log::info!("  targets:   {}", idx.num_targets());
    if lazy {
        log::info!("  mode:      lazy (BGZF virtual-position CIGAR)");
    }

    if let Some(outfile) = args.get_one::<String>("outfile") {
        idx.save(outfile)?;
        log::info!("  saved to {}", outfile);
    }

    Ok(())
}

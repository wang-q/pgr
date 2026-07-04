use clap::{Arg, ArgMatches, Command};
use std::fmt::Write;

/// Build the clap subcommand for consensus.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        Command::new("consensus")
            .about("Generates consensus sequences using POA")
            .after_help(
                r###"
Generates consensus sequences using POA (Partial Order Alignment) graph.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* POA Engine:
    * `--engine builtin` (default): Uses built-in Rust implementation.
    * `--engine spoa`: Forces use of external `spoa` command.
* Alignment Parameters:
    * Configurable via `--match`, `--mismatch`, `--gap-open`, `--gap-extend`, `--align-mode`.
    * Defaults: Global alignment; Match 5, Mismatch -4, GapOpen -8, GapExtend -6.
* Supports parallel processing for improved performance
    * Running in parallel mode with 1 reader, 1 writer and the corresponding number of workers
    * The order of output may be different from the original
* If outgroups are present, they are handled appropriately

Examples:
1. Generate consensus sequences from a block FA file:
   pgr fas consensus tests/fas/example.fas

2. Generate consensus sequences using built-in engine:
   pgr fas consensus tests/fas/example.fas --engine builtin

3. Generate consensus sequences with outgroups:
   pgr fas consensus tests/fas/example.fas --outgroup

4. Run in parallel with 4 threads:
   pgr fas consensus tests/fas/example.fas --parallel 4

5. Output results to a file:
   pgr fas consensus tests/fas/example.fas -o output.fas

"###,
            )
            .arg(crate::cmd_pgr::args::engine_arg(
                &["builtin", "spoa"],
                "builtin",
                "POA engine to use",
            ))
            .arg(
                Arg::new("align_mode")
                    .long("align-mode")
                    .value_parser(["local", "global", "semi_global"])
                    .default_value("global") // Default to global for fas consensus
                    .help("Alignment mode"),
            )
            .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
            .arg(
                Arg::new("consensus_name")
                    .long("consensus-name")
                    .num_args(1)
                    .default_value("consensus")
                    .help("Name of the consensus"),
            )
            .arg(crate::cmd_pgr::args::outgroup_arg())
            .arg(crate::cmd_pgr::args::parallel_arg())
            .arg(crate::cmd_pgr::args::outfile_arg()),
        true,
    )
}

/// Execute the consensus command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    pgr::libs::fmt::fas::run_pipeline(&mut writer, &infiles, parallel, |block| {
        proc_block(block, args)
    })
}

fn proc_block(block: &pgr::libs::fmt::fas::FasBlock, args: &ArgMatches) -> anyhow::Result<String> {
    let cname = args.get_one::<String>("consensus_name").unwrap();
    let has_outgroup = args.get_flag("outgroup");

    let engine = args.get_one::<String>("engine").unwrap();

    let params = crate::cmd_pgr::args::get_poa_params(args);
    let algorithm = args.get_one::<String>("align_mode").unwrap();

    // Map algorithm string to integer code (0=local, 1=global, 2=semi_global) for internal use/spoa
    const ALGO_LOCAL: i32 = 0;
    const ALGO_GLOBAL: i32 = 1;
    const ALGO_SEMI_GLOBAL: i32 = 2;
    let algo_code = match algorithm.as_str() {
        "local" => ALGO_LOCAL,
        "global" => ALGO_GLOBAL,
        "semi_global" => ALGO_SEMI_GLOBAL,
        _ => anyhow::bail!("unknown align_mode: {}", algorithm),
    };

    let mut seqs = vec![];

    let outgroup = if has_outgroup {
        block.entries.iter().last()
    } else {
        None
    };

    for entry in &block.entries {
        seqs.push(entry.seq().as_ref());
    }
    if outgroup.is_some() {
        seqs.pop(); // Remove the outgroup sequence
    }

    // Generate consensus sequence
    let mut cons = match engine.as_str() {
        "spoa" => pgr::libs::alignment::get_consensus_poa_external(
            &seqs,
            params.match_score,
            params.mismatch_score,
            params.gap_open,
            params.gap_extend,
            algo_code,
        )?,
        _ => pgr::libs::alignment::get_consensus_poa_builtin(
            &seqs,
            params.match_score,
            params.mismatch_score,
            params.gap_open,
            params.gap_extend,
            algo_code,
        )?,
    };
    cons = cons.replace('-', "");

    let mut range = match block.entries.first() {
        Some(e) => e.range().clone(),
        None => anyhow::bail!("empty block"),
    };

    let mut out_string = String::new();
    if range.is_valid() {
        *range.name_mut() = cname.to_string();
        writeln!(out_string, ">{}\n{}", range, cons)?;
    } else {
        writeln!(out_string, ">{}\n{}", cname, cons)?;
    }
    if let Some(og) = outgroup {
        out_string.push_str(&og.to_string());
    }

    // end of a block
    out_string.push('\n');

    Ok(out_string)
}

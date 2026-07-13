use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::fmt::fas::{consensus_block, run_pipeline, ConsensusOptions};

/// Build the clap subcommand for consensus.
pub fn make_subcommand() -> Command {
    crate::cmd_pgr::args::add_poa_args(
        Command::new("consensus")
            .about("Generates consensus sequences using POA")
            .after_help(
                r###"
Generates consensus sequences using POA (Partial Order Alignment).

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* POA engine:
    * `--engine builtin` (default): built-in Rust implementation.
    * `--engine spoa`: external `spoa` command.
* Alignment parameters:
    * `--match` (default: 5), `--mismatch` (default: -4)
    * `--gap-open` (default: -8), `--gap-extend` (default: -6)
    * `--align-mode` (default: global)
* `--consensus-name` sets the output header name (default: consensus)
* `--outgroup` preserves the last sequence as outgroup (excluded from consensus)
* Parallel mode (`-p`) may change output order

Examples:
1. Generate consensus sequences from a block FA file:
   pgr fas consensus tests/fas/example.fas

2. Generate consensus with outgroup:
   pgr fas consensus tests/fas/example.fas --outgroup

3. Run in parallel with 4 threads:
   pgr fas consensus tests/fas/example.fas --parallel 4

4. Output results to a file:
   pgr fas consensus tests/fas/example.fas -o output.fas

"###,
            )
            .arg(crate::cmd_pgr::args::engine_arg(
                &["builtin", "spoa"],
                "builtin",
                "POA engine to use",
            ))
            .arg(crate::cmd_pgr::args::align_mode_arg())
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
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    // Map algorithm string to integer code (0=local, 1=global, 2=semi_global) for internal use/spoa
    let algo_code = crate::cmd_pgr::args::get_align_mode_code(args)?;

    let opts = ConsensusOptions {
        cname: args.get_one::<String>("consensus_name").unwrap().clone(),
        has_outgroup: args.get_flag("outgroup"),
        engine: args.get_one::<String>("engine").unwrap().clone(),
        params: crate::cmd_pgr::args::get_poa_params(args),
        algo_code,
    };

    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    run_pipeline(&mut writer, &infiles, parallel, |block| {
        consensus_block(block, &opts)
    })
}

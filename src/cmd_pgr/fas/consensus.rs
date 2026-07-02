use clap::*;

use super::common;

// Create clap subcommand arguments
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
            .arg(
                Arg::new("engine")
                    .long("engine")
                    .value_parser(["builtin", "spoa"])
                    .default_value("builtin")
                    .help("POA engine to use"),
            )
            .arg(
                Arg::new("align_mode")
                    .long("align-mode")
                    .short('l')
                    .value_parser(["local", "global", "semi_global"])
                    .default_value("global") // Default to global for fas consensus
                    .help("Alignment mode"),
            )
            .arg(
                Arg::new("infiles")
                    .required(true)
                    .num_args(1..)
                    .index(1)
                    .help("Input block FA file(s) to process"),
            )
            .arg(
                Arg::new("cname")
                    .long("consensus-name")
                    .num_args(1)
                    .default_value("consensus")
                    .help("Name of the consensus"),
            )
            .arg(
                Arg::new("outgroup")
                    .long("outgroup")
                    .action(ArgAction::SetTrue)
                    .help("Indicates the presence of outgroups at the end of each block"),
            )
            .arg(crate::cmd_pgr::args::parallel_arg())
            .arg(crate::cmd_pgr::args::outfile_arg()),
        true,
    )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    common::run_pipeline(args, parallel, |block| proc_block(block, args))
}

fn proc_block(block: &pgr::libs::fmt::fas::FasBlock, args: &ArgMatches) -> anyhow::Result<String> {
    //----------------------------
    // Args
    //----------------------------
    let cname = args.get_one::<String>("cname").unwrap();
    let has_outgroup = args.get_flag("outgroup");

    let engine = args.get_one::<String>("engine").unwrap();

    let params = crate::cmd_pgr::args::get_poa_params(args);
    let algorithm = args.get_one::<String>("align_mode").unwrap();

    // Map algorithm string to integer code (0=local, 1=global, 2=semi_global) for internal use/spoa
    let algo_code = match algorithm.as_str() {
        "local" => 0,
        "global" => 1,
        "semi_global" => 2,
        _ => 1,
    };

    //----------------------------
    // Ops
    //----------------------------
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

    //----------------------------
    // Output
    //----------------------------
    let mut out_string = "".to_string();
    if range.is_valid() {
        *range.name_mut() = cname.to_string();
        out_string += format!(">{}\n{}\n", range, cons).as_ref();
    } else {
        out_string += format!(">{}\n{}\n", cname, cons).as_ref();
    }
    if let Some(og) = outgroup {
        out_string += og.to_string().as_ref();
    }

    // end of a block
    out_string += "\n";

    Ok(out_string)
}

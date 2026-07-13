use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use pgr::libs::fmt::fas::{refine_block, run_pipeline, RefineOptions};

/// Build the clap subcommand for refine.
pub fn make_subcommand() -> Command {
    Command::new("refine")
        .about("Realigns files with built-in or external programs and trims unwanted regions")
        .after_help(
            r###"
Realigns sequences in block FA files using built-in or external programs and trims unwanted regions.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Supported MSA programs (`--engine`):
    * `builtin` (default): built-in Rust POA implementation.
    * `clustalw`, `mafft`, `muscle`, `spoa`: external commands.
    * `none`: skip realigning (useful for trimming only).
* `--chop` trims head/tail indels (default: 0, disabled)
* `--quick` aligns only indel-adjacent regions (useful for .axt/.maf conversions)
    * `--indel-pad` enlarges indel regions in quick mode (default: 50)
    * `--fill` fills holes between indels in quick mode (default: 50)
* Parallel mode (`-p`) may change output order

Examples:
1. Realign block FA files using builtin (default):
   pgr fas refine tests/fas/refine.fas tests/fas/refine2.fas

2. Realign using mafft with 4 threads:
   pgr fas refine tests/fas/refine.fas --engine mafft --parallel 4

3. Quick alignment for files converted from pairwise alignments:
   pgr fas refine tests/fas/refine.fas --quick --parallel 4

4. Output results to a file:
   pgr fas refine tests/fas/refine.fas -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::engine_arg(
            &["builtin", "clustalw", "mafft", "muscle", "spoa", "none"],
            "builtin",
            "Aligning program (builtin/clustalw/mafft/muscle/spoa/none)",
        ))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(
            Arg::new("chop")
                .long("chop")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Chop head and tail indels"),
        )
        .arg(
            Arg::new("is_quick")
                .long("quick")
                .action(ArgAction::SetTrue)
                .help("Quick mode, only aligns indel adjacent regions"),
        )
        .arg(
            Arg::new("indel_pad")
                .long("indel-pad")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("50")
                .help("In quick mode, enlarge indel regions"),
        )
        .arg(
            Arg::new("fill")
                .long("fill")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("50")
                .help("In quick mode, fill holes between indel"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the refine command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let mut writer =
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?;

    let opts = RefineOptions {
        engine: args.get_one::<String>("engine").unwrap(),
        has_outgroup: args.get_flag("outgroup"),
        chop: *args.get_one::<usize>("chop").unwrap(),
        is_quick: args.get_flag("is_quick"),
        pad: *args.get_one::<usize>("indel_pad").unwrap(),
        fill: *args.get_one::<usize>("fill").unwrap(),
    };

    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    run_pipeline(&mut writer, &infiles, parallel, |block| {
        refine_block(block, &opts)
    })
}

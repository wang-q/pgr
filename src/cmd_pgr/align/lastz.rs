use clap::builder::PossibleValuesParser;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for lastz.
pub fn make_subcommand() -> Command {
    Command::new("lastz")
        .about("Wraps lastz alignment (Cactus style)")
        .after_help(format!(
            r###"
This command wraps lastz to perform alignments suitable for the Cactus RepeatMasking workflow.

It handles:
*   Parallel execution for multiple target files.
*   Directory recursion for target and query inputs.
*   Adding required modifiers: [multiple][nameparse=darkspace].
*   Setting output format to LAV.
*   Setting query depth threshold: --querydepth=keep,nowarn:N.
    N is the depth of coverage threshold (aligned bases / query length).
    When the threshold is exceeded, processing stops for that query/strand to save time.
    'keep' ensures alignments found *before* the threshold are reported (unlike default which discards all).
    'nowarn' suppresses warnings about exceeded depth.
    Note: Reported alignments are the first found, not necessarily optimal.
    Default depth 50 allows ~50x coverage.

{}
Examples:
1. Single target with set01:
   pgr align lastz target.fa query.fa --preset set01 -o lastz_out

2. Directory inputs:
   pgr align lastz target_dir/ query_dir/ --preset set03 -o lastz_out

3. Show parameters and matrix for set01:
   pgr align lastz --preset set01 --show-preset

"###,
            pgr::libs::lastz::preset_help()
        ))
        .arg(
            Arg::new("target")
                .required_unless_present("show_preset")
                .index(1)
                .help("Target FASTA file or directory"),
        )
        .arg(
            Arg::new("query")
                .index(2)
                .help("Query FASTA file or directory; omit for self-alignment"),
        )
        .arg(
            Arg::new("query_depth")
                .long("query-depth")
                .default_value("50")
                .value_parser(value_parser!(usize))
                .help("Query depth threshold"),
        )
        .arg(
            Arg::new("is_self")
                .long("self")
                .action(clap::ArgAction::SetTrue)
                .help("Self-alignment (query omitted or the same input as the target)"),
        )
        .arg(
            Arg::new("preset")
                .long("preset")
                .value_parser(PossibleValuesParser::new(pgr::libs::lastz::preset_names()))
                .help("Use predefined parameter sets (set01..set07)"),
        )
        .arg(
            Arg::new("show_preset")
                .long("show-preset")
                .action(ArgAction::SetTrue)
                .help("Display the configuration (parameters & matrix) for the selected preset and exit"),
        )
        .arg(
            Arg::new("lastz_args")
                .long("lastz-args")
                .help("Additional arguments passed directly to lastz (overrides preset)"),
        )
        .arg(crate::cmd_pgr::args::outdir_arg_with_default("lastz_out"))
        .arg(crate::cmd_pgr::args::parallel_arg_with_default("4"))
}
/// Execute the lastz command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let preset = args.get_one::<String>("preset");

    // Show preset and exit
    if args.get_flag("show_preset") {
        let preset_name = preset
            .ok_or_else(|| anyhow::anyhow!("--show-preset requires --preset to be specified."))?;
        let p = pgr::libs::lastz::find_preset(preset_name)
            .ok_or_else(|| anyhow::anyhow!("unknown preset: {}", preset_name))?;
        let stdout = std::io::stdout();
        let mut out = stdout.lock();
        writeln!(out, "Preset: {}", p.name)?;
        writeln!(out, "Description: {}", p.desc)?;
        writeln!(out, "Parameters: {}", p.params)?;
        if let Some(matrix) = p.matrix {
            writeln!(out, "\nMatrix Content:\n{}", matrix)?;
        }
        out.flush()?;
        return Ok(());
    }

    let arg_query = args.get_one::<String>("query");
    let arg_target = args.get_one::<String>("target").unwrap();
    let opt_depth = *args.get_one::<usize>("query_depth").unwrap();
    let opt_lastz_args = args.get_one::<String>("lastz_args");
    let opt_output = args.get_one::<String>("outdir").unwrap();
    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();
    let is_self = args.get_flag("is_self");
    if is_self {
        if let Some(q) = arg_query {
            anyhow::ensure!(
                q == arg_target,
                "--self expects the query to be the same input as the target \
                 (omit the query or pass the same path)"
            );
        }
    }
    let self_mode = is_self || arg_query.is_none();

    // Check if lastz is installed
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }

    let mut target_files = pgr::libs::fmt::fa::find_fasta_files(arg_target);
    target_files.sort();

    if target_files.is_empty() {
        anyhow::bail!("No target FASTA files found in {}", arg_target);
    }
    let query_files = if self_mode {
        target_files.clone()
    } else {
        let mut qf = pgr::libs::fmt::fa::find_fasta_files(arg_query.unwrap());
        qf.sort();
        qf
    };
    if query_files.is_empty() {
        anyhow::bail!(
            "No query FASTA files found in {}",
            if self_mode {
                arg_target
            } else {
                arg_query.unwrap()
            }
        );
    }

    // Common lastz arguments (query-depth, LAV format, preset params + matrix)
    // come from the shared builder; user overrides are appended afterwards.
    // The `_` prefix keeps the matrix temp file alive until lastz finishes.
    let (mut common_args, _temp_matrix_handle) =
        pgr::libs::lastz::build_common_args(preset.map(|s| s.as_str()), opt_depth)?;

    if let Some(args) = opt_lastz_args {
        for arg in args.split_whitespace() {
            common_args.push(arg.to_string());
        }
    }

    // Delegate the parallel orchestration to libs::lastz
    let opts = pgr::libs::lastz::RunLastzOptions {
        depth: opt_depth,
        is_self: self_mode,
        common_args,
        output_dir: opt_output.clone(),
        parallel: opt_parallel,
    };

    pgr::libs::lastz::run_lastz(target_files, query_files, opts)
}

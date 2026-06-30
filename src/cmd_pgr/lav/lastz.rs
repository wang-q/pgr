use clap::builder::PossibleValuesParser;
use clap::*;
use std::io::Write;
use tempfile::NamedTempFile;

// TODO: [multiple] on target
// TODO: unmask on t/q

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("lastz")
        .about("Wrapper for lastz alignment (Cactus style)")
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
    # Single target with set01
    pgr lav lastz target.fa query.fa --preset set01 -o lastz_out

    # Directory inputs
    pgr lav lastz target_dir/ query_dir/ --preset set03 -o lastz_out

    # Show parameters and matrix for set01
    pgr lav lastz --preset set01 --show-preset

"###,
            pgr::libs::lastz::preset_help()
        ))
        .arg(
            Arg::new("target")
                .required(true)
                .index(1)
                .help("Target FASTA file or directory"),
        )
        .arg(
            Arg::new("query")
                .required(true)
                .index(2)
                .help("Query FASTA file or directory"),
        )
        .arg(
            Arg::new("depth")
                .long("depth")
                .default_value("50")
                .value_parser(value_parser!(usize))
                .help("Query depth threshold"),
        )
        .arg(
            Arg::new("is_self")
                .long("self")
                .action(clap::ArgAction::SetTrue)
                .help("Self-alignment"),
        )
        .arg(
            Arg::new("preset")
                .long("preset")
                .short('s')
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
        .arg(crate::cmd_pgr::args::outdir_arg().default_value("lastz_out"))
        .arg(crate::cmd_pgr::args::parallel_arg().default_value("4"))
}

pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    let preset = matches.get_one::<String>("preset");

    // Show preset and exit
    if matches.get_flag("show_preset") {
        let preset_name = preset
            .ok_or_else(|| anyhow::anyhow!("--show-preset requires --preset to be specified."))?;
        let p = pgr::libs::lastz::find_preset(preset_name)
            .ok_or_else(|| anyhow::anyhow!("unknown preset: {}", preset_name))?;
        println!("Preset: {}", p.name);
        println!("Description: {}", p.desc);
        println!("Parameters: {}", p.params);
        if let Some(matrix) = p.matrix {
            println!("\nMatrix Content:\n{}", matrix);
        }
        return Ok(());
    }

    let arg_query = matches.get_one::<String>("query").unwrap();
    let arg_target = matches.get_one::<String>("target").unwrap();
    let opt_depth = *matches.get_one::<usize>("depth").unwrap();
    let opt_lastz_args = matches.get_one::<String>("lastz_args");
    let opt_output = matches.get_one::<String>("outdir").unwrap();
    let opt_parallel = *matches.get_one::<usize>("parallel").unwrap();
    let is_self = matches.get_flag("is_self");

    // Check if lastz is installed
    if which::which("lastz").is_err() {
        anyhow::bail!("lastz not found in PATH. Please install lastz first.");
    }

    // Expand files
    let mut query_files = pgr::libs::io::find_fasta_files(arg_query);
    query_files.sort();

    let mut target_files = pgr::libs::io::find_fasta_files(arg_target);
    target_files.sort();

    if query_files.is_empty() {
        anyhow::bail!("No query FASTA files found in {}", arg_query);
    }
    if target_files.is_empty() {
        anyhow::bail!("No target FASTA files found in {}", arg_target);
    }

    // Prepare matrix temp file if preset defines one (keep handle alive for the run)
    let mut _temp_matrix_handle: Option<NamedTempFile> = None;
    let mut matrix_path = String::new();

    if let Some(preset_name) = preset {
        let p = pgr::libs::lastz::find_preset(preset_name)
            .ok_or_else(|| anyhow::anyhow!("unknown preset: {}", preset_name))?;
        if let Some(matrix) = p.matrix {
            let mut t = NamedTempFile::new()?;
            t.write_all(matrix.as_bytes())?;
            matrix_path = t.path().to_string_lossy().to_string();
            _temp_matrix_handle = Some(t);
        }
    }

    // Build common args (depth, format, preset params, user args)
    let mut common_args = Vec::new();
    common_args.push(format!("--querydepth=keep,nowarn:{}", opt_depth));
    common_args.push("--format=lav".to_string());
    common_args.push("--markend".to_string());
    common_args.push("--ambiguous=iupac".to_string());

    if let Some(preset_name) = preset {
        let p = pgr::libs::lastz::find_preset(preset_name)
            .ok_or_else(|| anyhow::anyhow!("unknown preset: {}", preset_name))?;
        for arg in p.params.split_whitespace() {
            if !arg.starts_with("Q=") {
                common_args.push(arg.to_string());
            }
        }
        if !matrix_path.is_empty() {
            common_args.push(format!("Q={}", matrix_path));
        }
    }

    if let Some(args) = opt_lastz_args {
        for arg in args.split_whitespace() {
            common_args.push(arg.to_string());
        }
    }

    // Delegate the parallel orchestration to libs::lastz
    let opts = pgr::libs::lastz::RunLastzOptions {
        depth: opt_depth,
        is_self,
        common_args,
        output_dir: opt_output.clone(),
        parallel: opt_parallel,
    };

    pgr::libs::lastz::run_lastz(target_files, query_files, opts)
}

use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use std::io::Write;
/// Build the clap subcommand for to-dna.
pub fn make_subcommand() -> Command {
    Command::new("to-dna")
        .about("Converts ms output haplotypes (0/1) to DNA sequences (FASTA)")
        .arg(
            Arg::new("gc")
                .long("gc")
                .short('g')
                .num_args(1)
                .default_value("0.5")
                .value_parser(value_parser!(f64))
                .help("GC content ratio in ancestral sequence (0..1)"),
        )
        .arg(crate::cmd_pgr::args::seed_arg(
            None,
            Some('s'),
            "Random seed; default uses system time and PID",
        ))
        .arg(
            Arg::new("no_perturb")
                .long("no-perturb")
                .action(clap::ArgAction::SetTrue)
                .help("Disable positions micro-perturbation"),
        )
        .arg(
            Arg::new("verbose")
                .long("verbose")
                .short('v')
                .action(clap::ArgAction::SetTrue)
                .help("Print runtime information (paths and inputs)"),
        )
        .arg(
            Arg::new("doc")
                .long("doc")
                .action(clap::ArgAction::SetTrue)
                .help("Print the full documentation (markdown)"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
        .arg(
            Arg::new("infiles")
                .num_args(0..)
                .help("Input files with ms output; reads stdin when omitted"),
        )
        .after_help(
            r###"
Examples:
1. Pipe ms output to pgr ms to-dna:
   ms 10 1 -t 5 -r 0 1000 | pgr ms to-dna --gc 0.5 > out.fa

2. Read from file and write to output:
   pgr ms to-dna input.ms -o out.fa --seed 12345

3. Disable position perturbation (keep original ms positions):
   pgr ms to-dna input.ms --no-perturb

Output Format:
  FASTA format with single-line sequences.
  Headers: >[Lx_][Px_]Sx
    Lx: Batch/Replicate index (if multiple)
    Px: Population index (if multiple)
    Sx: Sample index
"###,
        )
}

/// Execute the to-dna command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    if args.get_flag("doc") {
        let stdout = std::io::stdout();
        let mut out = stdout.lock();
        writeln!(out, "{}", include_str!("../../../docs/ms-to-dna.md"))?;
        out.flush()?;
        return Ok(());
    }
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let gc = *args.get_one::<f64>("gc").unwrap();
    anyhow::ensure!(
        (0.0..=1.0).contains(&gc),
        "--gc must be in [0, 1], got {}",
        gc
    );
    let seed = args.get_one::<u64>("seed").copied();
    let no_perturb = args.get_flag("no_perturb");
    let verbose = args.get_flag("verbose");

    if verbose {
        let curdir = std::env::current_dir()?;
        let pgr = pgr::libs::io::current_exe_string()?;
        eprintln!("==> Paths");
        eprintln!("    \"pgr\"     = {}", pgr);
        eprintln!("    \"curdir\"  = {:?}", curdir);
    }

    if verbose {
        eprintln!("==> Inputs");
    }
    let files: Vec<String> = args
        .get_many::<String>("infiles")
        .map(|vals| vals.map(|s| s.to_string()).collect())
        .unwrap_or_default();
    let abs_files: Vec<String> = files
        .iter()
        .map(|f| intspan::absolute_path(f).map(|p| p.display().to_string()))
        .collect::<Result<_, _>>()?;
    if verbose {
        if abs_files.is_empty() {
            eprintln!("    [stdin]");
        } else {
            eprintln!("    files = {:?}", abs_files);
        }
    }

    let seed_final = seed.unwrap_or(pgr::libs::ms::system_seed());
    if verbose {
        eprintln!("==> Seed");
        eprintln!("    using = {}", seed_final);
    }

    // Writer
    let mut writer: Box<dyn Write> = Box::new(
        pgr::writer(outfile).with_context(|| format!("Failed to open writer for {}", outfile))?,
    );

    // Process inputs (stdin or files)
    if abs_files.is_empty() {
        pgr::libs::ms::convert_stream(
            pgr::reader("stdin").with_context(|| "Failed to open reader for stdin")?,
            gc,
            Some(seed_final),
            &mut writer,
            no_perturb,
        )?;
    } else {
        for path in abs_files {
            pgr::libs::ms::convert_stream(
                pgr::reader(&path)
                    .with_context(|| format!("Failed to open reader for {}", path))?,
                gc,
                Some(seed_final),
                &mut writer,
                no_perturb,
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

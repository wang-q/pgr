use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};

/// Build the clap subcommand for to-xlsx.
pub fn make_subcommand() -> Command {
    Command::new("to-xlsx")
        .about("Exports variations (substitutions/indels) to Excel")
        .after_help(
            r###"
Exports variations (substitutions/indels) to Excel.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'

Examples:
1. Export variations to an Excel file:
   pgr fas to-xlsx tests/fas/part1.fas -o variations.xlsx

2. Include indels and handle outgroup:
   pgr fas to-xlsx tests/fas/part1.fas --indel --outgroup

3. Filter variations by frequency (e.g., min 0.1, max 0.9):
   pgr fas to-xlsx tests/fas/part1.fas --min-freq 0.1 --max-freq 0.9

4. Omit singleton and complex variations:
   pgr fas to-xlsx tests/fas/part1.fas --no-single --no-complex

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(
            Arg::new("wrap")
                .long("wrap")
                .value_parser(value_parser!(u16))
                .num_args(1)
                .default_value("50")
                .help("Wrap length"),
        )
        .arg(
            Arg::new("indel")
                .long("indel")
                .action(ArgAction::SetTrue)
                .help("List indels"),
        )
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(
            Arg::new("no_single")
                .long("no-single")
                .action(ArgAction::SetTrue)
                .help("Omit singleton SNPs and indels"),
        )
        .arg(
            Arg::new("no_complex")
                .long("no-complex")
                .action(ArgAction::SetTrue)
                .help("Omit complex SNPs and indels"),
        )
        .arg(
            Arg::new("min_freq")
                .long("min-freq")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Minimal frequency"),
        )
        .arg(
            Arg::new("max_freq")
                .long("max-freq")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Maximal frequency"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg_with_default(
            "variations.xlsx",
        ))
}

/// Execute the to-xlsx command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let opt_wrap = *args.get_one::<u16>("wrap").unwrap();
    let is_indel = args.get_flag("indel");
    let is_outgroup = args.get_flag("outgroup");
    let is_nosingle = args.get_flag("no_single");
    let is_nocomplex = args.get_flag("no_complex");
    let opt_min = args.get_one::<f64>("min_freq").cloned();
    let opt_max = args.get_one::<f64>("max_freq").cloned();
    if let Some(v) = opt_min {
        anyhow::ensure!(
            v.is_finite() && (0.0..=1.0).contains(&v),
            "--min-freq must be in [0, 1]: {}",
            v
        );
    }
    if let Some(v) = opt_max {
        anyhow::ensure!(
            v.is_finite() && (0.0..=1.0).contains(&v),
            "--max-freq must be in [0, 1]: {}",
            v
        );
    }
    if let (Some(a), Some(b)) = (opt_min, opt_max) {
        anyhow::ensure!(a <= b, "--min-freq ({}) must be <= --max-freq ({})", a, b);
    }
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();

    pgr::libs::fas_xlsx::export_to_xlsx(
        &infiles,
        outfile,
        opt_wrap,
        is_indel,
        is_outgroup,
        is_nosingle,
        is_nocomplex,
        opt_min,
        opt_max,
    )
}

use clap::*;

// Create clap subcommand arguments
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
   pgr fas to-xlsx tests/fas/part1.fas --min 0.1 --max 0.9

4. Omit singleton and complex variations:
   pgr fas to-xlsx tests/fas/part1.fas --nosingle --nocomplex

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input block FA file(s) to process"),
        )
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
        .arg(
            Arg::new("outgroup")
                .long("outgroup")
                .action(ArgAction::SetTrue)
                .help("Indicates the presence of outgroups at the end of each block"),
        )
        .arg(
            Arg::new("nosingle")
                .long("nosingle")
                .action(ArgAction::SetTrue)
                .help("Omit singleton SNPs and indels"),
        )
        .arg(
            Arg::new("nocomplex")
                .long("nocomplex")
                .action(ArgAction::SetTrue)
                .help("Omit complex SNPs and indels"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Minimal frequency"),
        )
        .arg(
            Arg::new("max")
                .long("max")
                .value_parser(value_parser!(f64))
                .num_args(1)
                .help("Maximal frequency"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("variations.xlsx")
                .help("Output filename"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = args.get_one::<String>("outfile").unwrap();
    let opt_wrap = *args.get_one::<u16>("wrap").unwrap();
    let is_indel = args.get_flag("indel");
    let is_outgroup = args.get_flag("outgroup");
    let is_nosingle = args.get_flag("nosingle");
    let is_nocomplex = args.get_flag("nocomplex");
    let opt_min = args.get_one::<f64>("min").cloned();
    let opt_max = args.get_one::<f64>("max").cloned();
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

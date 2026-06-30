use clap::*;
use std::io::Write;

pub fn make_subcommand() -> Command {
    Command::new("transform")
        .about("Apply mathematical transformations to a matrix")
        .after_help(
            r###"
Transform matrix values element-wise.
Useful for converting similarity matrices to distance matrices.

Operations:
    * linear:     val = val * scale + offset
    * inv-linear: val = max - val
    * log:        val = -ln(val)
    * exp:        val = exp(-val)
    * square:     val = val * val
    * sqrt:       val = sqrt(val)

Normalization:
    If --normalize is set, values are normalized using diagonal elements before transformation:
    x_norm(i, j) = x(i, j) / sqrt(x(i, i) * x(j, j))

Examples:
    1. Convert Identity (0-100) to Distance (0-1):
       # Using linear: -0.01 * x + 1.0 = (100 - x) / 100
       pgr mat transform in.phy --op linear --scale -0.01 --offset 1.0

    2. Convert Identity (0-100) to Distance (0-100):
       pgr mat transform in.phy --op inv-linear --max 100

    3. Convert Similarity (0-1) to Distance (0-1):
       pgr mat transform in.phy --op inv-linear --max 1.0

    4. Log transformation with normalization (e.g. for probability):
       pgr mat transform in.phy --op log --normalize
"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .index(1)
                .help("Input PHYLIP matrix or pairwise TSV file"),
        )
        .arg(
            Arg::new("op")
                .long("op")
                .default_value("linear")
                .value_parser([
                    builder::PossibleValue::new("linear"),
                    builder::PossibleValue::new("inv-linear"),
                    builder::PossibleValue::new("log"),
                    builder::PossibleValue::new("exp"),
                    builder::PossibleValue::new("square"),
                    builder::PossibleValue::new("sqrt"),
                ])
                .help("Transformation operation"),
        )
        .arg(
            Arg::new("max")
                .long("max")
                .default_value("1.0")
                .value_parser(value_parser!(f32))
                .help("Maximum value for inv-linear"),
        )
        .arg(
            Arg::new("scale")
                .long("scale")
                .default_value("1.0")
                .value_parser(value_parser!(f32))
                .help("Scale factor for linear"),
        )
        .arg(
            Arg::new("offset")
                .long("offset")
                .default_value("0.0")
                .value_parser(value_parser!(f32))
                .help("Offset value for linear"),
        )
        .arg(
            Arg::new("normalize")
                .long("normalize")
                .action(ArgAction::SetTrue)
                .help("Normalize based on diagonal values"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .default_value("phylip")
                .value_parser([
                    builder::PossibleValue::new("phylip"),
                    builder::PossibleValue::new("pair"),
                ])
                .help("Input format"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let infile = args.get_one::<String>("infile").unwrap();
    let op = args.get_one::<String>("op").unwrap().as_str();
    let max_val = *args.get_one::<f32>("max").unwrap();
    let scale = *args.get_one::<f32>("scale").unwrap();
    let offset = *args.get_one::<f32>("offset").unwrap();
    let normalize = args.get_flag("normalize");
    let format = args.get_one::<String>("format").unwrap().as_str();
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap())?;

    //----------------------------
    // Load and Transform
    //----------------------------
    let matrix = if format == "pair" {
        pgr::libs::pairmat::NamedMatrix::from_pair_scores(infile, 0.0, 1.0)?
    } else {
        pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)?
    };

    let matrix =
        pgr::libs::pairmat::transform_matrix(&matrix, op, max_val, scale, offset, normalize)?;

    //----------------------------
    // Output
    //----------------------------
    let size = matrix.size();
    let names = matrix.get_names();

    writer.write_fmt(format_args!("{:>4}\n", size))?;
    for (i, name) in names.iter().enumerate() {
        writer.write_fmt(format_args!("{}", name))?;
        for j in 0..size {
            let val = matrix.get(i, j);
            writer.write_fmt(format_args!("\t{:.6}", val))?;
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(())
}

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
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
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
    let mut writer = pgr::writer(args.get_one::<String>("outfile").unwrap());

    //----------------------------
    // Load and Process
    //----------------------------
    let mut matrix = if format == "pair" {
        pgr::libs::pairmat::NamedMatrix::from_pair_scores(infile, 0.0, 1.0)
    } else {
        pgr::libs::pairmat::NamedMatrix::from_relaxed_phylip(infile)
    };

    // Get diagonals if needed
    // We always try to get diagonals because transform operations might apply to diagonal elements too
    // (e.g. inv-linear: max - x, if x=1.0 then 0.0)
    // NamedMatrix stores diagonals internally now.
    let diags = matrix.get_diags().cloned().unwrap_or_default();
    let size = matrix.size();
    // We clone names here to avoid borrowing matrix while we mutate it later
    let names: Vec<String> = matrix.get_names().iter().map(|n| n.to_string()).collect();

    // Warn if normalize is requested but diagonals are missing (all zero)
    if normalize {
        if diags.is_empty() {
            eprintln!("Warning: --normalize requested but no diagonal values found. Result will be Inf/NaN.");
        } else {
            let max_diag = diags.iter().fold(0.0f32, |a, &b| a.max(b));
            if max_diag == 0.0 {
                eprintln!("Warning: --normalize requested but all diagonal values are 0.0. Result will be Inf/NaN.");
            }
        }
    }

    // Transform values
    // NamedMatrix stores only upper triangle (or lower, logically symmetric).
    // We iterate i, j where i < j (upper triangle).
    // CondensedMatrix::set handles i != j.
    for i in 0..size {
        for j in (i + 1)..size {
            let mut val = matrix.get(i, j);

            // 1. Normalize
            if normalize {
                let d_i = diags[i];
                let d_j = diags[j];
                // Avoid division by zero
                if d_i > 1e-9 && d_j > 1e-9 {
                    val = val / (d_i * d_j).sqrt();
                } else {
                    // If diagonal is 0, similarity is undefined or 0?
                    // Assuming 0 if undefined.
                    val = 0.0;
                }
            }

            // 2. Transform
            val = match op {
                "linear" => val * scale + offset,
                "inv-linear" => max_val - val,
                "log" => {
                    if val > 0.0 {
                        -val.ln()
                    } else {
                        // -ln(0) = Inf. Use a large number?
                        // Or just let it be Inf?
                        // f32::INFINITY
                        1000.0 // Cap at reasonable max distance?
                    }
                }
                "exp" => (-val).exp(),
                "square" => val * val,
                "sqrt" => {
                    if val >= 0.0 {
                        val.sqrt()
                    } else {
                        0.0
                    }
                }
                _ => val,
            };

            matrix.set(i, j, val);
        }
    }

    //----------------------------
    // Output
    //----------------------------
    // We always output full matrix for now (to match other commands default behavior)
    // Or should we support --mode?
    // Let's stick to full matrix tab-separated, similar to `mat format`.

    writer.write_fmt(format_args!("{:>4}\n", size))?;
    for i in 0..size {
        writer.write_fmt(format_args!("{}", names[i]))?;
        for j in 0..size {
            let val = if i == j {
                // Handle diagonal specially
                // We need to calculate the transformed diagonal value here because
                // the main transformation loop only updated off-diagonal elements.
                let mut d = if !diags.is_empty() { diags[i] } else { 0.0 };
                if normalize {
                    // x_norm(i,i) = x(i,i) / sqrt(x(i,i)*x(i,i)) = 1.0
                    if d > 1e-9 {
                        d = 1.0;
                    } else {
                        d = 0.0;
                    }
                }

                match op {
                    "linear" => d * scale + offset,
                    "inv-linear" => max_val - d,
                    "log" => {
                        if d > 0.0 {
                            -d.ln()
                        } else {
                            0.0
                        }
                    } // Distance diag usually 0
                    "exp" => (-d).exp(),
                    "square" => d * d,
                    "sqrt" => {
                        if d >= 0.0 {
                            d.sqrt()
                        } else {
                            0.0
                        }
                    }
                    _ => d,
                }
            } else {
                matrix.get(i, j)
            };

            writer.write_fmt(format_args!("\t{:.6}", val))?;
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(())
}

use clap::{Arg, ArgMatches, Command};
use pgr::libs::chain::net::read_nets;
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::Path;

pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Split a net file into one file per chromosome")
        .arg(
            Arg::new("input")
                .help("Input net file")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("output_dir")
                .help("Output directory")
                .required(true)
                .index(2),
        )
}

pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("input").unwrap();
    let output_dir = args.get_one::<String>("output_dir").unwrap();

    let reader = pgr::reader(input_path)?;

    let chroms = read_nets(reader)?;

    fs::create_dir_all(output_dir)?;

    for chrom in chroms {
        let file_path = Path::new(output_dir).join(format!("{}.net", chrom.name));
        let mut file = BufWriter::new(File::create(file_path)?);
        chrom.write(&mut file)?;
    }

    Ok(())
}

use clap::{ArgMatches, Command};
use pgr::libs::chain::net::read_nets;
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::Path;
/// Build the clap subcommand for split.
pub fn make_subcommand() -> Command {
    Command::new("split")
        .about("Splits a net file into one file per chromosome")
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input net file",
        ))
        .arg(crate::cmd_pgr::args::outdir_arg_required())
}
/// Execute the split command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input_path = args.get_one::<String>("infile").unwrap();
    let output_dir = args.get_one::<String>("outdir").unwrap();

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

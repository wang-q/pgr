use anyhow::Context;
use clap::{ArgMatches, Command};
/// Build the clap subcommand for to-psl.
pub fn make_subcommand() -> Command {
    Command::new("to-psl")
        .about("Converts from axt to psl format")
        .after_help(
            r###"
Where tSizes and qSizes are tab-delimited files with <seqName> <size> columns.

Examples:
1. Convert axt to psl:
   pgr axt to-psl in.axt -t t.sizes -q q.sizes -o out.psl
"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg().help("Input AXT file. [stdin] for standard input"))
        .arg(crate::cmd_pgr::args::t_sizes_arg().required(true))
        .arg(crate::cmd_pgr::args::q_sizes_arg().required(true))
        .arg(crate::cmd_pgr::args::outfile_arg())
}
/// Execute the to-psl command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let input = crate::cmd_pgr::args::get_infile(args);
    let t_sizes_path = args.get_one::<String>("t_sizes").unwrap();
    let q_sizes_path = args.get_one::<String>("q_sizes").unwrap();
    let output = crate::cmd_pgr::args::get_outfile(args);

    let t_sizes = pgr::read_sizes::<usize>(t_sizes_path)?;
    let q_sizes = pgr::read_sizes::<usize>(q_sizes_path)?;

    let reader =
        pgr::reader(input).with_context(|| format!("Failed to open reader for {}", input))?;
    let mut writer =
        pgr::writer(output).with_context(|| format!("Failed to open writer for {}", output))?;

    pgr::libs::fmt::axt::axt_to_psl(reader, &mut writer, &t_sizes, &q_sizes)?;

    Ok(())
}

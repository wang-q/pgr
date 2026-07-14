pub mod masked;
pub mod range;
pub mod size;
pub mod some;
pub mod to_fa;

use clap::{ArgMatches, Command};
/// Build the clap subcommand for 2bit.
pub fn make_subcommand() -> Command {
    Command::new("2bit")
        .about("Manages 2bit files")
        .after_help(
            r###"Description:
    Manipulates 2bit binary sequence files.

Subcommand groups:
* Info: masked / size
* Subset: range / some
* Transform: to-fa

Notes:
* 2bit files are binary and require random access (seeking)
* Does not support stdin or gzipped inputs

Examples:
1. Convert FASTA to 2bit:
   pgr fa to-2bit in.fa -o out.2bit

2. Extract masked regions:
   pgr 2bit masked out.2bit -o masked.txt

3. Convert back to FASTA:
   pgr 2bit to-fa out.2bit -o out.fa
"###,
        )
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(masked::make_subcommand())
        .subcommand(range::make_subcommand())
        .subcommand(size::make_subcommand())
        .subcommand(some::make_subcommand())
        .subcommand(to_fa::make_subcommand())
}
/// Execute the 2bit command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("masked", sub_matches)) => masked::execute(sub_matches),
        Some(("range", sub_matches)) => range::execute(sub_matches),
        Some(("size", sub_matches)) => size::execute(sub_matches),
        Some(("some", sub_matches)) => some::execute(sub_matches),
        Some(("to-fa", sub_matches)) => to_fa::execute(sub_matches),
        _ => Ok(()),
    }
}

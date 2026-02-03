extern crate clap;
use clap::*;

mod cmd_fasr;

fn main() -> anyhow::Result<()> {
    let app = Command::new("fasr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("`fasr` operates block fasta files")
        .propagate_version(true)
        .arg_required_else_help(true)
        .color(ColorChoice::Auto)
        .subcommand(cmd_fasr::create::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* info: check
* creation: create
* records: (none)
* transform: (none)
* transmute: (none)

"###,
        );

    // Check which subcommand the user ran...
    match app.get_matches().subcommand() {
        // info
        // creation
        Some(("create", sub_matches)) => cmd_fasr::create::execute(sub_matches),
        // records
        // transform
        // transmute
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: fasr variation --indel
// TODO: fasr match
//  sparsemem -maxmatch -F -l %d -b -n -k 4 -threads 4 %s %s > %s
//  mummer -maxmatch -F -l %d -b -n %s %s > %s
//  $exe, $length, $genome, $query, $result

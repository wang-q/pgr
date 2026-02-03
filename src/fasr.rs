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
        .subcommand(cmd_fasr::check::make_subcommand())
        .subcommand(cmd_fasr::create::make_subcommand())
        .subcommand(cmd_fasr::filter::make_subcommand())
        .subcommand(cmd_fasr::pl_p2m::make_subcommand())
        .subcommand(cmd_fasr::stat::make_subcommand())
        .subcommand(cmd_fasr::variation::make_subcommand())
        .subcommand(cmd_fasr::vcf::make_subcommand())
        .subcommand(cmd_fasr::xlsx::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* info: check / stat
* creation: create
* records: (none)
* transform: filter
* transmute: pl-p2m
* variations: variation / vcf / xlsx

"###,
        );

    // Check which subcommand the user ran...
    match app.get_matches().subcommand() {
        // info
        Some(("check", sub_matches)) => cmd_fasr::check::execute(sub_matches),
        Some(("stat", sub_matches)) => cmd_fasr::stat::execute(sub_matches),
        // creation
        Some(("create", sub_matches)) => cmd_fasr::create::execute(sub_matches),
        // records
        // transform
        Some(("filter", sub_matches)) => cmd_fasr::filter::execute(sub_matches),
        // transmute
        Some(("pl-p2m", sub_matches)) => cmd_fasr::pl_p2m::execute(sub_matches),
        // variations
        Some(("variation", sub_matches)) => cmd_fasr::variation::execute(sub_matches),
        Some(("vcf", sub_matches)) => cmd_fasr::vcf::execute(sub_matches),
        Some(("xlsx", sub_matches)) => cmd_fasr::xlsx::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: fasr variation --indel
// TODO: fasr match
//  sparsemem -maxmatch -F -l %d -b -n -k 4 -threads 4 %s %s > %s
//  mummer -maxmatch -F -l %d -b -n %s %s > %s
//  $exe, $length, $genome, $query, $result

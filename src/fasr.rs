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
        .subcommand(cmd_fasr::axt2fas::make_subcommand())
        .subcommand(cmd_fasr::check::make_subcommand())
        .subcommand(cmd_fasr::concat::make_subcommand())
        .subcommand(cmd_fasr::consensus::make_subcommand())
        .subcommand(cmd_fasr::cover::make_subcommand())
        .subcommand(cmd_fasr::create::make_subcommand())
        .subcommand(cmd_fasr::filter::make_subcommand())
        .subcommand(cmd_fasr::join::make_subcommand())
        .subcommand(cmd_fasr::link::make_subcommand())
        .subcommand(cmd_fasr::maf2fas::make_subcommand())
        .subcommand(cmd_fasr::name::make_subcommand())
        .subcommand(cmd_fasr::pl_p2m::make_subcommand())
        .subcommand(cmd_fasr::refine::make_subcommand())
        .subcommand(cmd_fasr::replace::make_subcommand())
        .subcommand(cmd_fasr::separate::make_subcommand())
        .subcommand(cmd_fasr::slice::make_subcommand())
        .subcommand(cmd_fasr::split::make_subcommand())
        .subcommand(cmd_fasr::stat::make_subcommand())
        .subcommand(cmd_fasr::subset::make_subcommand())
        .subcommand(cmd_fasr::variation::make_subcommand())
        .subcommand(cmd_fasr::xlsx::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* info: check / cover / link / name / stat
* creation: axt2fas / maf2fas / create
* records: separate / split / subset
* transform: filter / replace / refine
* transmute: concat / consensus / join / pl-p2m / slice
* variations: variation / xlsx

"###,
        );

    // Check which subcommand the user ran...
    match app.get_matches().subcommand() {
        // info
        Some(("check", sub_matches)) => cmd_fasr::check::execute(sub_matches),
        Some(("cover", sub_matches)) => cmd_fasr::cover::execute(sub_matches),
        Some(("link", sub_matches)) => cmd_fasr::link::execute(sub_matches),
        Some(("name", sub_matches)) => cmd_fasr::name::execute(sub_matches),
        Some(("stat", sub_matches)) => cmd_fasr::stat::execute(sub_matches),
        // creation
        Some(("axt2fas", sub_matches)) => cmd_fasr::axt2fas::execute(sub_matches),
        Some(("maf2fas", sub_matches)) => cmd_fasr::maf2fas::execute(sub_matches),
        Some(("create", sub_matches)) => cmd_fasr::create::execute(sub_matches),
        // records
        Some(("separate", sub_matches)) => cmd_fasr::separate::execute(sub_matches),
        Some(("split", sub_matches)) => cmd_fasr::split::execute(sub_matches),
        Some(("subset", sub_matches)) => cmd_fasr::subset::execute(sub_matches),
        // transform
        Some(("filter", sub_matches)) => cmd_fasr::filter::execute(sub_matches),
        Some(("replace", sub_matches)) => cmd_fasr::replace::execute(sub_matches),
        Some(("refine", sub_matches)) => cmd_fasr::refine::execute(sub_matches),
        // transmute
        Some(("concat", sub_matches)) => cmd_fasr::concat::execute(sub_matches),
        Some(("consensus", sub_matches)) => cmd_fasr::consensus::execute(sub_matches),
        Some(("join", sub_matches)) => cmd_fasr::join::execute(sub_matches),
        Some(("pl-p2m", sub_matches)) => cmd_fasr::pl_p2m::execute(sub_matches),
        Some(("slice", sub_matches)) => cmd_fasr::slice::execute(sub_matches),
        // variations
        Some(("variation", sub_matches)) => cmd_fasr::variation::execute(sub_matches),
        Some(("xlsx", sub_matches)) => cmd_fasr::xlsx::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO: add more tools - vcf, match
// TODO: fasr variation --indel
// TODO: fasr match
//  sparsemem -maxmatch -F -l %d -b -n -k 4 -threads 4 %s %s > %s
//  mummer -maxmatch -F -l %d -b -n %s %s > %s
//  $exe, $length, $genome, $query, $result

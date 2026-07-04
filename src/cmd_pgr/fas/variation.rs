use clap::{ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for variation.
pub fn make_subcommand() -> Command {
    Command::new("variation")
        .about("Lists variations (substitutions)")
        .after_help(
            r###"
Lists variations (substitutions) from block FA files in TSV format.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Filter out complex variations: `tsv-filter -H --ne freq:-1`
* Filter out singletons: `tsv-filter -H --ne freq:1`

Examples:
1. List substitutions from block FA files:
   pgr fas variation tests/fas/part1.fas

2. Handle outgroup (last sequence) for polarization:
   pgr fas variation tests/fas/part1.fas --outgroup

3. Output results to a file:
   pgr fas variation tests/fas/part1.fas -o output.tsv

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the variation command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let has_outgroup = args.get_flag("outgroup");

    let field_names = vec![
        "#target",
        "chr",
        "chr_pos",
        "range",
        "pos",
        "tbase",
        "qbase",
        "bases",
        "mutant_to",
        "freq",
        "pattern",
        "obase",
    ];

    // Operating
    writer.write_all(format!("{}\n", field_names.join("\t")).as_ref())?;

    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::reader(infile)?;

        while let Ok(block) = pgr::libs::fmt::fas::next_fas_block(&mut reader) {
            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }

            // target range and sequence intspan
            let first = match block.entries.first() {
                Some(e) => e,
                None => continue,
            };
            let trange = first.range().clone();
            let t_ints_seq = pgr::libs::alignment::seq_intspan(first.seq());

            // pos, tbase, qbase, bases, mutant_to, freq, pattern, obase
            //   0,     1,     2,     3,         4,    5,       6,     7
            let seq_count = seqs.len();
            if has_outgroup && seq_count < 2 {
                anyhow::bail!(
                    "outgroup mode requires at least 2 sequences per block, got {}",
                    seq_count
                );
            }
            let subs = if has_outgroup {
                let mut unpolarized = pgr::libs::alignment::get_subs(&seqs[..(seq_count - 1)])?;
                pgr::libs::alignment::polarize_subs(&mut unpolarized, seqs[seq_count - 1])?;
                unpolarized
            } else {
                pgr::libs::alignment::get_subs(&seqs)?
            };

            for s in subs {
                let chr = trange.chr();

                let chr_pos = pgr::libs::alignment::align_to_chr(
                    &t_ints_seq,
                    s.pos,
                    trange.start,
                    trange.strand(),
                )?;
                let var_rg = format!("{}:{}", chr, chr_pos);

                writer.write_all(
                    format!("{}\t{}\t{}\t{}\t{}\n", trange, chr, chr_pos, var_rg, s,).as_ref(),
                )?;
            }
        }
    }

    Ok(())
}

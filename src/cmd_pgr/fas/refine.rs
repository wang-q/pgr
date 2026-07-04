use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::fmt::Write;

/// Build the clap subcommand for refine.
pub fn make_subcommand() -> Command {
    Command::new("refine")
        .about("Realigns files with built-in or external programs and trim unwanted regions")
        .after_help(
            r###"
Realigns sequences in block FA files using built-in or external programs (clustalw, mafft, muscle, spoa) and trims unwanted regions.

Notes:
* Supports both plain text and gzipped (.gz) files
* Reads from stdin if input file is 'stdin'
* Supported MSA programs (`--engine`):
    * `builtin` (default): Uses built-in Rust implementation (Partial Order Alignment).
    * `clustalw`: Uses external `clustalw` command.
    * `mafft`: Uses external `mafft` command.
    * `muscle`: Uses external `muscle` command.
    * `spoa`: Uses external `spoa` command (SIMD optimized).
    * `none`: Skips realigning (useful for trimming only).
* The `--quick` option aligns only indel-adjacent regions (useful for .axt/.maf conversions)
* Supports parallel processing. The output order may differ from the original

Examples:
1. Realign block FA files using builtin (default):
   pgr fas refine tests/fas/part1.fas tests/fas/part2.fas

2. Realign using mafft with 4 threads:
   pgr fas refine tests/fas/part1.fas --engine mafft --parallel 4

3. Quick alignment for files converted from pairwise alignments:
   pgr fas refine tests/fas/part1.fas --quick --parallel 4

4. Output results to a file:
   pgr fas refine tests/fas/part1.fas -o output.fas

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg("block FA"))
        .arg(crate::cmd_pgr::args::engine_arg(
            &["builtin", "clustalw", "mafft", "muscle", "spoa", "none"],
            "builtin",
            "Aligning program (builtin/clustalw/mafft/muscle/spoa/none)",
        ))
        .arg(crate::cmd_pgr::args::outgroup_arg())
        .arg(
            Arg::new("chop")
                .long("chop")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("0")
                .help("Chop head and tail indels"),
        )
        .arg(
            Arg::new("is_quick")
                .long("quick")
                .action(ArgAction::SetTrue)
                .help("Quick mode, only aligns indel adjacent regions"),
        )
        .arg(
            Arg::new("indel_pad")
                .long("indel-pad")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("50")
                .help("In quick mode, enlarge indel regions"),
        )
        .arg(
            Arg::new("fill")
                .long("fill")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("50")
                .help("In quick mode, fill holes between indel"),
        )
        .arg(crate::cmd_pgr::args::parallel_arg())
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the refine command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let parallel = *args.get_one::<usize>("parallel").unwrap();
    let mut writer = pgr::writer(crate::cmd_pgr::args::get_outfile(args))?;
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    pgr::libs::fmt::fas::run_pipeline(&mut writer, &infiles, parallel, |block| {
        proc_block(block, args)
    })
}

fn proc_block(block: &pgr::libs::fmt::fas::FasBlock, args: &ArgMatches) -> anyhow::Result<String> {
    let engine = args.get_one::<String>("engine").unwrap();
    let has_outgroup = args.get_flag("outgroup");
    let chop = *args.get_one::<usize>("chop").unwrap();
    let is_quick = args.get_flag("is_quick");
    let pad = *args.get_one::<usize>("indel_pad").unwrap();
    let fill = *args.get_one::<usize>("fill").unwrap();

    let n = block.entries.len();
    let mut seqs: Vec<String> = Vec::with_capacity(n);
    let mut ranges = Vec::with_capacity(n);
    for entry in &block.entries {
        seqs.push(String::from_utf8(entry.seq().to_vec())?);
        ranges.push(entry.range().clone());
    }

    let mut aligned = vec![];
    if engine.as_str() == "none" {
        aligned = seqs;
    } else {
        if is_quick {
            let pad_i32 = i32::try_from(pad)
                .map_err(|_| anyhow::anyhow!("--indel-pad {} exceeds i32 range", pad))?;
            let fill_i32 = i32::try_from(fill)
                .map_err(|_| anyhow::anyhow!("--fill {} exceeds i32 range", fill))?;
            aligned = pgr::libs::alignment::align_seqs_quick(&seqs, engine, pad_i32, fill_i32)?;
        } else {
            aligned = pgr::libs::alignment::align_seqs(&seqs, engine)?;
        }
    };

    pgr::libs::alignment::trim_pure_dash(&mut aligned);
    if has_outgroup {
        pgr::libs::alignment::trim_outgroup(&mut aligned)?;
        pgr::libs::alignment::trim_complex_indel(&mut aligned)?;
    }

    if chop > 0 {
        pgr::libs::alignment::trim_head_tail(&mut aligned, &mut ranges, chop);
    }

    let mut out_string = String::new();
    for (range, seq) in ranges.iter().zip(aligned) {
        writeln!(out_string, ">{}\n{}", range, seq)?;
    }

    // end of a block
    out_string.push('\n');

    Ok(out_string)
}

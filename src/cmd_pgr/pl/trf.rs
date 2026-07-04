use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use std::io::BufRead;

/// Build the clap subcommand for trf.
pub fn make_subcommand() -> Command {
    Command::new("trf")
        .about("Identifies tandem repeats in a genome")
        .after_help(
            r###"
This command identifies tandem repeats in a genome via `trf`.

* <infile> is path to fasta file, .fa.gz is supported. Cannot be stdin.

* All operations are running in a tempdir and no intermediate files are retained.

* External dependencies
    * trf
    * spanr

"###,
        )
        .arg(crate::cmd_pgr::args::infile_arg_required_with_help(
            "Input file to process",
        ))
        .arg(
            Arg::new("trf_match")
                .long("trf-match")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(usize))
                .help("TRF matching weight"),
        )
        .arg(
            Arg::new("trf_mismatch")
                .long("trf-mismatch")
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("TRF mismatching penalty"),
        )
        .arg(
            Arg::new("delta")
                .long("delta")
                .num_args(1)
                .default_value("7")
                .value_parser(value_parser!(usize))
                .help("Indel penalty"),
        )
        .arg(
            Arg::new("pm")
                .long("pm")
                .num_args(1)
                .default_value("80")
                .value_parser(value_parser!(usize))
                .help("Match probability"),
        )
        .arg(
            Arg::new("pi")
                .long("pi")
                .num_args(1)
                .default_value("10")
                .value_parser(value_parser!(usize))
                .help("Indel probability"),
        )
        .arg(crate::cmd_pgr::args::min_score_arg("50"))
        .arg(
            Arg::new("max_period")
                .long("max-period")
                .num_args(1)
                .default_value("2000")
                .value_parser(value_parser!(usize))
                .help("Maximum period size to report"),
        )
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the trf command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);

    let opt_trf_match = *args.get_one::<usize>("trf_match").unwrap();
    let opt_trf_mismatch = *args.get_one::<usize>("trf_mismatch").unwrap();
    let opt_delta = *args.get_one::<usize>("delta").unwrap();
    let opt_pm = *args.get_one::<usize>("pm").unwrap();
    let opt_pi = *args.get_one::<usize>("pi").unwrap();
    let opt_minscore = *args.get_one::<f64>("min_score").unwrap();
    if !opt_minscore.is_finite() || opt_minscore < 0.0 {
        anyhow::bail!("--minscore must be non-negative finite: {}", opt_minscore);
    }
    if opt_minscore.fract() != 0.0 {
        anyhow::bail!("--minscore must be an integer: {}", opt_minscore);
    }
    let opt_minscore_u = opt_minscore as usize;
    let opt_max_period = *args.get_one::<usize>("max_period").unwrap();

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_trf_")?;
    let pgr = ctx.pgr.clone();

    run_cmd!(info "==> Absolute paths")?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    ctx.enter()?;

    run_cmd!(info "==> Split by names")?;
    run_cmd!(
        ${pgr} fa split name ${abs_infile} -o .
    )?;

    run_cmd!(info "==> Process each chromosome")?;
    run_cmd!(
        ${pgr} fa size ${abs_infile} -o chr.sizes
    )?;
    let chrs = pgr::libs::io::read_names::<Vec<String>>("chr.sizes")?;

    let mut rg_files = vec![];
    for (i, chr) in chrs.iter().enumerate() {
        run_cmd!(
            trf ${chr}.fa ${opt_trf_match} ${opt_trf_mismatch} ${opt_delta} ${opt_pm} ${opt_pi} ${opt_minscore_u} ${opt_max_period} -d -h -ngs > trf.${i}.dat
        )?;

        // 198 229 12 2.7 12 90 0 50 34 46 3 15 1.62 CATTACCACCAC CATTAGCACCACCATTACCACCACCATCACCA ATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACG TTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAA
        // start end
        // period copy_number consensus_pattern_size
        // perc_matches perc_indels
        // alignment_score
        // perc_a perc_c perc_g perc_t
        // entropy
        // consensus_pattern
        // repeat_seq
        // 15 fields
        // The last 2 fields were introduced by -ngs
        // Matched with `pgr fa range mg1655.fa NC_000913:198-229`

        let reader = pgr::reader(&format!("trf.{}.dat", i))?;

        let rg_file = format!("trf.{}.rg", i);
        let mut writer = pgr::writer(&rg_file)?;
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_ascii_whitespace().collect();
            if fields.len() < 15 {
                log::debug!("skipping short TRF line: {}", line);
                continue;
            }

            let start = fields[0].parse::<usize>()?;
            let end = fields[1].parse::<usize>()?;

            writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
        }
        rg_files.push(rg_file);
    }

    run_cmd!(info "==> Outputs")?;
    run_cmd!(
        spanr cover $[rg_files] -o ${abs_outfile}
    )?;

    // Done
    ctx.leave()?;

    Ok(())
}

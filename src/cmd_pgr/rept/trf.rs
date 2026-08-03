use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use std::io::Write;

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
    if opt_minscore > usize::MAX as f64 {
        anyhow::bail!("--minscore too large: {}", opt_minscore);
    }
    let opt_minscore_u = opt_minscore as usize;
    let opt_max_period = *args.get_one::<usize>("max_period").unwrap();

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_rept_trf_")?;
    let pgr = ctx.pgr.clone();

    run_cmd!(info "==> Absolute paths")?;
    let abs_infile = ctx.abs_path(args.get_one::<String>("infile").unwrap())?;
    let abs_outfile = pgr::libs::pl::abs_path_or_stdout(outfile)?;

    let _cwd_guard = ctx.enter()?;

    run_cmd!(info "==> Split by names")?;
    run_cmd!(
        ${pgr} fa split name ${abs_infile} -o .
    )?;

    run_cmd!(info "==> Process each chromosome")?;
    run_cmd!(
        ${pgr} fa size ${abs_infile} -o chr.sizes
    )?;
    let chrs = pgr::libs::io::read_names::<Vec<String>>("chr.sizes")?;

    // `spanr cover` truncates dotted contig names (e.g. `NC_000913.1` ->
    // `1`) at the last '.', so map real names to dot-free placeholders and
    // restore them after the spanr pass.
    let mut name_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let mut safe_map: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    let safe_chrs: Vec<String> = chrs
        .iter()
        .map(|c| {
            let s = format!("c{}", name_map.len() + 1);
            name_map.insert(c.clone(), s.clone());
            safe_map.insert(s.clone(), c.clone());
            s
        })
        .collect();

    let mut rg_files = vec![];
    for (i, chr) in chrs.iter().enumerate() {
        // `fa split name` names each file after `sanitize_filename(name)`;
        // use the same sanitization so names with special chars resolve.
        let chr_file = pgr::libs::io::sanitize_filename(chr);
        run_cmd!(
            trf ${chr_file}.fa ${opt_trf_match} ${opt_trf_mismatch} ${opt_delta} ${opt_pm} ${opt_pi} ${opt_minscore_u} ${opt_max_period} -d -h -ngs > trf.${i}.dat
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

        let dat_file = format!("trf.{}.dat", i);
        let reader = pgr::reader(&dat_file)
            .with_context(|| format!("Failed to open reader for {}", dat_file))?;

        let rg_file = format!("trf.{}.rg", i);
        let mut writer = pgr::writer(&rg_file)
            .with_context(|| format!("Failed to open writer for {}", rg_file))?;
        pgr::libs::pl::parse_trf_output(reader, &safe_chrs[i], &mut writer)?;
        writer.flush()?;
        rg_files.push(rg_file);
    }

    run_cmd!(info "==> Outputs")?;
    if pgr::libs::pl::count_rg_lines(&rg_files)? == 0 {
        let empty = b"{}\n";
        if abs_outfile == "stdout" {
            use std::io::Write;
            std::io::stdout().write_all(empty)?;
        } else {
            std::fs::write(&abs_outfile, empty)?;
        }
        return Ok(());
    }
    let mut set: std::collections::BTreeMap<String, pgr::libs::ds::IntSpan> =
        std::collections::BTreeMap::new();
    for rg in &rg_files {
        let reader = pgr::reader(rg)?;
        for (chr, is) in pgr::libs::runlist::rg_to_set(reader)? {
            set.entry(chr).or_default().merge(&is);
        }
    }
    let json = pgr::libs::ds::intspan::set2json(&set);
    std::fs::write("out.json", serde_json::to_vec_pretty(&json)?)?;

    // Restore the real contig names in the runlist json.
    let mut val: serde_json::Value = serde_json::from_slice(&std::fs::read("out.json")?)?;
    if let Some(obj) = val.as_object_mut() {
        let old = std::mem::take(obj);
        for (k, v) in old {
            // Drop the empty marker `-` so the runlist stays clean.
            if v.as_str() == Some("-") {
                continue;
            }
            obj.insert(safe_map.get(&k).cloned().unwrap_or(k), v);
        }
    }
    let out_bytes = serde_json::to_vec_pretty(&val)?;
    if abs_outfile == "stdout" {
        use std::io::Write;
        let mut w = pgr::writer("stdout")?;
        w.write_all(&out_bytes)?;
        w.write_all(b"\n")?;
    } else {
        std::fs::write(&abs_outfile, out_bytes)?;
    }

    // Done

    Ok(())
}

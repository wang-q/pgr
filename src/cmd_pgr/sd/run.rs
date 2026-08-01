//! `pgr sd run` — end-to-end SD pipeline: search -> align -> cluster -> decompose -> cover.

use clap::{value_parser, Arg, ArgMatches, Command};
use cmd_lib::run_cmd;
use std::io::Write;

/// Build the clap subcommand for run.
pub fn make_subcommand() -> Command {
    Command::new("run")
        .about("Runs the full SD pipeline (search/align/cluster/decompose/cover)")
        .after_help(
            r###"
Runs the whole segmental duplication pipeline on one genome:
`pgr sd search` -> `pgr sd align` -> `pgr sd cluster` -> `pgr sd decompose`
-> `pgr sd cover`. The final CORE-annotated elementary BED is written to
`<outdir>/out.elem.bed`; intermediate files live in a temp workspace.

Examples:
1. Full pipeline:
   pgr sd run genome.fa -o sd_out/
"###,
        )
        .arg(
            Arg::new("genome")
                .index(1)
                .required(true)
                .help("Genome FASTA file"),
        )
        .arg(
            Arg::new("outdir")
                .long("outdir")
                .short('o')
                .default_value("sd_out")
                .help("Output directory"),
        )
        .arg(
            Arg::new("preset")
                .long("preset")
                .value_parser(clap::builder::PossibleValuesParser::new(
                    pgr::libs::lastz::preset_names(),
                ))
                .help("lastz parameter set for search"),
        )
        .arg(
            Arg::new("min_len")
                .long("min-len")
                .default_value("1000")
                .value_parser(value_parser!(u32))
                .help("Minimum SD block length in bp"),
        )
        .arg(
            Arg::new("min_identity")
                .long("min-identity")
                .default_value("0.90")
                .value_parser(value_parser!(f64))
                .help("Minimum SD block identity"),
        )
}
/// Execute the run command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let genome = args.get_one::<String>("genome").unwrap();
    let outdir = args.get_one::<String>("outdir").unwrap();
    let preset = args.get_one::<String>("preset").cloned();
    let min_len = *args.get_one::<u32>("min_len").unwrap();
    let min_identity = *args.get_one::<f64>("min_identity").unwrap();

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_sd_run_")?;
    let pgr = ctx.pgr.clone();
    let abs_genome = ctx.abs_path(genome)?;
    let abs_outdir = ctx.abs_path(outdir)?;
    let _cwd_guard = ctx.enter()?;

    let preset_args = preset.map(|p| format!("--preset {p}")).unwrap_or_default();
    run_cmd!(${pgr} sd search ${abs_genome} -o hits.psl ${preset_args} --min-len ${min_len} --min-identity ${min_identity})?;
    run_cmd!(${pgr} sd align ${abs_genome} hits.psl -o hits.paf)?;
    run_cmd!(${pgr} sd cluster ${abs_genome} hits.paf -o clusters)?;

    // Decompose each cluster and merge the elementary BEDs with global
    // set_id renumbering (each cluster restarts set_id at 1).
    let mut elems = String::new();
    let mut set_offset = 0u32;
    for fa in pgr::libs::io::list_files_ext("clusters", "fa") {
        let stem = pgr::libs::io::basename_or_err(&fa)?;
        let stem = stem.trim_end_matches(".fa");
        run_cmd!(${pgr} sd decompose ${fa} -o clusters/${stem}.elem.bed)?;
        let content = std::fs::read_to_string(format!("clusters/{stem}.elem.bed"))?;
        let mut cluster_max = 0u32;
        for line in content.lines() {
            if line.is_empty() {
                continue;
            }
            let mut f: Vec<&str> = line.split('\t').collect();
            let sid: u32 = f[4].parse()?;
            cluster_max = cluster_max.max(sid);
            elems.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                f[0],
                f[1],
                f[2],
                f[3],
                sid + set_offset,
                f[5],
                f[6],
                f[7]
            ));
        }
        set_offset += cluster_max;
    }
    std::fs::write("elems.bed", elems)?;

    std::fs::create_dir_all(&abs_outdir)?;
    let out = format!("{}/out.elem.bed", abs_outdir.trim_end_matches('/'));
    run_cmd!(${pgr} sd cover hits.paf elems.bed -o ${out})?;

    let mut w = pgr::writer("stdout")?;
    writeln!(w, "wrote {}", out)?;
    Ok(())
}

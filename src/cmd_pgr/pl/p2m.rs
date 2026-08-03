use clap::{ArgMatches, Command};
use cmd_lib::run_cmd;
use pgr::libs::pl::CwdGuard;
use std::collections::BTreeMap;
use std::{env, fs};

/// Build the clap subcommand for p2m.
pub fn make_subcommand() -> Command {
    Command::new("p2m")
        .about("Runs the pairwise-to-multiple alignments pipeline")
        .after_help(
            r###"
Pairwise to Multiple (p2m) Pipeline

This pipeline constructs a "core" Multiple Sequence Alignment (MSA) from multiple
pairwise alignment files (Block FASTA). It identifies the intersection of genomic
regions covered by all inputs (anchored to the reference) and stitches them together.

Key Features:
* Reference-Based: The first species of the first input file is used as the
  reference target. All inputs must be pairwise alignments against this target.
* Intersection Logic: Only regions present in ALL input files are retained.
  This results in a gap-free core alignment.
* Automation: Runs `fas cover` (range extraction), `runlist compare`
  (range intersection), `fas slice` (sequence extraction), and `fas join`
  (alignment merging).

Dependencies:
* `pgr`: This binary itself.

Notes:
* <infiles> can be plain or gzipped (.fas.gz) block fasta files.
* Input cannot be stdin.
* Output is written to a directory (default: `PL-p2m`), not stdout.

"###,
        )
        .arg(crate::cmd_pgr::args::infiles_arg_with_numargs(
            "Input files",
            2..,
        ))
        .arg(crate::cmd_pgr::args::outdir_arg_with_default("PL-p2m"))
}

/// Execute the p2m command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outdir = args.get_one::<String>("outdir").unwrap();
    fs::create_dir_all(outdir)?;

    let curdir = env::current_dir()?;
    let pgr = pgr::libs::io::current_exe_string()?;

    run_cmd!(echo "==> Paths")?;
    run_cmd!(echo "    \"pgr\"      = ${pgr}")?;
    run_cmd!(echo "    \"curdir\" = ${curdir:?}")?;
    run_cmd!(echo "    \"outdir\" = ${outdir}")?;

    // basename => [abs_path, .json, .slice.fas]
    let mut info_of: BTreeMap<String, Vec<String>> = BTreeMap::new();

    // Operating
    run_cmd!(echo "==> Basenames and absolute paths")?;
    for infile in args.get_many::<String>("infiles").unwrap() {
        let basename = pgr::libs::io::basename_or_err(infile)?;
        let absolute = pgr::libs::io::absolute_path(infile)?.display().to_string();

        info_of.insert(basename.to_string(), vec![absolute.to_string()]);
    }

    run_cmd!(echo "==> Switch to outdir")?;
    let _cwd_guard = CwdGuard::enter(outdir)?;

    run_cmd!(echo "==> pgr fas name - first")?;
    let target_name = first_target_name(&pgr, &info_of)?;
    run_cmd!(echo "    \"target_name\" = ${target_name}")?;

    run_cmd!(echo "==> pgr fas cover")?;
    for (basename, info) in info_of.iter_mut() {
        let infile = info
            .first()
            .ok_or_else(|| anyhow::anyhow!("no info for basename {}", basename))?;
        let outfile = format!("{}.json", basename);
        run_cmd!(${pgr} fas cover ${infile} --trim 10 --name ${target_name} -o ${outfile})?;

        info.push(outfile.to_string());
    }

    run_cmd!(echo "==> runlist compare")?;
    {
        let infiles: Vec<String> = info_of
            .iter()
            .map(|e| {
                e.1.get(1)
                    .ok_or_else(|| anyhow::anyhow!("missing cover output for {}", e.0))
                    .map(|s| s.to_string())
            })
            .collect::<Result<_, _>>()?;
        let first = pgr::libs::runlist::read_json(&infiles[0])?;
        let first = pgr::libs::runlist::json_to_sets(&first);
        let mut others = Vec::new();
        for f in &infiles[1..] {
            let json = pgr::libs::runlist::read_json(f)?;
            others.push(pgr::libs::runlist::json_to_set(&json));
        }
        let mut res = pgr::libs::runlist::compare_sets(
            &first,
            &others,
            pgr::libs::runlist::CompareOp::Intersect,
        );
        // Equivalent of the old `spanr span stdin --op excise -n 10`.
        for set in res.values_mut() {
            *set = pgr::libs::runlist::span_op(set, pgr::libs::runlist::SpanOp::Excise, 10);
        }
        pgr::libs::runlist::write_sets("intersect.json", &res)?;

        let mut all = infiles.clone();
        all.push("intersect.json".to_string());
        let merged = pgr::libs::runlist::merge_files(&all, false)?;
        pgr::libs::ds::intspan::write_json("merge.json", &merged)?;
    }

    run_cmd!(echo "==> pgr fas slice")?;
    for (basename, info) in info_of.iter_mut() {
        let infile = info
            .first()
            .ok_or_else(|| anyhow::anyhow!("no info for basename {}", basename))?;
        let outfile = format!("{}.slice.fas", basename);
        run_cmd!(${pgr} fas slice ${infile} --runlist intersect.json --name ${target_name} -o ${outfile})?;

        info.push(outfile.to_string());
    }

    run_cmd!(echo "==> pgr fas join")?;
    {
        let infiles: Vec<String> = info_of
            .iter()
            .map(|e| {
                e.1.get(2)
                    .ok_or_else(|| anyhow::anyhow!("missing slice output for {}", e.0))
                    .map(|s| s.to_string())
            })
            .collect::<Result<_, _>>()?;
        run_cmd!(
            ${pgr} fas join $[infiles] --name ${target_name} -o join.raw.fas
        )?;
    }

    run_cmd!(echo "==> pgr fas name && pgr fas subset")?;
    {
        run_cmd!(
            ${pgr} fas name join.raw.fas -o name.lst
        )?;
        run_cmd!(
            ${pgr} fas subset join.raw.fas --required name.lst -o join.subset.fas
        )?;
    }

    Ok(())
}

// Run `pgr fas name` on the first input file and return the first name line.
fn first_target_name(pgr: &str, info_of: &BTreeMap<String, Vec<String>>) -> anyhow::Result<String> {
    let infile = info_of
        .values()
        .next()
        .and_then(|v| v.first())
        .ok_or_else(|| anyhow::anyhow!("no cover input found"))?;
    run_cmd!(
        ${pgr} fas name ${infile} -o name.first.lst
    )?;
    let first_content = fs::read_to_string("name.first.lst")
        .map_err(|e| anyhow::anyhow!("failed to read name.first.lst: {}", e))?;
    first_content
        .lines()
        .next()
        .ok_or_else(|| anyhow::anyhow!("name.first.lst is empty"))
        .map(str::to_string)
}

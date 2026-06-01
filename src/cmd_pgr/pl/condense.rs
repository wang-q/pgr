use clap::*;
use cmd_lib::*;
use itertools::Itertools;
use std::env;
use std::fs;
use std::io::Write;
use tempfile::TempDir;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("condense")
        .about("Pipeline - condense subtrees based on taxonomy")
        .after_help(
            r###"
* <taxon.tsv> is a tab-separated file without header, containing at least 2 columns:
    node_name   taxonomic_term   [additional columns...]
* The first column is the node name (matching leaf labels in the Newick file)
* Use `--rank` to specify which column(s) to use for grouping (1-based index, default: 2)
* Can specify multiple `--rank` values to condense at multiple levels
* Monophyletic subtrees with the same taxonomic term will be condensed
* Condensed nodes are named as {term}___{count}

Examples:
1. Condense by species (2nd column):
   pgr pl condense -t taxon.tsv tree.nwk

2. Condense by genus (3rd column):
   pgr pl condense -t taxon.tsv -r 3 tree.nwk

3. Condense by multiple ranks:
   pgr pl condense -t taxon.tsv -r 2 -r 3 tree.nwk

4. Output mapping file:
   pgr pl condense -t taxon.tsv --map tree.nwk -o condensed.nwk

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input Newick filename. [stdin] for standard input"),
        )
        .arg(
            Arg::new("taxon")
                .long("taxon")
                .short('t')
                .num_args(1)
                .required(true)
                .help("Path to taxonomy TSV file"),
        )
        .arg(
            Arg::new("rank")
                .long("rank")
                .short('r')
                .num_args(1)
                .action(ArgAction::Append)
                .value_parser(value_parser!(usize))
                .help("Column index(es) to use for grouping (1-based, can be specified multiple times, default: 2)"),
        )
        .arg(
            Arg::new("map")
                .long("map")
                .action(ArgAction::SetTrue)
                .help("Write a map file `condensed.tsv`"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outfile = args.get_one::<String>("outfile").unwrap();
    let taxon_file = args.get_one::<String>("taxon").unwrap();

    let ranks: Vec<usize> = if args.contains_id("rank") {
        args.get_many::<usize>("rank").unwrap().copied().collect()
    } else {
        vec![2] // default to column 2
    };

    let curdir = env::current_dir()?;
    let exe = env::current_exe().unwrap().display().to_string();
    let tempdir = TempDir::new()?;
    let tempdir_str = tempdir.path().to_str().unwrap();
    let curdir_str = curdir.display().to_string();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"pgr\"     = ${exe}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir_str}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    //----------------------------
    // Operating
    //----------------------------
    run_cmd!(info "==> Absolute paths")?;
    let infile = args.get_one::<String>("infile").unwrap();
    let abs_infile = if infile == "stdin" {
        "stdin".to_string()
    } else {
        intspan::absolute_path(infile)
            .unwrap()
            .display()
            .to_string()
    };
    let abs_taxon = intspan::absolute_path(taxon_file)
        .unwrap()
        .display()
        .to_string();

    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> Start")?;
    run_cmd!(
        ${exe} nwk indent ${abs_infile} -o start.nwk
    )?;

    run_cmd!(info "==> Labels in the file")?;
    run_cmd!(
        ${exe} nwk label start.nwk -o labels.lst
    )?;

    run_cmd!(info "==> Create replace.tsv from taxon.tsv")?;
    let mut writer = pgr::writer("replace.tsv");
    for line in pgr::read_lines(&abs_taxon) {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        let node_name = parts[0];
        // Write all specified ranks as additional columns
        let mut row = vec![node_name.to_string()];
        for rank_col in &ranks {
            let rank_idx = rank_col.saturating_sub(1);
            if let Some(term) = parts.get(rank_idx) {
                row.push(newick_safe(term));
            }
        }
        if row.len() > 1 {
            writer.write_all(format!("{}\n", row.join("\t")).as_ref())?;
        }
    }
    writer.flush()?;

    run_cmd!(info "==> Add taxonomy info to the tree")?;
    run_cmd!(
        ${exe} nwk replace start.nwk replace.tsv --mode label -o commented.nwk
    )?;

    run_cmd!(info "==> Build groups")?;
    let mut groups = vec![];
    for line in pgr::read_lines("replace.tsv") {
        let parts: Vec<&str> = line.split('\t').collect();
        // Collect all taxonomic terms from columns 1+ (skip column 0 which is node_name)
        for i in 1..parts.len() {
            groups.push(parts[i].to_string());
        }
    }
    groups = groups.into_iter().unique().filter(|s| s.ne("NA")).collect();

    run_cmd!(info "==> Condensing")?;
    let mut cur_tree = "commented.nwk".to_string();
    let mut condensed = vec![];
    for group in groups.iter() {
        let labels: Vec<String> = run_fun!(
            ${exe} nwk label ${cur_tree} -n ${group} -M
        )
        .unwrap()
        .split('\n')
        .map(|s| s.to_string())
        .filter(|s| !s.is_empty())
        .collect();

        if labels.is_empty() {
            continue;
        }

        let new_label = format!("{}___{}", group, labels.len());

        labels
            .iter()
            .for_each(|e| condensed.push(format!("{}\t{}", e, new_label)));

        run_cmd!(
            ${exe} nwk subtree ${cur_tree} -n ${group} -M --condense ${new_label} -o condense.${group}.nwk
        )?;

        cur_tree = format!("condense.{}.nwk", group);
    }

    run_cmd!(info "==> Results")?;
    fs::copy(
        tempdir.path().join(cur_tree.clone()).to_str().unwrap(),
        "result.nwk",
    )?;

    let mut writer = pgr::writer("condensed.tsv");
    for line in condensed.iter() {
        writer.write_all(format!("{}\n", line).as_ref())?;
    }
    writer.flush()?;

    //----------------------------
    // Done
    //----------------------------
    if outfile == "stdout" {
        run_cmd!(cat ${cur_tree})?;
        env::set_current_dir(&curdir)?;
    } else {
        env::set_current_dir(&curdir)?;
        fs::copy(tempdir.path().join("result.nwk").to_str().unwrap(), outfile)?;
    }

    if args.get_flag("map") {
        fs::copy(
            tempdir.path().join("condensed.tsv").to_str().unwrap(),
            "condensed.tsv",
        )?;
    }

    Ok(())
}

fn newick_safe(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '(' | ')' | '[' | ']' | ',' | ':' | ';' | ' ' => '_',
            _ => c,
        })
        .collect()
}
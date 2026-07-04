use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use cmd_lib::{run_cmd, run_fun};
use itertools::Itertools;
use pgr::libs::phylo::tree::Tree;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;

/// Build the clap subcommand for condense.
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
   pgr pl condense --taxon taxon.tsv tree.nwk

2. Condense by genus (3rd column):
   pgr pl condense --taxon taxon.tsv --rank 3 tree.nwk

3. Condense by multiple ranks:
   pgr pl condense --taxon taxon.tsv --rank 2 --rank 3 tree.nwk

4. Output mapping file:
   pgr pl condense --taxon taxon.tsv --map tree.nwk -o condensed.nwk

"###,
        )
        .arg(
            crate::cmd_pgr::args::infile_arg_required_with_help(
                "Input Newick filename. [stdin] for standard input",
            ),
        )
        .arg(
            Arg::new("taxon")
                .long("taxon")
                .num_args(1)
                .required(true)
                .help("Path to taxonomy TSV file"),
        )
        .arg(
            Arg::new("rank")
                .long("rank")
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
        .arg(crate::cmd_pgr::args::outfile_arg())
}

/// Execute the condense command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outfile = crate::cmd_pgr::args::get_outfile(args);
    let taxon_file = args.get_one::<String>("taxon").unwrap();

    let ranks: Vec<usize> = if args.contains_id("rank") {
        args.get_many::<usize>("rank").unwrap().copied().collect()
    } else {
        vec![2] // default to column 2
    };
    for rank_col in &ranks {
        anyhow::ensure!(*rank_col >= 1, "--rank must be >= 1, got {}", rank_col);
    }

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_condense_")?;
    let exe = ctx.pgr.clone();

    // Operating
    run_cmd!(info "==> Absolute paths")?;
    let infile = args.get_one::<String>("infile").unwrap();
    let abs_infile = if infile == "stdin" {
        "stdin".to_string()
    } else {
        ctx.abs_path(infile)?
    };
    let abs_taxon = ctx.abs_path(taxon_file)?;

    ctx.enter()?;

    // Read tree leaf names for filtering
    let trees = Tree::from_file(&abs_infile)?;
    let leaf_names: BTreeSet<String> = if let Some(tree) = trees.first() {
        tree.get_leaf_names().into_iter().flatten().collect()
    } else {
        BTreeSet::new()
    };
    let leaf_count = leaf_names.len();
    run_cmd!(info "    Tree leaf count: ${leaf_count}")?;

    // Read taxonomy TSV
    run_cmd!(info "==> Read taxonomy TSV")?;

    // taxon_map: node_name -> Vec of terms (one per rank, None if column missing)
    let mut taxon_map: BTreeMap<String, Vec<Option<String>>> = BTreeMap::new();
    // groups: all unique terms for each rank
    let mut all_groups: Vec<Vec<String>> = vec![vec![]; ranks.len()];

    for line in pgr::read_lines(&abs_taxon)? {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        let node_name = parts[0].to_string();
        if !leaf_names.contains(&node_name) {
            continue;
        }
        let mut terms: Vec<Option<String>> = Vec::with_capacity(ranks.len());

        for (i, rank_col) in ranks.iter().enumerate() {
            let rank_idx = rank_col.saturating_sub(1);
            let term = parts.get(rank_idx).map(|s| pgr::libs::phylo::newick_safe(s));
            if let Some(t) = &term {
                all_groups[i].push(t.clone());
            }
            terms.push(term);
        }

        if terms.iter().any(|t| t.is_some()) {
            taxon_map.insert(node_name, terms);
        }
    }

    // Deduplicate groups and filter NA
    for groups in &mut all_groups {
        *groups = groups
            .clone()
            .into_iter()
            .sorted()
            .unique()
            .filter(|s| s.ne("NA"))
            .collect();
    }

    let taxon_count = taxon_map.len();
    run_cmd!(info "    Loaded ${taxon_count} taxonomy entries")?;
    for (i, rank) in ranks.iter().enumerate() {
        let rank_groups = all_groups[i].len();
        run_cmd!(info "    Rank ${rank}: ${rank_groups} groups")?;
    }

    // Start - copy input tree
    run_cmd!(info "==> Start")?;
    run_cmd!(
        ${exe} nwk indent ${abs_infile} -o start.nwk
    )?;

    // Condensing - process each rank
    run_cmd!(info "==> Condensing")?;
    let mut cur_tree = "start.nwk".to_string();
    let mut condensed: Vec<String> = vec![];
    let mut toggle = false;

    for (rank_idx, groups) in all_groups.iter().enumerate() {
        let rank_num = ranks[rank_idx];
        run_cmd!(info "    Processing rank ${rank_num}")?;

        for group in groups.iter() {
            // Find all original nodes that belong to this group at this rank
            let nodes_in_group: Vec<String> = taxon_map
                .iter()
                .filter(|(_, terms)| {
                    terms
                        .get(rank_idx)
                        .and_then(|t| t.as_deref())
                        .map(|t| t == group)
                        .unwrap_or(false)
                })
                .map(|(name, _)| name.clone())
                .collect();

            if nodes_in_group.len() < 2 {
                continue;
            }

            // Write node list to a reusable file
            let mut writer = pgr::writer("nodes.txt")?;
            for node in &nodes_in_group {
                writeln!(writer, "{}", node)?;
            }
            writer.flush()?;
            drop(writer);

            // Check if these nodes form a monophyletic group and get labels.
            // nwk label -M exits 0 with empty output for non-monophyletic groups;
            // a non-zero exit indicates a real error and is propagated via `?`.
            let labels_output = run_fun!(
                ${exe} nwk label ${cur_tree} -l nodes.txt -M
            )?;

            let labels: Vec<String> = labels_output
                .split('\n')
                .map(|s: &str| s.to_string())
                .filter(|s: &String| !s.is_empty())
                .collect();

            if labels.is_empty() {
                // Not monophyletic, skip
                continue;
            }

            let new_label = format!("{}___{}", group, nodes_in_group.len());

            // Record mapping: original node name -> condensed label
            for node in &nodes_in_group {
                condensed.push(format!("{}\t{}", node, new_label));
            }

            // Condense the subtree; use alternating output files to avoid clobbering
            toggle = !toggle;
            let new_tree = if toggle {
                "temp_a.nwk".to_string()
            } else {
                "temp_b.nwk".to_string()
            };
            run_cmd!(
                ${exe} nwk subtree ${cur_tree} -l nodes.txt -M --condense ${new_label} -o ${new_tree}
            )?;

            cur_tree = new_tree;
        }
    }

    // Results
    run_cmd!(info "==> Results")?;
    fs::copy(
        ctx.tempdir
            .path()
            .join(&cur_tree)
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("path is not utf-8"))?,
        "result.nwk",
    )?;

    let mut writer = pgr::writer("condensed.tsv")?;
    for line in condensed.iter() {
        writeln!(writer, "{}", line)?;
    }
    writer.flush()?;

    // Done
    if outfile == "stdout" {
        let result_content = fs::read_to_string("result.nwk")?;
        print!("{}", result_content);
        ctx.leave()?;
    } else {
        ctx.leave()?;
        fs::copy(
            ctx.tempdir
                .path()
                .join("result.nwk")
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("path is not utf-8"))?,
            outfile,
        )?;
    }

    if args.get_flag("map") {
        fs::copy(
            ctx.tempdir
                .path()
                .join("condensed.tsv")
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("path is not utf-8"))?,
            "condensed.tsv",
        )?;
    }

    Ok(())
}

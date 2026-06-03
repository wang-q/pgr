use clap::*;
use cmd_lib::*;
use itertools::Itertools;
use std::collections::BTreeMap;
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

    //----------------------------
    // Read taxonomy TSV
    //----------------------------
    run_cmd!(info "==> Read taxonomy TSV")?;
    
    // taxon_map: node_name -> Vec of terms (one per rank)
    let mut taxon_map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    // groups: all unique terms for each rank
    let mut all_groups: Vec<Vec<String>> = vec![vec![]; ranks.len()];
    
    for line in pgr::read_lines(&abs_taxon) {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        let node_name = parts[0].to_string();
        let mut terms = vec![];
        
        for (i, rank_col) in ranks.iter().enumerate() {
            let rank_idx = rank_col.saturating_sub(1);
            if let Some(term) = parts.get(rank_idx) {
                let term = newick_safe(term);
                terms.push(term.clone());
                all_groups[i].push(term);
            }
        }
        
        if !terms.is_empty() {
            taxon_map.insert(node_name, terms);
        }
    }
    
    // Deduplicate groups and filter NA
    for groups in &mut all_groups {
        *groups = groups.clone().into_iter().sorted().unique().filter(|s| s.ne("NA")).collect();
    }
    
    let taxon_count = taxon_map.len();
    run_cmd!(info "    Loaded ${taxon_count} taxonomy entries")?;
    for (i, rank) in ranks.iter().enumerate() {
        let rank_groups = all_groups[i].len();
        run_cmd!(info "    Rank ${rank}: ${rank_groups} groups")?;
    }

    //----------------------------
    // Start - copy input tree
    //----------------------------
    run_cmd!(info "==> Start")?;
    run_cmd!(
        ${exe} nwk indent ${abs_infile} -o start.nwk
    )?;

    //----------------------------
    // Condensing - process each rank
    //----------------------------
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
                .filter(|(_, terms)| terms.get(rank_idx).map(|t| t == group).unwrap_or(false))
                .map(|(name, _)| name.clone())
                .collect();

            if nodes_in_group.len() < 2 {
                continue;
            }

            // Write node list to a reusable file
            let mut writer = pgr::writer("nodes.txt");
            for node in &nodes_in_group {
                writer.write_all(format!("{}\n", node).as_ref())?;
            }
            writer.flush()?;
            drop(writer);

            // Check if these nodes form a monophyletic group and get labels
            let labels_result = run_fun!(
                ${exe} nwk label ${cur_tree} -f nodes.txt -M
            );
            
            let labels_output = match labels_result {
                Ok(output) => output,
                Err(_) => continue, // Not monophyletic or error
            };

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
            let new_tree = if toggle { "temp_a.nwk".to_string() } else { "temp_b.nwk".to_string() };
            run_cmd!(
                ${exe} nwk subtree ${cur_tree} -f nodes.txt -M --condense ${new_label} -o ${new_tree}
            )?;

            cur_tree = new_tree;
        }
    }

    //----------------------------
    // Results
    //----------------------------
    run_cmd!(info "==> Results")?;
    fs::copy(
        tempdir.path().join(&cur_tree).to_str().unwrap(),
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
        let result_content = fs::read_to_string("result.nwk")?;
        print!("{}", result_content);
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
            '(' | ')' | '[' | ']' | ',' | ':' | ';' | ' ' | '/' | '\\' => '_',
            _ => c,
        })
        .collect()
}

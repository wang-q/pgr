//! Shared helpers for `pgr nwk` subcommands.

use std::collections::{BTreeMap, BTreeSet};

use anyhow::anyhow;
use clap::ArgMatches;
use pgr::libs::phylo::node::{Node, NodeId};
use pgr::libs::phylo::tree::Tree;
use regex::RegexBuilder;

/// Parse a `--lca` argument value as two comma-separated names.
/// Returns `(&str, &str)` to avoid allocation; bails if the input does not
/// contain exactly one comma delimiting two non-empty names.
pub(crate) fn parse_lca_pair(lca: &str) -> anyhow::Result<(&str, &str)> {
    let mut parts = lca.splitn(2, ',');
    let first = parts.next().unwrap_or("");
    let last = parts.next().unwrap_or("");
    if lca.matches(',').count() != 1 || first.is_empty() || last.is_empty() {
        return Err(anyhow!(
            "--lca requires exactly two comma-separated names, got: {}",
            lca
        ));
    }
    Ok((first, last))
}

/// Format a node's label with extra columns (`dup`, `taxid`, `species`, `full`).
pub(crate) fn format_label_columns(node: &Node, name: &str, columns: &[String]) -> String {
    let mut out = String::from(name);
    if columns.is_empty() {
        return out;
    }
    for column in columns {
        match column.as_str() {
            "dup" => out.push_str(&format!("\t{}", name)),
            "taxid" => out.push_str(&format!(
                "\t{}",
                node.get_property("T").map(|s| s.as_str()).unwrap_or("")
            )),
            "species" => out.push_str(&format!(
                "\t{}",
                node.get_property("S").map(|s| s.as_str()).unwrap_or("")
            )),
            "full" => {
                let comment = node
                    .properties
                    .as_ref()
                    .filter(|p| !p.is_empty())
                    .map(|p| {
                        let pairs: Vec<String> = p
                            .iter()
                            .map(|(k, v)| {
                                if v.is_empty() {
                                    format!(":{}", k)
                                } else {
                                    format!(":{}={}", k, v)
                                }
                            })
                            .collect();
                        format!("[&&NHX{}]", pairs.join(""))
                    })
                    .unwrap_or_default();
                out.push_str(&format!("\t{}", comment));
            }
            _ => {}
        }
    }
    out
}

/// Returns IDs of named nodes matching the name selection rules from CLI args.
pub(crate) fn match_names(tree: &Tree, args: &ArgMatches) -> anyhow::Result<BTreeSet<NodeId>> {
    // IDs with names
    let id_of: BTreeMap<_, _> = tree.get_name_id();

    // all matched IDs
    let mut ids = BTreeSet::new();

    // ids supplied by --node
    if args.contains_id("node") {
        let names = args
            .get_many::<String>("node")
            .ok_or_else(|| anyhow::anyhow!("missing --node values"))?;
        for name in names {
            if let Some(id) = id_of.get(name) {
                ids.insert(*id);
            }
        }
    }

    // ids supplied by --name-list
    if args.contains_id("name_list") {
        let file = args
            .get_one::<String>("name_list")
            .ok_or_else(|| anyhow::anyhow!("missing --name-list value"))?;
        for name in pgr::libs::io::read_names::<Vec<String>>(file)?.iter() {
            if let Some(id) = id_of.get(name) {
                ids.insert(*id);
            }
        }
    }

    // ids matched with --regex
    if args.contains_id("regex") {
        let regexes = args
            .get_many::<String>("regex")
            .ok_or_else(|| anyhow::anyhow!("missing --regex values"))?;
        for regex in regexes {
            let re = RegexBuilder::new(regex).case_insensitive(true).build()?;
            for (name, id) in id_of.iter() {
                if re.is_match(name) {
                    ids.insert(*id);
                }
            }
        }
    }

    // Default is printing all named nodes
    let is_all =
        !(args.contains_id("node") || args.contains_id("name_list") || args.contains_id("regex"));

    if is_all {
        ids = id_of.values().cloned().collect();
    }

    // Include all descendants of internal nodes
    let is_descendants = if args.try_contains_id("descendants").is_ok() {
        args.get_flag("descendants")
    } else {
        false
    };

    if is_descendants {
        let internal_ids: Vec<NodeId> = ids
            .iter()
            .filter(|&&id| tree.get_node(id).map(|n| !n.is_leaf()).unwrap_or(false))
            .copied()
            .collect();
        for id in &internal_ids {
            let subtree_ids = tree.get_subtree(id);
            for sid in &subtree_ids {
                if let Some(n) = tree.get_node(*sid) {
                    if n.name.is_some() {
                        ids.insert(*sid);
                    }
                }
            }
        }
    }

    Ok(ids)
}

/// Returns IDs of nodes matching the position selection rules from CLI args.
pub(crate) fn match_positions(tree: &Tree, args: &ArgMatches) -> BTreeSet<NodeId> {
    let mut skip_internal = if args.try_contains_id("internal").is_ok() {
        args.get_flag("internal")
    } else {
        false
    };
    let skip_leaf = if args.try_contains_id("leaf").is_ok() {
        args.get_flag("leaf")
    } else {
        false
    };

    let is_monophyly = if args.try_contains_id("monophyly").is_ok() {
        args.get_flag("monophyly")
    } else {
        false
    };

    if is_monophyly {
        skip_internal = true;
    }

    // all matched IDs
    let mut ids = BTreeSet::new();

    let Some(root_id) = tree.get_root() else {
        return ids;
    };
    let preorder_ids = tree.preorder(&root_id);

    preorder_ids.iter().for_each(|id| {
        if let Some(node) = tree.get_node(*id) {
            if node.is_leaf() && !skip_leaf {
                ids.insert(*id);
            }
            if !node.is_leaf() && !skip_internal {
                ids.insert(*id);
            }
        }
    });

    ids
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_lca_pair_valid() {
        assert_eq!(parse_lca_pair("a,b").unwrap(), ("a", "b"));
        assert_eq!(parse_lca_pair("foo,bar").unwrap(), ("foo", "bar"));
    }

    #[test]
    fn parse_lca_pair_invalid() {
        assert!(parse_lca_pair("a").is_err());
        assert!(parse_lca_pair("a,b,c").is_err());
        assert!(parse_lca_pair(",b").is_err());
        assert!(parse_lca_pair("a,").is_err());
        assert!(parse_lca_pair("").is_err());
    }
}

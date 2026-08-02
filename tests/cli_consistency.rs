#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::PgrCmd;
use std::collections::BTreeMap;
use std::path::Path;

/// Subcommand names listed by clap in the `Commands:` section of a help text.
fn commands_from_help(help: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut in_commands = false;
    for line in help.lines() {
        if line.starts_with("Commands:") {
            in_commands = true;
            continue;
        }
        if in_commands {
            if line.trim().is_empty() || line.starts_with("Options:") {
                break;
            }
            // Command lines start with exactly two spaces; wrapped description
            // lines are indented deeper and are skipped.
            let Some(rest) = line.strip_prefix("  ") else {
                continue;
            };
            if rest.starts_with(' ') {
                continue;
            }
            let Some(name) = rest.split_whitespace().next() else {
                continue;
            };
            if name != "help" {
                out.push(name.to_string());
            }
        }
    }
    out
}

/// Top-level command names mentioned in the `after_help` bullet list.
fn after_help_commands(help: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in help.lines() {
        let trimmed = line.trim();
        let Some(rest) = trimmed.strip_prefix("* ") else {
            continue;
        };
        let Some(name) = rest.split_whitespace().next() else {
            continue;
        };
        if trimmed.ends_with(':') {
            continue; // group header, e.g. "* Sequences:"
        }
        if name != "help" {
            out.push(name.to_string());
        }
    }
    out
}

#[test]
fn after_help_matches_registered_top_level_commands() {
    let (stdout, stderr) = PgrCmd::new().args(&["--help"]).run();
    assert!(stderr.is_empty(), "pgr --help wrote to stderr: {}", stderr);
    let mut registered = commands_from_help(&stdout);
    let mut mentioned = after_help_commands(&stdout);
    registered.sort();
    mentioned.sort();
    assert_eq!(
        mentioned, registered,
        "after_help command list drifted from the registered top-level commands; \
         update src/pgr.rs after_help when adding/removing commands"
    );
}

#[test]
fn every_top_level_command_has_docs() {
    let (stdout, _) = PgrCmd::new().args(&["--help"]).run();
    for cmd in commands_from_help(&stdout) {
        // Historical naming: `pgr 2bit` is documented as docs/twobit.md.
        let file = if cmd == "2bit" {
            "twobit.md".to_string()
        } else {
            format!("{}.md", cmd)
        };
        assert!(
            Path::new("docs").join(&file).is_file(),
            "missing docs/{} for top-level command `{}`",
            file,
            cmd
        );
    }
}

#[test]
fn every_docs_file_maps_to_a_registered_command() {
    let (stdout, _) = PgrCmd::new().args(&["--help"]).run();
    let top: Vec<String> = commands_from_help(&stdout);

    let mut subcommands: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for cmd in &top {
        let (sub_stdout, _) = PgrCmd::new().args(&[cmd.as_str(), "--help"]).run();
        subcommands.insert(cmd.clone(), commands_from_help(&sub_stdout));
    }

    for entry in std::fs::read_dir("docs").expect("docs/ must exist") {
        let path = entry.expect("read_dir entry").path();
        let Some(name) = path.file_name().and_then(|n| n.to_str()).map(String::from) else {
            continue;
        };
        let Some(stem) = name.strip_suffix(".md") else {
            continue; // skip directories (formats/, benchmark/)
        };
        if stem == "usage_examples" {
            continue;
        }

        if stem == "twobit" {
            assert!(
                top.contains(&"2bit".to_string()),
                "docs/twobit.md has no registered `pgr 2bit` command"
            );
            continue;
        }
        if top.iter().any(|c| c == stem) {
            continue; // top-level command doc
        }

        // Leaf doc in <command>-<subcommand>.md form (e.g. align-pgi).
        let Some((cmd, sub)) = stem.split_once('-') else {
            panic!(
                "docs/{}.md matches neither a top-level command nor a \
                 <command>-<subcommand> pair",
                name
            );
        };
        let subs = subcommands.get(cmd).unwrap_or_else(|| {
            panic!(
                "docs/{}.md references unknown top-level command `{}`",
                name, cmd
            )
        });
        assert!(
            subs.iter().any(|s| s == sub),
            "docs/{}.md references subcommand `{}` not registered under `{}`",
            name,
            sub,
            cmd
        );
    }
}

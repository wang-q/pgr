//! Consistency between `notes/references/*.md` and the reference index
//! (`notes/project-understanding.md` §11), mirroring
//! `tests/cli_consistency.rs` for commands/docs.
//!
//! A note that is not indexed is effectively undiscoverable (the
//! hv.md/hypergen.md duplication was caused by exactly this drift), so a
//! new or renamed reference note must add its `[[file.md]]` row here.

use std::fs;
use std::path::Path;

fn notes_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("notes")
}

/// Lines of the section whose header starts with `header`, up to the next
/// `## ` section.
fn section_lines<'a>(pu: &'a str, header: &str) -> Vec<&'a str> {
    let lines: Vec<&str> = pu.lines().collect();
    let start = lines
        .iter()
        .position(|l| l.starts_with(header))
        .expect(header);
    let end = lines[start + 1..]
        .iter()
        .position(|l| l.starts_with("## "))
        .map(|i| start + 1 + i)
        .unwrap_or(lines.len());
    lines[start..end].to_vec()
}

/// All `.md` paths relative to `notes/` (recursive).
fn all_md_files() -> Vec<String> {
    fn walk(base: &Path, dir: &Path, out: &mut Vec<String>) {
        for entry in fs::read_dir(dir).unwrap().flatten() {
            let path = entry.path();
            if path.is_dir() {
                walk(base, &path, out);
            } else if path.extension().is_some_and(|e| e == "md") {
                out.push(
                    path.strip_prefix(base)
                        .unwrap()
                        .to_string_lossy()
                        .into_owned(),
                );
            }
        }
    }
    let mut out = Vec::new();
    let base = notes_dir();
    walk(&base, &base, &mut out);
    out
}

/// Every reference note must be indexed in the §11 table.
#[test]
fn every_reference_note_is_indexed() {
    let refs = notes_dir().join("references");
    let mut files: Vec<String> = fs::read_dir(&refs)
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".md"))
        .collect();
    files.sort();
    let pu = fs::read_to_string(notes_dir().join("project-understanding.md")).unwrap();
    let section = section_lines(&pu, "## 11");
    let missing: Vec<&String> = files
        .iter()
        .filter(|f| !section.iter().any(|l| l.contains(&format!("[[{f}]]"))))
        .collect();
    assert!(
        missing.is_empty(),
        "notes/references/*.md missing from project-understanding.md §11: {missing:?}\n\
         add a [[file]] row when creating or renaming a reference note"
    );
}

/// Every `[[x.md]]` link in the §11 table must resolve under `notes/`.
#[test]
fn every_s11_link_resolves() {
    let pu = fs::read_to_string(notes_dir().join("project-understanding.md")).unwrap();
    let section = section_lines(&pu, "## 11");
    let all = all_md_files();
    for line in section {
        let mut rest = line;
        while let Some(i) = rest.find("[[") {
            rest = &rest[i + 2..];
            let Some(j) = rest.find("]]") else {
                break;
            };
            let target = &rest[..j];
            if target.ends_with(".md") {
                assert!(
                    all.iter().any(|f| {
                        f == target || Path::new(f).file_name().is_some_and(|b| b == target)
                    }),
                    "§11 link [[{target}]] has no notes/ file"
                );
            }
            rest = &rest[j + 2..];
        }
    }
}

//! Flat clustering output formatting.
//!
//! Shared formatting logic for clustering algorithms that produce
//! `Vec<Vec<usize>>` results (DBSCAN, K-Medoids, MCL, Connected Components).

/// Sort and format flat clustering results (indices into `names`).
///
/// Members within each cluster are sorted alphabetically by name; clusters
/// are sorted by size (descending) then by first member name. `rep_fn`
/// selects the representative index for each cluster. For "pair" format,
/// returning `None` skips that cluster. For "cluster" format, the
/// representative is placed in the first column if one is returned.
pub fn format_flat_clusters<F>(
    clusters: &mut Vec<Vec<usize>>,
    names: &[String],
    format: &str,
    rep_fn: F,
) -> anyhow::Result<String>
where
    F: Fn(&[usize]) -> Option<usize>,
{
    // Sort members within each cluster alphabetically by name.
    for c in clusters.iter_mut() {
        c.sort_by_key(|&idx| &names[idx]);
    }
    // Sort clusters: size desc, then first member name.
    clusters.sort_by(|a, b| match b.len().cmp(&a.len()) {
        std::cmp::Ordering::Equal => names[a[0]].cmp(&names[b[0]]),
        other => other,
    });

    let mut out = String::new();
    match format {
        "cluster" => {
            for c in clusters {
                let rep_idx = rep_fn(c);
                let mut members: Vec<&str> = c.iter().map(|&idx| names[idx].as_str()).collect();
                if let Some(rep) = rep_idx {
                    // Move the representative to the first column.
                    if let Some(pos) = c.iter().position(|&idx| idx == rep) {
                        if pos > 0 && pos < members.len() {
                            let rep_name = members.remove(pos);
                            members.insert(0, rep_name);
                        }
                    }
                }
                out.push_str(&members.join("\t"));
                out.push('\n');
            }
        }
        "pair" => {
            for c in clusters.iter() {
                if let Some(rep_idx) = rep_fn(c) {
                    let rep_name = &names[rep_idx];
                    for &member_idx in c {
                        out.push_str(rep_name);
                        out.push('\t');
                        out.push_str(&names[member_idx]);
                        out.push('\n');
                    }
                }
            }
        }
        _ => anyhow::bail!("unsupported output format: {}", format),
    }
    Ok(out)
}

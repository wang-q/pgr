use anyhow::{anyhow, bail, Context, Result};
use indexmap::IndexMap;
use std::collections::HashMap;

/// Parsed NRPS modules: per-module domain lists and module metadata.
pub struct NrpsData {
    /// Domains per module, in insertion order.
    pub modules: IndexMap<String, Vec<HashMap<String, String>>>,
    /// Module metadata: id, color, prev.
    pub module_info: IndexMap<String, HashMap<String, String>>,
}

/// Parse NRPS TSV content into modules with computed domain positions.
pub fn parse_nrps(content: &str, default_color: &str) -> Result<NrpsData> {
    let mut modules: IndexMap<String, Vec<HashMap<String, String>>> = IndexMap::new();
    let mut current_module = String::from("");
    let mut current_color = default_color.to_string();
    let mut module_count = 1;

    // module info
    let mut module_info: IndexMap<String, HashMap<String, String>> = IndexMap::new();
    let mut prev_module = String::from("origin");

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields[0] == "Module" {
            // Get module name from second column or generate default name
            current_module = if fields.len() > 1 && !fields[1].is_empty() {
                fields[1].to_string()
            } else {
                format!("M{}", module_count)
            };
            module_count += 1;

            // Get color from third column or use default
            current_color = if fields.len() > 2 {
                fields[2].to_string()
            } else {
                default_color.to_string()
            };

            // Initialize new module vector and info
            modules.insert(current_module.clone(), Vec::new());
            let info = HashMap::from([
                ("id".to_string(), current_module.clone()),
                ("color".to_string(), current_color.clone()),
                ("prev".to_string(), prev_module.clone()),
            ]);
            module_info.insert(current_module.clone(), info);
            prev_module = current_module.clone();
            continue;
        }

        let domain_type = fields[0].to_string();
        let text = if fields.len() > 1 {
            let raw_text = fields[1].to_string();
            if raw_text.starts_with("D-") || raw_text.starts_with("L-") {
                let (prefix, rest) = raw_text.split_at(2);
                format!("{{\\scriptsize {}}}{}", prefix, rest)
            } else {
                raw_text
            }
        } else {
            String::new()
        };
        let color = if fields.len() > 2 {
            fields[2].to_string()
        } else {
            current_color.clone()
        };

        let (dx_before, dx_after) = match domain_type.as_str() {
            "A" => (0.4, 0.4),
            "C" | "E" | "CE" | "M" => (0.4, 0.4),
            "T" => (0.2, 0.2),
            "Te" | "R" => (0.3, 0.3),
            other => bail!("unknown domain type: {}", other),
        };

        let domain_id = if let Some(domains) = modules.get(&current_module) {
            format!("{}-{}", current_module, domains.len() + 1)
        } else {
            format!("{}-1", current_module)
        };

        let pos = if let Some(domains) = modules.get(&current_module) {
            if domains.is_empty() {
                0.0
            } else {
                let last_domain = domains
                    .last()
                    .ok_or_else(|| anyhow!("empty domains while computing pos"))?;
                let last_pos: f64 = last_domain
                    .get("pos")
                    .ok_or_else(|| anyhow!("missing pos in last domain"))?
                    .parse()
                    .context("failed to parse last domain pos")?;
                let last_dx_after: f64 = last_domain
                    .get("dx_after")
                    .ok_or_else(|| anyhow!("missing dx_after in last domain"))?
                    .parse()
                    .context("failed to parse last domain dx_after")?;
                last_pos + last_dx_after + dx_before
            }
        } else {
            0.0
        };

        let domain = HashMap::from([
            ("type".to_string(), domain_type),
            ("text".to_string(), text),
            ("color".to_string(), color),
            ("dx_before".to_string(), dx_before.to_string()),
            ("dx_after".to_string(), dx_after.to_string()),
            ("id".to_string(), domain_id),
            ("pos".to_string(), pos.to_string()),
        ]);

        if let Some(domains) = modules.get_mut(&current_module) {
            domains.push(domain);
        }
    }

    Ok(NrpsData {
        modules,
        module_info,
    })
}

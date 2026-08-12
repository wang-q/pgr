use anyhow::{anyhow, Result};

/// Replace the content between `begin_anchor` and `end_anchor` (inclusive of begin, exclusive of end)
/// in `template` with `replacement`. Anchors are matched as plain substring finds.
pub fn replace_section(
    template: &mut String,
    begin_anchor: &str,
    end_anchor: &str,
    replacement: &str,
) -> Result<()> {
    let begin = template
        .find(begin_anchor)
        .ok_or_else(|| anyhow!("template anchor {} not found", begin_anchor))?;
    let end = template
        .find(end_anchor)
        .ok_or_else(|| anyhow!("template anchor {} not found", end_anchor))?;
    template.replace_range(begin..end, replacement);
    Ok(())
}

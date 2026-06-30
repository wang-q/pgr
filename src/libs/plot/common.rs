use anyhow::{anyhow, Context, Result};
use std::io::Write;

/// Get a string value from a tera::Context, replacing the common
/// `context.get(k).unwrap().as_str().unwrap()` pattern with a friendly error.
pub fn context_get_str<'a>(context: &'a tera::Context, key: &str) -> Result<&'a str> {
    context
        .get(key)
        .ok_or_else(|| anyhow!("missing context key: {}", key))?
        .as_str()
        .ok_or_else(|| anyhow!("context key {} is not a string", key))
}

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

/// Render `template` with `context` via Tera and write the result to `writer`.
pub fn render_and_write<W: Write>(
    template: &str,
    context: &tera::Context,
    writer: &mut W,
) -> Result<()> {
    let mut tera = tera::Tera::default();
    tera.add_raw_templates(vec![("t", template)])
        .context("failed to register tera template")?;
    let rendered = tera
        .render("t", context)
        .context("failed to render tera template")?;
    writer.write_all(rendered.as_ref())?;
    Ok(())
}

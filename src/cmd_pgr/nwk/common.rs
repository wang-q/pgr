//! Shared helpers for `pgr nwk` subcommands.

use anyhow::anyhow;

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

//! Shared helpers for `pgr pl` pipeline subcommands.

use std::io::BufRead;

/// Read chromosome names from a `chr.sizes` file (lines of `<chr>\t<size>`).
pub fn read_chr_names(sizes_file: &str) -> anyhow::Result<Vec<String>> {
    let mut chrs: Vec<String> = Vec::new();
    for line in std::io::BufReader::new(std::fs::File::open(sizes_file)?).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 2 {
            chrs.push(fields[0].to_string());
        }
    }
    Ok(chrs)
}

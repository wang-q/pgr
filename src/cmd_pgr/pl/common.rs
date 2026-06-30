//! Shared helpers for `pgr pl` pipeline subcommands.

use cmd_lib::*;
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

/// Resolve `path` to an absolute path string. `stdout` is passed through as-is.
pub fn abs_path_or_stdout(path: &str) -> anyhow::Result<String> {
    if path == "stdout" {
        Ok(path.to_string())
    } else {
        Ok(intspan::absolute_path(path)?.display().to_string())
    }
}

/// Run `Profex -z genome` per chromosome and write `.rg` files.
///
/// For each chromosome, runs `Profex -z genome <sn>` writing `prof.<sn>.txt`,
/// then scans lines with `re_prof` capturing `start` and `end` (1-based inclusive
/// in output). If `min_depth` is set and the regex has a `depth` capture group,
/// entries with depth below the threshold are skipped. Returns the list of
/// `prof.<sn>.rg` file names.
#[allow(unused_variables)]
pub fn run_profex_per_chr(
    chrs: &[String],
    re_prof: &regex::Regex,
    min_depth: Option<usize>,
) -> anyhow::Result<Vec<String>> {
    let mut rg_files = vec![];
    for (i, chr) in chrs.iter().enumerate() {
        let sn = i + 1;
        run_cmd!(
            Profex -z genome ${sn} > prof.${sn}.txt
        )?;

        let reader = pgr::reader(&format!("prof.{}.txt", sn))?;

        let rg_file = format!("prof.{}.rg", sn);
        let mut writer = pgr::writer(&rg_file)?;

        for line in std::io::BufReader::new(reader)
            .lines()
            .map_while(Result::ok)
        {
            let Some(caps) = re_prof.captures(&line) else {
                continue;
            };

            if let Some(min_d) = min_depth {
                if let Some(depth_str) = caps.name("depth") {
                    let depth: usize = depth_str.as_str().parse()?;
                    if depth < min_d {
                        continue;
                    }
                }
            }

            let start = caps["start"].parse::<usize>()? + 1;
            let end = caps["end"].parse::<usize>()? + 1;

            writer.write_fmt(format_args!("{}:{}-{}\n", chr, start, end))?;
        }
        rg_files.push(rg_file);
    }
    Ok(rg_files)
}

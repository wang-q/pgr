//! FasBlock multiz-style merging: align multiple FAS files against a common
//! reference and emit merged multi-species blocks per genomic window.

mod banded_align;
mod merge;
#[cfg(test)]
mod tests;
mod windows;

pub use merge::merge_window;

use crate::libs::fmt::fas::{FasBlock, FasEntry};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FasMultizMode {
    Core,
    Union,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FasMultizGapModel {
    Constant,
    Medium,
    Loose,
}

#[derive(Clone, Debug)]
pub struct FasMultizConfig {
    pub ref_name: String,
    pub radius: usize,
    pub min_width: usize,
    pub mode: FasMultizMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_score: i32,
    pub gap_model: FasMultizGapModel,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
    pub score_matrix: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Window {
    pub chr: String,
    pub start: u64,
    pub end: u64,
}

fn find_ref_entry<'a>(block: &'a FasBlock, ref_name: &str) -> Option<&'a FasEntry> {
    block
        .entries
        .iter()
        .zip(block.names.iter())
        .find_map(|(entry, name)| if name == ref_name { Some(entry) } else { None })
}

fn ref_overlaps_window(entry: &FasEntry, window: &Window) -> bool {
    let range = entry.range();
    if range.chr() != &window.chr {
        return false;
    }
    let start = *range.start() as u64;
    let end = *range.end() as u64;
    start < window.end && end > window.start
}

pub fn merge_fas_files(
    ref_name: &str,
    infiles: &[impl AsRef<Path>],
    windows: &[Window],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Vec<FasBlock>> {
    let mut blocks_per_input: Vec<Vec<FasBlock>> = Vec::new();

    for infile in infiles {
        let infile_str = infile.as_ref().to_str().ok_or_else(|| {
            anyhow::anyhow!("path is not valid UTF-8: {}", infile.as_ref().display())
        })?;
        let mut reader = crate::reader(infile_str)?;
        let mut blocks = Vec::new();

        loop {
            match crate::libs::fmt::fas::next_fas_block(&mut reader) {
                Ok(block) => blocks.push(block),
                Err(e) => {
                    if e.kind() == std::io::ErrorKind::UnexpectedEof {
                        break;
                    } else {
                        return Err(e.into());
                    }
                }
            }
        }

        blocks_per_input.push(blocks);
    }

    if windows.is_empty() {
        return Ok(Vec::new());
    }

    let mut merged_blocks = Vec::new();
    for window in windows {
        if let Some(block) = merge_window(ref_name, window, &blocks_per_input, cfg)? {
            merged_blocks.push(block);
        }
    }

    Ok(merged_blocks)
}

pub fn merge_fas_files_auto_windows(
    ref_name: &str,
    infiles: &[impl AsRef<Path>],
    cfg: &FasMultizConfig,
) -> anyhow::Result<Vec<FasBlock>> {
    let mut blocks_per_input: Vec<Vec<FasBlock>> = Vec::new();

    for infile in infiles {
        let infile_str = infile.as_ref().to_str().ok_or_else(|| {
            anyhow::anyhow!("path is not valid UTF-8: {}", infile.as_ref().display())
        })?;
        let mut reader = crate::reader(infile_str)?;
        let mut blocks = Vec::new();

        loop {
            match crate::libs::fmt::fas::next_fas_block(&mut reader) {
                Ok(block) => blocks.push(block),
                Err(e) => {
                    if e.kind() == std::io::ErrorKind::UnexpectedEof {
                        break;
                    } else {
                        return Err(e.into());
                    }
                }
            }
        }

        blocks_per_input.push(blocks);
    }

    let windows = windows::derive_windows_from_blocks(ref_name, &blocks_per_input, cfg);
    if windows.is_empty() {
        return Ok(Vec::new());
    }

    let mut merged_blocks = Vec::new();
    for window in &windows {
        if let Some(block) = merge_window(ref_name, window, &blocks_per_input, cfg)? {
            merged_blocks.push(block);
        }
    }

    Ok(merged_blocks)
}

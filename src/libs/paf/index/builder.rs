//! Index construction: in-memory, multi-reader, and lazy BGZF modes.

use super::{build_trees, insert_record, PafIndex, PafMetadata};
use crate::libs::paf::parser::{parse_paf, parse_paf_line};
use coitrees::Interval;
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::sync::Mutex;

impl PafIndex {
    /// Build an index from a single PAF reader (in-memory CIGAR).
    pub fn build<R: BufRead>(reader: R) -> std::io::Result<Self> {
        let records = parse_paf(reader)?;
        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();

        for rec in &records {
            insert_record(rec, &mut names, &mut by_target, &mut by_query, None)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
        }

        let trees = build_trees(by_target);
        let reverse_trees = build_trees(by_query);
        Ok(PafIndex {
            names,
            trees,
            reverse_trees,
            lazy_source: None,
            lazy_source_path: None,
        })
    }

    /// Build an index from multiple PAF readers (in-memory CIGAR).
    pub fn build_multi<R: BufRead>(readers: Vec<R>) -> std::io::Result<Self> {
        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();

        for reader in readers {
            for rec in &parse_paf(reader)? {
                insert_record(rec, &mut names, &mut by_target, &mut by_query, None).map_err(
                    |e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                )?;
            }
        }

        let trees = build_trees(by_target);
        let reverse_trees = build_trees(by_query);
        Ok(PafIndex {
            names,
            trees,
            reverse_trees,
            lazy_source: None,
            lazy_source_path: None,
        })
    }

    /// Build an index from a file path, using lazy CIGAR loading for BGZF files.
    ///
    /// For `.gz` files that are BGZF-compressed: records the BGZF virtual position
    /// of each PAF line and stores `CigarStore::Lazy(vpos)`. CIGAR is fetched
    /// on-demand during queries, reducing memory at build time.
    ///
    /// For non-BGZF files (plain text, regular gzip): falls back to in-memory
    /// build (`CigarStore::Owned`).
    pub fn build_from_path(path: &str) -> anyhow::Result<Self> {
        if path == "stdin" {
            return Ok(Self::build(crate::libs::io::reader(path)?)?);
        }

        let p = std::path::Path::new(path);
        if p.extension() == Some(std::ffi::OsStr::new("gz")) && crate::is_bgzf(path) {
            Ok(Self::build_lazy_bgzf(path)?)
        } else {
            Ok(Self::build(crate::libs::io::reader(path)?)?)
        }
    }

    /// Lazy BGZF build: store virtual positions, fetch CIGAR on demand.
    pub(crate) fn build_lazy_bgzf(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = bgzf::io::Reader::new(file);

        let mut names = IndexMap::new();
        let mut by_target: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut by_query: HashMap<u32, Vec<Interval<PafMetadata>>> = HashMap::new();
        let mut line = String::new();

        loop {
            let vpos = u64::from(reader.virtual_position());
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end_matches('\n');
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let rec = match parse_paf_line(trimmed) {
                Ok(r) => r,
                Err(e) => {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("{}: {}", e, trimmed),
                    ));
                }
            };
            insert_record(&rec, &mut names, &mut by_target, &mut by_query, Some(vpos))
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
        }

        let trees = build_trees(by_target);
        let reverse_trees = build_trees(by_query);
        Ok(PafIndex {
            names,
            trees,
            reverse_trees,
            lazy_source: Some(Mutex::new(bgzf::io::Reader::new(File::open(path)?))),
            lazy_source_path: Some(path.to_string()),
        })
    }

    /// Reopen the lazy source file (used after `load` to restore lazy mode).
    pub(crate) fn reopen_lazy_source(&mut self) -> std::io::Result<()> {
        if let Some(ref path) = self.lazy_source_path {
            let file = File::open(path)?;
            self.lazy_source = Some(Mutex::new(bgzf::io::Reader::new(file)));
        }
        Ok(())
    }
}

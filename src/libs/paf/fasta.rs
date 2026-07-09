use crate::libs::loc;
use indexmap::IndexMap;
use noodles_core::Position;
use noodles_fasta as fasta;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::num::NonZeroUsize;

/// Per-file FASTA record LRU cache size.
const FASTA_LRU_SIZE: usize = 8;

/// Load TSV mapping genome_name -> bgzf_fasta_path.
/// Lines starting with '#' are comments; blank lines are skipped.
pub fn load_fasta_tsv(path: &str) -> anyhow::Result<IndexMap<String, String>> {
    let reader = crate::libs::io::reader(path)?;
    let mut map = IndexMap::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        anyhow::ensure!(
            fields.len() >= 2,
            "invalid TSV line '{line}': expected 2 tab-separated columns"
        );
        let name = fields[0].to_string();
        let fasta_path = fields[1].to_string();
        if map.insert(name.clone(), fasta_path).is_some() {
            anyhow::bail!("duplicate genome name in TSV: {name}");
        }
    }
    Ok(map)
}

/// Validate that every name in the PAF index is present in the TSV mapping.
pub fn validate_tsv_covers_index(
    seq_to_file: &IndexMap<String, String>,
    idx: &crate::libs::paf::index::PafIndex,
) -> anyhow::Result<()> {
    let mut missing: Vec<&str> = idx
        .names
        .keys()
        .filter(|n| !seq_to_file.contains_key(*n))
        .map(|n| n.as_str())
        .collect();
    missing.sort_unstable();
    if !missing.is_empty() {
        anyhow::bail!(
            "FASTA TSV is missing {} genome(s) present in PAF index: {}",
            missing.len(),
            missing.join(", ")
        );
    }
    Ok(())
}

/// Load FASTA TSV, validate it covers every name in `idx`, and build a [`FastaStore`].
///
/// Shared by `pgr paf` subcommands that need sequence access (to-fas, to-gfa,
/// to-maf, to-vcf). Combines [`load_fasta_tsv`], [`validate_tsv_covers_index`],
/// and [`FastaStore::new`] into one step.
pub fn prepare_store(
    tsv_path: &str,
    idx: &crate::libs::paf::index::PafIndex,
) -> anyhow::Result<FastaStore> {
    let seq_to_file = load_fasta_tsv(tsv_path)?;
    validate_tsv_covers_index(&seq_to_file, idx)?;
    FastaStore::new(&seq_to_file)
}

/// Load all sequences referenced by a TSV into a HashMap.
/// Returns an empty map if `tsv_path` is None.
/// Bails if the TSV loads 0 entries.
pub fn load_all_seqs(tsv_path: Option<&str>) -> anyhow::Result<HashMap<String, Vec<u8>>> {
    let Some(tsv) = tsv_path else {
        return Ok(HashMap::new());
    };
    let seq_to_file = load_fasta_tsv(tsv)?;
    if seq_to_file.is_empty() {
        anyhow::bail!("--fasta-tsv loaded 0 entries: {}", tsv);
    }
    let mut store = FastaStore::new(&seq_to_file)?;
    let mut map = HashMap::with_capacity(seq_to_file.len());
    for name in seq_to_file.keys() {
        map.insert(name.clone(), store.fetch_full(name)?);
    }
    Ok(map)
}

/// One opened BGZF FASTA file with its .loc index and a per-name record cache.
pub struct FastaEntry {
    reader: loc::Input,
    loc_of: IndexMap<String, (u64, usize)>,
    cache: lru::LruCache<String, fasta::Record>,
}

/// Manages multiple BGZF FASTA files keyed by file path, with a name -> file
/// mapping so multiple genome names can share one file (multi-chrom).
pub struct FastaStore {
    files: HashMap<String, FastaEntry>,
    name_to_file: HashMap<String, String>,
}

impl FastaStore {
    /// Create a FastaStore from a `seq_name → bgzf_fasta_path` mapping.
    pub fn new(seq_to_file: &IndexMap<String, String>) -> anyhow::Result<Self> {
        let mut files = HashMap::new();
        let name_to_file: HashMap<String, String> = seq_to_file
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();

        // Unique file paths
        let unique_paths: HashSet<&String> = seq_to_file.values().collect();
        for path in unique_paths {
            let loc_file = format!("{path}.loc");
            if !std::path::Path::new(&loc_file).is_file() {
                loc::create_loc(path, &loc_file, true)?;
            }
            let loc_of = loc::load_loc(&loc_file)?;
            let reader = loc::Input::Bgzf(
                noodles_bgzf::io::indexed_reader::Builder::default().build_from_path(path)?,
            );
            let cache = lru::LruCache::new(
                NonZeroUsize::new(FASTA_LRU_SIZE).expect("FASTA_LRU_SIZE must be non-zero"),
            );
            files.insert(
                path.clone(),
                FastaEntry {
                    reader,
                    loc_of,
                    cache,
                },
            );
        }

        Ok(Self {
            files,
            name_to_file,
        })
    }

    /// Fetch sequence [start, end) (0-based, half-open) and the total sequence
    /// length. Caches the underlying FASTA record keyed by `name`.
    pub fn fetch_range(
        &mut self,
        name: &str,
        start: i32,
        end: i32,
    ) -> anyhow::Result<(Vec<u8>, usize)> {
        let path = self
            .name_to_file
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("sequence '{name}' not in FASTA store"))?;
        let entry = self
            .files
            .get_mut(path)
            .ok_or_else(|| anyhow::anyhow!("file '{path}' not opened"))?;

        if !entry.cache.contains(name) {
            let record = loc::fetch_record(&mut entry.reader, &entry.loc_of, name)?;
            entry.cache.put(name.to_string(), record);
        }
        let record = entry
            .cache
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("cache miss after insert for '{name}'"))?;
        let total_len = record.sequence().len();

        // Validate coordinates before constructing noodles positions.
        if start < 0 {
            anyhow::bail!("start position {start} is negative for '{name}'");
        }
        if end < 0 {
            anyhow::bail!("end position {end} is negative for '{name}'");
        }
        if start >= end {
            anyhow::bail!("empty range [{start},{end}) for '{name}'");
        }
        if end as usize > total_len {
            anyhow::bail!("end position {end} exceeds sequence length {total_len} for '{name}'");
        }

        // noodles Position is 1-based inclusive; our coords are 0-based half-open.
        let start_pos = Position::new(start as usize + 1)
            .ok_or_else(|| anyhow::anyhow!("invalid start position {start}"))?;
        let end_pos = Position::new(end as usize)
            .ok_or_else(|| anyhow::anyhow!("invalid end position {end}"))?;
        let slice = record
            .sequence()
            .slice(start_pos..=end_pos)
            .ok_or_else(|| anyhow::anyhow!("slice [{start},{end}) out of range for '{name}'"))?;

        Ok((slice.as_ref().to_vec(), total_len))
    }

    /// Fetch the full sequence bytes for a name.
    pub fn fetch_full(&mut self, name: &str) -> anyhow::Result<Vec<u8>> {
        let path = self
            .name_to_file
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("sequence '{name}' not in FASTA store"))?;
        let entry = self
            .files
            .get_mut(path)
            .ok_or_else(|| anyhow::anyhow!("file '{path}' not opened"))?;

        if !entry.cache.contains(name) {
            let record = loc::fetch_record(&mut entry.reader, &entry.loc_of, name)?;
            entry.cache.put(name.to_string(), record);
        }
        let record = entry
            .cache
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("cache miss after insert for '{name}'"))?;
        Ok(record.sequence()[..].to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::IndexMap;
    use std::io::Write;

    /// Write a single-sequence BGZF FASTA file and build its .gzi index.
    /// Returns the path to the `.fa.gz` file.
    fn write_bgzf_fasta(dir: &std::path::Path, name: &str, seq: &str) -> String {
        let path = dir.join(format!("{name}.fa.gz"));
        let file = std::fs::File::create(&path).unwrap();
        let mut writer = noodles_bgzf::io::Writer::new(file);
        writeln!(writer, ">{name}").unwrap();
        writeln!(writer, "{seq}").unwrap();
        writer.flush().unwrap();
        drop(writer);
        crate::libs::fmt::fa::build_gzi_index(path.to_str().unwrap()).unwrap();
        path.to_string_lossy().to_string()
    }

    fn build_store(dir: &std::path::Path, name: &str, seq: &str) -> FastaStore {
        let gz = write_bgzf_fasta(dir, name, seq);
        let mut map = IndexMap::new();
        map.insert(name.to_string(), gz);
        FastaStore::new(&map).unwrap()
    }

    #[test]
    fn test_fetch_range_valid() {
        let dir = tempfile::tempdir().unwrap();
        let mut store = build_store(dir.path(), "chr1", "ACGTACGTAC");
        let (seq, len) = store.fetch_range("chr1", 2, 8).unwrap();
        assert_eq!(len, 10);
        assert_eq!(seq, b"GTACGT");
    }

    #[test]
    fn test_fetch_range_end_exceeds_length() {
        let dir = tempfile::tempdir().unwrap();
        let mut store = build_store(dir.path(), "chr1", "ACGTACGTAC");
        let err = store.fetch_range("chr1", 0, 100).unwrap_err();
        let msg = format!("{err}");
        assert!(
            msg.contains("end position 100 exceeds sequence length 10"),
            "{msg}"
        );
    }

    #[test]
    fn test_fetch_range_negative_start() {
        let dir = tempfile::tempdir().unwrap();
        let mut store = build_store(dir.path(), "chr1", "ACGTACGTAC");
        let err = store.fetch_range("chr1", -1, 5).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("start position -1 is negative"), "{msg}");
    }

    #[test]
    fn test_fetch_range_empty_range() {
        let dir = tempfile::tempdir().unwrap();
        let mut store = build_store(dir.path(), "chr1", "ACGTACGTAC");
        let err = store.fetch_range("chr1", 5, 5).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("empty range [5,5)"), "{msg}");
    }

    #[test]
    fn test_fetch_range_start_greater_than_end() {
        let dir = tempfile::tempdir().unwrap();
        let mut store = build_store(dir.path(), "chr1", "ACGTACGTAC");
        let err = store.fetch_range("chr1", 8, 2).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("empty range [8,2)"), "{msg}");
    }
}

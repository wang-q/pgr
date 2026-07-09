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
            let cache = lru::LruCache::new(NonZeroUsize::new(FASTA_LRU_SIZE).unwrap());
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

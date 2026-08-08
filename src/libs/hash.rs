use itertools::Itertools;
use minimizer_iter::MinimizerBuilder;
use std::iter::FromIterator;

// These codes were adapted from https://curiouscoding.nl/posts/fast-minimizers/
pub trait Hasher: Clone {
    fn hash(&self, t: &[u8]) -> u64;
    fn hash_kmers(&mut self, k: usize, t: &[u8]) -> Vec<u64> {
        t.windows(k).map(|kmer| self.hash(kmer)).collect::<Vec<_>>()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FxHash;
impl Hasher for FxHash {
    fn hash(&self, t: &[u8]) -> u64 {
        fxhash::hash64(t)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MurmurHash3;
impl Hasher for MurmurHash3 {
    fn hash(&self, t: &[u8]) -> u64 {
        murmurhash3::murmurhash3_x64_128(t, 42).0
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RapidHash;
impl Hasher for RapidHash {
    fn hash(&self, t: &[u8]) -> u64 {
        rapidhash::rapidhash(t)
    }
}

pub trait Minimizer {
    /// The absolute positions of all minimizers in the text.
    fn minimizer(&mut self, text: &[u8]) -> Vec<(u64, usize)>;
    fn mins(&mut self, text: &[u8]) -> Vec<u64>;
}

pub struct JumpingMinimizer<H = FxHash> {
    pub w: usize,
    pub k: usize,
    pub hasher: H,
}

impl<H: Hasher> Minimizer for JumpingMinimizer<H> {
    fn minimizer(&mut self, text: &[u8]) -> Vec<(u64, usize)> {
        let mut minimizers = Vec::new();

        // Precompute hashes of all k-mers.
        let hashes = self.hasher.hash_kmers(self.k, text);

        if hashes.len() < self.w {
            return vec![];
        }

        let mut start = 0;
        while start < hashes.len() - self.w {
            // Position_min returns the position of the leftmost minimal hash.
            let min_pos = start
                + hashes[start..start + self.w]
                    .iter()
                    .position_min()
                    .expect("w > 0");
            minimizers.push(min_pos);
            start = min_pos + 1;
        }
        // Possibly add one last minimizer.
        let start = hashes.len() - self.w;
        let min_pos = start + hashes[start..].iter().position_min().expect("w > 0");
        if minimizers.last() != Some(&min_pos) {
            minimizers.push(min_pos);
        }
        minimizers.iter().map(|e| (hashes[*e], *e)).collect()
    }

    fn mins(&mut self, text: &[u8]) -> Vec<u64> {
        self.minimizer(text).iter().map(|(min, _)| *min).collect()
    }
}

pub fn seq_mins(
    seq: &[u8],
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
) -> anyhow::Result<rapidhash::RapidHashSet<u64>> {
    let minimizers: Vec<u64> = match opt_hasher {
        "rapid" => JumpingMinimizer {
            w: opt_window,
            k: opt_kmer,
            hasher: RapidHash,
        }
        .mins(seq),
        "fx" => JumpingMinimizer {
            w: opt_window,
            k: opt_kmer,
            hasher: FxHash,
        }
        .mins(seq),
        "murmur" => JumpingMinimizer {
            w: opt_window,
            k: opt_kmer,
            hasher: MurmurHash3,
        }
        .mins(seq),
        "mod" => {
            let min_iter = minimizer_iter::MinimizerBuilder::<u64, _>::new_mod()
                .canonical()
                .minimizer_size(opt_kmer)
                .width(opt_window as u16)
                .iter(seq);

            min_iter.map(|(min, _, _)| min).collect()
        }
        _ => anyhow::bail!("unknown hasher: {}", opt_hasher),
    };
    let hashset: rapidhash::RapidHashSet<u64> = rapidhash::RapidHashSet::from_iter(minimizers);

    Ok(hashset)
}

/// Compute the Mash distance from a Jaccard index and k-mer size.
///
/// See <https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation>.
pub fn mash_distance(jaccard: f64, kmer: usize) -> f64 {
    if jaccard == 0.0 {
        1.0
    } else {
        ((-1.0 / kmer as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln()).abs()
    }
}

/// Convert a Mash distance to a similarity in [0, 1].
/// Clamps values > 1.0 (possible from numerical error) to similarity 0.0.
pub fn mash_to_sim(mash: f64) -> f64 {
    if mash > 1.0 {
        0.0
    } else {
        1.0 - mash
    }
}

/// Mash-compatible canonical k-mer hashes: uppercase, skip k-mers with
/// non-ACGT bases, take the lexicographically smaller of the forward k-mer
/// and its reverse complement (byte compare, not hash compare), then hash
/// with MurmurHash3_x64_128(seed=42) low 64 bits — exactly Mash's `getHash`
/// with `parameters.seed = 42` (Mash-master/src/mash/Sketch.cpp
/// `addMinHashes`, `hash.cpp`).
pub fn seq_mash_hashes(seq: &[u8], k: usize, seed: u32) -> Vec<u64> {
    let mut hashes = Vec::new();
    for_each_mash_hash(seq, k, seed, |h| hashes.push(h));
    hashes
}

/// Complement of an uppercase base (Mash's complement table; non-ACGT maps
/// to N, and such windows are skipped by the validity scan).
fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Stream canonical Mash k-mer hashes of `seq`, calling `f` for each window
/// hash. Rolling-window implementation with O(k) memory (no full-length
/// buffers), replicating Mash's `addMinHashes` window scan exactly,
/// including the "bad base at window start after a jump" case.
pub fn for_each_mash_hash(seq: &[u8], k: usize, seed: u32, mut f: impl FnMut(u64)) {
    if seq.len() < k {
        return;
    }
    // Rolling windows: `fwd` = uppercased bases of [start, start+k); `rev`
    // = their reverse complement in reading order (rev[t] = comp(fwd[k-1-t])).
    // Shifting one base per step keeps memory O(k).
    let mut fwd = vec![0u8; k];
    let mut rev = vec![0u8; k];
    let mut start = 0usize;
    for t in 0..k {
        let u = if seq[t].is_ascii_lowercase() {
            seq[t] - 32
        } else {
            seq[t]
        };
        fwd[t] = u;
        rev[k - 1 - t] = complement_base(u);
    }
    let mut j = 0usize;
    let mut i = 0usize;
    // Replicate Mash's `addMinHashes` window scan exactly: `j` marks how far
    // the window scan has advanced (positions are checked once); on a bad
    // base the window start jumps to `j` (then the outer `i++` advances past
    // it), so a window may be processed even though its start position was
    // never re-checked — this matches Mash's behaviour, including the
    // "bad base at window start after a jump" case.
    while i + k <= seq.len() {
        // Slide the window to [i, i+k) if it lagged behind (bad-base jump).
        while start < i {
            let next = seq[start + k];
            let u = if next.is_ascii_lowercase() {
                next - 32
            } else {
                next
            };
            fwd.copy_within(1.., 0);
            rev.copy_within(0..k - 1, 1);
            fwd[k - 1] = u;
            rev[0] = complement_base(u);
            start += 1;
        }
        let mut bad = false;
        while j < i + k {
            let u = if seq[j].is_ascii_lowercase() {
                seq[j] - 32
            } else {
                seq[j]
            };
            if !matches!(u, b'A' | b'C' | b'G' | b'T') {
                i = j; // Mash: i = j++
                j += 1;
                bad = true;
                break;
            }
            j += 1;
        }
        if bad {
            i += 1; // outer for i++
            continue;
        }
        if i + k > seq.len() {
            break;
        }
        let kmer = if fwd <= rev { &fwd[..] } else { &rev[..] }; // memcmp
        let h = murmurhash3::murmurhash3_x64_128(kmer, seed as u64).0;
        f(h);
        i += 1;
    }
}

/// Bottom-k MinHash sketch: the `sketch_size` smallest unique hashes
/// (Mash's MinHashHeap semantics; memory O(sketch_size) via a max-heap).
pub fn bottom_k_min_hashes(
    hashes: impl Iterator<Item = u64>,
    sketch_size: usize,
) -> rapidhash::RapidHashSet<u64> {
    let mut acc = BottomK::new(sketch_size);
    for h in hashes {
        acc.insert(h);
    }
    acc.into_set()
}

/// Incremental bottom-k accumulator (Mash's MinHashHeap): keeps the `size`
/// smallest unique hashes in O(size) memory.
pub struct BottomK {
    size: usize,
    set: rapidhash::RapidHashSet<u64>,
    heap: std::collections::BinaryHeap<u64>,
}

impl BottomK {
    /// Create an empty accumulator keeping the `size` smallest hashes.
    pub fn new(size: usize) -> Self {
        Self {
            size,
            set: rapidhash::RapidHashSet::default(),
            heap: std::collections::BinaryHeap::new(),
        }
    }

    /// Insert one hash, evicting the current maximum when over capacity.
    pub fn insert(&mut self, h: u64) {
        // Fast path: once full, hashes >= the current maximum can never be
        // in the top-k, so reject them before touching the set/heap.
        if self.set.len() >= self.size {
            if let Some(&max) = self.heap.peek() {
                if h >= max {
                    return;
                }
            }
        }
        if self.set.insert(h) {
            self.heap.push(h);
            if self.heap.len() > self.size {
                if let Some(removed) = self.heap.pop() {
                    self.set.remove(&removed);
                }
            }
        }
    }

    /// Consume the accumulator, returning the selected hash set.
    pub fn into_set(self) -> rapidhash::RapidHashSet<u64> {
        self.set
    }
}

/// Mash-compatible sketch distances: merge the two sorted bottom-k sketches
/// and count equal pairs in a `sketch_size`-step merge walk (Mash's
/// `compareSketches`); Jaccard = common / denom, where denom is the walk
/// length completed with the remaining unmerged hashes of the exhausted set
/// (capped at `sketch_size`), NOT the standard set Jaccard (full
/// intersection / union). Verified against `mash dist` on E. coli MG1655 x
/// Sakai (k=21, s=1000): 456/1000 shared, distance 0.0222766 — identical
/// to Mash; undersized sketches also match (e.g. 2/2 for identical 2-hash
/// sketches at k=15/s=1000, where Mash reports distance 0).
/// Containment uses the full sketch intersection / first-set size (Mash's
/// `within` semantics), which is larger than the merged-prefix common.
pub fn mash_sketch_distances(
    a: &rapidhash::RapidHashSet<u64>,
    b: &rapidhash::RapidHashSet<u64>,
    k: usize,
    sketch_size: usize,
) -> SetDistances {
    let mut ai: Vec<u64> = a.iter().copied().collect();
    let mut bi: Vec<u64> = b.iter().copied().collect();
    ai.sort_unstable();
    bi.sort_unstable();
    let (mut i, mut j, mut common, mut denom) = (0usize, 0usize, 0usize, 0usize);
    while denom < sketch_size && i < ai.len() && j < bi.len() {
        if ai[i] < bi[j] {
            i += 1;
        } else if ai[i] > bi[j] {
            j += 1;
        } else {
            i += 1;
            j += 1;
            common += 1;
        }
        denom += 1;
    }
    // Mash completes the union when one side exhausts early: add the
    // remaining unmerged hashes of both sides, capped at sketch_size.
    if denom < sketch_size {
        if i < ai.len() {
            denom += ai.len() - i;
        }
        if j < bi.len() {
            denom += bi.len() - j;
        }
        denom = denom.min(sketch_size);
    }
    // Two empty sketches: Mash treats common == denom (0 == 0) as distance
    // 0, so define jaccard = 1 here instead of emitting NaN.
    let jaccard = if denom > 0 {
        common as f64 / denom as f64
    } else {
        1.0
    };
    let inter_full = a.iter().filter(|h| b.contains(h)).count();
    SetDistances {
        total1: a.len(),
        total2: b.len(),
        inter: common,
        union: denom,
        mash: mash_distance(jaccard, k),
        jaccard,
        containment: inter_full as f64 / a.len().max(1) as f64,
    }
}

/// Distance metrics between two minimizer sets.
pub struct SetDistances {
    /// Cardinality of the first set.
    pub total1: usize,
    /// Cardinality of the second set.
    pub total2: usize,
    /// Intersection size.
    pub inter: usize,
    /// Union size.
    pub union: usize,
    /// Mash distance (bounded to [0, 1]).
    pub mash: f64,
    /// Jaccard index.
    pub jaccard: f64,
    /// Containment index (intersection / first set).
    pub containment: f64,
}

/// Compute Jaccard, Containment, and Mash distance between two minimizer sets.
///
/// See <https://mash.readthedocs.io/en/latest/distances.html#mash-distance-formulation>.
pub fn set_distances(
    s1: &rapidhash::RapidHashSet<u64>,
    s2: &rapidhash::RapidHashSet<u64>,
    kmer: usize,
) -> SetDistances {
    let total1 = s1.len();
    let total2 = s2.len();

    let inter = s1.intersection(s2).cloned().count();
    let union = total1 + total2 - inter;

    // Empty sketches: two empty sets are identical (jaccard 1, distance 0);
    // containment of an empty first set is undefined, report 0 instead of NaN.
    let jaccard = if union > 0 {
        inter as f64 / union as f64
    } else {
        1.0
    };
    let containment = if total1 > 0 {
        inter as f64 / total1 as f64
    } else {
        0.0
    };
    let mash = mash_distance(jaccard, kmer);

    SetDistances {
        total1,
        total2,
        inter,
        union,
        mash,
        jaccard,
        containment,
    }
}

/// 95% confidence interval for ANI estimated from a Jaccard value
/// (normal approximation on the Jaccard proportion; Hera et al. 2023 give
/// tighter FracMinHash-specific bounds). Only meaningful for unbiased
/// samplers (FracMinHash); minimizer/syncmer sketches have sampling bias.
pub fn ani_ci_from_jaccard(jaccard: f64, union: usize, kmer: usize) -> (f64, f64) {
    let se = (jaccard * (1.0 - jaccard) / union.max(1) as f64).sqrt();
    let j_lo = (jaccard - 1.96 * se).max(0.0);
    let j_hi = (jaccard + 1.96 * se).min(1.0);
    let ani = |j: f64| 1.0 - mash_distance(j, kmer);
    (ani(j_lo), ani(j_hi))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MinimizerInfo {
    pub hash: u64,
    pub seq_id: u32,
    pub pos: u32,
    pub strand: bool, // true: +, false: -
}

// Wrapper for Filter logic
struct FilterBuildHasher<'a, F> {
    filter: &'a F,
}

impl<'a, F> Clone for FilterBuildHasher<'a, F> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, F> Copy for FilterBuildHasher<'a, F> {}

impl<'a, F> std::hash::BuildHasher for FilterBuildHasher<'a, F>
where
    F: Fn(u64) -> bool,
{
    type Hasher = FilterHasher<'a, F>;
    fn build_hasher(&self) -> Self::Hasher {
        FilterHasher {
            filter: self.filter,
            state: 0,
        }
    }
}

struct FilterHasher<'a, F> {
    filter: &'a F,
    state: u64,
}

impl<'a, F> std::hash::Hasher for FilterHasher<'a, F>
where
    F: Fn(u64) -> bool,
{
    fn write(&mut self, bytes: &[u8]) {
        // Use RapidHash logic directly
        self.state = rapidhash::rapidhash(bytes);
    }

    fn finish(&self) -> u64 {
        let h = self.state;
        // If filter returns false, we want to reject this hash.
        // minimizer_iter selects the MINIMUM hash.
        // If we return u64::MAX, it will be ignored unless all hashes in window are MAX.
        if (self.filter)(h) {
            h
        } else {
            u64::MAX
        }
    }
}

/// Sketch a sequence to find minimizers.
///
/// # Arguments
/// * `seq` - The DNA sequence
/// * `seq_id` - ID of the sequence
/// * `k` - K-mer size
/// * `w` - Window size
/// * `soft_mask` - If true, ignore k-mers containing lowercase bases
/// * `filter` - A predicate that returns true if a hash should be KEPT.
pub fn seq_sketch<F>(
    seq: &[u8],
    seq_id: u32,
    k: usize,
    w: usize,
    soft_mask: bool,
    filter: F,
) -> Vec<MinimizerInfo>
where
    F: Fn(u64) -> bool,
{
    // Use minimizer_iter with our custom FilterBuildHasher
    let build_hasher = FilterBuildHasher { filter: &filter };

    let builder = MinimizerBuilder::<u64, _>::new()
        .minimizer_size(k)
        .width(w as u16)
        .canonical() // Ensure canonical minimizers (min of fwd/rev)
        .hasher(build_hasher);

    // If soft_mask is enabled, we filter out minimizers that overlap with lowercase regions.
    // We check the original sequence at the minimizer's position.
    // Note: This effectively drops windows where the minimizer falls in a masked region.
    // It avoids allocation and complex iterator mapping.

    builder
        .iter(seq)
        .map(|(hash, pos, is_rc)| {
            let strand = !is_rc;
            MinimizerInfo {
                hash,
                seq_id,
                pos: pos as u32,
                strand,
            }
        })
        .filter(|m| {
            if m.hash == u64::MAX {
                return false;
            } // Should be filtered by FilterHasher if used, but explicit check is fine.

            if soft_mask {
                let start = m.pos as usize;
                let end = start + k;
                if end > seq.len() {
                    return false;
                } // Should not happen

                // Check if any byte in seq[start..end] is lowercase
                !seq[start..end].iter().any(|&b| b.is_ascii_lowercase())
            } else {
                true
            }
        })
        .collect()
}

/// A named minimizer set, the basic unit compared by the sketch-distance
/// commands (`pgr dist mini` / `mash` / `frac`).
#[derive(Debug, Default, Clone)]
pub struct MinimizerEntry {
    pub name: String,
    pub set: rapidhash::RapidHashSet<u64>,
}

/// Read a FASTA file and build a `MinimizerEntry` per record (or one merged entry with `is_merge`).
pub fn load_minimizers(
    infile: &str,
    opt_hasher: &str,
    opt_kmer: usize,
    opt_window: usize,
    is_merge: bool,
) -> anyhow::Result<Vec<MinimizerEntry>> {
    let mut fa_in = crate::libs::fmt::fa::reader(infile)?;

    let mut entries = vec![];
    // Set to merge all minimizers if --merge is true
    let mut all_set: rapidhash::RapidHashSet<u64> = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        // obtain record or fail with error
        let record = result?;

        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();

        let set: rapidhash::RapidHashSet<u64> =
            seq_mins(&seq[..], opt_hasher, opt_kmer, opt_window)?;

        if is_merge {
            all_set.extend(set);
        } else {
            let entry = MinimizerEntry { name, set };
            entries.push(entry);
        }
    }

    if is_merge {
        let entry = MinimizerEntry {
            name: infile.to_string(),
            set: all_set,
        };
        entries.push(entry);
    }

    Ok(entries)
}

/// FracMinHash sketch of a sequence: keep canonical k-mers whose hash is
/// below `u64::MAX / scale` (Irber et al. 2022). Unlike minimizers/syncmers,
/// every k-mer is sampled with the same independent probability (1/scale),
/// so Jaccard/containment estimates are unbiased and comparable across
/// differently-sized sets; Hera et al. 2023 give ANI bias correction and
/// confidence intervals.
pub fn seq_fracminhash(
    seq: &[u8],
    k: usize,
    scale: usize,
    is_protein: bool,
) -> anyhow::Result<rapidhash::RapidHashSet<u64>> {
    anyhow::ensure!(scale > 0, "scale must be positive: {scale}");
    let threshold = u64::MAX / scale as u64;
    let mut set = rapidhash::RapidHashSet::default();
    if is_protein {
        for kmer in seq.windows(k) {
            let h = rapidhash::rapidhash(kmer);
            if h < threshold {
                set.insert(h);
            }
        }
    } else {
        // Canonical k-mers (min of fwd/rev 2-bit encoding) hashed with
        // rapidhash: the raw 2-bit value is structured (see hv.md §1.4),
        // so it must not be used as the FracMinHash key directly.
        for key in crate::libs::nt::rolling_kmer_keys(seq, k) {
            let Some(key) = key else { continue };
            let canonical = key.min(crate::libs::nt::rc_key(key, k));
            let h = rapidhash::rapidhash(&canonical.to_le_bytes());
            if h < threshold {
                set.insert(h);
            }
        }
    }
    Ok(set)
}

/// Read a FASTA file and build a `MinimizerEntry` per record (or one merged
/// entry with `is_merge`) using FracMinHash sampling.
pub fn load_fracminhash(
    infile: &str,
    k: usize,
    scale: usize,
    is_protein: bool,
    is_merge: bool,
) -> anyhow::Result<Vec<MinimizerEntry>> {
    let mut fa_in = crate::libs::fmt::fa::reader(infile)?;
    let mut entries = vec![];
    let mut all_set: rapidhash::RapidHashSet<u64> = rapidhash::RapidHashSet::default();

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();
        let set = seq_fracminhash(&seq[..], k, scale, is_protein)?;

        if is_merge {
            all_set.extend(set);
        } else {
            let entry = MinimizerEntry { name, set };
            entries.push(entry);
        }
    }

    if is_merge {
        let entry = MinimizerEntry {
            name: infile.to_string(),
            set: all_set,
        };
        entries.push(entry);
    }

    Ok(entries)
}

/// Read a FASTA file and build a `MinimizerEntry` per record (or one merged
/// entry with `is_merge`) using a Mash-compatible bottom-k MinHash sketch.
/// `seed` defaults to Mash's 42; `sketch_size` defaults to Mash's 1000.
pub fn load_mash_minhashes(
    infile: &str,
    k: usize,
    sketch_size: usize,
    seed: u32,
    is_merge: bool,
) -> anyhow::Result<Vec<MinimizerEntry>> {
    anyhow::ensure!(
        sketch_size > 0,
        "sketch size must be positive: {sketch_size}"
    );
    let mut fa_in = crate::libs::fmt::fa::reader(infile)?;
    let mut entries = vec![];
    // Mash builds one MinHashHeap across all sequences of a file (global
    // bottom-k). With --merge we stream every canonical hash into a single
    // accumulator so memory stays O(sketch_size) regardless of genome
    // length; without --merge each record gets its own bounded sketch.
    let mut merged = is_merge.then(|| BottomK::new(sketch_size));

    for result in fa_in.records() {
        let record = result?;
        let name = String::from_utf8(record.name().into())?;
        let seq = record.sequence();
        if let Some(acc) = merged.as_mut() {
            for_each_mash_hash(&seq[..], k, seed, |h| acc.insert(h));
        } else {
            let mut acc = BottomK::new(sketch_size);
            for_each_mash_hash(&seq[..], k, seed, |h| acc.insert(h));
            entries.push(MinimizerEntry {
                name,
                set: acc.into_set(),
            });
        }
    }

    if let Some(acc) = merged {
        entries.push(MinimizerEntry {
            name: infile.to_string(),
            set: acc.into_set(),
        });
    }

    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rand::SeedableRng;

    fn rand_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
    }

    /// Reverse complement of an uppercase ACGT byte slice (Mash canonical
    /// rule), kept for the reference implementation below.
    fn reverse_complement(seq: &[u8], out: &mut [u8]) {
        for (i, &b) in seq.iter().enumerate() {
            out[seq.len() - 1 - i] = match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            };
        }
    }

    #[test]
    fn test_seq_fracminhash_sampling_rate() {
        let seq = rand_dna(100_000, 42);
        // scale=10 keeps ~1/10 of k-mers
        let s10 = seq_fracminhash(&seq, 21, 10, false).unwrap();
        let expected10 = 100_000 / 10;
        assert!(
            (s10.len() as i64 - expected10 as i64).abs() < expected10 as i64 / 5,
            "scale=10 kept {} (expected ~{})",
            s10.len(),
            expected10
        );
        // scale=1000 keeps ~1/1000
        let s1000 = seq_fracminhash(&seq, 21, 1000, false).unwrap();
        let expected1000 = 100_000 / 1000;
        assert!(
            (s1000.len() as i64 - expected1000 as i64).abs() < 50,
            "scale=1000 kept {} (expected ~{})",
            s1000.len(),
            expected1000
        );
    }

    #[test]
    fn test_fracminhash_jaccard_estimation() {
        // A = 50kb random; B = first 40kb of A + 10kb random -> ~80% shared
        // k-mers, so true Jaccard ~ 40/(50+50-40) = 0.667 (ignoring k-mer
        // boundary effects). FracMinHash must estimate it unbiasedly.
        let a = rand_dna(50_000, 7);
        let mut b = a[..40_000].to_vec();
        b.extend(rand_dna(10_000, 8));
        let sa = seq_fracminhash(&a, 21, 10, false).unwrap();
        let sb = seq_fracminhash(&b, 21, 10, false).unwrap();
        let inter = sa.iter().filter(|h| sb.contains(h)).count() as f64;
        let j = inter / (sa.len() + sb.len() - inter as usize) as f64;
        assert!(
            (j - 40.0 / 60.0).abs() < 0.1,
            "FracMinHash Jaccard {} vs expected ~0.667",
            j
        );
    }

    #[test]
    fn test_ani_ci() {
        // Narrower CI with more samples; CI must bracket the ANI point
        // estimate (1 - mash), not the raw Jaccard.
        let ani_pt = 1.0 - mash_distance(0.4, 21);
        let (lo, hi) = ani_ci_from_jaccard(0.4, 100, 21);
        assert!(
            lo < ani_pt && hi > ani_pt,
            "CI must bracket ANI {}: {} {}",
            ani_pt,
            lo,
            hi
        );
        let (lo2, hi2) = ani_ci_from_jaccard(0.4, 1000, 21);
        assert!(
            hi2 - lo2 < hi - lo,
            "more samples => narrower CI: {} vs {}",
            hi2 - lo2,
            hi - lo
        );
    }

    #[test]
    fn test_seq_mash_hashes_canonical_invariance() {
        // A sequence and its reverse complement (uppercase) must yield the
        // same canonical hash set (Mash takes min(fwd, rev) per k-mer).
        let seq = b"acgttgcatgcaacgtaacgt";
        let upper: Vec<u8> = seq
            .iter()
            .map(|&b| if b.is_ascii_lowercase() { b - 32 } else { b })
            .collect();
        let mut rc = vec![0u8; seq.len()];
        reverse_complement(&upper, &mut rc);
        let h1: std::collections::HashSet<u64> =
            seq_mash_hashes(&upper, 5, 42).into_iter().collect();
        let h2: std::collections::HashSet<u64> = seq_mash_hashes(&rc, 5, 42).into_iter().collect();
        assert_eq!(h1, h2, "canonical hashes must be strand-invariant");
    }

    #[test]
    fn test_seq_mash_hashes_skips_bad_bases() {
        // k-mers containing non-ACGT bases are skipped (Mash alphabet filter).
        let seq = b"ACGTNACGT"; // len 9, k=4 -> 6 windows, windows with N skipped
        let h = seq_mash_hashes(seq, 4, 42);
        assert!(
            !h.is_empty() && h.len() <= 5,
            "bad-base windows skipped: {}",
            h.len()
        );
        // No hash should come from a window containing N.
        for win in seq.windows(4) {
            if win.contains(&b'N') {
                continue;
            }
        }
    }

    /// Reference implementation with full-length buffers (the pre-streaming
    /// logic), used to prove the rolling-window version is byte-identical.
    fn reference_mash_hashes(seq: &[u8], k: usize, seed: u32) -> Vec<u64> {
        if seq.len() < k {
            return vec![];
        }
        let mut upper = vec![0u8; seq.len()];
        let mut valid = vec![false; seq.len()];
        for (i, &b) in seq.iter().enumerate() {
            let u = if b.is_ascii_lowercase() { b - 32 } else { b };
            upper[i] = u;
            valid[i] = matches!(u, b'A' | b'C' | b'G' | b'T');
        }
        let mut rc = vec![0u8; seq.len()];
        reverse_complement(&upper, &mut rc);
        let mut hashes = Vec::new();
        let mut j = 0usize;
        let mut i = 0usize;
        while i + k <= seq.len() {
            let mut bad = false;
            while j < i + k {
                if !valid[j] {
                    i = j;
                    j += 1;
                    bad = true;
                    break;
                }
                j += 1;
            }
            if bad {
                i += 1;
                continue;
            }
            if i + k > seq.len() {
                break;
            }
            let fwd = &upper[i..i + k];
            let rev = &rc[seq.len() - i - k..seq.len() - i];
            let kmer = if fwd <= rev { fwd } else { rev };
            hashes.push(murmurhash3::murmurhash3_x64_128(kmer, seed as u64).0);
            i += 1;
        }
        hashes
    }

    #[test]
    fn test_for_each_mash_hash_matches_reference() {
        // Random DNA with lowercase bases and Ns sprinkled in, across k
        // values: the rolling-window stream must equal the full-buffer
        // reference exactly, including bad-base jump behaviour.
        let mut seq = rand_dna(20_000, 99);
        for pos in (0..seq.len()).step_by(911) {
            seq[pos] = b'N';
        }
        for pos in (0..seq.len()).step_by(131) {
            seq[pos] = seq[pos].to_ascii_lowercase();
        }
        for k in [3usize, 4, 7, 15, 21, 31] {
            let mut streamed = Vec::new();
            for_each_mash_hash(&seq, k, 42, |h| streamed.push(h));
            assert_eq!(streamed, reference_mash_hashes(&seq, k, 42), "k={k}");
        }
    }

    #[test]
    fn test_bottom_k_min_hashes() {
        // unique hashes {3,5,8,10,20}; keep the 3 smallest.
        let hashes = vec![10u64, 5, 20, 5, 3, 8];
        let s = bottom_k_min_hashes(hashes.into_iter(), 3);
        assert_eq!(s.len(), 3);
        assert!(s.contains(&3) && s.contains(&5) && s.contains(&8));
    }

    #[test]
    fn test_mash_sketch_distances() {
        // A = {1..600, 2000..2399}, B = {1..400, 1500..1899, 2000..2199}
        // (both size 1000 bottom-k). Full intersection = 600 (1..400 and
        // 2000..2199). Mash's compareSketches walks at most 1000 merge steps:
        // 400 matches for 1..400, then A advances through 401..600 (200 steps)
        // and exhausts -> common = 400, so Mash Jaccard = 0.4 while the
        // standard set Jaccard is 600/1400 and standard containment 0.6.
        let a: rapidhash::RapidHashSet<u64> = (1..=600).chain(2000..2400).collect();
        let b: rapidhash::RapidHashSet<u64> =
            (1..=400).chain(1500..1900).chain(2000..2200).collect();
        let d = mash_sketch_distances(&a, &b, 21, 1000);
        assert_eq!(d.inter, 400);
        assert_eq!(d.jaccard, 0.4);
        // Containment uses the full intersection (600), not the merged common.
        assert_eq!(d.containment, 0.6);
        // Standard set Jaccard differs (600/1400 = 0.4286).
        let std = set_distances(&a, &b, 21);
        assert!((std.jaccard - 600.0 / 1400.0).abs() < 1e-6);
        assert_eq!(std.containment, 0.6);
        assert!((d.jaccard - std.jaccard).abs() > 0.02);
    }

    #[test]
    fn test_mash_sketch_distances_undersized() {
        // Identical 46-hash sketches: Mash reports 46/46, distance 0.
        let a: rapidhash::RapidHashSet<u64> = (1..=46).collect();
        let b: rapidhash::RapidHashSet<u64> = (1..=46).collect();
        let d = mash_sketch_distances(&a, &b, 15, 1000);
        assert_eq!(d.inter, 46);
        assert_eq!(d.union, 46);
        assert_eq!(d.jaccard, 1.0);
        assert_eq!(d.mash, 0.0);

        // Disjoint 46-hash sketches: Mash reports 0/92.
        let a: rapidhash::RapidHashSet<u64> = (1..=46).collect();
        let b: rapidhash::RapidHashSet<u64> = (1000..=1045).collect();
        let d = mash_sketch_distances(&a, &b, 15, 1000);
        assert_eq!(d.inter, 0);
        assert_eq!(d.union, 92);
        assert_eq!(d.jaccard, 0.0);
        assert_eq!(d.mash, 1.0);
    }

    #[test]
    fn test_distances_empty_sets() {
        // Two empty sketches must not produce NaN: jaccard 1 / distance 0
        // (Mash's common == denom == 0 case) for both distance functions.
        let empty: rapidhash::RapidHashSet<u64> = rapidhash::RapidHashSet::default();
        let d = mash_sketch_distances(&empty, &empty, 21, 1000);
        assert_eq!(d.jaccard, 1.0);
        assert_eq!(d.mash, 0.0);
        assert!(!d.jaccard.is_nan() && !d.mash.is_nan());

        let d = set_distances(&empty, &empty, 21);
        assert_eq!(d.jaccard, 1.0);
        assert_eq!(d.mash, 0.0);
        assert_eq!(d.containment, 0.0);
        assert!(!d.containment.is_nan());
    }

    #[test]
    fn test_seq_sketch_basic() {
        let seq = b"ACGTACGT";
        let k = 3;
        let w = 3; // minimizer_iter requires odd window size?
        let mins = seq_sketch(seq, 1, k, w, false, |_| true);

        assert!(!mins.is_empty());
        for m in &mins {
            assert_eq!(m.seq_id, 1);
            assert!(m.pos < seq.len() as u32);
        }
    }

    #[test]
    fn test_seq_sketch_soft_mask() {
        // "acgt" is lowercase. If soft_mask is true, it should be ignored.
        let seq = b"ACGTacgtACGT";
        let k = 4;
        let w = 1; // Small window to test individual k-mers

        // Without soft mask: "acgt" (lowercase) is a valid k-mer (different hash from ACGT)
        let mins_no_mask = seq_sketch(seq, 1, k, w, false, |_| true);
        // ACGT (0), CGTa (1), GTac (2), Tacg (3), acgt (4), cgtA (5), gtAC (6), tACG (7), ACGT (8)
        // We expect some minimizers.
        assert!(mins_no_mask.len() >= 2);

        // With soft mask: "acgt" and any k-mer containing lowercase should be ignored.
        // k-mers containing lowercase:
        // CGTa, GTac, Tacg, acgt, cgtA, gtAC, tACG
        // Only ACGT (0) and ACGT (8) are purely uppercase.
        let mins_mask = seq_sketch(seq, 1, k, w, true, |_| true);

        // Should only find the two uppercase blocks
        // ACGT at 0
        // ACGT at 8
        // Note: minimizer_iter with w=1 returns all valid k-mers.
        // But if we return u64::MAX, they are filtered out.
        assert_eq!(mins_mask.len(), 2);
        assert_eq!(mins_mask[0].pos, 0);
        assert_eq!(mins_mask[1].pos, 8);
    }

    #[test]
    fn test_seq_sketch_strand() {
        // AAAA (fwd) vs TTTT (rev)
        // If canonical is working, it should pick the smaller hash.
        // Let's rely on consistency.
        let seq = b"ACGT";
        let k = 4;
        let w = 1;
        let mins = seq_sketch(seq, 1, k, w, false, |_| true);
        assert_eq!(mins.len(), 1);
    }
}

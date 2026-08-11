//! Tadpole-compatible contig assembly (contigMode).

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use crate::libs::fq::qual::{from_phred, to_phred};
use crate::libs::fq::tadpole::{
    argmax2, base_code, base_defined, number_to_base, second_highest_position, Kmer, TadpoleTable,
};
use crate::libs::nt::rev_comp;
use anyhow::Result;
use std::collections::HashSet;
use std::io::Write;

/// Result codes from `extendToRight` (ShaveObject).
const DEAD_END: i32 = 1;
const LOOP: i32 = 7;
const BAD_OWNER: i32 = 11;
const BAD_SEED: i32 = 12;
const F_BRANCH: i32 = 17;
const B_BRANCH: i32 = 18;
const D_BRANCH: i32 = 19;

/// Assembly options with tadpole.sh defaults.
#[derive(Debug, Clone)]
pub struct AssembleOptions {
    /// K-mer length.
    pub k: usize,
    /// Ignore k-mers below this error-free probability.
    pub min_prob: f32,
    /// Minimum k-mer depth to seed a contig.
    pub min_count_seed: usize,
    /// Minimum k-mer depth to continue an extension.
    pub min_count_extend: usize,
    /// Minimum added bases past the seed for a contig to be kept.
    pub min_extension: usize,
    /// Minimum contig length; `0` selects the tadpole auto value max(124, 2k).
    pub min_contig_len: usize,
    /// Minimum k-mer coverage for a contig.
    pub min_coverage: f32,
    /// Branch ratio at high depth (branchmult1).
    pub branch_mult1: f32,
    /// Branch ratio at low depth (branchmult2).
    pub branch_mult2: f32,
    /// Second-highest depth considered "low" (branchlower).
    pub branch_lower_const: usize,
    /// Number of seeding passes (contigPasses).
    pub contig_passes: usize,
    /// Seeding pass multiplier (contigPassMult).
    pub contig_pass_mult: f64,
    /// Merge parallel paths in the contig graph (Tadpole popbubbles).
    pub pop_bubbles: bool,
    /// Append `L:` links to unitig FASTA headers (BCALM format).
    pub emit_links: bool,
    /// Emit a GFA graph instead of FASTA.
    pub emit_gfa: bool,
}

impl Default for AssembleOptions {
    fn default() -> Self {
        Self {
            k: 31,
            min_prob: 0.5,
            min_count_seed: 3,
            min_count_extend: 2,
            min_extension: 2,
            min_contig_len: 0,
            min_coverage: 1.0,
            branch_mult1: 20.0,
            branch_mult2: 3.0,
            branch_lower_const: 3,
            contig_passes: 16,
            contig_pass_mult: 1.7,
            pop_bubbles: true,
            emit_links: false,
            emit_gfa: false,
        }
    }
}

impl AssembleOptions {
    fn resolved_min_contig_len(&self) -> usize {
        if self.min_contig_len > 0 {
            self.min_contig_len
        } else {
            (124).max(2 * self.k)
        }
    }
}

/// `Tadpole.isJunction(max, second)`: depth-ratio branch detection.
fn is_junction(max: u32, second: u32, opts: &AssembleOptions) -> bool {
    if second < 1
        || (second as f32) * opts.branch_mult1 < max as f32
        || (second <= opts.branch_lower_const as u32
            && (max as f32)
                >= (opts.min_count_extend as f32).max(second as f32 * opts.branch_mult2))
    {
        return false;
    }
    true
}

/// One assembled contig.
#[derive(Clone)]
struct Contig {
    bases: Vec<u8>,
    id: usize,
    coverage: f32,
    min_cov: usize,
    max_cov: usize,
    left_code: i32,
    right_code: i32,
    left_ratio: f32,
    right_ratio: f32,
    used: bool,
    associate: bool,
    flipped: bool,
    left_edges: Vec<EdgeRef>,
    right_edges: Vec<EdgeRef>,
}

/// Directed edge between two contigs (assemble.Edge).
#[derive(Clone)]
struct Edge {
    origin: usize,
    destination: usize,
    length: usize,
    /// bit 0: source connects on its right; bit 1: destination on its right.
    orientation: u8,
    depth: u32,
    bases: Vec<u8>,
}

impl Edge {
    fn dest_right(&self) -> bool {
        self.orientation & 2 == 2
    }

    fn flip_source(&mut self) {
        self.bases = rev_comp(&self.bases).collect();
        self.orientation ^= 1;
    }

    fn flip_dest(&mut self) {
        self.orientation ^= 2;
    }
}

/// Assembly statistics.
#[derive(Debug, Default, Clone)]
pub struct AssembleStats {
    pub reads_in: u64,
    pub bases_in: u64,
    pub contigs_built: u64,
    pub bases_built: u64,
    pub longest_contig: usize,
}

/// Assembles reads into contigs via the k-mer graph (tadpole contigMode).
///
/// Mirrors `Tadpole.process2(contigMode)`: canonical k-mer counting, then
/// multi-pass seeding with decreasing depth thresholds, bidirectional greedy
/// extension with ownership, and deterministic longest-first output.
pub fn assemble<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &AssembleOptions,
) -> Result<AssembleStats> {
    anyhow::ensure!(
        opts.k >= 1,
        "k-mer length must be at least 1, got {}",
        opts.k
    );

    let records = read_records(infiles)?;

    // Pass 2: count k-mers from the canonicalized (phred) qualities.
    let reads: Vec<(Vec<u8>, Vec<u8>)> = records
        .iter()
        .map(|r| {
            (
                r.sequence().to_vec(),
                to_phred(r.sequence(), r.quality_scores()),
            )
        })
        .collect();
    let table = TadpoleTable::build(&reads, opts.k, opts.min_prob);
    let bases_in: u64 = reads.iter().map(|(s, _)| s.len() as u64).sum();

    // Pass 3: multi-pass seeding and contig building (BuildThread.run).
    let mut claimed: HashSet<Kmer> = HashSet::new();
    let mut contigs: Vec<Contig> = Vec::new();
    let mut id_counter = 0usize;
    for i in (1..opts.contig_passes).rev() {
        let threshold = pass_threshold(opts, i);
        scan_table(
            &table,
            threshold,
            opts,
            &mut claimed,
            &mut contigs,
            &mut id_counter,
        );
    }
    scan_table(
        &table,
        opts.min_count_seed,
        opts,
        &mut claimed,
        &mut contigs,
        &mut id_counter,
    );

    // Contig graph + bubble popping (Tadpole.processContigs/popBubbles);
    // with --no-bubbles the pre-pop contigs are kept and only sorted and
    // renumbered.
    if opts.pop_bubbles {
        process_contigs(&mut contigs, &table, opts);
        pop_bubbles(&mut contigs, opts);
    } else {
        finalize_contigs(&mut contigs);
    }

    let mut stats = AssembleStats {
        reads_in: records.len() as u64,
        bases_in,
        ..AssembleStats::default()
    };
    let min_contig_len = opts.resolved_min_contig_len();
    for c in &contigs {
        if c.bases.len() >= min_contig_len {
            write_contig(out, c)?;
            stats.contigs_built += 1;
            stats.bases_built += c.bases.len() as u64;
            stats.longest_contig = stats.longest_contig.max(c.bases.len());
        }
    }
    Ok(stats)
}

/// Reads all records from 1 interleaved or 2 paired files, canonicalizing
/// qualities like BBTools (shared by the contig and unitig modes).
fn read_records(infiles: &[String]) -> Result<Vec<SeqRecord>> {
    let mut records: Vec<SeqRecord> = Vec::new();
    let mut reader1 = SeqReader::new(&infiles[0])?;
    let mut reader2 = if infiles.len() > 1 {
        Some(SeqReader::new(&infiles[1])?)
    } else {
        None
    };
    let mut rec = SeqRecord::new();
    loop {
        if !reader1.read_record(&mut rec)? {
            break;
        }
        canonicalize_quality(&mut rec);
        records.push(rec.clone());
        if let Some(r) = reader2.as_mut() {
            if !r.read_record(&mut rec)? {
                anyhow::bail!("unpaired trailing read in {}", infiles[0]);
            }
            canonicalize_quality(&mut rec);
            records.push(rec.clone());
        } else if !reader1.read_record(&mut rec)? {
            anyhow::bail!("unpaired trailing read in {}", infiles[0]);
        } else {
            canonicalize_quality(&mut rec);
            records.push(rec.clone());
        }
    }
    Ok(records)
}

/// One maximal unitig (non-branching path; BCALM `graph3` semantics).
#[derive(Clone)]
struct Unitig {
    bases: Vec<u8>,
    id: usize,
    coverage: f32,
    min_cov: usize,
    max_cov: usize,
    /// The k-mer path closes back on itself (a circular contig).
    circular: bool,
}

/// Assembles reads into maximal unitigs instead of seeded contigs.
///
/// BCALM-style compaction (`ograph.cpp` `graph3`): every solid k-mer
/// (count >= `min_count_seed`) compresses into its unique non-branching
/// path. A k-mer extends only while it has exactly one solid successor
/// whose own predecessor is also unique; parallel paths stay separate
/// (no bubble popping), and the result is independent of scan order.
pub fn assemble_unitigs<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &AssembleOptions,
) -> Result<AssembleStats> {
    anyhow::ensure!(
        opts.k >= 1,
        "k-mer length must be at least 1, got {}",
        opts.k
    );
    let records = read_records(infiles)?;
    let reads: Vec<(Vec<u8>, Vec<u8>)> = records
        .iter()
        .map(|r| {
            (
                r.sequence().to_vec(),
                to_phred(r.sequence(), r.quality_scores()),
            )
        })
        .collect();
    let table = TadpoleTable::build(&reads, opts.k, opts.min_prob);
    let bases_in: u64 = reads.iter().map(|(s, _)| s.len() as u64).sum();

    let mut unitigs = build_unitigs(&table, opts);
    unitigs.sort_by(unitig_cmp);
    let links = if opts.emit_links || opts.emit_gfa {
        compute_links(&unitigs, opts.k)
    } else {
        vec![Vec::new(); unitigs.len()]
    };
    let mut stats = AssembleStats {
        reads_in: records.len() as u64,
        bases_in,
        ..AssembleStats::default()
    };
    if opts.emit_gfa {
        writeln!(out, "H\tVN:Z:1.0\tks:i:{}", opts.k)?;
    }
    let min_len = opts.resolved_min_contig_len();
    for (i, u) in unitigs.iter_mut().enumerate() {
        u.id = i;
        if u.bases.len() >= min_len {
            if opts.emit_gfa {
                writeln!(out, "S\t{}\t{}", u.id, String::from_utf8_lossy(&u.bases))?;
                for l in &links[i] {
                    writeln!(
                        out,
                        "L\t{}\t{}\t{}\t{}\t{}M",
                        u.id,
                        if l.from_rc { '-' } else { '+' },
                        l.to,
                        if l.to_rc { '-' } else { '+' },
                        opts.k.saturating_sub(1),
                    )?;
                }
            } else {
                write_unitig(
                    out,
                    u,
                    if opts.emit_links {
                        Some(&links[i])
                    } else {
                        None
                    },
                )?;
            }
            stats.contigs_built += 1;
            stats.bases_built += u.bases.len() as u64;
            stats.longest_contig = stats.longest_contig.max(u.bases.len());
        }
    }
    Ok(stats)
}

/// One directed unitig link, starting at the owning unitig's right end.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Link {
    to: usize,
    from_rc: bool,
    to_rc: bool,
}

/// Computes links between unitigs sharing an endpoint (k-1)-mer (BCALM
/// LinkTigs semantics, simplified to actual-sequence matching).
///
/// Direction rule for the source's right-end (k-1)-mer `r` meeting the
/// target's end `a` (all in actual sequence):
/// * `a == r` on the target's left end -> `+`/`+` (3' -> 5');
/// * `a == rc(r)` on the target's left end -> `+`/`-`;
/// * 3'-3' or 5'-5' meetings are expressed on the reverse strand (`-`).
fn compute_links(unitigs: &[Unitig], k: usize) -> Vec<Vec<Link>> {
    if k < 2 || unitigs.is_empty() {
        return vec![Vec::new(); unitigs.len()];
    }
    // Endpoint (k-1)-mers indexed by canonical form; bool = right end.
    let mut idx: HashMap<Kmer, Vec<(usize, bool)>> = HashMap::new();
    for (i, u) in unitigs.iter().enumerate() {
        for right in [false, true] {
            let km = end_kmer1(&u.bases, k, right);
            idx.entry(km.canonical()).or_default().push((i, right));
        }
    }
    let mut links = vec![Vec::new(); unitigs.len()];
    for (i, u) in unitigs.iter().enumerate() {
        let r = end_kmer1(&u.bases, k, true);
        let rc_r = r.rc();
        let Some(candidates) = idx.get(&r.canonical()) else {
            continue;
        };
        for &(j, j_right) in candidates {
            if j == i {
                continue;
            }
            let actual = end_kmer1(&unitigs[j].bases, k, j_right);
            let (from_rc, to_rc) = if !j_right {
                if actual.cmp_bases(&r) == std::cmp::Ordering::Equal {
                    (false, false)
                } else if actual.cmp_bases(&rc_r) == std::cmp::Ordering::Equal {
                    (false, true)
                } else {
                    continue;
                }
            } else if actual.cmp_bases(&r) == std::cmp::Ordering::Equal {
                (true, true)
            } else if actual.cmp_bases(&rc_r) == std::cmp::Ordering::Equal {
                (true, false)
            } else {
                continue;
            };
            links[i].push(Link {
                to: j,
                from_rc,
                to_rc,
            });
        }
    }
    for l in &mut links {
        l.sort_unstable();
        l.dedup();
    }
    links
}

/// The (k-1)-mer at a unitig end as a `Kmer` (actual sequence).
fn end_kmer1(bases: &[u8], k: usize, right: bool) -> Kmer {
    let mut km = Kmer::new(k - 1);
    if right {
        for &b in &bases[bases.len() - (k - 1)..] {
            km.push_right(base_code(b));
        }
    } else {
        for &b in &bases[..k - 1] {
            km.push_right(base_code(b));
        }
    }
    km
}

/// Compresses every solid k-mer into its maximal unitig (order-independent).
fn build_unitigs(table: &TadpoleTable, opts: &AssembleOptions) -> Vec<Unitig> {
    let k = opts.k;
    let threshold = opts.min_count_seed as u32;
    let mut visited: HashSet<Kmer> = HashSet::new();
    let mut unitigs = Vec::new();
    for (seed, count) in table.sorted_entries().iter() {
        if *count < threshold || visited.contains(seed) {
            continue;
        }
        // `base_at(0)` is the 3' end (last base pushed); rebuild 5'->3'.
        let mut bb: Vec<u8> = (0..k)
            .map(|i| number_to_base(seed.base_at(k - 1 - i)))
            .collect();
        visited.insert(seed.clone());
        // Extend right while the path stays non-branching.
        let mut circular = false;
        let mut kmer = rightmost_kmer(&bb, k);
        while let Some(b) = unique_solid_out(&kmer, table, threshold) {
            let mut next = kmer.clone();
            next.push_right(b);
            let canon = next.canonical();
            if unique_solid_in(&next, table, threshold) != 1 {
                break;
            }
            if visited.contains(&canon) {
                circular = true;
                break;
            }
            bb.push(number_to_base(b));
            visited.insert(canon);
            kmer = next;
        }
        // Extend left by reverse-complementing and extending right.
        let mut rc: Vec<u8> = rev_comp(&bb).collect();
        let mut rkmer = rightmost_kmer(&rc, k);
        while let Some(b) = unique_solid_out(&rkmer, table, threshold) {
            let mut next = rkmer.clone();
            next.push_right(b);
            let canon = next.canonical();
            if unique_solid_in(&next, table, threshold) != 1 {
                break;
            }
            if visited.contains(&canon) {
                circular = true;
                break;
            }
            rc.push(number_to_base(b));
            visited.insert(canon);
            rkmer = next;
        }
        bb = rev_comp(&rc).collect();
        // Canonical orientation, like the contig mode.
        if !canonical(&bb) {
            bb = rev_comp(&bb).collect();
        }
        let (coverage, min_cov, max_cov) = calc_coverage(&bb, table, k);
        unitigs.push(Unitig {
            bases: bb,
            id: 0,
            coverage,
            min_cov,
            max_cov,
            circular,
        });
    }
    unitigs
}

/// The rightmost k-mer of `bb` (5'->3' order, pushed right).
fn rightmost_kmer(bb: &[u8], k: usize) -> Kmer {
    let mut kmer = Kmer::new(k);
    for &b in &bb[bb.len() - k..] {
        kmer.push_right(base_code(b));
    }
    kmer
}

/// The single solid successor of `kmer`, or None at a branch or dead end.
fn unique_solid_out(kmer: &Kmer, table: &TadpoleTable, threshold: u32) -> Option<u8> {
    let counts = table.fill_right_counts(kmer);
    let mut out = None;
    for b in 0..4u8 {
        if counts[b as usize] >= threshold {
            if out.is_some() {
                return None;
            }
            out = Some(b);
        }
    }
    out
}

/// Number of solid predecessors of `kmer`.
fn unique_solid_in(kmer: &Kmer, table: &TadpoleTable, threshold: u32) -> usize {
    let counts = table.fill_left_counts(kmer);
    (0..4).filter(|&b| counts[b] >= threshold).count()
}

/// Descending length / coverage / sequence / id order (shared shape with
/// `contig_cmp`, without the branch-code fields).
fn unitig_cmp(a: &Unitig, b: &Unitig) -> std::cmp::Ordering {
    match a.bases.len().cmp(&b.bases.len()).reverse() {
        std::cmp::Ordering::Equal => {}
        x => return x,
    }
    if a.coverage != b.coverage {
        return a.coverage.partial_cmp(&b.coverage).unwrap().reverse();
    }
    match a.bases.cmp(&b.bases).reverse() {
        std::cmp::Ordering::Equal => {}
        x => return x,
    }
    a.id.cmp(&b.id).reverse()
}

/// Writes one unitig in FASTA with the contig-mode header fields (no
/// left/right branch codes: unitigs have none).
fn write_unitig<W: Write>(w: &mut W, u: &Unitig, links: Option<&[Link]>) -> Result<()> {
    let (gc, hh, caga) = calc_scalars(&u.bases);
    write!(
        w,
        ">unitig_{},len={},cov={},gc={},min={},max={},hh={},caga={}",
        u.id,
        u.bases.len(),
        fmt_fixed(u.coverage as f64, 1),
        fmt_fixed(gc as f64, 3),
        u.min_cov,
        u.max_cov,
        fmt_fixed(hh as f64, 3),
        fmt_fixed(caga as f64, 3),
    )?;
    if let Some(links) = links {
        for l in links {
            write!(
                w,
                " L:{}:{}:{}",
                if l.from_rc { '-' } else { '+' },
                l.to,
                if l.to_rc { '-' } else { '+' },
            )?;
        }
    }
    if u.circular {
        write!(w, ",circular")?;
    }
    writeln!(w)?;
    for chunk in u.bases.chunks(70) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// Seeding threshold for pass `i` (Java `minCountSeedCurrent` formula).
fn pass_threshold(opts: &AssembleOptions, i: usize) -> usize {
    let t = (opts.min_count_seed as f64 * opts.contig_pass_mult.powi(i as i32) * 0.92 - 0.25)
        .floor() as i64;
    (opts.min_count_seed as i64 + i as i64)
        .max(t)
        .min(i32::MAX as i64) as usize
}

/// One seeding scan over all table k-mers (BuildThread.processNextTable).
#[allow(clippy::too_many_arguments)]
fn scan_table(
    table: &TadpoleTable,
    threshold: usize,
    opts: &AssembleOptions,
    claimed: &mut HashSet<Kmer>,
    contigs: &mut Vec<Contig>,
    id_counter: &mut usize,
) {
    // Deterministic scan order by canonical k-mer sequence (the BBTools
    // hash-table cell order is memory-dependent and not portable). The
    // sorted snapshot is cached in the table, so all 16 seeding passes
    // iterate it linearly instead of re-sorting the HashMap each pass.
    let entries = table.sorted_entries();
    for (kmer, count) in entries.iter() {
        if *count < threshold as u32 {
            continue;
        }
        if claimed.contains(kmer) {
            continue;
        }
        claimed.insert(kmer.clone());
        if let Some(c) = make_contig(kmer, table, opts, claimed) {
            let mut c = c;
            c.id = *id_counter;
            *id_counter += 1;
            contigs.push(c);
        }
    }
}

/// Builds one contig from a claimed seed (Tadpole2.makeContig).
fn make_contig(
    seed: &Kmer,
    table: &TadpoleTable,
    opts: &AssembleOptions,
    claimed: &mut HashSet<Kmer>,
) -> Option<Contig> {
    let k = opts.k;
    // `base_at(0)` is the 3' end (last base pushed); rebuild 5'->3'.
    let mut bb: Vec<u8> = (0..k)
        .map(|i| number_to_base(seed.base_at(k - 1 - i)))
        .collect();
    debug_assert_eq!(bb.len(), k);

    let (right_status, mut right_ratio) = extend_to_right(&mut bb, table, opts, claimed);
    match right_status {
        DEAD_END | LOOP => {}
        BAD_SEED => return None,
        _ => {
            if bb.len() == k {
                // A branch or ownership failure at the seed rejects the contig.
                return None;
            }
            match right_status {
                BAD_OWNER => return None,
                F_BRANCH | D_BRANCH => {
                    right_ratio = calc_ratio(&right_counts_of(bb.as_slice(), table, opts))
                }
                B_BRANCH => right_ratio = calc_ratio(&left_counts_of(bb.as_slice(), table, opts)),
                _ => return None,
            }
        }
    }

    // Extend the left end by reverse-complementing and extending right.
    let mut rc: Vec<u8> = rev_comp(&bb).collect();
    let (left_status, mut left_ratio) = extend_to_right(&mut rc, table, opts, claimed);
    match left_status {
        DEAD_END | LOOP => {}
        BAD_SEED => return None,
        _ => match left_status {
            BAD_OWNER => return None,
            F_BRANCH | D_BRANCH => {
                left_ratio = calc_ratio(&right_counts_of(rc.as_slice(), table, opts))
            }
            B_BRANCH => left_ratio = calc_ratio(&left_counts_of(rc.as_slice(), table, opts)),
            _ => return None,
        },
    }
    bb = rev_comp(&rc).collect();

    // With bubble popping enabled (the default), BBTools keeps every contig
    // of at least k+minExtension internally; the minContigLen filter applies
    // only at output time (short contigs still anchor graph edges).
    if bb.len() >= k + opts.min_extension {
        let (coverage, min_cov, max_cov) = calc_coverage(&bb, table, k);
        if coverage < opts.min_coverage {
            return None;
        }
        // Canonical orientation (Contig.canonical + rcomp).
        let (bases, left_code, right_code, left_ratio, right_ratio) = if canonical(&bb) {
            (bb, left_status, right_status, left_ratio, right_ratio)
        } else {
            (
                rev_comp(&bb).collect(),
                right_status,
                left_status,
                right_ratio,
                left_ratio,
            )
        };
        Some(Contig {
            bases,
            id: 0,
            coverage,
            min_cov,
            max_cov,
            left_code,
            right_code,
            left_ratio,
            right_ratio,
            used: false,
            associate: false,
            flipped: false,
            left_edges: Vec::new(),
            right_edges: Vec::new(),
        })
    } else {
        None
    }
}

/// Counts of the four right/left extensions of a k-mer at `bb`'s 3'/5' end.
fn right_counts_of(bb: &[u8], table: &TadpoleTable, opts: &AssembleOptions) -> [u32; 4] {
    let k = opts.k;
    let mut kmer = Kmer::new(k);
    for &b in &bb[bb.len() - k..] {
        kmer.push_right(base_code(b));
    }
    table.fill_right_counts(&kmer)
}

fn left_counts_of(bb: &[u8], table: &TadpoleTable, opts: &AssembleOptions) -> [u32; 4] {
    let k = opts.k;
    let mut kmer = Kmer::new(k);
    for &b in &bb[bb.len() - k..] {
        kmer.push_right(base_code(b));
    }
    table.fill_left_counts(&kmer)
}

/// `extendToRight` (contig mode): bidirectional-aware greedy extension.
///
/// Returns the exit status and, for branch exits, the branch ratio.
fn extend_to_right(
    bb: &mut Vec<u8>,
    table: &TadpoleTable,
    opts: &AssembleOptions,
    claimed: &mut HashSet<Kmer>,
) -> (i32, f32) {
    let k = opts.k;
    if bb.len() < k {
        return (BAD_SEED, 0.0);
    }
    // Rightmost k-mer of the current sequence.
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    for &b in &bb[bb.len() - k..] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
    }
    if len < k {
        return (BAD_SEED, 0.0);
    }
    if table.get_count(&kmer) < opts.min_count_seed as u32 {
        return (BAD_SEED, 0.0);
    }

    let mut left = table.fill_left_counts(&kmer);
    let mut left_max_pos = argmax2(&left, &mut 0);
    let mut left_max = left[left_max_pos];
    let left_second_pos = second_highest_position(&left);
    let left_second = left[left_second_pos];

    let mut right = table.fill_right_counts(&kmer);
    let mut right_max_pos = argmax2(&right, &mut 0);
    let mut right_max = right[right_max_pos];
    let right_second_pos = second_highest_position(&right);
    let right_second = right[right_second_pos];

    if right_max < opts.min_count_extend as u32 {
        return (DEAD_END, 0.0);
    }
    if is_junction(right_max, right_second, opts) {
        let d = is_junction(left_max, left_second, opts);
        return if d {
            (D_BRANCH, calc_ratio(&right))
        } else {
            (F_BRANCH, calc_ratio(&right))
        };
    }
    if is_junction(left_max, left_second, opts) {
        return (B_BRANCH, calc_ratio(&left));
    }

    // The seed was claimed by the caller (single-threaded ownership).
    let max_len = 1_000_000_000usize;
    while bb.len() < max_len {
        let b = right_max_pos as u8;
        let evicted = kmer.base_at(k - 1);
        kmer.push_right(b);

        left = table.fill_left_counts(&kmer);
        left_max_pos = argmax2(&left, &mut 0);
        left_max = left[left_max_pos];
        let left_second_pos = second_highest_position(&left);
        let left_second = left[left_second_pos];

        right = table.fill_right_counts(&kmer);
        right_max_pos = argmax2(&right, &mut 0);
        right_max = right[right_max_pos];
        let right_second_pos = second_highest_position(&right);
        let right_second = right[right_second_pos];

        let fbranch = is_junction(right_max, right_second, opts);
        let bbranch = is_junction(left_max, left_second, opts);
        let hbranch = left_max_pos != evicted as usize && opts.branch_mult1 > 0.0;
        if bbranch || hbranch {
            let ratio = if fbranch {
                calc_ratio(&right)
            } else {
                calc_ratio(&left)
            };
            return if fbranch {
                (D_BRANCH, ratio)
            } else {
                (B_BRANCH, ratio)
            };
        }

        bb.push(number_to_base(b));

        // Loop detection / ownership claim (single-thread id=0).
        let canonical = kmer.canonical();
        if claimed.contains(&canonical) {
            return if fbranch {
                (F_BRANCH, calc_ratio(&right))
            } else {
                (LOOP, 0.0)
            };
        }
        claimed.insert(canonical);

        if fbranch {
            return (F_BRANCH, calc_ratio(&right));
        }
        if right_max < opts.min_count_extend as u32 {
            return (DEAD_END, 0.0);
        }
    }
    (BAD_OWNER, 0.0)
}

/// `KmerTableSet.calcCoverage`: mean/min/max canonical k-mer counts.
fn calc_coverage(bases: &[u8], table: &TadpoleTable, k: usize) -> (f32, usize, usize) {
    if bases.len() < k {
        return (0.0, 0, 0);
    }
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let mut sum = 0u64;
    let mut max = 0usize;
    let mut min = usize::MAX;
    let mut kmers = 0usize;
    for &b in bases {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len >= k {
            let count = table.get_count(&kmer) as usize;
            sum += count as u64;
            max = max.max(count);
            min = min.min(count);
            kmers += 1;
        }
    }
    if sum == 0 {
        (0.0, 0, 0)
    } else {
        (sum as f32 / kmers as f32, min, max)
    }
}

/// `Contig.calcScalarsFast`: gc fraction plus dimer-based hh/caga.
fn calc_scalars(bases: &[u8]) -> (f32, f32, f32) {
    if bases.len() < 2 {
        return (0.0, 0.0, 0.0);
    }
    let mut counts = [0u64; 16];
    let mut prev_bad = 8u8; // "N" so the first dimer is skipped
    let mut prev_val = 0u8;
    let mut at_sum = 0u64;
    let mut gc_sum = 0u64;
    for &b in bases {
        let gcbit = b >> 1;
        at_sum += (!gcbit & 1) as u64;
        gc_sum += (gcbit & !(b >> 3) & 1) as u64;
        let mut val = (b & 6) >> 1;
        val ^= (val & 2) >> 1;
        let bad = b & 8;
        if (prev_bad | bad) == 0 {
            counts[((prev_val << 2) | val) as usize] += 1;
        }
        prev_val = val;
        prev_bad = bad;
    }
    let aa = counts[0b0000];
    let tt = counts[0b1111];
    let at = counts[0b0011];
    let ta = counts[0b1100];
    let cc = counts[0b0101];
    let gg = counts[0b1010];
    let cg = counts[0b0110];
    let gc = counts[0b1001];
    let ac = counts[0b0001];
    let tg = counts[0b1110];
    let ag = counts[0b0010];
    let ct = counts[0b0111];
    let tc = counts[0b1101];
    let ga = counts[0b1000];
    let gt = counts[0b1011];
    let ca = counts[0b0100];
    let hh = (aa + cc + gg + tt) as f32 / (aa + tt + at + ta + cc + gg + cg + gc).max(1) as f32;
    let caga = 0.5
        * (1.0
            + (ca as i64 + tg as i64 - ga as i64 - tc as i64) as f32
                / (ac + ag + ca + ga + tc + tg + ct + gt).max(1) as f32);
    let gc_frac = gc_sum as f32 / (at_sum + gc_sum).max(1) as f32;
    (gc_frac, hh, caga)
}

/// A contig is canonical iff its sequence <= its reverse complement.
fn canonical(bases: &[u8]) -> bool {
    let n = bases.len();
    for i in 0..n {
        let a = bases[i];
        let b = complement(bases[n - 1 - i]);
        if a < b {
            return true;
        }
        if b < a {
            return false;
        }
    }
    true
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b,
    }
}

/// `calcRatio`: highest / second-highest count, 99 when no second branch.
fn calc_ratio(counts: &[u32; 4]) -> f32 {
    let mut a = 0u32;
    let mut b = 0u32;
    for &x in counts {
        if x > a {
            b = a;
            a = x;
        } else if x > b {
            b = x;
        }
    }
    if b < 1 {
        99.0
    } else {
        a as f32 / b as f32
    }
}

/// `ContigLengthComparator` (descending): length, coverage, sequence, id.
fn contig_cmp(a: &Contig, b: &Contig) -> std::cmp::Ordering {
    match a.bases.len().cmp(&b.bases.len()).reverse() {
        std::cmp::Ordering::Equal => {}
        x => return x,
    }
    if a.coverage != b.coverage {
        return a.coverage.partial_cmp(&b.coverage).unwrap().reverse();
    }
    match a.bases.cmp(&b.bases).reverse() {
        std::cmp::Ordering::Equal => {}
        x => return x,
    }
    a.id.cmp(&b.id).reverse()
}

/// Writes one contig in FASTA (SHORT_NAMES header, 70-column wrap).
fn write_contig<W: Write>(w: &mut W, c: &Contig) -> Result<()> {
    let (gc, hh, caga) = calc_scalars(&c.bases);
    writeln!(
        w,
        ">contig_{},len={},cov={},gc={},min={},max={},hh={},caga={}",
        c.id,
        c.bases.len(),
        fmt_fixed(c.coverage as f64, 1),
        fmt_fixed(gc as f64, 3),
        c.min_cov,
        c.max_cov,
        fmt_fixed(hh as f64, 3),
        fmt_fixed(caga as f64, 3),
    )?;
    for chunk in c.bases.chunks(70) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// `ByteBuilder.append(double, decimals)`: half-up fixed-point formatting.
fn fmt_fixed(x: f64, decimals: usize) -> String {
    if x == x.trunc() {
        return format!("{}", x as i64);
    }
    if decimals < 1 {
        return format!("{}", (x + 0.5) as i64);
    }
    let neg = x < 0.0;
    let x = x.abs();
    let inv = 10f64.powi(-(decimals as i32));
    let x = x + 0.5 * inv;
    let upper = x as i64;
    let lower = ((x - upper as f64) * 10f64.powi(decimals as i32)) as i64;
    format!(
        "{}{}.{:0width$}",
        if neg { "-" } else { "" },
        upper,
        lower,
        width = decimals
    )
}

/// Applies the BBTools phred round-trip to a record's quality scores.
fn canonicalize_quality(rec: &mut SeqRecord) {
    if rec.quality_scores().is_empty() {
        return;
    }
    let seq = rec.sequence().to_vec();
    let raw = rec.quality_scores().to_vec();
    let phred = to_phred(&seq, &raw);
    rec.set_quality(from_phred(&phred));
}

/*--------------------------------------------------------------------*/
/*  Contig graph and bubble popping (Tadpole.processContigs)          */
/*--------------------------------------------------------------------*/

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

type EdgeRef = Rc<RefCell<Edge>>;

impl Contig {
    fn left_kmer(&self, k: usize) -> Kmer {
        let mut kmer = Kmer::new(k);
        for &b in &self.bases[..k] {
            kmer.push_right(base_code(b));
        }
        kmer
    }

    fn right_kmer(&self, k: usize) -> Kmer {
        let mut kmer = Kmer::new(k);
        let n = self.bases.len();
        for &b in &self.bases[n - k..] {
            kmer.push_right(base_code(b));
        }
        kmer
    }

    fn left_forward_branch(&self) -> bool {
        self.left_code == F_BRANCH
    }

    fn right_forward_branch(&self) -> bool {
        self.right_code == F_BRANCH
    }

    fn add_left_edge(&mut self, e: EdgeRef) {
        let (dest, orient, depth, len) = {
            let eb = e.borrow();
            (eb.destination, eb.orientation, eb.depth, eb.length)
        };
        if let Some(old) = self.get_left_edge(dest, Some(orient)) {
            let mut ob = old.borrow_mut();
            if depth >= ob.depth && (ob.depth == 1 || ob.length == len) {
                ob.bases = e.borrow().bases.clone();
                ob.length = len;
                ob.depth += depth;
                return;
            }
        }
        self.left_edges.push(e);
    }

    fn add_right_edge(&mut self, e: EdgeRef) {
        let (dest, orient, depth, len) = {
            let eb = e.borrow();
            (eb.destination, eb.orientation, eb.depth, eb.length)
        };
        if let Some(old) = self.get_right_edge(dest, Some(orient)) {
            let mut ob = old.borrow_mut();
            if depth >= ob.depth && (ob.depth == 1 || ob.length == len) {
                ob.bases = e.borrow().bases.clone();
                ob.length = len;
                ob.depth += depth;
                return;
            }
        }
        self.right_edges.push(e);
    }

    fn get_left_edge(&self, dest: usize, orientation: Option<u8>) -> Option<EdgeRef> {
        self.left_edges
            .iter()
            .find(|e| {
                let e = e.borrow();
                e.destination == dest
                    && (orientation.is_none() || orientation == Some(e.orientation))
            })
            .cloned()
    }

    fn get_right_edge(&self, dest: usize, orientation: Option<u8>) -> Option<EdgeRef> {
        self.right_edges
            .iter()
            .find(|e| {
                let e = e.borrow();
                e.destination == dest
                    && (orientation.is_none() || orientation == Some(e.orientation))
            })
            .cloned()
    }

    fn remove_edges_to(&mut self, dest: usize) {
        self.left_edges.retain(|e| e.borrow().destination != dest);
        self.right_edges.retain(|e| e.borrow().destination != dest);
    }

    fn flip(&mut self, inbound: Option<&[EdgeRef]>) {
        self.flipped = !self.flipped;
        self.bases = rev_comp(&self.bases).collect();
        std::mem::swap(&mut self.left_code, &mut self.right_code);
        std::mem::swap(&mut self.left_ratio, &mut self.right_ratio);
        std::mem::swap(&mut self.left_edges, &mut self.right_edges);
        for e in &self.left_edges {
            e.borrow_mut().flip_source();
        }
        for e in &self.right_edges {
            e.borrow_mut().flip_source();
        }
        if let Some(inbound) = inbound {
            for e in inbound {
                e.borrow_mut().flip_dest();
            }
        }
    }

    fn renumber(&mut self, new_id: usize, inbound: Option<&[EdgeRef]>) {
        if self.id == new_id {
            return;
        }
        for e in &self.left_edges {
            e.borrow_mut().origin = new_id;
        }
        for e in &self.right_edges {
            e.borrow_mut().origin = new_id;
        }
        if let Some(inbound) = inbound {
            for e in inbound {
                e.borrow_mut().destination = new_id;
            }
        }
        self.id = new_id;
    }
}

/// Clears a contig's edges and detaches them from live sources
/// (Contig.removeAllEdges); `inbound` is the dest-map entry for `id`.
fn remove_all_edges(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].left_edges.clear();
    contigs[id].right_edges.clear();
    if let Some(inbound) = inbound {
        for e in inbound {
            let (dest, origin) = {
                let eb = e.borrow();
                (eb.destination, eb.origin)
            };
            if dest == id && origin != id {
                let source = &mut contigs[origin];
                if !source.used && !source.associate {
                    source.remove_edges_to(id);
                }
            }
        }
    }
}

fn set_used(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].used = true;
    remove_all_edges(id, inbound, contigs);
}

fn set_associate(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].associate = true;
    remove_all_edges(id, inbound, contigs);
}

/// Builds the contig end-kmer ownership map and edges
/// (Tadpole.initializeContigs + ProcessContigThread).
fn process_contigs(contigs: &mut [Contig], table: &TadpoleTable, opts: &AssembleOptions) {
    let k = opts.k;
    let mut end_claims: HashMap<Kmer, usize> = HashMap::new();
    for (i, c) in contigs.iter().enumerate() {
        end_claims.entry(c.left_kmer(k).canonical()).or_insert(i);
        end_claims.entry(c.right_kmer(k).canonical()).or_insert(i);
    }
    for i in 0..contigs.len() {
        process_contig_left(i, contigs, table, opts, &end_claims);
        process_contig_right(i, contigs, table, opts, &end_claims);
    }
}

fn process_contig_left(
    c_id: usize,
    contigs: &mut [Contig],
    table: &TadpoleTable,
    opts: &AssembleOptions,
    end_claims: &HashMap<Kmer, usize>,
) {
    if contigs[c_id].left_code == DEAD_END {
        return;
    }
    let k = opts.k;
    let kmer0 = contigs[c_id].left_kmer(k);
    let left = table.fill_left_counts(&kmer0);
    let left_max_pos = argmax2(&left, &mut 0);
    let left_max = left[left_max_pos];
    let mut edges_to_add: Vec<EdgeRef> = Vec::new();
    for x in 0..4u8 {
        let count = left[x as usize];
        if count > 0 && is_junction(left_max, count, opts) {
            let mut kmer = kmer0.clone();
            kmer.push_left(x);
            // Tadpole1 (k <= 31) walks the left edge in reverse-complement
            // space (`processContigLeft` swaps kmer/rkmer into `exploreRight`);
            // Tadpole2 (k > 31) walks it in forward space.
            if opts.k <= 31 {
                kmer = kmer.rc();
            }
            let mut bb = vec![number_to_base(x)];
            let (target, last_length, last_orientation) =
                explore_right(&kmer, table, opts, end_claims, contigs, &mut bb);
            if let Some(target) = target {
                edges_to_add.push(Rc::new(RefCell::new(Edge {
                    origin: c_id,
                    destination: target,
                    length: last_length,
                    orientation: last_orientation,
                    depth: count,
                    bases: bb,
                })));
            }
        }
    }
    for e in edges_to_add {
        contigs[c_id].add_left_edge(e);
    }
}

fn process_contig_right(
    c_id: usize,
    contigs: &mut [Contig],
    table: &TadpoleTable,
    opts: &AssembleOptions,
    end_claims: &HashMap<Kmer, usize>,
) {
    if contigs[c_id].right_code == DEAD_END {
        return;
    }
    let k = opts.k;
    let kmer0 = contigs[c_id].right_kmer(k);
    let right = table.fill_right_counts(&kmer0);
    let right_max_pos = argmax2(&right, &mut 0);
    let right_max = right[right_max_pos];
    let mut edges_to_add: Vec<EdgeRef> = Vec::new();
    for x in 0..4u8 {
        let count = right[x as usize];
        if count > 0 && is_junction(right_max, count, opts) {
            let mut kmer = kmer0.clone();
            kmer.push_right(x);
            let mut bb = vec![number_to_base(x)];
            let (target, last_length, mut last_orientation) =
                explore_right(&kmer, table, opts, end_claims, contigs, &mut bb);
            if let Some(target) = target {
                last_orientation |= 1;
                edges_to_add.push(Rc::new(RefCell::new(Edge {
                    origin: c_id,
                    destination: target,
                    length: last_length,
                    orientation: last_orientation,
                    depth: count,
                    bases: bb,
                })));
            }
        }
    }
    for e in edges_to_add {
        contigs[c_id].add_right_edge(e);
    }
}

/// `ProcessContigThread.exploreRight`: walks from an end k-mer to the next
/// contig end; returns (destination contig, path length, destination-side
/// orientation bit).
fn explore_right(
    kmer0: &Kmer,
    table: &TadpoleTable,
    opts: &AssembleOptions,
    end_claims: &HashMap<Kmer, usize>,
    contigs: &[Contig],
    bb: &mut Vec<u8>,
) -> (Option<usize>, usize, u8) {
    let k = opts.k;
    let mut kmer = kmer0.clone();
    let mut length = 1usize;
    let mut owner: Option<usize> = None;
    while length < 500 {
        owner = end_claims.get(&kmer.canonical()).copied();
        if owner.is_some() {
            break;
        }
        let left = table.fill_left_counts(&kmer);
        let left_max_pos = argmax2(&left, &mut 0);
        let left_max = left[left_max_pos];
        let left_second_pos = second_highest_position(&left);
        let left_second = left[left_second_pos];
        if is_junction(left_max, left_second, opts) {
            return (None, length, 0);
        }
        let right = table.fill_right_counts(&kmer);
        let right_max_pos = argmax2(&right, &mut 0);
        let right_max = right[right_max_pos];
        let right_second_pos = second_highest_position(&right);
        let right_second = right[right_second_pos];
        if right_max < opts.min_count_extend as u32 {
            return (None, length, 0);
        }
        if is_junction(right_max, right_second, opts) {
            return (None, length, 0);
        }
        bb.push(number_to_base(right_max_pos as u8));
        kmer.push_right(right_max_pos as u8);
        length += 1;
    }
    if let Some(owner) = owner {
        // Orientation: 0 if the destination's left k-mer matches, 2 if its
        // right k-mer matches (canonical comparison, like Java Kmer.equals).
        let dest = &contigs[owner];
        let mut temp = dest.left_kmer(k);
        let orientation = if kmer_eq(&temp, &kmer) {
            0
        } else {
            temp = dest.right_kmer(k);
            if kmer_eq(&temp, &kmer) {
                2
            } else {
                debug_assert!(false, "exploreRight destination mismatch");
                return (None, length, 0);
            }
        };
        (Some(owner), length, orientation)
    } else {
        (None, length, 0)
    }
}

fn kmer_eq(a: &Kmer, b: &Kmer) -> bool {
    a.canonical().cmp_bases(&b.canonical()) == std::cmp::Ordering::Equal
}

/// BubblePopper over the contig graph (assemble.BubblePopper).
struct BubblePopper {
    contigs: Vec<Contig>,
    dest_map: HashMap<usize, Vec<EdgeRef>>,
    k: usize,
    min_len: usize,
    center: usize,
    dest: usize,
    last_mutual_dest: i64,
    last_mutual_dest_orientation: i64,
    expansions: usize,
    contigs_absorbed: usize,
}

impl BubblePopper {
    fn dest_to_edge_map(&self) -> HashMap<usize, Vec<EdgeRef>> {
        let mut map: HashMap<usize, Vec<EdgeRef>> = HashMap::new();
        for c in &self.contigs {
            if c.used || c.associate {
                continue;
            }
            for e in &c.left_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
            for e in &c.right_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
        }
        map
    }

    fn expand(&mut self, center_id: usize) -> usize {
        self.center = center_id;
        let mut count = 0;
        while self.expand_right_simple() {
            count += 1;
        }
        while self.contigs[center_id].right_forward_branch() && self.expand_right() {
            count += 1;
            while self.expand_right_simple() {
                count += 1;
            }
        }
        let left_ok = {
            let c = &self.contigs[center_id];
            (c.left_code != LOOP && c.left_code != DEAD_END && !c.left_edges.is_empty())
                || c.left_forward_branch()
        };
        if left_ok {
            let inbound = self.dest_map.get(&center_id).cloned();
            self.contigs[center_id].flip(inbound.as_deref());
            while self.expand_right_simple() {
                count += 1;
            }
            while self.contigs[center_id].right_forward_branch() && self.expand_right() {
                count += 1;
                while self.expand_right_simple() {
                    count += 1;
                }
            }
        }
        count
    }

    fn expand_right_simple(&mut self) -> bool {
        let center_id = self.center;
        let outbound = self.contigs[center_id].right_edges.clone();
        if outbound.is_empty() || self.contigs[center_id].right_code == LOOP || outbound.len() > 1 {
            return false;
        }
        let left_edge = outbound[0].clone();
        let dest_id = left_edge.borrow().destination;
        let dest_right = left_edge.borrow().dest_right();
        if self.contigs[dest_id].used || dest_id == center_id {
            return false;
        }
        let (outbound_right, right_code) = {
            let d = &self.contigs[dest_id];
            if dest_right {
                (d.right_edges.clone(), d.right_code)
            } else {
                (d.left_edges.clone(), d.left_code)
            }
        };
        if right_code == LOOP {
            return false;
        }
        if !outbound_right.is_empty() {
            if outbound_right.len() > 1 {
                return false;
            }
            if outbound_right[0].borrow().destination != center_id {
                return false;
            }
        }
        if self.count_inbound(center_id, true) > 1 {
            return false;
        }
        if self.count_inbound(dest_id, dest_right) > 1 {
            return false;
        }
        if dest_right {
            let inbound = self.dest_map.get(&dest_id).cloned();
            self.contigs[dest_id].flip(inbound.as_deref());
        }
        self.merge(center_id, dest_id, left_edge)
    }

    fn count_inbound(&self, id: usize, dest_right: bool) -> usize {
        self.dest_map
            .get(&id)
            .map(|v| {
                v.iter()
                    .filter(|e| e.borrow().dest_right() == dest_right)
                    .count()
            })
            .unwrap_or(0)
    }

    fn merge(&mut self, left_id: usize, right_id: usize, left_edge: EdgeRef) -> bool {
        let k = self.k;
        let original_left_len = self.contigs[left_id].bases.len();
        let mut bb: Vec<u8> = self.contigs[left_id].bases.clone();
        {
            let eb = left_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        bb.extend_from_slice(&self.contigs[right_id].bases);
        self.contigs[left_id].bases = bb;
        self.contigs[left_id].right_edges.clear();
        let right_right = self.contigs[right_id].right_edges.clone();
        if right_right.is_empty() {
            self.contigs[left_id].right_edges = Vec::new();
        } else {
            for e in &right_right {
                e.borrow_mut().origin = left_id;
            }
            self.contigs[left_id].right_edges = right_right;
        }
        self.redirect_edges(right_id, left_id, true);
        let inbound_right = self.dest_map.get(&right_id).cloned();
        set_used(right_id, inbound_right.as_deref(), &mut self.contigs);
        let right_len = self.contigs[right_id].bases.len();
        let (right_max_cov, right_min_cov, right_code, right_ratio, right_coverage) = {
            let r = &self.contigs[right_id];
            (
                r.max_cov,
                r.min_cov,
                r.right_code,
                r.right_ratio,
                r.coverage,
            )
        };
        {
            let left = &mut self.contigs[left_id];
            left.max_cov = left.max_cov.max(right_max_cov);
            left.min_cov = left.min_cov.min(right_min_cov);
            left.right_code = right_code;
            left.right_ratio = right_ratio;
            let coverage_sum = left.coverage as f64 * (original_left_len - k + 1) as f64
                + right_coverage as f64 * (right_len - k + 1) as f64;
            left.coverage = (coverage_sum / (left.bases.len() - k + 1) as f64) as f32;
        }
        if self.is_loop(left_id) {
            self.contigs[left_id].left_code = LOOP;
            self.contigs[left_id].right_code = LOOP;
            let inbound = self.dest_map.get(&left_id).cloned();
            remove_all_edges(left_id, inbound.as_deref(), &mut self.contigs);
        }
        self.expansions += 1;
        self.contigs_absorbed += 1;
        true
    }

    fn redirect_edges(&mut self, from: usize, to: usize, dest_right: bool) {
        if from == to {
            return;
        }
        let Some(inbound_from) = self.dest_map.remove(&from) else {
            return;
        };
        let mut inbound_to = self.dest_map.get(&to).cloned().unwrap_or_default();
        for e in &inbound_from {
            if e.borrow().dest_right() == dest_right {
                e.borrow_mut().destination = to;
                inbound_to.push(e.clone());
            }
        }
        if inbound_to.is_empty() {
            self.dest_map.remove(&to);
        } else {
            self.dest_map.insert(to, inbound_to);
        }
    }

    fn is_loop(&self, id: usize) -> bool {
        let c = &self.contigs[id];
        if c.left_code == LOOP && c.right_code == LOOP {
            return true;
        }
        if c.left_edges.len() != 1 || c.right_edges.len() != 1 {
            return false;
        }
        for e in &c.left_edges {
            let e = e.borrow();
            if e.destination != id || !e.dest_right() {
                return false;
            }
        }
        for e in &c.right_edges {
            let e = e.borrow();
            if e.destination != id || e.dest_right() {
                return false;
            }
        }
        if let Some(inbound) = self.dest_map.get(&id) {
            for e in inbound {
                if e.borrow().origin != id {
                    return false;
                }
            }
        }
        true
    }

    fn expand_right(&mut self) -> bool {
        let center_id = self.center;
        self.dest = usize::MAX;
        self.last_mutual_dest = -1;
        self.last_mutual_dest_orientation = -1;
        if !self.contigs[center_id].right_forward_branch()
            || self.contigs[center_id].right_edges.is_empty()
        {
            return false;
        }
        let outbound = self.contigs[center_id].right_edges.clone();
        let Some(left_mid_edge) = self.find_representative_mid_edge(&outbound) else {
            return false;
        };
        let mid_id = left_mid_edge.borrow().destination;
        if self.contigs[mid_id].bases.len() < self.min_len {
            return false;
        }
        let mutual_dest = self.find_mutual_dest(&outbound);
        let mutual_dest_orientation = self.last_mutual_dest_orientation;
        let mutual_dest_right = (mutual_dest_orientation & 2) == 2;
        if mutual_dest < 0 || mutual_dest_orientation < 0 {
            return false;
        }
        let dest_id = mutual_dest as usize;
        if self.contigs[dest_id].used || dest_id == center_id {
            return false;
        }
        if mutual_dest_right && !self.contigs[dest_id].right_forward_branch() {
            return false;
        }
        if !mutual_dest_right && !self.contigs[dest_id].left_forward_branch() {
            return false;
        }
        let dest_outbound = {
            let d = &self.contigs[dest_id];
            if mutual_dest_right {
                d.right_edges.clone()
            } else {
                d.left_edges.clone()
            }
        };
        if dest_outbound.is_empty() {
            return false;
        }
        let mutual_dest2 = self.find_mutual_dest(&dest_outbound);
        if mutual_dest2 < 0 || mutual_dest2 as usize != center_id {
            return false;
        }
        let Some(mid_nodes) = self.fetch_mid_nodes(&outbound, true) else {
            return false;
        };
        if !self.mid_nodes_concur(&mid_nodes) {
            return false;
        }
        if mutual_dest_right {
            let inbound = self.dest_map.get(&dest_id).cloned();
            self.contigs[dest_id].flip(inbound.as_deref());
        }
        let right_mid_edge = self.contigs[mid_id].get_right_edge(dest_id, Some(1));
        let Some(right_mid_edge) = right_mid_edge else {
            return false;
        };
        self.dest = dest_id;
        self.pop(
            center_id,
            dest_id,
            mid_id,
            left_mid_edge,
            right_mid_edge,
            &mid_nodes,
        )
    }

    fn find_representative_mid_edge(&self, edges: &[EdgeRef]) -> Option<EdgeRef> {
        let mut mid_edge: Option<EdgeRef> = None;
        let mut mid_len = 0usize;
        for e in edges {
            let c = &self.contigs[e.borrow().destination];
            let clen = c.bases.len();
            match &mid_edge {
                None => {
                    mid_edge = Some(e.clone());
                    mid_len = clen;
                }
                Some(me) => {
                    let me_depth = me.borrow().depth;
                    let e_depth = e.borrow().depth;
                    if clen >= self.min_len
                        && (mid_len < self.min_len
                            || e_depth > me_depth
                            || (e_depth == me_depth && clen > mid_len))
                    {
                        mid_edge = Some(e.clone());
                        mid_len = clen;
                    }
                }
            }
        }
        mid_edge
    }

    fn find_mutual_dest(&mut self, edges: &[EdgeRef]) -> i64 {
        self.last_mutual_dest = -2;
        self.last_mutual_dest_orientation = -1;
        for e in edges {
            let mid_id = e.borrow().destination;
            if mid_id == self.center {
                return -1;
            }
            let outbound = {
                let mid = &self.contigs[mid_id];
                if e.borrow().dest_right() {
                    mid.left_edges.clone()
                } else {
                    mid.right_edges.clone()
                }
            };
            for o in &outbound {
                let ob = o.borrow();
                if self.last_mutual_dest < 0 {
                    self.last_mutual_dest = ob.destination as i64;
                    self.last_mutual_dest_orientation = (ob.orientation & 2) as i64;
                } else if self.last_mutual_dest != ob.destination as i64
                    || self.last_mutual_dest_orientation != (ob.orientation & 2) as i64
                {
                    return -1;
                }
            }
        }
        self.last_mutual_dest
    }

    fn fetch_mid_nodes(
        &mut self,
        outbound: &[EdgeRef],
        flip_as_needed: bool,
    ) -> Option<Vec<usize>> {
        let mut mid_nodes: Vec<usize> = Vec::new();
        for e in outbound {
            let mid_id = e.borrow().destination;
            if mid_nodes.contains(&mid_id) {
                return None;
            }
            if self.contigs[mid_id].used {
                return None;
            }
            mid_nodes.push(mid_id);
            if flip_as_needed && e.borrow().dest_right() {
                let inbound = self.dest_map.get(&mid_id).cloned();
                self.contigs[mid_id].flip(inbound.as_deref());
            }
        }
        Some(mid_nodes)
    }

    fn mid_nodes_concur(&self, mid_nodes: &[usize]) -> bool {
        let center_id = self.center;
        let dest_id = self.dest;
        let mut left_dest: i64 = -1;
        let mut right_dest: i64 = -1;
        for &mid_id in mid_nodes {
            let c = &self.contigs[mid_id];
            if c.left_edges.is_empty() || c.right_edges.is_empty() {
                return false;
            }
            for e in &c.left_edges {
                let eb = e.borrow();
                if left_dest < 0 {
                    left_dest = eb.destination as i64;
                } else if left_dest != eb.destination as i64 {
                    return false;
                }
                if eb.origin == eb.destination {
                    return false;
                }
            }
            for e in &c.right_edges {
                let eb = e.borrow();
                if right_dest < 0 {
                    right_dest = eb.destination as i64;
                } else if right_dest != eb.destination as i64 {
                    return false;
                }
                if eb.origin == eb.destination {
                    return false;
                }
            }
            let incoming = self.dest_map.get(&mid_id);
            let Some(incoming) = incoming else {
                return false;
            };
            for e in incoming {
                let origin = e.borrow().origin;
                if origin != center_id && origin != dest_id {
                    return false;
                }
            }
        }
        if left_dest >= 0 && left_dest as usize != center_id {
            return false;
        }
        if right_dest >= 0 && right_dest as usize != dest_id {
            return false;
        }
        left_dest >= 0 && right_dest >= 0
    }

    fn pop(
        &mut self,
        left_id: usize,
        right_id: usize,
        mid_id: usize,
        left_mid_edge: EdgeRef,
        right_mid_edge: EdgeRef,
        mid_nodes: &[usize],
    ) -> bool {
        let k = self.k;
        let original_left_len = self.contigs[left_id].bases.len();
        let mut bb: Vec<u8> = self.contigs[left_id].bases.clone();
        {
            let eb = left_mid_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        {
            let mid = &self.contigs[mid_id];
            let lim = mid.bases.len() - k + 1;
            if k - 1 < lim {
                bb.extend_from_slice(&mid.bases[k - 1..lim]);
            }
        }
        {
            let eb = right_mid_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        bb.extend_from_slice(&self.contigs[right_id].bases);
        self.contigs[left_id].bases = bb;
        self.contigs[left_id].right_edges.clear();
        let right_right = self.contigs[right_id].right_edges.clone();
        if right_right.is_empty() {
            self.contigs[left_id].right_edges = Vec::new();
        } else {
            for e in &right_right {
                e.borrow_mut().origin = left_id;
            }
            self.contigs[left_id].right_edges = right_right;
        }
        self.redirect_edges(right_id, left_id, true);
        let inbound_right = self.dest_map.get(&right_id).cloned();
        set_used(right_id, inbound_right.as_deref(), &mut self.contigs);
        for &c in mid_nodes {
            let inbound = self.dest_map.get(&c).cloned();
            if c == mid_id {
                set_used(c, inbound.as_deref(), &mut self.contigs);
            } else {
                set_associate(c, inbound.as_deref(), &mut self.contigs);
            }
        }
        let right_len = self.contigs[right_id].bases.len();
        let (right_max_cov, right_min_cov, right_code, right_ratio, right_coverage) = {
            let r = &self.contigs[right_id];
            (
                r.max_cov,
                r.min_cov,
                r.right_code,
                r.right_ratio,
                r.coverage,
            )
        };
        let (mid_max_cov, mid_min_cov) = {
            let m = &self.contigs[mid_id];
            (m.max_cov, m.min_cov)
        };
        {
            let left = &mut self.contigs[left_id];
            left.max_cov = left.max_cov.max(right_max_cov).max(mid_max_cov);
            left.min_cov = left.min_cov.min(right_min_cov).min(mid_min_cov);
            left.right_code = right_code;
            left.right_ratio = right_ratio;
            let coverage_sum = left.coverage as f64 * (original_left_len - k + 1) as f64
                + right_coverage as f64 * (right_len - k + 1) as f64;
            left.coverage = (coverage_sum / (left.bases.len() - k + 1) as f64) as f32;
        }
        if self.is_loop(left_id) {
            self.contigs[left_id].left_code = LOOP;
            self.contigs[left_id].right_code = LOOP;
            let inbound = self.dest_map.get(&left_id).cloned();
            remove_all_edges(left_id, inbound.as_deref(), &mut self.contigs);
        }
        self.expansions += 1;
        self.contigs_absorbed += 1 + mid_nodes.len();
        true
    }

    fn remove_dead_edges(&self, c: &mut Contig) {
        c.left_edges.retain(|e| {
            let d = e.borrow().destination;
            let dc = &self.contigs[d];
            !(dc.used || dc.associate)
        });
        c.right_edges.retain(|e| {
            let d = e.borrow().destination;
            let dc = &self.contigs[d];
            !(dc.used || dc.associate)
        });
    }
}

/// `Tadpole.popBubbles`: one bubble-popping pass, then deterministic sort and
/// renumbering.
fn pop_bubbles(contigs: &mut Vec<Contig>, opts: &AssembleOptions) {
    let dest_map = {
        let mut map: HashMap<usize, Vec<EdgeRef>> = HashMap::new();
        for c in contigs.iter() {
            if c.used || c.associate {
                continue;
            }
            for e in &c.left_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
            for e in &c.right_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
        }
        map
    };
    let mut bp = BubblePopper {
        contigs: std::mem::take(contigs),
        dest_map,
        k: opts.k,
        min_len: 2 * opts.k - 1,
        center: 0,
        dest: usize::MAX,
        last_mutual_dest: -1,
        last_mutual_dest_orientation: -1,
        expansions: 0,
        contigs_absorbed: 0,
    };
    for i in 0..bp.contigs.len() {
        let c = &bp.contigs[i];
        if !c.used && (c.left_forward_branch() || c.right_forward_branch()) {
            bp.expand(i);
        }
    }
    let dest_map2 = bp.dest_to_edge_map();
    let mut temp: Vec<Contig> = Vec::new();
    for i in 0..bp.contigs.len() {
        if bp.contigs[i].used {
            continue;
        }
        let mut c = bp.contigs[i].clone();
        bp.remove_dead_edges(&mut c);
        temp.push(c);
    }
    temp.sort_by(contig_cmp);
    for (new_id, c) in temp.iter_mut().enumerate() {
        let inbound = dest_map2.get(&c.id).cloned();
        c.renumber(new_id, inbound.as_deref());
    }
    *contigs = temp;
}

/// Deterministic longest-first sort and renumbering for the no-bubbles path
/// (bubble popping performs the same step while also renumbering edges).
fn finalize_contigs(contigs: &mut [Contig]) {
    contigs.sort_by(contig_cmp);
    for (new_id, c) in contigs.iter_mut().enumerate() {
        c.renumber(new_id, None);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unitig(bases: &[u8]) -> Unitig {
        Unitig {
            bases: bases.to_vec(),
            id: 0,
            coverage: 0.0,
            min_cov: 0,
            max_cov: 0,
            circular: false,
        }
    }

    fn assert_link(links: &[Link], to: usize, from_rc: bool, to_rc: bool) {
        assert!(
            links
                .iter()
                .any(|l| l.to == to && l.from_rc == from_rc && l.to_rc == to_rc),
            "missing link to {to} ({from_rc},{to_rc}): {links:?}"
        );
    }

    /// Branching (3'->5' on the forward strand) and reverse-strand
    /// directions of shared endpoint (k-1)-mers.
    #[test]
    fn links_directions_branch_and_rc() {
        // S is a 30 bp random fragment (k = 31 -> k-1 = 30).
        let s: Vec<u8> = b"GCTAAAGACAATTACATAACATACACGTCAG"[..30].to_vec();
        assert_eq!(s.len(), 30);
        let poly_a: Vec<u8> = b"A".repeat(50);
        // Random filler fragments (poly-C/poly-G would share a canonical
        // (k-1)-mer: each is the other's reverse complement).
        let x1: Vec<u8> = b"TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAA".to_vec();
        let x2: Vec<u8> = b"CCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAG".to_vec();
        let s_rc: Vec<u8> = rev_comp(&s).collect();

        // Branch: U0's right end and U1/U2's left ends all share S.
        let uts = vec![
            unitig(&[&poly_a[..], &s[..]].concat()),
            unitig(&[&s[..], &x1[..]].concat()),
            unitig(&[&s[..], &x2[..]].concat()),
        ];
        let links = compute_links(&uts, 31);
        assert_eq!(links[0].len(), 2);
        assert_link(&links[0], 1, false, false);
        assert_link(&links[0], 2, false, false);
        assert!(links[1].is_empty() && links[2].is_empty());

        // Reverse: U0's right end is rc(S), U1's left end is S.
        let uts = vec![
            unitig(&[&poly_a[..], &s_rc[..]].concat()),
            unitig(&[&s[..], &x1[..]].concat()),
        ];
        let links = compute_links(&uts, 31);
        assert_link(&links[0], 1, false, true);

        // 3'-3': both right ends are S -> reverse-strand representation.
        let uts = vec![
            unitig(&[&poly_a[..], &s[..]].concat()),
            unitig(&[&x2[..], &s[..]].concat()),
        ];
        let links = compute_links(&uts, 31);
        assert_link(&links[0], 1, true, true);
    }
}

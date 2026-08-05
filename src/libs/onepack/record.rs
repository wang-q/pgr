//! `.1aln` record layer: schema records, trace points, and the GDB skeleton.
//!
//! Rides on top of the [`super::container::Reader`] and interprets its lines
//! as alignment records (`A`/`R`/`D`/`T`/`X`/`U`) and the GDB skeleton
//! (`g`/`S`/`G`/`C`). This is the P4 building block of the `.1aln` migration;
//! it is independent of any sequence source.

use anyhow::{bail, Result};

use super::container::{List, Reader};

/// A single `.1aln` alignment record, decoded from the `A`/`R`/`D`/`T`/`X`/`U`
/// line sequence defined by `alncode.c` `Read_Aln_Overlap`/`Read_Aln_Trace`.
#[derive(Debug, Clone)]
pub struct AlnRecord {
    /// Contig index of the reference (`a`) side.
    pub aread: i64,
    /// Reference interval start (0-based, contig-relative).
    pub abpos: i64,
    /// Reference interval end.
    pub aepos: i64,
    /// Contig index of the query (`b`) side.
    pub bread: i64,
    /// Query interval start (0-based, contig-relative).
    pub bbpos: i64,
    /// Query interval end.
    pub bepos: i64,
    /// Whether the query (`b`) side is reverse-complemented (`R` line).
    pub comp: bool,
    /// Number of differences (substitutions + indels, `D` line).
    pub diffs: i64,
    /// Trace point deltas in `b` (`T` line). `tpoints[0]` is relative to
    /// `bbpos`; each later entry is relative to the previous point.
    pub tpoints: Vec<i64>,
    /// Per-interval difference counts (`X` line), one per trace interval.
    pub tdiffs: Vec<i64>,
    /// TR alignment unit length (`U` line), 0 if absent.
    pub period: i64,
}

impl AlnRecord {
    /// The interleaved trace array used by the expansion code:
    /// `trace[2i] = tdiffs[i]`, `trace[2i+1] = tpoints[i]`.
    pub fn interleaved_trace(&self) -> Vec<i64> {
        let mut t = Vec::with_capacity(2 * self.tpoints.len());
        for i in 0..self.tpoints.len() {
            t.push(self.tdiffs[i]);
            t.push(self.tpoints[i]);
        }
        t
    }

    /// Number of trace points (`k`). `tlen == 2 * points`.
    pub fn num_points(&self) -> usize {
        self.tpoints.len()
    }
}

/// A GDB skeleton contig (mirrors `GDB_CONTIG`).
#[derive(Debug, Clone)]
pub struct Contig {
    /// Scaffold index this contig belongs to.
    pub scaf: usize,
    /// Offset of the contig within its scaffold.
    pub sbeg: i64,
    /// Contig length.
    pub clen: i64,
}

/// A GDB skeleton scaffold (mirrors `GDB_SCAFFOLD`).
#[derive(Debug, Clone)]
pub struct Scaffold {
    /// Scaffold id (from the `S` line).
    pub name: String,
    /// First contig index (inclusive).
    pub fctg: usize,
    /// Last contig index (exclusive).
    pub ectg: usize,
    /// Total scaffold length (sum of contigs + gaps).
    pub slen: i64,
}

/// The `.1aln` GDB skeleton (`g` object): a set of scaffolds, each containing
/// a contiguous run of contigs separated by optional gaps.
#[derive(Debug, Clone, Default)]
pub struct Skeleton {
    pub scaffolds: Vec<Scaffold>,
    pub contigs: Vec<Contig>,
}

impl Skeleton {
    /// Total number of contigs.
    pub fn ncontig(&self) -> usize {
        self.contigs.len()
    }

    /// Total number of scaffolds.
    pub fn nscaff(&self) -> usize {
        self.scaffolds.len()
    }

    /// The scaffold name containing contig `c`.
    pub fn scaffold_name(&self, c: usize) -> Option<&str> {
        let contig = self.contigs.get(c)?;
        self.scaffolds.get(contig.scaf).map(|s| s.name.as_str())
    }
}

/// A record-level reader over a `.1aln` file.
///
/// Wraps a [`Reader`], advancing past the `t`/skeleton header lines to the
/// first alignment, then yielding [`AlnRecord`]s one at a time.
pub struct AlnFile {
    reader: Reader,
    /// Trace point spacing (from the `t` line).
    pub tspace: i64,
    /// The GDB skeletons, in file order (one per source genome). An alignment's
    /// `aread`/`bread` index into `skeletons[0]`/`skeletons[1]` respectively
    /// (the same skeleton when a single genome is aligned to itself).
    pub skeletons: Vec<Skeleton>,
    /// The `A` line that started the current (as-yet-unread) record, if any.
    pending: Option<super::container::Line>,
}

impl AlnFile {
    /// Open a `.1aln` file and position at the first alignment record.
    pub fn open(path: &str) -> Result<AlnFile> {
        let mut reader = Reader::open(path)?;
        let mut tspace = 0i64;
        let mut skeletons = Vec::new();
        // Advance through header records to the first `A`, capturing the `t`
        // line's tspace and all `g` skeletons on the way (FastGA `open_Aln_Read`).
        let pending = loop {
            let Some(line) = reader.read_line()? else {
                bail!("no alignment records found in {path}");
            };
            match line.line_type {
                b't' => tspace = line.int(0),
                b'g' => skeletons.push(read_skeleton(&mut reader)?),
                b'A' => break line,
                _ => {}
            }
        };
        if tspace == 0 {
            bail!("no tspace (`t` line) found before first alignment/skeleton in {path}");
        }
        Ok(AlnFile {
            reader,
            tspace,
            skeletons,
            pending: Some(pending),
        })
    }

    /// Access the underlying counts (for `stat`).
    pub fn counts(&self) -> &[super::container::Counts] {
        self.reader.counts()
    }

    /// Access the reference entries (`<` lines).
    pub fn references(&self) -> &[super::container::Reference] {
        self.reader.references.as_slice()
    }

    /// Access the provenance entries.
    pub fn provenance(&self) -> &[super::container::Provenance] {
        self.reader.provenance.as_slice()
    }

    /// Read the next alignment record, or `None` at end of data.
    pub fn next_record(&mut self) -> Result<Option<AlnRecord>> {
        let a = match self.pending.take() {
            Some(l) => l,
            None => match self.reader.read_line()? {
                Some(l) => l,
                None => return Ok(None),
            },
        };
        if a.line_type != b'A' {
            bail!("expected `A` line, got `{}`", a.line_type as char);
        }
        let aread = a.int(0);
        let abpos = a.int(1);
        let aepos = a.int(2);
        let bread = a.int(3);
        let bbpos = a.int(4);
        let bepos = a.int(5);

        let mut comp = false;
        let mut diffs = 0i64;
        let mut tpoints = Vec::new();
        let mut tdiffs = Vec::new();
        let mut period = 0i64;
        // Read the `R`/`D`/`T`/`X`/`U` lines that follow; stop at the next `A`.
        loop {
            let Some(line) = self.reader.read_line()? else {
                break;
            };
            match line.line_type {
                b'R' => comp = true,
                b'D' => diffs = line.int(0),
                b'T' => tpoints = int_list(&line)?,
                b'X' => tdiffs = int_list(&line)?,
                b'U' => period = line.int(0),
                b'A' => {
                    self.pending = Some(line);
                    break;
                }
                _ => {}
            }
        }
        if tpoints.len() != tdiffs.len() {
            bail!(
                "T/X line length mismatch for record {aread}:{}..{}: T={} X={}",
                abpos,
                aepos,
                tpoints.len(),
                tdiffs.len()
            );
        }
        Ok(Some(AlnRecord {
            aread,
            abpos,
            aepos,
            bread,
            bbpos,
            bepos,
            comp,
            diffs,
            tpoints,
            tdiffs,
            period,
        }))
    }

    /// Number of `A` records announced in the footer counts.
    pub fn n_overlaps(&self) -> i64 {
        self.reader.counts()[b'A' as usize].count
    }
}

/// Extract the integer list from a `T`/`X` line.
fn int_list(line: &super::container::Line) -> Result<Vec<i64>> {
    match &line.list {
        Some(List::Ints(v)) => Ok(v.clone()),
        _ => bail!("line `{}` has no INT list", line.line_type as char),
    }
}

/// Read a GDB skeleton (`g` object) from the reader, which is positioned just
/// after the `g` line. Stops at the first non-`S`/`G`/`C` line (the next `A`),
/// which is held aside for the record reader via `unread_line`.
fn read_skeleton(reader: &mut Reader) -> Result<Skeleton> {
    let mut scaffolds = Vec::new();
    let mut contigs = Vec::new();
    let mut cur: Option<Scaffold> = None;
    let mut spos = 0i64; // running offset within the current scaffold
    let mut ncontig = 0usize;
    loop {
        let Some(line) = reader.read_line()? else {
            // End of data: close the final scaffold.
            break;
        };
        match line.line_type {
            b'S' => {
                if let Some(s) = cur.take() {
                    scaffolds.push(s);
                }
                let name = match &line.list {
                    Some(List::Bytes(b)) => String::from_utf8_lossy(b).to_string(),
                    _ => bail!("`S` line has no string list"),
                };
                cur = Some(Scaffold {
                    name,
                    fctg: ncontig,
                    ectg: ncontig,
                    slen: 0,
                });
                spos = 0;
            }
            b'G' => {
                spos += line.int(0);
            }
            b'C' => {
                let clen = line.int(0);
                let scaf = scaffolds.len();
                contigs.push(Contig {
                    scaf,
                    sbeg: spos,
                    clen,
                });
                spos += clen;
                ncontig += 1;
                if let Some(s) = cur.as_mut() {
                    s.ectg = ncontig;
                    s.slen = spos;
                }
            }
            _ => {
                // First non-skeleton line: retain it for the record loop.
                reader.unread_line(line);
                break;
            }
        }
    }
    if let Some(s) = cur.take() {
        scaffolds.push(s);
    }
    // Fix `fctg` (the incremental assignment above used a placeholder).
    let mut fctg = 0usize;
    for s in &mut scaffolds {
        s.fctg = fctg;
        fctg = s.ectg;
    }
    Ok(Skeleton { scaffolds, contigs })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn golden_path() -> String {
        env!("CARGO_MANIFEST_DIR").to_string() + "/tests/genome/mg1655-sakai.1aln"
    }

    #[test]
    fn reads_golden_header_and_skeleton() {
        let aln = AlnFile::open(&golden_path()).unwrap();
        assert_eq!(aln.tspace, 100);
        assert_eq!(aln.n_overlaps(), 700);
        // The golden file carries two GDB skeletons (one per source genome).
        assert_eq!(aln.skeletons.len(), 2);
        let g1 = &aln.skeletons[0];
        let g2 = &aln.skeletons[1];
        // Single-contig genomes: each side is one scaffold with one contig.
        for (i, g) in aln.skeletons.iter().enumerate() {
            eprintln!(
                "skeleton {i}: {} scaffolds, {} contigs",
                g.nscaff(),
                g.ncontig()
            );
            for s in &g.scaffolds {
                eprintln!(
                    "  scaffold '{}' fctg={} ectg={} slen={}",
                    s.name, s.fctg, s.ectg, s.slen
                );
            }
        }
        // MG1655 (NC_000913) is single-contig; Sakai has three sequences
        // (NC_002695 chromosome + NC_002127/NC_002128 plasmids).
        assert_eq!(g1.nscaff(), 1);
        assert_eq!(g1.ncontig(), 1);
        assert_eq!(g2.nscaff(), 3);
        assert_eq!(g2.ncontig(), 3);
        // Scaffold ids are the source FASTA sequence names.
        assert!(!g1.scaffold_name(0).unwrap().is_empty());
        assert!(!g2.scaffold_name(0).unwrap().is_empty());
        // References point at the two source genomes (count 1 and 2) plus a
        // count-3 entry for the intermediate-file working directory (cpath).
        let refs = aln.references();
        assert_eq!(refs.len(), 3);
        assert_eq!(refs[0].count, 1);
        assert_eq!(refs[1].count, 2);
        assert_eq!(refs[2].count, 3);
    }

    #[test]
    fn reads_all_golden_records() {
        let mut aln = AlnFile::open(&golden_path()).unwrap();
        let mut n = 0;
        let mut total_diffs = 0i64;
        while let Some(rec) = aln.next_record().unwrap() {
            n += 1;
            total_diffs += rec.diffs;
            // T and X lists always have equal length (validated internally).
            assert_eq!(rec.tpoints.len(), rec.tdiffs.len());
            assert!(!rec.tpoints.is_empty());
            assert!(!rec.tdiffs.is_empty());
        }
        assert_eq!(n, 700);
        assert!(total_diffs > 0);
    }

    #[test]
    fn first_record_has_expected_shape() {
        let mut aln = AlnFile::open(&golden_path()).unwrap();
        let rec = aln.next_record().unwrap().unwrap();
        // Contig indices are valid against the skeletons.
        let g1 = &aln.skeletons[0];
        let g2 = &aln.skeletons[1];
        assert!(rec.aread >= 0 && (rec.aread as usize) < g1.ncontig());
        assert!(rec.bread >= 0 && (rec.bread as usize) < g2.ncontig());
        assert!(rec.aepos > rec.abpos);
        assert!(rec.bepos > rec.bbpos);
        assert!(rec.diffs >= 0);
        // The interleaved trace is even length and twice the point count.
        assert_eq!(rec.interleaved_trace().len(), 2 * rec.num_points());
    }
}

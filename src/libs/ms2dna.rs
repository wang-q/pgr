use anyhow::Result;
use std::io::Write;
use crate::libs::ms::SimpleRng;

pub fn build_anc_seq(gc: f64, nsite: usize, rng: &mut SimpleRng) -> Vec<u8> {
    let mut seq = vec![b'A'; nsite];
    for i in 0..nsite {
        let r1 = rng.next_f64();
        let r2 = rng.next_f64();
        seq[i] = if r1 <= gc {
            if r2 <= 0.5 { b'G' } else { b'C' }
        } else {
            if r2 <= 0.5 { b'A' } else { b'T' }
        };
    }
    seq
}

pub fn map_positions(positions: &[f64], nsite: usize, rng: &mut SimpleRng) -> Vec<usize> {
    let segsites = positions.len();
    let mut used = vec![false; nsite];
    let mut map = Vec::with_capacity(segsites);
    let mut used_count = 0;
    for i in 0..segsites {
        if used_count >= nsite {
            break;
        }
        let pos = positions[i];
        let mut p = ((pos * (nsite as f64)) as isize).clamp(0, (nsite - 1) as isize) as usize;
        while used[p] {
            p = (rng.next_f64() * (nsite as f64)) as usize;
        }
        used[p] = true;
        used_count += 1;
        map.push(p);
    }
    map
}

pub fn build_mut_seq(seq_anc: &[u8], map: &[usize], gc: f64, rng: &mut SimpleRng, nsite: usize) -> Vec<u8> {
    let mut seq_mut = vec![b'x'; nsite];
    for &p in map {
        let r1 = rng.next_f64();
        let r2 = rng.next_f64();
        let anc = seq_anc[p];
        let nuc = if r1 <= gc {
            match anc {
                b'G' => b'C',
                b'C' => b'G',
                _ => if r2 <= 0.5 { b'G' } else { b'C' },
            }
        } else {
            match anc {
                b'A' => b'T',
                b'T' => b'A',
                _ => if r2 <= 0.5 { b'A' } else { b'T' },
            }
        };
        seq_mut[p] = nuc;
    }
    seq_mut
}

pub fn write_fasta(
    writer: &mut dyn Write,
    nsam: usize,
    nsite: usize,
    map: &[usize],
    seq_anc: &[u8],
    seq_mut: &[u8],
    haplotypes: &[Vec<u8>],
    howmany: usize,
    npop: usize,
    sample_sizes: Option<&[usize]>,
    sample_counter: usize,
) -> Result<()> {
    let mut sc = 0usize;
    let mut pc = 0usize;
    let mut pos_to_seg = vec![usize::MAX; nsite];
    for (seg_idx, &p) in map.iter().enumerate() {
        pos_to_seg[p] = seg_idx;
    }
    for i in 0..nsam {
        writer.write_all(b">")?;
        if howmany > 1 {
            writer.write_fmt(format_args!("L{}", sample_counter))?;
        }
        if npop > 1 {
            if nsam > 1 {
                writer.write_fmt(format_args!("_P{}", pc + 1))?;
            } else {
                writer.write_fmt(format_args!("P{}", pc + 1))?;
            }
        }
        sc += 1;
        if npop > 1 || howmany > 1 {
            writer.write_fmt(format_args!("_S{}\n", sc))?;
        } else {
            writer.write_fmt(format_args!("S{}\n", sc))?;
        }
        if npop > 1 {
            if let Some(sizes) = sample_sizes {
                if sc == sizes[pc] {
                    sc = 0;
                    pc += 1;
                }
            }
        }
        let hap = &haplotypes[i];
        for j in 0..nsite {
            let seg_idx = pos_to_seg[j];
            if seg_idx != usize::MAX && seg_idx < hap.len() {
                let derived = hap[seg_idx] == b'1';
                let base = if derived { seq_mut[j] } else { seq_anc[j] };
                writer.write_all(&[base])?;
            } else {
                writer.write_all(&[seq_anc[j]])?;
            }
        }
        writer.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::ms::SimpleRng;

    #[test]
    fn test_map_positions_unique_bounds() {
        let mut rng = SimpleRng::new(123);
        let positions = vec![0.0, 0.1, 0.5, 0.9];
        let nsite = 10;
        let map = map_positions(&positions, nsite, &mut rng);
        assert_eq!(map.len(), positions.len());
        for &p in &map {
            assert!(p < nsite);
        }
        let mut seen = vec![false; nsite];
        for &p in &map {
            assert!(!seen[p]);
            seen[p] = true;
        }
    }

    #[test]
    fn test_write_fasta_headers_simple() {
        let mut out = Vec::new();
        let nsam = 2;
        let nsite = 5;
        let map = vec![1, 3];
        let mut rng = SimpleRng::new(42);
        let seq_anc = build_anc_seq(0.5, nsite, &mut rng);
        let seq_mut = build_mut_seq(&seq_anc, &map, 0.5, &mut rng, nsite);
        let haplotypes = vec![b"10".to_vec(), b"01".to_vec()];
        write_fasta(
            &mut out,
            nsam,
            nsite,
            &map,
            &seq_anc,
            &seq_mut,
            &haplotypes,
            1,
            1,
            None,
            1,
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.starts_with(">S1\n"));
        assert!(s.contains("\n>S2\n"));
        // no hard wrapping for nsite=5
        let lines: Vec<&str> = s.lines().collect();
        assert!(lines[1].len() == nsite);
        assert!(lines[3].len() == nsite);
    }

    #[test]
    fn test_write_fasta_multipop_headers() {
        let mut out = Vec::new();
        let nsam = 3;
        let nsite = 4;
        let map = vec![1];
        let mut rng = SimpleRng::new(11);
        let seq_anc = build_anc_seq(0.5, nsite, &mut rng);
        let seq_mut = build_mut_seq(&seq_anc, &map, 0.5, &mut rng, nsite);
        let haplotypes = vec![b"1".to_vec(), b"0".to_vec(), b"1".to_vec()];
        write_fasta(
            &mut out,
            nsam,
            nsite,
            &map,
            &seq_anc,
            &seq_mut,
            &haplotypes,
            1,
            2,
            Some(&[2, 1]),
            1,
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        let headers: Vec<&str> = s.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headers.len(), 3);
        assert!(headers[0].contains("_P1_S1"));
        assert!(headers[1].contains("_P1_S2"));
        assert!(headers[2].contains("_P2_S1"));
    }

    #[test]
    fn test_write_fasta_multisample_label() {
        let mut out = Vec::new();
        let nsam = 2;
        let nsite = 6;
        let map = vec![2, 4];
        let mut rng = SimpleRng::new(7);
        let seq_anc = build_anc_seq(0.5, nsite, &mut rng);
        let seq_mut = build_mut_seq(&seq_anc, &map, 0.5, &mut rng, nsite);
        let haplotypes = vec![b"10".to_vec(), b"01".to_vec()];
        write_fasta(
            &mut out,
            nsam,
            nsite,
            &map,
            &seq_anc,
            &seq_mut,
            &haplotypes,
            2,  // howmany > 1 triggers L{sample_counter}
            1,
            None,
            1,  // sample_counter = 1
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        let headers: Vec<&str> = s.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headers.len(), 2);
        assert!(headers[0].starts_with(">L1_S1"));
        assert!(headers[1].starts_with(">L1_S2"));
        // sequences single-line length equals nsite
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), nsite);
        assert_eq!(lines[3].len(), nsite);
    }

    #[test]
    fn test_write_fasta_zero_segsites_ancestral_only() {
        let mut out = Vec::new();
        let nsam = 3;
        let nsite = 8;
        let map: Vec<usize> = vec![]; // no segregating sites
        let mut rng = SimpleRng::new(21);
        let seq_anc = build_anc_seq(0.5, nsite, &mut rng);
        // seq_mut content shouldn't matter when map is empty
        let seq_mut = vec![b'x'; nsite];
        let haplotypes = vec![b"".to_vec(), b"".to_vec(), b"".to_vec()];
        write_fasta(
            &mut out,
            nsam,
            nsite,
            &map,
            &seq_anc,
            &seq_mut,
            &haplotypes,
            1,
            1,
            None,
            1,
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        // for each sample: header + sequence
        assert_eq!(lines.len(), nsam * 2);
        // sequence lines equal to ancestral sequence
        for i in 0..nsam {
            let seq_line = lines[i * 2 + 1];
            assert_eq!(seq_line.len(), nsite);
            assert_eq!(seq_line.as_bytes(), &seq_anc);
        }
    }

    #[test]
    fn test_write_fasta_label_combo_lps() {
        let mut out = Vec::new();
        let nsam = 3;
        let nsite = 5;
        let map = vec![0, 2];
        let mut rng = SimpleRng::new(33);
        let seq_anc = build_anc_seq(0.5, nsite, &mut rng);
        let seq_mut = build_mut_seq(&seq_anc, &map, 0.5, &mut rng, nsite);
        let haplotypes = vec![b"10".to_vec(), b"01".to_vec(), b"00".to_vec()];
        write_fasta(
            &mut out,
            nsam,
            nsite,
            &map,
            &seq_anc,
            &seq_mut,
            &haplotypes,
            2,          // howmany > 1 -> L prefix
            2,          // npop > 1 -> P and _S formatting
            Some(&[2, 1]),
            1,          // sample_counter
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        let headers: Vec<&str> = s.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headers.len(), nsam);
        assert!(headers[0].starts_with(">L1_P1_S1"));
        assert!(headers[1].starts_with(">L1_P1_S2"));
        assert!(headers[2].starts_with(">L1_P2_S1"));
    }
}

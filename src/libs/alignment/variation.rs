use crate::libs::nt::NT_VAL;
use anyhow::{anyhow, bail};
use intspan::IntSpan;
use itertools::Itertools;
use std::cmp::min;
use std::collections::BTreeMap;
use std::fmt;

use super::coords::indel_intspan;

#[derive(Default, Clone, Debug)]
pub struct Substitution {
    pub pos: i32,
    pub tbase: String,
    pub qbase: String,
    pub bases: String,
    pub mutant_to: String,
    pub freq: i32,
    pub pattern: String,
    pub obase: String,
}

/// To string
impl fmt::Display for Substitution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.pos,
            self.tbase,
            self.qbase,
            self.bases,
            self.mutant_to,
            self.freq,
            self.pattern,
            self.obase,
        )?;
        Ok(())
    }
}

/// Returns unpolarized substitutions
///
/// ```
/// let seqs = vec![
///     //        *
///     b"AAAATTTTGG".as_ref(),
///     b"aaaatttttg".as_ref(),
///     b"AAAATTTTAG".as_ref(),
/// ];
/// let subs = pgr::libs::alignment::get_subs(&seqs[..2]).unwrap();
/// let sub = subs.first().unwrap();
/// assert_eq!(sub.pos, 9);
/// assert_eq!(sub.tbase, "G".to_string());
/// assert_eq!(sub.qbase, "T".to_string());
/// assert_eq!(sub.bases, "GT".to_string());
/// assert_eq!(sub.mutant_to, "G<->T".to_string());
/// assert_eq!(sub.freq, 1);
/// assert_eq!(sub.pattern, "10".to_string());
/// assert_eq!(sub.obase, "".to_string());
///
/// let seqs = vec![
///     //*   **     * *
///     b"TTAG--GCTGAGAAGC".as_ref(),
///     b"GTAGCCGCTGA-AGGC".as_ref(),
/// ];
/// let subs = pgr::libs::alignment::get_subs(&seqs).unwrap();
/// let sub = subs.first().unwrap();
/// assert_eq!(sub.pos, 1);
/// assert_eq!(sub.tbase, "T".to_string());
/// assert_eq!(sub.qbase, "G".to_string());
/// assert_eq!(sub.bases, "TG".to_string());
/// assert_eq!(sub.mutant_to, "T<->G".to_string());
/// assert_eq!(sub.freq, 1);
/// assert_eq!(sub.pattern, "10".to_string());
/// assert_eq!(sub.obase, "".to_string());
///
/// let sub = subs.get(1).unwrap();
/// assert_eq!(sub.pos, 14);
/// assert_eq!(sub.tbase, "A".to_string());
/// assert_eq!(sub.qbase, "G".to_string());
/// assert_eq!(sub.bases, "AG".to_string());
/// assert_eq!(sub.mutant_to, "A<->G".to_string());
/// assert_eq!(sub.freq, 1);
/// assert_eq!(sub.pattern, "10".to_string());
/// assert_eq!(sub.obase, "".to_string());
///
/// ```
pub fn get_subs(seqs: &[&[u8]]) -> anyhow::Result<Vec<Substitution>> {
    let seq_count = seqs.len();
    let length = seqs[0].len();

    // For each position, search for polymorphic sites
    let mut bases_of: BTreeMap<usize, Vec<u8>> = BTreeMap::new();
    #[allow(clippy::needless_range_loop)]
    for pos in 0..length {
        let mut column = vec![];
        for seq in seqs.iter().take(seq_count) {
            column.push(seq[pos].to_ascii_uppercase());
        }

        if column.iter().all(|e| NT_VAL[*e as usize] <= 3) {
            // comparable += 1;
            if column.iter().any(|e| *e != column[0]) {
                // difference += 1;
                bases_of.insert(pos, column);
            }
        }
    }

    let mut sites = vec![];
    for pos in bases_of.keys().sorted() {
        let bases = bases_of
            .get(pos)
            .ok_or_else(|| anyhow!("position {} not found in bases_of", pos))?;

        let tbase = bases
            .first()
            .ok_or_else(|| anyhow!("empty bases at position {}", pos))?;

        let class: Vec<_> = bases.iter().unique().collect();

        if class.len() < 2 {
            bail!("No subs found in pos {}", pos);
        }

        let (freq, pattern, qbase) = if class.len() > 2 {
            (-1, "unknown".to_string(), "".to_string())
        } else {
            let mut freq: i32 = 0;
            let mut pattern = "".to_string();
            for base in bases {
                if tbase != base {
                    freq += 1;
                    pattern += "0";
                } else {
                    pattern += "1";
                }
            }
            let qbase = bases
                .iter()
                .find(|e| *e != tbase)
                .ok_or_else(|| anyhow!("no variant base found at position {}", pos))?;

            (freq, pattern, String::from_utf8(vec![*qbase])?)
        };

        let tbase = String::from_utf8(vec![*tbase])?;
        let mutant_to = if class.len() > 2 {
            "Complex".to_string()
        } else {
            format!("{}<->{}", tbase, qbase)
        };

        // mask previous variables
        let bases = String::from_utf8(bases.clone())?;
        let obase = "".to_string();
        let sub = Substitution {
            pos: (pos + 1) as i32,
            tbase,
            qbase,
            bases,
            mutant_to,
            freq: min(freq, seq_count as i32 - freq),
            pattern,
            obase,
        };
        sites.push(sub);
    }

    Ok(sites)
}

/// Polarize substitutions
///
/// ```
/// let seqs = vec![
///     //        *
///     b"AAAATTTTGG".as_ref(),
///     b"AAAATTTTAG".as_ref(),
///     b"AAAATTTTAG".as_ref(),
/// ];
/// let mut subs = pgr::libs::alignment::get_subs(&seqs[0..2]).unwrap();
/// pgr::libs::alignment::polarize_subs(&mut subs, &seqs[2]).unwrap();
/// let sub = subs.first().unwrap();
/// assert_eq!(sub.pos, 9);
/// assert_eq!(sub.tbase, "G".to_string());
/// assert_eq!(sub.qbase, "A".to_string());
/// assert_eq!(sub.bases, "GA".to_string());
/// assert_eq!(sub.mutant_to, "A->G".to_string());
/// assert_eq!(sub.freq, 1);
/// assert_eq!(sub.pattern, "10".to_string());
/// assert_eq!(sub.obase, "A".to_string());
///
/// let seqs = vec![
///     //*   **     * *
///     b"TTAG--GCTGAGAAGC".as_ref(),
///     b"GTAGCCGCTGA-AGGC".as_ref(),
///     b"TTAGCCGCTGAGAGGC".as_ref(),
/// ];
/// let mut subs = pgr::libs::alignment::get_subs(&seqs[0..2]).unwrap();
/// pgr::libs::alignment::polarize_subs(&mut subs, &seqs[2]).unwrap();
/// let sub = subs.first().unwrap();
/// assert_eq!(sub.pos, 1);
/// assert_eq!(sub.tbase, "T".to_string());
/// assert_eq!(sub.qbase, "G".to_string());
/// assert_eq!(sub.bases, "TG".to_string());
/// assert_eq!(sub.mutant_to, "T->G".to_string());
/// assert_eq!(sub.freq, 1);
/// assert_eq!(sub.pattern, "01".to_string());
/// assert_eq!(sub.obase, "T".to_string());
///
/// let sub = subs.get(1).unwrap();
/// assert_eq!(sub.pos, 14);
/// assert_eq!(sub.tbase, "A".to_string());
/// assert_eq!(sub.qbase, "G".to_string());
/// assert_eq!(sub.bases, "AG".to_string());
/// assert_eq!(sub.mutant_to, "G->A".to_string());
/// assert_eq!(sub.freq, 1);
/// assert_eq!(sub.pattern, "10".to_string());
/// assert_eq!(sub.obase, "G".to_string());
///
/// ```
pub fn polarize_subs(subs: &mut Vec<Substitution>, og: &[u8]) -> anyhow::Result<()> {
    for sub in subs {
        let pos = sub.pos;
        let obase_u8 = og[(pos - 1) as usize].to_ascii_uppercase();
        let obase = String::from_utf8(vec![obase_u8])?;

        if sub.qbase.is_empty() {
            // complex ingroup bases
            sub.obase = obase.clone();
        } else if NT_VAL[obase_u8 as usize] <= 3 {
            if sub.bases.contains(&obase) {
                // can polarize subs
                // ingroup bases have 2 classes
                let mut mutant_to = "".to_string();
                let mut freq = 0;
                let mut pattern = "".to_string();
                for base in sub.bases.as_bytes() {
                    if *base == obase_u8 {
                        pattern += "0";
                    } else {
                        pattern += "1";
                        freq += 1;
                        mutant_to = format!("{}->{}", obase, String::from_utf8(vec![*base])?);
                    }
                }
                sub.mutant_to = mutant_to;
                sub.freq = freq;
                sub.pattern = pattern;
                sub.obase = obase.clone();
            } else {
                // outgroup base is not equal to any nts
                sub.mutant_to = "Complex".to_string();
                sub.freq = -1;
                sub.pattern = "unknown".to_string();
                sub.obase = obase.clone();
            }
        } else {
            // outgroup base is N
            sub.mutant_to = "Complex".to_string();
            sub.freq = -1;
            sub.pattern = "unknown".to_string();
            sub.obase = obase.clone();
        }
    }

    Ok(())
}

#[derive(Default, Clone, Debug)]
pub struct Indel {
    pub start: i32,
    pub end: i32,
    pub length: i32,
    pub seq: String,
    pub all_seqs: String,
    pub freq: i32,
    pub occurred: String,
    pub itype: String,
    pub og_seq: String,
}

/// To string for Indel
impl fmt::Display for Indel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.start,
            self.end,
            self.length,
            self.seq,
            self.all_seqs,
            self.freq,
            self.occurred,
            self.itype,
            self.og_seq,
        )
    }
}

/// Returns unpolarized indels
///
/// 'D': means deletion relative to target/first seq
///      target is ----
/// 'I': means insertion relative to target/first seq
///      target is AAAA
///
/// ```
/// let seqs = vec![
///     //        ****
///     b"AAAATTTTGGGG".as_ref(),
///     b"AAAATTTT----".as_ref(),
///     b"AAAATTTTGGGG".as_ref(),
/// ];
/// let indels = pgr::libs::alignment::get_indels(&seqs).unwrap();
/// let indel = indels.first().unwrap();
/// assert_eq!(indel.start, 9);
/// assert_eq!(indel.end, 12);
/// assert_eq!(indel.length, 4);
/// assert_eq!(indel.seq, "GGGG");
/// assert_eq!(indel.all_seqs, "GGGG|----|GGGG");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "101");
/// assert_eq!(indel.itype, "I");
///
/// let seqs = vec![
///     //****
///     b"----TTTTGGGG".as_ref(),
///     b"AAAATTTTGGGG".as_ref(),
///     b"----TTTTGGGG".as_ref(),
/// ];
/// let indels = pgr::libs::alignment::get_indels(&seqs).unwrap();
/// let indel = indels.first().unwrap();
/// assert_eq!(indel.start, 1);
/// assert_eq!(indel.end, 4);
/// assert_eq!(indel.length, 4);
/// assert_eq!(indel.seq, "AAAA");
/// assert_eq!(indel.all_seqs, "----|AAAA|----");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "101");
/// assert_eq!(indel.itype, "D");
///
/// let seqs = vec![
///     //*   **     * *
///     b"TTAG--GCTGAGAAGC".as_ref(),
///     b"GTAGCCGCTGA-AGGC".as_ref(),
/// ];
/// let indels = pgr::libs::alignment::get_indels(&seqs).unwrap();
/// let indel = indels.first().unwrap();
/// assert_eq!(indel.start, 5);
/// assert_eq!(indel.end, 6);
/// assert_eq!(indel.length, 2);
/// assert_eq!(indel.seq, "CC");
/// assert_eq!(indel.all_seqs, "--|CC");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "10");
/// assert_eq!(indel.itype, "D");
///
/// let indel = indels.get(1).unwrap();
/// assert_eq!(indel.start, 12);
/// assert_eq!(indel.end, 12);
/// assert_eq!(indel.length, 1);
/// assert_eq!(indel.seq, "G");
/// assert_eq!(indel.all_seqs, "G|-");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "10");
/// assert_eq!(indel.itype, "I");
///
/// ```
// cargo test --doc alignment::get_indels
pub fn get_indels(seqs: &[&[u8]]) -> anyhow::Result<Vec<Indel>> {
    let seq_count = seqs.len();

    // Find all indel regions
    let mut indel_set = IntSpan::new();
    for seq in seqs {
        let seq_indel_set = indel_intspan(seq);
        indel_set.merge(&seq_indel_set);
    }

    let mut sites = vec![];
    for (start, end) in indel_set.spans() {
        let indel_length = end - start + 1;

        // Extract subsequences for each sequence
        let mut indel_seqs = vec![];
        for seq in seqs {
            let subseq = seq[(start - 1) as usize..end as usize].to_vec();
            indel_seqs.push(String::from_utf8(subseq)?);
        }
        let all_seqs = indel_seqs.join("|");

        // Determine the indel type
        let uniq_indel_seqs = indel_seqs.iter().unique().collect::<Vec<_>>();
        // seqs with least '-' char wins
        let indel_seq = uniq_indel_seqs
            .iter()
            .min_by_key(|s| s.chars().filter(|c| *c == '-').count())
            .ok_or_else(|| anyhow!("no indel sequence found at {}..{}", start, end))?
            .to_string();

        let itype = if uniq_indel_seqs.len() < 2 {
            bail!("No indel found at position {}..{}", start, end);
        } else if uniq_indel_seqs.len() > 2 || indel_seq.contains('-') {
            "C".to_string() // Complex indel
        } else if indel_seqs[0] == indel_seq {
            "I".to_string() // Insertion
        } else {
            "D".to_string() // Deletion
        };

        // Calculate frequency and occurrence pattern
        let (freq, occurred) = if itype == "C" {
            (-1, "unknown".to_string())
        } else {
            let mut freq = 0;
            let mut occurred = String::new();
            for seq in &indel_seqs {
                if seq == &indel_seqs[0] {
                    freq += 1;
                    occurred.push('1');
                } else {
                    occurred.push('0');
                }
            }
            (freq.min(seq_count as i32 - freq), occurred)
        };

        // Add to sites
        sites.push(Indel {
            start,
            end,
            length: indel_length,
            seq: indel_seq,
            all_seqs,
            freq,
            occurred,
            itype,
            og_seq: "".to_string(),
        });
    }

    Ok(sites)
}

/// Polarize indels based on outgroup sequence
///
/// ```
/// let seqs = vec![
///     //        ****
///     b"AAAATTTTGGGG".as_ref(),
///     b"AAAATTTT----".as_ref(),
///     b"AAAATTTTGGGG".as_ref(),
/// ];
/// let mut indels = pgr::libs::alignment::get_indels(&seqs[0..2]).unwrap();
/// pgr::libs::alignment::polarize_indels(&mut indels, &seqs[2]).unwrap();
/// let indel = indels.first().unwrap();
/// assert_eq!(indel.start, 9);
/// assert_eq!(indel.end, 12);
/// assert_eq!(indel.length, 4);
/// assert_eq!(indel.seq, "GGGG");
/// assert_eq!(indel.all_seqs, "GGGG|----");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "01");
/// assert_eq!(indel.itype, "D");
/// assert_eq!(indel.og_seq, "GGGG");
///
/// let seqs = vec![
///     //  ****
///     b"----TTTTGGGG".as_ref(),
///     b"AAAATTTTGGGG".as_ref(),
///     b"----TTTTGGGG".as_ref(),
/// ];
/// let mut indels = pgr::libs::alignment::get_indels(&seqs[0..2]).unwrap();
/// pgr::libs::alignment::polarize_indels(&mut indels, &seqs[2]).unwrap();
/// let indel = indels.first().unwrap();
/// assert_eq!(indel.start, 1);
/// assert_eq!(indel.end, 4);
/// assert_eq!(indel.length, 4);
/// assert_eq!(indel.seq, "AAAA");
/// assert_eq!(indel.all_seqs, "----|AAAA");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "01");
/// assert_eq!(indel.itype, "I");
/// assert_eq!(indel.og_seq, "----");
///
/// let seqs = vec![
///     //*   **     * *
///     b"TTAG--GCTGAGAAGC".as_ref(),
///     b"GTAGCCGCTGA-AGGC".as_ref(),
///     b"GTAGCCGCTGA--GGC".as_ref(),
/// ];
/// let mut indels = pgr::libs::alignment::get_indels(&seqs[0..2]).unwrap();
/// pgr::libs::alignment::polarize_indels(&mut indels, &seqs[2]).unwrap();
/// let indel = indels.first().unwrap();
/// assert_eq!(indel.start, 5);
/// assert_eq!(indel.end, 6);
/// assert_eq!(indel.length, 2);
/// assert_eq!(indel.seq, "CC");
/// assert_eq!(indel.all_seqs, "--|CC");
/// assert_eq!(indel.freq, 1);
/// assert_eq!(indel.occurred, "10");
/// assert_eq!(indel.itype, "D");
/// assert_eq!(indel.og_seq, "CC");
///
/// let indel = indels.get(1).unwrap();
/// assert_eq!(indel.start, 12);
/// assert_eq!(indel.end, 12);
/// assert_eq!(indel.length, 1);
/// assert_eq!(indel.seq, "G");
/// assert_eq!(indel.all_seqs, "G|-");
/// assert_eq!(indel.freq, -1);
/// assert_eq!(indel.occurred, "unknown");
/// assert_eq!(indel.itype, "C");
/// assert_eq!(indel.og_seq, "-");
///
/// ```
// cargo test --doc alignment::polarize_indels
pub fn polarize_indels(indels: &mut Vec<Indel>, og: &[u8]) -> anyhow::Result<()> {
    let og_indel_set = indel_intspan(og);

    for indel in indels {
        let og_seq = og[(indel.start - 1) as usize..indel.end as usize].to_vec();
        let og_seq = String::from_utf8(og_seq)?;
        indel.og_seq = og_seq.clone();

        let indel_seqs: Vec<String> = indel.all_seqs.split('|').map(|s| s.to_string()).collect();

        // Unique indel sequences including outgroup
        let mut uniq_indel_seqs = indel_seqs.clone();
        uniq_indel_seqs.push(og_seq.clone());
        uniq_indel_seqs.sort();
        uniq_indel_seqs.dedup();

        // Find the sequence with the least gaps
        let indel_seq = uniq_indel_seqs
            .iter()
            .min_by_key(|s| s.chars().filter(|c| *c == '-').count())
            .ok_or_else(|| anyhow!("no indel sequence found at {}..{}", indel.start, indel.end))?
            .clone();

        if uniq_indel_seqs.len() < 2 {
            anyhow::bail!("No indel found at position {}..{}", indel.start, indel.end);
        } else if uniq_indel_seqs.len() > 2 || indel_seq.contains('-') {
            indel.itype = "C".to_string(); // Complex indel
        } else {
            // Check outgroup sequence and indel set
            let indel_set = IntSpan::from_pair(indel.start, indel.end);

            if !og_seq.contains('-') && indel_seq != og_seq {
                // Outgroup has no gaps and is different from the reference
                //    AAA
                //    A-A
                // og ACA
                indel.itype = "C".to_string(); // Complex indel
            } else {
                let intersect = og_indel_set.intersect(&indel_set);
                if !intersect.is_empty() {
                    let island = og_indel_set.find_islands_ints(&indel_set);
                    if island.equals(&indel_set) {
                        // Outgroup has the same indel
                        //    NNNN
                        //    N--N
                        // og N--N
                        indel.itype = "I".to_string(); // Insertion
                    } else {
                        // Outgroup has a different indel
                        //    NNNN
                        //    N-NN
                        // og N--N
                        // or
                        //    NNNN
                        //    N--N
                        // og N-NN
                        indel.itype = "C".to_string(); // Complex indel
                    }
                } else {
                    // Outgroup has no gaps in this region
                    //    NNNN
                    //    N--N
                    // og NNNN
                    indel.itype = "D".to_string(); // Deletion
                }
            }
        }

        // Update frequency and occurrence pattern
        if indel.itype == "C" {
            indel.freq = -1;
            indel.occurred = "unknown".to_string();
        } else {
            let mut freq = 0;
            let mut occurred = String::new();
            for seq in &indel_seqs {
                if seq == &og_seq {
                    occurred.push('0');
                } else {
                    occurred.push('1');
                    freq += 1;
                }
            }
            indel.freq = freq;
            indel.occurred = occurred;
        }
    }

    Ok(())
}

/// Collect substitutions, polarizing with outgroup if provided.
///
/// When `outgroup` is `Some`, the last element of `seqs` is treated as the
/// outgroup and used to polarize substitutions from the remaining ingroup sequences.
pub fn collect_subs(seqs: &[&[u8]], outgroup: Option<&[u8]>) -> anyhow::Result<Vec<Substitution>> {
    let ingroup_count = if outgroup.is_some() {
        seqs.len() - 1
    } else {
        seqs.len()
    };
    let mut subs = get_subs(&seqs[..ingroup_count])?;
    if let Some(og) = outgroup {
        polarize_subs(&mut subs, og)?;
    }
    Ok(subs)
}

/// Build deduplicated VCF alt allele chars from a substitution, excluding the ref base.
pub fn vcf_alt_bases(sub: &Substitution) -> Vec<char> {
    let ref_base = sub.tbase.chars().next();
    let mut alt_bases: Vec<char> = vec![];
    for b in sub.bases.chars() {
        if matches!(b, 'A' | 'C' | 'G' | 'T') && Some(b) != ref_base {
            alt_bases.push(b);
        }
    }
    alt_bases.into_iter().unique().collect()
}

/// Collect indels, polarizing with outgroup if provided.
///
/// When `outgroup` is `Some`, the last element of `seqs` is treated as the
/// outgroup and used to polarize indels from the remaining ingroup sequences.
pub fn collect_indels(seqs: &[&[u8]], outgroup: Option<&[u8]>) -> anyhow::Result<Vec<Indel>> {
    let ingroup_count = if outgroup.is_some() {
        seqs.len() - 1
    } else {
        seqs.len()
    };
    let mut indels = get_indels(&seqs[..ingroup_count])?;
    if let Some(og) = outgroup {
        polarize_indels(&mut indels, og)?;
    }
    Ok(indels)
}

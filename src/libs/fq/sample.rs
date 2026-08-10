//! Deterministic read subsampling to a target base count.

use crate::libs::fmt::seq::{SeqReader, SeqRecord};
use anyhow::{bail, Result};
use std::io::Write;

/// `2^-53`, the multiplier of BBTools `FastRandomXoshiro.nextDouble()`.
const TWO_POW_MINUS_53: f64 = f64::from_bits(0x3CA0000000000000);

/// xoshiro256+ with SplitMix64 seeding, output-compatible with BBTools
/// `shared.FastRandomXoshiro` (Java `nextLong`/`nextDouble` semantics).
struct Xoshiro {
    s: [u64; 4],
}

impl Xoshiro {
    /// Seeds the generator; four warm-up draws match the Java constructor.
    fn new(seed: u64) -> Self {
        let mut s = [0u64; 4];
        s[0] = seed;
        for i in 1..4 {
            s[i] = Self::mix(s[i - 1]);
        }
        if s == [0; 4] {
            s = [0x5DEECE66D, 0xB, 0xCCA, 0xF00];
        }
        let mut rng = Self { s };
        for _ in 0..4 {
            rng.next_long();
        }
        rng
    }

    /// SplitMix64 finalizer.
    fn mix(x: u64) -> u64 {
        let x = x.wrapping_add(0x9E3779B97F4A7C15);
        let x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        let x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
        x ^ (x >> 31)
    }

    /// xoshiro256+ next value.
    fn next_long(&mut self) -> u64 {
        let result = self.s[0].wrapping_add(self.s[3]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// `(nextLong() >>> 11) * 0x1.0p-53`, Java-compatible.
    fn next_double(&mut self) -> f64 {
        ((self.next_long() >> 11) as f64) * TWO_POW_MINUS_53
    }
}

/// Subsets reads so the kept output contains about `target_bases` bases,
/// reproducing BBTools `reformat.sh samplebasestarget=<n> sampleseed=<seed>`
/// (exact mode, no upsampling). Decisions are made per interleaved pair (both
/// mates are kept together) in input order; the weight of a pair is the sum of
/// its read lengths.
pub fn sample<W: Write>(infile: &str, out: &mut W, target_bases: i64, seed: u64) -> Result<()> {
    if infile == "stdin" {
        bail!("sampling requires a file input (two passes over the data)");
    }
    let total = count_bases(infile)?;
    let mut rng = Xoshiro::new(seed);
    let mut reader = SeqReader::new(infile)?;
    let mut rec1 = SeqRecord::new();
    let mut rec2 = SeqRecord::new();
    let mut remaining = total;
    let mut target = target_bases;
    loop {
        if !reader.read_record(&mut rec1)? {
            break;
        }
        let has2 = reader.read_record(&mut rec2)?;
        let len1 = rec1.sequence().len() as i64;
        let bases = len1
            + if has2 {
                rec2.sequence().len() as i64
            } else {
                0
            };
        // All remaining pairs are empty (0 bases) once `remaining` hits 0;
        // stop rather than divide by zero on a trailing empty record.
        if remaining == 0 {
            break;
        }
        let prob = target as f64 / remaining as f64;
        if rng.next_double() < prob {
            target -= bases;
            write_record(out, &rec1)?;
            if has2 {
                write_record(out, &rec2)?;
            }
        }
        remaining -= bases;
        if !has2 {
            break;
        }
    }
    Ok(())
}

/// Total sequence bases of `infile`.
fn count_bases(infile: &str) -> Result<i64> {
    let mut reader = SeqReader::new(infile)?;
    let mut rec = SeqRecord::new();
    let mut total = 0i64;
    while reader.read_record(&mut rec)? {
        total += rec.sequence().len() as i64;
    }
    Ok(total)
}

/// Writes a FASTQ record, preserving the `name comment` header layout.
fn write_record<W: Write>(w: &mut W, rec: &SeqRecord) -> anyhow::Result<()> {
    let comment = rec.comment();
    let header = if comment.is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), comment)
    };
    crate::libs::fmt::fq::write_fq(w, &header, rec.sequence(), rec.quality_scores())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xoshiro_matches_java_next_double_sequence() {
        // Expected values from BBTools FastRandomXoshiro(1) via a scratch Java
        // run; guards the PRNG port against silent divergence.
        let mut rng = Xoshiro::new(1);
        let v: Vec<u64> = (0..5).map(|_| rng.next_long()).collect();
        assert_eq!(
            v,
            [
                0xee26d5d2c9ae4a10,
                0x4c7da8dc5fbdd4dc,
                0x3a621704ee180c4f,
                0x8609f81561c8855d,
                0xc5af9593638e9d19,
            ]
        );
    }
}

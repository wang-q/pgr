//! vc Huffman codec + DNA 2-bit codec (FastGA `vc*` / `Compress_DNA`).
//!
//! The `VcCodec` builds a length-limited (12-bit) Huffman code from a byte
//! histogram (Larmore & Hirschberg, JACM 1990), optionally reserving one
//! escape symbol for bytes absent from the training set. `encode`/`decode`
//! operate on a 64-bit bit buffer. The `DNAcodec` is a fixed 2-bits-per-base
//! codec (little-endian packing, A/C/G/T = 0/1/2/3).

use anyhow::{bail, Result};

const HUFF_CUTOFF: usize = 12;

#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    Empty,
    Filled,
    CodedWith,
    CodedRead,
}

/// A Huffman (or DNA) compressor codec.
pub struct VcCodec {
    state: State,
    isbig: bool,
    codebits: [u16; 256],
    codelens: [u8; 256],
    lookup: [u8; 0x10000],
    esc_code: i32,
    esc_len: i32,
    hist: [u64; 256],
}

impl Default for VcCodec {
    fn default() -> Self {
        Self::new()
    }
}

impl VcCodec {
    /// Create an empty compressor with a zeroed histogram.
    pub fn new() -> Self {
        VcCodec {
            state: State::Empty,
            isbig: false, // x86 is little-endian
            codebits: [0; 256],
            codelens: [0; 256],
            lookup: [0; 0x10000],
            esc_code: -1,
            esc_len: 0,
            hist: [0; 256],
        }
    }

    /// Add byte frequencies from `bytes` to the histogram.
    pub fn add_to_table(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.hist[b as usize] += 1;
        }
        if self.state == State::Empty {
            self.state = State::Filled;
        }
    }

    /// Build the length-limited Huffman tables from the histogram.
    pub fn create_codec(&mut self, partial: bool) -> Result<()> {
        if self.state == State::CodedWith || self.state == State::CodedRead {
            bail!("compressor already has a codec");
        }
        if self.state == State::Empty {
            bail!("compressor has no byte distribution data");
        }

        let hist = &self.hist;
        let mut code = [0usize; 256];
        let mut leng = [0usize; 256];
        let mut bits = [0u16; 256];
        let mut ncode = 0usize;
        let mut ecode: i64 = -1;

        // Collect present bytes; if partial, reserve one absent byte as escape.
        for (i, &f) in hist.iter().enumerate() {
            if f > 0 {
                code[ncode] = i;
                ncode += 1;
            } else if ecode < 0 {
                ecode = i as i64;
                code[ncode] = i;
                ncode += 1;
            }
        }
        let dcode = 2 * ncode;
        if ecode < 0 {
            // No absent byte: no escape code.
        }
        let partial = partial && ecode >= 0;

        // Sort codes by histogram frequency (ascending).
        code[..ncode].sort_by_key(|&i| hist[i]);

        // Package–merge (Larmore–Hirschberg) length-limited Huffman.
        let mut matrix = vec![vec![0u8; dcode]; HUFF_CUTOFF];
        let mut count1 = vec![0u64; dcode];
        let mut count2 = vec![0u64; dcode];
        let mut countb = vec![0u64; ncode];
        for n in 0..ncode {
            count1[n] = hist[code[n]];
            countb[n] = hist[code[n]];
            leng[n] = 0;
        }
        let mut lcnt = &mut count1;
        let mut ccnt = &mut count2;
        let mut llen = ncode - 1;
        for level in (1..HUFF_CUTOFF).rev() {
            let mut j = 0usize;
            let mut k = 0usize;
            let mut n = 0usize;
            while j < ncode || k < llen {
                if k >= llen || (j < ncode && countb[j] <= lcnt[k] + lcnt[k + 1]) {
                    ccnt[n] = countb[j];
                    matrix[level][n] = 1;
                    j += 1;
                } else {
                    ccnt[n] = lcnt[k] + lcnt[k + 1];
                    matrix[level][n] = 0;
                    k += 2;
                }
                n += 1;
            }
            llen = n - 1;
            std::mem::swap(&mut lcnt, &mut ccnt);
        }

        // Back trace to recover code lengths.
        let mut span = 2 * (ncode - 1);
        for row in &matrix[1..HUFF_CUTOFF] {
            let mut j = 0usize;
            for &m in &row[..span] {
                if m != 0 {
                    leng[j] += 1;
                    j += 1;
                }
            }
            span = 2 * (span - j);
        }
        for l in leng[..span].iter_mut() {
            *l += 1;
        }

        // Build canonical codes.
        let mut llen = leng[0];
        let mut lbits = (1u16 << llen) - 1;
        bits[0] = lbits;
        for n in 1..ncode {
            while (lbits & 0x1) == 0 {
                lbits >>= 1;
                llen -= 1;
            }
            lbits -= 1;
            while llen < leng[n] {
                lbits = (lbits << 1) | 0x1;
                llen += 1;
            }
            bits[n] = lbits;
        }

        // Assign code lengths and bit values.
        for i in 0..256 {
            self.codelens[i] = 0;
            self.codebits[i] = 0;
        }
        for i in 0..ncode {
            self.codelens[code[i]] = leng[i] as u8;
            self.codebits[code[i]] = bits[i];
        }

        // Fill decoder lookup table: every 16-bit prefix maps to a symbol.
        for i in 0..256 {
            if self.codelens[i] > 0 {
                let base = (self.codebits[i] as usize) << (16 - self.codelens[i] as usize);
                let powr = 1usize << (16 - self.codelens[i] as usize);
                for j in 0..powr {
                    self.lookup[base + j] = i as u8;
                }
            }
        }

        if partial {
            self.esc_code = ecode as i32;
            self.esc_len = self.codelens[ecode as usize] as i32;
            self.codelens[ecode as usize] = 0;
        } else {
            self.esc_code = -1;
        }
        self.state = State::CodedWith;
        Ok(())
    }

    /// Serialize the codec into a byte blob (for the footer `;` line).
    pub fn serialize(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.push(self.isbig as u8);
        out.extend_from_slice(&self.esc_code.to_le_bytes());
        out.extend_from_slice(&self.esc_len.to_le_bytes());
        for i in 0..256 {
            out.push(self.codelens[i]);
            if self.codelens[i] > 0 || i as i32 == self.esc_code {
                out.extend_from_slice(&self.codebits[i].to_le_bytes());
            }
        }
        out
    }

    /// Deserialize a codec from a serialized blob (footer `;` line).
    pub fn deserialize(data: &[u8]) -> Result<VcCodec> {
        let mut v = VcCodec {
            state: State::CodedRead,
            isbig: false,
            codebits: [0; 256],
            codelens: [0; 256],
            lookup: [0; 0x10000],
            esc_code: -1,
            esc_len: 0,
            hist: [0; 256],
        };
        let inbig = data[0] != 0;
        let mut pos = 1;
        let rd_i32 = |p: &mut usize| -> i32 {
            let b: [u8; 4] = data[*p..*p + 4].try_into().unwrap();
            *p += 4;
            i32::from_le_bytes(b)
        };
        let rd_u16 = |p: &mut usize| -> u16 {
            let b: [u8; 2] = data[*p..*p + 2].try_into().unwrap();
            *p += 2;
            u16::from_le_bytes(b)
        };
        v.esc_code = rd_i32(&mut pos);
        v.esc_len = rd_i32(&mut pos);
        for i in 0..256 {
            v.codelens[i] = data[pos];
            pos += 1;
            if v.codelens[i] > 0 || i as i32 == v.esc_code {
                let b16 = rd_u16(&mut pos);
                v.codebits[i] = if inbig != v.isbig {
                    b16.swap_bytes()
                } else {
                    b16
                };
            }
        }
        if v.esc_code >= 0 {
            v.codelens[v.esc_code as usize] = v.esc_len as u8;
        }
        for i in 0..256 {
            if v.codelens[i] > 0 {
                let base = (v.codebits[i] as usize) << (16 - v.codelens[i] as usize);
                let powr = 1usize << (16 - v.codelens[i] as usize);
                for j in 0..powr {
                    v.lookup[base + j] = i as u8;
                }
            }
        }
        if v.esc_code >= 0 {
            v.codelens[v.esc_code as usize] = 0;
        }
        Ok(v)
    }

    /// Encode `ibytes`, returning `(compressed_bytes, bit_count)`.
    pub fn encode(&self, ibytes: &[u8]) -> (Vec<u8>, i64) {
        let ilen = ibytes.len();
        let ibits = (ilen * 8) as i64;

        let mut out = Vec::<u8>::new();
        let mut ocode: u64 = 0;
        let mut rem: i64 = 62;
        let mut tbits: i64 = 2;
        let mut k = 0usize;

        while k < ilen {
            let x = ibytes[k] as usize;
            let n = self.codelens[x] as i64;
            if n == 0 {
                let esc = self.esc_code;
                if esc < 0 {
                    // Unsolvable byte: fall back to plaintext (should not happen
                    // when the codec covers its input alphabet).
                    let mut v = Vec::with_capacity(1 + ilen);
                    v.push(0xff);
                    v.extend_from_slice(ibytes);
                    return (v, ibits + 8);
                }
                let elen = self.esc_len as i64;
                tbits += 8 + elen;
                if tbits > ibits {
                    break;
                }
                // escape symbol
                let c = self.codebits[esc as usize] as u64;
                ocode_put(elen, c, &mut ocode, &mut out, &mut rem);
                // raw byte
                let c = x as u64;
                ocode_put(8, c, &mut ocode, &mut out, &mut rem);
                k += 1;
            } else {
                tbits += n;
                if tbits > ibits {
                    break;
                }
                let c = self.codebits[x] as u64;
                ocode_put(n, c, &mut ocode, &mut out, &mut rem);
                k += 1;
            }
        }

        if k < ilen {
            // Plaintext fallback.
            let mut v = Vec::with_capacity(1 + ilen);
            v.push(0xff);
            v.extend_from_slice(ibytes);
            return (v, ibits + 8);
        }

        // Flush the partial word.
        let bcode = ocode.to_le_bytes();
        if self.isbig {
            let n = ((71 - rem) >> 3) as usize;
            out.extend_from_slice(&bcode[..n]);
        } else {
            let n = 7 - ((63 - rem) >> 3);
            for kk in ((n as usize)..=7).rev() {
                out.push(bcode[kk]);
            }
        }

        if tbits >= 64 && !self.isbig {
            out.swap(0, 7);
        }
        (out, tbits)
    }

    /// Decode `ilen` bits of `ibytes` back into bytes.
    pub fn decode(&self, ibytes: &[u8], ilen: i64) -> Vec<u8> {
        if ibytes[0] == 0xff {
            let olen = ((ilen >> 3) - 1) as usize;
            return ibytes[1..1 + olen].to_vec();
        }

        let mut buf = ibytes.to_vec();
        let inbig = buf[0] & 0x40 != 0;
        if !inbig && ilen >= 64 {
            buf.swap(7, 0);
        }
        if inbig != self.isbig {
            for chunk in buf.chunks_mut(8) {
                chunk.reverse();
            }
        }

        let mut il = ilen;
        let mut p = 0usize;
        let mut icode: u64 = if il < 64 {
            let mut v = 0u64;
            for (pp, &b) in buf[..((il + 7) / 8) as usize].iter().enumerate() {
                v |= (b as u64) << (56 - (pp * 8));
            }
            v
        } else {
            let v = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            p = 8;
            v
        };
        icode <<= 2;
        il -= 2;
        let mut rem = 62;
        if rem > il {
            rem = il;
        }
        let mut ncode: u64 = 0;
        let mut nem: i64 = 0;
        let mut o = Vec::new();

        while il > 0 {
            let c = self.lookup[(icode >> 48) as usize];
            if c as i32 == self.esc_code {
                // GET(esc_len)
                {
                    let n = self.esc_len as i64;
                    il -= n;
                    icode <<= n;
                    rem -= n;
                    while rem < 16 {
                        let z = 64 - rem;
                        icode |= ncode >> rem;
                        if nem > z {
                            nem -= z;
                            ncode <<= z;
                            rem = 64;
                            break;
                        } else {
                            rem += nem;
                            if rem >= il {
                                break;
                            } else if il - rem < 64 {
                                nem = il - rem;
                                ncode = 0;
                                for kk in (0..nem).step_by(8) {
                                    ncode |= (buf[p] as u64) << (56 - kk);
                                    p += 1;
                                }
                            } else {
                                ncode = u64::from_le_bytes(buf[p..p + 8].try_into().unwrap());
                                p += 8;
                                nem = 64;
                            }
                        }
                    }
                }
                let c = (icode >> 56) as u8;
                // GET(8)
                {
                    il -= 8;
                    icode <<= 8;
                    rem -= 8;
                    while rem < 16 {
                        let z = 64 - rem;
                        icode |= ncode >> rem;
                        if nem > z {
                            nem -= z;
                            ncode <<= z;
                            rem = 64;
                            break;
                        } else {
                            rem += nem;
                            if rem >= il {
                                break;
                            } else if il - rem < 64 {
                                nem = il - rem;
                                ncode = 0;
                                for kk in (0..nem).step_by(8) {
                                    ncode |= (buf[p] as u64) << (56 - kk);
                                    p += 1;
                                }
                            } else {
                                ncode = u64::from_le_bytes(buf[p..p + 8].try_into().unwrap());
                                p += 8;
                                nem = 64;
                            }
                        }
                    }
                }
                o.push(c);
            } else {
                let n = self.codelens[c as usize] as i64;
                il -= n;
                icode <<= n;
                rem -= n;
                while rem < 16 {
                    let z = 64 - rem;
                    icode |= ncode >> rem;
                    if nem > z {
                        nem -= z;
                        ncode <<= z;
                        rem = 64;
                        break;
                    } else {
                        rem += nem;
                        if rem >= il {
                            break;
                        } else if il - rem < 64 {
                            nem = il - rem;
                            ncode = 0;
                            for kk in (0..nem).step_by(8) {
                                ncode |= (buf[p] as u64) << (56 - kk);
                                p += 1;
                            }
                        } else {
                            ncode = u64::from_le_bytes(buf[p..p + 8].try_into().unwrap());
                            p += 8;
                            nem = 64;
                        }
                    }
                }
                o.push(c);
            }
        }
        o
    }
}

/// The bit-packing step (`OCODE` macro) of `vcEncode`.
fn ocode_put(l: i64, c: u64, ocode: &mut u64, out: &mut Vec<u8>, rem: &mut i64) {
    *rem -= l;
    if *rem <= 0 {
        *ocode |= c >> (-(*rem));
        out.extend_from_slice(&ocode.to_le_bytes());
        if *rem < 0 {
            *rem += 64;
            *ocode = c << *rem;
        } else {
            *rem = 64;
            *ocode = 0;
        }
    } else {
        *ocode |= c << *rem;
    }
}

/// The `Number[128]` table mapping ASCII bases to 2-bit codes.
const NUMBER: [u8; 128] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// Compress DNA into 2 bits per base (little-endian packing), `t` must hold
/// `(len+3)>>2` bytes. Returns the number of bits (2*len).
pub fn compress_dna(s: &[u8], t: &mut Vec<u8>) -> usize {
    let len = s.len();
    let mut i = 0usize;
    let lim = len.saturating_sub(3);
    while i < lim {
        let b = NUMBER[s[i] as usize & 0x7f]
            | (NUMBER[s[i + 1] as usize & 0x7f] << 2)
            | (NUMBER[s[i + 2] as usize & 0x7f] << 4)
            | (NUMBER[s[i + 3] as usize & 0x7f] << 6);
        t.push(b);
        i += 4;
    }
    match i.checked_sub(lim) {
        Some(0) => {
            let b = NUMBER[s[i] as usize & 0x7f]
                | (NUMBER[s[i + 1] as usize & 0x7f] << 2)
                | (NUMBER[s[i + 2] as usize & 0x7f] << 4);
            t.push(b);
        }
        Some(1) => {
            let b = NUMBER[s[i] as usize & 0x7f] | (NUMBER[s[i + 1] as usize & 0x7f] << 2);
            t.push(b);
        }
        Some(2) => t.push(NUMBER[s[i] as usize & 0x7f]),
        _ => {}
    }
    len << 1
}

/// Uncompress 2-bit DNA back into lowercase `a/c/g/t`. `s` holds the N
/// compressed bytes for `nbases` bases; `t` receives `nbases` bases.
pub fn uncompress_dna(s: &[u8], nbases: usize, t: &mut Vec<u8>) {
    let base = [b'a', b'c', b'g', b't'];
    let mut i = 0usize; // base output position
    let mut j = 0usize; // compressed byte index
    let lim = nbases.saturating_sub(3);
    while i < lim {
        let byte = s[j];
        j += 1;
        t.extend_from_slice(&[
            base[(byte & 0x3) as usize],
            base[((byte >> 2) & 0x3) as usize],
            base[((byte >> 4) & 0x3) as usize],
            base[((byte >> 6) & 0x3) as usize],
        ]);
        i += 4;
    }
    match i.checked_sub(lim) {
        Some(0) => {
            let byte = s[j];
            t.extend_from_slice(&[
                base[(byte & 0x3) as usize],
                base[((byte >> 2) & 0x3) as usize],
                base[((byte >> 4) & 0x3) as usize],
            ]);
        }
        Some(1) => {
            let byte = s[j];
            t.extend_from_slice(&[
                base[(byte & 0x3) as usize],
                base[((byte >> 2) & 0x3) as usize],
            ]);
        }
        Some(2) => t.push(base[(s[j] & 0x3) as usize]),
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_round_trip() {
        let s = b"ACGTACGTnnnnACGT";
        let mut comp = Vec::new();
        compress_dna(s, &mut comp);
        let mut out = Vec::new();
        uncompress_dna(&comp, s.len(), &mut out);
        assert_eq!(out, b"acgtacgtaaaaacgt");
    }

    #[test]
    fn huffman_round_trip() {
        let data = b"the quick brown fox jumps over the lazy dog";
        let mut vc = VcCodec::new();
        vc.add_to_table(data);
        vc.create_codec(true).unwrap();
        let (enc, bits) = vc.encode(data);
        let dec = vc.decode(&enc, bits);
        assert_eq!(dec, data);
    }

    #[test]
    fn huffman_serialize_round_trip() {
        let data = b"aaaaaaaaaaaaaaaaaaabbbbbccccd";
        let mut vc = VcCodec::new();
        vc.add_to_table(data);
        vc.create_codec(true).unwrap();
        let blob = vc.serialize();
        let vc2 = VcCodec::deserialize(&blob).unwrap();
        let (enc, bits) = vc.encode(data);
        let dec = vc2.decode(&enc, bits);
        assert_eq!(dec, data);
    }

    #[test]
    fn huffman_escape_handles_unseen_bytes() {
        // Train on a restricted alphabet, then encode bytes outside it (escaped).
        let train = b"aaaaaaaaaaaaaaaaaaabbbbbbbbbbbbbbbbbbbb";
        let data = b"aaaaXbbbbbbYccccc";
        let mut vc = VcCodec::new();
        vc.add_to_table(train);
        vc.create_codec(true).unwrap();
        let (enc, bits) = vc.encode(data);
        let dec = vc.decode(&enc, bits);
        assert_eq!(dec, data);
    }
}

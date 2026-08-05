//! LTF integer encoding (FastGA `intPut` / `intGet` / `ltfWrite` / `ltfRead`).
//!
//! A variable-length, sign-aware integer codec derived from htslib's CRAM
//! implementation. The first byte's high bits select the encoding:
//!
//! | first byte  | meaning                          | total bytes |
//! |-------------|----------------------------------|-------------|
//! | `0x40..0x7f`| 6-bit positive                    | 1           |
//! | `0xc0..0xff`| 8-bit negative (sign extended)   | 1           |
//! | `0x20..0x3f`| 13-bit positive                   | 2           |
//! | `0x00..0x07`| positive, low 3 bits = data bytes | 2..9        |
//! | `0x80..0x87`| negative, low 3 bits = data bytes | 2..9        |
//!
//! Multi-byte values are stored little-endian; negative ones sign-extend.

use anyhow::{bail, Result};

/// Number of bytes needed to encode `val` in LTF form.
#[allow(clippy::unreadable_literal, clippy::unusual_byte_groupings)]
pub fn int_size(val: i64) -> usize {
    if val >= 0 {
        if val & !0x3f == 0 {
            1
        } else if val & !0x1fff == 0 {
            2
        } else if val & !0xffff == 0 {
            3
        } else if val & !0xff_ffff == 0 {
            4
        } else if val & !0xff_ffff_ff == 0 {
            5
        } else if val & !0xff_ffff_ffff == 0 {
            6
        } else if val & !0xff_ff_ffff_ffff == 0 {
            7
        } else if val & !0xff_ff_ff_ffff_ffff == 0 {
            8
        } else {
            9
        }
    } else if !val & !0x3f == 0 {
        1
    } else if !val & !0xffff == 0 {
        3
    } else if !val & !0xff_ffff == 0 {
        4
    } else if !val & !0xff_ffff_ff == 0 {
        5
    } else if !val & !0xff_ffff_ffff == 0 {
        6
    } else if !val & !0xff_ff_ffff_ffff == 0 {
        7
    } else if !val & !0xff_ff_ff_ffff_ffff == 0 {
        8
    } else {
        9
    }
}

/// Append the LTF encoding of `val` to `out`, returning bytes written.
#[allow(clippy::unreadable_literal, clippy::unusual_byte_groupings)]
pub fn int_put(out: &mut Vec<u8>, val: i64) -> usize {
    if val >= 0 {
        if val & !0x3f == 0 {
            out.push((val as u8) | 0x40);
            1
        } else if val & !0x1fff == 0 {
            out.push(((val >> 8) as u8) | 0x20);
            out.push((val & 0xff) as u8);
            2
        } else if val & !0xffff == 0 {
            out.push(1);
            push_le(out, val, 2);
            3
        } else if val & !0xff_ffff == 0 {
            out.push(2);
            push_le(out, val, 3);
            4
        } else if val & !0xff_ffff_ff == 0 {
            out.push(3);
            push_le(out, val, 4);
            5
        } else if val & !0xff_ffff_ffff == 0 {
            out.push(4);
            push_le(out, val, 5);
            6
        } else if val & !0xff_ff_ffff_ffff == 0 {
            out.push(5);
            push_le(out, val, 6);
            7
        } else if val & !0xff_ff_ff_ffff_ffff == 0 {
            out.push(6);
            push_le(out, val, 7);
            8
        } else {
            out.push(7);
            push_le(out, val, 8);
            9
        }
    } else if !val & !0x3f == 0 {
        out.push((val as u8) | 0x40);
        1
    } else if !val & !0xffff == 0 {
        out.push(0x81);
        push_le(out, val, 2);
        3
    } else if !val & !0xff_ffff == 0 {
        out.push(0x82);
        push_le(out, val, 3);
        4
    } else if !val & !0xff_ffff_ff == 0 {
        out.push(0x83);
        push_le(out, val, 4);
        5
    } else if !val & !0xff_ffff_ffff == 0 {
        out.push(0x84);
        push_le(out, val, 5);
        6
    } else if !val & !0xff_ff_ffff_ffff == 0 {
        out.push(0x85);
        push_le(out, val, 6);
        7
    } else if !val & !0xff_ff_ff_ffff_ffff == 0 {
        out.push(0x86);
        push_le(out, val, 7);
        8
    } else {
        out.push(0x87);
        push_le(out, val, 8);
        9
    }
}

fn push_le(out: &mut Vec<u8>, val: i64, nbytes: usize) {
    let le = val.to_le_bytes();
    out.extend_from_slice(&le[..nbytes]);
}

/// Decode an LTF integer starting at `buf[pos]`, returning `(value, pos_after)`.
pub fn int_get(buf: &[u8], pos: usize) -> Result<(i64, usize)> {
    let u0 = buf[pos];
    match u0 >> 5 {
        2 | 3 => Ok(((u0 & 0x3f) as i64, pos + 1)),
        6 | 7 => Ok((u0 as i8 as i64, pos + 1)),
        1 => {
            let u1 = buf.get(pos + 1).copied().unwrap_or(0);
            Ok(((((u0 & 0x1f) as i64) << 8) | u1 as i64, pos + 2))
        }
        0 => {
            let h = (u0 & 0x07) as usize;
            if h == 0 {
                bail!("int packing error");
            }
            let nbytes = h + 1;
            let mut val: i64 = 0;
            for k in 0..nbytes {
                let b = buf.get(pos + 1 + k).copied().unwrap_or(0) as i64;
                val |= b << (8 * k);
            }
            Ok((val, pos + 1 + nbytes))
        }
        4 => {
            let h = (u0 & 0x07) as usize;
            if h == 0 {
                bail!("int packing error");
            }
            let nbytes = h + 1;
            let mut val: u64 = 0;
            for k in 0..nbytes {
                let b = buf.get(pos + 1 + k).copied().unwrap_or(0) as u64;
                val |= b << (8 * k);
            }
            // Negative multi-byte: all bits above the data field are forced to 1
            // (FastGA ORs the field with `0xffff...0000`), regardless of the
            // stored sign bit.
            let mask = if nbytes >= 8 {
                0
            } else {
                !0u64 << (8 * nbytes)
            };
            Ok(((val | mask) as i64, pos + 1 + nbytes))
        }
        _ => bail!("int packing error"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_various_values() {
        let vals = [
            0i64,
            1,
            63,
            64,
            8191,
            8192,
            65535,
            65536,
            1_000_000,
            46_416_520,
            -1,
            -63,
            -64,
            -65,
            -8192,
            -32768,
            -32769,
            i64::MAX,
            i64::MIN,
        ];
        for v in vals {
            let mut buf = Vec::new();
            int_put(&mut buf, v);
            let (got, _) = int_get(&buf, 0).unwrap();
            assert_eq!(got, v, "round trip failed for {v}");
        }
    }

    #[test]
    fn matches_fastga_byte_patterns() {
        // Single byte positive: 0x40 | val
        let mut buf = Vec::new();
        int_put(&mut buf, 5);
        assert_eq!(buf, vec![0x45]);
        // Single byte negative: val | 0x40 (sign extended low byte)
        let mut buf = Vec::new();
        int_put(&mut buf, -5);
        assert_eq!(buf, vec![0xfb]);
        // Two byte positive: (val>>8)|0x20, val&0xff
        let mut buf = Vec::new();
        int_put(&mut buf, 300);
        assert_eq!(buf, vec![0x21, 0x2c]);
    }
}

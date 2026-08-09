use criterion::{black_box, criterion_group, criterion_main, Criterion};
use flate2::{Compression, Decompress, FlushDecompress};
use std::io::Write;

// zlib-ng's C API requires a zero-initialized z_stream (zalloc/zfree == NULL
// selects the default allocator); the struct contains `extern "C" fn` fields,
// which trip the invalid_value lint even though null fn pointers are valid.
#[allow(invalid_value)]
fn zeroed_z_stream() -> libz_ng_sys::z_stream {
    unsafe { std::mem::MaybeUninit::<libz_ng_sys::z_stream>::zeroed().assume_init() }
}

/// One 64 KiB raw DEFLATE stream (a BGZF block payload), DNA-like data.
fn make_cdata() -> Vec<u8> {
    let bases = b"ACGT";
    let mut x = 0x1234_5678_9abc_def0u64;
    let mut raw = Vec::with_capacity(1 << 16);
    for _ in 0..(1 << 16) {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        raw.push(bases[(x >> 33) as usize & 3]);
    }
    let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), Compression::new(6));
    enc.write_all(&raw).unwrap();
    enc.finish().unwrap()
}

fn bench_inflate_block(c: &mut Criterion) {
    let cdata = make_cdata();
    assert!(cdata.len() < (1 << 16), "compressed block too large");
    let mut group = c.benchmark_group("bgzf_inflate/block");

    // flate2 (zlib-rs), inflater reused across blocks: CachedBgzfReader default.
    group.bench_function("flate2_zlibrs_reuse", |b| {
        let mut inf = Decompress::new(false);
        let mut out = vec![0u8; 1 << 16];
        b.iter(|| {
            inf.reset(false);
            let st = inf
                .decompress(&cdata, &mut out, FlushDecompress::Finish)
                .unwrap();
            assert_eq!(st, flate2::Status::StreamEnd);
            black_box(&out);
        });
    });

    // flate2 (zlib-rs), fresh inflater per block: noodles-bgzf behavior.
    group.bench_function("flate2_zlibrs_fresh", |b| {
        let mut out = vec![0u8; 1 << 16];
        b.iter(|| {
            let mut inf = Decompress::new(false);
            let st = inf
                .decompress(&cdata, &mut out, FlushDecompress::Finish)
                .unwrap();
            assert_eq!(st, flate2::Status::StreamEnd);
            black_box(&out);
        });
    });

    // linflate: full-buffer, SIMD match copy.
    group.bench_function("linflate", |b| {
        let mut out = vec![0u8; (1 << 16) + linflate::OVERWRITE_HEADROOM];
        b.iter(|| {
            let n = linflate::inflate_into(&cdata, &mut out).unwrap();
            black_box(&out[..n]);
        });
    });

    // libdeflater (libdeflate C), decompressor reused.
    group.bench_function("libdeflater", |b| {
        let mut dec = libdeflater::Decompressor::new();
        let mut out = vec![0u8; 1 << 16];
        b.iter(|| {
            let n = dec.deflate_decompress(&cdata, &mut out).unwrap();
            black_box(&out[..n]);
        });
    });

    // miniz_oxide: flate2 default backend, baseline for comparison.
    group.bench_function("miniz_oxide", |b| {
        b.iter(|| {
            let out = miniz_oxide::inflate::decompress_to_vec(&cdata).unwrap();
            black_box(out);
        });
    });

    // isal (Intel ISA-L) stateless inflate.
    group.bench_function("isal_stateless", |b| {
        use isal_sys::igzip_lib::{
            inflate_state, isal_inflate_init, isal_inflate_stateless, ISAL_DECOMP_OK, ISAL_DEFLATE,
        };
        let mut state = unsafe { std::mem::zeroed::<inflate_state>() };
        unsafe { isal_inflate_init(&mut state) };
        let mut out = vec![0u8; 1 << 16];
        b.iter(|| {
            unsafe {
                state.next_in = cdata.as_ptr().cast_mut();
                state.avail_in = cdata.len() as u32;
                state.next_out = out.as_mut_ptr();
                state.avail_out = out.len() as u32;
                state.crc_flag = ISAL_DEFLATE;
                let ret = isal_inflate_stateless(&mut state);
                assert_eq!(ret, ISAL_DECOMP_OK as i32);
            }
            black_box(&out);
        });
    });

    // libz-ng (zlib-ng C), stream reused via inflateReset.
    group.bench_function("libz_ng", |b| {
        let mut strm = zeroed_z_stream();
        let init_ret = unsafe {
            libz_ng_sys::inflateInit2_(
                &mut strm,
                -15,
                libz_ng_sys::zlibVersion(),
                std::mem::size_of::<libz_ng_sys::z_stream>() as i32,
            )
        };
        assert_eq!(init_ret, libz_ng_sys::Z_OK);
        let mut out = vec![0u8; 1 << 16];
        b.iter(|| {
            unsafe {
                strm.next_in = cdata.as_ptr().cast_mut();
                strm.avail_in = cdata.len() as u32;
                strm.next_out = out.as_mut_ptr();
                strm.avail_out = out.len() as u32;
                let ret = libz_ng_sys::inflate(&mut strm, libz_ng_sys::Z_FINISH);
                assert_eq!(ret, libz_ng_sys::Z_STREAM_END);
                libz_ng_sys::inflateReset(&mut strm);
            }
            black_box(&out);
        });
        unsafe {
            libz_ng_sys::inflateEnd(&mut strm);
        }
    });

    group.finish();
}

criterion_group!(benches, bench_inflate_block);
criterion_main!(benches);

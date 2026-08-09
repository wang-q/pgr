use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use flate2::{Decompress, FlushDecompress};
use isal_sys::igzip_lib::{
    inflate_state, isal_inflate_init, isal_inflate_stateless, ISAL_DECOMP_OK, ISAL_DEFLATE,
};
use noodles_bgzf as bgzf;
use pgr::libs::bgzf::{BlockInflater, CachedBgzfReader};
use pgr::libs::fmt::fa::build_gzi_index;
use std::io;
use std::io::{Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::sync::OnceLock;

/// Deterministic BGZF test set: ~50 MB FASTA compressed once, with a `.gzi`
/// index. Reused across benchmark runs via a stable temp path.
struct BgzfTestSet {
    path: std::path::PathBuf,
    /// One read per block (cold: every read misses the previous block).
    cold_reads: Vec<(u64, usize)>,
    /// 1000 reads inside a single block (warm: block reuse).
    warm_reads: Vec<(u64, usize)>,
    n_blocks: usize,
}

fn test_set() -> &'static BgzfTestSet {
    static SET: OnceLock<BgzfTestSet> = OnceLock::new();
    SET.get_or_init(build_test_set)
}

fn build_test_set() -> BgzfTestSet {
    let path = std::env::temp_dir().join("pgr_bgzf_bench_50m.fa.gz");
    if !path.is_file() {
        let fa = generate_fasta(50, 1_000_000);
        let mut writer =
            bgzf::io::Writer::new(std::fs::File::create(&path).expect("create bgzf bench file"));
        writer.write_all(&fa).expect("write bgzf bench data");
        writer.finish().expect("finish bgzf writer");
        build_gzi_index(path.to_str().unwrap()).expect("build gzi index");
    }
    let gzi_path = format!("{}.gzi", path.to_str().unwrap());
    assert!(std::path::Path::new(&gzi_path).is_file(), "missing .gzi");

    let n_blocks = 50 * 1_000_000 / (1 << 16) + 2;
    let mut x = 0x9e37_79b9_7f4a_7c15u64;
    let mut rng = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x >> 33
    };

    // Cold: one read per block, offset within the first 60 KB of the block.
    let n_cold = n_blocks.min(800);
    let cold_reads = (0..n_cold)
        .map(|b| {
            let off = b as u64 * (1 << 16) + rng() % 60_000;
            (off, 100)
        })
        .collect();

    // Warm: 1000 reads inside the middle block.
    let mid = (n_blocks / 2) as u64 * (1 << 16);
    let warm_reads = (0..1000).map(|_| (mid + rng() % 65_000, 100)).collect();

    BgzfTestSet {
        path,
        cold_reads,
        warm_reads,
        n_blocks,
    }
}

fn generate_fasta(n_seq: usize, seq_len: usize) -> Vec<u8> {
    let bases = b"ACGT";
    let mut x = 0x1234_5678_9abc_def0u64;
    let mut fa = Vec::with_capacity(n_seq * (seq_len + seq_len / 80 + 20));
    for s in 0..n_seq {
        fa.extend_from_slice(format!(">seq{s}\n").as_bytes());
        let mut col = 0usize;
        for _ in 0..seq_len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            fa.push(bases[(x >> 33) as usize & 3]);
            col += 1;
            if col == 80 {
                fa.push(b'\n');
                col = 0;
            }
        }
        if col != 0 {
            fa.push(b'\n');
        }
    }
    fa
}

fn replay_reads<R: Read + Seek>(reader: &mut R, reads: &[(u64, usize)]) {
    let mut buf = vec![0u8; 100];
    for &(off, len) in reads {
        reader.seek(SeekFrom::Start(off)).expect("seek");
        buf.resize(len, 0);
        reader.read_exact(&mut buf).expect("read");
        black_box(&buf);
    }
}

struct Flate2Inflater {
    decompress: Decompress,
}

impl Default for Flate2Inflater {
    fn default() -> Self {
        Self {
            decompress: Decompress::new(false),
        }
    }
}

impl BlockInflater for Flate2Inflater {
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()> {
        self.decompress.reset(false);
        let status = self
            .decompress
            .decompress(cdata, dst, FlushDecompress::Finish)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if status != flate2::Status::StreamEnd {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "incomplete inflate",
            ));
        }
        Ok(())
    }
}

fn bench_indexed_reader(c: &mut Criterion) {
    let set = test_set();
    let mut group = c.benchmark_group("bgzf/indexed_reader");

    group.bench_function(BenchmarkId::new("cold", set.n_blocks), |b| {
        b.iter(|| {
            let mut rdr = bgzf::io::indexed_reader::Builder::default()
                .build_from_path(&set.path)
                .unwrap();
            replay_reads(&mut rdr, &set.cold_reads);
        })
    });

    group.bench_function(BenchmarkId::new("warm", 1000), |b| {
        b.iter(|| {
            let mut rdr = bgzf::io::indexed_reader::Builder::default()
                .build_from_path(&set.path)
                .unwrap();
            replay_reads(&mut rdr, &set.warm_reads);
        })
    });

    group.finish();
}

struct LinflateInflater {
    buf: Vec<u8>,
}

impl BlockInflater for LinflateInflater {
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()> {
        self.buf.resize(dst.len() + linflate::OVERWRITE_HEADROOM, 0);
        let n = linflate::inflate_into(cdata, &mut self.buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        dst[..n].copy_from_slice(&self.buf[..n]);
        Ok(())
    }
}

struct IsalInflater {
    state: inflate_state,
}

impl BlockInflater for IsalInflater {
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()> {
        unsafe {
            self.state.next_in = cdata.as_ptr().cast_mut();
            self.state.avail_in = cdata.len() as u32;
            self.state.next_out = dst.as_mut_ptr();
            self.state.avail_out = dst.len() as u32;
            self.state.crc_flag = ISAL_DEFLATE;
            let ret = isal_inflate_stateless(&mut self.state);
            if ret != ISAL_DECOMP_OK as i32 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "isal inflate failed",
                ));
            }
        }
        Ok(())
    }
}

struct LibzNgInflater {
    strm: libz_ng_sys::z_stream,
}

impl BlockInflater for LibzNgInflater {
    fn inflate(&mut self, cdata: &[u8], dst: &mut [u8]) -> io::Result<()> {
        unsafe {
            let init_ret = libz_ng_sys::inflateInit2_(
                &mut self.strm,
                -15,
                libz_ng_sys::zlibVersion(),
                std::mem::size_of::<libz_ng_sys::z_stream>() as i32,
            );
            if init_ret != libz_ng_sys::Z_OK {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "zlib-ng init failed",
                ));
            }
            self.strm.next_in = cdata.as_ptr().cast_mut();
            self.strm.avail_in = cdata.len() as u32;
            self.strm.next_out = dst.as_mut_ptr();
            self.strm.avail_out = dst.len() as u32;
            let ret = libz_ng_sys::inflate(&mut self.strm, libz_ng_sys::Z_FINISH);
            if ret != libz_ng_sys::Z_STREAM_END {
                eprintln!(
                    "zlib-ng ret={ret} avail_in={} avail_out={} total_in={} total_out={} cdata_len={} dst_len={}",
                    self.strm.avail_in, self.strm.avail_out, self.strm.total_in, self.strm.total_out,
                    cdata.len(), dst.len()
                );
                libz_ng_sys::inflateEnd(&mut self.strm);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "zlib-ng inflate failed",
                ));
            }
            libz_ng_sys::inflateEnd(&mut self.strm);
        }
        Ok(())
    }
}

type InflaterFactory = fn() -> Box<dyn BlockInflater>;

// zlib-ng's C API requires a zero-initialized z_stream (zalloc/zfree == NULL
// selects the default allocator); the struct contains `extern "C" fn` fields,
// which trip the invalid_value lint even though null fn pointers are valid.
#[allow(invalid_value)]
fn zeroed_z_stream() -> libz_ng_sys::z_stream {
    unsafe { std::mem::MaybeUninit::<libz_ng_sys::z_stream>::zeroed().assume_init() }
}

fn inflater_variants() -> Vec<(&'static str, InflaterFactory)> {
    vec![
        ("flate2", || Box::<Flate2Inflater>::default()),
        ("linflate", || {
            Box::new(LinflateInflater { buf: Vec::new() })
        }),
        ("libdeflater", || {
            Box::<pgr::libs::bgzf::LibdeflaterInflater>::default()
        }),
        ("isal", || {
            let mut state = unsafe { std::mem::zeroed::<inflate_state>() };
            unsafe { isal_inflate_init(&mut state) };
            Box::new(IsalInflater { state })
        }),
        ("libz_ng", || {
            let mut strm = zeroed_z_stream();
            let ret = unsafe {
                libz_ng_sys::inflateInit2_(
                    &mut strm,
                    -15,
                    libz_ng_sys::zlibVersion(),
                    std::mem::size_of::<libz_ng_sys::z_stream>() as i32,
                )
            };
            assert_eq!(ret, libz_ng_sys::Z_OK);
            Box::new(LibzNgInflater { strm })
        }),
    ]
}

fn bench_cached_reader(c: &mut Criterion) {
    let set = test_set();
    let nz = NonZeroUsize::new(16).expect("non-zero");
    for (name, inflater) in inflater_variants() {
        let mut group = c.benchmark_group(format!("bgzf/cached_reader/{name}"));

        group.bench_function(BenchmarkId::new("cold", set.n_blocks), |b| {
            b.iter(|| {
                let mut rdr =
                    CachedBgzfReader::open_with_inflater(&set.path, nz, inflater()).unwrap();
                replay_reads(&mut rdr, &set.cold_reads);
            })
        });

        group.bench_function(BenchmarkId::new("warm", 1000), |b| {
            b.iter(|| {
                let mut rdr =
                    CachedBgzfReader::open_with_inflater(&set.path, nz, inflater()).unwrap();
                replay_reads(&mut rdr, &set.warm_reads);
            })
        });

        group.finish();
    }
}

fn bench_sequential_read(c: &mut Criterion) {
    let set = test_set();
    let mut group = c.benchmark_group("bgzf/sequential_read");

    group.bench_function("gz_reader", |b| {
        b.iter(|| {
            let mut rdr =
                pgr::libs::bgzf::GzReader::new(std::fs::File::open(&set.path).unwrap()).unwrap();
            let mut total = 0usize;
            let mut buf = vec![0u8; 1 << 16];
            loop {
                let n = rdr.read(&mut buf).unwrap();
                if n == 0 {
                    break;
                }
                total += n;
            }
            black_box(total);
        })
    });

    for workers in [1usize, 2, 4] {
        group.bench_function(BenchmarkId::new("parallel", workers), |b| {
            b.iter(|| {
                let mut rdr =
                    pgr::libs::bgzf::ParallelBgzfReader::open(&set.path, workers).unwrap();
                let mut total = 0usize;
                let mut buf = vec![0u8; 1 << 16];
                loop {
                    let n = rdr.read(&mut buf).unwrap();
                    if n == 0 {
                        break;
                    }
                    total += n;
                }
                black_box(total);
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_indexed_reader,
    bench_cached_reader,
    bench_sequential_read
);
criterion_main!(benches);
